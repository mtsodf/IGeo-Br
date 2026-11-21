"""Verificação: injeção em fratura de comprimento fixo (solução de Sneddon).

Caso: ``tests/sneddon_injection/sneddon_injection.geo``.

Uma fratura de meio comprimento L, em deformação plana, numa rocha elástica
praticamente impermeável, recebe fluido a pressão controlada no poço. Com a
fratura muito condutiva, a pressão dentro dela é uniforme, e a pressão
líquida é Dp = p_poço - Shmin. A abertura é dada pela solução de Sneddon para
uma fratura pressurizada (Sneddon, 1946; Sneddon & Lowengrub, 1969):

    w(x) = 4 (1 - nu^2) Dp / E * sqrt(L^2 - x^2)          (abertura total)

O iGeo modela 1/4 do domínio (simetria em x = 0 e y = 0) e calcula o
deslocamento da face y = 0, uy = w/2. O volume de fluido na fratura, no
quarto de modelo e com altura H, é

    V = H * Integral_0^L uy dx = pi (1 - nu^2) Dp L^2 H / (2 E).

São feitas duas simulações:
  * o caso do .geo (dx = 0,25 m, domínio 20L), comparado com Sneddon;
  * uma variante com malha grossa (dx = 0,5 m) e a pressão mantida por
    3 passos, que verifica a convergência com o refino da malha e que a
    abertura não muda enquanto a pressão é mantida.

O erro numérico vem dos elementos lineares e do contorno fixo a 10L do
poço, que deixam o modelo um pouco mais rígido que o meio infinito; ele
diminui com o refino da malha e com o tamanho do domínio.

Referências
  I. N. Sneddon (1946). The distribution of stress in the neighbourhood of
    a crack in an elastic solid. Proc. R. Soc. Lond. A 187, 229-260.
  I. N. Sneddon, M. Lowengrub (1969). Crack Problems in the Classical
    Theory of Elasticity. Wiley.

Executar (a partir da raiz do repositório):
    make -C iGeo_vrs2.0 test
    # ou, com o executável já compilado:
    python3 -m unittest discover -s tests -v

O executável usado é ``iGeo_vrs2.0/build/iGeo``; outro pode ser indicado na
variável de ambiente IGEO_BIN. Requer numpy.
"""

from __future__ import annotations

import functools
import os
import re
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np

RAIZ = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(RAIZ / "tools"))
from perfil_abertura_fratura import ler_vtk  # noqa: E402

CASO = RAIZ / "tests" / "sneddon_injection" / "sneddon_injection.geo"
POCO = "SneddonInj"
IGEO = Path(os.environ.get("IGEO_BIN", RAIZ / "iGeo_vrs2.0" / "build" / "iGeo"))
TEMPO_LIMITE = 600  # s

# Parâmetros do caso (iguais aos do .geo)
E = 1.0e7          # kPa
NU = 0.25
L = 5.0            # m (meio comprimento da fratura)
H = 1.0            # m (altura da fratura = espessura da camada)
SHMIN = 10000.0    # kPa
P_POCO = 11000.0   # kPa
DP = P_POCO - SHMIN
KGF_CM2_PARA_KPA = 98.0665

PASTAS_TEMPORARIAS: list[Path] = []


def uy_sneddon(x: np.ndarray, dp: float = DP) -> np.ndarray:
    """Deslocamento da face (meia abertura) da solução de Sneddon."""
    return 2.0 * (1.0 - NU**2) * dp / E * np.sqrt(np.clip(L**2 - x**2, 0.0, None))


def volume_sneddon(dp: float = DP) -> float:
    """Volume de fluido na fratura no quarto de modelo (m3)."""
    return np.pi * (1.0 - NU**2) * dp * L**2 * H / (2.0 * E)


def variante_do_caso(nx: int | None = None, passos: int = 1) -> str:
    """Texto do .geo com outro NX (mesmo domínio) e a pressão repetida em `passos` passos."""
    texto = CASO.read_text(encoding="utf-8")
    if nx is not None:
        texto, n = re.subn(r"\n\d+(  \d+  1  [\d.]+ +<- NX)", rf"\n{nx}\1", texto)
        assert n == 1, "linha NX não encontrada no .geo"
    tabela = "".join(f"{k}.0    {P_POCO:.1f}  -1000.0\n" for k in range(1, passos + 1))
    texto, n = re.subn(r"(#WELLPRESSURECONSTRAINT\n[^\n]*\n[^\n]*\n)(?:[^\n]*\n)*?<\n",
                       lambda m: m.group(1) + tabela + "<\n", texto)
    assert n == 1, "tabela de pressão não encontrada no .geo"
    return texto


def face_da_fratura(vtk: Path) -> dict[str, np.ndarray]:
    """Nós da face y = 0 no topo da camada, ordenados por x."""
    malha = ler_vtk(vtk)
    x, y, z = malha.pontos.T
    nos = np.flatnonzero((np.abs(y) < 1e-6) & (np.abs(z - z.max()) < 1e-6))
    nos = nos[np.argsort(x[nos])]
    return {
        "x": x[nos],
        "uy": malha.dados_ponto["displacement"][nos, 1],
        "p": malha.dados_ponto["Pp_kgf_cm2"][nos] * KGF_CM2_PARA_KPA,
    }


@functools.lru_cache(maxsize=None)
def simular(nx: int | None = None, passos: int = 1) -> tuple[dict[str, np.ndarray], ...]:
    """Roda o iGeo (uma vez por variante) e devolve a face da fratura em cada passo."""
    pasta = Path(tempfile.mkdtemp(prefix="igeo_sneddon_"))
    PASTAS_TEMPORARIAS.append(pasta)
    (pasta / CASO.name).write_text(variante_do_caso(nx, passos), encoding="utf-8")
    r = subprocess.run([str(IGEO), str(pasta), CASO.name], cwd=pasta,
                       capture_output=True, text=True, timeout=TEMPO_LIMITE)
    saida = r.stdout + r.stderr
    if r.returncode != 0 or "WORKFLOW FINISHED" not in saida:
        raise AssertionError(f"iGeo falhou (código {r.returncode}):\n{saida[-3000:]}")
    # <poço>_1.vtk é o estado inicial; <poço>_{k+1}.vtk é o fim do passo k.
    return tuple(face_da_fratura(pasta / f"{POCO}_{k + 1}.vtk") for k in range(1, passos + 1))


def erros(face: dict[str, np.ndarray]) -> dict[str, float]:
    """Erros relativos em relação a Sneddon: centro, perfil (L2) e volume."""
    x, uy = face["x"], face["uy"]
    trecho = x <= L + 1e-9
    ref = uy_sneddon(x[trecho])
    return {
        "centro": uy[0] / ref[0] - 1.0,
        "perfil_l2": float(np.sqrt(np.trapezoid((uy[trecho] - ref) ** 2, x[trecho])
                                   / np.trapezoid(ref**2, x[trecho]))),
        "volume": H * np.trapezoid(np.where(trecho, uy, 0.0), x) / volume_sneddon() - 1.0,
    }


def tearDownModule() -> None:
    for pasta in PASTAS_TEMPORARIAS:
        shutil.rmtree(pasta, ignore_errors=True)


@unittest.skipUnless(IGEO.is_file(), f"executável não encontrado: {IGEO} "
                     "(compile com 'make -C iGeo_vrs2.0' ou defina IGEO_BIN)")
class TesteInjecaoSneddon(unittest.TestCase):
    """O caso do .geo comparado com a solução de Sneddon."""

    @classmethod
    def setUpClass(cls) -> None:
        cls.face = simular()[0]
        cls.erros = erros(cls.face)
        cls.dentro = cls.face["x"] < L - 1e-9      # exclui o nó da ponta

    def test_pressao_uniforme_igual_a_do_poco(self):
        """No fim do 1º passo a fratura ainda termina de encher (o fluido que
        a abre vem do poço pela própria fratura): perto da ponta a pressão
        está ~0,5% abaixo da do poço. Tolerância de 1% de Dp."""
        p = self.face["p"][self.dentro]
        self.assertLess(np.ptp(p), 0.01 * DP, "pressão não uniforme na fratura")
        self.assertAlmostEqual(p.mean() - SHMIN, DP, delta=0.01 * DP)

    def test_abertura_no_centro(self):
        self.assertLess(abs(self.erros["centro"]), 0.02,
                        f"uy(0) difere de Sneddon em {self.erros['centro']:+.2%}")

    def test_perfil_da_abertura(self):
        self.assertLess(self.erros["perfil_l2"], 0.03,
                        f"erro L2 relativo do perfil = {self.erros['perfil_l2']:.2%}")

    def test_volume_da_fratura(self):
        self.assertLess(abs(self.erros["volume"]), 0.03,
                        f"volume difere de Sneddon em {self.erros['volume']:+.2%}")

    def test_face_fechada_fora_da_fratura(self):
        fora = self.face["x"] >= L - 1e-9          # inclui o nó da ponta
        self.assertLess(np.abs(self.face["uy"][fora]).max(), 1e-12)


@unittest.skipUnless(IGEO.is_file(), f"executável não encontrado: {IGEO}")
class TesteConvergenciaEPressaoMantida(unittest.TestCase):
    """Malha grossa (dx = 0,5 m) com a pressão mantida por 3 passos."""

    NX_GROSSO = 200      # dx = 0,5 m no domínio de 100 m (o caso usa NX = 400)
    PASSOS = 3

    @classmethod
    def setUpClass(cls) -> None:
        cls.grosso = simular(cls.NX_GROSSO, cls.PASSOS)
        cls.fino = simular()[0]

    def test_erro_diminui_com_o_refino(self):
        e_grosso, e_fino = erros(self.grosso[0]), erros(self.fino)
        for nome in ("centro", "perfil_l2", "volume"):
            with self.subTest(erro=nome):
                self.assertLess(abs(e_fino[nome]), abs(e_grosso[nome]),
                                f"{nome}: dx=0,5 m {e_grosso[nome]:+.2%}, "
                                f"dx=0,25 m {e_fino[nome]:+.2%}")

    def test_abertura_constante_com_pressao_mantida(self):
        """A abertura não pode mudar enquanto a pressão é mantida.

        Antes da correção do carregamento da face (UpdateHydraulicFracture),
        o resíduo (p + Syy médio do elemento) era reaplicado a cada passo e,
        no elemento da ponta, nunca se anulava: a abertura crescia passo a
        passo (+2,6% no 2º passo, +10% no 6º). Hoje o resíduo usa a força
        interna nodal K*u, que em equilíbrio é exatamente a carga aplicada.
        """
        uy0 = np.array([face["uy"][0] for face in self.grosso])
        variacao = np.abs(uy0 / uy0[0] - 1.0).max()
        self.assertLess(variacao, 1e-3, f"uy(0) nos passos: {uy0}")


if __name__ == "__main__":
    unittest.main(verbosity=2)

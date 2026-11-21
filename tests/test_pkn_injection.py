"""Verificação: fratura PKN (altura fixa) propagada por injeção a vazão constante.

Caso: ``tests/pkn_injection/pkn_injection.geo``.

Uma fratura vertical de altura H fixa, numa rocha elástica homogênea e
impermeável, recebe um fluido newtoniano a vazão constante. A solução de
referência é a de Nordgren (1972) para o modelo PKN (Perkins & Kern, 1961),
sem perda de fluido e sem tenacidade. Com q = vazão de uma asa:

    L(t)    = 0,68 [G q^3 / ((1-nu) mu H^4)]^(1/5) t^(4/5)
    w0(t)   = 2,5  [(1-nu) mu q^2 / (G H)]^(1/5) t^(1/5)    (abertura no poço)
    pnet(t) = G w0 / ((1-nu) H)

As hipóteses do PKN que o iGeo deve reproduzir, e que este teste verifica:

  * balanço de massa: o volume da fratura é o volume injetado;
  * a fratura não cresce para fora das camadas com #CDM_FRACTURE;
  * cada seção vertical é elíptica e está em deformação plana:
    w(z) = w_max sqrt(1 - (2z/H)^2) e pnet = G w_max / ((1-nu) H);
  * a vazão ao longo da fratura segue a lei cúbica com a abertura total:
    Q = -(dp/dx) Integral w(z)^3 dz / (12 mu);
  * a abertura diminui do poço para a ponta.

O comprimento, a abertura no poço e a pressão ainda diferem da solução de
Nordgren, porque o tratamento da ponta (critério de propagação e zona seca
à frente do fluido) não reproduz a ponta do PKN. Essa comparação está no
teste ``test_comparacao_com_nordgren``, marcado como falha esperada; quando
a ponta for corrigida, ele passa a ser reportado como sucesso inesperado.
Detalhes em ``docs/relatorio_pkn/README.md``.

Referências
  T. K. Perkins, L. R. Kern (1961). Widths of hydraulic fractures.
    J. Pet. Tech. 13(9), 937-949.
  R. P. Nordgren (1972). Propagation of a vertical hydraulic fracture.
    SPE Journal 12(4), 306-314.

Executar (a partir da raiz do repositório):
    ctest --test-dir build -L pkn --output-on-failure
    # ou, com o executável já compilado:
    python3 -m unittest -v tests/test_pkn_injection.py

O executável usado é ``iGeo_vrs2.0/build/iGeo``; outro pode ser indicado na
variável de ambiente IGEO_BIN. Requer numpy.
"""

from __future__ import annotations

import functools
import os
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

CASO = RAIZ / "tests" / "pkn_injection" / "pkn_injection.geo"
POCO = "PKN"
IGEO = Path(os.environ.get("IGEO_BIN", RAIZ / "iGeo_vrs2.0" / "build" / "iGeo"))
TEMPO_LIMITE = 1800  # s

# Parâmetros do caso (iguais aos do .geo)
E = 1.0e10             # Pa
NU = 0.25
G = E / (2.0 * (1.0 + NU))
H = 10.0               # m (altura da fratura)
Z_CENTRO = -1040.5     # m (centro da zona com fratura)
MU = 0.1               # Pa.s (100 cP)
SHMIN = 1.0e7          # Pa
Q_MODELO = 432.0 / 86400.0   # m3/s no quarto de modelo
Q_ASA = 2.0 * Q_MODELO       # m3/s numa asa (as duas metades em y)
DT = 10.0              # s (duração de cada passo)
PASSOS = 12
KGF_CM2_PARA_PA = 98066.5

PASTAS_TEMPORARIAS: list[Path] = []


def nordgren(t: float) -> dict[str, float]:
    """Comprimento, abertura no poço e pressão líquida no poço (Nordgren, 1972)."""
    L = 0.68 * (G * Q_ASA**3 / ((1 - NU) * MU * H**4)) ** 0.2 * t**0.8
    w0 = 2.5 * ((1 - NU) * MU * Q_ASA**2 / (G * H)) ** 0.2 * t**0.2
    return {"L": L, "w0": w0, "pnet": G * w0 / ((1 - NU) * H)}


def face_da_fratura(vtk: Path) -> dict[str, np.ndarray]:
    """Deslocamento uy e pressão líquida na face y = 0, em grade (z, x)."""
    malha = ler_vtk(vtk)
    x, y, z = malha.pontos.T
    sel = np.abs(y) < 1e-6
    xs = np.unique(np.round(x[sel], 6))
    zs = np.unique(np.round(z[sel], 6))
    ix = np.searchsorted(xs, np.round(x[sel], 6))
    iz = np.searchsorted(zs, np.round(z[sel], 6))
    uy = np.zeros((zs.size, xs.size))
    pnet = np.zeros_like(uy)
    uy[iz, ix] = malha.dados_ponto["displacement"][sel, 1]
    pnet[iz, ix] = malha.dados_ponto["Pp_kgf_cm2"][sel] * KGF_CM2_PARA_PA - SHMIN
    return {"x": xs, "z": zs, "uy": uy, "pnet": pnet}


@functools.lru_cache(maxsize=None)
def simular() -> tuple[dict[str, np.ndarray], ...]:
    """Roda o iGeo uma vez e devolve a face da fratura no fim de cada passo."""
    pasta = Path(tempfile.mkdtemp(prefix="igeo_pkn_"))
    PASTAS_TEMPORARIAS.append(pasta)
    shutil.copy(CASO, pasta / CASO.name)
    r = subprocess.run([str(IGEO), str(pasta), CASO.name], cwd=pasta,
                       capture_output=True, text=True, timeout=TEMPO_LIMITE)
    saida = r.stdout + r.stderr
    if r.returncode != 0 or "WORKFLOW FINISHED" not in saida:
        raise AssertionError(f"iGeo falhou (código {r.returncode}):\n{saida[-3000:]}")
    # <poço>_1.vtk é o estado inicial; <poço>_{k+1}.vtk é o fim do passo k.
    return tuple(face_da_fratura(pasta / f"{POCO}_{k + 1}.vtk") for k in range(1, PASSOS + 1))


def volume_modelo(face: dict[str, np.ndarray], x_max: float = np.inf) -> float:
    """Volume de fluido na fratura no quarto de modelo, entre x = 0 e x_max."""
    m = face["x"] <= x_max + 1e-9
    return float(np.trapezoid(np.trapezoid(face["uy"][:, m], face["x"][m], axis=1), face["z"]))


def linha_central(face: dict[str, np.ndarray]) -> int:
    return int(np.argmin(np.abs(face["z"] - Z_CENTRO)))


def frente_molhada(face: dict[str, np.ndarray]) -> float:
    """Posição da frente do fluido: primeiro nó com pressão líquida < 5% da do poço."""
    pn = face["pnet"][linha_central(face)]
    molhado = np.flatnonzero(pn > 0.05 * pn[0])
    return float(face["x"][molhado.max() + 1])


def tearDownModule() -> None:
    for pasta in PASTAS_TEMPORARIAS:
        shutil.rmtree(pasta, ignore_errors=True)


@unittest.skipUnless(IGEO.is_file(), f"executável não encontrado: {IGEO} "
                     "(compile com 'make -C iGeo_vrs2.0' ou defina IGEO_BIN)")
class TesteFraturaPKN(unittest.TestCase):
    """O caso do .geo comparado com as hipóteses do PKN e com Nordgren."""

    X_SECOES = (4.0, 8.0, 12.0)   # seções verticais longe do poço e da ponta (m)

    @classmethod
    def setUpClass(cls) -> None:
        cls.faces = simular()
        cls.final = cls.faces[-1]

    def test_volume_igual_ao_injetado(self):
        """Balanço de massa da fratura (a partir do 4º passo; o início é o
        rompimento da pré-fratura de um elemento)."""
        for k in range(4, PASSOS + 1):
            with self.subTest(passo=k):
                v, v_inj = volume_modelo(self.faces[k - 1]), Q_MODELO * DT * k
                self.assertLess(abs(v / v_inj - 1.0), 0.03,
                                f"volume {v:.4f} m3, injetado {v_inj:.4f} m3")

    def test_altura_contida(self):
        """Nenhuma abertura fora das camadas com #CDM_FRACTURE nem nas suas
        bordas (z = +-H/2), que são a frente da fratura e ficam presas.
        Antes da correção da liberação de nós (PropagateFracture), os nós da
        borda eram soltos quando um elemento novo abria, e a fratura crescia
        uma camada acima e outra abaixo."""
        f = self.final
        fora = np.abs(f["z"] - Z_CENTRO) >= H / 2 - 1e-6
        self.assertLess(np.abs(f["uy"][fora]).max(), 1e-12,
                        f"abertura máxima fora da zona: {2e3 * np.abs(f['uy'][fora]).max():.3f} mm")

    def test_secao_eliptica(self):
        f = self.final
        dentro = np.abs(f["z"] - Z_CENTRO) <= H / 2 + 1e-6
        z = f["z"][dentro] - Z_CENTRO
        for x in self.X_SECOES:
            with self.subTest(x=x):
                w = 2.0 * f["uy"][dentro, int(np.argmin(np.abs(f["x"] - x)))]
                ref = w.max() * np.sqrt(np.clip(1.0 - (2.0 * z / H) ** 2, 0.0, None))
                l2 = np.sqrt(np.trapezoid((w - ref) ** 2, z) / np.trapezoid(ref**2, z))
                self.assertLess(l2, 0.03, f"erro L2 em relação à elipse = {l2:.2%}")

    def test_rigidez_da_secao(self):
        """pnet / (G w_max / ((1-nu) H)) = 1 em deformação plana. A malha (10
        elementos na altura) deixa a seção ~5% mais rígida (verificação
        estática no relatório), daí o intervalo [1,00; 1,12]."""
        f = self.final
        jz = linha_central(f)
        for x in self.X_SECOES:
            with self.subTest(x=x):
                ix = int(np.argmin(np.abs(f["x"] - x)))
                razao = f["pnet"][jz, ix] / (G * 2.0 * f["uy"][jz, ix] / ((1 - NU) * H))
                self.assertGreater(razao, 1.00, f"razão = {razao:.3f}")
                self.assertLess(razao, 1.12, f"razão = {razao:.3f}")

    def test_lei_cubica_local(self):
        """Coeficiente C = Q 12 mu / (-dp/dx Integral w^3 dz) em três trechos,
        com Q (asa inteira) obtida do balanço de volume entre os dois últimos
        passos. Deve ser ~1 (0,93 com a abertura média de cada elemento);
        antes da correção da lei cúbica era ~0,2 (condutividade 4 vezes menor)."""
        a, b = self.faces[-2], self.final
        jz_dentro = np.abs(b["z"] - Z_CENTRO) <= H / 2 + 1e-6
        for xa, xb in ((4.0, 8.0), (8.0, 12.0), (12.0, 16.0)):
            with self.subTest(trecho=(xa, xb)):
                xc = 0.5 * (xa + xb)
                q = 2.0 * (Q_MODELO - (volume_modelo(b, xc) - volume_modelo(a, xc)) / DT)
                ia = int(np.argmin(np.abs(b["x"] - xa)))
                ib = int(np.argmin(np.abs(b["x"] - xb)))
                p = 0.5 * (a["pnet"] + b["pnet"])[jz_dentro]
                uy = 0.5 * (a["uy"] + b["uy"])[jz_dentro]
                aberto = uy[:, ia] > 0
                dpdx = (p[aberto, ib].mean() - p[aberto, ia].mean()) / (xb - xa)
                w = uy[:, ia] + uy[:, ib]          # média das duas seções, abertura total
                i3 = np.trapezoid(w**3, b["z"][jz_dentro])
                c = q * 12.0 * MU / (i3 * -dpdx)
                self.assertGreater(c, 0.85, f"C = {c:.3f}")
                self.assertLess(c, 1.05, f"C = {c:.3f}")

    def test_abertura_decresce_do_poco_para_a_ponta(self):
        f = self.final
        w = 2.0 * f["uy"][linha_central(f)]
        frente = int(np.searchsorted(f["x"], frente_molhada(f)))
        aumento = np.diff(w[:frente]) / w[:frente - 1]
        self.assertLess(aumento.max(), 0.01,
                        f"a abertura cresce {aumento.max():.1%} entre dois nós")

    @unittest.expectedFailure
    def test_comparacao_com_nordgren(self):
        """Falha esperada: o tratamento da ponta ainda não reproduz o PKN
        (em t = 120 s: frente do fluido ~-30%, w0 ~+25%, pnet ~+30%)."""
        f = self.final
        ref = nordgren(DT * PASSOS)
        jz = linha_central(f)
        obtido = {"L": frente_molhada(f), "w0": 2.0 * f["uy"][jz, 0], "pnet": f["pnet"][jz, 0]}
        for nome in ("L", "w0", "pnet"):
            with self.subTest(grandeza=nome):
                self.assertLess(abs(obtido[nome] / ref[nome] - 1.0), 0.10,
                                f"{nome}: iGeo {obtido[nome]:.4g}, Nordgren {ref[nome]:.4g}")


if __name__ == "__main__":
    unittest.main(verbosity=2)

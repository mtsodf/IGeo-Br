#!/usr/bin/env python3
"""Perfil da abertura da fratura hidráulica a partir das saídas VTK do iGeo.

O iGeo grava, a cada passo do workflow de injeção (#WF_WaterInj), um arquivo
VTK ASCII legado chamado <NomeDoPoço>_<passo>.vtk. A abertura da fratura está
no campo de célula ``FractureAperture_mm`` (em mm).

No modelo near-wellbore a fratura ocupa a primeira fileira de elementos da
malha (itheta = 1, veja HydraulicFractureModel.f90):

* malha BOX: elementos junto ao plano Y = 0; o comprimento da fratura é a
  coordenada X e a altura é Z;
* malha RAD: elementos junto ao plano theta = 0; o comprimento é o raio
  r = sqrt(X² + Y²) e a altura é Z.

O script identifica esse plano, organiza a abertura numa grade
(comprimento x profundidade) e gera uma figura com:

1. o perfil vertical da abertura junto ao poço (abertura x profundidade);
2. o perfil ao longo do comprimento, numa profundidade escolhida;
3. o mapa da abertura no plano da fratura, no último passo.

Os perfis mostram todos os passos lidos. Também imprime um resumo por passo
e, opcionalmente, exporta os dados do plano da fratura em CSV.

Requisitos: numpy e, para a figura, matplotlib (pip install numpy matplotlib).

Exemplos:
    # todos os passos de um caso (procura <poço>_<n>.vtk na pasta)
    python3 tools/perfil_abertura_fratura.py Simulation_cases --poco LinearFlow

    # arquivos explícitos, perfil horizontal em z = -2480 m, com CSV
    python3 tools/perfil_abertura_fratura.py caso/MLS3B_*.vtk --z -2480 \\
        --csv abertura.csv
"""

from __future__ import annotations

import argparse
import csv
import math
import re
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np

CAMPO_ABERTURA = "FractureAperture_mm"
TOL_FILEIRA = 1e-4  # m; veja identificar_plano
PADRAO_PASSO = re.compile(r"^(?P<poco>.+)_(?P<passo>\d+)\.vtk$", re.IGNORECASE)


# ---------------------------------------------------------------------------
# Leitura do VTK legado (ASCII, UNSTRUCTURED_GRID)
# ---------------------------------------------------------------------------
@dataclass
class MalhaVTK:
    pontos: np.ndarray            # (n_pontos, 3)
    celulas: list[np.ndarray]     # índices (base 0) dos nós de cada célula
    dados_celula: dict[str, np.ndarray]
    dados_ponto: dict[str, np.ndarray]


def _para_float(tokens: list[str]) -> np.ndarray:
    """Converte tokens em float; valores inválidos (ex.: '*****' do Fortran) viram NaN."""
    try:
        return np.array(tokens, dtype=float)
    except ValueError:
        saida = np.empty(len(tokens))
        for i, t in enumerate(tokens):
            try:
                saida[i] = float(t)
            except ValueError:
                saida[i] = math.nan
        return saida


def ler_vtk(caminho: Path) -> MalhaVTK:
    with open(caminho, encoding="latin-1") as f:
        cabecalho = [f.readline() for _ in range(3)]
        if "ASCII" not in cabecalho[2].upper():
            raise ValueError(f"{caminho}: apenas VTK ASCII é suportado")
        tokens = f.read().split()

    pontos = None
    celulas: list[np.ndarray] = []
    dados = {"CELL_DATA": {}, "POINT_DATA": {}}
    secao = None
    n_secao = 0
    i = 0
    n = len(tokens)
    while i < n:
        chave = tokens[i].upper()
        if chave == "POINTS":
            npts = int(tokens[i + 1])
            i += 3
            pontos = _para_float(tokens[i:i + 3 * npts]).reshape(npts, 3)
            i += 3 * npts
        elif chave == "CELLS":
            ncel, tamanho = int(tokens[i + 1]), int(tokens[i + 2])
            i += 3
            bloco = np.array(tokens[i:i + tamanho], dtype=np.int64)
            j = 0
            for _ in range(ncel):
                k = bloco[j]
                celulas.append(bloco[j + 1:j + 1 + k])
                j += k + 1
            i += tamanho
        elif chave == "CELL_TYPES":
            i += 2 + int(tokens[i + 1])
        elif chave in ("CELL_DATA", "POINT_DATA"):
            secao, n_secao = chave, int(tokens[i + 1])
            i += 2
        elif chave == "SCALARS":
            nome = tokens[i + 1]
            i += 3
            ncomp = 1
            if i < n and tokens[i].isdigit():
                ncomp = int(tokens[i])
                i += 1
            if i < n and tokens[i].upper() == "LOOKUP_TABLE":
                i += 2
            valores = _para_float(tokens[i:i + n_secao * ncomp])
            dados[secao][nome] = valores.reshape(n_secao, ncomp) if ncomp > 1 else valores
            i += n_secao * ncomp
        elif chave in ("VECTORS", "NORMALS"):
            nome = tokens[i + 1]
            i += 3
            dados[secao][nome] = _para_float(tokens[i:i + 3 * n_secao]).reshape(n_secao, 3)
            i += 3 * n_secao
        else:
            i += 1

    if pontos is None or not celulas:
        raise ValueError(f"{caminho}: malha (POINTS/CELLS) não encontrada")
    return MalhaVTK(pontos, celulas, dados["CELL_DATA"], dados["POINT_DATA"])


# ---------------------------------------------------------------------------
# Plano da fratura
# ---------------------------------------------------------------------------
@dataclass
class PlanoFratura:
    tipo_malha: str            # 'BOX' ou 'RAD'
    indices: np.ndarray        # células do plano da fratura
    s: np.ndarray              # comprimento (m) do centro de cada coluna da grade
    z: np.ndarray              # profundidade (m) do centro de cada linha (decrescente)
    s_bordas: np.ndarray       # (n_colunas, 2): início e fim de cada coluna (m)
    z_bordas: np.ndarray       # (n_linhas, 2): topo e base de cada linha (m)
    i_linha: np.ndarray        # linha da grade de cada célula do plano
    i_coluna: np.ndarray       # coluna da grade de cada célula do plano


def _agrupar(valores: np.ndarray, tol: float) -> tuple[np.ndarray, np.ndarray]:
    """Agrupa valores próximos (dentro de tol); devolve os centros e o índice de cada valor."""
    ordem = np.argsort(valores)
    centros: list[float] = []
    grupo = np.empty(len(valores), dtype=int)
    atual: list[float] = []
    for idx in ordem:
        v = valores[idx]
        if atual and v - atual[-1] > tol:
            centros.append(float(np.mean(atual)))
            atual = []
        atual.append(v)
        grupo[idx] = len(centros)
    centros.append(float(np.mean(atual)))
    return np.array(centros), grupo


def identificar_plano(malha: MalhaVTK, tipo: str = "auto") -> PlanoFratura:
    centroides = np.array([malha.pontos[c].mean(axis=0) for c in malha.celulas])
    cx, cy, cz = centroides.T

    if tipo == "auto":
        # A malha BOX tem um nó no eixo do poço (x = y = 0); a RAD tem o furo do poço.
        no_eixo = np.any(np.hypot(malha.pontos[:, 0], malha.pontos[:, 1]) < 1e-9)
        tipo = "BOX" if no_eixo else "RAD"
    tipo = tipo.upper()

    # Plano da fratura: na malha BOX, o plano Y = Y mínimo; na RAD, o plano
    # radial de menor ângulo, estimado nos nós de maior raio (onde o
    # arredondamento das coordenadas menos afeta o ângulo). As células da
    # fratura são as que têm uma face (4 nós) sobre esse plano. Os VTK têm
    # coordenadas com 5 casas decimais, daí a tolerância de 1e-4 m.
    xn, yn, z_no = malha.pontos.T
    if tipo == "BOX":
        dist_no = np.abs(yn - yn.min())
        s = cx
        s_no = xn
    else:
        r_no = np.hypot(xn, yn)
        externos = r_no >= 0.5 * r_no.max()
        theta0 = np.arctan2(yn[externos], xn[externos]).min()
        dist_no = np.abs(-np.sin(theta0) * xn + np.cos(theta0) * yn)
        s = np.hypot(cx, cy)
        s_no = r_no
    no_no_plano = dist_no <= TOL_FILEIRA
    indices = np.array([k for k, nos in enumerate(malha.celulas)
                        if np.count_nonzero(no_no_plano[nos]) >= 4])
    if indices.size == 0:
        raise ValueError("plano da fratura não encontrado; tente --malha box ou --malha rad")

    s_col, i_col = _agrupar(s[indices], TOL_FILEIRA)
    z_lin, i_lin = _agrupar(cz[indices], TOL_FILEIRA)
    # Linhas de cima para baixo (z decrescente), como numa seção vertical.
    ordem = np.argsort(-z_lin)
    novo = np.empty_like(ordem)
    novo[ordem] = np.arange(len(ordem))
    i_lin = novo[i_lin]

    # Extensão real de cada coluna e linha, a partir das coordenadas dos nós.
    s_bordas = np.tile([np.inf, -np.inf], (len(s_col), 1))
    z_bordas = np.tile([np.inf, -np.inf], (len(z_lin), 1))
    for k, cel in enumerate(indices):
        nos = malha.celulas[cel]
        c, l = i_col[k], i_lin[k]
        s_bordas[c] = [min(s_bordas[c, 0], s_no[nos].min()), max(s_bordas[c, 1], s_no[nos].max())]
        z_bordas[l] = [min(z_bordas[l, 0], z_no[nos].min()), max(z_bordas[l, 1], z_no[nos].max())]
    return PlanoFratura(tipo, indices, s_col, z_lin[ordem], s_bordas, z_bordas, i_lin, i_col)


def grade_abertura(malha: MalhaVTK, plano: PlanoFratura) -> np.ndarray:
    """Abertura (mm) no plano da fratura, como matriz (profundidade x comprimento)."""
    valores = malha.dados_celula[CAMPO_ABERTURA][plano.indices]
    grade = np.full((len(plano.z), len(plano.s)), np.nan)
    grade[plano.i_linha, plano.i_coluna] = valores
    return grade


# ---------------------------------------------------------------------------
# Entrada
# ---------------------------------------------------------------------------
def listar_arquivos(entradas: list[str], poco: str | None) -> list[tuple[int, Path]]:
    arquivos: dict[Path, int] = {}
    for entrada in entradas:
        p = Path(entrada)
        candidatos = sorted(p.glob("*.vtk")) if p.is_dir() else [p]
        for c in candidatos:
            m = PADRAO_PASSO.match(c.name)
            if not m:
                continue  # ex.: <poço>.vtk, só com a malha e dados estáticos
            if poco and m["poco"] != poco:
                continue
            arquivos[c] = int(m["passo"])

    if not arquivos:
        raise SystemExit("Nenhum arquivo <poço>_<passo>.vtk encontrado.")
    pocos = {PADRAO_PASSO.match(c.name)["poco"] for c in arquivos}
    if len(pocos) > 1:
        raise SystemExit(
            "Arquivos de mais de um poço encontrados: "
            + ", ".join(sorted(pocos)) + ". Use --poco para escolher."
        )
    return sorted(((passo, c) for c, passo in arquivos.items()))


# ---------------------------------------------------------------------------
# Resumo, CSV e figura
# ---------------------------------------------------------------------------
def extensao_aberta(bordas: np.ndarray, aberto: np.ndarray) -> float:
    """Extensão (m) do trecho aberto, das bordas das células abertas extremas."""
    if not aberto.any():
        return 0.0
    return float(bordas[aberto].max() - bordas[aberto].min())


def imprimir_resumo(passos, grades, plano: PlanoFratura) -> None:
    print(f"Malha {plano.tipo_malha}: plano da fratura com {len(plano.z)} linhas (Z) "
          f"x {len(plano.s)} colunas ({'X' if plano.tipo_malha == 'BOX' else 'r'})")
    print(f"{'passo':>6} {'abert. máx (mm)':>16} {'em s (m)':>10} {'em z (m)':>10} "
          f"{'altura (m)':>11} {'compr. (m)':>11} {'células abertas':>16}")
    for passo, g in zip(passos, grades):
        aberto = np.nan_to_num(g) > 0
        if aberto.any():
            il, ic = np.unravel_index(np.nanargmax(g), g.shape)
            maximo, s_max, z_max = g[il, ic], plano.s[ic], plano.z[il]
        else:
            maximo, s_max, z_max = 0.0, math.nan, math.nan
        altura = extensao_aberta(plano.z_bordas, aberto.any(axis=1))
        compr = extensao_aberta(plano.s_bordas, aberto.any(axis=0))
        print(f"{passo:>6} {maximo:>16.5f} {s_max:>10.3f} {z_max:>10.3f} "
              f"{altura:>11.3f} {compr:>11.3f} {int(aberto.sum()):>16}")


def exportar_csv(caminho: Path, passos, grades, plano: PlanoFratura) -> None:
    with open(caminho, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["passo", "s_m", "z_m", "abertura_mm"])
        for passo, g in zip(passos, grades):
            for il, z in enumerate(plano.z):
                for ic, s in enumerate(plano.s):
                    if not math.isnan(g[il, ic]):
                        w.writerow([passo, f"{s:.6f}", f"{z:.6f}", f"{g[il, ic]:.6f}"])
    print(f"CSV salvo em {caminho}")


def gerar_figura(caminho: Path, passos, grades, plano: PlanoFratura, z_alvo: float | None,
                 titulo: str, mostrar: bool, zoom: bool = True) -> None:
    import matplotlib
    if not mostrar:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    ultima = grades[-1]
    abertas = np.nan_to_num(np.array(grades)) > 0          # (passos, linhas, colunas)
    houve_abertura = bool(abertas.any())

    # Profundidade do perfil horizontal: a pedida, ou a da maior abertura final.
    if z_alvo is not None:
        il_alvo = int(np.argmin(np.abs(plano.z - z_alvo)))
    elif np.nan_to_num(ultima).max() > 0:
        il_alvo = int(np.unravel_index(np.nanargmax(ultima), ultima.shape)[0])
    else:
        il_alvo = 0
    nome_s = "X (m)" if plano.tipo_malha == "BOX" else "r (m)"

    # Bordas das células em ordem crescente, para os perfis em degraus.
    z_cresc = np.append(plano.z_bordas[::-1, 0], plano.z_bordas[0, 1])
    s_cresc = np.append(plano.s_bordas[:, 0], plano.s_bordas[-1, 1])

    # Faixa de Z mostrada: a região que abriu em algum passo, com margem.
    z_lim = (z_cresc[0], z_cresc[-1])
    if zoom and houve_abertura:
        linhas = abertas.any(axis=(0, 2))
        topo, base = plano.z_bordas[linhas, 1].max(), plano.z_bordas[linhas, 0].min()
        margem = max(0.5 * (topo - base), 2 * float(np.median(np.diff(z_cresc))))
        z_lim = (max(base - margem, z_cresc[0]), min(topo + margem, z_cresc[-1]))

    fig, (ax_v, ax_h, ax_m) = plt.subplots(1, 3, figsize=(16, 6), constrained_layout=True)
    cores = plt.cm.viridis(np.linspace(0, 1, len(passos)))

    for cor, passo, g in zip(cores, passos, grades):
        ax_v.stairs(np.nan_to_num(g[::-1, 0]), z_cresc, orientation="horizontal",
                    baseline=None, color=cor, label=f"passo {passo}")
        ax_h.stairs(np.nan_to_num(g[il_alvo, :]), s_cresc, baseline=None, color=cor,
                    label=f"passo {passo}")

    ax_v.set(xlabel="abertura (mm)", ylabel="Z (m)", ylim=z_lim,
             title=f"Perfil vertical junto ao poço\n(célula de {nome_s[0]} = "
                   f"{plano.s_bordas[0, 0]:.2f} a {plano.s_bordas[0, 1]:.2f} m)")
    ax_h.set(xlabel=nome_s, ylabel="abertura (mm)",
             title=f"Perfil ao longo do comprimento\n(célula de Z = "
                   f"{plano.z_bordas[il_alvo, 1]:.2f} a {plano.z_bordas[il_alvo, 0]:.2f} m)")
    for ax in (ax_v, ax_h):
        ax.grid(True, alpha=0.3)
    ax_v.set_xlim(left=0, right=None if houve_abertura else 1.0)
    ax_h.set_ylim(bottom=0, top=None if houve_abertura else 1.0)
    if len(passos) <= 12:
        ax_h.legend(fontsize=8)

    # Mapa no plano da fratura (último passo), com as bordas reais das células.
    z_desc = np.append(plano.z_bordas[:, 1], plano.z_bordas[-1, 0])
    vmax = float(np.nanmax(ultima)) if houve_abertura else 1.0
    vmax = vmax if vmax > 0 else 1.0
    mapa = ax_m.pcolormesh(s_cresc, z_desc, ultima, shading="flat", cmap="magma",
                           vmin=0, vmax=vmax)
    fig.colorbar(mapa, ax=ax_m, label="abertura (mm)")
    ax_m.axhline(plano.z[il_alvo], color="cyan", lw=0.8, ls="--")
    ax_m.set(xlabel=nome_s, ylabel="Z (m)", ylim=z_lim,
             title=f"Abertura no plano da fratura\n(passo {passos[-1]})")

    if not houve_abertura:
        titulo += "  (sem abertura: a fratura não abriu em nenhum passo)"
    fig.suptitle(titulo)
    fig.savefig(caminho, dpi=150)
    print(f"Figura salva em {caminho}")
    if mostrar:
        plt.show()


# ---------------------------------------------------------------------------
def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        description="Perfil da abertura da fratura a partir dos VTK de passo do iGeo.",
        formatter_class=argparse.RawDescriptionHelpFormatter, epilog=__doc__.split("Exemplos:")[1],
    )
    ap.add_argument("entradas", nargs="+",
                    help="arquivos <poço>_<passo>.vtk ou pastas que os contenham")
    ap.add_argument("--poco", help="nome do poço (prefixo dos arquivos), se houver mais de um")
    ap.add_argument("--malha", choices=["auto", "box", "rad"], default="auto",
                    help="tipo de malha (padrão: detectado pela geometria)")
    ap.add_argument("--z", type=float, help="profundidade (m) do perfil horizontal "
                    "(padrão: a da maior abertura no último passo)")
    ap.add_argument("--passos", help="passos a usar, ex.: '1,3,5' ou '2-6' (padrão: todos)")
    ap.add_argument("--saida", type=Path, help="arquivo da figura (padrão: "
                    "<pasta>/<poço>_perfil_abertura.png)")
    ap.add_argument("--csv", type=Path, help="exporta abertura do plano da fratura em CSV")
    ap.add_argument("--sem-figura", action="store_true", help="não gera a figura")
    ap.add_argument("--mostrar", action="store_true", help="abre a figura numa janela")
    ap.add_argument("--sem-zoom", action="store_true",
                    help="mostra toda a altura do modelo (padrão: só a região que abriu)")
    args = ap.parse_args(argv)

    arquivos = listar_arquivos(args.entradas, args.poco)
    if args.passos:
        escolhidos: set[int] = set()
        for parte in args.passos.split(","):
            a, _, b = parte.partition("-")
            escolhidos.update(range(int(a), int(b or a) + 1))
        arquivos = [(p, c) for p, c in arquivos if p in escolhidos]
        if not arquivos:
            raise SystemExit("Nenhum dos passos pedidos foi encontrado.")

    passos, grades, plano = [], [], None
    for passo, caminho in arquivos:
        malha = ler_vtk(caminho)
        if CAMPO_ABERTURA not in malha.dados_celula:
            print(f"Aviso: {caminho.name} não tem o campo {CAMPO_ABERTURA}; ignorado.",
                  file=sys.stderr)
            continue
        if plano is None:  # a malha é a mesma em todos os passos
            plano = identificar_plano(malha, args.malha)
        passos.append(passo)
        grades.append(grade_abertura(malha, plano))
    if not grades:
        raise SystemExit(f"Nenhum arquivo com o campo {CAMPO_ABERTURA}.")

    poco = PADRAO_PASSO.match(arquivos[0][1].name)["poco"]
    imprimir_resumo(passos, grades, plano)
    if args.csv:
        exportar_csv(args.csv, passos, grades, plano)
    if not args.sem_figura:
        saida = args.saida or arquivos[0][1].with_name(f"{poco}_perfil_abertura.png")
        gerar_figura(saida, passos, grades, plano, args.z,
                     f"iGeo — abertura da fratura: {poco}", args.mostrar,
                     zoom=not args.sem_zoom)
    return 0


if __name__ == "__main__":
    sys.exit(main())

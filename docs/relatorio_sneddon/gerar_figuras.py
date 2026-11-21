#!/usr/bin/env python3
"""Gera as figuras do relatório da verificação de Sneddon (docs/relatorio_sneddon).

Uso:
    python3 gerar_figuras.py <pasta_dados>

<pasta_dados> deve conter as simulações usadas no relatório:
    rel/antigo, rel/novo          caso do teste com a pressão mantida por 6 passos
                                  (versão anterior e versão corrigida do iGeo)
    conv/D50_dx0.5, conv/D50_dx0.25, conv/D50_dx0.125,
    tempo/D100_dx0.5, conv/D100_dx0.25, conv/D100_dx0.125, conv/D200_dx0.25
                                  estudo de convergência (versão corrigida)
    prop/antigo, prop/novo        caso de propagação M3_injection_test (opcional)

Cada figura é gravada em figuras/ como SVG (usada no README.md) e como PDF
(usada por `make pdf`, ver Makefile). Requer numpy e matplotlib.
"""

from __future__ import annotations

import os
import re
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.patches import FancyArrowPatch, Polygon, Rectangle  # noqa: E402

RAIZ = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(RAIZ / "tests"))
sys.path.insert(0, str(RAIZ / "tools"))
import test_sneddon_injection as t  # noqa: E402
from perfil_abertura_fratura import ler_vtk  # noqa: E402

SAIDA = Path(__file__).resolve().parent / "figuras"

# Cores (paleta categórica de referência, validada): série 1 e 2, tinta e grade.
NOVO = "#2a78d6"      # iGeo corrigido
ANTIGO = "#eb6834"    # iGeo antes da correção
TINTA = "#0b0b0b"
TINTA_2 = "#52514e"
GRADE = "#e4e3df"
SUPERFICIE = "#fcfcfb"
KGF = t.KGF_CM2_PARA_KPA

plt.rcParams.update({
    "figure.facecolor": SUPERFICIE, "axes.facecolor": SUPERFICIE,
    "savefig.facecolor": SUPERFICIE, "font.size": 10,
    "axes.edgecolor": TINTA_2, "axes.labelcolor": TINTA, "axes.titlecolor": TINTA,
    "xtick.color": TINTA_2, "ytick.color": TINTA_2, "text.color": TINTA,
    "axes.grid": True, "grid.color": GRADE, "grid.linewidth": 0.8,
    "axes.spines.top": False, "axes.spines.right": False,
    "legend.frameon": False, "lines.linewidth": 1.8, "svg.fonttype": "none",
    "pdf.fonttype": 42, "svg.hashsalt": "igeo-sneddon",
})


def salvar(fig, nome: str) -> None:
    SAIDA.mkdir(parents=True, exist_ok=True)
    # SVG para o README (GitHub) e PDF vetorial para o relatório em PDF (make pdf).
    # metadata sem data, para que figuras iguais gerem arquivos iguais.
    fig.savefig(SAIDA / nome, bbox_inches="tight", metadata={"Date": None})
    fig.savefig(SAIDA / nome.replace(".svg", ".pdf"), bbox_inches="tight",
                metadata={"CreationDate": None, "ModDate": None})
    previa = os.environ.get("PREVIA_PNG")   # pasta opcional para cópias PNG
    if previa:
        Path(previa).mkdir(parents=True, exist_ok=True)
        fig.savefig(Path(previa) / nome.replace(".svg", ".png"), bbox_inches="tight", dpi=110)
    plt.close(fig)
    print("figura:", SAIDA / nome)


def faces(pasta: Path, passos: int) -> list[dict]:
    return [t.face_da_fratura(pasta / f"{t.POCO}_{k + 1}.vtk") for k in range(1, passos + 1)]


# ---------------------------------------------------------------------------
def fig_esquema() -> None:
    """Domínio (1/4 do problema), condições de contorno e fratura."""
    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    ax.set_aspect("equal")
    ax.grid(False)
    for lado in ("left", "bottom"):
        ax.spines[lado].set_visible(False)
    ax.set_xticks([])
    ax.set_yticks([])
    lx, ly, lf = 10.0, 6.0, 3.0   # desenho fora de escala
    ax.add_patch(Rectangle((0, 0), lx, ly, facecolor="#f1f0ec", edgecolor=TINTA_2, lw=1.2))
    # Fratura (face y = 0, 0 <= x <= L) com pressão
    ax.plot([0, lf], [0, 0], color=NOVO, lw=4, solid_capstyle="butt")
    for x in np.linspace(0.3, lf - 0.3, 6):
        ax.add_patch(FancyArrowPatch((x, -0.75), (x, -0.08), arrowstyle="-|>",
                                     mutation_scale=9, color=NOVO, lw=1.2))
    ax.text(lf / 2, -1.05, "pressão do fluido p (uniforme)", ha="center", va="top", color=TINTA)
    # Apoios: simetria em x = 0 e em y = 0 fora da fratura; contornos externos fixos.
    def roletes(x0, y0, x1, y1, normal):
        for s in np.linspace(0.15, 0.85, 7):
            x, y = x0 + s * (x1 - x0), y0 + s * (y1 - y0)
            ax.plot(x + 0.18 * normal[0], y + 0.18 * normal[1], "o", ms=4.5,
                    mfc="none", mec=TINTA_2, mew=1)
    roletes(0, 0, 0, ly, (-1, 0))
    roletes(lf, 0, lx, 0, (0, -1))
    roletes(lx, 0, lx, ly, (1, 0))
    roletes(0, ly, lx, ly, (0, 1))
    ax.plot(0, 0, marker="o", ms=9, color=TINTA)
    ax.annotate("poço", (0, 0), (-1.5, 0.7), color=TINTA,
                arrowprops=dict(arrowstyle="-", color=TINTA_2, lw=0.8))
    ax.text(-0.45, ly / 2, "simetria: ux = 0", rotation=90, ha="right", va="center", color=TINTA_2)
    ax.text((lf + lx) / 2, -0.45, "simetria: uy = 0", ha="center", va="top", color=TINTA_2)
    ax.text(lx + 0.45, ly / 2, "ux = 0", rotation=90, ha="left", va="center", color=TINTA_2)
    ax.text(lx / 2, ly + 0.45, "uy = 0", ha="center", va="bottom", color=TINTA_2)
    ax.annotate("", (lf, 0.55), (0, 0.55), arrowprops=dict(arrowstyle="<->", color=TINTA, lw=0.9))
    ax.text(lf / 2, 0.7, "L = 5 m", ha="center", va="bottom")
    ax.text(lx * 0.62, ly * 0.55,
            "rocha elástica\nE = 10 GPa, ν = 0,25\nBiot ≈ 0, k = 1e-5 mD\n"
            "Shmin = 10 MPa (normal a y)\ndeformação plana (uz = 0)",
            ha="center", va="center", color=TINTA, linespacing=1.4)
    ax.annotate("", (lx, ly + 1.15), (0, ly + 1.15), arrowprops=dict(arrowstyle="<->", color=TINTA_2, lw=0.8))
    ax.text(lx / 2, ly + 1.3, "100 m (20 L)", ha="center", va="bottom", color=TINTA_2)
    ox, oy = 0.6, ly - 1.6
    ax.annotate("", (ox + 1.0, oy), (ox, oy), arrowprops=dict(arrowstyle="->", color=TINTA, lw=1))
    ax.annotate("", (ox, oy + 1.0), (ox, oy), arrowprops=dict(arrowstyle="->", color=TINTA, lw=1))
    ax.text(ox + 1.12, oy, "x", va="center")
    ax.text(ox, oy + 1.12, "y", ha="center", va="bottom")
    ax.set_xlim(-3.4, lx + 1.2)
    ax.set_ylim(-2.1, ly + 2.0)
    salvar(fig, "fig1_esquema.svg")


def fig_malha(pasta: Path) -> None:
    """Malha junto à fratura, vista em planta (topo da camada)."""
    malha = ler_vtk(pasta / f"{t.POCO}_1.vtk")
    zmax = malha.pontos[:, 2].max()
    fig, ax = plt.subplots(figsize=(7.2, 3.0))
    ax.grid(False)
    for cel in malha.celulas:
        topo = cel[np.abs(malha.pontos[cel, 2] - zmax) < 1e-6]
        xy = malha.pontos[topo, :2]
        if xy[:, 0].min() > 10.5 or xy[:, 1].min() > 2.6:
            continue
        centro = xy.mean(axis=0)
        ordem = np.argsort(np.arctan2(xy[:, 1] - centro[1], xy[:, 0] - centro[0]))
        frat = xy[:, 1].min() < 1e-9 and xy[:, 0].max() <= t.L + 1e-9
        ax.add_patch(Polygon(xy[ordem], closed=True, facecolor="#cde2fb" if frat else "none",
                             edgecolor=TINTA_2, lw=0.5))
    ax.plot([0, t.L], [0, 0], color=NOVO, lw=3, solid_capstyle="butt", label="fratura (20 elementos)")
    ax.set(xlim=(0, 10), ylim=(-0.05, 2.5), xlabel="x (m)", ylabel="y (m)", aspect="equal")
    ax.legend(loc="lower right", bbox_to_anchor=(1.0, 1.02), borderaxespad=0)
    salvar(fig, "fig2_malha.svg")


def fig_perfis(antigo: list[dict], novo: list[dict]) -> None:
    """Perfil da abertura: passo 1 e passo 6 (pressão mantida)."""
    x = np.linspace(0, t.L, 300)
    fig, axes = plt.subplots(1, 2, figsize=(9.6, 3.8), sharey=True)
    for ax, k, titulo in ((axes[0], 0, "(a) fim do 1º passo"),
                          (axes[1], 5, "(b) fim do 6º passo (pressão mantida)")):
        ax.plot(x, t.uy_sneddon(x) * 1e3, color=TINTA, lw=1.4, ls="--", label="Sneddon")
        for faces_v, cor, nome in ((antigo, ANTIGO, "iGeo antes da correção"),
                                   (novo, NOVO, "iGeo corrigido")):
            f = faces_v[k]
            m = f["x"] <= t.L + 1e-9
            ax.plot(f["x"][m], f["uy"][m] * 1e3, color=cor, marker="o", ms=4, label=nome)
        ax.set(title=titulo, xlabel="x (m)", xlim=(0, t.L + 0.2))
    axes[0].set_ylabel("uy na face = meia abertura (mm)")
    axes[0].legend(loc="lower left")
    salvar(fig, "fig3_perfis.svg")


def fig_deriva(antigo: list[dict], novo: list[dict]) -> None:
    """Abertura no centro / Sneddon ao longo dos passos com pressão mantida."""
    ref = t.uy_sneddon(np.array([0.0]))[0]
    passos = np.arange(1, len(antigo) + 1)
    fig, ax = plt.subplots(figsize=(6.4, 3.6))
    ax.plot([passos[0] - 0.3, passos[-1]], [1.0, 1.0], color=TINTA, lw=1.2, ls="--")
    ax.text(passos[-1] + 0.15, 1.0, "Sneddon", va="center", color=TINTA)
    for faces_v, cor, nome in ((antigo, ANTIGO, "antes da correção"), (novo, NOVO, "corrigido")):
        r = np.array([f["uy"][0] for f in faces_v]) / ref
        ax.plot(passos, r, color=cor, marker="o", ms=6, label=nome)
        ax.text(passos[-1] + 0.15, r[-1], f"{nome} ({r[-1] - 1:+.1%})", va="center", color=TINTA)
    ax.set(xlabel="passo de tempo (pressão no poço mantida em 11 MPa)",
           ylabel="uy(0) / Sneddon", xticks=passos, xlim=(0.7, passos[-1] + 2.6))
    ax.legend(loc="upper left")
    salvar(fig, "fig4_deriva.svg")


def fig_residuo(pasta: Path) -> None:
    """Resíduo p + Syy médio por elemento da face, no fim do 1º passo."""
    malha = ler_vtk(pasta / f"{t.POCO}_2.vtk")
    zmax = malha.pontos[:, 2].max()
    xs, res = [], []
    for k, cel in enumerate(malha.celulas):
        p = malha.pontos[cel]
        if p[:, 1].min() < 1e-9 and p[:, 0].max() <= t.L + 1e-9 and np.abs(p[:, 2] - zmax).min() < 1e-6:
            xs.append(p[:, 0].mean())
            res.append((malha.dados_celula["FracPressure"][k] + malha.dados_celula["SyyeffecAVG_kgf_cm2"][k]) * KGF)
    xs, res = np.array(xs), np.array(res)
    ordem = np.argsort(xs)
    xs, res = xs[ordem], res[ordem]
    fig, ax = plt.subplots(figsize=(6.8, 3.6))
    ax.bar(xs, res, width=0.2, color=ANTIGO)
    ax.axhline(0, color=TINTA_2, lw=0.8)
    junto = xs > 4.0
    longe = np.abs(res[~junto]).max()
    ax.text(0.15, 0.62 * res.max(),
            f"longe da ponta (x < 4 m): |resíduo| < {max(longe, 1):.0f} kPa\n"
            f"junto à ponta: oscila entre {res[junto].min():+.0f} e {res[junto].max():+.0f} kPa",
            color=TINTA, linespacing=1.5)
    ax.set(xlabel="centro do elemento da face (m)",
           ylabel="p + σyy médio do elemento (kPa)", xlim=(0, t.L))
    salvar(fig, "fig5_residuo.svg")


def fig_convergencia(dados: Path) -> None:
    """|erro| x dx para dois tamanhos de domínio (versão corrigida)."""
    series = {
        "domínio 10L (50 m)": [(0.5, dados / "conv/D50_dx0.5"), (0.25, dados / "conv/D50_dx0.25"),
                               (0.125, dados / "conv/D50_dx0.125")],
        "domínio 20L (100 m)": [(0.5, dados / "tempo/D100_dx0.5"), (0.25, dados / "conv/D100_dx0.25"),
                                (0.125, dados / "conv/D100_dx0.125")],
    }
    nomes = (("centro", "abertura no centro"), ("perfil_l2", "perfil (norma L2)"), ("volume", "volume"))
    fig, axes = plt.subplots(1, 3, figsize=(10.4, 3.4), sharey=True)
    estilos = (dict(ls="--", marker="s"), dict(ls="-", marker="o"))
    for (estilo, (rotulo, pontos)) in zip(estilos, series.items()):
        dx = np.array([p[0] for p in pontos])
        e = [t.erros(t.face_da_fratura(p[1] / f"{t.POCO}_2.vtk")) for p in pontos]
        for ax, (chave, titulo) in zip(axes, nomes):
            ax.plot(dx, [abs(x[chave]) * 100 for x in e], color=NOVO, ms=6, label=rotulo, **estilo)
            ax.set(title=titulo, xlabel="dx (m)", xscale="log", yscale="log")
    for ax in axes:
        ax.set_xticks([0.125, 0.25, 0.5], ["0,125", "0,25", "0,5"])
        ax.minorticks_off()
        ax.set_yticks([1, 2, 3, 5], ["1", "2", "3", "5"])
    axes[0].set_ylabel("|erro| em relação a Sneddon (%)")
    axes[2].legend(loc="lower right")
    salvar(fig, "fig6_convergencia.svg")


def fig_propagacao(dados: Path) -> None:
    """Caso de propagação M3_injection_test: pressão máxima e eventos de propagação."""
    padrao = re.compile(r"\s*CGM:\s+([\d.E+-]+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+\d+\s+\d+\s+([\d.E+-]+)")
    curvas = {}
    for v in ("antigo", "novo"):
        arq = dados / "prop" / v / "out.txt"
        if not arq.exists():
            return
        t_min, pmax, eventos = [], [], set()
        for linha in arq.read_text(errors="replace").splitlines():
            m = padrao.match(linha)
            if not m:
                continue
            if m.group(4) == "0":                 # passo de fluxo
                t_min.append(float(m.group(1)))
                pmax.append(float(m.group(6)) / 1000.0)
            elif m.group(5) == "2":               # iteração de propagação da fratura
                eventos.add(float(m.group(1)))
        curvas[v] = (np.array(t_min), np.array(pmax), np.array(sorted(eventos)))
    t_fim = min(curvas["antigo"][0].max(), curvas["novo"][0].max())
    fig, ax = plt.subplots(figsize=(7.2, 4.0))
    for i, (v, cor, nome) in enumerate((("antigo", ANTIGO, "antes da correção"),
                                        ("novo", NOVO, "corrigido"))):
        tm, p, ev = curvas[v]
        m = tm <= t_fim
        ax.plot(tm[m], p[m], color=cor, marker="o", ms=3, label=nome)
        ev = ev[ev <= t_fim]
        y = 24.8 - 0.6 * i
        ax.plot(ev, np.full_like(ev, y), ls="none", marker="|", ms=9, mew=1.6, color=cor)
        rotulo = "evento" if len(ev) == 1 else "eventos"
        ax.text(t_fim + 0.4, y, f"{len(ev)} {rotulo} de propagação", va="center", color=TINTA, fontsize=9)
    ax.set(xlabel="tempo (min)", ylabel="pressão máxima no modelo (MPa)",
           xlim=(-0.5, t_fim + 7.5), ylim=(23.7, 31.2))
    ax.legend(loc="center right")
    salvar(fig, "fig7_propagacao.svg")


def main() -> int:
    dados = Path(sys.argv[1])
    antigo = faces(dados / "rel/antigo", 6)
    novo = faces(dados / "rel/novo", 6)
    fig_esquema()
    fig_malha(dados / "rel/novo")
    fig_perfis(antigo, novo)
    fig_deriva(antigo, novo)
    fig_residuo(dados / "rel/antigo")
    fig_convergencia(dados)
    fig_propagacao(dados)
    return 0


if __name__ == "__main__":
    sys.exit(main())

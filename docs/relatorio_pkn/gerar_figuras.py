#!/usr/bin/env python3
"""Gera as figuras do relatório da verificação PKN (docs/relatorio_pkn).

Uso:
    python3 gerar_figuras.py <pasta_dados>

<pasta_dados> deve conter as simulações usadas no relatório, cada uma numa
pasta com o .geo e as saídas do iGeo:
    etapa_0          caso do teste com a versão anterior do iGeo (commit 7cc39f5)
    etapa_A          + balanço de massa da fratura (FluidFlowSolver)
    etapa_AB         + lei cúbica com a abertura total
    etapa_ABC        + liberação de nós na propagação
    etapa_final      + carga da face com a força interna K*u (versão atual)
    prop_gc0.001_dx2 versão atual, Gc = 0,001 J/m2
    prop_gc0.001_dx1 versão atual, Gc = 0,001 J/m2 e dx = 1 m (NX = 120)

Cada figura é gravada em figuras/ como SVG (usada no README.md) e como PDF
(usada por `make pdf`). Requer numpy e matplotlib.
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.patches import Ellipse, FancyArrowPatch, Rectangle  # noqa: E402

RAIZ = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(RAIZ / "tests"))
import test_pkn_injection as t  # noqa: E402

SAIDA = Path(__file__).resolve().parent / "figuras"

# Paleta categórica de referência (validada): séries 1 a 3, tinta e grade.
AZUL = "#2a78d6"       # versão atual
LARANJA = "#eb6834"    # versão anterior / etapa intermediária
VERDE_AGUA = "#1baf7a"
TINTA = "#0b0b0b"
TINTA_2 = "#52514e"
GRADE = "#e4e3df"
SUPERFICIE = "#fcfcfb"

plt.rcParams.update({
    "figure.facecolor": SUPERFICIE, "axes.facecolor": SUPERFICIE,
    "savefig.facecolor": SUPERFICIE, "font.size": 10,
    "axes.edgecolor": TINTA_2, "axes.labelcolor": TINTA, "axes.titlecolor": TINTA,
    "xtick.color": TINTA_2, "ytick.color": TINTA_2, "text.color": TINTA,
    "axes.grid": True, "grid.color": GRADE, "grid.linewidth": 0.8,
    "axes.spines.top": False, "axes.spines.right": False,
    "legend.frameon": False, "lines.linewidth": 1.8, "svg.fonttype": "none",
    "pdf.fonttype": 42, "svg.hashsalt": "igeo-pkn",
})

TEMPOS = t.DT * np.arange(1, t.PASSOS + 1)


def salvar(fig, nome: str) -> None:
    SAIDA.mkdir(parents=True, exist_ok=True)
    # metadata sem data, para que figuras iguais gerem arquivos iguais.
    fig.savefig(SAIDA / f"{nome}.svg", bbox_inches="tight", metadata={"Date": None})
    fig.savefig(SAIDA / f"{nome}.pdf", bbox_inches="tight",
                metadata={"CreationDate": None, "ModDate": None})
    plt.close(fig)
    print("figura:", SAIDA / f"{nome}.svg")


def faces(pasta: Path) -> list[dict]:
    return [t.face_da_fratura(pasta / f"{t.POCO}_{k + 1}.vtk") for k in range(1, t.PASSOS + 1)]


def comprimento_aberto(face: dict) -> float:
    """Posição do primeiro nó fechado na linha central (fim da região aberta)."""
    w = face["uy"][t.linha_central(face)]
    aberto = np.flatnonzero(w > 1e-9)
    return float(face["x"][aberto.max() + 1]) if aberto.size else 0.0


def no_poco(face: dict) -> tuple[float, float]:
    jz = t.linha_central(face)
    return 2.0 * face["uy"][jz, 0], face["pnet"][jz, 0]


def nordgren_perfil(x: np.ndarray, tempo: float) -> np.ndarray:
    """Perfil aproximado do PKN, w0 (1 - x/L)^(1/4)."""
    ref = t.nordgren(tempo)
    return ref["w0"] * np.clip(1.0 - x / ref["L"], 0.0, None) ** 0.25


def fig_esquema() -> None:
    fig, (a, b) = plt.subplots(1, 2, figsize=(9.0, 3.4), gridspec_kw={"width_ratios": [2.3, 1]})
    for ax in (a, b):
        ax.set_axis_off()
    # (a) face y = 0 em planta (x, z)
    a.add_patch(Rectangle((0, -3.2), 10, 6.4, facecolor="#f3f2ef", edgecolor=TINTA_2, lw=0.8))
    a.add_patch(Rectangle((0, -1.0), 6.2, 2.0, facecolor="#cde2fb", edgecolor=AZUL, lw=1.5))
    a.plot([0, 0], [-3.2, 3.2], color=TINTA, lw=2.2)
    a.text(0.15, 2.75, "poço (x = 0, simetria)", fontsize=9, color=TINTA)
    a.annotate("", xy=(6.2, -1.45), xytext=(0, -1.45),
               arrowprops=dict(arrowstyle="<->", color=TINTA, lw=1))
    a.text(3.1, -1.85, "L(t)", ha="center", fontsize=10)
    a.annotate("", xy=(6.6, 1.0), xytext=(6.6, -1.0), arrowprops=dict(arrowstyle="<->", color=TINTA, lw=1))
    a.text(6.8, 0, "H = 10 m\n(camadas com\n#CDM_FRACTURE)", va="center", fontsize=9)
    a.text(3.0, 0, "fratura", ha="center", va="center", fontsize=10, color=TINTA)
    a.text(5.0, -2.6, "camadas acima e abaixo: mesma rocha, sem propagação", ha="center", fontsize=8.5, color=TINTA_2)
    a.annotate("", xy=(1.4, 0.55), xytext=(0.05, 0.55),
               arrowprops=dict(arrowstyle="->", color=AZUL, lw=1.6))
    a.text(0.25, 0.68, "q/2", fontsize=9, color=TINTA)
    a.set_xlim(-0.3, 10.2)
    a.set_ylim(-3.4, 3.4)
    a.set_title("(a) face da fratura (y = 0)", fontsize=10, loc="left")
    # (b) seção vertical: elipse
    b.add_patch(Ellipse((0, 0), 1.2, 4.0, facecolor="#cde2fb", edgecolor=AZUL, lw=1.5))
    b.plot([0, 0], [-2.6, 2.6], color=TINTA_2, lw=0.8, ls=":")
    b.annotate("", xy=(0.6, -0.0), xytext=(-0.6, 0.0), arrowprops=dict(arrowstyle="<->", color=TINTA, lw=1))
    b.text(0.75, 0.05, "w(x)", fontsize=10, va="center")
    b.annotate("", xy=(-1.0, 2.0), xytext=(-1.0, -2.0), arrowprops=dict(arrowstyle="<->", color=TINTA, lw=1))
    b.text(-1.15, 0, "H", ha="right", va="center", fontsize=10)
    b.text(0, -2.75, "modelo: só y ≥ 0\n(simetria em y = 0)", ha="center", va="top", fontsize=8.5, color=TINTA_2)
    b.set_xlim(-1.6, 1.6)
    b.set_ylim(-3.4, 3.0)
    b.set_aspect("equal")
    b.set_title("(b) seção vertical", fontsize=10, loc="left")
    salvar(fig, "fig1_esquema")


def fig_antes(dados: Path) -> None:
    antes, atual = faces(dados / "etapa_0"), faces(dados / "etapa_final")
    fig, (a, b) = plt.subplots(1, 2, figsize=(9.0, 3.4))
    v_inj = t.Q_MODELO * TEMPOS
    a.semilogy(TEMPOS, [t.volume_modelo(f) for f in antes], "o-", color=LARANJA, ms=4, label="versão anterior")
    a.semilogy(TEMPOS, [t.volume_modelo(f) for f in atual], "o-", color=AZUL, ms=4, label="versão atual")
    a.semilogy(TEMPOS, v_inj, "--", color=TINTA, lw=1.2, label="volume injetado")
    a.set_xlabel("tempo (s)")
    a.set_ylabel("volume da fratura no modelo (m³)")
    a.set_title("(a) volume da fratura", fontsize=10, loc="left")
    a.legend(loc="center right")
    p_ref = [t.nordgren(s)["pnet"] / 1e6 for s in TEMPOS]
    b.semilogy(TEMPOS, [no_poco(f)[1] / 1e6 for f in antes], "o-", color=LARANJA, ms=4, label="versão anterior")
    b.semilogy(TEMPOS, [no_poco(f)[1] / 1e6 for f in atual], "o-", color=AZUL, ms=4, label="versão atual")
    b.semilogy(TEMPOS, p_ref, "--", color=TINTA, lw=1.2, label="Nordgren")
    b.set_xlabel("tempo (s)")
    b.set_ylabel("pressão líquida no poço (MPa)")
    b.set_title("(b) pressão líquida no poço", fontsize=10, loc="left")
    b.legend(loc="center right")
    fig.tight_layout()
    salvar(fig, "fig2_antes")


def fig_lei_cubica(dados: Path) -> None:
    """Coeficiente da lei cúbica local, antes e depois da correção."""
    def coeficientes(pasta: Path) -> list[float]:
        fs = faces(pasta)
        a, b = fs[-2], fs[-1]
        dentro = np.abs(b["z"] - t.Z_CENTRO) <= t.H / 2 + 1e-6
        out = []
        for xa, xb in ((4.0, 8.0), (8.0, 12.0), (12.0, 16.0)):
            xc = 0.5 * (xa + xb)
            q = 2.0 * (t.Q_MODELO - (t.volume_modelo(b, xc) - t.volume_modelo(a, xc)) / t.DT)
            ia, ib = (int(np.argmin(np.abs(b["x"] - v))) for v in (xa, xb))
            p = 0.5 * (a["pnet"] + b["pnet"])[dentro]
            uy = 0.5 * (a["uy"] + b["uy"])[dentro]
            ab = uy[:, ia] > 0
            dpdx = (p[ab, ib].mean() - p[ab, ia].mean()) / (xb - xa)
            w = uy[:, ia] + uy[:, ib]
            out.append(q * 12.0 * t.MU / (np.trapezoid(w**3, b["z"][dentro]) * -dpdx))
        return out
    antes, depois = coeficientes(dados / "etapa_A"), coeficientes(dados / "etapa_final")
    fig, ax = plt.subplots(figsize=(6.0, 3.2))
    pos = np.arange(3)
    ax.bar(pos - 0.19, antes, 0.36, color=LARANJA, label="meia abertura (antes)")
    ax.bar(pos + 0.19, depois, 0.36, color=AZUL, label="abertura total (versão atual)")
    ax.axhline(1.0, color=TINTA, lw=1.0, ls="--")
    ax.text(2.93, 1.03, "lei cúbica exata", ha="right", va="bottom", fontsize=8.5, color=TINTA_2)
    for x, v in zip(pos - 0.19, antes):
        ax.text(x, v + 0.03, f"{v:.2f}", ha="center", fontsize=8.5)
    for x, v in zip(pos + 0.19, depois):
        ax.text(x, v - 0.1, f"{v:.2f}", ha="center", fontsize=8.5, color="white")
    ax.set_xticks(pos, ["x = 4 a 8 m", "x = 8 a 12 m", "x = 12 a 16 m"])
    ax.set_ylabel("C = Q·12μ / (−dp/dx ∫w³dz)")
    ax.set_ylim(0, 1.25)
    ax.set_xlim(-0.6, 2.95)
    ax.grid(axis="x", visible=False)
    ax.legend(loc="upper left", ncol=2)
    salvar(fig, "fig3_lei_cubica")


def fig_altura(dados: Path) -> None:
    sem, com = faces(dados / "etapa_AB")[-1], faces(dados / "etapa_final")[-1]
    fig, ax = plt.subplots(figsize=(5.0, 3.6))
    x = 8.0
    for f, cor, rot in ((sem, LARANJA, "sem a correção dos nós"), (com, AZUL, "com a correção")):
        ix = int(np.argmin(np.abs(f["x"] - x)))
        sel = np.abs(f["z"] - t.Z_CENTRO) <= 8.0
        ax.plot(2e3 * f["uy"][sel, ix], f["z"][sel] - t.Z_CENTRO, "o-", color=cor, ms=4, label=rot)
    z = np.linspace(-t.H / 2, t.H / 2, 200)
    wmax = 2e3 * com["uy"][t.linha_central(com), int(np.argmin(np.abs(com["x"] - x)))]
    ax.plot(wmax * np.sqrt(1 - (2 * z / t.H) ** 2), z, "--", color=TINTA, lw=1.2, label="elipse (altura H)")
    ax.axhspan(-t.H / 2, t.H / 2, color="#f3f2ef", zorder=0)
    ax.text(1.0, 6.6, "fora das camadas com #CDM_FRACTURE", fontsize=8.5, color=TINTA_2)
    ax.set_xlabel("abertura w (mm), em x = 8 m")
    ax.set_ylabel("z − z_centro (m)")
    ax.legend(loc="center left", fontsize=8.5)
    salvar(fig, "fig4_altura")


def fig_carga(dados: Path) -> None:
    inc, atual = faces(dados / "etapa_ABC")[-1], faces(dados / "etapa_final")[-1]
    fig, (a, b) = plt.subplots(1, 2, figsize=(9.0, 3.4))
    for f, cor, rot in ((inc, LARANJA, "carga incremental (etapa C)"), (atual, AZUL, "carga com K·u (atual)")):
        jz = t.linha_central(f)
        a.plot(f["x"], 2e3 * f["uy"][jz], "o-", color=cor, ms=3.5, label=rot)
        frente = int(np.searchsorted(f["x"], t.frente_molhada(f)))
        razao = f["pnet"][jz, :frente - 1] / (t.G * 2.0 * f["uy"][jz, :frente - 1] / ((1 - t.NU) * t.H))
        b.plot(f["x"][:frente - 1], razao, "o-", color=cor, ms=3.5, label=rot)
    a.set_xlim(0, 44)
    a.set_xlabel("x (m)")
    a.set_ylabel("abertura no centro da altura (mm)")
    a.set_title("(a) perfil ao longo da fratura, t = 120 s", fontsize=10, loc="left")
    a.legend(loc="lower left", fontsize=8.5)
    b.axhline(1.0, color=TINTA, lw=1.0, ls="--")
    b.axhline(1 / 0.95, color=TINTA_2, lw=1.0, ls=":", label="verificação estática da malha (1/0,95)")
    b.set_xlabel("x (m)")
    b.set_ylabel("pnet / (G w / ((1−ν) H))")
    b.set_title("(b) relação pressão-abertura", fontsize=10, loc="left")
    b.legend(loc="upper left", fontsize=8.5)
    fig.tight_layout()
    salvar(fig, "fig5_carga_face")


def fig_tempo(dados: Path) -> None:
    fs = faces(dados / "etapa_final")
    ref = [t.nordgren(s) for s in TEMPOS]
    fig, axs = plt.subplots(2, 2, figsize=(9.0, 6.0))
    series = [
        ("frente do fluido L (m)", [t.frente_molhada(f) for f in fs], [r["L"] for r in ref]),
        ("abertura no poço w0 (mm)", [1e3 * no_poco(f)[0] for f in fs], [1e3 * r["w0"] for r in ref]),
        ("pressão líquida no poço (MPa)", [no_poco(f)[1] / 1e6 for f in fs], [r["pnet"] / 1e6 for r in ref]),
        ("volume no modelo (m³)", [t.volume_modelo(f) for f in fs], list(t.Q_MODELO * TEMPOS)),
    ]
    for ax, (rot, igeo, nord), letra in zip(axs.flat, series, "abcd"):
        ax.plot(TEMPOS, igeo, "o-", color=AZUL, ms=4, label="iGeo")
        ax.plot(TEMPOS, nord, "--", color=TINTA, lw=1.2,
                label="volume injetado" if letra == "d" else "Nordgren")
        ax.set_title(f"({letra}) {rot}", fontsize=10, loc="left")
        ax.set_xlabel("tempo (s)")
        ax.set_ylim(bottom=0)
        ax.legend(loc={"a": "upper left", "b": "lower right", "c": "upper right", "d": "upper left"}[letra],
                  fontsize=8.5)
    fig.tight_layout()
    salvar(fig, "fig6_comparacao_nordgren")


def fig_propagacao(dados: Path) -> None:
    fig, ax = plt.subplots(figsize=(6.4, 3.8))
    casos = (("etapa_final", AZUL, "Gc = 1 J/m², dx = 2 m (caso do teste)"),
             ("prop_gc0.001_dx2", LARANJA, "Gc = 0,001 J/m², dx = 2 m"),
             ("prop_gc0.001_dx1", VERDE_AGUA, "Gc = 0,001 J/m², dx = 1 m"))
    for pasta, cor, rot in casos:
        ax.plot(TEMPOS, [comprimento_aberto(f) for f in faces(dados / pasta)], "o-", color=cor, ms=4, label=rot)
    ax.plot(TEMPOS, [t.nordgren(s)["L"] for s in TEMPOS], "--", color=TINTA, lw=1.2, label="Nordgren")
    ax.set_xlabel("tempo (s)")
    ax.set_ylabel("comprimento aberto (m)")
    ax.set_ylim(bottom=0)
    ax.legend(loc="upper left", fontsize=8.5)
    salvar(fig, "fig7_propagacao")


def main() -> int:
    dados = Path(sys.argv[1])
    fig_esquema()
    fig_antes(dados)
    fig_lei_cubica(dados)
    fig_altura(dados)
    fig_carga(dados)
    fig_tempo(dados)
    fig_propagacao(dados)
    return 0


if __name__ == "__main__":
    sys.exit(main())

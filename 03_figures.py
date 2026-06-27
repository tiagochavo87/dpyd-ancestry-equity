"""
03_figures.py
Final publication composite (CHOL style), in English:
  figure_composite.png
    Panel A: DPYD actionable-variant frequency by ancestry (heatmap; includes
             Brazil/ABraOM column when step 05 has produced brazil_freq.csv).
    Panel B: risk carriers missed by the standard (Eurocentric) panel, with the
             out-of-panel risk variant labeled inside the missed segment.
Title, panel letters and caption included. No seaborn dependency.
The full set of alternative styles lives in viz_all.py.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.gridspec import GridSpec

import config as C

CAPT, MISS = "#0b7a75", "#d1495b"
POPS = ["AFR", "AMR", "EAS", "EUR", "SAS"]

plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11,
                     "axes.titleweight": "bold",
                     "axes.spines.top": False, "axes.spines.right": False})


def load():
    global POPS
    freq = pd.read_csv(C.ALLELE_FREQ_CSV, index_col=0)
    cov = pd.read_csv(C.PANEL_COVERAGE_CSV, index_col=0)
    var = pd.read_csv(C.VARIANTS_CSV, index_col=0)
    brazil = False
    try:  # coluna BR (ABraOM), se a etapa 05 rodou
        br = pd.read_csv(C.BRAZIL_FREQ_CSV, index_col=0)
        freq = freq.join(br["BR"])
        if "BR" not in POPS:
            POPS = POPS + ["BR"]
        brazil = True
    except Exception:
        pass
    order = freq[POPS].max(axis=1).sort_values().index
    return freq.loc[order], cov, var.reindex(order), brazil


def panel_freq(ax, freq, var):
    m = freq[POPS].values * 100
    im = ax.imshow(m, cmap="Reds", aspect="auto")
    ax.set_xticks(range(len(POPS))); ax.set_xticklabels(POPS)
    ax.set_yticks(range(len(freq)))
    ax.set_yticklabels([f"{rs}\n{var.loc[rs,'hgvs']}" for rs in freq.index],
                       fontsize=8)
    for i in range(m.shape[0]):
        for j in range(m.shape[1]):
            ax.text(j, i, f"{m[i,j]:.2f}", ha="center", va="center",
                    fontsize=7.5,
                    color="white" if m[i, j] > m.max() / 2 else "#333")
    cb = ax.figure.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cb.set_label("allele frequency (%)", fontsize=9)
    ax.set_xlabel("ancestry")
    ax.set_title("A", loc="left", fontsize=16)


def panel_equity(ax, cov, var):
    pops = [p for p in ["AFR", "AMR", "EAS", "EUR", "SAS"] if p in cov.index]
    c = cov.loc[pops]
    out = var[~var["in_eurocentric_panel"].astype(bool)]
    x = np.arange(len(pops))
    ax.bar(x, c["captados_euro"], color=CAPT, label="captured by standard panel")
    ax.bar(x, c["perdidos_euro"], bottom=c["captados_euro"], color=MISS,
           label="missed by standard panel")
    for i, p in enumerate(pops):
        tot = int(c.loc[p, "at_risk"]); miss = int(c.loc[p, "perdidos_euro"])
        cap = int(c.loc[p, "captados_euro"]); pct = c.loc[p, "pct_perdidos"]
        if tot > 0 and pd.notna(pct):
            ax.text(i, tot + 0.7, f"{pct:.0f}% missed", ha="center",
                    fontweight="bold", fontsize=9,
                    color=MISS if pct > 0 else "#444")
        if miss >= 5 and len(out):
            ax.text(i, cap + miss / 2, out["hgvs"].iloc[0], ha="center",
                    va="center", color="white", fontsize=8.5, fontweight="bold")
    ax.set_xticks(x); ax.set_xticklabels(pops)
    ax.set_ylim(0, max(c["at_risk"].max(), 1) * 1.28)
    ax.set_ylabel("individuals with reduced DPYD metabolism")
    ax.set_xlabel("superpopulation (1000 Genomes)")
    ax.legend(frameon=False, loc="upper center", fontsize=9)
    ax.set_title("B", loc="left", fontsize=16)


def main():
    freq, cov, var, brazil = load()
    fig = plt.figure(figsize=(15, 6))
    gs = GridSpec(1, 2, figure=fig, width_ratios=[1.15, 1], wspace=0.28)
    panel_freq(fig.add_subplot(gs[0, 0]), freq, var)
    panel_equity(fig.add_subplot(gs[0, 1]), cov, var)

    fig.suptitle(
        "DPYD and 5-FU toxicity: the standard panel does not protect everyone",
        fontsize=15, y=1.03)
    src = "1000 Genomes (GRCh38)"
    if brazil:
        src += " and ABraOM SABE-1171 (Brazil)"
    out = var[~var["in_eurocentric_panel"].astype(bool)]
    out_label = ", ".join(f"{r['hgvs']} ({rs})" for rs, r in out.iterrows())
    fig.text(0.5, -0.04,
             "A) Actionable DPYD variant frequency by ancestry. "
             "B) Risk carriers missed by the standard (Eurocentric) panel; "
             f"missed carriers are driven by the out-of-panel variant {out_label}. "
             f"Data: {src}; CPIC variants. Exploratory analysis.",
             ha="center", fontsize=9, color="#444", wrap=True)
    fig.savefig(C.FIGURES / "figure_composite.png", dpi=200,
                bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(">> figure_composite.png salvo"
          + (" (com coluna Brasil)" if brazil else " (sem Brasil ainda)"))


if __name__ == "__main__":
    main()

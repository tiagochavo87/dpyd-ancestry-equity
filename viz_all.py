"""
viz_all.py
Generates ALL figure options at once into figures/ (English labels).
Reads CSVs from results/. Uses only matplotlib/numpy/pandas.

Frequency:  freq_heatmap, freq_bubble, freq_slope, freq_forest,
            freq_radar, freq_radial, freq_stream, freq_panel_highlight
Equity:     equity_waffle, equity_lollipop, equity_bar (hero)
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import config as C

POPS = ["AFR", "AMR", "EAS", "EUR", "SAS"]
SUPERPOPS = ["AFR", "AMR", "EAS", "EUR", "SAS"]  # fixo (1000G), sem BR
CAPT, MISS = "#0b7a75", "#d1495b"
IN_PANEL, OUT_PANEL = "#33576b", "#d1495b"
FUNC_COLOR = {"No function": "#5c1a33", "Decreased function": "#d1495b"}

plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11,
                     "axes.titleweight": "bold",
                     "axes.spines.top": False, "axes.spines.right": False})


def load():
    global POPS
    freq = pd.read_csv(C.ALLELE_FREQ_CSV, index_col=0)
    var = pd.read_csv(C.VARIANTS_CSV, index_col=0)
    cov = pd.read_csv(C.PANEL_COVERAGE_CSV, index_col=0)
    try:
        ci = pd.read_csv(C.RESULTS / "allele_freq_ci.csv")
    except Exception:
        ci = None
    try:  # Brasil (ABraOM), se a etapa 05 tiver rodado
        br = pd.read_csv(C.BRAZIL_FREQ_CSV, index_col=0)
        freq = freq.join(br["BR"])
        if "BR" not in POPS:
            POPS = POPS + ["BR"]
    except Exception:
        pass
    order = (freq[POPS].max(axis=1)).sort_values().index
    return freq.loc[order], var.reindex(order), cov, ci


def save(fig, path):
    fig.savefig(path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def lab(var, rs):
    return f"{rs}\n{var.loc[rs,'hgvs']}"


def empty_pops(freq):
    return [p for p in POPS if freq[p].sum() == 0]


# ---------- FREQUENCY ----------
def freq_heatmap(freq, var, path):
    m = freq[POPS].values * 100
    fig, ax = plt.subplots(figsize=(8.5, 5.2))
    im = ax.imshow(m, cmap="Reds", aspect="auto")
    ax.set_xticks(range(len(POPS))); ax.set_xticklabels(POPS)
    ax.set_yticks(range(len(freq)))
    ax.set_yticklabels([lab(var, rs) for rs in freq.index], fontsize=8)
    for i in range(m.shape[0]):
        for j in range(m.shape[1]):
            ax.text(j, i, f"{m[i,j]:.2f}", ha="center", va="center",
                    color="white" if m[i, j] > m.max() / 2 else "#333",
                    fontsize=8)
    fig.colorbar(im, ax=ax, label="allele frequency (%)")
    ax.set_title("DPYD variant frequency by ancestry")
    save(fig, path)


def freq_bubble(freq, var, path):
    rsids = list(freq.index)
    fig, ax = plt.subplots(figsize=(8.5, 5.5))
    for yi, rs in enumerate(rsids):
        for xi, p in enumerate(POPS):
            f = freq.loc[rs, p] * 100
            if f > 0:
                ax.scatter(xi, yi, s=f * 130 + 25,
                           color=FUNC_COLOR[var.loc[rs, "function"]],
                           alpha=0.85, edgecolors="white", linewidths=1)
    ax.set_xticks(range(len(POPS))); ax.set_xticklabels(POPS)
    ax.set_yticks(range(len(rsids)))
    ax.set_yticklabels([lab(var, rs) for rs in rsids], fontsize=8)
    ax.set_xlim(-0.5, len(POPS) - 0.5); ax.set_ylim(-0.6, len(rsids) - 0.4)
    ax.set_title("DPYD variant frequency by ancestry (size = frequency)")
    for p in empty_pops(freq):
        xi = POPS.index(p)
        ax.axvspan(xi - 0.45, xi + 0.45, color="#f4f4f4", zorder=0)
        ax.text(xi, len(rsids) - 0.5, "no actionable\nvariant detected",
                ha="center", va="top", fontsize=7, color="#999", style="italic")
    for f, l in [(0.5, "0.5%"), (1.0, "1%"), (2.0, "2%")]:
        ax.scatter([], [], s=f * 130 + 25, color="#999", label=l)
    ax.legend(title="frequency", loc="upper right", frameon=False,
              labelspacing=1.4)
    save(fig, path)


def freq_panel_highlight(freq, var, path):
    rsids = list(freq.index)
    fig, ax = plt.subplots(figsize=(8.8, 5.5))
    for yi, rs in enumerate(rsids):
        inpanel = bool(var.loc[rs, "in_eurocentric_panel"])
        color = IN_PANEL if inpanel else OUT_PANEL
        for xi, p in enumerate(POPS):
            f = freq.loc[rs, p] * 100
            if f > 0:
                ax.scatter(xi, yi, s=f * 140 + 25, color=color, alpha=0.85,
                           edgecolors="white", linewidths=1)
        if not inpanel:
            ax.axhspan(yi - 0.45, yi + 0.45, color=OUT_PANEL, alpha=0.07)
            ax.text(len(POPS) - 0.4, yi, "outside panel", color=OUT_PANEL,
                    fontsize=9, fontweight="bold", va="center")
    ax.set_xticks(range(len(POPS))); ax.set_xticklabels(POPS)
    ax.set_yticks(range(len(rsids)))
    ax.set_yticklabels([lab(var, rs) for rs in rsids], fontsize=8)
    ax.set_xlim(-0.5, len(POPS) + 0.8); ax.set_ylim(-0.6, len(rsids) - 0.4)
    ax.set_title("What the standard (Eurocentric) panel misses")
    for p in empty_pops(freq):
        xi = POPS.index(p)
        ax.axvspan(xi - 0.45, xi + 0.45, color="#f4f4f4", zorder=0)
        ax.text(xi, len(rsids) - 0.5, "blind spot:\nnothing detected",
                ha="center", va="top", fontsize=7, color="#888", style="italic")
    h = [plt.Line2D([], [], marker="o", ls="", color=IN_PANEL,
                    label="in standard panel"),
         plt.Line2D([], [], marker="o", ls="", color=OUT_PANEL,
                    label="outside panel (escapes)")]
    ax.legend(handles=h, loc="upper right", frameon=False, fontsize=9)
    save(fig, path)


def freq_slope(freq, var, path):
    fig, ax = plt.subplots(figsize=(8.5, 6))
    x = range(len(POPS))
    cmap = plt.cm.tab10(np.linspace(0, 1, len(var)))
    for color, rs in zip(cmap, freq.index):
        y = [freq.loc[rs, p] * 100 for p in POPS]
        lw = 3.2 if not var.loc[rs, "in_eurocentric_panel"] else 1.8
        ax.plot(x, y, marker="o", lw=lw, color=color,
                label=lab(var, rs).replace("\n", " "))
    ax.set_xticks(list(x)); ax.set_xticklabels(POPS)
    ax.set_ylabel("allele frequency (%)")
    ax.set_title("Ancestry flip across DPYD variants")
    ax.legend(frameon=False, fontsize=8)
    save(fig, path)


def freq_forest(ci, var, path):
    if ci is None:
        return
    rsids = list(var.index)
    fig, axes = plt.subplots(1, len(rsids), figsize=(3 * len(rsids), 4.8),
                             sharex=True)
    for ax, rs in zip(axes, rsids):
        sub = ci[ci.rsid == rs].set_index("super_pop").reindex(SUPERPOPS)
        y = range(len(SUPERPOPS))
        f = np.nan_to_num(sub["freq"].values * 100)
        lo = np.clip((sub["freq"] - sub["ci_low"]).values * 100, 0, None)
        hi = np.clip((sub["ci_high"] - sub["freq"]).values * 100, 0, None)
        lo = np.nan_to_num(lo); hi = np.nan_to_num(hi)
        ax.errorbar(f, list(y), xerr=[lo, hi], fmt="o", color="#33576b",
                    ecolor="#aaa", capsize=3)
        ax.set_yticks(list(y))
        ax.set_yticklabels(SUPERPOPS if ax is axes[0] else [])
        ax.set_title(f"{rs}\n{var.loc[rs,'hgvs']}", fontsize=9)
        ax.set_xlabel("freq (%)")
    fig.suptitle("Allele frequency with 95% CI (Wilson) by ancestry", y=1.02)
    save(fig, path)


def freq_radar(freq, var, path):
    angles = np.linspace(0, 2 * np.pi, len(POPS), endpoint=False)
    ang = np.concatenate([angles, [angles[0]]])
    fig = plt.figure(figsize=(7.5, 7.5))
    ax = plt.subplot(111, polar=True)
    cmap = plt.cm.tab10(np.linspace(0, 1, len(var)))
    for c, rs in zip(cmap, freq.index):
        vals = [freq.loc[rs, p] * 100 for p in POPS]; vals += [vals[0]]
        lw = 3 if not var.loc[rs, "in_eurocentric_panel"] else 1.8
        ax.plot(ang, vals, color=c, lw=lw, label=var.loc[rs, "hgvs"])
        ax.fill(ang, vals, color=c, alpha=0.07)
    ax.set_xticks(angles); ax.set_xticklabels(POPS, fontweight="bold")
    ax.set_title("Ancestry signature (radius = freq %)", pad=25)
    ax.legend(loc="upper right", bbox_to_anchor=(1.32, 1.12),
              frameon=False, fontsize=8)
    save(fig, path)


def freq_radial(freq, var, path):
    nv = len(var)
    fig = plt.figure(figsize=(8.5, 8.5))
    ax = plt.subplot(111, polar=True)
    group = 2 * np.pi / len(POPS)
    width = group / (nv + 1.5)
    cmap = plt.cm.tab10(np.linspace(0, 1, nv))
    for pi, p in enumerate(POPS):
        for vi, rs in enumerate(freq.index):
            ax.bar(pi * group + (vi + 1) * width, freq.loc[rs, p] * 100,
                   width=width * 0.9, color=cmap[vi], edgecolor="white",
                   linewidth=0.3)
    ax.set_xticks([pi * group + group / 2 for pi in range(len(POPS))])
    ax.set_xticklabels(POPS, fontweight="bold")
    ax.set_theta_zero_location("N"); ax.set_theta_direction(-1)
    h = [plt.Rectangle((0, 0), 1, 1, color=cmap[i]) for i in range(nv)]
    ax.legend(h, [var.iloc[i]["hgvs"] for i in range(nv)], loc="upper right",
              bbox_to_anchor=(1.35, 1.1), frameon=False, fontsize=8)
    ax.set_title("DPYD frequencies by ancestry (radial)", pad=25)
    save(fig, path)


def freq_stream(freq, var, path):
    x = np.arange(len(POPS))
    data = np.array([[freq.loc[rs, p] * 100 for p in POPS] for rs in freq.index])
    baseline = -data.sum(0) / 2
    xs = np.linspace(0, len(POPS) - 1, 300)
    fig, ax = plt.subplots(figsize=(9.5, 5.5))
    cmap = plt.cm.viridis(np.linspace(0.1, 0.9, len(var)))
    cum = baseline.copy()
    for i, rs in enumerate(freq.index):
        lower = np.interp(xs, x, cum); cum = cum + data[i]
        upper = np.interp(xs, x, cum)
        ax.fill_between(xs, lower, upper, color=cmap[i], alpha=0.92,
                        label=var.loc[rs, "hgvs"], linewidth=0)
    ax.set_xticks(x); ax.set_xticklabels(POPS, fontweight="bold")
    ax.set_yticks([])
    for s in ["left", "right", "top"]:
        ax.spines[s].set_visible(False)
    ax.set_title("River of DPYD frequencies (pinches to zero at EAS)")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.08),
              ncol=len(var), frameon=False, fontsize=8)
    save(fig, path)


# ---------- EQUITY ----------
def equity_waffle(cov, path):
    pops = [p for p in ["AFR", "AMR", "EUR", "SAS"] if cov.loc[p, "at_risk"] > 0]
    ncols = 8
    fig, axes = plt.subplots(1, len(pops), figsize=(3.4 * len(pops), 5))
    if len(pops) == 1:
        axes = [axes]
    maxrows = max(int(np.ceil(cov.loc[p, "at_risk"] / ncols)) for p in pops)
    for ax, p in zip(axes, pops):
        at = int(cov.loc[p, "at_risk"]); missed = int(cov.loc[p, "perdidos_euro"])
        cols = [CAPT] * (at - missed) + [MISS] * missed
        for i, c in enumerate(cols):
            r, cc = divmod(i, ncols)
            ax.add_patch(plt.Rectangle((cc, -r), 0.88, 0.88, color=c))
        ax.set_xlim(-0.3, ncols); ax.set_ylim(-maxrows, 1.2)
        ax.set_aspect("equal"); ax.axis("off")
        ax.set_title(f"{p}\n{missed}/{at} missed "
                     f"({cov.loc[p,'pct_perdidos']:.0f}%)", fontsize=12)
    h = [plt.Rectangle((0, 0), 1, 1, color=CAPT, label="captured by panel"),
         plt.Rectangle((0, 0), 1, 1, color=MISS, label="missed by panel")]
    fig.legend(handles=h, loc="lower center", ncol=2, frameon=False)
    fig.suptitle("Risk carriers missed by the standard (Eurocentric) DPYD panel\n"
                 "each square = 1 individual", fontsize=14, y=1.04)
    save(fig, path)


def equity_lollipop(cov, path):
    d = cov[cov["pct_perdidos"].notna()].copy().sort_values("pct_perdidos")
    pops = list(d.index); vals = d["pct_perdidos"].values
    fig, ax = plt.subplots(figsize=(8.5, 4.5))
    y = range(len(pops))
    ax.hlines(list(y), 0, vals, color="#cccccc", lw=3)
    ax.scatter(vals, list(y), s=220,
               color=[MISS if v > 0 else CAPT for v in vals], zorder=3)
    for yi, v in zip(y, vals):
        ax.text(v + 2.5, yi, f"{v:.0f}%", va="center", fontweight="bold")
    ax.set_yticks(list(y)); ax.set_yticklabels(pops); ax.set_xlim(0, 100)
    ax.set_xlabel("% of risk carriers missed by the standard panel")
    ax.set_title("How exposed each ancestry is")
    save(fig, path)


def equity_bar(cov, var, path):
    pops = [p for p in POPS if p in cov.index]
    c = cov.loc[pops]
    out = var[~var["in_eurocentric_panel"].astype(bool)]
    out_label = ", ".join(f"{r['hgvs']} ({rs})" for rs, r in out.iterrows())
    x = np.arange(len(pops))
    fig, ax = plt.subplots(figsize=(9, 5.8))
    ax.bar(x, c["captados_euro"], color=CAPT, label="captured by standard panel")
    ax.bar(x, c["perdidos_euro"], bottom=c["captados_euro"], color=MISS,
           label="missed by standard panel")
    for i, p in enumerate(pops):
        tot = int(c.loc[p, "at_risk"]); miss = int(c.loc[p, "perdidos_euro"])
        cap = int(c.loc[p, "captados_euro"]); pct = c.loc[p, "pct_perdidos"]
        if tot > 0 and pd.notna(pct):
            ax.text(i, tot + 0.7, f"{pct:.0f}% missed", ha="center",
                    fontweight="bold", color=MISS if pct > 0 else "#444")
        if miss >= 5:  # rotula a variante responsavel dentro do segmento vermelho
            ax.text(i, cap + miss / 2, out["hgvs"].iloc[0],
                    ha="center", va="center", color="white", fontsize=9,
                    fontweight="bold")
    ax.set_xticks(x); ax.set_xticklabels(pops)
    ax.set_ylim(0, max(c["at_risk"].max(), 1) * 1.28)
    ax.set_ylabel("Individuals with reduced DPYD metabolism")
    ax.set_xlabel("Superpopulation (1000 Genomes)")
    ax.set_title("Risk carriers missed by the standard (Eurocentric) DPYD panel")
    ax.legend(frameon=False, loc="upper center")
    ax.text(0.5, -0.16,
            f"Missed carriers are driven by an out-of-panel variant: {out_label}",
            transform=ax.transAxes, ha="center", fontsize=8.5, color="#555")
    save(fig, path)


def main():
    freq, var, cov, ci = load()
    F = C.FIGURES
    jobs = [
        ("freq_by_ancestry", lambda: freq_heatmap(freq, var, F / "freq_by_ancestry.png")),
        ("freq_bubble", lambda: freq_bubble(freq, var, F / "freq_bubble.png")),
        ("freq_panel_highlight", lambda: freq_panel_highlight(freq, var, F / "freq_panel_highlight.png")),
        ("freq_slope", lambda: freq_slope(freq, var, F / "freq_slope.png")),
        ("freq_forest", lambda: freq_forest(ci, var, F / "freq_forest.png")),
        ("freq_radar", lambda: freq_radar(freq, var, F / "freq_radar.png")),
        ("freq_radial", lambda: freq_radial(freq, var, F / "freq_radial.png")),
        ("freq_stream", lambda: freq_stream(freq, var, F / "freq_stream.png")),
        ("equity_waffle", lambda: equity_waffle(cov, F / "equity_waffle.png")),
        ("equity_lollipop", lambda: equity_lollipop(cov, F / "equity_lollipop.png")),
        ("equity_bar", lambda: equity_bar(cov, var, F / "equity_bar.png")),
    ]
    ok = []
    for name, fn in jobs:
        try:
            fn(); ok.append(name)
        except Exception as e:
            print(f"   !! {name} falhou: {e}")
    print(">> Figuras geradas em figures/:")
    for n in ok:
        print("   ", n + ".png")


if __name__ == "__main__":
    main()

"""
07_gradient_ml.py
Cruza a ancestralidade continua (06) com o status de risco (01/02) para:

  1) figura de assinatura: PCA colorido por superpopulacao (ancestry_pca.png)
  2) curva de recall do painel eurocentrico ao longo do escore de ancestralidade
     africana (recall_by_ancestry.png) -> auditoria de justica do exame
  3) modelo de vies: preve P(perdido pelo painel) a partir SO da ancestralidade
     genomica (PCs, longe do DPYD); reporta AUC (bias_model.png)

Enquadramento honesto: o "modelo" nao prediz biologia nova; ele demonstra que a
FALHA do painel e previsivel a partir da ancestralidade, que e a definicao de
vies. Analise exploratoria sobre painel de referencia, nao pacientes.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score, roc_curve
from sklearn.model_selection import cross_val_predict

import config as C

CLINE = ["AFR", "AMR", "EUR", "SAS"]   # cline AFR-EUR; EAS fica de fora
POP_COLOR = {"AFR": "#d1495b", "AMR": "#edae49", "EAS": "#66a182",
             "EUR": "#33576b", "SAS": "#8d6a9f"}
MISS = "#d1495b"


def split_haplotypes(hap):
    h1 = hap.apply(lambda c: c.map(lambda s: np.nan if s == "." else int(str(s)[0])))
    h2 = hap.apply(lambda c: c.map(lambda s: np.nan if s == "." else int(str(s)[2])))
    return h1, h2


def activity(h1, h2, var, rsids):
    samples = h1.columns
    a1 = pd.Series(1.0, index=samples); a2 = pd.Series(1.0, index=samples)
    for rs in rsids:
        if rs not in h1.index:
            continue
        act = float(var.loc[rs, "activity"])
        a1[h1.loc[rs] == 1] = np.minimum(a1[h1.loc[rs] == 1], act)
        a2[h2.loc[rs] == 1] = np.minimum(a2[h2.loc[rs] == 1], act)
    return a1 + a2


def risk_status():
    """Devolve DataFrame por amostra: at_risk (painel completo) e missed (euro)."""
    hap = pd.read_csv(C.HAPLOTYPES_CSV, index_col=0)
    var = pd.read_csv(C.VARIANTS_CSV, index_col=0)
    h1, h2 = split_haplotypes(hap)
    full = activity(h1, h2, var, list(hap.index))
    euro = activity(h1, h2, var, C.EUROCENTRIC_PANEL)
    df = pd.DataFrame({"score_full": full, "score_euro": euro})
    df["at_risk"] = df["score_full"] < 2.0
    df["missed"] = df["at_risk"] & (df["score_euro"] >= 2.0)
    df.index.name = "sample"
    return df.reset_index()


def fig_pca(anc, path):
    fig, ax = plt.subplots(figsize=(7.5, 6.5))
    for pop, g in anc.groupby("super_pop"):
        ax.scatter(g["PC1"], g["PC2"], s=14, alpha=0.7,
                   color=POP_COLOR.get(pop, "#999"), label=pop)
    ax.set_xlabel("PC1"); ax.set_ylabel("PC2")
    ax.set_title("Genetic ancestry structure (1000 Genomes)")
    ax.legend(frameon=False, title="superpopulation")
    for s in ["top", "right"]:
        ax.spines[s].set_visible(False)
    fig.savefig(path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def fig_recall(d, path):
    """Recall do painel modelado ao longo do escore de ancestralidade (logistica
    + banda bootstrap), com os pontos observados por faixa (tamanho = n)."""
    sub = d[d["super_pop"].isin(CLINE) & d["at_risk"]].copy()
    x = sub["afr_score"].values
    y = sub["missed"].astype(int).values
    if len(np.unique(y)) < 2:
        print("   sem variacao em 'missed'; curva de recall pulada.")
        return
    grid = np.linspace(0, 1, 100).reshape(-1, 1)
    base = LogisticRegression().fit(x.reshape(-1, 1), y)
    recall = (1 - base.predict_proba(grid)[:, 1]) * 100

    rng = np.random.default_rng(0)
    boots = []
    for _ in range(300):
        idx = rng.integers(0, len(x), len(x))
        if len(np.unique(y[idx])) < 2:
            continue
        try:
            c = LogisticRegression().fit(x[idx].reshape(-1, 1), y[idx])
            boots.append((1 - c.predict_proba(grid)[:, 1]) * 100)
        except Exception:
            pass
    fig, ax = plt.subplots(figsize=(8, 5))
    if boots:
        boots = np.array(boots)
        ax.fill_between(grid.ravel(), np.percentile(boots, 2.5, axis=0),
                        np.percentile(boots, 97.5, axis=0),
                        color="#33576b", alpha=0.15, label="95% CI (bootstrap)")
    ax.plot(grid.ravel(), recall, color="#33576b", lw=2.5, label="logistic fit")

    bins = np.linspace(0, 1, 6)
    sub["bin"] = pd.cut(sub["afr_score"], bins, include_lowest=True)
    obs = sub.groupby("bin").apply(lambda g: pd.Series(
        {"recall": (1 - g["missed"].mean()) * 100, "n": len(g),
         "x": g["afr_score"].mean()})).dropna()
    ax.scatter(obs["x"], obs["recall"], s=obs["n"] * 6 + 25, color=MISS,
               zorder=3, alpha=0.85, label="observed (size = n)")

    ax.set_xlabel("African-ancestry score (genome-wide, AFR<->EUR cline)")
    ax.set_ylabel("standard panel recall (%)")
    ax.set_ylim(-5, 105); ax.set_xlim(-0.03, 1.03)
    ax.set_title("The standard panel's sensitivity falls as African ancestry rises")
    ax.legend(frameon=False, fontsize=9)
    for s in ["top", "right"]:
        ax.spines[s].set_visible(False)
    fig.savefig(path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def fig_bias_model(d, path):
    """Modelo: P(perdido pelo painel) a partir dos PCs (ancestralidade genomica)."""
    m = d[d["super_pop"].isin(CLINE)].copy()
    X = m[["PC1", "PC2"]].values
    y = m["missed"].astype(int).values
    if y.sum() < 5:
        print("   poucos casos 'missed' para o modelo; pulando.")
        return None
    clf = LogisticRegression(max_iter=1000, class_weight="balanced")
    proba = cross_val_predict(clf, X, y, cv=5, method="predict_proba")[:, 1]
    auc = roc_auc_score(y, proba)
    fpr, tpr, _ = roc_curve(y, proba)

    fig, axes = plt.subplots(1, 2, figsize=(13, 5.2))
    axes[0].plot(fpr, tpr, color="#33576b", lw=2, label=f"AUC = {auc:.2f}")
    axes[0].plot([0, 1], [0, 1], "--", color="#bbb")
    axes[0].set_xlabel("false positive rate"); axes[0].set_ylabel("true positive rate")
    axes[0].set_title("Predicting 'missed by panel' from ancestry alone")
    axes[0].legend(frameon=False)
    # P(missed) previsto vs escore de ancestralidade africana
    axes[1].scatter(m["afr_score"], proba, s=12, alpha=0.5, color=MISS)
    axes[1].set_xlabel("African-ancestry score")
    axes[1].set_ylabel("predicted P(missed by panel)")
    axes[1].set_title("Predicted miss-probability rises with African ancestry")
    for ax in axes:
        for s in ["top", "right"]:
            ax.spines[s].set_visible(False)
    fig.savefig(path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return auc


def main():
    anc = pd.read_csv(C.ANCESTRY_CSV)
    risk = risk_status()
    d = anc.merge(risk, on="sample", how="inner")
    print(f">> {len(d)} individuos com ancestralidade e status de risco.")

    fig_pca(anc, C.FIGURES / "ancestry_pca.png")
    fig_recall(d, C.FIGURES / "recall_by_ancestry.png")
    auc = fig_bias_model(d, C.FIGURES / "bias_model.png")

    # resumo numerico para o post
    cl = d[d["super_pop"].isin(CLINE) & d["at_risk"]]
    lo = cl[cl["afr_score"] < 0.2]["missed"].mean()
    hi = cl[cl["afr_score"] > 0.6]["missed"].mean()
    print(">> Resumo do gradiente (entre portadores de risco do cline):")
    if pd.notna(lo):
        print(f"   recall com baixa ancestralidade africana (<0.2): {(1-lo)*100:.0f}%")
    if pd.notna(hi):
        print(f"   recall com alta ancestralidade africana  (>0.6): {(1-hi)*100:.0f}%")
    if auc is not None:
        print(f"   AUC do modelo de vies (so ancestralidade): {auc:.2f}")
    print(">> Figuras: ancestry_pca, recall_by_ancestry, bias_model.")


if __name__ == "__main__":
    main()

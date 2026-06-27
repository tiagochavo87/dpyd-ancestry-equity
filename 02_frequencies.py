"""
02_frequencies.py  (versao com escore CPIC por haplotipo)
A partir dos haplotipos fasados:
  - frequencia alelica por superpopulacao + intervalo de confianca (Wilson)
  - escore de atividade DPYD pela regra CPIC (min por haplotipo, soma dos dois)
    -> trata homo, hetero e heterozigoto composto corretamente
  - cobertura: portadores de risco perdidos pelo painel eurocentrico, por
    superpopulacao
"""
import numpy as np
import pandas as pd

import config as C


def load():
    hap = pd.read_csv(C.HAPLOTYPES_CSV, index_col=0)         # variantes x amostras
    var = pd.read_csv(C.VARIANTS_CSV, index_col=0)
    panel = pd.read_csv(C.PANEL_TSV, sep="\t").set_index("sample")
    samples = [s for s in hap.columns if s in panel.index]
    return hap[samples], var, panel.loc[samples]


def split_haplotypes(hap):
    """Devolve hap1, hap2 (variantes x amostras) com 0/1 e NaN para faltantes."""
    def first(s):
        return np.nan if s == "." else int(str(s)[0])

    def second(s):
        return np.nan if s == "." else int(str(s)[2])

    hap1 = hap.apply(lambda col: col.map(first))
    hap2 = hap.apply(lambda col: col.map(second))
    return hap1, hap2


def wilson_ci(x, n, z=1.96):
    """Intervalo de Wilson para uma proporcao (x sucessos em n)."""
    if n == 0:
        return (np.nan, np.nan)
    p = x / n
    d = 1 + z**2 / n
    centre = (p + z**2 / (2 * n)) / d
    half = (z * np.sqrt(p * (1 - p) / n + z**2 / (4 * n**2))) / d
    return (max(0.0, centre - half), min(1.0, centre + half))


def allele_frequencies(hap1, hap2, panel):
    sp = panel["super_pop"]
    dosage = hap1.fillna(0) + hap2.fillna(0)
    rows = []
    for rsid in dosage.index:
        for pop in sorted(sp.unique()) + ["ALL"]:
            cols = sp.index if pop == "ALL" else sp.index[sp == pop]
            alt = int(dosage.loc[rsid, cols].sum())
            n = 2 * len(cols)
            lo, hi = wilson_ci(alt, n)
            rows.append({"rsid": rsid, "super_pop": pop, "freq": alt / n,
                         "ci_low": lo, "ci_high": hi, "n_alleles": n})
    freq = pd.DataFrame(rows)
    freq.to_csv(C.RESULTS / "allele_freq_ci.csv", index=False)   # long + IC
    # versao larga (so freq) para o heatmap das figuras
    wide = freq.pivot(index="rsid", columns="super_pop", values="freq")
    wide.to_csv(C.ALLELE_FREQ_CSV)
    return wide


def activity_score(hap1, hap2, var, rsids):
    """Escore CPIC: cada haplotipo recebe o MENOR valor de atividade entre as
    variantes que ele carrega (1.0 se nenhuma); o escore e a soma dos dois."""
    samples = hap1.columns
    h1 = pd.Series(1.0, index=samples)
    h2 = pd.Series(1.0, index=samples)
    for rsid in rsids:
        if rsid not in hap1.index:
            continue
        a = float(var.loc[rsid, "activity"])
        on1 = hap1.loc[rsid] == 1
        on2 = hap2.loc[rsid] == 1
        h1[on1] = np.minimum(h1[on1], a)
        h2[on2] = np.minimum(h2[on2], a)
    return h1 + h2


def phenotype(score):
    # CPIC: 2 = Normal; 1 a 1.5 = Intermediario; < 1 = Pobre
    return pd.cut(score, bins=[-0.1, 0.99, 1.99, 2.01],
                  labels=["Metabolizador pobre", "Metabolizador intermediario",
                          "Metabolizador normal"])


def main():
    hap, var, panel = load()
    hap1, hap2 = split_haplotypes(hap)

    print(">> Frequencias alelicas (com IC de Wilson):")
    wide = allele_frequencies(hap1, hap2, panel)
    pops_present = [p for p in ["AFR", "AMR", "EAS", "EUR", "SAS"]
                    if p in wide.columns]
    print((wide[pops_present] * 100).round(3).to_string())

    all_rsids = list(hap.index)
    score_full = activity_score(hap1, hap2, var, all_rsids)
    score_euro = activity_score(hap1, hap2, var, C.EUROCENTRIC_PANEL)

    df = pd.DataFrame({"super_pop": panel["super_pop"],
                       "score_full": score_full, "score_euro": score_euro})
    df["pheno_full"] = phenotype(df["score_full"])
    df["at_risk_full"] = df["score_full"] < 2.0
    df["at_risk_euro"] = df["score_euro"] < 2.0
    df["missed_by_euro"] = df["at_risk_full"] & (~df["at_risk_euro"])

    (df.groupby("super_pop")["pheno_full"]
       .value_counts(normalize=True).unstack().fillna(0)).to_csv(C.PHENOTYPE_CSV)

    cov = df.groupby("super_pop").agg(
        n=("at_risk_full", "size"),
        at_risk=("at_risk_full", "sum"),
        captados_euro=("at_risk_euro", "sum"),
        perdidos_euro=("missed_by_euro", "sum"),
    )
    cov["pct_perdidos"] = (cov["perdidos_euro"] /
                           cov["at_risk"].replace(0, np.nan) * 100).round(1)
    cov.to_csv(C.PANEL_COVERAGE_CSV)

    print("\n>> Cobertura do painel eurocentrico (portadores perdidos):")
    print(cov.to_string())
    print("\n>> Etapa 2 concluida (escore por haplotipo).")


if __name__ == "__main__":
    main()

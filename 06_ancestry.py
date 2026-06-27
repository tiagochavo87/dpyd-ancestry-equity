"""
06_ancestry.py
Inferencia de ancestralidade continua a partir do 1000 Genomes.

1) Extrai SNPs comuns (MAF>=ANC_MIN_MAF), espacados, de uma regiao de chr1
   LONGE do DPYD (para nao vazar o locus do gene no escore de ancestralidade).
2) PCA (numpy SVD) sobre o genotipo padronizado.
3) Eixo de ancestralidade africana continuo por individuo: projecao na direcao
   media(AFR) - media(nao-AFR) no espaco de genotipos, normalizada para [0,1].

Saidas: results/ancestry.csv (sample, super_pop, PC1, PC2, afr_score).
Observacao honesta: e um escore supervisionado baseado em PCA sobre um painel
de referencia, nao uma proporcao formal de ADMIXTURE (que pode ser adicionada).
"""
import sys

import numpy as np
import pandas as pd

import config as C


def build_matrix():
    """Extrai genotipos (dosagem 0/1/2) de SNPs comuns espacados.
    Cacheia a matriz para reruns instantaneos; filtra amostras do painel."""
    panel = pd.read_csv(C.PANEL_TSV, sep="\t").set_index("sample")
    cache = C.RESULTS / "ancestry_matrix.npz"

    if cache.exists():
        d = np.load(cache, allow_pickle=True)
        G = d["G"]; samples = [str(s) for s in d["samples"]]
        print(f">> Matriz carregada do cache ({G.shape[0]} x {G.shape[1]}).")
    else:
        import pysam
        vcf = pysam.VariantFile(C.VCF_URL_CHR1)
        samples = list(vcf.header.samples)
        chrom, start, end = C.ANC_CONTIG_REGION
        contig = None
        for c in (str(chrom), f"chr{chrom}"):
            if c in set(vcf.header.contigs):
                contig = c; break
        if contig is None:
            print("!! contig de ancestralidade nao encontrado", file=sys.stderr)
            sys.exit(1)

        cols, positions, last_pos = [], [], -10**9
        print(f">> Extraindo SNPs de {contig}:{start}-{end} ...")
        for rec in vcf.fetch(contig, start, end):
            if len(cols) >= C.ANC_N_SNPS:
                break
            if rec.alts is None or len(rec.alts) != 1:
                continue
            if len(rec.ref) != 1 or len(rec.alts[0]) != 1:
                continue  # so SNVs bialelicos
            if rec.pos - last_pos < C.ANC_MIN_SPACING:
                continue
            dos = []
            for s in samples:
                gt = rec.samples[s].get("GT")
                if gt is None or any(a is None for a in gt):
                    dos.append(np.nan)
                else:
                    dos.append(sum(1 for a in gt if a and a > 0))
            dos = np.array(dos, dtype=float)
            maf = np.nanmean(dos) / 2.0
            maf = min(maf, 1 - maf)
            if maf < C.ANC_MIN_MAF:
                continue
            cols.append(dos); positions.append(rec.pos); last_pos = rec.pos
            if len(cols) % 250 == 0:
                print(f"   {len(cols)} SNPs...")

        G = np.array(cols).T  # samples x SNPs
        col_mean = np.nanmean(G, axis=0)
        idx = np.where(np.isnan(G))
        G[idx] = np.take(col_mean, idx[1])
        np.savez(cache, G=G, samples=np.array(samples),
                 positions=np.array(positions))
        print(f">> Matriz extraida e cacheada ({G.shape[0]} x {G.shape[1]}).")

    # filtra para amostras presentes no painel (VCF traz amostras extras)
    keep = [i for i, s in enumerate(samples) if s in panel.index]
    G = G[keep]
    samples = [samples[i] for i in keep]
    superpop = panel.loc[samples, "super_pop"].values
    print(f">> {len(samples)} individuos apos cruzar com o painel.")
    return G, samples, superpop


def standardize(G):
    p = G.mean(0) / 2.0
    sd = np.sqrt(2 * p * (1 - p))
    sd[sd == 0] = 1.0
    return (G - G.mean(0)) / sd


def pca(Z, k=2):
    U, S, _ = np.linalg.svd(Z - Z.mean(0), full_matrices=False)
    return U[:, :k] * S[:k]


def african_axis(Z, superpop):
    """Eixo de ancestralidade africana ao longo do cline AFR<->EUR (polos limpos).
    EUR -> 0, AFR -> 1; miscigenados caem no intermediario. (EAS fica fora do
    cline e e tratado a parte na analise.)"""
    afr = superpop == "AFR"
    eur = superpop == "EUR"
    d = Z[afr].mean(0) - Z[eur].mean(0)
    raw = Z @ d
    lo, hi = raw[eur].mean(), raw[afr].mean()
    score = (raw - lo) / (hi - lo) if hi != lo else raw * 0
    return np.clip(score, 0, 1)


def main():
    G, samples, superpop = build_matrix()
    Z = standardize(G)
    pcs = pca(Z, 2)
    afr = african_axis(Z, superpop)
    out = pd.DataFrame({"sample": samples, "super_pop": superpop,
                        "PC1": pcs[:, 0], "PC2": pcs[:, 1], "afr_score": afr})
    out.to_csv(C.ANCESTRY_CSV, index=False)
    print(">> ancestry.csv salvo. Media do escore africano por grupo:")
    print(out.groupby("super_pop")["afr_score"].mean().round(3).to_string())


if __name__ == "__main__":
    main()

"""
01_download.py
1) Baixa o mapa amostra->populacao do 1000 Genomes.
2) Resolve cada rsID do painel para a posicao GRCh38 via API do Ensembl
   (nada de coordenada hardcoded; auto-verificavel).
3) Extrai do VCF remoto, casando por POSICAO (a coluna de ID do VCF do 1000G
   costuma nao trazer rsID), e salva a dosagem do alelo alternativo por amostra.
"""
import sys

import numpy as np
import pandas as pd
import requests

import config as C

ENSEMBL = "https://rest.ensembl.org/variation/human/{rsid}?content-type=application/json"


def download_panel():
    print(">> Baixando mapa de populacoes...")
    r = requests.get(C.PANEL_URL, timeout=120)
    r.raise_for_status()
    C.PANEL_TSV.write_text(r.text)
    panel = pd.read_csv(C.PANEL_TSV, sep="\t")
    print(f"   {len(panel)} amostras, superpops: "
          f"{sorted(panel['super_pop'].dropna().unique())}")


def resolve_position(rsid):
    """rsID -> (chrom, pos GRCh38, allele_string) via Ensembl REST."""
    r = requests.get(ENSEMBL.format(rsid=rsid), timeout=60)
    r.raise_for_status()
    for m in r.json().get("mappings", []):
        if m.get("assembly_name") == "GRCh38":
            return str(m["seq_region_name"]), int(m["start"]), m.get("allele_string")
    return None


def pick_contig(vcf, chrom):
    contigs = set(vcf.header.contigs)
    for c in (f"chr{chrom}", chrom, C.CONTIG, "chr1", "1"):
        if c in contigs:
            return c
    raise ValueError(f"Contig nao encontrado. VCF tem: {list(contigs)[:5]}")


def fetch_variants():
    import pysam

    print(">> Abrindo VCF remoto do 1000 Genomes...")
    vcf = pysam.VariantFile(C.VCF_URL_CHR1)
    samples = list(vcf.header.samples)

    rows, meta = {}, []
    for v in C.VARIANTS:
        rsid = v["rsid"]
        loc = resolve_position(rsid)
        if loc is None:
            print(f"   !! {rsid}: posicao GRCh38 nao resolvida, pulando")
            continue
        chrom, pos, alleles = loc
        contig = pick_contig(vcf, chrom)
        rec_found = None
        for rec in vcf.fetch(contig, pos - 1, pos + 1):
            if rec.pos == pos:
                rec_found = rec
                break
        if rec_found is None:
            print(f"   !! {rsid} ({chrom}:{pos}) nao encontrado no VCF")
            continue

        haps = []
        for sm in samples:
            gt = rec_found.samples[sm].get("GT")
            if gt is None or len(gt) < 2 or any(a is None for a in gt):
                haps.append(".")
            else:
                a = 1 if (gt[0] and gt[0] > 0) else 0
                b = 1 if (gt[1] and gt[1] > 0) else 0
                haps.append(f"{a}|{b}")
        rows[rsid] = haps
        meta.append({
            "rsid": rsid, "chrom": rec_found.chrom, "pos": rec_found.pos,
            "ref": rec_found.ref, "alt": ",".join(rec_found.alts or []),
            "ensembl_alleles": alleles,
        })
        print(f"   ok {rsid} -> {contig}:{pos} "
              f"({rec_found.ref}>{','.join(rec_found.alts or [])})")

    if not rows:
        print("!! Nenhuma variante extraida. Verifique conexao/URL.",
              file=sys.stderr)
        sys.exit(1)

    hap = pd.DataFrame(rows, index=samples).T
    hap.to_csv(C.HAPLOTYPES_CSV)

    vmeta = pd.DataFrame(meta).set_index("rsid")
    cfg = pd.DataFrame(C.VARIANTS).set_index("rsid")
    vmeta = vmeta.join(cfg[["star", "hgvs", "function", "activity"]])
    vmeta["in_eurocentric_panel"] = vmeta.index.isin(C.EUROCENTRIC_PANEL)
    vmeta.to_csv(C.VARIANTS_CSV)
    print(f">> {len(hap)} variantes x {hap.shape[1]} amostras salvas (haplotipos fasados).")


def main():
    download_panel()
    fetch_variants()
    print(">> Etapa 1 concluida.")


if __name__ == "__main__":
    main()

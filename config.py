"""
Config do projeto PGx + ancestralidade (foco DPYD / toxicidade ao 5-FU).
O painel de variantes fica aqui: para estender a outro gene (NUDT15, TPMT),
basta trocar GENE, REGION e VARIANTS.
"""
from pathlib import Path

GENE = "DPYD"

# --- Regiao GRCh38 do gene (janela generosa; variantes casadas por rsID) -----
# DPYD: cromossomo 1 (1p21.3), fita negativa, ~950 kb.
CONTIG = "chr1"            # alguns VCFs usam "1"; o codigo tenta os dois
REGION_START = 97_400_000
REGION_END = 98_600_000

# --- Fontes publicas (1000 Genomes, GRCh38, IGSR/EBI) ------------------------
# VCF bialelico, fasado, GRCh38, 2548 amostras, indexado com tabix (.tbi)
VCF_URL_CHR1 = (
    "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/"
    "1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/"
    "ALL.chr1.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz"
)
# Mapa amostra -> populacao / superpopulacao (sample, pop, super_pop, gender)
PANEL_URL = (
    "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/"
    "integrated_call_samples_v3.20130502.ALL.panel"
)

# --- Painel de variantes acionaveis do DPYD (CPIC, conjunto ampliado) ---------
# activity: valor de atividade do alelo variante (CPIC). Normal = 1.0.
# No function = 0.0 ; Decreased function = 0.5.
VARIANTS = [
    {"rsid": "rs3918290",   "star": "DPYD*2A", "hgvs": "c.1905+1G>A",
     "function": "No function",       "activity": 0.0},
    {"rsid": "rs55886062",  "star": "DPYD*13", "hgvs": "c.1679T>G",
     "function": "No function",       "activity": 0.0},
    {"rsid": "rs72549303",  "star": "DPYD*3",  "hgvs": "c.1898del",
     "function": "No function",       "activity": 0.0},
    {"rsid": "rs78060119",  "star": "DPYD*12", "hgvs": "c.1156G>T",
     "function": "No function",       "activity": 0.0},
    {"rsid": "rs67376798",  "star": "-",       "hgvs": "c.2846A>T",
     "function": "Decreased function", "activity": 0.5},
    {"rsid": "rs75017182",  "star": "HapB3",   "hgvs": "c.1129-5923C>G",
     "function": "Decreased function", "activity": 0.5},
    {"rsid": "rs115232898", "star": "-",       "hgvs": "c.557A>G",
     "function": "Decreased function", "activity": 0.5},
]

# "Painel eurocentrico" = as quatro variantes tier-1 do teste padrao.
EUROCENTRIC_PANEL = ["rs3918290", "rs55886062", "rs67376798", "rs75017182"]

# --- Caminhos ----------------------------------------------------------------
ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results"
FIGURES = ROOT / "figures"
for _p in (RESULTS, FIGURES):
    _p.mkdir(parents=True, exist_ok=True)

PANEL_TSV = RESULTS / "samples_panel.tsv"
HAPLOTYPES_CSV = RESULTS / "dpyd_haplotypes.csv"
VARIANTS_CSV = RESULTS / "dpyd_variants.csv"
ALLELE_FREQ_CSV = RESULTS / "allele_freq_by_superpop.csv"
PHENOTYPE_CSV = RESULTS / "phenotype_by_superpop.csv"
PANEL_COVERAGE_CSV = RESULTS / "panel_coverage_by_superpop.csv"

# --- Fonte brasileira: ABraOM (SABE-1171-WGS, hg38, USP) ---------------------
# Sem API aberta: baixe o VCF no formulario (uso academico, citar Naslavsky et
# al.) e coloque indexado (.tbi) em data/. O pipeline extrai a regiao do DPYD.
# Pagina: https://abraom.ib.usp.br/download
DATA = ROOT / "data"
DATA.mkdir(parents=True, exist_ok=True)
ABRAOM_TSV_TARGZ = DATA / "SABE1171.Abraom.clean.tar.gz"  # TSV anotado, hg38
BRAZIL_FREQ_CSV = RESULTS / "brazil_freq.csv"

# --- Ancestralidade (PCA + eixo africano) ------------------------------------
ANCESTRY_CSV = RESULTS / "ancestry.csv"
# SNPs comuns para PCA, amostrados de uma regiao de chr1 LONGE do DPYD
# (DPYD ~ 97-98 Mb; usamos uma janela distante para nao vazar o locus)
ANC_CONTIG_REGION = (1, 150000000, 200000000)  # chrom, start, end (GRCh38)
ANC_N_SNPS = 2500          # alvo de SNPs
ANC_MIN_MAF = 0.05         # apenas SNPs comuns
ANC_MIN_SPACING = 15000    # espacamento minimo (bp) para reduzir LD

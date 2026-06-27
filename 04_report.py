"""
04_report.py
Compila todos os resultados num unico relatorio de texto: results/RESULTADOS.txt
Le os CSVs gerados pelas etapas anteriores. Tolerante a arquivos ausentes.
"""
from datetime import datetime

import pandas as pd

import config as C

LINE = "=" * 70
SUB = "-" * 70


def _safe(fn):
    try:
        return fn()
    except Exception as e:
        return f"[indisponivel: {e}]"


def sec_panel():
    var = pd.read_csv(C.VARIANTS_CSV)
    cols = [c for c in ["rsid", "star", "hgvs", "function", "activity",
                        "in_eurocentric_panel"] if c in var.columns]
    return var[cols].to_string(index=False)


def sec_samples():
    panel = pd.read_csv(C.PANEL_TSV, sep="\t")
    vc = panel["super_pop"].value_counts().sort_index()
    vc.loc["TOTAL"] = vc.sum()
    return vc.to_string()


def sec_freq():
    df = pd.read_csv(C.RESULTS / "allele_freq_ci.csv")
    df = df[df["super_pop"] != "ALL"].copy()
    df["freq_%"] = (df["freq"] * 100).round(3)
    df["IC95%"] = ("[" + (df["ci_low"] * 100).round(3).astype(str) + ", "
                   + (df["ci_high"] * 100).round(3).astype(str) + "]")
    out = df.pivot(index="rsid", columns="super_pop", values="freq_%")
    return out.to_string()


def sec_freq_ci():
    df = pd.read_csv(C.RESULTS / "allele_freq_ci.csv")
    df = df[df["super_pop"] != "ALL"].copy()
    df["freq_%"] = (df["freq"] * 100).round(3)
    df["IC95_%"] = ("[" + (df["ci_low"] * 100).round(2).astype(str) + " - "
                    + (df["ci_high"] * 100).round(2).astype(str) + "]")
    return df[["rsid", "super_pop", "freq_%", "IC95_%", "n_alleles"]].to_string(index=False)


def sec_pheno():
    df = pd.read_csv(C.PHENOTYPE_CSV, index_col=0)
    return (df * 100).round(1).to_string()


def sec_cov():
    return pd.read_csv(C.PANEL_COVERAGE_CSV, index_col=0).to_string()


METODOLOGIA = """\
Exploracao da cobertura do teste de DPYD (toxicidade a fluoropirimidinas/5-FU)
entre ancestralidades, usando dados publicos.

- Dados: 1000 Genomes Project, GRCh38, release bialelico fasado (IGSR/EBI).
- Variantes e funcao alelica: diretriz CPIC para DPYD/fluoropirimidinas.
- Posicoes resolvidas via Ensembl REST (GRCh38); casamento por posicao.
- Escore de atividade DPYD pela regra CPIC: cada haplotipo recebe o MENOR
  valor de atividade entre as variantes que carrega (1.0 se nenhuma); o
  escore do gene e a soma dos dois haplotipos (trata homo, hetero e
  heterozigoto composto). Fenotipo: 2=Normal; 1 a 1.5=Intermediario; <1=Pobre.
- "Painel eurocentrico" = as 4 variantes tier-1 do teste padrao.
- Frequencias com intervalo de confianca de Wilson (95%)."""

RESSALVAS = """\
- Analise EXPLORATORIA, sem fim de inferencia clinica definitiva.
- 1000 Genomes e painel de referencia populacional, NAO uma coorte de
  pacientes; o trabalho mede cobertura sobre PORTADORES de variantes, nao
  desfechos clinicos medidos de toxicidade.
- Numero de portadores por grupo pode ser pequeno (ver n_alleles e IC).
- Frequencia baixa de variantes acionaveis em alguma ancestralidade pode
  refletir lacuna de pesquisa (vies de ascertainment) e nao ausencia de risco.
- A fase dos haplotipos vem do painel fasado do 1000 Genomes."""


def main():
    parts = [
        LINE,
        " RELATORIO DE RESULTADOS - DPYD x ANCESTRALIDADE",
        f" Projeto exploratorio | gerado em {datetime.now():%Y-%m-%d %H:%M}",
        LINE, "",
        "1. ESCOPO E METODOLOGIA", SUB, METODOLOGIA, "",
        "2. PAINEL DE VARIANTES (CPIC)", SUB, _safe(sec_panel), "",
        "3. AMOSTRAS POR SUPERPOPULACAO", SUB, _safe(sec_samples), "",
        "4. FREQUENCIA ALELICA POR SUPERPOPULACAO (%)", SUB, _safe(sec_freq), "",
        "   Detalhe com IC95% (Wilson):", _safe(sec_freq_ci), "",
        "5. DISTRIBUICAO DE FENOTIPOS DPYD POR SUPERPOPULACAO (%)", SUB,
        _safe(sec_pheno), "",
        "6. COBERTURA DO PAINEL EUROCENTRICO (EQUIDADE)", SUB, _safe(sec_cov), "",
        "   Leitura: 'perdidos_euro' = portadores de risco (metabolismo",
        "   reduzido) que o painel tier-1 padrao NAO sinalizaria.", "",
        "7. RESSALVAS", SUB, RESSALVAS, "", LINE,
    ]
    out = C.RESULTS / "RESULTADOS.txt"
    out.write_text("\n".join(parts), encoding="utf-8")
    print(f">> Relatorio salvo em {out}")


if __name__ == "__main__":
    main()

"""
05_brazil.py
Adiciona a populacao brasileira (ABraOM / SABE-1171-WGS, hg38).

O download do WGS-1171 e uma TABELA TSV anotada (dentro de um .tar.gz), nao um
VCF. Esta etapa le o TSV em streaming de dentro do .tar.gz (sem extrair o
arquivo gigante), filtra as posicoes GRCh38 do painel e grava brazil_freq.csv.

Coloque o arquivo (uso academico; citar Naslavsky et al.) em:
    data/SABE1171.Abraom.clean.tar.gz
Pagina: https://abraom.ib.usp.br/download

E fonte de FREQUENCIA (sem genotipos individuais): o Brasil entra na comparacao
de frequencias, nao na contagem individuo a individuo.
"""
import io
import sys
import tarfile
from pathlib import Path

import pandas as pd

import config as C


def find_source():
    """Procura o tar.gz/tsv do ABraOM em data/ e na raiz do projeto."""
    names = ["SABE1171.Abraom.clean.tar.gz", "SABE1171.Abraom.clean.tsv"]
    cands = [C.DATA / n for n in names] + [C.ROOT / n for n in names]
    cands += list(C.DATA.glob("*SABE*tar.gz")) + list(C.ROOT.glob("*SABE*tar.gz"))
    cands += list(C.DATA.glob("*SABE*tsv")) + list(C.ROOT.glob("*SABE*tsv"))
    for c in cands:
        if c.exists():
            return c
    return None


def open_stream(src):
    """Devolve um iteravel de linhas de texto do TSV (extrai do tar.gz se preciso)."""
    if str(src).endswith(".tar.gz") or str(src).endswith(".tgz"):
        tf = tarfile.open(src, "r:gz")
        member = next((m for m in tf.getmembers() if m.name.endswith(".tsv")), None)
        if member is None:
            print("!! Nenhum .tsv dentro do tar.gz", file=sys.stderr); sys.exit(1)
        return io.TextIOWrapper(tf.extractfile(member), encoding="utf-8"), tf
    return open(src, encoding="utf-8"), None


def colmap(header):
    """Mapeia nomes de coluna do ABraOM por correspondencia flexivel."""
    cols = {h.strip().lower(): i for i, h in enumerate(header)}

    def find(*opts):
        for o in opts:
            for name, i in cols.items():
                if o in name:
                    return i
        return None

    return {
        "chr": find("chr"),
        "pos": find("start position", "start", "position", "pos"),
        "ref": find("ref"),
        "alt": find("alt"),
        "rsid": find("dbsnp", "rsid", "rs"),
        "freq": find("frequency", "freq", " af"),
        "ac": find("allele count"),
        "an": find("allele number"),
    }


def main():
    src = find_source()
    if src is None:
        print("!! ABraOM nao encontrado (data/SABE1171.Abraom.clean.tar.gz).")
        print("   Etapa do Brasil pulada; o resto do pipeline segue normal.")
        return
    print(f">> Lendo ABraOM: {src.name}")

    var = pd.read_csv(C.VARIANTS_CSV)
    targets = {int(p): rs for rs, p in zip(var["rsid"], var["pos"])}
    alt_of = {rs: str(a).split(",")[0] for rs, a in zip(var["rsid"], var["alt"])}
    found = {}

    stream, tf = open_stream(src)
    try:
        header = stream.readline().rstrip("\n").split("\t")
        cm = colmap(header)
        if cm["chr"] is None or cm["pos"] is None:
            print(f"!! Cabecalho inesperado: {header}", file=sys.stderr); sys.exit(1)
        for line in stream:
            # atalho: so cromossomo 1
            if not (line.startswith("1\t") or line.startswith("chr1\t")):
                continue
            f = line.rstrip("\n").split("\t")
            try:
                pos = int(f[cm["pos"]])
            except (ValueError, IndexError):
                continue
            if pos not in targets:
                continue
            rsid = targets[pos]
            # frequencia: usa coluna Frequency; senao AC/AN
            af = None
            if cm["freq"] is not None and f[cm["freq"]] not in ("", "NA", "."):
                af = float(f[cm["freq"]])
                if af > 1:           # veio em %, normaliza
                    af /= 100.0
            elif cm["ac"] is not None and cm["an"] is not None:
                an = float(f[cm["an"]])
                if an > 0:
                    af = float(f[cm["ac"]]) / an
            if af is not None:
                # se ja achou esta variante, mantem a de maior freq (alelo certo)
                found[rsid] = max(found.get(rsid, 0.0), af)
    finally:
        stream.close()
        if tf is not None:
            tf.close()

    rows = [{"rsid": rs, "BR": found.get(rs, 0.0)} for rs in var["rsid"]]
    out = pd.DataFrame(rows).set_index("rsid")
    out.to_csv(C.BRAZIL_FREQ_CSV)
    for rs in var["rsid"]:
        tag = f"{found[rs]*100:.3f}%" if rs in found else "ausente (freq=0)"
        print(f"   {rs}: BR = {tag}")
    print(f">> Frequencias brasileiras salvas em {C.BRAZIL_FREQ_CSV}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
discover_primers.py

Recupera as sequencias dos primers do painel QIAseq 16S/ITS diretamente dos FASTQ,
sem depender do CLC Genomics.

Como funciona
-------------
Os primers do painel sao "phased": a QIAGEN insere de 0 a 11 bases antes do primer
para gerar diversidade de bases no inicio da read. Isso quebra qualquer corte por
posicao fixa, mas cria uma assinatura estatistica muito clara: o mesmo motivo
aparece repetido em varios deslocamentos.

O script:
  1. amostra reads de varias amostras;
  2. conta k-mers em todos os deslocamentos de 0 a MAXPHASE;
  3. agrupa k-mers parecidos e monta um consenso IUPAC por grupo;
  4. descarta grupos que sao apenas versoes deslocadas de outro grupo
     (criterio: um primer de verdade aparece tambem em deslocamento 0);
  5. estende o consenso ate o comprimento tipico de um primer;
  6. cruza os motivos de R1 com os de R2 e reporta os pares F/R observados,
     que sao as regioes do painel.

Saida: tabela no stdout + primers.fasta + primers.tsv

Uso
---
  python3 discover_primers.py \
      --dir /home/thiago.parente/data.thiago.parente/2026_tematicos/dados_brutos \
      --samples 24 --reads 50000 --out ./primers_descobertos

Nao requer internet nem bibliotecas externas (apenas Python 3.7+).
"""

import argparse
import gzip
import os
import random
import re
import sys
from collections import Counter

# ---------------------------------------------------------------- IUPAC

_SETS = {
    "A": "A", "C": "C", "G": "G", "T": "T",
    "AG": "R", "CT": "Y", "CG": "S", "AT": "W", "GT": "K", "AC": "M",
    "CGT": "B", "AGT": "D", "ACT": "H", "ACG": "V", "ACGT": "N",
}
TO_IUPAC = {frozenset(k): v for k, v in _SETS.items()}
EXPAND = {v: set(k) for k, v in _SETS.items()}


def iupac_code(bases):
    """Conjunto de bases -> letra IUPAC."""
    key = frozenset(b for b in bases if b in "ACGT")
    return TO_IUPAC.get(key, "N")


def mismatches(motif, seq, cap):
    """Numero de discordancias entre um motivo IUPAC e uma sequencia.
    Interrompe assim que passa de `cap` (retorna cap + 1)."""
    n = 0
    for m, s in zip(motif, seq):
        if s not in EXPAND.get(m, {m}):
            n += 1
            if n > cap:
                return cap + 1
    return n


def hamming(a, b, cap):
    n = 0
    for x, y in zip(a, b):
        if x != y:
            n += 1
            if n > cap:
                return cap + 1
    return n


def compativel(a, b, max_mm):
    """Dois consensos IUPAC batem se, posicao a posicao, os conjuntos de
    bases se cruzam. Tolera ate `max_mm` posicoes discordantes."""
    n = 0
    for x, y in zip(a, b):
        if not (EXPAND.get(x, {x}) & EXPAND.get(y, {y})):
            n += 1
            if n > max_mm:
                return False
    return True


def mesma_familia(a, b, min_comuns=50, fracao=0.8):
    """True se `a` e `b` sao o mesmo primer visto com deslocamentos
    diferentes. E o caso patologico deste metodo: o k-mer formado pela
    ultima base do bloco de fase mais o inicio do primer tambem e
    frequente, e vira um motivo fantasma.

    O criterio e o CONJUNTO DE READS, nao a semelhanca das sequencias.
    Comparar consensos IUPAC nao funciona: um consenso cheio de N, B, D e
    V e "compativel" com quase tudo, entao familias se fundem entre
    primers diferentes e primers reais somem da lista. Ja dois motivos que
    sao o mesmo primer deslocado casam nas MESMAS reads, sempre com a
    mesma diferenca de posicao; dois primers diferentes quase nunca
    aparecem na mesma read.

    `a` e `b` sao dicionarios com a chave "hits" (posicao por read, -1 se
    nao casou).
    """
    difs, comuns = Counter(), 0
    for ha, hb in zip(a["hits"], b["hits"]):
        if ha >= 0 and hb >= 0:
            comuns += 1
            difs[hb - ha] += 1
    if comuns < min_comuns:
        return False
    return difs.most_common(1)[0][1] >= fracao * comuns


def espalhamento(hist, minimo=0.02):
    """Quantos deslocamentos distintos concentram pelo menos `minimo` das
    reads. Um primer phased de verdade se espalha por 0..11; um motivo
    ancorado numa posicao fixa e artefato — sequencia interna que por
    acaso cai sempre no mesmo lugar da read."""
    total = sum(hist.values())
    if not total:
        return 0
    return sum(1 for n in hist.values() if n >= minimo * total)


# ---------------------------------------------------------------- FASTQ

FASTQ_EXT = (".fastq.gz", ".fq.gz", ".fastq", ".fq")

# Tokens que distinguem R1 de R2, na ordem em que sao testados.
MATE_TOKENS = (("_R1_", "_R2_"), ("_R1.", "_R2."), (".R1.", ".R2."),
               ("_1.", "_2."), ("_R1", "_R2"))


def all_fastqs(root):
    """Todos os FASTQ sob `root`, em qualquer profundidade.

    Dois detalhes que fazem diferenca em cluster:
      followlinks=True  — por padrao o os.walk NAO entra em diretorio que e
                          symlink, e em cluster e comum a area de dados
                          brutos ser um monte de links para outro filesystem.
      onerror           — por padrao o os.walk engole erros de permissao em
                          silencio, e o resultado vazio parece "nao tem
                          arquivo" quando na verdade e "nao pude olhar".
    """
    out, erros, vistos = [], [], set()

    def onerror(e):
        erros.append(e)

    for dirpath, dirnames, files in os.walk(root, onerror=onerror,
                                            followlinks=True):
        real = os.path.realpath(dirpath)
        if real in vistos:          # protege contra ciclo de symlinks
            dirnames[:] = []
            continue
        vistos.add(real)
        for fn in files:
            if fn.endswith(FASTQ_EXT):
                out.append(os.path.join(dirpath, fn))
    return sorted(out), erros


def diagnostico(root):
    """Linhas de diagnostico para quando nada e encontrado."""
    linhas = []
    if not os.path.exists(root):
        return ["O caminho nao existe."]
    if not os.path.isdir(root):
        return ["O caminho existe mas nao e um diretorio."]
    try:
        entradas = sorted(os.listdir(root))
    except OSError as e:
        return ["Nao consegui listar o diretorio: %s" % e]
    links = [e for e in entradas
             if os.path.islink(os.path.join(root, e))]
    linhas.append("%d entrada(s) no primeiro nivel, %d symlink(s)."
                  % (len(entradas), len(links)))
    for e in entradas[:3]:
        p = os.path.join(root, e)
        tipo = "dir" if os.path.isdir(p) else "arquivo"
        if os.path.islink(p):
            tipo += " -> " + os.path.realpath(p)
        linhas.append("  %-28s %s" % (e, tipo))
        if os.path.isdir(p):
            try:
                dentro = sorted(os.listdir(p))[:3]
                for d in dentro:
                    linhas.append("      %s" % d)
            except OSError as err:
                linhas.append("      (sem permissao: %s)" % err)
    return linhas


def mate_of(path):
    """Caminho do par R2 e o token usado. (None, None) se nao parecer R1."""
    d, base = os.path.split(path)
    for r1tok, r2tok in MATE_TOKENS:
        if r1tok in base:
            return os.path.join(d, base.replace(r1tok, r2tok, 1)), r1tok
    return None, None


def sample_of(base, r1tok):
    """Nome da amostra a partir do nome do arquivo.
    9_S277_L001_R1_001.fastq.gz -> 9"""
    m = re.match(r"(.+?)_S\d+_L\d+_R[12]", base)
    if m:
        return m.group(1)
    return base.split(r1tok)[0].rstrip("_.")


def find_pairs(root):
    """Localiza os pares R1/R2 em qualquer profundidade sob `root`."""
    todos, erros = all_fastqs(root)
    pairs, orfaos = [], 0
    for r1 in todos:
        r2, tok = mate_of(r1)
        if r2 is None:
            continue
        if os.path.exists(r2):
            pairs.append((sample_of(os.path.basename(r1), tok), r1, r2))
        else:
            orfaos += 1
    return pairs, todos, orfaos, erros


def head_seqs(path, nreads, width):
    """Le as primeiras `nreads` sequencias, truncadas em `width` bases."""
    out = []
    with gzip.open(path, "rt") as fh:
        for i, line in enumerate(fh):
            if i & 3 == 1:
                out.append(line.rstrip("\n")[:width].upper())
                if len(out) >= nreads:
                    break
    return out


# ---------------------------------------------------------------- nucleo

def count_kmers(reads, k, maxphase):
    c = Counter()
    for seq in reads:
        for off in range(maxphase + 1):
            kmer = seq[off:off + k]
            if len(kmer) == k and "N" not in kmer:
                c[kmer] += 1
    return c


def cluster(kmers, max_dist):
    """Agrupa k-mers por distancia de Hamming, do mais frequente ao menos.
    Retorna lista de (consenso_iupac, contagem_total)."""
    clusters = []  # (representante, [(kmer, contagem)])
    for kmer, n in kmers:
        placed = False
        for rep, members in clusters:
            if hamming(rep, kmer, max_dist) <= max_dist:
                members.append((kmer, n))
                placed = True
                break
        if not placed:
            clusters.append((kmer, [(kmer, n)]))

    out = []
    for _rep, members in clusters:
        total = sum(n for _, n in members)
        cons = []
        for pos in range(len(members[0][0])):
            per_base = Counter()
            for kmer, n in members:
                per_base[kmer[pos]] += n
            keep = {b for b, n in per_base.items() if n >= 0.10 * total}
            cons.append(iupac_code(keep) if keep else "N")
        out.append(("".join(cons), total))
    out.sort(key=lambda x: -x[1])
    return out


def locate(motif, reads, maxphase, max_mm):
    """Para cada read, acha o menor deslocamento em que o motivo casa.
    Retorna (lista de deslocamentos por read com -1 para nao-casou, histograma)."""
    hits = []
    hist = Counter()
    k = len(motif)
    for seq in reads:
        found = -1
        for off in range(maxphase + 1):
            if mismatches(motif, seq[off:off + k], max_mm) <= max_mm:
                found = off
                break
        hits.append(found)
        if found >= 0:
            hist[found] += 1
    return hits, hist


def extend(motif, reads, hits, length):
    """Estende o consenso para a direita usando as reads que casaram."""
    cols = [Counter() for _ in range(length)]
    used = 0
    for seq, off in zip(reads, hits):
        if off < 0:
            continue
        frag = seq[off:off + length]
        if len(frag) < length:
            continue
        used += 1
        for i, b in enumerate(frag):
            cols[i][b] += 1
    if used == 0:
        return motif
    cons = []
    for col in cols:
        keep = {b for b, n in col.items() if n >= 0.10 * used and b in "ACGT"}
        cons.append(iupac_code(keep) if keep else "N")
    return "".join(cons)


def discover(reads, k, maxphase, max_mm, top, primer_len, min_frac, label):
    n = len(reads)
    kmers = count_kmers(reads, k, maxphase).most_common(top)
    cands = cluster(kmers, max_dist=3)

    resolved = []
    for cons, _total in cands[: max(24, top // 2)]:
        hits, hist = locate(cons, reads, maxphase, max_mm)
        matched = sum(1 for h in hits if h >= 0)
        if matched < min_frac * n:
            continue
        full = extend(cons, reads, hits, primer_len)
        # A extensao entra no molde depois que o primer acaba, e ali o
        # consenso vira N. Aparar esses N devolve o comprimento real do
        # primer em vez de um numero arbitrario.
        full = full.rstrip("N")
        if len(full) < 12:
            full = cons
        media = sum(h for h in hits if h >= 0) / float(max(1, matched))
        resolved.append({"motif": cons, "primer": full, "n": matched,
                         "frac": matched / n, "hist": hist, "hits": hits,
                         "media": media, "esp": espalhamento(hist)})

    # Colapsa cada familia de deslocamentos num representante so.
    #
    # Duas formas de fantasma aparecem, e elas se distinguem:
    #   - deslocado a direita (primer sem as primeiras bases): casa nas
    #     MESMAS reads que o primer, so que mais adiante. Mesma contagem,
    #     posicao media maior.
    #   - estendido a esquerda (bloco de fase + primer): comeca mais cedo,
    #     mas so casa nas reads que tem aquele bloco — contagem bem menor.
    # Logo: fica quem comeca mais cedo, desde que nao perca reads.
    uniq = []
    for r in sorted(resolved, key=lambda x: -x["n"]):
        irmao = next((u for u in uniq if mesma_familia(r, u)), None)
        if irmao is None:
            uniq.append(r)
        # 0.98 e proposital: o estendido a esquerda perde poucas reads, mas
        # perde. Com limiar frouxo ele venceria e o primer sairia com uma ou
        # duas bases do bloco de fase grudadas — que depois fariam o cutadapt
        # ignorar justamente as reads com outro bloco.
        elif r["n"] >= 0.98 * irmao["n"] and r["media"] < irmao["media"] - 0.5:
            uniq[uniq.index(irmao)] = r

    sys.stderr.write("[%s] %d reads, %d motivos retidos\n" % (label, n, len(uniq)))
    return uniq


# ---------------------------------------------------------------- main

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dir", required=True, help="diretorio com os dados brutos")
    ap.add_argument("--samples", type=int, default=24, help="quantas amostras amostrar")
    ap.add_argument("--reads", type=int, default=50000, help="reads por amostra")
    ap.add_argument("--out", default="primers_descobertos", help="prefixo de saida")
    ap.add_argument("--k", type=int, default=18, help="tamanho do k-mer")
    ap.add_argument("--maxphase", type=int, default=13, help="deslocamento maximo do phasing")
    ap.add_argument("--primer-len", type=int, default=22, help="comprimento do primer a reportar")
    ap.add_argument("--max-mm", type=int, default=2, help="discordancias toleradas")
    ap.add_argument("--min-frac", type=float, default=0.01,
                    help="fracao minima de reads para um motivo ser retido")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--listar", action="store_true",
                    help="so lista os pares R1/R2 encontrados e sai")
    ap.add_argument("--filtrar", metavar="REGEX",
                    help="usa so as amostras cujo nome casa com a expressao. "
                         "Ex.: --filtrar '^Smart' para rodar so nos controles, "
                         "onde o construto sintetico tem sitio para TODOS os "
                         "primers do painel e serve de gabarito.")
    args = ap.parse_args()

    raiz = os.path.realpath(args.dir)
    pairs, todos, orfaos, erros = find_pairs(raiz)
    if not pairs:
        sys.stderr.write("\nNenhum par R1/R2 encontrado em:\n  %s\n\n" % raiz)
        sys.stderr.write("%d arquivo(s) FASTQ localizado(s).\n" % len(todos))
        for p in todos[:5]:
            sys.stderr.write("  %s\n" % os.path.relpath(p, raiz))
        if len(todos) > 5:
            sys.stderr.write("  ... e mais %d\n" % (len(todos) - 5))
        if todos:
            sys.stderr.write(
                "\nOs arquivos existem mas nao formam par. Verifique se os nomes\n"
                "usam um dos padroes reconhecidos: %s\n"
                % ", ".join(a for a, _ in MATE_TOKENS))
        else:
            sys.stderr.write("\nEstrutura encontrada:\n")
            for linha in diagnostico(raiz):
                sys.stderr.write("  %s\n" % linha)
        if erros:
            sys.stderr.write("\n%d erro(s) de acesso durante a varredura, "
                             "o primeiro deles:\n  %s\n" % (len(erros), erros[0]))
        sys.exit(1)
    sys.stderr.write("Encontrados %d pares R1/R2 em %d arquivos FASTQ"
                     % (len(pairs), len(todos)))
    sys.stderr.write((", %d R1 sem par.\n" % orfaos) if orfaos else ".\n")

    if args.filtrar:
        padrao = re.compile(args.filtrar)
        antes = len(pairs)
        pairs = [p for p in pairs if padrao.search(p[0])]
        sys.stderr.write("Filtro '%s': %d de %d amostras -> %s\n"
                         % (args.filtrar, len(pairs), antes,
                            ", ".join(s for s, _, _ in pairs[:8]) or "(nenhuma)"))
        if not pairs:
            sys.exit("Nenhuma amostra casou com o filtro.")
    if args.listar:
        for s, f1, _f2 in pairs[:20]:
            sys.stderr.write("  %-12s %s\n" % (s, os.path.relpath(f1, raiz)))
        sys.exit(0)

    random.seed(args.seed)
    chosen = random.sample(pairs, min(args.samples, len(pairs)))
    per = max(1, args.reads // len(chosen))

    r1_reads, r2_reads = [], []
    for sample, f1, f2 in chosen:
        a = head_seqs(f1, per, args.k + args.maxphase + args.primer_len)
        b = head_seqs(f2, per, args.k + args.maxphase + args.primer_len)
        m = min(len(a), len(b))
        r1_reads.extend(a[:m])
        r2_reads.extend(b[:m])
    sys.stderr.write("Amostradas %d amostras, %d pares de reads.\n"
                     % (len(chosen), len(r1_reads)))

    fwd = discover(r1_reads, args.k, args.maxphase, args.max_mm,
                   300, args.primer_len, args.min_frac, "R1")
    rev = discover(r2_reads, args.k, args.maxphase, args.max_mm,
                   300, args.primer_len, args.min_frac, "R2")

    # ---- relatorio dos motivos
    def report(title, items):
        print("\n%s" % title)
        print("%-26s %10s %8s %4s  %s"
              % ("PRIMER (consenso IUPAC)", "READS", "%", "ESP", "PHASE 0..n"))
        print("-" * 92)
        for r in items:
            hist = " ".join(str(r["hist"].get(o, 0) // 100) for o in range(args.maxphase + 1))
            print("%-26s %10d %7.2f%% %4d  %s"
                  % (r["primer"], r["n"], 100 * r["frac"], r["esp"], hist))

    report("PRIMERS FORWARD (R1)", fwd)
    report("PRIMERS REVERSE (R2)", rev)
    print("\n(PHASE: reads por deslocamento, em centenas.  ESP: quantos")
    print(" deslocamentos concentram >=2% das reads. Primer com bloco de fase")
    print(" espalha por varios; ESP=1 significa posicao fixa — o que pode ser")
    print(" artefato OU um primer do painel que simplesmente nao e phased.)")

    # ---- pareamento F x R: e o par que define a regiao, nao a read isolada
    #      varios primers do painel sao reverso-complementares entre si,
    #      entao atribuir regiao por R1 sozinho produz troca de regiao.
    f_idx = [r["hits"] for r in fwd]
    r_idx = [r["hits"] for r in rev]
    table = Counter()
    for i in range(len(r1_reads)):
        f = next((j for j, h in enumerate(f_idx) if h[i] >= 0), None)
        r = next((j for j, h in enumerate(r_idx) if h[i] >= 0), None)
        table[(f, r)] += 1

    print("\nPARES OBSERVADOS (F em R1  x  R em R2)")
    print("%-26s %-26s %10s %8s" % ("FORWARD", "REVERSE", "PARES", "%"))
    print("-" * 76)
    total = len(r1_reads)
    regions = []
    for (f, r), n in table.most_common(20):
        if n < 0.005 * total:
            continue
        fs = fwd[f]["primer"] if f is not None else "(nenhum)"
        rs = rev[r]["primer"] if r is not None else "(nenhum)"
        print("%-26s %-26s %10d %7.2f%%" % (fs, rs, n, 100 * n / total))
        if f is not None and r is not None:
            regions.append((fs, rs, n))

    # ---- arquivos de saida
    with open(args.out + ".tsv", "w") as fh:
        fh.write("region\tforward\treverse\tpairs\n")
        for i, (fs, rs, n) in enumerate(regions, 1):
            fh.write("regiao_%02d\t%s\t%s\t%d\n" % (i, fs, rs, n))
    with open(args.out + ".fasta", "w") as fh:
        for i, (fs, rs, n) in enumerate(regions, 1):
            fh.write(">regiao_%02d_F\n%s\n>regiao_%02d_R\n%s\n" % (i, fs, i, rs))

    print("\nEscrito: %s.tsv e %s.fasta  (%d regioes candidatas)"
          % (args.out, args.out, len(regions)))
    print("Confira se o numero de regioes bate com o esperado: 7 (6 x 16S + ITS).")


if __name__ == "__main__":
    main()

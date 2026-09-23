#!/usr/bin/env python3
"""
check_overlap.py — os pares de reads se sobrepoem? E, quando nao, o que sao?

    bin/check_overlap.py R1.fastq.gz R2.fastq.gz
    bin/check_overlap.py R1.fastq.gz R2.fastq.gz --fasta nao_funde.fasta

Por que isto existe: quando uma regiao funde mal, o DADA2 nao diz por que. Ele
so entrega o que sobrou, e o que sobrou ja esta censurado — o par que nao fundiu
nao aparece em lugar nenhum, com comprimento nenhum. Entao a pergunta "o
amplicon e longo demais ou as reads tem erro demais?" nao se responde pela
saida do pipeline: responde-se pelos FASTQ de entrada.

Este script tenta o alinhamento por conta propria, direto nas reads, e separa:

  - par que SOBREPOE: devolve o comprimento do amplicon inferido. Isso da a
    distribuicao NAO censurada, que e o que falta em toda essa discussao.
  - par que nao sobrepoe: o amplicon e mais longo que R1+R2 menos a
    sobreposicao minima. Nenhum truncLen alcanca isso.
  - par que sobrepoe mas com diferencas demais: e erro, nao comprimento.

E, para os que nao sobrepoem, agrupa as reads R1 identicas e escreve as mais
frequentes em FASTA. Isso e o que fecha o diagnostico: a sequencia vai para o
BLAST e diz o que estava sendo amplificado fora do alvo. Um par universal de
16S que tambem pega 18S de eucarioto produz exatamente este quadro — amplicon
longo demais, sempre, em qualquer parametro.

Sem dependencia externa: so a biblioteca padrao. O cluster nem sempre tem
vsearch, e instalar um pacote para responder uma pergunta de uma linha e caro.
"""
import argparse
import gzip
import sys
from collections import defaultdict

COMPLEMENTO = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(s):
    return s.translate(COMPLEMENTO)[::-1]


def le_fastq(caminho, limite):
    """Devolve ate `limite` sequencias. Aceita .gz ou texto."""
    abrir = gzip.open if caminho.endswith(".gz") else open
    seqs = []
    with abrir(caminho, "rt") as fh:
        for i, linha in enumerate(fh):
            if i % 4 == 1:
                seqs.append(linha.strip().upper())
                if len(seqs) >= limite:
                    break
    return seqs


def sobreposicao(r1, r2rc, min_ov, max_dif_frac, semente=16):
    """
    Procura o encaixe de r2rc na cauda de r1.

    Devolve (comprimento_do_amplicon, n_diferencas) ou None.

    A busca e por semente e nao por forca bruta: comparar todas as
    sobreposicoes possiveis de todos os pares custa O(n x L^2) e nao termina
    em tempo util em Python. Uma semente exata de 16 bases tirada do inicio de
    r2rc, procurada dentro de r1, ja da o deslocamento candidato; so ele e
    verificado. Como a semente pode cair em cima de um erro de sequenciamento,
    tentamos algumas posicoes diferentes antes de desistir.
    """
    l1, l2 = len(r1), len(r2rc)
    for inicio in (0, 20, 40, 60):
        if inicio + semente > l2:
            break
        chave = r2rc[inicio:inicio + semente]
        pos = r1.find(chave)
        while pos != -1:
            # r2rc[inicio] alinha com r1[pos]; entao r2rc[0] alinharia em
            # r1[pos-inicio]. Esse e o comeco da sobreposicao.
            comeco = pos - inicio
            if comeco >= 0:
                ov = l1 - comeco
                if ov >= min_ov and ov <= l2:
                    a, b = r1[comeco:], r2rc[:ov]
                    dif = sum(1 for x, y in zip(a, b) if x != y and x != "N" and y != "N")
                    if dif <= max(1, int(max_dif_frac * ov)):
                        return (l1 + l2 - ov, dif)
            pos = r1.find(chave, pos + 1)
    return None


def percentil(hist, p):
    total = sum(hist.values())
    if not total:
        return 0
    alvo = total * p / 100.0
    acum = 0
    for v in sorted(hist):
        acum += hist[v]
        if acum >= alvo:
            return v
    return max(hist)


def main():
    ap = argparse.ArgumentParser(
        description="Do the read pairs overlap? And if not, what are they?")
    ap.add_argument("r1")
    ap.add_argument("r2")
    ap.add_argument("--reads", type=int, default=20000,
                    help="how many pairs to examine (default 20000)")
    ap.add_argument("--min-overlap", type=int, default=12,
                    help="minimum overlap, as in DADA2 (default 12)")
    ap.add_argument("--max-diff", type=float, default=0.10,
                    help="mismatch fraction tolerated in the overlap "
                         "(default 0.10 — deliberately loose: the question "
                         "here is whether they CAN overlap at all)")
    ap.add_argument("--fasta", default=None,
                    help="write the most frequent non-overlapping R1 reads "
                         "here, as FASTA, ready to paste into BLAST")
    ap.add_argument("--top", type=int, default=10,
                    help="how many sequences to write (default 10)")
    ap.add_argument("--prefix", type=int, default=250,
                    help="how much of the R1 read to keep in the FASTA "
                         "(default 250)")
    ap.add_argument("--group", type=int, default=120,
                    help="prefix length used to group identical reads "
                         "(default 120; shorter tolerates sequencing error)")
    args = ap.parse_args()

    s1 = le_fastq(args.r1, args.reads)
    s2 = le_fastq(args.r2, args.reads)
    n = min(len(s1), len(s2))
    if not n:
        sys.exit("ERROR: no reads read. Check the paths.")

    hist = defaultdict(int)
    # Agrupar por um prefixo CURTO, e nao pela read inteira: duas reads do
    # mesmo organismo quase nunca sao identicas em 250 bases (erro de
    # sequenciamento basta para separa-las), e agrupar pela read inteira
    # devolveria dez sequencias de contagem 1, que nao dizem nada sobre
    # abundancia. O prefixo curto junta as variantes; guardamos uma read
    # inteira de cada grupo como representante para o BLAST.
    sem_encaixe = defaultdict(int)     # chave curta -> contagem
    representante = {}                 # chave curta -> read para o FASTA
    n_ok = n_dif = n_sem = 0
    for i in range(n):
        r1, r2rc = s1[i], revcomp(s2[i])
        res = sobreposicao(r1, r2rc, args.min_overlap, args.max_diff)
        if res:
            hist[res[0]] += 1
            n_ok += 1
        else:
            # segunda tentativa, tolerando muito mais diferenca: separa
            # "nao alinha" de "alinha mal"
            res2 = sobreposicao(r1, r2rc, args.min_overlap, 0.35)
            if res2:
                n_dif += 1
            else:
                n_sem += 1
                chave = r1[:args.group]
                sem_encaixe[chave] += 1
                representante.setdefault(chave, r1[:args.prefix])

    print("Pairs examined: %d   (from %s)" % (n, args.r1))
    print("")
    print("  overlap found, clean       %8d  %5.1f%%" % (n_ok, 100.0 * n_ok / n))
    print("  overlap found, noisy       %8d  %5.1f%%" % (n_dif, 100.0 * n_dif / n))
    print("  NO overlap at all          %8d  %5.1f%%" % (n_sem, 100.0 * n_sem / n))
    print("")

    if n_ok:
        print("Amplicon length among the pairs that DO overlap:")
        print("  p5 %d   p50 %d   p95 %d   max %d"
              % (percentil(hist, 5), percentil(hist, 50),
                 percentil(hist, 95), max(hist)))
        print("")
        print("  This is the UNCENSORED distribution — measured from the reads")
        print("  themselves, not from what survived DADA2.")
        print("")

    limite = len(s1[0]) + len(s2[0]) - args.min_overlap
    if n_sem:
        print("The %.1f%% with no overlap have an amplicon longer than %d bp"
              % (100.0 * n_sem / n, limite))
        print("(read lengths %d + %d minus the %d bp minimum overlap). No"
              % (len(s1[0]), len(s2[0]), args.min_overlap))
        print("truncLen reaches that: truncating only makes the reads shorter.")
        print("")

    if args.fasta and sem_encaixe:
        top = sorted(sem_encaixe.items(), key=lambda kv: -kv[1])[:args.top]
        with open(args.fasta, "w") as fh:
            for k, (chave, c) in enumerate(top, 1):
                fh.write(">nonoverlapping_%d count=%d pct=%.2f\n%s\n"
                         % (k, c, 100.0 * c / n, representante[chave]))
        print("Wrote %s — the %d most frequent non-overlapping R1 reads."
              % (args.fasta, len(top)))
        print("Paste them into BLAST (nt). If they come back as 18S, or as")
        print("anything that is not bacterial 16S, the loss is off-target")
        print("amplification and not a pipeline setting.")
        print("")
        print("  most frequent, alone: %.2f%% of all pairs"
              % (100.0 * top[0][1] / n))
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except BrokenPipeError:
        try:
            sys.stdout.close()
        finally:
            sys.exit(0)

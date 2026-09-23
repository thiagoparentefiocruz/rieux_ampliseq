#!/usr/bin/env python3
"""
check_overlap.py — mede o amplicon nas reads, antes de qualquer censura.

Dois modos.

1) UM PAR DE ARQUIVOS, para investigar uma regiao que fundiu mal:

       bin/check_overlap.py R1.fastq.gz R2.fastq.gz --fasta nao_funde.fasta

   Separa os pares em tres casos — sobrepoe limpo, sobrepoe com ruido, nao
   sobrepoe — e, para os que nao sobrepoem, escreve as R1 mais frequentes em
   FASTA, prontas para o BLAST dizer o que estava sendo amplificado fora do
   alvo.

2) O DIRETORIO DO SPLIT INTEIRO, para escolher truncLen com medida e nao com
   tabela:

       bin/check_overlap.py --split <projeto>/split/split \\
                            --primers <projeto>/primers.tsv \\
                            --out <projeto>/inserts.tsv

   Percorre regiao por regiao e escreve os percentis do inserto REAL.

Por que isto existe: quando uma regiao funde mal, o DADA2 nao diz por que. Ele
so entrega o que sobrou, e o que sobrou ja esta censurado — o par que nao
fundiu nao aparece na saida com comprimento nenhum. Estimar o inserto pela
mediana dos ASVs sobreviventes e, por isso, circular: um truncLen curto corta
a cauda longa, a mediana do que sobra desce, e a proxima estimativa confirma o
truncLen curto. Foi assim que o V7V9 do painel ficou com 214/213 — teto de 415
bp — para um amplicon de ~432 bp: os quatro controles fundiram ZERO read, e a
regiao inteira marcou 6,6%. Com 255/255 a mesma regiao foi a 99,9%.

A medida honesta vem dos FASTQ de entrada, que e o que este script faz.

Comprimento bruto e liquido: a read no split ainda carrega o bloco de fase e o
primer, porque o split roteia com `cutadapt --action=none` e nao corta. Quem
corta e o cutadapt do ampliseq, depois. Com --primers o script localiza o
primer em cada read e desconta fase e primer par a par, devolvendo o inserto
LIQUIDO — que e o comprimento com que o truncLen tem de ser comparado. Sem
--primers, os numeros sao brutos e nao servem para escolher truncLen.

Sem dependencia externa: so a biblioteca padrao.
"""
import argparse
import gzip
import os
import re
import sys
from collections import defaultdict

COMPLEMENTO = str.maketrans("ACGTNacgtn", "TGCANtgcan")

IUPAC = {
    "A": "A", "C": "C", "G": "G", "T": "T",
    "R": "[AG]", "Y": "[CT]", "S": "[GC]", "W": "[AT]",
    "K": "[GT]", "M": "[AC]", "B": "[CGT]", "D": "[AGT]",
    "H": "[ACT]", "V": "[ACG]", "N": "[ACGT]", "I": "[ACGT]",
}


def revcomp(s):
    return s.translate(COMPLEMENTO)[::-1]


def regex_primer(p):
    return re.compile("".join(IUPAC.get(c, c) for c in p.strip().upper()))


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
                    dif = sum(1 for x, y in zip(a, b)
                              if x != y and x != "N" and y != "N")
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


def desconto(read, rx, janela=45):
    """
    Quanto ha antes do fim do primer: bloco de fase + primer.

    O bloco de fase antecede o primer e varia de read para read, entao este
    desconto e por READ e nao uma constante. Devolve None quando o primer nao
    aparece na janela inicial — pair descartado da medida liquida, nunca
    'corrigido' por uma media.
    """
    if rx is None:
        return 0
    m = rx.search(read[:janela])
    return m.end() if m else None


def mede_par(r1, r2, min_ov, max_dif, rxF=None, rxR=None):
    """
    Devolve (bruto, liquido_ou_None, motivo).

    motivo: 'ok', 'ruido', 'sem_encaixe'
    """
    r2rc = revcomp(r2)
    res = sobreposicao(r1, r2rc, min_ov, max_dif)
    if not res:
        if sobreposicao(r1, r2rc, min_ov, 0.35):
            return (None, None, "ruido")
        return (None, None, "sem_encaixe")
    bruto = res[0]
    dF, dR = desconto(r1, rxF), desconto(r2, rxR)
    if dF is None or dR is None:
        return (bruto, None, "ok")
    return (bruto, bruto - dF - dR, "ok")


def primers_do_arquivo(caminho):
    """region -> (forward, reverse), do primers.tsv."""
    tab = {}
    with open(caminho) as fh:
        cab = None
        for linha in fh:
            if linha.startswith("#") or not linha.strip():
                continue
            c = linha.rstrip("\n").split("\t")
            if cab is None:
                cab = [x.strip().lower() for x in c]
                continue
            d = dict(zip(cab, c))
            reg = d.get("region") or d.get("regiao")
            fw = d.get("forward")
            rv = d.get("reverse")
            if reg and fw and rv:
                tab[reg.strip()] = (fw.strip(), rv.strip())
    return tab


# --------------------------------------------------------------------- modos
def modo_lote(args):
    """Percorre <split>/<REGIAO>/*_R1.fastq.gz e escreve os percentis."""
    if not os.path.isdir(args.split):
        sys.exit("ERROR: %s is not a directory" % args.split)
    primers = primers_do_arquivo(args.primers) if args.primers else {}
    if not primers:
        sys.stderr.write(
            "WARNING: no --primers table, so the lengths below are RAW "
            "(phasing block and primers included).\n"
            "         They are NOT comparable with truncLen.\n\n")

    rx_ctrl = re.compile(args.controls)
    regioes = sorted(d for d in os.listdir(args.split)
                     if os.path.isdir(os.path.join(args.split, d))
                     and d != "unknown" and not d.startswith("."))
    linhas = []
    sem_medida = []
    print("  %-7s %-8s %8s %9s %7s %7s %7s %7s"
          % ("region", "set", "pairs", "no_ovlap", "p50", "p75", "p90", "p95"))
    print("  " + "-" * 66)
    for reg in regioes:
        d = os.path.join(args.split, reg)
        arquivos = sorted(f for f in os.listdir(d) if "_R1" in f
                          and f.endswith(".fastq.gz"))
        if args.samples and len(arquivos) > args.samples:
            # ESPACADO, nao os primeiros N. Pegar os primeiros em ordem
            # alfabetica escolhe um subconjunto com cara de amostragem e que
            # nao e: numa corrida cujas amostras comecam com H e cujos
            # controles comecam com S, `[:6]` devolve seis amostras e zero
            # controles — e foi exatamente o que aconteceu, dando uma tabela
            # que descrevia so metade do experimento.
            ult = len(arquivos) - 1
            if args.samples == 1:
                idx = [0]
            else:
                # inclui o PRIMEIRO e o ULTIMO: controles costumam cair numa
                # das pontas da ordem alfabetica, e e exatamente a ponta que
                # nao pode faltar.
                idx = sorted(set(
                    int(round(i * ult / float(args.samples - 1)))
                    for i in range(args.samples)))
            arquivos = [arquivos[i] for i in idx]

        # hist[tipo] e n[tipo], com tipo em {"sample", "control"}
        hist = {"sample": defaultdict(int), "control": defaultdict(int)}
        n_tot = {"sample": 0, "control": 0}
        n_sem = {"sample": 0, "control": 0}
        rxF = rxR = None
        if reg in primers:
            rxF, rxR = (regex_primer(primers[reg][0]),
                        regex_primer(primers[reg][1]))
        for f1 in arquivos:
            c1 = os.path.join(d, f1)
            c2 = os.path.join(d, f1.replace("_R1", "_R2"))
            if not os.path.isfile(c2):
                continue
            nome = re.sub(r"_R1.*$", "", f1)
            tipo = "control" if rx_ctrl.search(nome) else "sample"
            s1 = le_fastq(c1, args.reads)
            s2 = le_fastq(c2, args.reads)
            for i in range(min(len(s1), len(s2))):
                n_tot[tipo] += 1
                bruto, liquido, motivo = mede_par(
                    s1[i], s2[i], args.min_overlap, args.max_diff, rxF, rxR)
                if motivo == "sem_encaixe":
                    n_sem[tipo] += 1
                elif liquido is not None and liquido > 0:
                    hist[tipo][liquido] += 1
                elif bruto is not None and not primers:
                    hist[tipo][bruto] += 1

        total = n_tot["sample"] + n_tot["control"]
        if not total:
            continue
        juntos = defaultdict(int)
        for tipo in ("sample", "control"):
            for v, c in hist[tipo].items():
                juntos[v] += c

        def linha(rotulo, h, nt, ns):
            if not nt:
                return None
            if not h:
                print("  %-7s %-8s %8d %8.1f%% %7s %7s %7s %7s"
                      % (reg, rotulo, nt, 100.0 * ns / nt, "-", "-", "-", "-"))
                return None
            p = [percentil(h, x) for x in (50, 75, 90, 95)]
            print("  %-7s %-8s %8d %8.1f%% %7d %7d %7d %7d"
                  % (reg, rotulo, nt, 100.0 * ns / nt, p[0], p[1], p[2], p[3]))
            return p

        linha("samples", hist["sample"], n_tot["sample"], n_sem["sample"])
        linha("controls", hist["control"], n_tot["control"], n_sem["control"])
        p = linha("ALL", juntos, total, n_sem["sample"] + n_sem["control"])
        print("")
        if p is None:
            sem_medida.append(reg)
            continue
        linhas.append((reg, total,
                       100.0 * (n_sem["sample"] + n_sem["control"]) / total)
                      + tuple(p))

    if sem_medida:
        print("  NOT MEASURED: %s" % " ".join(sem_medida))
        print("  No pair gave a usable length: either nothing overlapped, or")
        print("  the primer was not found in the reads (check that this")
        print("  region's row in the primers table matches these files).")
        print("  They are absent from the output file — no truncLen is")
        print("  invented for a region that was not measured.")
        print("")

    if primers:
        print("  Lengths are NET: the phasing block and the primer were located")
        print("  and subtracted read by read. Compare them with truncLenF +")
        print("  truncLenR - 12 directly.")
    print("  'no_ovlap' is the share of pairs whose amplicon exceeds what two")
    print("  reads can span. Those are invisible to every length above — and")
    print("  to DADA2. A high value there is off-target amplification, not a")
    print("  setting (check one sample with the single-pair mode and BLAST).")

    if args.out:
        pasta = os.path.dirname(os.path.abspath(args.out))
        if pasta:
            os.makedirs(pasta, exist_ok=True)
        with open(args.out, "w") as fh:
            fh.write("# Medido por check_overlap.py nas reads do split.\n"
                     "# Inserto LIQUIDO (sem fase, sem primer), em bp.\n"
                     if primers else
                     "# ATENCAO: medido SEM --primers: comprimento BRUTO.\n")
            fh.write("region\tpairs\tpct_no_overlap\tp50\tp75\tp90\tp95\n")
            for l in linhas:
                fh.write("%s\t%d\t%.2f\t%d\t%d\t%d\t%d\n" % l)
        print("\nWrote %s" % args.out)
    return 0


def modo_par(args):
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
    sem_encaixe = defaultdict(int)
    representante = {}
    n_ok = n_dif = n_sem = 0
    for i in range(n):
        r1 = s1[i]
        bruto, _, motivo = mede_par(r1, s2[i], args.min_overlap, args.max_diff)
        if motivo == "ok":
            hist[bruto] += 1
            n_ok += 1
        elif motivo == "ruido":
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
        print("  Careful: these lengths are RAW. On the output of a split done")
        print("  with cutadapt --action=none the reads still carry the phasing")
        print("  block and both primers, which ampliseq's own cutadapt removes")
        print("  later. Subtract that before comparing with ASV lengths, or the")
        print("  same molecule will look like two. The --split mode subtracts")
        print("  them read by read and reports the net insert.")
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


def main():
    ap = argparse.ArgumentParser(
        description="Measure the amplicon in the reads, before DADA2 censors it.")
    ap.add_argument("r1", nargs="?", help="R1 FASTQ (single-pair mode)")
    ap.add_argument("r2", nargs="?", help="R2 FASTQ (single-pair mode)")
    ap.add_argument("--split", default=None,
                    help="batch mode: the split/<REGION>/ directory tree")
    ap.add_argument("--primers", default=None,
                    help="primers.tsv — needed for NET lengths")
    ap.add_argument("--out", default=None,
                    help="batch mode: write the percentiles here as TSV")
    ap.add_argument("--samples", type=int, default=0,
                    help="batch mode: use at most N samples per region, "
                         "evenly spaced across the sorted list "
                         "(default 0 = all)")
    ap.add_argument("--controls", default="^[Ss]mart",
                    help="regex for control sample names (default ^[Ss]mart). "
                         "Controls are reported separately: a known community "
                         "and a real sample can fail for opposite reasons")
    ap.add_argument("--reads", type=int, default=20000,
                    help="pairs per file (default 20000)")
    ap.add_argument("--min-overlap", type=int, default=12,
                    help="minimum overlap, as in DADA2 (default 12)")
    ap.add_argument("--max-diff", type=float, default=0.10,
                    help="mismatch fraction tolerated in the overlap "
                         "(default 0.10 — deliberately loose: the question "
                         "here is whether they CAN overlap at all)")
    ap.add_argument("--fasta", default=None,
                    help="single-pair mode: write the most frequent "
                         "non-overlapping R1 reads here, as FASTA")
    ap.add_argument("--top", type=int, default=10,
                    help="how many sequences to write (default 10)")
    ap.add_argument("--prefix", type=int, default=250,
                    help="how much of the R1 read to keep in the FASTA")
    ap.add_argument("--group", type=int, default=120,
                    help="prefix length used to group identical reads "
                         "(default 120; shorter tolerates sequencing error)")
    args = ap.parse_args()

    if args.split:
        return modo_lote(args)
    if not (args.r1 and args.r2):
        ap.error("give R1 and R2, or --split DIR")
    return modo_par(args)


if __name__ == "__main__":
    try:
        sys.exit(main())
    except BrokenPipeError:
        try:
            sys.stdout.close()
        finally:
            sys.exit(0)

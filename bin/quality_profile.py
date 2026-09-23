#!/usr/bin/env python3
"""
quality_profile.py — escolhe truncLenF/truncLenR por regiao

Por que isto existe: o NextSeq 2000 emite qualidade em 4 bins (Q2, Q12,
Q24, Q40). O `--trunc_qmin` do ampliseq decide o corte pela mediana de
qualidade por ciclo — com 4 bins a mediana fica travada em Q40 ate
despencar de vez, e o corte automatico nao serve. Os valores tem que ser
definidos a mao, e este script os calcula.

Duas restricoes disputam o mesmo orcamento:

  QUALIDADE  — cortar cedo o bastante para nao arrastar bases ruins, que
               o DADA2 penaliza via erro esperado (maxEE).
  SOBREPOSICAO — cortar tarde o bastante para R1 e R2 ainda se
               sobreporem. sobreposicao = truncF + truncR - inserto.
               Sem ela o par nao funde e a regiao se perde inteira.

O script mede a qualidade real por ciclo, calcula o corte que o erro
esperado permite, e verifica se ele cabe no orcamento de sobreposicao.
Quando nao cabe, avisa e mostra o deficit em vez de devolver um numero
que so falharia horas depois.

Uso:
    python3 quality_profile.py <dir_split> [--amostras 12] [--reads 20000]

Compativel com Python 3.6.
"""

import argparse
import glob
import gzip
import os
import random
import sys
from collections import Counter

# Comprimento do inserto por regiao (amplicon menos os primers).
#
# ESTES VALORES SAO MEDIDOS, NAO DERIVADOS DE COORDENADAS.
#
# A versao anterior usava as coordenadas canonicas dos primers em E. coli, e
# isso nos custou caro: dava 271 para o V1V2 quando o real e ~310, e o truncLen
# calculado a partir dali deixava o teto de fusao ABAIXO da mediana da
# comunidade. Resultado: 98% dos pares do V1V2 nao fundiram, sem erro nenhum no
# log. E. coli nao e a regua — a regua e o que amplifica de fato.
#
# Os numeros abaixo sao a mediana do amplicon extraido do SILVA por PCR in
# silico com os primers do painel (bin/validate_sidle_regions.py), conferida
# contra a distribuicao de comprimento dos ASVs observados. As duas batem nas
# pontas em todas as seis regioes.
#
# Para OUTRO painel ou outra comunidade, meça de novo: --asv aponta para um
# asv_length.tsv de uma execucao anterior e usa a mediana observada, que e
# sempre melhor que qualquer tabela.
INSERTO = {
    "V1V2": 310,   # 27F  -> 338R
    "V2V3": 392,   # 104F -> 519R
    "V3V4": 421,   # 341F -> 806R
    "V4V5": 372,   # 515F -> 926R
    "V5V7": 368,   # 805F -> 1193R
    "V7V9": 377,   # 1100F-> 1492R
    "ITS1": None,  # comprimento muito variavel: nao truncar
}


def insertos_do_tsv(caminho, coluna):
    """region -> inserto, do inserts.tsv do check_overlap.py."""
    out = {}
    with open(caminho) as fh:
        cab = None
        for linha in fh:
            if linha.startswith("#") or not linha.strip():
                continue
            c = linha.rstrip("\n").split("\t")
            if cab is None:
                cab = [x.strip().lower() for x in c]
                if coluna.lower() not in cab:
                    sys.exit("ERROR: %s has no column '%s' (has: %s)"
                             % (caminho, coluna, ", ".join(cab)))
                continue
            d = dict(zip(cab, c))
            reg = (d.get("region") or "").strip()
            try:
                v = int(d[coluna.lower()])
            except (KeyError, ValueError):
                continue
            if reg and v > 0:
                out[reg] = v
    return out


def insertos_medidos(caminho):
    """Mediana do comprimento de ASV por regiao, de um asv_length.tsv."""
    from collections import defaultdict
    hist = defaultdict(lambda: defaultdict(int))
    with open(caminho) as fh:
        cab = fh.readline().rstrip("\n").split("\t")
        idx = {n: i for i, n in enumerate(cab)}
        creg = "region" if "region" in idx else "regiao"
        clen = "length" if "length" in idx else "comprimento"
        for linha in fh:
            c = linha.rstrip("\n").split("\t")
            if len(c) < len(cab):
                continue
            hist[c[idx[creg]]][int(c[idx[clen]])] += int(c[idx["n_asv"]])
    out = {}
    for reg, h in hist.items():
        total = sum(h.values())
        if not total:
            continue
        acum = 0
        for v in sorted(h):
            acum += h[v]
            if acum >= total / 2.0:
                out[reg] = v
                break
    return out

SOBREPOSICAO_MIN = 20
EE_MAX = 2.0            # erro esperado acumulado tolerado, padrao do DADA2
CICLOS = 300
FASE_MAX = 13           # bloco de fase observado nos dados (handbook diz 0-11)
MARGEM = 8              # folga para nao descartar read por ser curta demais

# Comprimento dos primers por regiao, para calcular quanto sobra da read
# depois que o cutadapt remove o bloco de fase e o primer.
PRIMER_LEN = {
    "V1V2": (20, 20), "V2V3": (17, 20), "V3V4": (17, 22),
    "V4V5": (20, 20), "V5V7": (22, 20), "V7V9": (15, 22),
}


def qualidades(caminho, nreads, maxciclos=300):
    """Soma de probabilidade de erro por ciclo, e contagem de bases."""
    soma = [0.0] * maxciclos
    n = [0] * maxciclos
    bins = Counter()
    lidas = 0
    with gzip.open(caminho, "rt") as fh:
        for i, linha in enumerate(fh):
            if i % 4 != 3:
                continue
            linha = linha.rstrip("\n")
            for c, ch in enumerate(linha[:maxciclos]):
                q = ord(ch) - 33
                soma[c] += 10.0 ** (-q / 10.0)
                n[c] += 1
                bins[q] += 1
            lidas += 1
            if lidas >= nreads:
                break
    return soma, n, bins, lidas


def corte_por_ee(soma, n, limite):
    """Ultimo ciclo em que o erro esperado acumulado ainda cabe no limite."""
    acumulado = 0.0
    for c in range(len(n)):
        if n[c] == 0:
            return c
        acumulado += soma[c] / n[c]
        if acumulado > limite:
            return c
    return len(n)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("split", help="the split/ directory (the one holding the regions)")
    ap.add_argument("--samples", type=int, default=12)
    ap.add_argument("--reads", type=int, default=20000)
    ap.add_argument("--ee", type=float, default=EE_MAX)
    ap.add_argument("--overlap", type=int, default=SOBREPOSICAO_MIN,
                    help="minimum acceptable overlap")
    ap.add_argument("--overlap-target", type=int, default=50,
                    help="overlap to aim for. Measured in the pilot: 50 merges better "
                         "than 89. Do not truncate at the maximum.")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--asv", default=None,
                    help="asv_length.tsv from an earlier run; the observed median per "
                         "region replaces the built-in table. CIRCULAR — see --inserts")
    ap.add_argument("--inserts", default=None,
                    help="inserts.tsv from check_overlap.py --split: the insert "
                         "measured in the READS, uncensored. Preferred over --asv")
    ap.add_argument("--insert-column", default="p90",
                    help="which percentile of the measured insert to aim at "
                         "(default p90: the median leaves half the community "
                         "above the ceiling)")
    ap.add_argument("--out", default=None,
                    help="write the chosen truncLen values as region_params.tsv")
    args = ap.parse_args()

    random.seed(args.seed)

    # --inserts tem precedencia sobre --asv de proposito.
    #
    # O --asv mede o inserto na mediana dos ASVs que SOBREVIVERAM, e isso e
    # circular: truncLen curto corta a cauda longa, a mediana do que sobra
    # desce, e a estimativa seguinte confirma o truncLen curto. Foi esse laco
    # que deixou o V7V9 com teto de 415 bp para um amplicon de ~432 — os
    # quatro controles fundiram ZERO read. O --inserts vem do
    # check_overlap.py, que mede nas reads de entrada, antes de qualquer
    # censura, e por isso enxerga a cauda que o --asv nunca ve.
    if args.inserts:
        medidos = insertos_do_tsv(args.inserts, args.insert_column)
        origem = "%s (%s)" % (args.inserts, args.insert_column)
    elif args.asv:
        medidos = insertos_medidos(args.asv)
        origem = "%s (median of surviving ASVs — may be censored)" % args.asv
    else:
        medidos = {}
        origem = None
    for reg, v in medidos.items():
        # ITS fica de fora: comprimento variavel, nao se trunca.
        if reg.upper().startswith("ITS"):
            continue
        INSERTO[reg] = v
    if origem:
        print("insert lengths taken from %s: %s\n"
              % (origem, ", ".join("%s=%d" % (r, medidos[r])
                                   for r in sorted(medidos))))
    regioes = sorted(d for d in os.listdir(args.split)
                     if os.path.isdir(os.path.join(args.split, d)) and d != "unknown")

    print("%-7s %6s %6s %7s %7s %9s  %s"
          % ("REGION", "truncF", "truncR", "insert", "overlap", "reads", "note"))
    print("-" * 82)

    linhas_cfg = []
    for reg in regioes:
        r1s = sorted(glob.glob(os.path.join(args.split, reg, "*_R1.fastq.gz")))
        r1s = [f for f in r1s if os.path.getsize(f) > 10000]
        if not r1s:
            print("%-7s %s" % (reg, "no files with data"))
            continue
        escolhidos = random.sample(r1s, min(args.samples, len(r1s)))
        por = max(1, args.reads // len(escolhidos))

        somaF = nF = somaR = nR = None
        binsT = Counter()
        total = 0
        for f1 in escolhidos:
            f2 = f1.replace("_R1.fastq.gz", "_R2.fastq.gz")
            for caminho, alvo in ((f1, "F"), (f2, "R")):
                s, n, b, lidas = qualidades(caminho, por)
                binsT += b
                if alvo == "F":
                    somaF = s if somaF is None else [a + x for a, x in zip(somaF, s)]
                    nF = n if nF is None else [a + x for a, x in zip(nF, n)]
                    total += lidas
                else:
                    somaR = s if somaR is None else [a + x for a, x in zip(somaR, s)]
                    nR = n if nR is None else [a + x for a, x in zip(nR, n)]

        cF = corte_por_ee(somaF, nF, args.ee)
        cR = corte_por_ee(somaR, nR, args.ee)

        inserto = INSERTO.get(reg)
        if inserto is None:
            # ITS nao tem comprimento fixo, entao nao se trunca — mas a regiao
            # ENTRA na tabela mesmo assim, com truncLen 0. Ficar de fora nao
            # significava "sem truncagem": significava sumir. O estagio `run`
            # percorre as regioes desta tabela, entao o ITS1 do painel
            # 16S/ITS simplesmente nao era analisado, sem uma linha de aviso.
            linhas_cfg.append((reg, 0, 0))
            print("%-7s %6s %6s %7s %7s %9d  %s"
                  % (reg, "-", "-", "variavel", "-", total,
                     "ITS: do not truncate (--illumina_pe_its)"))
            continue

        # O corte se aplica a read JA sem bloco de fase e sem primer, e o
        # DADA2 DESCARTA read mais curta que truncLen. Como o bloco de fase
        # varia de 0 a ~13 bases, o comprimento util varia entre reads da
        # mesma regiao: truncar no limite jogaria fora justamente as reads
        # de fase longa. O teto e calculado pelo pior caso, com margem.
        #
        # E nao ha custo em truncar abaixo do teto: o amplicon final sai da
        # fusao do par, cujo comprimento e o do inserto. Desde que a
        # sobreposicao feche, base extra na ponta da read e redundante.
        lenF, lenR = PRIMER_LEN.get(reg, (22, 22))
        tetoF = CICLOS - FASE_MAX - lenF - MARGEM
        tetoR = CICLOS - FASE_MAX - lenR - MARGEM
        cF, cR = min(cF, tetoF), min(cR, tetoR)
        teto = min(tetoF, tetoR)

        # Nao truncar no maximo: mirar uma SOBREPOSICAO-ALVO.
        #
        # Medido no piloto do Fabio (V3V4, 40 amostras): baixar o corte de
        # 257/257 para 240/235 — de 89 pb de sobreposicao para 50 — subiu a
        # fusao em +4,1 pontos percentuais na mediana, uniformemente.
        # Sobra maior nao ajuda: o DADA2 exige sobreposicao sem discordancia,
        # entao cada base extra ali e mais uma chance de erro derrubar o par,
        # e as bases extras sao as do fim da read, as piores.
        alvo = inserto + args.overlap_target
        if cF + cR > alvo:
            tF = min(tetoF, cF, (alvo + 1) // 2)
            tR = min(tetoR, cR, alvo // 2)
            falta = alvo - (tF + tR)
            if falta > 0:
                add = min(falta, tetoF - tF, cF - tF)
                tF += max(0, add)
                falta = alvo - (tF + tR)
            if falta > 0:
                add = min(falta, tetoR - tR, cR - tR)
                tR += max(0, add)
            cF, cR = tF, tR
        sobrep = cF + cR - inserto
        obs = ""
        if sobrep < args.overlap:
            falta = args.overlap - sobrep
            # tenta recuperar alongando ate o teto, se houver folga
            folga = (teto - cF) + (teto - cR)
            if folga >= falta:
                add = (falta + 1) // 2
                cF, cR = min(teto, cF + add), min(teto, cR + add)
                sobrep = cF + cR - inserto
                obs = "extended past the EE cut to close the merge"
            else:
                obs = "DOES NOT CLOSE: %d bp short even at the read limit" % falta
        linhas_cfg.append((reg, cF, cR))
        print("%-7s %6d %6d %7d %7d %9d  %s" % (reg, cF, cR, inserto, sobrep, total, obs))

    print("\nQuality bin distribution (confirms the NextSeq binning):")
    for q, c in sorted(binsT.items()):
        print("  Q%-3d %10d" % (q, c))

    print("\nampliseq parameters, per region:")
    for reg, cF, cR in linhas_cfg:
        if reg.upper().startswith("ITS"):
            print("  %-7s --illumina_pe_its --cut_its its1  (no trunclen)" % reg)
        else:
            print("  %-7s --trunclenf %d --trunclenr %d" % (reg, cF, cR))

    # Escrever o TSV, e nao so imprimir, e o que permite este passo ser um
    # ESTAGIO do wrapper em vez de uma consulta que alguem transcreve na mao.
    # Transcrever a mao foi exatamente como o truncLen do V1V2 entrou errado.
    if args.out:
        pasta = os.path.dirname(os.path.abspath(args.out))
        if pasta:
            os.makedirs(pasta, exist_ok=True)
        with open(args.out, "w") as fh:
            fh.write("# Gerado por quality_profile.py a partir das reads desta\n"
                     "# corrida. truncLen depende de comprimento e qualidade da\n"
                     "# corrida: NAO reaproveite esta tabela em outra.\n")
            fh.write("region\ttrunclenf\ttrunclenr\tdatabase\textra\n")
            for reg, cF, cR in linhas_cfg:
                if reg.upper().startswith("ITS"):
                    # ITS nao tem comprimento fixo: truncar corta ASV legitimo.
                    fh.write("%s\t0\t0\tunite\t--illumina_pe_its --cut_its its1\n"
                             % reg)
                else:
                    fh.write("%s\t%d\t%d\tsilva\t\n" % (reg, cF, cR))
        print("\nWrote %s" % args.out)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
perfil_qualidade.py — escolhe truncLenF/truncLenR por regiao

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
    python3 perfil_qualidade.py <dir_split> [--amostras 12] [--reads 20000]

Compativel com Python 3.6.
"""

import argparse
import glob
import gzip
import os
import random
import sys
from collections import Counter

# Comprimento do inserto por regiao (amplicon menos os primers), a partir
# das coordenadas canonicas em E. coli dos primers que recuperamos.
# Sobrescreva com --insertos se medir os valores reais nos seus dados.
INSERTO = {
    "V1V2": 271,   # 27F  -> 338R
    "V2V3": 375,   # 104F -> 519R
    "V3V4": 425,   # 341F -> 806R
    "V4V5": 371,   # 515F -> 926R
    "V5V7": 348,   # 805F -> 1193R
    "V7V9": 352,   # 1100F-> 1492R
    "ITS1": None,  # comprimento muito variavel: nao truncar
}

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
    ap.add_argument("split", help="diretorio split/ (o que tem as regioes dentro)")
    ap.add_argument("--amostras", type=int, default=12)
    ap.add_argument("--reads", type=int, default=20000)
    ap.add_argument("--ee", type=float, default=EE_MAX)
    ap.add_argument("--overlap", type=int, default=SOBREPOSICAO_MIN,
                    help="sobreposicao minima aceitavel")
    ap.add_argument("--overlap-alvo", type=int, default=50,
                    help="sobreposicao que se pretende obter. Medido no piloto: "
                         "50 funde melhor que 89. Nao truncar no maximo.")
    ap.add_argument("--seed", type=int, default=1)
    args = ap.parse_args()

    random.seed(args.seed)
    regioes = sorted(d for d in os.listdir(args.split)
                     if os.path.isdir(os.path.join(args.split, d)) and d != "unknown")

    print("%-7s %6s %6s %7s %7s %9s  %s"
          % ("REGIAO", "truncF", "truncR", "inserto", "sobrep", "reads", "observacao"))
    print("-" * 82)

    linhas_cfg = []
    for reg in regioes:
        r1s = sorted(glob.glob(os.path.join(args.split, reg, "*_R1.fastq.gz")))
        r1s = [f for f in r1s if os.path.getsize(f) > 10000]
        if not r1s:
            print("%-7s %s" % (reg, "sem arquivos com dado"))
            continue
        escolhidos = random.sample(r1s, min(args.amostras, len(r1s)))
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
            print("%-7s %6s %6s %7s %7s %9d  %s"
                  % (reg, "-", "-", "variavel", "-", total,
                     "ITS: nao truncar (--illumina_pe_its)"))
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
        alvo = inserto + args.overlap_alvo
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
                obs = "estendido alem do corte por EE para fechar o merge"
            else:
                obs = "NAO FECHA: faltam %d pb mesmo no limite da read" % falta
        linhas_cfg.append((reg, cF, cR))
        print("%-7s %6d %6d %7d %7d %9d  %s" % (reg, cF, cR, inserto, sobrep, total, obs))

    print("\nDistribuicao dos bins de qualidade (confirma o binning do NextSeq):")
    for q, c in sorted(binsT.items()):
        print("  Q%-3d %10d" % (q, c))

    print("\nParametros para o ampliseq, por regiao:")
    for reg, cF, cR in linhas_cfg:
        print("  %-7s --trunclenf %d --trunclenr %d" % (reg, cF, cR))
    print("\n  (ITS1: usar --illumina_pe_its --cut_its its1, sem trunclen)")


if __name__ == "__main__":
    main()

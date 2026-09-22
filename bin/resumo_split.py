#!/usr/bin/env python3
"""
resumo_split.py — metricas do split e samplesheets por regiao

Faz duas coisas depois que o split_regioes.sh roda:

  1. Le os relatorios JSON do cutadapt e monta a tabela de reads por
     amostra x regiao, com a fracao atribuida. O criterio de aceite do
     modulo e >90% de reads atribuidas por amostra; amostras abaixo disso
     sao listadas para inspecao.

  2. Escreve uma samplesheet por regiao no formato do nf-core/ampliseq
     (sample,fastq_1,fastq_2), pulando amostras com profundidade abaixo do
     piso. O piso se aplica POR REGIAO, nao sobre o total da amostra: uma
     amostra com 300 mil reads tem ~45 mil por regiao, e e esse numero que
     decide se ela sustenta uma analise de comunidade.

Uso:
    python3 resumo_split.py <dir_saida_do_split> [--minimo 1000]

Compativel com Python 3.6.
"""

import argparse
import csv
import glob
import gzip
import json
import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed


def contar_reads(caminho):
    """Numero de reads num FASTQ gzipado.

    Usa `zcat | wc -l` em vez do modulo gzip do Python: a descompressao em
    C e uma ordem de grandeza mais rapida que iterar linhas em Python, e
    aqui sao milhares de arquivos somando dezenas de GB. Se o zcat nao
    existir, cai no caminho puro-Python.
    """
    try:
        p1 = subprocess.Popen(["zcat", caminho], stdout=subprocess.PIPE,
                              stderr=subprocess.DEVNULL)
        p2 = subprocess.Popen(["wc", "-l"], stdin=p1.stdout,
                              stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
        p1.stdout.close()
        saida, _ = p2.communicate()
        if p2.returncode == 0:
            return int(saida.strip()) // 4
    except (OSError, ValueError):
        pass
    n = 0
    with gzip.open(caminho, "rb") as fh:
        for _ in fh:
            n += 1
    return n // 4


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("saida", help="diretorio usado no split_regioes.sh")
    ap.add_argument("--minimo", type=int, default=1000,
                    help="piso de reads por amostra x regiao (padrao 1000)")
    ap.add_argument("--jobs", type=int,
                    default=int(os.environ.get("SLURM_CPUS_PER_TASK", 8)),
                    help="processos paralelos para contar os FASTQ")
    ap.add_argument("--cabecalho", default="sampleID,forwardReads,reverseReads",
                    help="nomes das colunas da samplesheet. O nf-core/ampliseq "
                         "NAO usa o formato generico sample/fastq_1/fastq_2 — "
                         "ele tem os proprios nomes. Confira na sua versao com: "
                         "python3 -c \"import json,sys; "
                         "print(list(json.load(open(sys.argv[1]))"
                         "['items']['properties']))\" $WF/assets/schema_input.json")
    ap.add_argument("--colaboradores", metavar="DIR",
                    help="diretorio com <colaborador>/metadata.tsv (a saida do "
                         "organizar_colaboradores.py). Com isso as samplesheets "
                         "saem por colaborador e por regiao, ja com o ID "
                         "original no lugar do ID de sequenciamento.")
    args = ap.parse_args()

    split_dir = os.path.join(args.saida, "split")
    if not os.path.isdir(split_dir):
        sys.exit("Nao achei %s" % split_dir)

    regioes = sorted(d for d in os.listdir(split_dir)
                     if os.path.isdir(os.path.join(split_dir, d))
                     and d != "unknown")

    # ---- contagem por amostra x regiao, em paralelo
    alvos = []
    for reg in regioes + ["unknown"]:
        for f in glob.glob(os.path.join(split_dir, reg, "*_R1.fastq.gz")):
            amostra = os.path.basename(f)[:-len("_R1.fastq.gz")]
            alvos.append((amostra, reg, f))

    sys.stderr.write("Contando %d arquivos com %d processos...\n"
                     % (len(alvos), args.jobs))
    tabela = {}   # amostra -> {regiao: n}
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        futuros = {pool.submit(contar_reads, f): (a, r) for a, r, f in alvos}
        feitos = 0
        for fut in as_completed(futuros):
            amostra, reg = futuros[fut]
            tabela.setdefault(amostra, {})[reg] = fut.result()
            feitos += 1
            if feitos % 200 == 0:
                sys.stderr.write("  %d/%d\n" % (feitos, len(alvos)))

    if not tabela:
        sys.exit("Nenhum FASTQ encontrado em %s" % split_dir)

    # ---- relatorio
    cab = ["amostra"] + regioes + ["unknown", "total", "%atribuido"]
    largura = max(12, max(len(a) for a in tabela) + 1)
    print("".join(["%-*s" % (largura, cab[0])]
                  + ["%9s" % c for c in cab[1:]]))
    print("-" * (largura + 9 * (len(cab) - 1)))

    problemas, rasas = [], []
    linhas_csv = []
    for amostra in sorted(tabela):
        d = tabela[amostra]
        atrib = sum(d.get(r, 0) for r in regioes)
        desconhecido = d.get("unknown", 0)
        total = atrib + desconhecido
        frac = 100.0 * atrib / total if total else 0.0
        print("".join(["%-*s" % (largura, amostra)]
                      + ["%9d" % d.get(r, 0) for r in regioes]
                      + ["%9d" % desconhecido, "%9d" % total, "%8.1f%%" % frac]))
        linhas_csv.append([amostra] + [d.get(r, 0) for r in regioes]
                          + [desconhecido, total, round(frac, 2)])
        if frac < 90.0:
            problemas.append((amostra, frac))
        for r in regioes:
            if 0 < d.get(r, 0) < args.minimo:
                rasas.append((amostra, r, d[r]))

    with open(os.path.join(args.saida, "resumo_split.csv"), "w") as fh:
        w = csv.writer(fh)
        w.writerow(cab)
        w.writerows(linhas_csv)

    # ---- medianas por regiao
    print("\nMediana por regiao:")
    for r in regioes:
        vals = sorted(d.get(r, 0) for d in tabela.values())
        mediana = vals[len(vals) // 2] if vals else 0
        print("  %-10s %8d" % (r, mediana))

    # ---- criterio de aceite
    print("\n== Criterio de aceite ===========================================")
    if problemas:
        print("  %d amostra(s) com menos de 90%% atribuido:" % len(problemas))
        for a, f in sorted(problemas, key=lambda x: x[1])[:15]:
            print("     %-20s %.1f%%" % (a, f))
        print("  Investigue antes de seguir: primer incompleto, contaminacao")
        print("  de outra biblioteca, ou adaptador residual.")
    else:
        print("  Todas as amostras acima de 90% atribuido.")

    if rasas:
        print("\n  %d par(es) amostra x regiao abaixo de %d reads "
              "(excluidos das samplesheets):" % (len(rasas), args.minimo))
        for a, r, n in sorted(rasas, key=lambda x: x[2])[:15]:
            print("     %-20s %-10s %d" % (a, r, n))

    # ---- de qual colaborador e cada ID de sequenciamento
    dono_de, nome_de = {}, {}
    if args.colaboradores:
        for meta in glob.glob(os.path.join(args.colaboradores, "*", "metadata.tsv")):
            colaborador = os.path.basename(os.path.dirname(meta))
            with open(meta) as fh:
                cab = fh.readline().rstrip("\n").split("\t")
                iseq = cab.index("id_sequenciamento")
                isam = cab.index("sample")
                ictl = cab.index("controle") if "controle" in cab else None
                for linha in fh:
                    c = linha.rstrip("\n").split("\t")
                    if ictl is not None and c[ictl] == "sim":
                        continue      # controle entra em todos, tratado a parte
                    dono_de[c[iseq]] = colaborador
                    nome_de[c[iseq]] = c[isam]
        controles = [a for a in tabela if a.lower().startswith("smart")]
        print("\nColaboradores lidos: %d amostras mapeadas, %d controle(s)"
              % (len(dono_de), len(controles)))
    else:
        controles = []

    # ---- samplesheets
    print("\n== Samplesheets =================================================")
    ss_dir = os.path.join(args.saida, "samplesheets")

    colunas = [c.strip() for c in args.cabecalho.split(",")]

    def escrever(caminho, reg, amostras):
        n = 0
        with open(caminho, "w") as fh:
            fh.write("\t".join(colunas) + "\n")
            for amostra in amostras:
                if tabela.get(amostra, {}).get(reg, 0) < args.minimo:
                    continue
                r1 = os.path.abspath(os.path.join(split_dir, reg,
                                                  "%s_R1.fastq.gz" % amostra))
                r2 = r1.replace("_R1.fastq.gz", "_R2.fastq.gz")
                if not (os.path.exists(r1) and os.path.exists(r2)):
                    continue
                nome = nome_de.get(amostra, amostra)
                if not nome[:1].isalpha():
                    nome = "s" + nome
                fh.write("%s\t%s\t%s\n" % (nome, r1, r2))
                n += 1
        return n

    if dono_de:
        por_dono = {}
        for seq, colaborador in dono_de.items():
            por_dono.setdefault(colaborador, []).append(seq)
        for colaborador in sorted(por_dono):
            destino = os.path.join(ss_dir, colaborador)
            os.makedirs(destino, exist_ok=True)
            # os controles acompanham TODOS os colaboradores
            amostras = sorted(por_dono[colaborador]) + sorted(controles)
            linha = []
            for reg in regioes:
                n = escrever(os.path.join(destino, "samplesheet_%s.tsv" % reg),
                             reg, amostras)
                linha.append("%s=%d" % (reg, n))
            print("  %-10s %s" % (colaborador, "  ".join(linha)))
        print("\n  (cada samplesheet inclui os %d controle(s))" % len(controles))
    else:
        os.makedirs(ss_dir, exist_ok=True)
        for reg in regioes:
            caminho = os.path.join(ss_dir, "samplesheet_%s.tsv" % reg)
            n = escrever(caminho, reg, sorted(tabela))
            print("  %-10s %3d amostras  -> %s" % (reg, n, caminho))

    print("\nTabela completa em %s" % os.path.join(args.saida, "resumo_split.csv"))


if __name__ == "__main__":
    main()

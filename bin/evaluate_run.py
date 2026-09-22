#!/usr/bin/env python3
"""
evaluate_run.py — criterios de aceite de uma execucao do ampliseq

Le a saida de um `nextflow run nf-core/ampliseq` e responde as perguntas que
decidem se o resultado presta, sem precisar abrir PDF nem HTML:

  1. RETENCAO POR ETAPA. A cadeia entrada -> filtrado -> denoised -> merged
     -> nao-quimerico. Onde as reads somem diz o que ajustar: perda no merge
     e sobreposicao insuficiente (truncLen errado ou inserto maior que o
     suposto); perda em quimera alta demais sugere excesso de ciclos de PCR.

  2. CLASSIFICACAO POR NIVEL. Fracao de ASVs com atribuicao em cada rank.
     Reportar isso por projeto e obrigatorio: um numero global misturaria
     intestino bem resolvido com ambiental mal resolvido.

  3. O CONSTRUTO DO SMART CONTROL. O ASV dominante nos controles e o
     construto sintetico. Identifica-lo permite (a) remove-lo das amostras
     reais, onde aparece por carryover, e (b) usar o que sobra no controle
     como leitura de contaminacao.

Uso:
    python3 evaluate_run.py piloto/fabio_V3V4 [--controls '^Smart']

Compativel com Python 3.6.
"""

import argparse
import csv
import fnmatch
import os
import re
import sys


def achar(raiz, *padroes):
    """Primeiro arquivo que casa com um desses padroes (glob), em qualquer
    profundidade. Precisa ser glob e nao nome exato: o ampliseq acrescenta
    o identificador do banco ao nome da tabela de taxonomia, produzindo
    coisas como ASV_tax_species.silva_138.2.tsv."""
    for dirpath, _d, arquivos in os.walk(raiz):
        for p in padroes:
            casos = sorted(fnmatch.filter(arquivos, p))
            if casos:
                return os.path.join(dirpath, casos[0])
    return None


def ler_tsv(caminho):
    with open(caminho) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def coluna(campos, *chaves):
    """Nome real da coluna cujo nome contem alguma das chaves."""
    for k in chaves:
        for c in campos:
            if k.lower() == c.lower():
                return c
    for k in chaves:
        for c in campos:
            if k.lower() in c.lower():
                return c
    return None


def num(v):
    try:
        return float(str(v).replace(",", "").replace("%", ""))
    except (TypeError, ValueError):
        return None


def mediana(vals):
    v = sorted(x for x in vals if x is not None)
    return v[len(v) // 2] if v else None


def etapas(linhas):
    """Sequencia de colunas da cadeia, na ordem, com os nomes que existirem."""
    campos = list(linhas[0].keys())
    ordem = [
        ("entrada",      ("DADA2_input", "cutadapt_passing_filters", "input")),
        ("filtrado",     ("filtered",)),
        ("denoisedF",    ("denoisedF",)),
        ("denoisedR",    ("denoisedR",)),
        ("merged",       ("merged",)),
        ("nao-quimera",  ("nonchim", "non-chimeric")),
        ("pos-filtro-tax", ("filtered_tax_filter", "tax_filter")),
    ]
    saida = []
    for rotulo, chaves in ordem:
        c = coluna(campos, *chaves)
        if c:
            saida.append((rotulo, c))
    return saida


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("outdir")
    ap.add_argument("--controls", default="^Smart")
    ap.add_argument("--min-merge", type=float, default=70.0)
    ap.add_argument("--comparar", metavar="OUTRO_OUTDIR",
                    help="compara a taxa de fusao amostra a amostra com outra "
                         "execucao. Serve para decidir parametros (truncLen, "
                         "por exemplo) com numero em vez de intuicao.")
    args = ap.parse_args()

    # ------------------------------------------------ 1. retencao por etapa
    resumo = achar(args.outdir, "overall_summary.tsv")
    if not resumo:
        sys.exit("Nao achei overall_summary.tsv em %s" % args.outdir)
    linhas = ler_tsv(resumo)
    campos = list(linhas[0].keys())
    c_amostra = coluna(campos, "sample", "sampleID", "ID")
    passos = etapas(linhas)

    print("== 1. Retencao por etapa (mediana entre %d amostras) ============"
          % len(linhas))
    print("%-16s %12s %10s %10s" % ("ETAPA", "MEDIANA", "% ENTRADA", "% ETAPA ANT"))
    print("-" * 52)
    base = anterior = None
    for rotulo, col in passos:
        med = mediana([num(l.get(col)) for l in linhas])
        if med is None:
            continue
        if base is None:
            base = med
        pe = 100.0 * med / base if base else 0.0
        pa = 100.0 * med / anterior if anterior else 100.0
        print("%-16s %12.0f %9.1f%% %9.1f%%" % (rotulo, med, pe, pa))
        anterior = med

    # amostras problematicas no merge
    c_in = dict(passos).get("entrada")
    c_merge = dict(passos).get("merged")
    if c_in and c_merge:
        ruins = []
        for l in linhas:
            a, b = num(l.get(c_in)), num(l.get(c_merge))
            if a and b is not None and a > 0:
                pct = 100.0 * b / a
                if pct < args.min_merge:
                    ruins.append((l.get(c_amostra), pct))
        print()
        if ruins:
            print("  %d amostra(s) com fusao abaixo de %.0f%%:"
                  % (len(ruins), args.min_merge))
            for a, p in sorted(ruins, key=lambda x: x[1])[:12]:
                print("     %-28s %.1f%%" % (a, p))
            print("  Fusao baixa = sobreposicao insuficiente: reveja truncLen")
            print("  ou o comprimento do inserto assumido para esta regiao.")
        else:
            print("  Todas as amostras acima de %.0f%% de fusao." % args.min_merge)

    # ------------------------------------------- 1b. comparacao entre execucoes
    if args.comparar:
        outro = achar(args.comparar, "overall_summary.tsv")
        if not outro:
            print("\n  (nao achei overall_summary.tsv em %s)" % args.comparar)
        else:
            lb = ler_tsv(outro)
            cb = list(lb[0].keys())
            ca_b = coluna(cb, "sample", "sampleID", "ID")
            pb = dict(etapas(lb))

            def taxa(linhas_, c_amo, c_ent, c_mer):
                d = {}
                for l in linhas_:
                    a, b = num(l.get(c_ent)), num(l.get(c_mer))
                    if a and b is not None and a > 0:
                        d[l.get(c_amo)] = 100.0 * b / a
                return d

            A = taxa(linhas, c_amostra, c_in, c_merge)
            B = taxa(lb, ca_b, pb.get("entrada"), pb.get("merged"))
            comuns = sorted(set(A) & set(B))
            print("\n== 1b. Fusao: esta execucao x %s ==============" % args.comparar)
            if not comuns:
                print("  Nenhuma amostra em comum.")
            else:
                print("  %-28s %8s %8s %8s" % ("AMOSTRA", "ATUAL", "OUTRA", "DELTA"))
                deltas = []
                for a in sorted(comuns, key=lambda x: A[x]):
                    d = B[a] - A[a]
                    deltas.append(d)
                    print("  %-28s %7.1f%% %7.1f%% %+7.1f" % (a, A[a], B[a], d))
                med = mediana(deltas)
                print("\n  Delta mediano: %+.1f ponto(s) percentual(is)" % med)
                if med > 3:
                    print("  A outra execucao funde melhor — adote os parametros dela.")
                elif med < -3:
                    print("  Esta execucao funde melhor — mantenha os parametros atuais.")
                else:
                    print("  Empate tecnico: o parametro nao e a causa da perda.")

    # ------------------------------------------- 2. classificacao por nivel
    tax = achar(args.outdir, "ASV_tax_species*.tsv", "ASV_tax*.tsv", "*tax*species*.tsv")
    if tax:
        tl = ler_tsv(tax)
        ranks = [c for c in tl[0].keys()
                 if c.lower() in ("domain", "kingdom", "phylum", "class",
                                  "order", "family", "genus", "species")]
        print("\n== 2. Classificacao por nivel (%d ASVs) =========================" % len(tl))
        for r in ranks:
            n = sum(1 for l in tl
                    if l.get(r) and l[r].strip() not in ("", "NA", "NA_NA", "unclassified"))
            print("  %-9s %6d  %5.1f%%" % (r, n, 100.0 * n / len(tl)))
    else:
        print("\n== 2. Classificacao: tabela de taxonomia nao encontrada")

    # -------------------------------------- 3. construto do Smart Control
    tabela = achar(args.outdir, "ASV_table.tsv", "feature-table.tsv", "ASV_table*.tsv")
    if tabela:
        tl = ler_tsv(tabela)
        campos = list(tl[0].keys())
        c_id = campos[0]
        padrao = re.compile(args.controles)
        ctrl = [c for c in campos[1:] if padrao.search(c)]
        amostras = [c for c in campos[1:] if c not in ctrl]
        print("\n== 3. Smart Control: construto e carryover =====================")
        if not ctrl:
            print("  Nenhuma coluna casou com '%s' — controles nao entraram nesta execucao."
                  % args.controles)
        else:
            print("  Controles: %s" % ", ".join(ctrl))
            tot_ctrl = {c: sum(num(l.get(c)) or 0 for l in tl) for c in ctrl}
            ranking = []
            for l in tl:
                soma = sum(num(l.get(c)) or 0 for c in ctrl)
                ranking.append((soma, l[c_id], l))
            ranking.sort(reverse=True)
            print("\n  %-14s %12s %10s %14s" %
                  ("ASV", "nos controles", "% ctrl", "nas amostras"))
            base = sum(tot_ctrl.values()) or 1
            for soma, ident, l in ranking[:5]:
                nas = sum(num(l.get(c)) or 0 for c in amostras)
                print("  %-14s %12.0f %9.1f%% %14.0f"
                      % (ident[:14], soma, 100.0 * soma / base, nas))
            print("\n  O primeiro da lista e o construto sintetico. Remova-o das")
            print("  amostras antes de qualquer analise; o que sobra nos controles")
            print("  e contaminacao de reagente ou ambiente.")

    # ------------------------------- 4. cloroplasto e mitocondria
    # O ampliseq tem --exclude_taxa com default mitochondria,chloroplast, mas
    # o filtro nem sempre atua quando a taxonomia vem de banco customizado.
    # Conferir aqui evita levar sequencia de planta e de organela para dentro
    # da tabela de abundancia sem perceber.
    if tax and tabela:
        tl_tax = ler_tsv(tax)
        tl_tab = ler_tsv(tabela)
        campos_tab = list(tl_tab[0].keys())
        id_tab = campos_tab[0]
        amostras_tab = campos_tab[1:]
        id_tax = list(tl_tax[0].keys())[0]

        alvo = re.compile(r"chloroplast|mitochondri", re.I)
        marcados = {l[id_tax] for l in tl_tax
                    if alvo.search("\t".join(str(v) for v in l.values()))}

        total = 0.0
        removivel = 0.0
        for l in tl_tab:
            soma = sum(num(l.get(c)) or 0 for c in amostras_tab)
            total += soma
            if l[id_tab] in marcados:
                removivel += soma

        print("\n== 4. Cloroplasto e mitocondria ================================")
        print("  ASVs marcados : %d de %d" % (len(marcados), len(tl_tax)))
        print("  Reads neles   : %.0f de %.0f (%.2f%%)"
              % (removivel, total, 100.0 * removivel / total if total else 0))
        c_in_tf = coluna(list(linhas[0].keys()), "input_tax_filter")
        c_out_tf = coluna(list(linhas[0].keys()), "filtered_tax_filter")
        if c_in_tf and c_out_tf:
            ent = sum(num(l.get(c_in_tf)) or 0 for l in linhas)
            sai = sum(num(l.get(c_out_tf)) or 0 for l in linhas)
            print("  Filtro removeu: %.0f reads" % (ent - sai))
            # A ASV_table.tsv do DADA2 e a tabela ANTES do filtro; a filtrada
            # fica em qiime2/. Entao encontrar os ASVs marcados aqui e normal.
            # O que diz se o filtro atuou e a diferenca input/filtered_tax_filter.
            if marcados and ent == sai:
                print("  >>> O filtro do pipeline NAO atuou. Remova esses ASVs")
                print("      no downstream e reporte a fracao removida.")
            elif marcados and abs((ent - sai) - removivel) <= max(5, 0.02 * removivel):
                print("  Consistente: o filtro removeu exatamente esses ASVs.")
                print("  (a tabela do dada2/ e a de antes do filtro; use a de qiime2/)")

    # ------------------------------------------------------- onde olhar
    print("\n== Arquivos ====================================================")
    for nome in ("summary_report.html", "overall_summary.tsv",
                 "ASV_table*.tsv", "ASV_tax*.tsv"):
        p = achar(args.outdir, nome)
        if p:
            print("  %s" % p)


if __name__ == "__main__":
    main()

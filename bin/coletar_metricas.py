#!/usr/bin/env python3
"""
coletar_metricas.py — consolida as execucoes do ampliseq em final_reports/

Este script e a metade do contrato que fica do lado do cluster. Ele roda no
rieux, SO com biblioteca padrao — sem R, sem CRAN, sem versao de pacote para
brigar num ambiente que a gente nao controla — e escreve tabelas longas. O
desenho das figuras e a analise acontecem depois, no pacote R `aspp`, a partir
destes arquivos. Assim da para iterar no visual sem reprocessar nada.

E o mesmo papel que o awk/grep faz dentro do rieux_bacass: extrair metrica e
trabalho do wrapper, nao do pacote de visualizacao.

Saida (final_reports/)
----------------------
  reads_por_regiao.tsv        nome regiao amostra controle entrada filtrado
                              denoisedF denoisedR merged naochim
  classificacao.tsv           nome regiao rank n_asv n_classificado pct
  abundancia.tsv              nome regiao rank taxon reads_amostras rel_pct
                              n_amostras reads_controles
  comprimento_asv.tsv         nome regiao comprimento n_asv
  abundancia_por_amostra.tsv  nome regiao rank taxon amostra controle reads
                              pct_amostra
  prevalencia.tsv             nome regiao rank taxon n_presente n_total
                              pct_prevalencia pct_mediano pct_max amostra_max
                              n_controles_presente

As duas ultimas respondem uma pergunta que as agregadas nao respondem: um taxon
a 42% pode estar em todas as amostras ou empilhado numa so, e a conduta e
oposta nos dois casos.

Uso
---
    # layout do wrapper: <resultados>/<REGIAO>/
    coletar_metricas.py --resultados resultados/fabio --nome fabio \
                        --out resultados/fabio/final_reports

    # varios conjuntos de uma vez: <raiz>/resultados/<nome>/<REGIAO>/
    coletar_metricas.py --raiz .

    # so uma parte, e com foco num taxon
    coletar_metricas.py --resultados resultados/fabio --nome fabio \
                        --regioes ITS1 --foco 'antarctomyces|ochrolechia'

Compativel com Python 3.6, so biblioteca padrao.
"""

import argparse
import glob
import os
import re
import sys
from collections import defaultdict

RANKS = ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]
# O UNITE usa Kingdom; o SILVA que baixamos comeca em Domain. Aceitamos os dois.
ALIAS = {"Domain": "Kingdom"}

# diretorios que existem dentro da saida mas nao sao regioes
NAO_REGIAO = {"final_reports", "logs", "pipeline_info", "sidle", "multiqc"}

RANKS_DETALHE = ["Family", "Genus", "Species"]


# ------------------------------------------------------------------ leitura

def eh_controle(nome):
    return nome.lower().startswith("smart")


def ler_tsv(caminho):
    with open(caminho) as fh:
        cab = fh.readline().rstrip("\n").split("\t")
        linhas = [l.rstrip("\n").split("\t") for l in fh if l.strip()]
    return cab, linhas


def achar_tax(dir_dada2):
    """
    Prefere a tabela COM addSpecies, que e a que tem especie preenchida. O
    sufixo depois de ASV_tax varia com o id do banco (.user, .silva...), por
    isso o glob.
    """
    for padrao in ("ASV_tax_species*.tsv", "ASV_tax.*.tsv", "ASV_tax*.tsv"):
        achados = sorted(glob.glob(os.path.join(dir_dada2, padrao)))
        achados = [a for a in achados if "ITS" not in os.path.basename(a)]
        if achados:
            return achados[0]
    return None


def mediana(vals):
    if not vals:
        return 0.0
    v = sorted(vals)
    n = len(v)
    return v[n // 2] if n % 2 else 0.5 * (v[n // 2 - 1] + v[n // 2])


# ------------------------------------------------------------------ coleta

def retencao(dir_regiao, nome, regiao, saida):
    f = os.path.join(dir_regiao, "overall_summary.tsv")
    if not os.path.isfile(f):
        return 0
    cab, linhas = ler_tsv(f)
    idx = {n: i for i, n in enumerate(cab)}

    def pega(linha, campo):
        i = idx.get(campo)
        return linha[i] if i is not None and i < len(linha) else ""

    for l in linhas:
        amostra = l[0]
        saida.append([nome, regiao, amostra,
                      "sim" if eh_controle(amostra) else "nao",
                      pega(l, "DADA2_input"), pega(l, "filtered"),
                      pega(l, "denoisedF"), pega(l, "denoisedR"),
                      pega(l, "merged"), pega(l, "nonchim")])
    return len(linhas)


def comprimento(dir_regiao, nome, regiao, saida):
    f = os.path.join(dir_regiao, "dada2", "ASV_seqs.fasta")
    if not os.path.isfile(f):
        return 0
    hist = defaultdict(int)
    n = 0
    with open(f) as fh:
        for linha in fh:
            if linha.startswith(">"):
                if n:
                    hist[n] += 1
                n = 0
            else:
                n += len(linha.strip())
    if n:
        hist[n] += 1
    for c in sorted(hist):
        saida.append([nome, regiao, str(c), str(hist[c])])
    return len(hist)


def taxonomia(dir_regiao, nome, regiao, s_cls, s_abd, s_amo, s_prev,
              min_reads, rx_foco, focos):
    """
    Junta ASV_table.tsv (contagens) com a tabela de taxonomia e escreve as
    quatro saidas taxonomicas de uma passada so.
    """
    d = os.path.join(dir_regiao, "dada2")
    f_tab = os.path.join(d, "ASV_table.tsv")
    f_tax = achar_tax(d)
    if not os.path.isfile(f_tab) or not f_tax:
        return 0

    cab, linhas = ler_tsv(f_tab)
    amostras = cab[1:]
    n = len(amostras)
    ctrl = [eh_controle(a) for a in amostras]
    idx_am = [i for i, c in enumerate(ctrl) if not c]
    idx_ct = [i for i, c in enumerate(ctrl) if c]

    contagem = {}
    for l in linhas:
        vals = []
        for v in l[1:]:
            try:
                vals.append(float(v))
            except ValueError:
                vals.append(0.0)
        if len(vals) < n:
            vals += [0.0] * (n - len(vals))
        contagem[l[0]] = vals

    cab_t, linhas_t = ler_tsv(f_tax)
    pos = {}
    for i, campo in enumerate(cab_t):
        pos[ALIAS.get(campo, campo)] = i
    presentes = [r for r in RANKS if r in pos]
    if not presentes:
        return 0

    def valor(l, r):
        i = pos.get(r)
        if i is None:
            return ""
        v = l[i].strip() if i < len(l) else ""
        if not v or v.upper() == "NA" or "unidentified" in v.lower():
            return ""
        return v

    n_asv = 0
    classificado = defaultdict(int)
    acc = defaultdict(lambda: [0.0] * n)     # (rank, taxon) -> reads por amostra
    tot = dict((r, [0.0] * n) for r in presentes)

    for l in linhas_t:
        asv = l[0]
        if asv not in contagem:
            continue
        # Cloroplasto e mitocondria sao 16S de organela, nao da comunidade. O
        # ampliseq os remove na tabela do qiime2/, mas lemos a do dada2/, que e
        # anterior ao filtro — entao removemos aqui, com o mesmo criterio.
        # Cloroplasto e classificado DENTRO de Cyanobacteriota no SILVA: filtrar
        # por filo apagaria cianobacteria de verdade. O criterio e a Ordem.
        # Mitocondria e Familia, e isso importa: ela vive em Rickettsiales, a
        # mesma ordem de Anaplasmataceae — filtrar por ordem apagaria o alvo de
        # quem procura Anaplasma.
        if valor(l, "Order") == "Chloroplast" or valor(l, "Family") == "Mitochondria":
            continue
        n_asv += 1
        vals = contagem[asv]
        for r in presentes:
            v = valor(l, r)
            if not v:
                continue
            classificado[r] += 1
            # O SILVA guarda so o EPITETO na coluna Species ("acidifaciens").
            # Sem o genero a linha e ilegivel numa figura — e ambigua de
            # verdade: varios generos compartilham epiteto.
            if r == "Species":
                g = valor(l, "Genus")
                if g and not v.lower().startswith(g.lower()):
                    v = "%s %s" % (g, v)
            alvo = acc[(r, v)]
            alvo_tot = tot[r]
            for i, x in enumerate(vals):
                if x:
                    alvo[i] += x
                    alvo_tot[i] += x

    for r in presentes:
        pct = 100.0 * classificado[r] / n_asv if n_asv else 0.0
        s_cls.append([nome, regiao, r, str(n_asv), str(classificado[r]),
                      "%.1f" % pct])

    # --- agregado por regiao
    # A abundancia relativa e DENTRO da regiao: comparar contagem bruta entre
    # regioes seria comparar eficiencia de primer, nao composicao.
    soma_rank = defaultdict(float)
    for (r, _), vals in acc.items():
        soma_rank[r] += sum(vals[i] for i in idx_am)

    for (r, taxon), vals in sorted(acc.items()):
        reads_am = sum(vals[i] for i in idx_am)
        reads_ct = sum(vals[i] for i in idx_ct)
        n_pres = sum(1 for i in idx_am if vals[i] > 0)
        if reads_am > 0:
            rel = 100.0 * reads_am / soma_rank[r] if soma_rank[r] else 0.0
            s_abd.append([nome, regiao, r, taxon, "%.0f" % reads_am,
                          "%.6f" % rel, str(n_pres), "%.0f" % reads_ct])

        if r not in RANKS_DETALHE:
            continue

        # --- por amostra e prevalencia
        pcts, pct_max, amostra_max = [], 0.0, ""
        for i in idx_am:
            if vals[i] <= 0:
                continue
            den = tot[r][i]
            p = 100.0 * vals[i] / den if den else 0.0
            pcts.append(p)
            if p > pct_max:
                pct_max, amostra_max = p, amostras[i]
        n_ct = sum(1 for i in idx_ct if vals[i] > 0)
        if n_pres == 0 and n_ct == 0:
            continue

        s_prev.append([nome, regiao, r, taxon, str(n_pres), str(len(idx_am)),
                       "%.1f" % (100.0 * n_pres / len(idx_am) if idx_am else 0.0),
                       # mediana ENTRE AS AMOSTRAS EM QUE OCORRE. A mediana
                       # sobre todas seria zero para qualquer taxon esparso e
                       # esconderia exatamente o caso que interessa.
                       "%.3f" % mediana(pcts), "%.3f" % pct_max, amostra_max,
                       str(n_ct)])

        for i, v in enumerate(vals):
            if v <= 0 or v < min_reads:
                continue
            den = tot[r][i]
            s_amo.append([nome, regiao, r, taxon, amostras[i],
                          "sim" if ctrl[i] else "nao", "%.0f" % v,
                          "%.3f" % (100.0 * v / den if den else 0.0)])

        if rx_foco and rx_foco.search(taxon):
            det = [(amostras[i], ctrl[i], vals[i],
                    100.0 * vals[i] / tot[r][i] if tot[r][i] else 0.0)
                   for i in range(n) if vals[i] > 0]
            focos.append((nome, regiao, r, taxon, n_pres, len(idx_am),
                          len(idx_ct), det))

    return n_asv


# ------------------------------------------------------------------ principal

def regioes_em(dirbase):
    return sorted(d for d in os.listdir(dirbase)
                  if os.path.isdir(os.path.join(dirbase, d))
                  and d not in NAO_REGIAO and not d.startswith("."))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--resultados", help="diretorio com uma subpasta por regiao")
    ap.add_argument("--nome", default=None, help="rotulo do conjunto")
    ap.add_argument("--raiz", default=None,
                    help="alternativa: <raiz>/resultados/<nome>/<REGIAO>/")
    ap.add_argument("--out", default=None, help="diretorio de saida")
    ap.add_argument("--regioes", nargs="*", default=None)
    ap.add_argument("--min-reads", type=int, default=1)
    ap.add_argument("--foco", default=None,
                    help="regex; imprime o detalhe por amostra dos taxons que casarem")
    args = ap.parse_args()

    # --- que conjuntos processar: [(nome, dir)]
    conjuntos = []
    if args.resultados:
        nome = args.nome or os.path.basename(os.path.normpath(args.resultados))
        conjuntos.append((nome, args.resultados))
        saida_padrao = os.path.join(args.resultados, "final_reports")
    elif args.raiz:
        base = os.path.join(args.raiz, "resultados")
        if not os.path.isdir(base):
            sys.exit("ERRO: nao achei %s" % base)
        for d in sorted(os.listdir(base)):
            if os.path.isdir(os.path.join(base, d)):
                conjuntos.append((d, os.path.join(base, d)))
        saida_padrao = os.path.join(args.raiz, "final_reports")
    else:
        sys.exit("ERRO: informe --resultados ou --raiz")

    out = args.out or saida_padrao
    rx = re.compile(args.foco, re.I) if args.foco else None

    s_ret, s_cls, s_abd, s_len, s_amo, s_prev = [], [], [], [], [], []
    focos = []

    for nome, dirbase in conjuntos:
        regs = args.regioes if args.regioes else regioes_em(dirbase)
        for regiao in regs:
            d = os.path.join(dirbase, regiao)
            if not os.path.isdir(d):
                continue
            n1 = retencao(d, nome, regiao, s_ret)
            n2 = taxonomia(d, nome, regiao, s_cls, s_abd, s_amo, s_prev,
                           args.min_reads, rx, focos)
            comprimento(d, nome, regiao, s_len)
            marca = "" if n2 else "   (sem taxonomia — execucao incompleta)"
            print("  %-12s %-12s %3d amostras  %5d ASVs%s"
                  % (nome, regiao, n1, n2, marca))

    os.makedirs(out, exist_ok=True)

    def gravar(arquivo, cabecalho, linhas):
        caminho = os.path.join(out, arquivo)
        with open(caminho, "w") as fh:
            fh.write("\t".join(cabecalho) + "\n")
            for l in linhas:
                fh.write("\t".join(l) + "\n")
        print("  %-30s %8d linhas" % (arquivo, len(linhas)))

    print("\nGravando em %s/" % out)
    gravar("reads_por_regiao.tsv",
           ["nome", "regiao", "amostra", "controle", "entrada", "filtrado",
            "denoisedF", "denoisedR", "merged", "naochim"], s_ret)
    gravar("classificacao.tsv",
           ["nome", "regiao", "rank", "n_asv", "n_classificado", "pct"], s_cls)
    gravar("abundancia.tsv",
           ["nome", "regiao", "rank", "taxon", "reads_amostras", "rel_pct",
            "n_amostras", "reads_controles"], s_abd)
    gravar("comprimento_asv.tsv",
           ["nome", "regiao", "comprimento", "n_asv"], s_len)
    gravar("abundancia_por_amostra.tsv",
           ["nome", "regiao", "rank", "taxon", "amostra", "controle", "reads",
            "pct_amostra"], s_amo)
    gravar("prevalencia.tsv",
           ["nome", "regiao", "rank", "taxon", "n_presente", "n_total",
            "pct_prevalencia", "pct_mediano", "pct_max", "amostra_max",
            "n_controles_presente"], s_prev)

    if rx:
        print("\n" + "=" * 74)
        print("FOCO: /%s/" % args.foco)
        print("=" * 74)
        if not focos:
            print("  nenhum taxon casou com o padrao.")
        for nome, regiao, r, taxon, n_pres, n_total, n_ct, det in focos:
            print("\n%s · %s · %s · %s" % (nome, regiao, r, taxon))
            print("  em %d de %d amostras (%.0f%%)%s"
                  % (n_pres, n_total, 100.0 * n_pres / n_total if n_total else 0.0,
                     "" if not n_ct else "  ·  %d controles na corrida" % n_ct))
            for amostra, c, v, p in sorted(det, key=lambda x: -x[3]):
                print("    %-24s %10.0f reads  %7.2f%%%s"
                      % (amostra, v, p, "   <-- CONTROLE" if c else ""))

    print("\nNo R:  dados <- aspp::read_ampliseq_summary('%s')" % out)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
organize_project.py

Cruza a planilha de amostras com os dados em disco e monta, para cada
colaborador, um diretorio proprio com os FASTQ identificados pelo ID
original, mais a tabela de metadados correspondente.

DECISOES DE DESENHO, e o porque:

  Symlink, nao copia. Sao ~130 GB de FASTQ. Duplicar isso por colaborador
  nao traz beneficio nenhum: o cutadapt e o Nextflow leem atraves de link
  sem reclamar (ja verificado), e uma copia so cria a chance de as duas
  divergirem. Use --copiar se precisar mesmo de arquivos independentes.

  O nome do arquivo carrega os DOIS identificadores:
      <ID_original>__s<ID_sequenciamento>_R1.fastq.gz
  Renomear so para o ID original perderia dados: a Maria tem RM18, RM46,
  RM48 e RM55 duas vezes cada, e o Herbert tem "WT S" e "WT-S" que sao a
  mesma amostra escrita de dois jeitos. Sem o ID de sequenciamento no nome,
  um link sobrescreveria o outro em silencio.

  Nomes sao saneados para o formato que o nf-core/ampliseq exige (comeca
  com letra, so letras, digitos e underscore), mas o ID original intacto
  fica preservado na coluna id_original da tabela de metadados.

Uso:
    python3 organize_project.py planilha.csv dir_dados_brutos dir_saida
    python3 organize_project.py ... --executar     # cria os links
    python3 organize_project.py ... --copiar       # copia em vez de linkar

Sem --executar ele so relata; nada e criado.

Compativel com Python 3.6.
"""

import argparse
import csv
import os
import re
import shutil
import sys
from collections import Counter, defaultdict


def sanear(nome):
    """Nome aceito pelo ampliseq: comeca com letra, [A-Za-z0-9_]."""
    limpo = re.sub(r"[^A-Za-z0-9]+", "_", nome.strip()).strip("_")
    limpo = re.sub(r"_+", "_", limpo)
    if not limpo:
        limpo = "amostra"
    if not limpo[0].isalpha():
        limpo = "s" + limpo
    return limpo


def dono_base(valor):
    """'herbert - Fezes N Experimento 12' -> ('herbert', 'Fezes N Experimento 12')"""
    bruto = valor.strip()
    if " - " in bruto:
        dono, grupo = bruto.split(" - ", 1)
    else:
        dono, grupo = bruto, ""
    return dono.strip().lower(), grupo.strip()


def indexar_disco(raiz):
    """ID de sequenciamento -> (R1, R2). O diretorio do BaseSpace e
    <ID>_ds.<hash> e o arquivo <ID>_S<n>_L<n>_R1_001.fastq.gz."""
    achados = {}
    for dirpath, _dirs, arquivos in os.walk(raiz, followlinks=True):
        for a in arquivos:
            if "_R1_" not in a or not a.endswith((".fastq.gz", ".fq.gz")):
                continue
            r1 = os.path.join(dirpath, a)
            r2 = r1.replace("_R1_", "_R2_")
            if not os.path.exists(r2):
                continue
            ident = a.split("_S")[0]
            achados[ident] = (r1, r2)
    return achados


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("sample_table")
    ap.add_argument("raw_dir")
    ap.add_argument("out_dir")
    ap.add_argument("--apply", action="store_true",
                    help="actually create the links/copies (default: report only)")
    ap.add_argument("--copy", action="store_true",
                    help="copy the files instead of creating symlinks")
    ap.add_argument("--controls", default="^Smart",
                    help="regex of on-disk samples that are controls. They "
                         "nao estao na planilha e sao linkadas em TODAS as "
                         "pastas, marcadas na coluna 'controle'. Passe uma "
                         "string vazia para nao incluir nenhuma.")
    args = ap.parse_args()

    with open(args.sample_table, encoding="utf-8-sig") as fh:
        linhas = list(csv.DictReader(fh))
    cols = list(linhas[0].keys())
    C_DONO, C_ORIG, C_SEQ = cols[0], cols[1], cols[2]

    disco = indexar_disco(args.raw_dir)
    print("Sample table: %d rows | On disk: %d samples with an R1/R2 pair\n"
          % (len(linhas), len(disco)))

    # ---------------------------------------------------------- cruzamento
    por_dono = defaultdict(list)
    sem_dado = defaultdict(list)
    for r in linhas:
        dono, grupo = dono_base(r[C_DONO])
        seq = r[C_SEQ].strip()
        orig = r[C_ORIG].strip()
        if seq in disco:
            por_dono[dono].append((seq, orig, grupo))
        else:
            sem_dado[dono].append((seq, orig))

    em_disco_sem_planilha = sorted(
        set(disco) - {r[C_SEQ].strip() for r in linhas},
        key=lambda x: (not x.isdigit(), x))

    # Controles entram em TODAS as pastas. Cada execucao do pipeline precisa
    # deles por duas razoes independentes do lote: o construto sintetico tem
    # sitio para os 7 pares de primers, entao valida o split; e o ASV dele
    # identifica o carryover que aparece nas amostras reais.
    controles = []
    if args.controls:
        padrao = re.compile(args.controls)
        controles = sorted(s for s in em_disco_sem_planilha if padrao.search(s))

    # ---------------------------------------------------------- relatorio
    print("%-10s %6s %6s %s" % ("DONO", "COM", "SEM", "IDs sem dado"))
    print("-" * 74)
    for dono in sorted(set(list(por_dono) + list(sem_dado))):
        faltam = sem_dado.get(dono, [])
        com = len(por_dono.get(dono, []))
        alerta = ""
        if faltam and com + len(faltam) > 0:
            perda = 100.0 * len(faltam) / (com + len(faltam))
            if perda >= 25:
                alerta = "  <<< %.0f%% de perda" % perda
        print("%-10s %6d %6d %s%s"
              % (dono, com, len(faltam),
                 ", ".join(s for s, _ in faltam[:8]) + ("..." if len(faltam) > 8 else ""),
                 alerta))

    if em_disco_sem_planilha:
        print("\nOn disk with no row in the sample table: %s"
              % ", ".join(em_disco_sem_planilha))
    if controles:
        print("Treated as controls (they go into every folder): %s"
              % ", ".join(controles))
        ignorados = [s for s in em_disco_sem_planilha if s not in controles]
        if ignorados:
            print("Ignored (neither sample nor control): %s" % ", ".join(ignorados))

    # ------------------------------------------- replicatas (ID repetido)
    print("\n== Candidatos a replicata ======================================")
    print("Mesmo ID original no mesmo colaborador, ou nomes que so diferem")
    print("em pontuacao (ex.: 'WT S' e 'WT-S').\n")
    grupos = defaultdict(list)
    for dono, itens in por_dono.items():
        for seq, orig, _g in itens:
            grupos[(dono, sanear(orig).lower())].append((seq, orig))
    n_rep = 0
    for (dono, chave), membros in sorted(grupos.items()):
        if len(membros) > 1:
            n_rep += 1
            formas = sorted(set(o for _s, o in membros))
            marca = "  (grafias diferentes: %s)" % " | ".join(formas) if len(formas) > 1 else ""
            print("  %-9s %-28s -> %s%s"
                  % (dono, chave, ", ".join(s for s, _o in membros), marca))
    print("\n  %d grupo(s) de replicata." % n_rep)

    # ---------------------------------------------------------- execucao
    print("\n== Saida =======================================================")
    for dono, itens in sorted(por_dono.items()):
        destino = os.path.join(args.out_dir, dono, "raw")
        meta = os.path.join(args.out_dir, dono, "metadata.tsv")
        if args.apply:
            os.makedirs(destino, exist_ok=True)

        usados = Counter()
        linhas_meta = []
        fila = sorted(itens, key=lambda x: int(x[0]) if x[0].isdigit() else 0)
        fila += [(c, c, "controle") for c in controles]
        for seq, orig, grupo in fila:
            eh_controle = grupo == "controle"
            base = sanear(orig) if eh_controle else "%s__s%s" % (sanear(orig), seq)
            usados[base] += 1
            if usados[base] > 1:          # nao deveria ocorrer, mas garante
                base = "%s_%d" % (base, usados[base])
            r1, r2 = disco[seq]
            for origem, sufixo in ((r1, "R1"), (r2, "R2")):
                alvo = os.path.join(destino, "%s_%s.fastq.gz" % (base, sufixo))
                if not args.apply:
                    continue
                if os.path.lexists(alvo):
                    os.remove(alvo)
                if args.copy:
                    shutil.copy2(origem, alvo)
                else:
                    os.symlink(os.path.abspath(origem), alvo)
            chave = sanear(orig).lower()
            rep = (not eh_controle) and len(grupos[(dono, chave)]) > 1
            linhas_meta.append([base, seq, orig, dono, grupo or "-",
                                "yes" if rep else "no",
                                "yes" if eh_controle else "no",
                                os.path.abspath(os.path.join(destino, base + "_R1.fastq.gz")),
                                os.path.abspath(os.path.join(destino, base + "_R2.fastq.gz"))])

        if args.apply:
            with open(meta, "w") as fh:
                fh.write("sample\tseq_id\toriginal_id\tproject"
                         "\tgroup\treplicate\tcontrol\tfastq_1\tfastq_2\n")
                for l in linhas_meta:
                    fh.write("\t".join(l) + "\n")
        print("  %-10s %3d samples + %d control(s) -> %s"
              % (dono, len(itens), len(controles), destino))

    if not args.apply:
        print("\n(nothing was created — repeat with --apply)")


if __name__ == "__main__":
    main()

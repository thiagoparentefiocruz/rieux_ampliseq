#!/usr/bin/env python3
"""
make_samplesheet.py — samplesheet do ampliseq a partir do metadata.tsv

Para que serve
--------------
O ramo por regiao usa os samplesheets do split, um por regiao. O ramo Sidle
nao: ele recebe as reads INTEIRAS e faz o proprio roteamento, rodando cutadapt
uma vez por regiao com os primers do regions_multiregion.tsv. Alimentar o Sidle
com as reads ja divididas cortaria primer duas vezes e deslocaria as bordas.

O metadata.tsv que o organize_project.py grava ja tem o que falta:
`sample`, `fastq_1` e `fastq_2` apontando para as reads nao divididas. Este
script so projeta essas colunas no formato do ampliseq.

Verificacao que ele faz e' o motivo de existir como script
----------------------------------------------------------
Os FASTQ sao symlinks (sao ~130 GB; duplicar por colaborador nao traz
beneficio). Symlink quebra em silencio quando o alvo se move, e o Nextflow so
descobre isso na primeira tarefa, horas dentro da execucao. Aqui o custo de
checar e' um os.path.exists por arquivo, e a falha aparece antes de submeter.

Uso
---
    make_samplesheet.py colaboradores/renata/metadata.tsv \\
        --out split/samplesheets/renata/samplesheet_completo.tsv

    # varios de uma vez, um arquivo por colaborador
    make_samplesheet.py colaboradores/*/metadata.tsv --out-dir samplesheets/

Compativel com Python 3.6, so biblioteca padrao.
"""

import argparse
import os
import sys


def ler(caminho):
    with open(caminho) as fh:
        cab = fh.readline().rstrip("\n").split("\t")
        idx = {n: i for i, n in enumerate(cab)}
        for obrig in ("sample", "fastq_1", "fastq_2"):
            if obrig not in idx:
                sys.exit("ERROR: %s has no '%s' column — this script expects the "
                         "metadata.tsv written by organize_project.py"
                         % (caminho, obrig))
        linhas = []
        for l in fh:
            if not l.strip():
                continue
            c = l.rstrip("\n").split("\t")
            if len(c) < len(cab):
                continue
            linhas.append(dict((n, c[i]) for n, i in idx.items()))
    return linhas


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("metadata", nargs="+")
    ap.add_argument("--out", default=None,
                    help="output file (only with a single input metadata)")
    ap.add_argument("--out-dir", default=None,
                    help="directory; writes samplesheet_complete_<name>.tsv")
    ap.add_argument("--no-controls", action="store_true",
                    help="drop the Smart Controls. By default they are INCLUDED: the "
                         "synthetic construct carries a site for all seven primer "
                         "pairs, so it is what validates the per-region routing "
                         "that Sidle does on its own.")
    ap.add_argument("--no-check", action="store_true",
                    help="do not check that the FASTQs exist")
    args = ap.parse_args()

    if not args.out and not args.out_dir:
        sys.exit("ERROR: give --out or --out-dir")
    if args.out and len(args.metadata) > 1:
        sys.exit("ERROR: --out only works for one metadata; use --out-dir")

    total_falta = 0
    for caminho in args.metadata:
        linhas = ler(caminho)
        nome = os.path.basename(os.path.dirname(os.path.abspath(caminho)))

        escolhidas, faltando, controles = [], [], 0
        for r in linhas:
            eh_ctrl = r.get("control", r.get("controle", "no")).strip().lower() in ("yes", "sim")
            if eh_ctrl:
                controles += 1
                if args.no_controls:
                    continue
            r1, r2 = r["fastq_1"], r["fastq_2"]
            if not args.no_check and not (os.path.exists(r1) and os.path.exists(r2)):
                faltando.append((r["sample"], r1 if not os.path.exists(r1) else r2))
                continue
            escolhidas.append((r["sample"], r1, r2))

        if args.out:
            saida = args.out
        else:
            os.makedirs(args.out_dir, exist_ok=True)
            saida = os.path.join(args.out_dir,
                                 "samplesheet_complete_%s.tsv" % nome)
        pasta = os.path.dirname(os.path.abspath(saida))
        if pasta:
            os.makedirs(pasta, exist_ok=True)

        with open(saida, "w") as fh:
            fh.write("sampleID\tforwardReads\treverseReads\n")
            for s, r1, r2 in escolhidas:
                fh.write("%s\t%s\t%s\n" % (s, r1, r2))

        print("%-12s %3d samples (%d control%s) -> %s"
              % (nome, len(escolhidas), controles,
                 "" if controles == 1 else "s", saida))
        if faltando:
            total_falta += len(faltando)
            print("   %d missing FASTQ(s) — broken symlink or moved file:"
                  % len(faltando))
            for s, f in faltando[:5]:
                print("     %-22s %s" % (s, f))
            if len(faltando) > 5:
                print("     ... and %d more" % (len(faltando) - 5))

    if total_falta:
        print("\nWARNING: %d missing file(s). The matching samples were LEFT OUT"
              % total_falta)
        print("of the samplesheet — fix the links before running,")
        print("or the run goes ahead with fewer samples than you think.")
        sys.exit(1)


if __name__ == "__main__":
    main()

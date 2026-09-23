#!/usr/bin/env python3
"""
diagnose_merge.py — por que uma regiao perdeu reads, e em que passo.

    bin/diagnose_merge.py --reports <projeto>/final_reports \
                          --params  <projeto>/region_params.tsv

Uma regiao com fusao baixa tem poucas causas possiveis, e elas pedem
correcoes opostas. Olhar so a porcentagem final nao distingue:

  - o filtro de qualidade comeu as reads          -> truncar mais curto
  - a fusao nao fecha porque o truncLen e curto   -> truncar mais longo
  - a quimera comeu os ASVs                       -> nao e truncagem

O sinal que separa o segundo caso dos outros esta no comprimento dos ASVs.
O DADA2 so funde um par quando truncLenF + truncLenR - inserto >= minOverlap
(12 por padrao). Entao existe um TETO de comprimento fusionavel:

    teto = truncLenF + truncLenR - 12

Amplicon mais longo que isso nao aparece na saida — nao porque nao exista,
mas porque nao teve como fundir. O histograma de comprimento denuncia: em vez
de cair suave, ele bate numa parede e para. Uma pilha de ASVs encostada no
teto, com fusao baixa, e diagnostico de truncagem curta demais.

Cuidado de leitura que vale registrar: o asv_length.tsv so contem o que
fundiu. A distribuicao observada ja e a distribuicao CENSURADA pelo teto, e
por isso ela nunca vai "mostrar" o que ficou de fora. Quem estima o que
faltou e o PCR in-silico do validate_sidle_regions.py, contra o banco.
"""
import argparse
import os
import sys
from collections import defaultdict

SOBREPOSICAO_DADA2 = 12   # minOverlap padrao do DADA2
PERTO_DO_TETO = 3         # bp: o que conta como "encostado na parede"
LARGURA_ESTREITA = 12     # bp: p95-p5 abaixo disso nao e distribuicao, e fatia


def ler_tsv(caminho):
    """Linhas como dicionarios, ignorando comentarios e linhas curtas."""
    linhas = []
    with open(caminho) as fh:
        cab = None
        for linha in fh:
            if linha.startswith("#") or not linha.strip():
                continue
            c = linha.rstrip("\n").split("\t")
            if cab is None:
                cab = c
                continue
            if len(c) < len(cab):
                c += [""] * (len(cab) - len(c))
            linhas.append(dict(zip(cab, c)))
    return linhas


def percentil(hist, p):
    """hist: {comprimento: n}. Devolve o percentil p (0-100)."""
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
        description="Localize where each region lost its reads.")
    ap.add_argument("--reports", required=True,
                    help="final_reports/ directory written by collect_metrics.py")
    ap.add_argument("--params", default=None,
                    help="region_params.tsv (truncLen per region). Without it, "
                         "the merge ceiling cannot be computed.")
    ap.add_argument("--min-merge", type=float, default=85.0,
                    help="merge rate below which a region is reported (default 85)")
    ap.add_argument("--overlap", type=int, default=SOBREPOSICAO_DADA2,
                    help="DADA2 minOverlap (default 12)")
    args = ap.parse_args()

    f_reads = os.path.join(args.reports, "reads_per_region.tsv")
    f_len = os.path.join(args.reports, "asv_length.tsv")
    for f in (f_reads, f_len):
        if not os.path.isfile(f):
            sys.exit("ERROR: %s not found — run the 'collect' stage first." % f)

    # ---------------------------------------------------- retencao por passo
    passos = ["input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim"]
    soma = defaultdict(lambda: defaultdict(int))
    # Separar controle de amostra nao e detalhe de apresentacao: e o unico
    # corte que distingue "o parametro esta errado" de "a amostra tem outra
    # coisa dentro". O controle e uma comunidade conhecida; se ele funde e a
    # amostra nao, mexer em parametro nao vai resolver.
    por_tipo = {"sample": defaultdict(lambda: defaultdict(int)),
                "control": defaultdict(lambda: defaultdict(int))}
    amostras = defaultdict(int)
    controles = defaultdict(int)
    for l in ler_tsv(f_reads):
        r = l["region"]
        ctrl = (l.get("control", "no") or "no").strip().lower() in ("yes", "sim", "1", "true")
        if ctrl:
            controles[r] += 1
        else:
            amostras[r] += 1
        for p in passos:
            try:
                v = int(float(l.get(p) or 0))
            except ValueError:
                continue
            soma[r][p] += v
            por_tipo["control" if ctrl else "sample"][r][p] += v

    # ------------------------------------------------ histograma por regiao
    hist = defaultdict(lambda: defaultdict(int))
    for l in ler_tsv(f_len):
        try:
            hist[l["region"]][int(l["length"])] += int(l["n_asv"])
        except (ValueError, KeyError):
            pass

    # ------------------------------------------------------------ truncLens
    trunc = {}
    if args.params:
        if not os.path.isfile(args.params):
            sys.exit("ERROR: %s not found" % args.params)
        for l in ler_tsv(args.params):
            try:
                tF, tR = int(l["trunclenf"]), int(l["trunclenr"])
            except (ValueError, KeyError):
                continue
            if tF and tR:
                trunc[l["region"]] = (tF, tR)

    # ------------------------------------------------------------- relatorio
    def taxa_fusao(acum, r):
        """merged / denoised, ou None quando nao ha nada daquele tipo."""
        s = acum.get(r)
        if not s:
            return None
        den = min(s.get("denoisedF", 0), s.get("denoisedR", 0))
        if not den:
            return None
        return 100.0 * s.get("merged", 0) / den

    print("Read retention per region (% of the DADA2 input):\n")
    print("  %-7s %5s %5s %8s %8s %8s %8s %9s %9s"
          % ("region", "n", "ctrl", "filter", "denoise", "merge",
             "chimera", "merge_smp", "merge_ctl"))
    print("  " + "-" * 74)
    for r in sorted(soma):
        s = soma[r]
        ent = s["input"] or 1
        den = min(s["denoisedF"], s["denoisedR"])
        t_amo = taxa_fusao(por_tipo["sample"], r)
        t_ctl = taxa_fusao(por_tipo["control"], r)
        print("  %-7s %5d %5d %7.1f%% %7.1f%% %7.1f%% %7.1f%% %8s %9s"
              % (r, amostras[r], controles[r],
                 100.0 * s["filtered"] / ent,
                 100.0 * den / (s["filtered"] or 1),
                 100.0 * s["merged"] / (den or 1),
                 100.0 * s["nonchim"] / (s["merged"] or 1),
                 "-" if t_amo is None else "%.1f%%" % t_amo,
                 "-" if t_ctl is None else "%.1f%%" % t_ctl))

    print("\nEach column is the fraction that SURVIVED that step, not the")
    print("cumulative total. The last two split the merge between real")
    print("samples and controls — the cut that separates a wrong parameter")
    print("from a sample that contains something the control does not.")

    print("\n\nMerge ceiling vs observed ASV length:\n")
    print("  %-7s %7s %7s %8s %6s %6s %6s %6s %7s"
          % ("region", "truncF", "truncR", "ceiling",
             "p5", "p50", "p95", "max", "at_wall"))
    print("  " + "-" * 66)
    veredito = []
    for r in sorted(hist):
        h = hist[r]
        p5, p50 = percentil(h, 5), percentil(h, 50)
        p95, mx = percentil(h, 95), max(h) if h else 0
        largura = p95 - p5          # largura da distribuicao QUE SOBREVIVEU
        if r in trunc:
            tF, tR = trunc[r]
            teto = tF + tR - args.overlap
            total = sum(h.values()) or 1
            parede = sum(n for v, n in h.items() if v >= teto - PERTO_DO_TETO)
            pct_parede = 100.0 * parede / total
            print("  %-7s %7d %7d %8d %6d %6d %6d %6d %6.1f%%"
                  % (r, tF, tR, teto, p5, p50, p95, mx, pct_parede))
        else:
            teto = pct_parede = None
            print("  %-7s %7s %7s %8s %6d %6d %6d %6d %7s"
                  % (r, "-", "-", "-", p5, p50, p95, mx, "-"))

        s = soma.get(r, {})
        den = min(s.get("denoisedF", 0), s.get("denoisedR", 0))
        taxa = 100.0 * s.get("merged", 0) / (den or 1)
        if taxa >= args.min_merge:
            continue

        # O controle e uma comunidade conhecida e de comprimento conhecido.
        # Se ele funde e a amostra nao, o que difere nao e parametro: e o que
        # esta dentro da amostra. Este teste vem ANTES dos de comprimento
        # porque nenhum truncLen conserta conteudo.
        t_amo = taxa_fusao(por_tipo["sample"], r)
        t_ctl = taxa_fusao(por_tipo["control"], r)
        if (t_ctl is not None and t_amo is not None
                and t_ctl >= args.min_merge and t_ctl - t_amo >= 30.0):
            veredito.append((r, taxa,
                "NOT A PARAMETER — the controls merge and the samples do not. "
                "Controls %.1f%%, samples %.1f%%, same region, same truncLen, "
                "same run. A known community pairs normally here, so the "
                "pipeline is doing its job; what the samples carry is "
                "something whose amplicon these reads cannot span. The usual "
                "cause is off-target amplification: several universal 16S "
                "pairs also amplify eukaryotic 18S, which is far longer than "
                "the 16S product and therefore never merges at any truncLen. "
                "Raising truncLen will not recover it — confirm what the "
                "unmerged pairs are (vsearch --fastq_mergepairs reports WHY "
                "each pair failed) before treating this as a loss."
                % (t_ctl, t_amo)))
            continue

        if teto is None:
            veredito.append((r, taxa,
                             "merge %.1f%% — give --params to test the ceiling"
                             % taxa))
        elif mx >= teto - PERTO_DO_TETO or pct_parede >= 10.0:
            veredito.append((r, taxa,
                "CAPPED BY truncLen, with a pile-up at the wall. The histogram "
                "stops at the ceiling (%d bp) with %.1f%% of the ASVs against "
                "it, and only %.1f%% of the pairs merged. Anything longer than "
                "%d bp had no way to merge. Raise truncLenF+truncLenR (read "
                "length permitting) and re-run, or accept that this region "
                "only sees its short end."
                % (teto, pct_parede, taxa, teto)))
        elif largura <= LARGURA_ESTREITA:
            # Sem pilha na parede, mas o que sobrou e uma FATIA de comprimento,
            # nao uma distribuicao. Uma comunidade cujo comprimento e bimodal
            # some inteira acima do teto e nao deixa pilha nenhuma: some o
            # modo longo, fica o curto, e o histograma parece bem-comportado.
            # Ausencia de parede NAO e ausencia de censura.
            veredito.append((r, taxa,
                "INCONCLUSIVE from the lengths alone, and the shape is "
                "suspicious: only %.1f%% of the pairs merged, and everything "
                "that did merge sits in a %d bp window (p5 %d, p95 %d) well "
                "below the %d bp ceiling. A community whose lengths are "
                "bimodal disappears above the ceiling WITHOUT leaving a "
                "pile-up — the long mode vanishes and the short one looks "
                "tidy. Two ways to settle it: in-silico PCR against the "
                "reference (validate_sidle_regions.py), which gives the "
                "UNCENSORED expected distribution, and one test run of this "
                "region at the maximum truncLen the read length allows. If "
                "merge jumps, the ceiling was binding after all."
                % (taxa, largura, p5, p95, teto)))
        else:
            veredito.append((r, taxa,
                "merge %.1f%%, and the surviving lengths span %d bp (p5 %d, "
                "p95 %d, max %d) with room below the %d bp ceiling. "
                "Truncation does not look like the binding constraint — look "
                "at the filter and denoise columns above, at primer "
                "carryover, or at the region simply not amplifying here."
                % (taxa, largura, p5, p95, mx, teto)))

    if not veredito:
        print("\nEvery region merged above %.0f%%. Nothing to diagnose."
              % args.min_merge)
        return 0

    print("\n\n" + "=" * 70)
    print("REGIONS BELOW %.0f%% MERGE" % args.min_merge)
    print("=" * 70)
    for r, taxa, texto in sorted(veredito, key=lambda x: x[1]):
        print("\n%s  (merge %.1f%%)" % (r, taxa))
        linha = ""
        for palavra in texto.split():
            if len(linha) + len(palavra) + 1 > 68:
                print("  " + linha)
                linha = palavra
            else:
                linha = (linha + " " + palavra).strip()
        if linha:
            print("  " + linha)
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except BrokenPipeError:
        # sair de `| head` sem despejar traceback na cara de quem so queria
        # ver as primeiras linhas
        try:
            sys.stdout.close()
        finally:
            sys.exit(0)

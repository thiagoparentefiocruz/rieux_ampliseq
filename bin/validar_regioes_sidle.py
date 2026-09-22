#!/usr/bin/env python3
"""
validar_regioes_sidle.py — PCR in silico para checar os primers e escolher
                           o region_length do --multiregion

Por que este script existe
--------------------------
O `--multiregion` do ampliseq pede um TSV com quatro colunas por regiao:

    region   region_length   FW_primer   RV_primer

Os primers a gente ja tem (primers_painel.tsv, recuperados dos proprios FASTQ).
Mas eles foram recortados para um trabalho DIFERENTE: rotear read para regiao.
Ali, primer mais curto e mais permissivo e' seguro — no maximo classifica um
pouco a mais e a validacao empirica corrige.

No Sidle o mesmo primer faz outra coisa: ele RECORTA O BANCO DE REFERENCIA.
Um primer curto demais casa em lugar errado e produz uma regiao de referencia
deslocada — e o Sidle nao acusa erro, so reconstroi mal. O V7V9 e' o caso
explicito: a nota do primers_painel.tsv diz que ali o nucleo e' a maior
substring comum, nao o prefixo, entao a borda 5' dele e' incerta por
construcao.

E o `region_length` e' pior ainda, porque e' o mesmo mecanismo que ja nos
custou o V1V2: a documentacao diz que as sequencias sao cortadas nesse
comprimento e "shorter ones are omitted". Chutar esse numero e' repetir o erro
do truncLen — com a diferenca de que dessa vez a gente TEM a distribuicao de
comprimento observada.

O que o script faz
------------------
1. Extrai cada regiao do banco de referencia por PCR in silico, com os primers
   do painel (casamento IUPAC, tolerando ate N erros, nunca nas 3 ultimas bases
   do primer — que e' onde a polimerase de fato exige pareamento).
2. Compara a distribuicao de comprimento extraida da REFERENCIA com a
   distribuicao dos ASVs OBSERVADOS (comprimento_asv.tsv).
   - As duas medianas baterem e' a evidencia de que os primers e as bordas
     estao certos.
   - Uma diferenca aproximadamente constante denuncia deslocamento de borda, e
     diz de quantas bases.
   - Taxa de recuperacao baixa na referencia denuncia primer errado.
3. Escolhe `region_length` a partir dos dados, nao de chute, e diz quanto se
   perde com a escolha.
4. Grava o regions_multiregion.tsv pronto para o `--multiregion`.

Uso
---
    python3 validar_regioes_sidle.py \
        --primers primers_painel.tsv \
        --ref bancos/silva_nr99_v138.2_toSpecies_trainset.fa.gz \
        --asv metricas/comprimento_asv.tsv \
        --out regions_multiregion.tsv

    # opcoes uteis
        --n 20000          quantas sequencias de referencia amostrar
        --erros 1          erros tolerados no primer (0-2)
        --perda 5.0        % de ASV que se aceita perder ao cortar
        --colab renata     usar so os ASVs de um colaborador

Os ASVs observados devem estar SEM primer (que e' o caso: o cutadapt remove),
porque o script extrai a referencia sem primer tambem — e' o mesmo objeto.

So biblioteca padrao. Compativel com Python 3.6.
"""

import argparse
import gzip
import itertools
import os
import random
import re
import sys
from collections import defaultdict

# --------------------------------------------------------------- IUPAC

CLASSE = {
    "A": "A", "C": "C", "G": "G", "T": "T",
    "R": "AG", "Y": "CT", "S": "CG", "W": "AT", "K": "GT", "M": "AC",
    "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT",
}
COMPL = {"A": "T", "C": "G", "G": "C", "T": "A",
         "R": "Y", "Y": "R", "S": "S", "W": "W", "K": "M", "M": "K",
         "B": "V", "D": "H", "H": "D", "V": "B", "N": "N"}

# as 3 ultimas bases do primer nunca viram curinga: e' a extremidade 3', onde a
# polimerase exige pareamento de verdade. Tolerar erro ali seria modelar uma
# PCR que nao acontece.
FIXO_3P = 3


def revcomp(s):
    return "".join(COMPL.get(b, "N") for b in reversed(s))


def regex_primer(primer, erros, fixar="fim"):
    """
    Primer IUPAC -> regex, tolerando ate `erros` discordancias fora da
    extremidade 3'. Implementado por alternancia explicita de posicoes curinga:
    o `re` resolve isso em C, e a alternativa (varrer em Python) seria ordens de
    grandeza mais lenta num banco do tamanho do SILVA.

    `fixar` diz ONDE fica a extremidade 3' DENTRO DO PADRAO. Para o primer
    forward ela e' o fim do padrao; para o reverse, o padrao buscado e' o
    reverso-complemento, entao a extremidade 3' do primer e' o INICIO dele.
    Proteger o lado errado seria tolerar erro justamente onde a polimerase nao
    tolera, e exigir exatidao onde ela tolera — o contrario do que se quer.
    """
    base = [("[%s]" % CLASSE[b] if len(CLASSE[b]) > 1 else b) for b in primer]
    if fixar == "inicio":
        elegiveis = list(range(min(FIXO_3P, len(primer)), len(primer)))
    else:
        elegiveis = list(range(max(0, len(primer) - FIXO_3P)))
    variantes = set()
    for k in range(0, erros + 1):
        for combo in itertools.combinations(elegiveis, k):
            v = list(base)
            for i in combo:
                v[i] = "."
            variantes.add("".join(v))
    return re.compile("|".join(sorted(variantes)))


# --------------------------------------------------------------- IO

def abrir(caminho):
    with open(caminho, "rb") as fh:
        gz = fh.read(2) == b"\x1f\x8b"
    return gzip.open(caminho, "rt") if gz else open(caminho)


def amostrar_fasta(caminho, n, semente=1):
    """
    Amostragem por reservatorio. Pegar as n primeiras sequencias seria pior que
    inutil: o SILVA vem ordenado por taxonomia, entao as primeiras n sao um
    punhado de clados, nao o banco.
    """
    rng = random.Random(semente)
    reserva = []
    vistas = 0
    atual = []

    def guarda(seq):
        nonlocal vistas
        if not seq:
            return
        vistas += 1
        if len(reserva) < n:
            reserva.append(seq)
        else:
            j = rng.randrange(vistas)
            if j < n:
                reserva[j] = seq

    with abrir(caminho) as fh:
        for linha in fh:
            if linha.startswith(">"):
                guarda("".join(atual))
                atual = []
            else:
                atual.append(linha.strip().upper().replace("U", "T"))
        guarda("".join(atual))
    return reserva, vistas


def ler_primers(caminho):
    regioes = []
    with open(caminho) as fh:
        for linha in fh:
            if linha.startswith("#") or not linha.strip():
                continue
            campos = linha.rstrip("\n").split("\t")
            if campos[0] == "regiao":
                continue
            regioes.append((campos[0], campos[1].upper(), campos[2].upper()))
    return regioes


def ler_observado(caminho, colab=None, regioes_validas=None):
    """comprimento_asv.tsv -> {regiao: {comprimento: n_asv}}"""
    hist = defaultdict(lambda: defaultdict(int))
    with open(caminho) as fh:
        cab = fh.readline().rstrip("\n").split("\t")
        idx = {n: i for i, n in enumerate(cab)}
        for linha in fh:
            c = linha.rstrip("\n").split("\t")
            if not c or len(c) < len(cab):
                continue
            if colab and c[idx["colaborador"]] != colab:
                continue
            reg = c[idx["regiao"]]
            # so as regioes canonicas: as variantes de diagnostico (V1V2_t160_ruim)
            # nao podem entrar na distribuicao que decide o corte
            if regioes_validas and reg not in regioes_validas:
                continue
            hist[reg][int(c[idx["comprimento"]])] += int(c[idx["n_asv"]])
    return hist


# --------------------------------------------------------------- estatistica

def percentil(hist, p):
    """Percentil de um histograma {valor: peso}."""
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


def acima_de(hist, corte):
    """Fracao do peso com valor >= corte."""
    total = sum(hist.values())
    if not total:
        return 0.0
    return 100.0 * sum(w for v, w in hist.items() if v >= corte) / total


def resumo(hist):
    if not sum(hist.values()):
        return "—"
    return "p5=%d  med=%d  p95=%d" % (percentil(hist, 5), percentil(hist, 50),
                                      percentil(hist, 95))


# --------------------------------------------------------------- principal

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--primers", default="primers_painel.tsv")
    ap.add_argument("--ref", required=True,
                    help="fasta do banco de referencia (pode ser .gz)")
    ap.add_argument("--asv", default=None,
                    help="metricas/comprimento_asv.tsv (opcional, mas e' a "
                         "metade que valida)")
    ap.add_argument("--out", default="regions_multiregion.tsv")
    ap.add_argument("--n", type=int, default=20000)
    ap.add_argument("--erros", type=int, default=1, choices=[0, 1, 2])
    ap.add_argument("--perda", type=float, default=5.0,
                    help="%% de ASV que se aceita descartar no corte")
    ap.add_argument("--colab", default=None)
    ap.add_argument("--min-amplicon", type=int, default=80)
    ap.add_argument("--max-amplicon", type=int, default=1200)
    args = ap.parse_args()

    regioes = ler_primers(args.primers)
    nomes = set(r[0] for r in regioes)

    print("Amostrando o banco: %s" % args.ref)
    refs, total = amostrar_fasta(args.ref, args.n)
    print("  %d sequencias no banco, %d amostradas\n" % (total, len(refs)))

    obs = {}
    if args.asv and os.path.isfile(args.asv):
        obs = ler_observado(args.asv, args.colab, nomes)
    elif args.asv:
        print("AVISO: nao achei %s — sem a metade observada, o script so "
              "descreve a referencia e nao valida nada.\n" % args.asv)

    linhas_saida = []
    print("%-6s %6s   %-26s %-26s  %s" %
          ("regiao", "recup", "referencia (in silico)", "ASVs observados",
           "diagnostico"))
    print("-" * 104)

    for regiao, fw, rv in regioes:
        rv_busca = revcomp(rv)
        rx_f = regex_primer(fw, args.erros, "fim")
        rx_r = regex_primer(rv_busca, args.erros, "inicio")

        hist_ref = defaultdict(int)
        achados = 0
        for seq in refs:
            # Primers de regioes vizinhas do painel se sobrepoem no gene — o
            # 338R com o 341F, o 805F com o 806R (essa colisao ja nos apareceu
            # no roteamento das reads, com 11% delas). Entao pegar o primeiro
            # casamento forward e o primeiro reverse depois dele pode emendar
            # duas regioes num amplicon quimerico. Varremos TODOS os pares
            # validos e ficamos com o mais curto: o par espurio sempre atravessa
            # regiao inteira a mais, entao e' sempre o mais longo.
            melhor = None
            for mf in rx_f.finditer(seq):
                mr = rx_r.search(seq, mf.end())
                if not mr:
                    continue
                # entre os primers, SEM eles — que e' exatamente o que o
                # cutadapt deixa nos nossos ASVs
                comp = mr.start() - mf.end()
                if args.min_amplicon <= comp <= args.max_amplicon:
                    if melhor is None or comp < melhor:
                        melhor = comp
            if melhor is not None:
                hist_ref[melhor] += 1
                achados += 1

        recup = 100.0 * achados / len(refs) if refs else 0.0
        hist_obs = obs.get(regiao, {})

        # ---- diagnostico
        diag = []
        if regiao.startswith("ITS"):
            # esperado: o banco e' 16S, o ITS nao esta la. Nao e' defeito do
            # primer, e sinalizar como defeito mandaria caçar problema que nao
            # existe.
            diag.append("fora do escopo (banco 16S)")
        elif recup < 20:
            diag.append("RECUPERACAO BAIXA — primer suspeito")
        if hist_ref and hist_obs:
            d = percentil(hist_obs, 50) - percentil(hist_ref, 50)
            if abs(d) <= 5:
                diag.append("bordas batem")
            else:
                diag.append("DESLOCAMENTO de %+d nt" % d)
        elif not hist_obs:
            diag.append("sem ASV observado")

        # ---- escolha do region_length
        # O Sidle corta em region_length e DESCARTA o que for mais curto. Entao
        # o corte tem que caber tanto nos ASVs quanto na referencia: quem manda
        # e' o menor dos dois percentis.
        if hist_obs and hist_ref:
            alvo = min(percentil(hist_obs, args.perda),
                       percentil(hist_ref, args.perda))
        elif hist_ref:
            alvo = percentil(hist_ref, args.perda)
        elif hist_obs:
            alvo = percentil(hist_obs, args.perda)
        else:
            alvo = 0

        print("%-6s %5.1f%%   %-26s %-26s  %s"
              % (regiao, recup, resumo(hist_ref), resumo(hist_obs),
                 "; ".join(diag)))
        if alvo:
            print("%-6s %s region_length = %d  "
                  "(mantem %.1f%% dos ASVs, %.1f%% da referencia)"
                  % ("", " " * 5, alvo, acima_de(hist_obs, alvo),
                     acima_de(hist_ref, alvo)))
        linhas_saida.append((regiao, alvo, fw, rv, recup))

    # ---- arquivo do --multiregion
    # O Sidle e' 16S com regiao de referencia; ITS1 nao entra e continua trilha
    # separada.
    with open(args.out, "w") as fh:
        fh.write("region\tregion_length\tFW_primer\tRV_primer\n")
        n = 0
        for regiao, alvo, fw, rv, recup in linhas_saida:
            if regiao.startswith("ITS"):
                continue
            if not alvo:
                continue
            fh.write("%s\t%d\t%s\t%s\n" % (regiao, alvo, fw, rv))
            n += 1

    print("\nGravado %s com %d regioes (ITS excluido: o Sidle e' 16S)." %
          (args.out, n))
    print("\nComo ler isto:")
    print("  'bordas batem'  -> primer e recorte corretos; pode seguir.")
    print("  'DESLOCAMENTO'  -> a borda da referencia nao e' a dos ASVs. O sinal")
    print("                     diz o lado e o numero, quantas bases.")
    print("  'RECUPERACAO BAIXA' -> o primer nao encontra o sitio no banco;")
    print("                     rode de novo com --erros 2 antes de concluir.")


if __name__ == "__main__":
    main()

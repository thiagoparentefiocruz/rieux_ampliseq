#!/usr/bin/env python3
"""
validate_reports.py — o final_reports/ esta no formato que o aspp espera?

    bin/validate_reports.py <projeto>/final_reports

Por que isto existe: `final_reports/` e o contrato entre o wrapper e o pacote
R. Contrato que ninguem verifica nao e contrato, e erro de formato em TSV e
particularmente traicoeiro — um separador faltando desloca uma coluna inteira
e o R le numero como texto, ou pior, le o numero errado sem reclamar.

Verifica, arquivo por arquivo:

  - o cabecalho, nome por nome e na ordem;
  - o numero de campos de CADA linha (e aqui que mora o separador perdido);
  - o tipo de cada coluna: inteiro e inteiro, decimal e decimal;
  - o separador decimal, que tem de ser ponto — virgula decimal vinda de
    locale pt_BR e o jeito mais facil de entregar um arquivo que o R le como
    texto;
  - porcentagem dentro de 0-100;
  - 'yes'/'no' nas colunas de controle;
  - espaco sobrando nas pontas de um campo, que vira nivel de fator fantasma.

E entre arquivos:

  - as regioes sao as mesmas em todos;
  - toda amostra citada em abundance_per_sample e prevalence existe em
    reads_per_region;
  - o n_total do prevalence bate com o numero de amostras nao-controle
    daquela regiao.

Sai com status 1 se achar qualquer problema, para poder entrar em teste.
"""
import argparse
import os
import sys
from collections import defaultdict

INT, FLOAT, TXT, SIMNAO = "int", "float", "text", "yes/no"
TXT_NA = "text or NA"   # ausencia legitima, escrita como NA
INT_NA = "int or NA"

# Arquivos do nivel fino. Sao opcionais porque um final_reports/ gerado por
# uma versao anterior nao os tem — e isso nao o torna invalido.
OPCIONAIS = ("asv_table.tsv", "asv_taxonomy.tsv")

ESPERADO = {
    "reads_per_region.tsv": [
        ("project", TXT), ("region", TXT), ("sample", TXT), ("control", SIMNAO),
        ("input", INT), ("filtered", INT), ("denoisedF", INT),
        ("denoisedR", INT), ("merged", INT), ("nonchim", INT)],
    "classification.tsv": [
        ("project", TXT), ("region", TXT), ("rank", TXT),
        ("n_asv", INT), ("n_classified", INT), ("pct", FLOAT)],
    "abundance.tsv": [
        ("project", TXT), ("region", TXT), ("rank", TXT), ("taxon", TXT),
        ("reads_samples", INT), ("rel_pct", FLOAT), ("n_samples", INT),
        ("reads_controls", INT)],
    "asv_length.tsv": [
        ("project", TXT), ("region", TXT), ("length", INT), ("n_asv", INT)],
    "abundance_per_sample.tsv": [
        ("project", TXT), ("region", TXT), ("rank", TXT), ("taxon", TXT),
        ("sample", TXT), ("control", SIMNAO), ("reads", INT),
        ("pct_sample", FLOAT)],
    "asv_table.tsv": [
        ("project", TXT), ("region", TXT), ("asv_id", TXT), ("sample", TXT),
        ("control", SIMNAO), ("reads", INT)],
    "asv_taxonomy.tsv": [
        ("project", TXT), ("region", TXT), ("asv_id", TXT), ("length", INT_NA),
        ("Kingdom", TXT), ("Phylum", TXT), ("Class", TXT), ("Order", TXT),
        ("Family", TXT), ("Genus", TXT), ("Species", TXT)],
    "prevalence.tsv": [
        ("project", TXT), ("region", TXT), ("rank", TXT), ("taxon", TXT),
        ("n_present", INT), ("n_total", INT), ("pct_prevalence", FLOAT),
        ("pct_median", FLOAT), ("pct_max", FLOAT), ("sample_max", TXT_NA),
        ("n_controls_present", INT)],
}

# colunas cujo valor tem de caber em 0-100
PERCENTUAIS = {"pct", "rel_pct", "pct_sample", "pct_prevalence",
               "pct_median", "pct_max"}


class Erros(object):
    def __init__(self, limite=5):
        self.itens = []
        self.limite = limite
        self.n = 0

    def add(self, msg):
        self.n += 1
        if len(self.itens) < self.limite:
            self.itens.append(msg)

    def imprime(self, prefixo="      "):
        for m in self.itens:
            print(prefixo + m)
        if self.n > len(self.itens):
            print("%s... e mais %d" % (prefixo, self.n - len(self.itens)))


def confere_arquivo(caminho, colunas, err):
    """Devolve (linhas, dados) — dados como lista de dicts."""
    with open(caminho) as _fh:
        _cab = _fh.readline().rstrip("\n").rstrip("\r").split("\t")
    if _cab and _cab[-1] == "sequence" and "sequence" not in [c[0] for c in colunas]:
        # asv_taxonomy.tsv sem --no-sequences: coluna extra, legitima
        colunas = list(colunas) + [("sequence", TXT)]
    nomes = [c[0] for c in colunas]
    tipos = dict(colunas)
    dados = []
    with open(caminho, "rb") as fh:
        bruto = fh.read()
    if bruto.startswith(b"\xef\xbb\xbf"):
        err.add("arquivo comeca com BOM UTF-8 — o R le a primeira coluna com "
                "o nome sujo")
        bruto = bruto[3:]
    if b"\r\n" in bruto:
        err.add("fim de linha do Windows (CRLF) — o ultimo campo de cada "
                "linha vem com \\r colado")
    try:
        texto = bruto.decode("utf-8")
    except UnicodeDecodeError as e:
        err.add("nao e UTF-8 valido: %s" % e)
        texto = bruto.decode("utf-8", "replace")

    linhas = texto.rstrip("\n").split("\n")
    if not linhas or not linhas[0].strip():
        err.add("arquivo vazio")
        return 0, dados

    cab = linhas[0].rstrip("\r").split("\t")
    if cab != nomes:
        faltam = [n for n in nomes if n not in cab]
        sobram = [n for n in cab if n not in nomes]
        if faltam or sobram:
            err.add("cabecalho diferente do esperado. Faltam: %s. Sobram: %s"
                    % (faltam or "-", sobram or "-"))
        else:
            err.add("cabecalho na ordem errada: %s" % " ".join(cab))
        return len(linhas) - 1, dados

    for i, linha in enumerate(linhas[1:], start=2):
        campos = linha.rstrip("\r").split("\t")
        if len(campos) != len(nomes):
            err.add("linha %d tem %d campos, esperados %d — provavel "
                    "separador perdido: %.70s"
                    % (i, len(campos), len(nomes), linha))
            continue
        reg = {}
        for nome, valor in zip(nomes, campos):
            t = tipos[nome]
            if valor != valor.strip():
                err.add("linha %d, coluna '%s': espaco sobrando nas pontas "
                        "(%r)" % (i, nome, valor))
            v = valor.strip()
            if t == INT:
                if not (v.lstrip("-").isdigit()):
                    err.add("linha %d, coluna '%s': '%s' nao e inteiro"
                            % (i, nome, v))
                    continue
                reg[nome] = int(v)
            elif t == FLOAT:
                if "," in v:
                    err.add("linha %d, coluna '%s': '%s' usa virgula decimal "
                            "— o R le isso como texto" % (i, nome, v))
                    continue
                try:
                    reg[nome] = float(v)
                except ValueError:
                    err.add("linha %d, coluna '%s': '%s' nao e numero"
                            % (i, nome, v))
                    continue
                if nome in PERCENTUAIS and not (-0.01 <= reg[nome] <= 100.01):
                    err.add("linha %d, coluna '%s': %s fora de 0-100"
                            % (i, nome, v))
            elif t == INT_NA:
                if v != "NA" and not v.lstrip("-").isdigit():
                    err.add("linha %d, coluna '%s': '%s' nao e inteiro nem NA"
                            % (i, nome, v))
                    continue
                reg[nome] = None if v == "NA" else int(v)
            elif t == TXT_NA:
                if not v:
                    err.add("linha %d, coluna '%s': vazio. Ausencia se "
                            "escreve NA — campo vazio e indistinguivel de "
                            "separador perdido" % (i, nome))
                reg[nome] = v
            elif t == SIMNAO:
                if v not in ("yes", "no"):
                    err.add("linha %d, coluna '%s': '%s' nao e yes/no"
                            % (i, nome, v))
                reg[nome] = v
            else:
                if not v:
                    err.add("linha %d, coluna '%s': vazio" % (i, nome))
                reg[nome] = v
        dados.append(reg)
    return len(linhas) - 1, dados


def main():
    ap = argparse.ArgumentParser(
        description="Check that final_reports/ matches the contract with aspp.")
    ap.add_argument("reports", help="the final_reports/ directory")
    ap.add_argument("--max-errors", type=int, default=5,
                    help="how many examples to print per file (default 5)")
    args = ap.parse_args()

    if not os.path.isdir(args.reports):
        sys.exit("ERROR: %s is not a directory" % args.reports)

    print("Checking %s\n" % os.path.abspath(args.reports))
    print("  %-26s %9s %9s  %s" % ("file", "rows", "problems", "status"))
    print("  " + "-" * 60)

    tudo = {}
    total_erros = 0
    for arquivo in sorted(ESPERADO):
        caminho = os.path.join(args.reports, arquivo)
        err = Erros(args.max_errors)
        if not os.path.isfile(caminho):
            if arquivo in OPCIONAIS:
                print("  %-26s %9s %9s  absent (optional)"
                      % (arquivo, "-", "-"))
            else:
                print("  %-26s %9s %9s  MISSING" % (arquivo, "-", "-"))
                total_erros += 1
            continue
        n, dados = confere_arquivo(caminho, ESPERADO[arquivo], err)
        tudo[arquivo] = dados
        total_erros += err.n
        print("  %-26s %9d %9d  %s"
              % (arquivo, n, err.n, "ok" if err.n == 0 else "FAIL"))
        err.imprime()

    # ------------------------------------------------ coerencia entre arquivos
    print("\n  Cross-file checks")
    print("  " + "-" * 60)
    cruz = Erros(args.max_errors)
    notas = []

    # Coerencia de regioes, com a distincao que importa.
    #
    # As quatro tabelas taxonomicas tem de concordar ENTRE SI. Ja o
    # reads_per_region lista toda regiao que rodou, inclusive as que nao
    # produziram ASV nenhuma — e regiao sem ASV ausente das tabelas
    # taxonomicas nao e inconsistencia, e o unico resultado possivel.
    #
    # O que E problema: regiao que TEM ASV no asv_table e mesmo assim nao
    # aparece nas taxonomicas. Isso significa que a taxonomia nao foi
    # encontrada ou nao foi lida, e a regiao sai da analise em silencio.
    TAXONOMICAS = ("abundance.tsv", "abundance_per_sample.tsv",
                   "prevalence.tsv", "classification.tsv")
    regioes = {}
    for arquivo, dados in tudo.items():
        regioes[arquivo] = set(d.get("region") for d in dados) if dados else set()

    tax_regs = [regioes.get(a, set()) for a in TAXONOMICAS if a in tudo]
    if tax_regs:
        base_tax = set.union(*tax_regs)
        for arquivo in TAXONOMICAS:
            if arquivo not in tudo:
                continue
            faltando = base_tax - regioes[arquivo]
            if faltando:
                cruz.add("%s nao tem as regioes %s, que aparecem nas outras "
                         "tabelas taxonomicas" % (arquivo, " ".join(sorted(faltando))))

        com_asv = regioes.get("asv_table.tsv", set())
        rodaram = regioes.get("reads_per_region.tsv", set())
        RANKS_T = ("Kingdom", "Phylum", "Class", "Order",
                   "Family", "Genus", "Species")
        for r in sorted(rodaram - base_tax):
            if r not in com_asv:
                notas.append("regiao %s rodou e nao produziu ASV nenhuma "
                             "(ausente das tabelas taxonomicas, como esperado)"
                             % r)
                continue
            linhas_r = [d for d in tudo.get("asv_taxonomy.tsv", [])
                        if d.get("region") == r]
            com_tax = sum(1 for d in linhas_r
                          if any(d.get(k, "NA") != "NA" for k in RANKS_T))
            if com_tax:
                # A taxonomia EXISTE e mesmo assim a regiao nao entrou nas
                # tabelas agregadas: isso e falha de agregacao, nao do dado.
                cruz.add("regiao %s tem %d ASV COM taxonomia e nao aparece em "
                         "nenhuma tabela taxonomica — a agregacao perdeu a "
                         "regiao" % (r, com_tax))
            else:
                # Nenhuma ASV recebeu taxonomia. Acontece de verdade: no ITS,
                # o ITSx descarta ASV sem as ancoras conservadas, e se nao
                # sobrar nenhuma os passos seguintes nao rodam e o ampliseq
                # nem escreve a tabela de taxonomia.
                notas.append("regiao %s tem %d ASV e NENHUMA recebeu "
                             "taxonomia — fica fora das tabelas taxonomicas. "
                             "No ITS isso e esperado quando so ha controles: "
                             "o construto sintetico nao passa pelo ITSx"
                             % (r, len(linhas_r)))

    amostras = defaultdict(set)     # region -> samples
    nao_controle = defaultdict(set)
    for d in tudo.get("reads_per_region.tsv", []):
        amostras[d["region"]].add(d["sample"])
        if d.get("control") == "no":
            nao_controle[d["region"]].add(d["sample"])

    for arquivo, coluna in (("abundance_per_sample.tsv", "sample"),
                            ("prevalence.tsv", "sample_max")):
        for d in tudo.get(arquivo, []):
            s = d.get(coluna)
            if s == "NA":
                continue
            if s and s not in amostras.get(d["region"], set()):
                cruz.add("%s: '%s' (regiao %s) nao existe em "
                         "reads_per_region.tsv" % (arquivo, s, d["region"]))

    # NA em sample_max tem de significar exatamente "nao ocorre em amostra
    # nenhuma". Se aparecer com n_present > 0, o NA esta escondendo um erro.
    for d in tudo.get("prevalence.tsv", []):
        na = d.get("sample_max") == "NA"
        zero = d.get("n_present") == 0
        if na != zero:
            cruz.add("prevalence.tsv: taxon '%s' (%s) tem n_present=%s e "
                     "sample_max=%s — um contradiz o outro"
                     % (d.get("taxon"), d.get("region"),
                        d.get("n_present"), d.get("sample_max")))

    conhecidos = set((d["region"], d["asv_id"])
                     for d in tudo.get("asv_taxonomy.tsv", []))
    if conhecidos:
        for d in tudo.get("asv_table.tsv", []):
            if (d["region"], d["asv_id"]) not in conhecidos:
                cruz.add("asv_table.tsv: ASV '%s' (%s) nao aparece em "
                         "asv_taxonomy.tsv" % (d["asv_id"], d["region"]))
    for d in tudo.get("asv_table.tsv", []):
        if d["sample"] not in amostras.get(d["region"], set()):
            cruz.add("asv_table.tsv: '%s' (regiao %s) nao existe em "
                     "reads_per_region.tsv" % (d["sample"], d["region"]))

    vistos = set()
    for d in tudo.get("prevalence.tsv", []):
        chave = (d["region"], d.get("n_total"))
        if chave in vistos:
            continue
        vistos.add(chave)
        esperado = len(nao_controle.get(d["region"], set()))
        if esperado and d.get("n_total") != esperado:
            cruz.add("prevalence.tsv: regiao %s diz n_total=%s, mas "
                     "reads_per_region tem %d amostras nao-controle"
                     % (d["region"], d.get("n_total"), esperado))

    total_erros += cruz.n
    if cruz.n == 0:
        print("  regions, sample names and n_total are consistent   ok")
    else:
        cruz.imprime("      ")
    for n in notas:
        print("      note: %s" % n)

    print("")
    if total_erros == 0:
        print("  PASS — final_reports/ matches the contract.")
        print("  In R:  data <- aspp::read_ampliseq_summary('%s')"
              % os.path.abspath(args.reports))
        return 0
    print("  FAIL — %d problem(s). Do not feed this to aspp before fixing."
          % total_erros)
    return 1


if __name__ == "__main__":
    try:
        sys.exit(main())
    except BrokenPipeError:
        try:
            sys.stdout.close()
        finally:
            sys.exit(0)

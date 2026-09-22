#!/usr/bin/env python3
"""
prepare_unite.py — transforma o download do UNITE nos DOIS FASTA que o DADA2 usa

O UNITE e distribuido como TARBALL (.tgz) contendo varios FASTA, nao como um
FASTA comprimido. Quando o ampliseq baixa esse banco pelo caminho oficial
(--dada_ref_taxonomy unite-fungi) ele desempacota E REFORMATA sozinho, via
`bin/taxref_reformat_unite.sh`. Como pre-baixamos o arquivo e vamos passa-lo
com --dada_ref_tax_custom, NENHUM desses dois passos acontece.

O desempacotamento a versao anterior deste script ja resolvia. A reformatacao
nao — e era o buraco: o cabecalho cru do UNITE e

    >Gyroporus_purpurinus|KX389110|SH0879786.10FU|reps|k__Fungi;p__Basidio...

e o assignTaxonomy do DADA2 le o cabecalho como a propria string de taxonomia,
partindo no ';'. Sem reformatar, todo rank sai com o prefixo 'k__', 'p__'
colado, e ranks 'unidentified' entram como nome valido. Nao quebra: produz
taxonomia errada em silencio, que e pior.

Este script replica o taxref_reformat_unite.sh do ampliseq, byte a byte, e
produz os DOIS arquivos que ele produz:

    assignTaxonomy  -> DB_UNITE      (--dada_ref_tax_custom)
    addSpecies      -> DB_UNITE_SP   (--dada_ref_tax_custom_sp)

O segundo importa: o ampliseq NAO pula o addSpecies no ramo do UNITE. Ele
deriva um arquivo proprio, no formato ">ID Genero especie", do mesmo download.
Usar os dois reproduz o caminho oficial; pular o addSpecies seria uma escolha
metodologica diferente da do pipeline de referencia.

Por que Python e nao shell: os nos de computacao deste cluster nao tem `tar`
nem GNU sed garantidos no PATH. tarfile, gzip e re sao biblioteca padrao.

Uso:
    python3 prepare_unite.py [caminho/para/bancos.env]

E idempotente: reexecutar reformata a partir do arquivo bruto preservado em
DB_UNITE_BRUTO, nunca em cima de um arquivo ja reformatado.

Compativel com Python 3.6.
"""

import gzip
import os
import re
import shutil
import sys
import tarfile

# sem caminho de ninguem cravado: vem do RIEUX_PIPELINE_BASE, igual ao env.sh
PADRAO = os.path.join(os.environ.get("RIEUX_PIPELINE_BASE", "."), "bancos.env")


def ler_env(caminho):
    """Le as linhas `export VAR="valor"` de um arquivo de ambiente."""
    valores = {}
    with open(caminho) as fh:
        for linha in fh:
            m = re.match(r'\s*export\s+([A-Za-z_][A-Za-z0-9_]*)\s*=\s*(.*)', linha)
            if m:
                valores[m.group(1)] = m.group(2).strip().strip('"').strip("'")
    return valores


def comprimido(caminho):
    """
    Detecta gzip pelos BYTES MAGICOS, nao pela extensao.

    O UNITE entrega o tarball como .tgz, e ha downloads que chegam sem
    extensao nenhuma (o nosso veio como <uuid>.tgz). Decidir por sufixo faz
    o script ler bytes gzip como se fossem texto e concluir "formato
    desconhecido" — exatamente o que aconteceu aqui.
    """
    with open(caminho, "rb") as fh:
        return fh.read(2) == b"\x1f\x8b"


def abrir(caminho, modo="rt"):
    """Abre gzip ou texto puro de forma transparente."""
    abridor = gzip.open if comprimido(caminho) else open
    if "t" in modo:
        return abridor(caminho, modo, errors="replace")
    return abridor(caminho, modo)


def formato(caminho):
    """'fasta', 'tar' ou 'desconhecido', olhando o conteudo descomprimido."""
    abridor = gzip.open if comprimido(caminho) else open
    with abridor(caminho, "rb") as fh:
        inicio = fh.read(512)
    if inicio[:1] == b">":
        return "fasta", inicio
    # Bytes 257..261 sao "ustar" em qualquer tar moderno. Criterio muito mais
    # confiavel que chamar o `tar`, que pode nem existir na maquina.
    if inicio[257:262] == b"ustar":
        return "tar", inicio
    return "desconhecido", inicio


def seguro(membro, destino):
    """Impede que um membro do tar escreva fora do diretorio de destino."""
    alvo = os.path.realpath(os.path.join(destino, membro.name))
    return alvo.startswith(os.path.realpath(destino) + os.sep)


def escolher(fastas):
    """
    O UNITE publica variantes no mesmo tarball:
      sh_general_release_dynamic_<data>.fasta     <- padrao para DADA2
      sh_general_release_dynamic_s_<data>.fasta   <- inclui singletons
      developer/...                               <- campos extras
    Preferimos a general dynamic, sem singletons, fora de developer/.
    """
    def preferida(p):
        base = os.path.basename(p)
        return ("/developer/" not in p
                and "_s_" not in base
                and "general_release" in base)

    ideais = [p for p in fastas if preferida(p)]
    if ideais:
        return sorted(ideais)[0], None
    fora_dev = [p for p in fastas if "/developer/" not in p]
    candidatos = fora_dev or fastas
    maior = max(candidatos, key=os.path.getsize)
    return maior, "nenhum 'general_release' sem singletons; usando o maior FASTA"


# --------------------------------------------------------------------------
# Reformatacao — equivalente exato de bin/taxref_reformat_unite.sh
#
#   sed '/^>/s/;k__.*//'              -> corta uma taxonomia duplicada, se
#                                        houver. No release dynamic o
#                                        separador antes de k__ e '|', entao
#                                        nao casa; mantido por fidelidade.
#   sed '/^>/s/[a-z]__unidentified//g' -> esvazia ranks nao identificados
#   sed '/^>/s/[a-z]__//g'             -> remove os prefixos de rank
#   sed '/^>/s/ /_/g'                  -> espaco vira underscore
# --------------------------------------------------------------------------

def cab_assign(h):
    h = re.sub(r';k__.*$', '', h)
    h = re.sub(r'[a-z]__unidentified', '', h)
    h = re.sub(r'[a-z]__', '', h)
    return h.replace(' ', '_')


def cab_species(h):
    #   sed 's/>\\([^|]\\+\\)|\\([^|]\\+|[^|]\\+\\)|.*/>\\2 \\1/'  depois  s/_/ /g
    # Entrada : >Gyroporus_purpurinus|KX389110|SH0879786.10FU|reps|Fungi;...
    # Saida   : >KX389110|SH0879786.10FU Gyroporus purpurinus
    novo = re.sub(r'^>([^|]+)\|([^|]+\|[^|]+)\|.*$', r'>\2 \1', h)
    return novo.replace('_', ' ')


def reformatar(bruto, db_dir):
    """Gera os dois FASTA reformatados. Devolve (assign_gz, species_gz, n)."""
    assign = os.path.join(db_dir, "unite_assignTaxonomy.fasta.gz")
    especie = os.path.join(db_dir, "unite_addSpecies.fasta.gz")
    n = 0
    amostra = []
    with abrir(bruto, "rt") as ent, \
            gzip.open(assign, "wt") as sa, \
            gzip.open(especie, "wt") as se:
        for linha in ent:
            if linha.startswith(">"):
                h = linha.rstrip("\n")
                ha = cab_assign(h)
                hs = cab_species(ha)
                sa.write(ha + "\n")
                se.write(hs + "\n")
                n += 1
                if len(amostra) < 2:
                    amostra.append((h, ha, hs))
            else:
                sa.write(linha)
                se.write(linha)
    return assign, especie, n, amostra


def main():
    if any(a in ("-h", "--help") for a in sys.argv[1:]):
        # Sem argparse aqui de proposito (o script tem um argumento so), mas
        # uma ferramenta publica que nao responde a --help e uma porta fechada.
        print(__doc__.strip())
        return 0
    env_path = sys.argv[1] if len(sys.argv) > 1 else PADRAO
    if not os.path.isfile(env_path):
        sys.exit("ERROR: %s not found" % env_path)

    env = ler_env(env_path)
    base_dir = os.path.dirname(os.path.abspath(env_path))
    db_dir = os.path.join(base_dir, "db")

    # Fonte: o bruto preservado, se ja houver; senao o DB_UNITE atual. Isto e
    # o que torna a reexecucao segura — nunca reformatamos duas vezes.
    orig = env.get("DB_UNITE_BRUTO") or env.get("DB_UNITE_TARBALL") or env.get("DB_UNITE")
    if not orig or not os.path.isfile(orig):
        sys.exit("ERROR: DB_UNITE missing, or pointing at a file that does not exist")

    print("Fonte : %s" % orig)
    tipo, inicio = formato(orig)
    print("Formato: %s" % tipo)

    if tipo == "desconhecido":
        print("Primeiros bytes: %r" % inicio[:80], file=sys.stderr)
        sys.exit("ERROR: neither FASTA nor tar.")

    # ---- 1. desempacotar, se for tarball
    tarball = orig if tipo == "tar" else env.get("DB_UNITE_TARBALL", "")
    if tipo == "tar":
        destino = os.path.join(db_dir, "unite_extraido")
        os.makedirs(destino, exist_ok=True)
        print("\nDesempacotando em %s ..." % destino)
        with tarfile.open(orig, "r:gz") as tf:
            membros = [m for m in tf.getmembers() if seguro(m, destino)]
            tf.extractall(destino, members=membros)

        fastas = []
        for raiz, _dirs, arquivos in os.walk(destino):
            for a in arquivos:
                if a.endswith((".fasta", ".fa", ".fna")):
                    fastas.append(os.path.join(raiz, a))
        if not fastas:
            sys.exit("ERROR: no FASTA inside the tarball.")

        print("\nFASTA encontrados:")
        for f in sorted(fastas):
            print("  %-68s %6.1f MB"
                  % (os.path.relpath(f, destino), os.path.getsize(f) / 1e6))

        escolhido, aviso = escolher(fastas)
        if aviso:
            print("\nAVISO: %s" % aviso)
        print("\nEscolhido: %s" % os.path.basename(escolhido))

        bruto = os.path.join(db_dir,
                             re.sub(r"\.(fasta|fa|fna)$", "", os.path.basename(escolhido))
                             + ".fasta.gz")
        with open(escolhido, "rb") as ent, gzip.open(bruto, "wb") as sai:
            shutil.copyfileobj(ent, sai)
        print("Raw file written: %s (%.1f MB)" % (bruto, os.path.getsize(bruto) / 1e6))
    else:
        bruto = orig

    # ---- 2. conferir que o bruto ainda esta cru
    with abrir(bruto, "rt") as fh:
        primeiro = fh.readline().rstrip("\n")
    if "k__" not in primeiro:
        print("\nCabecalho: %s" % primeiro[:120], file=sys.stderr)
        sys.exit("ERROR: this file has no rank prefixes ('k__').\n"
                 "      Parece ja reformatado. Aponte DB_UNITE_BRUTO para o\n"
                 "      download original antes de rodar de novo.")

    # ---- 3. reformatar
    print("\nReformatando (equivalente a bin/taxref_reformat_unite.sh) ...")
    assign, especie, n, amostra = reformatar(bruto, db_dir)
    print("  %d sequencias" % n)
    for cru, ha, hs in amostra:
        print("\n  cru     : %s" % cru[:150])
        print("  assign  : %s" % ha[:150])
        print("  species : %s" % hs[:150])
    print("\n  %-52s %6.1f MB" % (os.path.basename(assign),
                                  os.path.getsize(assign) / 1e6))
    print("  %-52s %6.1f MB" % (os.path.basename(especie),
                                os.path.getsize(especie) / 1e6))

    # ---- 4. atualizar bancos.env, preservando o bruto
    shutil.copy(env_path, env_path + ".bak")
    with open(env_path) as fh:
        linhas = [l for l in fh
                  if not re.match(r'\s*export\s+DB_UNITE(_BRUTO|_SP|_TARBALL)?\s*=', l)]
    if tarball:
        linhas.append('export DB_UNITE_TARBALL="%s"\n' % tarball)
    linhas.append('export DB_UNITE_BRUTO="%s"\n' % bruto)
    linhas.append('export DB_UNITE="%s"\n' % assign)
    linhas.append('export DB_UNITE_SP="%s"\n' % especie)
    with open(env_path, "w") as fh:
        fh.writelines(linhas)

    print("\nbancos.env updated (backup at %s.bak):" % env_path)
    for chave, valor in sorted(ler_env(env_path).items()):
        if chave.startswith("DB_"):
            print("  %-20s %s" % (chave, valor))
    print("\nRecarregue o ambiente antes de rodar o ITS1:")
    print("  source <repo>/bin/env.sh")


if __name__ == "__main__":
    main()

#!/usr/bin/env bash
#
# rieux_ampliseq.sh — nf-core/ampliseq multi-regiao, por estagios
#
# Leva um painel de amplicons multi-regiao do FASTQ bruto ate as tabelas
# consolidadas, em seis estagios. Qualquer um pode ser rodado sozinho, e a
# execucao pode comecar ou parar em qualquer ponto.
#
#     organizar -> descobrir -> dividir -> perfilar -> rodar -> sidle -> consolidar
#
#   organizar   cruza a planilha de amostras com os FASTQ em disco e monta um
#               projeto por grupo de amostras. Opcional: quem ja tem os FASTQ
#               separados comeca em `descobrir`.
#   descobrir   recupera os primers do painel a partir das proprias reads.
#               Opcional: quem conhece os primers passa --primers.
#   dividir     roteia cada par de reads para a sua regiao, exigindo o PAR de
#               primers, e escreve uma samplesheet por regiao.
#   perfilar    mede a qualidade real por ciclo e escolhe truncLenF/R por
#               regiao. Vem DEPOIS de dividir porque o corte e por regiao.
#   rodar       executa o ampliseq uma vez por regiao.
#   sidle       executa o ramo multi-regiao (reconstrucao contra referencia).
#   consolidar  escreve final_reports/, o contrato com o pacote R `aspp`.
#
# Do zero, um comando:
#
#     rieux_ampliseq.sh --projeto renata --brutos brutos/renata
#
# Retomando no meio, quando o que veio antes ja existe:
#
#     rieux_ampliseq.sh --projeto renata --from rodar
#     rieux_ampliseq.sh --projeto renata --stage consolidar
#
# ---------------------------------------------------------------------------
# O DIRETORIO DO PROJETO
#
# Os estagios nao se comunicam por flags: eles se comunicam por um diretorio
# com layout conhecido. E isso que permite retomar em qualquer ponto sem
# reinformar o que veio antes.
#
#     <projeto>/
#       brutos/                    FASTQ (links), de `organizar`
#       metadata.tsv               amostra -> grupo, de `organizar`
#       primers.tsv                de `descobrir`, ou copiado de --primers
#       parametros_regioes.tsv     truncLen por regiao, de `perfilar`
#       split/                     reads roteadas, de `dividir`
#         samplesheets/            uma por regiao + a completa
#       <REGIAO>/                  uma execucao do ampliseq, de `rodar`
#       sidle/                     de `sidle`
#       final_reports/             de `consolidar`
#       logs/
#
# Qualquer um desses caminhos pode ser sobrescrito por flag. O layout e o
# padrao, nao uma imposicao.
#
# ---------------------------------------------------------------------------
# ONDE RODAR
#
# O driver do Nextflow roda no NO DE LOGIN, dentro de um `screen`. Ele fica
# ocioso esperando o SLURM e submete as tarefas; nao precisa de alocacao.
#
# As regioes rodam EM SEQUENCIA, de proposito: duas regioes do mesmo projeto so
# disputariam a mesma fila. Dois PROJETOS em paralelo e outra coisa — ai sao
# particoes diferentes, e ai ganha-se de verdade:
#
#     rieux_ampliseq.sh --projeto fabio    --from rodar                    # cpu
#     rieux_ampliseq.sh --projeto patricia --from rodar \
#                       --partition fat --work-dir exec/patricia
#
# O --work-dir nao e opcional nesse caso: o Nextflow guarda `.nextflow.log`,
# `.nextflow/history` e o cache do -resume no diretorio de lancamento, e dois
# drivers no mesmo diretorio destroem o historico um do outro. Nao corrompe
# resultado — o hash de cada tarefa continua unico no work/ — mas destroi a
# capacidade de depurar quando algo falha.
# ---------------------------------------------------------------------------

set -uo pipefail

AQUI="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Ctrl-C tem que parar TUDO, nao so o estagio da vez.
#
# Sem isto o SIGINT mata apenas o processo em primeiro plano; o laco do bash
# sobrevive e alegremente comeca o proximo. Ja aconteceu aqui: uma execucao
# interrompida "continuou sozinha".
trap 'echo; echo "interrompido pelo usuario."; exit 130' INT

ESTAGIOS=(organizar descobrir dividir perfilar rodar sidle consolidar)

# ------------------------------------------------------------------ padroes
PROJETO=""
RAIZ=""
BRUTOS=""
PLANILHA=""
PRIMERS=""
PARAMS=""
PARTICAO="cpu"
CONFIG="$AQUI/conf/rieux.config"
PIPE="${AMPLISEQ_HOME:-nf-core/ampliseq}"
DIR_EXEC="$PWD"
CONTROLES='^[Ss]mart'
MIN_AMOSTRAS=3
MIN_READS_REGIAO=1000
MULTIREGION=""
SIDLE_ENTRADA=""
SIDLE_REF="silva"
SIDLE_EXTRA=""
LOTE_SPLIT=10
REGIOES=()
DE=""; ATE=""; SO=""; PULAR=""
SIMULAR=0

ajuda() {
    sed -n '3,75p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'
    cat <<'FIM'

Opcoes
------
  --projeto NOME       rotulo do projeto: entra nos caminhos e nos logs.
                       Tambem e o diretorio de trabalho, salvo --outdir.
  --outdir DIR         raiz do projeto (padrao: ./<projeto>)

  Escolha de estagios (sem nenhuma delas, roda o fluxo todo):
  --stage NOME         roda SO este estagio
  --from NOME          deste estagio ate o fim
  --until NOME         do inicio ate este estagio
  --skip "A B"         pula estes

  Entradas (cada uma tem padrao dentro do diretorio do projeto):
  --brutos DIR         FASTQ brutos           (padrao <projeto>/brutos)
  --planilha ARQ       planilha de amostras   (so o estagio `organizar`)
  --primers ARQ        tabela regiao/forward/reverse. Passe se ja conhece os
                       primers do seu painel — pula o estagio `descobrir`.
  --params ARQ         tabela regiao/trunclenf/trunclenr/banco/extra. Passe se
                       ja tem truncLen escolhido — pula o estagio `perfilar`.
  --controles REGEX    nomes de amostra que sao controle (padrao: ^[Ss]mart)

  Execucao:
  --partition NOME     particao do SLURM (padrao: cpu)
  --regions "A B"      so estas regioes
  --work-dir DIR       diretorio de lancamento do Nextflow (padrao: o atual)
  --config ARQ         perfil do cluster
  --pipeline X         nome ou caminho do ampliseq (padrao: nf-core/ampliseq)
  --min-samples N      pula regiao com menos amostras que isto (padrao 3)
  --min-reads N        piso de reads por amostra POR REGIAO no split (padrao 1000)
  --lote N             amostras por tarefa do array de split (padrao 10)

  Ramo Sidle:
  --multiregion ARQ    regions_multiregion.tsv do validar_regioes_sidle.py
  --sidle-input ARQ    samplesheet das reads NAO divididas
  --sidle-ref NOME     banco do Sidle (padrao: silva)
  --sidle-extra "..."  argumentos extras so para este ramo

  --dry-run            so mostra os comandos
  -h, --help           esta ajuda

Nota sobre o estagio `organizar`
--------------------------------
Ele e o unico UM-PARA-MUITOS do fluxo: uma planilha com varios grupos de
amostras vira varios projetos. Por isso ele nao encadeia — roda, escreve os
projetos e imprime os comandos seguintes, um por projeto. Os outros estagios
sao todos um-para-um.

Nota sobre o ramo Sidle
-----------------------
Ele recebe as reads INTEIRAS, nao as divididas: faz o proprio roteamento,
rodando cutadapt uma vez por regiao com os primers do regions_multiregion.tsv.
Alimenta-lo com as reads ja divididas cortaria primer duas vezes e deslocaria
as bordas — que e exatamente o que o validar_regioes_sidle.py existe para
impedir.
FIM
}

# ------------------------------------------------------------------ opcoes
while [[ $# -gt 0 ]]; do
    case "$1" in
        --projeto|--name) PROJETO="${2:?}"; shift 2 ;;
        --projeto=*)    PROJETO="${1#*=}"; shift ;;
        --outdir)       RAIZ="${2:?}"; shift 2 ;;
        --outdir=*)     RAIZ="${1#*=}"; shift ;;
        --stage)        SO="${2:?}"; shift 2 ;;
        --from)         DE="${2:?}"; shift 2 ;;
        --until)        ATE="${2:?}"; shift 2 ;;
        --skip)         PULAR="${2:?}"; shift 2 ;;
        --brutos)       BRUTOS="${2:?}"; shift 2 ;;
        --planilha)     PLANILHA="${2:?}"; shift 2 ;;
        --primers)      PRIMERS="${2:?}"; shift 2 ;;
        --params)       PARAMS="${2:?}"; shift 2 ;;
        --controles)    CONTROLES="${2:?}"; shift 2 ;;
        --partition)    PARTICAO="${2:?}"; shift 2 ;;
        --regions)      read -r -a REGIOES <<< "${2:?}"; shift 2 ;;
        --work-dir)     DIR_EXEC="${2:?}"; shift 2 ;;
        --config)       CONFIG="${2:?}"; shift 2 ;;
        --pipeline)     PIPE="${2:?}"; shift 2 ;;
        --min-samples)  MIN_AMOSTRAS="${2:?}"; shift 2 ;;
        --min-reads)    MIN_READS_REGIAO="${2:?}"; shift 2 ;;
        --lote)         LOTE_SPLIT="${2:?}"; shift 2 ;;
        --multiregion)  MULTIREGION="${2:?}"; shift 2 ;;
        --sidle-input)  SIDLE_ENTRADA="${2:?}"; shift 2 ;;
        --sidle-ref)    SIDLE_REF="${2:?}"; shift 2 ;;
        --sidle-extra)  SIDLE_EXTRA="${2:?}"; shift 2 ;;
        --dry-run|--simular) SIMULAR=1; shift ;;
        -h|--help)      ajuda; exit 0 ;;
        *)              echo "ERRO: opcao desconhecida: $1" >&2; exit 1 ;;
    esac
done

[[ -n "$PROJETO" ]] || { echo "ERRO: --projeto e obrigatorio (use --help)" >&2; exit 1; }

# ------------------------------------------------- layout do projeto
[[ -n "$RAIZ" ]] || RAIZ="$PWD/$PROJETO"
[[ "$RAIZ" == /* ]] || RAIZ="$PWD/$RAIZ"

[[ -n "$BRUTOS" ]]  || BRUTOS="$RAIZ/brutos"
[[ -n "$PRIMERS" ]] || PRIMERS="$RAIZ/primers.tsv"
[[ -n "$PARAMS" ]]  || PARAMS="$RAIZ/parametros_regioes.tsv"
SPLIT="$RAIZ/split"
SS_DIR="$SPLIT/samplesheets"
LOGS="$RAIZ/logs"
[[ -n "$SIDLE_ENTRADA" ]] || SIDLE_ENTRADA="$SS_DIR/samplesheet_completo.tsv"

for v in BRUTOS PRIMERS PARAMS PLANILHA CONFIG DIR_EXEC MULTIREGION; do
    val="${!v}"
    [[ -z "$val" || "$val" == /* ]] || printf -v "$v" '%s' "$PWD/$val"
done

# ------------------------------------------------- quais estagios rodar
# O resultado sai por variavel, nao por stdout, de proposito. Com
# `i=$(indice X)` a funcao roda num SUBSHELL, e o `exit 1` dela encerraria
# apenas o subshell: um nome de estagio invalido passava batido e o script
# seguia com INI vazio. Bug real, pego no teste.
INDICE_RES=""
indice() {
    local alvo="$1" i
    for i in "${!ESTAGIOS[@]}"; do
        if [[ "${ESTAGIOS[$i]}" == "$alvo" ]]; then INDICE_RES="$i"; return 0; fi
    done
    echo "ERRO: estagio desconhecido: $alvo" >&2
    echo "      validos: ${ESTAGIOS[*]}" >&2
    return 1
}

if [[ -n "$SO" ]]; then
    indice "$SO" || exit 1; INI="$INDICE_RES"; FIM="$INDICE_RES"
else
    # Sem --from, o padrao comeca em `descobrir`, nao em `organizar`: organizar
    # e' um-para-muitos e exige uma planilha, entao entrar nele por acidente
    # seria surpresa. Quem quer a planilha pede por ela.
    indice "${DE:-descobrir}"   || exit 1; INI="$INDICE_RES"
    indice "${ATE:-consolidar}" || exit 1; FIM="$INDICE_RES"
fi
(( INI <= FIM )) || { echo "ERRO: --from vem depois de --until" >&2; exit 1; }

declare -A PULA=()
for p in $PULAR; do indice "$p" || exit 1; PULA["$p"]=1; done

rodar_estagio() {
    local nome="$1" i
    indice "$nome" || exit 1
    i="$INDICE_RES"
    (( i >= INI && i <= FIM )) || return 1
    [[ -n "${PULA[$nome]:-}" ]] && return 1
    return 0
}

exec_cmd() {
    if (( SIMULAR )); then
        printf '    '; printf '%q ' "$@"; echo
        return 0
    fi
    "$@"
}

precisa() {   # precisa <caminho> <descricao> <estagio que produz>
    [[ -e "$1" ]] && return 0
    # Em --dry-run o estagio anterior nao criou nada, entao exigir o arquivo
    # faria a simulacao do fluxo completo morrer sempre no segundo estagio —
    # justamente quando ela e mais util. Avisa e segue.
    if (( SIMULAR )); then
        echo "    (--dry-run: $2 nao existe ainda; viria de '$3')"
        return 0
    fi
    echo "ERRO: falta $2" >&2
    echo "      esperado em: $1" >&2
    echo "      produzido pelo estagio '$3' — rode-o antes, ou informe o" >&2
    echo "      caminho pela flag correspondente." >&2
    exit 1
}

mkdir -p "$RAIZ" "$LOGS"

echo "Projeto    : $PROJETO"
echo "Diretorio  : $RAIZ"
echo "Estagios   : ${ESTAGIOS[*]:$INI:$((FIM-INI+1))}${PULAR:+  (pulando: $PULAR)}"
echo "Particao   : $PARTICAO"
echo

# ------------------------------------------------- ambiente
if [[ -z "${DB_SILVA_GENERO:-}" && -r "$AQUI/bin/ambiente.sh" ]]; then
    # shellcheck disable=SC1091
    source "$AQUI/bin/ambiente.sh" || exit 1
fi

# =================================================================== organizar
if rodar_estagio organizar; then
    echo "== organizar ==================================================="
    [[ -n "$PLANILHA" ]] || { echo "ERRO: 'organizar' precisa de --planilha" >&2; exit 1; }
    precisa "$PLANILHA" "a planilha de amostras" organizar
    precisa "$BRUTOS" "o diretorio de FASTQ brutos (--brutos)" organizar
    exec_cmd python3 "$AQUI/bin/organizar_projeto.py" "$PLANILHA" "$BRUTOS" \
             "$RAIZ/projetos" --controles "$CONTROLES" --executar || exit 1
    echo
    echo "  'organizar' e um-para-muitos: ele escreveu um projeto por grupo em"
    echo "  $RAIZ/projetos/. Continue um por vez:"
    if (( ! SIMULAR )); then
        for d in "$RAIZ"/projetos/*/; do
            [[ -d "$d" ]] || continue
            n=$(basename "$d")
            echo "    rieux_ampliseq.sh --projeto $n --outdir $d --brutos $d/dados_brutos"
        done
    fi
    echo
    exit 0
fi

# =================================================================== descobrir
if rodar_estagio descobrir; then
    echo "== descobrir ==================================================="
    if [[ -s "$PRIMERS" ]]; then
        echo "  $PRIMERS ja existe — nada a fazer."
    else
        precisa "$BRUTOS" "o diretorio de FASTQ brutos (--brutos)" organizar
        exec_cmd python3 "$AQUI/bin/discover_primers.py" --dir "$BRUTOS" \
                 --out "$RAIZ/primers" || exit 1
        echo
        echo "  CONFIRA $RAIZ/primers.tsv antes de seguir. Primer recuperado"
        echo "  errado nao da erro: da resultado errado, em silencio."
    fi
    echo
fi

# ===================================================================== dividir
if rodar_estagio dividir; then
    echo "== dividir ====================================================="
    precisa "$PRIMERS" "a tabela de primers" descobrir
    precisa "$BRUTOS" "o diretorio de FASTQ brutos (--brutos)" organizar
    n_amostras=$(find -L "$BRUTOS" -name '*_R1*.fastq.gz' 2>/dev/null | wc -l)
    n_tarefas=$(( (n_amostras + LOTE_SPLIT - 1) / LOTE_SPLIT ))
    (( n_tarefas > 0 )) || { echo "ERRO: nenhum FASTQ R1 em $BRUTOS" >&2; exit 1; }
    echo "  $n_amostras amostras, lote $LOTE_SPLIT -> array 1-$n_tarefas"
    # --wait bloqueia ate o array terminar. O driver esta num screen no login,
    # entao bloquear e exatamente o comportamento desejado: os estagios
    # seguintes dependem deste.
    exec_cmd sbatch --wait --partition "$PARTICAO" \
             --array="1-${n_tarefas}%10" \
             --export="ALL,LOTE=$LOTE_SPLIT" \
             --output="$LOGS/split_%A_%a.log" \
             "$AQUI/bin/split_regioes.sh" "$BRUTOS" "$PRIMERS" "$SPLIT" || exit 1
    exec_cmd python3 "$AQUI/bin/resumo_split.py" "$SPLIT" \
             --minimo "$MIN_READS_REGIAO" || exit 1
    echo
fi

# ==================================================================== perfilar
if rodar_estagio perfilar; then
    echo "== perfilar ===================================================="
    if [[ -s "$PARAMS" ]]; then
        echo "  $PARAMS ja existe — nada a fazer."
    else
        precisa "$SPLIT/split" "as reads divididas por regiao" dividir
        exec_cmd python3 "$AQUI/bin/perfil_qualidade.py" "$SPLIT/split" \
                 --out "$PARAMS" || exit 1
    fi
    echo
fi

# ======================================================================= rodar
if rodar_estagio rodar; then
    echo "== rodar ======================================================="
    precisa "$PRIMERS" "a tabela de primers" descobrir
    precisa "$PARAMS" "a tabela de truncLen" perfilar
    precisa "$SS_DIR" "as samplesheets por regiao" dividir

    if command -v sinfo >/dev/null 2>&1; then
        sinfo -h -p "$PARTICAO" -o '%P' 2>/dev/null | grep -q . || {
            echo "ERRO: a particao '$PARTICAO' nao existe neste cluster" >&2
            sinfo -h -o '        %P  %D nos  %m MB' >&2; exit 1; }
    fi

    export RIEUX_QUEUE="$PARTICAO"
    case "$PARTICAO" in
        fat) export RIEUX_QUEUESIZE="${RIEUX_QUEUESIZE:-40}"
             export RIEUX_MEM_HIGH="${RIEUX_MEM_HIGH:-256}" ;;
        *)   export RIEUX_QUEUESIZE="${RIEUX_QUEUESIZE:-60}"
             export RIEUX_MEM_HIGH="${RIEUX_MEM_HIGH:-160}" ;;
    esac
    export RIEUX_QUEUE_TAX="${RIEUX_QUEUE_TAX:-fat}"

    mkdir -p "$DIR_EXEC"
    cd "$DIR_EXEC" || exit 1

    if [[ ${#REGIOES[@]} -eq 0 ]]; then
        mapfile -t REGIOES < <(awk -F'\t' '!/^#/ && !/^regiao\t/ && NF>=4 {print $1}' "$PARAMS")
    fi

    declare -a OK=() FALHOU=() PULADO=()
    for REG in "${REGIOES[@]}"; do
        SS="$SS_DIR/samplesheet_${REG}.tsv"
        [[ -s "$SS" ]] || { echo "[$REG] sem samplesheet — pulando"; PULADO+=("$REG"); continue; }
        N=$(( $(wc -l < "$SS") - 1 ))
        (( N >= MIN_AMOSTRAS )) || { echo "[$REG] so $N amostra(s) — pulando"; PULADO+=("$REG:$N"); continue; }

        read -r FW RV < <(awk -F'\t' -v r="$REG" '$1==r {print $2" "$3}' "$PRIMERS")
        [[ -n "${FW:-}" && -n "${RV:-}" ]] || { echo "[$REG] primers ausentes — pulando"; PULADO+=("$REG:primers"); continue; }

        read -r TF TR BANCO EXTRA < <(awk -F'\t' -v r="$REG" \
            '$1==r {printf "%s %s %s %s", $2, $3, $4, ($5==""?"-":$5)}' "$PARAMS")
        [[ "$EXTRA" == "-" ]] && EXTRA=""

        ARGS=(run "$PIPE" -profile singularity -c "$CONFIG" -resume
              --input "$SS" --outdir "$RAIZ/$REG"
              --FW_primer "$FW" --RV_primer "$RV"
              --illumina_novaseq --dada_min_boot 80)
        [[ "$TF" != "0" && "$TR" != "0" ]] && ARGS+=(--trunclenf "$TF" --trunclenr "$TR")

        if [[ "$BANCO" == "unite" ]]; then
            # DOIS arquivos, como no SILVA — e nao --skip_dada_addspecies.
            #
            # A armadilha central: --dada_ref_tax_custom NAO roda o fmtscript
            # que o --dada_ref_taxonomy roda. Quem passa o banco na mao passa o
            # banco JA formatado, ou o passa cru sem receber aviso nenhum. O
            # bin/preparar_unite.py existe para reproduzir aquele fmtscript, e
            # ele deriva DOIS FASTA: um para o assignTaxonomy e outro, no
            # formato ">ID Genero especie", para o addSpecies.
            if [[ -z "${DB_UNITE_SP:-}" || ! -r "${DB_UNITE_SP:-/dev/null}" ]]; then
                echo "[$REG] DB_UNITE_SP ausente — rode '$AQUI/bin/preparar_unite.py'" >&2
                PULADO+=("$REG:banco"); continue
            fi
            ARGS+=(--dada_ref_tax_custom "$DB_UNITE"
                   --dada_ref_tax_custom_sp "$DB_UNITE_SP"
                   --dada_ref_taxlevels "Kingdom,Phylum,Class,Order,Family,Genus,Species")
        else
            ARGS+=(--dada_ref_tax_custom "$DB_SILVA_GENERO"
                   --dada_ref_tax_custom_sp "$DB_SILVA_ESPECIE")
        fi
        # shellcheck disable=SC2206
        [[ -n "$EXTRA" ]] && ARGS+=($EXTRA)

        LOG="$LOGS/${PROJETO}_${REG}.log"
        echo "  [$REG] $N amostras | truncLen ${TF}/${TR} | banco $BANCO"
        if (( SIMULAR )); then
            printf '    nextflow'; printf ' %q' "${ARGS[@]}"; echo
            OK+=("$REG"); continue
        fi
        INICIO=$SECONDS
        if nextflow "${ARGS[@]}" 2>&1 | tee -a "$LOG"; then
            echo "    concluido em $(( (SECONDS-INICIO)/60 )) min"; OK+=("$REG")
        else
            echo "    FALHOU — veja $LOG"; FALHOU+=("$REG")
        fi
    done
    echo
    if (( SIMULAR )); then
        [[ ${#OK[@]} -gt 0 ]] && echo "  simuladas : ${OK[*]}"
    else
        [[ ${#OK[@]} -gt 0 ]] && echo "  concluidas: ${OK[*]}"
    fi
    [[ ${#PULADO[@]} -gt 0 ]] && echo "  puladas   : ${PULADO[*]}"
    [[ ${#FALHOU[@]} -gt 0 ]] && echo "  falharam  : ${FALHOU[*]}"
    echo
fi

# ======================================================================= sidle
if rodar_estagio sidle; then
    if [[ -z "$MULTIREGION" ]]; then
        if [[ -z "$SO" ]]; then
            echo "== sidle: pulado (sem --multiregion) ==========================="
            echo "  Gere o arquivo com:"
            echo "    $AQUI/bin/validar_regioes_sidle.py --primers $PRIMERS \\"
            echo "        --ref \"\$DB_SILVA_GENERO\" \\"
            echo "        --asv $RAIZ/final_reports/comprimento_asv.tsv \\"
            echo "        --out $RAIZ/regions_multiregion.tsv"
            echo
        else
            echo "ERRO: o estagio 'sidle' precisa de --multiregion" >&2; exit 1
        fi
    else
        echo "== sidle ======================================================="
        precisa "$MULTIREGION" "o regions_multiregion.tsv" sidle
        precisa "$SIDLE_ENTRADA" "a samplesheet das reads NAO divididas" dividir
        mkdir -p "$DIR_EXEC"; cd "$DIR_EXEC" || exit 1
        export RIEUX_QUEUE="$PARTICAO"
        # Nao passamos --trunclenf/--trunclenr aqui de proposito: eles sao
        # globais, e no multi-regiao cada regiao tem o seu. Quem governa o corte
        # e o region_length do proprio regions_multiregion.tsv, escolhido a
        # partir da distribuicao de comprimento observada.
        ARGS=(run "$PIPE" -profile singularity -c "$CONFIG" -resume
              --input "$SIDLE_ENTRADA" --multiregion "$MULTIREGION"
              --outdir "$RAIZ/sidle" --illumina_novaseq
              --sidle_ref_taxonomy "$SIDLE_REF")
        # shellcheck disable=SC2206
        [[ -n "$SIDLE_EXTRA" ]] && ARGS+=($SIDLE_EXTRA)
        LOG="$LOGS/${PROJETO}_sidle.log"
        echo "  $(( $(wc -l < "$MULTIREGION") - 1 )) regioes | banco $SIDLE_REF"
        if (( SIMULAR )); then
            printf '    nextflow'; printf ' %q' "${ARGS[@]}"; echo
        else
            nextflow "${ARGS[@]}" 2>&1 | tee -a "$LOG" || echo "  FALHOU — veja $LOG"
        fi
        echo
    fi
fi

# ================================================================== consolidar
if rodar_estagio consolidar; then
    echo "== consolidar =================================================="
    exec_cmd python3 "$AQUI/bin/coletar_metricas.py" \
             --resultados "$RAIZ" --nome "$PROJETO" \
             --out "$RAIZ/final_reports" \
             --primers "$PRIMERS" --controles "$CONTROLES" || exit 1
    echo
    echo "  no R:  dados <- aspp::read_ampliseq_summary('$RAIZ/final_reports')"
fi

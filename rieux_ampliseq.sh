#!/usr/bin/env bash
#
# rieux_ampliseq.sh — nf-core/ampliseq multi-regiao no cluster rieux
#
# Executa o ampliseq em cada regiao de um painel de amplicons multi-regiao
# (QIAseq 16S/ITS Pro Screening Panel), opcionalmente executa tambem o ramo
# multi-regiao com Sidle, e consolida tudo em final_reports/ — que e o que o
# pacote R `aspp` le.
#
#     rieux_ampliseq.sh --name fabio --input split/samplesheets/fabio
#
# O driver do Nextflow roda no NO DE LOGIN, dentro de um `screen`. Ele fica
# ocioso esperando o SLURM e submete as tarefas; nao precisa de alocacao.
#
#     screen -S fabio
#     rieux_ampliseq.sh --name fabio --input split/samplesheets/fabio
#     # Ctrl-A depois D
#
# As regioes rodam EM SEQUENCIA, de proposito: duas regioes do mesmo conjunto
# so disputariam a mesma fila, sem terminar antes. Dois CONJUNTOS em paralelo e
# outra coisa — ai sao particoes diferentes, e ai ganha-se de verdade:
#
#     rieux_ampliseq.sh --name fabio    --input .../fabio                    # cpu
#     rieux_ampliseq.sh --name patricia --input .../patricia \
#                       --partition fat --work-dir exec/patricia
#
# O --work-dir nao e opcional nesse caso. Ver DIRETORIO DE LANCAMENTO abaixo.
#
# E idempotente: cada regiao usa -resume e um outdir proprio, entao
# reexecutar depois de uma falha retoma de onde parou.
#
# ---------------------------------------------------------------------------
# DIRETORIO DE LANCAMENTO
#
# O Nextflow guarda estado no diretorio de onde foi chamado: `.nextflow.log`,
# `.nextflow/history` e o cache do `-resume` em `.nextflow/cache/`. Dois
# drivers rodando ao mesmo tempo no MESMO diretorio disputam esses arquivos — o
# segundo a subir rotaciona o `.nextflow.log` do primeiro, e o historico fica
# embaralhado. Nao corrompe resultado (o hash de cada tarefa continua unico no
# work/), mas destroi a capacidade de depurar quando algo falha.
#
# Por isso: execucoes simultaneas, cada uma no seu --work-dir.
#
# O padrao continua sendo o diretorio atual, e nao por acaso: mudar o diretorio
# de lancamento de uma execucao que ja rodou invalida o cache dela, e um
# -resume posterior refaria tudo do zero.
# ---------------------------------------------------------------------------

set -uo pipefail

AQUI="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Ctrl-C tem que parar TUDO, nao so a regiao da vez.
#
# Sem isto, o SIGINT mata apenas o nextflow em primeiro plano; o `for` do bash
# sobrevive e alegremente comeca a proxima regiao. Ja aconteceu aqui: o usuario
# interrompeu uma execucao e ela "continuou sozinha".
trap 'echo; echo "interrompido pelo usuario."; exit 130' INT

# ------------------------------------------------------------------ padroes
NOME=""
ENTRADA=""
SAIDA=""
PARTICAO="cpu"
CONFIG="$AQUI/conf/rieux.config"
PRIMERS="$AQUI/assets/primers_painel.tsv"
PARAMS="$AQUI/assets/parametros_regioes.tsv"
# Padrao e o nome no nf-core: o Nextflow baixa e versiona sozinho. Um caminho
# local so e necessario em cluster sem internet nos nos — e ai o usuario passa
# --pipeline. Cravar aqui o caminho de alguem quebra para todo mundo mais.
PIPE="${AMPLISEQ_HOME:-nf-core/ampliseq}"
DIR_EXEC="$PWD"
MIN_AMOSTRAS=3
MULTIREGION=""
SIDLE_ENTRADA=""
SIDLE_REF="silva"
SIDLE_EXTRA=""
SO_SIDLE=0
PULAR_SIDLE=0
PULAR_METRICAS=0
SIMULAR=0
REGIOES=()

ajuda() {
    sed -n '3,50p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'
    cat <<'FIM'

Opcoes
------
  --name NOME          rotulo do conjunto (colaborador, projeto, lote).
                       Entra nos caminhos de saida e nos logs.  [obrigatorio]
  --input DIR          diretorio com um samplesheet por regiao, no formato
                       samplesheet_<REGIAO>.tsv                 [obrigatorio]
  --outdir DIR         raiz da saida (padrao: ./resultados/<NOME>)
  --partition NOME     particao do SLURM (padrao: cpu)
  --regions "A B"      so estas regioes (padrao: todas as do --params)
  --work-dir DIR       diretorio de lancamento do Nextflow (padrao: o atual)

  --primers ARQ        tabela regiao/forward/reverse
  --params ARQ         tabela regiao/trunclenf/trunclenr/banco/extra
  --config ARQ         perfil do cluster
  --pipeline DIR       caminho do nf-core/ampliseq
  --min-samples N      pula regiao com menos amostras que isto (padrao: 3)

  --multiregion ARQ    liga o ramo Sidle. O arquivo e o regions_multiregion.tsv
                       gerado por bin/validar_regioes_sidle.py
  --sidle-input ARQ    samplesheet das reads NAO divididas (ver nota abaixo)
  --sidle-ref NOME     banco do Sidle (padrao: silva; ver --help do ampliseq)
  --sidle-extra "..."  argumentos extras so para o ramo Sidle
  --only-sidle         nao roda o ramo por regiao
  --skip-sidle         nao roda o ramo Sidle mesmo com --multiregion
  --skip-metrics       nao consolida final_reports/ ao terminar

  --dry-run            so mostra os comandos
  -h, --help           esta ajuda

Nota sobre --sidle-input
------------------------
O ramo Sidle recebe as reads INTEIRAS, nao as ja divididas por regiao: ele faz
o proprio roteamento, rodando cutadapt uma vez por regiao com os primers do
regions_multiregion.tsv. Alimentar o Sidle com as reads ja divididas cortaria
primer duas vezes e deslocaria as bordas — que e exatamente o que o
validar_regioes_sidle.py existe para impedir.

Saida
-----
  <outdir>/<REGIAO>/        uma execucao do ampliseq por regiao
  <outdir>/sidle/           o ramo multi-regiao, se ligado
  <outdir>/final_reports/   as tabelas consolidadas — o contrato com o `aspp`
FIM
}

# ------------------------------------------------------------------ opcoes
while [[ $# -gt 0 ]]; do
    case "$1" in
        --name)         NOME="${2:?--name precisa de um valor}"; shift 2 ;;
        --name=*)       NOME="${1#*=}"; shift ;;
        --input)        ENTRADA="${2:?--input precisa de um caminho}"; shift 2 ;;
        --input=*)      ENTRADA="${1#*=}"; shift ;;
        --outdir)       SAIDA="${2:?--outdir precisa de um caminho}"; shift 2 ;;
        --outdir=*)     SAIDA="${1#*=}"; shift ;;
        --partition)    PARTICAO="${2:?--partition precisa de um valor}"; shift 2 ;;
        --partition=*)  PARTICAO="${1#*=}"; shift ;;
        --regions)      read -r -a REGIOES <<< "${2:?}"; shift 2 ;;
        --regions=*)    read -r -a REGIOES <<< "${1#*=}"; shift ;;
        --work-dir)     DIR_EXEC="${2:?}"; shift 2 ;;
        --work-dir=*)   DIR_EXEC="${1#*=}"; shift ;;
        --primers)      PRIMERS="${2:?}"; shift 2 ;;
        --primers=*)    PRIMERS="${1#*=}"; shift ;;
        --params)       PARAMS="${2:?}"; shift 2 ;;
        --params=*)     PARAMS="${1#*=}"; shift ;;
        --config)       CONFIG="${2:?}"; shift 2 ;;
        --config=*)     CONFIG="${1#*=}"; shift ;;
        --pipeline)     PIPE="${2:?}"; shift 2 ;;
        --pipeline=*)   PIPE="${1#*=}"; shift ;;
        --min-samples)  MIN_AMOSTRAS="${2:?}"; shift 2 ;;
        --multiregion)  MULTIREGION="${2:?}"; shift 2 ;;
        --multiregion=*) MULTIREGION="${1#*=}"; shift ;;
        --sidle-input)  SIDLE_ENTRADA="${2:?}"; shift 2 ;;
        --sidle-ref)    SIDLE_REF="${2:?}"; shift 2 ;;
        --sidle-extra)  SIDLE_EXTRA="${2:?}"; shift 2 ;;
        --only-sidle)   SO_SIDLE=1; shift ;;
        --skip-sidle)   PULAR_SIDLE=1; shift ;;
        --skip-metrics) PULAR_METRICAS=1; shift ;;
        --dry-run|--simular) SIMULAR=1; shift ;;
        -h|--help)      ajuda; exit 0 ;;
        -*)             echo "ERRO: opcao desconhecida: $1" >&2; exit 1 ;;
        *)              REGIOES+=("$1"); shift ;;
    esac
done

[[ -n "$NOME" ]]    || { echo "ERRO: --name e obrigatorio (use --help)" >&2; exit 1; }
[[ -n "$ENTRADA" ]] || { echo "ERRO: --input e obrigatorio (use --help)" >&2; exit 1; }
[[ -n "$SAIDA" ]]   || SAIDA="$PWD/resultados/$NOME"

# caminhos relativos resolvem ANTES do cd, senao viram outro arquivo
for v in ENTRADA SAIDA DIR_EXEC CONFIG PRIMERS PARAMS MULTIREGION SIDLE_ENTRADA; do
    val="${!v}"
    [[ -z "$val" || "$val" == /* ]] || printf -v "$v" '%s' "$PWD/$val"
done

# ------------------------------------------------------------------ ambiente
# Carrega so se ainda nao foi: `source ambiente.sh` na mao continua valendo, e
# recarregar e inofensivo mas polui a tela.
if [[ -z "${DB_SILVA_GENERO:-}" && -r "$AQUI/bin/ambiente.sh" ]]; then
    # shellcheck disable=SC1091
    source "$AQUI/bin/ambiente.sh"
fi

# ------------------------------------------------------------------ checagens
for f in "$PRIMERS" "$PARAMS" "$CONFIG"; do
    [[ -r "$f" ]] || { echo "ERRO: nao achei $f" >&2; exit 1; }
done
[[ -d "$ENTRADA" ]] || { echo "ERRO: --input nao e um diretorio: $ENTRADA" >&2; exit 1; }
# so checamos existencia quando for caminho; nome de pipeline o Nextflow resolve
if [[ "$PIPE" == /* || "$PIPE" == ./* ]]; then
    [[ -d "$PIPE" ]] || { echo "ERRO: nao achei o ampliseq em $PIPE (use --pipeline)" >&2; exit 1; }
fi
command -v nextflow >/dev/null 2>&1 || {
    echo "ERRO: nextflow fora do PATH — rode 'source $AQUI/bin/ambiente.sh'" >&2; exit 1; }
[[ -n "${DB_SILVA_GENERO:-}" ]] || {
    echo "ERRO: DB_SILVA_GENERO vazia — rode 'source $AQUI/bin/ambiente.sh'" >&2; exit 1; }

# A particao pedida precisa existir e aceitar jobs. Descobrir isso agora custa
# um segundo; descobrir depois custa a execucao inteira falhando na primeira
# submissao, horas dentro do `screen`.
if command -v sinfo >/dev/null 2>&1; then
    if ! sinfo -h -p "$PARTICAO" -o '%P' 2>/dev/null | grep -q .; then
        echo "ERRO: a particao '$PARTICAO' nao existe neste cluster" >&2
        echo "      particoes disponiveis:" >&2
        sinfo -h -o '        %P  %D nos  %m MB' >&2
        exit 1
    fi
fi

# ------------------------------------------------------------------ particao
# As tres linhas que separavam o rieux.config do rieux_fat.config viraram isto.
export RIEUX_QUEUE="$PARTICAO"
case "$PARTICAO" in
    fat)
        # 2 nos, ~4 TB. Menos nos => fila menor, senao so se enfileira; mais
        # memoria => process_high pode crescer de verdade.
        export RIEUX_QUEUESIZE="${RIEUX_QUEUESIZE:-40}"
        export RIEUX_MEM_HIGH="${RIEUX_MEM_HIGH:-256}"
        export RIEUX_QUEUE_TAX="${RIEUX_QUEUE_TAX:-fat}"
        ;;
    *)
        export RIEUX_QUEUESIZE="${RIEUX_QUEUESIZE:-60}"
        export RIEUX_MEM_HIGH="${RIEUX_MEM_HIGH:-160}"
        # taxonomia com SILVA e o unico passo que realmente precisa da fat
        export RIEUX_QUEUE_TAX="${RIEUX_QUEUE_TAX:-fat}"
        ;;
esac

mkdir -p "$SAIDA" "$SAIDA/logs" "$DIR_EXEC"
cd "$DIR_EXEC" || { echo "ERRO: nao consegui entrar em $DIR_EXEC" >&2; exit 1; }

if [[ ${#REGIOES[@]} -eq 0 ]]; then
    mapfile -t REGIOES < <(awk -F'\t' '!/^#/ && !/^regiao\t/ && NF>=4 {print $1}' "$PARAMS")
fi

echo "Conjunto   : $NOME"
echo "Regioes    : ${REGIOES[*]}"
echo "Particao   : $PARTICAO  (fila $RIEUX_QUEUESIZE, process_high ${RIEUX_MEM_HIGH} GB)"
echo "Entrada    : $ENTRADA"
echo "Lancamento : $DIR_EXEC"
echo "Saida      : $SAIDA"
[[ -n "$MULTIREGION" ]] && echo "Sidle      : $MULTIREGION"
echo

# ------------------------------------------------------------------ por regiao
declare -a OK=() FALHOU=() PULADO=()

executar() {   # executar <rotulo> <log> <args...>
    local rotulo="$1" log="$2"; shift 2
    if (( SIMULAR )); then
        printf '    nextflow'; printf ' %q' "$@"; echo; echo
        return 0
    fi
    local inicio=$SECONDS
    if nextflow "$@" 2>&1 | tee -a "$log"; then
        echo "    concluido em $(( (SECONDS-inicio)/60 )) min"
        return 0
    fi
    echo "    FALHOU — veja $log"
    return 1
}

if (( ! SO_SIDLE )); then
for REG in "${REGIOES[@]}"; do
    SS="$ENTRADA/samplesheet_${REG}.tsv"
    if [[ ! -s "$SS" ]]; then
        echo "[$REG] sem samplesheet — pulando"; PULADO+=("$REG"); continue
    fi
    N=$(( $(wc -l < "$SS") - 1 ))
    if (( N < MIN_AMOSTRAS )); then
        echo "[$REG] so $N amostra(s) — pulando"; PULADO+=("$REG:$N"); continue
    fi

    read -r FW RV < <(awk -F'\t' -v r="$REG" '$1==r {print $2" "$3}' "$PRIMERS")
    if [[ -z "${FW:-}" || -z "${RV:-}" ]]; then
        echo "[$REG] primers ausentes em $PRIMERS — pulando"; PULADO+=("$REG:primers"); continue
    fi

    read -r TF TR BANCO EXTRA < <(awk -F'\t' -v r="$REG" \
        '$1==r {printf "%s %s %s %s", $2, $3, $4, ($5==""?"-":$5)}' "$PARAMS")
    [[ "$EXTRA" == "-" ]] && EXTRA=""

    ARGS=(run "$PIPE" -profile singularity -c "$CONFIG" -resume
          --input "$SS"
          --outdir "$SAIDA/$REG"
          --FW_primer "$FW"
          --RV_primer "$RV"
          --illumina_novaseq
          --dada_min_boot 80)

    [[ "$TF" != "0" && "$TR" != "0" ]] && ARGS+=(--trunclenf "$TF" --trunclenr "$TR")

    if [[ "$BANCO" == "unite" ]]; then
        # DOIS arquivos, como no SILVA — e nao --skip_dada_addspecies, que foi
        # o primeiro palpite aqui e estava errado.
        #
        # O bin/taxref_reformat_unite.sh do proprio ampliseq deriva DOIS FASTA
        # do mesmo download do UNITE: um para o assignTaxonomy e outro, no
        # formato ">ID Genero especie", para o addSpecies. O caminho oficial
        # (--dada_ref_taxonomy unite-fungi) usa os dois.
        #
        # E a armadilha central: --dada_ref_tax_custom NAO roda o fmtscript que
        # o --dada_ref_taxonomy roda. Quem passa o banco na mao passa o banco
        # JA formatado, ou o passa cru sem receber aviso nenhum.
        if [[ -z "${DB_UNITE_SP:-}" || ! -r "${DB_UNITE_SP:-/dev/null}" ]]; then
            echo "[$REG] DB_UNITE_SP ausente — rode 'python3 $AQUI/bin/preparar_unite.py'" >&2
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

    LOG="$SAIDA/logs/${NOME}_${REG}.log"
    echo "=== [$REG] $N amostras | truncLen ${TF}/${TR} | banco $BANCO"
    echo "    log: $LOG"
    if executar "$REG" "$LOG" "${ARGS[@]}"; then OK+=("$REG"); else FALHOU+=("$REG"); fi
    echo
done
fi

# ------------------------------------------------------------------ Sidle
if [[ -n "$MULTIREGION" ]] && (( ! PULAR_SIDLE )); then
    if [[ ! -r "$MULTIREGION" ]]; then
        echo "ERRO: nao achei $MULTIREGION" >&2
        FALHOU+=("sidle:arquivo")
    elif [[ -z "$SIDLE_ENTRADA" ]]; then
        echo "ERRO: o ramo Sidle precisa de --sidle-input (reads NAO divididas)." >&2
        echo "      Ele faz o proprio roteamento por regiao; alimentar com as" >&2
        echo "      reads ja divididas cortaria primer duas vezes." >&2
        FALHOU+=("sidle:entrada")
    else
        # Nao passamos --trunclenf/--trunclenr aqui de proposito: eles sao
        # globais, e no multi-regiao cada regiao tem o seu. Quem governa o corte
        # e o region_length do proprio regions_multiregion.tsv, escolhido a
        # partir da distribuicao observada pelo validar_regioes_sidle.py.
        ARGS=(run "$PIPE" -profile singularity -c "$CONFIG" -resume
              --input "$SIDLE_ENTRADA"
              --multiregion "$MULTIREGION"
              --outdir "$SAIDA/sidle"
              --illumina_novaseq
              --sidle_ref_taxonomy "$SIDLE_REF")
        # shellcheck disable=SC2206
        [[ -n "$SIDLE_EXTRA" ]] && ARGS+=($SIDLE_EXTRA)

        LOG="$SAIDA/logs/${NOME}_sidle.log"
        echo "=== [sidle] $(( $(wc -l < "$MULTIREGION") - 1 )) regioes | banco $SIDLE_REF"
        echo "    log: $LOG"
        if executar sidle "$LOG" "${ARGS[@]}"; then OK+=("sidle"); else FALHOU+=("sidle"); fi
        echo
    fi
fi

# ------------------------------------------------------------------ consolidacao
if (( ! PULAR_METRICAS )) && (( ! SIMULAR )) && [[ ${#OK[@]} -gt 0 ]]; then
    echo "Consolidando final_reports/ ..."
    python3 "$AQUI/bin/coletar_metricas.py" --resultados "$SAIDA" \
            --out "$SAIDA/final_reports" --nome "$NOME" || \
        echo "  (a consolidacao falhou; os resultados por regiao estao intactos)"
    echo
fi

# ------------------------------------------------------------------ resumo
echo "================================================================"
echo "Conjunto $NOME"
# Em --dry-run nada rodou: chamar de "concluida" uma regiao que so foi impressa
# e o tipo de mentira pequena que faz alguem achar que ja tem resultado.
if (( SIMULAR )); then
    [[ ${#OK[@]} -gt 0 ]] && echo "  simuladas : ${OK[*]}"
else
    [[ ${#OK[@]} -gt 0 ]] && echo "  concluidas: ${OK[*]}"
fi
[[ ${#PULADO[@]} -gt 0 ]] && echo "  puladas   : ${PULADO[*]}"
[[ ${#FALHOU[@]} -gt 0 ]] && echo "  falharam  : ${FALHOU[*]}"
echo
echo "  final_reports: $SAIDA/final_reports"
echo "  no R:  dados <- aspp::read_ampliseq_summary('$SAIDA/final_reports')"
[[ ${#FALHOU[@]} -gt 0 ]] && exit 1 || exit 0

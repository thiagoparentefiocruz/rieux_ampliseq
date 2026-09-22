#!/usr/bin/env bash
#SBATCH --job-name=split_regioes
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --output=logs/split_%A_%a.log
#
# split_regions.sh — separa as reads de cada amostra por regiao do painel
#
# O PROBLEMA QUE ELE RESOLVE
# Cada amostra tem UM par de FASTQ com os amplicons das 7 regioes
# misturados. Nenhum pipeline de 16S espera isso: todos assumem uma regiao
# por biblioteca. Este modulo le o par R1/R2 de uma amostra, decide a qual
# regiao cada PAR de reads pertence, e escreve em arquivos separados.
#
# COMO DECIDE
# Pelo par de primers, exigindo os dois lados: primer forward casando no R1
# E o reverse correspondente casando no R2 (--pair-adapters). Exigir os dois
# nao e preciosismo — o 805F (forward do V5V7) e quase o reverso-complemento
# do 806R (reverse do V3V4), e nos dados isso produz 10,9% de reads que, se
# atribuidas so pelo R1, iriam para a regiao errada.
#
# NAO CORTA NADA (--action=none). Apenas roteia. O corte de primer fica com
# o proprio ampliseq (--FW_primer/--RV_primer), que ja e testado, ja gera
# relatorio de QC e remove junto o bloco de fase, porque o cutadapt elimina
# o adaptador 5' e tudo que vem antes dele.
#
# SAIDA
#   split/<regiao>/<amostra>_R{1,2}.fastq.gz   reads daquela regiao
#   split/unknown/<amostra>_R{1,2}.fastq.gz    nao atribuidas
#   relatorios/<amostra>.cutadapt.json         metricas por amostra
#
# USO
#   mkdir -p logs
#   # 1. descobre quantos indices sao necessarios (roda fora do SLURM):
#   bash split_regions.sh <brutos> <primers.tsv> <saida>
#   # 2. submete com o --array que ele indicar:
#   sbatch --array=1-30%10 split_regions.sh <brutos> <primers.tsv> <saida>
#
# LOTE controla quantas amostras cada tarefa processa (padrao 10). Aumente
# se o QOS do cluster limitar jobs submetidos:
#   LOTE=30 sbatch --export=ALL,LOTE=30 --array=1-10 split_regions.sh ...
#
# E idempotente: amostra que ja tem relatorio JSON e pulada, entao dá para
# resubmeter depois de uma falha parcial sem refazer trabalho.

set -euo pipefail

BRUTOS="${1:?informe o diretorio dos dados brutos}"
TABELA="${2:?informe o primers_panel.tsv}"
SAIDA="${3:-split_saida}"

THREADS="${SLURM_CPUS_PER_TASK:-4}"
IDX="${SLURM_ARRAY_TASK_ID:-1}"

CACHE="${NXF_SINGULARITY_CACHEDIR:-${RIEUX_PIPELINE_BASE:-$HOME}/singularity}"

# ---------------------------------------------------------------- cutadapt
# O cluster nao tem cutadapt instalado, mas o cache do ampliseq tem a imagem.
# Reaproveitar garante que o split use exatamente a mesma versao que o
# pipeline usara depois.
IMG=$(find "$CACHE" -iname '*cutadapt*' \( -name '*.img' -o -name '*.sif' \) -print -quit 2>/dev/null || true)
if command -v cutadapt >/dev/null 2>&1; then
    CUTADAPT=(cutadapt)
elif [[ -n "$IMG" ]]; then
    # Os diretorios de amostra sao SYMLINKS para outro filesystem, e o home
    # tambem e link. O Singularity so enxerga o que esta montado: se a gente
    # bindar apenas o caminho aparente, o container abre o link e nao acha o
    # alvo. Entao resolvemos os caminhos de verdade e bindamos as raizes.
    declare -A RAIZES=()
    primeiro_r1=$(find -L "$BRUTOS" -name '*_R1_*.fastq.gz' -print -quit 2>/dev/null || true)
    for alvo in "$BRUTOS" "$SAIDA" "$TABELA" "$primeiro_r1" "$PWD"; do
        [[ -e "$alvo" ]] || continue
        real=$(readlink -f "$alvo")
        raiz="/$(echo "${real#/}" | cut -d/ -f1)"
        [[ -d "$raiz" ]] && RAIZES["$raiz"]=1
    done
    BINDS=()
    for r in "${!RAIZES[@]}"; do BINDS+=(-B "$r"); done
    echo "  binds do singularity: ${!RAIZES[*]}"
    CUTADAPT=(singularity exec "${BINDS[@]}" "$IMG" cutadapt)
else
    echo "ERRO: nem cutadapt no PATH nem imagem no cache ($CACHE)." >&2
    echo "Imagens disponiveis:" >&2
    find "$CACHE" \( -name '*.img' -o -name '*.sif' \) -printf '  %f\n' 2>/dev/null | head -20 >&2
    exit 1
fi

# ---------------------------------------------------- FASTA dos primers
# Dois arquivos, na MESMA ordem: --pair-adapters casa o n-esimo de -g com o
# n-esimo de -G. Ordem trocada = regiao trocada, em silencio.
#
# Os FASTA sao gerados num diretorio TEMPORARIO por tarefa, nao em $SAIDA.
# Com 40 tarefas do array comecando ao mesmo tempo, escrever o mesmo arquivo
# concorrentemente produz FASTA truncado ou intercalado — e o sintoma seria
# uma regiao sumindo em algumas amostras e nao em outras, o pior tipo de bug
# para diagnosticar depois.
mkdir -p "$SAIDA/relatorios" "$SAIDA/split"
TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT
F_FA="$TMP/forward.fasta"
R_FA="$TMP/reverse.fasta"

awk -F'\t' '!/^#/ && !/^regiao\t/ && NF>=3 {
        print ">"$1"\n"$2 > "'"$F_FA"'"
        print ">"$1"\n"$3 > "'"$R_FA"'"
        print $1
     }' "$TABELA" > "$TMP/regioes.txt"

mapfile -t REGIOES < "$TMP/regioes.txt"
[[ ${#REGIOES[@]} -gt 0 ]] || { echo "ERRO: nenhuma regiao lida de $TABELA" >&2; exit 1; }
for r in "${REGIOES[@]}" unknown; do mkdir -p "$SAIDA/split/$r"; done

# --------------------------------------------------- amostras desta tarefa
# Cada tarefa processa LOTE amostras em sequencia, em vez de uma. Motivo:
# o QOS do cluster limita jobs SUBMETIDOS por usuario, e um array conta cada
# indice como um job — 297 indices estouram o limite mesmo sem nada rodando.
# Com LOTE=10, 297 amostras cabem em 30 indices.
LOTE="${LOTE:-10}"

mapfile -t R1S < <(find -L "$BRUTOS" -name '*_R1_*.fastq.gz' | sort)
TOTAL=${#R1S[@]}
NTAREFAS=$(( (TOTAL + LOTE - 1) / LOTE ))

if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    echo "$TOTAL amostras, LOTE=$LOTE  ->  submeta com --array=1-$NTAREFAS"
    echo "  sbatch --array=1-$NTAREFAS%10 $0 $BRUTOS $TABELA $SAIDA"
    exit 0
fi

INICIO=$(( (IDX - 1) * LOTE ))
FIM=$(( INICIO + LOTE - 1 ))
(( FIM >= TOTAL )) && FIM=$(( TOTAL - 1 ))
(( INICIO >= TOTAL )) && { echo "tarefa $IDX nao tem amostras"; exit 0; }

echo "tarefa $IDX/$NTAREFAS: amostras $((INICIO+1)) a $((FIM+1)) de $TOTAL"

# -------------------------------------------------------------- executa
# -e 0.15   ~3 erros num primer de 20 nt, folga para a degenerescencia
# --no-indels  primer sintetico nao tem indel; permitir so gera falso match
# sem ancora (^): o bloco de fase desloca o primer de 0 a 11 bases
for (( k = INICIO; k <= FIM; k++ )); do
    R1="${R1S[$k]}"
    R2="${R1/_R1_/_R2_}"
    AMOSTRA=$(basename "$R1" | sed -E 's/_S[0-9]+_L[0-9]+_R1.*//')
    JSON="$SAIDA/relatorios/${AMOSTRA}.cutadapt.json"

    if [[ ! -f "$R2" ]]; then
        echo "  AVISO: $AMOSTRA sem par R2 — pulando" >&2
        continue
    fi
    # idempotente: reexecutar depois de uma falha parcial nao refaz o que ja
    # terminou. Apague o JSON da amostra para forcar o reprocessamento.
    if [[ -s "$JSON" ]]; then
        echo "  $AMOSTRA ja processada — pulando"
        continue
    fi

    echo "  $AMOSTRA"
    "${CUTADAPT[@]}" \
        -j "$THREADS" \
        -e 0.15 \
        --no-indels \
        --pair-adapters \
        --action=none \
        -g "file:$F_FA" \
        -G "file:$R_FA" \
        -o "$SAIDA/split/{name}/${AMOSTRA}_R1.fastq.gz" \
        -p "$SAIDA/split/{name}/${AMOSTRA}_R2.fastq.gz" \
        --json "$JSON" \
        "$R1" "$R2" > /dev/null
done

echo "tarefa $IDX concluida"

#!/usr/bin/env bash
#
# rieux_ampliseq.sh — nf-core/ampliseq for multi-region amplicon panels
#
# Takes a multi-region amplicon panel from raw FASTQ to consolidated tables, in
# seven stages. Any stage can be run on its own, and the chain can be entered or
# left at any point.
#
#     organize -> discover -> split -> profile -> run -> sidle -> collect
#
#   organize   cross a sample table with the FASTQs on disk, one project per
#              group. Optional: if your FASTQs are already separated, start at
#              `discover`.
#   discover   recover the panel's primers from the reads themselves.
#              Optional: if you know your primers, pass --primers.
#   split      route each read pair to its region, requiring the primer PAIR,
#              and write one samplesheet per region.
#   profile    measure real per-cycle quality and pick truncLenF/R per region.
#              Comes AFTER split, because truncation is chosen per region.
#   run        run ampliseq once per region.
#   sidle      run the multi-region branch (reconstruction against a reference).
#   collect    write final_reports/, the contract with the R package `aspp`.
#
# From scratch, one command:
#
#     rieux_ampliseq.sh --project renata --raw-dir raw/renata
#
# Re-entering in the middle, when what came before already exists:
#
#     rieux_ampliseq.sh --project renata --from run
#     rieux_ampliseq.sh --project renata --stage collect
#
# ---------------------------------------------------------------------------
# THE PROJECT DIRECTORY
#
# Stages do not talk to each other through flags: they talk through a directory
# with a known layout. That is what makes re-entry possible without re-stating
# everything that came before.
#
#     <project>/
#       raw/                   FASTQ (symlinks), from `organize`
#       metadata.tsv           sample -> group, from `organize`
#       primers.tsv            from `discover`, or copied from --primers
#       region_params.tsv      truncLen per region, from `profile`
#       split/                 routed reads, from `split`
#         samplesheets/        one per region, plus the complete one
#       <REGION>/              one ampliseq run, from `run`
#       sidle/                 from `sidle`
#       final_reports/         from `collect`
#       logs/
#
# Any of these paths can be overridden by a flag. The layout is the default,
# not a requirement.
#
# ---------------------------------------------------------------------------
# WHERE TO RUN IT
#
# The Nextflow driver runs on the LOGIN NODE, inside a `screen`. It sits idle
# waiting on SLURM and submits the tasks; it needs no allocation of its own.
#
# Regions run IN SEQUENCE, on purpose: two regions of the same project would
# only compete for the same queue. Two PROJECTS in parallel is a different
# matter — different partitions, and there the gain is real:
#
#     rieux_ampliseq.sh --project fabio    --from run                     # cpu
#     rieux_ampliseq.sh --project patricia --from run \
#                       --partition fat --work-dir exec/patricia
#
# --work-dir is not optional in that case: Nextflow keeps `.nextflow.log`,
# `.nextflow/history` and the -resume cache in the launch directory, and two
# drivers in the same directory shred each other's history. It does not corrupt
# results — each task's hash stays unique under work/ — but it destroys any
# chance of debugging a failure.
# ---------------------------------------------------------------------------

set -uo pipefail

AQUI="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Ctrl-C tem que parar TUDO, nao so o estagio da vez.
#
# Sem isto o SIGINT mata apenas o processo em primeiro plano; o laco do bash
# sobrevive e alegremente comeca o proximo. Ja aconteceu aqui: uma execucao
# interrompida "continuou sozinha".
trap 'echo; echo "interrupted by user."; exit 130' INT

ESTAGIOS=(organize discover split profile run sidle collect)

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
    cat <<'END'

Options
-------
  --project NAME       project label: goes into paths and logs. Also the
                       working directory, unless --outdir is given.
  --outdir DIR         project root (default: ./<project>)

  Stage selection (with none of these, the whole chain runs):
  --stage NAME         run ONLY this stage
  --from NAME          from this stage to the end
  --until NAME         from the start to this stage
  --skip "A B"         skip these

  Inputs (each defaults to a path inside the project directory):
  --raw-dir DIR        raw FASTQs              (default <project>/raw)
  --sample-table FILE  sample table            (stage `organize` only)
  --primers FILE       region/forward/reverse table. Pass it if you already
                       know your panel's primers — skips stage `discover`.
  --params FILE        region/trunclenf/trunclenr/database/extra table. Pass it
                       if you already chose truncLen — skips stage `profile`.
  --controls REGEX     sample names that are controls (default: ^[Ss]mart)

  Execution:
  --partition NAME     SLURM partition (default: cpu)
  --regions "A B"      only these regions
  --work-dir DIR       Nextflow launch directory (default: current)
  --config FILE        cluster profile
  --pipeline X         ampliseq name or path (default: nf-core/ampliseq)
  --min-samples N      skip a region with fewer samples than this (default 3)
  --min-reads N        floor of reads per sample PER REGION in split (default 1000)
  --batch-size N       samples per task of the split array (default 10)

  Sidle branch:
  --multiregion FILE   regions_multiregion.tsv from bin/validate_sidle_regions.py
  --sidle-input FILE   samplesheet of the UNDIVIDED reads
  --sidle-ref NAME     Sidle reference database (default: silva)
  --sidle-extra "..."  extra arguments for this branch only

  --dry-run            print the commands, run nothing
  -h, --help           this help

About the `organize` stage
--------------------------
It is the only ONE-TO-MANY stage: one sample table becomes several projects. So
it does not chain — it runs, writes the projects and prints the next command for
each one. Every other stage is one-to-one.

About the Sidle branch
----------------------
It is fed the UNDIVIDED reads, not the split ones: it does its own routing,
running cutadapt once per region with the primers from regions_multiregion.tsv.
Feeding it already-split reads would trim the primers twice and shift the
boundaries — exactly what validate_sidle_regions.py exists to prevent.
END
}

# ------------------------------------------------------------------ opcoes
while [[ $# -gt 0 ]]; do
    case "$1" in
        --project|--name) PROJETO="${2:?}"; shift 2 ;;
        --project=*)    PROJETO="${1#*=}"; shift ;;
        --outdir)       RAIZ="${2:?}"; shift 2 ;;
        --outdir=*)     RAIZ="${1#*=}"; shift ;;
        --stage)        SO="${2:?}"; shift 2 ;;
        --from)         DE="${2:?}"; shift 2 ;;
        --until)        ATE="${2:?}"; shift 2 ;;
        --skip)         PULAR="${2:?}"; shift 2 ;;
        --raw-dir)       BRUTOS="${2:?}"; shift 2 ;;
        --sample-table)     PLANILHA="${2:?}"; shift 2 ;;
        --primers)      PRIMERS="${2:?}"; shift 2 ;;
        --params)       PARAMS="${2:?}"; shift 2 ;;
        --controls)    CONTROLES="${2:?}"; shift 2 ;;
        --partition)    PARTICAO="${2:?}"; shift 2 ;;
        --regions)      read -r -a REGIOES <<< "${2:?}"; shift 2 ;;
        --work-dir)     DIR_EXEC="${2:?}"; shift 2 ;;
        --config)       CONFIG="${2:?}"; shift 2 ;;
        --pipeline)     PIPE="${2:?}"; shift 2 ;;
        --min-samples)  MIN_AMOSTRAS="${2:?}"; shift 2 ;;
        --min-reads)    MIN_READS_REGIAO="${2:?}"; shift 2 ;;
        --batch-size)         LOTE_SPLIT="${2:?}"; shift 2 ;;
        --multiregion)  MULTIREGION="${2:?}"; shift 2 ;;
        --sidle-input)  SIDLE_ENTRADA="${2:?}"; shift 2 ;;
        --sidle-ref)    SIDLE_REF="${2:?}"; shift 2 ;;
        --sidle-extra)  SIDLE_EXTRA="${2:?}"; shift 2 ;;
        --dry-run|--simular) SIMULAR=1; shift ;;
        -h|--help)      ajuda; exit 0 ;;
        *)              echo "ERROR: unknown option: $1" >&2; exit 1 ;;
    esac
done

[[ -n "$PROJETO" ]] || { echo "ERROR: --project is required (see --help)" >&2; exit 1; }

# ------------------------------------------------- layout do projeto
[[ -n "$RAIZ" ]] || RAIZ="$PWD/$PROJETO"
[[ "$RAIZ" == /* ]] || RAIZ="$PWD/$RAIZ"

[[ -n "$BRUTOS" ]]  || BRUTOS="$RAIZ/brutos"
[[ -n "$PRIMERS" ]] || PRIMERS="$RAIZ/primers.tsv"
[[ -n "$PARAMS" ]]  || PARAMS="$RAIZ/region_params.tsv"
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
    echo "ERROR: unknown stage: $alvo" >&2
    echo "      valid: ${ESTAGIOS[*]}" >&2
    return 1
}

if [[ -n "$SO" ]]; then
    indice "$SO" || exit 1; INI="$INDICE_RES"; FIM="$INDICE_RES"
else
    # Sem --from, o padrao comeca em `descobrir`, nao em `organizar`: organizar
    # e' um-para-muitos e exige uma planilha, entao entrar nele por acidente
    # seria surpresa. Quem quer a planilha pede por ela.
    indice "${DE:-discover}"   || exit 1; INI="$INDICE_RES"
    indice "${ATE:-collect}" || exit 1; FIM="$INDICE_RES"
fi
(( INI <= FIM )) || { echo "ERROR: --from comes after --until" >&2; exit 1; }

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
        echo "    (--dry-run: $2 does not exist yet; it would come from '$3')"
        return 0
    fi
    echo "ERROR: missing $2" >&2
    echo "      expected at: $1" >&2
    echo "      produced by stage '$3' — run it first, or give the path" >&2
    echo "      with the matching flag." >&2
    exit 1
}

mkdir -p "$RAIZ" "$LOGS"

echo "Project    : $PROJETO"
echo "Directory  : $RAIZ"
echo "Stages     : ${ESTAGIOS[*]:$INI:$((FIM-INI+1))}${PULAR:+  (skipping: $PULAR)}"
echo "Partition  : $PARTICAO"
echo

# ------------------------------------------------- ambiente
if [[ -z "${DB_SILVA_GENERO:-}" && -r "$AQUI/bin/env.sh" ]]; then
    # shellcheck disable=SC1091
    source "$AQUI/bin/env.sh" || exit 1
fi

# =================================================================== organizar
if rodar_estagio organize; then
    echo "== organize ==================================================="
    [[ -n "$PLANILHA" ]] || { echo "ERROR: stage 'organize' needs --sample-table" >&2; exit 1; }
    precisa "$PLANILHA" "the sample table" organize
    precisa "$BRUTOS" "the raw FASTQ directory (--raw-dir)" organize
    exec_cmd python3 "$AQUI/bin/organize_project.py" "$PLANILHA" "$BRUTOS" \
             "$RAIZ/projects" --controls "$CONTROLES" --apply || exit 1
    echo
    echo "  'organize' is one-to-many: it wrote one project per group under"
    echo "  $RAIZ/projects/. Continue one at a time:"
    if (( ! SIMULAR )); then
        for d in "$RAIZ"/projects/*/; do
            [[ -d "$d" ]] || continue
            n=$(basename "$d")
            echo "    rieux_ampliseq.sh --project $n --outdir $d --raw-dir $d/dados_brutos"
        done
    fi
    echo
    exit 0
fi

# =================================================================== descobrir
if rodar_estagio discover; then
    echo "== discover ==================================================="
    if [[ -s "$PRIMERS" ]]; then
        echo "  $PRIMERS already exists — nothing to do."
    else
        precisa "$BRUTOS" "the raw FASTQ directory (--raw-dir)" organize
        exec_cmd python3 "$AQUI/bin/discover_primers.py" --dir "$BRUTOS" \
                 --out "$RAIZ/primers" || exit 1
        echo
        echo "  CHECK $RAIZ/primers.tsv before going on. A wrongly recovered"
        echo "  primer does not raise an error: it gives a wrong result, silently."
    fi
    echo
fi

# ===================================================================== dividir
if rodar_estagio split; then
    echo "== split ====================================================="
    precisa "$PRIMERS" "the primers table" discover
    precisa "$BRUTOS" "the raw FASTQ directory (--raw-dir)" organize
    n_amostras=$(find -L "$BRUTOS" -name '*_R1*.fastq.gz' 2>/dev/null | wc -l)
    n_tarefas=$(( (n_amostras + LOTE_SPLIT - 1) / LOTE_SPLIT ))
    (( n_tarefas > 0 )) || { echo "ERROR: no R1 FASTQ found in $BRUTOS" >&2; exit 1; }
    echo "  $n_amostras samples, batch $LOTE_SPLIT -> array 1-$n_tarefas"
    # --wait bloqueia ate o array terminar. O driver esta num screen no login,
    # entao bloquear e exatamente o comportamento desejado: os estagios
    # seguintes dependem deste.
    exec_cmd sbatch --wait --partition "$PARTICAO" \
             --array="1-${n_tarefas}%10" \
             --export="ALL,LOTE=$LOTE_SPLIT" \
             --output="$LOGS/split_%A_%a.log" \
             "$AQUI/bin/split_regions.sh" "$BRUTOS" "$PRIMERS" "$SPLIT" || exit 1
    exec_cmd python3 "$AQUI/bin/split_summary.py" "$SPLIT" \
             --minimum "$MIN_READS_REGIAO" || exit 1
    echo
fi

# ==================================================================== perfilar
if rodar_estagio profile; then
    echo "== profile ===================================================="
    if [[ -s "$PARAMS" ]]; then
        echo "  $PARAMS already exists — nothing to do."
    else
        precisa "$SPLIT/split" "the reads routed per region" split
        exec_cmd python3 "$AQUI/bin/quality_profile.py" "$SPLIT/split" \
                 --out "$PARAMS" || exit 1
    fi
    echo
fi

# ======================================================================= rodar
if rodar_estagio run; then
    echo "== run ======================================================="
    precisa "$PRIMERS" "the primers table" discover
    precisa "$PARAMS" "the truncLen table" profile
    precisa "$SS_DIR" "the per-region samplesheets" split

    if command -v sinfo >/dev/null 2>&1; then
        sinfo -h -p "$PARTICAO" -o '%P' 2>/dev/null | grep -q . || {
            echo "ERROR: partition '$PARTICAO' does not exist on this cluster" >&2
            sinfo -h -o '        %P  %D nodes  %m MB' >&2; exit 1; }
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
        [[ -s "$SS" ]] || { echo "  [$REG] no samplesheet — skipping"; PULADO+=("$REG"); continue; }
        N=$(( $(wc -l < "$SS") - 1 ))
        (( N >= MIN_AMOSTRAS )) || { echo "  [$REG] only $N sample(s) — skipping"; PULADO+=("$REG:$N"); continue; }

        read -r FW RV < <(awk -F'\t' -v r="$REG" '$1==r {print $2" "$3}' "$PRIMERS")
        [[ -n "${FW:-}" && -n "${RV:-}" ]] || { echo "  [$REG] primers missing — skipping"; PULADO+=("$REG:primers"); continue; }

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
            # bin/prepare_unite.py existe para reproduzir aquele fmtscript, e
            # ele deriva DOIS FASTA: um para o assignTaxonomy e outro, no
            # formato ">ID Genero especie", para o addSpecies.
            if [[ -z "${DB_UNITE_SP:-}" || ! -r "${DB_UNITE_SP:-/dev/null}" ]]; then
                echo "  [$REG] DB_UNITE_SP missing — run '$AQUI/bin/prepare_unite.py'" >&2
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
        echo "  [$REG] $N samples | truncLen ${TF}/${TR} | db $BANCO"
        if (( SIMULAR )); then
            printf '    nextflow'; printf ' %q' "${ARGS[@]}"; echo
            OK+=("$REG"); continue
        fi
        INICIO=$SECONDS
        if nextflow "${ARGS[@]}" 2>&1 | tee -a "$LOG"; then
            echo "    done in $(( (SECONDS-INICIO)/60 )) min"; OK+=("$REG")
        else
            echo "    FAILED — see $LOG"; FALHOU+=("$REG")
        fi
    done
    echo
    if (( SIMULAR )); then
        [[ ${#OK[@]} -gt 0 ]] && echo "  simulated : ${OK[*]}"
    else
        [[ ${#OK[@]} -gt 0 ]] && echo "  done      : ${OK[*]}"
    fi
    [[ ${#PULADO[@]} -gt 0 ]] && echo "  skipped   : ${PULADO[*]}"
    [[ ${#FALHOU[@]} -gt 0 ]] && echo "  failed    : ${FALHOU[*]}"
    echo
fi

# ======================================================================= sidle
if rodar_estagio sidle; then
    if [[ -z "$MULTIREGION" ]]; then
        if [[ -z "$SO" ]]; then
            echo "== sidle: skipped (no --multiregion) ==========================="
            echo "  Generate it with:"
            echo "    $AQUI/bin/validate_sidle_regions.py --primers $PRIMERS \\"
            echo "        --ref \"\$DB_SILVA_GENERO\" \\"
            echo "        --asv $RAIZ/final_reports/asv_length.tsv \\"
            echo "        --out $RAIZ/regions_multiregion.tsv"
            echo
        else
            echo "ERROR: stage 'sidle' needs --multiregion" >&2; exit 1
        fi
    else
        echo "== sidle ======================================================="
        precisa "$MULTIREGION" "the regions_multiregion.tsv" sidle
        precisa "$SIDLE_ENTRADA" "the samplesheet of the UNDIVIDED reads" split
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
        echo "  $(( $(wc -l < "$MULTIREGION") - 1 )) regions | db $SIDLE_REF"
        if (( SIMULAR )); then
            printf '    nextflow'; printf ' %q' "${ARGS[@]}"; echo
        else
            nextflow "${ARGS[@]}" 2>&1 | tee -a "$LOG" || echo "  FAILED — see $LOG"
        fi
        echo
    fi
fi

# ================================================================== consolidar
if rodar_estagio collect; then
    echo "== collect =================================================="
    exec_cmd python3 "$AQUI/bin/collect_metrics.py" \
             --results "$RAIZ" --name "$PROJETO" \
             --out "$RAIZ/final_reports" \
             --primers "$PRIMERS" --controls "$CONTROLES" || exit 1

fi

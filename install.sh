#!/usr/bin/env bash
#
# install.sh — install rieux_ampliseq as a native command
#
#   ./install.sh --mode command      # launcher in ~/.local/bin
#   ./install.sh --mode source       # shell function in your rc file
#   ./install.sh --mode command --base /path/to/pipeline
#   ./install.sh --uninstall
#
# Same two modes as the other rieux_* tools. Two differences, both deliberate:
#
#   - `source` mode defines a shell FUNCTION that calls the script, instead of
#     sourcing the script itself. This wrapper is an executable that uses
#     `exit`, so sourcing it would run it immediately and could close your
#     shell. The function gives the same result honestly.
#
#   - it also records RIEUX_PIPELINE_BASE and checks the prerequisites, because
#     this tool needs a pipeline installation (containers, NXF_HOME, reference
#     databases) and finding that out halfway through a run is expensive.

set -uo pipefail

SCRIPT_NAME="rieux_ampliseq.sh"
COMMAND_NAME="rieux_ampliseq"
DEST_DIR="$HOME/.local/bin"
REPO_DIR="$(cd "$(dirname "$0")" && pwd)"
SCRIPT_PATH="$REPO_DIR/$SCRIPT_NAME"
CONF="$HOME/.rieux_ampliseq.conf"
MARCA="# Added by rieux_ampliseq installer"

usage() {
    cat <<EOF
Usage:
  ./install.sh --mode command [--base DIR]
  ./install.sh --mode source  [--base DIR]
  ./install.sh --uninstall

Modes:
  command   Install a runnable command in ~/.local/bin
  source    Add a shell function to your rc file

Options:
  --base DIR   your pipeline installation (containers, NXF_HOME, bancos.env).
               Asked interactively if omitted; saved to ~/.rieux_ampliseq.conf
  --ampliseq DIR  a LOCAL copy of nf-core/ampliseq. Strongly recommended:
               without it Nextflow pulls whatever is current on GitHub, which
               may require a newer Nextflow than your cluster module has.
  --uninstall  remove the command and the rc lines (keeps the config file)
EOF
}

MODE=""
BASE=""
AMPLISEQ=""
UNINSTALL=0
while [[ $# -gt 0 ]]; do
    case "$1" in
        --mode)      MODE="${2:?--mode needs a value}"; shift 2 ;;
        --base)      BASE="${2:?--base needs a path}"; shift 2 ;;
        --ampliseq)  AMPLISEQ="${2:?--ampliseq needs a path}"; shift 2 ;;
        --base=*)     BASE="${1#*=}"; shift ;;
        --ampliseq=*) AMPLISEQ="${1#*=}"; shift ;;
        --uninstall) UNINSTALL=1; shift ;;
        -h|--help)   usage; exit 0 ;;
        *)           usage; exit 1 ;;
    esac
done

if [[ ! -f "$SCRIPT_PATH" ]]; then
    echo "Error: $SCRIPT_NAME not found in repository root." >&2
    exit 1
fi

detect_shell_rc() {
    if [[ "${SHELL:-}" == *"zsh"* ]]; then
        echo "$HOME/.zshrc"
    elif [[ "${SHELL:-}" == *"bash"* ]]; then
        echo "$HOME/.bashrc"
    else
        echo ""
    fi
}

# Every line this installer adds carries the marker, so uninstall removes
# exactly what was added and a second install does not pile up duplicates.
add_line() {   # add_line <rc> <line>
    local rc="$1" linha="$2"
    [[ -n "$rc" ]] || return 1
    if grep -Fq "$linha" "$rc" 2>/dev/null; then
        echo "Already present in $rc"
        return 0
    fi
    { echo ""; echo "$MARCA"; echo "$linha"; } >> "$rc"
    echo "Added to $rc"
}

ensure_path() {
    local rc="$1"
    mkdir -p "$DEST_DIR"
    if [[ ":$PATH:" == *":$DEST_DIR:"* ]]; then
        echo "$DEST_DIR is already in your PATH."
        return 0
    fi
    if [[ -n "$rc" ]]; then
        add_line "$rc" "export PATH=\"$DEST_DIR:\$PATH\""
    else
        echo "Could not detect bash or zsh automatically."
        echo "Add this line to your shell configuration:"
        echo "export PATH=\"$DEST_DIR:\$PATH\""
    fi
}

# ------------------------------------------------------------------ uninstall
if (( UNINSTALL )); then
    rc="$(detect_shell_rc)"
    rm -f "$DEST_DIR/$COMMAND_NAME" && echo "Removed $DEST_DIR/$COMMAND_NAME"
    if [[ -n "$rc" && -f "$rc" ]] && grep -Fq "$MARCA" "$rc"; then
        cp "$rc" "$rc.rieux_ampliseq.bak"
        # remove a linha marcada E a marca; deixar so a marca seria lixo
        awk -v m="$MARCA" '
            $0 == m {pula=1; next}
            pula    {pula=0; next}
            {print}' "$rc" > "$rc.tmp" && mv "$rc.tmp" "$rc"
        echo "Cleaned $rc (backup: $rc.rieux_ampliseq.bak)"
    fi
    echo "$CONF was NOT removed — it holds your settings."
    exit 0
fi

[[ "$MODE" == "command" || "$MODE" == "source" ]] || { usage; exit 1; }

chmod +x "$SCRIPT_PATH" "$REPO_DIR"/bin/* 2>/dev/null
RC="$(detect_shell_rc)"

# --------------------------------------------------------------------- modes
if [[ "$MODE" == "command" ]]; then
    mkdir -p "$DEST_DIR"
    cat > "$DEST_DIR/$COMMAND_NAME" <<EOF
#!/usr/bin/env bash
exec "$SCRIPT_PATH" "\$@"
EOF
    chmod +x "$DEST_DIR/$COMMAND_NAME"
    echo "Installed command:"
    echo "  $DEST_DIR/$COMMAND_NAME"
    ensure_path "$RC"
else
    if [[ -z "$RC" ]]; then
        echo "Could not detect bash or zsh automatically."
        echo "Add this line to your shell configuration:"
        echo "$COMMAND_NAME() { \"$SCRIPT_PATH\" \"\$@\"; }"
    else
        add_line "$RC" "$COMMAND_NAME() { \"$SCRIPT_PATH\" \"\$@\"; }"
    fi
fi

# ---------------------------------------------------------- pipeline base dir
# O arquivo e lido UMA vez, aqui, antes de qualquer reescrita. Ele e regravado
# do zero mais abaixo; ler depois disso perderia o que ja estava la — foi assim
# que um `--base` explicito apagava um AMPLISEQ_HOME gravado antes.
CONF_BASE=""
CONF_AMPLISEQ=""
if [[ -r "$CONF" ]]; then
    # shellcheck disable=SC1090
    source "$CONF"
    CONF_BASE="${RIEUX_PIPELINE_BASE:-}"
    CONF_AMPLISEQ="${AMPLISEQ_HOME:-}"
fi
if [[ -z "$BASE" && -n "$CONF_BASE" ]]; then
    BASE="$CONF_BASE"
    echo "Reusing RIEUX_PIPELINE_BASE from $CONF"
fi
if [[ -z "$AMPLISEQ" && -n "$CONF_AMPLISEQ" ]]; then
    AMPLISEQ="$CONF_AMPLISEQ"
fi
if [[ -z "$BASE" ]]; then
    echo ""
    echo "Where does your pipeline installation live? That directory holds"
    echo "nextflow_home/, singularity/ and bancos.env. It is NOT the raw data."
    if [[ -t 0 ]]; then
        read -r -p "RIEUX_PIPELINE_BASE: " BASE
    else
        echo "(no terminal: skipped. Re-run with --base DIR when you know it.)"
    fi
fi
[[ -n "$BASE" ]]     && BASE="${BASE/#\~/$HOME}"
[[ -n "$AMPLISEQ" ]] && AMPLISEQ="${AMPLISEQ/#\~/$HOME}"

# O arquivo e escrito de uma vez so, com tudo que sabemos agora. Acrescentar
# linha a linha e o que produz arquivos com a mesma chave duas vezes.
if [[ -n "$BASE" || -n "$AMPLISEQ" ]]; then
    : > "$CONF"
    [[ -n "$BASE" ]]     && printf 'RIEUX_PIPELINE_BASE=%s\n' "$BASE"     >> "$CONF"
    [[ -n "$AMPLISEQ" ]] && printf 'AMPLISEQ_HOME=%s\n'       "$AMPLISEQ" >> "$CONF"
    echo "Wrote $CONF"
fi
if [[ -n "$BASE" && ! -d "$BASE" ]]; then
    echo "WARNING: $BASE does not exist yet."
fi
if [[ -n "$AMPLISEQ" ]]; then
    echo "Local ampliseq pinned: $AMPLISEQ"
    [[ -d "$AMPLISEQ" ]] || echo "WARNING: $AMPLISEQ does not exist yet."
else
    echo ""
    echo "NOTE: no local ampliseq pinned (--ampliseq DIR). Without one, Nextflow"
    echo "      pulls whatever is current on GitHub, which may require a newer"
    echo "      Nextflow than your cluster module has."
fi

# ----------------------------------------------------------- prerequisites
echo ""
echo "Prerequisites:"
missing=0
check() {   # check <command> <req|opt> <note>
    if command -v "$1" >/dev/null 2>&1; then
        printf '  %-12s %s\n' "$1" "ok"
    else
        printf '  %-12s %s\n' "$1" "MISSING — $3"
        [[ "$2" == "req" ]] && missing=1
    fi
}
check python3    req "needed by every helper in bin/"
check nextflow   req "load your cluster module, or install it"
check sbatch     opt "only needed to submit to SLURM"
if command -v singularity >/dev/null 2>&1 || command -v apptainer >/dev/null 2>&1; then
    printf '  %-12s %s\n' "singularity" "ok"
else
    printf '  %-12s %s\n' "singularity" "MISSING — or apptainer; needed for the containers"
fi
if [[ -n "$BASE" && -r "$BASE/bancos.env" ]]; then
    printf '  %-12s %s\n' "bancos.env" "ok"
else
    printf '  %-12s %s\n' "bancos.env" "MISSING — run bin/prepare_unite.py, or set the DB_* variables"
fi

echo ""
(( missing )) && echo "Installed, but something required is missing (see above)." \
              || echo "Installed."
echo "Reload your shell with:  source ${RC:-your rc file}"
echo "Then run:                $COMMAND_NAME --help"

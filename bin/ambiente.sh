# ambiente.sh — carregue antes de qualquer chamada ao Nextflow.
#
#     source ambiente.sh
#
# Estas variaveis precisam existir no ambiente do PROCESSO que conduz o
# pipeline (o "driver"), nao dentro do nextflow.config: o driver le a
# config, entao o que estiver la dentro chega tarde demais para governar
# como a config e lida.
#
# ONDE RODAR O DRIVER: no NO DE LOGIN, dentro de um `screen`.
#   - o login tem internet, entao o Nextflow resolve plugins sem drama;
#   - o driver precisa poder chamar sbatch/squeue, o que nem sempre
#     funciona de dentro de um job;
#   - ele e um processo ocioso, so espera o SLURM — nao pesa no login.
# As TAREFAS e que rodam nos nos de computacao, e essas sim nao tem
# internet: por isso os containers e os bancos foram pre-baixados.

# Raiz da instalacao: containers, NXF_HOME e bancos de referencia. Ela e
# especifica de cada usuario, entao NAO pode estar cravada aqui. Defina uma vez:
#
#     echo 'RIEUX_PIPELINE_BASE=/caminho/para/pipeline' > ~/.rieux_ampliseq.conf
#
# ou exporte RIEUX_PIPELINE_BASE no seu ~/.bashrc.
if [[ -z "${RIEUX_PIPELINE_BASE:-}" && -r "$HOME/.rieux_ampliseq.conf" ]]; then
    # shellcheck disable=SC1090
    source "$HOME/.rieux_ampliseq.conf"
fi
BASE="${RIEUX_PIPELINE_BASE:-}"
if [[ -z "$BASE" || ! -d "$BASE" ]]; then
    echo "ERRO: RIEUX_PIPELINE_BASE nao definida (ou aponta para nada)." >&2
    echo "      Crie ~/.rieux_ampliseq.conf com uma linha:" >&2
    echo "        RIEUX_PIPELINE_BASE=/caminho/para/pipeline" >&2
    echo "      Esse diretorio guarda nextflow_home/, singularity/ e bancos.env." >&2
    return 1 2>/dev/null || exit 1
fi

export NXF_HOME="$BASE/nextflow_home"
export NXF_SINGULARITY_CACHEDIR="$BASE/singularity"

# Parser antigo: o Nextflow 26.04 tornou o estrito padrao e ele rejeita a
# interpolacao de `manifest` e `validation` no nextflow.config dos
# pipelines nf-core da geracao 24.x.
export NXF_SYNTAX_PARSER=v1

# Os nos de computacao nao tem o terminfo do xterm-256color, e o log enche
# de "tput: unknown terminal". Mas NAO da para forcar TERM=dumb sempre: um
# terminal dumb nao tem capacidade de limpar a tela, e o `screen` se recusa
# a iniciar com "Clear screen capability required".
# Entao so trocamos quando o TERM atual nao existe na base de terminfo.
if ! infocmp "${TERM:-dumb}" >/dev/null 2>&1; then
    export TERM=dumb
fi

# NXF_OFFLINE fica DESLIGADO no driver, porque o login tem internet.
# Ligue apenas se for conduzir o pipeline de dentro de um job:
#   export NXF_OFFLINE=true
# Mesmo com internet, a versao do plugin segue fixada no rieux.config —
# isso e reprodutibilidade, nao contorno: sem fixar, uma execucao daqui a
# seis meses pode pegar um nf-schema diferente e se comportar diferente.

# Uma unica versao de Nextflow para todo mundo. O conda instalou uma
# 26.04.6 no env_nfcore que pode sombrear o modulo; o modulo e o que a
# instituicao mantem e o que existe nos nos de computacao.
module load nextflow 2>/dev/null || true

# O modulo de Nextflow deste cluster aponta o curl para um cacert.pem que
# nao existe (caminho com nome de env errado), e toda chamada HTTPS do
# driver reclama "error adding trust anchors". Como o driver precisa de
# internet para resolver plugins, vale apontar para um bundle valido.
if ! curl -sS --max-time 5 https://api.github.com -o /dev/null 2>/dev/null; then
    for ca in /etc/pki/tls/certs/ca-bundle.crt \
              /etc/ssl/certs/ca-certificates.crt \
              /etc/pki/ca-trust/extracted/pem/tls-ca-bundle.pem; do
        if [[ -r "$ca" ]]; then
            export CURL_CA_BUNDLE="$ca"
            export SSL_CERT_FILE="$ca"
            echo "  (cacert do modulo quebrado; usando $ca)"
            break
        fi
    done
fi

# Caminhos dos bancos de referencia (DB_SILVA_GENERO, DB_SILVA_ESPECIE,
# DB_UNITE), gravados pelo instalador.
# shellcheck disable=SC1090
[[ -r "$BASE/bancos.env" ]] && source "$BASE/bancos.env"

echo "ambiente carregado:"
echo "  nextflow : $(command -v nextflow || echo AUSENTE) ($(nextflow -v 2>/dev/null | grep -oE '[0-9]+\.[0-9]+\.[0-9]+' | head -1))"
echo "  NXF_HOME : $NXF_HOME"
echo "  plugins  : $(ls "$NXF_HOME/plugins" 2>/dev/null | tr '\n' ' ')"
echo "  offline  : ${NXF_OFFLINE:-nao (driver no login)}"

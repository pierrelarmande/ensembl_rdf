#!/bin/bash
#SBATCH --job-name=ensembl_rdf
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#
# Convert one species per array task.
#
#   mkdir -p logs
#   sbatch --array=1-$(python3 bin/species_config.py species ensembl_rdf/config/species_agrold.yaml | wc -l) \
#          bin/convert_slurm.sh ensembl_rdf/config/species_agrold.yaml /path/to/workdir
#
# Each task downloads the core tables of its species, converts them and gzips
# the Turtle files, all inside WORKDIR/<species>_core_*. Tasks are independent,
# so a failed one can be resubmitted alone with --array=<n>.
#
# Memory: the converter loads every table it needs in RAM, about 35-40x the
# size of the compressed dumps (Arabidopsis: 44 MB of dumps -> 1.6 GB). 16 GB
# covers every plant species, bread wheat included; lower it to 8G if your
# cluster is tight, or raise it for a genome larger than wheat.
set -euo pipefail

CONFIG=${1:?usage: $0 CONFIG_YAML [WORKDIR]}
WORKDIR=${2:-$PWD}
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
TASK=${SLURM_ARRAY_TASK_ID:-1}

# The task changes directory, so make the paths absolute first
[ -f "$CONFIG" ] || { echo "Error: no such config file: $CONFIG" >&2; exit 1; }
CONFIG=$(cd "$(dirname "$CONFIG")" && pwd)/$(basename "$CONFIG")
mkdir -p "$WORKDIR"
WORKDIR=$(cd "$WORKDIR" && pwd)

# Cluster environment. MODULES names the environment modules to load, so a
# different site only has to override it:
#   MODULES="raptor2 python/3.11" sbatch ... bin/convert_slurm.sh ...
# Set it empty to load nothing.
MODULES=${MODULES-"bioinfo-trop raptor2/2.0.16"}
if [ -n "$MODULES" ]; then
    # `module` is a shell function, not a binary, and a batch shell may not
    # have sourced it yet.
    if ! type module >/dev/null 2>&1 && [ -n "${MODULESHOME:-}" ] \
       && [ -f "$MODULESHOME/init/bash" ]; then
        . "$MODULESHOME/init/bash"
    fi
    if type module >/dev/null 2>&1; then
        for m in $MODULES; do
            module load "$m" || echo "Warning: could not load module $m" >&2
        done
    else
        echo "Warning: no module command; skipping $MODULES" >&2
    fi
fi
# A conda environment or an extracted package works just as well:
#   export PATH="$HOME/envs/rapper/bin:$PATH"

command -v rapper >/dev/null || echo "Warning: rapper not in PATH; Turtle will not be validated" >&2
python3 -c "import yaml" 2>/dev/null || {
    echo "Error: PyYAML not installed (pip install -r requirements.txt)" >&2; exit 1; }

SPECIES=$(python3 "$SCRIPT_DIR/species_config.py" species "$CONFIG" | sed -n "${TASK}p")
if [ -z "$SPECIES" ]; then
    echo "Error: no species at index $TASK in $CONFIG" >&2
    exit 1
fi

echo "[$(date)] task $TASK: $SPECIES"
cd "$WORKDIR"

# Download this species only; already downloaded tables are overwritten, so a
# resubmitted task restarts cleanly.
python3 "$SCRIPT_DIR/download_files.py" -f "$CONFIG" -s "$SPECIES"

# Convert, then rapper + gzip. Model, base URI and entity selection come from
# the config file.
bash "$SCRIPT_DIR/convert.sh" -f "$CONFIG" -s "$SPECIES"

echo "[$(date)] task $TASK: $SPECIES done"
du -sh "${SPECIES}"_core_* 2>/dev/null || true

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
#          bin/convert_slurm.sh ensembl_rdf/config/species_agrold.yaml /path/to/workdir [/path/to/project_dir]
#
# Submit from the repository root. Each task downloads the core tables of its
# species, converts them and gzips the Turtle files, all inside
# WORKDIR/<species>_core_*. Tasks are independent, so a failed one can be
# resubmitted alone with --array=<n>.
#
# Give a third argument to also copy the results — *.ttl.gz, conversion.json,
# triple_counts.tsv, xref_report.tsv — from WORKDIR into
# PROJECT_DIR/<species>_core_*/ once the species converts successfully. The
# downloaded MySQL dumps are not copied, since they are large and
# re-downloadable; set COPY_RAW_TABLES=1 to copy them too.
#
# This script locates its own bin/ directory on its own (Slurm copies the
# batch script to a spool directory, so $0 cannot be trusted for that); if it
# still fails, set ENSEMBL_RDF_BIN to bin/'s absolute path. If `module load`
# fails here but works in an interactive shell, your site likely needs its
# init script sourced explicitly in a batch shell — point SITE_INIT at it.
#
# Memory: the converter loads every table it needs in RAM, about 35-40x the
# size of the compressed dumps (Arabidopsis: 44 MB of dumps -> 1.6 GB). 16 GB
# covers every plant species, bread wheat included; lower it to 8G if your
# cluster is tight, or raise it for a genome larger than wheat.
set -euo pipefail

CONFIG=${1:?usage: $0 CONFIG_YAML [WORKDIR] [PROJECT_DIR]}
WORKDIR=${2:-$PWD}
PROJECT_DIR=${3:-}
TASK=${SLURM_ARRAY_TASK_ID:-1}

# Locate this script's own directory (bin/), where species_config.py,
# download_files.py and convert.sh live. Under sbatch, Slurm copies the batch
# script into its spool directory and runs it from there, so $0 resolves to
# something like /var/spool/slurmd/jobNNNN/slurm_script, not to bin/ — dirname
# "$0" is then useless. Try, in order: an explicit override, the plain dirname
# (correct outside sbatch, e.g. when testing this script directly), what
# `scontrol` knows the job was submitted as, and the standard repo layout
# relative to the directory sbatch was run from.
find_script_dir() {
    local candidate
    if [ -n "${ENSEMBL_RDF_BIN:-}" ] && [ -f "$ENSEMBL_RDF_BIN/species_config.py" ]; then
        echo "$ENSEMBL_RDF_BIN"; return
    fi
    candidate=$(cd "$(dirname "$0")" 2>/dev/null && pwd) || candidate=""
    if [ -n "$candidate" ] && [ -f "$candidate/species_config.py" ]; then
        echo "$candidate"; return
    fi
    if [ -n "${SLURM_JOB_ID:-}" ] && command -v scontrol >/dev/null; then
        candidate=$(scontrol show job "$SLURM_JOB_ID" 2>/dev/null \
            | grep -oE 'Command=\S+' | cut -d= -f2 | xargs -r dirname 2>/dev/null)
        if [ -n "$candidate" ] && [ -f "$candidate/species_config.py" ]; then
            echo "$candidate"; return
        fi
    fi
    if [ -n "${SLURM_SUBMIT_DIR:-}" ] && [ -f "$SLURM_SUBMIT_DIR/bin/species_config.py" ]; then
        echo "$SLURM_SUBMIT_DIR/bin"; return
    fi
    return 1
}
SCRIPT_DIR=$(find_script_dir) || {
    echo "Error: cannot find this script's directory (bin/) — Slurm copies the" >&2
    echo "       batch script to a spool directory, so \$0 is not reliable here." >&2
    echo "       Set ENSEMBL_RDF_BIN to the absolute path of bin/ and resubmit," >&2
    echo "       e.g.: sbatch --export=ALL,ENSEMBL_RDF_BIN=\$PWD/bin ..." >&2
    exit 1
}

# The task changes directory, so make the paths absolute first
[ -f "$CONFIG" ] || { echo "Error: no such config file: $CONFIG" >&2; exit 1; }
CONFIG=$(cd "$(dirname "$CONFIG")" && pwd)/$(basename "$CONFIG")
mkdir -p "$WORKDIR"
WORKDIR=$(cd "$WORKDIR" && pwd)
if [ -n "$PROJECT_DIR" ]; then
    mkdir -p "$PROJECT_DIR"
    PROJECT_DIR=$(cd "$PROJECT_DIR" && pwd)
fi

# Cluster environment. MODULES names the environment modules to load, so a
# different site only has to override it:
#   MODULES="raptor2 python/3.11" sbatch ... bin/convert_slurm.sh ...
# Set it empty to load nothing.
MODULES=${MODULES-"bioinfo-itrop raptor2/2.0.16"}
if [ -n "$MODULES" ]; then
    # `module` is a shell function, not a binary: a non-login batch shell may
    # not have sourced the site's init script, in which case the command
    # either does not exist, or runs but with a MODULEPATH too narrow to find
    # a module that a login shell sees fine ("Unable to locate a modulefile").
    # Point SITE_INIT at your site's script if none of the common ones work.
    for f in "${SITE_INIT:-}" /etc/profile.d/modules.sh /etc/profile.d/lmod.sh \
             "${MODULESHOME:-}/init/bash" /usr/share/lmod/lmod/init/bash \
             /usr/share/Modules/init/bash; do
        [ -n "$f" ] && [ -f "$f" ] && . "$f" && break
    done
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

if [ -n "$PROJECT_DIR" ]; then
    for species_dir in "${SPECIES}"_core_*; do
        [ -d "$species_dir" ] || continue
        dest="$PROJECT_DIR/$species_dir"
        mkdir -p "$dest"
        if command -v rsync >/dev/null; then
            patterns=(--include='*.ttl.gz' --include='conversion.json' \
                      --include='triple_counts.tsv' --include='xref_report.tsv')
            [ "${COPY_RAW_TABLES:-0}" = "1" ] && patterns+=(--include='*.txt.gz')
            rsync -a "${patterns[@]}" --exclude='*' "$species_dir/" "$dest/"
        else
            for pat in '*.ttl.gz' conversion.json triple_counts.tsv xref_report.tsv; do
                cp -f "$species_dir"/$pat "$dest/" 2>/dev/null || true
            done
            if [ "${COPY_RAW_TABLES:-0}" = "1" ]; then
                cp -f "$species_dir"/*.txt.gz "$dest/" 2>/dev/null || true
            fi
        fi
        echo "[$(date)] task $TASK: copied results to $dest"
    done
fi

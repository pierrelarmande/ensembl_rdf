#!/bin/bash
#SBATCH --job-name=ensembl_rdf_rapper
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#
# Validate and normalize with rapper the Turtle a prior run already produced
# but could not pass through it (e.g. rapper's module was not found) — one
# species per array task, without re-downloading or re-converting anything.
#
#   mkdir -p logs
#   sbatch --array=1-$(python3 bin/species_config.py species ensembl_rdf/config/species_agrold.yaml | wc -l) \
#          bin/rapper_postprocess.sh ensembl_rdf/config/species_agrold.yaml /path/to/workdir
#
# Submit from the repository root. Afterwards, `python3 bin/run.py CONFIG -o
# WORKDIR` rebuilds manifest.json with the triple counts this fills in —
# every species is already "done" (conversion.json and *.ttl.gz both exist),
# so that call only reads what is on disk, touching nothing else.
#
# Same bin/ resolution and module-loading caveats as convert_slurm.sh
# (duplicated here too, deliberately — see that script's header for
# ENSEMBL_RDF_BIN, SITE_INIT, MODULES).
set -euo pipefail

CONFIG=${1:?usage: $0 CONFIG_YAML WORKDIR}
WORKDIR=${2:?usage: $0 CONFIG_YAML WORKDIR}
TASK=${SLURM_ARRAY_TASK_ID:-1}

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

[ -f "$CONFIG" ] || { echo "Error: no such config file: $CONFIG" >&2; exit 1; }
CONFIG=$(cd "$(dirname "$CONFIG")" && pwd)/$(basename "$CONFIG")
[ -d "$WORKDIR" ] || { echo "Error: no such directory: $WORKDIR" >&2; exit 1; }
WORKDIR=$(cd "$WORKDIR" && pwd)

MODULES=${MODULES-"bioinfo-itrop raptor2/2.0.16"}
if [ -n "$MODULES" ]; then
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
command -v rapper >/dev/null || {
    echo "Error: rapper still not found in PATH after loading modules ($MODULES)." >&2
    echo "       Set SITE_INIT to your site's module init script, or MODULES to" >&2
    echo "       the right module name(s), and resubmit." >&2
    exit 1
}

SPECIES=$(python3 "$SCRIPT_DIR/species_config.py" species "$CONFIG" | sed -n "${TASK}p")
[ -n "$SPECIES" ] || { echo "Error: no species at index $TASK in $CONFIG" >&2; exit 1; }

cd "$WORKDIR"
shopt -s nullglob
species_dirs=("${SPECIES}"_core_*)
shopt -u nullglob
if [ "${#species_dirs[@]}" -eq 0 ]; then
    echo "Error: no ${SPECIES}_core_* directory in $WORKDIR" >&2
    exit 1
fi

# Same split-at-large-file + validate + record-triple-count logic as
# convert.sh's process_turtle_file/record_triples, duplicated (not sourced)
# for the same reason as the rest of this script's plumbing.
SPLIT_THRESHOLD=20000000

record_triples() {
    local file=$1 stderr=$2
    local n
    n=$(grep -oE 'Parsing returned [0-9]+ triples' "$stderr" | grep -oE '[0-9]+' \
        | awk '{s+=$1} END{print s+0}')
    printf '%s\t%s\n' "$file" "$n" >> triple_counts.tsv
}

process_turtle_file() {
    local file=$1
    echo "Processing $file..." >&2
    local line_count=$(wc -l < "$file")
    echo "$file contains $line_count lines" >&2

    if [ "$line_count" -gt "$SPLIT_THRESHOLD" ]; then
        echo "File is large, splitting into chunks of $SPLIT_THRESHOLD lines" >&2
        local tmp_dir="tmp_split_$(basename "$file" .ttl)"
        mkdir -p "$tmp_dir"
        awk -v threshold="$SPLIT_THRESHOLD" -v outdir="$tmp_dir" '
        BEGIN { chunk_num = 0; line_count = 0; current_file = sprintf("%s/chunk_%02d", outdir, chunk_num) }
        { print $0 > current_file; line_count++
          if ($0 ~ /\.$/ && line_count >= threshold) {
              close(current_file); chunk_num++
              current_file = sprintf("%s/chunk_%02d", outdir, chunk_num); line_count = 0 } }
        END { close(current_file) }
        ' "$file"
        > "${file}.processed"
        grep -E "^@prefix|^@base" "$file" > "${tmp_dir}/prefixes.ttl"
        for chunk in "$tmp_dir"/chunk_*; do
            echo "Processing chunk $chunk" >&2
            cat "${tmp_dir}/prefixes.ttl" "$chunk" > "${chunk}.with_prefix"
            rapper -i turtle -o turtle "${chunk}.with_prefix" 2> "${chunk}.log" > "${chunk}.processed"
            cat "${chunk}.log" >> "${tmp_dir}/rapper.log"
            cat "${chunk}.log" >&2
            if [ "$chunk" = "$tmp_dir/chunk_00" ]; then
                cat "${chunk}.processed" >> "${file}.processed"
            else
                grep -v -E "^@prefix|^@base" "${chunk}.processed" >> "${file}.processed"
            fi
        done
        mv "${file}.processed" "$file"
        record_triples "$file" "${tmp_dir}/rapper.log"
        rm -rf "$tmp_dir"
    else
        rapper -i turtle -o turtle "$file" 2> "${file}.rapper.log" > "${file}.rapper.ttl"
        cat "${file}.rapper.log" >&2
        mv "${file}.rapper.ttl" "$file"
        record_triples "$file" "${file}.rapper.log"
        rm -f "${file}.rapper.log"
    fi
}

for species_dir in "${species_dirs[@]}"; do
    echo "[$(date)] $species_dir"
    cd "$WORKDIR/$species_dir"
    : > triple_counts.tsv
    shopt -s nullglob
    gz_files=(*.ttl.gz)
    shopt -u nullglob
    for gz in "${gz_files[@]}"; do
        f="${gz%.gz}"
        gunzip -f "$gz"
        process_turtle_file "$f"
        gzip -f "$f"
    done
    cd "$WORKDIR"
done
echo "[$(date)] task $TASK: $SPECIES done"

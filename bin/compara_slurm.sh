#!/bin/bash
#SBATCH --job-name=ensembl_rdf_compara
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --time=2-00:00:00
#
# Download an Ensembl Compara dump and emit pairwise orthology/paralogy
# relations (bin/compara_orthology.py) for the species listed in a species
# configuration. A single task, not an array: Compara is one dump for the
# whole division, not per species.
#
#   mkdir -p logs
#   sbatch bin/compara_slurm.sh ensembl_rdf/config/species_agrold.yaml /path/to/workdir [/path/to/project_dir]
#
# Submit from the repository root. Give a third argument to also copy
# orthology.ttl.gz and triple_counts.tsv into PROJECT_DIR/<compara_dir>/ once
# done (see bin/convert_slurm.sh's header — same idea, one file instead of a
# whole species directory).
#
# This script resolves its own bin/ directory and loads cluster modules the
# same way bin/convert_slurm.sh does (duplicated rather than shared between
# the two, deliberately: a change to one must not risk the other while a
# cluster run is in flight). See that script's header for ENSEMBL_RDF_BIN,
# SITE_INIT and MODULES.
#
# Sizing: homology.txt.gz + homology_member.txt.gz alone were 65 GB
# compressed for the whole Ensembl Plants 63 division (every species pair
# Ensembl computed, not just the ones asked for) — download time dominates.
# A resubmission skips files already downloaded (download_compara.py only
# writes a file under its final name once a fetch fully succeeds, so a file
# present there is already complete; --force redownloads anyway). Memory is
# unmeasured for a real run: it holds, for species of interest, every
# relevant homology_id and gene_member_id-to-stable_id mapping, which could
# be a few GB for 51 species or well more; 64G is a starting point —
# compara_orthology.py logs row counts and memory as the small tables load,
# before the run reaches the two large files, to judge a first real run.
set -euo pipefail

CONFIG=${1:?usage: $0 CONFIG_YAML [WORKDIR] [PROJECT_DIR]}
WORKDIR=${2:-$PWD}
PROJECT_DIR=${3:-}

# Locate this script's own directory (bin/) — see convert_slurm.sh's header:
# Slurm copies the batch script to a spool directory and runs it from there,
# so $0 cannot be trusted under sbatch.
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
mkdir -p "$WORKDIR"
WORKDIR=$(cd "$WORKDIR" && pwd)
if [ -n "$PROJECT_DIR" ]; then
    mkdir -p "$PROJECT_DIR"
    PROJECT_DIR=$(cd "$PROJECT_DIR" && pwd)
fi

# Cluster environment — same as convert_slurm.sh.
MODULES=${MODULES-"bioinfo-trop raptor2/2.0.16"}
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

HAVE_RAPPER=1
command -v rapper >/dev/null || {
    HAVE_RAPPER=0
    echo "Warning: rapper not in PATH; orthology.ttl will not be validated or normalized" >&2
}
python3 -c "import yaml" 2>/dev/null || {
    echo "Error: PyYAML not installed (pip install -r requirements.txt)" >&2; exit 1; }

URL=$(python3 "$SCRIPT_DIR/species_config.py" url "$CONFIG")
[ -n "$URL" ] || { echo "Error: $CONFIG has no url" >&2; exit 1; }

cd "$WORKDIR"
echo "[$(date)] downloading the Compara dump from $URL"
python3 "$SCRIPT_DIR/download_compara.py" "$URL"

COMPARA_DIR=$(ls -d *_compara_* 2>/dev/null | head -1) || true
[ -n "$COMPARA_DIR" ] || { echo "Error: no *_compara_* directory after download" >&2; exit 1; }

echo "[$(date)] generating orthology/paralogy relations from $COMPARA_DIR"
python3 "$SCRIPT_DIR/compara_orthology.py" "$COMPARA_DIR" -f "$CONFIG"

ORTHOLOGY_FILE="$COMPARA_DIR/orthology.ttl"
: > "$COMPARA_DIR/triple_counts.tsv"
if [ "$HAVE_RAPPER" -eq 1 ]; then
    # Same split-at-large-file + validate + record-triple-count logic as
    # convert.sh's process_turtle_file/record_triples, duplicated (not
    # sourced) for the same reason as the rest of this script's plumbing.
    SPLIT_THRESHOLD=20000000
    line_count=$(wc -l < "$ORTHOLOGY_FILE")
    echo "$ORTHOLOGY_FILE contains $line_count lines" >&2
    if [ "$line_count" -gt "$SPLIT_THRESHOLD" ]; then
        tmp_dir="tmp_split_orthology"
        mkdir -p "$tmp_dir"
        awk -v threshold="$SPLIT_THRESHOLD" -v outdir="$tmp_dir" '
        BEGIN { chunk_num = 0; line_count = 0; current_file = sprintf("%s/chunk_%02d", outdir, chunk_num) }
        { print $0 > current_file; line_count++
          if ($0 ~ /\.$/ && line_count >= threshold) {
              close(current_file); chunk_num++
              current_file = sprintf("%s/chunk_%02d", outdir, chunk_num); line_count = 0 } }
        END { close(current_file) }
        ' "$ORTHOLOGY_FILE"
        > "${ORTHOLOGY_FILE}.processed"
        grep -E "^@prefix|^@base" "$ORTHOLOGY_FILE" > "${tmp_dir}/prefixes.ttl"
        for chunk in "$tmp_dir"/chunk_*; do
            echo "Processing chunk $chunk" >&2
            cat "${tmp_dir}/prefixes.ttl" "$chunk" > "${chunk}.with_prefix"
            rapper -i turtle -o turtle "${chunk}.with_prefix" 2> "${chunk}.log" > "${chunk}.processed"
            cat "${chunk}.log" >> "${tmp_dir}/rapper.log"
            cat "${chunk}.log" >&2
            if [ "$chunk" = "$tmp_dir/chunk_00" ]; then
                cat "${chunk}.processed" >> "${ORTHOLOGY_FILE}.processed"
            else
                grep -v -E "^@prefix|^@base" "${chunk}.processed" >> "${ORTHOLOGY_FILE}.processed"
            fi
        done
        mv "${ORTHOLOGY_FILE}.processed" "$ORTHOLOGY_FILE"
        n=$(grep -oE 'Parsing returned [0-9]+ triples' "${tmp_dir}/rapper.log" | grep -oE '[0-9]+' | awk '{s+=$1} END{print s+0}')
        printf 'orthology.ttl\t%s\n' "$n" >> "$COMPARA_DIR/triple_counts.tsv"
        rm -rf "$tmp_dir"
    else
        rapper -i turtle -o turtle "$ORTHOLOGY_FILE" 2> "${ORTHOLOGY_FILE}.rapper.log" > "${ORTHOLOGY_FILE}.rapper.ttl"
        cat "${ORTHOLOGY_FILE}.rapper.log" >&2
        mv "${ORTHOLOGY_FILE}.rapper.ttl" "$ORTHOLOGY_FILE"
        n=$(grep -oE 'Parsing returned [0-9]+ triples' "${ORTHOLOGY_FILE}.rapper.log" | grep -oE '[0-9]+' | awk '{s+=$1} END{print s+0}')
        printf 'orthology.ttl\t%s\n' "$n" >> "$COMPARA_DIR/triple_counts.tsv"
        rm -f "${ORTHOLOGY_FILE}.rapper.log"
    fi
fi
gzip -f "$ORTHOLOGY_FILE"

echo "[$(date)] done"
du -sh "$COMPARA_DIR/orthology.ttl.gz" 2>/dev/null || true

if [ -n "$PROJECT_DIR" ]; then
    dest="$PROJECT_DIR/$COMPARA_DIR"
    mkdir -p "$dest"
    if command -v rsync >/dev/null; then
        rsync -a --include='orthology.ttl.gz' --include='triple_counts.tsv' --exclude='*' \
            "$COMPARA_DIR/" "$dest/"
    else
        cp -f "$COMPARA_DIR/orthology.ttl.gz" "$dest/" 2>/dev/null || true
        cp -f "$COMPARA_DIR/triple_counts.tsv" "$dest/" 2>/dev/null || true
    fi
    echo "[$(date)] copied results to $dest"
fi

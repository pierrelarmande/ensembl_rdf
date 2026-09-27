#!/usr/bin/bash
# set -euo pipefail
THIS_DIR=$(cd $(dirname $0); pwd)
SCRIPT_DIR=$THIS_DIR/../ensembl_rdf
CONFIG_DIR=$SCRIPT_DIR/config

# 分割処理の最大行数
SPLIT_THRESHOLD=20000000

SPECIES=()
SPECIES_FILE=""
TYPES=()
EXCLUDE=()
BASE_URI=""
TERMS_URI=""
MODEL=""

usage() {
    cat >&2 <<EOT
Usage: $0 [-s SPECIES ...] [-f CONFIG_YAML] [-t ENTITY ...] [-x ENTITY ...]
          [-b BASE_URI] [-u TERMS_URI] [-m MODEL] [dir ...]
Convert Ensembl core MySQL dumps to RDF. Species are resolved to the
<species>_core_* directories of the current directory.
  -s SPECIES  production name (e.g. arabidopsis_thaliana); may be repeated
  -f FILE     YAML file with \`species\`, \`entities\`, \`exclude\`, \`model\`,
              \`base_uri\` (see ensembl_rdf/config/species.yaml)
  -t ENTITY   only output this entity type; may be repeated (default: all of
              gene transcript translation exon exon_transcript xref chromosome)
  -x ENTITY   do not output this entity type; may be repeated
  -b URI      base of the resource URIs (default: http://rdf.ebi.ac.uk, e.g.
              -b http://purl.agrold.org gives http://purl.agrold.org/resource/...)
  -u URI      namespace of the model vocabulary, overriding the profile
  -m MODEL    vocabulary profile: a name in ensembl_rdf/config/models/
              (ensembl, agrold) or a YAML file (default: ensembl)
  dir         core database directory (as downloaded by download_files.py)
EOT
    exit 1
}

while getopts "s:f:t:x:b:u:m:h" opt; do
  case "$opt" in
    s) SPECIES+=("$OPTARG") ;;
    f) SPECIES_FILE=$OPTARG ;;
    t) TYPES+=("$OPTARG") ;;
    x) EXCLUDE+=("$OPTARG") ;;
    b) BASE_URI=$OPTARG ;;
    u) TERMS_URI=$OPTARG ;;
    m) MODEL=$OPTARG ;;
    *) usage ;;
  esac
done
shift $((OPTIND-1))

DIRS=("$@")
if [ -n "$SPECIES_FILE" ]; then
    # -s (or directories) selects the species; the file then only supplies options
    if [ "${#SPECIES[@]}" -eq 0 ] && [ "${#DIRS[@]}" -eq 0 ]; then
        while IFS= read -r d; do
            DIRS+=("$d")
        done < <(python3 "$THIS_DIR/species_config.py" dirs "$SPECIES_FILE" *_core_*)
    fi
    for e in $(python3 "$THIS_DIR/species_config.py" entities "$SPECIES_FILE"); do TYPES+=("$e"); done
    for e in $(python3 "$THIS_DIR/species_config.py" exclude "$SPECIES_FILE"); do EXCLUDE+=("$e"); done
    [ -z "$BASE_URI" ] && BASE_URI=$(python3 "$THIS_DIR/species_config.py" base_uri "$SPECIES_FILE")
    [ -z "$TERMS_URI" ] && TERMS_URI=$(python3 "$THIS_DIR/species_config.py" terms_uri "$SPECIES_FILE")
    [ -z "$MODEL" ] && MODEL=$(python3 "$THIS_DIR/species_config.py" model "$SPECIES_FILE")
fi

for sp in "${SPECIES[@]+"${SPECIES[@]}"}"; do
    matched=$(ls -d ${sp}_core_* 2>/dev/null || true)
    if [ -z "$matched" ]; then
        echo "Warning: no directory matches ${sp}_core_*, skipping." >&2
        continue
    fi
    for d in $matched; do DIRS+=("$d"); done
done

if [ "${#DIRS[@]}" -eq 0 ]; then
    usage
fi

# Options for the converter, and the resolved entity list for the rapper loop
CONV_OPTS=()
[ "${#TYPES[@]}" -gt 0 ] && CONV_OPTS+=(-t "${TYPES[@]}")
[ "${#EXCLUDE[@]}" -gt 0 ] && CONV_OPTS+=(-x "${EXCLUDE[@]}")
[ -n "$BASE_URI" ] && CONV_OPTS+=(-b "$BASE_URI")
[ -n "$TERMS_URI" ] && CONV_OPTS+=(--terms-uri "$TERMS_URI")
[ -n "$MODEL" ] && CONV_OPTS+=(-m "$MODEL")
SELECTED=$(python3 "$SCRIPT_DIR/rdf_converter_ensembl_db.py" --list-targets \
    "$CONFIG_DIR/dbinfo.json" . "${CONV_OPTS[@]+"${CONV_OPTS[@]}"}")
echo "Entities: $SELECTED" >&2

# Record how many triples rapper parsed, for the run manifest (bin/run.py)
record_triples() {
    local file=$1 stderr=$2
    local n
    n=$(grep -oE 'Parsing returned [0-9]+ triples' "$stderr" | grep -oE '[0-9]+' \
        | awk '{s+=$1} END{print s+0}')
    printf '%s\t%s\n' "$file" "$n" >> triple_counts.tsv
}

# Turtle ファイルを分割して rapper で処理する関数
process_turtle_file() {
    local file=$1
    echo "Processing $file..." >&2

    # ファイルの行数を取得
    local line_count=$(wc -l < "$file")
    echo "$file contains $line_count lines" >&2

    if [ "$line_count" -gt "$SPLIT_THRESHOLD" ]; then
        echo "File is large, splitting into chunks of $SPLIT_THRESHOLD lines" >&2

        # 一時ディレクトリを作成
        local tmp_dir="tmp_split_$(basename "$file" .ttl)"
        mkdir -p "$tmp_dir"

        # ファイルを分割（トリプルの境界で）
        awk -v threshold="$SPLIT_THRESHOLD" -v outdir="$tmp_dir" '
        BEGIN {
            chunk_num = 0
            line_count = 0
            current_file = sprintf("%s/chunk_%02d", outdir, chunk_num)
        }
        {
            print $0 > current_file
            line_count++
            # `.` で終わる行かつ閾値を超えた場合に次のファイルへ
            if ($0 ~ /\.$/ && line_count >= threshold) {
                close(current_file)
                chunk_num++
                current_file = sprintf("%s/chunk_%02d", outdir, chunk_num)
                line_count = 0
            }
        }
        END {
            close(current_file)
        }
        ' "$file"


        # 分割されたファイルを処理して結合
        > "${file}.processed"

        # プレフィックス部分を保存（通常ファイルの先頭部分）
        grep -E "^@prefix|^@base" "$file" > "${tmp_dir}/prefixes.ttl"

        for chunk in "$tmp_dir"/chunk_*; do
            echo "Processing chunk $chunk" >&2

            # プレフィックスをチャンクの先頭に追加して rapper で処理
            cat "${tmp_dir}/prefixes.ttl" "$chunk" > "${chunk}.with_prefix"
            rapper -i turtle -o turtle "${chunk}.with_prefix" 2> "${chunk}.log" > "${chunk}.processed"
            cat "${chunk}.log" >> "${tmp_dir}/rapper.log"
            cat "${chunk}.log" >&2

            # プレフィックス部分を除去して結合（最初のチャンクを除く）
            if [ "$chunk" = "$tmp_dir/chunk_00" ]; then
                cat "${chunk}.processed" >> "${file}.processed"
            else
                grep -v -E "^@prefix|^@base" "${chunk}.processed" >> "${file}.processed"
            fi
        done

        # 元のファイルを置き換え
        mv "${file}.processed" "$file"
        record_triples "$file" "${tmp_dir}/rapper.log"

        # 一時ディレクトリを削除
        rm -rf "$tmp_dir"
    else
        # サイズが閾値以下なら通常処理
        rapper -i turtle -o turtle "$file" 2> "${file}.rapper.log" > "${file}.rapper.ttl"
        cat "${file}.rapper.log" >&2
        mv "${file}.rapper.ttl" "$file"
        record_triples "$file" "${file}.rapper.log"
        rm -f "${file}.rapper.log"
    fi
}

for d in "${DIRS[@]}"; do
    (
        set -euo pipefail
        if [ ! -d "$d" ]; then
            echo "Warning: '$d' is not a directory or does not exist, skipping." >&2
            exit 1
        fi

        echo "$d" >&2
        cd "$d"
        : > triple_counts.tsv

        python3 "$SCRIPT_DIR/rdf_converter_ensembl_db.py" "$CONFIG_DIR/dbinfo.json" . \
            "${CONV_OPTS[@]+"${CONV_OPTS[@]}"}"

        for f in $SELECTED; do
            if [ -f "$f.ttl" ]; then
                process_turtle_file "$f.ttl"
                gzip -f "$f.ttl"
            else
                echo "Warning: $f.ttl not found" >&2
            fi
        done
    )
done

#!/usr/bin/bash
set -euo pipefail
SCRIPT_DIR=$(cd $(dirname $0); pwd)
CONFIG_DIR=$SCRIPT_DIR/../config

# 分割処理の最大行数
SPLIT_THRESHOLD=20000000

usage() {
    cat >&2 <<EOT
Usage: $0 [-s SPECIES ...] [-f CONFIG_YAML] [-e ENTITY ...] [-x ENTITY ...] [-b BASE_URI] [-m MODEL] [dir ...]
Convert Ensembl core MySQL dumps to RDF. Species are resolved to the
<species>_core_* directories of the current directory.
  -s SPECIES  production name (e.g. arabidopsis_thaliana); may be repeated
  -f FILE     YAML file with \`species\`, \`entities\`, \`exclude\` lists (see config/species.yaml)
  -e ENTITY   only output this entity type; may be repeated (default: all of
              gene transcript translation exon exon_transcript xref chromosome)
  -x ENTITY   do not output this entity type; may be repeated
  -b URI      base of the resource URIs (default: http://rdf.ebi.ac.uk, e.g.
              -b http://purl.agrold.org gives http://purl.agrold.org/resource/ensembl/...)
  -t URI      namespace of the terms: vocabulary (default: http://rdf.ebi.ac.uk/terms/ensembl/)
  -m MODEL    vocabulary profile: a name in config/models/ (ensembl, agrold) or
              a YAML file (default: ensembl)
  dir         core database directory (as downloaded by download_files.py)
EOT
    exit 1
}

species=()
species_file=""
entities=()
exclude=()
base_uri=""
terms_uri=""
model=""
while getopts "s:f:e:x:b:t:m:h" opt; do
    case $opt in
        s) species+=("$OPTARG") ;;
        f) species_file=$OPTARG ;;
        e) entities+=("$OPTARG") ;;
        x) exclude+=("$OPTARG") ;;
        b) base_uri=$OPTARG ;;
        t) terms_uri=$OPTARG ;;
        m) model=$OPTARG ;;
        *) usage ;;
    esac
done
shift $((OPTIND - 1))

dirs=("$@")
if [ -n "$species_file" ]; then
    # -s (or directories) selects the species; the file then only supplies options
    if [ "${#species[@]}" -eq 0 ] && [ "${#dirs[@]}" -eq 0 ]; then
        while IFS= read -r d; do
            dirs+=("$d")
        done < <(python3 "$SCRIPT_DIR/species_config.py" dirs "$species_file" *_core_*)
    fi
    for e in $(python3 "$SCRIPT_DIR/species_config.py" entities "$species_file"); do entities+=("$e"); done
    for e in $(python3 "$SCRIPT_DIR/species_config.py" exclude "$species_file"); do exclude+=("$e"); done
    # command line options take precedence over the YAML file
    [ -z "$base_uri" ] && base_uri=$(python3 "$SCRIPT_DIR/species_config.py" base_uri "$species_file")
    [ -z "$terms_uri" ] && terms_uri=$(python3 "$SCRIPT_DIR/species_config.py" terms_uri "$species_file")
    [ -z "$model" ] && model=$(python3 "$SCRIPT_DIR/species_config.py" model "$species_file")
fi

# Converter options: entity selection and URIs
conv_opts=()
[ -n "$base_uri" ] && conv_opts+=(-b "$base_uri")
[ -n "$terms_uri" ] && conv_opts+=(--terms-uri "$terms_uri")
[ -n "$model" ] && conv_opts+=(-m "$model")

# Entity selection, used for the rapper/gzip loop too
entity_opts=()
if [ "${#entities[@]}" -gt 0 ]; then entity_opts+=(-e "${entities[@]}"); fi
if [ "${#exclude[@]}" -gt 0 ]; then entity_opts+=(-x "${exclude[@]}"); fi
selected=$(python3 "$SCRIPT_DIR/rdf_converter_ensembl_db.py" --list-entities "$CONFIG_DIR/dbinfo.json" "${entity_opts[@]+"${entity_opts[@]}"}")
echo "Entities: $selected" >&2
conv_opts+=("${entity_opts[@]+"${entity_opts[@]}"}")
for sp in "${species[@]+"${species[@]}"}"; do
    matched=$(ls -d ${sp}_core_* 2>/dev/null || true)
    if [ -z "$matched" ]; then
        echo "Warning: no directory matches ${sp}_core_*, skipping." >&2
        continue
    fi
    for d in $matched; do dirs+=("$d"); done
done

if [ "${#dirs[@]}" -eq 0 ]; then
    usage
fi

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

        # ファイルを分割
        split -l "$SPLIT_THRESHOLD" "$file" "$tmp_dir/chunk_"

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
            if [ "$chunk" = "$tmp_dir/chunk_aa" ]; then
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

for d in "${dirs[@]}"; do
    if [ ! -d "$d" ]; then
        echo "Warning: '$d' is not a directory or does not exist, skipping." >&2
        continue
    fi

    echo "$d" >&2
    cd "$d"
    : > triple_counts.tsv
    python3 "$SCRIPT_DIR/rdf_converter_ensembl_db.py" "$CONFIG_DIR/dbinfo.json" "${conv_opts[@]+"${conv_opts[@]}"}"
    #echo "Validating turtle files..."
    for f in $selected; do
        if [ -f "$f.ttl" ]; then
            process_turtle_file "$f.ttl"
            gzip -f "$f.ttl"
        else
            echo "Warning: $f.ttl not found" >&2
        fi
    done
    cd ..
done

#!/usr/bin/bash
set -euo pipefail
SCRIPT_DIR=$(cd $(dirname $0); pwd)
CONFIG_DIR=$SCRIPT_DIR/../config

# 分割処理の最大行数
SPLIT_THRESHOLD=20000000

usage() {
    cat >&2 <<EOT
Usage: $0 [-s SPECIES ...] [-f SPECIES_YAML] [dir ...]
Convert Ensembl core MySQL dumps to RDF. Species are resolved to the
<species>_core_* directories of the current directory.
  -s SPECIES  production name (e.g. arabidopsis_thaliana); may be repeated
  -f FILE     YAML file with a \`species\` list (see config/species.yaml)
  dir         core database directory (as downloaded by download_files.py)
EOT
    exit 1
}

species=()
species_file=""
while getopts "s:f:h" opt; do
    case $opt in
        s) species+=("$OPTARG") ;;
        f) species_file=$OPTARG ;;
        *) usage ;;
    esac
done
shift $((OPTIND - 1))

dirs=("$@")
if [ -n "$species_file" ]; then
    while IFS= read -r d; do
        dirs+=("$d")
    done < <(python3 "$SCRIPT_DIR/species_config.py" "$species_file" *_core_*)
fi
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
            rapper -i turtle -o turtle "${chunk}.with_prefix" > "${chunk}.processed"

            # プレフィックス部分を除去して結合（最初のチャンクを除く）
            if [ "$chunk" = "$tmp_dir/chunk_aa" ]; then
                cat "${chunk}.processed" >> "${file}.processed"
            else
                grep -v -E "^@prefix|^@base" "${chunk}.processed" >> "${file}.processed"
            fi
        done

        # 元のファイルを置き換え
        mv "${file}.processed" "$file"

        # 一時ディレクトリを削除
        rm -rf "$tmp_dir"
    else
        # サイズが閾値以下なら通常処理
        rapper -i turtle -o turtle "$file" > "${file}.rapper.ttl"
        mv "${file}.rapper.ttl" "$file"
    fi
}

for d in "${dirs[@]}"; do
    if [ ! -d "$d" ]; then
        echo "Warning: '$d' is not a directory or does not exist, skipping." >&2
        continue
    fi

    echo "$d" >&2
    cd "$d"
    python3 $SCRIPT_DIR/rdf_converter_ensembl_db.py $CONFIG_DIR/dbinfo.json
    #echo "Validating turtle files..."
    for f in gene transcript translation exon exon_transcript xref; do
        if [ -f "$f.ttl" ]; then
            process_turtle_file "$f.ttl"
            gzip -f "$f.ttl"
        else
            echo "Warning: $f.ttl not found" >&2
        fi
    done
    cd ..
done

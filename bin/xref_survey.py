"""Survey the cross-reference sources used by a set of Ensembl core databases.

Reports, for every external database actually referenced, how many species use
it, how many xrefs it holds, an example accession and whether
config/external_db_url.tsv already gives it a URL. Use it to decide which
sources are worth mapping before converting a whole division.

It streams xref.txt.gz (and external_db.txt.gz), so it needs neither much
memory nor a full conversion; `download_files.py --dbinfo ensembl_rdf/config/xref_dbinfo.json`
fetches just those two tables.

Usage: xref_survey.py [DIR ...]          # default: every *_core_* of the cwd
       xref_survey.py --missing-only     # only the sources without a URL
       xref_survey.py --tsv              # lines ready for external_db_url.tsv
"""
import os
import re
import csv
import gzip
import glob
import sys
import argparse
from collections import defaultdict

BASE_DIR = os.path.dirname(os.path.abspath(__file__)) + "/../"
URL_TABLE = BASE_DIR + "ensembl_rdf/config/external_db_url.tsv"


def load_url_table(path):
    """external_db_id -> URL (empty string when the id is listed but unmapped)."""
    urls = {}
    with open(path) as f:
        for row in csv.reader(f, delimiter="\t"):
            if row:
                urls[row[0]] = row[1] if len(row) > 1 else ""
    return urls


def read_external_db(directory):
    """external_db_id -> db_name"""
    names = {}
    path = os.path.join(directory, "external_db.txt.gz")
    if not os.path.exists(path):
        return names
    with gzip.open(path, "rt") as f:
        for line in f:
            fields = line.rstrip("\n").split("\t")
            if len(fields) > 1:
                names[fields[0]] = fields[1]
    return names


def survey(directories):
    """Aggregate xref counts per external_db_id over the given core directories."""
    counts = defaultdict(int)
    species = defaultdict(set)
    examples = {}
    names = {}
    conflicts = defaultdict(set)

    for directory in directories:
        path = os.path.join(directory, "xref.txt.gz")
        if not os.path.exists(path):
            print(f"Warning: no xref.txt.gz in {directory}, skipping", file=sys.stderr)
            continue
        name = re.sub(r"_core_.*", "", os.path.basename(directory.rstrip("/")))
        print(f"Reading {directory}", file=sys.stderr)
        for db_id, db_name in read_external_db(directory).items():
            conflicts[db_id].add(db_name)
            names[db_id] = db_name
        with gzip.open(path, "rt") as f:
            for line in f:
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 3:
                    continue
                db_id, acc = fields[1], fields[2]
                counts[db_id] += 1
                species[db_id].add(name)
                if db_id not in examples and acc not in ("", "\\N"):
                    examples[db_id] = acc

    return counts, species, examples, names, conflicts


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("directories", nargs="*", help="core database directories (default: ./*_core_*)")
    parser.add_argument("--missing-only", action="store_true", help="only sources without a URL")
    parser.add_argument("--tsv", action="store_true",
                        help="print the ids missing from external_db_url.tsv as lines to append "
                             "(id, empty URL, empty prefix), commented with the db name and an example")
    parser.add_argument("--url-table", default=URL_TABLE, help="default: ensembl_rdf/config/external_db_url.tsv")
    args = parser.parse_args()

    directories = args.directories or sorted(glob.glob("*_core_*"))
    if not directories:
        sys.exit("Error: no core directory given and none found in the current directory")

    urls = load_url_table(args.url_table)
    counts, species, examples, names, conflicts = survey(directories)

    rows = []
    for db_id, count in counts.items():
        url = urls.get(db_id, "")
        if args.missing_only and url:
            continue
        rows.append((db_id, names.get(db_id, "?"), len(species[db_id]), count,
                     examples.get(db_id, ""), url,
                     "" if db_id in urls else "id absent from the URL table"))
    rows.sort(key=lambda r: -r[3])

    if args.tsv:
        # Only ids absent from the table: the others are already listed and must
        # be filled in place, not appended.
        new = [r for r in rows if r[6]]
        for db_id, name, n_species, count, example, url, note in new:
            print(f"# {name}: {count} xrefs in {n_species} species, e.g. {example}")
            print(f"{db_id}\t\t")
        listed = [r for r in rows if not r[6] and not r[5]]
        print(f"{len(new)} id(s) to append; {len(listed)} already listed in "
              f"{args.url_table} with no URL, fill those in place", file=sys.stderr)
    else:
        print("external_db_id\tdb_name\tspecies\txrefs\texample\turl\tnote")
        for row in rows:
            print("\t".join(str(v) for v in row))

    ambiguous = {k: v for k, v in conflicts.items() if len(v) > 1}
    if ambiguous:
        print(f"Warning: {len(ambiguous)} external_db_id have different names across species: "
              f"{dict(list(ambiguous.items())[:5])}", file=sys.stderr)


if __name__ == "__main__":
    main()

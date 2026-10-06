"""Download the tables needed for pairwise orthology from an Ensembl Compara dump.

Ensembl Compara is one dump per division (not per species), under a single
*_compara_* directory, e.g. ensembl_compara_plants_63_116. This fetches the
tables listed in a dbinfo file from it, reusing download_files.py's HTTP
source, retry logic and .part-then-rename downloads.

Two dbinfo files cover two different uses:
  ensembl_rdf/config/dbinfo_compara.json          gene tree clusters (small)
  ensembl_rdf/config/dbinfo_compara_homology.json pairwise orthology (large:
                                                   homology_member alone can
                                                   exceed 50 GB for a division)
--dbinfo defaults to the homology one, the main reason to run this at all.

Usage:
    download_compara.py https://ftp.ebi.ac.uk/pub/ensemblgenomes/plants/current/mysql/
"""
import os
import sys
import json
import argparse

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from download_files import HttpSource, retry, log

CONFIG_DIR = os.path.dirname(os.path.abspath(__file__)) + "/../ensembl_rdf/config/"


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("url", help="base URL of the mysql dump directory (https://)")
    parser.add_argument("--dbinfo", default=CONFIG_DIR + "dbinfo_compara_homology.json",
                        help="default: ensembl_rdf/config/dbinfo_compara_homology.json")
    parser.add_argument("-o", "--output", default=".",
                        help="where the *_compara_* directory goes (default: .)")
    parser.add_argument("--force", action="store_true",
                        help="redownload files that already exist (default: skip them — "
                             "a complete file only exists under its final name once the "
                             "download that produced it succeeded, since every fetch writes "
                             "to a .part file renamed only on completion)")
    args = parser.parse_args()

    with open(args.dbinfo, "r") as f:
        dbinfo = json.load(f)
    wanted = {v["filename"] for v in dbinfo.values()}

    source = HttpSource(args.url)
    candidates = [d for d in retry("listing the release", lambda: source.list())
                  if "_compara_" in d]
    if not candidates:
        sys.exit(f"Error: no *_compara_* directory under {args.url}")
    if len(candidates) > 1:
        log(f"Warning: several *_compara_* directories found, using the first: {candidates}")
    compara_dir = candidates[0]

    local_dir = os.path.join(args.output, compara_dir)
    os.makedirs(local_dir, exist_ok=True)
    files = [os.path.basename(f) for f in retry(f"listing {compara_dir}",
                                                 lambda: source.list(compara_dir))]
    for file in files:
        if file not in wanted:
            continue
        path = os.path.join(local_dir, file)
        if not args.force and os.path.exists(path):
            log(f"Skipping (already downloaded): {compara_dir}/{file}")
            continue
        log(f"Downloading: {compara_dir}/{file}")
        retry(f"downloading {path}", lambda: source.fetch(compara_dir + "/" + file, path + ".part"))
        os.replace(path + ".part", path)
    log(f"Done: {local_dir}")


if __name__ == "__main__":
    main()

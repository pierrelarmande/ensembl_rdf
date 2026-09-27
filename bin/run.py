"""Run the whole conversion for a config file: download, convert, manifest.

One command for a pipeline to call, whatever orchestrates it:

    run.py ensembl_rdf/config/species_agrold.yaml -o /path/to/workdir

Each species is independent: its core tables are downloaded, converted, checked
by rapper and gzipped, and a species that is already done is skipped, so an
interrupted run resumes where it stopped. The run writes manifest.json next to
the species directories, listing what was produced with its triple counts, for
whatever loads the result afterwards.

Exit status is 0 when every species succeeded, 1 otherwise; the manifest and
the summary name the ones that failed.
"""
import os
import re
import sys
import json
import glob
import shutil
import argparse
import datetime
import subprocess
import concurrent.futures

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from species_config import load_config

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
BASE_DIR = SCRIPT_DIR + "/../"
MANIFEST = "manifest.json"


def log(msg):
    print(f'[{datetime.datetime.now().isoformat(timespec="seconds")}] {msg}', file=sys.stderr)


def tool_version():
    """Commit this converter runs from, so a manifest says what produced it."""
    try:
        out = subprocess.run(["git", "-C", BASE_DIR, "describe", "--always", "--dirty"],
                             capture_output=True, text=True, timeout=10)
        return out.stdout.strip() or "unknown"
    except Exception:
        return "unknown"


def core_dir_of(species, workdir):
    """The <species>_core_* directory of this species, if it was downloaded."""
    matches = [d for d in glob.glob(os.path.join(workdir, species + "_core_*")) if os.path.isdir(d)]
    return sorted(matches)[-1] if matches else None


def is_done(directory):
    """A species is done when its conversion report sits next to gzipped Turtle."""
    return (directory
            and os.path.exists(os.path.join(directory, "conversion.json"))
            and glob.glob(os.path.join(directory, "*.ttl.gz")))


def run(command, cwd):
    subprocess.run(command, cwd=cwd, check=True)


def convert_species(species, config_path, workdir, skip_download):
    """Download and convert one species; returns its manifest entry."""
    entry = {"species": species, "status": "failed"}
    try:
        if not skip_download:
            log(f"{species}: downloading")
            run([sys.executable, SCRIPT_DIR + "/download_files.py",
                 "-f", config_path, "-s", species], cwd=workdir)

        directory = core_dir_of(species, workdir)
        if directory is None:
            raise RuntimeError(f"no {species}_core_* directory (download failed or wrong name?)")
        entry["core_db"] = os.path.basename(directory)

        log(f"{species}: converting")
        run(["bash", SCRIPT_DIR + "/convert.sh", "-f", config_path, "-s", species], cwd=workdir)

        entry.update(collect(directory))
        entry["status"] = "done"
        log(f"{species}: done, {entry['triples']} triples")
    except Exception as e:
        entry["error"] = f"{type(e).__name__}: {e}"
        log(f"{species}: FAILED, {entry['error']}")
    return entry


def collect(directory):
    """What a converted species directory holds: files, sizes and triple counts."""
    with open(os.path.join(directory, "conversion.json")) as f:
        report = json.load(f)

    counts = {}
    counts_file = os.path.join(directory, "triple_counts.tsv")
    if os.path.exists(counts_file):
        for line in open(counts_file):
            name, _, count = line.rstrip("\n").partition("\t")
            if count.isdigit():
                counts[name] = int(count)

    files, total = [], 0
    for path in sorted(glob.glob(os.path.join(directory, "*.ttl.gz"))):
        name = os.path.basename(path)[:-3]  # without .gz
        triples = counts.get(name)
        files.append({
            "file": os.path.basename(path),
            "bytes": os.path.getsize(path),
            "triples": triples,
            "statements_written": report["statements"].get(name),
        })
        total += triples or 0

    return {
        "taxonomy_ids": report["taxonomy_ids"],
        "production_names": report["species"],
        "ensembl_version": report["ensembl_version"],
        "entities": report["entities"],
        "files": files,
        "triples": total,
        "converted_at": report["converted_at"],
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("config", help="species YAML (url, species, model, base_uri, entities)")
    parser.add_argument("-o", "--workdir", default=".", help="where the species directories go (default: .)")
    parser.add_argument("-s", "--species", nargs="+", metavar="NAME", action="extend", default=[],
                        help="only these species of the config (default: all of them)")
    parser.add_argument("-j", "--jobs", type=int, default=1,
                        help="species to convert at once; each needs up to 4 GB of RAM (default: 1)")
    parser.add_argument("--skip-download", action="store_true",
                        help="use the tables already downloaded")
    parser.add_argument("--force", action="store_true",
                        help="convert species that are already done")
    parser.add_argument("--stop-on-error", action="store_true",
                        help="stop at the first failure instead of carrying on")
    parser.add_argument("--dry-run", action="store_true",
                        help="report what would be done and exit")
    args = parser.parse_args()

    config_path = os.path.abspath(args.config)
    conf = load_config(config_path)
    workdir = os.path.abspath(args.workdir)
    os.makedirs(workdir, exist_ok=True)

    species = args.species or conf["species"]
    unknown = [s for s in species if s not in conf["species"]]
    if unknown:
        sys.exit(f"Error: {unknown} not listed in {args.config}")

    todo, skipped = [], []
    for name in species:
        if not args.force and is_done(core_dir_of(name, workdir)):
            skipped.append(name)
        else:
            todo.append(name)

    log(f"{len(species)} species: {len(todo)} to convert, {len(skipped)} already done")
    if args.dry_run:
        print("\n".join(todo))
        return 0
    if shutil.which("rapper") is None:
        sys.exit("Error: rapper (raptor) not found in PATH")

    entries = []
    if args.jobs > 1 and not args.stop_on_error:
        with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as pool:
            futures = {pool.submit(convert_species, name, config_path, workdir,
                                   args.skip_download): name for name in todo}
            for future in concurrent.futures.as_completed(futures):
                entries.append(future.result())
    else:
        for name in todo:
            entry = convert_species(name, config_path, workdir, args.skip_download)
            entries.append(entry)
            if entry["status"] == "failed" and args.stop_on_error:
                break

    # Species already done keep their place in the manifest
    for name in skipped:
        directory = core_dir_of(name, workdir)
        entry = {"species": name, "status": "done", "core_db": os.path.basename(directory)}
        try:
            entry.update(collect(directory))
        except Exception as e:
            entry = {"species": name, "status": "unreadable", "error": str(e)}
        entries.append(entry)

    entries.sort(key=lambda e: e["species"])
    failed = [e["species"] for e in entries if e["status"] != "done"]
    manifest = {
        "generated_at": datetime.datetime.now().isoformat(timespec="seconds"),
        "tool": {"name": "ensembl_rdf", "version": tool_version()},
        "config": {
            "path": config_path,
            "url": conf["url"],
            "model": conf["model"] or "ensembl",
            "base_uri": conf["base_uri"],
            "entities": conf["entities"],
            "exclude": conf["exclude"],
        },
        "totals": {
            "species": len(entries),
            "done": len(entries) - len(failed),
            "failed": len(failed),
            "triples": sum(e.get("triples", 0) for e in entries),
            "bytes": sum(f["bytes"] for e in entries for f in e.get("files", [])),
        },
        "species": entries,
    }
    manifest_path = os.path.join(workdir, MANIFEST)
    with open(manifest_path, "w") as f:
        json.dump(manifest, f, indent=2)

    t = manifest["totals"]
    log(f"{t['done']}/{t['species']} species, {t['triples']} triples, "
        f"{t['bytes'] / 1024 ** 3:.2f} GB -> {manifest_path}")
    if failed:
        log(f"failed: {' '.join(failed)}")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())

"""Pairwise orthology/paralogy relations from an Ensembl Compara dump.

For every ENSEMBL_ORTHOLOGUES comparison between two species both listed in a
species configuration, emits one statement per gene pair and per direction:

    gene:A  sio:SIO_000558  gene:B .   # is_orthologous_to
    gene:A  sio:SIO_000630  gene:B .   # is_paralogous_to (within_species_paralog, ...)

`homology_member` alone can exceed 50 GB for a whole division (Ensembl
Plants: ~60 GB compressed), covering every species pair Ensembl computed —
most of no interest here. So `homology.txt.gz` and `homology_member.txt.gz`
are streamed line by line, never loaded; `genome_db`, `species_set`,
`method_link` and `method_link_species_set` are small enough to load fully,
and are used first to narrow the scan to mlss ids between species of
interest, and `gene_member` is streamed once, kept only for species of
interest, to resolve a Compara gene_member_id to the stable_id used
elsewhere in this repo (so the output links directly to the per-species
RDF ensembl_rdf/rdf_converter_ensembl_db.py produces).

Usage:
    download_compara.py https://ftp.ebi.ac.uk/pub/ensemblgenomes/plants/current/mysql/
    compara_orthology.py <compara_dir> -f ensembl_rdf/config/species_agrold.yaml
"""
import os
import sys
import gzip
import argparse
import datetime

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)) + "/../ensembl_rdf")
from utils import iri_escape, percent_encode, turtle_escape, memory_usage_mb  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from rdf_converter_ensembl_db import load_model, Genome2turtle  # noqa: E402
from species_config import load_config  # noqa: E402

# homology.description -> relation kind; anything else is logged once and skipped
HOMOLOGY_KIND = {
    "ortholog_one2one": "ortholog",
    "ortholog_one2many": "ortholog",
    "ortholog_many2many": "ortholog",
    "within_species_paralog": "paralog",
    "other_paralog": "paralog",
    "gene_split": "paralog",
}
ORTHOLOGY_MODEL_KEYS = ["is_orthologous_to", "is_paralogous_to"]


def log(msg):
    print(f'[{datetime.datetime.now().isoformat(timespec="seconds")}] {msg}', file=sys.stderr)


def open_table(compara_dir, table):
    path = os.path.join(compara_dir, table + ".txt.gz")
    if not os.path.exists(path):
        sys.exit(f"Error: {path} not found (see download_compara.py)")
    return gzip.open(path, "rt")


def load_small_table(compara_dir, table, key_idx, value_idxs, as_list=False):
    """Fully load a small lookup table: {key: [values...]} or {key: [[values...], ...]}."""
    result = {}
    with open_table(compara_dir, table) as f:
        for line in f:
            fields = line.rstrip("\n").split("\t")
            key = fields[key_idx]
            values = [fields[i] for i in value_idxs]
            if as_list:
                result.setdefault(key, []).append(values)
            else:
                result[key] = values
    return result


def relevant_mlss_ids(compara_dir, species):
    """{mlss_id: True} for every ENSEMBL_ORTHOLOGUES comparison entirely within `species`."""
    method_link = load_small_table(compara_dir, "method_link", 0, [1])
    orthologues_ids = {k for k, v in method_link.items() if v[0] == "ENSEMBL_ORTHOLOGUES"}
    if not orthologues_ids:
        sys.exit("Error: no ENSEMBL_ORTHOLOGUES entry in method_link.txt.gz")

    genome_db = load_small_table(compara_dir, "genome_db", 0, [2])
    in_scope_genome_db = {k for k, v in genome_db.items() if v[0] in species}

    species_set = load_small_table(compara_dir, "species_set", 0, [1], as_list=True)
    in_scope_species_set = set()
    for ss_id, members in species_set.items():
        genome_db_ids = [m[0] for m in members]
        if genome_db_ids and all(g in in_scope_genome_db for g in genome_db_ids):
            in_scope_species_set.add(ss_id)

    mlss = load_small_table(compara_dir, "method_link_species_set", 0, [1, 2])
    relevant = {mlss_id for mlss_id, (method_link_id, species_set_id) in mlss.items()
                if method_link_id in orthologues_ids and species_set_id in in_scope_species_set}
    log(f"{len(relevant)} method_link_species_set ids between species of interest "
        f"(of {len(mlss)} total)")
    return relevant


def load_gene_members(compara_dir, species, genome_db_name_by_id):
    """Stream gene_member.txt.gz, keeping only genes of the species of interest.

    Returns {gene_member_id: stable_id}.
    """
    in_scope_genome_db_ids = {gid for gid, name in genome_db_name_by_id.items() if name in species}
    kept = {}
    n = 0
    with open_table(compara_dir, "gene_member") as f:
        for line in f:
            n += 1
            fields = line.rstrip("\n").split("\t")
            # gene_member: gene_member_id, stable_id, version, source_name, taxon_id, genome_db_id, ...
            if fields[5] in in_scope_genome_db_ids:
                kept[fields[0]] = fields[1]
    log(f"gene_member: kept {len(kept)} of {n} rows (species of interest)")
    return kept


def load_homology_kinds(compara_dir, relevant_mlss):
    """Stream homology.txt.gz, keeping only rows of a relevant mlss. {homology_id: kind}."""
    kept = {}
    unknown = set()
    n = 0
    with open_table(compara_dir, "homology") as f:
        for line in f:
            n += 1
            fields = line.rstrip("\n").split("\t")
            # homology: homology_id, method_link_species_set_id, description, ...
            if fields[1] not in relevant_mlss:
                continue
            description = fields[2]
            kind = HOMOLOGY_KIND.get(description)
            if kind is None:
                if description not in unknown:
                    log(f"Warning: unknown homology description `{description}`, skipped")
                    unknown.add(description)
                continue
            kept[fields[0]] = kind
    log(f"homology: kept {len(kept)} of {n} rows (relevant mlss and known description)")
    return kept


def gene_uri(stable_id, resources, escape):
    """A gene's URI, as the model profile's `resources.gene` section builds it
    (see Genome2turtle.uri in ensembl_rdf/rdf_converter_ensembl_db.py — kept in
    sync with it so orthology triples reference the same URIs the per-species
    conversion emits).
    """
    gene = resources["gene"]
    if gene["prefix"]:
        return gene["prefix"] + ":" + escape(stable_id)
    return "<" + gene["uri"] + iri_escape(stable_id) + ">"


def resolve_resources(model, base_uri):
    """model["resources"], with {base} filled in — same logic as Genome2turtle.__init__."""
    resources = {}
    for key in Genome2turtle.resource_keys:
        conf = model["resources"].get(key) or {}
        resources[key] = {
            "prefix": conf.get("prefix", ""),
            "uri": conf.get("uri", "{base}/resource/").replace("{base}", base_uri),
        }
    return resources


def convert(compara_dir, species, model, base_uri, output_path,
            include_orthologs=True, include_paralogs=True, source_note=""):
    terms = model["terms"]
    missing = [k for k in ORTHOLOGY_MODEL_KEYS if k not in terms]
    if missing:
        sys.exit(f"Error: model profile does not define {missing}")
    predicate = {"ortholog": terms["is_orthologous_to"], "paralog": terms["is_paralogous_to"]}
    wanted_kinds = {k for k, on in (("ortholog", include_orthologs),
                                   ("paralog", include_paralogs)) if on}

    resources = resolve_resources(model, base_uri)
    escape = percent_encode if model["options"].get("percent_encode_ids", True) else turtle_escape

    genome_db = load_small_table(compara_dir, "genome_db", 0, [2])
    genome_db_name_by_id = {k: v[0] for k, v in genome_db.items()}

    relevant_mlss = relevant_mlss_ids(compara_dir, species)
    gene_members = load_gene_members(compara_dir, species, genome_db_name_by_id)
    homology_kinds = load_homology_kinds(compara_dir, relevant_mlss)
    homology_kinds = {h: k for h, k in homology_kinds.items() if k in wanted_kinds}
    log(f"memory after the small/medium tables: {memory_usage_mb() or 0:.0f} MB")

    # homology_member has exactly two rows per homology_id, not necessarily
    # adjacent, so a pair is completed when the second row of a known
    # homology_id is seen.
    pending = {}
    emitted = 0
    incomplete_other_species = 0
    n = 0
    with open(output_path, "w") as out:
        if source_note:
            print(f"# {source_note}", file=out)
        prefixes = {"sio:": "<http://semanticscience.org/resource/>"}
        if resources["gene"]["prefix"]:
            prefixes[resources["gene"]["prefix"] + ":"] = "<" + resources["gene"]["uri"] + ">"
        for p, u in prefixes.items():
            print(f"@prefix {p} {u} .", file=out)
        print(file=out)

        with open_table(compara_dir, "homology_member") as f:
            for line in f:
                n += 1
                fields = line.rstrip("\n").split("\t")
                # homology_member: homology_id, gene_member_id, seq_member_id, ...
                homology_id, gene_member_id = fields[0], fields[1]
                kind = homology_kinds.get(homology_id)
                if kind is None:
                    continue
                first = pending.pop(homology_id, None)
                if first is None:
                    pending[homology_id] = gene_member_id
                    continue
                stable_a = gene_members.get(first)
                stable_b = gene_members.get(gene_member_id)
                if stable_a is None or stable_b is None:
                    # one member turned out to be outside the species of interest
                    incomplete_other_species += 1
                    continue
                uri_a = gene_uri(stable_a, resources, escape)
                uri_b = gene_uri(stable_b, resources, escape)
                pred = predicate[kind]
                print(uri_a, pred, uri_b, ".", file=out)
                print(uri_b, pred, uri_a, ".", file=out)
                emitted += 1
                if emitted % 1_000_000 == 0:
                    log(f"... {emitted} pairs emitted ({n} homology_member rows read)")

    log(f"homology_member: {n} rows read, {emitted} pairs emitted "
        f"({2 * emitted} statements), {len(pending)} pairs left incomplete "
        f"(should be 0), {incomplete_other_species} pairs skipped (member outside scope)")
    log(f"Done: {output_path}")


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("compara_dir", help="the downloaded *_compara_* directory")
    parser.add_argument("-f", "--species-file", required=True,
                        help="species YAML (species list, model, base_uri — "
                             "see ensembl_rdf/config/species_agrold.yaml)")
    parser.add_argument("-m", "--model", help="override the species file's model")
    parser.add_argument("-b", "--base-uri", help="override the species file's base_uri")
    parser.add_argument("-o", "--output", help="output file (default: COMPARA_DIR/orthology.ttl)")
    parser.add_argument("--orthologs-only", action="store_true")
    parser.add_argument("--paralogs-only", action="store_true")
    args = parser.parse_args()
    if args.orthologs_only and args.paralogs_only:
        sys.exit("Error: --orthologs-only and --paralogs-only are mutually exclusive")

    conf = load_config(args.species_file)
    species = set(conf["species"])
    if not species:
        sys.exit(f"Error: {args.species_file} lists no species")

    model_name = args.model or conf["model"] or "ensembl"
    models_dir = os.path.join(Genome2turtle.base_dir, "config", "models")
    model_path = model_name if os.path.exists(model_name) else \
        os.path.join(models_dir, model_name + ".yaml")
    if not os.path.exists(model_path):
        sys.exit(f"Error: unknown model '{model_name}'")
    model = load_model(model_path)
    base_uri = (args.base_uri or conf["base_uri"] or Genome2turtle.default_base_uri)

    output_path = args.output or os.path.join(args.compara_dir, "orthology.ttl")
    log(f"Species: {len(species)} ({', '.join(sorted(species)[:5])}{'...' if len(species) > 5 else ''})")
    log(f"Model: {model_path}")
    log(f"Base URI: {base_uri}")

    convert(args.compara_dir, species, model, base_uri, output_path,
            include_orthologs=not args.paralogs_only,
            include_paralogs=not args.orthologs_only,
            source_note=f"Derived from {os.path.basename(args.compara_dir)} "
                        f"(ENSEMBL_ORTHOLOGUES), {len(species)} species")


if __name__ == "__main__":
    main()

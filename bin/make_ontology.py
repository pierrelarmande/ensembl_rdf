"""Derive the ontology of a vocabulary profile from ontology/ensembl_ontology.ttl.

The source ontology defines the classes and properties of the Ensembl model in
its own namespace (`terms:`). A profile (config/models/*.yaml) may rename them
and publish them elsewhere, e.g. terms:EnsemblGene becomes
agrold_vocabulary:Gene; this script rewrites the ontology accordingly, so the
vocabulary the converter emits is actually published somewhere.

Terms the profile maps outside its own vocabulary (to SIO, SO, ...) are dropped:
they are defined by the ontology they come from. Each renamed term keeps an
owl:equivalentClass / owl:equivalentProperty link to the Ensembl term it
derives from.

Usage: make_ontology.py MODEL [-o OUTPUT]
"""
import os
import re
import sys
import argparse

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from rdf_converter_ensembl_db import Ensembl2turtle, load_model

BASE_DIR = os.path.dirname(os.path.abspath(__file__)) + "/../"
SOURCE_ONTOLOGY = BASE_DIR + "ontology/ensembl_ontology.ttl"
ENSEMBL_MODEL = Ensembl2turtle.models_dir + "ensembl.yaml"


def local_name(term, vocabulary_prefix):
    """Local part of a term written in the vocabulary namespace, else None."""
    prefix, _, name = term.partition(":")
    return name if prefix == vocabulary_prefix and name else None


def build_renaming(source, target):
    """Map each source local name to its target local name, or None if dropped."""
    src_prefix = source.get("vocabulary_prefix") or "terms"
    dst_prefix = target.get("vocabulary_prefix") or "terms"
    renaming = {}
    for key, src_term in source["terms"].items():
        src_name = local_name(src_term, src_prefix)
        if src_name is None:
            continue  # not defined by this ontology (so:, sio:, dcterms:, ...)
        renaming[src_name] = local_name(target["terms"].get(key, ""), dst_prefix)
    return renaming


def words_of(name):
    """Turn a term name into the words of its label: EnsemblGene -> "Ensembl gene"."""
    words = re.sub(r"(?<!^)(?=[A-Z])", " ", name).replace("_", " ").split()
    return " ".join([words[0]] + [w.lower() for w in words[1:]])


def block_subject(block):
    match = re.match(r"\s*:([A-Za-z_][A-Za-z0-9_]*)", block)
    return match.group(1) if match else None


def rewrite(text, renaming, source_uri, target_uri, target_prefix):
    blocks = text.split("\n\n")
    kept, dropped = [], []
    for block in blocks:
        subject = block_subject(block)
        if subject is not None and subject in renaming and renaming[subject] is None:
            dropped.append(subject)
            continue
        kept.append(block)
    text = "\n\n".join(kept)

    # Rename the terms the profile renames, longest first to avoid partial hits
    for old in sorted(renaming, key=len, reverse=True):
        new = renaming[old]
        if new is None or new == old:
            continue
        text = re.sub(r":%s\b" % re.escape(old), ":" + new, text)
        # Labels follow the term name: :EnsemblGene "Ensembl gene" -> :Gene "gene"
        old_words, new_words = words_of(old), words_of(new)
        if old_words != new_words:
            text = text.replace('"%s"' % old_words, '"%s"' % new_words)

    text = text.replace("@prefix : <%s> ." % source_uri, "@prefix : <%s> ." % target_uri)
    text = text.replace("<%s> a owl:Ontology ." % source_uri, "<%s> a owl:Ontology ." % target_uri)

    # Traceability: link each renamed term to the Ensembl term it derives from
    equivalences = ""
    for old in sorted(renaming):
        new = renaming[old]
        if new is None or new == old:
            continue
        predicate = "owl:equivalentClass" if old[0].isupper() else "owl:equivalentProperty"
        equivalences += ":%s %s <%s%s> .\n" % (new, predicate, source_uri, old)
    if equivalences:
        text = text.rstrip("\n") + "\n\n\n# Terms of the Ensembl vocabulary these are derived from\n" + equivalences

    return text, dropped


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("model", help="vocabulary profile: a name in config/models/ or a YAML file")
    parser.add_argument("-o", "--output", help="output file (default: ontology/<model>_ontology.ttl)")
    parser.add_argument("--source", default=SOURCE_ONTOLOGY,
                        help="ontology to derive from (default: ontology/ensembl_ontology.ttl)")
    args = parser.parse_args()

    model_path = args.model
    name = os.path.basename(args.model)
    if not os.path.exists(model_path):
        model_path = Ensembl2turtle.models_dir + args.model + ".yaml"
        if not os.path.exists(model_path):
            sys.exit(f"Error: unknown model '{args.model}'")
    else:
        name = re.sub(r"\.ya?ml$", "", name)

    source = load_model(ENSEMBL_MODEL)
    target = load_model(model_path)
    source_uri = source.get("vocabulary_uri") or Ensembl2turtle.default_terms_uri
    target_uri = target.get("vocabulary_uri") or Ensembl2turtle.default_terms_uri
    target_prefix = target.get("vocabulary_prefix") or "terms"

    renaming = build_renaming(source, target)
    with open(args.source, "r") as f:
        text = f.read()
    text, dropped = rewrite(text, renaming, source_uri, target_uri, target_prefix)

    output = args.output or BASE_DIR + "ontology/" + name + "_ontology.ttl"
    with open(output, "w") as f:
        f.write(text)
    renamed = {k: v for k, v in renaming.items() if v and v != k}
    print(f"{output}: {len(renamed)} terms renamed, {len(dropped)} dropped "
          f"({', '.join(dropped) if dropped else 'none'})", file=sys.stderr)


if __name__ == "__main__":
    main()

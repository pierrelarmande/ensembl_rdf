"""Species selection shared by download_files.py and convert.sh.

A species is given by its production name (directory prefix before `_core_`,
e.g. `arabidopsis_thaliana`); a regular expression is accepted too.
"""
import re
import sys


def load_species_file(path):
    """Return (url, [species]) from a YAML file with keys `url` (optional) and `species`."""
    try:
        import yaml
    except ImportError:
        sys.exit("Error: PyYAML is required to read a species file (pip install pyyaml)")
    with open(path, "r") as f:
        conf = yaml.safe_load(f) or {}
    species = conf.get("species") or []
    if isinstance(species, str):
        species = [species]
    return conf.get("url"), [str(s) for s in species]


def match_core_dir(species_patterns, dirname):
    """True if `dirname` is a core DB directory of one of the given species."""
    if "_core_" not in dirname:
        return False
    if not species_patterns:
        return True
    return any(re.match(p + "_core_", dirname) for p in species_patterns)


if __name__ == "__main__":
    # Usage: species_config.py FILE [DIR ...]
    # Prints the matching core directories among DIR (or the species list if none given).
    url, species = load_species_file(sys.argv[1])
    dirs = sys.argv[2:]
    if dirs:
        print("\n".join(d for d in dirs if match_core_dir(species, d.rstrip("/"))))
    else:
        print("\n".join(species))

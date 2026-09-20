"""Species selection shared by download_files.py and convert.sh.

A species is given by its production name (directory prefix before `_core_`,
e.g. `arabidopsis_thaliana`); a regular expression is accepted too.
"""
import re
import sys


def _as_list(value):
    if value is None:
        return []
    if isinstance(value, str):
        return [value]
    return [str(v) for v in value]


def load_config(path):
    """Return the YAML config as a dict with keys url, species, entities, exclude."""
    try:
        import yaml
    except ImportError:
        sys.exit("Error: PyYAML is required to read a species file (pip install pyyaml)")
    with open(path, "r") as f:
        conf = yaml.safe_load(f) or {}
    return {
        "url": conf.get("url"),
        "base_uri": conf.get("base_uri"),
        "terms_uri": conf.get("terms_uri"),
        "species": _as_list(conf.get("species")),
        "entities": _as_list(conf.get("entities")),
        "exclude": _as_list(conf.get("exclude")),
    }


def load_species_file(path):
    """Return (url, [species]) from a YAML file with keys `url` (optional) and `species`."""
    conf = load_config(path)
    return conf["url"], conf["species"]


def match_core_dir(species_patterns, dirname):
    """True if `dirname` is a core DB directory of one of the given species."""
    if "_core_" not in dirname:
        return False
    if not species_patterns:
        return True
    return any(re.match(p + "_core_", dirname) for p in species_patterns)


if __name__ == "__main__":
    # Usage: species_config.py dirs FILE [DIR ...]   -> core directories of the file's species among DIR
    #        species_config.py entities FILE         -> the file's `entities` list (space separated)
    #        species_config.py exclude FILE          -> the file's `exclude` list (space separated)
    #        species_config.py base_uri|terms_uri FILE -> that value (empty if unset)
    if len(sys.argv) < 3 or sys.argv[1] not in ("dirs", "entities", "exclude", "base_uri", "terms_uri"):
        sys.exit(__doc__)
    conf = load_config(sys.argv[2])
    if sys.argv[1] == "dirs":
        print("\n".join(d for d in sys.argv[3:] if match_core_dir(conf["species"], d.rstrip("/"))))
    elif sys.argv[1] in ("base_uri", "terms_uri"):
        print(conf[sys.argv[1]] or "")
    else:
        print(" ".join(conf[sys.argv[1]]))

import re
import sys
import urllib.parse
from datetime import datetime

def quote_str(string):
    # Enclose a string with double-quotation marks.
    # Double-quotation marks within the input string are escaped with backslashes.
    return "\"" + string.replace("\"", "\\\"") + "\""


def percent_encode(string):
    #return re.sub(r"([()])", r"\\\1", string)
    return urllib.parse.quote(string, safe="")


def turtle_escape(string):
    """Escape a stable ID for the local part of a Turtle prefixed name.

    Unlike percent_encode, this leaves the identifier itself untouched — the
    IRI still holds `transcript-3'rps12-E1` — and only escapes what the Turtle
    grammar requires: the characters of PN_LOCAL_ESC, and a trailing dot, which
    would otherwise end the statement.
    """
    string = re.sub(r"([~!$&'()*+,;=/?#@%])", r"\\\1", string)
    return re.sub(r"\.$", r"\\.", string)


def iri_escape(string):
    """Percent-encode only what an IRI written as <...> cannot hold."""
    return re.sub(r'[<>"{}|^`\\\s]', lambda m: "%%%02X" % ord(m.group()), string)


def memory_usage_mb():
    """Resident memory of this process, or None when it cannot be read.

    psutil gives it on every platform; without it `resource` is close enough on
    Linux and macOS, and a missing figure is no reason to fail a conversion.
    """
    try:
        import psutil
        return psutil.Process().memory_info().rss / (1024 * 1024)
    except ImportError:
        pass
    try:
        import resource
        peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        # Linux reports kilobytes, macOS and the BSDs bytes
        return peak / 1024 if sys.platform.startswith("linux") else peak / (1024 * 1024)
    except Exception:
        return None


def strand2faldo(s):
    if s == "1":
        return "faldo:ForwardStrandPosition"
    elif s == "-1":
        return "faldo:ReverseStrandPosition"
    else:
        print(f"Error: Invalid argument \"{s}\" for strand2faldo", file=sys.stderr)
        sys.exit(1)


def log_time(text):
    dt_now = datetime.now()
    print(f"[{dt_now}] {text}", file=sys.stderr)
    return


class Bnode:
    def __init__(self):
        self.properties = []

    def add(self, tpl):
        # `tpl` is a tuple of strings (e.g. ("rdf:type", "owl:Class"))
        self.properties.append(tpl)

    def serialize(self, level=1):
        s = "[\n"
        indent = "    " * level
        for tpl in self.properties:
            s += indent + tpl[0] + " " + tpl[1] + " ;\n"
        s += "    " * (level-1) + "]"
        return s

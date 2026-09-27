import os
import sys
import json
import time
import datetime
import argparse
from ftplib import FTP
from html.parser import HTMLParser
from urllib.parse import urlparse, urljoin
from urllib.request import urlopen

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from species_config import load_species_file, match_core_dir


CONFIG_DIR = os.path.dirname(os.path.abspath(__file__)) + "/../config/"


def log(msg):
    dt_now = datetime.datetime.now()
    print(f'[{dt_now}] {msg}', file=sys.stderr)


def retry(what, call, attempts=5, delay=10):
    """Run `call`, retrying on network errors with a growing delay.

    Long unattended runs hit transient failures (the server closing the
    connection, a timeout); one failure should not lose a whole job.
    """
    for attempt in range(1, attempts + 1):
        try:
            return call()
        except Exception as e:
            if attempt == attempts:
                raise
            wait = delay * attempt
            log(f'Warning: {what} failed ({type(e).__name__}: {e}); '
                f'retry {attempt}/{attempts - 1} in {wait}s')
            time.sleep(wait)


class LinkParser(HTMLParser):
    """Collect href values of an HTTP directory listing."""
    def __init__(self):
        super().__init__()
        self.links = []

    def handle_starttag(self, tag, attrs):
        if tag == "a":
            for name, value in attrs:
                if name == "href" and value:
                    self.links.append(value)


class HttpSource:
    def __init__(self, base_url):
        self.base_url = base_url if base_url.endswith("/") else base_url + "/"

    def list(self, path=""):
        url = urljoin(self.base_url, path)
        with urlopen(url) as res:
            parser = LinkParser()
            parser.feed(res.read().decode("utf-8", errors="replace"))
        # Keep relative entries only (skip parent dir, sort links, absolute URLs)
        entries = [l.rstrip("/") for l in parser.links
                   if not l.startswith(("/", "?", "http")) and l not in ("../", "./")]
        return entries

    def fetch(self, path, dest):
        url = urljoin(self.base_url, path)
        with urlopen(url) as res, open(dest, 'wb') as f:
            while True:
                chunk = res.read(1024 * 1024)
                if not chunk:
                    break
                f.write(chunk)

    def close(self):
        return


class FtpSource:
    def __init__(self, host, base_dir):
        self.host = host
        self.base_dir = base_dir
        self.connect()

    def connect(self):
        self.ftp = FTP(self.host)
        self.ftp.login()

    def list(self, path=""):
        target = self.base_dir + "/" + path if path else self.base_dir
        try:
            return self.ftp.nlst(target)
        except Exception:
            self.connect()  # the control connection may have timed out
            return self.ftp.nlst(target)

    def fetch(self, path, dest):
        with open(dest, 'wb') as f:
            self.ftp.retrbinary('RETR %s' % (self.base_dir + "/" + path), f.write)

    def close(self):
        self.ftp.quit()


def make_source(args):
    # Backward compatible form: `download_files.py ftp.ensembl.org /pub/current_mysql/`
    if args.directory:
        return FtpSource(args.url, args.directory.rstrip("/"))
    parsed = urlparse(args.url)
    if parsed.scheme in ("http", "https"):
        return HttpSource(args.url)
    elif parsed.scheme == "ftp":
        return FtpSource(parsed.netloc, parsed.path.rstrip("/"))
    else:
        sys.exit(f"Error: unsupported URL scheme '{parsed.scheme}' (use http(s):// or ftp://)")


def download_files(source, directory, dbs):
    files = [os.path.basename(f) for f in retry(f'listing {directory}',
                                                lambda: source.list(directory))]
    os.makedirs(directory, exist_ok=True)
    for file in files:
        if file in dbs:
            path = directory + "/" + file
            log(f'Downloading: {path}')
            retry(f'downloading {path}',
                  lambda: source.fetch(path, path + ".part"))
            os.replace(path + ".part", path)


def process_directory(source, dbs, species_patterns):
    subdirectories = [os.path.basename(d) for d in retry('listing the release',
                                                         lambda: source.list())]
    found = False
    for subdirectory in subdirectories:
        if not match_core_dir(species_patterns, subdirectory):
            continue
        found = True
        log(subdirectory)
        download_files(source, subdirectory, dbs)
    if not found:
        log(f"Warning: no *_core_* directory matched {species_patterns}")


def main():
    parser = argparse.ArgumentParser(
        description="Download Ensembl core MySQL dumps listed in config/dbinfo.json.",
        epilog="Examples:\n"
               "  %(prog)s https://ftp.ebi.ac.uk/pub/ensemblgenomes/plants/current/mysql/\n"
               "  %(prog)s ftp.ensembl.org /pub/current_mysql/\n"
               "  %(prog)s https://ftp.ebi.ac.uk/pub/ensemblgenomes/plants/current/mysql/ -s arabidopsis_thaliana oryza_sativa\n"
               "  %(prog)s -f config/species.yaml",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("url", nargs="?", help="base URL of the mysql dump directory (https:// or ftp://), or an FTP host; "
                                               "may be omitted if the species file has a `url` key")
    parser.add_argument("directory", nargs="?", help="FTP directory (legacy form, with an FTP host as first argument)")
    parser.add_argument("-s", "--species", nargs="+", default=[], metavar="NAME",
                        help="only download these species (production name, e.g. arabidopsis_thaliana; a regex is accepted)")
    parser.add_argument("-f", "--species-file", metavar="FILE",
                        help="YAML file with a `species` list and optionally a `url` (see config/species.yaml)")
    parser.add_argument("--dbinfo", default=CONFIG_DIR + "dbinfo.json", help="dbinfo.json (default: config/dbinfo.json)")
    args = parser.parse_args()

    species = list(args.species)
    if args.species_file:
        url, file_species = load_species_file(args.species_file)
        # -s selects the species; the file then only supplies the URL
        if not species:
            species = file_species
        if not args.url:
            args.url = url
    if not args.url:
        parser.error("no URL given (either on the command line or as `url` in the species file)")

    with open(args.dbinfo, "r") as f:
        dbinfo = json.load(f)
    dbs = [dbinfo[k]["filename"] for k in dbinfo]

    source = make_source(args)
    process_directory(source, dbs, species)
    source.close()


if __name__ == '__main__':
    main()

# Ensembl RDF converter by DBCLS
## Data download
Give the base URL (https:// or ftp://) of the MySQL dump directory. Only the tables listed in `config/dbinfo.json` are downloaded, for every `*_core_*` database.
```
# Ensembl Plants
$ python3 /path/to/ensembl_rdf/bin/download_files.py https://ftp.ebi.ac.uk/pub/ensemblgenomes/plants/current/mysql/
# Ensembl (vertebrates)
$ python3 /path/to/ensembl_rdf/bin/download_files.py https://ftp.ensembl.org/pub/current_mysql/
```
The legacy form `download_files.py ftp.ensembl.org /pub/current_mysql/` (FTP host + directory) is still accepted.

### Selecting species
By default every `*_core_*` database is downloaded. To restrict to some species, give their production names (the directory name before `_core_`; a regular expression is accepted) with `-s`, or list them in a YAML file given with `-f` (see [config/species.yaml](config/species.yaml), which may also hold the `url`):
```
$ python3 /path/to/ensembl_rdf/bin/download_files.py https://ftp.ebi.ac.uk/pub/ensemblgenomes/plants/current/mysql/ -s arabidopsis_thaliana oryza_sativa
$ python3 /path/to/ensembl_rdf/bin/download_files.py -f /path/to/ensembl_rdf/config/species.yaml
```
Reading a species file requires [PyYAML](https://pypi.org/project/PyYAML/).

## RDF conversion
In the directory where `download_files.py` was executed, run `convert.sh` with the directories to convert, or the same species selection (`-s` / `-f`) as above:
```
$ bash /path/to/ensembl_rdf/bin/convert.sh *_core_*
$ bash /path/to/ensembl_rdf/bin/convert.sh -s arabidopsis_thaliana
$ bash /path/to/ensembl_rdf/bin/convert.sh -f /path/to/ensembl_rdf/config/species.yaml
```
### Selecting entity types
Each entity type is written to its own Turtle file: `gene`, `transcript`, `translation`, `exon`, `exon_transcript`, `xref`. All are output by default; restrict with `-e` (only these) or `-x` (all but these), on the command line or as `entities` / `exclude` lists in the YAML file. Only the tables needed by the selected types are loaded.
```
$ bash /path/to/ensembl_rdf/bin/convert.sh -s oryza_sativa -x exon -x exon_transcript
$ bash /path/to/ensembl_rdf/bin/convert.sh -s oryza_sativa -e gene -e transcript
```
The converter itself accepts the same selection: `rdf_converter_ensembl_db.py config/dbinfo.json -x exon exon_transcript`.

`convert.sh` requires [rapper](https://librdf.org/raptor/rapper.html) (Raptor RDF Syntax Library) to normalize the Turtle files.

## RDF schema
[RDF-config](https://github.com/dbcls/rdf-config/blob/master/config/ensembl/model.yaml)

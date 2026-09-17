# Ensembl RDF converter by DBCLS
## Data download
Give the base URL (https:// or ftp://) of the MySQL dump directory. Only the tables listed in `config/dbinfo.json` are downloaded, for every `*_core_*` database.
```
# Ensembl Plants
$ python3 /path/to/ensembl_rdf/bin/download_files.py https://ftp.ebi.ac.uk/pub/ensemblgenomes/plants/current/mysql/
# Ensembl (vertebrates)
$ python3 /path/to/ensembl_rdf/bin/download_files.py https://ftp.ensembl.org/pub/current_mysql/
```
Use `-s PATTERN [PATTERN ...]` to restrict to species directories matching a regex, e.g. `-s '^arabidopsis_thaliana_' '^oryza_sativa_'`.
The legacy form `download_files.py ftp.ensembl.org /pub/current_mysql/` (FTP host + directory) is still accepted.
## RDF conversion
In the directory where `download_files.py` were executed, run `convert.sh`.
```
$ bash /path/to/ensembl_rdf/bin/convert.sh *_core_*
```

## RDF schema
[RDF-config](https://github.com/dbcls/rdf-config/blob/master/config/ensembl/model.yaml)

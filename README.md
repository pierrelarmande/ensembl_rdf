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
Each entity type is written to its own Turtle file: `gene`, `transcript`, `translation`, `exon`, `exon_transcript`, `xref`, `chromosome` (the regions the other entities sit on — only those actually referenced, so the contigs and scaffolds of a `seq_region` table are left out; it is written last for that reason). All are output by default; restrict with `-e` (only these) or `-x` (all but these), on the command line or as `entities` / `exclude` lists in the YAML file. Only the tables needed by the selected types are loaded. Excluding `exon` also excludes `exon_transcript`, which only links transcripts to exons.
```
$ bash /path/to/ensembl_rdf/bin/convert.sh -s oryza_sativa -x exon
$ bash /path/to/ensembl_rdf/bin/convert.sh -s oryza_sativa -e gene -e transcript
```
The converter itself accepts the same selection: `rdf_converter_ensembl_db.py config/dbinfo.json -x exon`.

### Base URI
Resource URIs are built on `http://rdf.ebi.ac.uk` by default (`http://rdf.ebi.ac.uk/resource/ensembl/ENSG...`, `.../resource/ensembl.transcript/...`, chromosomes `.../resource/ensembl/116/oryza_sativa/IRGSP-1.0/1`). Change the base with `-b` or `base_uri` in the YAML file:
```
$ bash /path/to/ensembl_rdf/bin/convert.sh -s oryza_sativa -b http://purl.agrold.org
```
gives `http://purl.agrold.org/resource/ensembl/...`. The vocabulary of the model itself (classes and properties defined in [ontology/ensembl_ontology.ttl](ontology/ensembl_ontology.ttl)) is set by the model profile, or by `-t` / `terms_uri`; the ontology file must then be republished with the same namespace.

### Vocabulary profile
The classes and properties used for the model are not hard-coded: a profile in [config/models/](config/models/) maps each element to a term, selected with `-m` / `model`.
```
$ bash /path/to/ensembl_rdf/bin/convert.sh -s oryza_sativa -m ensembl   # default, Ensembl RDF as published by the EBI
$ bash /path/to/ensembl_rdf/bin/convert.sh -s oryza_sativa -m agrold -b http://purl.agrold.org
```
`agrold` names the classes after the feature (`Gene`, `Transcript`, `Protein`, `Exon`) and keeps the Ensembl properties, everything published under `http://purl.agrold.org/vocabulary/`. Pass a YAML file instead of a name to use your own profile.

A profile also sets the URIs of the resources themselves, so they need not mention `ensembl`:
```yaml
resources:
  gene:       {prefix: gene,       uri: "{base}/resource/"}            # gene:AT1G01010
  transcript: {prefix: transcript, uri: "{base}/resource/transcript/"}
  protein:    {prefix: protein,    uri: "{base}/resource/protein/"}
  exon:       {prefix: exon,       uri: "{base}/resource/exon/"}
  chromosome: {uri: "{base}/resource/chromosome/"}
```
`{base}` is replaced by `-b`; `prefix` is the Turtle prefix (omit it to write full IRIs).

The shape of a region URI and the form of the FALDO locations are profile options:
```yaml
options:
  chromosome_pattern: "{taxon}/{assembly}/{chromosome}"   # {version}, {production_name} also available
  faldo_named_regions: true    # <chromosome>:<begin>-<end>:<strand> instead of blank nodes
  biotype_as_literal: false    # true writes has_biotype "protein_coding"
```
With `faldo_named_regions`, regions and positions are named resources that features sharing coordinates share, as in the AgroLD graph; the `ensembl` profile keeps the blank nodes the EBI publishes.

The ontology of a profile is derived from [ontology/ensembl_ontology.ttl](ontology/ensembl_ontology.ttl) — renamed terms, profile namespace, and an `owl:equivalentClass` / `owl:equivalentProperty` link back to the Ensembl term each one derives from:
```
$ python3 /path/to/ensembl_rdf/bin/make_ontology.py agrold      # -> ontology/agrold_ontology.ttl
```
Regenerate it whenever the profile changes, and publish it at the namespace the profile declares.

`convert.sh` requires [rapper](https://librdf.org/raptor/rapper.html) (Raptor RDF Syntax Library) to normalize the Turtle files.

### Species of a division
[config/species_agrold.yaml](config/species_agrold.yaml) lists the 51 plant species of AgroLD, matched to Ensembl Plants by NCBI taxon ID against [config/species_EnsemblPlants_63.txt](config/species_EnsemblPlants_63.txt) (the release's own species table, `species_EnsemblPlants.txt` on the FTP). Where several genomes share a taxon ID — 18 bread wheat cultivars, 7 barleys, 10 indica rices — the reference genome is listed; matching is anchored on `<name>_core_`, so `triticum_aestivum` does not pull in `triticum_aestivum_cadenza`.

## Running on a cluster
[bin/convert_slurm.sh](bin/convert_slurm.sh) converts one species per array task — download, conversion, `rapper` and `gzip`:
```
$ mkdir -p logs
$ sbatch --array=1-$(python3 bin/species_config.py species config/species_agrold.yaml | wc -l) \
         bin/convert_slurm.sh config/species_agrold.yaml /path/to/workdir
```
Tasks are independent, so a failed one is resubmitted alone with `--array=<n>`. Edit the `#SBATCH` header and the commented `module load` / `conda activate` lines for your site.

Sizing, measured on Ensembl Plants 63: the 51 species of `species_agrold.yaml` amount to 1.3 GB of compressed dumps. The converter loads the tables it needs in RAM, roughly 35-40x their compressed size (Arabidopsis: 44 MB of dumps, 1.6 GB of RAM, 39 s, 560 MB of Turtle before `rapper` and `gzip`, 28 MB after); 16 GB per task covers every plant species, bread wheat included.

## Cross-reference sources
Cross-references are turned into `rdfs:seeAlso` links using [config/external_db_url.tsv](config/external_db_url.tsv): one line per `external_db_id` (numeric, not the name), giving the URL the accession is appended to and, optionally, a prefix to strip from it. Sources with no URL are reported in `xref_report.tsv` next to each converted database, and no link is emitted for them.

Before converting a whole division, survey which sources its species actually use — it only needs `xref.txt.gz` and `external_db.txt.gz`:
```
$ python3 /path/to/ensembl_rdf/bin/download_files.py URL --dbinfo /path/to/ensembl_rdf/config/xref_dbinfo.json
$ python3 /path/to/ensembl_rdf/bin/xref_survey.py                  # every *_core_* of the cwd
$ python3 /path/to/ensembl_rdf/bin/xref_survey.py --missing-only   # only the sources with no URL
$ python3 /path/to/ensembl_rdf/bin/xref_survey.py --missing-only --tsv   # ids to append to the table
```
The report gives, per source, the number of species using it, the number of xrefs, an example accession and its current URL, so the sources worth mapping come first. Ids already listed with an empty URL are filled in place; ids the table does not know are printed by `--tsv`, ready to append.

## RDF schema
[RDF-config](https://github.com/dbcls/rdf-config/blob/master/config/ensembl/model.yaml)

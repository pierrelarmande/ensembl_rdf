# Ensembl RDF converter by DBCLS

## Requirements
Python 3.8 or later with [PyYAML](https://pypi.org/project/PyYAML/), and [rapper](https://librdf.org/raptor/rapper.html) from the Raptor RDF Syntax Library, which `convert.sh` uses to check and normalize the Turtle files:
```
$ pip install -r requirements.txt
$ apt install raptor2-utils        # Debian/Ubuntu; brew install raptor on macOS
$ conda create -n rapper -c conda-forge 'raptor=2'      # or with conda
```
Beware of a homonym: on **bioconda**, `raptor` is SeqAn's sequence pre-filter (3.x) and provides no `rapper`; the RDF library is the **conda-forge** `raptor` (2.x). Check with `rapper --version`.

`rapper` validates and normalizes the Turtle but does not produce it: without it the conversion still writes valid files, only unnormalized, unchecked, and with no triple count in the manifest. A warning says so.

`psutil` is optional: it reports resident memory while the tables load, and the converter falls back to the standard library without it.

## Data download
Give the base URL (https:// or ftp://) of the MySQL dump directory. Only the tables listed in `ensembl_rdf/config/dbinfo.json` are downloaded, for every `*_core_*` database.
```
# Ensembl Plants
$ python3 /path/to/ensembl_rdf/bin/download_files.py https://ftp.ebi.ac.uk/pub/ensemblgenomes/plants/current/mysql/
# Ensembl (vertebrates)
$ python3 /path/to/ensembl_rdf/bin/download_files.py https://ftp.ensembl.org/pub/current_mysql/
```
The legacy form `download_files.py ftp.ensembl.org /pub/current_mysql/` (FTP host + directory) is still accepted. Network failures are retried, and each file is downloaded to a `.part` renamed on completion, so an interrupted run leaves nothing truncated.

### Selecting species
By default every `*_core_*` database is downloaded. To restrict to some species, give their production names (the directory name before `_core_`; a regular expression is accepted) with `-s`, or list them in a YAML file given with `-f` (see [ensembl_rdf/config/species.yaml](ensembl_rdf/config/species.yaml), which may also hold the `url`, the model and the base URI). An explicit `-s` selects the species and the file then only supplies the options.
```
$ python3 .../bin/download_files.py https://ftp.ebi.ac.uk/pub/ensemblgenomes/plants/current/mysql/ -s arabidopsis_thaliana oryza_sativa
$ python3 .../bin/download_files.py -f /path/to/ensembl_rdf/ensembl_rdf/config/species_agrold.yaml
```
[ensembl_rdf/config/species_agrold.yaml](ensembl_rdf/config/species_agrold.yaml) lists the 51 plant species of AgroLD, matched to Ensembl Plants by NCBI taxon ID against [ensembl_rdf/config/species_EnsemblPlants_63.txt](ensembl_rdf/config/species_EnsemblPlants_63.txt). Where several genomes share a taxon ID — 18 bread wheat cultivars, 7 barleys, 10 indica rices — the reference genome is listed; matching is anchored on `<name>_core_`, so `triticum_aestivum` does not pull in `triticum_aestivum_cadenza`.

## RDF conversion
In the directory where `download_files.py` was executed, run `convert.sh` with the directories to convert, or the same species selection (`-s` / `-f`) as above:
```
$ bash /path/to/ensembl_rdf/bin/convert.sh *_core_*
$ bash /path/to/ensembl_rdf/bin/convert.sh -s arabidopsis_thaliana
$ bash /path/to/ensembl_rdf/bin/convert.sh -f /path/to/ensembl_rdf/ensembl_rdf/config/species_agrold.yaml
```

### Selecting entity types
Each entity type is written to its own Turtle file: `gene`, `transcript`, `translation`, `exon`, `exon_transcript`, `xref`, `chromosome` (the regions the other entities sit on — only those actually referenced, so the contigs and scaffolds of a `seq_region` table are left out; it is written last for that reason). All are output by default; restrict with `-t` (only these) or `-x` (all but these). Only the tables the selected types need are loaded, so converting part of the model costs part of the memory. Excluding `exon` also excludes `exon_transcript`, which only links transcripts to exons.
```
$ bash /path/to/ensembl_rdf/bin/convert.sh -t gene -t transcript *_core_*
$ bash /path/to/ensembl_rdf/bin/convert.sh -x exon *_core_*
```
The python script takes the same selection, and only one `-t` is needed for several types:
```
$ python3 ensembl_rdf/rdf_converter_ensembl_db.py ensembl_rdf/config/dbinfo.json path/to/input/ -t gene transcript
```

### Base URI
Resource URIs are built on `http://rdf.ebi.ac.uk` by default. Change the base with `-b` or `base_uri` in the YAML file:
```
$ bash /path/to/ensembl_rdf/bin/convert.sh -s oryza_sativa -b http://purl.agrold.org
```

### Vocabulary profile
The classes and properties used for the model are not hard-coded: a profile in [ensembl_rdf/config/models/](ensembl_rdf/config/models/) maps each element to a term, selected with `-m` / `model`.
```
$ bash .../bin/convert.sh -s oryza_sativa -m ensembl   # default, Ensembl RDF as published by the EBI
$ bash .../bin/convert.sh -s oryza_sativa -m agrold -b http://purl.agrold.org
```
`agrold` names the classes after the feature (`Gene`, `Transcript`, `Protein`, `Exon`) and keeps the Ensembl properties, everything published under `http://purl.agrold.org/vocabulary/`. A profile also decides the resource URIs, so they need not mention `ensembl`, and a few things the two models express differently:
```yaml
options:
  percent_encode_ids: false    # keep the stable ID in the URI, escaping only for Turtle
  in_taxon_object: taxonomy    # the NCBI taxon, not the Ensembl genome URI
  chromosome_pattern: "{taxon}/{assembly}/{chromosome}"
  faldo_named_regions: true    # <chromosome>:<begin>-<end>:<strand> instead of blank nodes
  biotype_as_literal: false    # true writes has_biotype "protein_coding"
```
The ontology of a profile is derived from [ensembl_rdf/ontology/ensembl_ontology.ttl](ensembl_rdf/ontology/ensembl_ontology.ttl) — renamed terms, profile namespace, an `owl:equivalentClass` back to the Ensembl term each one derives from, and a declaration for the terms the profile introduces:
```
$ python3 /path/to/ensembl_rdf/bin/make_ontology.py agrold   # -> ensembl_rdf/ontology/agrold_ontology.ttl
```
Regenerate it whenever the profile changes, and publish it at the namespace the profile declares.

## Running a whole configuration
[bin/run.py](bin/run.py) is the single entry point a pipeline calls: it downloads, converts, checks with `rapper`, gzips and writes a manifest.
```
$ python3 /path/to/ensembl_rdf/bin/run.py ensembl_rdf/config/species_agrold.yaml -o /path/to/workdir
```
A species already converted is skipped, so an interrupted run resumes where it stopped (`--force` converts it again). `-s` restricts the run to some species, `-j` converts several at once (each needs up to 4 GB), `--skip-download` reuses the tables already there, `--dry-run` reports what would be done. The exit status is 0 only when every species succeeded; failures are named in the summary and in the manifest, and the run carries on unless `--stop-on-error` is given.

`manifest.json`, written next to the species directories, says what was produced — the converter's version, the configuration, and for each species its taxon, the core database and the files with their size and triple count. Triple counts come from `rapper`, which already parses every file, so they cost nothing; `statements_written` is lower because one statement may carry several objects (`skos:altLabel "a", "b"`).

## Running on a cluster
[bin/convert_slurm.sh](bin/convert_slurm.sh) converts one species per array task — download, conversion, `rapper` and `gzip`:
```
$ mkdir -p logs
$ sbatch --array=1-$(python3 bin/species_config.py species ensembl_rdf/config/species_agrold.yaml | wc -l) \
         bin/convert_slurm.sh ensembl_rdf/config/species_agrold.yaml /path/to/workdir
```
Tasks are independent, so a failed one is resubmitted alone with `--array=<n>`. Edit the `#SBATCH` header and the commented `module load` / `conda activate` lines for your site.

Sizing, measured on Ensembl Plants 63: the 51 species of `species_agrold.yaml` amount to 1.3 GB of compressed dumps. The converter loads the tables it needs in RAM, roughly 35-40x their compressed size (Arabidopsis: 44 MB of dumps, 1.6 GB of RAM, 39 s, 560 MB of Turtle before `rapper` and `gzip`, 28 MB after); 16 GB per task covers every plant species, bread wheat included.

## Cross-reference sources
Cross-references are turned into `rdfs:seeAlso` links using [ensembl_rdf/config/external_db_url.tsv](ensembl_rdf/config/external_db_url.tsv): one line per `external_db_id` (numeric, not the name), giving the URL the accession is appended to and, optionally, a prefix to strip from it. Sources with no URL are reported in `xref_report.tsv` next to each converted database, and no link is emitted for them.

Before converting a whole division, survey which sources its species actually use — it only needs `xref.txt.gz` and `external_db.txt.gz`:
```
$ python3 .../bin/download_files.py URL --dbinfo /path/to/ensembl_rdf/ensembl_rdf/config/xref_dbinfo.json
$ python3 .../bin/xref_survey.py                  # every *_core_* of the cwd
$ python3 .../bin/xref_survey.py --missing-only   # only the sources with no URL
$ python3 .../bin/xref_survey.py --missing-only --tsv   # ids to append to the table
```
The report gives, per source, the number of species using it, the number of xrefs, an example accession and its current URL, so the sources worth mapping come first. Ids already listed with an empty URL are filled in place; ids the table does not know are printed by `--tsv`.

## RDF schema
The model is described in the [rdf-config](https://github.com/dbcls/rdf-config) format, one entity per Turtle file with an example value for each property:

- [rdf-config/model.yaml](rdf-config/model.yaml) — the Ensembl model, as the `ensembl` profile emits it (upstream version: [dbcls/rdf-config](https://github.com/dbcls/rdf-config/blob/master/config/ensembl/model.yaml))
- [rdf-config/agrold/model.yaml](rdf-config/agrold/model.yaml) — the same model as the `agrold` profile emits it: AgroLD classes, Ensembl properties, `purl.agrold.org` URIs and named FALDO regions

Each comes with its `prefix.yaml`. Update them along with the converter: they are the reference for anyone querying the result.

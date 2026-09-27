import gzip
import sys
import argparse
import os
import json
import re
import datetime
import urllib.parse

input_dir = "./"


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


def quote(string):
    return "\"" + string.replace("\"", "\\\"") + "\""


def escape(string):
    """Escape a stable ID for use as the local part of a Turtle prefixed name.

    Characters of PN_LOCAL_ESC (Turtle grammar) are backslash-escaped, e.g. the
    apostrophe of plastid exon IDs like `transcript-3'rps12-E1`. A trailing dot
    is escaped too, as it would otherwise end the statement.
    """
    string = re.sub(r"([~!$&'()*+,;=/?#@%])", r"\\\1", string)
    return re.sub(r"\.$", r"\\.", string)


def iri_escape(string):
    """Percent-encode characters not allowed inside an IRI written as <...>."""
    return re.sub(r'[<>"{}|^`\\\s]', lambda m: "%%%02X" % ord(m.group()), string)


def load_model(path):
    """Load a vocabulary profile (config/models/*.yaml): prefixes, terms, options."""
    try:
        import yaml
    except ImportError:
        sys.exit("Error: PyYAML is required to read the model profile "
                 f"{path} (pip install -r requirements.txt)")
    with open(path, "r") as f:
        conf = yaml.safe_load(f) or {}
    vocabulary = conf.get("vocabulary") or {}
    resources = conf.get("resources") or {}
    missing_res = [k for k in Ensembl2turtle.resource_keys if k not in resources]
    if missing_res:
        sys.exit(f"Error: model {path} does not define resources {missing_res}")
    model = {
        "resources": resources,
        "prefixes": conf.get("prefixes") or {},
        "terms": conf.get("terms") or {},
        "options": conf.get("options") or {},
        "vocabulary_prefix": vocabulary.get("prefix", "terms"),
        "vocabulary_uri": vocabulary.get("uri"),
    }
    missing = [k for k in Ensembl2turtle.model_keys if k not in model["terms"]]
    if missing:
        sys.exit(f"Error: model {path} does not define {missing}")
    return model


def strand2faldo(s):
    if s == "1":
        return "faldo:ForwardStrandPosition"
    elif s == "-1":
        return "faldo:ReverseStrandPosition"
    else:
        print(f"Error: Invalid argument \"{s}\" for strand2faldo", file=sys.stderr)
        sys.exit(1)


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


class Ensembl2turtle:
    # Base of the resource URIs (genes, transcripts, proteins, exons, chromosomes)
    # and namespace of the terms:/ontology vocabulary. Both can be overridden.
    default_base_uri = "http://rdf.ebi.ac.uk"
    default_terms_uri = "http://rdf.ebi.ac.uk/terms/ensembl/"

    common_prefixes = [
        ['rdf:', '<http://www.w3.org/1999/02/22-rdf-syntax-ns#>'],
        ['rdfs:', '<http://www.w3.org/2000/01/rdf-schema#>'],
        ['faldo:', '<http://biohackathon.org/resource/faldo#>'],
        ['obo:', '<http://purl.obolibrary.org/obo/>'],
        ['so:', '<http://purl.obolibrary.org/obo/so#>'],
        ['dc:', '<http://purl.org/dc/elements/1.1/>'],
        ['dcterms:', '<http://purl.org/dc/terms/>'],
        ['owl:', '<http://www.w3.org/2002/07/owl#>'],
        ['ensgloss:', '<http://ensembl.org/glossary/>'],
        ['ensi:', '<http://identifiers.org/ensembl/>'],
        ['taxonomy:', '<http://identifiers.org/taxonomy/>'],
        ['uniprot:', '<http://purl.uniprot.org/uniprot/>'],
        ['refseq:', '<http://identifiers.org/refseq/>'],
        ['sio:', '<http://semanticscience.org/resource/>'],
        ['skos:', '<http://www.w3.org/2004/02/skos/core#>']
    ]

    # Keys are `code` of the attrib_type table to be used as transcript flags
    transcript_flags = {
        "gencode_basic": {
            "GENCODE basic": "ensgloss:ENSGLOSSARY_0000020",
            "1": "ensgloss:ENSGLOSSARY_0000020"
        },
        "appris": {
            "principal1": "ensgloss:ENSGLOSSARY_0000013",
            "principal2": "ensgloss:ENSGLOSSARY_0000014",
            "principal3": "ensgloss:ENSGLOSSARY_0000015",
            "principal4": "ensgloss:ENSGLOSSARY_0000016",
            "principal5": "ensgloss:ENSGLOSSARY_0000017",
            "alternative1": "ensgloss:ENSGLOSSARY_0000018",
            "alternative2": "ensgloss:ENSGLOSSARY_0000019"
        },
        "TSL": {
            "tsl1": "ensgloss:ENSGLOSSARY_0000006",
            "tsl2": "ensgloss:ENSGLOSSARY_0000007",
            "tsl3": "ensgloss:ENSGLOSSARY_0000008",
            "tsl4": "ensgloss:ENSGLOSSARY_0000009",
            "tsl5": "ensgloss:ENSGLOSSARY_0000010",
            "tslNA": "ensgloss:ENSGLOSSARY_0000011"
        },
        "is_canonical": {"1": "ensgloss:ENSGLOSSARY_0000023"},
        # "remark": {"MANE_select": "ensgloss:ENSGLOSSARY_0000365"},
        "MANE_Select": {},
        "MANE_Plus_Clinical": {},
        "mRNA_start_NF": {"1": "ensgloss:ENSGLOSSARY_0000021"},
        "mRNA_end_NF": {"1": "ensgloss:ENSGLOSSARY_0000022"},
        "cds_start_NF": {"1": "ensgloss:ENSGLOSSARY_0000021"},
        "cds_end_NF": {"1": "ensgloss:ENSGLOSSARY_0000022"}
    }

    # `code` of attrib_type entries that appear on transcripts but are not flags.
    # Listing them keeps the warning below meaningful: an attribute that is not
    # here has not been looked at yet, and may well deserve to be converted.
    non_flag_attributes = {
        # sequence checksums
        "md5_cdna", "md5_cds", "sha512t24u_cdna", "sha512t24u_cds",
        # secondary structure annotations
        "miRNA",   # Micro RNA
        "ncRNA",   # Structure
    }

    hco_chr_names = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10",
                     "11", "12", "13", "14", "15", "16", "17", "18",
                     "19", "20", "21", "22", "X", "Y", "MT"]

    base_dir = os.path.dirname(os.path.abspath(__file__)) + "/../"
    models_dir = base_dir + "config/models/"

    # Entity types that can be output (one Turtle file each), in output order,
    # with the tables each one needs. `meta`, `seq_region` and `coord_system`
    # are always loaded.
    entity_tables = {
        "gene": ["gene", "xref", "external_synonym"],
        "transcript": ["transcript", "transcript_attrib", "attrib_type", "xref", "gene", "translation"],
        "translation": ["translation", "transcript"],
        "exon": ["exon"],
        "exon_transcript": ["exon_transcript", "transcript", "exon"],
        "xref": ["gene", "transcript", "translation", "xref", "object_xref", "external_db"],
        "chromosome": ["seq_region", "coord_system"]
    }
    entities = list(entity_tables.keys())

    # URIs of the resources themselves; `uri` may use {base}, set by --base-uri
    resource_keys = ["gene", "transcript", "protein", "exon", "chromosome"]

    # Terms whose object is a literal, not a resource (used to declare them
    # as owl:DatatypeProperty when deriving a profile's ontology)
    literal_terms = {"label", "description", "identifier", "alt_label", "has_version",
                     "in_assembly", "in_schema_number", "ordered_exon_rank"}

    # Every element of the model a profile must map (config/models/*.yaml)
    model_keys = [
        "gene_class", "transcript_class", "protein_class", "exon_class", "chromosome_class",
        "exon_so_class", "ordered_exon_class", "ordered_list_item_class",
        "versioned_transcript_class",
        "label", "description", "identifier", "alt_label", "see_also",
        "in_taxon", "location",
        "has_biotype", "part_of", "transcribed_from", "translates_to",
        "translation_of", "has_exon", "has_ordered_exon",
        "ordered_exon_refers_to", "ordered_exon_rank",
        "has_transcript_flag", "has_versioned_transcript", "has_version",
        "has_counterpart",
    ]

    def __init__(self, input_dbinfo_file, entities=None, base_uri=None, terms_uri=None,
                 model=None):
        self.entities = [e for e in Ensembl2turtle.entities if entities is None or e in entities]
        self.base_uri = (base_uri or Ensembl2turtle.default_base_uri).rstrip("/")
        self.terms_uri = terms_uri
        self.model = model or {"prefixes": {}, "terms": {}, "options": {}, "resources": {},
                               "vocabulary_prefix": "terms", "vocabulary_uri": None}
        self.t = self.model["terms"]
        self.opt = self.model["options"]
        # The vocabulary namespace holds the classes and properties of the model
        # itself (terms: in the Ensembl model). A profile may rename it and move
        # it; --terms-uri still wins.
        self.vocabulary_prefix = self.model.get("vocabulary_prefix") or "terms"
        if self.terms_uri is None:
            self.terms_uri = self.model.get("vocabulary_uri") or Ensembl2turtle.default_terms_uri
        # Resource namespaces: {prefix, uri} per resource type, from the profile
        self.res = {}
        for key in Ensembl2turtle.resource_keys:
            conf = self.model["resources"].get(key) or {}
            self.res[key] = {
                "prefix": conf.get("prefix", ""),
                "uri": conf.get("uri", "{base}/resource/").replace("{base}", self.base_uri),
            }
        self.prefixes = Ensembl2turtle.common_prefixes + [
            [self.vocabulary_prefix + ':', '<' + self.terms_uri + '>'],
        ] + [[r["prefix"] + ':', '<' + r["uri"] + '>']
             for key, r in self.res.items() if r["prefix"]
        ] + [[k + ':', '<' + v + '>'] for k, v in self.model["prefixes"].items()]
        self.dbinfo = self.load_dbinfo(input_dbinfo_file)
        self.dbs = self.load_dbs()
        # self.taxonomy_id = self.get_taxonomy_id()
        self.ensembl_version = self.get_ensembl_version()
        # self.production_name = self.get_production_name()
        self.species_id2taxonomy_id = self.get_species_id2taxonomy_id()
        # print(self.species_id2taxonomy_id)
        self.species_id2production_name = self.get_species_id2production_name()
        self.xref_url_dic = {}
        self.xref_prefix_dic = {}
        self.init_xref_url_dic()
        self.xrefed_dbs = {"Gene": {}, "Transcript": {}, "Translation": {}}
        self.not_xrefed_dbs = {"Gene": {}, "Transcript": {}, "Translation": {}}
        self.output_file = sys.stdout
        self.biotype_url_dic = {}
        self.init_biotype_url_dic()
        # Regions the converted features sit on; `chromosome` describes those
        self.referenced_seq_regions = set()
        self.statement_count = 0
        self.statements_per_file = {}

    def init_biotype_url_dic(self):
        biotype_url_dic_tsv = "ontology/biotype_url.tsv"
        with open(Ensembl2turtle.base_dir+biotype_url_dic_tsv, "r") as input_table:
            line = input_table.readline()
            while (line):
                line = line.rstrip('\n')
                sep_line = line.split('\t')
                if sep_line[1] != "":
                    self.biotype_url_dic[sep_line[0]] = re.sub(
                        r"^terms:", self.vocabulary_prefix + ":", sep_line[1])
                line = input_table.readline()
        return

    def init_xref_url_dic(self):
        xref_url_dic_tsv = "config/external_db_url.tsv"

        with open(Ensembl2turtle.base_dir+xref_url_dic_tsv, "r") as input_table:
            line = input_table.readline()
            while (line):
                line = line.rstrip('\n')
                sep_line = line.split('\t')
                if sep_line[1] != "":
                    self.xref_url_dic[sep_line[0]] = sep_line[1]
                if sep_line[2] != "":
                    self.xref_prefix_dic[sep_line[0]] = sep_line[2]
                line = input_table.readline()
        return

    def triple(self, s, p, o):
        print(s, p, o, ".", file=self.output_file)
        self.statement_count += 1
        return

    def get_ensembl_version(self):
        ensembl_version = [v[2] for k, v in self.dbs["meta"].items() if v[1] == 'schema_version']
        return ensembl_version[0]

    def get_taxonomy_id(self):
        taxonomy_ids = [v[2] for k, v in self.dbs["meta"].items() if v[1] == 'species.taxonomy_id']
        if len(taxonomy_ids) >= 2:
            print("Error: `meta` table has multiple taxonomy_id. This seems to be multi-species database.", file=sys.stderr)
            print("taxonomy_ids: ", taxonomy_ids, file=sys.stderr)
            sys.exit(1)

        return taxonomy_ids[0]

    def get_production_name(self):
        production_names = [v[2] for k, v in self.dbs["meta"].items() if v[1] == 'species.production_name']
        return production_names[0]

    def get_species_id2production_name(self):
        return {v[0]: v[2] for v in self.dbs["meta"].values() if v[1] == 'species.production_name'}

    def get_species_id2taxonomy_id(self):
        return {v[0]: v[2] for v in self.dbs["meta"].values() if v[1] == 'species.taxonomy_id'}

    def load_dbinfo(self, input_dbinfo_file):
        dbinfo_dict = {}
        with open(input_dbinfo_file, 'r') as input_dbinfo:
            dbinfo_dict = json.load(input_dbinfo)
        return dbinfo_dict

    def load_db(self, db):
        dt_now = datetime.datetime.now()
        print(f"[{dt_now}] Loading DB: {db}", file=sys.stderr)
        dic = {}
        table_file = self.dbinfo[db]["filename"]
        key_indices = self.dbinfo[db]["key_indices"]
        val_indices = [v["index"] for v in self.dbinfo[db]["values"]]
        with gzip.open(table_file, 'rt') as input_table:
            line = input_table.readline()
            while (line):
                line = line.rstrip('\n')
                sep_line = line.split('\t')
                # e.g. transcript_attrib: [transcript_id, attrib_type_id, value]
                key_list = [sep_line[i] for i in key_indices]
                if len(key_list) >= 2:
                    key = tuple(key_list)
                else:
                    key = key_list[0]
                vals = [sep_line[i] for i in val_indices]
                if self.dbinfo[db].get("list", False):
                    if key in dic:
                        dic[key].append(vals)
                    else:
                        dic[key] = [vals]
                else:
                    dic[key] = vals

                line = input_table.readline()

        return dic

    def load_dbs(self):
        needed = {"meta", "seq_region", "coord_system"}
        for entity in self.entities:
            needed.update(Ensembl2turtle.entity_tables[entity])
        db_dics = {}
        for db in self.dbinfo:
            if db not in needed:
                continue
            db_dics[db] = self.load_db(db)
            memory = memory_usage_mb()
            if memory is not None:
                print(f'memory: {memory:.1f} MB', file=sys.stderr)

        return db_dics

    def rdfize_gene(self):
        gene = self.dbs["gene"]
        xref = self.dbs["xref"]
        seq_region = self.dbs["seq_region"]
        external_synonym = self.dbs["external_synonym"]
        f = open("gene.ttl", mode="w")
        self.output_file = f
        self.output_prefixes()
        for id in gene:
            sbj = self.uri("gene", gene[id][6])
            xref_id = gene[id][4]
            seq_region_id = gene[id][7]

            self.triple(sbj, "a", self.t["gene_class"])
            self.output_biotype(sbj, gene[id][0])
            if xref_id == "\\N":
                label = gene[id][6]  # Substitute ID for label
            else:
                label = xref[xref_id][2]
            self.triple(sbj, self.t["label"], quote(label))
            description = gene[id][5]
            if description == "\\N":
                description = ""
            self.triple(sbj, self.t["description"], quote(description))
            self.triple(sbj, self.t["identifier"], quote(gene[id][6]))
            self.triple(sbj, self.t["in_taxon"], "taxonomy:"+self.seq_region_id_to_taxonomy_id(seq_region_id))

            # synonym
            synonyms = ", ".join([quote(v[0]) for v in external_synonym.get(xref_id, [])])
            if len(synonyms) >= 1:
                self.triple(sbj, self.t["alt_label"], synonyms)

            # location
            chromosome_urls = self.seq_region_id_to_chr(seq_region_id)
            self.output_location(sbj, gene[id][1], gene[id][2], gene[id][3], chromosome_urls)
            for chromosome_url in chromosome_urls:
                self.triple(sbj, self.t["part_of"], "<" + chromosome_url + ">")
        self.output_file = sys.stdout
        f.close()
        return

    def rdfize_transcript(self):
        transcript = self.dbs["transcript"]
        transcript_attrib = self.dbs["transcript_attrib"]
        xref = self.dbs["xref"]
        gene = self.dbs["gene"]
        translation = self.dbs["translation"]
        attrib_type = self.dbs["attrib_type"]
        unknown_flags = set()
        f = open("transcript.ttl", mode="w")
        self.output_file = f
        self.output_prefixes()
        for id in transcript:
            stable_id = transcript[id][7]
            sbj = self.uri("transcript", stable_id)
            xref_id = transcript[id][4]

            self.triple(sbj, "a", self.t["transcript_class"])
            self.output_biotype(sbj, transcript[id][5])
            if xref_id == "\\N":
                label = stable_id  # Substitute ID for label
            else:
                label = xref[xref_id][2]
            self.triple(sbj, self.t["label"], quote(label))
            self.triple(sbj, self.t["identifier"], quote(stable_id))
            self.triple(sbj, self.t["transcribed_from"], self.uri("gene", gene[transcript[id][0]][6]))
            translates_to = transcript[id][6]
            if translates_to != "\\N":
                self.triple(sbj, self.t["translates_to"], self.uri("protein", translation[translates_to][1]))

            # location
            chromosome_urls = self.seq_region_id_to_chr(transcript[id][8])
            self.output_location(sbj, transcript[id][1], transcript[id][2], transcript[id][3],
                                 chromosome_urls)

            # flag
            attribs = transcript_attrib.get(id, [])
            flag_dic = Ensembl2turtle.transcript_flags
            for attrib in attribs:
                attrib_code = attrib_type[attrib[0]][0]
                if attrib_code in flag_dic:
                    attrib_val = attrib[1]
                    if attrib_code == "TSL":
                        comment = ""
                        match = re.search(r'\((.*?)\)', attrib_val)
                        attrib_val = re.sub(r" .*", "", attrib[1])
                        if match:
                            comment = match.group(1)
                            statement = "<" + self.res["transcript"]["uri"] + "#_" + iri_escape(stable_id) + "-has_transcript_flag-"+attrib_val+">"
                            self.triple(statement, "a", "rdf:Statement")
                            self.triple(statement, "rdf:subject", sbj)
                            self.triple(statement, "rdf:predicate", self.t["has_transcript_flag"])
                            self.triple(statement, "rdf:object", flag_dic[attrib_code][attrib_val])
                            self.triple(statement, "rdfs:comment", quote(comment))
                    # elif attrib_code == "remark":
                    #     if attrib_val != "MANE_select":
                        #     continue
                        # else:
                    ## For the MANE transcripts, triples with a versioned transcript is added.
                    elif attrib_code == "MANE_Select" or attrib_code == "MANE_Plus_Clinical":
                        if attrib_code == "MANE_Select":
                            ensgloss_term = "ensgloss:ENSGLOSSARY_0000365"
                        else:
                            ensgloss_term = "ensgloss:ENSGLOSSARY_0000375"
                        version = transcript[id][9]
                        versioned_id = escape(stable_id) + "." + version
                        versioned_sbj = self.uri("transcript", stable_id + "." + version)
                        self.triple(sbj, self.t["has_transcript_flag"], ensgloss_term)
                        self.triple(sbj, self.t["has_versioned_transcript"], versioned_sbj)
                        self.triple(versioned_sbj, "a", self.t["versioned_transcript_class"])
                        self.triple(versioned_sbj, self.t["has_version"], version)
                        self.triple(versioned_sbj, self.t["identifier"], quote(versioned_id))
                        self.triple(versioned_sbj, self.t["has_transcript_flag"], ensgloss_term)
                        counterpart = "refseq:" + attrib_val
                        self.triple(versioned_sbj, self.t["has_counterpart"], counterpart)
                        self.triple(re.sub(r"\.[0-9]+$", "", counterpart), self.t["has_versioned_transcript"], counterpart)
                        continue
                    # self.triple(sbj, "terms:has_transcript_flag", flag_dic[attrib_code][attrib_val])
                    try:
                        self.triple(sbj, self.t["has_transcript_flag"], flag_dic[attrib_code][attrib_val])
                    except KeyError as e:
                        print(f"Warning: KeyError: {e}; {sbj} {attrib_code} {attrib_val}", file=sys.stderr)
                elif attrib_code not in Ensembl2turtle.non_flag_attributes:
                    if attrib[0] not in unknown_flags:
                        print(f"Warning: Attribute not treated as flag: {attrib[0]} {attrib_code}", file=sys.stderr)
                        unknown_flags.add(attrib[0])
        self.output_file = sys.stdout
        f.close()
        return

    def uri(self, kind, stable_id):
        """URI of a resource, as a prefixed name when the profile declares a prefix."""
        prefix = self.res[kind]["prefix"]
        if prefix:
            return prefix + ":" + escape(stable_id)
        return "<" + self.res[kind]["uri"] + iri_escape(stable_id) + ">"

    def output_biotype(self, sbj, biotype):
        """Emit the biotype, as an ontology term or as a literal (see the model profile)."""
        if self.opt.get("biotype_as_literal"):
            self.triple(sbj, self.t["has_biotype"], quote(biotype))
            return
        if biotype not in self.biotype_url_dic:
            print(f'Warning: Unknown biotype `{biotype}`', file=sys.stderr)
            return
        if self.opt.get("biotype_as_type", True):
            self.triple(sbj, "a", self.biotype_url_dic[biotype])
        self.triple(sbj, self.t["has_biotype"], self.biotype_url_dic[biotype])

    def seq_region_id_to_taxonomy_id(self, seq_region_id):
        seq_region = self.dbs["seq_region"]
        coord_system = self.dbs["coord_system"]
        coord_system_id = seq_region[seq_region_id][1]
        species_id = coord_system[coord_system_id][0]
        return self.species_id2taxonomy_id[species_id]

    def seq_region_id_to_production_name(self, seq_region_id):
        seq_region = self.dbs["seq_region"]
        coord_system = self.dbs["coord_system"]
        coord_system_id = seq_region[seq_region_id][1]
        species_id = coord_system[coord_system_id][0]
        return self.species_id2production_name[species_id]

    def seq_region_id_to_chr(self, seq_region_id):
        """Bare URIs of the region a feature sits on (several for human, with HCO).

        The shape comes from the model profile, e.g.
        `{version}/{production_name}/{assembly}/{chromosome}` for Ensembl or
        `{taxon}/{assembly}/{chromosome}` for AgroLD.
        """
        seq_region = self.dbs["seq_region"]
        coord_system = self.dbs["coord_system"]
        chromosome_name = seq_region[seq_region_id][0]
        coord_system_id = seq_region[seq_region_id][1]
        taxonomy_id = self.seq_region_id_to_taxonomy_id(seq_region_id)
        fields = {
            "version": self.ensembl_version,
            "production_name": self.seq_region_id_to_production_name(seq_region_id),
            "assembly": coord_system[coord_system_id][2],  # e.g. "GRCm38", "TAIR10"
            "chromosome": chromosome_name,
            "taxon": taxonomy_id,
        }
        # LRG regions have no assembly
        pattern = self.opt.get("chromosome_pattern", "{version}/{production_name}/{assembly}/{chromosome}")
        if coord_system[coord_system_id][1] == "lrg":
            pattern = self.opt.get("chromosome_lrg_pattern", "{version}/{production_name}/{chromosome}")
        chromosome_urls = [self.res["chromosome"]["uri"] + pattern.format(**fields)]

        if taxonomy_id == "9606" and chromosome_name in Ensembl2turtle.hco_chr_names:
            chromosome_urls.append("http://identifiers.org/hco/" + chromosome_name + "/" + fields["assembly"])

        self.referenced_seq_regions.add(seq_region_id)
        return chromosome_urls

    def create_location_str(self, beg, end, strand, chromosome_urls, level=1):
        """Blank node describing the location, as in the Ensembl RDF."""
        loc = Bnode()
        loc_beg = Bnode()
        loc_end = Bnode()

        loc_beg.add(("a", "faldo:ExactPosition"))
        loc_beg.add(("a", strand2faldo(strand)))
        loc_beg.add(("faldo:position", beg))
        loc_end.add(("a", "faldo:ExactPosition"))
        loc_end.add(("a", strand2faldo(strand)))
        loc_end.add(("faldo:position", end))

        for chromosome_url in chromosome_urls:
            loc_beg.add(("faldo:reference", "<" + chromosome_url + ">"))
            loc_end.add(("faldo:reference", "<" + chromosome_url + ">"))

        loc.add(("a", "faldo:Region"))
        loc.add(("faldo:begin", loc_beg.serialize(level=level+1)))
        loc.add(("faldo:end", loc_end.serialize(level=level+1)))
        return loc.serialize()

    def output_location(self, sbj, beg, end, strand, chromosome_urls):
        """Attach a location to `sbj`.

        With `faldo_named_regions`, the region and its two positions are named
        resources shaped after the region they sit on
        (<chromosome>:<begin>-<end>:<strand>), so that features sharing
        coordinates share them; otherwise they are blank nodes, as in the
        Ensembl RDF.
        """
        if not self.opt.get("faldo_named_regions"):
            self.triple(sbj, self.t["location"], self.create_location_str(beg, end, strand, chromosome_urls))
            return

        chromosome_url = chromosome_urls[0]
        region = "<" + chromosome_url + ":" + beg + "-" + end + ":" + strand + ">"
        self.triple(sbj, self.t["location"], region)
        begin_uri = "<" + chromosome_url + ":" + beg + ":" + strand + ">"
        end_uri = "<" + chromosome_url + ":" + end + ":" + strand + ">"
        self.triple(region, "a", "faldo:Region")
        self.triple(region, "faldo:begin", begin_uri)
        self.triple(region, "faldo:end", end_uri)
        for position_uri, position in ((begin_uri, beg), (end_uri, end)):
            self.triple(position_uri, "a", "faldo:ExactPosition")
            self.triple(position_uri, "a", strand2faldo(strand))
            self.triple(position_uri, "faldo:position", position)
            for url in chromosome_urls:
                self.triple(position_uri, "faldo:reference", "<" + url + ">")

    def rdfize_translation(self):
        transcript = self.dbs["transcript"]
        translation = self.dbs["translation"]
        f = open("translation.ttl", mode="w")
        self.output_file = f
        self.output_prefixes()
        for id in translation:
            sbj = self.uri("protein", translation[id][1])

            self.triple(sbj, "a", self.t["protein_class"])
            self.triple(sbj, self.t["identifier"], quote(translation[id][1]))
            self.triple(sbj, self.t["translation_of"], self.uri("transcript", transcript[translation[id][0]][7]))
        self.output_file = sys.stdout
        f.close()
        return

    def rdfize_exon(self):
        exon = self.dbs["exon"]
        f = open("exon.ttl", mode="w")
        self.output_file = f
        self.output_prefixes()
        for id in exon:
            sbj = self.uri("exon", exon[id][3])

            self.triple(sbj, "a", self.t["exon_class"])
            self.triple(sbj, "a", self.t["exon_so_class"])
            self.triple(sbj, self.t["identifier"], quote(exon[id][3]))

            # location
            chromosome_urls = self.seq_region_id_to_chr(exon[id][4])
            self.output_location(sbj, exon[id][0], exon[id][1], exon[id][2], chromosome_urls)
        self.output_file = sys.stdout
        f.close()
        return

    def rdfize_exon_transcript(self):
        exon_transcript = self.dbs["exon_transcript"]
        transcript = self.dbs["transcript"]
        exon = self.dbs["exon"]
        f = open("exon_transcript.ttl", mode="w")
        self.output_file = f
        self.output_prefixes()
        for id in exon_transcript:
            exon_id = id[0]
            transcript_id = id[1]
            exon_stable_id = exon[exon_id][3]
            transcript_stable_id = transcript[transcript_id][7]
            rank = exon_transcript[id][0]
            ordered_exon_uri = "<"+self.res["transcript"]["uri"]+iri_escape(transcript_stable_id)+"#Exon_"+rank+">"
            exon_uri = self.uri("exon", exon_stable_id)
            transcript_uri = self.uri("transcript", transcript_stable_id)

            self.triple(ordered_exon_uri, "a", self.t["ordered_exon_class"])
            self.triple(ordered_exon_uri, "a", self.t["ordered_list_item_class"])
            self.triple(ordered_exon_uri, self.t["ordered_exon_refers_to"], exon_uri)
            self.triple(ordered_exon_uri, self.t["ordered_exon_rank"], rank)

            self.triple(transcript_uri, self.t["has_exon"], exon_uri)
            self.triple(transcript_uri, self.t["has_ordered_exon"], ordered_exon_uri)
        self.output_file = sys.stdout
        f.close()
        return

    def rdfize_chromosome(self):
        """Describe the regions the converted features sit on.

        Only those actually referenced are described, so the thousands of
        contigs and scaffolds of a seq_region table are left out. When no other
        entity ran, every region of a chromosome coordinate system is taken.
        """
        seq_region = self.dbs["seq_region"]
        coord_system = self.dbs["coord_system"]
        ids = self.referenced_seq_regions
        if not ids:
            ids = {i for i in seq_region if coord_system[seq_region[i][1]][1] == "chromosome"}
        f = open("chromosome.ttl", mode="w")
        self.output_file = f
        self.output_prefixes()
        for seq_region_id in sorted(ids):
            coord_system_id = seq_region[seq_region_id][1]
            name = seq_region[seq_region_id][0]
            sbj = "<" + self.seq_region_id_to_chr(seq_region_id)[0] + ">"

            if coord_system[coord_system_id][1] == "chromosome":
                self.triple(sbj, "a", self.t["chromosome_class"])
            if self.t.get("chromosome_so_class"):
                self.triple(sbj, "a", self.t["chromosome_so_class"])
            self.triple(sbj, self.t["label"], quote(name))
            self.triple(sbj, self.t["identifier"], quote(name))
            self.triple(sbj, self.t["in_taxon"],
                        "taxonomy:" + self.seq_region_id_to_taxonomy_id(seq_region_id))
            if self.t.get("in_assembly"):
                self.triple(sbj, self.t["in_assembly"], quote(coord_system[coord_system_id][2]))
            if self.t.get("in_schema_number"):
                self.triple(sbj, self.t["in_schema_number"], quote(self.ensembl_version))
        self.output_file = sys.stdout
        f.close()
        return

    def rdfize_xref(self):
        gene = self.dbs["gene"]
        transcript = self.dbs["transcript"]
        translation = self.dbs["translation"]
        xref = self.dbs["xref"]
        object_xref = self.dbs["object_xref"]
        external_db = self.dbs["external_db"]
        f = open("xref.ttl", mode="w")
        self.output_file = f
        self.output_prefixes()
        for id in object_xref:
            xref_id = object_xref[id][2]
            subject_id = object_xref[id][0]
            subject_type = object_xref[id][1]
            if subject_type == "Gene":
                subject_url = self.uri("gene", gene[subject_id][6])
            elif subject_type == "Transcript":
                subject_url = self.uri("transcript", transcript[subject_id][7])
            elif subject_type == "Translation":
                subject_url = self.uri("protein", translation[subject_id][1])
            else:
                continue
            # xref_node = Bnode()
            # xref_node.add(("terms:id_of", quote(external_db[xref[xref_id][0]][0])))  # FIXME
            # xref_node.add(("dcterms:identifier", quote(xref[xref_id][1])))
            # self.triple(subject_url, "rdfs:seeAlso", xref_node.serialize())
            external_db_id = xref[xref_id][0]
            external_db_code = external_db[external_db_id][0]
            if self.xref_url_dic.get(external_db_id, "") != "":
                dbprimary_acc = xref[xref_id][1]
                if external_db_id in self.xref_prefix_dic:
                    dbprimary_acc = dbprimary_acc.replace(self.xref_prefix_dic[external_db_id], "")
                xref_url = self.xref_url_dic[external_db_id] + urllib.parse.quote(dbprimary_acc, safe=":/")
                self.triple(subject_url, self.t["see_also"], "<"+xref_url+">")
                if external_db_code not in self.xrefed_dbs[subject_type]:
                    self.xrefed_dbs[subject_type][external_db_code] = [subject_url, xref_url, 0]
                self.xrefed_dbs[subject_type][external_db_code][2] += 1
            else:
                if external_db_code not in self.not_xrefed_dbs[subject_type]:
                    self.not_xrefed_dbs[subject_type][external_db_code] = [subject_url, xref[xref_id][1], 0]
                self.not_xrefed_dbs[subject_type][external_db_code][2] += 1
        self.output_file = sys.stdout
        f.close()
        return

    def output_prefixes(self):
        for prefix in self.prefixes:
            self.triple("@prefix", prefix[0], prefix[1])
        self.statement_count -= len(self.prefixes)  # header lines, not statements
        return

    def output_turtle(self):
        rdfizers = {
            "gene": self.rdfize_gene,
            "transcript": self.rdfize_transcript,
            "translation": self.rdfize_translation,
            "exon": self.rdfize_exon,
            "exon_transcript": self.rdfize_exon_transcript,
            "xref": self.rdfize_xref,
            "chromosome": self.rdfize_chromosome
        }
        for entity in self.entities:
            dt_now = datetime.datetime.now()
            print(f"[{dt_now}] Output turtle: {entity}", file=sys.stderr)
            before = self.statement_count
            rdfizers[entity]()
            self.statements_per_file[entity + ".ttl"] = self.statement_count - before

        if "xref" in self.entities:
            self.output_xref_report()

        self.output_conversion_report()
        dt_now = datetime.datetime.now()
        print(f"[{dt_now}] Done.", file=sys.stderr)

    def output_conversion_report(self):
        """Describe this conversion for the run manifest (see bin/run.py)."""
        report = {
            "species": sorted(self.species_id2production_name.values()),
            "taxonomy_ids": sorted(self.species_id2taxonomy_id.values()),
            "ensembl_version": self.ensembl_version,
            "base_uri": self.base_uri,
            "vocabulary_uri": self.terms_uri,
            "entities": self.entities,
            # Statements written per file. A statement may carry several objects
            # (`skos:altLabel "a", "b"`), so rapper counts more triples than this.
            "statements": self.statements_per_file,
            "converted_at": datetime.datetime.now().isoformat(timespec="seconds"),
        }
        with open("conversion.json", "w") as f:
            json.dump(report, f, indent=2, sort_keys=True)

    def output_xref_report(self):
        with open("xref_report.tsv", "w") as f:
            cwd = os.getcwd()
            dir_prod_name = re.sub(r"(.*/)|(_core_[^/]+$)", "", cwd)
            for subject_type, dbs in self.xrefed_dbs.items():
                for db in dbs:
                    print(dir_prod_name, "xref", subject_type, db,
                          dbs[db][0], dbs[db][1], dbs[db][2], sep="\t", file=f)
            for subject_type, dbs in self.not_xrefed_dbs.items():
                for db in dbs:
                    print(dir_prod_name, "unknown", subject_type, db,
                          dbs[db][0], dbs[db][1], dbs[db][2], sep="\t", file=f)


def select_entities(include=None, exclude=None):
    """Return the entity types to process, in output order.

    `include` restricts to these types (default: all), `exclude` removes some.
    Both accept names separated by commas or spaces. Unknown names are an error.
    """
    def parse(names):
        result = []
        for n in names or []:
            result += [x for x in re.split(r"[,\s]+", n) if x]
        unknown = [x for x in result if x not in Ensembl2turtle.entities]
        if unknown:
            sys.exit(f"Error: unknown entity type(s) {unknown}; choose among {Ensembl2turtle.entities}")
        return result

    selected = parse(include) or list(Ensembl2turtle.entities)
    excluded = parse(exclude)
    result = [e for e in Ensembl2turtle.entities if e in selected and e not in excluded]
    # exon_transcript only links transcripts to exons, so it is meaningless without them
    if "exon_transcript" in result and "exon" not in result:
        print("Warning: exon is not selected, exon_transcript is skipped too", file=sys.stderr)
        result.remove("exon_transcript")
    return result


def main():
    parser = argparse.ArgumentParser(
        description="Convert the Ensembl core MySQL dumps of the current directory to Turtle files "
                    "(one per entity type: " + ", ".join(Ensembl2turtle.entities) + ").")
    parser.add_argument("dbinfo", help="config/dbinfo.json")
    parser.add_argument("-e", "--entities", nargs="+", metavar="ENTITY",
                        action="extend", default=[],
                        help="only output these entity types (default: all); repeatable")
    parser.add_argument("-x", "--exclude", nargs="+", metavar="ENTITY",
                        action="extend", default=[],
                        help="do not output these entity types (e.g. -x exon); repeatable")
    parser.add_argument("-b", "--base-uri", metavar="URI", default=Ensembl2turtle.default_base_uri,
                        help="base of the resource URIs, e.g. http://purl.agrold.org gives "
                             "http://purl.agrold.org/resource/ensembl/... (default: %(default)s)")
    parser.add_argument("-t", "--terms-uri", metavar="URI", default=None,
                        help="namespace of the model vocabulary, overriding the profile "
                             "(default: " + Ensembl2turtle.default_terms_uri + ")")
    parser.add_argument("-m", "--model", metavar="NAME_OR_FILE", default="ensembl",
                        help="vocabulary profile: a name in config/models/ (ensembl, agrold) "
                             "or a YAML file (default: %(default)s)")
    parser.add_argument("--list-entities", action="store_true",
                        help="print the selected entity types and exit (used by convert.sh)")
    args = parser.parse_args()

    entities = select_entities(args.entities, args.exclude)
    if args.list_entities:
        print(" ".join(entities))
        return
    print(f"Entities: {' '.join(entities)}", file=sys.stderr)
    model_path = args.model
    if not os.path.exists(model_path):
        model_path = Ensembl2turtle.models_dir + args.model + ".yaml"
        if not os.path.exists(model_path):
            available = sorted(f[:-5] for f in os.listdir(Ensembl2turtle.models_dir) if f.endswith(".yaml"))
            sys.exit(f"Error: unknown model '{args.model}'; choose among {available} or give a YAML file")
    model = load_model(model_path)

    print(f"Model: {model_path}", file=sys.stderr)
    print(f"Base URI: {args.base_uri}", file=sys.stderr)
    converter = Ensembl2turtle(args.dbinfo, entities, args.base_uri, args.terms_uri, model)
    converter.output_turtle()


if __name__ == '__main__':
    main()

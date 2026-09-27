import os
import sys
import json
import gzip
from utils import log_time, memory_usage_mb

class Ensembl2turtle:
    base_dir = os.path.dirname(os.path.abspath(__file__))

    def __init__(self, input_dbinfo_file, input_data_dir, needed_tables=None):
        self.input_data_dir = input_data_dir
        self.dbinfo = self.load_dbinfo(input_dbinfo_file)
        # Tables to read; the rest of dbinfo is left alone, which keeps a
        # conversion of part of the model out of the memory the whole one needs.
        self.needed_tables = needed_tables
        self.dbs = self.load_dbs()
        self.output_file = sys.stdout
        self.prefixes = [
            ["rdf:", "<http://www.w3.org/1999/02/22-rdf-syntax-ns#>"],
            ["rdfs:", "<http://www.w3.org/2000/01/rdf-schema#>"],
            ["dcterms:", "<http://purl.org/dc/terms/>"],
            ["obo:", "<http://purl.obolibrary.org/obo/>"],
            ["sio:", "<http://semanticscience.org/resource/>"],
            ["taxonomy:", "<http://identifiers.org/taxonomy/>"],
            ["so:", "<http://purl.obolibrary.org/obo/so#>"],
            ["dc:", "<http://purl.org/dc/elements/1.1/>"],
            ["owl:", "<http://www.w3.org/2002/07/owl#>"],
            ["ensg:", "<http://rdf.ebi.ac.uk/resource/ensembl/>"],
            ["terms:", "<http://rdf.ebi.ac.uk/terms/ensembl/>"],
            ["ensp:", "<http://rdf.ebi.ac.uk/resource/ensembl.protein/>"]
        ]

    def triple(self, s, p, o):
        print(s, p, o, ".", file=self.output_file)
        return

    def load_dbinfo(self, input_dbinfo_file):
        dbinfo_dict = {}
        with open(input_dbinfo_file, "r") as input_dbinfo:
            dbinfo_dict = json.load(input_dbinfo)
        return dbinfo_dict

    def load_db(self, db):
        log_time(f"Loading DB: {db}")
        dic = {}
        table_file = os.path.join(self.input_data_dir, self.dbinfo[db]["filename"])
        key_indices = self.dbinfo[db]["key_indices"]
        val_indices = [v["index"] for v in self.dbinfo[db]["values"]]
        with gzip.open(table_file, "rt") as input_table:
            line = input_table.readline()
            while (line):
                line = line.rstrip("\n")
                sep_line = line.split("\t")
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
        db_dics = {}
        for db in self.dbinfo:
            if self.needed_tables is not None and db not in self.needed_tables:
                continue
            db_dics[db] = self.load_db(db)
            memory = memory_usage_mb()
            if memory is not None:
                print(f"memory: {memory:.1f} MB", file=sys.stderr)

        return db_dics

    def output_prefixes(self):
        for prefix in self.prefixes:
            self.triple("@prefix", prefix[0], prefix[1])
        return


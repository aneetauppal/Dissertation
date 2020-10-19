#!/usr/bin/env python3

"""
Authors: Aaron Trautman, Aneeta Uppal, Steven Blanchard 
"""

import xml.etree.ElementTree as ET
import csv
import gzip
from zipfile import ZipFile
import logging
#from parse_it import handle

class Drugbank():
    """
    Parser class for parsing a drugbank xml file
    """

    def __init__(self, node_file, edge_file, db2uni, uni2ncbi, namespace="drugbank"):
        self.namespace = namespace
        self.conv_ID = generate_conversion_dict(db2uni,uni2ncbi)

        self.node_file = node_file
        self.node_file.write(
            "source_id:ID|name:string|url:string|synonyms:string[]|:LABEL\n"
            )
        
        self.edge_file = edge_file
        self.edge_file.write(
            ":START_ID|:TYPE|source:string|:END_ID\n"
        )

    def parse(self, xmlFile):
        """
        parsing function that parses an xml file
        """
        nodes = {}
        tree = ET.parse(xmlFile) # Memory intensive.. Could be optimized
        root = tree.getroot() # 
        for child in root: # 
            node = Node()
            node.parse(child)
            self.node_file.write(str(node))
            self.node_file.write("\n")
            self.__writeEdges(node)
    
    def __writeEdges(self, node):
        for d in node.mapsTo:
            new_id = d["outID"].strip("BE").lstrip("0")
            if new_id in self.conv_ID.keys():
                edge = "|".join([node.db_id, d["desc"],"drugbank",self.conv_ID[new_id]])
                self.edge_file.write(edge)
                self.edge_file.write("\n") 
        
class Node():

    __ns = "drugbank"
    __baseUrl = "https://www.drugbank.ca/drugs/{nID}"
    def __init__(self):
        self.db_id = None
        self.d_type = None
        self.name = None
        self.synonyms = set()
        self.mapsTo = []
        self.labels = set()
        self.labels.add(type(self).__ns)
    
    def parse(self, child):
        """
        Parse each xml child (a drug record)
        """
        xmlns = "{http://www.drugbank.ca}"
        self.d_type = child.attrib['type']
        self.labels.add("Drug")
        # Get primary ID
        for elem in child.findall(f"{xmlns}drugbank-id"):
            if 'primary' in elem.keys():
                self.db_id = elem.text
        # Get name
        self.name = child.find(f"{xmlns}name").text
        # Get Synonyms
        for elem in child.findall(f"{xmlns}synonyms"):
            for e in elem:
                if e.attrib["language"] == "english":
                    self.synonyms.add(e.text)
        # Get Drug MoA
        for elem in child.find(f"{xmlns}targets"):
            d2g = {}
            d2g["outID"] = elem.find(f"{xmlns}id").text
            d2g["desc"] = ""
            self.mapsTo.append(d2g)
     
    def __str__(self):
        """
            Redefining builtin str function to join records
        """
        return "|".join([
            self.db_id,
            self.name,
            self.__baseUrl.format(nID=self.db_id),
            ";".join(self.synonyms),
            ";".join(self.labels)
        ])

def generate_conversion_dict(db2uni,uni2ncbi):
    """
    Generates a conversion dictionary for converting drugbank
    protein IDs to ncbi gene IDs
    """
    ncbi_prefix = "nih.nlm.ncbi.gene.id:"
    d2u = {}
    u2n = {}
    d2n = {}
    # Drugbank to uniprot
    with handle(db2uni) as i1:
        reader = csv.DictReader(i1)
        for row in reader:
            d2u[row["ID"]] = row["UniProt ID"]
    # Uniprot to ncbi
    # Drugbank doesn't currently export their files with headers
    # If you have issues with ncbi gene ids, this might be somewhere to look.
    # Grab files for uniprot: human, rat, and mouse
    for f in uni2ncbi:
        with handle(f) as i2:
            lines = i2.readlines()
            for row in lines:
                u2n[row.split("\t")[0]] = ncbi_prefix+row.split("\t")[2]
    for k,v in d2u.items():
        if v in u2n.keys():
            d2n[k] = u2n[v]
        
    return d2n

def handle(filehandle):
    """
    File handling method for mix of files.
    """
    logging.info("handling {}".format(filehandle))
    if filehandle.endswith(".gz"):
        return gzip.open(filehandle, mode="rt")
    elif filehandle.endswith(".zip"):
        with ZipFile(filehandle) as myzip:
            records = myzip.infolist()
            if len(records) == 1:
                return myzip.open(records[0])
            logging.error("zip fail")
    return open(filehandle, "r")

if __name__ == "__main__":
    uniprot_base = "/Users/aneetauppal/Graduate_PhD/DissertationResearch/Knowledgebase/Drugbank/uniprot/{}.tab.gz"
    uniprot_files = [
        uniprot_base.format("HUMAN_9606"),
        uniprot_base.format("MOUSE_10090"),
        uniprot_base.format("RAT_10116")
        ]
    parser = Drugbank(open("/Users/aneetauppal/drugbankresults/drugbank.nodes.csv","w"),open("/Users/aneetauppal/drugbankresults/drugbank.edges.csv","w"),"/Users/aneetauppal/Graduate_PhD/DissertationResearch/Knowledgebase/Drugbank/drugbank2/drugbank_all_target_polypeptide_ids.csv/all.csv" ,uniprot_files)
    parser.parse(handle("/Users/aneetauppal/Graduate_PhD/DissertationResearch/Knowledgebase/Drugbank/drugbank2/drugbank_all_full_database.xml.zip"))
    
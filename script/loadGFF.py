import sys
import os
import csv
import gzip
import datetime
import numpy as np
import pandas as pd

DataFrame = pd.DataFrame  
from intervaltree import Interval, IntervalTree


##Load GenomicInterval ##

def loadAnnotation(annotation):
    intervals = {}
    for gene in annotation.gene_annot:
        if annotation.gene_annot[gene].nbExons() == 0:
            continue
        intervals = addInterval(annotation.gene_annot[gene],intervals)
        for exon in annotation.gene_annot[gene].allExons():
            intervals = addInterval(exon,intervals)
    return intervals

def addInterval(genomic_interval,intervals):
    key = str(genomic_interval.chr)+"@"+str(genomic_interval.strand)

    if key not in intervals:
        new_tree = IntervalTree()
        intervals[key] = new_tree

    intervals[key].addi(genomic_interval.start,genomic_interval.end,genomic_interval)
    return intervals


##Load Annotation ##
class AnnotationClass:
    def __init__(self,gene_annot,exon_annot):
        self.gene_annot = gene_annot
        self.exon_annot = exon_annot

class GeneClass:
    def __init__(self,chr,strand,id,symbol,biotype,gff_start,gff_end):
        self.chr = chr
        self.strand = strand
        self.id = id
        self.symbol = symbol
        self.biotype = biotype
        self.gff_start = gff_start
        self.gff_end = gff_end
        self.exon_hash = {}
        self.all_exon = []
        self.start = 999999999999999
        self.end = 0

    def addExon(self,new_exon):
        exon_key = self._getExonKey(new_exon)
        if exon_key in self.exon_hash:
            exon = self.exon_hash[exon_key]
            for transcript in new_exon.allTranscript():
                exon.addTranscript(transcript)
            return exon
        else:
            self.exon_hash[exon_key] = new_exon
            return new_exon

    def nbExons(self):
        return len(self.exon_hash)

    def allExons(self):
        for key in self.exon_hash:
            self.all_exon.append(self.exon_hash[key])
        return self.all_exon

    def length(self):
        return(int(self.gff_end)-int(self.gff_start))

    def sortedExons(self):
        return

    def _getExonKey(self,exon):
        return str(exon.start)+","+str(exon.end)

class ExonClass:
    def __init__(self,chr,strand,start,end,transcript):
        self.chr = chr
        self.strand = strand
        self.start = start
        self.end = end
        self.gene = None
        self.transcript = transcript

    def addTranscript(self,transcript_id):
        try:
            self.transcript.append(transcript_id)
        except:
            self.transcript = [self.transcript]
            self.transcript.append(transcript_id)

    def removeTranscript(self,transcript_id):
        self.transcript.remove(transcript_id)

    def allTranscript(self):
        if type(self.transcript) is list:
            return self.transcript
        else:
            return [self.transcript]

def loadFromGTF(gtf_file):
    gtf_it = gffFileIterator(gtf_file,'gtf')
    genes = {}
    exons = []
    for index, row in gtf_it.iterrows():
        if row["feature"] != "exon":
            continue
        else:
            gene_id = row["attribute_hash"]["gene_id"]
            if gene_id not in genes:
                try:
                    gene_symbol = row["attribute_hash"]["Name"]
                except:
                    gene_symbol = "Unknown"
                gene = GeneClass(row["chr"],row["strand"],gene_id,gene_symbol,row["attribute_hash"]["biotype"],int(row["start"])-1,int(row["end"]))
                genes[gene_id] = gene
            else:
                gene = genes[gene_id]
            exon = ExonClass(row["chr"],row["strand"],int(row["start"])-1,int(row["end"]),row["attribute_hash"]["transcript_id"])
            exon.gene = gene
            genes[gene_id].addExon(exon)
            exons.append(exon)
    for gene_key in genes:
        for exon in genes[gene_key].allExons():
            if exon.start < genes[gene_key].start:
                genes[gene_key].start = exon.start
            if exon.end > genes[gene_key].end:
                genes[gene_key].end = exon.end

    annot = AnnotationClass(genes,exons)
    return annot


def loadFromGFF(gff_file):
    gff_it = gffFileIterator(gff_file,'gff3')

    id_to_parents = {}
    genes = {}
    exons = []
    nb_genes = 0
    nb_exons = 0

    for index, row in gff_it.iterrows():
        try:
            id = row["attribute_hash"]["ID"]
            id = parseEnsemblID(id)
            id = getAtomicGeneID(id)
        except:
            id = None

        try :
            parent = row["attribute_hash"]["Parent"]
            parent = parseEnsemblID(parent)
            parent = getAtomicGeneID(parent)
        except :
            parent = None


        if row["feature"] in ["transcript","mRNA"] and parent is not None and id is not None:
            id_to_parents[id] = parent

        if row["feature"] in ["gene","pseudogene","ncRNA_gene"]:
            gene_id = id
            try:
                gene_symbol = row["attribute_hash"]["Name"]
            except:
                gene_symbol = "Unknown"
            gene_biotype = row["attribute_hash"]["biotype"]

            if gene_id not in genes:
                gene = GeneClass(row["chr"],row["strand"],gene_id,gene_symbol,gene_biotype,int(row["start"])-1,int(row["end"]))
                nb_genes += 1
                genes[gene_id] = gene

        if row["feature"] == "exon":
            exon = ExonClass(row["chr"],row["strand"],int(row["start"])-1,int(row["end"]),parent)
            exons.append(exon)

    for exon in exons:
        transcript = exon.allTranscript()
        if len(transcript) > 0:
            transcript_id = transcript[0]
            gene_id = id_to_parents[transcript_id]
            gene = genes[gene_id]
            exon.gene = gene
            genes[gene_id].addExon(exon)
            nb_exons += 1

    for gene_key in genes:
        if genes[gene_key].nbExons == 0:
            exon = ExonClass(genes[gene_key].chr,genes[gene_key].strand,genes[gene_key].gff_start,genes[gene_key].gff_end,[])
            exon.gene = genes[gene_key]
            exons.append(exon)
            genes[gene_key].start = genes[gene_key].gff_start
            genes[gene_key].end = genes[gene_key].gff_end
        else:
            for exon in genes[gene_key].allExons():
                if exon.start < genes[gene_key].start:
                    genes[gene_key].start = exon.start
                if exon.end > genes[gene_key].end:
                    genes[gene_key].end = exon.end

    annot = AnnotationClass(genes,exons)
    return annot

def gffFileIterator(file,type):
    return getFileIterator(file,type)

def getFileIterator(file,type,parsing_method=None):
    if type is not None and parsing_method is None:
        if type == 'gff3' or type == 'gff2' or type == "gtf":
            header_regex = '#'
            parsing_method = lambda x : parseGFFLine(x,type)
        else:
            sys.exit("Undefined format type")
    parsed_file_raw = []
    with open(file) as fh:
        count_line = 0
        for line in fh:
            count_line += 1
            if header_regex is not None:
                if line.startswith(header_regex):
                    continue
            parsed_line = parsing_method(line)
            parsed_line_count = parsed_line + (count_line,)
            parsed_file_raw.append(parsed_line_count)
    parsed_file = pd.DataFrame(parsed_file_raw, columns = ["chr","source","feature","start","end","score","strand","frame","attribute_hash","line"])
    return parsed_file


def parseGFFLine(line,type):
    if type == "gff3":
        attribute_split = lambda x: x.split("=")
    elif type == "gff2" or type == "gtf" :
        attribute_split = lambda x: x.split('"')
    else:
        sys.exit("Unknown GFF format (must be gff3, gff2 or gtf)")

    chr,source,feature,start,end,score,strand,frame,attribute = line.split("\t")
    attribute_hash = {}
    if attribute is not None :
        attribute_tab = attribute.split(";")
        for attr in attribute_tab:
            try:
                k = attribute_split(attr)[0].strip()
                v = attribute_split(attr)[1].strip()
                attribute_hash[k] = v
            except:
                continue

    conv_attr = {'gene_name':'Name','gene_type':'biotype','gene_biotype':'biotype'}
    for attr in conv_attr:
        if attr in attribute_hash and conv_attr[attr] not in attribute_hash :
            attribute_hash[conv_attr[attr]] = attribute_hash.pop(attr)

    return chr,source,feature,start,end,score,strand,frame,attribute_hash

def getAtomicGeneID(base_id):
    atomic = base_id.split(".")[0]
    try :
        version = base_id.split(".")[1]
        return atomic
    except :
        return base_id

def parseEnsemblID(base_id):
    type = base_id.split(":")[0]
    try :
        id = base_id.split(":")[1]
        return id
    except :
        return base_id

import sys
import os
import csv
import gzip
import datetime
import numpy as np
import pandas as pd
import script.loadGFF as lGFF

DataFrame = pd.DataFrame

from intervaltree import Interval, IntervalTree

def loadFormerTable(path: str, header=0, dtype=None, **kwargs) -> DataFrame:
    #File handler
    options = ""
    if ".tsv" in path:
        df = pd.read_csv(path, dtype=dtype, sep="\\t",engine='python')
    elif ".csv" in path:
        df = pd.read_csv(path, dtype=dtype, sep=",",engine='python')
    else:
        sys.exit("Invalid format of input file. Should be a TSV or CSV file (can be gzipped)")
    return df

class GenomicIntervalClass:
    def __init__(self,chr,start,end):
        self.chr = chr
        self.start = start
        self.end = end
        self.strand = "+"

def fetchByRegion(genomic_interval,intervals):
    key = str(genomic_interval.chr)+"@"+str(genomic_interval.strand)
    try:
        tree = intervals[key]
        start = genomic_interval.start
        end = genomic_interval.end
        query_result = tree.overlap(start,end+1)
    except:
        query_result = set()
    return query_result

def selectBestCanditate(query,result):
    exonic = 0
    intronic = 1
    current_gene = None
    gene_overlap_length = 0

    for res in result:
        if isinstance(res.data,lGFF.GeneClass):
            if exonic == 0:
                current_gene = res.data
        elif isinstance(res.data,lGFF.ExonClass):
            gene_overlap = min(res.data.gene.end,query.end) - max(res.data.gene.start,query.start) + 1
            if gene_overlap > gene_overlap_length:
                gene_overlap_length = gene_overlap
                intronic = 1 if query.start < res.data.start or query.end > res.data.end else 0
                current_gene = res.data.gene
            elif gene_overlap == gene_overlap_length:
                if(res.data.gene.length() > current_gene.length()):
                    intronic = 1 if query.start < res.data.start or query.end > res.data.end else 0
                    current_gene = res.data.gene
            exonic = 1
    return exonic, intronic, current_gene

def reverseStrand(strand):
    if strand == "+":
        return "-"
    elif strand == "-":
        return "+"


def addGFFAnnotation(bam_annot,gff_loaded,output_file,is_stranded):
    GFF_data_raw = []
    output_index = ["gene_id","gene_symbol","gene_biotype","gene_strand","as_gene_id","as_gene_symbol","as_gene_biotype","as_gene_strand","is_exonic","is_intronic"]
    former_table = loadFormerTable(bam_annot)
    for index, row in former_table.iterrows():
        chromosome = row["chromosome"]
        start = row["start"]
        end = row["end"]
        query = GenomicIntervalClass(chromosome,start,end)
        if is_stranded:
            query.strand = row["strand"]
            fwd_res = fetchByRegion(query,gff_loaded)

            query.strand = reverseStrand(row["strand"])
            rev_res = fetchByRegion(query,gff_loaded)
        else:
            query.strand = "+"
            fwd_res = fetchByRegion(query,gff_loaded)
            query.strand = "-"
            fwd_res.update(fetchByRegion(query,gff_loaded))
            rev_res = set()

        fid = fsymbol = fbiotype = fstrand = rid = rsymbol = rbiotype = rstrand = "NA"
        is_exonic = is_intronic = 0
        for strand in ["forward","reverse"]:
            result = fwd_res if strand == "forward" else rev_res
            if result == set():
                candidate = None
            else:
                exonic,intronic,candidate = selectBestCanditate(query,result)
                if strand == "forward":
                    fid = candidate.id
                    fsymbol = candidate.symbol
                    fbiotype = candidate.biotype
                    fstrand = candidate.strand
                    is_exonic = exonic
                    is_intronic = intronic
                else:
                    rid = candidate.id
                    rsymbol = candidate.symbol
                    rbiotype = candidate.biotype
                    rstrand = candidate.strand
        new_line_list = [fid,fsymbol,fbiotype,fstrand,rid,rsymbol,rbiotype,rstrand,str(is_exonic),str(is_intronic)]
        GFF_data_raw.append(new_line_list)
    GFF_data = pd.DataFrame(GFF_data_raw,columns = output_index)
    GFF_data.to_csv(output_file,sep="\t",index = False)

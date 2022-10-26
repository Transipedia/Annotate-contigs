import sys
import os
import csv
import gzip
import datetime
import numpy as np
import pandas as pd

DataFrame = pd.core.frame.DataFrame

import pysam

def extractFromCIGAR(cigar,strand):
    clipped_5p = 0
    clipped_3p = 0
    ref_aln_len = 0
    query_aln_len = 0
    nb_match = 0
    nb_ins = 0
    nb_del = 0
    nb_splice = 0

    if len(cigar)>1:
        if strand in ["+","."]:
            clipped_5p = cigar[0][1] if cigar[0][0] in [4,5] else 0
            clipped_3p = cigar[len(cigar)-1][1] if cigar[len(cigar)-1][0] in [4,5] else 0
        elif strand == "-":
            clipped_5p = cigar[len(cigar)-1][1] if cigar[len(cigar)-1][0] in [4,5] else 0
            clipped_3p = cigar[0][1] if cigar[0][0] in [4,5] else 0
    is_clipped_5p = 1 if clipped_5p != 0 else 0
    is_clipped_3p = 1 if clipped_5p != 0 else 0

    for tuple in cigar:
        if tuple[0] in [0,7,8]: #Match
            ref_aln_len += tuple[1]
            query_aln_len += tuple[1]
            nb_match += tuple[1]
        elif tuple[0] == 1: #Insertion
            query_aln_len += tuple[1]
            nb_ins += 1
        elif tuple[0] == 2: #Deletion
            ref_aln_len += tuple[1]
            nb_del += 1
        elif tuple[0] == 3: #Splice
            ref_aln_len += tuple[1]
            nb_splice += 1

    return nb_match,nb_ins,nb_del,nb_splice,ref_aln_len,query_aln_len,clipped_5p,clipped_3p,is_clipped_5p,is_clipped_3p

def getStrandFromFlag(is_stranded,read):
    strand = "."
    if is_stranded is True:
        if read.is_reverse:
            strand = "-"
        else:
            strand = "+"
    return strand

def addBAMAnnotation(bam,output_file,is_stranded,unique_id_col,genome):
    bam_handle = pysam.AlignmentFile(bam)
    output_index = [unique_id_col,"mapped_to","chromosome","start","end","strand","cigar","nb_ins","nb_del","nb_splice","nb_snv","clipped_5p","clipped_3p","query_cover","alignment_identity","nb_hit","nb_mismatch"]
    BAM_data_raw = []

    count = 0
    for read in bam_handle:
        count += 1
        if read.is_secondary or read.is_supplementary:
            continue
        else: # Select only primary alignments
            tag = read.query_name
            line_in_sam = count
            strand = getStrandFromFlag(is_stranded,read)
            cigar_tuple = read.cigartuples
            cigar = read.cigarstring
            nb_match,nb_ins,nb_del,nb_splice,ref_aln_len,query_aln_len,clipped_5p,clipped_3p,is_clipped_5p,is_clipped_3p = extractFromCIGAR(cigar_tuple,strand)
            chromosome = read.reference_name
            start = read.reference_start - 1
            end = start + ref_aln_len - 1
            query_cover = query_aln_len/read.query_length
            for sup_tag in read.tags:
                if sup_tag[0] =="NH" :
                    nb_hit = sup_tag[1]
                if sup_tag[0] == "NM" :
                    nb_mismatch = sup_tag[1]
            nb_snv = nb_mismatch - nb_ins - nb_del
            alignment_identity = (nb_match - nb_mismatch) / query_aln_len
            new_line_list = [tag,genome,chromosome,str(start),str(end),strand,cigar,str(nb_ins),str(nb_del),str(nb_splice),str(nb_snv),str(clipped_5p),str(clipped_3p),str(query_cover),str(alignment_identity),str(nb_hit),str(nb_mismatch)]
            BAM_data_raw.append(new_line_list)

    BAM_data = pd.DataFrame(BAM_data_raw,columns = output_index)
    BAM_data.to_csv(output_file,sep="\t",index = False)

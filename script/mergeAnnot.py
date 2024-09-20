import sys
import os
import csv
import gzip
import datetime
import numpy as np
import pandas as pd

DataFrame = pd.DataFrame

def loadTable(path: str, dtype=str, **kwargs) -> DataFrame:
    #File handler
    options = ""
    if ".tsv" in path:
        df = pd.read_csv(path, dtype=dtype, sep="\\t",engine='python',**kwargs)
    elif ".csv" in path:
        df = pd.read_csv(path, dtype=dtype, sep=",",engine='python',**kwargs)
    else:
        sys.exit("Invalid format of input file. Should be a TSV or CSV file (can be gzipped)")
    return df


def mergeAll(input_table,bam_and_gff,blast,chim,unique_id_col,column_to_keep,output_file):
    #Load BAM/GFF annotated table
    merged = mergeBAM_GFF(input_table,bam_and_gff,unique_id_col,column_to_keep,output_file)

    #Add Chimeric informations ("junc_type","seg1_cj","seg2_cj")
    chim_loaded = loadTable(chim)
    merged=pd.merge(left=merged,right=chim_loaded,how="outer",on=unique_id_col)
    # We set set off the circular information if the sequence is multimapped.
    merged['nb_hit'] = merged['nb_hit'].fillna("0")
    merged['is_circ'] = np.where(pd.isna(merged["reference_chim"]),"0",merged['is_circ'])
    merged['is_chimeric'] = np.where(pd.isna(merged["reference_chim"]),"0",merged['is_chimeric'])
    # Fill 'mapped_to' column and delete the chimeric reference column
    merged['mapped_to'] = np.where(pd.isna(merged['mapped_to']),merged["reference_chim"],merged["mapped_to"])
    merged.drop('reference_chim',axis=1,inplace=True)


    #Add (optionnal) supplementary mapping informations
    for blast_ali in blast:
        blast_loaded = loadTable(blast_ali)
        merged = pd.merge(left = merged,right = blast_loaded, how = "outer", on = unique_id_col)
    merged.to_csv(output_file,sep="\t", index = False,na_rep="NA")


def mergeBAM_GFF(input_table,bam_and_gff,unique_id_col,column_to_keep,output_file):
    #Load selected columns from initial file
    if column_to_keep == "all":
        base_file = loadTable(input_table)
    else:
        base_file = loadTable(input_table,usecols=column_to_keep)
    #Load BAM/GFF annotation
    bam_gff = loadTable(bam_and_gff)
    #Add BAM/GFF annotation to corresponding contigs
    merged = pd.merge(left = base_file,right = bam_gff, how = "outer", on = unique_id_col)
    merged.fillna(value={"mapped_to":"None"})
    return merged

import sys
import os
import csv
import gzip
import datetime
import numpy as np
import pandas as pd

DataFrame = pd.core.frame.DataFrame

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





def mergeAll(input_table,bam_and_gff,blast,chim, contamination, unique_id_col,column_to_keep,output_file):
    print("===== mergeAll() INPUT FILES =====")
    print("Base file:", input_table)
    print("BAM and GFF file:", bam_and_gff)
    print("Blast files:", blast)
    print("Chim file:", chim)
    print("Contamination hits:", contamination)

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

   #Add (optional) contamination (bacteria, virus, fungi) information
    for cont in contamination: 
        cont_loaded = loadTable(cont)
        merged = pd.merge(left = merged, right = cont_loaded, how = "outer", on = unique_id_col)
    merged.to_csv(output_file, sep="\t", index= False,na_rep="NA")



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



if __name__ == "__main__":
    if len(sys.argv) < 10:
        print("Usage: python script.py <input_table> <bam_and_gff> <blast> <chim> <contamination> <unique_id_col> <column_to_keep> <output_file>")
        sys.exit(1)

    input_table = sys.argv[1]
    bam_and_gff = sys.argv[2]
    blast = sys.argv[3].split(',')
    chim = sys.argv[4]
    contamination = sys.argv[5].split(',')
    unique_id_col = sys.argv[8]
    column_to_keep = sys.argv[9].split(',')
    output_file = sys.argv[10]

    # Check if input files exist
    for file in [input_table, bam_and_gff, chim] + blast + contamination:
        if not os.path.exists(file):
            print(f"Error: File {file} does not exist.")
            sys.exit(1)

    # Run the merging function
    mergeAll(input_table, bam_and_gff, blast, chim, contamination, unique_id_col, column_to_keep, output_file)
    print(f"Merging completed. Output saved to {output_file}.")


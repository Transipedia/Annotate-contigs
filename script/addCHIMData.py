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
    df = pd.read_csv(path, dtype=dtype, sep="\\t",engine='python',header=None, **kwargs)
    return df

def extractChimFromFile(input_file,output_file,unique_id_col,reference):
    output_index = [unique_id_col,"is_chimeric","is_circ","seg1_cj","seg2_cj","reference_chim"]
    CHIM_data_raw = []
    if (os.stat(input_file).st_size == 0): #If no Chimeric junction, print empty temporary file
        CHIM_data = pd.DataFrame(CHIM_data_raw,columns = output_index)
        CHIM_data.to_csv(output_file,sep="\t",index = False) #Temporary file chim_annotation_{ref}.tsv
    else :
        base_chim = loadTable(input_file) #Load Chimeric.out.junction output from STAR
        for index, row in base_chim.iterrows(): #Iter through chimeric file
            # From STAR Manual (5 Chimeric and circular alignments) : https://physiology.med.cornell.edu/faculty/skrabanek/lab/angsd/lecture_notes/STARmanual.pdf
            # 1: chromosome of the donor / 2: first base of the intron of the donor (1-based) / 3: strand of the donor
            # 4: chromosome of the acceptor / 5: first base of the intron of the acceptor (1-based) / 6: strand of the acceptor
            # 7: junction type / 8: repeat length to the left of the junction / 9: repeat length to the right of the junction
            # 10: read name / 11: first base of the first segment (on the + strand) / 12: CIGAR of the first segment / 13: first base of the second segment /14: CIGAR of the second segment

            tag = row[9]
            chimeric = 1
            circ = 0
            chr1 = row[0]
            strand1 = row[2]
            chr2 = row[3]
            strand2 = row[5]
            tag = row[9]
            start1 = row[1]
            start2 = row[4]

            #Test for circular RNA
            if chr1 == chr2 and strand1 == strand2: #If same chromosome and strand
                if abs(int(row[1]) - int(row[4])) <= 30000: #If roughly on the same gene (30000 being the mean size of a gene)
                    #Depending on the strand, the order of fragment isn't the same
                    if strand1 == "+":
                        if int(row[1]) > int(row[4]):
                            circ = 1
                    elif strand1 == "-":
                        if int(row[4]) > int(row[1]):
                            circ = 1

            seg_1 = ":".join([chr1,start1,strand1])
            seg_2 = ":".join([chr2,start2,strand2])
            new_line = [tag,chimeric,circ,seg_1,seg_2,reference+"(chimeric)"]
            CHIM_data_raw.append(new_line)
        CHIM_data = pd.DataFrame(CHIM_data_raw,columns = output_index)
        CHIM_data.to_csv(output_file,sep="\t",index = False) #Temporary file chim_annotation_{ref}.tsv

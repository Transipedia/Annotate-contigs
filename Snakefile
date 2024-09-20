"""
@created: 29/09/21
@modified: 28/10/22
@main_dev: Antoine Laine
    I2BC(SSFA)

Adds general annotation information to any type of sequence. Mainly designed around kamrat/dekupl output formats, and inspired by dekupl-annotation process.
"""
#!/bin/python3 
configfile: "config.json"
import sys
import re
import os
import csv
import gzip
import datetime
import numpy as np
import pandas as pd
import json 
DataFrame = pd.core.frame.DataFrame

import script.addBAMData as aBD
import script.loadGFF as lGFF
import script.addGFFData as aGFFD
import script.addCHIMData as aCD
import script.mergeAnnot as mA

shell.executable('bash')

with open('config.json', 'r') as config_file:
        config = json.load(config_file)

### ========== Variables ========== ###
SEQUENCE_FILE = config['input_file']
SEQUENCE_COLNAME = config['sequence_col'] if 'sequence_col' in config else "contig"
UNIQUE_ID_COL = config['id_col'] if 'id_col' in config else "tag"
OUTPUT_DIR = config['output_dir'] if 'output_dir' in config else "./output"
KEEP_COLUMN = config['keep_col'] if 'keep_col' in config else "all"
PRESET = config['preset'] if 'preset' in config else "map-ont" 
LIB_TYPE = config['library_type'] if 'library_type' in config else "rf"
STRAND = True if LIB_TYPE in ["rf","fr","stranded"] else False
MAP_TO = config['map_to']
SUPP_MAP_TO = config['supp_map_to'] if 'supp_map_to' in config else []
REFERENCE = config['reference']
ANNOTATION = config['annotation']
INDEX_STAR = config['star_index'] 
INDEX_MINIMAP2 = config['minimap2_index']
MODE = config['mode']
MAX_THREADS = config['threads']
TMP_FOLDER = OUTPUT_DIR + "/tmp"
blast_dict = {}
for REF in SUPP_MAP_TO:
    blast_dict[REF] = config["supp_map_to_fasta"]


### ========== Rules ========== ###         

rule all:
    input:
       expand(OUTPUT_DIR + "/{ref}_tmp/annotation.gtf", ref=MAP_TO),
       OUTPUT_DIR + "/query.fa.gz",
       expand(OUTPUT_DIR + "/STAR_{ref}/Aligned.out.sam", ref=MAP_TO),
       expand(OUTPUT_DIR + "/Minimap2_{ref}/Alignement_minimap2.sam", ref=MAP_TO),
       OUTPUT_DIR + "/gff_annotation_minimap2_" + MAP_TO[0] + ".tsv",
       OUTPUT_DIR + "/merged_annotation.tsv"
        
if MODE == "index":
	rule organize_reference:
               input:
                  fasta = expand(REFERENCE),
                  gff = expand(ANNOTATION)
               output:
                  fasta = expand("{ref}_index/reference.fa.gz", ref=MAP_TO),
                  gff = expand("{ref}_index/annotation.gtf.gz", ref=MAP_TO)
               shell:
                  "cp {input.fasta} {output.cp_fasta}"
                  "cp {input.gff} {output.cp_gff}"
	rule build_star_index:
                input:
                  unzip_fasta = expand("{output}/{ref}_tmp/reference.fa",output=OUTPUT_DIR,ref=MAP_TO),
                  unzip_gff = expand("{output}/{ref}_tmp/annotation.gtf",output=OUTPUT_DIR,ref=MAP_TO)
                output:
                  star_sa = expand("{ref}_index/star/SA",ref=MAP_TO)
                params:
                  star = expand("{ref}_index/star",ref=MAP_TO)
                threads: MAX_THREADS
                log: 
                  log_file = expand(OUTPUT_DIR + "/LOGS/{ref}_index/star/{ref}.log",ref=MAP_TO)
                shell:
                   "STAR --runThreadN {threads} --runMode genomeGenerate --genomeDir {params.star} --genomeFastaFiles {input.unzip_fasta} --genomeSAindexNbases 13 -- sjdbGTFfile {input.unzip_gff} > {log.log_file} 2>&1"

	rule build_minimap2_index: 
                input: 
                  expand("{output}/{ref}_tmp/reference.fa",output=OUTPUT_DIR,ref=MAP_TO)
                output: 
                  expand("{ref}_index/minimap2/{ref}.mmi",ref=MAP_TO)
                threads: MAX_THREADS
                log : 
                  log_file = expand(OUTPUT_DIR + "/LOGS/{ref}_index/minimap2/{ref}.log",ref=MAP_TO)
                params: PRESET
                shell:
                   "minimap2 -x {params} -k 14 -d {output} {input} > {log.log_file} 2>&1 "

rule copy_reference:
    input:
        fasta = REFERENCE,
        gff = ANNOTATION
    output:
        cp_fasta = expand(OUTPUT_DIR + "/{ref}_tmp/reference.fa.gz",ref=MAP_TO),
        cp_gff = expand(OUTPUT_DIR + "/{ref}_tmp/annotation.gtf.gz",ref=MAP_TO)
    shell:
        "cp {input.fasta} {output.cp_fasta} ;"
        "cp {input.gff} {output.cp_gff}"

rule gunzip_reference:
    input:
        fasta = OUTPUT_DIR + "/{ref}_tmp/reference.fa.gz",
        gff = OUTPUT_DIR + "/{ref}_tmp/annotation.gtf.gz"
    output:
        unzip_fasta = OUTPUT_DIR + "/{ref}_tmp/reference.fa",
        unzip_gff = OUTPUT_DIR + "/{ref}_tmp/annotation.gtf"
    shell:
        "gunzip -c {input.fasta} > {output.unzip_fasta} ;"
        "gunzip -c {input.gff} > {output.unzip_gff}"


if SEQUENCE_FILE.endswith(".gz"):
    PLAIN_SEQUENCE_FILE = SEQUENCE_FILE[:-3]
    rule unzip_input:
        input:
            base_file = SEQUENCE_FILE
        output:
            plain_base_file = temp(PLAIN_SEQUENCE_FILE)
        shell:
            "gunzip -c {input.base_file} > {output.plain_base_file}"
else:
    PLAIN_SEQUENCE_FILE = SEQUENCE_FILE

rule create_fasta:
    input:
        base_file = PLAIN_SEQUENCE_FILE
    params:
        seq_col = SEQUENCE_COLNAME,
        id_col = UNIQUE_ID_COL
    output:
        fasta = temp(OUTPUT_DIR + "/query.fa"),
        fasta_gz = OUTPUT_DIR + "/query.fa.gz"
    shell:
        "awk -f script/generate_fasta.awk c1={params.id_col} c2={params.seq_col} < {input.base_file}> {output.fasta} ;"
        "gzip -c {output.fasta} > {output.fasta_gz}"

rule split_query :
     input:
        fasta = OUTPUT_DIR + "/query.fa"
     output:
        short = OUTPUT_DIR + "/query_lt_200.fa",
        long = OUTPUT_DIR + "/query_gt_200.fa"
     shell :
        " python script/categorize_sequences.py -fasta {input.fasta} -short {output.short} -long {output.long} -threshold 200 "




if MODE == "index":
	rule primary_alignment:
         input:
            star_sa = "{ref}_index/star/SA",
            query_fasta = OUTPUT_DIR + "/query_lt_200.fa"
         params:
            star_index = "{ref}_index/star",
            out_path = OUTPUT_DIR + "/STAR_{ref}/",
            mode = MODE
         threads: MAX_THREADS
         output:
            chimeric = OUTPUT_DIR + "/STAR_{ref}/Chimeric.out.junction",
            bam = OUTPUT_DIR + "/STAR_{ref}/Aligned.out.sam",
            unmapped = OUTPUT_DIR + "/STAR_{ref}/Unmapped.txt",
            chim_tags = OUTPUT_DIR + "/STAR_{ref}/Chim_tags.txt"
         log :
            log_file = OUTPUT_DIR + "/LOGS/STARalignment_{ref}.log"
         shell:           
            """
            STAR --genomeDir {params.star_index} --runThreadN {threads} --readFilesIn {input.query_fasta} --outFileNamePrefix {params.out_path} --outSAMattributes NH NM nM HI AS --outFilterMultimapNmax 100 --winAnchorMultimapNmax 200 --outReadsUnmapped Fastx  --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMismatchNmax 2 --alignIntronMax 100000 --alignSJoverhangMin 10 --chimSegmentMin 12 --chimJunctionOverhangMin 12 --alignSJDBoverhangMin 10 --chimSegmentReadGapMax 3 --alignSJstitchMismatchNmax 5 -1 5 5 > {log.log_file} 2>&1
            awk '{{print ">"$10}}' {params.out_path}Chimeric.out.junction > {output.chim_tags} 
            set +o pipefail; grep -n -A1 -f {output.chim_tags} {params.out_path}Unmapped.out.mate1 | sed -n 's/^\\([0-9]\\{{1,\\}}\\).*/\\1d/p' | sed -f - {params.out_path}Unmapped.out.mate1 > {output.unmapped}
            """

	rule minimap2_alignement: 
         input: 
            query_fasta = OUTPUT_DIR + "/query_gt_200.fa",
            minimap2_index = "{ref}_index/minimap2/{ref}.mmi"
         output: 
            aln = OUTPUT_DIR + "/Minimap2_{ref}/Alignement_minimap2.sam"
         threads: MAX_THREADS
         params : 
            preset = PRESET,
            mode = MODE
         log:
            log_file = OUTPUT_DIR + "/LOGS/runMinimap2alignment_{ref}.log"
         shell:     
            """  
             minimap2 -ax {params.preset} --MD -k 14 --sam-hit-only {input.minimap2_index} {input.query_fasta} > {output.aln} 2>> {log.log_file} 
            """


if MODE == "table": 
	rule primary_alignment:
         input:
            star_index = INDEX_STAR,
            query_fasta = OUTPUT_DIR + "/query_lt_200.fa"
         params:
            out_path = OUTPUT_DIR + "/STAR_{ref}/",
            mode = MODE
         threads: MAX_THREADS
         output:
            chimeric = OUTPUT_DIR + "/STAR_{ref}/Chimeric.out.junction",
            bam = OUTPUT_DIR + "/STAR_{ref}/Aligned.out.sam",
            unmapped = OUTPUT_DIR + "/STAR_{ref}/Unmapped.txt",
            chim_tags = OUTPUT_DIR + "/STAR_{ref}/Chim_tags.txt"
         log :
            log_file = OUTPUT_DIR + "/LOGS/STARalignment_{ref}.log"
         shell:
            """           
            STAR --genomeDir {input.star_index} --runThreadN {threads} --readFilesIn {input.query_fasta} --outFileNamePrefix {params.out_path} --outSAMattributes NH NM nM HI AS --outFilterMultimapNmax 100 --winAnchorMultimapNmax 200 --outReadsUnmapped Fastx  --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMismatchNmax 2 --alignIntronMax 100000 --alignSJoverhangMin 10 --chimSegmentMin 12 --chimJunctionOverhangMin 12 --alignSJDBoverhangMin 10 --chimSegmentReadGapMax 3 --alignSJstitchMismatchNmax 5 -1 5 5  
            awk '{{print ">"$10}}' {params.out_path}Chimeric.out.junction > {output.chim_tags}
            set +o pipefail; grep -n -A1 -f {output.chim_tags} {params.out_path}Unmapped.out.mate1 | sed -n 's/^\\([0-9]\\{{1,\\}}\\).*/\\1d/p' | sed -f - {params.out_path}Unmapped.out.mate1 > {output.unmapped}            
            """

	rule minimap2_alignement: 
         input: 
            query_fasta = OUTPUT_DIR + "/query_gt_200.fa",
            minimap2_index = INDEX_MINIMAP2
         output: 
            aln = OUTPUT_DIR + "/Minimap2_{ref}/Alignement_minimap2.sam"
         threads: MAX_THREADS
         params : 
            preset = PRESET,
            mode = MODE
         log:
            log_file = OUTPUT_DIR + "/LOGS/runMinimap2alignment_{ref}.log"
         shell:     
            """  
             minimap2 -ax {params.preset} --MD -k 14 --sam-hit-only {input.minimap2_index} {input.query_fasta} > {output.aln} 2>> {log.log_file}
           """





if len(INDEX_STAR) > 1:
    for i in range(len(MAP_TO)-1):
        rule: #secondary_aligmnments
            input:
                star_sa = INDEX[i+1] + "/star/SA",
                query_fasta = OUTPUT_DIR + "/STAR_" + MAP_TO[i] +"/Unmapped.txt"
            params:
                star_index = INDEX[i+1] + "/star",
                out_path = OUTPUT_DIR + "/STAR_" + MAP_TO[i+1] +"/",
                threads = MAX_THREADS
            threads: MAX_THREADS
            output:
                chimeric = OUTPUT_DIR + "/STAR_" + MAP_TO[i+1] + "/Chimeric.out.junction",
                bam = OUTPUT_DIR + "/STAR_" + MAP_TO[i+1] +"/Aligned.out.sam",
                unmapped = OUTPUT_DIR + "/STAR_" + MAP_TO[i+1] +"/Unmapped.txt",
                chim_tags = OUTPUT_DIR + "/STAR_" + MAP_TO[i+1] + "/Chim_tags.txt"
            log :
                OUTPUT_DIR + "/LOGS/runSTARalignment"+ MAP_TO[i+1] +".log"
            shell:
                """
                STAR --genomeDir {params.star_index} --runThreadN {params.threads} --readFilesIn {input.query_fasta} --outFileNamePrefix {params.out_path} --outSAMattributes NH NM nM HI AS --outFilterMultimapNmax 100 --winAnchorMultimapNmax 200 --outReadsUnmapped Fastx --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMismatchNmax 2 --alignIntronMax 100000 --alignSJoverhangMin 10 --chimSegmentMin 12 --chimJunctionOverhangMin 12 --alignSJDBoverhangMin 10 --chimSegmentReadGapMax 3 --alignSJstitchMismatchNmax 5 -1 5 5 > {log.log_file} 2>&1
                awk '{{print \">\"$10}}' {params.out_path}Chimeric.out.junction > {output.chim_tags}
                set +o pipefail; grep -n -A1 -f {output.chim_tags} {params.out_path}Unmapped.out.mate1 | sed -n 's/^\\([0-9]\\{{1,\\}}\\).*/\\1d/p' | sed -f - {params.out_path}Unmapped.out.mate1 > {output.unmapped}
                """

rule remove_STAR_duplicates:
    input:
        bam = OUTPUT_DIR + "/STAR_" + MAP_TO[0] + "/Aligned.out.sam"
    output:
        bam_fixed = OUTPUT_DIR + "/STAR_" + MAP_TO[0] + "/Aligned-fixed.out.sam",
        dedup = OUTPUT_DIR + "/STAR_" + MAP_TO[0] + "/Aligned.out.sam.remove_dedup",
        count_by_query = OUTPUT_DIR + "/STAR_" + MAP_TO[0] + "/count_by_query.txt"
    log :
        OUTPUT_DIR + "/LOGS/remove_STAR_duplicates.log"
    shell:        
        "cat {input.bam} | awk '{{if ($1$3$4$6 != prev) {{print}}; prev=$1$3$4$6}}' > {output.dedup} ;"
        "awk '$1 !~ \"^@\" {{print $1}}' {output.dedup} | uniq -c | awk '{{print $2\"\\t\"$1}}' > {output.count_by_query} ;"
        "awk 'NR==FNR {{nh[$1]=$2;next}} {{if ($1 !~ \"^@\") {{print $1\"\\t\"$2\"\\t\"$3\"\\t\"$4\"\\t\"$5\"\\t\"$6\"\\t\"$7\"\\t\"$8\"\\t\"$9\"\\t\"$10\"\\t\"$11\"\\tNH:i:\"nh[$1]\"\\t\"$13\"\\t\"$14\"\\t\"$15\"\\t\"$16}} else {{print}}}}' {output.count_by_query} {output.dedup} > {output.bam_fixed}" 

rule add_hit_tag:
     input: 
        aln = OUTPUT_DIR + "/Minimap2_{ref}/Alignement_minimap2.sam", 
     output: 
        aln_hit = OUTPUT_DIR + "/Minimap2_{ref}/Alignement_minimap2_add_hit.sam",
     shell: 
        " python script/add_hits.py {input} {output}"

rule add_bam_data:
    input:
        bam_fixed = OUTPUT_DIR + "/STAR_" + MAP_TO[0] + "/Aligned-fixed.out.sam",
        sam_minimap2 = OUTPUT_DIR + "/Minimap2_" + MAP_TO[0] + "/Alignement_minimap2_add_hit.sam"
    params:
        strand = STRAND,
        id_col = UNIQUE_ID_COL,
        reference = MAP_TO[0]
    output:
        bam_annot = temp(OUTPUT_DIR + "/bam_annotation_star_" + MAP_TO[0] + ".tsv"),
        minimap2_annot = temp(OUTPUT_DIR + "/bam_annotation_minimap2_" + MAP_TO[0] + ".tsv")
    run:        
        aBD.addBAMAnnotation(input.bam_fixed, output.bam_annot, params.strand, params.id_col, params.reference)
        aBD.addBAMAnnotation(input.sam_minimap2, output.minimap2_annot, params.strand, params.id_col, params.reference)

rule add_gff_data:
    input:
        unzip_gff = OUTPUT_DIR + "/" + MAP_TO[0] + "_tmp/annotation.gtf",
        bam_annot = OUTPUT_DIR + "/bam_annotation_star_" + MAP_TO[0] + ".tsv",
        minimap2_annot = OUTPUT_DIR + "/bam_annotation_minimap2_" + MAP_TO[0] + ".tsv"
    output: 
        gff_annot = temp(OUTPUT_DIR + "/gff_annotation_star_" + MAP_TO[0] + ".tsv"),
        gff_minimap2_annot = temp(OUTPUT_DIR + "/gff_annotation_minimap2_" + MAP_TO[0] + ".tsv")
    params:
        strand = STRAND
    run:
        annot = lGFF.loadFromGTF(input.unzip_gff)
        annot_interval = lGFF.loadAnnotation(annot)
        aGFFD.addGFFAnnotation(input.bam_annot,annot_interval,output.gff_annot,params.strand)
        aGFFD.addGFFAnnotation(input.minimap2_annot,annot_interval,output.gff_minimap2_annot,params.strand)

#Generating Chimeric file from minimap2 

rule create_chimeric_file:
    input:
        sam = OUTPUT_DIR + "/Minimap2_" + MAP_TO[0] +  "/Alignement_minimap2_add_hit.sam"
    output: 
        chimeric_tmp = expand(OUTPUT_DIR + "/Minimap2_" + MAP_TO[0] + "/Chimerics_minimap2_tmp.out"),
        chimeric = OUTPUT_DIR + "/Minimap2_" + MAP_TO[0] + "/Chimerics_minimap2.out" 
    shell: 
        "grep 'SA:' {input} > {output.chimeric_tmp} ;"
        "python script/add_chimeric.py {output.chimeric_tmp} | sed '1~2d' > {output.chimeric} "


rule add_chimeric_data:
    input:
        chim_file = OUTPUT_DIR + "/STAR_" + MAP_TO[0] + "/Chimeric.out.junction",
        chimeric = OUTPUT_DIR + "/Minimap2_" + MAP_TO[0] + "/Chimerics_minimap2.out"
    output:
        chim_annot = temp(OUTPUT_DIR + "/chim_annotation_star_" + MAP_TO[0] + ".tsv"),
        chim_annot_minimap2 = temp(OUTPUT_DIR + "/chim_annotation_minimap2_" + MAP_TO[0] + ".tsv")
    params:
        id_col = UNIQUE_ID_COL,
        reference = MAP_TO[0]
    run:
        aCD.extractChimFromFile(input.chim_file,output.chim_annot,params.id_col,params.reference)
        aCD.extractChimFromFile(input.chimeric,output.chim_annot_minimap2,params.id_col,params.reference)
        
  


if SUPP_MAP_TO:
    for i in range(len(SUPP_MAP_TO)):
        rule: #supplementaryAlignment
            input:
                query_fasta = OUTPUT_DIR + "/query.fa",
                index_fasta = blast_dict[SUPP_MAP_TO[i]]
            output:
                blast_out = OUTPUT_DIR + "/blast/" + SUPP_MAP_TO[i] + ".tsv"
            params:
                ref = SUPP_MAP_TO[i],
                id_col = UNIQUE_ID_COL
            log: 
                log_file = OUTPUT_DIR + "/LOGS/runBLAST.log"
            run:
                shell("makeblastdb -in {input.index_fasta} -dbtype nucl >> {log.log_file} 2>&1")
                shell("echo -e \"{params.id_col}\t{params.ref}\" > {output.blast_out}")
                shell("blastn -query {input.query_fasta} -db {input.index_fasta} -word_size 11 -max_hsps 1 -max_target_seqs 1 -evalue 1e-3 -outfmt \"6 qseqid sseqid\" 1>>{output.blast_out} 2>>{log.log_file}")


rule merge_files:
     input: 
         bam_star= OUTPUT_DIR + "/bam_annotation_star_" + MAP_TO[0] + ".tsv",
         bam_minimap2=OUTPUT_DIR + "/bam_annotation_minimap2_" + MAP_TO[0] + ".tsv",
         gff_star= OUTPUT_DIR + "/gff_annotation_star_" + MAP_TO[0] + ".tsv",
         gff_minimap2=OUTPUT_DIR + "/gff_annotation_minimap2_" + MAP_TO[0] + ".tsv",
         chim_star= OUTPUT_DIR + "/chim_annotation_star_" + MAP_TO[0] + ".tsv",
         chim_minimap2=OUTPUT_DIR + "/chim_annotation_minimap2_" + MAP_TO[0] + ".tsv"
     output: 
         bam_merged =  temp(OUTPUT_DIR + "/bam_annotation_" + MAP_TO[0] + ".tsv"),
         gff_merged =  temp(OUTPUT_DIR + "/gff_annotation_" + MAP_TO[0] + ".tsv"), 
         chim_merged = temp(OUTPUT_DIR + "/chim_annotation_" + MAP_TO[0] + ".tsv"),
     shell: 
         "cat {input.bam_star} > {output.bam_merged} ;"
         "tail -n +2 {input.bam_minimap2} >> {output.bam_merged} ;"  
         "cat {input.gff_star} > {output.gff_merged} ;"
         "tail -n +2 {input.gff_minimap2} >> {output.gff_merged} ;" 
         "cat {input.chim_star} > {output.chim_merged} ;"
         "tail -n +2 {input.chim_minimap2} >> {output.chim_merged} ;" 


rule merge_annot:
    input:
        bam_annot = expand(OUTPUT_DIR + "/bam_annotation_{ref}.tsv", ref=MAP_TO),
        gff_annot = expand( OUTPUT_DIR + "/gff_annotation_{ref}.tsv", ref=MAP_TO),
        chim_annot = expand( OUTPUT_DIR + "/chim_annotation_{ref}.tsv" ,ref=MAP_TO),
        blast_annot = expand( OUTPUT_DIR + "/blast/{ref}.tsv", ref=SUPP_MAP_TO) if SUPP_MAP_TO else [],
        base_file = SEQUENCE_FILE
    output:
        merged_annot = OUTPUT_DIR + "/merged_annotation.tsv",
        bam_and_gff = temp(OUTPUT_DIR + "/tmp.tsv"),
        bam_all = temp(OUTPUT_DIR + "r/bam.tsv"),
        gff_all = temp(OUTPUT_DIR + "/gff.tsv"),
        chim_all = temp(OUTPUT_DIR + "/chim.tsv"),
    params:
        id_col = UNIQUE_ID_COL,
        keep_col = KEEP_COLUMN
    run:
        shell( "rm -f {output.bam_all} {output.gff_all}" )
        for i in range(len(MAP_TO)):
            shell(r'cat {bam_annot} >> {bam_all}'.format(bam_annot=input.bam_annot[i],bam_all=output.bam_all ))
            shell("cat {gff_annot} >> {gff_all}".format(gff_annot=input.gff_annot[i],gff_all=output.gff_all ))
            shell("cat {chim_annot} >> {chim_all}".format(chim_annot=input.chim_annot[i],chim_all=output.chim_all ))
        shell(r"paste -d$'\t' {output.bam_all} {output.gff_all} > {output.bam_and_gff}")
        mA.mergeAll( input.base_file,output.bam_and_gff,input.blast_annot,output.chim_all,params.id_col,params.keep_col,output.merged_annot )

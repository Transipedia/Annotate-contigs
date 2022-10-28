"""
@created: 29/09/21
@modified: 28/10/22
@main_dev: Antoine Laine
    I2BC(SSFA)

Adds general annotation information to any type of sequence. Mainly designed around kamrat/dekupl output formats, and inspired by dekupl-annotation process.
"""
import sys
import os
import csv
import gzip
import datetime
import numpy as np
import pandas as pd

DataFrame = pd.core.frame.DataFrame

import script.addBAMData as aBD
import script.loadGFF as lGFF
import script.addGFFData as aGFFD
import script.addCHIMData as aCD
import script.mergeAnnot as mA

shell.executable('bash')

### ========== Variables ========== ###
SEQUENCE_FILE = config['input_file']
SEQUENCE_COLNAME = config['sequence_col'] if 'sequence_col' in config else "contig"
UNIQUE_ID_COL = config['id_col'] if 'id_col' in config else "tag"
OUTPUT_DIR = config['output_dir'] if 'output_dir' in config else "./output"
KEEP_COLUMN = config['keep_col'] if 'keep_col' in config else "all"
LIB_TYPE = config['library_type'] if 'library_type' in config else "rf"
STRAND = True if LIB_TYPE in ["rf","fr","stranded"] else False

MAP_TO = config['map_to']
SUPP_MAP_TO = config['supp_map_to'] if 'supp_map_to' in config else []

index_dict = {}
for REF in MAP_TO:
    ref_info=[]
    if REF + "_index" in config:
        ref_info.append(True)
        ref_info.append(config[REF + "_index"])
        ref_info.append("")
        ref_info.append("")
    else :
        ref_info.append(False)
        ref_info.append(REF + "_index")
        ref_info.append(config[REF + "_fasta"])
        ref_info.append(config[REF + "_gff"])

    ref_info.append(ref_info[1] + "/reference.fa.gz")
    ref_info.append(ref_info[1] + "/annotation.gtf.gz")

    index_dict[REF] = ref_info

blast_dict = {}
for REF in SUPP_MAP_TO:
    blast_dict[REF] = config[REF+"_fasta"]

MAX_THREADS = config['threads'] if 'threads' in config else 1
TMP_FOLDER = OUTPUT_DIR + "/tmp"
LOG_FOLDER = OUTPUT_DIR + "/logs"

### ========== Functions ========== ###
def current_date():
    return datetime.datetime.now().strftime("%d %B %Y, %H:%M:%S")

def start_log(log_file, rule_name):
    with open(log_file, "w") as f:
        f.write("******\n")
        f.write("start of rule " + rule_name + " : " + current_date() + "\n")

def end_log(log_file, rule_name):
    with open(log_file, "a") as f:
        f.write("\nend of rule " + rule_name + " : " + current_date() + "\n")
        f.write("******\n");

def add_log(log_file, finished_step):
    with open(log_file, "a") as f:
        f.write("\nDone: " + finished_step + " : " + current_date() + "\n");

def index_to_create(ref_dict):
    list_ref = []
    for ref in ref_dict:
        if ref_dict[ref][0]==False:
            list_ref.append(ref)
    return list_ref

def select_unmapped(ref,output_dir):
    return output_dir + "/STAR_" + ref + "/Unmapped.out.mate1"
### ========== Rules ========== ###

#rule all:
#    input: OUTPUT_DIR + "/human_tmp/reference.fa"

rule all:
    input: OUTPUT_DIR + "/merged_annotation.tsv"

BUILD = index_to_create(index_dict)
INDEX = [index_dict[MAP_TO[i]][1] for i in range(len(MAP_TO))]

if BUILD:
    rule organize_reference:
        input:
            fasta = [index_dict[ref][2] for ref in BUILD],
            gff = [index_dict[ref][3] for ref in BUILD]
        output:
            fasta = expand("{ref}_index/reference.fa.gz",ref=BUILD),
            gff = expand("{ref}_index/annotation.gtf.gz",ref=BUILD)
        run:
            for i in range(len(MAP_TO)):
                shell("cp {fasta} {cp_fasta}".format(fasta=input.fasta[i],cp_fasta=output.fasta[i]))
                shell("cp {gff} {cp_gff}".format(gff=input.gff[i],cp_gff=output.gff[i]))

rule copy_reference:
    input:
        fasta = [index_dict[ref][4] for ref in MAP_TO],
        gff = [index_dict[ref][5] for ref in MAP_TO]
    output:
        fasta = temp(expand("{output}/{ref}_tmp/reference.fa.gz",output=OUTPUT_DIR,ref=MAP_TO)),
        gff = temp(expand("{output}/{ref}_tmp/annotation.gtf.gz",output=OUTPUT_DIR,ref=MAP_TO))
    run:
        for i in range(len(MAP_TO)):
            shell("cp {fasta} {cp_fasta}".format(fasta=input.fasta[i],cp_fasta=output.fasta[i]))
            shell("cp {gff} {cp_gff}".format(gff=input.gff[i],cp_gff=output.gff[i]))

rule gunzip_reference:
    input:
        fasta = OUTPUT_DIR + "/{ref}_tmp/reference.fa.gz",
        gff = OUTPUT_DIR + "/{ref}_tmp/annotation.gtf.gz"
    output:
        unzip_fasta = temp(OUTPUT_DIR + "/{ref}_tmp/reference.fa"),
        unzip_gff = temp(OUTPUT_DIR + "/{ref}_tmp/annotation.gtf")
    run:
        shell("gunzip -c {input.fasta} > {output.unzip_fasta}")
        shell("gunzip -c {input.gff} > {output.unzip_gff}")

if BUILD:
    rule build_star_index:
        input:
            unzip_fasta = expand("{output}/{ref}_tmp/reference.fa",output=OUTPUT_DIR,ref=BUILD),
            unzip_gff = expand("{output}/{ref}_tmp/annotation.gtf",output=OUTPUT_DIR,ref=BUILD)
        output:
            star_sa = expand("{ref}_index/star/SA",ref=BUILD)
        params:
            star = expand("{ref}_index/star",ref=BUILD),
            threads = MAX_THREADS
        threads: MAX_THREADS
        run:
            for i in range(len(BUILD)):
                shell("STAR --runThreadN {threads} --runMode genomeGenerate --genomeDir {star} --genomeFastaFiles {unzip_fasta} --genomeSAindexNbases 14 -- sjdbGTFfile {unzip_gff}".format(threads=params.threads,star=params.star[i],unzip_fasta=input.unzip_fasta[i],unzip_gff=input.unzip_gff[i]))

if SEQUENCE_FILE.endswith(".gz"):
    PLAIN_SEQUENCE_FILE = SEQUENCE_FILE[:-3]
    rule unzip_input:
        input:
            base_file = SEQUENCE_FILE
        output:
            plain_base_file = temp(PLAIN_SEQUENCE_FILE)
        run:
            shell("gunzip -c {input.base_file} > {output.plain_base_file}")
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
    run:
        shell("awk -f /annot/script/generate_fasta.awk c1={params.id_col} c2={params.seq_col} < {input.base_file}> {output.fasta}")
        shell("gzip -c {output.fasta} > {output.fasta_gz}")


rule primary_alignment:
    input:
        star_sa = INDEX[0] + "/star/SA",
        query_fasta = OUTPUT_DIR + "/query.fa.gz"
    params:
        star_index = INDEX[0] + "/star",
        out_path = OUTPUT_DIR + "/STAR_" + MAP_TO[0] + "/" ,
        threads = MAX_THREADS
    threads: MAX_THREADS
    output:
        chimeric = OUTPUT_DIR + "/STAR_" + MAP_TO[0] + "/Chimeric.out.junction",
        bam = OUTPUT_DIR + "/STAR_" + MAP_TO[0] + "/Aligned.out.sam",
        unmapped = OUTPUT_DIR + "/STAR_" + MAP_TO[0] + "/Unmapped.txt",
        chim_tags = OUTPUT_DIR + "/STAR_" + MAP_TO[0] + "/Chim_tags.txt"
    log :
        LOG_FOLDER + "/runSTARalignment"+ MAP_TO[0] +".log"
    run:
        start_log(log[0],"STAR Alignment")
        shell("STARlong --genomeDir {params.star_index} --runThreadN {params.threads} --readFilesIn {input.query_fasta} --readFilesCommand gunzip -c --outFileNamePrefix {params.out_path} --outSAMattributes NH NM nM HI AS --outFilterMultimapNmax 100 --winAnchorMultimapNmax 200 --outReadsUnmapped Fastx  --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMismatchNmax 2 --alignIntronMax 100000 --alignSJoverhangMin 10 --chimSegmentMin 12 --chimJunctionOverhangMin 12 --alignSJDBoverhangMin 10 --chimSegmentReadGapMax 3 --alignSJstitchMismatchNmax 5 -1 5 5")
        shell("awk '{{print \">\"$10}}' {params.out_path}Chimeric.out.junction > {output.chim_tags}")
        shell("set +o pipefail; grep -n -A1 -f {output.chim_tags} {params.out_path}Unmapped.out.mate1 |sed -n 's/^\([0-9]\{{1,\}}\).*/\\1d/p' | sed -f - {params.out_path}Unmapped.out.mate1 > {output.unmapped}")
        end_log(log[0],"STAR Alignment")


if len(INDEX) > 1:
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
                LOG_FOLDER + "/runSTARalignment"+ MAP_TO[i+1] +".log"
            run:
                start_log(log[0],"STAR Alignment")
                shell("STARlong --genomeDir {params.star_index} --runThreadN {params.threads} --readFilesIn {input.query_fasta} --outFileNamePrefix {params.out_path} --outSAMattributes NH NM nM HI AS --outFilterMultimapNmax 100 --winAnchorMultimapNmax 200 --outReadsUnmapped Fastx --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMismatchNmax 2 --alignIntronMax 100000 --alignSJoverhangMin 10 --chimSegmentMin 12 --chimJunctionOverhangMin 12 --alignSJDBoverhangMin 10 --chimSegmentReadGapMax 3 --alignSJstitchMismatchNmax 5 -1 5 5")
                shell("awk '{{print \">\"$10}}' {params.out_path}Chimeric.out.junction > {output.chim_tags}")
                shell("set +o pipefail; grep -n -A1 -f {output.chim_tags} {params.out_path}Unmapped.out.mate1 |sed -n 's/^\([0-9]\{{1,\}}\).*/\\1d/p' | sed -f - {params.out_path}Unmapped.out.mate1 > {output.unmapped}")
                end_log(log[0],"STAR Alignment")


rule add_bam_data:
    input:
        bam = OUTPUT_DIR + "/STAR_{ref}/Aligned.out.sam"
    params:
        strand = STRAND,
        id_col = UNIQUE_ID_COL,
        reference = "{ref}"
    output:
        bam_annot = temp(OUTPUT_DIR + "/bam_annotation_{ref}.tsv")
    log :
        LOG_FOLDER + "/add_bam_data.log"
    run:
        start_log(log[0],"add_bam_data")
        aBD.addBAMAnnotation(input.bam,output.bam_annot,params.strand,params.id_col,params.reference)
        end_log(log[0],"add_bam_data")


rule add_gff_data:
    input:
        unzip_gff = OUTPUT_DIR + "/{ref}_tmp/annotation.gtf",
        bam_annot = OUTPUT_DIR + "/bam_annotation_{ref}.tsv"
    output:
        gff_annot = temp(OUTPUT_DIR + "/gff_annotation_{ref}.tsv")
    params:
        strand = STRAND
    log :
        LOG_FOLDER + "/add_gff_data.log"
    run:
        start_log(log[0],"add_gff_data")
        annot = lGFF.loadFromGTF(input.unzip_gff)
        add_log(log[0],"Loaded GFF")
        annot_interval = lGFF.loadAnnotation(annot)
        add_log(log[0],"Loaded Interval Tree")
        aGFFD.addGFFAnnotation(input.bam_annot,annot_interval,output.gff_annot,params.strand)
        end_log(log[0],"add_gff_data")

rule add_chimeric_data:
    input:
        chim_file = OUTPUT_DIR + "/STAR_{ref}/Chimeric.out.junction"
    output:
        chim_annot = temp(OUTPUT_DIR + "/chim_annotation_{ref}.tsv")
    params:
        id_col = UNIQUE_ID_COL,
        reference = "{ref}"
    log :
        LOG_FOLDER + "/add_chimeric_data.log"
    run:
        start_log(log[0],"add_chimeric_data")
        aCD.extractChimFromFile(input.chim_file,output.chim_annot,params.id_col,params.reference)
        end_log(log[0],"add_chimeric_data")

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
            run:
                shell("makeblastdb -in {input.index_fasta} -dbtype nucl 2>/dev/null")
                shell("echo -e \"{params.id_col}\t{params.ref}\" > {output.blast_out}")
                shell("blastn -query {input.query_fasta} -db {input.index_fasta} -word_size 11 -max_hsps 1 -max_target_seqs 1 -evalue 1e-3 -outfmt \"6 qseqid sseqid\" 1>>{output.blast_out} 2>/dev/null")

rule merge_annot:
    input:
        bam_annot = expand("{output_dir}/bam_annotation_{ref}.tsv",output_dir=OUTPUT_DIR,ref=MAP_TO),
        gff_annot = expand("{output_dir}/gff_annotation_{ref}.tsv",output_dir=OUTPUT_DIR,ref=MAP_TO),
        chim_annot = expand("{output_dir}/chim_annotation_{ref}.tsv",output_dir=OUTPUT_DIR,ref=MAP_TO),
        blast_annot = expand("{output_dir}/blast/{ref}.tsv",output_dir=OUTPUT_DIR,ref=SUPP_MAP_TO) if SUPP_MAP_TO else [],
        base_file = SEQUENCE_FILE
    output:
        merged_annot = OUTPUT_DIR + "/merged_annotation.tsv",
        bam_and_gff = temp(OUTPUT_DIR + "/tmp.tsv"),
        bam_all = temp(OUTPUT_DIR + "/bam.tsv"),
        gff_all = temp(OUTPUT_DIR + "/gff.tsv"),
        chim_all = temp(OUTPUT_DIR + "/chim.tsv"),
    params:
        id_col = UNIQUE_ID_COL,
        keep_col = KEEP_COLUMN
    log:
        LOG_FOLDER + "/mergeAnnot.log"
    run:
        shell("rm -f {output.bam_all} {output.gff_all}")
        for i in range(len(MAP_TO)):
            shell("cat {bam_annot} >> {bam_all}".format(bam_annot=input.bam_annot[i],bam_all=output.bam_all))
            shell("cat {gff_annot} >> {gff_all}".format(gff_annot=input.gff_annot[i],gff_all=output.gff_all))
            shell("cat {chim_annot} >> {chim_all}".format(chim_annot=input.chim_annot[i],chim_all=output.chim_all))
        shell("paste -d$'\t' {output.bam_all} {output.gff_all} > {output.bam_and_gff}")
        mA.mergeAll(input.base_file,output.bam_and_gff,input.blast_annot,output.chim_all,params.id_col,params.keep_col,output.merged_annot)

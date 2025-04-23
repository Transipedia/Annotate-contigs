# Annotate-Contigs:

This tool takes a k-mer/contig table as input and produces an annotation file with mapping information for each k-mer/contig including position, intron/exon/intergenic location, gene name, CIGAR etc. Should work with any reference genome and annoattion as input.
In this pipeline, we use two alignment tools, STAR and Minimap2. STAR is used for aligning sequences that are 200 bases or shorter, while Minimap2 is used for aligning sequences greater than or equal 200 bases.

This tool was initially thought as a less restrictive alternative to [DEkupl-annotation](https://github.com/Transipedia/dekupl-annotation). So many of its aspects are similar.

- [Usage](#usage)
- [Installation](#installation)
    - [Option 1: With singularity](#option-1-with-singularity)
    - [Option 2: From source](#option-2-from-source)
- [Configuration](#configuration)
- [Output files](#output-files)
- [Ontology](#ontology)



## Usage

-   In order to run the tool, you will need at least 4 specific files 
    -   Input file : A tsv/csv file with at least 2 named columns. One contains the sequences you want to annotate, and each sequence must have a unique identifier. Typical input files are Dekupl-run/Kamrat outputs (example available in data/input_table.tsv).
    -   Config file : As the pipeline is designed with snakemake, any run requires a configuration file. See below for specifications of available parameters.
    -   Genome and annotation files : Associated fasta and gtf files of an organism (gz). Typically downloaded from [Ensembl](https://www.ensembl.org/index.html) or [Gencode](https://www.gencodegenes.org/) websites.(examples available in data/)

## Installation

### Clone the Repository

To download this project, use the following command:  
 
    ```
    git clone https://github.com/Transipedia/Annotate-contigs && cd Annotate-contigs
    ```

We recommand to use [singularity](https://singularity.lbl.gov/) to use the tool, or using the manual installation.


### Option 1: With singularity

- **Step 1: Upload Singularity Image**
You can upload the singularity image directly from the [link](https://zenodo.org/records/13789508/files/annotatecontig.sif?download=1)
    ```
    wget https://zenodo.org/records/13789508/files/annotatecontig.sif?download=1 -O annotatecontig.sif
    ```
- **Step 2: Create your configuration file**
This tool is designed to work with Snakemake, which means that all user inputs must be defined in a configuration file (`config.json`). You can find an example of this configuration file in the repository. A comprehensive list of all parameters is provided in the following section.


- **Step 3: Run with mounted volumes**
It is advised to mount certain volumes (input/output directories). By default, a Singularity image cannot access external data. To fix this, you need to mount your directories as volumes. 
Using the parameter `-B /store:/store` tells Singularity to reference your store directory when mentioned (notably in your configuration file). It is recommended that all your input files be located in the `/store` directory.

    ```
    singularity -v run -B /home:/home annotatecontig.sif -s ./Snakefile --configfile ./config.json --cores $nb_cores  
    ```

### Option 2: From conda 

- **Step 1: Install dependancies**. Before using the tool, install the dependencies. You can install them manually using the conda environnement file [annotatecontig.yml](https://zenodo.org/records/13789508/files/annotatecontig.yml?download=1) :

    ```
    conda env create -f annotatecontig.yml
    conda activate AnnotateContig

    ```

- **Step 2: Edit config file & run with Snakemake.**    

    ```
    snakemake --configfile config.json --cores $nb_cores

    ```
## Configuration : 

inside the `config.json` you find all the parameters: 

### Mandatory parameters :

- **mode**: can be either "index" or "table".

  "index": Use this when running the pipeline for the first time to build the indexes.
  "table": Use this if the STAR and Minimap2 indexes already exist, and you only need to generate the table.

- **input_file**: Path to the file containing sequences to annotate (supports tsv/csv, gzipped or uncompressed). Example in (data/input_file.tsv)

- **map_to**: Name of the organism to which the tool will  map your sequences.

- **reference**: Path to the fasta.gz file used to build the index for the specified organism. Only required in "index" mode, otherwise, you put "". Example for test [reference](https://zenodo.org/records/13820050/files/reference.fa.gz?download=1)

- **annotation**: Path to the gtf.gz file used to build the index for the specified organism. Required in "index" and "table" mode. Example for test [annotation](https://zenodo.org/records/13820050/files/part_annotation.gtf.gz?download=1)

- **preset**: (Default : "splice") Adjusts internal parameters of Minimap2 (e.g., k-mer size, scoring schemes, alignment heuristics) to optimize performance and accuracy for specific data types, you can find other presets [here](https://lh3.github.io/minimap2/minimap2.html#8.). 

- **minimap2_index**: Path to the pre-built Minimap2 index for the organism, if previously created. if "index" mode, add ""

- **star_index**: Path to the pre-built STAR index for the organism, if previously created. if "index" mode, add ""

    **The preset used for index building must be consistent. Some presets may not provide information about chimeric reads. In such cases, you may need to build the index again using a different preset.**

### *About the GTF
Only the "exon" features of the GTF file will be used. In order for the program to run properly, the mandatory attributes (column 9) are : "gene_id", "transcript_id", "gene_type".

### Optionnal parameters (and default values) :
- **sequence_col**: (Default :"contig"). Name of the column in input file containing the sequences to annotate.

- **id_col**: (Default:"tag"). Name of the column in input file containing the unique identifier of the sequence.

- **output_dir**: (Default:"./output"). Path to where the results will be generated.

- **keep_col**: (Default:"all"). Either "all" or a list of column names you want to keep from the input file.

- **library_type**: (Default:"rf"). Strandedness: "rf", "fr" or "unstranded".

- **supp_map_to**: (Default:[""]). List of supplementary reference names you want to map your sequences to, with no further information (using blast).

- **supp_map_to_fasta**: For each reference in `supp_map_to`, path to its fasta sequence. An exemple of a typical fasta file you could use (Human repeats from Dfam) is available in data/.

- **contamination** : (Default: False). Toggle contaminant detection.

- **database** : If contamination is set to True, specify the path to the BLAST-generated databases here. 

### *About supplementary alignment

Any amount of supplementary alignment columns can be added to the output. For each supplementary reference provided, a single column will be added at the end of the output file specifying where the annotated sequence was aligned on this reference.  

Example: with a reference of human repeats provided in this repository (data/human_repeat_ref.fasta):
- "supp_map_to":["HumanRepeats"],
- "supp_map_to_fasta" : ["/home/Documents/Annotate-contigs/data/human_repeat_ref.fasta"],

Example with Multiple References:
If you have two supplementary references, e.g., human repeats and viral elements, the configuration would look like this:

- "supp_map_to": ["HumanRepeats", "ViralElements"],
- "supp_map_to_fasta": ["/home/Documents/Annotate-contigs/data/human_repeat_ref.fasta","/home/Documents/Annotate-contigs/data/viral_elements_ref.fasta" ]


## Output file

- Table `merged_annotation.tsv`, summarizing for each contig, its location on the genome (if it's aligned), the sequence alignment informations, and other optionnal alignment informations.

N.B : 
You will also find some intermediate files in the output folder, specifically query_lt_200.fa and query_gt_200.fa. 
- If **query_lt_200.fa** is empty, it means that all the sequences in your query have a length of less than 200 bases, so you will have empty output files from STAR. 
- If **query_gt_200.fa** is empty, it means that all the sequences in your query have a length greater than 200 bases, so you will have empty output files from Minimap2.

## Annotated values

| Term       | Type        | Description         |
| ----------------------- | ----- | ----------- |
| mapped_to               | Str   | Reference to which the sequence was aligned                                               |
| chromosome              | Str   | Chromosome                                                                                         |
| start                   | Int   | Beginning of the alignment on the reference                                                         |
| end                     | Int   | End of the alignment on the reference                                                              |
| strand                  | Char  | Strand of the alignment (+/-). set to "." in unstranded data.                                      |
| cigar                   | Str   | CIGAR string from the SAM alignment.                                                               |
| nb_insertion            | Int   | Number of insertions in the alignment (infered from cigar)                                                 |
| nb_deletion             | Int   | Number of deletions in the alignment (infered from cigar)                                                  |
| nb_splice               | Int   | Number of splices in the alignment (infered from cigar)                                                    |
| nb_snv                  | Int   | Number of SNV in the contigs (computed as the number of mismatches minus indels)                   |
| clipped_3p              | Int   | Number of clipped bases (soft/hard) from 3prim contig                                              |
| clipped_5p              | Int   | Number of clipped bases (soft/hard) from 5prim contig                                              |
| query_cover             | Float | Fraction of the query that have been aligned to the reference                                      |
| alignment_identity      | Float | Fraction of exact match over the query alignment length (splices do not count)                     |
| nb_hit                  | Int   | Number of alignment given for the contig (NH field)                                                |
| nb_mismatches           | Int   | Number of mismatches in the alignment (NM field)                                                   |                         |
| gene_id                 | Str   | Overlapping gene ID (from GTF ID field)                   |
| gene_symbol             | Str   | Overlapping gene symbol (from GTF Name field)            |
| gene_biotype            | Str   | Overlapping gene biotype (from GTF biotype field)        |
| gene_strand             | Char  | Overlapping gene strand (+/-)                           |
| as_gene_id              | Str   | Overlapping antisense gene ID (from GFF ID field). Defined only when working with stranded datas.              |
| as_gene_symbol          | Str   | Overlapping antisense gene symbol (from GFF Name field). Defined only when working with stranded datas.        |
| as_gene_strand          | Char  | Overlapping antisense gene strand (+/-). Defined only when working with stranded datas.                        |
| as_gene_biotype         | Str   | Overlapping antisense gene biotype (from GFF biotype field). Defined only when working with stranded datas.    |
| is_exonic               | Bool  | Overlap between an exon and the contig. Same strand if working with stranded datas, both strand otherwise.   |
| is_intronic             | Bool  | Overlap between an intron and the contig. Same strand if working with stranded datas, both strand otherwise. |
| is_chimeric             | Bool  | The contig contains a chimeric junction                                                             |
| is_circ                 | Bool  | The chimeric junction behaves like a circular RNA.                                                        |
| seg1_cj                 | Str   | First segment of the chimeric junction |
| seg2_cj                 | Str   | Second segment of the chimeric junction |


<details>
</summary> Example usage (click to expand) </summary>

After installing the tool, in the data folder, we have the directory containing the files we'll need.  
Assuming that you have already generated the minimap2 and STAR indexes,  
retrieved the `annotation.gtf.gz` file,  
and also downloaded and decompressed the BLAST database.

### Config file should be like this:

```json
{
  "mode" : "table",
  "input_file": "Annotate-contigs/data/Pancreas_data_example.tsv",
  "sequence_col": "contig",
  "id_col": "tag",
  "map_to" : ["human"],
  "supp_map_to" : ["repeats"],
  "preset" : "splice",
  "supp_map_to_fasta" : ["Annotate-contigs/data/human_repeat_ref.fasta"],
  "output_dir": "Example",
  "keep_col": ["contig","tag"],
  "threads": 5,
  "library_type": "unstranded",
  "reference": "",
  "annotation" : "annotation.gtf.gz",
  "star_index" : "human_index/star",
  "minimap2_index" : "human_index/minimap2/human.mmi",
  "contamination": true,
  "database": "Annotate-contigs/databases/db_bacteria_virus_fungi"
}

```

Then, execute this command line :


```bash
singularity -v run -B /home:/home annotatecontig.sif -s ./Snakefile --configfile ./config.json --cores 4 
```

At the end of the execution, you wil have a folder named "Example" containing the following structure: 

```bash
.
|-- LOGS  
|-- Minimap2_human  
|-- STAR_human  
|-- blast  
|-- contamination  
|-- human_tmp  
|-- merged_annotation.tsv  
|-- query.fa.gz  
|-- query_gt_200.fa  
`-- query_lt_200.fa  

```

The final file, which contains the global annotation table, is merged_annotation.tsv.
It is located in the output folder inside data/.

Note: The last two columns are optional. In this example, we chose to include repeat elements and detect contaminations.

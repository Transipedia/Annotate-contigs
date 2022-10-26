This tool takes a k-mer/contig table as input and produces an annotation file with mapping information for each k-mer/contig including position, intron/exon/intergenic location, gene name, CIGAR etc. Should work with any reference genome and annoattion as input.

This tool was initially thought as a less restrictive alternative to [DEkupl-annotation](https://github.com/Transipedia/dekupl-annotation). So many of its aspects are similar.

- [Usage](#usage)
- [Installation](#installation)
    - [Option 1: With singularity](#option-1-with-singularity)
    - [Option 2: With Docker](#option-2-with-docker)
    - [Option 3: From source (not recommended)](#option-3-from-source-not-recommended)
- [Configuration](#configuration)
- [Output files](#output-files)
- [Ontology](#ontology)

## Usage

-   In order to run the tool, you will need at least 4 specific files (more if you intend to use multiple organisms):
    -   Input file : A tsv/csv file with at least 2 named columns. One contains the sequences you want to annotate, and each sequence must have a unique identifier. Typical input files are Dekupl-run/Kamrat outputs.
    -   Config file : As the pipeline is designed with snakemake, any run requires a cinfiguration file. See below for specifications of available parameters.
    -   Genome and annotation files : Associated fasta and gtf files of an organism, gz.

## Installation

We recommand tu use [singularity](https://singularity.lbl.gov/) to install the tool. Alternatives are Docker and manual installation.

### Option 1: With singularity

- **Step 1: Build Singularity image**
You can create a singularity container from the docker image.
    ```
    singularity build general-annot.simg docker://alaine1/gen-annot:latest
    ```

- **Step 2: Create your configuration file**
The tool was designed with Snakemake, so any user input as to be specified in a configuration file (`config.json`). An example of configuration file is available in the repository. The exhaustive list of all parameters is available in the next section.


- **Step 3: Run with mounted volumes**
It's advised to mount some volumes (input/output directories). Natively, a singularity image isn't able to access external data. To fix this, you have to mount your your directories as volumes.
Using the parameter "-B /store:/store" will indicate singularity to reference your "store" directory if it is mentionned (in your configuration file, notably).
    ```
    singularity -v run -B /store:/store ./general-annot.simg --configfile config.json --cores 4  
    ```


### Option 2: With Docker

- **Step 1: Retrieve the docker image.**
    ```
    docker pull alaine1/gen-annot:latest
    ```
- **Step 2: Run**.
    You might need to mount some volumes:
    - Your `config.json` to `/annot/my-config.json`
    - Any folder referenced in your `config.json`
    The following command is an example working with the `config.json` available in the repository.
    ```
    sudo docker run -v ${PWD}/config.json:/annot/config.json -v /store:/store alaine1/gen-annot:latest
    ```


### Option 3: From source (not recommended)

- **Step 1: Install dependancies**. Before using the tool, install the dependencies (You can use conda):
    - **Required**: bash (version >= 4.3.46), blast (version >= 2.5.0+), STAR (version >= 2.5.3), Python (version >= 3.6.3), Snakemake (version >= 3.12.0)
    - Python packages (using pip install): intervaltree-3.1.0, numpy-1.18.1, pandas-1.1.4, pysam-0.16.0.1, python-dateutil-2.8.2, pytz-2021.1, sortedcontainers-2.4.0

- **Step 2: Install the tool from Github repository**
    ```
    git clone https://github.com/aLaine1/general-annotation && cd dekupl-annotation
    ```
- **Step 3: Edit config file & run with Snakemake.**    
    ```
    snakemake
    ```

## Configuration

Your `config.json` should be the only file you have to interact with in order to run the annotation.

### Mandatory parameters :
- **input_file**: Path to the file containing sequences to annotate. (TSV/CSV, gz or not)
- **map_to**: A list of species/organism on which the tool will try to map your sequences. Mapping is sequential and substractive, meaning if a sequence is mapped on the first organism of the list, we won't try to map it on the second, etc...

- EITHER **[organism]_fasta & [organism]_gff** : Links to fasta.gz and gtf.gz* (respectively) used to build the index of said organism, if it's the first time you use the tool with this organism
- OR **[organism]_index**: Path to built index of said organism, if you already used the tool once with this organism.


### *About the GTF
In the way it's designed, only the "exon" features of the GTF file will be used by this tool. In order for the program to annotate properly, the mandatory attributes (column 9) are : "gene_id", "transcript_id", "gene_type".

### Optionnal parameters (and default values) :
- **sequence_col**: (Default :"contig"). Name of the column of your file containing the sequences to annotate.
- **id_col**: (Default:"tag"). Name of the column of your file containing the unique identifier of the sequence.
- **output_dir**: (Default:"./output"). Path to where the results will be generated.
- **keep_col**: (Default:"all"). Either "all" or a list of column names you want to keep from the input file.
- **library_type**: (Default:"rf"). Specify if your data is "stranded", "rf", "fr" or "unstranded".
- **supp_map_to**: (Default:"None"). A list of supplementary references you want to map your sequences to, with no further informations (using blast).
- **[reference]_fasta**: For each reference in `supp_map_to`, path to its fasta sequence.

## Output file

- Table `merged_annotation.tsv`, summarizing for each contig, its location on the genome (if it's aligned), the sequence alignment informations, and other optionnal alignment informations.




## Ontology

WIP

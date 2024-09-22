from Bio import SeqIO
import argparse

parser = argparse.ArgumentParser(description="Script to classify sequences into two categories based on contig length (less than 200 base and greater than 200 base).")
parser.add_argument('-fasta', type=str, help="Path to the input FASTA file")
parser.add_argument('-short', type=str, help="Path to the output FASTA file for contigs shorter than 200 base ")
parser.add_argument('-long', type=str, help="Path to the output FASTA file for contigs longer than 200 base")
parser.add_argument('-threshold', type=int, help="Length threshold for categorizing contigs")


args = parser.parse_args()
input_file = args.fasta
output_short = args.short 
output_long  = args.long
threshold = args.threshold 


with open(output_short, "w") as short_out, open(output_long, "w") as long_out:
    for record in SeqIO.parse(input_file, "fasta"):
        if len(record.seq) < threshold:
            SeqIO.write(record, short_out, "fasta")
        else:
            SeqIO.write(record, long_out, "fasta")

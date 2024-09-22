import sys
import re

def parse_sam_line(sam_line):
    parts = sam_line.strip().split('\t')
    read_id = parts[0]
    flag = int(parts[1])
    chr_donor = parts[2]
    pos_donor = int(parts[3])
    cigar_donor = parts[5]

    strand_donor = '-' if flag & 16 else '+'
    
    sa_tag = re.search(r'SA:Z:(\S+)', sam_line)
    if not sa_tag:
        return None
    
    sa_fields = sa_tag.group(1).split(',')
    chr_acceptor = sa_fields[0]
    pos_acceptor = int(sa_fields[1])
    strand_acceptor = sa_fields[2]
    cigar_acceptor = sa_fields[3]
    
    return f"{chr_donor}\t{pos_donor}\t{strand_donor}\t{chr_acceptor}\t{pos_acceptor}\t{strand_acceptor}\t1\t*\t*\t{read_id}\t{pos_donor}\t{cigar_donor}\t{pos_acceptor}\t{cigar_acceptor}"

def main():
    if len(sys.argv) != 2:
        print("Usage: python sam_to_chimeric_format.py <input_sam_file>")
        return
    
    input_sam_file = sys.argv[1]
    
    with open(input_sam_file, 'r') as f:
        for line in f:
            if line.startswith('@'):
                continue  # skip header lines
            chimeric_output = parse_sam_line(line)
            if chimeric_output:
                print(chimeric_output)

if __name__ == "__main__":
    main()

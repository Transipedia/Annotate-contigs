import sys

def count_hits_and_update_sam(input_sam, output_sam):
    # Dictionary to store counts of hits for each read ID
    read_counts = {}

    # Count the number of hits for each read ID
    with open(input_sam, 'r') as f:
        for line in f:
            if line.startswith('@'):
                continue  
            fields = line.split('\t')
            read_id = fields[0]
            if read_id in read_counts:
                read_counts[read_id] += 1
            else:
                read_counts[read_id] = 1

    # Rewrite the SAM file with NH:i: tag
    with open(input_sam, 'r') as f, open(output_sam, 'w') as out_f:
        for line in f:
            if line.startswith('@'):
                out_f.write(line)
                continue
            fields = line.strip().split('\t')
            read_id = fields[0]
            nh_tag = f'NH:i:{read_counts[read_id]}'
            fields.append(nh_tag)
            out_f.write('\t'.join(fields) + '\n')

if __name__ == '__main__':
    if len(sys.argv) != 3:
        sys.exit(1)

    input_sam = sys.argv[1]
    output_sam = sys.argv[2]
    count_hits_and_update_sam(input_sam, output_sam)

import sys
import gzip
from Bio import SeqIO

def lowercase_regions(seq_record):
    start = None
    for i, nucleotide in enumerate(seq_record.seq):
        if nucleotide.islower() and start is None:
            start = i
        elif nucleotide.isupper() and start is not None:
            yield (start, i)
            start = None
    if start is not None:
        yield (start, len(seq_record.seq))

def process_fasta_to_bed(input_file, output_file):
    with open(output_file, 'w') as out:
        open_func = gzip.open if input_file.endswith('.gz') else open
        with open_func(input_file, 'rt') as handle:
            for record in SeqIO.parse(handle, 'fasta'):
                for start, end in lowercase_regions(record):
                    out.write(f"{record.id}\t{start}\t{end}\n")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: softmasked_to_bed.py <softmasked_fasta(.gz)> <output_bed>")
        sys.exit(1)
    input_fasta = sys.argv[1]
    output_bed = sys.argv[2]
    process_fasta_to_bed(input_fasta, output_bed)

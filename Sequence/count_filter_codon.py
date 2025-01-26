from Bio import SeqIO
from Bio.Seq import Seq
import sys

def filter_genes_with_stop_codons(dna_file, output_file):
    with open(dna_file, 'r') as infile, open(output_file, 'w') as outfile:
        for dna_record in SeqIO.parse(infile, "fasta"):
            dna_sequence = str(dna_record.seq)

            # Translate the DNA sequence
            translated_seq = str(Seq(dna_sequence).translate())

            # Count the number of stop codons
            stop_codon_count = translated_seq.count('*')

            # If there are more than 2 stop codons, skip this gene
            if stop_codon_count > 2:
                print("Gene {} contains more than 2 stop codons. Deleting.".format(dna_record.id))
            else:
                # Write the sequence to the output file
                SeqIO.write(dna_record, outfile, "fasta")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python filter_genes_with_stop_codons.py <dna_file> <output_file>")
        sys.exit(1)

    dna_file = sys.argv[1]
    output_file = sys.argv[2]

    filter_genes_with_stop_codons(dna_file, output_file)
    print("Processed file saved as {}".format(output_file))


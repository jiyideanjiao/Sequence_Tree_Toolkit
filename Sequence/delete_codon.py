from Bio import SeqIO
from Bio.Seq import Seq
import sys

def remove_stop_codons(dna_file, output_file):
    with open(dna_file, 'r') as infile, open(output_file, 'w') as outfile:
        for dna_record in SeqIO.parse(infile, "fasta"):
            dna_sequence = str(dna_record.seq)

            # Trim the DNA sequence to the nearest codon if necessary
            if len(dna_sequence) % 3 != 0:
                print("Sequence {} length is not a multiple of 3. Trimming.".format(dna_record.id))
                dna_sequence = dna_sequence[:-(len(dna_sequence) % 3)]

            # Translate the DNA sequence
            translated_seq = str(Seq(dna_sequence).translate())

            # Find positions of stop codons in the translated sequence
            stop_codon_positions = [i for i, aa in enumerate(translated_seq) if aa == '*']

            # Remove corresponding codons in the DNA sequence
            for pos in stop_codon_positions[::-1]:  # Reverse order to maintain index correctness
                dna_sequence = dna_sequence[:pos*3] + dna_sequence[(pos*3)+3:]

            dna_record.seq = Seq(dna_sequence)  # Update the sequence
            SeqIO.write(dna_record, outfile, "fasta")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python remove_stop_codons.py <dna_file> <output_file>")
        sys.exit(1)

    dna_file = sys.argv[1]
    output_file = sys.argv[2]

    remove_stop_codons(dna_file, output_file)
    print("Processed file saved as {}".format(output_file))

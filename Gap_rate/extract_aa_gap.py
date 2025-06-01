#!/usr/bin/env python3

import sys

def parse_fasta(filename):
    """
    A simple FASTA parser yielding (sequence_id, sequence) tuples.
    """
    name = None
    seq_chunks = []
    
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue  # Skip empty lines
            if line.startswith(">"):
                # If we already have a sequence stored, yield it
                if name:
                    yield name, "".join(seq_chunks)
                # Start a new sequence
                name = line[1:].split()[0]  # Use the first token after '>'
                seq_chunks = []
            else:
                # Accumulate sequence lines
                seq_chunks.append(line)
    
    # After the loop, yield the last sequence if it exists
    if name:
        yield name, "".join(seq_chunks)


def main():
    if len(sys.argv) != 3:
        print("Usage: python count_seq_length_aa_gap.py <input_fasta> <output_csv>")
        sys.exit(1)

    input_fasta = sys.argv[1]
    output_csv = sys.argv[2]

    # Define gap characters (extend this set if your alignment uses other symbols)
    GAP_CHARS = {"-"}  # Add other gap symbols if necessary, e.g., {"-", ".", "~"}

    with open(output_csv, 'w') as out_handle:
        # Removed the header line
        # out_handle.write("species_name,length,aa,gap\n")

        # Parse input FASTA and compute values
        for seq_id, seq in parse_fasta(input_fasta):
            sequence_length = len(seq)
            gap_count = sum(seq.count(gc) for gc in GAP_CHARS)
            aa_length = sequence_length - gap_count  # Alternatively: len(seq.replace("-", ""))
            out_handle.write("{},{},{},{}\n".format(seq_id, sequence_length, aa_length, gap_count))


if __name__ == "__main__":
    main()

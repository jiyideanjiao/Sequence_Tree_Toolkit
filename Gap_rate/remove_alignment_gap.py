#!/usr/bin/env python

import sys

def remove_gaps(input_file, output_file):
    with open(input_file, 'r') as fin, open(output_file, 'w') as fout:
        current_header = ''
        current_sequence = ''
        for line in fin:
            line = line.rstrip('\n')
            if line.startswith('>'):
                if current_header:
                    # Write the previous header and cleaned sequence
                    fout.write(current_header + '\n')
                    fout.write(current_sequence + '\n')
                current_header = line
                current_sequence = ''
            else:
                # Remove '-' characters and append to the sequence
                current_sequence += line.replace('-', '')
        # Write the last sequence
        if current_header:
            fout.write(current_header + '\n')
            fout.write(current_sequence + '\n')

def main():
    if len(sys.argv) != 3:
        print "Usage: python remove_gaps.py input.fasta output.fasta"
        sys.exit(1)
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    remove_gaps(input_file, output_file)

if __name__ == "__main__":
    main()

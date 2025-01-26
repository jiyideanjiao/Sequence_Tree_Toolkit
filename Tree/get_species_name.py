#!/usr/bin/env python

import sys

def extract_species_names(input_file, output_file):
    with open(input_file, 'r') as fin, open(output_file, 'w') as fout:
        for line in fin:
            line = line.strip()
            if line.startswith('>'):
                species_name = line[1:].strip()
                fout.write(species_name + '\n')

def main():
    if len(sys.argv) != 3:
        print "Usage: python get_species_name.py input.fa output.txt"
        sys.exit(1)
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    extract_species_names(input_file, output_file)

if __name__ == "__main__":
    main()

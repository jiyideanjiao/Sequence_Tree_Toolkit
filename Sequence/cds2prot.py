#!/usr/bin/env python2
# -*- coding: utf-8 -*-

from __future__ import print_function
import sys
import argparse

TABLES = {
    'standard': {
        'TTT':'F','TTC':'F','TTA':'L','TTG':'L',
        'TCT':'S','TCC':'S','TCA':'S','TCG':'S',
        'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*',
        'TGT':'C','TGC':'C','TGA':'*','TGG':'W',
        'CTT':'L','CTC':'L','CTA':'L','CTG':'L',
        'CCT':'P','CCC':'P','CCA':'P','CCG':'P',
        'CAT':'H','CAC':'H','CAA':'Q','CAG':'Q',
        'CGT':'R','CGC':'R','CGA':'R','CGG':'R',
        'ATT':'I','ATC':'I','ATA':'I','ATG':'M',
        'ACT':'T','ACC':'T','ACA':'T','ACG':'T',
        'AAT':'N','AAC':'N','AAA':'K','AAG':'K',
        'AGT':'S','AGC':'S','AGA':'R','AGG':'R',
        'GTT':'V','GTC':'V','GTA':'V','GTG':'V',
        'GCT':'A','GCC':'A','GCA':'A','GCG':'A',
        'GAT':'D','GAC':'D','GAA':'E','GAG':'E',
        'GGT':'G','GGC':'G','GGA':'G','GGG':'G',
    },
    # Vertebrate mitochondrial: TGA=W, AGA/AGG=*, ATA=M
    'vertmt': {
        'TTT':'F','TTC':'F','TTA':'L','TTG':'L',
        'TCT':'S','TCC':'S','TCA':'S','TCG':'S',
        'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*',
        'TGT':'C','TGC':'C','TGA':'W','TGG':'W',
        'CTT':'L','CTC':'L','CTA':'L','CTG':'L',
        'CCT':'P','CCC':'P','CCA':'P','CCG':'P',
        'CAT':'H','CAC':'H','CAA':'Q','CAG':'Q',
        'CGT':'R','CGC':'R','CGA':'R','CGG':'R',
        'ATT':'I','ATC':'I','ATA':'M','ATG':'M',
        'ACT':'T','ACC':'T','ACA':'T','ACG':'T',
        'AAT':'N','AAC':'N','AAA':'K','AAG':'K',
        'AGT':'S','AGC':'S','AGA':'*','AGG':'*',
        'GTT':'V','GTC':'V','GTA':'V','GTG':'V',
        'GCT':'A','GCC':'A','GCA':'A','GCG':'A',
        'GAT':'D','GAC':'D','GAA':'E','GAG':'E',
        'GGT':'G','GGC':'G','GGA':'G','GGG':'G',
    },
    # Invertebrate mitochondrial: TGA=W, AGA/AGG=S, ATA=M
    'invertmt': {
        'TTT':'F','TTC':'F','TTA':'L','TTG':'L',
        'TCT':'S','TCC':'S','TCA':'S','TCG':'S',
        'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*',
        'TGT':'C','TGC':'C','TGA':'W','TGG':'W',
        'CTT':'L','CTC':'L','CTA':'L','CTG':'L',
        'CCT':'P','CCC':'P','CCA':'P','CCG':'P',
        'CAT':'H','CAC':'H','CAA':'Q','CAG':'Q',
        'CGT':'R','CGC':'R','CGA':'R','CGG':'R',
        'ATT':'I','ATC':'I','ATA':'M','ATG':'M',
        'ACT':'T','ACC':'T','ACA':'T','ACG':'T',
        'AAT':'N','AAC':'N','AAA':'K','AAG':'K',
        'AGT':'S','AGC':'S','AGA':'S','AGG':'S',
        'GTT':'V','GTC':'V','GTA':'V','GTG':'V',
        'GCT':'A','GCC':'A','GCA':'A','GCG':'A',
        'GAT':'D','GAC':'D','GAA':'E','GAG':'E',
        'GGT':'G','GGC':'G','GGA':'G','GGG':'G',
    }
}

def parse_args():
    p = argparse.ArgumentParser(
        description="Translate CDS FASTA to protein (Python 2.7 compatible)."
    )
    p.add_argument('fastas', nargs='+', help='Input CDS FASTA file(s)')
    p.add_argument('--table', '-g', default='standard',
                   choices=['standard', 'vertmt', 'invertmt'],
                   help='Genetic code (default: standard)')
    p.add_argument('--trim_at_stop', action='store_true',
                   help='Truncate at first stop (*)')
    p.add_argument('--drop_terminal_stop', action='store_true',
                   help='Drop a single terminal * if present')
    p.add_argument('--frame', type=int, default=1, choices=[1,2,3],
                   help='Reading frame (default: 1)')
    p.add_argument('--width', type=int, default=80,
                   help='Wrap width (default: 80)')
    p.add_argument('--clean_id', action='store_true',
                   help='Keep only first token in header')
    return p.parse_args()

def read_fastas(paths):
    order = []
    seqs = {}
    for path in paths:
        with open(path, 'r') as fh:
            hdr = None
            for raw in fh:
                line = raw.rstrip('\n\r')
                if not line:
                    continue
                if line.startswith('>'):
                    hdr = line
                    order.append(hdr)
                    if hdr not in seqs:
                        seqs[hdr] = []
                else:
                    line = line.upper().replace('U', 'T')
                    line = ''.join(line.split())
                    if hdr is not None:
                        seqs[hdr].append(line)
    # join
    for k in seqs:
        seqs[k] = ''.join(seqs[k])
    return order, seqs

def translate(seq, code, trim_at_stop):
    """Translate DNA (ACGT + IUPAC) -> protein.
       Any codon containing chars outside A/C/G/T => 'X'."""
    L = len(seq)
    prot = []
    i = 0
    while i + 2 < L:
        codon = seq[i:i+3]
        i += 3
        if not codon or len(codon) < 3:
            break
        if _has_non_acgt(codon):
            aa = 'X'
        else:
            aa = code.get(codon, 'X')  # ACGT-only but unexpected => X
        if trim_at_stop and aa == '*':
            break
        prot.append(aa)
    return ''.join(prot)

def _has_non_acgt(codon):
    for ch in codon:
        if ch not in ('A','C','G','T'):
            return True
    return False

def wrap_seq(s, w):
    if w <= 0:
        return s
    out = []
    for i in range(0, len(s), w):
        out.append(s[i:i+w])
    return '\n'.join(out)

def main():
    args = parse_args()
    code = TABLES[args.table]
    if args.width < 1:
        sys.stderr.write('--width must be >= 1\n')
        sys.exit(2)

    order, seqs = read_fastas(args.fastas)

    for hdr in order:
        dna = seqs.get(hdr, '')
        if not dna:
            continue

        # frame offset
        if args.frame > 1:
            if len(dna) < (args.frame - 1):
                sys.stderr.write('%s: sequence shorter than frame offset; skipping\n' % hdr)
                continue
            dna = dna[args.frame - 1:]

        L = len(dna)
        rem = L % 3
        if rem != 0:
            sys.stderr.write('%s length %d not multiple of 3; trailing %d nt ignored\n' % (hdr, L, rem))

        prot = translate(dna, code, args.trim_at_stop)

        if (not args.trim_at_stop) and args.drop_terminal_stop and prot.endswith('*'):
            prot = prot[:-1]

        out_hdr = hdr
        if args.clean_id:
            h = hdr[1:] if hdr.startswith('>') else hdr
            h = h.split()[0]
            out_hdr = '>' + h

        print(out_hdr)
        print(wrap_seq(prot, args.width))

if __name__ == '__main__':
    main()

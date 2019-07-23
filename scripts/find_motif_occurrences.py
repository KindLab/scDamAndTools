#!/usr/bin/env python3

"""
# Find motif occurences

Finds motif occurences in FASTA or gzipped FASTA file. Outputs the found
positions in BED format to STDOUT.

Note, does not handle revcomp well. But that doesn't matter if
your motif is itself palindromic.
"""


import os
import sys
import gzip
from collections import deque

from support.fastaiter import fasta_to_seqiter, buffered_seqlineiter

BUFSIZE = 1000

def fh_from_filename(filename):
    if filename == "-":
        return sys.stdin.buffer
    elif os.path.splitext(filename)[1] == ".gz":
        return gzip.open(filename, 'rb')
    else:
        return open(filename, 'rb')


def find_motif_occurences(outfh, filename, motif):
    motif = motif.lower().encode("ascii")
    l = len(motif)
    npos = 0

    infh = fh_from_filename(filename)
    for deflines, seqlineiter in fasta_to_seqiter(infh):
        defline = (b"\n".join(deflines)).decode("ascii")
        assert defline.startswith(">")
        chrom = defline[1:].split()[0]

        seqit = buffered_seqlineiter(seqlineiter)

        #currseq = seqit.read(BUFSIZE).lower()
        currseq = b""
        i = 0
        while True:
            currseq += seqit.read(BUFSIZE - len(currseq)).lower()
            lcurrseq = len(currseq)
            if lcurrseq < l:
                break
            ii = 0
            while ii < lcurrseq:
                try:
                    j = currseq[ii:].index(motif)
                    outfh.write("%s\t%d\t%d\n" % (chrom, i + ii + j, i + ii + j + l))
                except ValueError:
                    break
                ii += j + l

            i += lcurrseq - (l - 1)
            currseq = currseq[lcurrseq - (l - 1):]

    return


if __name__ == "__main__":
    import argparse

    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('filename', nargs=1, help="Genome FASTA file")
    ap.add_argument('motif', nargs=1, help="Motif to be found")
    args = ap.parse_args()

    outfh = sys.stdout

    find_motif_occurences(outfh, args.filename[0], args.motif[0])

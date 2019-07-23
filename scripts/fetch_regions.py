#!/usr/bin/env python3

"""
# Generate _in silico_ reads based on motif occurrences

Based on the motif occurrences in the genome and the desired read length,
generate _in silico_ reads.
"""

import os
import sys
import logging
from collections import namedtuple

import pysam


log = logging.getLogger(__name__)
logsh = logging.StreamHandler()
logsh.setLevel(logging.DEBUG)
logfmt = logging.Formatter('[%(asctime)s] [%(levelname)s] %(message)s')
logsh.setFormatter(logfmt)
log.addHandler(logsh)


FASTASEQFMT = ">{chrom}_{start}_{end}_{strand}\n{seq}\n"

Region = namedtuple("region", "chrom start end")


def open_file(fn, mode='rt'):
    if os.path.splitext(fn)[1].lower() == ".gz":
        import gzip
        return gzip.open(fn, mode)
    else:
        return open(fn, mode)


def bediter(fh):
    for line in fh:
        fields = line.rstrip().split("\t")
        yield Region(fields[0], int(fields[1]), int(fields[2]))


_REVCOMPTABLE = str.maketrans("ACGT", "TGCA")
def revcomp(s):
    return s.translate(_REVCOMPTABLE)[::-1]


def write_regions(outfh, bedfile, fastafile, length):
    fasta = pysam.FastaFile(fastafile)
    for region in bediter(open_file(bedfile)):
        seqplus = fasta.fetch(region.chrom, region.start, region.start + length)
        formatargs = {"chrom": region.chrom, "start": region.start, "end": region.end, "strand": "plus", "seq": seqplus}
        outfh.write(FASTASEQFMT.format(**formatargs))

        if (region.end - length) > 0:
            seqminus = revcomp(fasta.fetch(region.chrom, region.end - length, region.end))
            formatargs = {"chrom": region.chrom, "start": region.start, "end": region.end, "strand": "minus", "seq": seqminus}
            outfh.write(FASTASEQFMT.format(**formatargs))

    return


def main():
    import argparse
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('-l', '--length', type=int, required=True)
    ap.add_argument('--verbose', '-v', action='count', required=False, default=0)
    ap.add_argument('--quiet', '-q', action="store_true", required=False)
    ap.add_argument('bedfile')
    ap.add_argument('fastafile')

    args = ap.parse_args()

    if args.quiet:
        log.setLevel(logging.CRITICAL)
    else:
        log.setLevel([logging.WARNING, logging.INFO, logging.DEBUG][min(2, args.verbose)])

    assert os.access(args.bedfile, os.R_OK)
    assert os.access(args.fastafile, os.R_OK)

    outfh = sys.stdout
    write_regions(outfh, args.bedfile, args.fastafile, args.length)

    return


if __name__ == "__main__":
    from signal import signal, SIGPIPE, SIG_DFL
    signal(SIGPIPE, SIG_DFL)

    try:
        main()
    except IOError as e:
        if e.errno != 32:
            raise
    except KeyboardInterrupt as e:
        pass

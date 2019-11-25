#!/usr/bin/env python3

"""
# Generating a transcript count table from CELseq2 alignments

This script takes an alignment file (BAM) as input and generates a table of
UMI-unique transcript counts. Whether an alignment is considered to overlap a
feature is determined by the choice of overlap mode. The available modes are
modeled after the HTseq modes (https://htseq.readthedocs.io/en/release_0.11.1/count.html).

"""

import os
import sys
import logging
from collections import defaultdict, Counter, OrderedDict
from itertools import groupby, chain

import numpy as np
import pandas as pd
import h5py
import pysam
from tqdm import tqdm

from support.stepvector import StepVector
from support.gtfparser import GTFParser


log = logging.getLogger(__name__)
logsh = logging.StreamHandler()
logsh.setLevel(logging.DEBUG)
logfmt = logging.Formatter('[%(asctime)s] [%(levelname)s] %(message)s')
logsh.setFormatter(logfmt)
log.addHandler(logsh)


_MODE_CHOICES = {"intersection-strict", "intersection-nonempty", "union"}
_MAX_HAMMING_DIST_DERIVS = 1
_CHUNKSIZE_MAX = 256 * 1024


def read_transcriptome_from_gtf(gtffile, feature, attr_key):
    gtf = GTFParser(gtffile)
    transcriptome = defaultdict(lambda: defaultdict(lambda: StepVector(set)))  # chrom: strand: stepvector
    keys = set()

    chrom_strand_key = lambda f: (f.iv.chrom, f.iv.strand)
    attr_groupby_key = lambda f: f.attr.get(attr_key)

    gtf = filter(lambda f: f.feature == feature, gtf)
    gtf = chain.from_iterable(map(lambda t: t[1], groupby(gtf, key=chrom_strand_key)))
    grouped_gtf = groupby(gtf, attr_groupby_key)
    progressbar = tqdm(grouped_gtf, desc="Reading GTF", unit="features", smoothing=0.15, disable=None)
    for k, f_iter in grouped_gtf:
        keys.add(k)
        for f in f_iter:
            progressbar.update()
            chrom = f.iv.chrom
            strand = f.iv.strand
            start = f.iv.start
            end = f.iv.end
            transcriptome[chrom][strand].add_value(start, end, {k})

    progressbar.close()

    return transcriptome, keys


def hamming_dist(s1, s2):
    dist = abs(len(s1) - len(s2))
    for char1, char2 in zip(s1, s2):
        if char1 != char2:
            dist += 1

    return dist


def get_umi_from_read(read):
    readname_fields = read.qname.split(":")
    assert readname_fields[-2] == "UMI"
    return readname_fields[-1]


def count_reads_with_umis(bamfns, transcriptome, min_mapq, mode):
    if mode not in _MODE_CHOICES:
        raise ValueError("Invalid mode: %s" % mode)

    is_valid_read = lambda read: not read.is_secondary

    stats = OrderedDict([("total_reads", 0), ("valid_reads", 0), ("ambiguous_reads", 0), ("no_feature_reads", 0), ("unmapped_reads", 0), ("lowmapq_reads", 0), ("unique_counts", 0)])

    # open the BAM files
    bfs = [pysam.Samfile(bamfn, 'rb') for bamfn in bamfns]
    # concatenate them into a single stream
    bf = chain(*map(lambda bf: filter(is_valid_read, bf), bfs))

    umi_counts = defaultdict(lambda: Counter())
    for read in tqdm(bf, desc="Reading BAM files", unit="reads", smoothing=0.15, disable=None):
        stats["total_reads"] += 1
        if read.is_unmapped:
            stats["unmapped_reads"] += 1
            continue
        elif read.mapq < min_mapq:
            stats["lowmapq_reads"] += 1
            continue
        chrom = read.reference_name
        strand = "+" if not read.is_reverse else "-"

        sv = transcriptome[chrom][strand]

        if mode == "intersection-strict":
            overlapping = None
            for block in read.get_blocks():
                # block is (start, end) reference coordinates
                blockstart, blockend = block
                for stepstart, stepend, stepval in sv[blockstart:blockend]:
                    if overlapping is None:
                        overlapping = stepval.copy()
                    else:
                        overlapping = overlapping.intersection(stepval)
            assert overlapping is not None
        elif mode == "intersection-nonempty":
            # as HTseq overlap_mode == "intersection-nonempty"
            overlapping = None
            for block in read.get_blocks():
                # block is (start, end) reference coordinates
                blockstart, blockend = block
                for stepstart, stepend, stepval in sv[blockstart:blockend]:
                    if stepval:
                        if overlapping is None:
                            overlapping = stepval.copy()
                        else:
                            overlapping = overlapping.intersection(stepval)
        elif mode == "union":
            # as HTseq overlap_mode == "union"
            overlapping = set()
            for block in read.get_blocks():
                # block is (start, end) reference coordinates
                blockstart, blockend = block
                for stepstart, stepend, stepval in sv[blockstart:blockend]:
                    overlapping = overlapping.union(stepval)
        else:
            raise ValueError("Invalid mode: %s" % mode)

        if overlapping is None or len(overlapping) == 0:
            stats["no_feature_reads"] += 1
        elif len(overlapping) > 1:
            stats["ambiguous_reads"] += 1
        else:
            assert len(overlapping) == 1
            stats["valid_reads"] += 1
            k = next(iter(overlapping))
            umi_counts[k][get_umi_from_read(read)] += 1

    return umi_counts, stats


def reduce_umi_sequence_derivatives(counts):
    reduced_counts = {}
    num_unique_counts = 0
    for k in counts:
        k_counts = counts[k].copy()
        picked = set()
        unpicked = set(k_counts.keys())
        while unpicked:
            seq = k_counts.most_common(1)[0][0]
            unpicked.remove(seq)
            del k_counts[seq]
            picked.add(seq)
            derivs = set()
            for seq2 in unpicked:
                dist = hamming_dist(seq, seq2)
                assert dist > 0
                if dist <= _MAX_HAMMING_DIST_DERIVS:
                    derivs.add(seq2)

            for seq2 in derivs:
                unpicked.remove(seq2)
                del k_counts[seq2]
        reduced_counts[k] = len(picked)
        num_unique_counts += len(picked)

    return reduced_counts, num_unique_counts


def main():
    import argparse

    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--verbose', '-v', action='count', default=0, required=False)
    ap.add_argument('--quiet', '-q', action='store_true', default=False, required=False)
    ap.add_argument('--gtf', '-g', required=True, help="GTF file")
    ap.add_argument('--mode', '-m', choices=_MODE_CHOICES, default="intersection-strict", required=False, help="Overlap mode for aligning reads to features. Based on HTseq.")
    ap.add_argument('--min-mapq', type=int, default=0, required=False, help="Set a minimum mapping quality as a threshold for using alignments. Default is 0.")
    ap.add_argument('--outfile', '-o', required=False, metavar="HDF5_FILE", help="Output file name for the resulting counts.")
    ap.add_argument('--save-stats', action="store_true", default=False, required=False, help="If set, a file containing the count statistics is generated. Default is False.")
    ap.add_argument('infiles', nargs="+", metavar="bamfile", help="Input BAM file(s).")

    args = ap.parse_args()

    if args.quiet:
        log.setLevel(logging.CRITICAL)
    else:
        log.setLevel([logging.WARNING, logging.INFO, logging.DEBUG][min(2, args.verbose)])

    for fn in args.infiles + [args.gtf]:
        if not os.access(fn, os.R_OK):
            ap.error("Can't read file: %s" % fn)

    if args.outfile:
        if os.access(args.outfile, os.F_OK):
            ap.error("File exists: %s" % args.outfile)

    if (args.min_mapq < 0):
        ap.error("Invalid argument")

    log.debug("Reading transcriptome from GTF file: %s" % args.gtf)
    transcriptome, keys = read_transcriptome_from_gtf(args.gtf, feature="exon", attr_key="gene_id")

    log.debug("Counting reads")
    umi_counts, stats = count_reads_with_umis(args.infiles, transcriptome, args.min_mapq, args.mode)


    log.debug("Reducing UMI sequence derivatives")
    counts, num_unique_counts = reduce_umi_sequence_derivatives(umi_counts)
    stats["unique_counts"] += num_unique_counts
    log.debug("Stats: %s" % "; ".join(["%s: %d" % (k, v) for (k, v) in stats.items()]))
    if args.save_stats:
        var = list(stats.keys())
        val = [stats[k] for k in var]
        stats = pd.DataFrame({'sample': np.repeat(args.outfile.split('/')[-1], len(var)), 'type': var , 'counts': val})
        stats.to_csv(args.outfile.replace('.hdf5', '.stats.tsv'), index=False, header=False, sep="\t")

    log.debug("Writing counts")
    if args.outfile is None:
        outfh = sys.stdout
        for k in sorted(keys):
            if k not in counts:
                n = 0
            else:
                n = counts[k]

            outfh.write("%s\t%d\n" % (k, n))
    else:
        with h5py.File(args.outfile, 'w') as h5file:
            d = np.array([counts.get(k, 0) for k in sorted(keys)]).astype(int)
            chunks = (min(d.shape[0], _CHUNKSIZE_MAX / d.itemsize), )
            h5file.create_dataset("counts", data=d, chunks=chunks, compression=9)

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

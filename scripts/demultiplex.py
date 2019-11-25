#!/usr/bin/env python3

"""
# Demultiplexing of scDamID&T-seq data

Matches barcodes by pre-generating a sequence-to-barcode lookup table, for up
to "m" mismatches.

Additionally, allows demultiplexing of paired-end data (currently, only using a
barcode found in R1).

Additionally, multiple read-layouts are supported. In practice this is useful
for reads that are derived from adapters with different UMI lengths. The
lengths of the barcodes need not be equal, either.
The sequences of the UMIs will be append to the (first field of the) read name.
If UMIs are in the read (or the barcodes are not immediately at the 5'-end),
the barcode file must follow the format:
```
Barcodename <tab> 3-ATCG
[...]
```

where "ATCG" is the barcode sequence, "Barcodename" is the name of the barcode
(will used to name output files, see {name}) and 3 is the length of the UMI (or
otherwise the offset of the barcode). The dash between the offset and the
barcode sequence is mandatory.
"""

import os
import sys
from collections import defaultdict, namedtuple, OrderedDict
import gzip
from functools import reduce
from operator import itemgetter
from itertools import product, repeat, combinations, chain, groupby, accumulate
import logging

log = logging.getLogger(__name__)
logsh = logging.StreamHandler()
logsh.setLevel(logging.DEBUG)
logfmt = logging.Formatter('[%(asctime)s] [%(levelname)s] %(message)s')
logsh.setFormatter(logfmt)
log.addHandler(logsh)


BASES = "ACGT"
_PE_READNAMES = ["R1", "R2"]

FASTQRead = namedtuple("FASTQRead", ("name", "seq", "qual"))

def parse_fastq(infh):
    while True:
        try:
            seqnameline = next(infh)
            assert seqnameline.startswith("@")
            seqname = seqnameline.rstrip()

            seqline = next(infh)
            seq = seqline.rstrip()

            plusline = next(infh)
            assert plusline == "+\n"

            qualline = next(infh)
            qual = qualline.rstrip()

            yield FASTQRead(seqname, seq, qual)
        except StopIteration:
            break

def parse_prefix(s):
    fields = s.split("-")
    if not fields[0].isdigit():
        # no UMI at start, offset is implicitly 0
        fields = ["0"] + fields

    if fields[-1].isdigit():
        # no seq at end, seq is implicitly the empty string
        fields = fields + [""]

    digit_fields = fields[0::2]
    bc_fields = fields[1::2]
    assert len(digit_fields) == len(bc_fields)

    if not all(field.isdigit() for field in digit_fields):
        raise ValueError("Invalid prefix string: %s" % s)

    if not all((field.strip(BASES) == "") for field in bc_fields):
        raise ValueError("Invalid prefix string: %s" % s)

    digits = [int(x) for x in digit_fields]
    lengths = [len(bc) for bc in bc_fields]
    offsets = [sum(digits[:i + 1]) + sum(lengths[:i]) for i in range(len(digits))]
    return tuple(zip(offsets, bc_fields))


def read_prefixes(fn):
    d = dict()
    for line in open(fn, 'r'):
        fields = line.rstrip().split("\t")
        name, s = fields
        prefix = parse_prefix(s)

        assert prefix not in d, "Duplicate prefix found"
        d[prefix] = name

    return d


def find_breaks(prefix_table):
    breaks = sorted(reduce(
        set.union,
        (
            set(chain.from_iterable([(o, o + len(s)) for (o, s) in k]))
            for k in prefix_table
        ),
    ))

    return breaks


def break_prefix(prefix, breaks):
    broken = [None] * (len(breaks) - 1)
    for i, (start, end) in enumerate(zip(breaks[:-1], breaks[1:])):
        for offset, seq in prefix:
            l = len(seq)
            seqoffset_start = max(offset, start) - offset
            seqoffset_end = min(offset + l, max(offset, end)) - offset
            if seqoffset_start < seqoffset_end:
                assert broken[i] is None
                broken[i] = seq[seqoffset_start:seqoffset_end]

    return broken


def match_broken(bk, segment_sets):
    assert len(bk) == len(segment_sets)
    candidates = reduce(
        set.intersection,
        (
            segment_sets[i].get(k, segment_sets[i][None])
            for (i, k) in enumerate(bk) if k is not None
        )
    )
    return candidates


def break_seq(seq, breaks):
    return [seq[a:b] for (a, b) in zip(breaks[:-1], breaks[1:])]


def enumerate_variants(seq, mismatches=1, match_read_wildcard=False):
    if match_read_wildcard:
        bases = BASES + "N"
    else:
        bases = BASES

    l = len(seq)

    for ps in combinations(range(l), mismatches):
        for ps_bases in product(*repeat(bases, mismatches)):
            # build up sequence
            prev = 0
            derivative = ""
            for p, b in zip(ps, ps_bases):
                derivative += seq[prev:p] + b
                prev = p + 1
            derivative += seq[prev:]

            yield derivative


def prefix2pseudoseq(prefix):
    return "".join((seq for (offset, seq) in prefix))


def pseudoseq2prefix(pseq, prefix):
    l = len(prefix2pseudoseq(prefix))
    assert l == len(pseq)
    ol = [(o, len(seq)) for (o, seq) in prefix]
    oo = [0] + list(accumulate([l for (_, l) in ol]))
    assert oo[-1] == l
    return tuple([
        (o, pseq[a:b])
        for ((o, _), (a, b)) in zip(
            ol,
            zip(oo[:-1], oo[1:]),
        )
    ])


def parse_umi(readseq, prefix):
    ol = [(o, len(seq)) for (o, seq) in prefix]
    oe = list(map(sum, ol))
    return "".join(
        readseq[a:b]
        for (a, b) in zip(
            [0] + oe,
            [o for (o, _) in prefix] + [sum(ol[-1])],
        )
    )


def paste_umi_to_readname(name, umiseq):
    fields = name.split(" ")
    field0 = fields[0] + ":UMI:" + umiseq
    return " ".join([field0] + fields[1:])


def demultiplex(iters, segment_sets, breaks, prefix_table, i2prefix, orig_prefixes_lut, match_outfhs, unmatched_outfhs=None, ambiguous_outfhs=None, infofh=None, keep=False):
    stats = OrderedDict([("matching_perfect", 0), ("matching_mismatch", 0), ("unmatching", 0), ("ambiguous", 0)])
    matchingstats = {bc: [0, 0] for bc in match_outfhs}

    readiter = zip(*iters)

    for iseq, reads in enumerate(readiter):
        if iseq % 100000 == 0:
            log.debug("%d reads done" % iseq)

        seq = reads[0].seq  # take R1
        segments = break_seq(seq, breaks)
        candidates = match_broken(segments, segment_sets)

        prefixes = {i2prefix[i] for i in candidates}
        barcodes = reduce(set.union, (prefix_table[prefix] for prefix in prefixes), set())

        if len(candidates) > 1:
            # for sure ambiguous
            # test assumption
            assert len(barcodes) > 1

        if len(barcodes) == 0:
            assert len(candidates) == 0

        if not any(barcodes):
            stats["unmatching"] += 1
            if unmatched_outfhs is not None:
                for unmatched_outfh, read in zip(unmatched_outfhs, reads):
                    unmatched_outfh.write("%s\n%s\n+\n%s\n" % (read.name, read.seq, read.qual))
        elif len(barcodes) > 1:
            stats["ambiguous"] += 1
            if ambiguous_outfhs is not None:
                for ambiguous_outfh, read in zip(ambiguous_outfhs, reads):
                    ambiguous_outfh.write("%s\n%s\n+\n%s\n" % (read.name, read.seq, read.qual))
        else:
            assert len(prefixes) == 1
            bcname = next(iter(barcodes))
            prefix = next(iter(prefixes))
            is_perfect = orig_prefixes_lut[bcname] == prefix

            if is_perfect:
                stats["matching_perfect"] += 1
            else:
                stats["matching_mismatch"] += 1
            matchingstats[bcname][1 - int(is_perfect)] += 1

            umiseq = parse_umi(seq, prefix)
            for iread, (match_outfh, read) in enumerate(zip(match_outfhs[bcname], reads)):
                readname = paste_umi_to_readname(read.name, umiseq)

                readseq = read.seq
                readqual = read.qual

                if (iread == 0) and (not keep):
                    seq_offset = prefix[-1][0] + len(prefix[-1][1])
                    readseq = readseq[seq_offset:]
                    readqual = readqual[seq_offset:]

                match_outfh.write("%s\n%s\n+\n%s\n" % (readname, readseq, readqual))

    log.info("Stats: %s" % "; ".join([("%s: %d" % (k, d)) for (k, d) in sorted(stats.items())]))

    if infofh is not None:
        infofh.write("\n".join([("%s\t%d" % (k, d)) for (k, d) in sorted(stats.items())]))
        infofh.write("\n")
        infofh.write("\n".join([("%s\t%d\t%d" % (k, d_perfect, d_mismatch)) for (k, (d_perfect, d_mismatch)) in sorted(matchingstats.items())]))
        infofh.write("\n")
        infofh.flush()

    return


def main():
    import argparse
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('bcfile', help="TAB separated file with <bcname> and <bcseq> columns")
    ap.add_argument('infiles', help="FASTQ file to be analyzed. Use '-' to read from STDIN. For paired-end input, two files can be listed. In that case, outfmt needs to also contain {readname} where R1 and R2 will be substituted.", nargs="+")
    ap.add_argument('-o', '--outfmt', required=True, help="Format for output files. Should contain {name} (and {readname} for PE input). If the extension ends with '.gz', files that are written are automatically gzipped.")
    ap.add_argument('-m', '--mismatches', required=False, type=int, default=1, help="The number of mismatches that are allowed in the barcode sequence. Default is 1.")
    ap.add_argument('-k', '--keep', required=False, default=False, action="store_true", help="If True, retains barcode sequence in the read. Default is False.")
    ap.add_argument('--verify-no-ambiguous', action="store_true", required=False, default=False, help="If True, verifies whether two or more provided barcode sequences are inherently ambiguous given the number of allowed mismatches. Default is False.")
    ap.add_argument('--unmatched-outfile', required=False, default=None, metavar="FILE", help="File in which to store reads that don't match any provided barcode sequence. Should contain {readname} for PE input")
    ap.add_argument('--ambiguous-outfile', required=False, default=None, metavar="FILE", help="File in which to store reads that match multiple barcode sequences. Should contain {readname} for PE input")
    ap.add_argument('--infofile', required=False, default=None, metavar="FILE", help="File in which to report demultiplexing stats.")
    ap.add_argument('--verbose', '-v', action='count', required=False, default=0)
    ap.add_argument('--quiet', '-q', action="store_true", required=False, default=False)

    args = ap.parse_args()

    if len(args.infiles) <= 0 or len(args.infiles) > 2:
        ap.error("Single-end input: 1 infile; Paired-end input: 2 infiles")

    if sum(map(lambda x: x == "-", args.infiles)) > 1:
        ap.error("Only one input file can be STDIN")

    if args.mismatches < 0:
        ap.error("Invalid number of mismatches")

    if "{name}" not in args.outfmt:
        ap.error("Invalid outfmt: %s" % args.outfmt)

    is_paired_end = len(args.infiles) == 2

    if args.quiet:
        log.setLevel(logging.CRITICAL)
    else:
        log.setLevel([logging.WARNING, logging.INFO, logging.DEBUG][min(2, args.verbose)])

    outdir = "/".join(args.outfmt.split("/")[:-1])
    if not os.path.isdir(outdir):
        os.mkdir(outdir)
    if args.unmatched_outfile is not None:
        outdir = "/".join(args.unmatched_outfile.split("/")[:-1])
        if not os.path.isdir(outdir):
            os.mkdir(outdir)
    if args.ambiguous_outfile is not None:
        outdir = "/".join(args.ambiguous_outfile.split("/")[:-1])
        if not os.path.isdir(outdir):
            os.mkdir(outdir)

    def _open_fn(fn, writeable=False):
        if fn == "-":
            return sys.stdin

        if writeable:
            mode = "wt"
        else:
            mode = "rt"

        if os.path.splitext(fn)[1].lower() == ".gz":
            fh = gzip.open(fn, mode, compresslevel=6)
        else:
            fh = open(fn, mode)

        return fh

    infhs = list(map(_open_fn, args.infiles))

    if args.unmatched_outfile is not None:
        if is_paired_end:
            if "{readname}" not in args.unmatched_outfile:
                ap.error("Need {readname} in outfile format")
            unmatched_outfns = [args.unmatched_outfile.format(readname=readname) for readname in _PE_READNAMES]
        else:
            unmatched_outfns = [args.unmatched_outfile]

        unmatched_outfhs = list(map(lambda fn: _open_fn(fn, writeable=True), unmatched_outfns))
    else:
        unmatched_outfhs = None

    if args.ambiguous_outfile is not None:
        if is_paired_end:
            if "{readname}" not in args.ambiguous_outfile:
                ap.error("Need {readname} in outfile format")
            ambiguous_outfns = [args.ambiguous_outfile.format(readname=readname) for readname in _PE_READNAMES]
        else:
            ambiguous_outfns = [args.ambiguous_outfile]

        ambiguous_outfhs = list(map(lambda fn: _open_fn(fn, writeable=True), ambiguous_outfns))
    else:
        ambiguous_outfhs = None

    if args.infofile is not None:
        infofh = open(args.infofile, 'wt')
    else:
        infofh = None

    log.debug("Reading prefixes from file: %s" % args.bcfile)
    orig_prefixes = read_prefixes(args.bcfile)
    orig_prefixes_lut = {v: k for (k, v) in orig_prefixes.items()}

    log.debug("Enumerating prefix derivatives within %d Hamming distance" % args.mismatches)
    prefix_table = defaultdict(lambda: set())
    for prefix, bcname in orig_prefixes.items():
        pseq = prefix2pseudoseq(prefix)
        for deriv_seq in enumerate_variants(pseq, args.mismatches):
            deriv_prefix = pseudoseq2prefix(deriv_seq, prefix)
            prefix_table[deriv_prefix].add(bcname)
    log.debug("Length of prefix table: %d" % len(prefix_table))

    i2prefix = {i: prefix for (i, prefix) in enumerate(sorted(prefix_table))}

    breaks = find_breaks(prefix_table)
    log.debug("Set of breaks: %s" % breaks)

    # ---
    bpfx = list(map(lambda k: break_prefix(k, breaks), sorted(prefix_table)))
    # now we built a dict for every segment (between two breaks)
    # prefixes that have None for a segment are included in all matches for that segment
    segment_sets = []
    for i in range(len(breaks) - 1):
        wildcard_idxs = set(
            map(
                itemgetter(0),
                filter(
                    lambda i_v: i_v[1] is None,
                    enumerate(map(itemgetter(i), bpfx))
                ),
            ),
        )
        segment_idxs = dict([
            (k, set(map(itemgetter(0), v)))
            for (k, v) in groupby(
                sorted(
                    filter(
                        lambda i_v: i_v[1] is not None,
                        enumerate(map(itemgetter(i), bpfx))
                    ),
                    key=itemgetter(1),
                ),
                key=itemgetter(1),
            )
        ])
        segment_idxs[None] = wildcard_idxs
        for k in segment_idxs:
            segment_idxs[k] |= wildcard_idxs
        segment_sets.append(segment_idxs)

    # check for, and list, ambiguities
    ambiguity_table = dict((
        (i, match_broken(bk, segment_sets))
        for (i, bk) in enumerate(bpfx)
    ))

    assert len(ambiguity_table) == len(prefix_table)

    has_ambiguous = False
    for i, prefix in enumerate(sorted(prefix_table)):
        matches = ambiguity_table[i]
        assert i in matches, "prefix %s doesn't match itself(?!?!)" % prefix
        bcmatches = reduce(set.union, (prefix_table[i2prefix[j]] for j in matches))
        if len(bcmatches) > 1:
            has_ambiguous = True
            log.info(
                "Ambiguous barcode found: %s (%s) matches %s" % (
                    prefix,
                    ", ".join(sorted(prefix_table[prefix])),
                    " and ".join([str((bcname, str(orig_prefixes_lut[bcname]))) for bcname in bcmatches]),
                )
            )

    if has_ambiguous and args.verify_no_ambiguous:
        ap.error("One or more ambiguous prefix(es) found")

    log.debug("Reading prefixes from file: %s [DONE]" % args.bcfile)

    log.debug("Setting up outfiles")
    outfhs = {}
    for prefix, bcname in orig_prefixes.items():
        if is_paired_end:
            outfns = [args.outfmt.format(name=bcname, readname=readname) for readname in _PE_READNAMES]
        else:
            outfns = [args.outfmt.format(name=bcname)]

        assert bcname not in outfhs
        outfhs[bcname] = [_open_fn(outfn, writeable=True) for outfn in outfns]
    log.debug("Setting up outfiles [DONE]")

    log.debug("Demultiplexing")
    fastqiters = list(map(parse_fastq, infhs))

    demultiplex(fastqiters, segment_sets, breaks, prefix_table, i2prefix, orig_prefixes_lut, outfhs, unmatched_outfhs, ambiguous_outfhs, infofh, keep=args.keep)
    log.debug("Demultiplexing [DONE]")

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

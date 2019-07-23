#!/usr/bin/env python3

"""
# Write GATC position array

Generates a table in HDF5 format with all GATC positions genome-wide. For each
chromosome, the first entry is the position of the first GATC occurrence on the
plus strand, subsequent values represent the distance to the next GATC. Only
positions on the plus strand are indicated, GATC positions on the minus strand
can be inferred.
"""

import os
import gzip
from collections import defaultdict

import numpy as np
import h5py


def get_fh(fn, mode='rt'):
    if os.path.splitext(fn)[1].lower() == ".gz":
        fh = gzip.open(fn, mode)
    else:
        fh = open(fn, mode)

    return fh


def read_pos_file(infile, element_length=None):
    infh = get_fh(infile)

    mappable_pos = defaultdict(lambda: set())
    fielditer = (line.rstrip().split("\t") for line in infh)
    for fields in fielditer:
        chrom = fields[0]

        start = int(fields[1])
        end = int(fields[2])

        if element_length is None:
            element_length = end - start
        else:
            assert end - start == element_length

        mappable_pos[chrom].add(start)

    return mappable_pos, element_length


def write_pos_arrays(outfn, mappable_pos, element_length):
    h5file = h5py.File(outfn, 'w')
    h5file.attrs["element_length"] = element_length
    for chrom in sorted(mappable_pos):
        A = np.array(sorted(mappable_pos[chrom]), dtype="uint32")
        B = np.r_[A[0], np.diff(A)]
        h5file.create_dataset(chrom, data=B, compression=9)
        h5file.flush()

    h5file.close()

    return


def main():
    import argparse

    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('-l', '--element-length', default=None, required=False, type=int, help="Length of the motif (e.g. 4 for GATC motif). Default is to infer motif length from input.")
    ap.add_argument('pos_file', help="BED input file with GATC positions.")
    ap.add_argument('--outfile', required=True, help="Output file name for HDF5 posarray output.")

    args = ap.parse_args()

    mappable_pos, element_length = read_pos_file(args.pos_file, args.element_length)
    write_pos_arrays(args.outfile, mappable_pos, element_length)

    return


if __name__ == "__main__":
    main()

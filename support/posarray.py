import numpy as np
import h5py


def read_hdf5(fn, chroms=None):
    h5file = h5py.File(fn, 'r')
    if chroms is None:
        return {chrom: h5file[chrom][:] for chrom in h5file}
    else:
        return {chrom: h5file[chrom][:] for chrom in chroms}


def read_posfile(fn):
    f = h5py.File(fn, 'r')
    pos = {}
    for chrom in f:
        pos[chrom] = f[chrom][:].cumsum()

    return pos


# this needs some thought/work... not all functions take elements as *args, some as list/iterable (e.g. 'sum')
def chromwise(fn, *args, **kwargs):
    aslist = False
    if len(args) == 1:
        items = args[0]
        if type(items) != dict:
            assert type(items) == list
            aslist = True
        else:
            items = args
    else:
        items = args

    chroms = items[0].keys()

    if aslist:
        return {chrom: fn([item[chrom] for item in items], **kwargs) for chrom in chroms}
    else:
        return {chrom: fn(*(arg[chrom] for arg in args), **kwargs) for chrom in chroms}


def flatten(v):
    return np.clip(v, 0, 1).astype("int8")


def project_to_bins(c, pos, chromsizes, binsize):
    assert binsize > 0

    As = {}
    for chrom in c:
        assert chrom in chromsizes
        assert chrom in pos

        dtype = c[chrom].dtype

        p = pos[chrom][:]
        binnedpos = p // binsize

        A = np.zeros((int(np.ceil(chromsizes[chrom] / binsize)), 2), dtype=dtype)
        np.add.at(A, binnedpos, c[chrom])

        As[chrom] = A

    return As


def extract_genomic_signal(positions, pc, width):
    v = next(iter(pc.values()))
    shapesuffix = v.shape[1:]
    dtype = v.dtype

    assert width % 2 == 1

    l = len(positions)
    Z = np.empty((l, width) + shapesuffix, dtype=dtype)
    Zw = np.ones(l, dtype="bool")

    Wh = width // 2
    for i, item in enumerate(positions):
        chrom = item[0]
        c = item[1]

        a = c - Wh
        b = c + Wh + 1

        if (a < 0) or (b >= pc[chrom].shape[0]):
            Zw[i] = False
        else:
            Z[i] = pc[chrom][a:b]

    return Z, Zw

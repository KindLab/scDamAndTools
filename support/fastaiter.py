import io
from itertools import groupby, repeat


isdefline = lambda line: line.startswith(b">")


def fasta_to_seqiter(fh):
    g = groupby(fh, isdefline)
    deflines = []
    for key, it in g:
        if not deflines:
            assert key
            deflines = list(it)
        else:
            yield (deflines, it)
            deflines = []


class RawReaderFromIter(io.RawIOBase):
    def __init__(self, it):
        self.buf = None
        assert hasattr(it, '__next__')
        self.it = it

    def readable(self):
        return True

    def readinto(self, b):
        try:
            l = len(b)
            chunk = self.buf or next(self.it).rstrip()
            o, self.buf = chunk[:l], chunk[l:]
            b[:len(o)] = o
            return len(o)
        except StopIteration:
            return 0


def buffered_seqlineiter(it, buffer_size=io.DEFAULT_BUFFER_SIZE):
    r = RawReaderFromIter(it)
    return io.BufferedReader(r)


if __name__ == "__main__":
    fn = "/home/ircheese/Bioinformatics/references/Homo_sapiens.GRCh38.dna.toplevel.fa.gz"
    import gzip
    fh = gzip.open(fn, 'rb')
    seqiter = fasta_to_seqiter(fh)
    for deflines, seqlineiter in seqiter:
        assert len(deflines) == 1
        defline = deflines[0]
        if defline.split()[-1] != "REF":
            continue
        break

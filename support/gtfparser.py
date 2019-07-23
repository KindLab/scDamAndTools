import sys
import os
import re
from collections import namedtuple

from attrsplitter import parse_GTF_attribute_string


_RE_GTF_META_COMMENT = re.compile("##\s*(\S+)\s+(\S*)")

GenomicInterval = namedtuple("GenomicInterval", "chrom start end strand")


class GenomicFeature(namedtuple("GenomicFeature", "feature source iv score frame attrstr")):
    _lazy_attr = None
    _lazy_name = None

    def _compute_attr_and_name(self):
        (attr, name) = parse_GTF_attribute_string(self.attrstr.encode("ascii"), True)
        self._lazy_attr = attr
        self._lazy_name = name

        return

    @property
    def attr(self):
        if self._lazy_attr is None:
            self._compute_attr_and_name()
        return self._lazy_attr

    @property
    def name(self):
        if self._lazy_name is None:
            self._compute_attr_and_name()
        return self._lazy_name


def get_fh(fn):
    if hasattr(fn, 'read'):
        fh = fn
    else:
        if fn == "-":
            fh = sys.stdin
        elif os.path.splitext(fn)[1].lower() == ".gz":
            import gzip
            fh = gzip.open(fn, 'rt')
        else:
            fh = open(fn, 'rt')

    return fh


class GTFParser():
    def __init__(self, filename, end_included=True):
        self.fh = get_fh(filename)
        self.end_included = end_included
        self.metadata = {}

    def __iter__(self):
        for line in self.fh:
            if line == "\n":
                continue
            if line.startswith('#'):
                if line.startswith("##"):
                    mo = _RE_GTF_META_COMMENT.match(line)
                    if mo:
                        self.metadata[mo.group(1)] = mo.group(2)
                continue

            (seqname, source, feature, start, end, score,
                strand, frame, attributestr) = line.split("\t", 8)

            if self.end_included:
                iv = GenomicInterval(seqname, int(start)-1, int(end), strand)
            else:
                iv = GenomicInterval(seqname, int(start)-1, int(end)-1, strand)

            if score != ".":
                score = float(score)

            if frame != ".":
                frame = int(frame)

            f = GenomicFeature(feature=feature, iv=iv,
                               source=source, score=score,
                               frame=frame, attrstr=attributestr)

            yield f

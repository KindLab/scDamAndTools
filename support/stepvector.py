from sortedcontainers import SortedDict


class StepVector():
    @classmethod
    def sliced(cls, other, start, end):
        newobj = cls(other.datatype, _tree=other._t, _bounds=(start, end))
        return newobj

    def __init__(self, datatype, _tree=None, _bounds=None):
        self.datatype = datatype

        if _tree is not None:
            self._t = _tree
        else:
            self._t = SortedDict()

        if _bounds is not None:
            self._bounds = _bounds
        else:
            self._bounds = (None, None)  # set upon slicing/subsetting

    def __getitem__(self, key):
        if type(key) == slice:
            if (key.step is not None) and (key.step != 1):
                raise ValueError("Invalid step value")

            start = key.start
            end = key.stop

            if self._bounds[0] is not None:
                if start is None:
                    start = self._bounds[0]
                else:
                    if start < self._bounds[0]:
                        raise ValueError("Start out of bounds")
            if self._bounds[1] is not None:
                if end is None:
                    end = self._bounds[1]
                else:
                    if end > self._bounds[1]:
                        raise ValueError("End out of bounds")

            return self.sliced(self, start, end)
        else:
            assert type(key) == int

            if self._bounds[0] is not None:
                if key < self._bounds[0]:
                    raise ValueError("Key out of bounds")
            if self._bounds[1] is not None:
                if key >= self._bounds[0]:
                    raise ValueError("Key out of bounds")

            if self._t:
                try:
                    prevkey = self._floor_key(key)
                    return self._t[prevkey]
                except KeyError:
                    # no item smaller than or equal to key
                    return self.datatype()
            else:
                # empty tree
                return self.datatype()

    def __setitem__(self, key, value):
        if type(key) == slice:
            start = key.start
            end = key.stop
        else:
            assert type(key) == int
            start = key
            end = key + 1

        assert start is not None
        assert end is not None

        assert type(value) == self.datatype
        assert end >= start

        if start == end:
            return

        # check next val
        if self._t:
            try:
                nkey = self._floor_key(end, bisect="right")
                nvalue = self._t[nkey]
            except KeyError:
                nkey = None
                nvalue = None
        else:
            # empty tree
            nkey = None
            nvalue = None

        # check prev val
        if self._t:
            try:
                pkey = self._floor_key(start)
                pvalue = self._t[pkey]
            except KeyError:
                pkey = None
                pvalue = None
        else:
            pkey = None
            pvalue = None

        # remove intermediate steps if any
        if self._t:
            a = self._t.bisect_left(start)
            b = self._t.bisect(end)
            assert a <= b
            del self._t.iloc[a:b]

        # set an end marker if necessary
        if nkey is None:
            self._t[end] = self.datatype()
        elif nvalue != value:
            self._t[end] = nvalue

        # set a start marker if necessary
        if pkey is None or pvalue != value:
            self._t[start] = value

    def __iter__(self):
        start, end = self._bounds

        if not self._t:
            # empty tree
            if start is None or end is None:
                raise StopIteration  # FIXME: can't figure out a better thing to do if only one is None
            else:
                if start < end:
                    yield (start, end, self.datatype())
                raise StopIteration

        if start is None:
            a = 0
        else:
            a = max(0, self._bisect_right(start) - 1)

        if end is None:
            b = len(self._t)
        else:
            b = self._bisect_right(end)

        assert b >= a
        if a == b:
            if a is None:
                start = self._t[a]
            if b is None:
                end = self._t[b]

            if start < end:
                yield (start, end, self.datatype())

            raise StopIteration

        it = self._t.islice(a, b)

        currkey = next(it)
        currvalue = self._t[currkey]
        if start is not None:
            currkey = max(start, currkey)
            if start < currkey:
                yield (start, currkey, self.datatype())

        prevkey, prevvalue = currkey, currvalue
        for currkey in it:
            currvalue = self._t[currkey]
            yield (prevkey, currkey, prevvalue)
            prevkey = currkey
            prevvalue = currvalue

        if end is not None:
            if currkey < end:
                yield (currkey, end, prevvalue)

    def add_value(self, start, end, value):
        assert type(value) == self.datatype

        # can't modify self while iterating over values; will change the tree, and thus fuck up iteration
        items = list(self[start:end])

        for a, b, x in items:
            if self.datatype == set:
                y = x.copy()
                y.update(value)
            else:
                y = x + value

            self[a:b] = y

    def _bisect_left(self, key):
        return self._t.bisect_left(key)

    def _bisect_right(self, key):
        return self._t.bisect_right(key)

    def _floor_key(self, key, bisect="left"):
        """
        Returns the greatest key less than or equal to key
        """

        if bisect == "right":
            p = self._bisect_right(key)
        else:
            p = self._bisect_left(key)

        if p == 0:
            raise KeyError
        else:
            return self._t.iloc[p - 1]

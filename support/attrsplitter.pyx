import sys


cdef char * _DEFAULT_SPLIT = ";"
cdef char * _DEFAULT_QUOTES = "\"\'"
cdef int _DEFAULT_NUMQUOTES = 2
cdef char * _DEFAULT_ESCAPE = "\\"

cpdef inline list quotesafe_split(char * s, char * split=_DEFAULT_SPLIT, char * quotes=_DEFAULT_QUOTES, int numquotechars=_DEFAULT_NUMQUOTES, char * escape=_DEFAULT_ESCAPE):
    cdef list l = []
    cdef int i = 0
    cdef int j = 0
    cdef int p = 0
    cdef bint in_quote = False
    cdef char split_c = split[0]
    cdef char escape_c = escape[0]

    while s[i] != 0:
        for j in range(numquotechars):
            if s[i] == quotes[j]:
                if (i > 0) and (s[i - 1] != escape_c):
                    in_quote = not in_quote
                    break
        else:
            if (not in_quote) and (s[i] == split_c):
                if (i == 0) or (s[i - 1] != escape_c):
                    #l.append(s[p:i] )
                    l.append((p, i))
                    p = i + 1
        i += 1

    #l.append(s[p:])
    l.append((p, i))

    if in_quote:
        raise ValueError("Quote error")

    return l


cpdef parse_GTF_attribute_string(bytes attrstr, bint extra_return_first_value=False):
    cdef char * s = attrstr
    cdef str first_val = "_unnamed_"

    cdef int i
    cdef int a
    cdef int b
    cdef bytes attr
    cdef bytes key
    cdef bytes val
    cdef int quotechar
    cdef str strkey
    cdef str strval

    cdef dict d = {}

    for (i, (a, b)) in enumerate(quotesafe_split(s)):
        attr = attrstr[a:b]

        attr = attr.strip()
        if not attr:
            continue

        key, val = attr.split(maxsplit=1)

        for quotechar in _DEFAULT_QUOTES:
            if val[0] == quotechar:
                assert val[-1] == quotechar
                val = val[1:-1]

        strkey = key.decode("ascii")
        strval = val.decode("ascii")

        if strkey in d:
            if type(d[strkey]) == list:
                d[strkey].append(strval)
            else:
                d[strkey] = [d[strkey], strval]
        else:
            d[strkey] = strval

        if extra_return_first_value and i == 0:
            first_val = strval

    if extra_return_first_value:
        return (d, first_val)

    else:
        return d

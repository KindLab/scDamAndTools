#!/usr/bin/awk -f

BEGIN {
	t = 0;

	L = length(PREFIX);
	for (i = 0; i < L; i++) {
		QUALPREFIX = QUALPREFIX "I";
	}
}
{
	t = t + 1;

	if ((t % 4) == 1) {
		seqname = $0;
	} else if ((t % 4) == 2) {
		seq = $0;
	} else if ((t % 4) == 3) {
		plusline = $0;
	} else if ((t % 4) == 0) {
		qual = $0;
	}

	if ((t % 4) == 0) {
		seq = PREFIX seq;
		qual = QUALPREFIX qual;

		printf("%s\n%s\n%s\n%s\n", seqname, seq, plusline, qual);
	}
}

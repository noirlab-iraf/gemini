# Copyright(c) 2004-2006 Association of Universities for Research in Astronomy, Inc.
#
# GT_DICSTR - Invert strdic.  Given a dictionary string and an index,
# return the text (or the empty string)

procedure gt_dicstr (index, match, maxlen, dict)

int	index			# I The index whose text we require
char	match[ARB]		# O Where to save the text
int	maxlen			# I Length of match
char	dict[ARB]		# I The dictionary

int	idict, count, imatch
char	sep

begin
	match[1] = EOS
	idict = 1
	imatch = 1
	count = 1

	if (dict[idict] == EOS)
	    return
	sep = dict[idict]
	idict = idict + 1

	while (dict[idict] != EOS) {
	    if (count == index) {
		if (dict[idict] == sep || imatch == maxlen)
		    break
		match[imatch] = dict[idict]
		imatch = imatch + 1
	    } else {
		if (dict[idict] == sep)
		    count = count + 1
	    }
	    idict = idict + 1
	}

	match[imatch] = EOS
	return
end


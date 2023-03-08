# Copyright(c) 2004-2006 Association of Universities for Research in Astronomy, Inc.
#

include	<ctotok.h>
include	<error.h>


# GT_GARGCOMMA - A simple wrapper over gt_sscan to split a string using
# commas

# No support for quotes - "A,B" will parse into two tokens: "A and B"

# Returns true if data were found (false at end of list)
# Returns true with a null string between repeated commas 
# Returns false with a null string if no more data available


define	COMMA	','			# Could make this a param...


bool procedure gt_gargcomma (scn, outstr, maxlen)

pointer	scn				# IO The scan data
char	outstr[ARB]			# O Destination for text
int	maxlen				# I maximum size of outstr

char	ch
int	iout
bool	found

begin
	found = false
	iout = 1

	while (iout < maxlen) {
	    call gt_scanc (scn, ch)
	    if (ch == EOS)
		break
	    found = true
	    if (ch == COMMA)
		break
	    outstr[iout] = ch
	    iout = iout+1
	}

	outstr[iout] = EOS

	return found
end

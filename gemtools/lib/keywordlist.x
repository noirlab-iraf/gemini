# Copyright(c) 2004-2006 Association of Universities for Research in Astronomy, Inc.
#

include	<error.h>
include	"../pkg/gemextn/gemerrors.h"


# GT_KEYWORD_LIST -- Convert a comma-separated list of text options to an
# array of the corresponding indices via a dictionary.  The number of
# indices is one less than the first value in the array (this saves
# having to carry round a separate length).

# So, for example, if two values are found (10 and 17), the contents
# of "options" are [3,10,17].  If no values are present, the array
# "options" is [1].

# The error message template should contain a "%s" to receive the
# appropriate value (the unmatched word)

procedure gt_keyword_list (param, options, maxopt, dictionary, errmsg)

char	param[ARB]		# List as comma separated string
int	options[ARB]		# List of indices found, with length
int	maxopt			# Size of options
char	dictionary[ARB]		# Dictionary string
char	errmsg[ARB]		# Error message template for bad value

pointer	sp, word, match, msg, scn
int	count, index
int	strdic()
bool	ok
bool	gt_gargcomma()
pointer	gt_sscan()

begin
	call smark (sp)
	call salloc (word, SZ_LINE, TY_CHAR)
	call salloc (match, SZ_LINE, TY_CHAR)

	count = 1
	scn = gt_sscan (param)

	repeat {
	    ok = gt_gargcomma (scn, Memc[word], SZ_LINE)
	    if (ok) {
		if (count >= maxopt) {
		    call error (GEM_ERR, "Too many keywords in list")
		}
		count = count+1
		index = strdic (Memc[word], Memc[match], SZ_LINE, dictionary)
		if (index <= 0) {
		    call salloc (msg, SZ_LINE, TY_CHAR)
		    call sprintf (Memc[msg], SZ_LINE, errmsg)
			call pargstr (Memc[word])
		    iferr (call error (GEM_ERR, Memc[msg])) {
			call sfree (sp)
			call gt_free_scan (scn)
			call erract (EA_ERROR)
		    }
		}
		options[count] = index
	    }
	} until (! ok)
	options[1] = count

	call sfree (sp)
	call gt_free_scan (scn)
end

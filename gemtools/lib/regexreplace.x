# Copyright(c) 2004-2006 Association of Universities for Research in Astronomy, Inc.
#

include	<error.h>
include	<pattern.h>
include	"../pkg/gemextn/gemerrors.h"


# Perform a string substitution that, hopefully, follows the same
# principles as the substitution in image specifications.

# gt_regex_replace - do a substitution
# gt_parse_replace - splits the specification into the four components
# gt_apply_replace - applies the four components to a string

# The specification has the form before%remove%insert%after where the
# % characters separate either regular expressions (before, remove and
# after) or text (insert).  The % characters can be escaped, if
# required.

# Substitution matches "before", then replaces "remove" with "insert".
# The "after" expression is an additional (optional) condition - if
# present, the substitution only occurs if "after" matches the text
# immediately after "remove".

# In addition, to make prepending directories easy in gemextn, the
# special character '^' in "before" identifies the start of the
# string. So gt_regex_replace ("abc", 999, destn, "^%%pre%") will give
# destn="preabc".


define	CH_EDIT		'%'	# See fio/fntgfn.x
define	CH_ESC		'\\'	# Allow escaping of % signs
define  CH_START        '^'     # Matches start of string in "begin"


# GT_REGEX_REPLACE -- Do substitution with regexps
# May raise exceptions

procedure gt_regex_replace (oldstring, len, newstring, spec)

char	oldstring[ARB]		# I Initial string
int	len			# I Length of output buffer
char	newstring[len]		# O Destination for substituted text
char	spec[ARB]		# I Describes the substitution

pointer	sp, before, remove, insert, after

errchk	gt_parse_replace, gt_apply_replace

begin
	call smark (sp)
	call salloc (before, SZ_FNAME, TY_CHAR)
	call salloc (remove, SZ_FNAME, TY_CHAR)
	call salloc (insert, SZ_FNAME, TY_CHAR)
	call salloc (after, SZ_FNAME, TY_CHAR)

	iferr {
	    call gt_parse_replace (spec, SZ_FNAME, before, remove, insert, 
		after)
	    call gt_apply_replace (oldstring, len, newstring, Memc[before],
		Memc[remove], Memc[insert], Memc[after])
	} then {
	    call sfree (sp)
	    call erract (EA_ERROR)
	}

	call sfree (sp)
end


# GT_PARSE_REPLACE -- Split before%remove%insert%after into separate buffers

procedure gt_parse_replace (replace, len, before, remove, insert, after)

char	replace[ARB]		# I xxx%foo%bar%yyy style string
int	len			# I Buffer sizes
pointer	before			# O Starting pattern (xxx)
pointer	remove			# O Text to remove from within pattern (foo)
pointer	insert			# O Text to replace in pattern (bar)
pointer	after			# O Ending pattern (yyy)

int	i, iin, iout, idest
int	strlen()
pointer	sp, destn
bool	escape

begin
	call smark (sp)
	call salloc (destn, 4, TY_POINTER)
	Memi[destn] = before
	Memi[destn+1] = remove
	Memi[destn+2] = insert
	Memi[destn+3] = after

	idest = 0
	iout = 0
	escape = false

	for (iin = 1; iin <= strlen (replace); iin = iin+1) {
	    if (replace[iin] == CH_EDIT) {
		Memc[Memi[destn+idest]+iout] = EOS
		idest = idest+1
		if (idest > 3) {
		    call sfree (sp)
		    call error (RPL_ERR,
			"Too many % signs in replace expression")
		}
		iout = 0
	    } else {
		if (iout+1 >= len) {
		    call sfree (sp)
		    call error (RPL_ERR, "Replace expression too long")
		} else {
		    if (replace[iin] == CH_ESC &&
			replace[iin+1] == CH_EDIT) {
			# hop over escape and treat % as ordinary character
			iin = iin+1
		    }
		    Memc[Memi[destn+idest]+iout] = replace[iin]
		    iout = iout+1
		}
	    }
	}

	for (i = idest; i < 4; i = i+1) {
	    Memc[Memi[destn+i]+iout] = EOS
	    iout = 0
	}

	call sfree (sp)
end


# GT_APPLY_REPLACE -- Do substitution

procedure gt_apply_replace (oldstring, len, newstring, before, remove, insert,
    after)

char	oldstring[ARB]		# I String to be modified
int	len			# I Size of newstring
char	newstring[len]		# O Modified string
char	before[ARB]		# I Pattern before subsn
char	remove[ARB]		# I Pattern to remove
char	insert[ARB]		# I Text to replace remove
char	after[ARB]		# I Pattern to match after

pointer	sp, pbefore, pafter, premove
int	size, count, iin, length, iout, skip, peek, more, i
bool	nobefore, noafter, noremove, firstormatch, debug
int	patmake(), patmatch(), pat_amatch(), strlen()

begin
        debug = false

        if (debug) {
            call eprintf ("gt_apply_replace %s,%s,%s,%s,%s\n")
            call pargstr (oldstring)
            call pargstr (before)
            call pargstr (remove)
            call pargstr (insert)
            call pargstr (after)
        }

	call smark (sp)
	call salloc (pbefore, SZ_LINE, TY_CHAR)
	call salloc (premove, SZ_LINE, TY_CHAR)
	call salloc (pafter, SZ_LINE, TY_CHAR)

	nobefore = before[1] == EOS
	noafter = after[1] == EOS
	noremove = remove[1] == EOS

        if (! nobefore && CH_START == before[1]) {
            # Drop initial marker
            size = patmake (before[2], Memc[pbefore], SZ_LINE)
        } else {
            size = patmake (before, Memc[pbefore], SZ_LINE)
        }
	size = patmake (remove, Memc[premove], SZ_LINE)
	size = patmake (after, Memc[pafter], SZ_LINE)
	length = strlen (oldstring)

	iin = 1
	iout = 1
	firstormatch = true

	while (iin <= length) {

	    # First, match and copy "before"
	    # If undefined, step through character by character,
	    # unless we're just starting or have just matched
	    if (nobefore) {
		if (firstormatch) {
		    firstormatch = false
		    count = 1
		} else {
		    count = 2
		}
	    } else {
                # Handle special start character
                if (CH_START == before[1]) {
                    if (1 == iin && firstormatch) {
                        if (before[2] == EOS) {
                            if (debug)
                                call eprintf ("single start\n")
                            count = 0
                        } else {
                            if (debug)
                                call eprintf ("complex start\n")
                            count = patmatch (oldstring[iin], Memc[pbefore])
                            if (count == 0)
                                count = length-iin+2
                        }
                    } else {
                        count = length-iin+2
                    }
                } else {
                    count = patmatch (oldstring[iin], Memc[pbefore])
                    if (count == 0)
                        count = length-iin+2
                }
	    }

	    for (i = 1; i < count; i = i+1) {

                if (debug)
                    call eprintf ("copying initial match\n")

		if (iout >= len) {
		    call sfree (sp)
		    call error (RPL_ERR, "Pattern too large")
		}
		newstring[iout] = oldstring[iin]
		iout = iout+1
		iin = iin+1
	    }

	    # next, see if "remove" is present
	    if (iin <= length) {

		if (noremove) {
		    skip = 0
		} else {
		    skip = pat_amatch (oldstring[iin], 1, Memc[premove])
		}

		if (noremove || skip > 0) {

		    # before replacing, check that "after" exists
		    peek = iin+skip
		    if (noafter) {
			more = 0
		    } else {
			more = pat_amatch (oldstring[peek], 1, Memc[pafter])
		    }

		    # ok, so copy across inserted text
		    if (noafter || more > 0) {

                        if (debug)
                            call eprintf ("copying replacement\n")

			for (i = 1; i <= strlen (insert); i = i+1) {
			    if (iout >= len) {
				call sfree (sp)
				call error (RPL_ERR, "Pattern too large")
			    }
			    newstring[iout] = insert[i]
			    iout = iout+1
			}
			iin = peek

			# avoid retesting if noremove
			firstormatch = ! noremove
		    }
		}
	    }
	}

	if (iout >= len) {
	    call sfree (sp)
	    call error (RPL_ERR, "Pattern too large")
	}
	newstring[iout] = EOS

	call sfree (sp)
end

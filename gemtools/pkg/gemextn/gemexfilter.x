# Copyright(c) 2004-2006 Association of Universities for Research in Astronomy, Inc.
#

include	<mach.h>
include	<error.h>
include	"gemextn.h"
include	"../../lib/imexplode.h"


# See gemextensions.x
#
# This file contains:
#
#                 gx_filter(...) - check the file against index, name, 
#                                  and version params
#
# Support routines:
#  bool gx_in_range(desc, value) - is value in the range described?
# bool gx_has_match(list, value) - is value matched by a pattern in the list?


# GX_FILTER -- GX_FILTER a file
# Throws failure exceptions from gx_verify via gx_check_and_out

procedure gx_filter (imx, gxn)

pointer	imx			# I The image to output
pointer	gxn			# IO Task parameters

bool	debug
bool	gx_in_range(), gx_has_match(), g_whitespace()

begin
	debug = false

	if (debug) {
	    call eprintf ("filter: ")
	    call imx_debug (imx)
	}

	if (! g_whitespace (GXN_INDEX(gxn)) &&
	    ! gx_in_range (GXN_INDEX(gxn), IMX_CLINDEX(imx)))
	    return

	if (! g_whitespace (GXN_EXTNAME(gxn)) &&
	    ! gx_has_match (GXN_EXTNAME(gxn), Memc[IMX_EXTNAME(imx)]))
	    return

	if (! g_whitespace (GXN_EXTVER(gxn)) &&
	    ! gx_in_range (GXN_EXTVER(gxn), IMX_EXTVERSION(imx)))
	    return

	call gx_check_and_out (imx, gxn)
end


# GX_IN_RANGE -- Check range against what's defined
# Returns true if within range, false if outside or no range

bool procedure gx_in_range (desc, value)

char	desc[ARB]			# Description of the range
int	value				# The value to check

pointer	sp, range
int	nvalues
int	ldecode_ranges()
bool	ok
bool	lis_in_range(), g_whitespace()

begin
	# g_whitespace needed for undefined ranges, 
	# which don't accept NO_INDEX
	ok = false
	if (! g_whitespace (desc)) {
	    call smark (sp)
	    call salloc (range, 3 * MAX_RANGES, TY_INT)
	    if (ERR !=
		ldecode_ranges (desc, Memi[range], MAX_RANGES, nvalues)) {
		ok = lis_in_range (Memi[range], value) || 0 == nvalues
	    }
	    call sfree (sp)
	}

	return ok
end


# GX_HAS_MATCH -- Check for presence in comma separated list of regexps
# Returns true if in list, false if missing or no list

bool procedure gx_has_match (list, value)

char	list[ARB]			# Comma separated list of patterns
char	value[ARB]			# Value to look for

bool	havepattern, havematch
pointer	sp, word, pattern, scn
bool	gt_gargcomma(), g_whitespace()
int	patmatch(), patmake()
pointer	gt_sscan()

begin
	havematch = false
	call smark (sp)
	call salloc (word, SZ_LINE, TY_CHAR)
	call salloc (pattern, SZ_LINE, TY_CHAR)
	scn = gt_sscan (list)

	havepattern = gt_gargcomma (scn, Memc[word], SZ_LINE)
	while (havepattern && ! havematch) {
	    if (! g_whitespace (Memc[word])) {
		if (ERR != patmake (Memc[word], Memc[pattern], SZ_LINE))
		    havematch = patmatch (value, Memc[pattern]) > 0
	    }
	    havepattern = gt_gargcomma (scn, Memc[word], SZ_LINE)
	}
	call gt_free_scan (scn)
	call sfree (sp)

	return havematch
end

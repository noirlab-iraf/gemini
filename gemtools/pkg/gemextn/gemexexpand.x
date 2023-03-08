# Copyright(c) 2004-2006 Association of Universities for Research in Astronomy, Inc.
#

include	<error.h>
include	<syserr.h>
include	"../../lib/imexplode.h"
include	"gemextn.h"


# See gemextensions.x
#
# This file contains:
#
#               gx_expand(...) - expand extensions by reading data from disk
#
# Support routines:
# bool gx_range_defined(range) - does range define at least one value?
#   bool gx_list_defined(list) - does list contain at least one value?


# GX_EXPAND -- Expand a file by accessing the data
# Throws failure exceptions via gx_check_and_out

procedure gx_expand (imx, gxn)

pointer	imx			# I The image to output
pointer	gxn			# IO Task parameters

pointer	sp, name, copy, mef, errmsg
pointer	imx_copy(), imx_mef_open(), mefgeti()
int	i, ver
bool	hasindex, hasver, hasname, special, debug
bool	gx_range_defined(), gx_list_defined()
bool	imx_h_clindex(), imx_h_extname(), imx_h_extver(), streq()
errchk	mefgstr, mefgeti

begin
	debug = false

	if (debug) {
	    call eprintf ("expand: ")
	    call imx_debug (imx)
	}

	call smark (sp)
	call salloc (name, SZ_FNAME, TY_CHAR)
	call salloc (errmsg, SZ_FNAME, TY_CHAR)
	hasindex = gx_range_defined (GXN_INDEX(gxn))
	hasver = gx_range_defined (GXN_EXTVER(gxn))
	hasname = gx_list_defined (GXN_EXTNAME(gxn))
	# special case - display everything available
	special = ! (hasindex || hasver || hasname) && 1 == GXN_NARGS(gxn)

	copy = imx_copy (imx)
	call imx_d_extname (copy)
	call imx_d_extver (copy)
	call imx_d_clindex (copy)
	ifnoerr (mef = imx_mef_open (copy, READ_ONLY, 0)) {

	    # try all indices in turn
	    for (i = 0; true; i = i+1) {

		if (debug) {
		    call eprintf ("index: %d\n")
			call pargi (i)
		}

		call imx_d_extname (copy)
		call imx_d_extver (copy)
		call imx_d_clindex (copy)

		if (debug)
		    call eprintf ("mef_rdhdr_gn\n")

		iferr (call imx_rdhdr_gn (mef, i)) {
		    if (debug)
			call eprintf ("imx_rdhdr_gn error\n")
		    if (0 == i)
			next
		    else
			break
		}

		call imx_s_clindex (copy, i)

		iferr {
		    if (debug)
			call eprintf ("mefgstr\n")
		    call mefgstr (mef, EXTNAME, Memc[name], SZ_FNAME)
		    call imx_s_extname (copy, Memc[name])
		} then { }

		iferr {
		    if (debug)
			call eprintf ("mefgeti\n")
		    ver = mefgeti (mef, EXTVERSION)
		    call imx_s_extver (copy, ver)
		} then { }

		if (imx_h_clindex (imx) && 
		    IMX_CLINDEX(imx) != IMX_CLINDEX(copy)) {
		    if (debug)
			call eprintf ("index given, doesn't match\n")
		    next
		}

		if (imx_h_extname (imx) &&
		    (! imx_h_extname (copy) ||
		     (! streq (Memc[IMX_EXTNAME(imx)],
			       Memc[IMX_EXTNAME(copy)])))) {
		    if (debug)
			call eprintf ("name given, doesn't match %s %s\n")
			    call pargstr (Memc[IMX_EXTNAME(imx)])
			    call pargstr (Memc[IMX_EXTNAME(copy)])
		    next
		}

		if (imx_h_extver (imx) &&
		    (! imx_h_extver (copy) ||
		     (IMX_EXTVERSION(imx) != IMX_EXTVERSION(copy)))) {
		    if (debug)
			call eprintf ("version given, doesn't match\n")
		    next
		}

		if (! special && ! hasindex)
		    call imx_d_clindex (copy)

		if (! special && ! hasname)
		    call imx_d_extname (copy)

		if (! special && ! hasver)
		    call imx_d_extver (copy)

		iferr (call gx_filter (copy, gxn))
		    call gx_report_fail (gxn)
	    }

	    call mefclose (mef)
	}

	call imx_free (copy)
	call sfree (sp)
end


# GX_RANGE_DEFINED -- Is range defined?
# Return true if at least one value defined in range

bool procedure gx_range_defined(range)

char	range[ARB]		# I The range string to test

pointer	sp, rng
bool	ok
bool	g_whitespace()
int	nvalues
int	ldecode_ranges()

begin
	if (g_whitespace (range)) {
	    ok = false
	} else {
	    call smark (sp)
	    call salloc (rng, 3 * MAX_RANGES, TY_INT)
	    if (ERR != ldecode_ranges (range, Memi[rng], MAX_RANGES, nvalues))
		ok = 0 != nvalues
	    call sfree (sp)
	}
	return ok
end


# GX_LIST_DEFINED -- Is a list defined?
# Return true if at least one value in the list is non-null

bool procedure gx_list_defined(list)

char	list[ARB]		# I The list to test

pointer	sp, word, scn
bool	ok
bool	gt_gargcomma(), g_whitespace()
pointer	gt_sscan()

begin
	if (g_whitespace (list))
	    return false

	ok = false
	call smark (sp)
	call salloc (word, SZ_LINE, TY_CHAR)
	Memc[word] = EOS

	scn = gt_sscan (list)
	while (!ok && Memc[word] == EOS) {
	    ok = gt_gargcomma (scn, Memc[word], SZ_LINE)
	}

	call gt_free_scan (scn)
	call sfree (sp)
	return ok
end


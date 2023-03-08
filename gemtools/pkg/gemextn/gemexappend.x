# Copyright(c) 2004-2006 Association of Universities for Research in Astronomy, Inc.
#

include	<mach.h>
include	<error.h>
include	"gemerrors.h"
include	"gemextn.h"


# See gemextensions.x
#
# This file contains:
#
#               gx_append(...) - append combinations of 
#                                extension params to file
#
# Support routines:
#    bool gx_append_index(...) - append indices
# bool gx_append_name_ver(...) - append name and version


# GX_APPEND -- Create files by appending extension information
# May throw fail exception via gx_check_and_out and other exceptions via
# gx_append_index and gx_append_name_ver

procedure gx_append (imx, gxn)

pointer	imx			# I The image to output
pointer	gxn			# IO Task parameters

bool	done
bool	imx_h_any_extension(), gx_append_index(), gx_append_name_ver()

errchk	gx_append_index, gx_append_name_ver, gx_check_and_out

begin
	if (imx_h_any_extension (imx)) {
	    call gx_check_and_out (imx, gxn)
	} else {
	    done = gx_append_index (imx, gxn)
	    if (! done)
		done = gx_append_name_ver (imx, gxn)
	}

	# TODO - if ! done here then we could issue an error about
	# having neither index nor name defined
end


# GX_APPEND_INDEX -- GX_APPEND the index ranges
# Returns true if some index values were tried

bool procedure gx_append_index (imx, gxn)

pointer	imx			# I The image to output
pointer	gxn			# IO Task parameters

pointer	sp, range
int	nvalues, i
int	lget_next_number(), ldecode_ranges()
bool	done
bool	lis_in_range(), g_whitespace()

begin
	done = false
	call smark (sp)
	call salloc (range, 3 * MAX_RANGES, TY_INT)

	if (! g_whitespace (GXN_INDEX(gxn))) {
	    if (ERR == ldecode_ranges (GXN_INDEX(gxn), Memi[range], MAX_RANGES,
		    nvalues)) {
		call sfree (sp)
		call error (GEM_ERR, "Bad range syntax for index.")
	    }

	    if (lis_in_range (Memi[range], MAX_INT)) {
		call sfree (sp)
		call error (GEM_ERR, "Open range for index")
	    }

	    i = -1
	    while (EOF != lget_next_number (Memi[range], i)) {
		done = true
		call imx_s_clindex (imx, i)
		iferr (call gx_check_and_out (imx, gxn)) {
		    call gx_report_fail (gxn)
		}
	    }
	}

	call sfree (sp)
	return done
end


# GX_APPEND_NAME_VER GX_APPEND name and, possibly, version
# Returns true if some specification was attempted

bool procedure gx_append_name_ver (imx, gxn)

pointer	imx			# I The image to output
pointer	gxn			# IO Task parameters

pointer	sp, range, word, scn
int	i, nvalues
int	ldecode_ranges(), lget_next_number()
bool	done, havename, havever
bool	gt_gargcomma(), lis_in_range(), g_whitespace()
pointer	gt_sscan()

begin
	done = false;
	call smark (sp)
	call salloc (word, SZ_LINE, TY_CHAR)
	call salloc (range, 3 * MAX_RANGES, TY_INT)
	scn = gt_sscan(GXN_EXTNAME(gxn))

	havename = gt_gargcomma (scn, Memc[word], SZ_LINE)
	while (havename) {
	    if (! g_whitespace (Memc[word])) {
		havever = false
		done = true
		# should we check for wildcards here?
		# i think not - let other routines complain if
		# a character is illegal
		call imx_s_ext_name (imx, Memc[word])
		if (! g_whitespace (GXN_EXTVER(gxn))) {
		    if (ERR == ldecode_ranges (GXN_EXTVER(gxn), Memi[range],
			    MAX_RANGES, nvalues)) {
			call gt_free_scan (scn)
			call sfree (sp)
			call error (GEM_ERR, "Bad range syntax for extversion")
		    }

		    if (lis_in_range (Memi[range], MAX_INT)) {
			call gt_free_scan (scn)
			call sfree (sp)
			call error (GEM_ERR, "Open range for extversion")
		    }

		    i = -1
		    while (EOF != lget_next_number (Memi[range], i)) {
			havever = true
			call imx_s_extver (imx, i)
			iferr (call gx_check_and_out (imx, gxn))
			    call gx_report_fail (gxn)
		    }
		}
		if (! havever) {
		    iferr (call gx_check_and_out (imx, gxn))
			call gx_report_fail (gxn)
		}
	    }
	    havename = gt_gargcomma (scn, Memc[word], SZ_LINE)
	}

	call gt_free_scan (scn)
	call sfree (sp)
	return done
end

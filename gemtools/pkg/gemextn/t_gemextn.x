# Copyright(c) 2004-2013 Association of Universities for Research in Astronomy, Inc.
#

include	<error.h>
include	<syserr.h>
include	<ctotok.h>
include	"gemerrors.h"
include	"gemextn.h"
include "../gemlog/glog.h"


# This file contains:
#
#                 t_gemextn() - the main task procedure
#    gx_raise_fail(file, msg) - raise failure exception
# gx_raise_fail_imx(imx, msg) - raise failure exception
#         gx_report_fail(gxn) - report failure
# 
# Support routines:
#          gx_open_null(file) - open file or return NULL
#         gx_parse_check(...) - parse the check parameter
#       gx_check_depends(...) - check dependncies between options
#        gx_inconsistent(...) - check dependncies between restricted options
#            gx_contains(...) - search array for value
#    gx_loop_over_inputs(...) - main loop over all files
#
# See also:
# gemexchkout.x - file spec checking and output
# gemexexpand.x - reading extensions from files on disk
# gemexfilter.x - filtering specs against index, extname and extver
# gemexappend.x - appending extensions to plain file names
# gemexverify.x - checking access conditions etc
# gemexomit.x   - omitting file spec components


define	CHECK_IMG	"img"
define	CHECK_EXT	"ext"
define	CHECK_AUTO	"auto"
define	EQUALS		"="
define	CHECK_SEPN	";"
define	COMMA		","	# Argument separator in lists

define	PROCESS		"|none|expand|filter|append|"
define	PROC_NONE	1	# No processing
define	PROC_EXPAND	2	# Expansion using disk file
define	PROC_FILTER	3	# Filter given text
define	PROC_APPEND	4	# Generate textual template

define	done_		99	# Error exit for main loop

define	FAILURE		"FAIL: %s - %s"

define	SZ_OPTIONS	10	# Size of option list


# T_GEMEXTN -- Expand a list of files into a list of image
# extensions, with verification, counting and renaming.

procedure t_gemextn ()

pointer	images			# I List of file names
pointer	check			# I Checks to be performed
# GXN_IMG_CHECK			# I Image check array
# GXN_EXT_CHECK			# I Extension check array
# GXN_AUTO_CHECK		# I Automatic check array
pointer	process			# I Processing to be done
# GXN_INDEX			# I Range list of extension indexes
# GXN_EXTNAME			# I Pattern for extension names
# GXN_EXTVER			# I Range list of extension versions
# GXN_IKPARAMS			# I Image kernel parameters
# GXN_OMIT			# I Format control
# GXN_REPLACE			# I Late modification
# GXN_FDOUT			# I Output stream for results
# GXN_NARGS			# I Number of command line arguments
# GXN_FAILCOUNT			# O Output count of failed checks
# GXN_COUNT			# O Output count of generated extensions

int	iprocess
pointer	gxn, sp, junk, omit, outfile, logfile
bool	verbose, debug
int	gx_open_null(), strdic(), clgeti()
pointer	gxn_alloc()
bool	clgetb()

# GEMLOG variables and functions
pointer	gl, pp, op, msgstr
bool	g_whitespace()
pointer	glogopen()

int	btoi()
pointer clopset()

errchk	gt_keyword_list(), gx_parse_check()

begin
        debug = false

        if (debug)
            call eprintf ("starting gemextn\n")

	call clputi ("status", GEM_ERR)
	call clputi ("count", 0)
	call clputi ("fail_count", 0)

	call smark (sp)

	call salloc (images, SZ_LINE, TY_CHAR)
	call salloc (check, SZ_LINE, TY_CHAR)
	call salloc (process, SZ_LINE, TY_CHAR)
	call salloc (omit, SZ_LINE, TY_CHAR)
	call salloc (outfile, SZ_LINE, TY_CHAR)
	call salloc (logfile, SZ_LINE, TY_CHAR)

	gxn = gxn_alloc ()

	call salloc (junk, SZ_LINE, TY_CHAR)
	call salloc (msgstr, SZ_LINE, TY_CHAR)

	# Task parameters
	call clgstr ("inimages", Memc[images], SZ_LINE)
	call clgstr ("check", Memc[check], SZ_LINE)
	call clgstr ("process", Memc[process], SZ_LINE)
	call clgstr ("index", GXN_INDEX(gxn), SZ_LINE)
	call clgstr ("extname", GXN_EXTNAME(gxn), SZ_LINE)
	call clgstr ("extver", GXN_EXTVER(gxn), SZ_LINE)
	call clgstr ("ikparams", GXN_IKPARAMS(gxn), SZ_LINE)
	call clgstr ("omit", Memc[omit], SZ_LINE)
	call clgstr ("replace", GXN_REPLACE(gxn), SZ_LINE)
	call clgstr ("outfile", Memc[outfile], SZ_LINE)
	call clgstr ("logfile", Memc[logfile], SZ_LINE)
	verbose = clgetb ("verbose")
	GXN_NARGS(gxn) = clgeti ("$nargs")
	

	# Do not force the use of a logfile, as gemextn is also useful as
	# a standalone tool.  Ignore GEMLOG support if logfile is empty.
	
	gl = NULL
	if ( ! g_whitespace (Memc[logfile]) ) {
	    
	    call opalloc (op)
	    GXN_P_OP(gxn) = op
	    
	    # Get GLOGPARS pset pointer  (closed in glogopen)
	    pp = clopset ("glogpars")
	    
	    # Open the logfile and assign the GL structure
	    OP_FL_APPEND(op) = YES
	    OP_FORCE_APPEND(op) = NO
	    OP_VERBOSE(op) = btoi (verbose)
	    iferr (gl = glogopen(Memc[logfile], "gemextn", "", pp, op)) {
		gl = NULL
		call opfree (op)
		call gt_log_error (gl, NULL)
		goto done_
	    }
	}
	GXN_P_GL(gxn) = gl

	# Open output stream
	GXN_FDOUT(gxn) = gx_open_null (Memc[outfile])

	# Parse/check list style parameters
	iferr {
	    iprocess = strdic (Memc[process], Memc[junk], SZ_LINE, PROCESS)
	    if (iprocess == NO_MATCH)
		call error (GEM_ERR, "Process value not recognised")
	    call gt_keyword_list (Memc[omit], GXN_OMIT(gxn), SZ_OPTIONS,
		OMIT, "Omit value not recognised (%s)")
	    call gx_parse_check (Memc[check], GXN_IMG_CHECK(gxn),
		GXN_EXT_CHECK(gxn), GXN_AUTO_CHECK(gxn), SZ_OPTIONS)
	} then {
	    call gt_log_error (gl, op)
	    goto done_
	}

	# Zero counters
	GXN_COUNT(gxn) = 0
	GXN_FAILCOUNT(gxn) = 0

	# Process images
	iferr (call gx_loop_over_inputs (Memc[images], iprocess, gxn)) {
	    call gt_log_error (gl, op)
	    goto done_
	}

	# Update counts and status
        if (debug)
            call eprintf ("updating counts\n")
	call clputi ("fail_count", GXN_FAILCOUNT(gxn))
	call clputi ("count", GXN_COUNT(gxn))
	if (GXN_FAILCOUNT(gxn) == 0)
	    call clputi ("status", 0)
	    
done_
        if (debug)
            call eprintf ("done\n")
	if (NULL != GXN_FDOUT(gxn))
	    call close (GXN_FDOUT(gxn))
	if (NULL != gl) {
	    call gl_close (gl)
	    call opfree (op)
	}
	call gxn_free (gxn)
	call sfree (sp)
        if (debug)
            call eprintf ("end of gemextn\n")
end


# GT_LOG_ERROR -- Note the error *somewhere* and set the status

procedure gt_log_error (gl, op)

pointer	gl			# I The gemlog pointer
pointer	op			# I Gemlog output pointer (used if gl != NULL)

pointer	sp, msg
int	status, junk

int	errget(), glogprint()

errchk	glogprint()

begin
	call smark (sp)
	call salloc (msg, SZ_LINE, TY_CHAR)

	status = errget (Memc[msg], SZ_LINE)
	call clputi ("status", status)
	    
	if (NULL != gl) {
	    OP_ERRNO(op) = 0
	    ifnoerr {
		junk = glogprint (gl, STAT_LEVEL, G_ERR_LOG, Memc[msg], op)
	    } then {
		call sfree (sp)
		return
	    }	
	}

	call fprintf (STDERR, "GEMEXTN ERROR %d %s\n")
	    call pargi (status)
	    call pargstr (Memc[msg])
	call sfree (sp)
	return
end


# GXN_ALLOC -- Allocate memory for the gxn structure

pointer procedure gxn_alloc ()

pointer	gxn

begin
	call malloc (gxn, LEN_GXN, TY_STRUCT)
	call malloc (GXN_P_IMG_CHECK(gxn), SZ_OPTIONS, TY_INT)
	call malloc (GXN_P_EXT_CHECK(gxn), SZ_OPTIONS, TY_INT)
	call malloc (GXN_P_AUTO_CHECK(gxn), SZ_OPTIONS, TY_INT)
	call malloc (GXN_P_INDEX(gxn), SZ_LINE, TY_CHAR)
	call malloc (GXN_P_EXTNAME(gxn), SZ_LINE, TY_CHAR)
	call malloc (GXN_P_EXTVER(gxn), SZ_LINE, TY_CHAR)
	call malloc (GXN_P_IKPARAMS(gxn), SZ_LINE, TY_CHAR)
	call malloc (GXN_P_OMIT(gxn), SZ_OPTIONS, TY_INT)
	call malloc (GXN_P_REPLACE(gxn), SZ_LINE, TY_CHAR)
	return gxn
end


# GXN_FREE -- Free memory allocated to the gxn structure

procedure gxn_free (gxn)

pointer	gxn				# IO The structure to free

begin
	call mfree (GXN_P_REPLACE(gxn), TY_CHAR)
	call mfree (GXN_P_OMIT(gxn), TY_INT)
	call mfree (GXN_P_IKPARAMS(gxn), TY_CHAR)
	call mfree (GXN_P_EXTVER(gxn), TY_CHAR)
	call mfree (GXN_P_EXTNAME(gxn), TY_CHAR)
	call mfree (GXN_P_INDEX(gxn), TY_CHAR)
	call mfree (GXN_P_AUTO_CHECK(gxn), TY_INT)
	call mfree (GXN_P_EXT_CHECK(gxn), TY_INT)
	call mfree (GXN_P_IMG_CHECK(gxn), TY_INT)
	call mfree (gxn, TY_STRUCT)

	gxn = NULL
end


# GX_OPEN_NULL -- Open a stream or return NULL

int procedure gx_open_null (file)

char	file[ARB]		# I The file description (or "")

int	open()

begin
	if (EOS == file[1]) {
	    return NULL
	} else {
	    return open (file, APPEND, TEXT_FILE)
	}
end


# GX_PARSE_CHECK -- Do the tricky parsing for the check parameter.
# Return the final results in the same format defined by gt_keyword_list.

procedure gx_parse_check (param, imgopts, extopts, autopts, maxopt)

char	param[ARB]		# I Check param (see help page for format)
int	imgopts[ARB]		# O List of image options found
int	extopts[ARB]		# O List of extension options found
int	autopts[ARB]		# O List of automatic options found
int	maxopt			# I Size of option arrays

pointer	sp, imgparam, extparam, autoparam, word, destn
int	token
bool	streq()

errchk	gt_keyword_list, gx_check_depends

begin
	call smark (sp)
	call salloc (word, SZ_LINE, TY_CHAR)
	call salloc (imgparam, SZ_LINE, TY_CHAR)
	call salloc (extparam, SZ_LINE, TY_CHAR)
	call salloc (autoparam, SZ_LINE, TY_CHAR)

	Memc[imgparam] = EOS
	Memc[extparam] = EOS
	Memc[autoparam] = EOS
	destn = autoparam
	token = TOK_UNKNOWN
	call sscan (param)

	while (token != TOK_EOS) {
	    call gargtok (token, Memc[word], SZ_LINE)

	    switch (token) {

	    case TOK_IDENTIFIER:

		# for img or ext, lookahead one token to see if we
		# have assignment.  if so, then set destn and
		# continue.  otherwise, add to current destn

		if (streq (Memc[word], CHECK_IMG)) {
		    call gargtok (token, Memc[word], SZ_LINE)
		    if (token == TOK_OPERATOR && streq (Memc[word], EQUALS)) {
			destn = imgparam
		    } else {
			call sfree (sp)
			call error (GEM_ERR,
			    "'img' not followed by '=' in check")
		    }
		} else if (streq (Memc[word], CHECK_EXT)) {
		    call gargtok (token, Memc[word], SZ_LINE)
		    if (token == TOK_OPERATOR && streq (Memc[word], EQUALS)) {
			destn = extparam
		    } else {
			call sfree (sp)
			call error (GEM_ERR,
			    "'ext' not followed by '=' in check")
		    }
		} else if (streq (Memc[word], CHECK_AUTO)) {
		    call gargtok (token, Memc[word], SZ_LINE)
		    if (token == TOK_OPERATOR && streq (Memc[word], EQUALS)) {
			destn = autoparam
		    } else {
			call sfree (sp)
			call error (GEM_ERR,
			    "'auto' not followed by '=' in check")
		    }
		} else {
		    call strcat (Memc[word], Memc[destn], SZ_LINE)
		    call strcat (COMMA, Memc[destn], SZ_LINE)
		}

	    case TOK_PUNCTUATION:

		# if we have a semicolon, revert back to default (img)
		# destination.  otherwise ignore commas and flag other
		# punctuation as an error (equals are slurped above).

		if (streq (Memc[word], CHECK_SEPN)) {
		    destn = autoparam;
		} else if (streq (Memc[word], COMMA)) {
		    # nothing
		} else {
		    call sfree (sp)
		    call error (GEM_ERR, "Unexpected punctuation in check")
		}

	    case TOK_EOS:

		# this is ok - will exit on while test

	    default:
		call sfree (sp)
		call error (GEM_ERR, "Syntax error in check")
	    }
	}

	# the above constructs separate strings for the different
	# test target.  now we can parse them in the usual way and
	# check for inconsistencies.

	iferr {
	    call gt_keyword_list (Memc[imgparam], imgopts, maxopt, IMG_CHECKS,
		"Image check parameter not recognised (%s)")
	    call gt_keyword_list (Memc[extparam], extopts, maxopt, EXT_CHECKS,
		"Extension check parameter not recognised (%s)")
	    call gt_keyword_list (Memc[autoparam], autopts, maxopt, 
		AUTO_CHECKS, "Automatic check parameter not recognised (%s)")

	    call gx_check_depends (imgopts, extopts, autopts)

	} then {
	    call sfree (sp)
	    call erract (EA_ERROR)
	}

	call sfree (sp)
end


# CHECK_DEPENDS -- Check inconsistent dependencies between all the
# different combinations of check parameters.
# Raises error on inconsistencies

procedure gx_check_depends (imgopts, extopts, autopts)

int	imgopts[ARB]		# I List of image options
int	extopts[ARB]		# I List of extension options
int	autopts[ARB]		# I List of automatic options

# Array sizes used in the tables below
# Convention is set, value, comparison set
# So AUT_AB_IMG are the img dependencies when auto=absent

define	AUT_AB_AUTO	6
define	AUT_AB_IMG	3
define	AUT_AB_EXT	4
define	AUT_IM_AUTO	1
define	AUT_IM_EXT	1
define	IMG_AB_AUTO	6
define	IMG_AB_IMG	3
define	IMG_AB_EXT	5
define	EXT_AB_AUTO	6
define	EXT_AB_EXT	4
define	EXT_IM_AUTO	1
define	EXT_IM_EXT	1

int	autabsent_badauto[AUT_AB_AUTO], autabsent_badimg[AUT_AB_IMG]
int	autabsent_badext[AUT_AB_EXT]
int	autimage_badauto[AUT_IM_AUTO], autimage_badext[AUT_IM_EXT]
int	imgabsent_badauto[IMG_AB_AUTO], imgabsent_badimg[IMG_AB_IMG]
int	imgabsent_badext[IMG_AB_EXT]
int	extabsent_badauto[EXT_AB_AUTO], extabsent_badext[EXT_AB_EXT]
int	extimage_badauto[EXT_IM_AUTO], extimage_badext[EXT_IM_EXT]

# here are the tables of values inconsistent with particular checks -
# array lengths must correspond to the definitions above

data	autabsent_badauto / AUTO_EMPTY, AUTO_EXISTS, AUTO_IMAGE, AUTO_MEF, AUTO_TABLE, AUTO_WRITE /
data	autabsent_badimg / IMG_EXISTS, IMG_MEF, IMG_WRITE /
data	autabsent_badext / EXT_EMPTY, EXT_EXISTS, EXT_IMAGE, EXT_TABLE /

data	autimage_badauto / AUTO_TABLE /
data	autimage_badext / EXT_TABLE /

data	imgabsent_badauto / AUTO_EMPTY, AUTO_EXISTS, AUTO_IMAGE, AUTO_MEF, AUTO_TABLE, AUTO_WRITE /
data	imgabsent_badimg / IMG_EXISTS, IMG_MEF, IMG_WRITE /
data	imgabsent_badext / EXT_ABSENT, EXT_EMPTY, EXT_EXISTS, EXT_IMAGE, EXT_TABLE /

data	extabsent_badauto / AUTO_EMPTY, AUTO_EXISTS, AUTO_IMAGE, AUTO_MEF, AUTO_TABLE, AUTO_WRITE /
data	extabsent_badext / EXT_EMPTY, EXT_EXISTS, EXT_IMAGE, EXT_TABLE /

data	extimage_badauto / AUTO_TABLE /
data	extimage_badext / EXT_TABLE /

errchk	gx_inconsistent
bool	gx_contains()

begin
	# since these are rather complex, allow a future work-round -
	# the option "force" for auto omits these checks
	if (gx_contains (autopts, AUTO_FORCE))
	    return

	call gx_inconsistent (AUTO_ABSENT, autopts, autabsent_badauto,
	    AUT_AB_AUTO, autopts, AUTO_CHECKS,
	    "check: auto=%s inconsistent with auto=absent")
	call gx_inconsistent (AUTO_ABSENT, autopts, autabsent_badimg,
	    AUT_AB_IMG, imgopts, IMG_CHECKS,
	    "check: img=%s inconsistent with auto=absent")
	call gx_inconsistent (AUTO_ABSENT, autopts, autabsent_badext,
	    AUT_AB_EXT, extopts, EXT_CHECKS,
	    "check: ext=%s inconsistent with auto=absent")

	call gx_inconsistent (AUTO_IMAGE, autopts, autimage_badauto,
	    AUT_IM_AUTO, autopts, AUTO_CHECKS,
	    "check: auto=%s inconsistent with auto=image")
	call gx_inconsistent (AUTO_IMAGE, autopts, autimage_badext,
	    AUT_IM_EXT, extopts, AUTO_CHECKS,
	    "check: ext=%s inconsistent with auto=image")

	call gx_inconsistent (IMG_ABSENT, imgopts, imgabsent_badauto,
	    IMG_AB_AUTO, autopts, AUTO_CHECKS,
	    "check: auto=%s inconsistent with img=absent")
	call gx_inconsistent (IMG_ABSENT, imgopts, imgabsent_badimg,
	    IMG_AB_IMG, imgopts, IMG_CHECKS,
	    "check: img=%s inconsistent with img=absent")
	call gx_inconsistent (IMG_ABSENT, imgopts, imgabsent_badext,
	    IMG_AB_EXT, extopts, EXT_CHECKS,
	    "check: ext=%s inconsistent with img=absent")

	call gx_inconsistent (EXT_ABSENT, extopts, extabsent_badauto,
	    EXT_AB_AUTO, autopts, AUTO_CHECKS,
	    "check: auto=%s inconsistent with ext=absent")
	call gx_inconsistent (EXT_ABSENT, extopts, extabsent_badext,
	    EXT_AB_EXT, extopts, EXT_CHECKS,
	    "check: ext=%s inconsistent with ext=absent")

	call gx_inconsistent (EXT_IMAGE, extopts, extimage_badauto,
	    EXT_IM_AUTO, autopts, AUTO_CHECKS,
	    "check: auto=%s inconsistent with ext=image")
	call gx_inconsistent (EXT_IMAGE, extopts, extimage_badext,
	    EXT_IM_EXT, extopts, EXT_CHECKS,
	    "check: ext=%s inconsistent with ext=image")
end


# GX_INCONSISTENT -- Flag an error for inconsistent check parameters.
# This processes a single set of dependencies, between one value and
# the options selected for a particular target (image, extension or
# automatic)
# Raises an error on inconsistencies
# msg must contain a "%s" where the conflicting value will be inserted

procedure gx_inconsistent (value, valopts, errors, length, erropts, dict, msg)

int	value			# I Value that is inconsistent with others
int	valopts[ARB]		# I Array in which value may be present
int	errors[ARB]		# I Array of values inconsistent with value
int	length			# I Number of error values
int	erropts[ARB]		# I Array in which error values may be present
char	dict[ARB]		# I The original dictionary, for error msg
char	msg[ARB]		# I Error message template

int	i
bool	gx_contains()
pointer	sp, word, errmsg

begin
	if (gx_contains (valopts, value)) {
	    for (i = 1; i < length; i = i+1) {
		if (gx_contains (erropts, errors[i])) {
		    call smark (sp)
		    call salloc (word, SZ_LINE, TY_CHAR)
		    call salloc (errmsg, SZ_LINE, TY_CHAR)
		    call gt_dicstr (errors[i], Memc[word], SZ_LINE, dict)
		    call sprintf (Memc[errmsg], SZ_LINE, msg)
			call pargstr (Memc[word])
		    iferr (call error (GEM_ERR, Memc[errmsg])) {
			call sfree (sp)
			call erract (EA_ERROR)
		    } else {
			call sfree (sp)
		    }
		}
	    }
	}
end


# GX_CONTAINS - search for int in array of values (first element is length)
# Returns true if list contains value

bool procedure gx_contains (list, value)

int	list[ARB]			# I List to search
int	value			# I Value to look for

int	i

begin
	for (i = 2; i <= list[1]; i = i+1) {
	    if (list[i] == value)
		return true
	}
	return false
end


# GX_LOOP_OVER_INPUTS -- The main loop over the input list

procedure gx_loop_over_inputs (images, process, gxn)

char	images[ARB]		# I List of file names
int	process			# I Process case
pointer	gxn			# IO Task parameters

pointer	sp, imt, image, imx
pointer	gemtopen(), imx_explod()
int	gemtgetim()
bool	debug

errchk	gemtopen(), gemtgetim()

begin
	debug = false
	call smark (sp)
	call salloc (image, SZ_FNAME, TY_CHAR)

	if (debug) {
	    call eprintf ("opening %s\n")
	    call pargstr (images)
	}
	imt = gemtopen (images)
        if (debug)
            call eprintf ("opened\n")

	while (gemtgetim (imt, Memc[image], SZ_FNAME) != EOF) {

            if (debug) {
                call eprintf ("processing %s\n")
                call pargstr (Memc[image])
            }

	    imx = imx_explod (Memc[image])
	    iferr {
		switch (process) {
		case PROC_NONE:
		    call gx_check_and_out (imx, gxn)
		case PROC_EXPAND:
		    call gx_expand (imx, gxn)
		case PROC_FILTER:
		    call gx_filter (imx, gxn)
		case PROC_APPEND:
		    call gx_append (imx, gxn)
		default:
		    call error (GEM_ERR, "Error in switch: loopoverinputs")
		}
	    } then {
                if (debug)
                    call eprintf ("error handling\n")
		call imx_free (imx)
		call gx_report_fail (gxn)
	    } else {
                if (debug)
                    call eprintf ("loop free\n")
		call imx_free (imx)
	    }
	}

        if (debug)
            call eprintf ("loop closing\n")
	call gemtclose (imt)
	call sfree (sp)
        if (debug)
            call eprintf ("loop done\n")
end


# GX_RAISE_FAIL -- Generate an error for later reporting

procedure gx_raise_fail (file, msg)

char	file[ARB]		# I The file spec that failed
char	msg[ARB]		# I Description of error

pointer	sp, text
int	length
int	strlen()

begin
	call smark (sp)
	length = strlen (file) + strlen (msg) + strlen (FAILURE) + 1
	call salloc (text, length, TY_CHAR)
	call sprintf (Memc[text], length, FAILURE)
	    call pargstr (file)
	    call pargstr (msg)
	iferr (call error (GEMEX_FAIL, Memc[text])) {
	    call sfree (sp)
	    call erract (EA_ERROR)
	} else {
	    call sfree (sp)
	}
end


# GX_RAISE_FAIL_IMX -- Generate an error for later reporting (imx arg)

procedure gx_raise_fail_imx (imx, msg)

pointer	imx			# I The file spec that failed
char	msg[ARB]		# I Description of error

pointer	sp, file

begin
	call smark (sp)
	call salloc (file, SZ_FNAME, TY_CHAR)
	call imx_condense (imx, Memc[file], SZ_FNAME)
	iferr (call gx_raise_fail (Memc[file], msg)) {
	    call sfree (sp)
	    call erract (EA_ERROR)
	} else {
	    call sfree (sp)
	}
end


# GX_REPORT_FAIL -- Notify the user about a failure

procedure gx_report_fail (gxn)

pointer	gxn			# IO Task parameters

pointer	sp, msg, gl, op
int	code
int	errcode(), errget()

int	glogprint()

begin
	gl = GXN_P_GL(gxn)
	op = GXN_P_OP(gxn)
	
	if (GEMEX_FAIL != errcode ()) {
	    call erract (EA_ERROR)
	} else {
	    GXN_FAILCOUNT(gxn) = GXN_FAILCOUNT(gxn) + 1
	    call smark (sp)
	    call salloc (msg, SZ_LINE, TY_CHAR)
	    code = errget (Memc[msg], SZ_LINE)
	    
	    if (NULL != gl) {		# GEMLOG mode
		OP_ERRNO(op) = 0
		code = glogprint (gl, STAT_LEVEL, G_ERR_LOG, Memc[msg], op)
	    } else {			# No GEMLOG mode  ("logfile" == "")
		call fprintf (STDERR, "%s\n")
		    call pargstr (Memc[msg])
		call flush (STDERR)
	    }
	    
	    call sfree (sp)
	}
end

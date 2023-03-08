# Copyright(c) 2004-2013 Association of Universities for Research in Astronomy, Inc.
#

include	<mach.h>
include	<error.h>
include	<syserr.h>
include	<ctotok.h>
include	<mef.h>
include	"../pkg/gemextn/gemerrors.h"
include	"imexplode.h"


# Separate a file specification into small pieces.  This doesn't claim
# to be as fancy as the real thing.  In particular, extension name and
# version identification may differ from fits routines.

# The handling of image kernel parameters is not perfect.
# - possible incorrect handling of abbreviated extname or extversion
# - possible rearrangement and expansion of parameters
# For example:
# imx_condense (imx_explod ("foo[bar=gtr,extver=3]")) 
#    == "foo[extversion=3,bar=gtr]"

# This file contains:
# 
#                  pointer imx_explod(image) - parse image and generate imx
#                       pointer imx_copy(imx) - copy an imx structure
#                               imx_free(imx) - free an imx and associated data
#            imx_condense(imx, value, length) - imx to string (all components)
#           imx_full_file(imx, value, length) - imx to string (file only)
#     pointer imx_mef_open(imx, acmode, oldp) - open using meflib
# pointer imx_imio_open(imx, acmode, hdr_arg) - open using meflib
#                              imx_debug(imx) - print details to stdout
#                            imx_detail(imx) - more debug details
#
# Support routines:
#                         pointer imx_alloc() - allocate imx
#                   imx_copy_string(p, value) - allocate and copy a string
#                         imx_path_split(...) - parse file description
#                       imx_parse_params(...) - parse kernel parameters
#                        bool imx_abrevn(...) - test for abbreviation
#            imx_intcat(value, destn, maxlen) - gx_append an int to a string
#                       imx_rdhdr_gn(mef, gn) - simplify mef_rdhdr_gn iface
#                              imx_rdhdr(...) - simplify mef_rdhdr iface

# See also:
# imexplode.h - details of imx structure
# imxdelete.x - routines to delete components
# imxhas.x    - routines to test for presence of components
# imxset.x    - routines to set components


define	FILE_SEP	"/"		# TODO - need cross-platform value?
define	CLVAR_SEP	"$"
define	EXTN_SEP	"."
define	CLUSTER_SEP	"/"
define	IKPARAMS_SEP	","
define	IKPARAMS_EQ	"="
define	OPEN_SQ		"["
define	CLOSE_SQ	"]"

define	MIN_MATCH	4

define	mef_open_	10
define	imio_open_	20

# IMX_ALLOC -- Allocate a structure with no values
# Returns a pointer to the structure

pointer procedure imx_alloc ()

pointer	imx

begin
	call malloc (imx, LEN_IMEXPLODE, TY_STRUCT)
	IMX_PATH(imx) = NULL
	IMX_FILE(imx) = NULL
	IMX_EXTENSION(imx) = NULL
	IMX_EXTNAME(imx) = NULL
	IMX_EXTVERSION(imx) = NULL
	IMX_IKPARAMS(imx) = NULL
	IMX_SECTION(imx) = NULL
	IMX_CLSIZE(imx) = NO_INDEX
	IMX_CLINDEX(imx) = NO_INDEX

	return imx
end


# IMX_COPY_STRING -- Allocate and copy memory for a string

procedure imx_copy_string (p, value)

pointer	p			# O Where the string is stored
char	value[ARB]		# I The string to copy

int	length
int	strlen()

begin
	length = strlen (value)
	call malloc (p, length, TY_CHAR)
	call strcpy (value, Memc[p], length)
end


# IMX_COPY -- Copy an imexplode structure
# Returns a distinct copy

pointer procedure imx_copy (imx)

pointer	imx			# I The structure to copy

pointer	copy
pointer	imx_alloc()

begin
	copy = imx_alloc ()

	if (NULL != IMX_PATH(imx))
	    call imx_s_path (copy, Memc[IMX_PATH(imx)])
	if (NULL != IMX_FILE(imx))
	    call imx_s_file (copy, Memc[IMX_FILE(imx)])
	if (NULL != IMX_EXTENSION(imx))
	    call imx_s_extension (copy, Memc[IMX_EXTENSION(imx)])
	call imx_s_clindex (copy, IMX_CLINDEX(imx))
	call imx_s_clsize (copy, IMX_CLSIZE(imx))
	if (NULL != IMX_EXTNAME(imx))
	    call imx_s_extname (copy, Memc[IMX_EXTNAME(imx)])
	call imx_s_extver (copy, IMX_EXTVERSION(imx))
	if (NULL != IMX_IKPARAMS(imx))
	    call imx_s_ikparams (copy, Memc[IMX_IKPARAMS(imx)])
	if (NULL != IMX_SECTION(imx))
	    call imx_s_section (copy, Memc[IMX_SECTION(imx)])

	return copy
end


# IMX_FREE -- Free memory related to an imexplode structure

procedure imx_free (imx)

pointer	imx			# I The pointer to the structure

begin
	if (imx != NULL) {
	    call imx_d_path (imx)
	    call imx_d_file (imx)
	    call imx_d_extension (imx)
	    call imx_d_extname (imx)
	    call imx_d_extver (imx)
	    call imx_d_ikparams (imx)
	    call imx_d_section (imx)
	    call mfree (imx, TY_STRUCT)
	    imx = NULL
	}
end


# IMX_CONDENSE -- Generate a file specification from the structure

procedure imx_condense (imx, value, length)

pointer	imx			# The structure to convert
char	value[ARB]		# Where to store the string
int	length			# Allocated size of string

bool	debug
bool	imx_h_clindex()
bool	imx_h_extname(), imx_h_extver(), imx_h_ikparams(), imx_h_section()

begin
	debug = false

	if (debug) {
	    call eprintf ("condense:\n")
	    call imx_detail (imx)
	}

	value[1] = EOS
	call imx_full_file (imx, value, length)

	if (imx_h_clindex (imx)) {
	    
	    call strcat (OPEN_SQ, value, length)
	    call imx_intcat (IMX_CLINDEX(imx), value, length)

	    # ugly hack due to strange size problem
	    if (IMX_CLSIZE(imx) > 0 && IMX_CLSIZE(imx) < MAX_INT) {
		call strcat (CLUSTER_SEP, value, length)
		call imx_intcat (IMX_CLSIZE(imx), value, length)
	    }

	    # Include kernel params here if no name
	    if (! imx_h_extname (imx)) {
		if (imx_h_extver (imx)) {
		    call strcat (IKPARAMS_SEP, value, length)
		    call strcat (EXTVERSION, value, length)
		    call strcat (IKPARAMS_EQ, value, length)
		    call imx_intcat (IMX_EXTVERSION(imx), value, length)
		}		
		if (imx_h_ikparams (imx)) {
		    call strcat (IKPARAMS_SEP, value, length)
		    call strcat (Memc[IMX_IKPARAMS(imx)], value, length)
		}
	    }
	    call strcat (CLOSE_SQ, value, length)

	}

	if (imx_h_extname (imx)) {

	    call strcat (OPEN_SQ, value, length)
	    call strcat (Memc[IMX_EXTNAME(imx)], value, length)
	    if (imx_h_extver (imx) || imx_h_ikparams (imx)) {
		call strcat (IKPARAMS_SEP, value, length)
	    }

	    if (imx_h_extver (imx)) {
		if (! imx_h_extname (imx)) {
		    call strcat (EXTVERSION, value, length)
		    call strcat (IKPARAMS_EQ, value, length)
		}
		call imx_intcat (IMX_EXTVERSION(imx), value, length)
		if (imx_h_ikparams (imx)) {
		    call strcat (IKPARAMS_SEP, value, length)
		}
	    }

	    if (imx_h_ikparams (imx)) {
		call strcat (Memc[IMX_IKPARAMS(imx)], value, length)
	    }
	    call strcat (CLOSE_SQ, value, length)

	}

	if (! imx_h_clindex (imx) && ! imx_h_extname (imx)  &&
	    (imx_h_extver (imx) || imx_h_ikparams (imx))) {

	    if (debug) call eprintf ("no index or name\n")

	    call strcat (OPEN_SQ, value, length)
	    if (imx_h_extver (imx)) {
		call strcat (EXTVERSION, value, length)
		call strcat (IKPARAMS_EQ, value, length)
		call imx_intcat (IMX_EXTVERSION(imx), value, length)
		if (imx_h_ikparams (imx)) {
		    call strcat (IKPARAMS_SEP, value, length)
		}
	    }		
	    if (imx_h_ikparams (imx)) {
		call strcat (Memc[IMX_IKPARAMS(imx)], value, length)
	    }
	    call strcat (CLOSE_SQ, value, length)

	}


	if (imx_h_section (imx)) {
	    call strcat (Memc[IMX_SECTION(imx)], value, length)
	}
end


# IMX_FULL_FILE -- Generate the file (without extensions) from the structure

procedure imx_full_file (imx, value, length)

pointer	imx			# I The structure to convert
char	value[ARB]		# O Where to store the string
int	length			# I Allocated size of string

bool	imx_h_path(), imx_h_extension()

begin
	value[1] = EOS
	if (imx_h_path (imx)) {
	    call strcat (Memc[IMX_PATH(imx)], value, length)
	}
	call strcat (Memc[IMX_FILE(imx)], value, length)
	if (imx_h_extension (imx)) {
	    call strcat (EXTN_SEP, value, length)
	    call strcat (Memc[IMX_EXTENSION(imx)], value, length)
	}
end


# IMX_EXPLOD -- Parse a file descriptor in string form, returning the 
# imx structure
# Returns imx structure or NULL

pointer procedure imx_explod (image)

char	image[ARB]		# I The file to parse

bool	debug
pointer	sp, imx, imgfile, imgkernel, section, path, file, extension
pointer	ikparams, extname
int	clindex, clsize, extversion
pointer	imx_alloc()

begin
	debug = false

	imx = NULL
	call smark (sp)

	call salloc (imgfile, SZ_FNAME, TY_CHAR)
	call salloc (imgkernel, SZ_FNAME, TY_CHAR)
	call salloc (section, SZ_FNAME, TY_CHAR)
	call salloc (path, SZ_FNAME, TY_CHAR)
	call salloc (file, SZ_FNAME, TY_CHAR)
	call salloc (extension, SZ_FNAME, TY_CHAR)
	call salloc (ikparams, SZ_FNAME, TY_CHAR)
	call salloc (extname, SZ_FNAME, TY_CHAR)

	call imparse (image, Memc[imgfile], SZ_FNAME, Memc[imgkernel],
	    SZ_FNAME, Memc[section], SZ_FNAME, clindex, clsize)
	if (debug) {
	    call eprintf ("imparse (%s)\n")
		call pargstr (image)
	    call eprintf ("file: %s\n")
		call pargstr (Memc[imgfile])
	    call eprintf ("kernel: %s\n")
		call pargstr (Memc[imgkernel])
	    call eprintf ("section: %s\n")
		call pargstr (Memc[section])
	    call eprintf ("index, size: %d %d \n")
		call pargi (clindex)
		call pargi (clsize)
	}
	call imx_path_split (Memc[imgfile], Memc[path], Memc[file],
	    Memc[extension], SZ_FNAME)
	call imx_parse_params (Memc[imgkernel], Memc[extname], extversion,
	    Memc[ikparams], clindex, clsize, SZ_FNAME)

	imx = imx_alloc ()
	call imx_s_path (imx, Memc[path])
	call imx_s_file (imx, Memc[file])
	call imx_s_extension (imx, Memc[extension])
	call imx_s_clindex (imx, clindex)
	call imx_s_clsize (imx, clsize)
	call imx_s_extname (imx, Memc[extname])
	call imx_s_extver (imx, extversion)
	call imx_s_ikparams (imx, Memc[ikparams])
	call imx_s_section (imx, Memc[section])

	call sfree (sp)
	return imx
end


# IMX_PATH_SPLIT -- Split a file into its components

procedure imx_path_split (all, path, file, extension, maxlen)

char	all[ARB]		# I Full file name
char	path[ARB]		# O Path name
char	file[ARB]		# O File name
char	extension[ARB]		# O Extension name
int	maxlen			# I Size of path, file and extensions

int	slash, dot, length, i
int	strldx(), strlen

begin
	path[1] = EOS
	file[1] = EOS
	extension[1] = EOS
	length = strlen (all)

	# handle dev/pix and dev$pix in the same way
	slash = max(strldx (FILE_SEP[1], all), strldx (CLVAR_SEP[1], all)) 
	dot = strldx (EXTN_SEP[1], all)
	if (dot < slash) 
	    dot = 0		# dot in path, not file
	if (dot == 0)
	    dot = length+1	# add dot after
	if (slash > maxlen || length - dot > maxlen)
	    return

	for (i = 1; i <= length+1; i = i+1) {
	    if (i <= slash) {
		path[i] = all[i]
		if (i == slash) 
		    path[i+1] = EOS
	    } else if (i < dot) {
		file[i-slash] = all[i]
	    } else if (i == dot) {
		file[i-slash] = EOS
	    } else if (i > dot && i <= length) {
		extension[i-dot] = all[i]
	    } else if (i > dot) {
		extension[i-dot] = EOS
	    }
	}
end


# IMX_PARSE_PARAMS --- A possibly incomplete parsing of image kernel 
#                      parameters.
# We might want to call imio/iki/fxf/fxfksection instead (if we can
# get hold of a fits structure pointer)

procedure imx_parse_params (all, extname, extversion, ikparams, clindex, 
    clsize, maxlen)

char	all[ARB]		# I Kernel params as string
char	extname[ARB]		# O Destination for name
int	extversion		# O Destination for version
char	ikparams[ARB]		# O Remaining params
int	clindex			# IO Extension index
int	clsize			# IO Extension size
int	maxlen			# I Length of extname and ikparams

int	length, ignored, token, start
int	strlen(), ctoi()
bool	first, haveparam, debug
bool	streq(), imx_abrevn()
pointer	sp, name, peek, errmsg

begin
	debug = false

	extname[1] = EOS
	ikparams[1] = EOS
	extversion = NO_INDEX
	length = strlen (all)
	if (length < 1)
	    return

	if (debug)
	    call eprintf ("imx_parse_params: %s\n")
		call pargstr (all)

	call smark (sp)
	call salloc (name, length+1, TY_CHAR)
	call salloc (peek, length+1, TY_CHAR)
	Memc[name] = EOS

	call sscan (all)
	call gargtok (token, Memc[peek], length)
	if (! streq (Memc[peek], OPEN_SQ)) {
	    call sfree (sp)
	    call error (IMX_ERR, "Image kernel without '['")
	}

	first = true
	haveparam = false

	repeat {
	    call gargtok (token, Memc[peek], length)
	    if (Memc[peek] == EOS)
		break

	    if (debug)
		call eprintf ("token %s\n")
		    call pargstr (Memc[peek])

	    switch (token) {

	    case TOK_PUNCTUATION:
		# if we have a name, then it's either an isolated
		# parameter or extname
		if (EOS != Memc[name]) {
		    if (first) {
			call strcpy (Memc[name], extname, maxlen)
		    } else {
			if (haveparam) {
			    call strcat (IKPARAMS_SEP, ikparams, maxlen)
			}
			call strcat (Memc[name], ikparams, maxlen)
			haveparam = true
		    }
		}
		Memc[name] = EOS
		first = false

	    case TOK_OPERATOR:
		# we should have either "name=value" or "name+/-"
		# unless the name corresponds to extname or
		# extver.

		# extname candidate - only use if we don't already
		# have a value
		if (imx_abrevn (Memc[name], EXTNAME, MIN_MATCH) &&
		    streq (Memc[peek], IKPARAMS_EQ) && EOS == extname[1]) {
		    call gargtok (token, Memc[peek], length)
		    if (token != TOK_PUNCTUATION) {
			call strcpy (Memc[peek], extname, maxlen)
			if (debug) {
			    call eprintf ("extname found %s\n")
				call pargstr (extname)
			} 
		    }

		# extversion candidate - if the value turns out to
		# not be numeric then we'll save it as a parameter
		} else if (imx_abrevn (Memc[name], EXTVERSION, MIN_MATCH) &&
		    streq (Memc[peek], IKPARAMS_EQ) &&
		    extversion == NO_INDEX) {
		    call gargtok (token, Memc[peek], length)
		    if (token != TOK_PUNCTUATION) {
			if (token == TOK_NUMBER) {
			    start = 1
			    ignored = ctoi (Memc[peek], start, extversion)
			    if (debug) {
				call eprintf ("extversion found %d\n")
				    call pargi (extversion)
			    }
			} else {
			    if (haveparam) {
				call strcat (IKPARAMS_SEP, ikparams, maxlen)
			    }
			    call strcat (Memc[name], ikparams, maxlen)
			    if (debug) {
				call eprintf ("param %s (not version)\n")
				    call pargstr (Memc[name])
			    }
			    call strcat (IKPARAMS_EQ, ikparams, maxlen)
			    call strcat (Memc[peek], ikparams, maxlen)
			    if (debug) {
				call eprintf ("value? %s\n")
				    call pargstr (Memc[peek])
			    }
			    haveparam = true
			}
		    }

		# the simple parameter case
		} else {
		    if (haveparam) {
			call strcat (IKPARAMS_SEP, ikparams, maxlen)
		    }
		    call strcat (Memc[name], ikparams, maxlen)
		    if (debug) {
			call eprintf ("param %s\n")
			    call pargstr (Memc[name])
		    }
		    call strcat (Memc[peek], ikparams, maxlen)
		    if (streq (Memc[peek], IKPARAMS_EQ)) {
			call gargtok (token, Memc[peek], length)
			if (token != TOK_PUNCTUATION) {
			    call strcat (Memc[peek], ikparams, maxlen)
			    if (debug) {
				call eprintf ("value? %s\n")
				    call pargstr (Memc[peek])
			    }
			}
		    } else {
			if (debug) {
			    call eprintf ("+/-? %s\n")
				call pargstr (Memc[peek])
			}
		    }
		    haveparam = true
		}

		Memc[name] = EOS
		first = false

	    case TOK_NUMBER:
		# this should be extversion, or index - assume index if
		# it comes before a name
		if (EOS == extname[1]) {
		    start = 1
		    ignored = ctoi (Memc[peek], start, clindex)
		    if (debug) {
			call eprintf ("index found in kernel params %d\n")
			    call pargi (clindex)
		    }
		} else {
		    start = 1
		    ignored = ctoi (Memc[peek], start, extversion)
		    if (debug) {
			call eprintf ("extversion found %d (no keyword)\n")
			    call pargi (extversion)
		    }
		}

	    default:
		# treat *everything* else a an identifier
		# (again, leave complaints to the FITS layer)

		# save as name, and interpret using next token
		# unless, of course, we already have a name - then
		# it's an error (can't just throw stuff away)
		if (EOS == Memc[name]) {
		    call strcpy (Memc[peek], Memc[name], length)
		} else {
		    call salloc (errmsg, SZ_LINE, TY_CHAR)
		    call sprintf (Memc[errmsg], SZ_LINE,
			"Syntax error in kernel parameters: %s")
			call pargstr (all)
		    iferr (call error (IMX_ERR, Memc[errmsg])) {
			call sfree (sp)
			call erract (EA_ERROR)
		    }
		}
	    }

	} until (streq (Memc[peek], CLOSE_SQ))

	if (Memc[peek] == EOS) {
	    call sfree (sp)
	    call error (IMX_ERR, "Image kernel without ']'")
	}

	call sfree (sp)
end


# IMX_ABREVN -- Check whether a string is an abbreviation of a certain text

bool procedure imx_abrevn (str, text, minlen)

char	str[ARB]		# I The possible abbreviation
char	text[ARB]		# I The full text
int	minlen			# I Minimum length required for match

int	i, length
int	strlen()

begin
	length = strlen (str)
	if (length < minlen)
	    return false
	if (length > strlen (text))
	    return false
	for (i = 1; i < length; i = i+1) {
	    if (str[i] != text[i])
		return false
	}
	return true
end


# IMX_INTCAT -- GX_APPEND an int

procedure imx_intcat (value, destn, maxlen)

int	value			# I What to gx_append
char	destn[ARB]		# IO Appended to
int	maxlen			# I maximum length of destn

pointer	sp, text

begin
	call smark (sp)
	call salloc (text, SZ_LINE, TY_CHAR)
	call itoc (value, Memc[text], SZ_LINE)
	call strcat (Memc[text], destn, maxlen)
	call sfree (sp)
end


# IMX_MEF_OPEN -- Open a mef structure from meflib
# This tries adding .fit and .fits extensions if no extension supplied
# It attempts to locate the extension given by index/name/version if given
# Raises errors on failure or inconsistency in index/name/version

pointer	procedure imx_mef_open (imx, acmode, oldp)

pointer	imx			# I The structure to open
int	acmode			# I Meflib access mode
pointer	oldp			# IO Meflib reserved value

pointer	sp, mef, errmsg, extname
int	code, index, version
bool	debug
pointer	errget(), imx_mef_file()
bool	imx_h_clindex(), imx_h_extname(), imx_h_extver()
int	strcmp(), mef_geti()

begin
	debug = false

	call smark (sp)
	call salloc (extname, SZ_LINE, TY_CHAR)
	call salloc (errmsg, SZ_LINE, TY_CHAR)

	if (debug) {
	    call eprintf ("imx_mef_open\n")
	    call imx_debug (imx)
	}

	ifnoerr (mef = imx_mef_file (imx, acmode, oldp)) {

	    # Try to match/select extension information

	    if (imx_h_clindex (imx)) {

		if (debug) {
		    call eprintf ("selecting index %d\n")
			call pargi (IMX_CLINDEX(imx))
		}

		iferr (call imx_rdhdr_gn (mef, IMX_CLINDEX(imx))) {
		    if (debug)
			call eprintf ("index failed\n")
		    call sfree (sp)
		    call mefclose (mef)
		    call erract (EA_ERROR)
		}

		if (imx_h_extname (imx)) {
		    call mef_gstr (mef, "EXTNAME", Memc[extname], SZ_LINE)
		    if (0 != strcmp (Memc[IMX_EXTNAME(imx)], Memc[extname])) { 
			if (debug) {
			    call eprintf ("inconsistent names: '%s' '%s'\n")
				call pargstr (Memc[IMX_EXTNAME(imx)])
				call pargstr (Memc[extname])
			}
			call sfree (sp)
			call mefclose (mef)
			call error (GEM_ERR,
			    "Index inconsistent with extension name")
		    } else {
			if (debug)
			    call eprintf ("extension name ok\n")
		    }
		}

		if (imx_h_extver (imx)) {
		    version = mef_geti (mef, "EXTVER")
		    if (IMX_EXTVERSION(imx) != version) {
			if (debug) {
			    call eprintf ("inconsistent versions: %d %d\n")
				call pargi (IMX_EXTVERSION(imx))
				call pargi (version)
			}
			call sfree (sp)
			call mefclose (mef)
			call error (GEM_ERR,
			    "Index inconsistent with extension version")
		    } else {
			if (debug)
			    call eprintf ("extension index ok\n")
		    }
		}

	    } else if (imx_h_extname (imx)) {

		index = -2 	# avoids bug in meflib (MEF_CGROUP)
		version = INDEFL
		if (imx_h_extver (imx))
		    version = IMX_EXTVERSION(imx)

		if (debug) {
		    call eprintf ("selecting name/version: '%s', %d\n")
			call pargstr (Memc[IMX_EXTNAME(imx)])
			call pargi (IMX_EXTVERSION(imx))
		}

		iferr (call imx_rdhdr (mef, index, Memc[IMX_EXTNAME(imx)], 
		    version)) {
		    if (debug) {
			code = errget (Memc[errmsg], SZ_FNAME)
			call eprintf ("name/version failed: %d '%s'\n")
			    call pargi (code)
			    call pargstr (Memc[errmsg])
		    }
		    call sfree (sp)
		    call mefclose (mef)
		    call erract (EA_ERROR)
		}

		index = MEF_CGROUP(mef)

		if (debug) {
		    call eprintf ("found index %d\n")
			call pargi (index)
		}
	    }

	    if (debug)
		call eprintf ("imx_mef_open ok\n")

	    call sfree (sp)
	    return mef

	} else {
	    call sfree (sp)
	    call erract (EA_ERROR)
	}
end


# IMX_MEF_FILE -- Open a mef structure from meflib
# This tries adding .fit and .fits extensions if no extension supplied
# Extension information is ignored
# Raises errors on failure

pointer	procedure imx_mef_file (imx, acmode, oldp)

pointer	imx			# I The structure to open
int	acmode			# I Meflib access mode
pointer	oldp			# IO Meflib reserved value

pointer	sp, image, section, mef, errmsg, extname, extver
int	code, index
bool	debug, ok
pointer	mef_open(), errget()
bool	imx_h_extension()

begin
	debug = false
	ok = false

	call smark (sp)
	call salloc (image, SZ_FNAME, TY_CHAR)
	call salloc (errmsg, SZ_LINE, TY_CHAR)

	# Note that we save pointers to strings, so don't use get/set
	# routines when temporarily deleting values below (which would
	# deallocate the values)

	# cannot specify section, so save it
	section = IMX_SECTION(imx); IMX_SECTION(imx) = NULL
	# handle index and extension after opening (see below)
	index = IMX_CLINDEX(imx); call imx_d_clindex (imx)
	extname = IMX_EXTNAME(imx); IMX_EXTNAME(imx) = NULL
	extver = IMX_EXTVERSION(imx); call imx_d_extver (imx)

	call imx_condense (imx, Memc[image], SZ_FNAME)

	if (debug) {
	    call eprintf ("imx_mef_file: '%s'\n")
		call pargstr (Memc[image])
	}

	iferr (mef = mef_open (Memc[image], acmode, oldp)) {

	    if (debug) {
		code = errget (Memc[errmsg], SZ_FNAME)
		call eprintf ("imx_mef_file: mef_open: %s '%s'\n")
		    call pargi (code)
		    call pargstr (Memc[errmsg])
	    }

	    # Try with different fits extensions
	    if (! ok && ! imx_h_extension (imx)) {
		call imx_s_extension (imx, "fit")
		call imx_condense (imx, Memc[image], SZ_FNAME)
		ifnoerr (mef = mef_open (Memc[image], acmode, oldp)) {
		    ok = true
		}
		call imx_d_extension (imx)
	    }
	    if (code == SYS_FOPEN && ! imx_h_extension (imx)) {
		call imx_s_extension (imx, "fits")
		call imx_condense (imx, Memc[image], SZ_FNAME)
		ifnoerr (mef = mef_open (Memc[image], acmode, oldp)) {
		    ok = true
		}
		call imx_d_extension (imx)
	    }

	} else {
	    ok = true
	}

	IMX_SECTION(imx) = section
	call imx_s_clindex (imx, index)
	IMX_EXTNAME(imx) = extname
	call imx_s_extver (imx, extver)

	call sfree (sp)

	if (debug) {
	    call eprintf ("imx_mef_file: %b\n")
		call pargb (ok)
	}

	if (! ok) {
	    call erract (EA_ERROR)
	} else {
	    return mef
	}
end


# IMX_IMIO_OPEN -- Open an im structure from imio
# This tries adding .fit and .fits extensions if no extension supplied
# If an section/mef extension or other kernel information is supplied, 
# it is used (note that imio will not open binary tables)
# Raises errors on failure

pointer	procedure imx_imio_open (imx, acmode, hdr_arg)

pointer	imx			# I The structure to open
int	acmode			# I Meflib access mode
int	hdr_arg			# IO See sys/imio/immap

pointer	sp, image, im, errmsg
int	code
pointer	immap(), errget()
bool	imx_h_extension()
bool	debug

begin
	debug = false

	call smark (sp)
	call salloc (image, SZ_FNAME, TY_CHAR)
	call salloc (errmsg, SZ_LINE, TY_CHAR)

	call imx_condense (imx, Memc[image], SZ_FNAME)

	if (debug) {
	    call eprintf ("imx_imio_open: '%s'\n")
		call pargstr (Memc[image])
	}

	iferr (im = immap (Memc[image], acmode, hdr_arg)) {

	    code = errget (Memc[errmsg], SZ_FNAME)

	    if (debug) {
		call eprintf ("imx_imio_open: immap: %s '%s'\n")
		    call pargi (code)
		    call pargstr (Memc[errmsg])
	    }

	    # Try with different fits extensions
	    if (code == SYS_FOPEN && ! imx_h_extension (imx)) {
		call imx_s_extension (imx, "fit")
		call imx_condense (imx, Memc[image], SZ_FNAME)
		iferr (im = immap (Memc[image], acmode, hdr_arg)) {
		    call imx_d_extension (imx)
		    code = errget (Memc[errmsg], SZ_FNAME)
		} else {
		    call imx_d_extension (imx)
		    goto imio_open_
		}
	    }
	    if (code == SYS_FOPEN && ! imx_h_extension (imx)) {
		call imx_s_extension (imx, "fits")
		call imx_condense (imx, Memc[image], SZ_FNAME)
		iferr (im = immap (Memc[image], acmode, hdr_arg)) {
		    call imx_d_extension (imx)
		    code = errget (Memc[errmsg], SZ_FNAME)
		} else {
		    call imx_d_extension (imx)
		    goto imio_open_
		}
	    }

	    # Return on open failure
	    call sfree (sp)
	    call erract (EA_ERROR)
	}

imio_open_

	if (debug)
	    call eprintf ("imx_imio_open end\n")

	call sfree (sp)
	return im
end


# IMX_DEBUG -- Print the file details to stdout (saves messing round with
# memory allocations in debug code)

procedure imx_debug (imx)

pointer	imx			# I The structure to display

pointer	sp, file

begin
	call smark (sp)
	call salloc (file, SZ_FNAME, TY_CHAR)
	call imx_condense (imx, Memc[file], SZ_FNAME)
	call eprintf ("imx: '%s'\n")
	    call pargstr (Memc[file])
	call sfree (sp)
end


# IMX_DETAIL -- Print the expansion details to stdout (for debug)

procedure imx_detail (imx)

pointer	imx			# I The structure to display

bool	imx_h_path(), imx_h_extension(), imx_h_clindex(), imx_h_clsize()
bool	imx_h_extname(), imx_h_extver(), imx_h_ikparams(), imx_h_section()

begin
	if (NULL != imx) {
	    call eprintf ("imx details...\n")
	    if (imx_h_path (imx)) {
		call eprintf ("imx path: '%s'\n")
		    call pargstr (Memc[IMX_PATH(imx)])
	    } else {
		call eprintf ("no path\n")
	    }
	    call eprintf ("imx file: '%s'\n")
		call pargstr (Memc[IMX_FILE(imx)])
	    if (imx_h_extension (imx)) {
		call eprintf ("imx extension: '%s'\n")
		    call pargstr (Memc[IMX_EXTENSION(imx)])
	    } else {
		call eprintf ("no extension\n")
	    }
	    if (imx_h_clindex (imx)) {
		call eprintf ("imx index: '%d'\n")
		    call pargi (IMX_CLINDEX(imx))
	    } else {
		call eprintf ("no index\n")
	    }
	    if (imx_h_clsize (imx)) {
		call eprintf ("imx size: '%d'\n")
		    call pargi (IMX_CLSIZE(imx))
	    } else {
		call eprintf ("no size\n")
	    }
	    if (imx_h_extname (imx)) {
		call eprintf ("imx extname: '%s'\n")
		    call pargstr (Memc[IMX_EXTNAME(imx)])
	    } else {
		call eprintf ("no extname\n")
	    }
	    if (imx_h_extver (imx)) {
		call eprintf ("imx version: '%d'\n")
		    call pargi (IMX_EXTVERSION(imx))
	    } else {
		call eprintf ("no version\n")
	    }
	    if (imx_h_ikparams (imx)) {
		call eprintf ("imx ikparams: '%s'\n")
		    call pargstr (Memc[IMX_IKPARAMS(imx)])
	    } else {
		call eprintf ("no ikparams\n")
	    }
	    if (imx_h_section (imx)) {
		call eprintf ("imx section: '%s'\n")
		    call pargstr (Memc[IMX_SECTION(imx)])
	    } else {
		call eprintf ("no section\n")
	    }
	} else {
	    call eprintf ("NULL imx pointer in imx_detail\n")
	}
end


# IMX_RDHDR_GN -- Call mef_rdhdr_gn and throw an error on failure or EOF
# (EOF is only an error if you were expecting otherwise, of course)

procedure imx_rdhdr_gn (mef, gn)

pointer	mef			# I Mef descriptior
int	gn			# I Index to read

pointer	sp, msg
int	status, err
bool	debug
int	mef_rdhdr_gn(), errget()

begin
	debug = false
	if (debug) {
	    call eprintf ("calling mef_rdhdr_gn %d\n")
		call pargi (gn)
	}

	iferr (status = mef_rdhdr_gn(mef, gn)) {
	    if (debug) {
		call smark (sp)
		call salloc (msg, SZ_LINE, TY_CHAR)
		err = errget (Memc[msg], SZ_LINE)
		call eprintf ("Error: %d '%s'\n")
		    call pargi (err)
		    call pargstr (Memc[msg])
		call sfree (sp)
	    }
	    call erract(EA_ERROR)
	}
	    
	if (debug) {
	    call eprintf ("mef_rdhdr_gn status %d\n")
		call pargi (status)
	}

	if (ERR == status) {
	    call smark (sp)
	    call salloc (msg, SZ_LINE, TY_CHAR)
	    call sprintf (Memc[msg], SZ_LINE, "Error opening group %d")
		call pargi (gn)
	    iferr (call error (GEM_ERR, Memc[msg])) {
		call sfree (sp)
		call erract (EA_ERROR)
	    }
	}

	if (EOF == status) {
	    call smark (sp)
	    call salloc (msg, SZ_LINE, TY_CHAR)
	    call sprintf (Memc[msg], SZ_LINE, "EOF opening group %d")
		call pargi (gn)
	    iferr (call error (GEM_ERR, Memc[msg])) {
		call sfree (sp)
		call erract (EA_ERROR)
	    }
	}
end

# IMX_RDHDR -- Call mef_rdhdr and throw an error on failure or EOF
# (EOF is only an error if you were expecting otherwise, of course)

procedure imx_rdhdr (mef, group, extname, extver)

pointer	mef			# I Mef descriptior
int	group			# I Index to read
char	extname[ARB]		# I Extension name to read
int	extver			# I Extension version to read

pointer	sp, msg
int	status
bool	debug
int	mef_rdhdr()

errchk	mef_rdhdr

begin
	debug = false
	if (debug) {
	    call eprintf ("calling mef_rdhdr %d %s %d\n")
		call pargi (group)
		call pargstr (extname)
		call pargi (extver)
	}

	status = mef_rdhdr (mef, group, extname, extver)
	if (debug) {
	    call eprintf ("mef_rdhdr status %d\n")
		call pargi (status)
	}

	if (ERR == status) {
	    call smark (sp)
	    call salloc (msg, SZ_LINE, TY_CHAR)
	    call sprintf (Memc[msg], SZ_LINE, "Error opening %d/%s/%d")
		call pargi (group)
		call pargstr (extname)
		call pargi (extver)
	    iferr (call error (GEM_ERR, Memc[msg])) {
		call sfree (sp)
		call erract (EA_ERROR)
	    }
	}

	if (EOF == status) {
	    call smark (sp)
	    call salloc (msg, SZ_LINE, TY_CHAR)
	    call sprintf (Memc[msg], SZ_LINE, "EOF opening %d/%s/%d")
		call pargi (group)
		call pargstr (extname)
		call pargi (extver)
	    iferr (call error (GEM_ERR, Memc[msg])) {
		call sfree (sp)
		call erract (EA_ERROR)
	    }
	}
end


# Copyright(c) 2004-2013 Association of Universities for Research in Astronomy, Inc.
#

include	<imhdr.h>
include	<error.h>
include	<mef.h>
include	"gemerrors.h"
include	"../../lib/imexplode.h"
include	"gemextn.h"


# See gemextensions.x
#
# This file contains:
#
#              gx_verify(...) - verify an image spec against requirements
#
# Support routines:
#         gx_auto_verify(...) - checks against image or extension
# gx_im_verify(imx, imgcheck) - checks against the image alone
#      gx_test_map(imx, mode) - attempt to map file
#           gx_im_absent(imx) - test for img=absent
#           gx_im_exists(imx) - test for img=exists
#             gx_any_mef(imx) - test for img/ext=mef
#     bool gx_verify_mef(imx) - gx_verify mef format
#           gx_any_write(imx) - test for img/ext=write
# gx_ex_verify(imx, extcheck) - checks against extensions
#           gx_ex_absent(imx) - test for ext=absent
#            gx_ex_empty(imx) - test for ext=empty
#           gx_ex_exists(imx) - test for ext=exists
#  gx_ex_xtension(image, imx) - test for ext=image/table
#       gx_read_xtension(...) - read extension with meflib


define	XTENSION	"XTENSION"
define	EXTEND		"EXTEND"
define	IMAGE		"IMAGE"
define	BINTABLE	"BINTABLE"


# GX_VERIFY -- GX_VERIFY access to image, extension, etc
# Throws an exception if a test fails

procedure gx_verify (imx, gxn)

pointer	imx			# I The image specification to gx_verify
pointer	gxn			# I Task parameters

bool	imx_h_any_extension()

errchk	gx_raise_fail_imx, gx_auto_verify, gx_im_verify, gx_ex_verify

begin
	if (1 < GXN_AUTO_CHECK(gxn))
	    call gx_auto_verify (imx, GXN_AUTO_CHECK(gxn))

	else if (1 < GXN_IMG_CHECK(gxn))
	    call gx_im_verify (imx, GXN_IMG_CHECK(gxn))

	else if (1 < GXN_EXT_CHECK(gxn))
	    if (! imx_h_any_extension (imx))
		call gx_raise_fail_imx (imx, 
		    "Extensions check without extension")
	    else
	        call gx_ex_verify (imx, GXN_EXT_CHECK(gxn))
end


# GX_AUTO_VERIFY -- GX_VERIFY image related properties
# Throws an exception if a test fails

procedure gx_auto_verify (imx, autocheck)

pointer	imx			# I File spec to check
int	autocheck[ARB]		# I Checks for the image

int	i
bool	with_extn
bool	imx_h_any_extension()

errchk	gx_raise_fail_imx, gx_ex_absent, gx_im_absent, gx_ex_empty
errchk	gx_ex_exists, gx_im_exists, gx_ex_xtension, gx_any_mef, gx_any_write

begin
	with_extn = imx_h_any_extension (imx)

	iferr {
	    for (i = 2; i <= autocheck[1]; i = i+1) {

		switch (autocheck[i]) {

		case AUTO_ABSENT:
		    if (with_extn)
			call gx_ex_absent (imx)
		    else
			call gx_im_absent (imx)

		case AUTO_EMPTY:
		    if (with_extn) {
			call gx_ex_exists (imx)
			call gx_ex_empty (imx)
		    }

		case AUTO_EXISTS:
		    if (with_extn)
			call gx_ex_exists (imx)
		    else
			call gx_im_exists (imx)

		case AUTO_IMAGE:
		    if (with_extn) {
			call gx_ex_exists (imx)
			call gx_ex_xtension (true, imx)
		    }

		case AUTO_MEF:
		    if (with_extn)
			call gx_ex_exists (imx)
		    else
			call gx_im_exists (imx)
		    call gx_any_mef (imx)

		case AUTO_TABLE:
		    if (with_extn)
			call gx_ex_xtension (false, imx)

		case AUTO_WRITE:
		    if (with_extn)
			call gx_ex_exists (imx)
		    else
			call gx_im_exists (imx)
		    call gx_any_write (imx)

		case AUTO_FORCE:
		    # not used here (see gx_check_depends in t_gemextn.x)

		default:
		    call error (GEM_ERR, "Missing case in auto_verify")
		}
	    }
	} then {
	    call erract (EA_ERROR)
	}
end


# GX_IM_VERIFY -- GX_VERIFY image related properties
# Throws an exception if a test fails

procedure gx_im_verify (imx, imgcheck)

pointer	imx			# I Full file spec for error msg
int	imgcheck[ARB]		# I Checks for the image

int	i
pointer	copy
bool	ok
pointer	imx_copy()
bool	imx_h_any_extension()

errchk	gx_im_absent, gx_im_exists, gx_any_write, gx_any_mef

begin
	copy = NULL
	ok = true

	# remove extension information for testing
	if (imx_h_any_extension (imx)) {
	    copy = imx_copy (imx)
	    call imx_d_any_extension (imx)
	}

	iferr {

	    for (i = 2; i <= imgcheck[1]; i = i+1) {

		switch (imgcheck[i]) {

		case IMG_ABSENT:
		    call gx_im_absent (imx)

		case IMG_EXISTS:
		    call gx_im_exists (imx)

		case IMG_WRITE:
		    call gx_im_exists (imx)
		    call gx_any_write (imx)

		case IMG_MEF:
		    call gx_im_exists (imx)
		    call gx_any_mef (imx)

		default:
		    call error (GEM_ERR, "Missing case in img_verify")
		}
	    }

	} then {ok = false}

	if (NULL != copy) {
	    call imx_free (imx)
	    imx = copy
	}

	if (! ok)
	    call erract (EA_ERROR)
end


# GX_TEST_MAP -- Test for access
# Returns true if access possible

bool procedure gx_test_map (imx, mode)

pointer	imx			# I The file specification
int	mode			# I Access mode

pointer	mef, im
pointer	imx_mef_open(), imx_imio_open()

begin
	iferr {
	    mef = imx_mef_open (imx, mode, 0)
	} then {
	    iferr {
		im = imx_imio_open (imx, mode, 0)
	    } then {
		return false
	    } else {
		call imunmap (im)
		return true
	    }
	} else {
	    call mef_close (mef)
	    return true
	}
end


# GX_IM_ABSENT -- Test for img=absent (no file)
# Complains if access to any kind of file is possible (change help docs)

procedure gx_im_absent (imx)

pointer	imx			# I The file specification

pointer	name, sp, imx2, mef, im
bool	exists
int	access()
pointer	imx_mef_open(), imx_explod(), imx_imio_open()

begin
	call smark (sp)
	call salloc (name, SZ_LINE, TY_CHAR)
	call imx_full_file (imx, Memc[name], SZ_LINE)
	exists = (YES == access (Memc[name], 0, 0))
	if (! exists) {
	    # try opening with implicit .fits etc
	    imx2 = imx_explod (Memc[name])
	    ifnoerr (mef = imx_mef_open (imx2, READ_ONLY, 0)) {
		exists = true
		call mef_close (mef)
	    }
	    ifnoerr (im = imx_imio_open (imx2, READ_ONLY, 0)) {
		exists = true
		call mef_close (im)
	    }
	    call imx_free (imx2)
	}
	call sfree (sp)
	if (exists)
	    call gx_raise_fail_imx (imx, "File exists")
end


# IMG_EXISTS -- Test for img=exists (can be opened)
# Assumes no extension info in imx

procedure gx_im_exists (imx)

pointer	imx			# I The spec to gx_verify

bool	gx_test_map()

begin
	if (! gx_test_map (imx, READ_ONLY)) {
	    IMX_CLINDEX(imx) = 0
	    if (! gx_test_map (imx, READ_ONLY)) {
		call imx_d_clindex (imx)
		call gx_raise_fail_imx (imx, "File not readable or absent")
	    } else {
		call imx_d_clindex (imx)
	    }
	}
end


# GXANY_MEF -- Test for img=mef (EXTEND = true)
# Assumes no extension info in imx

procedure gx_any_mef (imx)

pointer	imx			# I The spec to gx_verify

pointer	copy
pointer	imx_copy()
bool	gx_verify_mef()

begin
	copy = imx_copy (imx)
	call imx_d_any_extension (copy)
	if (! gx_verify_mef (copy)) {
	    IMX_CLINDEX(copy) = 0
	    if (! gx_verify_mef (copy)) {
		call imx_free (copy)
		call gx_raise_fail_imx (imx, "File not MEF")
	    } else {
		call imx_free (copy)
	    }
	}
end


# GX_VERIFY_MEF -- Test for MEF format (local routine)
# Returns true if EXTEND true

bool procedure gx_verify_mef (imx)

pointer	imx			# I The spec to gx_verify

pointer	mef
pointer	imx_mef_open()
bool	ok
bool	mef_getb()

begin
	ok = false
	ifnoerr (mef = imx_mef_open (imx, READ_ONLY, 0)) {
	    if (0 != MEF_HDRP(mef)) {
		iferr (ok = mef_getb (mef, EXTEND)) {
		    ok = false
		}
	    }
	    call mef_close (mef)
	}
	return ok
end


# ANY_WRITE -- Test for img/ext=write (can be opened for writing)

procedure gx_any_write (imx)

pointer	imx			# I The spec to gx_verify

pointer	copy
pointer	imx_copy()
bool	gx_test_map()

begin
	copy = imx_copy (imx)
	call imx_d_any_extension (copy)
	if (! gx_test_map (copy, WRITE_ONLY) &&
	    ! gx_test_map (copy, READ_WRITE)) {
	    IMX_CLINDEX(copy) = 0
	    if (! gx_test_map (copy, WRITE_ONLY) &&
		! gx_test_map (copy, READ_WRITE)) {
		call imx_free (copy)
		call gx_raise_fail_imx (imx, "File not writable")
	    } else {
		call imx_free (copy)
	    }
	}
end


# GX_EX_VERIFY -- GX_VERIFY extension related properties
# Throws an exception if a test fails

procedure gx_ex_verify (imx, extcheck)

pointer	imx			# I The spec to gx_verify
int	extcheck[ARB]		# I Checks for the extensions

int	i

errchk	gx_ex_absent, gx_ex_empty, gx_ex_exists, gx_ex_xtension

begin
	for (i = 2; i <= extcheck[1]; i = i+1) {

	    switch (extcheck[i]) {

	    case EXT_ABSENT:
		call gx_ex_absent (imx)

	    case EXT_EMPTY:
		call gx_ex_exists (imx)
		call gx_ex_empty (imx)

	    case EXT_EXISTS:
		call gx_ex_exists (imx)

	    case EXT_IMAGE:
		call gx_ex_exists (imx)
		call gx_ex_xtension (true, imx)

	    case EXT_TABLE:
		call gx_ex_xtension (false, imx)

	    default:
		call error (GEM_ERR, "Missing case in ext_verify")
	    }
	}
end


# GX_EX_ABSENT -- Test for ext=absent (no extension)

procedure gx_ex_absent (imx)

pointer	imx			# I The spec to gx_verify

bool	gx_test_map()

begin
	if (gx_test_map (imx, READ_ONLY))
	    call gx_raise_fail_imx (imx, "Extension exists")
end


# GX_EX_EMPTY -- Test for ext=empty (no data)
# Assumes image data (should we test tables are empty too?)

procedure gx_ex_empty (imx)

pointer	imx			# I The spec to gx_verify

pointer	im, sp, file
pointer	immap()
bool	ok

begin
	call smark (sp)
	call salloc (file, SZ_LINE, TY_CHAR)
	call imx_condense (imx, Memc[file], SZ_LINE)
	iferr (im = immap (file, READ_ONLY, 0)) {
	    call sfree (sp)
	    call gx_raise_fail_imx (imx, "Extension not readable or absent")
	} else {
	    call sfree (sp)
	    ok = 0 == IM_NDIM(im)
	    call imunmap (im)
	    if (! ok)
		call gx_raise_fail_imx (imx, "Extension contains data")
	}
end


# EXT_EXISTS -- Test for ext=exists (can be opened)

procedure gx_ex_exists (imx)

pointer	imx			# I The spec to gx_verify

bool	gx_test_map()

begin
	if (! gx_test_map (imx, READ_ONLY))
	    call gx_raise_fail_imx (imx, "Extension not readable or absent")
end


# GX_EX_XTENSION -- Test for ext=image or img=table (XTENSION values)

procedure gx_ex_xtension (image, imx)

bool	image			# I True for image, otherwise table
pointer	imx			# I The spec to gx_verify

pointer	sp, xtension
bool	ok, debug
bool	streq()

begin
	debug = false

	call smark (sp)
	ok = false
	call salloc (xtension, SZ_LINE, TY_CHAR)
	call gx_read_xtension (imx, Memc[xtension], SZ_LINE)

	if (debug) {
	    call eprintf ("XTENSION: %s\n")
		call pargstr (Memc[xtension])
	}

	if (image) {
	    ok = streq (Memc[xtension], IMAGE)
	} else {
	    ok = streq (Memc[xtension], BINTABLE)
	}
	call sfree (sp)
	if (! ok) {
	    if (image) {
		call gx_raise_fail_imx (imx, 
		    "File not image (XTENSION != IMAGE)")
	    } else {
		call gx_raise_fail_imx (imx,
		    "File not table (XTENSION != BINTABLE)")
	    }
	}
end


# GX_READ_XTENSION -- Read the value of the XTENSION keyword in the FITS
# header.

procedure gx_read_xtension (imx, value, length)

pointer	imx			# I The file spec to read from
char	value[ARB]		# O The value of XTENSION (or empty string)
int	length			# I Length of value

pointer	mef
pointer	imx_mef_open()

begin
	value[1] = EOS
	ifnoerr (mef = imx_mef_open (imx, READ_ONLY, 0)) {
	    if (0 != MEF_HDRP(mef))
		iferr (call mef_gstr (mef, XTENSION, value, length))
		    value[1] = EOS
	    call mef_close (mef)
	}
end


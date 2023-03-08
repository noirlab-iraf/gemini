# Copyright(c) 2004-2013 Association of Universities for Research in Astronomy, Inc.
#

include	<error.h>
include	"gemextn.h"
include	"../../lib/imexplode.h"


# See gemextensions.x
#
# This file contains:
#
#   gx_check_and_out(...) - check and output the file spec


# GX_CHECK_AND_OUT -- Check and output the file
# May throw exceptions received from gx_verify

procedure gx_check_and_out (imx, gxn)

pointer	imx			# I The image to output
pointer	gxn			# IO Task parameters

pointer	sp, image, imxikpadded, imxreplaced, imgreplaced, params
pointer	imx_explod(), imx_copy()
bool	debug

begin
	debug = false

	if (debug) {
	    call eprintf ("check_and_out\n")
	    call imx_debug (imx)
	}

	call smark (sp)
	call salloc (image, SZ_FNAME, TY_CHAR)
	call salloc (imgreplaced, SZ_FNAME, TY_CHAR)
	call salloc (params, SZ_LINE, TY_CHAR)

	# add image kernel params
	imxikpadded = imx_copy (imx)
	Memc[params] = EOS
	call strcat (Memc[IMX_IKPARAMS[imxikpadded]], Memc[params], SZ_LINE)
	call strcat (GXN_IKPARAMS(gxn), Memc[params], SZ_LINE)
	call imx_s_ikparams (imxikpadded, Memc[params])
	call imx_condense (imxikpadded, Memc[image], SZ_FNAME)
	if (debug) {
	    call eprintf ("name before substitution: %s\n")
	    call pargstr (Memc[image])
	    call imx_detail (imxikpadded)
	}
	call imx_free (imxikpadded)

	# do substitution
	call gt_regex_replace (Memc[image], SZ_FNAME, Memc[imgreplaced],
	    GXN_REPLACE(gxn))
	imxreplaced = imx_explod (Memc[imgreplaced])
	if (debug) {
	    call eprintf ("name after substitution: %s\n")
	    call pargstr (Memc[imgreplaced])
	}
	# Not really needed, but forces consistent format after substitution
	call imx_condense (imxreplaced, Memc[imgreplaced], SZ_FNAME)
	if (debug) {
	    call eprintf ("name after substitution + parsing: %s\n")
	    call pargstr (Memc[imgreplaced])
	    call imx_detail (imxreplaced)
	}

	# drop components as required
	call gx_omit_comps (imxreplaced, GXN_OMIT(gxn))
	if (debug) {
	    call eprintf ("name after dropping components:\n")
	    call imx_detail (imxreplaced)
	}

	# gx_verify
	# if this is replaced with iferr(...) then the program fails
	# to work (unlearn and run with an unknown file).
	iferr {
	    call gx_verify (imxreplaced, gxn)
	} then {
	    call imx_free (imxreplaced)
	    call sfree (sp)
	    call erract (EA_ERROR)
	}
        if (debug)
            call eprintf ("verified\n")

	# if we get this far, gx_verify succeeded and we're done
	call imx_condense (imxreplaced, Memc[imgreplaced], SZ_FNAME)
	if (NULL != GXN_FDOUT(gxn)) {
            if (debug) {
                call eprintf ("output: %s\n")
                call pargstr (Memc[imgreplaced])
            }
	    call fprintf (GXN_FDOUT(gxn), "%s\n")
		call pargstr (Memc[imgreplaced])
            if (debug)
                call eprintf ("flushing\n")
	    call flush (GXN_FDOUT(gxn))
            if (debug)
                call eprintf ("flushed\n")
	}
	GXN_COUNT(gxn) = GXN_COUNT(gxn) + 1

        if (debug)
            call eprintf ("freeing\n")
	call imx_free (imxreplaced)
	call sfree (sp)
        if (debug)
            call eprintf ("checked and outted\n")
end


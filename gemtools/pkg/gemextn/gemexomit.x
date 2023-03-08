# Copyright(c) 2004-2006 Association of Universities for Research in Astronomy, Inc.
#

include	"gemextn.h"


# See gemextensions.x
#
# This file contains:
#
# gx_omit_comps(imx, omit) - omit components of file spec


# GX_OMIT_COMPS -- Omit components of the file

procedure gx_omit_comps (imx, omit)

pointer	imx			# I The image to output
int	omit[ARB]		# I Omissions for output

int	i

begin
	for (i = 2; i <= omit[1]; i = i+1) {
	    switch (omit[i]) {
	    case OMIT_PATH:	call imx_d_path (imx)
	    case OMIT_EXTENSION: call imx_d_extension (imx)
	    case OMIT_INDEX: 	call imx_d_clindex (imx)
	    case OMIT_NAME:	call imx_d_extname (imx)
	    case OMIT_VERSION:	call imx_d_extver (imx)
	    case OMIT_PARAMS:	call imx_d_ikparams (imx)
	    case OMIT_KERNEL:	call imx_d_any_ikparams (imx)
	    case OMIT_SECTION:	call imx_d_section (imx)
	    }
	}
end


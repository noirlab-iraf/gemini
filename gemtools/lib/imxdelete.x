# Copyright(c) 2004-2006 Association of Universities for Research in Astronomy, Inc.
#

include	"imexplode.h"

# See imexplode.x for full documentation

# This file contains:
#
#        imx_d_path(imx) - delete the file path
#        imx_d_file(imx) - delete the file name (without path or extension)
#   imx_d_extension(imx) - delete the file extension
#     imx_d_clindex(imx) - delete the cluster index
#     imx_d_extname(imx) - delete the image extension name
#      imx_d_extver(imx) - delete the image extension version
#    imx_d_ikparams(imx) - delete other image kernel parameters
# imx_d_any_ikparamsimx) - delete the name, version and other kernel parameters
#     imx_d_section(imx) - delete the image section
#
# Note that these routines (should) manage memory allocation
# correctly.  Don't delete components by clearing the structure values
# directly.

# Warning - not sure how to handle cluster size.  Currently persistent
# if index is deleted, but not displayed.

procedure imx_d_path (imx)
pointer	imx
begin
	if (NULL != IMX_PATH(imx)) {
	    call mfree (IMX_PATH(imx), TY_CHAR)
	    IMX_PATH(imx) = NULL
	}
end

procedure imx_d_file (imx)
pointer	imx
begin
	if (NULL != IMX_FILE(imx)) {
	    call mfree (IMX_FILE(imx), TY_CHAR)
	    IMX_FILE(imx) = NULL
	}
end

procedure imx_d_extension (imx)
pointer	imx
begin
	if (NULL != IMX_EXTENSION(imx)) {
	    call mfree (IMX_EXTENSION(imx), TY_CHAR)
	    IMX_EXTENSION(imx) = NULL
	}
end

procedure imx_d_clindex (imx)
pointer	imx
begin
	IMX_CLINDEX(imx) = NO_INDEX
end

procedure imx_d_extname (imx)
pointer	imx
begin
	if (NULL != IMX_EXTNAME(imx)) {
	    call mfree (IMX_EXTNAME(imx), TY_CHAR)
	    IMX_EXTNAME(imx) = NULL
	}
end

procedure imx_d_extver (imx)
pointer	imx
begin
	IMX_EXTVERSION(imx) = NO_INDEX
end

procedure imx_d_ikparams (imx)
pointer	imx
begin
	if (NULL != IMX_IKPARAMS(imx)) {
	    call mfree (IMX_IKPARAMS(imx), TY_CHAR)
	    IMX_IKPARAMS(imx) = NULL
	}
end

procedure imx_d_any_ikparams (imx)
pointer	imx
begin
	call imx_d_extname (imx)
	call imx_d_extver (imx)
	call imx_d_ikparams (imx)
end

procedure imx_d_any_extension (imx)
pointer	imx
begin
	call imx_d_extname (imx)
	call imx_d_extver (imx)
	call imx_d_clindex (imx)
end

procedure imx_d_section (imx)
pointer	imx
begin
	if (NULL != IMX_SECTION(imx)) {
	    call mfree (IMX_SECTION(imx), TY_CHAR)
	    IMX_SECTION(imx) = NULL
	}
end

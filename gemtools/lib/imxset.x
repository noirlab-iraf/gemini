# Copyright(c) 2004-2006 Association of Universities for Research in Astronomy, Inc.
#

include	"imexplode.h"

# See imexplode.x for full documentation

# This file contains:
#
#           imx_s_path(imx, path) - set file path
#           imx_s_file(imx, file) - set file name
# imx_s_extension(imx, extension) - set file extension
#     imx_s_clindex(imx, clindex) - set cluster index
#       imx_s_clsize(imx, clsize) - set cluster size
#     imx_s_extname(imx, extname) - set the image extension name
#       imx_s_extver(imx, extver) - set the image extension version
#   imx_s_ikparams(imx, ikparams) - set other image kernel parameters
#     imx_s_section(imx, section) - set the image section
#
# Note that these routines (should) manage memory allocation
# correctly.


procedure imx_s_path (imx, path)
pointer	imx
char	path[ARB]
begin
	call imx_d_path (imx)
	call imx_copy_string (IMX_PATH(imx), path)
end

procedure imx_s_file (imx, file)
pointer	imx
char	file[ARB]
begin
	call imx_d_file (imx)
	call imx_copy_string (IMX_FILE(imx), file)
end

procedure imx_s_extension (imx, extension)
pointer	imx
char	extension[ARB]
begin
	call imx_d_extension (imx)
	call imx_copy_string (IMX_EXTENSION(imx), extension)
end

procedure imx_s_clindex (imx, clindex)
pointer	imx
int	clindex
begin
	IMX_CLINDEX(imx) = clindex
end

procedure imx_s_clsize (imx, clsize)
pointer	imx
int	clsize
begin
	IMX_CLSIZE(imx) = clsize
end

procedure imx_s_extname (imx, extname)
pointer	imx
char	extname[ARB]
begin
	call imx_d_extname (imx)
	call imx_copy_string (IMX_EXTNAME(imx), extname)
end

procedure imx_s_extver (imx, extver)
pointer	imx
int	extver
begin
	IMX_EXTVERSION(imx) = extver
end

procedure imx_s_ikparams (imx, ikparams)
pointer	imx
char	ikparams[ARB]
begin
	call imx_d_ikparams (imx)
	call imx_copy_string (IMX_IKPARAMS(imx), ikparams)
end

procedure imx_s_section (imx, section)
pointer	imx
char	section[ARB]
begin
	call imx_d_section (imx)
	call imx_copy_string (IMX_SECTION(imx), section)
end

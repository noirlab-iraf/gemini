# Copyright(c) 2004-2006 Association of Universities for Research in Astronomy, Inc.
#

include	"imexplode.h"

# See imexplode.x for full documentation

# This file contains:
#
#          bool imx_h_path(imx) - has file path?
#     bool imx_h_extension(imx) - has file extension?
#       bool imx_h_clindex(imx) - has image cluster index?
#        bool imx_h_clsize(imx) - has image cluster size?
#       bool imx_h_extname(imx) - has image extension name?
#        bool imx_h_extver(imx) - has image extension version?
#      bool imx_h_ikparams(imx) - has any other image kernel parameters?
#  bool imx_h_any_ikparams(imx) - has extension name, version or other param?
# bool imx_h_any_extension(imx) - has index, extension name or version?
#       bool imx_h_section(imx) - has image section?
#
# Tests return true if the associated component(s) is present

bool procedure imx_h_path (imx)
pointer	imx
begin
	return NULL != IMX_PATH(imx) && EOS != Memc[IMX_PATH(imx)]
end

bool procedure imx_h_extension (imx)
pointer	imx
begin
	return NULL != IMX_EXTENSION(imx) && EOS != Memc[IMX_EXTENSION(imx)]
end

bool procedure imx_h_clindex (imx)
pointer	imx
begin
	return IMX_CLINDEX(imx) > NO_INDEX
end

bool procedure imx_h_clsize (imx)
pointer	imx
begin
	return IMX_CLSIZE(imx) > NO_INDEX
end

bool procedure imx_h_extname(imx)
pointer	imx
begin
	return NULL != IMX_EXTNAME(imx) && EOS != Memc[IMX_EXTNAME(imx)]
end

bool procedure imx_h_extver (imx)
pointer	imx
begin
	return IMX_EXTVERSION(imx) > NO_INDEX
end

bool procedure imx_h_ikparams (imx)
pointer	imx
begin
	return NULL != IMX_IKPARAMS(imx) && EOS != Memc[IMX_IKPARAMS(imx)]
end

bool procedure imx_h_any_ikparams (imx)
pointer	imx
bool	imx_h_extname(), imx_h_ikparams(), imx_h_extver()
begin
	return imx_h_extname (imx) || imx_h_extver (imx) ||
	imx_h_ikparams (imx)
end

bool procedure imx_h_any_extension(imx)
pointer	imx
bool	imx_h_extname(), imx_h_clindex(), imx_h_extver()
begin
	return imx_h_extname (imx) || imx_h_extver (imx) || imx_h_clindex (imx)
end

bool procedure imx_h_section(imx)
pointer	imx
begin
	return NULL != IMX_SECTION(imx) && EOS != Memc[IMX_SECTION(imx)]
end


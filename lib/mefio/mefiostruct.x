# Copyright(c) 2004-2005 Association of Universities for Research in Astronomy, Inc.

# ... include statements ...
#
include "mefio.h"

# This file contains: Routines to access the MEFIO structures
#	
#         mesnumexts(...) - get number of extensions, including PHU
#

#
# MESNUMEXTS -- return number of extensions member, ME_NUMEXTS

int procedure mesnumexts (mep)
pointer mep

begin

return ME_NUMEXTS(mep)

end

procedure i2ddims(i2dp, x, y)
pointer i2dp
int x,y
begin
	
	x = I2D_XW(i2dp)
	y = I2D_YW(i2dp)

end

# Copyright(c) 2004-2006 Association of Universities for Research in Astronomy, Inc.
#

include	"scan.h"

# SCANC -- Return the next character from the scanned input.

procedure gt_scanc (scn, cval)

pointer	scn			# IO The scan data
char    cval			# O Th echaracter read

begin
        cval = GT_SC_SCANBUF(scn, GT_SC_IP(scn))
        if (cval != EOS)
            GT_SC_IP(scn) = GT_SC_IP(scn) + 1
end

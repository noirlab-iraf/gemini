# Copyright(c) 2004-2006 Association of Universities for Research in Astronomy, Inc.
#

include "scan.h"

# GT_RESET_SCAN -- Initialize the scan common at the start of a scan.
# May also be called by the user to rescan a line, following a
# conversion failure.

procedure gt_reset_scan (scn)

pointer	scn			# IO The scan data

begin
        GT_SC_IP(scn) = 1
        GT_SC_NTOKENS(scn) = 0
        GT_SC_STOPSCAN(scn) = false
end

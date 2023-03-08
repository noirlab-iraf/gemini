# Copyright(c) 2004-2006 Association of Universities for Research in Astronomy, Inc.
#

include "scan.h"

# GT_FREE_SCAN -- Delete memory allocated for scan

procedure gt_free_scan (scn)

pointer	scn			# IO The scan data

begin
	if (NULL != scn) {
	    call mfree (GT_P_SC_SCANBUF(scn), TY_CHAR)
	    call mfree (scn, TY_STRUCT)
	    scn = NULL
	}
end

# Copyright(c) 2004-2006 Association of Universities for Research in Astronomy, Inc.
#

include "scan.h"

# GT_SSCAN -- Begin a scan from a string.  Only the first newline terminated
# line in the string buffer will be scanned.

# Returns a structure that must be passed to related scanner routines

pointer procedure gt_sscan (str)

char    str[ARB]		# I The string to scan

int     len, ip, op
pointer	scn
int	strlen()

begin
	len = strlen (str)
	call malloc (scn, LEN_SSCAN, TY_STRUCT)
	call malloc (GT_P_SC_SCANBUF(scn), len, TY_CHAR)

        op = 1
        for (ip=1;  str[ip] != EOS && str[ip] != '\n';  ip=ip+1) {
            GT_SC_SCANBUF(scn, op) = str[ip]
            op = op + 1
        }

        GT_SC_SCANBUF(scn, op) = EOS
        call gt_reset_scan (scn)	# initialize scan

	return scn
end

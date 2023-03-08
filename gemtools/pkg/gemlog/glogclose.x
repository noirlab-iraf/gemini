# Copyright(c) 2004-2006 Association of Universities for Research in Astronomy, Inc.

# Author: Kathleen Labrie 25-June-2004

include "glog.h"

.help
.nf
GLOGCLOSE - Close a Gemini log file block
 	status = glogclose ( gl, fl_success )
	
	status		: Exit status code  [return value, (int)]
	gl		: GL Structure  [input,(GL)]
	fl_success	: Was task successful? (YES|NO)  [input, (int)]
.fi
.endhelp

#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
#--------------------------------------------------------------------------

# GLOGCLOSE - Close a Gemini log file block
#  	status = glogclose ( gl, fl_success )
#
#	status		: Exit status code  [return value, (int)]
# 	gl		: GL Structure  [input,(GL)]
# 	fl_success	: Was task successful? (YES|NO)  [input, (int)]

int procedure glogclose (gl, fl_success)

pointer	gl			#I GL Structure 
int	fl_success		#I Was task successful? (YES|NO)

int	status

# Other variable
int	tmpstatus
pointer	sp, tmpstr

# Gemini functions
int	glw_vis(), glw_stat(), glw_title(), glw_tag()

# IRAF functions
int	errget()

errchk	glw_vis(), glw_stat()

begin
	status = 0

	# Allocate stack memory
	call smark (sp)
	call salloc (tmpstr, SZ_LINE, TY_CHAR)

	# Write formatted log entries ...
	iferr {
	    # Write the exit status of the task
	    status = status + glw_vis (gl, G_SHRT_DASH)
	    status = status + glw_stat (gl, fl_success , STAT_LEVEL)

	    # Write closing log comment
	    status = status + glw_title (gl, STAT_LEVEL, END_TAG)
	} then {
	    status = errget (Memc[tmpstr], SZ_LINE)
	    # postpone error call. will try to tag and close cleanly first.
	} else {
	    if (status != 0)
		call sprintf(Memc[tmpstr], SZ_LINE,
		    "1 Error closing the logfile.")
	    # postpone error call. will try to tag and close cleanly first.
	}

	# Tag the end of the parent task's log entries
	iferr ( tmpstatus = glw_tag (gl, END_TAG) ) {
	    # previous error (status != 0) takes precedence
	    if (status == 0)
	        status = errget (Memc[tmpstr], SZ_LINE)
	    # postpone error call.  will try to close cleanly first
	}

	# Close the log file.
	# No matter what happen in gl_close(), gl will be freed up.
	iferr ( call gl_close(gl) ) {
	    # previous error (status != 0) takes precedence
	    if (status == 0)
	        status = errget (Memc[tmpstr], SZ_LINE)
	    gl = NULL
	    call sfree (sp)
	    call error (status, Memc[tmpstr])
	} else {
	    if (status != 0) {		# from previous error that was postponed
	        gl = NULL
	        call sfree (sp)
		call error (status, Memc[tmpstr])
	    }
	}

	# Free stack memory
	call sfree (sp)

	return (status)
end

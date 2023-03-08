# Copyright(c) 2004-2006 Association of Universities for Research in Astronomy, Inc.

# Author: Kathleen Labrie 29-June-2004

include "glog.h"
include "gemerrors.h"
include <time.h>

.help
.nf
GLOGPRINT - Print to a Gemini log file  (API for t_glprint)
	status = glogprint( gl, loglevel, type, str, op )
	
	status		: Exit status  [return value, (int)]
	gl		: GL Structure  [input, (GL)]
	loglevel	: Level of the log entry [input, (int)]
	type		: Type of log entry [input, (int)]
	str		: String attached to the log entry [input, (string)]
	op		: OP (API options) Structure [input, (OP)]
.fi
.endhelp

#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
#--------------------------------------------------------------------------

# GLOGPRINT - Print to a Gemini log file  (API for t_glprint)
# 	status = glogprint( gl, loglevel, type, str, op )
# 	
# 	status		: Exit status  [return value, (int)]
# 	gl		: GL Structure  [input, (GL)]
# 	loglevel	: Level of the log entry [input, (int)]
# 	type		: Type of log entry [input, (int)]
# 	str		: String attached to the log entry [input, (string)]
# 	op		: OP (API options) Structure [input, (OP)]

int procedure glogprint( gl, loglevel, type, str, op )

pointer	gl		#I GL Structure
int	loglevel	#I Level of the log entry (eg. STAT_LEVEL, VIS_LEVEL...)
int	type		#I Type of log entry (eg. G_STR_LOG, G_ERR_LOG, etc)
char	str[ARB]	#I String to be written or name of file
pointer	op		#I/O OP structure

int	status		#O Exit status

# Other variable
char	errmsg[SZ_LINE]
int	tmpverbose
long	ltime
pointer	sp, timestr, timedate

# Gemini functions
int	glw_vis(), glw_err(), glw_file(), glw_fork(), glw_str(), glw_warn()

# IRAF functions
int	errget(), strcmp()
long	clktime()

begin
	status = 0

	# Decide whether this entry should be written at all, before formatting.
	# (glw_str() does this check but there is no point formatting
	#  a string we know won't be written)
	
	switch (loglevel) {
	case ENG_LEVEL:
	    if (GL_REQENG(gl) == NO)
		return (status)
	case SCI_LEVEL:
	    if (GL_REQSCI(gl) == NO)
		return (status)
	case STAT_LEVEL:
	    if (GL_REQSTAT(gl) == NO)
		return (status)
	case VIS_LEVEL:
	    if (GL_REQVIS(gl) == NO)
		return (status)
	case TSK_LEVEL:
	    if (GL_REQTSK(gl) == NO)
		return (status)
	default:
	    status = G_INTERNAL_ERROR
	    call sprintf (errmsg, SZ_LINE, "Unrecognized log level.")
	    call error (status, errmsg)
	}

	# Okay - If we made it to this point, we know that the line should be 
	# written to the logfile; we are not working for nothing.

	# Allocate stack memory
	call smark (sp)
	call salloc (timestr, SZ_LINE, TY_CHAR)
	call salloc (timedate, SZ_LINE, TY_CHAR)

	# Get the time and convert to long string.  Then format string.
	ltime = clktime(0)
	call cnvtime (ltime, Memc[timestr], SZ_TIME)
	call sprintf( Memc[timedate], SZ_LINE, "at [%s]" )
	    call pargstr( Memc[timestr] )

	# Call the appropriate glwrite routine
	iferr {
	    switch (type) {
	    case G_ERR_LOG:
	        # Store user's verbose preference ERROR must be written to 
		# screen, unless the logfile is dev$null, in which case the 
		# user's verbose setting is kept.
		if ( strcmp ("dev$null", GL_LOGFILE(gl)) != 0 ) {
		    tmpverbose = GL_VERBOSE(gl)
		    GL_VERBOSE(gl) = YES
		}
		status = status + glw_vis (gl, G_LONG_DASH)
		status = status + 
		    glw_err (gl, -(OP_ERRNO(op)), Memc[timedate], loglevel)
		status = status + glw_err (gl, OP_ERRNO(op), str, loglevel)
		status = status + glw_vis (gl, G_LONG_DASH)
		if ( strcmp ("dev$null", GL_LOGFILE(gl)) != 0 )
		    GL_VERBOSE(gl) = tmpverbose	    
	    case G_FILE_LOG:
        	status = glw_file (gl, str, loglevel)
	    case G_FORK_LOG:
        	status = glw_fork (gl, OP_CHILD(op), OP_FORK(op), loglevel)
	    case G_STR_LOG:
        	status = glw_str (gl, str, loglevel)
	    case G_VIS_LOG:
        	status = glw_vis (gl, OP_VISTYPE(op))
	    case G_WARN_LOG:
        	status = status + 
		    glw_warn (gl, -(OP_ERRNO(op)), Memc[timedate], loglevel)
		status = status + glw_warn (gl, OP_ERRNO(op), str, loglevel)
	    default:
        	status = G_INTERNAL_ERROR
		call sprintf (errmsg, SZ_LINE, "Unrecognized log type.")
		call error (status, errmsg)
	    }
	} then {
	    status = errget (errmsg, SZ_LINE)
	    OP_STATUS(op) = status
	    call sfree (sp)
	    call error (status, errmsg)
	}

	# Free memory
	call sfree (sp)

	return (status)
end

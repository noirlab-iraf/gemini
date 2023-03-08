# Copyright(c) 2004-2006 Association of Universities for Research in Astronomy, Inc.

# Author: Kathleen Labrie 22-June-2004

include "glog.h"
include "gemerrors.h"
include <time.h>

.help
.nf
GLOGINIT - Set up a Gemini log file  (API for t_glinit)
	gl = gloginit( logfile, curtask, curpack, paramstr, pp, op )
	
	gl		: GL Structure  [return value, (GL)]
	logfile		: Name of the logfile  [input, (string)]
	curtask		: Name of the current task  [input, (string)]
	curpack		: Name of the current task's package  [input, (string)]
	paramstr	: Formatted list of task's param/value pairs [ ", (") ]
	pp		: glogpars pset pointer  [input, (pointer)]
	op		: OP (API options) Structure [input, (OP)]
.fi
.endhelp

#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
#--------------------------------------------------------------------------

# GLOGINIT - Set up a Gemini log file  (API for t_glinit)
# 	gl = gloginit( logfile, curtask, curpack, paramstr, pp, op )
# 	
# 	gl		: GL Structure  [return value, (GL)]
# 	logfile		: Name of the logfile  [input, (string)]
# 	curtask		: Name of the current task  [input, (string)]
# 	curpack		: Name of the current task's package  [input, (string)]
# 	paramstr	: Formatted list of task's param/value pairs [ ", (") ]
#	pp		: glogpars pset pointer  [input, (pointer)]
# 	op		: OP (API options) Structure [input, (OP)]

pointer procedure gloginit ( logfile, curtask, curpack, paramstr, pp, op )

char	logfile[ARB]		#I Name of the logfile
char	curtask[ARB]		#I Name of the current task
char	curpack[ARB]		#I Name of the current task's package
char	paramstr[ARB]		#I Formatted list of task's param/value pairs
pointer	pp			#I glogpars pset pointer
pointer	op			#I/O OP structure

pointer gl			#O GL structure, return value

int	status
long	ltime
pointer sp, tmpstr, tskdate, timestr

# Gemini functions
int	glw_param(), glw_tag(), glw_title(), glw_vis(), glw_warn()
bool	g_whitespace()
pointer	glogopen()

# IRAF functions
int	errget()
long	clktime()

errchk	glw_tag(), glw_title(), glw_vis(), glw_warn(), glw_param()

begin
	status = 0

	# Allocate stack memory
	call smark (sp)
	call salloc (tmpstr, SZ_FNAME, TY_CHAR)

	# Open the logfile and assign the GL structure  (retrieve glogpars)
	iferr ( gl = glogopen (logfile, curtask, curpack, pp, op) ) {
	    gl = NULL
	    OP_STATUS(op) = errget (Memc[tmpstr], SZ_LINE)
	    call sfree (sp)
	    call error (OP_STATUS(op), Memc[tmpstr])
	}

	# Tag the beginning of the parent task's log entries
	iferr {
	    status = status + glw_tag (gl, BEGIN_TAG)
	    status = status + glw_title (gl, STAT_LEVEL, BEGIN_TAG)
	    status = status + glw_vis (gl, G_EMPTY)
	} then {
	    OP_STATUS(op) = errget (Memc[tmpstr], SZ_LINE)
	    call gl_close (gl)
	    call sfree (sp)
	    call error (OP_STATUS(op), Memc[tmpstr])
	} else {
	    if (status != 0) {
 		OP_STATUS(op) = 1
 		call sprintf (Memc[tmpstr],SZ_LINE,
		    "1 Error initializing the logfile.")
 		call gl_close (gl)
		call sfree (sp)
 		call error (OP_STATUS(op), Memc[tmpstr])
	    }
	}

	# If using default log name, write warning.
	# Whether or not a default is used has been determined in glogopen()
	
	if (OP_DEFLOG(op) == YES) {

	    # Time string
	    call malloc (timestr, SZ_TIME, TY_CHAR)
	    call malloc (tskdate, SZ_LINE, TY_CHAR)
	    ltime = clktime(0)
	    call cnvtime (ltime, Memc[timestr], SZ_TIME)
	    call sprintf (Memc[tskdate], SZ_LINE, "at [%s]")
 		call pargstr (Memc[timestr])

	    # Message string
	    call sprintf( Memc[tmpstr], SZ_LINE,
 		"\"%s.logfile\" and \"%s.logfile\" empty.\nUsing default name: \"%s\"")
 		call pargstr (GL_CURTASK(gl))
 		call pargstr (GL_CURPACK(gl))
 		call pargstr (GL_LOGFILE(gl))

	    # Write warning to log
	    iferr {
 		status = status +
 		    glw_warn (gl, -(G_USING_DEFAULT), Memc[tskdate], STAT_LEVEL)
 		status = status +
 		    glw_warn (gl, G_USING_DEFAULT, Memc[tmpstr], STAT_LEVEL)
 		status = status + glw_vis (gl, G_EMPTY)
	    } then {
 		OP_STATUS(op) = errget (Memc[tmpstr], SZ_LINE)
 		call gl_close (gl)
 		call mfree (timestr, TY_CHAR)
 		call mfree (tskdate, TY_CHAR)
 		call sfree (sp)
 		call error (OP_STATUS(op), Memc[tmpstr])
	    } else {
 		if (status != 0) {
 		    OP_STATUS(op) = 2
 		    call sprintf (Memc[tmpstr], SZ_LINE,
 			"2 Error initializing the logfile.")
 		    call gl_close (gl)
 		    call mfree (timestr, TY_CHAR)
 		    call mfree (tskdate, TY_CHAR)
 		    call sfree (sp)
		    call error (OP_STATUS(op), Memc[tmpstr])
 		}
	    }

	    # we do not need those anymore
	    call mfree (timestr, TY_CHAR)
	    call mfree (tskdate, TY_CHAR)
	}

	# Write the parameters string, if non-empty
	if ( ! g_whitespace (paramstr) ) {
	    iferr {
        	status = status + glw_vis (gl, G_SHRT_DASH)
 		status = status + glw_param (gl, paramstr)
 		status = status + glw_vis (gl, G_SHRT_DASH)
	    } then {
 		OP_STATUS(op) = errget (Memc[tmpstr], SZ_LINE)
 		call gl_close (gl)
 		call sfree (sp)
 		call error (OP_STATUS(op), Memc[tmpstr])
	    } else {
 		if (status != 0) {
 		    OP_STATUS(op) = 3
 		    call sprintf(Memc[tmpstr],SZ_LINE,
		        "3 Error initializing the logfile.")
 		    call gl_close (gl)
 		    call sfree (sp)
 		    call error (OP_STATUS(op), Memc[tmpstr])
 		}
	    }
	}

	# Free stack memory
	call sfree (sp)

	return (gl)
end

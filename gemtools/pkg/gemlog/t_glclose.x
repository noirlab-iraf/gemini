# Copyright(c) 2004-2006 Association of Universities for Research in Astronomy, Inc.

# Author: Kathleen Labrie  3-May-2004

# GLCLOSE -- Close a Gemini log file

include "glog.h"

procedure t_glclose ()

#char	logfile[SZ_FNAME]		# Name of the log file
#char	curtask[SZ_FNAME]		# Name of the parent task
#bool	fl_success			# Was task successful?
#bool	verbose				# Verbose?
#int	status				# Exit status (0=good)

# Local variables for task parameters
char	l_logfile[SZ_FNAME], l_curtask[SZ_FNAME]
int	l_fl_success, l_verbose
int	l_status

# Other variables
pointer	gl, pp, op, sp, tmpstr

# Gemini functions
int 	glogclose()
pointer	glogopen()

# IRAF functions
pointer	clopset()
bool	clgetb()
int	btoi(), errget()

begin
	l_status = 0

	# Get task parameter values
	call clgstr ("logfile", l_logfile, SZ_FNAME)
	call clgstr ("curtask", l_curtask, SZ_FNAME)
	l_fl_success = btoi (clgetb ("fl_success"))
	l_verbose = btoi (clgetb ("verbose"))

	# Get GLOGPARS pset pointer  (closed in glogopen() )
	pp = clopset ("glogpars")

	# Allocate stack memory
	call smark (sp)
	call salloc (tmpstr, SZ_LINE, TY_CHAR)

	# Allocate structure memory ('gl' is allocated in glogopen() )
	call opalloc (op)

	# Set required log options, then open file
	OP_FL_APPEND(op) = YES
	OP_FORCE_APPEND(op) = YES	 #Log file *must* already exist.
	OP_VERBOSE(op) = l_verbose
	gl = NULL
	call strcpy ("", Memc[tmpstr], SZ_LINE)

	iferr (gl = glogopen (l_logfile, l_curtask, Memc[tmpstr], pp, op)) {
	   gl = NULL
	   l_status = errget (Memc[tmpstr], SZ_LINE)
	   call printf ("GLOGCLOSE ERROR: %d %s.\n")
	       call pargi (l_status)
	       call pargstr (Memc[tmpstr])
	   call clputi ("status", l_status)
	   if (pp != NULL)
	       call clcpset (pp)
	   call opfree (op)
	   call sfree (sp)
	   return
	}

	# Write closing comments, and close logfile
	iferr ( l_status = glogclose (gl, l_fl_success) ) {
	    gl = NULL		# glogclose should have freed gl already.
	    l_status = errget (Memc[tmpstr], SZ_LINE)
	    call printf ("GLOGCLOSE ERROR: %d %s.\n")
		call pargi (l_status)
		call pargstr (Memc[tmpstr])
	    call clputi ("status", l_status)
	}

	# Free memory
	call opfree (op)
	call sfree (sp)

	call clputi ("status", l_status)
	return
end

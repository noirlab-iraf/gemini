# Copyright(c) 2004-2006 Association of Universities for Research in Astronomy, Inc.

# Author: Kathleen Labrie  26-August-2004

# GLFIX -- Fix a Gemini log file by adding missing EOE tags.

include "glog.h"

procedure t_glfix ()

#char	logfile[SZ_FNAME]		# Name of the log file
#bool	verbose				# Verbose?
#int	status				# Exit status (0=good)

# Local variables for task parameters
char	l_logfile[SZ_FNAME]
bool	l_verbose
int	l_status

# Other variables
char	curtask[SZ_FNAME]
pointer	gl, op, pp, sp, tmpstr

# Gemini functions
int	glogfix()
pointer	glogopen()

# IRAF functions
int	errget(), btoi()
bool	clgetb()
pointer	clopset()

begin
	l_status = 0
	call strcpy ("glogfix", curtask, SZ_FNAME)
	
	# Get task parameter values
	call clgstr ("logfile", l_logfile, SZ_FNAME)
	l_verbose = clgetb ("verbose")
	
	# Get GLOGPARS pset pointer   (closed in glogopen() )
	pp = clopset ("glogpars")

	# Allocate stack memory
	call smark (sp)
	call salloc (tmpstr, SZ_LINE, TY_CHAR)
	
	# Allocate structure memory
	call opalloc (op)
	
	# Set required log options, then open file
	OP_FL_APPEND(op) = YES
	OP_FORCE_APPEND(op) = YES	# Log file *must* already exist.
	OP_VERBOSE(op) = btoi (l_verbose)
	gl = NULL
	call strcpy ("", Memc[tmpstr], SZ_LINE)
	
	iferr (gl = glogopen (l_logfile, curtask, Memc[tmpstr], pp, op)) {
	    gl = NULL
	    l_status = errget (Memc[tmpstr], SZ_LINE)
	    call printf ("GLOGFIX ERROR: %d %s.\n")
	        call pargi (l_status)
		call pargstr (Memc[tmpstr])
	    call clputi ("status", l_status)
	    if (pp != NULL)
	        call clcpset (pp)
	    call opfree (op)
	    call sfree (sp)
	    return
	}
	
	# Clear unused memory
	call opfree (op)

	# Fix log
	iferr ( l_status = glogfix (gl) ) {
	    l_status = errget (Memc[tmpstr], SZ_LINE)
	    call printf ("GLOGFIX ERROR: %d %s.\n")
	        call pargi (l_status)
		call pargstr (Memc[tmpstr])
	    call clputi ("status", l_status)
	}
	
	# Close logfile  (nothing fancy, just use gl_close())
	iferr ( call gl_close (gl) ) {
	    gl = NULL
	    l_status = errget (Memc[tmpstr], SZ_LINE)
	    call printf ("GLOGFIX ERROR: %d %s.\n")
	        call pargi (l_status)
		call pargstr (Memc[tmpstr])
	}
	
	# Free remaining memory
	call sfree (sp)
	
	# Exiting...
	call clputi ("status", l_status)
	return

end

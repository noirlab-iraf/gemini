# Copyright(c) 2004-2006 Association of Universities for Research in Astronomy, Inc.

# Author: Kathleen Labrie  3-May-2004

# GLINIT -- Set up a Gemini log file 

include "glog.h"

procedure t_glinit ()

#char	logfile[SZ_FNAME]		# Name of the log file
#char	curtask[SZ_FNAME]		# Name of the current task
#char	curpack[SZ_FNAME]		# Name of the task's package
#char	paramstr[G_SZ_PARAMSTR]	# User formatted list of the task's
					# param/value pairs
#bool	fl_append			# Append to logfile?
#bool	verbose				# Verbose?
#int	status				# Exit status (0=good)

# Local variables for task parameters
char	l_logfile[SZ_FNAME], l_curtask[SZ_FNAME], l_curpack[SZ_FNAME]
char	l_paramstr[G_SZ_PARAMSTR]
bool	l_fl_append, l_verbose
int	l_status

# Other variables
int	sz_paramstr, ind, indnx, len
pointer	newpstr, workpstr, tmpstr
pointer	gl, pp, sp, op

# Gemini functions
pointer	gloginit()
bool	g_whitespace()

# IRAF functions
pointer clopset()
bool	clgetb()
int	btoi(), strlen(), strsearch(), errget()

begin
	l_status = 0

	# Get task parameter values
	call clgstr ("gloginit.logfile", l_logfile, SZ_FNAME)
	call clgstr ("curtask", l_curtask, SZ_FNAME)
	call clgstr ("curpack", l_curpack, SZ_FNAME)
	call clgstr ("paramstr", l_paramstr, G_SZ_PARAMSTR)
	l_fl_append = clgetb ("fl_append")
	l_verbose = clgetb ("verbose")

	# Get GLOGPARS pset pointer  (closed in gloginit) 
	pp = clopset ("glogpars")

	# Allocate stack memory
	call smark (sp)
	call salloc (tmpstr, SZ_FNAME, TY_CHAR)
	sz_paramstr = strlen (l_paramstr) + 1
	call salloc (newpstr, sz_paramstr, TY_CHAR)
	call salloc (workpstr, sz_paramstr, TY_CHAR)

	# Allocate structure memory  ('gl' is allocated in gloginit() )
	call opalloc (op)

	# Fix the CL version of paramstr.
	# This means replacing the "\n" by true newline characters, '\n'
	#  (CL does not play nice with '\n')

	if ( ! g_whitespace( l_paramstr ) ) {
	    ind = 1
	    indnx = -1
	    call strcpy ("", Memc[newpstr], sz_paramstr)
	    while ( (indnx != 0) && ( ind < sz_paramstr) ) {
 		call strcpy (l_paramstr[ind], Memc[workpstr], sz_paramstr)
        	indnx = strsearch (Memc[workpstr], "\\n")
        	if ( indnx != 0 ) {
		    len = indnx - 3
 		    call strcpy (l_paramstr[ind], Memc[workpstr], len)
 		    call strcat (Memc[workpstr], Memc[newpstr], sz_paramstr)
 		    call strcat ("\n", Memc[newpstr],  sz_paramstr)
 		} else {
 		    len = strlen (l_paramstr) - ind +1
 		    call strcpy (l_paramstr[ind], Memc[workpstr], len)
 		    call strcat (Memc[workpstr], Memc[newpstr], sz_paramstr)
 		}
 		ind = ind + indnx - 1
	    }
	    len = strlen (Memc[newpstr])
	    call strcpy (Memc[newpstr], l_paramstr, len)
	}

	# Set log options required to initialize logfile, then initialize.
	OP_FL_APPEND(op) = btoi (l_fl_append)
	OP_FORCE_APPEND(op) = NO
	OP_VERBOSE(op) = btoi (l_verbose)
	gl = NULL
	iferr (gl = gloginit(l_logfile, l_curtask, l_curpack,l_paramstr,pp,op)){
	    gl = NULL
	    l_status = errget (Memc[tmpstr], SZ_LINE)
	    call printf ("GLOGINIT ERROR: %d %s.\n")
		call pargi (l_status)
		call pargstr (Memc[tmpstr])
	    call clputi ("status", l_status)
	    if (pp != NULL)
		call clcpset (pp)
	    call opfree (op)
	    call sfree (sp)
	    return
	}

	# Store the logfile name used in the 'logfile' parameter of this task
	# This way, if the logfile is set in here, the calling task can 
	# retrieve it.
	
	call clpstr ("logfile",GL_LOGFILE(gl))

	# Close the log file
	iferr ( call gl_close(gl) ) {
	    l_status = errget (Memc[tmpstr], SZ_LINE)
	    call printf ("GLOGINIT ERROR: %d %s.\n")
		call pargi (l_status)
		call pargstr (Memc[tmpstr])
	}

	# Free memory
	call opfree (op)
	call sfree (sp)

	call clputi ("status", l_status)
	return
end

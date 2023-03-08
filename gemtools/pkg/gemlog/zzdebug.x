# Copyright(c) 2004-2006 Association of Universities for Research in Astronomy, Inc.

# Test/example routines for the GEMLOG interface and functions

task	demoapi		= t_demoapi

include "glog.h"


# DEMOAPI -- Demo of the GEMLOG API
#
# The parameter file of DEMOAPI is:
#
# logfile,s,a,,,,"Logfile name"
# glogpars,pset,h,"",,,"Logging preferences"
# verbose,b,h,yes,,,"Verbose"
# status,i,h,0,,,"Exit status (0=good)"
#
# Note: All GEMLOG-enabled task must have the GLOGPARS pset as one of the
# parameters
#
# To compile and run DEMOAPI:
#   cl> gemini
#   ge> gemtools
#   ge> cd gemtools$pkg/gemlog
#   ge> mkpkg -p gemini zzdebug
#   ge> cd
#   ge> task demoapi="gemtools$pkg/gemlog/zzdebug.e"
#   ge> demoapi

procedure t_demoapi ()

# Local variables for task parameters
char	l_logfile[SZ_FNAME]
int	l_status
bool	l_verbose
pointer glpset

# Others
char	paramstr[G_SZ_PARAMSTR]
int	fd, success
pointer	gl, op, sp, tmpstr, scratch

# Gemini function and GEMLOG application interface
bool	g_whitespace()
pointer gloginit()
int	glogclose(), glogprint()

# IRAF functions
int	errget(), open(), btoi(), access()
bool	clgetb()
pointer clopset()

errchk	glogprint()

begin
	l_status = 0
	success = G_SUCCESS		# G_SUCCESS is defined in 'glog.h'

	#Get task parameter values
	#   Here only the value of 'logfile' needs to be retrieved.
	
	call clgstr ("logfile", l_logfile, SZ_FNAME)
	l_verbose = clgetb ("verbose")


	#Get GLOGPARS pset pointer
	#   The pset values are actually queried for in gloginit().
	#   The pset will be closed by gloginit().  No further manipulation of
	#   the pset pointer is required in the current routine.
	
	glpset = clopset ("glogpars")


	# Allocate memory
	#    The memory for the GL structure is allocated by gloginit()
	#    The memory for the OP structure must be allocated here since 
	#    we will need it to send options to gloginit().
	
	call smark (sp)
	call salloc (tmpstr, SZ_LINE, TY_CHAR)
	call salloc (scratch, SZ_FNAME, TY_CHAR)
	call opalloc (op)

        
	# Build the parameter/value pairs string
	#   The string 'paramstr' will be passed to gloginit().  gloginit()
	#   will parse the lines and write the parameter information to the log
	#   file.
	#   First make sure you start with an empty 'paramstr'.  This is not
	#   critical but it will avoid bad surprises.
	#   Then, append to 'paramstr' the name of the parameter and its value.
	
	call strcpy ("", paramstr, G_SZ_PARAMSTR)
	call glogpstring (paramstr, "logfile", l_logfile)
	call glogpb (paramstr, "verbose", l_verbose)
	call glogpi (paramstr, "status", l_status)
     
     
	# Set log options required to initialize logfile, then initialize.
	#   Setting the OP_FL_APPEND(op) to yes means that the log file will be
	#   open in 'append' mode, or it will be created if it does not already
	#   exists.
	#   The gloginit() routine takes care of allocating memory for the GL 
	#   structure.  It also initialize it with the information found in 
	#   GLOGPARS and in the OP structure.  The log file is open, or 
	#   created, and the file descriptor is saved in the GL structure.
	#   Essential log entries are written to the log file, as well as the
	#   list of parameters and their values, 'paramstr'.
	
	OP_FL_APPEND(op) = YES
	OP_VERBOSE(op) = btoi (l_verbose)
	gl = NULL
	iferr (gl = gloginit (l_logfile, "demoapi", "gemtools", paramstr, glpset, op)){
	    gl = NULL		# memory already been freed in gloginit
	    l_status = errget (Memc[tmpstr], SZ_LINE)
	    call printf ("DEMOAPI ERROR: %d %s.\n")
		call pargi (l_status)
		call pargstr (Memc[tmpstr])
	    call opfree (op)
	    call sfree (sp)
	    call clputi ("status", l_status)
	    return
	} else {		# This block is optional
	    if (OP_STATUS(op) != 0) {
 		l_status = l_status + errget (Memc[tmpstr], SZ_LINE)
 		call printf ("DEMOAPI WARNING: %d %s.\n")
 		    call pargi (OP_STATUS(op))
 		    call pargstr (Memc[tmpstr])
 		OP_STATUS(op) = 0
	    }
	}
	call flush (GL_FD(gl))
    
    
	# Print stuff to the logfile
	#   Each type of entries defined in glogprint() is demonstrated.  The
	#   steps directly relevant to the log entry are found between the
	#   ' #--- ' delimiters.
	
	iferr {
	
	    # A 'visual' aid entry
	    l_status = glogprint (gl, VIS_LEVEL, G_STR_LOG, 
        	"Below is a long line of dashes.", op)
	    #---
	    OP_VISTYPE(op) = G_LONG_DASH
	    l_status = glogprint (gl, VIS_LEVEL, G_VIS_LOG, "", op)
	    #---



	    # A 'string' entry
	    #---
	    l_status = glogprint (gl, ENG_LEVEL, G_STR_LOG, 
		"This is a normal string.",op)
	    #---



	    # A 'fork' entries
	    l_status = glogprint (gl, STAT_LEVEL, G_STR_LOG, 
        	"Below are a fork-to and a fork-back-to messages.", op)
	    #---
	    call strcpy ("greattask", OP_CHILD(op), SZ_FNAME)
	    OP_FORK(op) = G_FORWARD
	    l_status = glogprint (gl, STAT_LEVEL, G_FORK_LOG, "", op)
	    OP_FORK(op) = G_BACKWARD
	    l_status = glogprint (gl, STAT_LEVEL, G_FORK_LOG, "", op)
	    #---



	    # A 'warning' entry
	    l_status = glogprint (gl, STAT_LEVEL, G_STR_LOG, 
		"Below is a warning.", op)
	    #---
	    OP_ERRNO(op) = G_INTERNAL_ERROR
	    l_status = glogprint (gl, STAT_LEVEL, G_WARN_LOG, 
		"Watch your back!", op)
	    #---
	    
	    

	    # An 'error' entry
	    l_status = glogprint (gl, STAT_LEVEL, G_STR_LOG, 
		"Below is an error.", op)
	    #---
	    OP_ERRNO(op) = G_OP_UNRECOGNIZED
	    l_status = glogprint (gl, STAT_LEVEL, G_ERR_LOG,"Just kidding!", op)
	    #---
	    
	    
	    
	    # A 'file' entry
	    #    This mode is used to write the content of a file to 'logfile'.
	    #    Here we write the content of the file named 'scratch'.
	    
	    call mktemp ( "scratch", Memc[scratch], SZ_FNAME )
	    fd = open (Memc[scratch], NEW_FILE, TEXT_FILE)
	    call fprintf (fd, "Twinkle, twinkle, little star,\n")
	    call fprintf (fd, "How I wonder what you are.\n")
	    call close (fd)
	    l_status = glogprint (gl, SCI_LEVEL, G_STR_LOG, 
		"Below we add a file.", op)
	    #---
	    l_status = glogprint (gl, SCI_LEVEL, G_FILE_LOG, Memc[scratch], op)
	    #---
	    call delete (Memc[scratch])


	} then {
	    l_status = errget (Memc[tmpstr], SZ_LINE)
	    call printf ("DEMOAPI ERROR: %d %s\n")
		call pargi (l_status)	
		call pargstr (Memc[tmpstr])
	    call flush (GL_FD(gl))
	    if ( access(Memc[scratch], 0, 0) == YES )
	        call delete (Memc[scratch])
	}



	# Close the logfile
	#     The memory for gl is freed in glogclose (whether glogclose is
	#     successful or not).
	if (l_status != 0)
	    success = G_FAILURE
	iferr ( l_status = glogclose (gl, success) ) {
	    l_status = errget (Memc[tmpstr], SZ_LINE)
	    call printf ("DEMOAPI ERROR: %d %s\n")
		call pargi (l_status)
		call pargstr (Memc[tmpstr])
	}
    
	# Free memory
	call opfree (op)
	call sfree (sp)

	call clputi ("status", l_status)
	return
end

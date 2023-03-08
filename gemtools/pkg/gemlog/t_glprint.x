# Copyright(c) 2004-2006 Association of Universities for Research in Astronomy, Inc.

# Author: Kathleen Labrie 4-May,2004

# T_GLPRINT -- Print to a Gemini log file

include "glog.h"

procedure t_glprint ()

#char	logfile[SZ_FNAME]		# Name of the log file
#char	curtask[SZ_FNAME]		# Name of the parent task
#char	loglevel[G_SZ_SSTR]		# Level of the log entry
#char	type[G_SZ_SSTR]			# Type of log entry
#char	str[SZ_LINE]			# String attached to the log entry
#char	fork[G_SZ_SSTR]			# Fork to or back from child process?
#char	child[SZ_FNAME]			# Name of child process
#char	vistype[G_SZ_SSTR]		# Type of visual enhancement
#int	errno				# Gemini error code
#bool	verbose				# Verbose?
#int	status				# Exit status (0=good)

char	l_logfile[SZ_FNAME], l_curtask[SZ_FNAME], l_loglevel[G_SZ_SSTR]
char	l_type[G_SZ_SSTR], l_str[SZ_LINE]
char	l_fork[G_SZ_SSTR], l_child[SZ_FNAME]
char	l_vistype[G_SZ_SSTR]
bool	l_verbose
int	l_errno
int	l_status
pointer	pp

int	level, typecode, junk
pointer	gl, op, tmpstr
pointer	sp, typedict, forkdict, leveldict, visdict

# GEMLOG functions
pointer	glogopen()
int	glogprint()

# IRAF functions
pointer	clopset()
bool	clgetb()
int	btoi(), clgeti(), strcmp(), errget(), clgwrd()

begin
	l_status = 0

	# Allocate stack memory
	call smark(sp)
	call salloc (tmpstr, SZ_LINE, TY_CHAR)

	# Allocate static memory
	call malloc (leveldict, SZ_LINE, TY_CHAR)
	call malloc (typedict, SZ_LINE, TY_CHAR)
	call malloc (forkdict, SZ_LINE, TY_CHAR)
	call malloc (visdict, SZ_LINE, TY_CHAR)

	# Allocate structure memory ('gl' is allocated in glogopen() )
	call opalloc (op)

	# Get task parameter values
	call clgstr ("logfile", l_logfile, SZ_FNAME)
	call clgstr ("curtask", l_curtask, SZ_FNAME)
	call clgstr ("loglevel.p_min", Memc[leveldict], SZ_LINE)
	junk = clgwrd ("loglevel", l_loglevel, G_SZ_SSTR, Memc[leveldict])
	call clgstr ("type.p_min", Memc[typedict], SZ_LINE)
	junk = clgwrd ("type",l_type, G_SZ_SSTR, Memc[typedict])
	call clgstr ("str", l_str, SZ_LINE)
	call clgstr ("fork.p_min", Memc[forkdict], SZ_LINE)
	junk = clgwrd ("fork", l_fork, G_SZ_SSTR, Memc[forkdict])
	call clgstr ("child", l_child, SZ_FNAME)
	call clgstr ("vistype.p_min", Memc[visdict], SZ_LINE)
	junk = clgwrd ("vistype", l_vistype, G_SZ_SSTR, Memc[visdict])
	l_errno = clgeti ("errno")
	l_verbose = clgetb ("verbose")

	# Free static memory (dictionary memory)
	call mfree (leveldict, TY_CHAR)
	call mfree (typedict, TY_CHAR)
	call mfree (forkdict, TY_CHAR)
	call mfree (visdict, TY_CHAR)

	# Get GLOGPARS pset pointer   (closed in glogopen() )
	pp = clopset ("glogpars")

	# Set log options required for opening log
	OP_FL_APPEND(op) = YES
	OP_FORCE_APPEND(op) = YES  #Log file *must* already exist
	OP_VERBOSE(op) = btoi(l_verbose)
	gl = NULL
	call strcpy ("", Memc[tmpstr], SZ_LINE)

	# Open logfile and set GL structure
	iferr (gl = glogopen (l_logfile, l_curtask, Memc[tmpstr], pp, op)) {
	    gl = NULL
	    l_status = errget (Memc[tmpstr], SZ_LINE )
	    call printf ("GLOGPRINT ERROR: %d %s.\n")
		call pargi (l_status)
		call pargstr (Memc[tmpstr])
	    call clputi ("status", l_status)
	    if (pp != NULL)
		call clcpset (pp)
	    call opfree (op)
	    call sfree (sp)
	    return
	}

	# Decipher loglevel
	if ( strcmp ("engineering", l_loglevel) == 0 )
	    level = ENG_LEVEL
	else if ( strcmp ("science", l_loglevel) == 0 )
	    level = SCI_LEVEL
	else if ( strcmp ("status", l_loglevel) == 0 )
	    level = STAT_LEVEL
	else if ( strcmp ("visual", l_loglevel) == 0 )
	    level = VIS_LEVEL
	else if ( strcmp ("task", l_loglevel) == 0 )
	    level = TSK_LEVEL
	else if ( strcmp ("none", l_loglevel) == 0 )
	    level = NO_LEVEL

	# Decipher type, and assign options to OP structure
	if ( strcmp ("error", l_type) == 0 ) {
	    typecode = G_ERR_LOG
	    OP_ERRNO(op) = l_errno

	} else if ( strcmp ("file", l_type) == 0 ) {
	    typecode = G_FILE_LOG

	} else if ( strcmp ("fork", l_type) == 0 ) {
	    typecode = G_FORK_LOG
	    call strcpy (l_child, OP_CHILD(op), SZ_FNAME)
	    if ( strcmp ("forward", l_fork) == 0 )
		OP_FORK(op) = G_FORWARD
	    else
		OP_FORK(op) = G_BACKWARD

	} else if ( strcmp ("string", l_type) == 0 ) {
	    typecode = G_STR_LOG

	} else if ( strcmp ("visual", l_type) == 0 ) {
	    typecode = G_VIS_LOG
	    if ( strcmp ("empty", l_vistype) == 0 )
        	OP_VISTYPE(op) = G_EMPTY
	    else if ( strcmp ("longdash", l_vistype) == 0 )
        	OP_VISTYPE(op) = G_LONG_DASH
	    else if ( strcmp ("shortdash", l_vistype) == 0 )
        	OP_VISTYPE(op) = G_SHRT_DASH

	} else if ( strcmp ("warning", l_type) == 0 ) {
	    typecode = G_WARN_LOG
	    OP_ERRNO(op) = l_errno
	}    

	# Write to log
	iferr ( l_status = glogprint (gl, level, typecode, l_str, op) ) {
	    l_status = errget (Memc[tmpstr], SZ_LINE)
	    call printf ("GLOGPRINT ERROR: %d %s.\n")
		call pargi (l_status)
		call pargstr (Memc[tmpstr])
	}


	# Close the log file
	iferr ( call gl_close(gl) ) {
	    l_status = errget (Memc[tmpstr], SZ_LINE)
	    call printf ("GLOGPRINT ERROR: %d %s.\n")
		call pargi (l_status)
		call pargstr (Memc[tmpstr])
	}

	# Free memory
	call opfree (op)
	call sfree (sp)

	call clputi ("status", l_status)
	return
end

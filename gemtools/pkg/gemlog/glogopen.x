# Copyright(c) 2004-2006 Association of Universities for Research in Astronomy, Inc.

# Author: Kathleen Labrie 12-July-2004

include "glog.h"

.help
.nf
GLOGOPEN - Open a Gemini log file.
	gl = glogopen( logfile, curtask, curpack, pp, op )
	
	gl		: GL Structure  [return value, (GL)]
	logfile		: Name of the logfile  [input, (string)]
	curtask		: Name of the current task  [input, (string)]
	curpack		: Name of the current task's package  [input, (string)]
	pp		: glogpars pset pointer  [input, (pointer)]
	op		: OP (API options) Structure  [input, (OP)]

GLOGOPEN() does not initialize a log block.  It simply opens the log file.  
It is used in particular to append to an already initialized log block.  Use 
GLOGINIT() instead if the log block has not been initialized.
.fi
.endhelp

#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
#--------------------------------------------------------------------------

# GLOGOPEN - Open a Gemini log file.
# 	gl = glogopen( logfile, curtask, curpack, pp, op )
# 	
# 	gl		: GL Structure  [return value, (GL)]
# 	logfile		: Name of the logfile  [input, (string)]
# 	curtask		: Name of the current task  [input, (string)]
# 	curpack		: Name of the current task's package  [input, (string)]
# 	pp		: glogpars pset pointer  [input, (pointer)]
# 	op		: OP (API options) Structure  [input/output, (OP)]

pointer procedure glogopen (logfile, curtask, curpack, pp, op)

char	logfile[ARB]		#I Name of the logfile
char	curtask[ARB]		#I Name of the current task
char	curpack[ARB]		#I Name of the current task's package
pointer	pp			#I glogpars pset pointer
pointer	op			#I/O OP structure

pointer gl			#O GL structure, return value

int	status
int	l_fl_stat_level, l_fl_tsk_level, l_fl_sci_level
int	l_fl_eng_level, l_fl_vis_level

int	fileaccess
pointer	sp, tmpstr

# Gemini functions
bool	g_whitespace()
pointer	gl_open()

# IRAF functions
int	btoi(), errget(), access()
bool	clgpsetb()

begin
	# Initialize
	status = 0
	OP_DEFLOG(op) = NO
	OP_STATUS(op) = status

	# Input arguments ...

	# Convert task and package names to lower case
	call strlwr (curtask)
	call strlwr (curpack)

	# Look into glogpars pset if 'logfile' empty
	if ( g_whitespace (logfile) )
	    call clgpset (pp, "logfile", logfile, SZ_FNAME)

	# Get user's preferences from glogpars pset ('pp').  Then close pset.
	l_fl_stat_level = btoi (clgpsetb (pp, "fl_stat_level"))
	l_fl_tsk_level  = btoi (clgpsetb (pp, "fl_tsk_level"))
	l_fl_sci_level  = btoi (clgpsetb (pp, "fl_sci_level"))
	l_fl_eng_level  = btoi (clgpsetb (pp, "fl_eng_level"))
	l_fl_vis_level  = btoi (clgpsetb (pp, "fl_vis_level"))
	call clcpset (pp)
	pp = NULL

	# ... done with input arguments

	# Allocate stack memory
	call smark (sp)
	call salloc (tmpstr, SZ_FNAME, TY_CHAR)

	# Assign log file name if left empty
	if ( g_whitespace (logfile) ) {
	    call sprintf (Memc[tmpstr], SZ_FNAME, "%s.logfile")
		call pargstr (curpack)
	    call clgstr (Memc[tmpstr], logfile, SZ_FNAME)

	    if ( g_whitespace(logfile) ) {

	        # Still empty. Using default name. Will write warning to log.
 		call sprintf (logfile, SZ_FNAME, "%s.log")
 		    call pargstr (curpack)
 		OP_DEFLOG(op) = YES

	    }
	}

	# Open/create the log file with correct mode.  
	# Allocate and initialize the GL structure, 'gl'.
	
	gl = NULL
	if (OP_FL_APPEND(op) == YES) {
	    if (OP_FORCE_APPEND(op) == YES)
		fileaccess = APPEND
	    else if ( access (logfile, 0, 0) == NO )
		fileaccess = NEW_FILE
	    else
		fileaccess = APPEND
	} else
	    fileaccess = NEW_FILE

	iferr ( gl = gl_open (logfile, fileaccess, status) ) {
	    gl = NULL
	    OP_STATUS(op) = errget (Memc[tmpstr], SZ_LINE)
	    call sfree (sp)
	    call error (OP_STATUS(op), Memc[tmpstr])
	}

	# Fill up the remaining of the GL structure, 'gl'.
	GL_REQSTAT(gl) = l_fl_stat_level
	GL_REQSCI(gl)  = l_fl_sci_level
	GL_REQENG(gl)  = l_fl_eng_level
	GL_REQVIS(gl)  = l_fl_vis_level
	GL_REQTSK(gl)  = l_fl_tsk_level
	GL_VERBOSE(gl) = OP_VERBOSE(op)
	call strcpy (curtask, GL_CURTASK(gl), SZ_FNAME)
	call strcpy (curpack, GL_CURPACK(gl), SZ_FNAME)

	# Free memory
	call sfree (sp)

	return (gl)
end

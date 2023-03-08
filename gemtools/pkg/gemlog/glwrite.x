# Copyright(c) 2004-2006 Association of Universities for Research in Astronomy, Inc.

# Author: Kathleen Labrie   3-May-2004

include "glog.h"
include "gemerrors.h"
include <time.h>

define	G_SZ_DATE	10	#Use in: glw_tag()
define	G_SZ_TIME 	8	#Use in: glw_tag()

.help
.nf
This file contains procedures to format and write a log entry to a Gemini log 
file.

      glw_str - write a string to log, if level requested
      glw_err - write an error message
     glw_file - write the content of a file
     glw_fork - write process info when forking to/from another task
    glw_param - write task parameter/value pairs
     glw_stat - write the exit status of a task
      glw_tag - write begin/end tag
    glw_title - write log opening/closing information entry
      glw_vis - write visual/log-readability improvement entries
     glw_warn - write a warning message

GLW_STR -- Write a string to logfile if the log entry level if that level has 
           been requested by the user.
	status = glw_str( gl, str, level )
	
	status		: Status code		 [return value, (int)]
	gl		: GL structure		 [input, (GL)]
	str		: String to write	 [input, (string)]
	level		: Level of the log entry [input, (int)]

GLW_ERR -- Write an error message to the Gemini log file.
	status = gwl_err( gl, errno, errmsg, level )
	
	status		: Status code		[return value, (int)]
	gl		: GL structure		[input, (GL)]
	errno		: Error code		[input, (int)]
	usrmsg		: User's error message	[input, (char[])]

    If "errno" is not NULL, GLW_ERR will use it the prepend to "usrmsg" the 
    Gemini error message corresponding to "errno".

GLW_FILE -- Write the content of a file to the Gemini log file.  This is used
            to format and write outputs from IRAF tasks.
	status = glw_file( gl, fname, level )
	
	status		: Status code		   [return value, (int)]
	gl		: GL structure		   [input, (GL)]
	fname		: File to read and transfer to logs [input, (char[])]
	level		: Level of the log entries [input, (int)]

    The "level" is applied to the entire content of the file.

GLW_FORK -- Write process info when forking to/from another task.
	status = glw_fork( gl, child, direction, level )
	
	status		: Status code			[return value, (int)]
	gl		: GL structure (parent is GL_CURTASK)	[input, (GL)]
	child		: Name of the child process		[input, (GL)]
	direction	: Forking to or back from child process	[input, (int)]
	level		: Level of the log entry		[input, (int)]

    Direction definitions: Forking to -> G_FORWARD, forking back -> G_BACKWARD
    Suggested log level: STAT_LEVEL.

GLW_PARAM -- Write a task's parameter/value pairs to the Gemini log file.
	status = glw_param( gl, kvsbuf )
	
	status		: Status code		[return value, (int)]
	gl		: GL Structure		[input, (GL)]
	kvsbuf		: String buffer for key/value entries [input, (char[])]

    To write multiple key/value entries at once, each key/value entry must be 
    delimited by newline characters, '\n'.

GLW_STAT -- Write the exit status of a task.
	status = glw_stat( gl, statno, level )
	
	status		: Status code		[return value, (int)]
	gl		: GL structure		[input, (GL)]
	statno		: State code [G_SUCCESS|G_FAILURE  [input, (int)]
	level		: Level of the entry	[input, (int)]

    Suggested log level: STAT_LEVEL

GLW_TAG -- Write the begin/end tag used by GLOGEXTRACT.
	status = glw_tag( gl, tagmode )
	
	status		: Status code		[return value, (int)]
	gl		: GL structure		[input, (GL)]
	tagmode		: Type of tag (BEGIN_TAG|END_TAG) [input, (int)]

GLW_TITLE -- Write log opening/closing information entry
	status = glw_title( gl, level, mode )
	
	status		: Status code		[return value, (int)]
	gl		: GL structure		[input, (GL)]
	level		: Level of the entry	[input, (int)]
	mode		: Type of entry (BEGIN_TAG|END_TAG) [input, (int)]

    Suggested log level: STAT_LEVEL

GLW_VIS -- Write visual/log-readability improvement entries
	status = glw_vis( gl, type )
	
	status		: Status code		[return value, (int)]
	gl		: GL structure		[input, (GL)]
	type		: Type of string	[input, (int)]

    Type of strings are: G_EMPTY, G_LONG_DASH, G_SHRT_DASH, for an empty line,
    a long series of '-' characters, and a short series (20) of '-' characters,
    respectively.

GLW_WARN -- Write a warning message to the Gemini Log file.
        status = glw_warn( gl, errno, errmsg, level )
        
        status          : Status code   [return value, (int)]
        gl              : GL Structure  [input, (GL)]
        errno           : Error code    [input, (int)]
        usrmsg          : User's error message [input, (char[])]
        level           : Level of the entry   [input, (int)]

    If "errno" is not NULL, GLW_ERR will use it to prepend the official error
    message to "usrmsg".

.fi
.endhelp

#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
#--------------------------------------------------------------------------


#GLW_STR -- Write a string to logfile if the log entry level if that level has 
#           been requested by the user.
#	status = glw_str( gl, str, level )
#	
#	status		: Status code		 [return value, (int)]
#	gl		: GL structure		 [input, (GL)]
#	str		: String to write	 [input, (string)]
#	level		: Level of the log entry [input, (int)]

#### MallocDebug sometimes detects a memory leak in one or more of the pargstr()

int procedure glw_str (gl, str, level)

pointer gl              #I GL structure
char    str[ARB]        #I String to write
int     level           #I Level of the log entry

int	status

# Other variables
char	msg[SZ_LINE]
pointer	sp, curtask

# Note to developer: If you change the strings below, change them also
#                    in glchk.x
string	NO_LEVEL_STR	"    "
string	STAT_LEVEL_STR  "STAT"
string	SCI_LEVEL_STR   "SCI "
string	ENG_LEVEL_STR   "ENG "
string	VIS_LEVEL_STR   "VIS "
string	TSK_LEVEL_STR   "TSK "

begin
	status = 0

	#Allocate stack memory
	call smark (sp)
	call salloc (curtask, SZ_FNAME, TY_CHAR)

	#Get uppercase version of GL_CURTASK(gl)
	call strcpy (GL_CURTASK(gl), Memc[curtask], SZ_FNAME)
	call strupr (Memc[curtask])

	#Switch to appropriate log level
	switch (level) {
	case IGNORE_LEVEL:
	    if ( (GL_FD(gl) == STDOUT) || (GL_FD(gl) == STDERR) ) {
		if (GL_VERBOSE(gl) == YES)
		    call fprintf (GL_FD(gl), "%s\n")
			call pargstr (str)
	    } else {
		if (GL_VERBOSE(gl) == YES)
        	    call printf ("%s\n")
        		call pargstr ( str )
		call fprintf (GL_FD(gl), "%s\n")
        	    call pargstr ( str )
	    }

	case NO_LEVEL:
	    if ( (GL_FD(gl) == STDOUT) || (GL_FD(gl) == STDERR) ) {
		if (GL_VERBOSE(gl) == YES)
		    call fprintf (GL_FD(gl), "%s %s\n")
        		call pargstr (NO_LEVEL_STR)
			call pargstr (str)
	    } else {
		if (GL_VERBOSE(gl) == YES)
        	    call printf ("%s %s\n")
        		call pargstr (NO_LEVEL_STR)
        		call pargstr (str)
		call fprintf (GL_FD(gl), "%s %s\n")
        	    call pargstr (NO_LEVEL_STR)
        	    call pargstr (str)
	    }

	case STAT_LEVEL:
	    if (GL_REQSTAT(gl) == YES) {
	        if ( (GL_FD(gl) == STDOUT) || (GL_FD(gl) == STDERR) ) {
		    if (GL_VERBOSE(gl) == YES)
		        call fprintf (GL_FD(gl), "%s %s\n")
		            call pargstr (Memc[curtask])
			    call pargstr (str)
		} else {
        	    if (GL_VERBOSE(gl) == YES)
        		call printf ("%s %s\n")
                	    call pargstr (Memc[curtask])
                	    call pargstr (str)
        	    call fprintf (GL_FD(gl), "%s %s %s\n")
        		call pargstr (STAT_LEVEL_STR)
        		call pargstr (Memc[curtask])
        		call pargstr (str)
		}
	    }
	case SCI_LEVEL:
	    if (GL_REQSCI(gl) == YES) {
	        if ( (GL_FD(gl) == STDOUT) || (GL_FD(gl) == STDERR) ) {
		    if (GL_VERBOSE(gl) == YES)
		        call fprintf (GL_FD(gl), "%s %s\n")
		            call pargstr (Memc[curtask])
			    call pargstr (str)
		} else {
        	    if (GL_VERBOSE(gl) == YES)
        		call printf ("%s %s\n")
                	    call pargstr (Memc[curtask])
                	    call pargstr (str)
        	    call fprintf (GL_FD(gl), "%s %s %s\n")
        		call pargstr (SCI_LEVEL_STR)
        		call pargstr (Memc[curtask])
        		call pargstr (str)
		}
	    }
	case ENG_LEVEL:
	    if (GL_REQENG(gl) == YES) {
	        if ( (GL_FD(gl) == STDOUT) || (GL_FD(gl) == STDERR) ) {
		    if (GL_VERBOSE(gl) == YES)
		        call fprintf (GL_FD(gl), "%s %s\n")
		            call pargstr (Memc[curtask])
			    call pargstr (str)
		} else {
        	    if (GL_VERBOSE(gl) == YES)
        		call printf ("%s %s\n")
                	    call pargstr (Memc[curtask])
                	    call pargstr (str)
        	    call fprintf (GL_FD(gl), "%s %s %s\n")
        		call pargstr (ENG_LEVEL_STR)
        		call pargstr (Memc[curtask])
        		call pargstr (str)
		}
	    }
	case VIS_LEVEL:
	    if (GL_REQVIS(gl) == YES) {
	        if ( (GL_FD(gl) == STDOUT) || (GL_FD(gl) == STDERR) ) {
		    if (GL_VERBOSE(gl) == YES)
		        call fprintf (GL_FD(gl), "%s %s\n")
		            call pargstr (Memc[curtask])
			    call pargstr (str)
		} else {
        	    if (GL_VERBOSE(gl) == YES)
        		call printf ("%s %s\n")
                	    call pargstr (Memc[curtask])
                	    call pargstr (str)
        	    call fprintf (GL_FD(gl), "%s %s %s\n")
        		call pargstr (VIS_LEVEL_STR)
        		call pargstr (Memc[curtask])
        		call pargstr (str)
		}
	    }
	case TSK_LEVEL:
	    if (GL_REQTSK(gl) == YES) {
	        if ( (GL_FD(gl) == STDOUT) || (GL_FD(gl) == STDERR) ) {
		    if (GL_VERBOSE(gl) == YES)
		        call fprintf (GL_FD(gl), "%s %s\n")
		            call pargstr (Memc[curtask])
			    call pargstr (str)
		} else {
         	    if (GL_VERBOSE(gl) == YES)
        		call printf ("%s %s\n")
                	    call pargstr (Memc[curtask])
                	    call pargstr (str)
        	    call fprintf (GL_FD(gl), "%s %s %s\n")
        		call pargstr (TSK_LEVEL_STR)
        		call pargstr (Memc[curtask])
        		call pargstr (str)
		}
	    }
	default:
	    status = G_INTERNAL_ERROR
	    call sprintf (msg, SZ_LINE, "Unrecognized entry level in glw_str()")
	    call sfree (sp)
	    call error (status, msg)
	}

	#Free stack memory
	call sfree (sp)

	return (status)
end

#--------------------------------------------------------------------------

# GLW_ERR -- Write an error message to the Gemini log file.
# 	status = gwl_err( gl, errno, errmsg, level )
# 	
# 	status		: Status code		[return value, (int)]
# 	gl		: GL structure		[input, (GL)]
# 	errno		: Error code		[input, (int)]
# 	usrmsg		: User's error message	[input, (char[])]
# 
#     If "errno" is not NULL, GLW_ERR will use it the prepend to "usrmsg" the 
#     Gemini error message corresponding to "errno".

int procedure glw_err (gl, errno, usrmsg, level)

pointer gl              #I GL structure
int     errno           #I Error code (optional - set it to NULL if not used)
char    usrmsg[ARB]     #I Error message (optional)
int     level           #I Level of the log entry

int     status

# Other variables
char	token, msg[SZ_LINE], curtask[SZ_FNAME]
int     i, nstr, strptr[G_MAX_LINES]
int	tmpstatus
pointer errmsg, tmpstr, sp

# Gemini functions
int	glw_str(), g_splitstr()
bool	g_whitespace()

# IRAF functions
int	errget()

errchk	glw_str()

begin
	status = 0
	tmpstatus = 0

	# Allocate stack memory
	call smark (sp)
	call salloc (errmsg, SZ_LINE, TY_CHAR)
	call salloc (tmpstr, SZ_LINE, TY_CHAR)

	# Initialize
	token = '\n'    

	# Fetch and write the Gemini error message
	if ((errno >= 99) && (errno < 500)) {
	    ifnoerr ( call gemerrmsg (errno, Memc[tmpstr]) ) {
        	call sprintf( Memc[errmsg], SZ_LINE, "ERROR: %d %s")
        	    call pargi (errno)
        	    call pargstr (Memc[tmpstr])
        	tmpstatus = glw_str (gl, Memc[errmsg], level)
		status = status + tmpstatus
	    } else {
		tmpstatus = errget (msg, SZ_LINE)
		if (tmpstatus == G_INPUT_ERROR) {
		    call strcpy (GL_CURTASK(gl), curtask, SZ_FNAME)
		    call strupr (curtask)
		    call printf ("%s WARNING: %d %s\n")
			call pargstr (curtask)
			call pargi (tmpstatus)
			call pargstr (msg)
		} else {
		    status = errget (msg, SZ_LINE)
		    call sfree (sp)
		    call error (status, msg)
		}
	    }
	} else if (errno == NULL) {
	    errno = 1
	} else if (errno < 0) {     # Use errno but ignore Gemini error msg
	    errno = -errno
	}

	# Write the user's error message.
	if ( ! g_whitespace (usrmsg) ) {
	
	    # Split strings on the newline character, '\n'
	    nstr = g_splitstr (usrmsg, token, strptr)

	    for (i = 1; i <= nstr; i = i+1) {
        	call sprintf (Memc[errmsg], SZ_LINE, "ERROR: %d %s")
        	    call pargi (errno)
        	    call pargstr (usrmsg[strptr[i]])
        	tmpstatus = glw_str (gl, Memc[errmsg], level)
		status = status + tmpstatus
	    }
	}

	# Free memory
	call sfree (sp)

	return (status)
end

#--------------------------------------------------------------------------

# GLW_FILE -- Write the content of a file to the Gemini log file.  This is used
#             to format and write outputs from IRAF tasks.
# 	status = glw_file( gl, fname, level )
# 	
# 	status		: Status code		   [return value, (int)]
# 	gl		: GL structure		   [input, (GL)]
# 	fname		: File to read and transfer to logs [input, (char[])]
# 	level		: Level of the log entries [input, (int)]
# 
#     The "level" is applied to the entire content of the file.

int procedure glw_file (gl, fname, level)

pointer gl              #I GL structure
char    fname[ARB]      #I File to read and transfer to logfile
int     level           #I Level of the log entry

int     status

# Other variables
char	msg[SZ_LINE]
int	fd, len, tmpstatus
pointer	sp, linebuf

# Gemini functions
int	glw_str()

# IRAF functions
int	open(), getline(), strlen(), errget()

errchk	glw_str(), open(), getline()

begin
	status = 0
	tmpstatus = 0

	# Open file to read  (there is an errchk on open() )
	# Note to self: what if nothing was sent to the file?  In the case of
	# a redirection, the file would not exist.  Is it really cause for
	# an error?  I think glogprint and t_glprint should take care of that.
	# Anyway, this is something to keep in mind.
	
	fd = open (fname, READ_ONLY, TEXT_FILE)

	# Allocate stack memory
	call smark (sp)
	call salloc (linebuf, SZ_LINE, TY_CHAR)

	# Read file, then write to Gemini log file (one line at a time)
	iferr {
	    while ( getline(fd, Memc[linebuf]) != EOF ) {

        	# Remove newline character
		len = strlen (Memc[linebuf])
		Memc[linebuf+(len-1)] = EOS

		# Write to Gemini log file
		tmpstatus = glw_str (gl, Memc[linebuf], level)
		status = status + tmpstatus
	     }
	} then {
	    status = errget (msg, SZ_LINE)
	    call close (fd)
	    call sfree (sp)
	    call error (status, msg)
	}

	# Close file
	call close (fd)

	# Free stack memory
	call sfree(sp)

	return(status)
end

#--------------------------------------------------------------------------

# GLW_FORK -- Write process info when forking to/from another task.
# 	status = glw_fork( gl, child, direction, level )
# 	
# 	status		: Status code			[return value, (int)]
# 	gl		: GL structure (parent is GL_CURTASK)	[input, (GL)]
# 	child		: Name of the child process		[input, (GL)]
# 	direction	: Forking to or back from child process	[input, (int)]
# 	level		: Level of the log entry		[input, (int)]
# 
#     Direction definitions: Forking to -> G_FORWARD, forking back -> G_BACKWARD
#     Suggested log level: STAT_LEVEL.

int procedure glw_fork (gl, child, direction, level)

pointer gl              #I GL Structure
char    child[ARB]      #I Name of the child process
int     direction       #I Forking to or back from child process
                        #  [G_FORWARD|G_BACKWARD]
int     level           #I Level of the log entry  (suggested level: STAT)

int	status

# Other variables
char	msg[SZ_LINE], curtask[SZ_FNAME]
pointer	sp, str

# Gemini functions
int	glw_str()

errchk	glw_str()

begin
	status = 0

	# Allocate stack memory
	call smark (sp)
	call salloc (str, SZ_LINE, TY_CHAR)

	# Convert to upper case
	call strcpy (GL_CURTASK(gl), curtask, SZ_FNAME)
	call strupr (curtask)
	call strupr (child)

	# Select direction
	switch (direction) {
	case G_FORWARD:
	    call sprintf (Memc[str], SZ_LINE, "FORK -- Forking to %s ...")
        	call pargstr (child)
	case G_BACKWARD:
	    call sprintf (Memc[str], SZ_LINE, "FORK -- Returning to %s ...")
        	call pargstr (curtask)
	default:
	    status = G_INTERNAL_ERROR
	    call sprintf (msg, SZ_LINE, "Unrecognized fork direction")
	    call sfree (sp)
	    call error (status, msg)
	}

	# Write to Gemini log file
	status = glw_str (gl, Memc[str], level)

	# Free memory
	call sfree (sp)

	return (status)
end

#--------------------------------------------------------------------------

# GLW_PARAM -- Write a task's parameter/value pairs to the Gemini log file.
# 	status = glw_param( gl, kvsbuf )
# 	
# 	status		: Status code		[return value, (int)]
# 	gl		: GL Structure		[input, (GL)]
# 	kvsbuf		: String buffer for key/value entries [input, (char[])]
# 
#     To write multiple key/value entries at once, each key/value entry must be 
#     delimited by newline characters, '\n'.

int procedure glw_param (gl, kvsbuf)

pointer gl              #I GL structure
char    kvsbuf[ARB]     #I String of key/value pair entries

int	status

# Other variables
char	token
int	i, nstr, kvsptr[G_MAX_PARAM]
int	tmpstatus
pointer	sp, str

# Gemini functions
int	glw_str(), g_splitstr()

errchk	glw_str()

begin
	# Set defaults and initial values
	status = 0
	token = '\n'

	# Allocate stack memory
	call smark (sp)
	call salloc (str, SZ_LINE, TY_CHAR)

	# Split at newlines
	nstr = g_splitstr (kvsbuf, token, kvsptr)

	# Write title
	call sprintf (Memc[str], SZ_LINE, "Input Parameters:")
	tmpstatus = glw_str (gl, Memc[str], TSK_LEVEL)
	status = status + tmpstatus

	# Write parameter/value entries
	for (i = 1; i <= nstr; i = i+1) {
	    call sprintf (Memc[str], SZ_LINE, "     %s")
        	call pargstr (kvsbuf[kvsptr[i]])
	    tmpstatus = glw_str (gl, Memc[str], TSK_LEVEL)
	    status = status + tmpstatus
	}

	# Free stack memory
	call sfree (sp)

	return (status)
end

#--------------------------------------------------------------------------

# GLW_STAT -- Write the exit status of a task.
# 	status = glw_stat( gl, statno, level )
# 	
# 	status		: Status code		[return value, (int)]
# 	gl		: GL structure		[input, (GL)]
# 	statno		: State code [G_SUCCESS|G_FAILURE  [input, (int)]
# 	level		: Level of the entry	[input, (int)]
# 
#     Suggested log level: STAT_LEVEL

int procedure glw_stat (gl, statno, level)

pointer gl              #I GL Structure
int     statno          #I State code
int     level           #I Level of the log entry

int	status

# Other variables
char	msg[SZ_LINE]
pointer	sp, str

# Gemini functions
int	glw_str()

errchk	glw_str()

begin
	status = 0

	# Allocate stack memory
	call smark (sp)
	call salloc (str, SZ_LINE, TY_CHAR)

	# Select exit status and format entry
	switch (statno) {
	case G_SUCCESS:
	    call sprintf (Memc[str], SZ_LINE, "Exit status: SUCCESS")
	case G_FAILURE:
	    call sprintf (Memc[str], SZ_LINE, "Exit status: FAILURE")
	default:
	    status = G_INTERNAL_ERROR
	    call sprintf (msg, SZ_LINE, "Unrecognized exit status code.\n")
	    call sfree (sp)
	    return (status)
	}

	# Write to Gemini log file
	status = glw_str (gl, Memc[str], level)

	# Free stack memory
	call sfree (sp)

	return (status)
end

#--------------------------------------------------------------------------
 
# GLW_TAG -- Write the begin/end tag used by GLOGEXTRACT.
# 	status = glw_tag( gl, tagmode )
# 	
# 	status		: Status code		[return value, (int)]
# 	gl		: GL structure		[input, (GL)]
# 	tagmode		: Type of tag (BEGIN_TAG|END_TAG) [input, (int)]

int procedure glw_tag (gl, tagmode)

pointer gl		#I GL structure
int	tagmode		#I Tag mode (BEGIN_TAG | END_TAG)

int	status

# Other variables
char	msg[SZ_LINE]
int	tm[LEN_TMSTRUCT], tmpverbose
long	ltime
pointer	sp, curtask, datestr, timestr, tag, tagstr

# Gemini functions
int	glw_str()

# IRAF functions
long	clktime()

errchk	glw_str()

begin
	status=0

	# Allocate stack memory for strings
	call smark (sp)
	call salloc (datestr, G_SZ_DATE, TY_CHAR)
	call salloc (timestr, G_SZ_TIME, TY_CHAR)
	call salloc (tag, G_SZ_LABEL, TY_CHAR)
	call salloc (tagstr, SZ_LINE, TY_CHAR)
	call salloc (curtask, SZ_FNAME, TY_CHAR)

	# Define date and time strings
	ltime = clktime (0)
	call brktime (ltime, tm)

	call sprintf (Memc[datestr], G_SZ_DATE, "%04d-%02d-%02d")
	    call pargi (TM_YEAR(tm))
	    call pargi (TM_MONTH(tm))
	    call pargi (TM_MDAY(tm))
	call sprintf (Memc[timestr], G_SZ_TIME, "%02d:%02d:%02d")
	    call pargi (TM_HOUR(tm))
	    call pargi (TM_MIN(tm))
	    call pargi (TM_SEC(tm))

	# Define the tag label
	if (tagmode == BEGIN_TAG) 
	    call strcpy ("BOE ",Memc[tag], G_SZ_LABEL)
	else if (tagmode == END_TAG)
	    call strcpy ("EOE ",Memc[tag], G_SZ_LABEL)
	else {
	    status = G_INTERNAL_ERROR
	    call sprintf (msg, SZ_LINE, "Unrecognized tag mode in glw_tag()")
	    call sfree (sp)
	    call error (status, msg)
	}

	# Get uppercase version of GL_CURTASK(gl)
	call strcpy (GL_CURTASK(gl), Memc[curtask], SZ_FNAME)
	call strupr (Memc[curtask])

	# Create tag string
	call sprintf (Memc[tagstr], SZ_LINE, "%s %s %sT%s")
	    call pargstr (Memc[tag])
	    call pargstr (Memc[curtask])
	    call pargstr (Memc[datestr])
	    call pargstr (Memc[timestr])

	# Write tag string (verbose disabled)
	tmpverbose = GL_VERBOSE(gl)
	GL_VERBOSE(gl) = NO
	status = glw_str (gl, Memc[tagstr], IGNORE_LEVEL) 
	GL_VERBOSE(gl) = tmpverbose

	# Free stack memory
	call sfree (sp)

	return (status)
end

#--------------------------------------------------------------------------

# GLW_TITLE -- Write log opening/closing information entry
# 	status = glw_title( gl, level, mode )
# 	
# 	status		: Status code		[return value, (int)]
# 	gl		: GL structure		[input, (GL)]
# 	level		: Level of the entry	[input, (int)]
# 	mode		: Type of entry (BEGIN_TAG|END_TAG) [input, (int)]
# 
#     Suggested log level: STAT_LEVEL

int procedure glw_title (gl, level, mode)

pointer gl              #I GL structure
int     level           #I Level of the log entry (suggested level: STAT)
int     mode            #I Beginning or end of log entries (BEGIN_TAG|END_TAG)

int	status

# Other variables
char	msg[SZ_LINE]
long	ltime
pointer	sp, timestr, str

# Gemini functions
int	glw_str(), glw_vis()

# IRAF functions
int	errget()
long	clktime()

errchk	glw_str(), glw_vis()

begin
	status = 0
	call strcpy ("", msg, SZ_LINE)

	# Allocate stack memory
	call smark (sp)
	call salloc (str, SZ_LINE, TY_CHAR)
	call salloc (timestr, SZ_TIME, TY_CHAR)

	# Get the time and convert to long string
	ltime = clktime(0)
	call cnvtime (ltime, Memc[timestr], SZ_TIME)    

	# Write title string to log
	iferr {
	    if (mode == BEGIN_TAG) {
        	call sprintf (Memc[str], SZ_LINE, "Log opened at [%s]")
        	    call pargstr (Memc[timestr])
		status = glw_vis (gl, G_EMPTY)
		status = glw_vis (gl, G_LONG_DASH)
		status = glw_str (gl, Memc[str], level)
	    } else if (mode == END_TAG) {
        	call sprintf (Memc[str], SZ_LINE, "Log closed at [%s]")
        	    call pargstr (Memc[timestr])
		status = glw_str (gl, Memc[str], level)
		status = glw_vis (gl, G_LONG_DASH)
		status = glw_vis (gl, G_EMPTY)
	    } else {
        	status = G_INTERNAL_ERROR
        	call sprintf (msg, SZ_LINE,
		    "Unrecognized tag mode in glw_title()")
	    }
	} then
	    status = errget (msg, SZ_LINE)

	# Free stack memory
	call sfree (sp)

	if (status >= G_INTERNAL_ERROR)
	    call error (status, msg)

	return (status)
end

#--------------------------------------------------------------------------

# GLW_VIS -- Write visual/log-readability improvement entries
# 	status = glw_vis( gl, type )
# 	
# 	status		: Status code		[return value, (int)]
# 	gl		: GL structure		[input, (GL)]
# 	type		: Type of string	[input, (int)]
# 
#     Type of strings are: G_EMPTY, G_LONG_DASH, G_SHRT_DASH, for an empty line,
#     a long series of '-' characters, and a short series (20) of '-' characters,
#     respectively.

int procedure glw_vis (gl, type)

pointer gl              #I GL structure
int     type            #I Type of string (G_EMPTY,G_LONG_DASH,G_SHRT_DASH)

int	status

# Other variables
char	msg[SZ_LINE]
int	i, len
pointer	sp, str

# Gemini functions
int	glw_str()

# IRAF functions
int	strlen()

errchk	glw_str()

begin
	status = 0

	# Allocate stack memory
	call smark (sp)
	call salloc (str, SZ_LINE, TY_CHAR)

	# Format entry
	call strcpy ("", Memc[str], SZ_LINE)    
	switch (type) {
	    case G_EMPTY:
	    ;
	case G_LONG_DASH:
	    len = 5 + strlen (GL_CURTASK(gl)) + 1
	    for (i = 1; i <= 80-len; i = i+1) {
        	call strcat ("-", Memc[str], SZ_LINE)
	    }
	case G_SHRT_DASH:
	    for (i = 1; i <= 20; i = i+1) {
        	call strcat ("-", Memc[str], SZ_LINE)
	    }
	default:
	    status = G_INTERNAL_ERROR
	    call sprintf (msg, SZ_LINE, 
	        "Unrecognized type of readability improvement")
	    call sfree (sp)
	    call error (status, msg)
	}

	# Write the string to log
	status = glw_str (gl, Memc[str], VIS_LEVEL)

	# Free stack memory
	call sfree (sp)

	return (status)
end

#--------------------------------------------------------------------------

# GLW_WARN -- Write a warning message to the Gemini Log file.
#         status = glw_warn( gl, errno, errmsg, level )
#         
#         status          : Status code   [return value, (int)]
#         gl              : GL Structure  [input, (GL)]
#         errno           : Error code    [input, (int)]
#         usrmsg          : User's error message [input, (char[])]
#         level           : Level of the entry   [input, (int)]
# 
#     If "errno" is not NULL, GLW_ERR will use it to prepend the official error
#     message to "usrmsg".

int procedure glw_warn (gl, errno, usrmsg, level)

pointer gl              #I GL structure
int     errno           #I Error code (optional - set it to NULL if not used)
char    usrmsg[ARB]     #I Error message (optional)
int     level           #I Level of the log entry

int	status

# Other variables
char	token, msg[SZ_LINE], curtask[SZ_FNAME]
int	i, nstr, strptr[G_MAX_LINES]
int 	tmpstatus
pointer	sp, errmsg, tmpstr

# Gemini functions
int	glw_str(), g_splitstr()
bool	g_whitespace()

# IRAF functions
int	errget()

errchk	glw_str()

begin
	status = 0
	tmpstatus = 0

	# Allocate stack memory
	call smark (sp)
	call salloc (errmsg, SZ_LINE, TY_CHAR)
	call salloc (tmpstr, SZ_LINE, TY_CHAR)

	# Initialize
	token = '\n'

	# Fetch Gemini error message
	if ((errno >= 99) && (errno < 500)) {
	    ifnoerr ( call gemerrmsg (errno, Memc[tmpstr]) ) {
		call sprintf( Memc[errmsg], SZ_LINE, "WARNING: %d %s")
		    call pargi (errno)
		    call pargstr (Memc[tmpstr])
        	tmpstatus = glw_str (gl, Memc[errmsg], level)
		status = status + tmpstatus
	    } else {
		tmpstatus = errget( msg, SZ_LINE )
		if (tmpstatus == G_INPUT_ERROR) {
		    call strcpy (GL_CURTASK(gl), curtask, SZ_FNAME)
		    call strupr (curtask)
		    call printf ("%s WARNING: %d %s\n")
			call pargstr (curtask)
			call pargi (tmpstatus)
			call pargstr (msg)
		} else {
		    status = errget (msg, SZ_LINE)
		    call sfree (sp)
		    call error (status, msg)
		}
	    }
	} else if (errno == NULL) {
	    errno = 1
	} else if (errno < 0) {  # Use errno but ignore Gemini error msg
	    errno = -errno
	}

	# Write the user's warning message.
	if ( ! g_whitespace (usrmsg) ) {
	
	    # Split strings on the newline character, '\n'
	    nstr = g_splitstr (usrmsg, token, strptr)

	    for (i = 1; i <= nstr; i = i+1) {
        	call sprintf (Memc[errmsg], SZ_LINE, "WARNING: %d %s")
        	    call pargi (errno)
        	    call pargstr (usrmsg[strptr[i]])
        	tmpstatus = glw_str (gl, Memc[errmsg], level)
		status = status + tmpstatus
	    }
	}

	# Free stack memory
	call sfree (sp)

	return (status)
end

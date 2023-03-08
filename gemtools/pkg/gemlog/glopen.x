# Copyright(c) 2004-2006 Association of Universities for Research in Astronomy, Inc.

# Author: Kathleen Labrie  3-May-2004

include "glog.h"
include "gemerrors.h"
include <error.h>

.help
.nf
GL_OPEN -- Create a Gemini log file structure and open the logfile
	gl = gl_open ( logfile, acmode, status )
	
	gl		: GL Structure  [return value, (GL)]
	logfile		: Name of the logfile  [input, (string)]
	acmode		: File access mode (e.g. APPEND)  [input, (int)]
	status		: Status  [output, (int)]
.fi
.endhelp

#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
#--------------------------------------------------------------------------

# GL_OPEN -- Create a Gemini log file structure and open the logfile
# 	gl = gl_open ( logfile, acmode, status )
# 	
# 	gl		: GL Structure  [return value, (GL)]
# 	logfile		: Name of the logfile  [input, (string)]
# 	acmode		: File access mode (e.g. APPEND)  [input, (int)]
# 	status		: Status  [output, (int)]

pointer procedure gl_open ( logfile, acmode, status )

char	logfile[ARB]		#I Name of the logfile
int	acmode			#I Logfile access mode
int	status			#I Status code

pointer gl

char	msg[SZ_LINE]

int 	access(), open(), strcmp()

begin
	status = 0
	
	# Allocate the GL structure
	call glalloc (gl)

	# Fill the string structure elements
	call strcpy (logfile, GL_LOGFILE(gl), SZ_FNAME)
	call strcpy ("INDEF", GL_CURPACK(gl), SZ_FNAME)
	call strcpy ("INDEF", GL_CURTASK(gl), SZ_FNAME)

	# Set GL defaults
	GL_FD(gl) = NULL
	GL_REQSTAT(gl) = YES
	GL_REQSCI(gl)  = YES
	GL_REQENG(gl)  = YES
	GL_REQVIS(gl)  = YES
	GL_REQTSK(gl)  = YES
	GL_VERBOSE(gl) = YES

	if ( strcmp("STDOUT", GL_LOGFILE(gl)) == 0 )
	    GL_FD(gl) = STDOUT
	else if ( strcmp("STDERR", GL_LOGFILE(gl)) == 0 )
	    GL_FD(gl) = STDERR
	else if ( strcmp("dev$null", GL_LOGFILE(gl)) == 0 )
	    GL_FD(gl) = open (GL_LOGFILE(gl), APPEND, TEXT_FILE)
	else {
	    # Determine the accessbility of the logfile
	    if (acmode == NEW_FILE) {
		if ( access (GL_LOGFILE(gl), 0, 0) == YES ) {
 		    status = G_FILE_EXISTS
  		    call sprintf (msg, SZ_LINE, 
  			"Cannot create new file, %s.  Already exists.")
  			call pargstr ( GL_LOGFILE(gl) )
        	    call gl_close (gl)
 		    call error (status, msg )
		}
	    }
	    else if ( access (GL_LOGFILE(gl), acmode, TEXT_FILE) == NO ) {
		status = G_FILE_NOT_ACCESSIBLE
		call sprintf (msg, SZ_LINE, "Unable to access file '%s'.")
  		    call pargstr (GL_LOGFILE(gl))
		call gl_close (gl)
		call error (status, msg)
	    }

	    # If we get here, then file is accessible or can be created.
	    iferr ( GL_FD(gl) = open (GL_LOGFILE(gl), acmode, TEXT_FILE) ) {
	        GL_FD(gl) = NULL
		call gl_close (gl)
		call erract (EA_ERROR)
	    }
	}

	# Return the structure
	return (gl)
end

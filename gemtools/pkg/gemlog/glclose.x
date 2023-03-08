# Copyright(c) 2004-2006 Association of Universities for Research in Astronomy, Inc.

# Author: Kathleen Labrie  3-May-2004

include "glog.h"

.help
.nf
GL_CLOSE -- Close the logfile and free the memory associated with the GL
            structure.
	gl_close( gl )

	gl	: GL Structure  [input (GL)]
.fi
.endhelp

#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
#--------------------------------------------------------------------------

# GL_CLOSE -- Close the logfile and free the memory associated with the GL
#             structure.
# 	gl_close( gl )
# 
# 	gl	: GL Structure  [input (GL)]


procedure gl_close (gl)

pointer gl			#I GL Structure

int	l_status
char	msg[SZ_LINE]

int	errget()

begin
	l_status = 0

	# Close the log file
	if ( GL_FD(gl) != NULL ) {
	    call flush (GL_FD(gl))
	    iferr ( call close (GL_FD(gl)) ) {
		l_status = errget (msg, SZ_LINE)
	    }
	    ## ??? catch error ##
	    ## If there is an error I still want to free the memory, and
	    ## then returns with a catchable error.
	    ## What's the best way to do that?
	    ## How's the current implementation?
	}

	# Free GL structure
	call glfree (gl)

	if (l_status != 0)
	    call error (l_status, msg)

	return
end

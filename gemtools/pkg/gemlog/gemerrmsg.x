# Copyright(c) 2004-2006 Association of Universities for Research in Astronomy, Inc.

include "gemerrors.h"

.help
.nf
GEMERRMSG -- Given a Gemini error code, get Gemini error message.  The error
              message file is './gemerrmsg'
	call gemerrmsg ( errno, errmsg )
	
	errno		: Gemini error code	[input, (int)]
	errmsg		: Gemini error message	[output, (char[])]

.fi
.endhelp

# Kathleen Labrie   3-May-2004

procedure gemerrmsg ( errno, errmsg )

int	errno			#I Error code
char	errmsg[SZ_LINE]		#O Error message associated with errno

int 	fd, code, i, status
char	msg[SZ_LINE], fullname[SZ_FNAME], buf[SZ_FNAME]
bool	found

string	MSG_FILE	"pkg/gemlog/gemerrmsg"  # relative to gemtools$

int	access(), open(), fscan(), envfind(), errget()

errchk	gargi()

begin
	# Initialize variables
	status = 0 
	found = FALSE
	call strcpy ("", errmsg, SZ_LINE)
	call strcpy ("", fullname, SZ_FNAME)
	call strcpy ("", buf, SZ_FNAME)
	call strcpy ("", msg, SZ_LINE)

	# Generate the OS pathname of the "gemerrmsg" file
	if ( envfind ("gemtools", buf, SZ_FNAME) > 0 ) {
	    call strcat (MSG_FILE, buf, SZ_FNAME)
	    call fpathname (buf, fullname, SZ_FNAME)
	} else {
	    call sprintf (msg, SZ_LINE, 
	        "'gemtools' environment variable not defined'")
	    call error ( G_INTERNAL_ERROR, msg )
	}

	# Open the 'gemerrmsg' file
	if ( access ( fullname, READ_ONLY, TEXT_FILE) == YES ) {
	    iferr ( fd = open ( fullname, READ_ONLY, TEXT_FILE ) ) {
		status = errget (msg, SZ_LINE)
		call error (status, msg)
	    }
	} else {
	    call sprintf (msg, SZ_LINE, "'gemerrmsg' not found")
	    call error (G_INTERNAL_ERROR, msg)
	}

	# Find the correct line (based on code), and retrieve error msg.
	while (fscan (fd) != EOF) {
	    ifnoerr {
        	call gargi (code)
		call gargstr (msg, SZ_LINE)
	    } then {
		if (code == errno) {
		    i = 0
		    repeat
			i = i+1
		    until (msg[i] != 32)
		    call strcpy (msg[i], errmsg, SZ_LINE)
		    found = TRUE
		    break
		}
	    } else {
		call sprintf (msg, SZ_LINE, "Error reading 'gemerrmsg'")
		call error (G_INTERNAL_ERROR, msg)
	    }
	}

	# The code was not found
	if (found == FALSE) {
	    call sprintf (msg, SZ_LINE, "Unrecognized Gemini error code (%d)")
		call pargi(errno)
	    call error (G_INPUT_ERROR, msg)
	}
	
	# Close up, we're done
	call close (fd)

	return
end

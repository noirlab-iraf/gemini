# Copyright(c) 2004-2009 Association of Universities for Research in Astronomy, Inc.

include "glog.h"
include <fset.h>
include <fio.h>

# GLOGFIX -- Fix a Gemini log file by adding missing EOE tags.
#
#	status = glogfix ( gl )
#
#	status		: Exit status code		[return value, (int)]
#	gl		: GL structure			[input, (GL)]

int procedure glogfix( gl )

pointer	gl	#I GL structure

int	status

# Other variables
char	tag[LEN_LEVEL_STR]
int	nboe, i, mode
pointer	sp, buffer, tmpcurtask

# Gemini functions
int	glw_tag()

# IRAF functions
int	getline(), strcmp(), open()
int	fstati()

begin
	status = 0

	# Allocate stack memory
	call smark (sp)
	call salloc (buffer, SZ_LINE, TY_CHAR)
	call salloc (tmpcurtask, SZ_FNAME, TY_CHAR)

	mode = fstati (GL_FD(gl), F_MODE)

	if ((mode != READ_ONLY) && (mode != READ_WRITE)) {
	    call close (GL_FD(gl))
	    GL_FD(gl) = open (GL_LOGFILE(gl), READ_WRITE, TEXT_FILE)
	} else {
	    call seek (GL_FD(gl), BOFL)
	}
	
	nboe = 0
	while ( getline(GL_FD(gl), Memc[buffer]) != EOF ) {
	    call sscan (Memc[buffer])
	        call gargwrd (tag, LEN_LEVEL_STR)
	    if ( strcmp ("BOE", tag) == 0 )
	        nboe = nboe + 1
	    else if ( strcmp ("EOE", tag) == 0 )
	        nboe = nboe - 1
	    else
	    	next
	}
	
	if ( nboe != 0 ) {	# logfile needs fixing!
	    call strcpy (GL_CURTASK(gl), Memc[tmpcurtask], SZ_FNAME)
	    call strcpy ("glogfix", GL_CURTASK(gl), SZ_FNAME)
	    for (i = 1; i <= nboe; i = i+1) {
	        status = glw_tag (gl, END_TAG)
	    }
	    call strcpy (Memc[tmpcurtask], GL_CURTASK(gl), SZ_FNAME)
	}

	# Close the file and reopen it with the access mode it had a the top.
	if ((mode != READ_ONLY) && (mode != READ_WRITE)) {
	    call close (GL_FD(gl))
	    GL_FD(gl) = open (GL_LOGFILE(gl), mode, TEXT_FILE)
	}

	# Free memory
	call sfree (sp)
	
	return (status)
end

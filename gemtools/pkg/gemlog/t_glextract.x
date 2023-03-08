# Copyright(c) 2004-2006 Association of Universities for Research in Astronomy, Inc.

# Author: Kathleen Labrie  3-May-2004

# GLEXTRACT -- Extract log entries from a Gemini log file

include "glog.h"
include <mach.h>

define	GETLAST	-2

procedure t_glextract ()

#char	logfile[SZ_FNAME]	# Name of the log file

#char	taskname[SZ_FNAME]	# Retrieve logs for this task
#char	block[SZ_LINE]		# Which block(s)? (e.g. last, 2...)
#char	ltime[SZ_LINE]		# Lower limit on the time range
#char	utime[SZ_LINE]		# Upper limit on the time range
#bool	fl_child		# Include logs for the child processes?
#int	nchild			# Maximum number of subprocess levels
#bool	fl_stat_level		# Retrieve state level logs?
#bool	fl_tsk_level		# Retrieve task level logs?
#bool	fl_sci_level		# Retrieve science level logs?
#bool	fl_eng_level		# Retrieve engineering level logs?
#bool	fl_vis_level		# Retrieve visual (readability) level logs?
#int	status			# Exit status (0=good)

# Local variables for task parameters
char	l_logfile[SZ_FNAME]
char	l_taskname[SZ_FNAME], l_block[SZ_LINE]
char	l_ltime[SZ_LINE], l_utime[SZ_LINE]
int	l_fl_child,l_nchild
int	l_fl_stat_level, l_fl_sci_level, l_fl_eng_level, l_fl_vis_level
int	l_fl_tsk_level
int	l_status

# Other variables
int	i, nlines, blkwanted, tindex, lastnlines
bool	wantlastblk
pointer	gl, sl, sp, blkbuf, tmpstr, lastblk

# Gemini functions
int	glr_blk(), prs_blk(), prs_time()
bool	g_whitespace()
pointer	gl_open()

# IRAF functions
bool	clgetb()
int	btoi(), stridx(), clgeti(), errget()

errchk	prs_blk(), prs_time()

begin
	l_status=0

	# Get task parameter values
	call clgstr ("logfile", l_logfile, SZ_FNAME)
	call clgstr ("taskname", l_taskname, SZ_FNAME)
	call clgstr ("block", l_block, SZ_LINE)
	call clgstr ("ltime", l_ltime, SZ_LINE)
	call clgstr ("utime", l_utime, SZ_LINE)
	l_fl_child = btoi (clgetb ("fl_child"))
	l_nchild = clgeti ("nchild")
	l_fl_stat_level = btoi (clgetb ("fl_stat_level"))
	l_fl_tsk_level  = btoi (clgetb ("fl_tsk_level"))
	l_fl_sci_level  = btoi (clgetb ("fl_sci_level"))
	l_fl_eng_level  = btoi (clgetb ("fl_eng_level"))
	l_fl_vis_level  = btoi (clgetb ("fl_vis_level"))

	# Allocate stack memory
	call smark (sp)
	call salloc (blkbuf, SZ_LINE*G_MAX_ENTRY_PER_BLK, TY_CHAR)
	call salloc (tmpstr, SZ_LINE, TY_CHAR)

	# Allocate memory for SL Structure
	call slalloc (sl)

	# Initialize
	SL_LTIME(sl) = 0
	SL_UTIME(sl) = 0
	SL_NCHILD(sl) = MAX_INT
	wantlastblk = FALSE

	# Open the log file for reading only
	iferr ( gl = gl_open (l_logfile, READ_ONLY, l_status) ) {
	    gl = NULL
	    l_status = errget (Memc[tmpstr], SZ_LINE)
	    call printf ("GLOGEXTRACT ERROR: %d %s.\n")
		call pargi (l_status)
		call pargstr (Memc[tmpstr])
	    call clputi ("status", l_status)
	    call slfree (sl)
	    call sfree (sp)
	    return
	}

	# Fill the GL structure
	GL_REQSTAT(gl) = l_fl_stat_level
	GL_REQSCI(gl)  = l_fl_sci_level
	GL_REQENG(gl)  = l_fl_eng_level
	GL_REQVIS(gl)  = l_fl_vis_level
	GL_REQTSK(gl)  = l_fl_tsk_level
	GL_VERBOSE(gl) = YES	# We want glr_blk() to print to STDOUT

	# Fill the SL structure  (includes prs_blk() and prs_time())
	call strupr (l_taskname)
	call strcpy (l_taskname, SL_TSKNAME(sl), SZ_FNAME)
	SL_NBLK(sl) = 0
	SL_CHILD(sl) = l_fl_child
	if (l_nchild != INDEFI)
	    SL_NCHILD(sl) = l_nchild

	# Parse 'block, 'ltime', and 'utime'
	iferr {
	
	    # "block"  (Create block position array, and sort - ascending)
	    for (i = 1; i <= G_MAX_BLK; i = i+1)
        	SL_BLKS(sl,i) = G_INDEF		#Get all (default)
	    SL_BPOS(sl) = 0			#Set counter to zero

	    if ( ! g_whitespace(l_block) )
		l_status = l_status + 
		    prs_blk (l_block, SL_NBLK(sl), SL_BLKS(sl,1))
	    else
		SL_NBLK(sl) = G_MAX_BLK		#Get all the good blocks

	    # "ltime/utime" (Convert ltime and utime to integer, if defined)
	    if ( ! g_whitespace (l_ltime) )
        	l_status = l_status + prs_time (l_ltime, SL_LTIME(sl))
	    if ( ! g_whitespace (l_utime) ) {
        	l_status = l_status + prs_time (l_utime, SL_UTIME(sl))

		# Check to see if a time has specified (look for the 'T')
        	# YYYY-MM-DDT -> will be interpreted as YYYY-MM-DDT00:00:00
        	# If time not specified (no 'T'), add 24 hours - 1sec.  
		# This way the user can ask for whole days. Eg. : 
		#	ltime=2004-04-19, utime=2004-04-19  => all logs for the
		#					       19th
		#	ltime=2004-04-19, utime=2004-04-21  => logs from 19 to 
		#					       21, incl.
		
		tindex = stridx ("T",l_utime)
		if (tindex == 0)
		    SL_UTIME(sl) = SL_UTIME(sl) + (24*3600) - 1
	    }
	} then {
	    l_status = errget (Memc[tmpstr], SZ_LINE)
	    call printf ("GLOGEXTRACT ERROR: %d %s\n")
		call pargi (l_status)
		call pargstr (Memc[tmpstr])
	    call gl_close (gl)
	    call clputi ("status", l_status)
	    call slfree (sl)
	    call sfree (sp)
	    return
	} else {
	    if (l_status != 0) {
		call printf ("GLOGEXTRACT ERROR: 1 Error extracting the log entries.\n")
		call gl_close (gl)
		call clputi ("status", l_status)
		call slfree (sp)
		call sfree (sp)
		return
	    }
	}

	# Check if 'last' block was requested
	if (SL_BLKS(sl, SL_NBLK(sl)) == G_MAX_BLK) {
	   wantlastblk = TRUE
	   
	   # Set up to let the loop go to EOF, so that the very last good 
	   # block will be retrieved.
	   
	   SL_BLKS(sl, SL_NBLK(sl)) = GETLAST
	   SL_NBLK(sl) = G_MAX_BLK	# this will let the loop go to EOF
	   call malloc (lastblk, SZ_LINE*G_MAX_ENTRY_PER_BLK, TY_CHAR)
	}

	# Loop through SL_NBLK(sl)
	for (i = 1; i <= SL_NBLK(sl); i = i+1) {
	    blkwanted = SL_BLKS(sl, i)
	    if (blkwanted == GETLAST) {
	        # don't print anything until the last blk
		blkwanted = G_INDEF
		GL_VERBOSE(gl) = NO
	    }
	    iferr (l_status = glr_blk(gl, sl, blkwanted, nlines, Memc[blkbuf])){
		l_status = errget (Memc[tmpstr], SZ_LINE)
		call printf ("GLOGEXTRACT ERROR: %d %s.\n")
		    call pargi (l_status)
		    call pargstr (Memc[tmpstr])
		call gl_close (gl)
		call clputi ("status", l_status)
		if (wantlastblk)
		    call mfree (lastblk, TY_CHAR)
		call slfree (sp)
		call sfree (sp)
		return
	    }
	    if (nlines == 0) {	# No more blocks found, we're almost done
		break
	    }
	    if (wantlastblk) {  # Keep a copy of the last good block
		call bufcpyc (Memc[blkbuf], Memc[lastblk], 
		    SZ_LINE*G_MAX_ENTRY_PER_BLK)
		lastnlines = nlines
	    }
	}

	if (wantlastblk) {
	    # Reset GL_VERBOSE(gl) to yes and print last block.
	    # (It was set no when the for-loop reached the GETLAST block.)
	    
	    GL_VERBOSE(gl) = YES
	    call printlast (lastnlines, Memc[lastblk])
	    call mfree (lastblk, TY_CHAR)
	}

	# Close the log file
	iferr ( call gl_close(gl) ) {
	    l_status = errget (Memc[tmpstr], SZ_LINE)
	    call printf ("GLOGEXTRACT ERROR: %d %s.\n")
		call pargi (l_status)
		call pargstr (Memc[tmpstr])
	}

	# Free memory
	call slfree(sl)
	call sfree(sp)

	return
end

##################################
# Utilities
procedure bufcpyc (inbuffer, outbuffer, size)

char	inbuffer[size]
char	outbuffer[size]
int	size

int	i

begin
	for (i = 1; i <= size; i = i+1)
	    outbuffer[i] = inbuffer[i]

	return
end

procedure printlast (nlines, strbuf)

int	nlines
char 	strbuf[SZ_LINE,nlines]

int	i

begin
	for (i = 1; i <= nlines; i = i+1) {
	    call printf ("%s")
		call pargstr (strbuf[1,i])
	}

	return
end

# Copyright(c) 2004-2006 Association of Universities for Research in Astronomy, Inc.

# Author: Kathleen Labrie  3-May-2004

include "glog.h"
include "gemerrors.h"
include <ctype.h>
include <time.h>

.help
.nf
This file contains procedures used to parse GEMLOG user inputs.

    prs_blk - Parse list of blocks to retrieve (glogextract.block)
   prs_time - Parse time string (glogextract.ltime; glogextract.utime)

PRS_BLK -- Parse the user's list of blocks to retrieve.
	status = prs_blk( blkstr, nblk, blks )
	
	status		: Exit status code		[return value, (int)]
	blkstr		: Block parameter string	[input, (char[])]
	nblk		: Number of blocks to retrieve	[output, (int)]
	blks		: Parsed list of blocks		[output, (int[])]

PRS_TIME -- Parse user's lower/upper limit to the valid time range.
	status = prs_time( timestr, ltime )
	
	status		: Exit status code		[return value, (int)]
	timestr		: Time string (from user)	[input, (char[])]
	ltime		: Seconds since 1-1-1980	[output, (long)]

.fi
.endhelp

define 	SECPHR	3600
define	SECPMIN	60
define	SECPSEC	1

#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
#--------------------------------------------------------------------------

# PRS_BLK -- Parse the user's list of blocks to retrieve.
# 	status = prs_blk( blkstr, nblk, blks )
# 	
# 	status		: Exit status code		[return value, (int)]
# 	blkstr		: Block parameter string	[input, (char[])]
# 	nblk		: Number of blocks to retrieve	[output, (int)]
# 	blks		: Parsed list of blocks		[output, (int[])]

int procedure prs_blk( blkstr, nblk, blks )

char	blkstr[ARB]             #I Block parameter string
int	nblk                    #O Number of blocks to retrieve
int	blks[ARB]               #O Numerical list of blocks to retrieve

int	status

char	msg[SZ_LINE]
int	ip, tmpip, lblkstr, len, tmpblk, b
int	tblks[G_MAX_BLK], tnblk
bool	f_range, f_valid

int	ctoi(), strlen(), strncmp()

begin
	# Initialize
	status = 0
	
	tnblk = 0
	lblkstr = strlen(blkstr)
	f_range = FALSE
	f_valid = FALSE

	# Go through the length of blkstr
	ip = 1
	while (ip <= lblkstr) {
	    f_valid = FALSE				# Reset flag

	    # Skip white spaces and commas
	    while ( IS_WHITE( blkstr[ip] ) || blkstr[ip] == ',' )
        	ip = ip + 1

	    # Decode
	    if ( IS_DIGIT( blkstr[ip] ) ) {
        	f_valid = TRUE              		# Valid value
        	len = ctoi (blkstr, ip, tmpblk)    	# len is added to ip
	    } else if ( blkstr[ip] == '-' ) {
        	ip = ip + 1
        	f_range = TRUE              	    # Next value is end of range
	    } else {
        	tmpip = ip
        	while ( IS_ALPHA( blkstr[tmpip] ) )	# 'first' or 'last' ?
        	    tmpip = tmpip + 1
        	len = tmpip - ip			# len is added to ip
        	if ( strncmp( blkstr[ip], "first", len ) == 0) {
        	    f_valid = TRUE
        	    tmpblk = 1
        	} else if ( strncmp( blkstr[ip], "last", len ) == 0) {
        	    f_valid = TRUE
        	    tmpblk = G_MAX_BLK		    # not the real last position
        	} else {
		    status = G_INPUT_ERROR
		    call sprintf (msg, SZ_LINE, "Invalid block string (%s)")
			call pargstr(blkstr[ip])
		    call error (status, msg)
        	}
        	ip = ip + len
	    }

	    # Assign positions to blks
	    if ( f_valid ) {
        	if ( f_range ) {
        	    for (b = tblks[tnblk] + 1; b <= tmpblk; b = b+1) {
                	tnblk = tnblk + 1
                	tblks[tnblk] = b
        	    }
        	    f_range = FALSE         		#Reset flag
        	} else {
        	    tnblk = tnblk + 1
        	    tblks[tnblk] = tmpblk
        	}
	    }

	}

	# Sort block position vector
	call asrti (tblks, blks, tnblk)
 
	# Remove duplicates
	if (tnblk > 1) {
	    # Copy blks back to work vector, tblks
	    call amovi (blks, tblks, tnblk)

	    # Take care of the first element
	    blks[1] = tblks[1]
	    nblk = 1

	    # Check the others
	    for (b = 2; b <= tnblk; b = b+1) {
        	if (tblks[b] != tblks[b-1]) {
        	    nblk = nblk + 1
        	    blks[nblk] = tblks[b]
        	}
	    }
	} else
	    nblk = tnblk

	return (status)
end

#--------------------------------------------------------------------------

# PRS_TIME -- Parse user's lower/upper limit to the valid time range.
# 	status = prs_time( timestr, ltime )
# 	
# 	status		: Exit status code		[return value, (int)]
# 	timestr		: Time string (from user)	[input, (char[])]
# 	ltime		: Seconds since 1-1-1980	[output, (long)]

int procedure prs_time( timestr, ltime )

char	timestr[ARB]		#I Time string (from user)
long	ltime                   #O Seconds since 1-1-1980

int     status

char	msg[SZ_LINE]
int	index, i, n, ip, nchar
int	tm[LEN_TMSTRUCT], timecnv[3], p_buf[3]
long	now
real	rval
pointer	sp, strbuf, dstr, tstr

int	g_splitstr()
bool	g_whitespace()

int	stridx(), dtm_ltime(), strlen(), ctor()
long	clktime()

begin
	# Initialize
	status = 0
	ltime = 0

	# Allocate stack memory
	call smark (sp)
	call salloc (strbuf, SZ_LINE, TY_CHAR)
	call salloc (dstr, SZ_LINE, TY_CHAR)
	call salloc (tstr, SZ_LINE, TY_CHAR)

	# Start parsing ...

	# First try the normal FITS format
	status = dtm_ltime (timestr, ltime)

	if (status == ERR) {	# Not a date or dateTtime FITS format
	    #Reset status; the string might still be a valid gemlog time format.
	    status = 0

	    #Split date and time at the 'T'.  If no 'T', then has to be a date.
	    index = stridx ("T",timestr)
	    if (index != 0) {
        	call sscan (timestr)
        	    call gargstr (Memc[dstr], index-1)
        	    call gargstr (Memc[strbuf], 1)
        	    call gargstr (Memc[tstr], SZ_LINE)
	    } else {
        	call strcpy (timestr, Memc[dstr], SZ_LINE)
        	call strcpy ("", Memc[tstr], SZ_LINE)
	    }

	    # Deal with the date.
	    # Is there a date string and is it a valid format?
	    
	    if ( ! g_whitespace (Memc[dstr]) ) {   	#yes, there is a date
        	status = dtm_ltime (Memc[dstr], ltime)
        	if (status == ERR) {       	#not a valid date => INPUT ERROR
		    status = G_INPUT_ERROR
        	    call sprintf (msg, SZ_LINE, "Invalid date string (%s)")
                	call pargstr (Memc[dstr])
        	    call sfree (sp)
        	    call error (status, msg)
        	}
	    } else {				# there is no date, use today
        	now = clktime(0)
        	call brktime (now,tm)
        	call sprintf (Memc[dstr], 10, "%04d-%02d-%02d")
        	    call pargi (TM_YEAR(tm))
        	    call pargi (TM_MONTH(tm))
        	    call pargi (TM_MDAY(tm))
        	status = dtm_ltime (Memc[dstr], ltime)
        	if (status == ERR) {		# INTERNAL ERROR !!!
		    status = G_INTERNAL_ERROR
        	    call sprintf (msg, SZ_LINE, 
			"Internal error - Invalid date string (%s)")
                	call pargstr (Memc[dstr])
        	    call sfree (sp)
        	    call error (status, msg)
        	}
	    }

	    # Deal with the time
	    # Is there a time string, and is it a valid format?
	    
	    if ( ! g_whitespace (Memc[tstr]) ) {   # yes, time string defined
        	timecnv[1] = SECPHR
        	timecnv[2] = SECPMIN
        	timecnv[3] = SECPSEC
        	call strcpy (Memc[tstr], Memc[strbuf], SZ_LINE)
        	n = g_splitstr (Memc[strbuf], ":", p_buf)
		if (n > 3) {
		    status = G_INPUT_ERROR
		    call sprintf (msg, SZ_LINE, "Invalid time string (%s)")
			call pargstr (Memc[tstr])
		    call sfree (sp)
		    call error (status, msg)
		}
        	for (i = 1; i <= n; i = i+1) {
        	    ip = 1
        	    nchar = ctor (Memc[strbuf+p_buf[i]-1], ip, rval)
        	    if ( nchar != strlen (Memc[strbuf+p_buf[i]-1]) ) {
			status = G_INPUT_ERROR
                	call sprintf (msg, SZ_LINE, "Invalid time string (%s)")
                	    call pargstr (Memc[tstr])
                	call sfree (sp)
                	call error (status, msg)
        	    }
        	    ltime = ltime + int(rval*timecnv[i])
        	}
	    } 
	    #else : there is no time string and we are done.
	    #      Side effect: If user enters 'T' only, => 0 hours today.

	}

	# Free stack memory
	call sfree (sp)

	return (status)
end

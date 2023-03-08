# Copyright(c) 2004-2006 Association of Universities for Research in Astronomy, Inc.

# Author: Kathleen Labrie  3-May-2004

include "glog.h"
include <mach.h>

.help
.nf
This file contains procedures used to select entries in a Gemini log file. 
These procedures are used by GLOGEXTRACT.

     blkchk - compare requested log block position to current block position
     lvlchk - compare entry's log level to requested levels
    taskchk - compare entry's task tag to the requested task name
    timechk - compare BOE time tag to the time range requested

BLKCHK -- Compare requested log block position position to current block 
          position.  The position is absolute and is counted from the 
	  beginning of the log file; it does not take in consideration the 
          other selection criteria.
	okay = blkchk( wanted, curblk )
	
	okay		: Found a request block? [return value, (bool)]
	wanted		: Block position to search for  [input, (int)]
	curblk		: Current block position	[input, (int)]

    If no specific block is requested, 'wanted' should be equal to 'G_INDEF'.

LVLCHK -- Compare an entry's log level to the list of requested levels.
	okay = lvlchk( gl, line )
	
	okay		: Retrieve this entry?	[return value, (bool)]
	gl		: GL structure		[input, (GL)]
	line		: Log entry		[input, (char[])]

TASKCHK -- Compare an entry's task tag to the requested task name.
	okay = taskchk( sl, line )
	
	okay		: Retrieve this entry?	    [return value, (bool)]
	sl		: SL (selection) structure  [input, (SL)]
	line		: Log entry		    [input, (char[])]		:

TIMECHK -- Compare the time tag of the log block to the time range requested.
	okay = timechk( sl, line )
	
	okay		: Current block within time range? [return value, (bool)]
	sl		: SL (selection) structure	[input, (SL)]
	line		: Log entry			[input, (char[])]
.fi
.endhelp

#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
#--------------------------------------------------------------------------

# BLKCHK -- Compare requested log block position number to current block position
# 	okay = procedure blkchk( wanted, curblk )
# 	
# 	okay		: Found a request block? [return value, (bool)]
# 	wanted		: Block position to search for  [input, (int)]
# 	curblk		: Current block position	[input, (int)]

bool procedure blkchk (wanted, curblk)

int     wanted          	#I Block number to search for
int     curblk          	#I Current block count

bool    okay

begin
	okay = FALSE

	if (wanted == G_INDEF)    # No specific block requested
	    okay = TRUE
	else if (wanted == curblk)
	    okay = TRUE
	else
	    okay = FALSE

	return (okay)
end

#--------------------------------------------------------------------------

# LVLCHK -- Compare an entry's log level to the list of requested levels
# 	okay = procedure lvlchk( gl, line )
# 	
# 	okay		: Retrieve this entry?	[return value, (bool)]
# 	gl		: GL structure		[input, (GL)]
# 	line		: Log entry		[input, (char[])]

bool procedure lvlchk (gl, line)

pointer gl              	#I GL structure (for level preferences)
char    line[ARB]       	#I Line from a GEMLOG file

bool	okay

char	level[LEN_LEVEL_STR]

# Note to developer: If you change the strings below, change them also
#                    in glwrite.x.
#                    Also, note that unlike for glwrite, here there should be
#                    no extra whitespace at the end of the strings, since 
#                    sscan does not preserve the whitespaces.
string	STAT_LEVEL_STR  "STAT"
string	SCI_LEVEL_STR   "SCI"
string	ENG_LEVEL_STR   "ENG"
string	VIS_LEVEL_STR   "VIS"
string	TSK_LEVEL_STR   "TSK"

int	strcmp()

errchk	sscan()

begin
	okay = FALSE

	# Get the entry's log level  (errchk on sscan)
	call sscan (line)
	    call gargwrd (level, LEN_LEVEL_STR)

	# Check level
	if ( strcmp (level, STAT_LEVEL_STR) == 0 ) {
	    if (GL_REQSTAT(gl) == YES)
        	okay = TRUE
	} else if ( strcmp (level, SCI_LEVEL_STR) == 0 ) {
	    if (GL_REQSCI(gl) == YES)
        	okay = TRUE
	} else if ( strcmp (level, ENG_LEVEL_STR) == 0 ) {
	    if (GL_REQENG(gl) == YES)
        	okay = TRUE
	} else if ( strcmp (level, VIS_LEVEL_STR) == 0 ) {
	    if (GL_REQVIS(gl) == YES)
        	okay = TRUE
	} else if ( strcmp (level, TSK_LEVEL_STR) == 0 ) {
	    if (GL_REQTSK(gl) == YES)
        	okay = TRUE
	}

	return (okay)
end

#--------------------------------------------------------------------------

# TASKCHK -- Compare an entry's task tag to the requested task name.
# 	okay = taskchk( sl, line )
# 	
# 	okay		: Retrieve this entry?	    [return value, (bool)]
# 	sl		: SL (selection) structure  [input, (SL)]
# 	line		: Log entry		    [input, (char[])]

bool procedure taskchk (sl, line)

pointer	sl			#I SL (selection) structure
char	line[SZ_LINE]		#I Line from a GEMLOG file

bool	okay

char	taskname[SZ_FNAME], level[LEN_LEVEL_STR]

bool	g_whitespace()
int	strcmp()

errchk	sscan()

begin
	okay = FALSE

	if ( g_whitespace (SL_TSKNAME(sl)) )	#Task selection not requested
	    okay = TRUE
	else {					#Task selection requested
	    #Get task name from line  (errchk on sscan)
	    call sscan (line)
        	call gargwrd (level, LEN_LEVEL_STR)
        	call gargwrd (taskname, SZ_FNAME)

	    if ( strcmp (taskname, SL_TSKNAME(sl)) == 0 )
        	okay = TRUE
	}

	return (okay)
end

#--------------------------------------------------------------------------

# TIMECHK -- Compare the time tag of the log block to the time range requested.
# 	okay = timechk( sl, line )
# 	
# 	okay		: Current block within time range? [return value, (bool)]
# 	sl		: SL (selection) structure	[input, (SL)]
# 	line		: Log entry			[input, (char[])]

bool procedure timechk (sl, line)

pointer	sl			#I SL (selection) structure
char	line[SZ_LINE]		#I BOE line of a GEMLOG file

bool	okay

char	level[LEN_LEVEL_STR], taskname[SZ_FNAME], datestr[SZ_LINE]
int	junk, ltime

int	dtm_ltime()

errchk	sscan()

begin
	okay = FALSE

	if ( (SL_LTIME(sl) == 0) && (SL_UTIME(sl) == 0) ) #no time range request
	    okay = TRUE
	else {
	    # If the upper limit not set, set it to a large number
	    if (SL_UTIME(sl) == 0)
        	SL_UTIME(sl) = MAX_INT

	    #Get time out of 'line'  (errchk on sscan)
	    call sscan (line)
        	call gargwrd (level, LEN_LEVEL_STR)
        	call gargwrd (taskname, SZ_FNAME)
        	call gargwrd (datestr, SZ_LINE)

	    junk = dtm_ltime (datestr, ltime)

	    if ( (ltime >= SL_LTIME(sl)) && (ltime <= SL_UTIME(sl)) )
	       okay = TRUE
	}

	return (okay)
end

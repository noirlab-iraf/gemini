# Copyright(c) 2004-2006 Association of Universities for Research in Astronomy, Inc.

# Author: Kathleen Labrie  3-May-2004

include "glog.h"

.help
.nf
GLR_BLK -- Read a Gemini log file.  Find and retrieve block based on 
           selection criteria.  If verbose is set, print entries to STDOUT.
	status = glr_blk( gl, sl, blkwanted, nlines, blk )
	
	status		: Exit status code		[return value, (int)]
	gl		: GL structure			[input, (GL)]
	sl		: SL (selection) structure	[input, (SL)]
	blkwanted	: Block desired (blk position)	[input, (int)]
	nlines		: Number of lines retrieved	[output, (int)]
	blk[SZ_LINE,ARB]: Block retrieved (array of strings) [output, (char[,])]
.fi
.endhelp

#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
#--------------------------------------------------------------------------

# GLR_BLK -- Read a Gemini log file.  Find and retrieve block based on 
#            selection criteria.  If verbose is set, print entries to STDOUT.
# 	status = glr_blk( gl, sl, blkwanted, nlines, blk )
# 	
# 	status		: Exit status code		[return value, (int)]
# 	gl		: GL structure			[input, (GL)]
# 	sl		: SL (selection) structure	[input, (SL)]
# 	blkwanted	: Block desired (blk position)	[input, (int)]
# 	nlines		: Number of lines retrieved	[output, (int)]
# 	blk[SZ_LINE,ARB]: Block retrieved (array of strings) [output, (char[,])]

int procedure glr_blk( gl, sl, blkwanted, nlines, blk )

pointer	gl		#I GL structure
pointer sl		#I SL (selection) structure
int	blkwanted	#I Block desired
int	nlines		#O Number of lines in the block
char	blk[SZ_LINE,ARB]	#O Block retrieved (array of strings)

int	status

# Other variables
char	level[LEN_LEVEL_STR], msg[SZ_LINE]
int	nboe, nvis
bool	goodblk, okay
pointer	sp, buffer, visbuf
string	VIS_LEVEL_STR	"VIS"

# Gemini functions
bool	taskchk(), blkchk(), lvlchk(), timechk(), g_whitespace()

# IRAF functions
int 	getline(), strcmp(), errget()

errchk	getline(), sscan(), taskchk(), timechk()

begin
	status = 0

	# Allocate stack memory
	call smark (sp)
	call salloc (buffer, SZ_LINE, TY_CHAR)
	call salloc (visbuf, SZ_LINE, TY_CHAR)

	# Initialize
	nlines = 0
	nvis = 0
	goodblk = FALSE
	call strcpy ("", Memc[visbuf], SZ_LINE)

	iferr {
	    while ( getline(GL_FD(gl), Memc[buffer]) != EOF ) {
        	call sscan (Memc[buffer])
		    call gargwrd (level, LEN_LEVEL_STR)
		if ( strcmp ("BOE",level) == 0) {
        	    if (timechk (sl, Memc[buffer]) && 
		        taskchk (sl, Memc[buffer])) {
			
	        	# We have found a block of possible interest.
			SL_BPOS(sl) = SL_BPOS(sl) + 1

			# Is it the right one?
	        	if ( blkchk (blkwanted, SL_BPOS(sl)) )
			    goodblk=TRUE
		    }
		    if (goodblk)
	        	nboe = nboe+1
		}
		else if (strcmp ("EOE",level) == 0) {
		    if (goodblk) {
	        	nboe = nboe - 1
	        	if (nboe == 0) {  
			    # we've reached the end of good block, we're done
			    goodblk = FALSE
			    break
			}
		    }
		}	
		else if (goodblk) {
		
		    # If we are in a child process section of the logs and 
		    # child processes were not requested...
		    
		    if ( (SL_CHILD(sl) == NO) && (nboe > 1) ) {
	        	# And if current line not owned by the requested task...
			# Otherwise, continue: althought it is a child it is 
			# also a requested task.
			
	        	if ( (g_whitespace (SL_TSKNAME(sl))) || 
			     (!taskchk (sl, Memc[buffer])) ) {
			    next		# Skip to next line
			}	    
		    } else if ( (SL_CHILD(sl) == YES) && 
		        (nboe-1 > SL_NCHILD(sl)) ) {
			
	        	#We want just SL_NCHILD(sl) levels of child processes
			#But if child is the requested task, keep it anyway.

			if ( (g_whitespace (SL_TSKNAME(sl))) ||
			     (!taskchk (sl, Memc[buffer])) ) {
			    next		# Skip to next line
			}
		    }
		    
		    # ... done with this child stuff

		    # Line is owned by a valid process.  Check level.
		    if ( lvlchk (gl, Memc[buffer]) ) {
		    
	        	# If it is a VIS entry decide if it is okay to write it.
			# The idea is to avoid two consecutive VIS of the same 
			# type.  (Does not look good.)
			
	        	if ( strcmp (VIS_LEVEL_STR, level) == 0 ) {
			    if ( strcmp (Memc[buffer], Memc[visbuf]) != 0 ) {
				nvis = nvis+1
				call strcpy (Memc[buffer], Memc[visbuf],
				    SZ_LINE)
				okay = TRUE
			    } else {
				okay = FALSE
			    }
			} else {
			    nvis = 0
			    call strcpy ("", Memc[visbuf], SZ_LINE)
			    okay = TRUE
			}
			
			#And avoid more than 3 successive VIS statements
			if (okay && (nvis <= 3)) {
       			    nlines = nlines + 1
       			    call strcpy (Memc[buffer], blk[1,nlines], SZ_LINE)
			    if (GL_VERBOSE(gl) == YES) {
				call printf ("%s")
				    call pargstr (blk[1,nlines])
			    }
			}
		    }
		}
	    }
	} then {
	    # A syserror was caught in getline(), sscan(), taskchk(), 
	    # or timechk().
	    
	    status = errget (msg, SZ_LINE)
	    call sfree (sp)
	    call error (status, msg)
	}

	# Free memory
	call sfree (sp)

	return(status)
end

# Copyright(c) 2004-2009 Association of Universities for Research in Astronomy, Inc.

include <time.h>
include "mefio.h"
include "../../gemtools/pkg/gemlog/glog.h"

# This file contains:
#
#	log_init(..) initialized glog system, setting up globals in a "glogcommon"
#   log_close() closes the glog system
#   me_error(..) wraps error for additional log output or checking
#   me_errreport(..) accepts just error value, and knows associated string
#   task_error(..) doesn't actually throw an error (for tast level functions)
#	log_err(..) logs an error and throws
#	log_warn(..) logs a warning
#	log_info(..) logs informational messages
#
#	i2dp i2dmalloc(xwidth, ywidth)
#	i2dfree(i2dp)
#	i2dset(i2dp, x,y,value)
#	i2dget(i2dp, x,y)
#
#	paryget(pary, index)
#	paryset(pary, index, val)
#	iaryget(iary, index)
#	iaryset(iary, index, val)	

# function below are used by the above (for bad historical
# reasons...should 	be cleaned up)

#	pset( pary, index, val)
#	pget( pary, index)
#	iset( iary, index, val)
#   iget( iary, index)

# ISSUES:
#  the iaryxxx and paryxxx functions are 0 relative, for good reason, but
#  this is not doubt confusing from an IRAF standard practice issue since
#  native SPP arrays are 1 relative. i2dxxx functions are 1 relative.

#
# TODO:: find which of these is called and cull this down... some of it makes
#  TODO:: little sense and is the way it is as these were the first routines 
#  TODO:: I wrote to play with memory.  Most of them are not used!

####### Centralized Error Reporting
# PROCEDUREs me_error(..) and me_errreport(..)
# central error handling

procedure me_errreport(errval)
int errval
# local
char errstr[SZ_LINE]
begin
    switch (errval){
        case MEERR_NOFILE:
            call strcpy(errstr,"MEERR_NOFILE: file not found", SZ_LINE)            

        case MEERR_BADSTRUCT: 
            call strcpy(errstr,"MEERR_BADSTRUCT: structure malformed", SZ_LINE)                   

        case MEERR_OUTOFRANGE:
            call strcpy(errstr, "MEERR_OUTOFRANGE: index out of range", SZ_LINE)                    

        case MEERR_INVAL:
            call strcpy(errstr, "MEERR_INVAL: invalid internal argument", SZ_LINE)                    

        case MEERR_GDRELATE:
            call strcpy(errstr, "MEER_GDRELATE: could not obtain basic frame relation table", SZ_LINE)
        
        default:
            call strcpy(errstr, "MEERR: Internal Error", SZ_LINE)

#       case MEERR_
#           errstr=""
    }
    call me_error(errval, errstr)
end

procedure me_error(errval, errmsg)
int errval
char errmsg[ARB]
begin
    # output to log
    call error(errval, errmsg) 
    
end

procedure log_init(curtask, curpkg)
char curtask[ARB]
char curpkg[ARB]
# LOG RELATED

int success, l_status
bool l_verbose
pointer glpset, sp
pointer op, gl 
pointer tmpstr
bool	clgetb()
pointer clopset()
int	errget(), btoi()
char	l_logfile[SZ_FNAME]

pointer glogopen()

include "glogcommon.h"

begin
	l_status = 0
	success = G_SUCCESS		# G_SUCCESS is defined in 'glog.h'

	#Get task parameter values
	#   Here only the value of 'logfile' needs to be retrieved.
	
	call clgstr ("logfile", l_logfile, SZ_FNAME)
	l_verbose = clgetb ("verbose")

	#Get GLOGPARS pset pointer
	#   The pset values are actually queried for in gloginit().
	#   The pset will be closed by gloginit().  No further manipulation of
	#   the pset pointer is required in the current routine.
	
	glpset = clopset ("glogpars")

	# Allocate memory
	#    The memory for the GL structure is allocated by gloginit()
	#    The memory for the OP structure must be allocated here since 
	#    we will need it to send options to gloginit().
	
	call smark (sp)
	call salloc (tmpstr, SZ_LINE, TY_CHAR)
	call opalloc (op)

#   GLOGPARS and in the OP structure.  The log file is open, or 
	#   created, and the file descriptor is saved in the GL structure.
	#   Essential log entries are written to the log file, as well as the
	#   list of parameters and their values, 'paramstr'.
	
	OP_FL_APPEND(op) = YES
    OP_FORCE_APPEND(op) = NO
	OP_VERBOSE(op) = btoi (l_verbose)
    OP_ERRNO(op) = 1
	gl = NULL
	iferr (gl = glogopen (l_logfile, curtask, curpkg, glpset, op)){
	    gl = NULL		# memory already been freed in glogopen
	    l_status = errget (Memc[tmpstr], SZ_LINE)
	    call printf ("ERR: Can't open logfile:%d %s.\n")
		call pargi (l_status)
		call pargstr (Memc[tmpstr])
	    call opfree (op)
	    call sfree (sp)
	    call clputi ("status", l_status)
	    return
	} else {		# This block is optional
	    if (OP_STATUS(op) != 0) {
 		l_status = l_status + errget (Memc[tmpstr], SZ_LINE)
 		call printf ("GLOGOPEN WARNING: %d %s.\n")
 		    call pargi (OP_STATUS(op))
 		    call pargstr (Memc[tmpstr])
 		OP_STATUS(op) = 0
	    }
	}
	call flush (GL_FD(gl))
    
    # set glog common globals
    g_gl = gl
    call strcpy (l_logfile, g_logfile, SZ_FNAME)
    g_op = op
end

procedure log_close()
include "glogcommon.h"
begin
    call gl_close(g_gl)
	call opfree (g_op)
end

procedure task_error(errval, errmsg)
int errval
char errmsg[ARB]
# -- locals
char outmsg[SZ_LINE]
int l_status
int glogprint()
# -- glog globals
include "glogcommon.h"
begin
#    call printf("TASK_ERR: %d, %s\n")
#    call pargi(errval)
#    call pargstr(errmsg)
#    call flush(STDOUT)
# we don't want tasks to create errors on the top level so to allow CL to recover
    call sprintf(outmsg, SZ_LINE, "(%d), %s")
    call pargi(errval)
    call pargstr(errmsg)
#    call flush(STDOUT)
    l_status = glogprint (g_gl, ENG_LEVEL, G_ERR_LOG, outmsg, g_op)
    call clputi("status", 1) #DOC: gemini tasks MUST have status parameters
end

procedure log_info(logmsg)
char logmsg[ARB]
# get global glog parameters
int l_status

#functions
int glogprint()

include "glogcommon.h"
begin
    l_status = glogprint (g_gl, ENG_LEVEL, G_STR_LOG, logmsg, g_op)
end

procedure log_warn(logmsg)
char logmsg[ARB]
# get global glog parameters
int l_status

#functions
int glogprint()

include "glogcommon.h"
begin
    l_status = glogprint (g_gl, STAT_LEVEL, G_WARN_LOG, logmsg, g_op)
end

procedure log_err(errnum, logmsg)
int errnum
char logmsg[ARB]

# get global glog parameters
include "glogcommon.h"
begin
    call task_error(errnum, logmsg)
end

#### I2D PROCEDURES (see mefio.h)
#
#PROCEDURE i2d
pointer procedure i2dmalloc(xwidth, ywidth)
int xwidth, ywidth
# -- locals
pointer retval
begin
	call calloc(retval, LEN_I2D + ((xwidth+1) * (ywidth+1)), TY_INT)
	I2D_XW(retval) =  xwidth
	I2D_YW(retval) =  ywidth
    I2D_REALNROWS(retval) = ywidth
	return retval
end

procedure i2dfree(i2dp)
pointer i2dp
begin
	call mfree(i2dp, TY_STRUCT)
	i2dp = 0
end

procedure i2dset(i2dp, x, y, val)
pointer i2dp
int x,y,val
#--- local
int lx,ly
int index
begin
    # note: XW and WY refer to the body of the table
    # there is also allocated a hidden column and row
    # used for type labeling the table/array
    # the default behavior is to ignore these rows
    # because they contain meta-information...
    # therefore to index the correct data one dimensionally
    # involves adding to the width (I2D_XW(i2dp)+1)
	lx = x
	ly = y
	index = LEN_I2D + (ly * (I2D_XW(i2dp)+1)) + lx;
	call iset(i2dp, index+1, val)
end

int procedure i2dget(i2dp, x, y)
pointer i2dp
int x,y
#-- locals
int lx,ly
int index
#-- procs
int iget()
begin
    # note: there is an extra row for metadata (one int for rowtype)
    #  at the 0 index.
	lx = x
	ly = y
	index = LEN_I2D + (ly * (I2D_XW(i2dp)+1)) + lx;
	
	return iget(i2dp, index+1)
end

procedure i2dprint(i2dp)
pointer i2dp
#--locals
int i, j
int i2dget()
begin
	call printf("i2dprint(): array struct at 0x%x I2D_XW=%d I2D_YW=%d\n")
	call pargi(i2dp)
	call pargi(I2D_XW(i2dp))
	call pargi(I2D_YW(i2dp))
	call flush(STDOUT)
	
	for (i = 0; i <= I2D_YW(i2dp); i = i+1)
	{
		for (j = 0 ; j <= I2D_XW(i2dp); j = j+1)
		{
			call printf("%d ")
			#call pargi(7)
			call pargi(i2dget(i2dp, j, i))
			call flush(STDOUT)
		}
		call printf("\n")
		call flush(STDOUT)
	}
	call printf("\n")
	call flush(STDOUT)
end
#
#### END OF I2D PROCEDURES

#PROCEDURE parymalloc()
# dynamically allocated pointer array
pointer procedure parymalloc(siz)
int siz
#-- locals
pointer retval
begin
	call calloc(retval, siz, TY_STRUCT)
	return retval
end


# PROCEDURE iset()
# set member of dynamically allocated pointer array
procedure iset( iary, index, val)
pointer iary
int index
int val

begin
	Memi[iary+index-1]=val
	return
end


#PROCEDURE iget()
# get member of dynamically allocation pointer array
int procedure iget( iary, index)
pointer iary
int index

int reti
begin
	reti = Memi[iary + index - 1]
	return reti
end


#PROCEDURE paryset()
procedure paryset( pary, index, val)
pointer pary
int index
pointer val
#--- locals below
begin
	call pset(Memi[pary],index,val)
	return;
end

#PROCEDURE paryget()
# get member of dynamically allocation pointer array
pointer procedure paryget(pary, index)
pointer pary
int index
#----
pointer pget()
begin
	return pget(Memi[pary], index)
#NOTE: the code below is NOT equiv to the above, and does not work
#      why?
#	ary = Memi[pary]
#	return pget(ary, index)	
end

# START OF REDUNDANT PARYGET PROCS (THESE ARE NOT USED IN MEFIO...)
# NOTE::!!!!!!!!!!!!!!!!!!!!!!
# the extra paryget functions were/are for testing (see gtest task)
# and remain now because I may need to use them again although
# I've solved the original mystery as to why these versions are
# not equivalent in SPP.

#PROCEDURE paryget2()
# get member of dynamically allocation pointer array
pointer procedure paryget2(pary, index)
pointer pary
int index
#----
int ary
pointer pget()
begin
#	return pget(Memi[pary], index)
#NOTE: the code below is NOT equiv to the above, and does not work
#      why?
	ary = Memi[pary]
	return pget(ary, index)	
end
#PROCEDURE paryget3()
# get member of dynamically allocation pointer array
pointer procedure paryget3(pary, index)
pointer pary
int index
#----
pointer ary
pointer pget()
begin
#	return pget(Memi[pary], index)
#NOTE: the code below is NOT equiv to the above, and does not work
#      why?
	ary = Memi[pary]
	return pget(ary , index)	
end

# !!!!!!! see note above !!!!!!
# END OF EXTRA (TEST) PARYGET CALLS


#PROCEDURE iaryget()
# get member of dynamically allocation int array
int procedure iaryget(iary, index)
pointer iary
int index
#----
pointer iget()
begin
	return iget(Memi[iary], index)
#NOTE: see paryget notes
end

#PROCEDURE iaryset()
procedure iaryset( iary, index, val)
pointer iary
int index
int val
#--- locals below
begin
	call iset(Memi[iary],index,val)
	return
end


#
#
#
#
### THOSE BELOW THIS LINE ARE NOT USED (I THINK, TODO:: check and clean up)
#
#
#
#


# PROCEDURE pset()
# set member of dynamically allocated pointer array
procedure pset( pary, index, val)

pointer pary[ARB]
int index
pointer val

begin
	pary[index] = val
	return
end


#PROCEDURE pget(..)
# get member of dynamically allocation pointer array
pointer procedure pget( pary, index)

pointer pary[ARB]
int index

begin
	return pary[index]
end

#PROCEDURE memdump(..)

procedure memdump( pvoid, bytesize)
pointer pvoid
int		bytesize

int i 

begin
	
	call eprintf("pvoid = 0x%x, dumplen=%d")
	call pargi(pvoid)
	call pargi(bytesize)
	for (i = 0; i < bytesize; i = i+1)
	{
		call eprintf("Memc[P2C(0x%x + %2d)] = %10d (0x%10x) \n")
		call pargi(pvoid)
		call pargi(i)
		call pargi(Memc[P2C(pvoid+i)]);
		call pargi(Memc[P2C(pvoid+i)]);
	}	
	
end

procedure printop(op)
pointer op

begin

    call printf("OP_FL_APPEND = %d\n")
    call pargi(OP_FL_APPEND(op))
    call printf("OP_FORCE_APPEND = %d\n")
    call pargi(OP_FORCE_APPEND(op))
    call printf("OP_DEFLOG = %d\n")
    call pargi(OP_DEFLOG(op))
    call printf("OP_VERBOSE = %d\n")
    call pargi(OP_VERBOSE(op))
    call printf("OP_STATUS = %d\n")
    call pargi(OP_STATUS(op))
    call printf("OP_VISTYPE = %d\n")
    call pargi(OP_VISTYPE(op))
    call printf("OP_ERRNO = %d\n")
    call pargi(OP_ERRNO(op))
    call printf("OP_FORK = %d\n")
    call pargi(OP_FORK(op))
    call printf("OP_CHILD_P = %d\n")
    call pargi(OP_CHILD_P(op))
    
    call flush(STDOUT)
    
end

bool procedure mefexists(mefname, fthrow)
char mefname[ARB] # MEF name to check... IMIO style, does not need file extension
bool fthrow		  # flag which if true causes error to be thrown if not found#
# ---------
char tmpname[SZ_FNAME]
pointer tmpimp
bool ldebug
int strlen()
pointer immap()
begin
	
	ldebug = false
	
	if (ldebug) {
		call printf("existsmef(..)\n")
		call flush(STDOUT)
	}
	# check 
	if (strlen(mefname) > (SZ_FNAME - 4)) {
		call log_info(mefname);
		call log_err(MEERR_BUFFERERR, "mef name TOO LONG (existmef BUFFERERR)");		
	}
	if (ldebug) {
		call printf("existsmef(%s):bufside\n")
		call pargstr(mefname)
		call flush(STDOUT)
	}
	
	call strcpy(mefname, tmpname, SZ_FNAME)
	call strcat("[0]", tmpname, SZ_FNAME)
	
    if (ldebug) {
        call printf("testname = %s\n")
        call pargstr(tmpname)
        call flush(STDOUT)
    }
	iferr (	tmpimp = immap(tmpname, READ_ONLY, NULL) )	{
		# doesn't exist
		if (fthrow) {
		# throw based on flag
			call sprintf(tmpname, SZ_FNAME, "existsmef(): %s doesn't exist according to IMIO (not a MEF?)")
			call pargstr(mefname)
			call log_err(MEERR_NOTEXIST, tmpname);
		}
		return false
	} else {
		call imunmap(tmpimp)
		#call imflush()
		return true
	}
end

bool procedure encmp(key, name)
char key[ARB]
char name[ARB]
#---
bool streq()
char upkey[SZ_LINE]
char lwkey[SZ_LINE]
begin
    call strcpy(key,upkey, SZ_LINE)
    call strcpy(key,lwkey, SZ_LINE)
    call strupr(upkey)
    call strlwr(lwkey)
    
    if (streq(upkey,name) || streq(lwkey,name)) {
        return true;
    }
    else
    {
        return false;
    }
end

procedure mtstamp(mefn, taskn)
char mefn[ARB]  # mef file name to stamp (PHU)
char taskn[ARB] # taskname to use for timestamp name

# --- local vars
pointer mep
long tmpl,tmpu #for local and UT time respectively
pointer pphu
char tmpstr[SZ_LINE]
char aboutstr[SZ_LINE]
int tm[LEN_TMSTRUCT]
bool ldebug
# --- functions
long clktime(),lsttogmt()
pointer megphu()
pointer memap()

begin
    
    ldebug = false
    
    # NOTE: TLM gets in via mimcopy, so I'm not stamping TLM-GEMINI from here
    mep = memap(mefn)
    
    if (mep != 0){
        tmpl = 0
        while (tmpl < 400000) {
            tmpl = clktime(tmpl)
            tmpu = lsttogmt(tmpl)
#call cnvtime(tmpl, tmpstr, SZ_LINE)
            call brktime(tmpu, tm)
            call sprintf(tmpstr,SZ_LINE, "%4d-%02d-%02dT%02d:%02d:%02d")
            call pargi(TM_YEAR(tm))
            call pargi(TM_MONTH(tm))
            call pargi(TM_MDAY(tm))
            call pargi(TM_HOUR(tm))
            call pargi(TM_MIN(tm))
            call pargi(TM_SEC(tm))
            if (ldebug) {
                call printf("-*-%s-*-\n")
                call pargstr(tmpstr)
                call flush(STDOUT)
            }
        }
        
        pphu = megphu(mep)
        #call imastr(im2, "GEM-TLM", tmpstr) 
        call sprintf(aboutstr,SZ_LINE, "UT Time Stamp for %s")
        call pargstr(taskn)
        
        if (pphu != 0) {
            # this is the task time stamp
            call imastrc(pphu, taskn, tmpstr, aboutstr)
            # this is the GEM-TLM stamp
            call imastrc(pphu, "GEM-TLM", tmpstr, "UT Last modification with GEMINI")
        }
        
        call meunmap(mep)
    }
        
end

# setNEXTEND() will set the NEXTEND header keyword to the actual number of 
#  extensions in the file, which is good because CL has no way to easily get
#  at the actual number, so mis-set NEXTEND values are likely to screw up 
#  loops in cl scripts.

procedure setNEXTEND(mefn)
char mefn[ARB]  # mef file name to stamp (PHU)

pointer mep,pphu
int nextend

pointer memap(), megphu()
begin

# map the mef
mep = memap (mefn)    
    
# get actual number of extensions from struct (includes PHU)
# note: NEXTEND is not to include PHU
nextend = ME_NUMEXTS(mep) - 1

# set the NEXTEND header in the PHU
pphu = megphu(mep)
if (pphu != 0) {
        call imakic(pphu, "NEXTEND", nextend, "Number of extensions" )
    }

# unmap the mef
    call meunmap(mep)
    
end

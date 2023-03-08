# Copyright(c) 2004-2005 Association of Universities for Research in Astronomy, Inc.

#
# not hard coded any more define MDFEXTNAME	"MINIMAL.TAB"
#

# This file contains:
#         gdpropagate(...) - 
#
# Support routines:
#         gdpropagate(...) - Short description
#
# See also:
#         gdpropagate(...) - Short description

#
# ... define statements ...
#

# GDPROPAGATE --
#	This routine checks to see if key files have been propagated based
# 	on Gemini MEF Data Type.  (WARNING: Data Type Currently
# 	Ignored... propagates 	MDF, aka MINIMAL.TAB extensions)
# NOTE:: TODO:: I want to put PHU propagation in here!
include "mefio.h"
include <error.h>


procedure gdpropagate (source, dest, tname)
char source[ARB]	# source MEF
char dest[ARB]		# dest MEF
char tname[ARB] 	# type of propagation (default="AUTO" PHU/MDF) 
# locals below #
char mdf_ext[SZ_LINE]
pointer smep	# for mefio memap
pointer smef, dmef # for fitsutil manipulations
int	mdfind		# mdf index (aka mdf "group" number)
bool ldebug
pointer mef_open()
int meindbnv()
pointer memap()
bool strneq()
char tmpstr[SZ_LINE]
int errc
bool verbose,ftmp
bool clgetb()
int errget()
pointer immap()
pointer tmpimp
char tmpimname[SZ_FNAME]
char tmpimname2[SZ_FNAME]
bool mefexists()
bool f_mefis
errchk immap(), imunmap(), memap(),meindvnv(),mef_open(), mef_copy_extn(), mef_close(), meunmap() 

begin
	
	ldebug = false # true

    if (ldebug) {
        call printf("gdpropagate(%s,%s,%s)\n")
        call pargstr(source)
        call pargstr(dest)
        call pargstr(tname)
        call flush(STDOUT)
    }
    

	# MDF propagation... done by default
    call clgstr("mdf_ext", mdf_ext, SZ_LINE)
    verbose= clgetb("verbose")
# memap
    
    iferr(smep = memap(source)){
        errc = errget(tmpstr, SZ_LINE)
        call task_err(errc, tmpstr)
        call sprintf(tmpstr, SZ_LINE, "gdpropagate(): cannot copy MDF,\"%s\"\nIs this is a correct reference image (not a constant)")
        call pargstr(source)
        call log_err(MEERR_MEMAPFAIL, tmpstr)
        call erract(EA_ERROR)
    }
   
    if (ldebug){
        call printf("smep = %x\n")
        call pargi(smep)
        call flush(STDOUT)
        call meprint(smep)
    }
# getindex by name and ver (get just one (NOTE::is that ok?))
	
    if (ldebug) {
	    call printf("before meindbnv\n")
	    call flush(STDOUT)
    }
	
	mdfind = meindbnv(smep, mdf_ext, -1)
	call meunmap(smep)

    if (ldebug) {
        call printf("after - meindvnv: mdfind=%d\n")
	    call pargi(mdfind)
        call flush(STDOUT)
    }

	if (mdfind > 0)
	{
		if (verbose) {
			call sprintf(tmpstr, SZ_LINE, "Going to copy MDF from %s[%d]")
			call pargstr(source)
			call pargi(mdfind)
			if (ldebug){
				call printf(tmpstr)
				call flush(STDOUT)
			}
            call log_info(tmpstr)
		}
	}
	else
	{
        if (verbose){
		    call log_info("No MDF present in source");
			if (ldebug){
				call printf("debug: No MDF present in source\n")
				call flush(STDOUT)
			}
        }
	}

	# use fitsutil open dest
	f_mefis = mefexists(dest,false)
    
	if(ldebug) {
		call printf("smef=%x\ndmef=%x, f_exist=%b\n")
		call pargi(smef)
		call pargi(dmef)
        call pargb(f_mefis)
		call flush(STDOUT)
	}
    # use fitsutil to copy PHU
	# PHU propagation...

#if the file doesn't exist
#	ftmp = mefexists(source, false)
	if (!f_mefis) {
		call sprintf(tmpimname, SZ_FNAME, "%s[0]")
		call pargstr(source)
		call sprintf(tmpimname2, SZ_FNAME, "%s[APPEND]")
        call pargstr(dest)
        call mimcopy(tmpimname, tmpimname2)
#		if (ldebug) call printf("does not exist\n")		
#		call mef_copy_extn(smef, dmef, 0)
#		if (ldebug) call printf("copied extn\n")				
	}
	else {
		if (ldebug) call printf("prop destination does exist\n")
	}

	call flush(STDOUT)

    # use fitsutil to copy source to dest
	if (mdfind > 0) {
	    smef = mef_open(source, READ_ONLY, 0)
	    dmef = mef_open(dest, APPEND, 0)
	    
        call mef_copy_extn (smef, dmef, mdfind)
       	# use fitsutil open source

	    # use fitsutil close dest
	    call mef_close(dmef)
	
        # use fitsutil close source
	    call mef_close(smef)
    }

end

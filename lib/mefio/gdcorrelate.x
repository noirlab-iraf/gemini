# Copyright(c) 2004-2005 Association of Universities for Research in Astronomy, Inc.

#
include "mefio.h"
#

# This file contains:
#
#         gdrelate(...) - Short description
#
# Support routines:
#         gdrelate(...) - Short description
#
# See also:
#         gdrelate(...) - Short description

#
# ... define statements ...
#

# GDRELATE -- Description

pointer procedure gdrelate (mefpary,nomefs, name)
pointer mefpary	# array of mefpointers
int nomefs		# number of mefs to relate
char name[ARB] 	# extnam (e.g. SCI)
#locals
pointer corary
int rowi,coli
int invalidrow
pointer fmrelate()
int i2dget()

bool ldebug

char tmpstr[SZ_LINE]
int  tmpi
int  errget()

errchk me_errreport()
begin
    ldebug = false
    
    if(ldebug) { 
        call printf("gdrelate(%x, %d, %s)\n")
        call pargi(mefpary)
        call pargi(nomefs)
        call pargstr(name)
        call flush(STDOUT)
    }
#
    
    iferr (	corary = fmrelate(mefpary, nomefs, name) ) {
        tmpi = errget(tmpstr, SZ_LINE)
        call error(tmpi, tmpstr)
    }

	invalidrow = I2D_ALLVALID
    # mark rows with missing frames (element will have 0
	for (rowi = 1; rowi <= I2D_NROWS(corary); rowi=rowi+1)
	{
		for (coli =1 ; coli <= I2D_NCOLS(corary); coli=coli+1)
		{
            if (ldebug) {
                call printf("(%d,%d)=%d")
                call pargi(coli)
                call pargi(rowi)
                call pargi(i2dget(corary,coli,rowi))
                call flush(STDOUT)
            }
            # note: I2DROWVALID == 0, i2d array should be initted to 0s
			if (i2dget(corary, coli,rowi) == 0) {
                if (ldebug) {
                    call printf("| wrote INVAL |")
                    call flush(STDOUT)
                }
                call i2dset(corary, 0, rowi, I2DROWINVAL)
			}
            
            if (ldebug){
                call printf("  (%d,%d)=%d\n")
                call pargi(coli)
                call pargi(rowi)
                call pargi(i2dget(corary,coli,rowi))
                call flush(STDOUT)
            }

		}
	}
    
    #call i2dprint(corary)
	
# now we have the raw corary, it may have 0's where related frames
# were not found... we need to prune out such entries.
# INFO:: pruning rows in corary with 0s is a matter of starting at
# -INFO:: extver row #1 and truncating the list at the first row with a 
# -INFO:: zero in it.
	return corary
end

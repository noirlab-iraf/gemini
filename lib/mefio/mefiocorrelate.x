# Copyright(c) 2004-2005 Association of Universities for Research in Astronomy, Inc.

# ... include statements ...
include <mefio.h>
# This file contains:
#
#         fmrelate(...) - Frame relation function

# ROUTINE_NAME -- Description

pointer procedure fmrelate(mefpary,nomefs, name)
pointer mefpary	# array of mefpointers
int nomefs		# number of mefs to relate
char name[ARB] 	# extnam (e.g. SCI)

# -- local
int i,j,k
int noexts0
pointer mep0, mep
pointer corary
int paryget()
pointer i2dmalloc()
int meindbnv()
int	index
char tmpstr[SZ_LINE]
bool locdebug
bool nonzeroext	# flag for checking if at least one non-zero extension exists
int indbnv
# functions
int memaxver()
int maxv
bool streq()

include <mefiocommon.h>

begin
	locdebug = false
    if (locdebug) {
        call printf("fmrelate(): entered\nmefpary=%x\nnomefs=%d\nextname=%s\n")
        call pargi(mefpary)
        call pargi(nomefs)
        call pargstr(name)
        call flush(STDOUT)
    }
	corary = 0
	noexts0=1;
	for (i = 1; i <= nomefs; i= i+1)
	{
		mep0 = paryget(mefpary, i-1)
        if (locdebug) {
            call printf("fmrelate(): paryget[%d]=%x\n")
            call pargi(i-1)
            call pargi(mep0)
            call flush(STDOUT)
        }
          
        if (mep0 != NULL) {
            if (locdebug){
                indbnv = meindbnv(mep0, name, i)
                call printf("fmrelate(): memaxver(%x,%s)=%d <%d>\n")
                call pargi(mep0)                
                call pargstr(name)
                call pargi(maxv)
                call pargi(indbnv)
                call flush(STDOUT)                
            }
            maxv = memaxver(mep0, name)
		    noexts0 = max(noexts0, maxv)
        }
	}
	
	if (locdebug)
	{	
		call printf("fmrelate(): number of MEFs = %d\n")
		call pargi(nomefs)

		call printf("fmrelate(): num %s exts in most populated mef = %d\n")
		call pargstr(name)
		call pargi(noexts0)
		call flush(STDOUT)
	}
	# allocate the 2d dynamic int array using i2dxxx
	corary = i2dmalloc(nomefs, noexts0)
	
	if (locdebug)
	{
		call printf("fmrelate(): allocated correlation array\n")
		call flush(STDOUT)
	}
	
	# loop over all the mefs
	for (i = 1; i <= nomefs; i= i+1)
	{
        nonzeroext=false
		mep = paryget(mefpary, i-1)
		for (j = 1; j <= noexts0; j=j+1) {
			if (mep == NULL) {
                # if the mep is null, the correlation frame is FM_SCALEROPEAND
                #NOTE: why support mep of NULL?  because gemexpr, 
				# for example, has numeric operands sometimes, so things that
				# match up it's easier to give these a column of special values than 
				# not support them.
				index = FM_SCALEROPERAND
			}
			else
			{
				if (locdebug)
				{
					call printf("fmrelate():calling meindbnv(0x%x,%s,%d)\n")
					call pargi(mep)
					call pargstr(name)
					call pargi(j)
					call flush(STDOUT)
				}
                
                #note: the conditional tests if this operand is a special optype
                #      which means to use a special extname other than that handed in
                if ((ME_OPTYPE(mep) == OPT_NORMAL) || (ME_OPTYPE(mep) == OPT_NULL)) {
                    if (locdebug){
                        call printf("fmrelate(): OPTYPE == OPT_NORMAL\n")
                        call flush(STDOUT)
                    }
				    index = meindbnv(mep, name, j)
                    if ((index == 0) && (j == 1)) { 
                        # ext not found, EXTVER == 1, then check for unnumbered extension
                        index = meindbnv(mep, name, INDEFI)
                        if (index != 0) {
                            call sprintf(tmpstr, SZ_LINE, "Inferred EXTVER=1 for %s[%d] (EXTNAME=%s)")
                            call pargstr(ME_FILENAME(mep))
                            call pargi(index)
                            call pargstr(name)
                            call log_warn(tmpstr)
                        }                       
                    }
                } else {
                    index = meindbnv(mep, opt_sci_ext, j)
                }
			}
			call i2dset(corary, i, j, index)
            if (index != 0) nonzeroext=true
			
			if (locdebug)
			{
				call printf("fmrelate():corary[%d,%d] = %d\n")
				call pargi(i)
				call pargi(j)
				call pargi(index)
				call flush(STDOUT)
			}
		}
        # handle special case... no extensions matched, then we do allow
# mefs with unnamed and unnumbered extensions if we are looking for SCI
        if ((mep!= NULL) && (nonzeroext == false) && streq(opt_sci_ext, name)) {
            for (j = 1; j <= noexts0; j=j+1) {
                call sprintf(tmpstr, SZ_LINE, "Inferring EXTNAME=(%s) and EXTVERs for %s")
                call pargstr(opt_sci_ext)
                call pargstr(ME_FILENAME(mep))
                call log_warn(tmpstr)
                if (locdebug) {
                    call printf("extname=-%s-%d\n")
                    call pargstr(Memc[paryget(ME_NAMES(mep),j)])
                    call pargi(paryget(ME_NAMES(mep),j))
                    call flush(STDOUT)
                }
             	if (streq(Memc[paryget(ME_NAMES(mep),j)], "")) {
				    call i2dset(corary, i,j,j)
                } else {
                    # if one of these extensions IS named... zero them all, not a legal inference
                    for (k = 1; k <= noexts0; k = k+1) {
                        call i2dset(corary, i,k,0)       
                    }
                    call error(MEERR_BADOPERAND, "Inference FAILED-some extensions have EXTNAME")
                    break
                }
            }
        }
	}
	
    if (locdebug) {
        call i2dprint(corary)
    }
	return corary
end

# this just call i2dfree, but it's good to call it to allow any special
# fmfree related cleanup if it's nec in the future
procedure fmfree(corary)
pointer corary
begin
	call i2dfree(corary)
	corary = 0
end
	
#TODO:: fmaryfree(), i2dfree()

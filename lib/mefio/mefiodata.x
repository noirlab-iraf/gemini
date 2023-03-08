# Copyright(c) 2004-2005 Association of Universities for Research in Astronomy, Inc.

include <mefio.h>

# mefiodata.gx
#  (this file is a generic file (you may be reading this in an instance, eg. mefiodatai.x))
#
# These routines thinly wrap IMIO or TABLES calls, currently with
#  single line thunks to the libraries in question.  The purpose
#  of these routines however is to act as a mediation for access
#  to data which is in the same routine.  So future work (TODO::)
#  involves mediating file access contention (e.g. simultaneous
#  get by next line and put by next line)... such access is not
#  going to be simultaneous in a non-threaded system like IRAF
#  and yet there are still possible issues since reads or writes
#  to one extension may conflict with reads or writes to another
#  extension if it is in the same MEF but not otherwise.
#
# This layer exists now in order to allow development against the
#  final interface before contention management is implimented

# non generic data access

#
pointer procedure imgnlb(ep, buffer, v)
pointer ep,buffer
long v[ARB]

#
# ... variables declaration ...
#
pointer imgnli()
begin
	return imgnli(EXT_EXTP(ep), buffer, v)
end

pointer procedure impnlb(ep, buffer, v)
pointer ep,buffer
long v[ARB]

#
# ... variables declaration ...
#
pointer impnli()
begin
	return impnli(EXT_EXTP(ep), buffer, v)
end

int procedure megftype(mextp, parmname)
pointer mextp
char parmname[ARB]
# ---
char extname[SZ_FNAME]
pointer immap()
int imgftype()
begin
    # note: I'd like to have a function for this lazy loading
    if (ME_PHU(mextp)==0) {
		call sprintf(extname, SZ_LINE, "%s[0]")
		call pargstr(Memc[ME_PFILENAME(mextp)])
        ME_PHU(mextp)= immap(extname, READ_ONLY, 0)
	}		
    return imgftype(ME_PHU(mextp), parmname)
end


# GENERIC FUNCTIONS, expanded in-file
# For Most functions, cover these Types
#

bool procedure megetb(mextp, parmname) 
pointer mextp
char parmname[ARB]
#---
char extname[SZ_FNAME]
bool imgetb()
pointer immap()
begin
    # note: I'd like to have a function for this lazy loading
	if (ME_PHU(mextp)==0) {
		call sprintf(extname, SZ_LINE, "%s[0]")
		call pargstr(Memc[ME_PFILENAME(mextp)])
        ME_PHU(mextp)= immap(extname, READ_ONLY, 0)
	}		
    
    return imgetb(ME_PHU(mextp),parmname)
end

short procedure megets(mextp, parmname) 
pointer mextp
char parmname[ARB]
#---
char extname[SZ_FNAME]
short imgets()
pointer immap()
begin
    # note: I'd like to have a function for this lazy loading
	if (ME_PHU(mextp)==0) {
		call sprintf(extname, SZ_LINE, "%s[0]")
		call pargstr(Memc[ME_PFILENAME(mextp)])
        ME_PHU(mextp)= immap(extname, READ_ONLY, 0)
	}		
    
    return imgets(ME_PHU(mextp),parmname)
end

int procedure megeti(mextp, parmname) 
pointer mextp
char parmname[ARB]
#---
char extname[SZ_FNAME]
int imgeti()
pointer immap()
begin
    # note: I'd like to have a function for this lazy loading
	if (ME_PHU(mextp)==0) {
		call sprintf(extname, SZ_LINE, "%s[0]")
		call pargstr(Memc[ME_PFILENAME(mextp)])
        ME_PHU(mextp)= immap(extname, READ_ONLY, 0)
	}		
    
    return imgeti(ME_PHU(mextp),parmname)
end

long procedure megetl(mextp, parmname) 
pointer mextp
char parmname[ARB]
#---
char extname[SZ_FNAME]
long imgetl()
pointer immap()
begin
    # note: I'd like to have a function for this lazy loading
	if (ME_PHU(mextp)==0) {
		call sprintf(extname, SZ_LINE, "%s[0]")
		call pargstr(Memc[ME_PFILENAME(mextp)])
        ME_PHU(mextp)= immap(extname, READ_ONLY, 0)
	}		
    
    return imgetl(ME_PHU(mextp),parmname)
end

real procedure megetr(mextp, parmname) 
pointer mextp
char parmname[ARB]
#---
char extname[SZ_FNAME]
real imgetr()
pointer immap()
begin
    # note: I'd like to have a function for this lazy loading
	if (ME_PHU(mextp)==0) {
		call sprintf(extname, SZ_LINE, "%s[0]")
		call pargstr(Memc[ME_PFILENAME(mextp)])
        ME_PHU(mextp)= immap(extname, READ_ONLY, 0)
	}		
    
    return imgetr(ME_PHU(mextp),parmname)
end

double procedure megetd(mextp, parmname) 
pointer mextp
char parmname[ARB]
#---
char extname[SZ_FNAME]
double imgetd()
pointer immap()
begin
    # note: I'd like to have a function for this lazy loading
	if (ME_PHU(mextp)==0) {
		call sprintf(extname, SZ_LINE, "%s[0]")
		call pargstr(Memc[ME_PFILENAME(mextp)])
        ME_PHU(mextp)= immap(extname, READ_ONLY, 0)
	}		
    
    return imgetd(ME_PHU(mextp),parmname)
end





pointer procedure mignls(ep, buffer, v)
pointer ep,buffer
long v[ARB]

#
# ... variables declaration ...
#
pointer imgnls()
bool ldebug
begin
	ldebug = false
	 
		if (ldebug)
		{
			call printf("at mignl$t()\n")
			call flush(STDOUT)
		}

	
	return imgnls(EXT_EXTP(ep), buffer, v)

end

pointer procedure mipnls(ep, buffer, v)
pointer ep,buffer
long v[ARB]

#
# ... variables declaration ...
#
pointer impnls()
bool ldebug
begin
	ldebug = false
	 
		if (ldebug)
		{
			call printf("in mipnl$t()\n")
			call flush(STDOUT)
		}

	
	return impnls(EXT_EXTP(ep), buffer, v)

end

# migl1$t -- get line from a 1D image

pointer procedure migl1s(ep)
pointer ep
#
# ... variables declaration ...
#
pointer imgl1s()

begin
	
	return imgl1s(EXT_EXTP(ep))

end


# miplt$t -- put line into a 1D image

pointer procedure mipl1s(ep)
pointer ep
#
# ... variables declaration ...
#
pointer impl1s()

begin
	
	return impl1s(EXT_EXTP(ep))

end

# migl2$t -- get line from a 2D image

pointer procedure migl2s(ep, line)
pointer ep
int line
#
# ... variables declaration ...
#
pointer imgl2s()

begin
	
	return imgl2s(EXT_EXTP(ep), line)

end

# mipl2$t -- put line into a 2D image

pointer procedure mipl2s(ep, line)
pointer ep
int line
#
# ... variables declaration ...
#
pointer impl2s()

begin
	
	return impl2s(EXT_EXTP(ep), line)

end

# migl3$t -- get line from a 3D image

pointer procedure migl3s(ep, line, band)
pointer ep
int line, band
#
# ... variables declaration ...
#
pointer imgl3s()

begin
	
	return imgl3s(EXT_EXTP(ep), line, band)

end


# mipl3$t -- put line into a 3D image

pointer procedure mipl3s(ep, line, band)
pointer ep
int line, band
#
# ... variables declaration ...
#
pointer impl3s()

begin

	return impl3s(EXT_EXTP(ep), line, band)

end

# NOTES:
# TODO:: these are not all tested... a related TODO:: note exists in the 
# gtest task t_testmefio.x code (that's where testing code should go)





pointer procedure mignli(ep, buffer, v)
pointer ep,buffer
long v[ARB]

#
# ... variables declaration ...
#
pointer imgnli()
bool ldebug
begin
	ldebug = false
	 
		if (ldebug)
		{
			call printf("at mignl$t()\n")
			call flush(STDOUT)
		}

	
	return imgnli(EXT_EXTP(ep), buffer, v)

end

pointer procedure mipnli(ep, buffer, v)
pointer ep,buffer
long v[ARB]

#
# ... variables declaration ...
#
pointer impnli()
bool ldebug
begin
	ldebug = false
	 
		if (ldebug)
		{
			call printf("in mipnl$t()\n")
			call flush(STDOUT)
		}

	
	return impnli(EXT_EXTP(ep), buffer, v)

end

# migl1$t -- get line from a 1D image

pointer procedure migl1i(ep)
pointer ep
#
# ... variables declaration ...
#
pointer imgl1i()

begin
	
	return imgl1i(EXT_EXTP(ep))

end


# miplt$t -- put line into a 1D image

pointer procedure mipl1i(ep)
pointer ep
#
# ... variables declaration ...
#
pointer impl1i()

begin
	
	return impl1i(EXT_EXTP(ep))

end

# migl2$t -- get line from a 2D image

pointer procedure migl2i(ep, line)
pointer ep
int line
#
# ... variables declaration ...
#
pointer imgl2i()

begin
	
	return imgl2i(EXT_EXTP(ep), line)

end

# mipl2$t -- put line into a 2D image

pointer procedure mipl2i(ep, line)
pointer ep
int line
#
# ... variables declaration ...
#
pointer impl2i()

begin
	
	return impl2i(EXT_EXTP(ep), line)

end

# migl3$t -- get line from a 3D image

pointer procedure migl3i(ep, line, band)
pointer ep
int line, band
#
# ... variables declaration ...
#
pointer imgl3i()

begin
	
	return imgl3i(EXT_EXTP(ep), line, band)

end


# mipl3$t -- put line into a 3D image

pointer procedure mipl3i(ep, line, band)
pointer ep
int line, band
#
# ... variables declaration ...
#
pointer impl3i()

begin

	return impl3i(EXT_EXTP(ep), line, band)

end

# NOTES:
# TODO:: these are not all tested... a related TODO:: note exists in the 
# gtest task t_testmefio.x code (that's where testing code should go)





pointer procedure mignll(ep, buffer, v)
pointer ep,buffer
long v[ARB]

#
# ... variables declaration ...
#
pointer imgnll()
bool ldebug
begin
	ldebug = false
	 
		if (ldebug)
		{
			call printf("at mignl$t()\n")
			call flush(STDOUT)
		}

	
	return imgnll(EXT_EXTP(ep), buffer, v)

end

pointer procedure mipnll(ep, buffer, v)
pointer ep,buffer
long v[ARB]

#
# ... variables declaration ...
#
pointer impnll()
bool ldebug
begin
	ldebug = false
	 
		if (ldebug)
		{
			call printf("in mipnl$t()\n")
			call flush(STDOUT)
		}

	
	return impnll(EXT_EXTP(ep), buffer, v)

end

# migl1$t -- get line from a 1D image

pointer procedure migl1l(ep)
pointer ep
#
# ... variables declaration ...
#
pointer imgl1l()

begin
	
	return imgl1l(EXT_EXTP(ep))

end


# miplt$t -- put line into a 1D image

pointer procedure mipl1l(ep)
pointer ep
#
# ... variables declaration ...
#
pointer impl1l()

begin
	
	return impl1l(EXT_EXTP(ep))

end

# migl2$t -- get line from a 2D image

pointer procedure migl2l(ep, line)
pointer ep
int line
#
# ... variables declaration ...
#
pointer imgl2l()

begin
	
	return imgl2l(EXT_EXTP(ep), line)

end

# mipl2$t -- put line into a 2D image

pointer procedure mipl2l(ep, line)
pointer ep
int line
#
# ... variables declaration ...
#
pointer impl2l()

begin
	
	return impl2l(EXT_EXTP(ep), line)

end

# migl3$t -- get line from a 3D image

pointer procedure migl3l(ep, line, band)
pointer ep
int line, band
#
# ... variables declaration ...
#
pointer imgl3l()

begin
	
	return imgl3l(EXT_EXTP(ep), line, band)

end


# mipl3$t -- put line into a 3D image

pointer procedure mipl3l(ep, line, band)
pointer ep
int line, band
#
# ... variables declaration ...
#
pointer impl3l()

begin

	return impl3l(EXT_EXTP(ep), line, band)

end

# NOTES:
# TODO:: these are not all tested... a related TODO:: note exists in the 
# gtest task t_testmefio.x code (that's where testing code should go)





pointer procedure mignlr(ep, buffer, v)
pointer ep,buffer
long v[ARB]

#
# ... variables declaration ...
#
pointer imgnlr()
bool ldebug
begin
	ldebug = false
	 
		if (ldebug)
		{
			call printf("at mignl$t()\n")
			call flush(STDOUT)
		}

	
	return imgnlr(EXT_EXTP(ep), buffer, v)

end

pointer procedure mipnlr(ep, buffer, v)
pointer ep,buffer
long v[ARB]

#
# ... variables declaration ...
#
pointer impnlr()
bool ldebug
begin
	ldebug = false
	 
		if (ldebug)
		{
			call printf("in mipnl$t()\n")
			call flush(STDOUT)
		}

	
	return impnlr(EXT_EXTP(ep), buffer, v)

end

# migl1$t -- get line from a 1D image

pointer procedure migl1r(ep)
pointer ep
#
# ... variables declaration ...
#
pointer imgl1r()

begin
	
	return imgl1r(EXT_EXTP(ep))

end


# miplt$t -- put line into a 1D image

pointer procedure mipl1r(ep)
pointer ep
#
# ... variables declaration ...
#
pointer impl1r()

begin
	
	return impl1r(EXT_EXTP(ep))

end

# migl2$t -- get line from a 2D image

pointer procedure migl2r(ep, line)
pointer ep
int line
#
# ... variables declaration ...
#
pointer imgl2r()

begin
	
	return imgl2r(EXT_EXTP(ep), line)

end

# mipl2$t -- put line into a 2D image

pointer procedure mipl2r(ep, line)
pointer ep
int line
#
# ... variables declaration ...
#
pointer impl2r()

begin
	
	return impl2r(EXT_EXTP(ep), line)

end

# migl3$t -- get line from a 3D image

pointer procedure migl3r(ep, line, band)
pointer ep
int line, band
#
# ... variables declaration ...
#
pointer imgl3r()

begin
	
	return imgl3r(EXT_EXTP(ep), line, band)

end


# mipl3$t -- put line into a 3D image

pointer procedure mipl3r(ep, line, band)
pointer ep
int line, band
#
# ... variables declaration ...
#
pointer impl3r()

begin

	return impl3r(EXT_EXTP(ep), line, band)

end

# NOTES:
# TODO:: these are not all tested... a related TODO:: note exists in the 
# gtest task t_testmefio.x code (that's where testing code should go)





pointer procedure mignld(ep, buffer, v)
pointer ep,buffer
long v[ARB]

#
# ... variables declaration ...
#
pointer imgnld()
bool ldebug
begin
	ldebug = false
	 
		if (ldebug)
		{
			call printf("at mignl$t()\n")
			call flush(STDOUT)
		}

	
	return imgnld(EXT_EXTP(ep), buffer, v)

end

pointer procedure mipnld(ep, buffer, v)
pointer ep,buffer
long v[ARB]

#
# ... variables declaration ...
#
pointer impnld()
bool ldebug
begin
	ldebug = false
	 
		if (ldebug)
		{
			call printf("in mipnl$t()\n")
			call flush(STDOUT)
		}

	
	return impnld(EXT_EXTP(ep), buffer, v)

end

# migl1$t -- get line from a 1D image

pointer procedure migl1d(ep)
pointer ep
#
# ... variables declaration ...
#
pointer imgl1d()

begin
	
	return imgl1d(EXT_EXTP(ep))

end


# miplt$t -- put line into a 1D image

pointer procedure mipl1d(ep)
pointer ep
#
# ... variables declaration ...
#
pointer impl1d()

begin
	
	return impl1d(EXT_EXTP(ep))

end

# migl2$t -- get line from a 2D image

pointer procedure migl2d(ep, line)
pointer ep
int line
#
# ... variables declaration ...
#
pointer imgl2d()

begin
	
	return imgl2d(EXT_EXTP(ep), line)

end

# mipl2$t -- put line into a 2D image

pointer procedure mipl2d(ep, line)
pointer ep
int line
#
# ... variables declaration ...
#
pointer impl2d()

begin
	
	return impl2d(EXT_EXTP(ep), line)

end

# migl3$t -- get line from a 3D image

pointer procedure migl3d(ep, line, band)
pointer ep
int line, band
#
# ... variables declaration ...
#
pointer imgl3d()

begin
	
	return imgl3d(EXT_EXTP(ep), line, band)

end


# mipl3$t -- put line into a 3D image

pointer procedure mipl3d(ep, line, band)
pointer ep
int line, band
#
# ... variables declaration ...
#
pointer impl3d()

begin

	return impl3d(EXT_EXTP(ep), line, band)

end

# NOTES:
# TODO:: these are not all tested... a related TODO:: note exists in the 
# gtest task t_testmefio.x code (that's where testing code should go)



# this is the generic for for silrd

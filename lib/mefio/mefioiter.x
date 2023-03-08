# Copyright(c) 2004-2009 Association of Universities for Research in Astronomy, Inc.

include <mefio.h>
include <imhdr.h>

# This file contains:
#
#         megetep(...) - Short description
#
# Support routines:
#         megetep(...) - Short description
#
# See also:
#         megetep(...) - Short description

#
# TODO:
#  (*) Need(?) to refcount writes seperate from reads (just one refcount right	now)
#  (*) purge the 'p' for pointer from functionnames 
#

# MEGETEXTP -- Description

pointer procedure megetep (mep, exti, mode)
pointer mep	# pointer to mextp
int exti
int mode
# -- local vars --
int ettype
pointer ep	# return value, pointer to ext_ struct 
pointer extp 	# extension pointer as returned from immap or equiv for other
char 	extname[SZ_LINE]
char	tmpstr[SZ_LINE]
# --- proc defs below ---
pointer paryget()
pointer immap()
int	iget()
bool ldebug
# --- other ---)
include "mefiocommon.h"

begin
    
    ldebug = debug
    
	ep = 0  # clear in case of null return

    
    if (ldebug) {
        call printf("megetep %s %d %d\n")
        call pargstr(ME_FILENAME(mep))
        call pargi(exti)
        call pargi(mode)
        call flush(STDOUT)
    }
    
	if (( exti >= ME_NUMEXTS(mep))  || (exti <= 0))
	{
		call sprintf (tmpstr, SZ_LINE, "megetep(): extension index %d out of range")
		call pargi(exti)
        call log_err(MEERR_OUTOFRANGE, tmpstr)
        call error(MEERR_OUTOFRANGE, "megetep failed")
		# OBSOLETE, PREVIOUS BEHAVIOR
        # OBSOLETE  I've decided to return null to allow the idiom of 
		# OBSOLETE looping over a fixed array of extp... not sure I like this though
		# OBSOLETE  TODO:: consider if this is an error of WARNING level
		# return NULL
	}
	ettype = ME_TYPES(mep, exti)
	
	# this switch allows one call to map any kind of MEFIO supported
	#  extension type (only planned to be IMAGE and TABLE)
	switch(iget(ME_TYPES(mep), exti))
	{
		case ET_IMAGE:
#		call eprintf("megetep(): ET_IMAGE\n")
		ep = paryget(ME_EPS(mep), exti)

		# handle OPENMODE "ANY_MODE"
		if (mode == ANY_MODE){
		#default is read only
			if (ep == NULL){
				mode = READ_ONLY
			} else {
				mode = EXT_OPENMODE(ep)
			}
		}
        
# I have not thought of a decent way to ensure getting the right mode
# if it's requested differently... so just warn about conflict
        if (ep != NULL) {
            if (mode != EXT_OPENMODE(ep)){
                call sprintf(tmpstr,SZ_LINE,"Mode request (%d) mismatchs current (%d)")
                call pargi(mode)
                call pargi(EXT_OPENMODE(ep))
                call log_warn(tmpstr)
                return NULL;
            }                            
        }
	
		if (ep == NULL) {
		# this means it's not mapped yet
		# NOTE:: SECTIONS:: handle image sections here in future
			call sprintf(extname, SZ_LINE, "%s[%d]")
			call pargstr(Memc[ME_PFILENAME(mep)])
			call pargi(exti)
			
			# add Section... it's EOS only string when "NULL", 
			# ... so this strcat should be safe
            # !! Commented out, no section support for now by group decision
			# call strcat(extname, ME_SECSTR(mep), SZ_LINE)
			
			
			if (ldebug)
			{
				call eprintf("megetep(): extension identifier = %s\n")
				call pargstr(extname)
			}

			iferr( extp = immap (extname, mode, NULL))
			{
				call eprintf("megetep(): immap(%s, %d, NULL) failed\n")
				call pargstr(ME_FILENAME(mep))
				call pargi(mode)
				ep = NULL # doc:: immap returns NULL on error
				return NULL
			}
			
			# if I'm here then I have an image extension mapped
			# create the ep (EXT_xxx) structure
			call calloc (ep, LEN_EXT, TY_STRUCT)
			EXT_MST(ep) 	= MST_EXT
			EXT_MEXTP(ep) 	= mep
			EXT_INDEX(ep) 	= exti
			EXT_EXTTYPE(ep)	= ET_IMAGE
			EXT_EXTP(ep) 	= extp
			EXT_OPENMODE(ep)= mode
			# put this pointer in the MEXT_ structure
			call paryset( ME_EPS(mep), exti, ep)
		}	
		
		case ET_TABLE:
        call log_err(MEWRN_UNKNOWN_EXTTYPE,"megetep(): ET_TABLE type not supported yet (UNDER DEVELOPMENT)")
		
		case ET_OTHER:
		call log_err(MEWRN_UNKNOWN_EXTTYPE,"megetep(): ET_OTHER type, type not supported by MEFIO")
		
		default:
		call log_err(MEWRN_UNKNOWN_EXTTYPE,"megetep(): apparently corrupt extension type!(bad mext pointer)")
	}
	
	#if not null, increment refcount
	if (ep != NULL)
	{
		call iset( ME_EXTCOUNTS(mep), exti, iget(ME_EXTCOUNTS(mep), exti) + 1)
	}
    
    if (ldebug) {
        call printf("megetep: about to return %d\n")
        call pargi(ep)
        call flush(STDOUT)
    }
	return ep

end

#MEGNV: MEFIO Get By extName and extVer
# this routines checks to see if EXTNAME matches name,
#  THEN if ver == EXTVER, the appropriate pointer is mapped and returned

pointer procedure megbnv(mep, name, ver, mode)
pointer mep
char	name[ARB]
int ver
int mode

# --- local vars ---
int i
char uname[SZ_LINE]
char lname[SZ_LINE]

# --- proc definitions ---
pointer megetep()

pointer paryget()
int iget()
bool streq()
int strcmp()

bool ldebug

begin

    ldebug =  false
call strcpy(name, uname, SZ_LINE)  
call strcpy(name, lname, SZ_LINE)
    
    call strlwr(lname);
    call strupr(uname);
	
	for (i = 1; i < ME_NUMEXTS(mep); i = i + 1)
	{
#		if (strcmp(Memc[paryget(ME_NAMES(mep),i)], name) == 0 )
            if (ldebug)
            {
                call printf("checking if %s can match %s or %s\n")
                call pargstr(Memc[paryget(ME_NAMES(mep),i)])
                call pargstr(lname)
                call pargstr(uname)
                call flush(STDOUT)
            }
		if ( streq(Memc[paryget(ME_NAMES(mep),i)], uname) || streq(Memc[paryget(ME_NAMES(mep),i)], lname))
		{
            if (ldebug)
            {
                call printf("%s matched (%s,%s)\n")
                call pargstr(Memc[paryget(ME_NAMES(mep),i)])
                call pargstr(lname)
                call pargstr(uname)
                call flush(STDOUT)
            }

            
			if (ver == iget(ME_EXTVERS(mep), i))
			{
				return megetep(mep, i, mode)
			}
		}
	}

#if I'm returning from here I've got nothing.  Question: should this throw
# an error? Answer: I don't thinnk so... let the caller loop for NULL
# they want.
	
	return NULL
end

# memaxver(mep, name, ver)
# MEfio return MAX extVER, by extname

int procedure memaxver(mep, name)
pointer mep
char	name[ARB]

# --- local vars ---
int i
int ver
bool nonenamed

# --- proc definitions ---
pointer paryget()
int iget()
bool streq()
int retval # return value
int max()
bool encmp()
begin
	
    retval = 0 
    
	for (i = 1; i < ME_NUMEXTS(mep); i = i + 1)	{
		if (encmp(Memc[paryget(ME_NAMES(mep),i)], name)) {
			ver = iget(ME_EXTVERS(mep), i)
            retval = max(retval, ver)
        }
	}
    # note: to support raw images
    if (retval == 0) {
        nonenamed = true
        for (i = 1; i < ME_NUMEXTS(mep); i = i + 1)	{
		    if (!streq(Memc[paryget(ME_NAMES(mep),i)], "")) {
                nonenamed = false
            }
	    }
        if (nonenamed == true) {
            retval = ME_NUMEXTS(mep) - 1
        }
    }
    if (false) {
        call printf("memaxver()=%d\n")
        call pargi(retval)
        call flush(STDOUT)
    }
    
	return retval
end


#MEGNV: MEFIO Get INDex By extName and extVer
# this routines checks to see if EXTNAME matches name,
#  THEN if ver == EXTVER (0 means EXTVER not PRESENT)
#  Note: ver == -1 means, first extension with given name

int procedure meindbnv(mep, name, ver)
pointer mep
char	name[ARB]
int ver

# --- local vars ---
int i
int val
pointer enam

# --- proc definitions ---
pointer paryget()
int iget()
bool streq()
bool ldebug
bool encmp()

begin
    ldebug = false 
    if (ldebug){
        call printf ("meindbnv(%x,%s,%d) [ME_NAMES(mep) == %x\n")
        call pargi(mep)
        call pargstr(name)
        call pargi(ver)
        call pargi(ME_NAMES(mep))
        call flush(STDOUT)
    }
	if (mep == NULL)
	{
		call error(MEERR_INVAL, "meindbnv(): ERROR: MEF pointer is NULL.")
	}
		
	
	for (i = 1; i < ME_NUMEXTS(mep); i = i + 1)
	{
        enam = paryget(ME_NAMES(mep),i)
        if (ldebug) {
            call printf("Memc[paryget(ME_NAMES(mep),i)] = %x,%s\n")
            call pargi(enam)
            if (enam == NULL) {
                call pargstr("NULL")
            }
            else {
                call pargstr(Memc[enam])
            }
        }
		if (encmp(Memc[paryget(ME_NAMES(mep),i)], name))
		{
            if (false){
                val = iget(ME_EXTVERS(mep), i)

                call printf("%s %d ver=%d\n")
                call pargstr(Memc[paryget(ME_NAMES(mep),i)])
                call pargi(val)
                call pargi(ver)

                if (IS_INDEF(val)) {
                    call printf("indef!!!!\n")
                }
                call flush(STDOUT)
            }
			if (ver == -1 || ver == iget(ME_EXTVERS(mep), i)) {
				return i
			}
            if (ver == 0 && IS_INDEFI(iget(ME_EXTVERS(mep), i))) {
                return i
            }
                
		}
	}

#if I'm returning from here I've got nothing.  Question: should this throw
# an error? Answer: I don't think so... let the caller loop for NULL
# they want.  WARNING:: Caveat: although callers should know 0 is the PHU... 
#   it's possible that if they pass it on without checking that it will 
#   appear to work with strange/bad results.
	
	return 0
end

# PROCEDURE: COUNTBN
# COUNT By Name... get number of extensions by given extension name

int procedure countbn(mep, name)
pointer mep
char name[ARB]
# -- locals
int extver,ind
int meindbnv()
begin
	extver = 1
	while (TRUE)
	{
		ind = meindbnv(mep, name, extver)
		if (ind == 0)
		{
			break
		}
		extver = extver + 1
	}
	extver = extver - 1
	return extver
end


# PROCEDURE: MEHCT()
# get the highest common type, e.g. the type of the highest-typed
# SCI extension 

int procedure mehct(mep)
pointer mep
# -- locals
int extver,ind
pointer exp
pointer imp2
int ftype, tftype
int imgeti()
pointer megbnv()
pointer memap(), immap(), imunmap()

begin
#3    imp = immap("svqa.fits[SCI,1]", READ_WRITE, NULL)
#   ftype = imgeti(imp, "i_pixtype")
#
#    call printf("type = %d %d\n")
#    call pargi(imp)
#    call pargi(ftype)
#    call flush(STDOUT)
#    call imunmap(imp)
#    return ftype
    
    ftype = 0
    
	extver = 1
	while (TRUE)
	{
		exp = megbnv(mep, "SCI", extver, READ_ONLY)
		if (exp == 0)
		{
			break
		}
        tftype = imgeti(EXT_EXTP(exp), "i_pixtype")
        if ((tftype == 11) || (tftype == 12)) {
           tftype = 4
        }
        if (false) {
            call printf("type[SCI,%d] = %d\n")
            call pargi(extver)
            call pargi(tftype)
        }
        extver = extver + 1
        if (tftype > ftype) {
            ftype = tftype
        }
        call meepunmap(exp)
	}
	extver = extver - 1
    
	return tftype
end



# MESETITER: MEFIO Set Iterator
procedure mesetiter(mep, exti)
pointer mep	# pointer to MEXT_ struct
int exti	# extension index to make current

# --- local vars below

begin
	# check the range
	#  0 is allowed as the "before all" starting point
	if (exti >= ME_NUMEXTS(mep) || exti < 0) {
		call eprintf("recieved out of range index, %d,  for interator ")
		call pargi(exti)
		call error(MEERR_OUTOFRANGE, "mesetiter(): out of range index")
	}
	
	ME_CUREXTI(mep)= exti

end

# MEGNEP: MEFIO Get Next Extension (Pointer)

pointer procedure megnep(mep, mode)
pointer mep
int mode
# -- local decl
pointer megntep()
begin

	return megntep(mep, ET_ANY, mode)
	
end

# MEGNNEP: MEFIO Get Next Type of Extension (Pointer)
# WARNING:: you cannot call MEGNNEP with two different
# WARNING::: types in an interleaved way, you must reset
# WARNING::: the iterator, i.e. mesetiter(mep, 0)

pointer procedure megntep(mep, type, mode )
pointer mep	# pointer to ME_ struct
int type	# type of extension desired
int mode 	# r/w mode  e.g. READ_ONLY
# -- local vars
int nexti
int pget()
pointer ep	# for return value
int extype
pointer megetep()

begin
	
	ep = NULL
	
	while (true)
	{
		nexti = ME_CUREXTI(mep) + 1
		
		if (nexti >= ME_NUMEXTS(mep) || (nexti <= 0))
		{
# then we are past the end
# TODO:: review this logic, it causes megntep and megnep to wrap
# TODO:::around to the start of the script after returning this NULL
			
# note: this is silent, there is no error, simply return NULL as EOF
			ME_CUREXTI(mep) = 0
			return NULL		
		}
		extype= pget(ME_TYPES(mep), nexti)
		if ((type == ET_ANY) || (type == extype))
		{
			ep=megetep(mep, nexti, mode)
			ME_CUREXTI(mep) = nexti	
			return ep
		}
	}
	return NULL; # don't think it can get here!
end



# these go to mefioutil.x
# IDECR: decrement member of int array pointer

procedure idecr(iptr, index)
pointer iptr
int index
# locals
int ival
int iget()

begin

	ival = iget(iptr, index)

	ival = ival-1

	# as I'm using this for reference counting, I don't want negative counts
	if (ival < 0)
	{
		ival = 0
		call eprintf("idecr(): Attempt to Decrement Intary Member to < 0, ignored\n")
	}
	call iset(iptr, index, ival)	
		
end

# IINCR: increment member of int array pointer
procedure iincr(iptr, index)
pointer iptr
int index
# locals
int ival
int iget()
begin

	ival = iget(iptr, index)

	ival = ival+1

	call iset(iptr, index, ival)	
		
end

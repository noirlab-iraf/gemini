# Copyright(c) 2004-2005 Association of Universities for Research in Astronomy, Inc.

include <mef.h>		#from fitsutil

include "mefio.h"	#libmefio main include

# This file contains:
#
#       memap(char filename[ARB])	- Map multiple extention file
#		meprint(pointer mextp)		- Print out information in mextp
#		meunmap(pointer mextp)		- Free the mextp datastructure 
#
# Support routines:
#         memap(...) - Short description
#
# See also:
#         memap(...) - Short description


# NOTES:: TODO:: ################################################
# () ME_NAMES should be ME_EXTNAMES
# ###############################################################

# MEMAP -- Open MEF and load some information about the extensions

pointer procedure memap(mefname)
# ... variables declaration ...
char mefname[ARB]

pointer mef
pointer mextp,extnms
pointer tp 	#temporary pointer

int 	retval,i,numexts, offset, ival
int	len, siz
int	nextend		#hold the NEXTEND value from PHU

# proc declarations
pointer 	pget()

char	cluster[SZ_FNAME]
char	section[SZ_FNAME]
char	tmpstr[SZ_LINE]

# proc includes
include "mefiocommon.h"
bool debughere

int strlen()
bool streq(), strne()
pointer	mef_open()
int mef_rdhdr_gn(), mefgeti(), imaccess()

begin
	debughere = false
	
	
	# initialize variable that need it
	nextend = 0	
	# ---

	if (debughere) 
	{
		call eprintf ("entering memap\n") 
	}
	
	mef = 0;
	
	call imgcluster(mefname, cluster, SZ_FNAME)
	call imgsection(mefname, section, SZ_FNAME)
    
    if ( strne(section,"") ) {
        # note: sections have some internal support, but 
        #  there are unresolved issue regarding how they
        #  should work, so currently, they are prohibited
        #  although you will see them "handled" in other
        #  parts of the code.  The prohibition is at this
        #  point as a gatekeeper.
        call error(MEERR_NOSECTIONS, "memap(): Image sections are not supported")
    }

	iferr (mef = mef_open(cluster, READ_ONLY, 0))
	{
		if (debughere)
		{	
			call eprintf ("memap(): error: mef_open failed for %s\ncluster=%s\tsection\n")
			call pargstr(mefname)
			call pargstr(cluster)
			call pargstr(section)
		}
        call sprintf(tmpstr, SZ_LINE, "mefio: cannot open MEF \"%s\"\n%s")
        call pargstr(cluster)
        if (imaccess(cluster, 0) == YES){
            call pargstr("(image does exist)")
        } else {
            call pargstr("(image does NOT exist)")
        }
		call error(MEERR_NOFILE,tmpstr)
	}
	
	if (debughere) 
	{
		call eprintf ("memap(): mef_open is fine mef pointer = 0x%x\n")
		call pargi(mef)
	}

# get info from PHU if present	
# note:: this PHU retrieval does not work if placed below the loop below after
#       the extensions have been counted (seems to not like reseeking
#       after EOF)
	
	retval = mef_rdhdr_gn(mef, 0) #get PHU
	
	iferr (nextend = mefgeti(mef, "NEXTEND"))
	{
		nextend = 0;
		if (debughere)
		{
			call error(MEERR_NONEXTEND, "memap(): NEXTEND not set in PHU (or other err retrieving it)")
		}
	}
	
	if (debughere)
	{
		call eprintf ("memap(): PHU NEXTEND = %d\n")
		call pargi(nextend)
	}

# rifle through file to physically count number of extensions
# (NEXTEND should represent this but might not be be accurate)
	
	i = 1;

	if (debughere)
	{
		call printf("memap(): pre-count extensions\n")
		call flush(STDOUT)
	}
	while (true) #note: we use break to exit the loop	
	{
		retval =  mef_rdhdr_gn(mef, i)
		if (retval == EOF ) 
		{
			break;
		}
		if (debughere)
		{
			call eprintf("ext#%d MEF_EXTTYPE=%s  MEF_EXTNAME=%s  MEF_EXTVER=%d\n")
			call pargi(i)
			call pargstr(MEF_EXTTYPE(mef))
			call pargstr(MEF_EXTNAME(mef))
			call pargi(MEF_EXTVER(mef))
		}
		
		i = i+1
	}

	numexts = i
	
	if (debughere)
	{
		call eprintf ("memap(): number of extensions (excluding PHU) = %d\n")
		call pargi(numexts)
	}

# KLUDGE: close and open file to clear the EOF (please tell me
#  there is a better way, please? I must have missed something, but
#  this works for now, I did try other things)
	call mef_close(mef)
	iferr (mef = mef_open(cluster, READ_ONLY, 0))
	{
		call eprintf ("memap(): error: re-mef_open failed for %s\n")
		call pargstr(mefname)
		call error(MEERR_NOFILE,"mefio: cannot open file")
	}
	
	
# START: allocate and fill the mext structure
	
# allocate mext structure (mextp = mefio multiple extension pointer)
	call calloc(mextp, LEN_MEXT, TY_STRUCT)
	if (debughere) {
		call eprintf ("memap(): calloc'ed mef (0x%x)\n")
		call pargi(mextp)
	}

# set MEFIO filename and MEFIO structure type type in mextp
# todo: I need to validate this name and also handle the image section parsing!)
	siz = SZ_FNAME # strlen(mefname)+1
	
	call calloc(tp, siz, TY_CHAR)
	ME_PFILENAME(mextp) = tp
	call imgcluster(mefname, Memc[ME_PFILENAME(mextp)], siz)
	call calloc(tp, siz, TY_CHAR)
	ME_PSECSTR(mextp) = tp
	call imgsection(mefname, Memc[ME_PSECSTR(mextp)], siz)
	
	ME_MST(mextp) = MST_MEXT

# set NEXTEND and numexts, and current extension (should already be 0
# from calloc, but... better safe than sorry!)
	ME_NUMEXTS(mextp) = numexts
	ME_NEXTEND(mextp) = nextend
	ME_CUREXTI(mextp)  = 0

# calloc substructures/arrays
	# extension names
	call calloc(extnms, numexts-1, TY_INT)
	ME_NAMES(mextp) = extnms
    if (false) # DEBUG
    {
        call printf("\n\nmemap: numexts = %d\n\n")
        call pargi(numexts)
        call flush(STDOUT)
    }
	# extension type array
	call calloc(tp, numexts-1, TY_INT)
	ME_TYPES(mextp) = tp
	# extension pointer array
	call calloc(tp, numexts-1, TY_INT)
	ME_EPS(mextp) = tp
# calloc pointer to extension EXTVERs array
	call calloc(tp, numexts-1, TY_INT)
	ME_EXTVERS(mextp) = tp
# calloc pointer to extension reference counter array
	call calloc(tp, numexts-1, TY_INT)
	ME_EXTCOUNTS(mextp) = tp
	

# MEXT EXTENSION LOOP: fill structure with extension information
	for (i = 1; i < numexts; i = i + 1)
	{
		#call eprintf("top o' loop #%d of %d\n");
		#call pargi(i)
		#call pargi(numexts-1)
		
		offset = i - 1
		retval = mef_rdhdr_gn(mef, i)

		if (debughere)
		{
			#call eprintf("memap(): retval=%d EOF=%d\n")
			#call pargi(retval)
			#call pargi(EOF)
			call printf("memap(): getting values from ext#%d(%s) of %d\n")
			call pargi(i)
			call pargstr(MEF_EXTNAME(mef))
			call pargi(numexts-1)
            call flush(STDOUT)
		}

# ME_NAMES allocate string for name, put it in the mext structure

		len =  strlen(MEF_EXTNAME(mef))
		if (len > 0)
		{
			ME_NUMNAMEDEXTS(mextp) = ME_NUMNAMEDEXTS(mextp) + 1
		}
		call calloc(tp, len+1, TY_CHAR)
		call strcpy(MEF_EXTNAME(mef), Memc[tp], len+1)
		#call pset(Memi[extnms], i, tp)
		call paryset(extnms, i, tp)
		
#test code: ref/deref issues: sets refcounts to see how they look elsewhere
#  the following two lines are equivalent
#		Memi[ME_EXTCOUNTS(mextp)+offset] = i;
#		call iset(ME_EXTCOUNTS(mextp), i, i)
		
# ME_EXTVERS
        
		Memi[ME_EXTVERS(mextp)+offset] = MEF_EXTVER(mef)

# ME_TYPES
		if (streq(MEF_EXTTYPE(mef), "IMAGE"))
		{
			if (debughere)
			{
				call eprintf("exttype is ET_IMAGE\n")
			}	
			Memi[ME_TYPES(mextp)+offset] = ET_IMAGE
		} 
		else if (streq(MEF_EXTTYPE(mef), "TABLE"))
		{
			if (debughere)
			{
				call eprintf("exttype is ET_TABLE\n")
			}	
			Memi[ME_TYPES(mextp)+offset] = ET_TABLE
		}
		else
		{
			if (debughere)
			{
				call eprintf("exttype is ET_OTHER\n")
			}	
			Memi[ME_TYPES(mextp)+offset] = ET_OTHER
		}

# Clear extp pointer

		
		#if (debughere)
		#{
		#	call eprintf("tp=0x%x\nstr=%s\t(len=%d) offset=%d extnms=0x%x\n")
		#	call pargi(tp)
		#	call pargstr(MEF_EXTNAME(mef))
		#	call pargi(len)
		#	call pargi(offset)
		#	call pargi(extnms)
		#}
# NOTE: ME_EPS is left empty (0's) as extensions are lazy loaded

		if (debughere)
		{
			ival = pget(Memi[extnms],i)
			call eprintf("memap(): loaded info for %s@%d:%s - %s\n")
			call pargstr(mefname)
			call pargi(i)
            call pargstr(Memc[tp])
			call pargstr(Memc[ival])
		}
	}
	
	if (debughere)
	{
		call printf("Number Of Named Extensions %d of %d\n")
		call pargi(ME_NUMNAMEDEXTS(mextp))
		call pargi(ME_NUMEXTS(mextp))
		call flush(STDOUT)
	}
#close mef file
	call mef_close(mef)
    
	return mextp
end

# MEEPUNMAP and MEEPFUNMAP
# TODO:: should we check type of structure? ME_ v EXT_?

procedure meepunmap(ep)
pointer ep

begin
	call meepfunmap(ep,false)
end

procedure meepfunmap(ep, force)
pointer ep		# extension pointer (EXT_ structure)
bool	force 	# used to override reference counting
# -- locals
pointer mep
int iget()

include "mefiocommon.h"

begin
	

	if (debug)
	{
		call eprintf("meepunmap(): ep=0x%x\n")
		call pargi(ep)
	}	
	
	mep = EXT_MEXTP(ep)
	
	if (mep == NULL)
	{
		call eprintf("meepunmap(): NO MEP POINTER!  Dangling EXT_xxx structure\n")
		call error(MEERR_BADSTRUCT, "meepunmap(): NULL mep pointer error")
	}
	
	# reduce the reference count
	if (force)
	{
		call iset(ME_REFCOUNTS(mep), EXT_INDEX(ep), 0)
	}
	else
	{
		call idecr(ME_REFCOUNTS(mep), EXT_INDEX(ep))
	}

	# free the structure
	# ONLY FREE IF REFCOUNTS is 0!!! this structure is shared
	# TODO:: reevaluate sharing this structure...
	if (iget(ME_REFCOUNTS(mep), EXT_INDEX(ep)) == 0)
	{
		# we free this imio provided pointer
		call imunmap(EXT_EXTP(ep))
		# we nullify the pointer to this struct from the mep ME_EPS array
		call paryset(ME_EPS(mep), EXT_INDEX(ep), NULL)
		
		# we free the extp structure
		call mfree(ep, TY_STRUCT)
	}
end


# PROCEDURE meunmap(mextp)
# this frees the mextp structure, and release any mapped extensions

# TODO:: behavior when extensions are still mapped (free or error)

# TODO:: do heap checking to check this doesn't leak... how do you do
#  this in SPP on linux?

procedure meunmap(mextp)

pointer mextp	# the pointer to the MEFIO MEXT structure pointer

#-- local vars (parameters above)
int i, offset
int exttype
pointer extp
char tmpstr[SZ_LINE]
#-- called procs
pointer paryget()
#-- includes/common/debug
include "mefiocommon.h"
bool dbgunmap
int iget()

begin

	dbgunmap = debug
	
	if (mextp == NULL)
	{ 
# why check?  meunmap sets the structure pointer to NULL, so \
#	people can call in with null, if they have called unmap with \
#	the same pointer previously, benign error
		call eprintf("meunmap(): pointer NULL\n")
		return
	}
    
    if ( ME_PHU(mextp) != 0) {
        call imunmap(ME_PHU(mextp))
        ME_PHU(mextp) = 0
    }
# free extension related memory
	for (i = 1; i < ME_NUMEXTS(mextp); i = i + 1)
	{
		offset = i -1
		extp = NULL
		exttype = Memi[ME_TYPES(mextp)+offset]
		extp = paryget(ME_EPS(mextp) , i)
		if (extp != NULL)
		{
# TODO:: see above... right now it's a warning condition to call this
# while extensions are still mapped (they should be unmapped by now)
			call sprintf(tmpstr, SZ_LINE,"meunmap(): WARNING: extension #%d still has %d references when freed")
			call pargi(i)
			call pargi(iget(ME_REFCOUNTS(mextp), i))
            call log_warn(tmpstr)
            
			switch(exttype)
			{
				case ET_OTHER:
				   call log_err(MEERR_BADTYPE, "meunmap(): error: ext pointer to ET_OTHER EXT TYPE\n") 

				case ET_IMAGE:
				if (dbgunmap) {
					call sprintf(tmpstr, SZ_LINE, "meunmap(): Unmapping ext#%d ET_IMAGE\n")
					call pargi(i)
                    call log_info(tmpstr)
                }
				call meepunmap(extp)

				case ET_TABLE:
				if (dbgunmap) {
					call sprintf(tmpstr, SZ_LINE,"meunmap(): Unmapping ext#%d ET_TABLE\n")
					call pargi(i)
                    call log_info(tmpstr)
				}
				default:
                if (dbgunmap) {
                    call log_err(MEERR_BADTYPE, "meunmap(): error: ext pointer to out of range EXT TYPE\n");
			    }
            }
			call paryset(ME_EPS(mextp), i, NULL)
		}
# free string containing image EXTNAME
		call mfree(Memi[ME_NAMES(mextp)+offset], TY_CHAR);
		# shouldn't be necessary, but can help when client code keeps freed pointers
		# clear the ref counts
		call iset(ME_REFCOUNTS(mextp),i, 0)
	}

	if (dbgunmap)
	{
		call eprintf("meunmap(): after loop\n")
	}
	
# free/clean everything at the top level
if (dbgunmap) { call eprintf("meunmap(): before ME_PFILENAME free\n") }	
    call mfree(ME_PFILENAME(mextp), TY_CHAR) # filename string
if (dbgunmap) { call eprintf("meunmap(): before ME_PSECSTR free\n") }	
	call mfree(ME_PSECSTR(mextp), TY_CHAR)
if (dbgunmap) { call eprintf("meunmap(): before ME_NAMES free\n") }
	call mfree(ME_NAMES(mextp), TY_INT)    # array of pointers to EXTNAME seperately freed above
if (dbgunmap) { call eprintf("meunmap(): before ME_TYPES free\n") }
	call mfree(ME_TYPES(mextp), TY_INT)    # the array of type integers
if (dbgunmap) { call eprintf("meunmap(): before ME_EXT free\n") }
	call mfree(ME_EPS(mextp), TY_INT)     # EXTP pointer array (extensions are freed/unmapped above
if (dbgunmap) { call eprintf("meunmap(): before ME_EXTVERS free\n") }
	call mfree(ME_EXTVERS(mextp), TY_INT)  # EXTVERS array 
if (dbgunmap) { call eprintf("meunmap(): before ME_EXTCOUNTS free\n") }
	call mfree(ME_EXTCOUNTS(mextp), TY_INT)  # EXTCOUNTS array
	
	if (dbgunmap)
	{
		call eprintf("meunmap(): after freeing sub-buffers\n")
	}
	
# free the structure itself, now empty (reviewing the code above
#   couldn't hurt though!)
	call mfree(mextp, TY_STRUCT)
	if (dbgunmap)
	{
		call eprintf("meunmap(): after freeing base mext\n")
	}


# now set mextp to NULL to render it dead
	mextp = NULL
end

# PROCEDURE meprint(mextp)
# dumps information about the mext structure pointed at

procedure meprint(mextp)

pointer mextp	# the pointer to the MEFIO MEXT structure pointer

#-- local vars (parameters above)
int i, offset
int exttype
pointer extp

pointer paryget()
int iget()

begin
	
	if (mextp == NULL)
	{
		call eprintf("meprint(): mext pointer NULL\n")
		return;
	}
	
	call eprintf("MEFIO MEXT Structure for %s\n")
	call pargstr(Memc[ME_PFILENAME(mextp)])
	call eprintf("Number of Extensions: %d\n")
	call pargi(ME_NUMEXTS(mextp))
	call eprintf("NEXTEND from PHU:     %d\n")
	call pargi(ME_NEXTEND(mextp))
	call eprintf("Current Extension:    %d\n")
	call pargi(ME_CUREXTI(mextp))
	call eprintf("Extensions Info\n---------------\n")
	for (i = 1; i < ME_NUMEXTS(mextp); i = i + 1)
	{
		offset = i -1
		extp = NULL
		exttype = Memi[ME_TYPES(mextp)+offset]
#		extp = pget(Memi[ME_EPS(mextp)], i) this is equiv: to paryget call below
		extp = paryget(ME_EPS(mextp) , i)
		
		call eprintf("#%d (EXTNAME=%s EXTVER=%d type=%s mapped now?=%s (%d)\n")
		call pargi(i)
		call pargstr(Memc[Memi[ME_NAMES(mextp)+offset]]);
		call pargi(Memi[ME_EXTVERS(mextp)+offset])
		#this switch provides extension type as a string
		switch(exttype)
		{
			case ET_OTHER:
			   call pargstr("ET_OTHER")
			case ET_IMAGE:
			   call pargstr("ET_IMAGE")
			case ET_TABLE:
			   call pargstr("ET_TABLE")
			default:
			   call pargstr("BAD TYPE")
		}
		#this conditional provides the YES/NO for mapped
		if (extp == NULL)
		{
			call pargstr("NO")
		}   
		else
		{
			call pargstr("YES")
		}
		# equiv to below: call pargi(Memi[ME_EXTCOUNTS(mextp)+offset])
		call pargi(iget(ME_EXTCOUNTS(mextp), i));
	}

	
	call flush(STDOUT)
end

# function provides lazy loading of files
# returns either already mapped PHU or maps... no ref counting
# goes away when mef is closed
pointer procedure megphu(mextp)
pointer mextp

# local
char extname[SZ_LINE]
pointer immap()
begin
    if (ME_PHU(mextp)==0) {
		call sprintf(extname, SZ_LINE, "%s[0]")
		call pargstr(Memc[ME_PFILENAME(mextp)])
        ME_PHU(mextp)= immap(extname, READ_WRITE, 0)
	}
    
    return ME_PHU(mextp)
end	


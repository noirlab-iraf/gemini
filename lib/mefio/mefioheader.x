# Copyright(c) 2005 Association of Universities for Research in Astronomy, Inc.

#
# ... include statements ...
#
include <mefio.h>

# This file contains:
#

# todo: put argument list in
#         meaheadT(...) - Add header value (T is the numeric type character)
#         mepheadT(...) - Put header value (T is the numeric type character)
#  		  megheadT(...) - Get header value (T is the numeric type character)
#         meahstr(...) - Add header string value
#         mephstr(...) - Put header string value
#  		  meghstr(...) - Get header string value
#	      medhead(pointer mep, int ext, char keword[ARB]) - Delete Header Value

# ROUTINE_NAME -- Description

procedure meaheads (mep, ext, keyword, value )
pointer mep
int ext
char keyword[ARB]
short value

#
# static variables declaration ...
#
bool ldebug

pointer ep
pointer imp

pointer megphu()
pointer megetep()
begin
    ldebug = false
    ep = 0
    imp = 0
    if (ldebug) {
        call printf("planning to make %s = %f\n")
        call pargstr(keyword)
        call pargd(value)
        call flush(STDOUT)
	}
	if (ext == 0){
        if (ldebug) {
            call printf("PHU header modification...\n")
            call flush(STDOUT)
        }
        # PHU is special case
        ep = 0
        imp = megphu(mep)
    }
    else {
        if (ldebug) {
            call printf("extension %d modification...\n")
            call pargi(ext)
            call flush(STDOUT)
        }
            
        ep = megetep(mep, ext, READ_WRITE)
        imp = EXT_EXTP(ep)    
        if (ep == 0) {
            call error(MEERR_EXTNOTAVAILABLE, "MEFIO: Can't Access Extension")
        }
	}
    
    if (ldebug) {
        call printf("imp pointer is %x\n")
        call pargi(imp)
        call flush(STDOUT)
    }
    
	call imadds(imp, keyword, value)
    
    if (ldebug) {
        call printf("%s = %d\n")
        call pargstr(keyword)
        call pargs(value)
        call flush(STDOUT)
	}
    
    if (ext != 0) {
	    call  meepunmap(ep)
    }

end

procedure meaheadi (mep, ext, keyword, value )
pointer mep
int ext
char keyword[ARB]
int value

#
# static variables declaration ...
#
bool ldebug

pointer ep
pointer imp

pointer megphu()
pointer megetep()
begin
    ldebug = false
    ep = 0
    imp = 0
    if (ldebug) {
        call printf("planning to make %s = %f\n")
        call pargstr(keyword)
        call pargd(value)
        call flush(STDOUT)
	}
	if (ext == 0){
        if (ldebug) {
            call printf("PHU header modification...\n")
            call flush(STDOUT)
        }
        # PHU is special case
        ep = 0
        imp = megphu(mep)
    }
    else {
        if (ldebug) {
            call printf("extension %d modification...\n")
            call pargi(ext)
            call flush(STDOUT)
        }
            
        ep = megetep(mep, ext, READ_WRITE)
        imp = EXT_EXTP(ep)    
        if (ep == 0) {
            call error(MEERR_EXTNOTAVAILABLE, "MEFIO: Can't Access Extension")
        }
	}
    
    if (ldebug) {
        call printf("imp pointer is %x\n")
        call pargi(imp)
        call flush(STDOUT)
    }
    
	call imaddi(imp, keyword, value)
    
    if (ldebug) {
        call printf("%s = %d\n")
        call pargstr(keyword)
        call pargi(value)
        call flush(STDOUT)
	}
    
    if (ext != 0) {
	    call  meepunmap(ep)
    }

end

procedure meaheadl (mep, ext, keyword, value )
pointer mep
int ext
char keyword[ARB]
long value

#
# static variables declaration ...
#
bool ldebug

pointer ep
pointer imp

pointer megphu()
pointer megetep()
begin
    ldebug = false
    ep = 0
    imp = 0
    if (ldebug) {
        call printf("planning to make %s = %f\n")
        call pargstr(keyword)
        call pargd(value)
        call flush(STDOUT)
	}
	if (ext == 0){
        if (ldebug) {
            call printf("PHU header modification...\n")
            call flush(STDOUT)
        }
        # PHU is special case
        ep = 0
        imp = megphu(mep)
    }
    else {
        if (ldebug) {
            call printf("extension %d modification...\n")
            call pargi(ext)
            call flush(STDOUT)
        }
            
        ep = megetep(mep, ext, READ_WRITE)
        imp = EXT_EXTP(ep)    
        if (ep == 0) {
            call error(MEERR_EXTNOTAVAILABLE, "MEFIO: Can't Access Extension")
        }
	}
    
    if (ldebug) {
        call printf("imp pointer is %x\n")
        call pargi(imp)
        call flush(STDOUT)
    }
    
	call imaddl(imp, keyword, value)
    
    if (ldebug) {
        call printf("%s = %d\n")
        call pargstr(keyword)
        call pargl(value)
        call flush(STDOUT)
	}
    
    if (ext != 0) {
	    call  meepunmap(ep)
    }

end

procedure meaheadr (mep, ext, keyword, value )
pointer mep
int ext
char keyword[ARB]
real value

#
# static variables declaration ...
#
bool ldebug

pointer ep
pointer imp

pointer megphu()
pointer megetep()
begin
    ldebug = false
    ep = 0
    imp = 0
    if (ldebug) {
        call printf("planning to make %s = %f\n")
        call pargstr(keyword)
        call pargd(value)
        call flush(STDOUT)
	}
	if (ext == 0){
        if (ldebug) {
            call printf("PHU header modification...\n")
            call flush(STDOUT)
        }
        # PHU is special case
        ep = 0
        imp = megphu(mep)
    }
    else {
        if (ldebug) {
            call printf("extension %d modification...\n")
            call pargi(ext)
            call flush(STDOUT)
        }
            
        ep = megetep(mep, ext, READ_WRITE)
        imp = EXT_EXTP(ep)    
        if (ep == 0) {
            call error(MEERR_EXTNOTAVAILABLE, "MEFIO: Can't Access Extension")
        }
	}
    
    if (ldebug) {
        call printf("imp pointer is %x\n")
        call pargi(imp)
        call flush(STDOUT)
    }
    
	call imaddr(imp, keyword, value)
    
    if (ldebug) {
        call printf("%s = %d\n")
        call pargstr(keyword)
        call pargr(value)
        call flush(STDOUT)
	}
    
    if (ext != 0) {
	    call  meepunmap(ep)
    }

end

procedure meaheadd (mep, ext, keyword, value )
pointer mep
int ext
char keyword[ARB]
double value

#
# static variables declaration ...
#
bool ldebug

pointer ep
pointer imp

pointer megphu()
pointer megetep()
begin
    ldebug = false
    ep = 0
    imp = 0
    if (ldebug) {
        call printf("planning to make %s = %f\n")
        call pargstr(keyword)
        call pargd(value)
        call flush(STDOUT)
	}
	if (ext == 0){
        if (ldebug) {
            call printf("PHU header modification...\n")
            call flush(STDOUT)
        }
        # PHU is special case
        ep = 0
        imp = megphu(mep)
    }
    else {
        if (ldebug) {
            call printf("extension %d modification...\n")
            call pargi(ext)
            call flush(STDOUT)
        }
            
        ep = megetep(mep, ext, READ_WRITE)
        imp = EXT_EXTP(ep)    
        if (ep == 0) {
            call error(MEERR_EXTNOTAVAILABLE, "MEFIO: Can't Access Extension")
        }
	}
    
    if (ldebug) {
        call printf("imp pointer is %x\n")
        call pargi(imp)
        call flush(STDOUT)
    }
    
	call imaddd(imp, keyword, value)
    
    if (ldebug) {
        call printf("%s = %d\n")
        call pargstr(keyword)
        call pargd(value)
        call flush(STDOUT)
	}
    
    if (ext != 0) {
	    call  meepunmap(ep)
    }

end



procedure mepheads (mep, ext, keyword, value )
pointer mep
int ext
char keyword[ARB]
short value

#
# static variables declaration ...
#
pointer ep
pointer imp

pointer megphu()
pointer megetep()
begin
    ep = 0
    imp = 0
    
	if (ext == 0){
        # PHU is special case
        ep = 0
        imp = megphu(mep)
    }
    else {
        ep = megetep(mep, ext, READ_WRITE)
        imp = EXT_EXTP(ep)
        if (ep == 0) {
            call error(MEERR_EXTNOTAVAILABLE, "MEFIO: Can't Access Extension")
        }
    }
	
	call imputs(imp, keyword, value)
	
    if (ext != 0) {
	    call  meepunmap(ep)
    }

end

procedure mepheadi (mep, ext, keyword, value )
pointer mep
int ext
char keyword[ARB]
int value

#
# static variables declaration ...
#
pointer ep
pointer imp

pointer megphu()
pointer megetep()
begin
    ep = 0
    imp = 0
    
	if (ext == 0){
        # PHU is special case
        ep = 0
        imp = megphu(mep)
    }
    else {
        ep = megetep(mep, ext, READ_WRITE)
        imp = EXT_EXTP(ep)
        if (ep == 0) {
            call error(MEERR_EXTNOTAVAILABLE, "MEFIO: Can't Access Extension")
        }
    }
	
	call imputi(imp, keyword, value)
	
    if (ext != 0) {
	    call  meepunmap(ep)
    }

end

procedure mepheadl (mep, ext, keyword, value )
pointer mep
int ext
char keyword[ARB]
long value

#
# static variables declaration ...
#
pointer ep
pointer imp

pointer megphu()
pointer megetep()
begin
    ep = 0
    imp = 0
    
	if (ext == 0){
        # PHU is special case
        ep = 0
        imp = megphu(mep)
    }
    else {
        ep = megetep(mep, ext, READ_WRITE)
        imp = EXT_EXTP(ep)
        if (ep == 0) {
            call error(MEERR_EXTNOTAVAILABLE, "MEFIO: Can't Access Extension")
        }
    }
	
	call imputl(imp, keyword, value)
	
    if (ext != 0) {
	    call  meepunmap(ep)
    }

end

procedure mepheadr (mep, ext, keyword, value )
pointer mep
int ext
char keyword[ARB]
real value

#
# static variables declaration ...
#
pointer ep
pointer imp

pointer megphu()
pointer megetep()
begin
    ep = 0
    imp = 0
    
	if (ext == 0){
        # PHU is special case
        ep = 0
        imp = megphu(mep)
    }
    else {
        ep = megetep(mep, ext, READ_WRITE)
        imp = EXT_EXTP(ep)
        if (ep == 0) {
            call error(MEERR_EXTNOTAVAILABLE, "MEFIO: Can't Access Extension")
        }
    }
	
	call imputr(imp, keyword, value)
	
    if (ext != 0) {
	    call  meepunmap(ep)
    }

end

procedure mepheadd (mep, ext, keyword, value )
pointer mep
int ext
char keyword[ARB]
double value

#
# static variables declaration ...
#
pointer ep
pointer imp

pointer megphu()
pointer megetep()
begin
    ep = 0
    imp = 0
    
	if (ext == 0){
        # PHU is special case
        ep = 0
        imp = megphu(mep)
    }
    else {
        ep = megetep(mep, ext, READ_WRITE)
        imp = EXT_EXTP(ep)
        if (ep == 0) {
            call error(MEERR_EXTNOTAVAILABLE, "MEFIO: Can't Access Extension")
        }
    }
	
	call imputd(imp, keyword, value)
	
    if (ext != 0) {
	    call  meepunmap(ep)
    }

end


# ------ megheadT: mefio get header value of type

short procedure megheads (mep, ext, keyword)
pointer mep
int ext
char keyword[ARB]
short value

#
# static variables declaration ...
#
pointer ep  #ext struct pointer
pointer imp #image pointer
pointer megetep()
short imgets()
int imaccf()
pointer megphu() # get PHU
errchk imgets()
begin
    ep = 0
    imp = 0
    
	if (ext == 0){
        # PHU is special case
        ep = 0
        imp = megphu(mep)
    }
    else {
        ep = megetep(mep, ext, READ_WRITE)
        imp = EXT_EXTP(ep)
    }
    
    if (imp == 0) {
        call error(MEERR_EXTNOTAVAILABLE, "MEFIO:  Can't Access Extension")
    }
	
    if (imaccf(imp,keyword) == NO ) {
        call meepunmap(ep)
        call error(MEERR_VALNOTFOUND, "MEFIO: header keyword/value didn't exist");
    }
        
    value = imgets(imp, keyword)

    if (ext != 0) {
	    call  meepunmap(ep)
    }

	return value
end

int procedure megheadi (mep, ext, keyword)
pointer mep
int ext
char keyword[ARB]
int value

#
# static variables declaration ...
#
pointer ep  #ext struct pointer
pointer imp #image pointer
pointer megetep()
int imgeti()
int imaccf()
pointer megphu() # get PHU
errchk imgeti()
begin
    ep = 0
    imp = 0
    
	if (ext == 0){
        # PHU is special case
        ep = 0
        imp = megphu(mep)
    }
    else {
        ep = megetep(mep, ext, READ_WRITE)
        imp = EXT_EXTP(ep)
    }
    
    if (imp == 0) {
        call error(MEERR_EXTNOTAVAILABLE, "MEFIO:  Can't Access Extension")
    }
	
    if (imaccf(imp,keyword) == NO ) {
        call meepunmap(ep)
        call error(MEERR_VALNOTFOUND, "MEFIO: header keyword/value didn't exist");
    }
        
    value = imgeti(imp, keyword)

    if (ext != 0) {
	    call  meepunmap(ep)
    }

	return value
end

long procedure megheadl (mep, ext, keyword)
pointer mep
int ext
char keyword[ARB]
long value

#
# static variables declaration ...
#
pointer ep  #ext struct pointer
pointer imp #image pointer
pointer megetep()
long imgetl()
int imaccf()
pointer megphu() # get PHU
errchk imgetl()
begin
    ep = 0
    imp = 0
    
	if (ext == 0){
        # PHU is special case
        ep = 0
        imp = megphu(mep)
    }
    else {
        ep = megetep(mep, ext, READ_WRITE)
        imp = EXT_EXTP(ep)
    }
    
    if (imp == 0) {
        call error(MEERR_EXTNOTAVAILABLE, "MEFIO:  Can't Access Extension")
    }
	
    if (imaccf(imp,keyword) == NO ) {
        call meepunmap(ep)
        call error(MEERR_VALNOTFOUND, "MEFIO: header keyword/value didn't exist");
    }
        
    value = imgetl(imp, keyword)

    if (ext != 0) {
	    call  meepunmap(ep)
    }

	return value
end

real procedure megheadr (mep, ext, keyword)
pointer mep
int ext
char keyword[ARB]
real value

#
# static variables declaration ...
#
pointer ep  #ext struct pointer
pointer imp #image pointer
pointer megetep()
real imgetr()
int imaccf()
pointer megphu() # get PHU
errchk imgetr()
begin
    ep = 0
    imp = 0
    
	if (ext == 0){
        # PHU is special case
        ep = 0
        imp = megphu(mep)
    }
    else {
        ep = megetep(mep, ext, READ_WRITE)
        imp = EXT_EXTP(ep)
    }
    
    if (imp == 0) {
        call error(MEERR_EXTNOTAVAILABLE, "MEFIO:  Can't Access Extension")
    }
	
    if (imaccf(imp,keyword) == NO ) {
        call meepunmap(ep)
        call error(MEERR_VALNOTFOUND, "MEFIO: header keyword/value didn't exist");
    }
        
    value = imgetr(imp, keyword)

    if (ext != 0) {
	    call  meepunmap(ep)
    }

	return value
end

double procedure megheadd (mep, ext, keyword)
pointer mep
int ext
char keyword[ARB]
double value

#
# static variables declaration ...
#
pointer ep  #ext struct pointer
pointer imp #image pointer
pointer megetep()
double imgetd()
int imaccf()
pointer megphu() # get PHU
errchk imgetd()
begin
    ep = 0
    imp = 0
    
	if (ext == 0){
        # PHU is special case
        ep = 0
        imp = megphu(mep)
    }
    else {
        ep = megetep(mep, ext, READ_WRITE)
        imp = EXT_EXTP(ep)
    }
    
    if (imp == 0) {
        call error(MEERR_EXTNOTAVAILABLE, "MEFIO:  Can't Access Extension")
    }
	
    if (imaccf(imp,keyword) == NO ) {
        call meepunmap(ep)
        call error(MEERR_VALNOTFOUND, "MEFIO: header keyword/value didn't exist");
    }
        
    value = imgetd(imp, keyword)

    if (ext != 0) {
	    call  meepunmap(ep)
    }

	return value
end


# STRING VERSIONS
procedure meahstr (mep, ext, keyword, value )
pointer mep
int ext
char keyword[ARB]
char value[ARB]

#
# static variables declaration ...
#
bool ldebug

pointer ep
pointer imp #image pointer
pointer megetep()
int imaccf()
pointer megphu() # get PHU
errchk imgeti()

begin

    ldebug = false
    ep = 0
    imp = 0
    if (ldebug) {
        call printf("planning to make %s = %s\n")
        call pargstr(keyword)
        call pargstr(value)
        call flush(STDOUT)
	}
	if (ext == 0){
        if (ldebug) {
            call printf("PHU string header modification...\n")
            call flush(STDOUT)
        }
        # PHU is special case
        ep = 0
        imp = megphu(mep)
    }
    else {
        if (ldebug) {
            call printf("extension %d modification...\n")
            call pargi(ext)
            call flush(STDOUT)
        }
            
        ep = megetep(mep, ext, READ_WRITE)
        imp = EXT_EXTP(ep)    
        if (ep == 0) {
            call error(MEERR_EXTNOTAVAILABLE, "MEFIO: meahstr(): Can't Access Extension")
        }
	}
        	
	call imastr(imp, keyword, value)
	
    if (ep != 0) {
	    call meepunmap(ep)
    }

end

procedure mephstr (mep, ext, keyword, value )
pointer mep
int ext
char keyword[ARB]
char value[ARB]

#
# static variables declaration ...
#
pointer ep
pointer megetep()
begin

	ep = megetep(mep, ext, READ_WRITE)
    
    if (ep == 0) {
        call error(MEERR_EXTNOTAVAILABLE, "MEFIO: Can't Access Extension")
    }
	
	call impstr(EXT_EXTP(ep), keyword, value)
	
	call  meepunmap(ep)

end

procedure meghstr (mep, ext, keyword, outbuf, maxch)
pointer mep
int ext
char keyword[ARB]
char outbuf[maxch]
int maxch
#
# static variables declaration ...
#
pointer ep
pointer megetep()

begin

	ep = megetep(mep, ext, READ_WRITE)
    
    if (ep == 0) {
        call error(MEERR_EXTNOTAVAILABLE, "MEFIO: Can't Access Extension")
    }
	
	call imgstr(EXT_EXTP(ep), keyword, outbuf, maxch)
    
#    call printf("outbuf = %s\nlen=%d\n")
#    call pargstr(outbuf)
#    call pargi(maxch)
#    call flush(STDOUT)
    
	call  meepunmap(ep)

end

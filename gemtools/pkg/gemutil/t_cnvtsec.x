# Copyright(c) 2005-2006 Association of Universities for Research in Astronomy, Inc.
#
# CNVTSEC -- Convert a FITS format date+time string to a number of seconds 
#            since 00:00:00 01-Jan-80, and print it.

include <time.h>

procedure t_cnvtsec ()

#char    dateobs[SZ_TIME+SZ_DATE]       # Value of the DATE-OBS keyword
#char    timeobs[SZ_TIME]               # Value of the TIME-OBS keyword
#int     status                         # Return status

# Local variables for task paramters
char    l_dateobs[SZ_TIME+SZ_DATE]
char    l_timeobs[SZ_TIME]
int     l_status

# Other variables
int     index, junk
long    ltime
double  dtime, fracsec
pointer dtstr, ndtstr, sp

# Gemini functions
bool    g_whitespace()

# IRAF functions
int     dtm_ltime(), stridx(), ctod()

begin
    l_status = 0
    ltime = 0
    
    # Allocate stack memory
    call smark (sp)
    call salloc (dtstr, SZ_TIME+SZ_DATE, TY_CHAR)
    call salloc (ndtstr, SZ_TIME+SZ_DATE, TY_CHAR)
    
    # Get task parameter values
    call clgstr ("dateobs", l_dateobs, SZ_TIME+SZ_DATE)
    call clgstr ("timeobs", l_timeobs, SZ_TIME)
    
    # Build string if necessary
    index = stridx ("T", l_dateobs)
    if ( index == 0 ) {      # time not in dateobs, get it from timeobs
        if ( g_whitespace(l_timeobs) ) {
            l_status = 121
            call printf ("CNVTSEC ERROR: %d No time information.\n")
                call pargi (l_status)
            call clputi ("status", l_status)
            return
        }
        call sprintf ( Memc[dtstr], SZ_TIME+SZ_DATE, "%sT%s" )
            call pargstr (l_dateobs)
            call pargstr (l_timeobs)
    } else {
        call sprintf ( Memc[dtstr], SZ_TIME+SZ_DATE, "%s" )
            call pargstr (l_dateobs)  
    } 

    # dtm_ltime ignores fractions of seconds.  To keep those fractions
    # of seconds, we extract and remove them from dtstr, then add them
    # back to the total integer number of seconds outputed by dtm_ltime
    
    index = stridx (".", Memc[dtstr])
    if (index != 0) {
        call strcpy (Memc[dtstr], Memc[ndtstr], index-1)
        junk = ctod(Memc[dtstr], index, fracsec)
    } else {
        call strcpy (Memc[dtstr], Memc[ndtstr], SZ_TIME+SZ_DATE)
        fracsec=0.
    }
        
    
    # Convert to seconds
    l_status = dtm_ltime(Memc[ndtstr], ltime)
    dtime = double(ltime)+fracsec
    
    call clputi ("status", l_status)
    call printf ("%.4f\n")
        call pargd (dtime)
    
    call sfree(sp)
    return    
end

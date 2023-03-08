# Copyright(c) 2007 Association of Universities for Research in Astronomy, Inc.

procedure getfakeUT ()

# Generate the fake UT date used to name Gemini data.
# At Gemini, the transit time should be 14:00:00 local time.
# For GN, that corresponds to midnight UT so the name is not faked, but for
# GS, a transit of 14hr is totally artificial.
#
# Before transit, UT of last night
# After transit, UT of coming night
#
# Note that the transit time is not hardcoded and the code should continue
# to work if Gemini's policy changes the transit time in the future
#
# Original author:  Kathleen Labrie
# Version:  2007.03.29 KL

char    transit     {"14:00:00", prompt="Local time of UT transit"}
char    fakeUT      {"", prompt="Output fake UT string"}
bool    verbose     {no, prompt="Verbose"}
int     status      {0, prompt="Exit status (0=good)"}

begin

    char    l_transit = ""
    char    l_fakeUT = ""
    bool    l_verbose
    
    char    locstr, utstr
    char    tmplocdate, tmputdate
    int     junk, idx
    int     ilocdate, iutdate
    real    ftransit, floctime, futtime
    int     day, month, year   
    int     dymnth[12]= 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31


    junk = fscan (transit, l_transit)
    l_verbose = verbose

    cache ("gemdate")
    
    # Convert transit string to a floating point value
    ftransit = real(l_transit)
    
    # Get local and UT date/time strings (FITS format: YYYY-MM-DDThh:mm:ss)
    gemdate (zone="local")
    locstr = gemdate.outdate
    gemdate (zone="UT")
    utstr = gemdate.outdate
    
    # Parse local date/time strings.  Date becomes an integer (YYYYMMDD),
    # and time becomes a floating point value (hh.hhhh)
    
    idx = stridx ("T", locstr)
    tmplocdate = substr (locstr, 1, idx-1)
    print (tmplocdate) | translit ("STDIN", "-", "", delete+, collapse-) | \
        scan (ilocdate)
    floctime = real(substr(locstr, idx+1, strlen(locstr)))
    
    # Parse UT date/time strings.  Date becomes an integer (YYYYMMDD),
    # and time becomes a floating point value (hh.hhhh)
    
    idx = stridx ("T", utstr)
    tmputdate = substr (utstr, 1, idx-1)
    print (tmputdate) | translit ("STDIN", "-", "", delete+, collapse-) | \
        scan (iutdate)
    futtime = real(substr(utstr, idx+1, strlen(locstr)))
    
    # Now we have everything we need.  So let's generate the fake UT date.    
    if (floctime < ftransit) {
    
        if (iutdate == ilocdate)
            l_fakeUT = str(iutdate)
        else  {                 # UT has changed before transit => fake the UT
            year = int(iutdate / 1.e4)
            month = int((iutdate - year*1.e4) / 1.e2)
            day = int(iutdate - year*1.e4 - month*1.e2)
            
            if (day == 1) {
                if (month == 1) {
                    year = year - 1
                    month = 12
                } else
                    month = month - 1
                
                if ((year % 4) == 0) dymnth[2] = dymnth[2] + 1

                day = dymnth[month]
                
                l_fakeUT = str(int((year*1e4)+(month*1e2)+day))
            } else
                l_fakeUT = str(iutdate - 1)

        }
            
    } else {    # flocaltime >= ftransit
    
        if (iutdate == ilocdate) { # UT has not changed yet; transit reached => fake the UT
            year = int(iutdate / 1.e4)
            month = int((iutdate - year*1.e4) / 1.e2)
            day = int(iutdate - year*1.e4 - month*1.e2)

            if ((year % 4) == 0) dymnth[2] = dymnth[2] + 1
            if (day == dymnth[month]) {
                if (month == 12) {
                    year = year + 1
                    month = 1
                } else
                    month = month + 1
                
                day = 1
                
                l_fakeUT = str(int((year*1e4)+(month*1e2)+day))
            } else            
                l_fakeUT = str(iutdate + 1)
                
        } else
            l_fakeUT = str(iutdate)
            
    }
    
    # Save result to output parameter
    fakeUT = l_fakeUT
    
end

# Copyright(c) 2005-2006 Association of Universities for Research in Astronomy, Inc.

procedure gemdate ()

# Get FITS format date string

char    zone        {"UT", enum="UT|local", prompt = "Time zone"}
bool    verbose     {no, prompt = "Verbose"}
char    outdate     {"", prompt = "Output formatted date"}
int     status      {0, prompt = "Exit status (0=good)"}

begin
    
    char    l_zone
    bool    l_verbose
    
    char    zoneopt, fmt
    struct  sdate
    
    status = 0
    
    l_zone = zone
    l_verbose = verbose
    
    fmt = "+%Y-%m-%dT%H:%M:%S"
    if (l_zone == "UT")
        zoneopt = "-u"
    else
        zoneopt = ""

    date (zoneopt, fmt) | scan(sdate)
    outdate = sdate
    
    if (outdate == "")
        status = 99
        
    if (l_verbose)
        print (sdate)    

end

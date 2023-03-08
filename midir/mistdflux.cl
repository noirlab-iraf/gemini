# Copyright(c) 2007-2011 Association of Universities for Research in Astronomy, Inc.

procedure mistdflux (stdname, filters)

char    stdname     {prompt="Name of Cohen or TIMMI2 standard"}
char    filters     {prompt="Filter names, or 'all' for all filters"}
bool    fl_matchall {no, prompt="All name matches?"}
char    atmos       {"none", enum="none|normalized|extinction", prompt="Model atmosphere for output in-band flux density"}
char    instrument  {"trecs", enum="trecs|michelle", prompt="Which instrument? trecs|michelle"}
char    logfile     {"", prompt="Name of logfile"}
bool    verbose     {no, prompt="Verbose"}
int     status      {0, prompt="Exit status"}

begin
    char    l_stdname = ""
    char    l_filters = ""
    char    l_atmos = ""
    char    l_instrument = ""
    char    l_logfile = ""
    bool    l_fl_matchall, l_verbose

    char    paramstr, errmsg
    char    tmpmatch, tmpstring
    char    object, units, filtername[20], chline, frmt
    char    filtinput[20]
    int     junk, i, j, index, len, idx
    int     nfound, entry, nfilters, nreqfilt, nfiltinput
    int     reqfilt[20]
    real    fnu[20]
    bool    matchfilt
    struct  scanline

    cache ("gloginit")

    junk = fscan (stdname, l_stdname)
    l_filters = filters     # Want to keep leading/trailing 'space' characters
    l_fl_matchall = fl_matchall
    junk = fscan (atmos, l_atmos)
    junk = fscan (instrument, l_instrument)
    junk = fscan (logfile, l_logfile)
    l_verbose = verbose

    # Initialize
    status = 0
    
    # Create temp file names
    tmpmatch = mktemp("tmpmatch")
    
    # Create the list of parameter/value pairs. One pair per line.
    # All lines combined into one string.  Line delimiter is '\n'.
    paramstr  = "stdname        = "//stdname.p_value//"\n"
    paramstr += "filters        = "//filters.p_value//"\n"
    paramstr += "fl_matchall    = "//fl_matchall.p_value//"\n"
    paramstr += "atmos          = "//atmos.p_value//"\n"
    paramstr += "instrument     = "//instrument.p_value//"\n"
    paramstr += "logfile        = "//logfile.p_value//"\n"
    paramstr += "verbose        = "//verbose.p_value
    
    # Assign a logfile name if not given.  Open logfile and start log.
    # Write parameter/value pairs ("paramstr") to log.
    gloginit (l_logfile, "mistdflux", "midir", paramstr, fl_append=yes,
        verbose=l_verbose)
    if (gloginit.status != 0) {
        status = gloginit.status
        goto exitnow
    }
    l_logfile = gloginit.logfile.p_value
    

    # More initialization based on the instrument    
    if (l_instrument == "trecs") {
        nfilters = 20
        glogprint (l_logfile, "mistdflux", "science", type="string",
            str="  Instrument = T-ReCS", verbose=l_verbose)
        filtername[1]  = " N (broad 10um)  "
        filtername[2]  = " Q (Broad 20.8um)"
        filtername[3]  = "Qshort 17.65um   "
        filtername[4]  = "Qa 18.30um       "
        filtername[5]  = "Qb 24.56um       "
        filtername[6]  = "Si-1 7.73um      "
        filtername[7]  = "Si-2 8.74um      "
        filtername[8]  = "Si-3 9.69um      "
        filtername[9]  = "Si-4 10.38um     "
        filtername[10] = "Si-5 11.66um     "
        filtername[11] = "Si-6 12.5um      "
        filtername[12] = "[ArIII] 8.99um   "
        filtername[13] = "[NeII] 12.81um   "
        filtername[14] = "[NeII]cont13.10um"
        filtername[15] = "[SIV] 10.52um    "
        filtername[16] = "PAH 8.6um        "
        filtername[17] = "PAH 11.3um       "
        filtername[18] = " K (2.2um)       "
        filtername[19] = " L (3.4um)       "
        filtername[20] = " M (4.6um)       "
    } else if (l_instrument == "michelle") {
        nfilters = 10
        glogprint (l_logfile, "mistdflux", "science", type="string",
            str="  Instrument = Michelle", verbose=l_verbose)
        filtername[1]  = " N (broad 10um)  "
        filtername[2]  = " Nprime          "
        filtername[3]  = "Qa 18.50um       "
        filtername[4]  = " Q (Broad 20um)  "
        filtername[5]  = "Si-1 7.9um       "
        filtername[6]  = "Si-2 8.8um       "
        filtername[7]  = "Si-3 9.7um       "
        filtername[8]  = "Si-4 10.3um      "
        filtername[9]  = "Si-5 11.6um      "
        filtername[10] = "Si-6 12.5um      "
    } else {    # should not happen because of the 'enum' in the param defn.
        errmsg = "Unrecognized instrument name ("//l_instrument//")"
        status = 99
        glogprint (l_logfile, "mistdflux", "engineering", type="error",
            errno=status, str=errmsg, verbose=yes)
        goto clean
    }

    
    # Parse the parameter 'filters' and set the list of requested filters
    if (strlwr(l_filters) == "all") {
        for (i=1; i<=nfilters; i=i+1) {
            reqfilt[i] = i
        }
        nreqfilt = nfilters
    } else {
        # Parse l_filters (can be a comma delimited list)
        nfiltinput = 0
        tmpstring = l_filters
        len = strlen (tmpstring)
        while (len > 0) {
            nfiltinput = nfiltinput+1
            idx = stridx (",", tmpstring)
            if (idx > 0) {
                filtinput[nfiltinput] = substr(tmpstring, 1, idx-1)
                tmpstring = substr(tmpstring, idx+1, len)
                len = strlen(tmpstring)
            } else {
                filtinput[nfiltinput] = substr(tmpstring, 1, len)
                len = 0
            }
        }
        if (nfiltinput == 0) {
            errmsg = "No filter specified"
            status = 121
            glogprint (l_logfile, "mistdflux", "status", type="error",
                errno=status, str=errmsg, verbose=yes)
            goto clean
        }
    
    
        nreqfilt = 0
        for (i=1; i<=nfiltinput; i=i+1) {
            matchfilt = no
            for (j=1; j<=nfilters; j=j+1) {
                if (strstr(strlwr(filtinput[i]), strlwr(filtername[j])) > 0) {
                    nreqfilt = nreqfilt + 1
                    reqfilt[nreqfilt] = j
                    matchfilt = yes
                }
            }
            if (matchfilt == no) {
                errmsg = "No match found for '"//filtinput[i]
                status = 121
                glogprint (l_logfile, "mistdflux", "status", type="warning",
                    errno=status, str=errmsg, verbose=yes)
            }
        }        
        if (nreqfilt == 0) {
            errmsg = "None of the filters requested are valid."
            status = 121
            glogprint (l_logfile, "mistdflux", "status", type="error",
                errno=status, str=errmsg, verbose=yes)
            glogprint (l_logfile, "mistdflux", "status", type="string",
                str="    Try setting 'filter' to 'all' to get a list of valid \
                filter names", verbose=yes)
            goto clean
        }
    }


    # Set the file name for the flux data based on instrument and 'atmos' 
    if (l_atmos == "none") {
        if (l_instrument == "trecs")
            chline = "midir$data/trecs/all.list"
        else    # michelle
            chline = "midir$data/michelle/all.list"
    } else if (l_atmos == "normalized") {
        if (l_instrument == "trecs")
            chline = "midir$data/trecs/all.list.atmos2p0"
        else    # michelle
            chline = "midir$data/michelle/allatmos.list"
    } else if (l_atmos == "extinction") {
        if (l_instrument == "trecs")
            chline = "midir$data/trecs/all.list.altatmos2p0"
        else
            chline = "midir$data/michelle/allaltatmos.list"
    } else {    # should not happen with 'enum' in parameter definition
        errmsg = "Invalid value for parameter 'atmos' ("//l_atmos//")"
        status = 99
        glogprint (l_logfile, "mistdflux", "engineering", type="error",
            errno=status, str=errmsg, verbose=yes)
        goto clean
    }

    # Ensure that 'chline' exists
    if (access(chline) == no) {
        errmsg = "Package installation problem.  The flux density file was \
            not found ("//chline//")"
        status = 99
        glogprint (l_logfile, "mistdflux", "status", type="error",
            errno=status, str=errmsg, verbose=yes)
        goto clean
    }


    # Search for stdname in the flux file (chline)
    nfound = 0
    match ("{"//l_stdname//"}", chline, stop=no, print=no, meta=yes, \
        > tmpmatch)
    wc (tmpmatch) | scan (nfound)       # nfound should be > 1
    if (nfound == 0) {
        errmsg = "No standard star name matches the string '"//l_stdname//"'"
        status = 121
        glogprint (l_logfile, "mistdflux", "status", type="error",
            errno=status, str=errmsg, verbose=yes)
        goto clean
    } else
        if (l_fl_matchall == no)    nfound = 1

    if (l_verbose == no) print ("")        
    for (entry=1; entry<=nfound; entry=entry+1) {
        fields (tmpmatch, "1", lines=str(entry)) | scan (object)
        for (i=1; i<=nfilters; i=i+1) {
            fields (tmpmatch, str(i+1), lines=str(entry)) | scan (fnu[i])
        }
        if (stridx ("*", object) > 0) {
            units = "mJy"
            object = substr(object, 1, strlen(object)-1)
        } else
            units = "Jy"
        
        for (i=1; i<=nreqfilt; i=i+1) {
            frmt = "  %-12s flux density (%s) = "
            if ((fnu[reqfilt[i]] < 1) || (fnu[reqfilt[i]] >= 1e4)) {
                frmt = frmt // "%.4e"
            } else if (fnu[reqfilt[i]] > 1000.) {
                frmt = frmt // "%.0f."
            } else {
                frmt = frmt // "%5f"
            }
            frmt = frmt // " %s\n"
            printf (frmt, object, filtername[reqfilt[i]], fnu[reqfilt[i]],
                units) | scan (scanline)
            glogprint (l_logfile, "mistdflux", "science", type="string",
                str=scanline, verbose=yes)
        }
    }
    if (l_verbose == no) print ("")
    delete (tmpmatch, verify-, >& "dev$null")

clean:
    delete (tmpmatch, verify-, >& "dev$null")
        
    if (status == 0)
        glogclose (l_logfile, "mistdflux", fl_success=yes, verbose=l_verbose)
    else
        glogclose (l_logfile, "mistdflux", fl_success=no, verbose=l_verbose)
    

exitnow:
    ;

end

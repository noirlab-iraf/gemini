# Copyright(c) 2006-2011 Association of Universities for Research in Astronomy, Inc.

procedure mcheckheader (inimages)

# This script checks the headers of Michelle raw data files, and if the 
# ENDTIME keyword is not present in the primary header it inserts an 
# approximate value.  Other header keywords are inserted as well to 
# get something that passes CADC muster.
#
#
# Version:  October 6, 2005 KV wrote the original script
#           Developers see CVS for recent updates


char    inimages    {prompt="Input Michell image(s)"}       # OLDP-1-input-primary-single-prefix=s
char    rawpath     {"", prompt="Path for in raw images"}   # OLDP-4
bool    update      {yes, prompt="Update headers?"}         # OLDP-4
char    logfile     {"", prompt="Logfile"}                  # OLDP-1
bool    verbose     {yes, prompt="Verbose"}                 # OLDP-4
int     status      {0, prompt="Exit status (0=good)"}      # OLDP-4
struct  *scanfile   {"", prompt="Internal use only"}        # OLDP-4

begin

    char    l_inimages = ""
    char    l_rawpath = ""
    char    l_logfile = ""
    bool    l_verbose, l_update

    char    in[200], out[200], filename
    char    tmplog, tmpinimg
    char    paramstr, time, newtime, errmsg, lastchar
    int     junk, i, k, l, m, n
    int     hr, min, sec, duration
    int     nimages, maximages
    real    start
    struct  l_struct


    # Open temp files
    tmplog = mktemp("tmplog")
    tmpinimg = mktemp("tmpinimg")

    junk = fscan (inimages, l_inimages)
    junk = fscan (rawpath, l_rawpath)
    l_update = update
    junk = fscan (logfile, l_logfile)
    l_verbose = verbose

    nimages = 0
    maximages = 200
    status = 0
    
    cache ("gemdate")
    
    # Create the list of parameter/value pairs.  One pair per line.
    # All lines combined into one string.  Line delimiter is '\n'.
    paramstr =  "inimages       = "//inimages.p_value//"\n"
    paramstr += "rawpath        = "//rawpath.p_value//"\n"
    paramstr += "update         = "//update.p_value//"\n"
    paramstr += "logfile        = "//logfile.p_value//"\n"
    paramstr += "verbose        = "//verbose.p_value

    # Assign a logfile name if not given. Open logfile and start log.
    # Write parameter/value pairs ("paramstr") to log.
    gloginit (l_logfile, "mcheckheader", "midir", paramstr, fl_append+,
        verbose=l_verbose)
    if (gloginit.status != 0) {
        status = gloginit.status
        goto exitnow
    }
    l_logfile = gloginit.logfile.p_value

    #----------------------------------------------------------------------
    # Check the rawpath name for a final /
    lastchar = substr (l_rawpath, strlen(l_rawpath), strlen(l_rawpath))
    if ((l_rawpath != "") && (lastchar != "$") && (lastchar != "/"))
        l_rawpath = l_rawpath//"/"
    #-----------------------------------------------------------------------
    # Load up arrays of input name lists

    gemextn (l_inimages, check="exists,mef", process="none", index="",
        extname="", extversion="", ikparam="", omit="kernel,exten",
        replace="^%%" // l_rawpath // "%", outfile=tmpinimg, logfile=l_logfile,
        verbose=l_verbose)
    nimages = gemextn.count
    
    # Check that the list file exists and contains filenames
    if ((gemextn.fail_count > 0) || (nimages == 0) || \
        (nimages > maximages)) {
        
        if (gemextn.fail_count > 0) {
            errmsg = gemextn.fail_count//" images were not found."
            status = 101
        } else if (nimages == 0) {
            errmsg = "No input images defined."
            status = 121
        } else if (nimages > maximages) {
            errmsg = "Maximum number of input images ("//str(maximages)//") \
                has been exceeded."
            status = 121
        }
        
        glogprint (l_logfile, "mcheckheader", "status", type="error", 
            errno=status, str=errmsg, verbose+)
        goto clean
    } else {
        scanfile = tmpinimg
        i = 0
        while (fscan (scanfile, filename) != EOF) {
            i += 1
            in[i] = filename
        }
        scanfile = ""
        if (i != nimages) {
            status = 99
            errmsg = "Error while counting the input images."
            glogprint (l_logfile, "mcheckheader", "status", type="error",
                errno=status, str=errmsg, verbose+)
        }
    }
    delete (tmpinimg, ver-, >& "dev$null")

    # Starting task

    i = 1
    while (i <= nimages) {
        imgets (in[i]//"[0]", "ENDTIME", >& "dev$null")
        if (imgets.value == "0") {
            imgets (in[i]//"[0]", "STARTTIM", >& "dev$null")
            if (imgets.value == "0") {
                status = 131
                errmsg = "File "//in[i]//" has no STARTTIM value."
                glogprint (l_logfile, "mcheckheader", "status", type="error",
                    errno=status, str=errmsg, verbose+)
                goto nextimage
            } else {
                k = fscanf (imgets.value, "%f", start)
                l = 1
                while (l > 0) {
                    if (!imaccess(in[i]//"["//str(l)//"]")) {
                        n = l-1
                        l = -50
                    }
                    l = l+1
                }
                # Estimate the end time from the start time, 45 seconds per 
                # extention. This assumes 45 seconds per NOD position, a crude 
                # general value
                start = start + (n * 45. / 3600.)
                if (start > 24.0) start = start-24.
                
                # we should not be taking observations at UT 0 hours, but just 
                # to be safe
                
                gemhedit (in[i]//"[0]", "ENDTIME", str(start), 
                    "UT at observation end", delete-, >& tmplog)
                start = n * 45.
                gemhedit (in[i]//"[0]", "OBSTIME", str(start), 
                    "[s] Duration of observation", delete-, >>& tmplog)
                glogprint (l_logfile, "mcheckheader", "task", type="file",
                    str=tmplog, verbose=l_verbose)
                delete (tmplog, ver-, >& "dev$null")

                imgets (in[i]//"[0]", "OBSSTART", >& "dev$null")
                time = substr (imgets.value, 12, 19)
                k = fscanf (time,"%2d:%2d:%2d", hr, min, sec)
                duration = n * 45
                glogprint (l_logfile, "mcheckheader", "engineering",
                    type="string", str="              duration = "//\
                    str(duration), verbose=l_verbose)
                k = duration / 3600
                duration = duration - k*3600
                l = duration / 60
                duration = duration - l*60
                glogprint (l_logfile, "mcheckheader", "engineering",
                    type="string", str="              start (HMS) = \
                    "//str(hr)//":"//str(min)//":"//str(sec), verbose=l_verbose)
                glogprint (l_logfile, "mcheckheader", "engineering",
                    type="string", str="              duration (HMS) = \
                    "//str(k)//":"//str(l)//":"//str(duration),
                    verbose=l_verbose)
                sec = sec + duration
                if (sec > 60) {
                    m = sec / 60
                    sec = sec - 60*m
                    l = l + m
                }
                min = min+l
                if (min > 60) {
                    m = min / 60
                    min = min - m*60
                    k = k + m
                }
                hr = hr + k
                # I am not checking for hr crossing 24:00:00; changing the date 
                # would be too complicated.  For science observations we should 
                # never have to worry about this.
                #
                
                newtime = substr (imgets.value, 1, 11)//\
                    str(hr)//":"//str(min)//":"//str(sec)//"Z" 
                             
                glogprint (l_logfile, "mcheckheader", "engineering",
                    type="string", str="              time="//str(time),
                    verbose=l_verbose)
                glogprint (l_logfile, "mcheckheader", "engineering",
                    type="string", str="              finish (HMS) = "//\
                    str(hr)//":"//str(min)//":"//str(sec), verbose=l_verbose)    
                glogprint (l_logfile, "mcheckheader", "engineering",
                    type="string", str="              newtime= "//str(newtime),
                    verbose=l_verbose)
                
                gemhedit (in[i]//"[0]", "OBSEND", newtime, "obsEnd", 
                    delete-, >& tmplog)
                gemhedit (in[i]//"[0]", "OBSCLASS", "unknown", "Observe class",
                    delete-, >>& tmplog)
                glogprint (l_logfile, "mcheckheader", "task", type="file",
                    str=tmplog, verbose=l_verbose)
                delete (tmplog, ver-, >& "dev$null")
            }
        }
        gemdate ()
        glogprint (l_logfile, "mcheckheader", "visual", type="visual",
            vistype="empty", verbose=l_verbose)
        glogprint (l_logfile, "mcheckheader", "status", type="string",
            str="Finished with image number "//i//" at ["//gemdate.outdate//"]",
            verbose=l_verbose)
        glogprint (l_logfile, "mcheckheader", "visual", type="visual",
            vistype="shortdash", verbose=l_verbose)
	
        # jump to here if there is a problem

nextimage:
        i = i+1
    }
clean:
    scanfile = ""
    delete (tmpinimg, ver-, >& "dev$null")

    if (status == 0) {
        glogclose (l_logfile, "mcheckheader", fl_success+, verbose=l_verbose)
    } else {
        glogclose (l_logfile, "mcheckheader", fl_success-, verbose=l_verbose)
    }	

exitnow:
    ;

end

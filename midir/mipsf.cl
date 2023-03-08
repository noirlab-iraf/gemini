# Copyright(c) 2007 Association of Universities for Research in Astronomy, Inc.

procedure mipsf (inimage)

char    inimage         {prompt="Michelle or T-ReCS image name"}
char    rawpath         {"", prompt="Path for raw images"}
bool    fl_stacked      {yes, prompt="Stacked Michelle or T-ReCS image"}
bool    fl_display      {yes, prompt="Current displayed image"}
char    logfile         {"", prompt="Logfile name"}
bool    verbose         {yes, prompt="Verbose"}
int     status          {0, prompt="Exit status"}
struct  *scanfile       {"", prompt="Internal use only"}

begin
    char    l_inimage = ""
    char    l_rawpath = ""
    char    l_logfile = ""
    bool    l_fl_stacked, l_fl_display, l_verbose

    char    paramstr, lastchar, errmsg
    char    inimg, instrument, filter1, filter2, filtername
    char    tmpinimg, tmplog, tmpimage
    char    tmpstring
    real    airmass
    real    xcoo, ycoo, r, mag, flux, sky, peak, ecc, pa
    real    beta, encl, mof, direct, gauss
    real    fwhm, arcfwhm, normpeak, strehl
    int     junk, nfields, star, line
    struct  scanline

    cache ("gloginit", "gemextn", "mireduce", "fparse")
    
    junk = fscan (inimage, l_inimage)
    junk = fscan (rawpath, l_rawpath)
    l_fl_stacked = fl_stacked
    l_fl_display = fl_display
    junk = fscan (logfile, l_logfile)
    l_verbose = verbose
    
    # Initialize
    status = 0
    
    # Create temp file names
    tmpinimg = mktemp ("tmpinimg")
    tmplog = mktemp ("tmplog")
    
    # Create the list of parameter/value pairs.  One pair per line.
    # All lines combined into one string.  Line delimiter is '\n'.
    paramstr  = "inimage        = "//inimage.p_value//"\n"
    paramstr += "rawpath        = "//rawpath.p_value//"\n"
    paramstr += "fl_stacked     = "//fl_stacked.p_value//"\n"
    paramstr += "fl_display     = "//fl_display.p_value//"\n"
    paramstr += "logfile        = "//logfile.p_value//"\n"
    paramstr += "verbose        = "//verbose.p_value
    
    # Assign a logfile name if not given.  Open logfile and start log.
    # Write parameter/value pairs ("paramstr") to log.
    gloginit (l_logfile, "mipsf", "midir", paramstr, fl_append+,
        verbose=l_verbose)
    if (gloginit.status != 0) {
        status = gloginit.status
        goto exitnow
    }
    l_logfile = gloginit.logfile.p_value

    # Add the trailing slash to rawpath, if missing
    lastchar = substr (l_rawpath, strlen(l_rawpath), strlen(l_rawpath))
    if ((l_rawpath != "") && (lastchar != "$") && (lastchar != "/"))
        l_rawpath = l_rawpath//"/"

    gemextn (l_inimage, check="exists,mef", process="none", index="",
        extname="", extversion="", ikparam="", omit="kernel,exten",
        replace="^%%"//l_rawpath//"%", outfile=tmpinimg, logfile=l_logfile,
        verbose=l_verbose)
    if (gemextn.count == 1) {
        fields (tmpinimg, "1", lines="1", quit=no, print=no) | scan (inimg)
    } else {
        errmsg = "Image not found, not MEF, or more than one image given"
        status = 121
        glogprint (l_logfile, "mipsf", "status", type="error", errno=status,
            str=errmsg, verbose=yes)
        goto clean
    }
    delete (tmpinimg, ver-, >& "dev$null")

    if (l_fl_stacked == no) {
        # mireduce does not handle inimg name with a path correctly.
        # I (K.Labrie) don't have time to fix that now, so I will tweak
        # mipsf to send info to mireduce that will actually work.
        
        fparse (inimg, verbose=no)
        # input to mireduce -> fparse.root
        # rawpath for mireduce -> fparse.directory
        
        # Okay, now we are ready to call mireduce
        glogprint (l_logfile, "mipsf", "task", type="string",
            str="Reducing the Michelle/T-ReCS image, please wait....",
            verbose=l_verbose)

        glogprint (l_logfile, "mipsf", "status", type="fork", fork="forward",
            child="mireduce", verbose=l_verbose)
        tmpimage = mktemp ("tmpimage")
        mireduce (fparse.root, outimage=tmpimage, rawpath=fparse.directory, outpref="j",
            fl_background-, fl_view-, fl_flat-, fl_display-, verbose-,
            logfile=l_logfile)
        status = mireduce.status
        glogprint (l_logfile, "mipsf", "status", type="fork", fork="backward",
            child="mireduce", verbose=l_verbose)
        glogprint (l_logfile, "mipsf", "visual", type="visual",
            vistype="empty", verbose=l_verbose)

        if (status != 0) {
            errmsg = "MIREDUCE failed for image "//inimg//".  Exiting."
            glogprint (l_logfile, "mipsf", "status", type="error",
                errno=status, str=errmsg, verbose=yes)
            goto clean
        }
    } else {
        tmpimage = inimg
    }

    # Determine which instrument this data is from
    hselect (tmpimage//"[0]", "INSTRUME", yes) | scan (instrument)
    hselect (tmpimage//"[0]", "AIRMASS", yes) | scan (airmass)
    
    # Get filter name from PHU
    if (strlwr(instrument) == "michelle") {
        hselect (tmpimage//"[0]", "FILTER", yes) | scan (filter1)
        filtername = filter1        
    } else if (instrument == "TReCS") {
        # Get filter name from PHU
        hselect (tmpimage//"[0]", "FILTER1", yes) | scan (filter1)
        hselect (tmpimage//"[0]", "FILTER2", yes) | scan (filter2)
        if (substr(filter1, 1, 4) == "Open")
            filtername = filter2
        else
            filtername = filter1
    } else {
        errmsg = "Data from unsupported instrument ("//instrument//")"
        status = 121
        glogprint (l_logfile, "mipsf", "status", type="error", errno=status,
            str=errmsg, verbose=yes)
        goto clean
    }

    if (l_fl_display) {
        display (tmpimage//"[1]", 1, zrange+, zscale-, nsample=76400,
            ztrans="linear", >& "dev$null")
    }

    glogprint (l_logfile, "mipsf", "status", type="fork", fork="forward",
        child="imexam", verbose=l_verbose)
    print ("")
    print ("    Entering IMEXAM: please use 'r' on the star or stars to ")
    print ("    get the radial profile, then q to quit.")
    print ("")

    imexam (tmpimage//"[1]", 1, logfile=tmplog, keeplog+, > "dev$null")
    glogprint (l_logfile, "mipsf", "status", type="fork", fork="backward",
        child="imexam", verbose=l_verbose)

    scanfile = tmplog
    line = 0
    star = 0
    while (fscan (scanfile, scanline) != EOF) {
        line = line + 1
        if (substr(scanline, 1, 1) == '#') {
            next
        }
        print (scanline) | count ("STDIN") | scan (junk, nfields, junk)
        if ((nfields != 14) && (nfields != 15)) {
            next
        }

        star = star + 1
        print (scanline) | fields ("STDIN", "1,2", lines="1", quit=no,
            print=no) | scan (xcoo, ycoo)
        print (scanline) | fields ("STDIN", "5,6,7,8,9", lines="1", quit=no,
            print=no) | scan (r, mag, flux, sky, peak)
        print (scanline) | fields ("STDIN", "10,11", lines="1", quit=no,
            print=no) | scan (ecc, pa)
            
        if (ecc < 0.005) pa = 0.
        
        tmpstring = str(r)//" "//str(mag)//" "//str(flux)//" "//\
                    str(sky)//" "//str(peak)//" "//str(ecc)//" "//\
                    str(pa)
                    
        if (nfields == 15) {
            print (scanline) | fields ("STDIN", "12", lines="1", quit=no,
                print=no) | scan (beta)
            print (scanline) | fields ("STDIN", "13,14,15", lines="1",
                quit=no, print=no) | scan (encl, mof, direct)
            tmpstring = tmpstring//" "//str(beta)//" "//str(encl)//" "//\
                        str(mof)//" "//str(direct)
        } else if (nfields == 14) {
            print (scanline) | fields ("STDIN", "12,13,14", lines="1",
                quit=no, print=no) | scan (encl, gauss, direct)
            tmpstring = tmpstring//" "//str(encl)//" "//str(gauss)//" "//\
                        str(direct)
        } else {
            errmsg = "Problem with results from IMEXAM"
            status = 99
            glogprint (l_logfile, "mipsf", "status", type="error",
                errno=status, str=errmsg, verbose=yes)
            goto clean
        }
        
        glogprint (l_logfile, "mipsf", "engineering", type="string",
            str=tmpstring, verbose=l_verbose)
        
        # Set the FWHM from the IMEXAM output
        fwhm = direct
        
        glogprint (l_logfile, "mipsf", "visual", type="visual",
            vistype="empty", verbose=l_verbose)
        printf ("Position of star #%d: %.2f %.2f\n", star, xcoo, ycoo) |\
            scan (scanline)
        glogprint (l_logfile, "mipsf", "science", type="string",
            str=scanline, verbose=l_verbose)
        printf ("\tFWHM value: %.2f pixels for filter %s.\n",
            fwhm, filtername) | scan (scanline)
        glogprint (l_logfile, "mipsf", "science", type="string",
            str=scanline, verbose=l_verbose)
        
        if (strlwr(instrument) == "michelle") arcfwhm = 0.1005 * fwhm
        else                                  arcfwhm = 0.08964 * fwhm # TReCS
        
        printf ("\tThis is nominally %.3f arc seconds.\n", arcfwhm) |\
            scan (scanline)
        glogprint (l_logfile, "mipsf", "science", type="string",
            str=scanline, verbose=l_verbose)
        
#        arcfwhm = arcfwhm / (airmass**0.6)
        printf ("\tAirmass is: %.3f.\n", airmass) | scan (scanline)
        glogprint (l_logfile, "mipsf", "science", type="string",
            str=scanline, verbose=l_verbose)
        
#        printf ("\tCorrected to zenith the FWHM would be %.3f arc seconds.\n",
#            arcfwhm) | scan (scanline)
#
# We should not correct the Michelle/T-ReCS FWHM to zenith the same way as 
# in the optical.  I have commented this part of the script out.  [Kevin Volk
# 3 May 2008].
# 
       glogprint (l_logfile, "mipsf", "science", type="string",
            str=scanline, verbose=l_verbose)     
        glogprint (l_logfile, "mipsf", "visual", type="visual",
            vistype="empty", verbose=l_verbose)
        
        printf ("\tTotal counts: %.4e\n", flux) | scan (scanline)
        glogprint (l_logfile, "mipsf", "science", type="string",
            str=scanline, verbose=l_verbose) 
        
        normpeak = peak / flux
        printf ("\tNormalized peak: %.5f\n", normpeak) | scan (scanline)
        glogprint (l_logfile, "mipsf", "science", type="string",
            str=scanline, verbose=l_verbose)
                
        strehl = 0.
        
        if (strlwr(instrument) == "michelle") {            
            if ((filtername == 'I79B10') || (filtername == 'IP79B10')) {
                strehl = normpeak / 0.1979890
            } else if ((filtername == 'I88B10') || (filtername == 'IP88B10')) {
                strehl = normpeak / 0.1517450
            } else if ((filtername == 'I97B10') || (filtername == 'IP97B10')) {
                strehl = normpeak / 0.1250020
            } else if ((filtername == 'I103B10') || (filtername == 'IP103B10')) {
                strehl = normpeak / 0.1109290
            } else if ((filtername == 'I116B9') || (filtername == 'IP116B9')) {
                strehl = normpeak / 0.0875688
            } else if (filtername == 'I125B9' || filtername == 'IP125B9') {
                strehl = normpeak / 0.0754795
            } else if ((filtername == 'I112B21') || (filtername == 'IP112B21')) {
                strehl = normpeak / 0.0938997
            } else if ((filtername == 'I185B9') || (filtername == 'IP185B9')) {
                strehl = normpeak / 0.0354149
            } else {
                errmsg = "Unsupported filter name ("//filtername//")"
                status = 123
                glogprint (l_logfile, "mipsf", "status", type="warning",
                    errno=status, str=errmsg, verbose=l_verbose)
            }
        } else {    # TReCS
            if (filtername == 'Si1-7.9um') {
                strehl = normpeak / 0.1300672
            } else if (filtername == 'Si2-8.8um') {
                strehl = normpeak / 0.1082258
            } else if (filtername == 'Si3-9.7um') {
                strehl = normpeak / 0.0872771
            } else if (filtername == 'Si4-10.4um') {
                strehl = normpeak / 0.0760649
            } else if (filtername == 'Si5-11.7um') {
                strehl = normpeak / 0.0601185
            } else if (filtername == 'Si6-12.3um') {
                strehl = normpeak / 0.0526666
            } else if (filtername == 'Qa-18.3um') {
                strehl = normpeak / 0.02253
                # The "theoretical" value from Francois is 0.0202076 but 
                # that gives Strehls greater than 1 so it is not correct.  
                # I (K.Volk) just assumed a maximum strehl of about 0.93 
                # to get this above value
            } else {
                errmsg = "Unsupported filter name ("//filtername//")"
                status = 123
                glogprint (l_logfile, "mipsf", "status", type="warning",
                    errno=status, str=errmsg, verbose=l_verbose)
            }
        }
        
        glogprint (l_logfile, "mipsf", "visual", type="visual",
            vistype="empty", verbose=l_verbose)
        
        if (strehl > 0.) {
            printf ("\t** Approximate Strehl: %.3f **\n", strehl) | \
                scan (scanline)
            glogprint (l_logfile, "mipsf", "science", type="string",
                str=scanline, verbose=l_verbose)
            glogprint (l_logfile, "mipsf", "visual", type="visual",
                vistype="empty", verbose=l_verbose)
        }
        
    }


clean:
    scanfile = ""
    delete (tmpinimg, ver-, >& "dev$null")
    delete (tmplog, ver-, >& "dev$null")

    if (fl_stacked == no)
        imdelete (tmpimage, ver-, >& "dev$null")

    if (status == 0)
        glogclose (l_logfile, "mipsf", fl_success+, verbose=l_verbose)
    else
        glogclose (l_logfile, "mipsf", fl_success-, verbose=l_verbose)

exitnow:
    ;

end

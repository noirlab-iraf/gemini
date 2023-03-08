# Copyright(c) 2002-2011 Association of Universities for Research in Astronomy, Inc.

procedure gemoffsetlist (infiles, reffile, distance, age, targetlist, offsetlist) 

# Version Aug 26, 2002  MT,JJ,IJ v1.4 release
#         Aug 12, 2003  KL IRAF2.12 - hedit: addonly-
#	  Dec 2, 2003   AC Tidied, distance and age tested together

char    infiles     {prompt = "Input images"}
char    reffile     {prompt = "Reference file"}

real    distance    {INDEF, prompt = "Radius (arcsec) from ref for target list"}
bool    fl_nearer   {yes, prompt = "Select files nearer than the distance limit?"}
int     direction   {3, min=1, max=3, prompt = "Axes to measure direction along (1=x/ra, 2=y/dec, 3=radial)"}
real    age         {INDEF, prompt = "Time difference (s) for target list"}
bool    fl_younger  {yes, prompt = "Select files younger than the age limit?"}
bool    fl_noref    {no, prompt = "Exclude the reference file from infiles?"}
char    wcs_source  {"phu", prompt = "Source of WCS?", enum = "none|phu|direct|inherit"}

char    targetlist  {prompt = "Output for images near to reference"}
char    offsetlist  {prompt = "Output for images far from reference"}

char    key_xoff    {"XOFFSET", prompt = "PHU keyword for instrument x offset (in arcsec)"}
char    key_yoff    {"YOFFSET", prompt = "PHU keyword for instrument y offset (in arcsec)"}
char    key_date    {"DATE-OBS", prompt = "Header keyword for the date"}
char    key_time    {"UT", prompt = "Header keyword for the time"}

char    logfile     {"", prompt = "Logfile"}
bool    verbose     {no, prompt = "Verbose output?"}
bool    force       {no, prompt = "Force use with earlier IRAF versions?"}

int     status      {0, prompt = "Exit status (0=good)"}
int     count       {0, prompt = "Number of files copied to targetlist"}

struct	*scanfile   {"", prompt = "Internal use only"} 

begin
    char    l_infiles = ""
    char    l_reffile = ""
    real    l_distance, l_age
    int     l_direction
    bool    l_fl_younger, l_fl_nearer, l_fl_noref
    char    l_wcs_source = ""
    char    l_targetlist = ""
    char    l_offsetlist = ""
    char    l_key_xoff = ""
    char    l_key_yoff = ""
    char    l_key_date = ""
    char    l_key_time = ""
    char    l_logfile = ""
    bool    l_verbose
    bool    l_force

    int     junk, nbad, total, countin, countout, idx
    struct  sdate
    char    tmpinfiles, tmprefandin, tmpoffsets, tmpref
    char    tmpcoordin
    char    img, ref, reason, imgwcs, refwcs, msg, values, dtlog
    char    inherit, check, index, noreftxt, datestr, timestr, keystr
    bool    use_wcs, first, time_ok, dist_ok, debug, use_inherit, use_phu
    bool    versionok
    real    xoff, yoff, toff, pixscale, xoff_ref, yoff_ref, toff_ref
    real    dx, dy, dt, r, raref, decref, crval1, crval2, cd11, cd22, fjunk

    cache ("gemextn", "keypar")

    debug = no
    status = 1
    count = 0
    total = 0
    versionok = yes

    junk = fscan (  infiles, l_infiles)
    junk = fscan (  reffile, l_reffile)
    l_distance =     distance
    l_direction =   direction
    l_age =         age
    l_fl_younger =  fl_younger
    l_fl_nearer =   fl_nearer
    l_fl_noref =    fl_noref
    junk = fscan (  wcs_source, l_wcs_source)
    junk = fscan (  targetlist, l_targetlist)
    junk = fscan (  offsetlist, l_offsetlist)
    junk = fscan (  key_xoff, l_key_xoff)
    junk = fscan (  key_yoff, l_key_yoff)
    junk = fscan (  key_date, l_key_date)
    junk = fscan (  key_time, l_key_time)
    junk = fscan (  logfile, l_logfile)
    l_verbose =     verbose
    l_force =       force

    use_wcs = l_wcs_source != "none"
    use_phu = l_wcs_source == "phu"
	use_inherit = l_wcs_source == "inherit"   
    
    tmpinfiles = mktemp ("tmpinfiles")
    tmprefandin = mktemp ("tmprefandin")
    tmpoffsets = mktemp ("tmpoffsets")
    tmpref = mktemp ("tmpref")
    tmpcoordin = mktemp ("tmpcoordin")

    if ("" == l_targetlist) l_targetlist = "STDOUT"
    if ("" == l_offsetlist) l_offsetlist = "dev$null"


    # Start logging to file

    if (l_logfile == "") {
        l_logfile = "gemoffsetlist.log"
        printlog ("WARNING - GEMOFFSETLIST: logfile is empty.", \
            l_logfile, verbose+) 
        printlog ("          Using default file " // l_logfile // ".", \
            l_logfile, verbose+) 
    }
    date | scan (sdate) 
    printlog ("----------------------------------------------------\
        --------------------------", l_logfile, l_verbose) 
    printlog ("GEMOFFSETLIST -- " // sdate, l_logfile, l_verbose) 
    printlog (" ", l_logfile, l_verbose) 


    if (defpar ("release")) versionok = (no == (release <= "2.12.2"))
    else versionok = no

    if (no == versionok && no == force) {
        printlog ("ERROR - GEMOFFSETLIST: This task will not work with \
            IRAF versions earlier", l_logfile, verbose+)
        printlog ("                       than 2.12.2.", \
            l_logfile, verbose+)
        printlog ("                       To enable this task anyway, \
            use force+.", l_logfile, verbose+)
        goto clean
    }
    if (l_force) {
        printlog ("WARNING - GEMOFFSETLIST: Version test ignored", \
            l_logfile, verbose+)
    }

    # Check that target and offset lists do not exist
    if (debug) print ("checking destination files")
    if (access (l_targetlist) && "dev$null" != l_targetlist \
        && strstr ("STD", strupr (l_targetlist)) != 1) {
        printlog ("ERROR - GEMOFFSETLIST: " // l_targetlist \
            // " already exists.", l_logfile, verbose+)
        goto clean
    }	
    if (access (l_offsetlist) && "dev$null" != l_offsetlist \
        && strstr ("STD", strupr (l_offsetlist)) != 1) {
        printlog ("ERROR - GEMOFFSETLIST: " // l_offsetlist \
            // " already exists.", l_logfile, verbose+) 
        goto clean
    }


    # Check that input files exist
    if (debug) print ("checking input")
    gemextn (l_infiles, check="exists,image", process="none", \
        index="", extname="", extversion="", ikparams="", \
        omit="extension", replace="", outfile=tmpinfiles, \
        logfile="", glogpars="", verbose=l_verbose)
    if (gemextn.fail_count != 0 || gemextn.count == 0) {
        printlog ("ERROR - GEMOFFSETLIST: Missing input files", \
            l_logfile, verbose+) 
        goto clean
    }
	

    # Check for single reference file
    if (debug) print ("checking ref")
    gemextn (l_reffile, check="exists,image", process="none",
        index="", extname="", extversion="", ikparams="",
        omit="extension", replace="", outfile=tmpref, logfile="",
        glogpars="", verbose=l_verbose)
    if (0 != gemextn.fail_count || 1 != gemextn.count) {
        printlog ("ERROR - GEMOFFSETLIST: Missing reference file", \
            l_logfile, verbose+) 
        goto clean
    }
    scanfile = tmpref
    junk = fscan (scanfile, ref)
    if (debug) print ("ref: " // ref)


    # Combine in a single file, reference first
    files ("@" // tmpref // ",@" // tmpinfiles, sort-, > tmprefandin)
	

    # Check that distance info available, if needed

    if (use_wcs) {
        if (use_inherit) {
            inherit = "inherit"
            check = "exists,image"
        } else {
            inherit = ""
            check = "exists"
        }
        if (use_phu)    index = "0"
        else            index = "1-"

        if (debug) {
            print ("inherit: " // inherit)
            print ("check: " // check)
            print ("index: " // index)
        }
    }


    count (tmpinfiles) | scan (total)


    if (INDEF != l_distance) {
        if (debug) print ("checking distance info")

        first = yes
        nbad = 0 
        scanfile = tmprefandin
        while (fscan (scanfile, img) != EOF) {

            if (use_wcs) {

                # get appropriate mef
                if (debug) print ("expanding " // img)
                gemextn (img, proc="expand", index=index, extname="", \
                    extver="", ikparams=inherit, check=check, \
                    omit="name,version", replace="", outfile="STDOUT", \
                    logfile="dev$null", glogpars="", verbose-) | \
                    scan (imgwcs)

                if (0 == gemextn.count) {

                    printlog ("WARNING - GEMOFFSETLIST: Couldn't find the \
                        wcs extension in " // img, l_logfile, verbose+) 
                    nbad += 1

                } else {

                    if (debug) print ("check wcs: " // imgwcs)

                    crval1 = 0.0; crval2 = 0.0
                    hselect (imgwcs, "CRVAL1,CRVAL2", yes) \
                        | scan (crval1, crval2)
                    if (2 != nscan ()) {
                        nbad += 1
                        if (debug) print ("CRVAL1/2 missing")
                    } else if (0.0 == crval1 || 0.0 == crval2) {
                        nbad += 1
                        if (debug) print ("CRVAL1/2 zero")
                    } else {
                        if (debug) print ("CRVAL ok")
                        raref = crval1
                        decref = crval2
                    }

                }
            }
            first = no
        }

        if (nbad != 0) {
            printlog ("WARNING - GEMOFFSETLIST: " // nbad \
                // " files do not have specified WCS information.", \
                l_logfile, verbose+) 
            printlog ("                         Using stored header \
                offsets instead.", l_logfile, verbose+) 
            use_wcs = no
        }


        # Check again, without WCS
        if (no == use_wcs) {

            nbad = 0
            first = yes
            scanfile = tmprefandin

            while (fscan (scanfile, img) != EOF) {
                if (debug) print ("check: " // img)
                hselect (img // "[0]", l_key_xoff, yes) | scan (fjunk)
                if (1 != nscan ()) {
                    nbad += 1
                } else {
                    if (first) raref = fjunk
                    hselect (img // "[0]", l_key_yoff, yes) | scan (fjunk)
                    if (1 != nscan ()) {
                        nbad += 1
                    } else {
                        if (first) decref = fjunk
                    }
                } 
                if (debug && nbad > 0) print ("nbad: " // nbad)
                first = no
            }

            if (nbad != 0) {
                printlog ("ERROR - GEMOFFSETLIST: " // nbad \
                    // " files do not have instrument offset information.",
                    l_logfile, verbose+) 
                goto clean
            }
        }

    } else {

        # No distance calc
        raref = 0; decref = 0
    }


    # Check the time info, if needed

    if (INDEF != l_age) {
        if (debug) {
            print ("checking date info: " // l_key_date)
            print ("checking time info: " // l_key_time)
        }

        nbad = 0
        scanfile = tmprefandin
        while (fscan (scanfile, img) != EOF) {
            hselect (img//"[0]", l_key_date, yes) | scan (keystr)
            if ( nscan () != 1 ) {  #key_date not found, we're in trouble
                nbad += 1
            } else {                #key_date found, check for time now
                idx = stridx ("T", keystr)
                if (idx == 0) {   # key_date does not contain time
                    hselect (img//"[0]", l_key_time, yes) | scan (keystr)
                    if (nscan() != 1)  nbad += 1    #key_time not found
                }
            }
        }	
        if (nbad != 0) {
            printlog ("ERROR - GEMOFFSETLIST: " // nbad \
                // " files do not have date/time in header.", \
                l_logfile, verbose+) 
            goto clean
        }	
    }


    # Get pixel scale, if required
    if (INDEF != l_distance && use_wcs) {
        if (debug) print ("getting pixel scale")
        hselect (ref // "[0]", "PIXSCALE", yes) | scan (pixscale)
        if (1 == nscan()) {
            if (debug) print ("from header " // pixscale)
        } else {
            hselect (ref // "[0]", "CD1_1", yes) | scan (cd11)
            hselect (ref // "[0]", "CD2_2", yes) | scan (cd22)
            pixscale = 3600 * sqrt (cd11*cd11 + cd22 * cd22)
        }
    }


    if (debug) print ("verification done\n")
    if (no == use_wcs) {
        printlog ("GEMOFFSETLIST: Using " // l_key_xoff // " and " \
            // l_key_yoff // " from header.", l_logfile, verbose+)
    } else {
        printf ("%12.8f %12.8f\n", raref, decref, > tmpcoordin)
        if (debug) type (tmpcoordin)
    }

    # Process each entry, generating a file with four columns:
    # File name, x and y offsets and time offset

    scanfile = tmprefandin
    while (fscan (scanfile, img) != EOF) {
        if (debug) print ("processing: " // img)

        xoff = 0.0
        yoff = 0.0
        if (INDEF != l_distance) {

            if (use_wcs) {

                # get first image in the mef
                gemextn (img, proc="expand", index=index, extname="", \
                    extver="", ikparams=inherit, check=check, \
                    omit="name,version", replace="", outfile="STDOUT", \
                    logfile="dev$null", glogpars="", verbose-) | \
                scan (imgwcs)
                if (0 == gemextn.count) imgwcs = img
                if (debug) print ("wcs from: " // imgwcs)

                if (debug) {
                    print ("wcsctran:")
                    wcsctran (tmpcoordin, "STDOUT", imgwcs, \
                        "world", "logical", columns = "1 2", \
                        units = "", formats = "", min_sigdigit = 7, \
                        verbose-)
                }
                wcsctran (tmpcoordin, "STDOUT", imgwcs, \
                    "world", "logical", columns = "1 2", \
                    units = "", formats = "", min_sigdigit = 7, \
                    verbose-) | scan (xoff, yoff)

                xoff = (-xoff) * pixscale
                yoff = (-yoff) * pixscale

            } else {

                hselect (img // "[0]", l_key_xoff, yes) | scan (xoff)
                hselect (img // "[0]", l_key_yoff, yes) | scan (yoff)

            }
        }

        toff = 0.0
        if (INDEF != l_age) {
            datestr="" ; timestr=""
#            hselect (img // "[0]", l_key_date, yes) | scan (datestr)
#            hselect (img // "[0]", l_key_time, yes) | scan (timestr)
            keypar(img // "[0]", l_key_date)
            datestr = keypar.value
            keypar(img // "[0]", l_key_time)
            timestr = keypar.value
            cnvtsec (datestr, timestr) | scan (toff)
            if (cnvtsec.status != 0) {
                status = cnvtsec.status
                printlog ("GEMOFFSETLIST ERROR Unable to retrieve the \
                    date/time of the observation.", l_logfile, verbose+)
                printlog ("GEMOFFSETLIST ERROR Check header of image "//\
                    img, l_logfile, verbose+)
                goto clean
            }
        }

        printf("%s %.2f %.2f %.10f\n", img, xoff, yoff, toff) | scan (line)
        if (debug) print (line)
        print (line, >> tmpoffsets)
    }


    # Select results
    if (debug) print ("selecting")
    first = yes
    scanfile = tmpoffsets
    countin = 0
    countout = 0
    while (fscan (scanfile, img, xoff, yoff, toff) != EOF) {

        if (debug) print (img)
        reason = ""

        if (first) {

            # Store/drop reference values
            xoff_ref = xoff
            yoff_ref = yoff
            toff_ref = toff
            first = no

            if (l_verbose) {
                msg = "GEMOFFSETLIST: matching " // img
                if (INDEF != l_age) {
                    if (l_fl_younger) msg = msg // "; age <= " // l_age
                    else              msg = msg // "; age > " // l_age
                }
                if (INDEF != l_distance) {
                    if (l_fl_nearer) msg = msg // "; dist <= " // l_distance
                    else             msg = msg // "; dist > " // l_distance
                }
                printlog (msg, l_logfile, l_verbose)
                printlog ("               ref age: " // toff \
                    // "; ref x: " // xoff // "; ref y: " // yoff, \
                    l_logfile, l_verbose)
            }
            dx = 0
            dy = 0

        } else if (no == l_fl_noref || (ref != img)) {

            dx = xoff - xoff_ref
            dy = yoff - yoff_ref
            if (1 == l_direction) r = abs (dx)
            else if (2 == l_direction) r = abs (dy)
            else if (3 == l_direction) r = sqrt (dx*dx+dy*dy) 
            dt = abs (toff - toff_ref)

            if (debug) print (dx // ", " // dy // ", " // dt)

            time_ok = yes
            if (INDEF != l_age) {
                if (l_fl_younger) time_ok = dt <= l_age
                else              time_ok = dt > l_age
                if (time_ok) reason = "time ok"
                else	 reason = "time not ok"
            }

            dist_ok = yes
            if (INDEF != l_distance) {
                if (l_fl_nearer) dist_ok = r <= l_distance
                else	     dist_ok = r > l_distance
                if ("" != reason) reason = reason // ", "
                if (dist_ok) {
                    reason = reason // "distance ok"
                } else {
                    reason = reason // "distance not ok"
                }
            }

            if (dist_ok && time_ok) {
                print (img, >> l_targetlist)
                countin = countin + 1

                if (l_verbose) {
                    msg = "GEMOFFSETLIST: including " // img // ": " // reason
                    printlog (msg, l_logfile, l_verbose)
                }

            } else {
                print (img, >> l_offsetlist)
                countout = countout + 1

                if (l_verbose) {
                    msg = "GEMOFFSETLIST: excluding " // img // ": " // reason
                    printlog (msg, l_logfile, l_verbose)
                }

            }

            printf ("%.1f\n", dt) | scan (dtlog)
            printlog ("         delta age: " // dtlog \
                // " (s); delta x: " // dx // " ; delta y: " // dy \
                // "; distance: " // r // " (arcsec)", \
                l_logfile, l_verbose)

        } else {

            # Have ref and fl_noref+
            dx = 0
            dy = 0

        }

        # These header values are used as a return value for 
        # some tasks (nssdist) so that spatial offsets are
        # still available after wavelength cal.

        # There was some logic here that made this
        # conditional on the keys NOT being GXOFF etc.  This
        # was removed to support gnirs, but may have introduced
        # some subtle bug with niri (the whole thing is an
        # ugly hack...)

        if (use_wcs && "" != l_key_xoff && "" != l_key_yoff) {
            dx = (int (dx*100) ) /100.
            dy = (int (dy*100) ) /100.
            gemhedit (img // "[0]", "GOFFREF", ref, \
                "GEMOFFSETLIST spatial reference image") 
            gemhedit (img // "[0]", l_key_xoff, dx, \
                "x offset from reference image (arcsec)") 
            gemhedit (img // "[0]", l_key_yoff, dy, \
                "y offset from reference image (arcsec)") 
        }
    }

    if (debug)
        print ("divided " // total // " into " // countin // "/" // countout)
    noreftxt = ""
    if (l_fl_noref) noreftxt = " (reference excluded from list)"
    printlog ("GEMOFFSETLIST: Divided " // total // " into " // countin \
        // "/" // countout // noreftxt, l_logfile, l_verbose)

    status = 0	# Success
    count = countin

clean:
    scanfile = ""
    delete (tmpinfiles // "," // tmprefandin // "," // tmpoffsets \
        // "," // tmpref // "," // tmpcoordin, verify-, >>& "dev$null") 
    if (0 == status)
        printlog ("GEMOFFSETLIST Exit status - SUCCESS", l_logfile, l_verbose)
    else
        printlog ("GEMOFFSETLIST Exit status - FAILURE", l_logfile, verbose+)
    printlog ("-------------------------------------------------\
        -------------------------------", l_logfile, l_verbose)

end

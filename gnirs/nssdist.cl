# Copyright(c) 2000-2012 Association of Universities for Research in Astronomy, Inc.

procedure nssdist (inimages) 

# Measure the S-distortion of GNIRS/NIRI spectra
# This task establishes the calibration using stepped images of a star 
# along a slit, or perhaps data from a "pinhole" slit
#
# Version: Sept 20, 2002 JJ v1.4 release
#          Aug 19, 2003  KL IRAF2.12 - new parameter, addonly, in hedit
#          Oct 29, 2003  KL moved from niri to gnirs package

char    inimages    {prompt = "Input GNIRS/NIRI spectra (pinhole or star spectra)"}
char    outsuffix   {"_sdist", prompt = "Output spectrum suffix (when combining multiple inputs)"}
real    pixscale    {1., prompt = "Pixel scale in arcsec/pixel"}
int     dispaxis    {1, min = 1, max = 2, prompt = "Dispersion axis"}
char    database    {"", prompt = "Directory for files containing feature data"}
real    firstcoord  {0., min = 0., prompt = "spatial coord of star in first input image"}
char    coordlist   {"", prompt = "Co-ordinate list for pinhole spectra"}
char    aptable     {"gnirs$data/apertures.fits", prompt = "Table of aperture data"}
bool    fl_inter    {no, prompt = "Examine identifications interactively"}
bool    fl_dbwrite  {yes, prompt = "Write results to database"}
char    section     {"default", prompt = "Image section for running identify"}
int     nsum        {30, min = 1, prompt = "Number of lines or columns to sum"}
char    ftype       {"emission", min = "emission|absorption", prompt = "Feature type"}
real    fwidth      {10., min = 2, prompt = "Feature width in pixels"}
real    cradius     {10., min = 2, prompt = "Centering radius in pixels"}
real    threshold   {50., prompt = "Feature threshold for centering"}
real    minsep      {10., prompt = "Minimum pixel separation for features"}
real    match       {-6., prompt = "Coordinate list matching limit, <0 pixels, >0 user"}
char    function    {"chebyshev", min = "legendre|chebyshev|spline1|spline3", prompt = "Coordinate fitting function"}
int     order       {2, min = 1, prompt = "Order of coordinate fitting function"}
char    sample      {"", prompt = "Coordinate sample regions"}
int     niterate    {3, min = 0, prompt = "Rejection iterations"}
real    low_reject  {5., min = 0, prompt = "Lower rejection sigma"}
real    high_reject {5., min = 0, prompt = "Upper rejection sigma"}
real    grow        {0., min = 0, prompt = "Rejection growing radius"}
bool    refit       {yes, prompt = "Refit coordinate function when running reidentify"}
int     step        {10, prompt = "Steps in lines or columns for reidentification"}
bool    trace       {no, prompt = "Use fit from previous step rather than central aperture"}
int     nlost       {10, min = 0, prompt = "Maximum number of lost features"}
char    aiddebug    {"", prompt = "Debug parameter for aidpars"}
char    logfile     {"", prompt = "Logfile"}
bool    verbose     {yes, prompt = "Verbose"}
bool    debug       {no, prompt = "Very verbose"}
bool    force       {no, prompt = "Force use with earlier IRAF versions?"}

int     status      {0, prompt = "Exit status (0=Good)"}
struct  *scanfile   {"", prompt = "Internal use"}

begin
        
    char    l_inimages = ""
    char    l_outsuffix = ""
    real    l_pixscale
    int     l_dispaxis
    char    l_database = ""
    real    l_firstcoord
    char    l_coordlist = ""
    char    l_aptable = ""
    bool    l_fl_inter
    bool    l_fl_dbwrite
    char    l_section = ""
    int     l_nsum
    char    l_ftype = ""
    real    l_fwidth
    real    l_cradius
    real    l_threshold
    real    l_minsep
    real    l_match
    char    l_function = ""
    int     l_order
    char    l_sample = ""
    int     l_niterate
    real    l_low_reject
    real    l_high_reject
    real    l_grow
    bool    l_refit
    int     l_step
    bool    l_trace
    int     l_nlost
    char    l_aiddebug = ""
    char    l_logfile = ""
    bool    l_verbose
    bool    l_debug
    bool    l_force

    char    l_key_xoff = ""
    char    l_key_yoff = ""
    char    l_key_pixscale = ""
    char    l_key_dispaxis = ""
    char    l_sci_ext = ""
    char    l_var_ext = ""
    char    l_dq_ext = ""
    char    l_key_arrayid = ""
    char    l_key_prism = ""
    char    l_key_decker = ""
    char    l_key_fpmask = ""
    char    l_key_exptime = ""
    char    l_key_ron = ""
    char    l_key_gain = ""
    char    l_key_instrument = ""

    char    l_units, l_outspectra, l_prefix, l_intemp, l_outtemp, tmpinlist
    char    tmpphu, l_outsec, keyfound, s_inter, s_dbwrite, tmpname, destname
    char    l_temp, img, imgout, sec, badhdr, phu, tmpcoordlist, tmphead
    char    tmpshiftcoord, tmpfakedin, offkey, arrayid, prism, decker, fpmask
    char    aprow, search, secn, subsec, unused 
    bool    l_fl_trimmed, usetable, finalwarn, usecentre, l_flsetarray
    bool    l_fl_specsec, fl_newout, f2, intdbg 
    real    l_crval, l_cdelt, l_crpix, apshift, shift, coord, yoff
    int     ii, nlamps, nout, l_spataxis, n_lines, junk, nfiles, nbad, count
    int     nin, minfound, nver, version, nbegin, nindef
    struct  l_struct

    tmpinlist = mktemp ("tmpinl") 
    tmpphu = mktemp ("tmpphu") 
    tmphead = mktemp ("tmphead") 
    tmpcoordlist = mktemp ("tmpcoordlist") 

    # Initialize parameters
    intdbg = no
    status = 1
    fl_newout = no # new output when combining several files
    finalwarn = no
    usecentre = no
    f2 = no
    unused = "(unused)"

    # Cache parameter files
    cache ("niri", "keypar", "gemextn", "gemcombine", "gemdate") 

    # Set the local variables

    junk = fscan (inimages, l_inimages)
    junk = fscan (outsuffix, l_outsuffix)
    l_pixscale = pixscale
    l_dispaxis = dispaxis
    junk = fscan (database, l_database)
    l_firstcoord = firstcoord
    junk = fscan (coordlist, l_coordlist)
    junk = fscan (aptable, l_aptable)
    l_fl_inter = fl_inter
    l_fl_dbwrite = fl_dbwrite
    l_section = section # may include spaces, like "first line"
    l_nsum = nsum
    junk = fscan (ftype, l_ftype)
    l_fwidth = fwidth
    l_cradius = cradius
    l_threshold = threshold
    l_minsep = minsep
    l_match = match
    junk = fscan (function, l_function)
    l_order = order
    junk = fscan (sample, l_sample)
    l_niterate = niterate
    l_low_reject = low_reject
    l_high_reject = high_reject
    l_grow = grow
    l_refit = refit
    l_step = step
    l_trace = trace
    l_nlost = nlost
    junk = fscan (aiddebug, l_aiddebug)
    junk = fscan (logfile, l_logfile)
    l_verbose = verbose 
    l_debug = debug
    l_force = force

    if (l_verbose == no)
        l_aiddebug = ""

    while (strlen (l_section) > 0 && " " == substr (l_section, 1, 1)) {
        l_section = substr (l_section, 2, strlen (l_section))
    }
    if (intdbg) print ("section: '" // l_section // "'")

    badhdr = ""
    junk = fscan (  nsheaders.key_xoff, l_key_xoff)
    if ("" == l_key_xoff) badhdr = badhdr + " key_xoff"
    junk = fscan (  nsheaders.key_yoff, l_key_yoff)
    if ("" == l_key_yoff) badhdr = badhdr + " key_yoff"
    junk = fscan (  nsheaders.key_pixscale, l_key_pixscale)
    if ("" == l_key_pixscale) badhdr = badhdr + " key_pixscale"
    junk = fscan (  nsheaders.key_dispaxis, l_key_dispaxis)
    if ("" == l_key_dispaxis) badhdr = badhdr + " key_dispaxis"
    junk = fscan (  nsheaders.sci_ext, l_sci_ext)
    if ("" == l_sci_ext) badhdr = badhdr + " sci_ext"
    junk = fscan (  nsheaders.var_ext, l_var_ext)
    if ("" == l_var_ext) badhdr = badhdr + " var_ext"
    junk = fscan (  nsheaders.dq_ext, l_dq_ext)
    if ("" == l_dq_ext) badhdr = badhdr + " dq_ext"
    junk = fscan (  nsheaders.key_arrayid, l_key_arrayid)
    if ("" == l_key_arrayid) badhdr = badhdr + " key_arrayid"
    junk = fscan (  nsheaders.key_prism, l_key_prism)
    if ("" == l_key_prism) badhdr = badhdr + " key_prism"
    junk = fscan (  nsheaders.key_decker, l_key_decker)
    if ("" == l_key_decker) badhdr = badhdr + " key_decker"
    junk = fscan (  nsheaders.key_fpmask, l_key_fpmask)
    if ("" == l_key_fpmask) badhdr = badhdr + " key_fpmask"
    junk = fscan (  nsheaders.key_exptime, l_key_exptime)
    if ("" == l_key_exptime) badhdr = badhdr + " key_exptime"
    junk = fscan (  nsheaders.key_ron, l_key_ron)
    if ("" == l_key_ron) badhdr = badhdr + " key_ron"
    junk = fscan (  nsheaders.key_gain, l_key_gain)
    if ("" == l_key_gain) badhdr = badhdr + " key_gain"
    junk = fscan (  nsheaders.key_instrument, l_key_instrument)
    if ("" == l_key_instrument) badhdr = badhdr + " key_instrument"

    if (l_logfile == "") {
        junk = fscan (gnirs.logfile, l_logfile)
        if (l_logfile == "") {
            l_logfile = "gnirs.log"
            printlog ("WARNING - NSSDIST: Both nssdist.logfile and \
                gnirs.logfile are", l_logfile, verbose+) 
            printlog ("                   undefined. Using " \ 
                // l_logfile, l_logfile, verbose+) 
        }
    }

    # Start logging
    printlog ("-----------------------------------------------------------\
        -------------------", l_logfile, l_verbose) 
    date | scan (l_struct) 
    printlog ("NSSDIST -- " // l_struct, l_logfile, l_verbose) 
    printlog ("", l_logfile, l_verbose) 

    if (l_database == "") {
        junk = fscan (gnirs.database, l_database)
        if (l_database == "") {
            l_database = "database"
            printlog ("WARNING - NSSDIST: Both nssdist.database and \
                gnirs.database are", l_logfile, verbose+)
            printlog ("                   undefined. Using " // l_database, \
                l_logfile, verbose+)
        }
    }

    # Remove any trailing slashes from the database name
    if (strlen (l_database) > 1 && \
        substr (l_database, strlen (l_database), strlen (l_database)) \
        == "/") {
        l_database = substr (l_database, 1, strlen (l_database) - 1)
    }

    if ("" != badhdr) {
        printlog ("ERROR - NSSDIST: Parameter(s) missing from \
            nsheaders: " // badhdr, l_logfile, verbose+) 
        goto clean
    }

    if (l_fl_inter)
        s_inter = "YES"
    else
        s_inter = "NO"

    if (l_fl_inter && l_fl_dbwrite)
        s_dbwrite = "yes"
    else if (no ==l_fl_inter && l_fl_dbwrite)
        s_dbwrite = "YES"
    else
        s_dbwrite = "NO"

    # Expand and verify input
    if (intdbg) print ("expansion and verification of input")
    gemextn (l_inimages, check="mef,exists", process="none", index="", \
        extname="", extversion="", ikparams="", omit="extension, kernel", \
        replace="", outfile=tmpinlist, logfile="", glogpars="",
        verbose=l_verbose)
    if (0 != gemextn.fail_count || 0 == gemextn.count) {
        printlog ("ERROR - NSSDIST: Problems with input.", l_logfile, \
            verbose+)
        goto clean
    }
    nfiles = gemextn.count
    fl_newout = nfiles > 1

    if (intdbg) print ("nfiles: " // nfiles)

    # Decide on output file
    head (tmpinlist, nlines=1) | scan (l_temp)
    if (fl_newout) {
        l_outspectra = l_temp // l_outsuffix
        printlog ("NSSDIST: Using output file " // l_outspectra, l_logfile, \
            l_verbose) 
        if (imaccess (l_outspectra) ) {
            printlog ("WARNING - NSSDIST: Output file " // l_outspectra \
                // " already exists.", l_logfile, verbose+)
            status = 2 # magic number to avoid deleting file
            goto clean
        }
    } else {
        l_outspectra = l_temp
    }

    # Check PHUs
    gemextn ("@" // tmpinlist, check="mef,exists", process="append", \
        index="0", extname="", extversion="", ikparams="", \
        omit="extension,kernel,", replace="", outfile=tmpphu, 
        logfile="", glogpars="", verbose=l_verbose)
    if (nfiles != gemextn.count || 0 != gemextn.fail_count) {
        printlog ("ERROR - NSSDIST: Problems with PHUs.", l_logfile, verbose+)
        goto clean
    }

    nbad = 0
    scanfile = tmpphu
    while (fscan (scanfile, img) != EOF) {

        if (intdbg) print ("file: " // img)

        keypar (img, "NSSDIST", >& "dev$null")
        if (keypar.found) {
            printlog ("WARNING - NSSDIST: Image " // img // " has been \
                run", l_logfile, verbose+) 
            printlog ("                   through NSSDIST before.", \
                l_logfile, verbose+) 
            if (l_fl_dbwrite && ! fl_newout) {
                printlog ("WARNING - NSSDIST: Existing database entries \
                    will be overwritten.", l_logfile, verbose+) 
            }
        }

        keyfound=""
        hselect(img, "*PREPAR*", yes) | scan(keyfound)
        if (keyfound == "") {
            printlog ("ERROR - NSSDIST: Image " // img \
                // " not PREPAREd.", l_logfile, verbose+)
            nbad += 1
        }

        keypar (img, l_key_instrument, silent+)
        if (intdbg) print (keypar.value)
        if (keypar.value == "F2" || keypar.value == "Flam")
            f2 = yes
    }

    if (nbad > 0) {
        printlog ("ERROR - NSSDIST: " // nbad // " image(s) \
            have not been run through *PREPARE", l_logfile, verbose+) 
        goto clean
    }

    # Check that all images have the same number of SCI extensions.
    # Compare against the first input image
    if (intdbg) print ("checking version numbers")
    head (tmpinlist, nlines=1) | scan (img)
    gemextn (img, check="exists", process="expand", index="", \
        extname=l_sci_ext, extversion="1-", ikparams="", omit="", \
        replace="", outfile="dev$null", logfile="", glogpars="",
        verbose-)
    if (0 != gemextn.fail_count || 0 == gemextn.count) {
        printlog ("ERROR - NSSDIST: Bad science data in " // img // ".",
            l_logfile, verbose+)
        goto clean
    }

    # Assume they're 1..n - this is checked by the extver range below
    nver = int (gemextn.count)

    scanfile = tmpinlist
    while (fscan (scanfile, img) != EOF) {
        gemextn (img, check="exists", process="expand", index="", \
            extname=l_sci_ext, extversion="1-" // nver, \
            ikparams="", omit="", replace="", outfile="dev$null", \
            logfile="", glogpars="", verbose-)
        if (0 != gemextn.fail_count || nver != gemextn.count) {
            printlog ("ERROR - NSSDIST: Bad or missing science data \
                in " // img // ".", l_logfile, verbose+)
            goto clean
        }
    }

    if (fl_newout) {
        printlog ("NSSDIST: Combining input files:", l_logfile, l_verbose) 
        if (l_verbose) type (tmpinlist) 
        type (tmpinlist, >> l_logfile) 
    } else {
        head (tmpinlist, nlines=1) | scan (img)
        printlog ("NSSDIST: Using single input file " // img, l_logfile, \
            l_verbose) 
    }

    # TODO - I haven't cleaned the code related to combination
    # apart from reformatting and some very basic fixes.

    # Find offsets for coordlist, combine input using gemcombine
    if (fl_newout) {
        if (intdbg) print ("combine")

        # Get offsets into headers
        head (tmpinlist, nlines=1) | scan (img)
        phu = img // "[0]"
        keypar (phu, l_key_dispaxis) 
        if (keypar.found)
            l_dispaxis = int (keypar.value)
        if (intdbg) print ("dispaxis: " // l_dispaxis)

        gemoffsetlist ("@"//tmpinlist, l_temp, distance=99999., age=INDEF, \
            targetlist="dev$null", offsetlist="dev$null", fl_nearer+, \
            direction=3-l_dispaxis, fl_younger+, fl_noref-, \
            wcs_source="phu", key_xoff=l_key_xoff, key_yoff=l_key_yoff, \
            key_date="", key_time="", logfile=l_logfile, verbose=l_debug, \
            force=l_force)
            # don't hide possible error message otherwise we'll fail silently
            # on old iraf versions. >& "dev$null") 

        if (gemoffsetlist.status != 0) {
            printlog ("ERROR - NSSDIST: Problem determining offsets \
                using GEMOFFSETLIST.", l_logfile, verbose+) 
            goto clean
        }

        # Check that coordlist not specified
        if (l_coordlist != "") {
            printlog ("ERROR - NSSDIST: coordlist given for multiple \
                input files.", l_logfile, verbose+)
            printlog ("                 Multiple input files imply a \
                single emission source in each.", l_logfile, verbose+)
            goto clean
        }

        # Create coordlist
        l_coordlist = tmpcoordlist
        scanfile = tmpinlist
        while (fscan (scanfile, img) != EOF) {

            phu = img // "[0]"

            keypar (phu, l_key_dispaxis) 
            if (keypar.found)
                l_dispaxis = int (keypar.value) 

            if (l_dispaxis == 1)
                offkey = l_key_yoff
            else
                offkey = l_key_xoff

            keypar (phu, offkey) 
            if (keypar.found) {
                yoff = real (keypar.value) 
                if (intdbg) print ("yoff: " // yoff)
                keypar (phu, l_key_pixscale, silent+) 
                if (keypar.found) {
                    yoff = l_firstcoord - (yoff/real (keypar.value) ) 
                    print (yoff, >> l_coordlist)
                    printlog ("NSSDIST: Constructing coordlist: " \
                        // img // ", " // yoff // " (" // keypar.value \
                        // ")", l_logfile, l_verbose)
                } else {
                    yoff = l_firstcoord - (yoff/l_pixscale) 
                    print (yoff, >> l_coordlist) 
                    printlog ("NSSDIST: Constructing coordlist: " \
                        // img // ", " // yoff, l_logfile, l_verbose)
                }
            } else {
                printlog ("ERROR - NSSDIST: Can't determine the offset \
                    for file " // img, l_logfile, verbose+) 
                goto clean
            }
        }

        # Combine input without rejection
        gemcombine ("@" // tmpinlist, l_outspectra, title = "", \
            combine = "average", reject = "none", offsets = "none", \
            masktype = "goodvalue", maskvalue=0., \
            scale = "none", zero = "none", weight = "none", \
            statsec = "[*,*]", expname = l_key_exptime, \
            lthreshold = INDEF, hthreshold = INDEF, nlow = 1, nhigh = 1, \
            nkeep = 1, mclip+, lsigma = 3, hsigma = 3, \
            key_ron = l_key_ron, key_gain = l_key_gain, ron = 0, \
            gain = 1, snoise = 0, sigscale = 0.1, pclip = -0.5, \
            grow = 0, bpmfile = "", nrejfile = "", sci_ext = l_sci_ext, \
            var_ext = l_var_ext, dq_ext = l_dq_ext, fl_vardq-, \
            logfile = l_logfile, verbose-) 
            # too complex a task to hide error messages - if things crash we
            # get no feedback. logfile = l_logfile, verbose-, >& "dev$null") 
        if (gemcombine.status != 0) {
            printlog ("ERROR - NSSDIST: Problem combining input \
                images using GEMCOMBINE.", l_logfile, verbose+) 
            goto clean
        }

    } else {

        if (l_coordlist == "" || !access (l_coordlist)) {
            l_coordlist = tmpcoordlist
            if (!f2) {
                # The slit edge coordinate information in the header is used
                # for FLAMINGOS-2 MOS data
                printlog ("WARNING - NSSDIST: Only 1 input image without \
                    coordlist.", l_logfile, verbose+) 
                printlog ("                   Will user image centre as \
                    single coordinate.", l_logfile, verbose+) 
                usecentre = yes
            }
        }
    }

    # All spectra to be processed are now combined in l_outspectra
    imgout = l_outspectra
    phu = l_outspectra // "[0]"

    if (intdbg) print ("axis logic")
    keypar (phu, l_key_dispaxis) 
    if (keypar.found)
        l_dispaxis = int (keypar.value) 
    if (intdbg) print ("dispaxis: " // l_dispaxis) 

    # Set dispaxis to the other axis, set section to match that
    l_spataxis = 1
    longslit.dispaxis = 1
    if (l_dispaxis == 1) {
        longslit.dispaxis = 2
        l_spataxis = 2
    }
    if ("default" == l_section) {
        if (1 == l_dispaxis)
            l_section = "middle column"
        else
            l_section = "middle line"
    }

    # Read header info for aperture info
    prism = ""
    if (l_key_prism != unused) {
        keypar (phu, l_key_prism, silent+)
        if (keypar.found) {
            prism = keypar.value
        } else {
            printlog ("WARNING - NSSDIST: No " // l_key_prism // " in " \
                // phu // ".", l_logfile, verbose+)
        }
    }
    if (intdbg) print ("prism: " // prism)

    decker = ""
    if (l_key_decker != unused) {
        keypar (phu, l_key_decker, silent+)
        if (keypar.found) {
            decker = keypar.value
        } else {
            printlog ("WARNING - NSSDIST: No " // l_key_decker // " in " \
                // phu // ".", l_logfile, verbose+)
        }
    }
    if (intdbg) print ("decker: " // decker)

    fpmask = ""
    if (l_key_fpmask != unused) {
        keypar (phu, l_key_fpmask, silent+)
        if (keypar.found) {
            fpmask = keypar.value
        } else {
            printlog ("WARNING - NSSDIST: No " // l_key_fpmask // " in " \
                // phu // ".", l_logfile, verbose+)
        }
    }
    if (intdbg) print ("fpmask: " // fpmask)
    
    arrayid = ""
    if (l_key_arrayid != unused) {
        keypar (phu, l_key_arrayid, silent+)
        if (keypar.found) {
            arrayid = keypar.value
        } else {
            printlog ("WARNING - NSSDIST: No " // l_key_arrayid // " in " \
                // phu // ".", l_logfile, verbose+)
        }
    }
    if (intdbg) print ("arrayid: " // arrayid)

    # Do the work
    for (version = 1; version <= nver; version = version+1) {

        # Pixel coords
        l_crpix = 1
        l_crval = 1
        l_cdelt = 1

        sec = "[" // l_sci_ext // "," // version // "]"

        # Generate default coordlist, if required
        if (usecentre) {
            coord = INDEF
            if (intdbg)
            hselect (img // sec, "i_naxis" // (3-l_dispaxis), yes)
            hselect (img // sec, "i_naxis" // (3-l_dispaxis), yes) \
                | scan (coord)
            if (isindef (coord)) {
                printlog ("ERROR - NSSDIST: cannot generate default \
                    coordlist (no image dimension)!", l_logfile, verbose+)
                goto clean
            }
            delete (l_coordlist, verify-, >& "dev$null")
            coord = (coord + 1) / 2.0
            printf ("NSSDIST: Using central coordinate %4.1f\n", coord) \
                | scan (l_struct)
            printlog (l_struct, l_logfile, l_verbose)
            print (coord, > l_coordlist)
        }

        # Use the slit edge coordinate information in the header for
        # FLAMINGOS-2 MOS data
        if (f2) {
            delete (l_coordlist, verify-, >& "dev$null")
            keypar (img // sec, "LSLTEDGE", silent+)
            if (keypar.found) {
                print (keypar.value, >> l_coordlist)
            } else {
                printlog ("ERROR - NSSDIST: No slit edge coordinate \
                    information in input.", l_logfile, verbose+)
                goto clean
            }
            keypar (img // sec, "RSLTEDGE", silent+)
            if (keypar.found) {
                print (keypar.value, >> l_coordlist)
            } else {
                printlog ("ERROR - NSSDIST: No slit edge coordinate \
                    information in input.", l_logfile, verbose+)
                goto clean
            }
        }

        # Check for enough entries in coordlist
        n_lines = 0
        count (l_coordlist) | scan (n_lines) 
        if (n_lines < 1) {
            printlog ("ERROR - NSSDIST: at least 1 position is needed", 
                l_logfile, verbose+) 
            printlog ("                 in coordlist.", l_logfile, verbose+) 
            goto clean
        }

        # Update coordlist to include shift for this particular version
        aprow = mktemp ("aprow")//".fits"
        search = "prism == '" // prism // "' && " \
            // "slit == '" // fpmask // "' && " \
            // "decker == '" // decker // "' && " \
            // "arrayid == '" // arrayid // "' && " \
            // "order == " // version
        if (intdbg) print (search)
        tselect (l_aptable, aprow, search)

        if (intdbg) tprint (aprow)
        tinfo (table = aprow, >>& "dev$null")
        usetable = (1 == tinfo.nrows)
        if (intdbg) print (tinfo.nrows // ", " // usetable)

        if (no == usetable) {
            if (1 < nver) {
                printlog ("WARNING - NSSDIST: No aperture parameters \
                    in " // l_aptable, l_logfile, verbose+)
                printlog ("                   for " // search, \
                    l_logfile, verbose+)
                printlog ("                   so no shift applied for \
                    each order.", l_logfile, verbose+)
            }
            apshift = 0.0
            secn = l_section
            subsec = "[*,*]"
        } else {
            tprint (aprow, prparam-, prdata+, \
                showrow-, showhdr-, showunits-, \
                col="apshift,apline,sdsample", rows=1, \
                pwidth=160) | scan (apshift,secn,subsec)
            if (1 == l_dispaxis) {
                secn = "column " // secn
                subsec = "[" // subsec // ",*]"
            } else {
                secn = "line " // secn
                subsec = "[*," // subsec // "]"
            }
            printlog ("NSSDIST: Using the following values from " \
                // l_aptable, l_logfile, l_verbose)
            printlog ("         where " // search // ":", \
                l_logfile, l_verbose)
            printlog ("         aperture at (column apline): " \
                // secn, l_logfile, l_verbose)
            printlog ("         shift (column apshift):      " \
                // apshift, l_logfile, l_verbose)
            # not using subsec at the moment
        }
        delete (aprow, verify-, >>& "dev$null")

        tmpshiftcoord = mktemp ("tmpshiftcoord")
        if (n_lines > 0) {
            scanfile = l_coordlist
            while (fscan (scanfile, shift) != EOF) {
                print ((shift + apshift), >> tmpshiftcoord)
            }
        } else {
            print (1, >> tmpshiftcoord)
        }

        # Write offsets to logfile
        printlog (" ", l_logfile, l_verbose) 
        printlog ("NSSDIST: Coordinate list:", l_logfile, l_verbose)
        type (tmpshiftcoord, >> l_logfile)
        if (l_verbose) type (tmpshiftcoord) 

        # Delete database entry if fl_dbwrite=yes
        if (l_fl_dbwrite) {
            delete (l_database // "/id" // imgout, verify-, >& "dev$null") 
        }

        # This seems to be necessary to lose the coordinates from the
        # original frame (originally only if !l_fl_inter)

        # autoidentify think that it is a rotated spectrum - but it isn't
        # so need to clean them. Clean the rest as well. 
        gemhedit (imgout // sec, "LTV1", "", "", delete+)
        gemhedit (imgout // sec, "LTV2", "", "", delete+)
        gemhedit (imgout // sec, "LTM1_1", "", "", delete+) 
        gemhedit (imgout // sec, "LTM2_2", "", "", delete+) 
        gemhedit (imgout // sec, "WAT0_001", "", "", delete+) 
        gemhedit (imgout // sec, "WCSDIM", "", "", delete+) 

        printlog ("NSSDIST: Running IDENTIFY", l_logfile, l_verbose) 

        if (l_fl_inter) {

            specred.logfile = l_logfile
            print ("") 
            print ("NSSDIST: Mark features to use") 
            print ("   Accept coordinate or assign new coordinate to \
                each feature") 
            print ("   m - mark feature ") 
            print ("   ? - help, q - exit ")

            if (intdbg) print ("section: " // secn)

            identify (imgout // sec, l_crval, l_cdelt, \
                section = secn, database = l_database, \
                coordlist = tmpshiftcoord, units = "", nsum = l_nsum, \
                match = l_match, maxfeatures = n_lines, zwidth = 100, \
                ftype = l_ftype, fwidth = l_fwidth, cradius = l_cradius, \
                threshold = l_threshold, minsep = l_minsep, \
                function = l_function, order = l_order, \
                sample = l_sample, niterate = l_niterate, \
                low_reject = l_low_reject, high_reject = l_high_reject, \
                grow = l_grow, autowrite = l_fl_dbwrite, \
                graphics = "stdgraph", cursor = "", \
                aidpars.debug = l_aiddebug, aidpars.crpix = l_crpix)
            specred.logfile = ""

        } else {

            # All new code - identifies the n_lines strongest lines
            # Old code didn't work with small numbers of traces
            specred.logfile = l_logfile
            tmpfakedin = mktemp ("tmpfakedin")
            print ("e\nq\n", > tmpfakedin)
            identify (imgout // sec, l_crval, l_cdelt, \
                section = secn, database = l_database, \
                coordlist = tmpshiftcoord, units = "", nsum = l_nsum, \
                match = l_match, maxfeatures = n_lines, zwidth = 100, \
                ftype = l_ftype, fwidth = l_fwidth, cradius = l_cradius, \
                threshold = l_threshold, minsep = l_minsep, \
                function = l_function, order = l_order, \
                sample = l_sample, niterate = l_niterate, \
                low_reject = l_low_reject, high_reject = l_high_reject, \
                grow = l_grow, autowrite = l_fl_dbwrite, \
                graphics = "stdgraph", cursor = tmpfakedin, \
                aidpars.debug = l_aiddebug, aidpars.crpix = l_crpix, \
                >G "dev$null")
            delete (tmpfakedin, verify-, >& "dev$null")
            specred.logfile = ""
        }

        if (l_fl_dbwrite) {

            # Change name to reflect individual extensions
            destname = l_database // "/id" // imgout // "_" // l_sci_ext \
                // "_" // version // "_"
            delete (destname, verify-, >& "dev$null")

            tmpname = l_database // "/id" // imgout

            if (no ==access (tmpname)) {

                # Currently only a warning, as sometimes impossible to 
                # trace x-dispersed orders
                printlog ("WARNING - NSSDIST: No s-distortion solution \
                    found", l_logfile, verbose+) 

            } else {

                # Check for INDEF user coords
                match ("INDEF", tmpname, meta-, stop-, print-) \
                    | count | scan (nindef)
                match ("begin", tmpname, meta-, stop-, print-) \
                    | count | scan (nbegin)
                if (nindef >= nbegin) {
                    if (l_fl_inter) {
                        printlog ("ERROR - NSSDIST: Undefined user \
                            coordinates in database (" // tmpname \
                            // "),", l_logfile, verbose+)
                        printlog ("                 you need to enter a \
                            pixel value when marking the line.", \
                            l_logfile, verbose+)
                        goto clean
                    } else {
                        printlog ("WARNING - NSSDIST: Undefined user \
                            coordinates in database (" // tmpname \
                            // "),", l_logfile, verbose+)
                        printlog ("                   values in coordlist \
                            or aptable or firstcoord may be bad.", \
                            l_logfile, verbose+)
                        finalwarn = yes
                    }
                }

                printlog ("NSSDIST: Running REIDENTIFY", l_logfile, l_verbose)

                reidentify (imgout // sec, imgout // sec, \
                    coordlist = tmpshiftcoord, interactive = l_fl_inter, \
                    section = secn, newaps = yes, refit = l_refit, \
                    trace = l_trace, step = l_step, nsum = l_nsum, \
                    shift = "0.", search = 0., nlost = l_nlost, \
                    cradius = l_cradius, threshold = l_threshold, \
                    addfeatures = no, match = 6, maxfeatures = 250, \
                    minsep = l_minsep, override = l_fl_dbwrite, \
                    database = l_database, verbose = l_verbose, \
                    logfile = l_logfile, plotfile = "", \
                    graphics = "stdgraph", cursor = "") 

                # Rename database entry
                rename (tmpname, destname, field="all")
            }
        }
        delete (tmpshiftcoord, verify-, >>& "dev$null")
    } # end for-loop

    #----------------------------------------------------------------------
    # Update header
    gemdate () 

    if (fl_newout) {
        ii = 1
        scanfile = tmpinlist
        while (fscan (scanfile, img) != EOF) {
            printf ("%-8s= \'%-18s\' / %-s\n", "NSSDIS0" // str (ii), img,
                "NSSDIST input image", >> tmphead)
            ii += 1
        }
        # put all the new stuff in the header
        mkheader (imgout // "[0]", tmphead, append+, verbose-) 

    }
    gemhedit (imgout // "[0]", "NSSDIST", gemdate.outdate, 
        "UT Time stamp for NSSDIST", delete-) 
    gemhedit (imgout // "[0]", "GEM-TLM", gemdate.outdate, 
        "UT Last modification with GEMINI", delete-) 

    # Remove REFSPEC1 since this otherwise refers to pinholes not wavelength
    gemhedit (imgout // sec, "REFSPEC1", "", "", delete+) 

    # Completed successfully
    status = 0

    #----------------------------------------------------------------------
    # Clean up
clean:
    scanfile = ""
    delete (tmpinlist, verify-, >>& "dev$null") 
    delete (tmphead // "," // tmpcoordlist // "," // tmpphu, \
        verify-, >& "dev$null") 
    if (status != 0 && status != 2 && fl_newout)
    imdelete (l_outspectra // ".fits", verify-, >& "dev$null") 

    printlog ("", l_logfile, l_verbose) 
    if (status == 0) {
        if (finalwarn) {
            printlog ("WARNING - NSSDIST: Check messages above about \
                undefined values.", l_logfile, verbose+)
        }
        printlog ("NSSDIST Exit status: good.", l_logfile, l_verbose) 
    } else {
        printlog ("NSSDIST Exit status: failed.", l_logfile, l_verbose) 
    }
    printlog ("------------------------------------------------------------\
        ------------------", l_logfile, l_verbose) 

end


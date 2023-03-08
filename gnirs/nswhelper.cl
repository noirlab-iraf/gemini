# Copyright(c) 2001-2015 Association of Universities for Research in Astronomy, Inc.

procedure nswhelper (lampspectra) 

# Wavelength calibration of GNIRS/NIRI arclamp spectra
#
# Version: Sept 20, 2002 JJ v1.4 release
#          Aug 19, 2003  KL IRAF2.12 - new parameter, addonly, in hedit
#          Oct 29, 2003  Kl moved from niri to gnirs package
#
# Planned: automatic linelist selection based on header keyword;
#          take out nsappwave call if added to nscut for all spectra;

char    lampspectra {prompt = "Input GNIRS/NIRI lamp spectra"}              # OLDP-1-input
char    outspectra  {"", prompt = "Output spectra"}                         # OLDP-1-output-single
char    outprefix   {"w", prompt = "Prefix for output spectra"}             # OLDP-4
char    reference   {"", prompt = "Reference wavelength image"}             # OLDP-3
real    crval       {INDEF, prompt = "Approximate wavelength at coordinate reference pixel"} # OLDP-3
real    cdelt       {INDEF, prompt = "Approximate dispersion"}              # OLDP-3
real    crpix       {INDEF, prompt = "Coordinate reference pixel"}          # OLDP-3
int     dispaxis    {1, min = 1, max = 2, prompt = "Dispersion axis if not defined in the header"} # OLDP-3

char    database    {"", prompt = "Directory for files containing feature data"} # OLDP-3
char    coordlist   {"gnirs$data/argon.dat", prompt = "User coordinate list, line list"} # OLDP-3
bool    fl_inter    {no, prompt = "Examine identifications interactively"}  # OLDP-4

char    section     {"default", prompt = "Image section for running identify"}  # OLDP-2
int     nsum        {10, min = 1, prompt = "Number of lines or columns to sum"} # OLDP-2
char    ftype       {"emission", min = "emission|absorption", prompt = "Feature type"} # OLDP-2
real    fwidth      {4., min = 2, prompt = "Feature width in pixels"}       # OLDP-2
real    cradius     {5., min = 2, prompt = "Centering radius in pixels"}    # OLDP-2
real    threshold   {0., prompt = "Feature threshold for centering"}        # OLDP-2
real    minsep      {2., prompt = "Minimum pixel separation for features"}  # OLDP-2
real    match       {-6., prompt = "Coordinate list matching limit, <0 pixels, >0 user"} # OLDP-2
char    function    {"chebyshev", min = "legendre|chebyshev|spline1|spline3", prompt = "Coordinate fitting function"} # OLDP-2
int     order       {4, min = 1, prompt = "Order of coordinate fitting function"} # OLDP-2
char    sample      {"*", prompt = "Coordinate sample regions"}             # OLDP-2
int     niterate    {10, min = 0, prompt = "Rejection iterations"}          # OLDP-2
real    low_reject  {3., min = 0, prompt = "Lower rejection sigma"}         # OLDP-2
real    high_reject {3., min = 0, prompt = "Upper rejection sigma"}         # OLDP-2
real    grow        {0., min = 0, prompt = "Rejection growing radius"}      # OLDP-2
bool    refit       {yes, prompt = "Refit coordinate function when running reidentify"} # OLDP-2
int     step        {10, prompt = "Steps in lines or columns for reidentification"} # OLDP-2
bool    trace       {no, prompt = "Use fit from previous step rather than central aperture"} # OLDP-2
int     nlost       {3, min = 0, prompt = "Maximum number of lost features"}    # OLDP-2
char    aiddebug    {"", prompt = "Debug parameter for aidpars"}            # OLDP-3
real    fmatch      {0.2, prompt = "Acceptable fraction of unmatched reference lines?"} # OLDP-3
int     nfound      {6, prompt = "Miniumum number of identified lines required?"}
real    sigma       {0.05, prompt = "Sigma of line centering (pixels)"}
real    rms         {0.1, prompt = "RMS goal (fwidths)"}
char    nsappwavedb {"gnirs$data/nsappwave.fits", prompt = "Database for nsappwave"} # OLDP-2
bool    fl_overwrite    {yes, prompt = "Overwrite existing database entries"}   # OLDP-2
char    logfile     {"", prompt = "Logfile"}                                # OLDP-1
bool    verbose     {yes, prompt = "Verbose"}                               # OLDP-2
bool    debug       {no, prompt = "Very verbose"}                           # OLDP-3
int     status      {0, prompt = "Exit status (0=Good, 1=Fail, 2=Some bad cal)"} # OLDP-4
struct  *scanfile1  {"", prompt = "Internal use"}                           # OLDP-4
struct  *scanfile2  {"", prompt = "Internal use"}                           # OLDP-4

begin

    char    l_inimages = ""
    char    l_outspectra = ""
    char    l_outprefix = ""
    char    l_reference = ""
    real    l_crval
    real    l_cdelt
    real    l_crpix
    int     l_dispaxis
    char    l_database = ""
    char    l_coordlist = ""
    bool    l_fl_inter
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
    real    l_fmatch
    int     l_nfound
    real    l_sigma
    real    l_rms
    char    l_nsappwavedb = ""
    bool    l_fl_overwrite
    char    l_logfile = ""
    bool    l_verbose
    bool    l_debug

    char    l_sci_ext = ""
    char    l_key_dispaxis = ""
    char    l_key_order = ""
    char    l_key_cradius = ""

    char    tmpinlist, tmpoutlist, tmpin, s_inter, tmpphu, tmpxd, tmpfile1
    char    l_temp, img, phu, tmpref, keyfound, instrument, imgin, imgout, sec
    char    unixsec, badhdr, name, secn, ref, fullref, pixel_axis
    int     nbad, nfiles, count, junk, nver, version, axis, nfail
    int     worder, prevworder, arc_lines_centre
    bool    fl_runbefore, intdbg, fl_reference, do_autoid, auto_find_section
    real    radius
    struct  l_date

    # Define temporary files

    tmpin = mktemp ("tmpin") 
    tmpphu = mktemp ("tmpphu") 
    tmpinlist = mktemp ("tmpinl") 
    tmpoutlist = mktemp ("tmpout") 
    tmpfile1 = mktemp ("tmpfile1") 
    tmpxd = mktemp ("tmpxd")
    tmpref = ""

    # Initialize parameters
    nfail = 0
    status = 1
    intdbg = no
    auto_find_section = no

    # Cache parameter files
    cache ("niri", "nsappwave", "nswavelength", "gemextn", "gemdate") 

    # Set the local variables

    junk = fscan (  lampspectra, l_inimages)
    junk = fscan (  outspectra, l_outspectra)
    junk = fscan (  outprefix, l_outprefix)
    junk = fscan (  reference, l_reference)
    l_crval =       crval
    l_cdelt =       cdelt
    l_crpix =       crpix
    l_dispaxis =    dispaxis
    junk = fscan (  database, l_database)
    junk = fscan (  coordlist, l_coordlist)
    l_fl_inter =    fl_inter
    junk = fscan (  section, l_section)
    l_nsum =        nsum
    junk = fscan (  ftype, l_ftype)
    l_fwidth =      fwidth
    l_cradius =     cradius
    l_threshold =   threshold
    l_minsep =      minsep
    l_match =       match
    junk = fscan (  function, l_function)
    l_order =       order
    junk = fscan (  sample, l_sample)
    l_niterate =    niterate
    l_low_reject =  low_reject
    l_high_reject = high_reject
    l_grow =        grow
    l_refit =       refit
    l_step =        step
    l_trace =       trace
    l_nlost =       nlost
    junk = fscan (  aiddebug, l_aiddebug)
    l_fmatch =      fmatch
    l_nfound =      nfound
    l_sigma =       sigma
    l_rms =         rms
    junk = fscan (  nsappwavedb, l_nsappwavedb)
    l_fl_overwrite =    fl_overwrite
    junk = fscan (  logfile, l_logfile)
    l_verbose =     verbose
    l_debug =       debug

    badhdr = ""
    junk = fscan (nsheaders.sci_ext, l_sci_ext)
    if ("" == l_sci_ext) badhdr = badhdr + " sci_ext"
    junk = fscan (nsheaders.key_dispaxis, l_key_dispaxis)
    if ("" == l_key_dispaxis) badhdr = badhdr + " key_dispaxis"
    junk = fscan (nsheaders.key_order, l_key_order)
    if ("" == l_key_order) badhdr = badhdr + " key_order"
    junk = fscan (nsheaders.key_cradius, l_key_cradius)
    if ("" == l_key_cradius) badhdr = badhdr + " key_cradius"

    if (no == l_verbose) l_aiddebug = ""
    if (intdbg) print ("verbose/debug: " // l_verbose // " " // l_aiddebug)

    if (l_logfile == "") {
        junk = fscan (gnirs.logfile, l_logfile)
        if (l_logfile == "") {
            l_logfile = "gnirs.log"
            printlog ("WARNING - NSWHELPER: Both nswavelength.logfile \
                and gnirs.logfile", l_logfile, verbose+) 
            printlog ("                         are empty. Using "\
                // l_logfile, l_logfile, verbose+) 
        }
    }

    s_inter = str (l_fl_inter) 
    if (no ==l_fl_inter) s_inter = "NO"

    # Write to the logfile

    printlog ("--------------------------------------------------------\
        -----------------------", l_logfile, l_verbose) 
    date | scan (l_date) 
    printlog ("NSWHELPER -- " // l_date, l_logfile, l_verbose) 
    printlog ("", l_logfile, l_verbose) 

    if (l_database == "") {
        junk = fscan (gnirs.database, l_database)
        if (l_database == "") {
            l_database = "database"
            printlog ("WARNING - NSWHELPER: Both nswhelper.database \
                and gnirs.database are", l_logfile, verbose+)
            printlog ("                      undefined. Using " \
                // l_database, l_logfile, verbose+)
        }
    }

    # Remove any trailing slashes from the database name
    if (strlen (l_database) > 1 && \
        substr (l_database, strlen (l_database), strlen (l_database)) \
        == "/") {
        l_database = substr (l_database, 1, strlen (l_database) - 1)
    }

    # Check access to coordinate file

    if (l_coordlist == "") {
        printlog ("ERROR - NSWHELPER: Coordinate list not defined", \
            l_logfile, verbose+) 
        goto clean
    }
    if (no ==access (l_coordlist) ) {
        printlog ("ERROR - NSWHELPER: Coordinate list " \
            // l_coordlist // " not found", l_logfile, verbose+) 
        goto clean
    }

    if ("" != badhdr) {
        printlog ("ERROR - NSWHELPER: Parameter(s) missing from \
            nsheaders: " // badhdr, l_logfile, verbose+) 
        goto clean
    }

    # Expand and verify input

    if (intdbg) print ("expansion and verification of input")

    gemextn (l_inimages, check="mef,exists", process="none", index="", \
        extname="", extversion="", ikparams="", omit="extension, kernel", \
        replace="", outfile=tmpfile1, logfile="", glogpars="",
        verbose=l_verbose)
    if (0 != gemextn.fail_count || 0 == gemextn.count) {
        printlog ("ERROR - NSWHELPER: Problems with input.", \
            l_logfile, verbose+)
        goto clean
    }
    nfiles = gemextn.count

    # Check PHUs

    gemextn ("@" // tmpfile1, check="mef,exists", process="append", \
        index="0", extname="", extversion="", ikparams="", \
        omit="extension,kernel,", replace="", outfile=tmpphu, 
        logfile="", glogpars="", verbose=l_verbose)
    if (nfiles != gemextn.count || 0 != gemextn.fail_count) {
        printlog ("ERROR - NSWHELPER: Problems with PHUs.", \
            l_logfile, verbose+)
        goto clean
    }

    scanfile1 = tmpphu
    scanfile2 = tmpfile1
    nbad = 0
    nfiles = 0
    while (fscan (scanfile1, phu) != EOF) {

        junk = fscan (scanfile2, img)
        if (intdbg) print ("file: " // phu)
        fl_runbefore = no

        keypar (phu, "NSWAVELE", >& "dev$null")
        if (keypar.found) {
            printlog ("WARNING - NSWHELPER: Image " // img // " has \
                been run", l_logfile, verbose+) 
            printlog ("                   through NSWAVELENGTH before.", \
                l_logfile, verbose+) 
            fl_run_before = yes
        }

        keyfound = ""
        hselect(phu, "*PREPAR*", yes) | scan(keyfound)
        if (keyfound == "") {
            printlog ("WARNING - NSWHELPER: Image " // img // " has \
                not been PREPAREd", l_logfile, verbose+) 
            nbad = nbad + 1
        }

        if (no == fl_runbefore) {
            print (img, >> tmpinlist)
            nfiles = nfiles + 1
        }

       hselect (phu,"INSTRUME", expr=yes) | scan (instrument)

    }

    # Check for empty file list

    if (nfiles == 0) {
        printlog ("ERROR - NSWHELPER: No input images to process.", 
            l_logfile, verbose+) 
        goto clean
    }

    # Exit if problems found with input files

    if (nbad > 0) {
        printlog ("ERROR - NSWHELPER: " // nbad // " bad image(s)", \
            l_logfile, verbose+) 
        goto clean
    }

    printlog ("Using input files:", l_logfile, l_verbose) 
    if (l_verbose) type (tmpinlist) 
    type (tmpinlist, >> l_logfile) 
    delete (tmpfile1, verify-) 

    # Expand and verify output

    if (intdbg) print ("expansion and verification of output")

    nbad = 0
    if (l_outspectra != "") {
        # If outspectra is supplied by the user, it takes precedence over
        # outprefix, so check if outspectra already exists
        gemextn (l_outspectra, check="absent", process="none", index="", \
            extname="", extversion="", ikparams="", omit="extension, kernel", \
            replace="", outfile=tmpoutlist, logfile="", glogpars="", verbose-)
        if (gemextn.fail_count != 0) {
            nbad += 1
        }
    } else {
        # Check if "outprefix + lampspectra" already exists
        gemextn ("%^%" // l_outprefix // "%" // "@" // tmpinlist,
            check="absent", process="none", index="", extname="", \
            extversion="", ikparams="", omit="extension, kernel", replace="", \
            outfile=tmpoutlist, logfile="", glogpars="", verbose-)
        if (gemextn.fail_count != 0) {
            nbad += 1
        }
    }

    if (nbad > 0) {
        printlog ("ERROR - NSWHELPER: Existing or incorrectly formatted \
            output files.", l_logfile, verbose+) 
        goto clean
    }

    printlog ("Using output files:", l_logfile, l_verbose) 
    if (l_verbose) type (tmpoutlist) 
    type (tmpoutlist, >> l_logfile) 

    # Check number of input and output images

    if (nfiles != gemextn.count) {
        printlog ("ERROR - NSWHELPER: Different number of input and \
            output files", l_logfile, verbose+) 
        goto clean
    }

    # Check that all images have the same number of SCI
    # extensions.

    # Compare against the first input image

    if (intdbg) print ("checking version numbers")

    head (tmpinlist, nlines=1) | scan (img)
    gemextn (img, check="exists", process="expand", index="", \
        extname=l_sci_ext, extversion="1-", ikparams="", omit="", \
        replace="", outfile="dev$null", logfile="", glogpars="",
        verbose-)
    if (0 != gemextn.fail_count || 0 == gemextn.count) {
        printlog ("ERROR - NSWHELPER: Bad science data in " // img \
            // ".", l_logfile, verbose+)
        goto clean
    }
    # Assume they're 1..n - this is checked by the extver range below
    nver = int (gemextn.count)

    scanfile1 = tmpinlist
    while (fscan (scanfile1, img) != EOF) {
        gemextn (img, check="exists", process="expand", index="", \
            extname=l_sci_ext, extversion="1-" // nver, \
            ikparams="", omit="", replace="", outfile="dev$null", \
            logfile="", glogpars="", verbose-)
        if (0 != gemextn.fail_count || nver != gemextn.count) {
            printlog ("ERROR - NSWHELPER: Bad or missing science data \
                in " // img // ".", l_logfile, verbose+)
            goto clean
        }
    }

    # Make the file list
    joinlines (tmpinlist, tmpoutlist, output=tmpin, delim=" ", missing=" ", \
        maxchars=161, shortest+, verbose-)

    # Check reference exists, if given

    fl_reference = l_reference != ""
    if (fl_reference) {
        gemextn (l_reference, check="exists", process="none", index="", \
            extname="", extversion="", ikparams="", omit="", replace="", \
            outfile="dev$null", logfile="", glogpars="", verbose-)
        if (0 != gemextn.fail_count || 1 != gemextn.count) {
            printlog ("ERROR - NSWHELPER: Bad or missing reference \
                data in " // l_reference // ".", l_logfile, verbose+)
            goto clean
        }
    }

    # Main loop

    scanfile1 = tmpin
    while (fscan (scanfile1, imgin, imgout) != EOF) {

        keypar (imgin // "[0]", l_key_dispaxis, silent+) 
        if (keypar.found)
            axis = int (keypar.value) 
        else
            axis = l_dispaxis

        if (intdbg) print ("dispaxis: " // axis)
        if (l_section != "") l_section = section # may contain spaces
        if (l_section == "default" || l_section == "") {
            if (axis == 1) {
                secn = "middle line"
            } else if (axis == 2) {
                secn = "middle column"
            }
        } else if (l_section == "auto") {
            # Automatically determine the section for autoidentify and
            # reidentify (this option is normally used for XD data that has
            # already had the s-distortion correction applied)
            secn = l_section
            auto_find_section = yes
        } else {
            secn = l_section
        }
        if (intdbg) print ("section: " // secn)

        # Run nsappwave if necessary

        keypar (imgin // "[0]", "NSAPPWAV", >& "dev$null") 
        if (no == keypar.found) {
            if (intdbg) print ("nsappwave with " // l_nsappwavedb)
            nsappwave (imgin, outspectra = imgout, outprefix = "", \
                nsappwavedb = l_nsappwavedb, crval = l_crval, \
                cdelt = l_cdelt, logfile = l_logfile, \
                verbose = l_verbose, debug = l_debug) 

            if (nsappwave.status > 0) {
                printlog ("ERROR - NSWHELPER: nsappwave returned \
                    with error", l_logfile, verbose+) 
                goto clean
            } 

        } else {
            if (intdbg) print ("already nsappwaved")
            copy (imgin // ".fits", imgout // ".fits", verbose+)
        }

        tmpref = ""
        prevworder = -1

        for (version = 1; version <= nver; version = version + 1) {

            sec = "[" // l_sci_ext // "," // version // "]"
            printlog ("Lamp image: " // imgin // sec, l_logfile, l_verbose)

            # Get the approximate wavelength info from the header
            radius = indef
            hselect (imgout // sec, "CRPIX" // axis // ",CRVAL" // axis \
                // ",CD" // axis // "_" // axis // "," // l_key_cradius, \
                yes) | scan (l_crpix, l_crval, l_cdelt, radius)
            if (no == isindef (l_cradius)) radius = l_cradius
            if (intdbg) {
                hselect (imgout // sec, \
                    "$I,CRPIX,CRVAL,CD"// axis // "_" // axis, yes)
                if (no == isindef (l_crpix))
                    print ("crpix  " // l_crpix)
                if (no == isindef (l_crval))
                    print ("crval  " // l_crval)
                if (no == isindef (l_cdelt))
                    print ("cd     " // l_cdelt)
                if (no == isindef (radius))
                    print ("radius " // radius)
            }

            if (isindef (cradius)) {
                printlog ("ERROR - NSWHELPER: No cradius value " \
                    // imgin // " (extver " // version // ").", \
                    l_logfile, verbose+)
                goto clean
            }

            longslit.dispaxis = axis

            do_autoid = yes
            ref = imgout
            fullref = imgout // sec

            worder = version

            if (instrument == "GNIRS") {
                keypar (imgout // sec, l_key_order, silent+)
                if (keypar.found)
                    worder = int (keypar.value)
            }

            if (intdbg) print ("worder: " // worder // ", " // prevworder)

            if (tmpref != "" && worder != prevworder) {
                if (intdbg) print ("deleting old ifu reference copy")
                imdelete (tmpver, verify-, >& "dev$null")
                tmpref = ""
            }

            if (fl_reference) {

                # Rename file from stored (per extension) format

                if (intdbg) print ("copying in reference")
                name = l_database // "/id" // l_reference // "_" \
                    // l_sci_ext // "_" // version // "_"
                if (access (name)) {
                    copy (name, l_database // "/id" // l_reference)
                    do_autoid = no
                    ref = l_reference
                    fullref = l_reference // sec
                } else {
                    printlog ("WARNING - NSWHELPER: No reference " \
                        // l_reference // " (extver " // version // ").", \
                        l_logfile, verbose+)
                }

            } else if (worder == prevworder) {

                # This is all a bit of a hack, unfortunately.

                # IRAF uses a single database file for all extensions
                # but Gemini doesn't. So we copy database files around
                # as necessary. That's fine, but what happens when we
                # want to use the previous order of the same file as
                # a reference? Answer - we can't, because the database
                # info has been copied into an extension specific
                # record. So we dream up a new image - a copy of the
                # image we're calibrating - and construct a database
                # entry so that we can use that file as reference.
                # since it has a different name we avoid the conflicts
                # and things (should) work...

                if (tmpref == "") {
                    if (intdbg) print ("copying image for ifu-reference")
                    tmpref = mktemp ("ifu-reference-")
                    copy (imgout // ".fits", tmpref // ".fits")
                }

                if (intdbg) print ("copying in previous extn")
                name = l_database // "/id" // imgout // "_" \
                    // l_sci_ext // "_" // (version-1) // "_"
                if (access (name)) {
                    ref = tmpref
                    fullref = tmpref // "[" // l_sci_ext // "," \
                        // (version-1) // "]"
                    if (intdbg)
                    print ("copy " // name // " -> " // ref)
                    sed ("-e", "s/" // imgout // "/" // tmpref // "/g", \
                        name, > l_database // "/id" // ref)
                    do_autoid = no
                } else {
                    printlog ("WARNING - NSWHELPER: previous \
                        extension missing " // imgout \
                        // " (extver " // version // ").", \
                        l_logfile, verbose+)
                }

            }

            if (do_autoid) {

                if (intdbg) print ("autoidentify")

                if (auto_find_section) {
                    # Determine the approximate centre of the arc lines. First,
                    # project the extension along the spectral axis, which sums
                    # the arc lines to create a 1D image
                    improject (imgout // sec, tmpxd, projaxis=axis, \
                        average-, highcut=0., lowcut=0., pixtype="real", \
                        verbose-)

                    # Replace pixels with a value greater than 100 (randomly
                    # selected threshold that should correspond only to the
                    # summed arc lines; the background pixels in XD data at
                    # this stage of the reduction are approximately 0) with the
                    # value of the pixel coordinate in the 1D image, so that
                    # the mean of the resulting image corresponds to the
                    # central pixel coordinate of the arc lines
                    if (axis == 1)
                        pixel_axis = "J"
                    else
                        pixel_axis = "I"

                    imexpr ("(a > 100) ? " // pixel_axis // " : 0", 
                        tmpxd // "[0,overwrite+]", tmpxd, dims="auto", \
                        intype="auto", outtype="auto", refim="auto", \
                        bwidth=0, btype="nearest", bpixval=0., \
                        rangecheck=yes, verbose=yes, exprdb="none")

                    imstat (tmpxd, fields="mean", lower=1, \
                        upper=INDEF, nclip=0, lsigma=3., usigma=3., \
                        binwidth=0.1, format=no, cache=no) | \
                        scan (arc_lines_centre)

                    imdelete (tmpxd, verify-, >& "dev$null")

                    if (axis == 1)
                        secn = "[*," // arc_lines_centre // "]"
                    else
                        secn = "[" // arc_lines_centre // ",*]"

                    if (intdbg) {
                        print ("Automatically determined section: " // secn)
                    }
                }

                autoidentify (imgout // sec, l_crval, l_cdelt, \
                    coordlist = l_coordlist, units = "Angstrom", \
                    interactive = s_inter, aidpars = "", section = secn, \
                    nsum = l_nsum, ftype = l_ftype, fwidth = l_fwidth, \
                    cradius = radius, threshold = l_threshold, \
                    minsep = l_minsep, match = l_match, \
                    function = l_function, order = l_order, \
                    sample = l_sample, niterate = l_niterate, \
                    low_reject = l_low_reject, \
                    high_reject = l_high_reject, \
                    grow = l_grow, dbwrite = "YES", \
                    overwrite = l_fl_overwrite, database = l_database, \
                    verbose = l_verbose, logfile = l_logfile, \
                    plotfile = "", graphics = "stdgraph", cursor = "", \
                    aidpars.fmatch = l_fmatch, \
                    aidpars.nfound = l_nfound, \
                    aidpars.debug = l_aiddebug, aidpars.cddir = "sign", \
                    aidpars.crpix = l_crpix, aidpars.sigma = l_sigma, \
                    aidpars.rms = l_rms)

                if (intdbg) print (aidpars.fmatch)

            }

            if (intdbg) print ("reidentify")

            name = l_database // "/id" // ref
            if (intdbg) print ("checking " // name)

            if (no == access (name)) {

                printlog ("WARNING - NSWHELPER: no calibration for " \
                    // ref // " (extver " // version // ").", \
                    l_logfile, verbose+)
                nfail = nfail + 1

            } else {

                if (intdbg) print ("refit: " // l_refit)

                reidentify (reference=fullref, images = imgout // sec, \
                    coordli = l_coordlist, interactive = s_inter, \
                    section = secn, newaps = yes, refit = l_refit, \
                    trace = l_trace, step = l_step, nsum = l_nsum, \
                    shift = "0.", search = 0., nlost = l_nlost, \
                    cradius = radius, threshold = l_threshold, \
                    addfeatures = no, match = l_match, maxfeatures = 250, \
                    minsep = l_minsep, override = l_fl_overwrite, \
                    database = l_database, verbose = l_verbose, \
                    logfile = l_logfile, plotfile = "", \
                    graphics = "stdgraph", cursor = "") 

                # Rename the idfile so fitcoords can find it
                name = l_database // "/id" // imgout // "_" \
                    // l_sci_ext // "_" // version // "_"

                if (l_fl_overwrite) {
                    if (intdbg) print ("deleting " // name)
                    delete (name, verify-, >& "dev$null")
                }

                if (intdbg) print ("moving")
                rename (l_database // "/id" // imgout, name, field="all")

                if (no == do_autoid) {
                    if (intdbg) print ("tidying away reference")
                    delete (l_database // "/id" // ref, verify-, \
                    >& "dev$null")
                }
            }
            prevworder = worder
        } # end of for-loop

        # update header
        gemdate ()
        gemhedit (imgout // "[0]", "NSWAVELE", gemdate.outdate, \
            "UT Time stamp for NSWAVELENGTH", delete-) 
        gemhedit (imgout // "[0]", "GEM-TLM", gemdate.outdate, \
            "UT Last modification with GEMINI", delete-) 

    }

    # Completed successfully
    if (0 == nfail) status = 0
    else status = 2 # magic value - failures during cal

    # Clean up
clean:
    scanfile1 = ""
    scanfile2 = ""
    delete (tmpin, verify-, >& "dev$null") 
    delete (tmpphu, verify-, >& "dev$null") 
    delete (tmpinlist, verify-, >& "dev$null") 
    delete (tmpoutlist, verify-, >& "dev$null") 
    delete (tmpfile1, verify-, >& "dev$null") 
    imdelete (tmpref, verify-, >& "dev$null")

    printlog ("", l_logfile, l_verbose) 
    if (status == 0) 
        printlog ("NSWHELPER Exit status good", l_logfile, l_verbose) 
    if (status == 2) 
        printlog ("WARNING - NSWHELPER Some calibration failed", \
            l_logfile, verbose+) 
    printlog ("---------------------------------------------------\
        ----------------------------", l_logfile, l_verbose) 

end


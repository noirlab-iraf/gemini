# Copyright(c) 2001-2015 Association of Universities for Research in Astronomy, Inc.

procedure nswavelength (lampspectra) 

# Wavelength calibration of GNIRS/NIRI arclamp spectra

char    lampspectra    {prompt = "Input GNIRS/NIRI lamp spectra"}                        # OLDP-1-primary-single-prefix=w
char    outspectra     {"", prompt = "Output spectra"}                                   # OLDP-1-output
char    outprefix      {"w", prompt = "Prefix for output spectra"}                       # OLDP-4
real    crval          {INDEF, prompt = "Approximate wavelength at coordinate reference pixel"} # OLDP-3
real    cdelt          {INDEF, prompt = "Approximate dispersion"}                        # OLDP-3
real    crpix          {INDEF, prompt = "Coordinate reference pixel"}                    # OLDP-3
int     dispaxis       {1, min = 1, max = 2, prompt = "Dispersion axis if not defined in the header"} # OLDP-3

char    database       {"", prompt = "Directory for files containing feature data"}      # OLDP-3
char    coordlist      {"gnirs$data/lowresargon.dat", prompt = "User coordinate list, line list"} # OLDP-3
bool    fl_inter       {no, prompt = "Examine identifications interactively"}            # OLDP-4
char    nsappwavedb    {"gnirs$data/nsappwave.fits", prompt = "Database for nsappwave"}  # OLDP-2
bool    fl_median      {yes, prompt = "Median filter XD arc (instead of transform)?"}    # OLDP-2
char    sdist          {"", prompt = "Aperture file for XD data (for nssdist)"}          # OLDP-2
int     sdorder        {4, min = 1, prompt = "Order of nssdist fitting function"}        # OLDP-2
int     xorder         {2, min = 1, prompt = "X order of nsfitcoords fitting function"}  # OLDP-2
int     yorder         {2, min = 1, prompt = "Y order of nsfitcoords fitting function"}  # OLDP-2
char    aptable        {"gnirs$data/apertures.fits", prompt = "Table of aperture data for nssdist"} # OLDP-3

char    section        {"default", prompt = "Image section for running identify"}        # OLDP-2
int     nsum           {10, min = 1, prompt = "Number of lines or columns to sum"}       # OLDP-2
char    ftype          {"emission", min = "emission|absorption", prompt = "Feature type"} # OLDP-2
real    fwidth         {4., min = 2, prompt = "Feature width in pixels"}                 # OLDP-2
real    cradius        {5., min = 2, prompt = "Centering radius in pixels"}              # OLDP-2
real    threshold      {100., prompt = "Feature threshold for centering"}                # OLDP-2
real    minsep         {2., prompt = "Minimum pixel separation for features"}            # OLDP-2
real    match          {-6., prompt = "Coordinate list matching limit, <0 pixels, >0 user"} # OLDP-2
char    function       {"chebyshev", min = "legendre|chebyshev|spline1|spline3", prompt = "Coordinate fitting function"} # OLDP-2
int     order          {4, min = 1, prompt = "Order of coordinate fitting function"}     # OLDP-2
char    sample         {"*", prompt = "Coordinate sample regions"}                       # OLDP-2
int     niterate       {10, min = 0, prompt = "Rejection iterations"}                    # OLDP-2
real    low_reject     {3., min = 0, prompt = "Lower rejection sigma"}                   # OLDP-2
real    high_reject    {3., min = 0, prompt = "Upper rejection sigma"}                   # OLDP-2
real    grow           {0., min = 0, prompt = "Rejection growing radius"}                # OLDP-2
bool    refit          {yes, prompt = "Refit coordinate function when running reidentify"} # OLDP-2
int     step           {10, prompt = "Steps in lines or columns for reidentification"}   # OLDP-2
bool    trace          {no, prompt = "Use fit from previous step rather than central aperture"} # OLDP-2
int     nlost          {3, min = 0, prompt = "Maximum number of lost features"}          # OLDP-2
bool    fl_overwrite   {yes, prompt = "Overwrite existing database entries"}             # OLDP-2
char    aiddebug       {"", prompt = "Debug parameter for aidpars"}                      # OLDP-3
real    fmatch         {0.2, prompt = "Acceptable fraction of unmatched reference lines?"} # OLDP-3
int     nfound         {6, prompt = "Miniumum number of identified lines required?"}     # OLDP-3
real    sigma          {0.05, prompt = "Sigma of line centering (pixels)"}               # OLDP-3
real    rms            {0.1, prompt = "RMS goal (fwidths)"}                              # OLDP-3

char    logfile        {"", prompt = "Logfile"}                                          # OLDP-1
bool    verbose        {yes, prompt = "Verbose"}                                         # OLDP-2
bool    debug          {no, prompt = "Very verbose"}                                     # OLDP-3
int     status         {0, prompt = "Exit status (0=Good, 1=Fail, 2=Some bad cal)"}      # OLDP-4

struct    *scanfile    {prompt = "Internal use only"}                                    # OLDP-4

begin

    char    l_inimages = ""
    char    l_outspectra = ""
    char    l_outprefix = ""
    real    l_crval
    real    l_cdelt
    real    l_crpix
    int     l_dispaxis
    char    l_database = ""
    char    l_coordlist = ""
    bool    l_fl_inter
    bool    l_fl_median
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
    char    l_aptable = ""
    bool    l_fl_overwrite
    char    l_sdist = ""
    int     l_sdorder
    int     l_xorder
    int     l_yorder
    char    l_logfile = ""
    bool    l_verbose
    bool    l_debug

    char    l_sci_ext = ""
    char    l_key_slit = ""
    char    l_key_dispaxis = ""

    int     junk, nbad, naxis, nver, version, nxx, nyy, nspatial
    bool    intdbg, haveifu, havexd, havemos, xd_is_straight
    struct  s_date, line
    char    badhdr, inphu, sdistphu, l_reference, secn, tmpstripe, tmpphu, imgin
    char    keyfound, tmpsci, tmpmedian, tmpout, key_value

    # Initialize parameters
    status = 1
    intdbg = no
    haveifu = no
    havexd = no
    havemos = no
    l_reference = ""

    tmpstripe = mktemp ("tmpstripe")
    tmpphu = mktemp ("tmpphu")
    tmpsci = mktemp ("tmpsci")
    tmpout = mktemp ("tmpout")

    # Cache parameter files
    cache ("gemextn", "nswhelper", "keypar", "nssdist", "nsfitcoords", \
        "nstransform") 

    # Set the local variables

    junk = fscan (lampspectra, l_inimages)
    junk = fscan (outspectra, l_outspectra)
    junk = fscan (outprefix, l_outprefix)
    l_crval = crval
    l_cdelt = cdelt
    l_crpix = crpix
    l_dispaxis = dispaxis
    junk = fscan (database, l_database)
    junk = fscan (coordlist, l_coordlist)
    l_fl_inter = fl_inter
    l_fl_median = fl_median
    junk = fscan (section, l_section)
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
    l_fmatch = fmatch
    l_nfound = nfound
    l_sigma = sigma
    l_rms = rms
    junk = fscan (nsappwavedb, l_nsappwavedb)
    l_fl_overwrite = fl_overwrite
    junk = fscan (sdist, l_sdist)
    l_sdorder = sdorder
    l_xorder = xorder
    l_yorder = yorder
    junk = fscan (aptable, l_aptable)
    junk = fscan (logfile, l_logfile)
    l_verbose = verbose
    l_debug = debug

    badhdr = ""
    junk = fscan (nsheaders.sci_ext, l_sci_ext)
    if ("" == l_sci_ext) badhdr = badhdr + " sci_ext"
    junk = fscan (nsheaders.key_slit, l_key_slit)
    if ("" == l_key_slit) badhdr = badhdr + " key_slit"
    junk = fscan (nsheaders.key_dispaxis, l_key_dispaxis)
    if ("" == l_key_dispaxis) badhdr = badhdr + " key_dispaxis"

    if (l_logfile == "") {
        junk = fscan (gnirs.logfile, l_logfile)
        if (l_logfile == "") {
            l_logfile = "gnirs.log"
            printlog ("WARNING - NSWAVELENGTH: Both nswavelength.logfile \
                and gnirs.logfile", l_logfile, verbose+) 
            printlog ("                         are empty. Using "\
                // l_logfile, l_logfile, verbose+) 
        }
    }

    # Write to the logfile

    printlog ("--------------------------------------------------------\
    -----------------------", l_logfile, l_verbose) 
    date | scan (s_date) 
    printlog ("NSWAVELENGTH -- " // s_date, l_logfile, l_verbose) 
    printlog ("", l_logfile, l_verbose) 

    # Check for *PREPARE
    gemextn (inimages=l_inimages, check="exists,mef", process="expand", \
        index="0", extname="", extversion="", ikparams="", omit="", \
        replace="", outfile=tmpphu, logfile="dev$null", glogpars="", 
        verbose-)
    nbad = gemextn.fail_count
    scanfile = tmpphu
    while (fscan (scanfile, imgin) != EOF) {
        keyfound=""
        hselect(imgin, "*PREPAR*", yes) | scan(keyfound)
        if (keyfound == "") {
            printlog ("ERROR - NSWAVELENGTH: Image " // imgin \
                // " not PREPAREd.", l_logfile, verbose+)
            nbad += 1
        }
    }
    if (nbad > 0) {
        printlog ("ERROR - NSWAVELENGTH: " // nbad // " image(s) \
            have not been run through *PREPARE", l_logfile, verbose+) 
        goto clean
    }

    # Check for 2D
    gemextn (inimages=l_inimages, check="exists,mef", process="expand", \
        index="", extname=l_sci_ext, extversion="1-", ikparams="", \
        omit="", replace="", outfile=tmpsci, logfile="dev$null",
        glogpars="", verbose-)
    nbad = gemextn.fail_count
    scanfile = tmpsci
    while (fscan (scanfile, imgin) != EOF) {
        hselect (imgin, "NAXIS", yes) | scan (naxis)
        if (nscan() != 1) {
            printlog ("WARNING - NSWAVELENGTH: Data in " // imgin \
                // " missing NAXIS.", l_logfile, verbose+)
            nbad = nbad + 1
        } else if (naxis != 2) {
            printlog ("WARNING - NSWAVELENGTH: Data in " // imgin \
                // " are not 2-D.", l_logfile, verbose+)
            nbad = nbad + 1
        }
    }
    if (nbad > 0) {
        printlog ("ERROR - NSWAVELENGTH: " // nbad // " image(s) \
            are not 2-D.", l_logfile, verbose+) 
        goto clean
    }

    nbad = 0
    if (l_outspectra != "") {
        # If outspectra is supplied by the user, it takes precedence over
        # outprefix, so check if outspectra already exists
        gemextn (l_outspectra, check="absent", process="none", index="", \
            extname="", extversion="", ikparams="", omit="extension, kernel", \
            replace="", outfile="dev$null", logfile="", glogpars="", verbose-)
        if (gemextn.fail_count != 0) {
            nbad += 1
        }
    } else {
        # Check if "outprefix + lampspectra" already exists
        gemextn ("%^%" // l_outprefix // "%" // "@" // tmpphu,
            check="absent", process="none", index="", extname="", \
            extversion="", ikparams="", omit="extension, kernel", replace="", \
            outfile="dev$null", logfile="", glogpars="", verbose-)
        if (gemextn.fail_count != 0) {
            nbad += 1
        }
    }

    if (nbad > 0) {
        printlog ("ERROR - NSWAVELENGTH: Existing or incorrectly formatted \
            output files.", l_logfile, verbose+) 
        goto clean
    }

    # Get example input file phu
    gemextn (l_inimages, check="exists,mef", process="expand", index="0", \
        extname="", extversion="", ikparams="", omit="", replace="", \
        outfile="STDOUT", logfile="dev$null", glogpars="", verbose-) | \
        scan (inphu)
    if (gemextn.fail_count != 0 || gemextn.count == 0) {
        printlog ("ERROR - NSWAVELENGTH: no PHU in input.", l_logfile, \
            verbose+)
        goto clean
    }

    # Detect IFU or MOS
    keypar (inphu, l_key_slit, silent+)
    if (keypar.found) {
        key_value = keypar.value
        if (substr(key_value, 1, 3) == "mos") {
            havemos = yes
        }
        if (strstr("IFU", strupr(key_value)) > 0) {
            haveifu = yes
        }
    }
    if (intdbg) print (l_key_slit // ": " // key_value)
    if (intdbg) print ("ifu: " // haveifu)
    if (intdbg) print ("mos: " // havemos)

    # Detect XD (multiple extensions, but not IFU or MOS)
    gemextn (inimages=l_inimages, check="exists,mef", process="none", \
        index="", extname="", extversion="", ikparams="", omit="", \
        replace="", outfile="STDOUT", logfile="dev$null", glogpars="",
        verbose-) | scan (line)
    gemextn (inimages=line, check="exists,mef", process="expand", \
        index="", extname=l_sci_ext, extversion="1-", ikparams="", \
        omit="", replace="", outfile="dev$null", logfile="dev$null",
        glogpars="", verbose-)
    havexd = gemextn.count > 1 && haveifu == no && havemos == no
    if (intdbg) print ("xd: " // line // ": " // havexd)
    
    # If the XD data have already had the s-distortion correction applied,
    # treat the data as longslit
    xd_is_straight = no
    if (havexd) {
        keypar (line // "[0]", "NSTRANSF", silent+)
        if (keypar.found) {
            havexd = no
            xd_is_straight = yes
        }
    }

    # If longslit, straightened XD or IFU, call nswhelper directly
    if (havexd == no && havemos == no) {
        if (haveifu) {
            printlog ("NSWAVELENGTH: IFU data detected. Calling \
                nswhelper to process data.", l_logfile, l_verbose) 
        } else if (xd_is_straight) {
            printlog ("NSWAVELENGTH: Straightened XD data detected. Calling \
                nswhelper to process data.", l_logfile, l_verbose) 
        } else {
            printlog ("NSWAVELENGTH: LS data detected. Calling \
                nswhelper to process data.", l_logfile, l_verbose) 
        }
        if (intdbg) print ("calling nswhelper with " // l_nsappwavedb)

        nswhelper (lampspectra=l_inimages, outspectra=l_outspectra, \
            outprefix=l_outprefix, reference=l_reference, crval=l_crval, \
            cdelt=l_cdelt, crpix=l_crpix, dispaxis=l_dispaxis, \
            database=l_database, coordlist=l_coordlist, \
            fl_inter=l_fl_inter, section=l_section, nsum=l_nsum, \
            ftype=l_ftype, fwidth=l_fwidth, cradius=l_cradius, \
            threshold=l_threshold, minsep=l_minsep, match=l_match, \
            function=l_function, order=l_order, sample=l_sample, \
            niterate=l_niterate, low_reject=l_low_reject, \
            high_reject=l_high_reject, grow=l_grow, refit=l_refit, \
            step=l_step, trace=l_trace, nlost=l_nlost, \
            aiddebug=l_aiddebug, fmatch=l_fmatch, nfound=l_nfound, \
            sigma=l_sigma, rms=l_rms, nsappwavedb=l_nsappwavedb, \
            fl_overwrite=l_fl_overwrite, logfile=l_logfile, \
            verbose=l_verbose, debug=l_debug)

        if (nswhelper.status != 0)
            goto clean

    } else {

        # Process XD or MOS data
        printlog ("NSWAVELENGTH: XD / MOS data detected.", l_logfile, \
            l_verbose)

        # If fl_median+, don't need to mess with sdist or transform
        if (l_fl_median) {
            # First lamp image used as reference.
            gemextn (inimages=l_inimages, check="exists,mef", process="none", \
                index="", extname="", extversion="", ikparams="", \
                omit="extension,index,kernel,section", replace="", \
                outfile="STDOUT", logfile="dev$null", verbose-) | \
                scan (l_reference)
                
            if (intdbg) print ("Medianing XD data, reference = " // l_reference)
            # Copy reference image to tmpstripe and set dispaxis
            copy (l_reference // ".fits", tmpstripe // ".fits", verbose-)
            keypar (tmpstripe // "[0]", l_key_dispaxis, silent+)
            if (keypar.found)
                l_dispaxis = int (keypar.value)
            # (else uses l_dispaxis value from input parameter)

            # Get number of science extensions
            gemextn (l_reference, check="exists", process="expand", index="", \
                extname=l_sci_ext, extversion="1-", ikparams="", omit="", \
                replace="", outfile="dev$null", logfile="", verbose-)
            if ((gemextn.fail_count != 0) || (gemextn.count == 0)) {
                printlog ("ERROR - NSWAVELENGTH: Bad science data in " //\
                    l_reference// ".", l_logfile, verbose+)
                goto clean
            }
            # Assume they're 1..n 
            nver = int (gemextn.count)

            # Median each science extension
            for (version = 1; version <= nver; version = version + 1) {
                tmpmedian = mktemp ("tmpmedian")
                imcopy (tmpstripe // "[" // l_sci_ext // \
                    "," // version // "]", tmpmedian, >& "dev$null")

                # Get size of spatial axis
                hselect (tmpmedian, "i_naxis" // (3 - l_dispaxis), yes) \
                    | scan (nspatial)

                if (l_dispaxis == 1) {
                    nxx = 1
                    nyy = nspatial
                } else {
                    nxx = nspatial
                    nyy = 1
                }

                if (intdbg)
                    print ("version " // version // ": " // "striping " // \
                        nxx // "," // nyy)
                median (tmpmedian, tmpmedian, nxx, nyy, boundary = "wrap", 
                    zlo = 0.001, zhi = INDEF, constant = 0, >& "dev$null")
                imcopy (tmpmedian, 
                    tmpstripe//"["//l_sci_ext//","//version//",overwrite]",
                    >& "dev$null")
                imdelete (tmpmedian, verify-, >& "dev$null")
            }

            # Now tmpstripe contains arc with median value across full extent 
            # of the spatial direction, but without rectifying; for this to 
            # work assumes that more than 1/2 the pixels in the x direction are
            # illuminated proceed to nswhelper

        } else {

            # Verify input data

            gemextn (l_sdist, check="exists,mef", process="expand", \
                index="0", extname="", extversion="", ikparams="", omit="", \
                replace="", outfile="STDOUT", logfile="", glogpars="", \
                verbose-) | scan (sdistphu)
            if (gemextn.fail_count != 0 || gemextn.count != 1) {
                printlog ("ERROR - NSWAVELENGTH: Bad sdist file.", l_logfile, \
                    verbose+)
                goto clean
            }

            keypar (sdistphu, "NSSDIST")
            if (keypar.found) {
                printlog ("NSWAVELENGTH: " // l_sdist \
                    // " already processed by nssdist", l_logfile, l_verbose) 
            } else {
                printlog ("NSWAVELENGTH: Tracing reference aperture.", \
                    l_logfile, l_verbose)

                # Trace

                keypar (sdistphu, l_key_dispaxis) 
                if (keypar.found)
                    l_dispaxis = int (keypar.value) 

                if (l_dispaxis == 1)
                    secn = "first column"
                else
                    secn = "first line"

                nssdist (l_sdist, outsuffix="_sdist", pixscale=1, \
                    database=l_database, firstcoord=0, coordlist="", \
                    aptable=l_aptable, fl_inter=l_fl_inter, fl_dbwrite=yes, \
                    section=secn, nsum=l_nsum, ftype="emission", \
                    fwidth=10, cradius=10, threshold=50, minsep=10, match=-6, \
                    function="chebyshev", order=l_sdorder, sample="", \
                    niterate=3, low_reject=5, high_reject=5, grow=0, refit+, \
                    step=10, trace-, nlost=10, aiddebug="", \
                    logfile=l_logfile, verbose=l_verbose)
                if (nssdist.status != 0)
                    goto clean

                gemextn (inimages=l_sdist, check="exists", process="none", \
                    index="", extname="", extversion="", ikparams="", \
                    omit="index,kernel,section", replace="", \
                    outfile="STDOUT", logfile="", glogpars="",
                    verbose-) | scan (l_sdist)
            }

            # Generate reference stripes

            # First image of input list used as reference.
            gemextn (inimages=l_inimages, check="exists,mef", \
                process="none", index="", extname="", extversion="", \
                ikparams="", omit="index,kernel,section", replace="", \
                outfile="STDOUT", logfile="dev$null", glogpars="", \
                verbose-) | scan (l_reference)

            keypar (inphu, "NSFITCOO", silent+)
            if (keypar.found) {
                printlog ("NSWAVELENGTH: input images already processed by \
                    nsfitcoords", l_logfile, l_verbose)
                copy (l_reference, tmpout // ".fits", verbose-, >& "dev$null")

            } else {
                printlog ("NSWAVELENGTH: Determine the s-distortion \
                    transformation for the arc.", l_logfile, l_verbose)

                nsfitcoords (inimages=l_reference, outspectra=tmpout, \
                    outprefix="", lamptransf="", sdisttransf=l_sdist, \
                    dispaxis=l_dispaxis, database=l_database, \
                    fl_inter=l_fl_inter, fl_align-, function="chebyshev", \
                    lxorder=l_xorder, lyorder=l_yorder, sxorder=l_xorder, \
                    syorder=l_yorder, pixscale=1, logfile=l_logfile, \
                    verbose=l_verbose, debug=no, force=no)
                if (nsfitcoords.status != 0)
                    goto clean
            }

            keypar (inphu, "NSTRANSF", silent+)
            if (keypar.found) {
                printlog ("NSWAVELENGTH: input images already processed by \
                    nstransform", l_logfile, l_verbose)
                copy (tmpout // ".fits", tmpstripe // ".fits", verbose-, \
                    >& "dev$null")

            } else {
                printlog ("NSWAVELENGTH: Generating s-distortion corrected \
                    arc.", l_logfile, l_verbose)

                nstransform (inimages=tmpout, outspectra=tmpstripe, \
                    dispaxis=l_dispaxis, database=l_database, \
                    fl_stripe+, interptype="poly3", xlog-, ylog-,\
                    pixscale=1, logfile=l_logfile, verbose=l_verbose, debug=no)
                if (nstransform.status != 0)
                    goto clean
            }
        }

        # Wavelength cal reference

        printlog ("NSWAVELENGTH: Calibrating straightened arc.", \
            l_logfile, l_verbose)

        nswhelper (lampspectra=tmpstripe, outspectra="", \
            outprefix="w", reference="", crval=l_crval, \
            cdelt=l_cdelt, crpix=l_crpix, dispaxis=l_dispaxis, \
            database=l_database, coordlist=l_coordlist, \
            fl_inter=l_fl_inter, section=l_section, nsum=l_nsum, \
            ftype=l_ftype, fwidth=l_fwidth, cradius=l_cradius, \
            threshold=l_threshold, minsep=l_minsep, match=l_match, \
            function=l_function, order=l_order, sample=l_sample, \
            niterate=l_niterate, low_reject=l_low_reject, \
            high_reject=l_high_reject, grow=l_grow, refit=l_refit, \
            step=l_step, trace=l_trace, nlost=l_nlost, \
            aiddebug=l_aiddebug, fmatch=l_fmatch, nfound=l_nfound, \
            sigma=l_sigma, rms=l_rms, nsappwavedb=l_nsappwavedb, \
            fl_overwrite=l_fl_overwrite, logfile=l_logfile, \
            verbose=l_verbose, debug=l_debug)

        if (nswhelper.status != 0) {
            if (nswhelper.status == 2) {
                printlog ("WARNING - NSWAVELENGTH: Reference is not \
                    complete (continuing).", l_logfile, verbose+)
            } else 
                goto clean
        }

        # Wavelength cal data

        printlog ("NSWAVELENGTH: Calibrating the input arcs and creating \
            appropriate database", \
            l_logfile, l_verbose)
        printlog ("              files using the fit determined from the \
            straightened arc.", l_logfile, l_verbose)

        # This call to nswhelper should always be non-interactive and the
        # coordinate function should not be refit. If the reference parameter
        # in nswhelper is set to the straightened arc and the refit parameter
        # in nswhelper is set to no, the fit determined from the straightened
        # arc is copied to the database files associated with the
        # un-straightened arc (i.e., the input to nswavelength).
        nswhelper (lampspectra=l_inimages, outspectra=l_outspectra, \
            outprefix=l_outprefix, reference="w" // tmpstripe, crval=l_crval, \
            cdelt=l_cdelt, crpix=l_crpix, dispaxis=l_dispaxis, \
            database=l_database, coordlist=l_coordlist, fl_inter-, \
            section=l_section, nsum=l_nsum, ftype=l_ftype, fwidth=l_fwidth, \
            cradius=l_cradius, threshold=l_threshold, minsep=l_minsep, \
            match=l_match, function=l_function, order=l_order, \
            sample=l_sample, niterate=l_niterate, low_reject=l_low_reject, \
            high_reject=l_high_reject, grow=l_grow, refit-, step=l_step, \
            trace=l_trace, nlost=l_nlost, aiddebug=l_aiddebug, \
            fmatch=l_fmatch, nfound=l_nfound, sigma=l_sigma, rms=l_rms, \
            nsappwavedb=l_nsappwavedb, fl_overwrite=l_fl_overwrite, \
            logfile=l_logfile, verbose=l_verbose, debug=l_debug)
        status = nswhelper.status    # Special value 2 possible
        goto clean                   # Includes exit on success
    }

    status = 0

    # Clean up
clean:
    imdelete (tmpstripe // ", w" // tmpstripe // "," // tmpout, verify-, \
        >& "dev$null")
    delete (tmpphu // "," // tmpsci, verify-, >& "dev$null")
    scanfile = ""

    printlog ("", l_logfile, l_verbose) 
    if (status == 0) 
        printlog ("NSWAVELENGTH Exit status good", l_logfile, l_verbose) 
        printlog ("---------------------------------------------------------\
            ----------------------", l_logfile, l_verbose) 

end

# Copyright(c) 2001-2012 Association of Universities for Research in Astronomy, Inc.

procedure nstelluric (inimages, cal)

# Written May 25 2000: Marianne Takamiya
#         Jun 16 2000: MT major debug
#         Jul 12 2000: MT allow input and cal files to have different 
#                      extension names
#         14 Oct 2002, major rewrite and cleanup, JJ
#         19 Aug 2003, KL  IRAF2.12 - new parameter, addonly, in hedit
#         29 Oct 2003, KL  moved from niri to gnirs package
#          

char    inimages    {prompt = "Input spectra"}                              # OLDP-1-primary-single-prefix=a
char    cal         {prompt = "Telluric calibration spectrum"}              # OLDP-1-input
char    extract     {"", prompt = "Optional extracted spectrum"}            # OLDP-3
char    outspectra  {"", prompt = "Output spectra"}                         # OLDP-1-output
char    outprefix   {"a", prompt = "Output prefix"}                         # OLDP-4
int     dispaxis    {1, min=1, max=2, prompt="Dispersion axis if not present in header"}    # OLDP-3
char    fitfunction {"spline3", prompt = "Fitting function for normalizing (fit1d)"}    # OLDP-2
int     fitorder    {10, prompt="Order of fitting function (fit1d)"}        # OLDP-2
bool    fl_xcorr    {yes, prompt="Cross correlate for shift? (telluric)"}   # OLDP-3
bool    fl_tweakrms {yes, prompt="Tweak to minimize RMS? (telluric)"}       # OLDP-3
char    sample      {"*", prompt = "Sample ranges for fits (fit1d, telluric)"}  # OLDP-3
real    threshold   {0., prompt = "Threshold for calibration (telluric)"}   # OLDP-3
int     lag         {20, prompt="Cross correlation lag (pixels) (telluric)"}# OLDP-3
real    shift       {0., prompt = "Initial shift of calibration spectrum (pixels) (telluric)"}    # OLDP-3
real    scale       {1., prompt = "Initial scale factor multiplying airmass ratio (telluric)"}    # OLDP-3
real    dshift      {0.5, prompt = "Initial shift search step (telluric)"}  # OLDP-3
real    dscale      {0.1, prompt = "Initial scale factor search step (telluric)"}    # OLDP-3
real    low_reject  {1.5, prompt = "Low rejection in sigma of fit (fit1d)"} # OLDP-3
real    high_reject {3., prompt = "Low rejection in sigma of fit (fit1d)"}   # OLDP-3
int     niterate    {2, prompt = "Number of rejection iterations (fit1d)"}  # OLDP-3
real    grow        {1.0, prompt = "Rejection growing radius in pixels (fit1d)"}    # OLDP-3
int     ifuextn     {12, prompt = "Extension version to use for IFU"}       # OLDP-3
bool    fl_inter    {no, prompt = "Run task interactively?"}                # OLDP-4
char    database    {"", prompt = "Directory for database files for nsextract"} # OLDP-4
char    logfile     {"", prompt = "Logfile"}                                # OLDP-4
bool    verbose     {yes, prompt="Run task verbosely"}                      # OLDP-2
int     status      {0, prompt="Exit status of task"}                       # OLDP-4
struct  *scanfile1  {"", prompt = "Internal use only"}                      # OLDP-4
struct  *scanfile2  {"", prompt = "Internal use only"}                      # OLDP-4


begin
    char    l_inimages = ""
    char    l_cal = ""
    char    l_extract = ""
    char    l_outspectra = ""
    char    l_outprefix = ""
    int     l_dispaxis
    char    l_fitfunction = ""
    int     l_fitorder
    bool    l_fl_xcorr
    bool    l_fl_tweakrms
    char    l_sample = ""
    real    l_threshold
    int     l_lag
    real    l_shift
    real    l_scale
    real    l_dshift
    real    l_dscale
    real    l_low_reject
    real    l_high_reject
    int     l_niterate
    real    l_grow
    int     l_ifuextn
    bool    l_fl_inter
    char    l_database = ""
    char    l_logfile = ""
    bool    l_verbose

    char    l_sci_ext = ""
    char    l_var_ext = ""
    char    l_dq_ext = ""
    char    l_key_airmass = ""
    char    l_key_dispaxis = ""
    char    l_key_mode = ""

    char    imagesdel, filesdel, imgi, imgo, tmpimg0, caldel, imgin, phu
    int     ii, inputnum, outnum, junk, ref_dispaxis, zero
    real    scaler, shiftr
    char    shifts, scales, badhdr, sci, sciapp, var, varapp, dq, dqapp
    bool    madenewoutput, debug, haveifu, ifuproc
    real    sciairmass, calairmass, scaleexp
    struct  sdate
    int     nfiles, nbad, nversion, version, n, i
    char    airmass, keyfound
    char    tmpsci1d
    char    tmpimout, tmpout, tmpscinorm, tmpcalnorm
    char    tmpcalshift, tmpcalshiftvar, tmpcalshiftdq, tmpcalfit
    char    tmpcalscale, tmpcalscalevar, tmphead
    char    tmpinimages, tmpoutspectra, tmptwodsci, tmptwodvar
    char    basissci, basisvar
    bool    twod, haveextract


    junk = fscan (  inimages, l_inimages)
    junk = fscan (  cal, l_cal)
    junk = fscan (  extract, l_extract)
    junk = fscan (  outspectra, l_outspectra)
    junk = fscan (  outprefix, l_outprefix)
    ref_dispaxis =  dispaxis
    junk = fscan (  fitfunction, l_fitfunction)
    l_fitorder =    fitorder
    l_fl_xcorr =    fl_xcorr
    l_fl_tweakrms = fl_tweakrms
    junk = fscan (  sample, l_sample)
    l_threshold =   threshold
    l_lag =     lag
    l_shift =   shift
    l_scale =   scale
    l_dshift =  dshift
    l_dscale =  dscale
    l_low_reject =  low_reject
    l_high_reject = high_reject
    l_niterate =    niterate
    l_grow =        grow
    l_ifuextn =     ifuextn
    l_fl_inter =    fl_inter
    junk = fscan (  database, l_database)
    junk = fscan (  logfile, l_logfile)
    l_verbose =     verbose


    badhdr = ""
    junk = fscan (  nsheaders.sci_ext, l_sci_ext)
    if ("" == l_sci_ext) badhdr = badhdr + " sci_ext"
    junk = fscan (  nsheaders.var_ext, l_var_ext)
    if ("" == l_var_ext) badhdr = badhdr + " var_ext"
    junk = fscan (  nsheaders.dq_ext, l_dq_ext)
    if ("" == l_dq_ext) badhdr = badhdr + " dq_ext"
    junk = fscan (  nsheaders.key_dispaxis, l_key_dispaxis)
    if ("" == l_key_dispaxis) badhdr = badhdr + " key_dispaxis"
    junk = fscan (  nsheaders.key_airmass, l_key_airmass)
    if ("" == l_key_airmass) badhdr = badhdr + " key_airmass"
    junk = fscan (  nsheaders.key_mode, l_key_mode)
    if ("" == l_key_mode) badhdr = badhdr + " key_mode"


    status = 1
    cache ("gimverify", "gemextn", "keypar", "nsextract", "gemdate")
    imagesdel = ""
    filesdel = ""
    madenewoutput = no
    debug = no

    if (l_logfile == "") {
        junk = fscan (gnirs.logfile, l_logfile)
        if (l_logfile == "") {
            l_logfile = "gnirs.log"
            printlog ("WARNING - NSTELLURIC: Both nstelluric.logfile \
                and gnirs.logfile are", l_logfile, verbose+)
            printlog ("                      undefined.  Using " \
                // l_logfile, l_logfile, verbose+)
        }
    }

    if ("" != badhdr) {
        printlog ("ERROR - NSTELLURIC: Parameter(s) missing from \
            nsheaders: " // badhdr, l_logfile, verbose+) 
        goto clean
    }

    # define tmp file/image names
    tmpout = mktemp ("tmpout")

    tmphead = mktemp ("tmphead")
    tmpinimages = mktemp ("tmpinimages")
    tmpoutspectra = mktemp ("tmpoutspectra")
    filesdel = tmphead // "," // tmpinimages // "," // tmpoutspectra

    tmpcalnorm = mktemp ("tmpcalnorm")
    tmpcalfit = mktemp ("tmpcalfit")
    caldel = tmpcalnorm // "," // tmpcalfit

    tmpsci1d = mktemp ("tmpsci1d")
    tmptwodsci = mktemp ("tmptwodsci")
    tmptwodvar = mktemp ("tmptwodvar")



    date | scan (sdate)
    printlog ("-------------------------------------------------------\
        -------------------------", l_logfile, l_verbose)
    printlog ("NSTELLURIC -- " // sdate, l_logfile, l_verbose)
    printlog ("", l_logfile, l_verbose)
    printlog ("   inimages = " // l_inimages, l_logfile, l_verbose)
    printlog ("        cal = " // l_cal, l_logfile, l_verbose)
    printlog (" outspectra = " // l_outspectra, l_logfile, l_verbose)
    printlog ("    sci_ext = " // l_sci_ext, l_logfile, l_verbose)
    printlog ("    var_ext = " // l_var_ext, l_logfile, l_verbose)
    printlog ("     dq_ext = " // l_dq_ext, l_logfile, l_verbose)
    printlog ("", l_logfile, l_verbose)


    if (l_outspectra == "" && l_outprefix == "") {
        printlog ("ERROR - NSTELLURIC: parameters outspectra and \
            outprefix are empty", l_logfile, verbose+)
        goto clean
    }
    if (l_inimages == "") {
        printlog ("ERROR - NSTELLURIC: parameter inimages is empty.", \
            l_logfile, verbose+)
        goto clean
    }
    if (l_cal == "") {
        printlog ("ERROR - NSTELLURIC: parameter cal is empty.", \
            l_logfile, verbose+)
        goto clean
    }



    # Check input files

    gemextn (l_inimages, check="mef,exists", process="none", index="", \
        extname="", extversion="", ikparams="", omit="extension, kernel", \
        replace="", outfile=tmpinimages, logfile="", glogpars="",
        verbose=l_verbose)
    if (0 != gemextn.fail_count || 0 == gemextn.count) {
        printlog ("ERROR - NSTELLURIC:  Problems with input.", \
            l_logfile, verbose+)
        goto clean
    }
    nfiles = gemextn.count


    nbad = 0
    scanfile1 = tmpinimages
    while (fscan (scanfile1, imgin) != EOF) {
        phu = imgin // "[0]"
        keyfound = ""
        hselect (phu, "*PREPAR*", yes) | scan (keyfound)
        if (keyfound == "") {
            printlog ("ERROR - NSTELLURIC: Image " // imgin \
                // " not PREPAREd.", l_logfile, verbose+) 
            nbad += 1
        }
    }

    if (nbad > 0) {
        printlog ("ERROR - NSTELLURIC: " // nbad // " image(s) \
            have not been run through *PREPARE", l_logfile, verbose+) 
        goto clean
    }


    # Check output files

    gemextn (l_outspectra, check="absent", process="none", index="", \
        extname="", extversion="", ikparams="", omit="extension, kernel", \
        replace="", outfile=tmpoutspectra, logfile="", glogpars="",
        verbose=l_verbose)
    if (0 != gemextn.fail_count) {
        printlog ("ERROR - NSTELLURIC:  Problems with output.", \
            l_logfile, verbose+)
        goto clean
    }
    if (0 == gemextn.count) {
        gemextn ("%^%" // l_outprefix // "%" // "@" // tmpinimages, \
            check="absent", process="none", index="", \
            extname="", extversion="", ikparams="", omit="kernel,exten", \
            replace="", outfile=tmpoutspectra, logfile="", glogpars="",
            verbose=l_verbose)
        if (gemextn.fail_count != 0 || gemextn.count == 0) {
            printlog ("ERROR - NSTELLURIC: Bad or incorrectly formatted \
                output files", l_logfile, verbose+) 
            goto clean
        }
    }

    if (nfiles != gemextn.count) {
        printlog ("ERROR - NSTELLURIC:  Number of input and output \
            files doesn't match.",l_logfile, verbose+)
        goto clean
    }


    # Generate output files!

    scanfile1 = tmpinimages
    scanfile2 = tmpoutspectra
    while (fscan (scanfile1, imgin) != EOF) {
        junk = fscan (scanfile2, imgo)
        imcopy (imgin//"[0]", imgo, verbose-)
    }
    madenewoutput = yes


    # Check threshold parameter

    if (no == isindef (l_threshold) && l_threshold > 1) {
        printlog ("WARNING - NSTELLURIC: threshold parameter larger than \
            expected for normalised data (see help).", l_logfile, verbose+)
    }


    # Check cal spectrum

    gemextn (l_cal, check="exists,mef", process="none", index="", \
        extname="", extversion="", ikparams="", omit="extension, kernel", \
        replace="", outfile="STDOUT", logfile="", glogpars="",
        verbose-) | scan (l_cal)
    if (gemextn.count != 1 || gemextn.fail_count != 0) {
        printlog ("ERROR - NSTELLURIC: Calibration spectrum " // l_cal \
            // " not found.", l_logfile, verbose+)
        goto clean
    }
    keypar (l_cal // "[" // l_sci_ext // "]", "i_naxis")
    if (no == keypar.found || int (keypar.value) != 1) {
        printlog ("ERROR - NSTELLURIC: Input calibration spectrum is \
            not 1-dimensional.", l_logfile, verbose+)
        printlog ("                    Use NSEXTRACT to extract a 1-d \
            spectrum first.", l_logfile, verbose+)
        goto clean
    }
    keypar (l_cal // "[0]", l_key_airmass, silent+)
    if (no == keypar.found) keypar.value = "1.0"
    calairmass = real (keypar.value)


    # Number of extension versions

    gemextn (l_cal, check="", process="expand", index="", \
        extname=l_sci_ext, extversion="1-", ikparams="", \
        omit="extension, kernel", replace="", outfile="dev$null", \
        logfile="", glogpars="", verbose=l_verbose)
    nversion = gemextn.count
    if (nversion == 0) {
        printlog ("ERROR - NSTELLURIC: Input calibration spectrum has \
            no science data.", l_logfile, verbose+)
        goto clean
    }


    # Check optional extracted spectrum

    if ("" == l_extract) {
        haveextract = no
    } else {
        gemextn (l_extract, check="exists,mef", process="none", index="", \
            extname="", extversion="", ikparams="", omit="", \
            replace="", outfile="STDOUT", logfile="", glogpars="",
            verbose-) | scan (l_extract)
        if (gemextn.count != 1 || gemextn.fail_count != 0) {
            printlog ("ERROR - NSTELLURIC: Extracted spectrum " \
                // l_extract // " not found.", l_logfile, verbose+)
            goto clean
        }
        haveextract = yes
        keypar (l_extract // "[" // l_sci_ext // "]", "i_naxis")
        if (no == keypar.found || int (keypar.value) != 1) {
            printlog ("ERROR - NSTELLURIC: Extracted spectrum is \
                not 1-dimensional.", l_logfile, verbose+)
            printlog ("                    Use NSEXTRACT to extract a 1-d \
                spectrum first.", l_logfile, verbose+)
            goto clean
        }
        gemextn (l_extract, check="", process="expand", index="", \
            extname=l_sci_ext, extversion="1-", ikparams="", \
            omit="", replace="", outfile="dev$null", \
            logfile="", glogpars="", verbose-)
        if (no == (nversion == gemextn.count)) {
            printlog ("ERROR - NSTELLURIC: Extracted spectrum has \
                wrong number of science extensions.", l_logfile, verbose+)
            goto clean
        }
    }


    for (version = 1; version <= nversion; version = version + 1) {

        sci = "[" // l_sci_ext // "," // version // "]"
        sciapp = "[" // l_sci_ext // "," // version // ",append]"
        var = "[" // l_var_ext // "," // version // "]"
        varapp = "[" // l_var_ext // "," // version // ",append]"

        gemhedit (l_cal // sci, "AIRMASS", calairmass, "", delete-)


        # Normalize for input into telluric

        if (version == 1) {
            imcopy (l_cal // "[0]", tmpcalnorm // "[0, overwrite]", \
                verbose-, >& "dev$null")
        }

        l_dispaxis = ref_dispaxis
        keypar (l_cal // sci, l_key_dispaxis, silent+)
        if (keypar.found) l_dispaxis = int (keypar.value)

        if (fl_inter)
            printlog ("NSTELLURIC: Fitting standard for normalization.", \
                l_logfile, l_verbose)
        imdelete (tmpcalfit, verify-, >& "dev$null")
        fit1d (l_cal // sci, tmpcalfit, "fit",
            axis = l_dispaxis, interac = l_fl_inter, sample = l_sample, \
            naverage = 1, function = l_fitfunction, order = l_fitorder, \
            low_rej = l_low_reject, high_rej = l_high_reject, \
            grow = l_grow, niter = l_niterate, \
            cursor = "", graphics = "stdgraph")
        imexpr ("a/b", tmpcalnorm // sciapp, l_cal // sci, tmpcalfit, \
            outtype = "real", verbose-)
        imexpr ("a/(b*b)", tmpcalnorm // varapp, l_cal // var, \
            tmpcalfit, outtype = "real", verbose-)

    }


    # Run telluric. If dispaxis=1, spectrum runs along rows.

    shifts = "0"
    scales = "0"

    scanfile1 = tmpinimages
    scanfile2 = tmpoutspectra
    while (fscan (scanfile1, imgin) != EOF) {
        junk = fscan (scanfile2, imgo)

        # tmp images used within this loop
        delete (tmpsci1d, verify-, >& "dev$null")
        tmpsci1d = mktemp (imgin // "-")

        keypar (imgin // "[0]", l_key_mode, silent+)
        haveifu = keypar.found && ("IFU" == strupr (keypar.value))

        # Copy the airmass header keyword to the [SCI] extension

        keypar (imgin // "[0]", l_key_airmass, silent+)
        if (no == keypar.found) keypar.value = "1.0"
        gemhedit (imgin // "[" // l_sci_ext // "]", "AIRMASS", \
            real (keypar.value), "", delete-)

        if (1 == nversion) {
            sci = "[" // l_sci_ext // "]"
            sciapp = "[" // l_sci_ext // ",append]"
            var = "[" // l_var_ext // "]"
            varapp = "[" // l_var_ext // ",append]"
        } else {
            sci = "[" // l_sci_ext // ",1]"
            sciapp = "[" // l_sci_ext // ",1,append]"
            var = "[" // l_var_ext // ",1]"
            varapp = "[" // l_var_ext // ",1,append]"
        }

        # Extract data if necessary

        # if input is 2d, run nsextract; normalize input for telluric
        # changed to use local logfile (prev deleted tmp.log)

        if (debug) print ("naxis check: " // imgin // sci)
        keypar (imgin // sci, "i_naxis", silent+)
        if (debug) print ("naxis check: " // keypar.found)
        if (debug && keypar.found) print ("naxis check: " // keypar.value)

        twod = keypar.found && int (keypar.value) == 2

        if (twod) {

            if (haveextract) {
                printlog ("NSTELLURIC: Using " // l_extract \
                    // " as extracted spectrum.", l_logfile, l_verbose)
                imgi = l_extract
            } else {

                # This never worked before, as far as I can tell
                # (had to fix bugs present in original).

                printlog ("NSTELLURIC: Calling NSEXTRACT to extract 1D \
                    spectrum.", l_logfile, l_verbose)

                nsextract (imgin, outspectra=tmpsci1d, outprefix="", \
                    dispaxis=l_dispaxis, database=l_database, line=500, \
                    nsum=10, ylevel=0.1, upper=10, lower=-10, \
                    background="average", fl_vardq+, fl_addvar-, \
                    fl_skylines-, fl_inter-, fl_apall+, fl_trace-, \
                    aptable="gnirs$data/apertures.fits", fl_usetabap-, \
                    fl_flipped+, fl_project-, fl_findneg-, bgsample="*", \
                    trace="", tr_nsum=10, tr_step=10, tr_nlost=3, \
                    tr_function="legendre", tr_order=5, tr_sample="*", \
                    tr_naver=1, tr_niter=0, tr_lowrej=3., tr_highrej=3., \
                    tr_grow=0, weights="variance", logfile=l_logfile, \
                    verbose=l_verbose)

                if (nsextract.status != 0) {
                    printlog ("ERROR - NSTELLURIC: Error in nsextract.", \
                        l_logfile, verbose+)
                    goto clean
                }

                imgi = tmpsci1d
            }
        } else {
            imgi = imgin
        }


        # Add an extra iteration for IFU to do single measurement

        if (haveifu) zero = 0
        else         zero = 1

        for (i = zero; i <= nversion; i = i + 1) {

            ifuproc = (haveifu && 0 == i)
            if (ifuproc) {
                gemextn (imgi, check="exists", process="expand", \
                    index="", extname=l_sci_ext, extversion=l_ifuext, \
                    ikparams="", omit="", replace="", outfile="dev$null", \
                    logfile="", glogpars="", verbose-)
                if (0 == gemextn.count) {
                    printlog ("ERROR - NSTELLURIC: ifuextn not available \
                        in " // imgi, l_logfile, verbose+)
                    goto clean
                }
                version = l_ifuextn
                printlog ("NSTELLURIC: Using extension " // l_ifuextn \
                    // " of " // imgi // " for reference.", \
                    l_logfile, verbose+)

                keyfound = ""
                hselect (imgo // "[0]", "NSTRANSF", yes) | scan (keyfound)
                if (keyfound == "") {
                    printlog ("WARNING - NSTELLURIC: Image " // imgo \
                        // " not NSTRANSFORMed.", l_logfile, verbose+)
                }
            } else {
                version = i
            }

            # tmp images used within this loop
            if (ifuproc || no == haveifu) {
                tmpimg0 = mktemp ("tmpimg0")
                tmpimout = mktemp ("tmpimout")
                tmpscinorm = mktemp ("tmpscinorm")
                tmpcalshift = mktemp ("tmpcalshift")
                tmpcalshiftvar = mktemp ("tmpcalshiftvar")
                tmpcalshiftdq = mktemp ("tmpcalshiftdq")
                tmpcalscale = mktemp ("tmpcalscale")
                tmpcalscalevar = mktemp ("tmpcalscalevar")

                imagesdel = tmpimg0 // "," // tmpimout // "," \
                    // tmpcalshift // "," \
                    // tmpscinorm // "," // tmpcalshift // "," \
                    // tmpcalshiftvar // "," // tmpcalshiftdq // "," \
                    // tmpcalscale // "," // tmpcalscalevar
            }

            if (1 == nversion) {
                sci = "[" // l_sci_ext // "]"
                sciapp = "[" // l_sci_ext // ",append]"
                var = "[" // l_var_ext // "]"
                varapp = "[" // l_var_ext // ",append]"
                var = "[" // l_var_ext // "]"
                varapp = "[" // l_var_ext // ",append]"
                dq = "[" // l_dq_ext // "]"
                dqapp = "[" // l_dq_ext // ",append]"
            } else {
                sci = "[" // l_sci_ext // "," // version // "]"
                sciapp = "[" // l_sci_ext // "," // version // ",append]"
                var = "[" // l_var_ext // "," // version // "]"
                varapp = "[" // l_var_ext // "," // version // ",append]"
                dq = "[" // l_dq_ext // "," // version // "]"
                dqapp = "[" // l_dq_ext // "," // version // ",append]"
            }

            if (ifuproc || no == haveifu) {

                if (no == imaccess (imgi // sci)) {
                    printlog ("WARNING - NSTELLURIC: No data for " \
                        // imgi // " extension " // version, \
                        l_logfile, verbose+)
                    next
                }

                l_dispaxis = ref_dispaxis
                keypar (imgi // sci, l_key_dispaxis, silent+)
                if (keypar.found) l_dispaxis = int (keypar.value)

                if (fl_inter)
                    printlog ("NSTELLURIC: Fitting obs for \
                        normalization.", l_logfile, l_verbose)
                fit1d (imgi // sci, tmpscinorm, "ratio", \
                    axis = l_dispaxis, interac = l_fl_inter, \
                    sample = l_sample, naverage = 1, \
                    function = l_fitfunction, order = l_fitorder, \
                    low_rej = l_low_reject, high_rej = l_high_reject, \
                    grow = l_grow, niter = l_niterate, \
                    cursor = "", graphics = "stdgraph")

                keypar (imgin // sci, l_key_airmass, silent+)
                if (no == keypar.found)
                    keypar (imgin // "[0]", l_key_airmass, silent+)
                sciairmass = real (keypar.value)
                gemhedit (tmpscinorm, l_key_airmass, sciairmass, "", delete-)

                # Get shift and scale factor using telluric

                if (fl_inter)
                    printlog ("NSTELLURIC: Calling telluric.", \
                        l_logfile, l_verbose)
                telluric (tmpscinorm, tmpimout, tmpcalnorm // sci, \
                    ignorea = yes, xcorr = l_fl_xcorr, \
                    tweakrm = l_fl_tweakrms, \
                    interac = l_fl_inter, sample = l_sample, \
                    threshold = l_threshold, lag = l_lag, \
                    shift = l_shift, scale = l_scale, \
                    dshift = l_dshift, dscale = l_dscale, \
                    offset = 1., smooth = 1, cursor = "", > tmpout)

                match ("shift", tmpout, stop-, print_file_n-, metach-) | \
                    scan (shifts, shifts, shifts, shifts, scales, scales, \
                    scales)

                delete (tmpout, verify-)
                shiftr = real (shifts)
                scaler = real (scales)

                scaleexp = (sciairmass/calairmass) * scaler

                printlog ("Calibration spectrum shifted by " // shiftr \
                    // " pixels.", l_logfile, l_verbose)
                printlog ("Airmass ratio = " // (sciairmass/calairmass), \
                    l_logfile, l_verbose)
                printlog ("Exponent scale factor = " // scaler, \
                    l_logfile, l_verbose)

                l_dispaxis = ref_dispaxis
                keypar (imgin // sci, l_key_dispaxis, silent+)
                if (keypar.found) l_dispaxis = int (keypar.value)

                # Need to generate appropriate 2d calibration spectrum?
                if (twod) {

                    printlog ("NSTELLURIC: Expanding to 2D (please \
                        wait).", l_logfile, l_verbose)

                    imdelete (tmptwodsci, verify-, >& "dev$null")
                    imdelete (tmptwodvar, verify-, >& "dev$null")
                    tmptwodsci = mktemp ("tmptwodsci")
                    tmptwodvar = mktemp ("tmptwodvar")
                    imstack (tmpcalnorm // sci, tmptwodsci)
                    imstack (tmpcalnorm // var, tmptwodvar)
                    hselect (imgin // sci, "i_naxis"//(3-l_dispaxis), yes) \
                        | scan (n)
                    blkrep (tmptwodsci, tmptwodsci, 1, n)
                    blkrep (tmptwodvar, tmptwodvar, 1, n)
                    if (2 == l_dispaxis) {
                        imtranspose (tmptwodsci, tmptwodsci)
                        imtranspose (tmptwodvar, tmptwodvar)
                    }
                    basissci = tmptwodsci
                    basisvar = tmptwodvar
                } else {
                    basissci = tmpcalnorm // sci
                    basisvar = tmpcalnorm // var
                }

                # Shift and scale the cal spectrum as determined
                # by telluric.  In order to shift images by
                # fractional pixels, the image has to be real.

                if (l_dispaxis == 1) {
                    imshift (basissci, tmpcalshift, shiftr, 0., \
                        shifts_file = "", interp = "linear", \
                        boundary = "nearest", constant = 0.)
                        imshift (basisvar, tmpcalshiftvar, shiftr, 0., \
                        shifts_file = "", interp = "linear", \
                        boundary = "nearest", constant = 0.)
                } else {
                    imshift (basissci, tmpcalshift, 0., shiftr, \
                        shifts_file = "", interp = "linear", \
                        boundary = "nearest", constant = 0.)
                        imshift (basisvar, tmpcalshiftvar, 0., shiftr, \
                        shifts_file = "", interp = "linear", \
                        boundary = "nearest", constant = 0.)
                }


                imexpr("a>0 ? max(a,0.00001)**b : 1.0", \
                    tmpcalscale, tmpcalshift, \
                    scaleexp, outtype="real",verbose-)
                    imexpr("a>0 ? (max(a,0.00001)**(2*(b-1)))*b*c : 1.0", \
                    tmpcalscalevar, tmpcalshift, scaleexp, \
                    tmpcalshiftvar, outtype="real", verbose-)
            }


            if (no == ifuproc) {

                printlog ("NSTELLURIC: Processing extension " \
                    // version, l_logfile, l_verbose)

                # Make the output spectrum (except when doing 
                # initial IFU processing)
                imexpr ("a/b", imgo // sciapp, imgin // sci, tmpcalscale, \
                    outtype = "real", verbose-)
                imexpr ("(c/(b*b))+(a*a*d/(b*b*b*b))", \
                    imgo // varapp, imgin // sci, tmpcalscale, \
                    imgin // var, tmpcalscalevar, outtype = "real", \
                    verbose-)
                imcopy (imgin // dq, imgo // dqapp, verbose-)

            }

            # delete temp files
            if (no == haveifu || i == nversion)
                imdelete (imagesdel, verify-, >& "dev$null")
                
        } # end for-loop


        # Add header info (needs much more!)
        gemdate ()
        gemhedit (imgo // "[0]", "NSTELLUR", gemdate.outdate, \
            "UT Time stamp for NSTELLURIC", delete-)
        gemhedit (imgo // "[0]", "GEM-TLM", gemdate.outdate, \
            "UT Last modification with GEMINI", delete-)
        printf ("%-8s= \'%-18s\' / %-s\n", "NSTELCAL", l_cal, \
            "Calibration file for NSTELLURIC", >> tmphead)
        printf ("%-8s= %20.3f / %-s\n", "NSTELSHF", shiftr, \
            "Applied shift (pixels)", >> tmphead)
        printf ("%-8s= %20.3f / %-s\n", "NSTELSCL", scaler, \
            "Scale factor for airmass ratio", >> tmphead)
        printf ("%-8s= \'%-18s\' / %-s\n", "NSTELFUN", l_fitfunction, \
            "Normalization fitting function", >> tmphead)
        printf ("%-8s= %20.0f / %-s\n", "NSTELORD", l_fitorder, \
            "Normalization fitting function order", >> tmphead)
        printf ("%-8s= \'%-18b\' / %-s\n", "NSTELXCO", l_fl_xcorr, \
            "Use cross-correlation for shift?", >> tmphead)
        printf ("%-8s= \'%-18b\' / %-s\n", "NSTELTWK", l_fl_tweakrms, \
            "Tweak to minimize RMS?", >> tmphead)
        printf ("%-8s= \'%-18s\' / %-s\n", "NSTELSAM", l_sample, \
            "Sample ranges for fits", >> tmphead)
        printf ("%-8s= %20.2f / %-s\n", "NSTELTHR", l_threshold, \
            "Threshold for calibration", >> tmphead)
        printf ("%-8s= %20.0f / %-s\n", "NSTELLAG", l_lag, \
            "Cross correlation lag (pixels)", >> tmphead)
        printf ("%-8s= %20.3f / %-s\n", "NSTELISH", l_shift, \
            "Initial shift (pixels)", >> tmphead)
        printf ("%-8s= %20.3f / %-s\n", "NSTELISC", l_scale, \
            "Initial scale factor for airmass ratio", >> tmphead)
        printf ("%-8s= %20.3f / %-s\n", "NSTELDSH", l_dshift, \
            "Initial shift search step (pixels)", >> tmphead)
        printf ("%-8s= %20.3f / %-s\n", "NSTELDSC", l_dscale, \
            "Initial scale factor search step", >> tmphead)
        mkheader (imgo // "[0]", tmphead, append+, verbose-)

        # Delete tmp files
        delete (tmphead, verify-, >& "dev$null")
    }


    status = 0


clean:
    if (status == 1 && madenewoutput) {
        delete ("@" // tmpoutspectra, verify-, >& "dev$null")
    }
    imdelete (tmpsci1d, verify-, >& "dev$null")
    imdelete (imagesdel, verify-, >& "dev$null")
    imdelete (caldel, verify-, >& "dev$null")
    delete (filesdel, verify-, >& "dev$null")
    imdelete (tmptwodsci, verify-, >& "dev$null")
    imdelete (tmptwodvar, verify-, >& "dev$null")

    scanfile1 = ""
    scanfile2 = ""

    printlog ("", l_logfile, l_verbose)
    if (status == 0)
        printlog ("NSTELLURIC: Exit status good", l_logfile, l_verbose)
    printlog ("----------------------------------------------------\
        ----------------------------", l_logfile, l_verbose)

end

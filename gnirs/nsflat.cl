# Copyright(c) 2001-2012 Association of Universities for Research in Astronomy, Inc.

procedure nsflat (lampson) 

# Derive GNIRS/NIRI bad pixel mask and spectroscopic flat field image
#
# Version  Sept 20, 2002 JJ v1.4 release
#          Aug 19, 2003  KL IRAF 2.12 - new/modified parameters
#                             hedit: addonly
#                             imcombine: headers,bpmasks,expmasks,outlimits
#                                 rejmask->rejmasks, plfile->nrejmasks
#                             imstat: nclip,lsigma,usigma,cache
#          Oct 29, 2003  KL moved from niri to gnirs package
#          Nov 9, 2003   AC Initial tidying, with gemextn
#          Nov 11, 2003  AC Testing and associated mods

char    lampson     {prompt = "Input flat field images"}                      # OLDP-1-primary-combine-suffix=_flat
char    darks       {"", prompt = "Input dark images (equal exposure)"}           # OLDP-1-input
char    flatfile    {"", prompt = "Output flat field image"}                      # OLDP-1-output
char    darkfile    {"", prompt = "File to save combined darks, if generated?"}   # OLDP-3
bool    fl_corner   {no, prompt = "Zero darks corners, if specified in MDF?"}     # OLDP-2
bool    fl_save_darks   {no, prompt = "Save cut darks, if generated?"}            # OLDP-3
char    flattitle   {"default", prompt = "Title for output flat image"}           # OLDP-3
char    bpmtitle    {"default", prompt = "Title for output bad pixel mask"}       # OLDP-3
char    bpmfile     {"default", prompt = "Name of output bad pixel mask PL file"} # OLDP-1

char    process     {"auto", enum = "auto|fit|median|trace", prompt = "Smoothing process"} # OLDP-2

char    statsec     {"default", prompt = "Image section for statistics (or 'MDF' to use MDF region)"} # OLDP-2
char    fitsec      {"default", prompt = "Image section for fitting (or 'MDF' to use MDF region)"} # OLDP-2
real    thr_flo     {0.35, prompt = "Lower bad-pixel threshold (fraction of peak) for flats"} # OLDP-2
real    thr_fup     {1.25, prompt = "Upper bad-pixel threshold (fraction of peak) for flats"} # OLDP-2
real    thr_dlo     {-20., prompt = "Lower bad-pixel threshold (ADU) for darks"}  # OLDP-2
real    thr_dup     {100., prompt = "Upper bad-pixel threshold (ADU) for darks"}  # OLDP-2
bool    fl_inter    {no, prompt = "Process data interactively?"}                  # OLDP-4
bool    fl_range    {no, prompt = "Set ranges interactively (if fl_inter selected)?"} # OLDP-4
bool    fl_fixbad   {yes, prompt = "Fix bad pixels in the output flat field image?"} # OLDP-2
real    fixvalue    {1.0, prompt = "Replace bad pixels in output flat with this value"} # OLDP-2
char    function    {"spline3", prompt = "Fitting function for illumination pattern"} # OLDP-2
int     order       {-1, min = -1, prompt = "Order of polynomial or number of spline pieces (-1 for default)"} # OLDP-2
char    normstat    {"midpt", enum = "mean|midpt", prompt = "Statistic to use to fine-tune normalization."} # OLDP-2
char    combtype    {"default", enum = "default|average|median", prompt = "Type of combine operation"} # OLDP-2
char    rejtype     {"ccdclip", enum = "none|minmax|ccdclip|crreject|sigclip|avsigclip|pclip", prompt = "Type of rejection"} # OLDP-2
char    masktype    {"goodvalue", enum="none|goodvalue", prompt="Mask type"}      # OLDP-3
real    maskvalue   {0., prompt="Mask value"}
char    scale       {"none", prompt = "Image scaling"}                            # OLDP-2
char    zero        {"none", prompt = "Image zeropoint offset"}                   # OLDP-2
char    weight      {"none", prompt = "Image weights"}                            # OLDP-2
real    lthreshold  {INDEF, prompt = "Lower threshold for input pixels"}          # OLDP-2
real    hthreshold  {INDEF, prompt = "Upper threshold for input pixels"}          # OLDP-2
int     nlow        {1, min = 0, prompt = "minmax: Number of low pixels to reject"} # OLDP-2
int     nhigh       {1, min = 0, prompt = "minmax: Number of high pixels to reject"} # OLDP-2
int     nkeep       {0, prompt = "Minimum to keep or maximum to reject"}          # OLDP-2
bool    mclip       {yes, prompt = "Use median in sigma clipping algorithms?"}    # OLDP-2
real    lsigma      {3., min = 0., prompt = "Lower sigma clipping factor"}        # OLDP-2
real    hsigma      {3., min = 0., prompt = "Upper sigma clipping factor"}        # OLDP-2
char    snoise      {"0.0", prompt = "ccdclip: Sensitivity noise (electrons)"}    # OLDP-2
real    sigscale    {0.1, min = 0., prompt = "Tolerance for sigma clipping scaling correction"} # OLDP-2
real    pclip       {-0.5, prompt = "pclip: Percentile clipping parameter"}       # OLDP-2
real    grow        {0.0, min = 0., prompt = "Radius (pixels) for neighbor rejection"} # OLDP-2
int     box_width   {20, prompt = "Median box width (ie spatial size)"}           # OLDP-2
int     box_length  {1, prompt = "Median box length (ie in dispersion direction)"} # OLDP-2
char    trace       {"", prompt = "Name of reference to trace"}                   # OLDP-2
char    traceproc   {"none", enum = "none|left|right|centre|smooth", prompt = "Processing for trace image"}
real    threshold   {100, prompt = "Detection threshold for aperture"}            # OLDP-3
char    aptable     {"gnirs$data/apertures.fits", prompt = "Table of aperture data"} # OLDP-3
char    database    {"", prompt = "Directory for files containing trace data"}    # OLDP-4
int     apsum       {10, prompt = "Number of dispersion lines to sum or median in fit."} # OLDP-2
int     tr_step     {10, prompt = "Tracing step"}                                 # OLDP-3
int     tr_nlost    {3, prompt = "Number of consecutive times profile can be lost"} # OLDP-3
char    tr_function {"legendre", prompt = "S-distortion fitting function"}        # OLDP-3
int     tr_order    {5, prompt = "S-distortion fitting order"}                    # OLDP-3
int     tr_naver    {1, prompt = "Trace average or median"}                       # OLDP-3
int     tr_niter    {0, prompt = "Trace rejection iterations"}                    # OLDP-3
real    tr_lowrej   {3., prompt = "S-distortion lower rejection threshold"}       # OLDP-3
real    tr_highrej  {3., prompt = "S-distortion upper rejection threshold"}       # OLDP-3
real    tr_grow     {0., prompt = "Trace rejection growing radius"}               # OLDP-3
int     ap_lower    {-30, prompt = "Default lower aperture relative to trace (for apall)"} # OLDP-3
int     ap_upper    {30, prompt = "Default upper aperture relative to trace (for apall)"} # OLDP-3

bool    fl_vardq    {yes, prompt = "Create output variance and data quality frames?"} # OLDP-2
char    logfile     {"", prompt = "Logfile name"}                                 # OLDP-1
bool    verbose     {yes, prompt = "Verbose output?"}                             # OLDP-2
int     status      {0, prompt = "Exit status (0=good)"}                          # OLDP-4
struct  *scanin1    {"", prompt = "Internal use only"}                            # OLDP-4
struct  *scanin2    {"", prompt = "Internal use only"}                            # OLDP-4

begin
    char    l_logfile = ""
    char    l_flattitle_test = ""
    char    l_bpmtitle_test = ""
    char    l_darkfile = ""
    char    l_statsec = ""
    char    l_fitsec = ""
    char    l_lampson = ""
    char    l_darks = ""
    char    l_flatfile = ""
    char    l_bpmfile = ""
    char    l_process = ""
    char    l_trace = ""
    char    l_traceproc = ""
    char    l_aptable = ""
    char    l_database = ""
    char    l_tr_function = ""
    char    l_sci_ext = ""
    char    l_var_ext = ""
    char    l_dq_ext = ""
    char    l_key_instrument = ""
    char    l_key_section = ""
    char    l_key_exptime = ""
    char    l_key_arrayid = ""
    char    l_key_filter = ""
    char    l_key_gain = ""
    char    l_key_ron = ""
    char    l_key_sat = ""
    char    l_key_nonlinear = ""
    char    l_key_cut_section = ""
    char    l_key_dispaxis = ""
    char    l_key_prism = ""
    char    l_key_decker = ""
    char    l_key_fpmask = ""
    char    l_key_mode = ""

    bool    l_fl_save_darks, l_fl_corner
    int     l_box_width, l_box_length, l_apsum
    char    tmpbpmfile, l_flattitle, l_bpmtitle, phu
    char    filter, instrument, specsec
    char    l_combtype, l_rejtype, l_masktype, l_scale, l_zero, l_weight
    char    l_normstat
    char    l_snoise, l_title, l_function
    real    ron, gain, sat, l_lthreshold, l_hthreshold, l_lsigma, l_hsigma
    real    l_grow, l_sigscale, l_pclip, dark1, dark2, dqlimit, dqbase
    real    l_maskvalue
    real    lower, upper, l_sigthresh, nonlinear, l_fixvalue, l_thresh
    real    satsci, nonlsci
    int     l_nlow, l_nhigh, l_nkeep, l_order, nfiles[2]
    bool    fl_dark, l_fl_fixbad, l_fl_vardq
    bool    l_verbose, l_mclip, l_fl_inter, l_fl_range, find
    int     l_tr_step, l_tr_nlost, l_tr_order, l_tr_naver, l_tr_niter
    int     l_ap_lower, l_ap_upper
    real    l_tr_grow
    real    l_tr_lowrej, l_tr_highrej
    int     warning, warninglimit
    bool    first, iscut, do_norm
    char    sgain, sron
    real    l_thr_dlo, l_thr_dup, l_thr_flo, l_thr_fup
    char    infiles[2], darksci, darkvar, darkdq
    char    img, imgfits, imgdq, imgsci, imgvar, imgout, dqsum, suf
    char    tmphead, aprow, tmpbpm, line
    char    tmpsci, tmpdq, combout[2], tmpout, tmpdarks, tmpmdfs, mdfdark
    char    tmptrace, tmptrace2, tmpapalllog
    char    dsub, esub, mask[3], ratio, tmpmdf, firstimg, dqmask, xsub
    char    project, projfit, section, mdfsec, statused, fitused
    real    norm, stddev, temp, mean
    int     i, nbad, len, dum, xsize, ysize, junk, maxin, dash
    int     korder, allotherorder, nver, version, dispaxis
    int     nx, ny, x1, x2, y1, y2, imdf, wx, wy, width, xover, corner
    int     slitwidth
    char    maskexpr
    real    slitsize_mx, slitsize_my, x_ccd, y_ccd, pixscale, s2n
    char    badhdr, appver, mdftemplate, tracefile, trproc
    char    tmpsubsec, tmpnorm, tmpfullnorm, imgmdf, tmpdiv, tmpdarklist
    struct  sdate, sline
    bool    cut, cutdarks, mdfdarks, dispaxis_copy
    char    arrayid
    int     apline, florder, apthresh, apupper, aplower, ipeak
    char    apsample, tracename, sdsample, bgsample, obsmode
    char    fpmask, decker, prism, search, title[2], secn, statline
    char    logfilesave, databasesave, keyfound, previmg, root
    real    shiftx, shifty, exptime
    int     specno
    char    tmpsec, tmpstr, tmpfsec, tmpssec
    int     s1,s2,f1,f2
    bool    warn
    bool    nscutfound

    char    dflt = "default"
    bool    debug = no

    cache ("keypar", "gemextn", "tinfo", "gemcombine", "specred", "gemdate", \
        "nsmdfhelper")

    # Default values

    korder = 45
    allotherorder = 20
    filter = "none"
    status = 1 # Error by default
    cutdarks = no
    title[1] = "flats"
    title[2] = "darks"
    exptime = INDEF
    dqlimit = 0.25
    obsmode = ""
    warninglimit = 10

    # Prompt for values

    junk = fscan (logfile, l_logfile)
    l_fl_vardq = fl_vardq
    l_flattitle = flattitle
    junk = fscan (l_flattitle, l_flattitle_test)
    l_bpmtitle = bpmtitle
    junk = fscan (l_bpmtitle, l_bpmtitle_test)
    junk = fscan (darkfile, l_darkfile)
    l_fl_inter = fl_inter
    l_fl_range = fl_range
    l_thr_flo = thr_flo; l_thr_fup = thr_fup
    l_thr_dlo = thr_dlo; l_thr_dup = thr_dup
    l_fl_fixbad = fl_fixbad
    junk = fscan (statsec, l_statsec)
    junk = fscan (fitsec, l_fitsec)
    junk = fscan (lampson, l_lampson)
    junk = fscan (darks, l_darks)
    l_fl_corner = fl_corner
    l_fl_save_darks = fl_save_darks
    junk = fscan (flatfile, l_flatfile)
    junk = fscan (bpmfile, l_bpmfile)
    junk = fscan (process, l_process)
    junk = fscan (trace, l_trace)
    junk = fscan (traceproc, l_traceproc)
    l_thresh = threshold
    junk = fscan (aptable, l_aptable)
    junk = fscan (database, l_database)
    junk = fscan (tr_function, l_tr_function)
    l_tr_step = tr_step
    l_tr_nlost = tr_nlost
    l_tr_order = tr_order
    l_tr_naver = tr_naver
    l_tr_niter = tr_niter
    l_tr_lowrej = tr_lowrej
    l_tr_highrej = tr_highrej
    l_tr_grow = tr_grow
    l_ap_lower = ap_lower
    l_ap_upper = ap_upper
    l_verbose = verbose
    l_box_width = box_width
    l_box_length = box_length
    l_apsum = apsum
    l_traceproc = strlwr (l_traceproc)

    badhdr = ""
    junk = fscan (nsheaders.sci_ext, l_sci_ext)
    if ("" == l_sci_ext) badhdr = badhdr + " sci_ext"
    junk = fscan (nsheaders.var_ext, l_var_ext)
    if ("" == l_var_ext) badhdr = badhdr + " var_ext"
    junk = fscan (nsheaders.dq_ext, l_dq_ext)
    if ("" == l_dq_ext) badhdr = badhdr + " dq_ext"
    junk = fscan (nsheaders.key_instrument, l_key_instrument)
    if ("" == l_key_instrument) badhdr = badhdr + " key_instrument"
    junk = fscan (nsheaders.key_section, l_key_section)
    if ("" == l_key_section) badhdr = badhdr + " key_section"
    junk = fscan (nsheaders.key_exptime, l_key_exptime)
    if ("" == l_key_exptime) badhdr = badhdr + " key_exptime"
    junk = fscan (nsheaders.key_filter, l_key_filter)
    if ("" == l_key_filter) badhdr = badhdr + " key_filter"
    junk = fscan (nsheaders.key_arrayid, l_key_arrayid)
    if ("" == l_key_arrayid) badhdr = badhdr + " key_arrayid"
    junk = fscan (nsheaders.key_gain, l_key_gain)
    if ("" == l_key_gain) badhdr = badhdr + " key_gain"
    junk = fscan (nsheaders.key_ron, l_key_ron)
    if ("" == l_key_ron) badhdr = badhdr + " key_ron"
    junk = fscan (nsheaders.key_sat, l_key_sat)
    if ("" == l_key_sat) badhdr = badhdr + " key_sat"
    junk = fscan (nsheaders.key_nonlinear, l_key_nonlinear)
    if ("" == l_key_nonlinear) badhdr = badhdr + " key_nonlinear"
    junk = fscan (nsheaders.key_cut_section, l_key_cut_section)
    if ("" == l_key_cut_section) badhdr = badhdr + " key_cut_section"
    junk = fscan (nsheaders.key_dispaxis, l_key_dispaxis)
    if ("" == l_key_dispaxis) badhdr = badhdr + " key_dispaxis"
    junk = fscan (nsheaders.key_prism, l_key_prism)
    if ("" == l_key_prism) badhdr = badhdr + " key_prism"
    junk = fscan (nsheaders.key_decker, l_key_decker)
    if ("" == l_key_decker) badhdr = badhdr + " key_decker"
    junk = fscan (nsheaders.key_fpmask, l_key_fpmask)
    if ("" == l_key_fpmask) badhdr = badhdr + " key_fpmask"
    junk = fscan (nsheaders.key_mode, l_key_mode)
    if ("" == l_key_mode) badhdr = badhdr + " key_mode"

    # Temporary files

    tmpout = mktemp ("tmpout")
    infiles[1] = mktemp ("tmponlist") 
    tmpdarklist = mktemp ("tmpdarklist") 
    tmphead = mktemp ("tmphead") 
    tmpmdfs = mktemp ("tmpmdfs") 
    tmpmdf = mktemp ("tmpmdf") 
    tmpdarks = mktemp ("tmpdarks") 
    combout[1] = mktemp ("tmpon") 
    combout[2] = mktemp ("tmpdark") 
    dsub = ""
    esub = ""
    dqsum = ""
    tmpbpmfile = ""
    tmpfullnorm = ""
    mask[1] = ""
    mask[2] = ""
    mask[3] = ""
    tmpsci = ""
    tmpdq = ""
    tmpbpm = ""
    ratio = ""
    project = ""
    projfit = ""
    tmpnorm = ""
    tmpsubsec = ""
    tmptrace = ""
    aprow = ""
    tmpapalllog = mktemp ("tmpapalllog")

    infiles[2] = tmpdarklist
    databasesave = specred.database
    logfilesave = specred.logfile

    # Start log and check for blanks / replace / validate input

    if (debug) print ("validating input")

    if (l_logfile == "") {
        junk = fscan (gnirs.logfile, l_logfile)
        if (l_logfile == "") {
            l_logfile = "gnirs.log"
            printlog ("WARNING - NSFLAT: Both nsflat.logfile and \
                gnirs.logfile are", l_logfile, verbose+) 
            printlog ("                  undefined.  Using " // l_logfile,
                l_logfile, verbose+)
        }
    }

    date | scan (sdate) 
    printlog ("----------------------------------------------------------\
        ---------------------", l_logfile, l_verbose) 
    printlog ("NSFLAT -- " // sdate, l_logfile, l_verbose) 
    printlog (" ", l_logfile, l_verbose) 

    if ("" != badhdr) {
        printlog ("ERROR - NSFLAT: Parameter(s) missing from \
            nsheaders: " // badhdr, l_logfile, verbose+) 
        goto clean
    }
 
    if (l_database == "") {
        junk = fscan (gnirs.database, l_database)
        if (l_database == "") {
            l_database = "database"
            printlog ("WARNING - NSFLAT: Both nsflat.database and \
                gnirs.database are", l_logfile, verbose+)
            printlog ("                  undefined.  Using " // l_database, \
                l_logfile, verbose+)
        }
    }

    # Remove any trailing slashes from the database name
    if (strlen (l_database) > 1 && \
        substr (l_database, strlen (l_database), strlen (l_database)) \
        == "/") {
        l_database = substr (l_database, 1, strlen (l_database) - 1)
    }

    specred.database = l_database
    specred.logfile = l_logfile


    if ("" != l_flatfile) {
        gemextn (l_flatfile, check="absent", process="none", \
            index="", extname="", extversion="", ikparams="", \
            omit="kernel,exten,section", replace="", outfile=tmpout, \
            logfile="", glogpars="", verbose=l_verbose)
        if (0 != gemextn.fail_count || 1 != gemextn.count) {
            printlog ("ERROR - NSFLAT: Output file " // l_flatfile \
                // " already exists.", l_logfile, verbose+) 
            status = 2 # Special value - l_flatfile deleted if status == 1
            goto clean
        }
        scanin1 = tmpout
        junk = fscan (scanin1, l_flatfile)
        if (debug) print ("output file: " // l_flatfile)
    }

    if ("" == l_flattitle_test || dflt == l_flattitle_test)
        l_flattitle = "FLAT FIELD IMAGE from gemini.gnirs.nsflat"
    if ("" == l_bpmtitle_test || dflt == l_bpmtitle_test)
        l_bpmtitle = "BAD PIXEL MASK from gemini.gnirs.nsflat"
    if ("" == l_key_exptime) {
        printlog ("WARNING - NSFLAT: key_exptime not set, cannot check \
            exposure times.", l_logfile, verbose+) 
    }
    if (no == l_fl_inter) {
        if ( (l_thr_flo > l_thr_fup) || (l_thr_dlo > l_thr_dup) ) {
            printlog ("ERROR - NSFLAT: Lower threshold must be lower \
                than the upper threshold.", l_logfile, verbose+) 
            goto clean
        }
    }
    if (l_fl_fixbad && !l_fl_vardq) {
        printlog ("WARNING - NSFLAT: Cannot fix bad pixels without bad \
            pixel mask, but", l_logfile, verbose+) 
        printlog ("                  fl_vardq=no.  Setting fl_fixbad=no.",
            l_logfile, verbose+) 
        l_fl_fixbad = no
    } 
    fl_dark = ! ("" == l_darks)
    if (no == fl_dark) {
        printlog ("WARNING - NSFLAT: No darks supplied.", \
            l_logfile, verbose+) 
    }

    if (fl_dark && "" != l_darkfile) {
        delete (tmpout, verify-, >& "dev$null")
        gemextn (l_darkfile, check="absent", process="none", \
            index="", extname="", extversion="", ikparams="", \
            omit="kernel,exten,section", replace="", outfile=tmpout, \
            logfile="", glogpars="", verbose=l_verbose)
        if (0 != gemextn.fail_count || 1 != gemextn.count) {
            printlog ("ERROR - NSFLAT: Output file " // l_darkfile \
                // " already exists.", l_logfile, verbose+) 
            goto clean
        }
        scanin1 = tmpout
        junk = fscan (scanin1, l_darkfile)
        if (debug) print ("output file: " // l_darkfile)
        combout[2] = l_darkfile
    }

    if (no == l_fl_vardq) {
        printlog ("WARNING - NSFLAT: Variance and quality data not \
            processed, so", l_logfile, verbose+) 
        printlog ("                  statistics will be noiser and \
            less complete.", l_logfile, verbose+) 
    }


    # Values below not processed for spaces or otherwise verified

    l_combtype = combtype; l_rejtype = rejtype; l_normstat = normstat
    l_masktype = masktype
    l_maskvalue = maskvalue
    l_scale = scale; l_zero = zero; l_weight = weight
    l_lthreshold = lthreshold; l_hthreshold = hthreshold
    l_nlow = nlow; l_nhigh = nhigh; l_nkeep = nkeep
    l_lsigma = lsigma; l_hsigma = hsigma; l_mclip = mclip
    l_snoise = snoise; l_sigscale = sigscale; l_pclip = pclip
    l_grow = grow
    l_fixvalue = fixvalue; l_order = order; l_function = function

    # Input info to logfile

    printlog ("lampson = " // l_lampson, l_logfile, l_verbose) 
    if (fl_dark) printlog ("darks = " // l_darks, l_logfile, l_verbose) 
    if (no == l_fl_inter) {
        if (fl_dark) printlog ("thr_dlo = " // l_thr_dlo, l_logfile, l_verbose)
        if (fl_dark) printlog ("thr_dup = " // l_thr_dup, l_logfile, l_verbose)
    }
    printlog ("fitting function = " // l_function, l_logfile, l_verbose) 
    printlog ("fit order = " // l_order, l_logfile, l_verbose) 
    printlog ("fitting section = " // l_fitsec, l_logfile, l_verbose) 
    printlog ("Normalization statistic = " // l_normstat, l_logfile, l_verbose)
    printlog ("Statistics region = " // l_statsec, l_logfile, l_verbose) 

    # Check files exist, are MEF, and have been through nprepare 
    # and nscut (and set order if required)

    if (debug) print ("running gemextn on file lists")

    gemextn (l_lampson, check="exists,mef", process="none", \
        index="", extname="", extversion="", ikparams="", \
        omit="kernel,exten,section", replace="", outfile=infiles[1], \
        logfile="", glogpars="", verbose=l_verbose)
    nbad = gemextn.fail_count
    if (0 == gemextn.count) {
        printlog ("ERROR - NSFLAT: No input images found for " \
            // l_lampson, l_logfile, verbose+) 
        goto clean
    }
    if (fl_dark) {
        gemextn (l_darks, check="exists,mef", process="none", \
            index="", extname="", extversion="", ikparams="", \
            omit="kernel,exten,section", replace="", outfile=infiles[2], \
            logfile="", glogpars="", verbose=l_verbose)
        nbad = nbad + gemextn.fail_count
        if (0 == gemextn.count) {
            printlog ("ERROR - NSFLAT: No input images found for " \
                // l_darks, l_logfile, verbose+) 
            goto clean
        }
    }

    # Check for previous processing

    if (debug) print ("check for nscut, nprepare")

    # this bit required because NIRI does not conform when it comes to MDFs
    scanin1 = infiles[1]
    junk = fscan (scanin1, img)
    instrument = ""
    hselect (img // "[0]", l_key_instrument, yes) | scan(instrument)

    maxin = 1
    if (fl_dark) {maxin = 2}

    for (i = 1; i <= maxin; i += 1) {
        scanin1 = infiles[i]

        first = yes

        while (fscan (scanin1, img) != EOF) {

            if (debug) print ("img: " // img)

            # Messy checking of cut and mdf status so that we can
            # handle darks from a different configuration

            keypar (img // "[0]", "NSCUT", >& "dev$null")
            nscutfound = keypar.found

            if (i == 1 && first) {

                # store first image as default status
                iscut = nscutfound
                mdftemplate = img
                if (no == iscut) {
                    printlog ("WARNING - NSFLAT: Input images not \
                        extracted with NSCUT.", l_logfile, l_verbose) 
                }

                # For NIRI, grab the specsec to apply to the darks
                if (instrument == "NIRI") {
                    specsec = ""
                    for (specno=1; specno<=3; specno=specno+1) {
                        keypar (img//"[0]", l_key_section//str(specno),
                            >& "dev$null")
                        if (keypar.found) {
                            specsec=keypar.value
                            break
                        }
                    }
                }

            } else if (i == 1 && (iscut != nscutfound)) {

                # check other images are consistent
                printlog ("ERROR - NSFLAT: Input images " // mdftemplate // \
                    " and " // img // " have inconsistent NSCUT histories.", \
                    l_logfile, l_verbose) 
                goto clean

            } else if (i == 2) {
                if (instrument == "NIRI") {
                    gemhedit (img//"[0]", l_key_section//str(specno), specsec,
                        "Region of image containing spectrum #"//str(specno),
                        delete-)
                } 
                if (first) {

                    # use first dark as default status
                    cutdarks = ! nscutfound
                    if (no == cutdarks) {
                        printlog ("WARNING - NSFLAT: Dark images already \
                            extracted with NSCUT.", l_logfile, l_verbose) 
                    }
                    if (cutdarks) {
                        gemextn (img, proc="expand", check="exists,mef", \
                            index="", extver="", extname="MDF", ikparams="", \
                            omit="", replace="", outfile="dev$null", \
                            logfile="dev$null", glogpars="", verbose-)
                        mdfdarks = !(gemextn.count == 0)
                    }

                } else if (cutdarks == nscutfound) {

                    # check other daks have same scut status
                    printlog ("ERROR - NSFLAT: Darks images have inconsistent \
                        NSCUT histories.", l_logfile, l_verbose) 
                    goto clean

                } else if (cutdarks) {

                    # check other darks have same mdf status
                    gemextn (img, proc="expand", check="exists,mef", \
                        index="", extver="", extname="MDF", ikparams="", \
                        omit="", replace="", outfile="dev$null", \
                        logfile="dev$null", glogpars="", verbose-)
                    if (mdfdarks == (gemextn.count == 0)) {
                        printlog ("ERROR - NSFLAT: Darks images have \
                            inconsistent NSPREPARE (MDF) histories.", 
                            l_logfile, l_verbose) 
                        goto clean
                    }
                }
            }

            keyfound = ""
            hselect (img // "[0]", "*PREPAR*", yes) | scan(keyfound)
            if (keyfound == "") {
                printlog ("ERROR - NSFLAT: Image " // img // " not PREPAREd.",
                    l_logfile, verbose+) 
                nbad += 1
                next
            }

            # Set default order using first file if param = -1
            if (i == 1 && first && -1 == l_order) {

                if (debug) print ("defaults from first image")

                keypar (img // "[0]", l_key_filter, >& "dev$null")
                if (keypar.found && "0" != keypar.value)
                filter = substr (keypar.value, 1, 1) 
                if (filter == "K" || filter == "Ks") {
                    l_order = korder
                } else {
                    l_order = allotherorder
                }
            }
            first = no
        }
    }

    # Process darks if required

    if (cutdarks) {

        if (debug) print ("process darks")

        if (mdfdarks) {

            # need index of MDF
            gemextn (mdftemplate, proc="expand", index="0-", \
                extname="MDF", extver="", check="exists", \
                omit="kernel", ikparams="", replace="", \
                outfile="STDOUT", logfile="dev$null", glogpars="",
                verbose-) | scan (sline)
            if (1 != gemextn.count || 0 != gemextn.fail_count) {
                printlog ("ERROR - NSFLAT: bad MDF in " // mdftemplate, \
                    l_logfile, verbose+) 
                goto clean
            }
            if (debug) print (sline)
            if (debug) print (substr (sline, stridx ("[", sline) + 1, \
                stridx ("]", sline) - 1))
            imdf = int (substr (sline, stridx ("[", sline) + 1, \
                stridx ("]", sline) - 1))

            shiftx = INDEF; shifty = INDEF
            hselect (mdftemplate // "[0]", "MDF_XSHF", yes) | scan (shiftx)
            hselect (mdftemplate // "[0]", "MDF_YSHF", yes) | scan (shifty)

            scanin1 = infiles[2]
            while (fscan (scanin1, img) != EOF) {
                mdfdark = mktemp (img // "-")
                print (mdfdark, >> tmpmdfs)
                if (debug) print (img // " --> " // mdfdark)
                copy (img // ".fits", mdfdark // ".fits")
                fxinsert (input = mdftemplate // ".fits[" // imdf // "]", \
                    output=mdfdark // ".fits[0]", groups="", verbose=debug)
                if (no == isindef (shiftx)) {
                    gemhedit (mdfdark // "[0]", "MDF_XSHF", shiftx, "", \
                        delete-)
                }
                if (no == isindef (shifty)) {
                    gemhedit (mdfdark // "[0]", "MDF_YSHF", shifty, "", \
                        delete-)
                }
            }
            infiles[2] = tmpmdfs
        }

        if (iscut) {

            if (l_fl_save_darks) {
                gemextn ("c@" // infiles[2], proc="none", check="absent", \
                    omit="", replace="", ikparams="", outfile=tmpdarks, \
                    logfile="", glogpars="", verbose=l_verbose)
                if (gemextn.fail_count > 0) {
                    printlog ("ERROR - NSFLAT: can't save darks - already \
                        exist.", l_logfile, verbose+) 
                    goto clean
                }
            } else {
                scanin1 = infiles[2]
                while (fscan (scanin1, img) != EOF) {
                    dash = strldx ("-", img)
                    if (dash > 1)
                        root = substr (img, 1, dash-1)
                    else
                        root = img
                    img = mktemp (root // "-")
                    print (img, >> tmpdarks)
                }
            }
            # nscut.section is defined by nsheaders ... don't specify it here
            nscut (inimages="@"//infiles[2], outspectra="@"//tmpdarks, \
                outprefix="", fl_corner = l_fl_corner, logfile=l_logfile, \
                verbose=l_verbose) 
            if (nscut.status != 0) {
                printlog ("ERROR - NSFLAT: problem cutting darks.", \
                    l_logfile, verbose+) 
                goto clean
            }
            infiles[2] = tmpdarks
        }
    }

    # Check that all images have the same number of SCI
    # extensions, and that corresponding var and dq data exist.

    # Compare against the first input image

    head (infiles[1], nlines=1) | scan (firstimg)
    gemextn (firstimg, check="exists", process="expand", index="", \
        extname=l_sci_ext, extversion="1-", ikparams="", omit="", \
        replace="", outfile="dev$null", logfile="", glogpars="",
        verbose-)
    if (0 != gemextn.fail_count || 0 == gemextn.count) {
        printlog ("ERROR - NSFLAT: Bad science data in " // firstimg // ".",
            l_logfile, verbose+)
        goto clean
    }
    # Assume they're 1..n - this is checked by the extver range below
    nver = gemextn.count

    for (i = 1; i <= maxin; i += 1) {
        scanin1 = infiles[i]
        while (fscan (scanin1, img) != EOF) {
            gemextn (img, check="exists", process="expand", index="", \
                extname=l_sci_ext, extversion="1-" // nver, \
                ikparams="", omit="", replace="", outfile="dev$null", \
                logfile="", glogpars="", verbose-)
            if (0 != gemextn.fail_count || nver != gemextn.count) {
                printlog ("ERROR - NSFLAT: Bad or missing science data \
                    in " // img // ".", l_logfile, verbose+)
                goto clean
            }

# Don't check for these because imcombine can generate info based on 
# the variation in the images given

#   if (l_fl_vardq) {
#       gemextn (img, check="exists", process="expand", index="", \
#           extname=l_var_ext, extversion="1-" // nver, \
#           ikparams="", omit="", replace="", outfile="dev$null", \
#           errfile="STDERR")
#       if (0 != gemextn.fail_count || nver != gemextn.count) {
#           printlog ("ERROR - NSFLAT: Bad or missing variance \
#               data in " // img // ".", l_logfile, verbose+)
#           goto clean
#       }
#       gemextn (img, check="exists", process="expand", index="", \
#           extname=l_dq_ext, extversion="1-" // nver, \
#           ikparams="", omit="", replace="", outfile="dev$null", \
#           errfile="STDERR")
#       if (0 != gemextn.fail_count || nver != gemextn.count) {
#           printlog ("ERROR - NSFLAT: Bad or missing quality \
#               data in " // img // ".", l_logfile, verbose+)
#           goto clean
#       }
#   }

        }
    }

    if (nbad > 0) {
        printlog ("ERROR - NSFLAT: " // nbad // " image(s) either do \
            not exist, are not MEF files, or", l_logfile, verbose+) 
        printlog ("                have not been run through *PREPARE \
            and NSCUT.", l_logfile, verbose+) 
        goto clean
    }

    # More checking and main process

    if (debug) print ("further checking and processing: " // maxin)

    scanin1 = infiles[1]
    junk = fscan (scanin1, img)
    scanin1 = ""
    hselect (img//"[0]", l_key_instrument, yes) | scan (instrument)
    if (dflt == l_statsec) {
        if (instrument == "NIRI")
            l_statsec = "[*,*]"
        else
            l_statsec = "MDF"
    }
    if (dflt == l_fitsec) {
        if (instrument == "NIRI")
            l_fitsec = "[*,*]"
        else
            l_fitsec = "MDF"
    }

    first = yes
    warning = 0

    # Data then darks
    for (i = 1; i <= maxin; i += 1) {

        scanin1 = infiles[i]
        nfiles[i] = 0

        while (fscan (scanin1, img) != EOF) {
            nfiles[i] = nfiles[i] + 1
            imgfits = img // ".fits"
            if (debug) print ("img: " // img // " (" // i // ")")

            # Construct default names from first file, if necessary,
            # and copy PHU to output

            if (first) {

                # Simpler than original, which wasn't consistent (imho)
                if ("" == l_flatfile) l_flatfile = img // "_flat"
                if (l_fl_vardq && ("" == l_bpmfile || dflt == l_bpmfile)) {
                    l_bpmfile = l_flatfile + "_bpm.pl"
                }
                printlog ("flatfile = " // l_flatfile, l_logfile, l_verbose) 
                printlog ("bpmfile = " // l_bpmfile, l_logfile, l_verbose) 
                if (imaccess (l_flatfile) ) {
                    printlog ("ERROR - NSFLAT: Output file " // l_flatfile // \
                        " already exists.", l_logfile, verbose+) 
                    goto clean
                }
                if (imaccess (l_bpmfile) ) {
                    printlog ("ERROR - NSFLAT: Output BPM " // l_bpmfile \
                        // " already exists.", l_logfile, verbose+) 
                    goto clean
                }

                # Create output file
                imcopy (imgfits // "[0]", l_flatfile, verbose-) 

                dispaxis = 1
                phu = img // "[0]"
                keypar (phu, l_key_dispaxis, silent+)
                if (keypar.found) {
                    dispaxis = int (keypar.value)
                } else {
                    printlog ("WARNING - NSFLAT: Can't determine dispersion \
                        axis, assuming " // dispaxis, l_logfile, l_verbose)
                }
                if (debug) print ("dispaxis: " // dispaxis)

                # Get MDF size

                if ("MDF" == l_fitsec || "MDF" == l_statsec) {

                    if (iscut) {

                        mdfsec = ""

                    } else {

                        tinfo (img // ".fits[MDF]", ttout-)
                        if (tinfo.nrows > 1) {
                            printlog ("ERROR - NSFLAT: Multiple MDF \
                                with uncut data; can't use MDF \
                                fit/statsec.", l_logfile, l_verbose)
                            goto clean
                        }

                        phu = img // "[0]"
                        imgsci = img // "[" // l_sci_ext // ",1]"

                        keypar (phu, "PIXSCALE", silent+)
                        if (no == keypar.found) {
                            printlog ("ERROR - NSFLAT: No PIXSCALE in " \
                                // phu // ".", l_logfile, verbose+)
                            goto clean
                        }
                        pixscale = real (keypar.value)
                        keypar (imgsci, "i_naxis1", silent+)
                        if (no == keypar.found) {
                            printlog ("ERROR - NSFLAT: No i_naxis1 in " \
                                // imgsci // ".", l_logfile, l_verbose)
                            goto clean
                        }
                        nx = int (keypar.value)
                        keypar (imgsci, "i_naxis2", silent+)
                        if (no == keypar.found) {
                            printlog ("ERROR - NSFLAT: No i_naxis2 in " \
                                // imgsci // ".", l_logfile, l_verbose)
                            goto clean
                        }
                        ny = int (keypar.value)

                        gemextn (img, check="exists,mef", proc="append", 
                            index="", extname="MDF", extver="", ikparams="",
                            omit="", replace="", outfile="dev$null", 
                            logfile="", glogpars="", verbose=l_verbose)
                        if (1 != gemextn.count || 0 != gemextn.fail_count) {
                            printlog ("ERROR - NSFLAT: No MDF extension.", 
                                l_logfile, verbose+)
                            goto clean
                        }

                        tprint (img // ".fits[MDF]", prparam-, prdata+, \
                            showrow-, showhdr-, showunits-, \
                            col="x_ccd,y_ccd,slitsize_mx,slitsize_my", \
                            rows=1, pwidth=160) \
                            | scan (x_ccd,y_ccd,slitsize_mx,slitsize_my)

                        slitsize_mx = slitsize_mx / pixscale
                        slitsize_my = slitsize_my / pixscale


                        # Copy subsection directly
                        if (1 == dispaxis) {
                            x1 = 1
                            x2 = nx
                            y1 = y_ccd - 0.5 * slitsize_my
                            y2 = y_ccd + 0.5 * slitsize_my
                        } else {
                            y1 = 1
                            y2 = ny
                            x1 = x_ccd - 0.5 * slitsize_mx
                            x2 = x_ccd + 0.5 * slitsize_mx
                        }
                        mdfsec = "[" // x1 // ":" // x2 // "," \
                            // y1 // ":" // y2 // "]"
                        if (debug) print ("mdf rectangle: " // mdfsec)
                    }
                }

                if ("MDF" == l_fitsec)
                    fitused = mdfsec
                else
                    fitused = l_fitsec

                if ("MDF" == l_statsec)
                    statused = mdfsec
                else
                    statused = l_statsec
            }

            # Check exposure time from frame to frame (inc darks)
            keypar (imgfits // "[0]", l_key_exptime, >& "dev$null")
            if (no == keypar.found) {
                printlog ("ERROR - NSFLAT: missing exposure time in " // img,
                    l_logfile, verbose+) 
                goto clean
            } else if (first) {
                exptime = real (keypar.value)
            } else if (warning < warninglimit) {
                if (abs (real (keypar.value) - exptime) > 0.1 && \
                    no == ("exposure" == strlwr (l_scale))) {
                    
                    printlog ("WARNING - NSFLAT: exposure time changed",
                        l_logfile, verbose+)
                    printf ("                  %-25s: %6.2f\n", \
                        previmg, exptime) | scan (sline)
                    printlog (sline, l_logfile, verbose+)
                    printf ("                  %-25s: %6.2f\n", \
                        img, real (keypar.value)) | scan (sline)
                    printlog (sline, l_logfile, verbose+)
                    warning = warning + 1
                }
                exptime = real (keypar.value)
            }

            # Check gain, readnoise, saturation and non-linear values

            keypar (imgfits // "[0]", l_key_gain, >& "dev$null") 
            if (no == keypar.found) {
                printlog ("ERROR - NSFLAT: Could not get gain from \
                    header of image " // img, l_logfile, verbose+) 
                goto clean
            } else if (first) {
                gain = real (keypar.value) 
                sgain = keypar.value
            } else if (warning < warninglimit) {
                if (abs (real (keypar.value) - gain) > 0.5) {
                    printlog ("WARNING - NSFLAT: gain values changed", \
                        l_logfile, verbose+) 
                    printf ("                  %-25s: %6.2f\n", \
                        previmg, gain) | scan (sline)
                    printlog (sline, l_logfile, verbose+)
                    printf ("                  %-25s: %6.2f\n", \
                        img, real (keypar.value)) | scan (sline)
                    printlog (sline, l_logfile, verbose+)
                    warning = warning + 1
                }
                gain = real (keypar.value) 
                sgain = keypar.value
            }

            keypar (imgfits // "[0]", l_key_ron, >& "dev$null") 
            if (no == keypar.found) {
                printlog ("ERROR - NSFLAT: Could not get read noise \
                    from header of image " // img, l_logfile, verbose+) 
                goto clean
            } else if (first) {
                ron = real (keypar.value) 
                sron = keypar.value
            } else if (warning < warninglimit) {
                if (abs (real (keypar.value) - ron) > 0.5) {
                    printlog ("WARNING - NSFLAT: read noise values \
                        changed", l_logfile, verbose+) 
                    printf ("                  %-25s: %6.2f\n", \
                        previmg, ron) | scan (sline)
                    printlog (sline, l_logfile, verbose+)
                    printf ("                  %-25s: %6.2f\n", \
                        img, real (keypar.value)) | scan (sline)
                    printlog (sline, l_logfile, verbose+)
                    warning = warning + 1
                }
                ron = real (keypar.value) 
                sron = keypar.value
            }

            keypar (imgfits // "[0]", l_key_sat, >& "dev$null") 
            if (no == keypar.found) {
                printlog ("ERROR - NSFLAT: Could not get saturation \
                    level from header of image " // img, l_logfile, verbose+) 
                goto clean
            } else if (first) {
                sat = real (keypar.value) 
            } else if (warning < warninglimit) {
                if (abs (real (keypar.value) - sat) > 1.0) {
                    printlog ("WARNING - NSFLAT: saturation values changed",
                        l_logfile, verbose+) 
                    printf ("                  %-25s: %6.0f\n", \
                        previmg, sat) | scan (sline)
                    printlog (sline, l_logfile, verbose+)
                    printf ("                  %-25s: %6.0f\n", \
                        img, real (keypar.value)) | scan (sline)
                    printlog (sline, l_logfile, verbose+)
                    warning = warning + 1
                }
                sat = real (keypar.value) 
            }

            keypar (imgfits // "[0]", l_key_nonlinear, >& "dev$null") 
            if (no == keypar.found) {
                printlog ("ERROR - NSFLAT: Could not get non-linear \
                    level from header of image " // img, l_logfile, verbose+) 
                goto clean
            } else if (first) {
                nonlinear = real (keypar.value) 
            } else if (warning < warninglimit) {
                if (abs (real (keypar.value) - nonlinear) > 1.0) {
                    printlog ("WARNING - NSFLAT: nonlinear values changed",
                        l_logfile, verbose+) 
                    printf ("                  %-25s: %6.0f\n", \
                        previmg, nonlinear) | scan (sline)
                    printlog (sline, l_logfile, verbose+)
                    printf ("                  %-25s: %6.0f\n", \
                        img, real (keypar.value)) | scan (sline)
                    printlog (sline, l_logfile, verbose+)
                    warning = warning + 1
                }
                nonlinear = real (keypar.value) 
            }

            if (warning == warninglimit) {
                printlog ("WARNING - NSFLAT: further header value \
                    warnings supressed.", l_logfile, verbose+) 
                warning = warning + 1
            }

            first = no
            previmg = img

            if (i == 1) {
                nonlsci = nonlinear
                satsci = sat
            }

        } # For all files

        # Combine images (with special cases for small numbers)
        if (debug) print ("combining")
        if (1 == nfiles[i]) {
            printlog ("WARNING - NSFLAT: only one input image (" \
                // title[i] // ").", l_logfile, verbose+) 
            imcopy (imgfits // "[0]", combout[i] // "[0,overwrite]", \
                verbose-, >& "dev$null")

            # Need to get index that corresponds to mdf, because
            # fxinsert doesn't handle extname(!)
            delete (tmpmdf, verify-, >& "dev$null")
            gemextn (imgfits, check="exists,mef", proc="expand", \
                index="1-", extname="MDF", extver="", ikparams="", \
                omit="kernel", replace="", outfile=tmpmdf, \
                logfile="dev$null", glogpars="", verbose-)
            if (debug) type (tmpmdf)
            if (1 == gemextn.count) {
                scanin2 = tmpmdf
                junk = fscan (scanin2, imgmdf)
                if (debug) print (junk)
                if (debug) print ("fxinsert " // imgmdf)
                fxinsert (imgmdf, combout[i] // ".fits[1]", group="", verbose-)
            }

            for (version = 1; version <= nver; version = version+1) {
                imgsci = imgfits // "[" // l_sci_ext // "," // version // "]"
                imgvar = imgfits // "[" // l_var_ext // "," // version // "]"
                imgdq = imgfits // "[" // l_dq_ext // "," // version // "]"
                imcopy (imgsci, combout[i] // ".fits[" // l_sci_ext \
                    // "," // version // ",append]", verbose-) 
                if (imaccess (imgdq) && imaccess (imgvar)) {
                    imcopy (imgvar, combout[i] // ".fits[" // l_var_ext \
                        // "," // version // ",append]", verbose-) 
                    imcopy (imgdq, combout[i] // ".fits[" // l_dq_ext \
                        // "," // version // ",append]", verbose-) 
                } else if (l_fl_vardq) {
                    l_fl_vardq = no
                    l_fl_fixbad = no
                    fl_dark = no
                    printlog ("WARNING - NSFLAT: only one input image \
                        without VAR and DQ planes.", l_logfile, verbose+) 
                    printlog ("                  Proceeding with \
                        fl_vardq=no and fl_fixbad=no", l_logfile, verbose+)
                }
            }

        } else {
            if (2 == nfiles[i]) {
                printlog ("WARNING - NSFLAT: only two images to \
                    combine, turning off rejection.", l_logfile, verbose+) 
                l_rejtype = "none"
            } else if (4 >= nfiles[i]) {
                printlog ("WARNING - NSFLAT: four or less images to combine.",
                    l_logfile, verbose+) 
            }

            if (l_combtype == dflt) l_combtype = "average"

            if (debug) print ("gemcombine")

            # Note that we have fl_vardq+ here, whatever is chosen
            # in the parameter list.  This is necessary for the task
            # logic, which uses the bad pixel mask from combination
            # whatever the value of the task parameter.

            # This is a carry-over from the original code, which 
            # did not handle multiple extensions, and which called
            # imcombine directly.  The call below reproduces that
            # behaviour.

            printlog ("NSFLAT: Combining data.", l_logfile, l_verbose)

            gemcombine ("@" // infiles[i], combout[i], \
                logfile = l_logfile, combine = l_combtype, \
                reject = l_rejtype, offsets = "none", \
                masktype=l_masktype, maskvalue=l_maskvalue, \
                scale = l_scale, zero = l_zero, weight = l_weight, \
                statsec = statused, expname = l_key_exptime, \
                lthreshold = l_lthreshold, hthreshold = l_hthreshold, \
                nlow = l_nlow, nhigh = l_nhigh, nkeep = l_nkeep, \
                mclip = l_mclip, lsigma = l_lsigma, hsigma = l_hsigma, \
                key_ron = l_key_ron, key_gain = l_key_gain, ron = ron, \
                gain = gain, snoise = l_snoise, sigscale = l_sigscale, \
                pclip = l_pclip, grow = l_grow, bpmfile = "", \
                nrejfile = "", sci_ext = l_sci_ext, var_ext = l_var_ext, \
                dq_ext = l_dq_ext, fl_vardq+, verbose-)

            if (debug) print ("gemcombine done")

            if (0 != gemcombine.status) {
                printlog ("ERROR - NSFLAT: Error in gemcombine.", \
                    l_logfile, verbose+) 
                goto clean
            }
        }

        if (debug) print ("header mods")

    } # And do darks after onlamps

    # Decide processing type

    #Modified to remove process=trace option (Jul2006, BR)
    if ("auto" == l_process) {
        hselect (mdftemplate // "[0]", l_key_mode, yes) | scan (obsmode)
        obsmode = strupr (obsmode)
        if ("IFU" == obsmode)
            l_process = "fit"
        else if ("XD" == obsmode)
            l_process = "fit"
        else
            l_process = "fit"
        printlog ("NSFLAT: process type: " // l_process, l_logfile, l_verbose)
    }
    if ("trace" == l_process) {
        printlog ("WARNING - NSFLAT: Use of process=trace is NOT recommended. \
            We recommend process=fit.", l_logfile, verbose+)
    }

    if ("fit" == l_process) {
        # a statused along the dispersion larger than fitused
        # does not make sense for 'fit'.  Check, correct and issue warning
        
        # KL 2007.03.19 - What a convoluted piece of code!  From what I
        #       make out, 'statused' and 'fitused' are set earlier.  If they
        #       are 'MDF', then they are already set to 'mdfsec', and the
        #       code below will not change that.
        #
        #       Changes made to help clarity (sigh!) but mostly to make sure
        #       it won't crash under PyRAF.
        #
        #       TODO:  really, this need to be rewritten.  actually, the whole
        #           script need to be rewritten.

        warn = no
        if ((l_fitsec != "MDF") && (l_statsec == "MDF")) {
            warn = yes
            statused = fitused
        } else if ((l_fitsec == "MDF") && (l_statsec == "MDF")) {
            # KL 2007.03.19 - fitused and statused already set to 'mdfsec'.  
            #       If data already cut, 'mdfsec' is an empty string and the
            #       whole image is used for statistics and fitting.
            #       Just continue.
            ;
        } else {
            if (dispaxis == 1)
                tmpfsec = substr (fitused, 2, stridx(",",fitused)-1)
            else
                tmpfsec = substr (fitused, stridx(",",fitused)+1, 
                    strlen(fitused)-1)

            if (tmpfsec == "*")
                statused = l_statsec
            else {
                if (dispaxis == 1)
                    tmpssec = substr (l_statsec, 2, stridx(",",l_statsec)-1)
                else
                    tmpssec = substr (l_statsec, stridx(",",l_statsec)+1,
                        strlen(l_statsec)-1)

                if (tmpssec == "*") {
                    warn = yes
                    tmpssec = tmpfsec
                } else {
                    junk = fscanf (tmpssec, "%d:%d", s1, s2)
                    junk = fscanf (tmpfsec, "%d:%d", f1, f2)
                    if (s1 < f1) {
                        warn = yes
                        s1 = f1
                    }
                    if (s2 > f2) {
                        warn = yes
                        s2 = f2
                    }
                    if (warn) { tmpssec = s1//":"//s2 }
                }
                if (warn) {
                    if (dispaxis == 1) {
                        printf ("[%s%s\n", tmpssec, substr (l_statsec,
                            stridx(",",l_statsec), strlen(l_statsec))) | \
                            scan(statused)
                    } else {
                        printf ("%s%s]\n", substr (l_statsec, 1, 
                            stridx(",",l_statsec)), tmpssec) | scan(statused)
                    }
                }
            }
        }
        if (warn) {
            print ("WARNING - NSFLAT: statsec made to fit within fitsec \
                (along the dispersion axis)")
            print ("WARNING - NSFLAT: New statsec = "//statused)
        }
    }

    if ("trace" == l_process && l_fl_fixbad) {
        printlog ("WARNING - NSFLAT: Cannot have fl_fixbad with \
            trace processing.", l_logfile, l_verbose)
        printlog ("                  Setting fl_fixbad='no'", 
            l_logfile, l_verbose)
        l_fl_fixbad = no
    }
    if ("trace" == l_process && "" == l_trace) {
        printlog ("WARNING - NSFLAT: Tracing flatfield data directly.", \
            l_logfile, l_verbose)
    }
    
    # Print some basic help for the interactive mode
    if (l_fl_inter) {
        if ("fit" == l_process) {
            printf ("\nCommon Interactive Commands:\n")
            printf ("   'd'  -  Delete point nearest to cursor\n")
            printf ("   'f'  -  Fit data and redraw\n")
            printf ("   'q'  -  Exit the interactive curve fitting\n")
            printf ("\n")
            printf ("   '?'  -  Print complete list of options\n\n")
        } else if ("trace" == l_process) {
            printf ("\nCommon Interactive Commands:\n")
            printf ("  Defining the Aperture:\n")
            printf ("     'd'  -  Delete aperture\n")
            printf ("     'l'  -  Define lower bound for the aperture\n")
            printf ("     'm'  -  Mark new aperture\n")
            printf ("     'u'  -  Define upper bound for the aperture\n")
            printf ("\n")
            printf ("     'q'  -  Done\n")
            printf ("\n")
            printf ("  Fitting:\n")
            printf ("     'd'  -  Delete point nearest to cursor\n")
            printf ("     'f'  -  Fit data and redraw\n")
            printf ("     'u'  -  Un-delete point nearest to cursor\n")
            printf ("\n")
            printf ("     'q'  -  Done\n\n")
        }
    }

    # Read header info for tracing
    if ("trace" == l_process) {

        if ("" == l_trace)
            tracefile = combout[1]
        else
            tracefile = l_trace

        phu = tracefile // "[0]"

        keypar (phu, l_key_prism, silent+)
        if (no == keypar.found) {
            printlog ("ERROR - NSFLAT: No " // l_key_prism // " in " \
                // phu // ".", l_logfile, l_verbose)
            goto clean
        }
        prism = keypar.value
        if (debug) print ("prism: " // prism)

        keypar (phu, l_key_decker, silent+)
        if (no == keypar.found) {
            printlog ("ERROR - NSFLAT: No " // l_key_decker // " in " \
                // phu // ".", l_logfile, l_verbose)
            goto clean
        }
        decker = keypar.value
        if (debug) print ("decker: " // decker)

        keypar (phu, l_key_fpmask, silent+)
        if (no == keypar.found) {
            printlog ("ERROR - NSFLAT: No " // l_key_fpmask // " in " \
                // phu // ".", l_logfile, l_verbose)
            goto clean
        }
        fpmask = keypar.value
        if (debug) print ("fpmask: " // fpmask)
        
        keypar (phu, l_key_arrayid, silent+)
        if (no == keypar.found) {
            printlog ("ERROR - NSFLAT: No " // l_key_arrayid // " in " \
                // phu // ".", l_logfile, l_verbose)
            goto clean
        }
        arrayid = keypar.value
        if (debug) print ("arrayid: " // arrayid)

    }

    # Now further process the generated files - subtract and
    # normalize flats, combine bad pixel masks.

    # This is done extension by extension

    if (debug) print ("process by version: " // nver)
    printlog ("", l_logfile, l_verbose)
    printlog (" Extn            Stats.        Lamps         Darks         \
        Flat", l_logfile, l_verbose)
    printlog (" no.             section     Mean   S/N    Mean   S/N    \
        Mean   S/N", l_logfile, l_verbose)

    for (version = 1; version <= nver; version = version+1) {

        if (debug) print ("combining flats etc: " // version)
        statline = ""
        do_norm = yes

        if (version > 1)
            imdelete (dsub // "," // esub // "," // dqsum // "," \
                // tmpbpmfile // "," // mask[1] // "," // mask[2] \
                // "," // mask[3] // "," // tmpsci // "," // tmpdq \
                // "," // ratio // "," // project // "," // projfit \
                // "," // tmpnorm // "," // tmpsubsec // "," // tmpbpm, \
                verify-, >& "dev$null") 

        dsub = mktemp ("tmpdsub") 
        esub = mktemp ("tmpesub") 
        dqsum = mktemp ("tmpdqsum")
        tmpbpmfile = mktemp ("tmpbpmfile") 
        mask[1] = mktemp ("tmpmask1") 
        mask[2] = mktemp ("tmpmask2") 
        mask[3] = mktemp ("tmpmask3") 
        tmpsci = mktemp ("tmpsci") 
        tmpdq = mktemp ("tmpdq") 
        tmpbpm = mktemp ("tmpbpm") 
        ratio = mktemp ("tmpratio")
        project = mktemp ("tmpproject") 
        projfit = mktemp ("tmpprojfit") 
        tmpnorm = mktemp ("tmpnorm")
        tmpsubsec = mktemp ("tmpsubsec")

        if (nver > 1) appver = "," // version
        else appver = ""
        imgsci = combout[1] // "[" // l_sci_ext // appver // "]"
        if (l_fl_vardq) {
            imgvar = combout[1] // "[" // l_var_ext // appver // "]"
            imgdq = combout[1] // "[" // l_dq_ext // appver // "]"
            imstat (imgdq, fields="mean", lower=INDEF, \
                upper=INDEF, nclip=0, lsigma=INDEF, usigma=INDEF, \
                binwidth=0.1, format-, cache-) | scan (dqbase)
        }

        if (debug) {
            print ("imgsci stats")
            imstat (imgsci)
        }

        if (fl_dark) {

            # Subtracting darks from lampson

            if (debug) print ("subtracting darks from lampson")

            if (fl_dark) {
                darksci = combout[2] // "[" // l_sci_ext // appver // "]"
                darkvar = combout[2] // "[" // l_var_ext // appver // "]"
                darkdq = combout[2] // "[" // l_dq_ext // appver // "]"
                if (debug) {
                    print ("darkdq stats:")
                    imstat (darkdq)
                }
            }

            # Limits either from plot or params
            if (l_fl_inter && l_fl_range) {
                imhistogram (imgsci, z1 = -100, z2 = 500, \
                    binwidth = INDEF, nbins = 512, autoscale+, \
                    top_closed-, hist_type = "normal", listout-, \
                    plot_type = "line", logy+, device = "stdgraph") 
                print ("Dark histogram") 
                print ("Hit any key to mark lower limit for good pixels") 
                dum = fscan (gcur, lower) 
                print ("Hit any key to mark upper limit for good pixels") 
                dum = fscan (gcur, upper) 
                print (" ") 
                if (lower > upper) {
                    temp = upper
                    upper = lower
                    lower = temp
                }
            } else {
                lower = l_thr_dlo
                upper = l_thr_dup
            }

            if (debug) print ("now do stats")

            # Set DQ using old DQ and thresholds with new data
            if (debug) {
                print ("mask from ranges : " // lower // ", " // upper)
                print (mask[2] // ", " // darksci // ", " // l_fl_vardq)
            }
            if (l_fl_vardq) {
                imexpr ("((a<" // lower // ") || (a>" // upper \
                    // ") || (b>0)) ? 1 : 0", mask[2], darksci, darkdq, \
                    outtype = "ushort", verbose-) 
                if (debug) {
                    print ("clipping dark, mask stats:")
                    imstat (mask[2])
                }
            } else {
                imexpr ("((a<" // lower // ") || (a>" // upper // ")) \
                    ? 1 : 0", mask[2], darksci, outtype = "ushort", verbose-) 
            }
            if (debug) print ("dq test for " // mask[2] // fitused)
            if (debug) 
                imstat (mask[2] // fitused, fields="mean", lower=INDEF, \
                    upper=INDEF, nclip=0, lsigma=INDEF, usigma=INDEF, \
                    binwidth=0.1, format+, cache-)
            # Upper limit excludes mask value (8) for XD data
            imstat (mask[2] // fitused, fields="mean", lower=INDEF, \
                upper=INDEF, nclip=0, lsigma=INDEF, usigma=INDEF, \
                binwidth=0.1, format-, cache-) | scan (temp)
            if (temp > dqlimit) {
                printlog ("WARNING - NSFLAT: DQ for darks is poor \
                    (thresh: " // lower // " to " // upper // ")", \
                    l_logfile, verbose+)
            }

            # Save for logging
            dark1 = lower; dark2 = upper


            if (debug) print ("now do subtraction")

            # Do the subtraction
            imarith (imgsci, "-", darksci, dsub, pixtype="real", verbose-)

            # Update data quality
            if (l_fl_vardq) {
                if (debug) print (imgsci)
                imexpr ("(a>" // satsci // ") ? 4 : ((a>" // nonlsci // \
                    ") ? 2 : 0)", mask[3], imgsci, outtype = "ushort", \
                    verbose-)

                imarith (imgvar, "+", darkvar, esub, pixtype = "real", 
                    verbose-) 

                if (debug) {
                    print ("imgdqs stats:")
                    imstat (imgdq)
                }

                addmasks (imgdq // "," // darkdq // "," // mask[3], \
                    dqsum // ".fits", "im1 || im2 || im3", flags = " ")
            }

        } else {

            # No darks to subtract
            
            imcopy (imgsci, dsub, verbose-) 
            if (l_fl_vardq) {
                imexpr ("(a>" // satsci // ") ? 4 : ((a>" // nonlsci // \
                    ") ? 2 : 0)", mask[3], imgsci, outtype = "ushort", \
                    verbose-)
                addmasks (imgdq // "," // mask[3], dqsum // ".fits", \
                    "im1 || im2", flags = " ")
                imcopy (imgvar, esub, verbose-)
            }
        }

        # Check that less than 10% of dq (approx, depending on values) is bad

        if (l_fl_vardq) {
            imstat (mask[3], fields="mean", lower=INDEF, upper=INDEF, \
                nclip=0, lsigma=INDEF, usigma=INDEF, binwidth=0.1, \
                format-, cache-) | scan (temp)
            if (temp > dqlimit) {
                printlog ("WARNING - NSFLAT: DQ for lampson is poor.", \
                    l_logfile, verbose+)
                printlog ("                  (saturation=" // satsci \
                    // "; nonlinear=" // nonlsci // ")", l_logfile, verbose+)
            }
        }

        # OK, this is where we do the real work.
        # We need to generate a smoothed version of the data for 
        # normalization - then we divide the flat obs by that.
        # The smoothed data can be either a fit along the dispersion
        # direction or a median filter.

        if (debug) print ("interp/project")

        if (l_fl_vardq) {
            if (debug) print ("fixpix: " // dsub // ", " // dqsum)
            if (debug) imhead (dsub)
            if (debug) imhead (dqsum)

            # This call to fixpix can sometimes cause the output image (dsub)
            # to have different values, depending on whether IRAF or PyRAF 
            # was used. The dsub and dqsum created by IRAF and PyRAF are 
            # identical before running fixpix. The dsub file output by fixpix 
            # using IRAF has values != 0 in the "bad pixel" areas and should 
            # have values = 0. This is where the problem occurs. Changing the 
            # mask from .fits to .pl using the imarith call below solves the 
            # problem. EH

            imarith (dqsum//".fits", "*", 1, dqsum//".pl", pixtype="ushort", \
                verbose-)

            proto.fixpix (dsub, dqsum // ".pl", linterp="INDEF", \
                cinterp="INDEF", verbose-, pixels-, >& "dev$null")

            # fixpix fixes pixels in the masked out corner area giving values
            # to non-values.  (i.e. the corners are zeros and should remain 
            # zeros). Here we re-mask the corners.
            #
            # Warning: before this re-mask was implemented the fit1d done later
            # was actually done on an average that *included* the garbage 
            # corners!  This new implementation should be scientifically 
            # correct.  I doubt the previous one was correct at all.  KL

            # Remember than NIRI doesn't have MDFs or corners

            if (instrument == "GNIRS") {
                nsmdfhelper (combout[1] // ".fits[MDF]", version, \
                    combout[1] // "[" // l_sci_ext // ",1]", \
                    area="spectrum", logfile=l_logfile, logname="NSFLAT", \
                    verbose=debug)

                if (no == (0 == nsmdfhelper.status)) {
                    printlog ("ERROR - NSFLAT: Could not read MDF.", \
                        l_logfile, verbose+)
                    goto clean
                }

                corner = 0
                xover = 0
                slitwidth = 0
                corner = nsmdfhelper.corner
                xover = nsmdfhelper.ixovershoot
                slitwidth = nsmdfhelper.slitwidth
                x1 = nsmdfhelper.ixlo
                x2 = nsmdfhelper.ixhi
                nx = x2 - x1
                keypar (dsub, "i_naxis2", silent+)

                if (no == keypar.found) {
                    printlog ("ERROR - NSFLAT: No i_naxis2 in " // dsub // \
                        ".", l_logfile, verbose+)
                    goto clean
                } else {
                    ny = int (keypar.value)
                }

                if (l_fl_corner && 0 != corner && (nx + xover - corner) > 0) {
                    xsub = mktemp("tmpxsub")
                    maskexpr = "(I < " // corner // "*J/" // ny // " || " // \
                       "I - " // xover // " > " // slitwidth // "+" // \
                       corner // "*J/" // ny // ")"
                    if (debug) print (maskexpr)
                    imexpr (maskexpr // "? 0 : a", xsub, dsub, verbose=debug)
                    imdelete (dsub, verify-, >& "dev$null")
                    dsub = xsub
                } else {
                    if (debug) print ("no corner: " // l_fl_corner // ", " \
                        // corner // ", " // nx)
                }
            }
            imdelete (dqsum // ".pl", verify-, >& "dev$null")
        }

        # Print statistics on combined data
        for (i = 1; i <= maxin; i += 1) {
            secn = "[" // l_sci_ext // "," // version // "]"
            if (l_fl_vardq) {
                mimstat (combout[i] // secn // statused, imasks = combout[i] \
                    // "[" // l_dq_ext // "," // version // "]" \
                    // statused, omasks = "", fields = "midpt,stddev", 
                    lower = INDEF, upper = satsci, nclip = 0, lsigma = INDEF, \
                    usigma = INDEF, binwidth = 0.1, format-, cache-) \
                    | scan (temp, stddev) 
                mimstat (combout[i] // secn // statused, imasks = combout[i] \
                    // "[" // l_dq_ext // "," // version // "]" \
                    // statused, omasks = "", fields = "mean", 
                    lower = (temp-4*stddev), upper = (temp+4*stddev), \
                    nclip = 0, lsigma = INDEF, usigma = INDEF, 
                    binwidth = 0.1, format-, cache-) | scan (temp) 
            } else {
                imstat (combout[i] // secn // statused, \
                    fields = "midpt,stddev", lower = INDEF, upper = satsci, \
                    nclip = 0, lsigma = INDEF, usigma = INDEF, \
                    binwidth = 0.1, format-, cache-) | scan (temp, stddev) 
                imstat (combout[i] // secn // statused, fields = "mean", \
                    lower = (temp-4*stddev), upper = (temp+4*stddev), \
                    nclip = 0, lsigma = INDEF, usigma = INDEF, 
                    binwidth = 0.1, format-, cache-) | scan (temp) 
            }
            s2n = ((gain*temp)/(gain*temp+ron*ron)**0.5)*(nfiles[i]**0.5)

            if (1 == i) {
                printf ("%4d %20s  \n", version, statused) | scan (sline)
                statline = statline // substr (sline, 1, strlen (sline) - 1)
            }
            printf ("%7.1f %5.1f  \n", temp, s2n)  | scan (sline)
            statline = statline // substr (sline, 1, strlen (sline) - 1)
        }

        # Skip o/p if no darks
        if (1 == maxin) statline = statline // "             "

        # Get subsection required

        if ("" != fitused) {
            if (debug) print ("fit subsection: " // fitused)
            tmpsubsec = mktemp("tmpsubsec")
            imcopy (dsub // fitused, tmpsubsec, verbose-)
        } else {
            if (debug) print ("fit full image")
            tmpsubsec = dsub
        }

        if ("fit" == l_process) {

            # Normalization from fitting

            # Project across dispersion direction
            if (debug) print ("project")
            improject (tmpsubsec, project, 3-dispaxis, average+, \
                highcut = satsci, lowcut = 0, verbose-)

            # Fit a curve
            if (debug) print ("fit")
            # Note: axis is 1 here no matter the instrument because 'project'
            # is 1-D.
            fit1d (project, projfit, "fit", axis=1, inter = l_fl_inter, \
                sample = "*", naverage = 1, function = l_function, \
                order = l_order, low_reject = 3., high_reject = 3., \
                niterate = 3, grow = 0.) 

            # Generate the nomralizn image
            
            if (1 == dispaxis) {
                if (debug) print ("expand")
                imstack (projfit, tmpnorm)
                keypar (dsub, "i_naxis2", silent+)
                if (no == keypar.found) {
                    printlog ("ERROR - NSFLAT: No i_naxis2 in " \
                        // dsub // ".", l_logfile, l_verbose)
                    goto clean
                }
                ny = int (keypar.value)
                blkrep (tmpnorm, tmpnorm, 1, ny)
            } else {
                if (debug) print ("rotate + expand")
                imstack (projfit, tmpnorm)
                imtranspose (tmpnorm, tmpnorm)

                keypar (dsub, "i_naxis1", silent+)
                if (no == keypar.found) {
                    printlog ("ERROR - NSFLAT: No i_naxis1 in " \
                        // dsub // ".", l_logfile, l_verbose)
                    goto clean
                }
                nx = int (keypar.value)
                blkrep (tmpnorm, tmpnorm, nx, 1)
            }

        } else if ("median" == l_process) {

            # Normalization from a median filter

            if (1 == dispaxis) {
                wx = l_box_length
                wy = l_box_width
            } else {
                wx = l_box_width
                wy = l_box_length
            }
            if (debug) print ("median: " // wx // "x" // wy)
            median (tmpsubsec, tmpnorm, xwin=wx, ywin=wy, zloreject=0, \
                zhireject=satsci, boundary="reflect")

        } else if ("trace" == l_process) {

            # Normalization by fitting along the traced aperture

            if ("" == l_trace)
                tracename = l_database // "/ap" // tmpsubsec
            else
                tracename = l_database // "/ap" // l_trace // "_" \
                    // l_sci_ext // "_" // version // "_"
            if (debug) print (tracename)

            # Read params from aptable
            aprow = mktemp ("aprow")//".fits"
            search = "prism == '" // prism // "' && " \
                // "slit == '" // fpmask // "' && " \
                // "decker == '" // decker // "' && " \
                // "arrayid == '" // arrayid // "' && " \
                // "order == " // version
            if (debug) print (search)
            tselect (l_aptable, aprow, search)

            apline = 1
            apthresh = 1
            florder = l_order
            apsample = "*"
            apupper = l_ap_upper
            aplower = l_ap_lower
            sdsample = "*"
            bgsample = "*"
            trproc = l_traceproc

            tinfo (table = aprow, ttout-, >& "dev$null")
            if (1 != tinfo.nrows) {

                printlog ("WARNING - NSFLAT: No aperture parameters \
                    in " // l_aptable, l_logfile, verbose+)
                printlog ("                  for " // search, \
                    l_logfile, verbose+) 
                printlog ("                  - flat quality may be \
                    poor", l_logfile, verbose+) 
                if (l_fl_inter || access (tracename)) {
                    printlog ("                (run with fl_inter+ or \
                        update/change the aperture table)", \
                        l_logfile, verbose+) 
                }

            } else {

                apsample = ""
                if (debug) 
                tprint (aprow, prparam-, prdata+, showrow-, showhdr-,
                    showunits-, col="apline,apthresh,florder,apsample,apupper,\
                    aplower,sdsample,bgsample,trproc", rows=1, pwidth=160)
                tprint (aprow, prparam-, prdata+, showrow-, showhdr-,
                    showunits-, 
                    col="apline,apthresh,florder,apsample,apupper,aplower,\
                    sdsample,bgsample,trproc", rows=1, pwidth=160) \
                    | scan (apline, apthresh, florder, apsample, apupper, \
                    aplower, sdsample, bgsample, trproc)

                if ("" == apsample) {
                    printlog ("ERROR - NSFLAT: No aperture info in MDF.", \
                        l_logfile, verbose+)
                    goto clean
                }

            }

            delete (aprow, verify-, >& "dev$null")

            if ("" == l_trace) {
                tracename = tmpsubsec
                width = apupper - aplower
                if (debug) print ("width " // width)
            } else {
                tracename = l_trace // "[" // l_sci_ext // "," // version // \
                    "]"
                width = 5
            }
            if (debug) print (tracename)

            gemhedit (tmpsubsec, l_key_dispaxis, str (dispaxis), "", delete-)

            # additional processing for trace file
            if (no == ("none" == trproc)) {
                imdelete (tmptrace, verify-, >& "dev$null")
                tmptrace = mktemp ("tmptrace")
                if ("right" == trproc) {
                    if (1 == dispaxis) {
                        imexpr ("a*J", tmptrace, tracename, dims="auto", \
                            intype="auto", outtype="auto", refim="auto", \
                            bwidth=0, btype="nearest", bpixval=0, \
                            rangecheck+, verbose-, exprdb="none")
                    } else {
                        imexpr ("a*I", tmptrace, tracename, dims="auto", \
                            intype="auto", outtype="auto", refim="auto", \
                            bwidth=0, btype="nearest", bpixval=0, \
                            rangecheck+, verbose-, exprdb="none")
                    }
                } else if ("left" == trproc) {
                    if (1 == dispaxis) {
                        hselect (tracename, "naxis2", yes) | scan (ipeak)
                        imexpr ("a*(" // ipeak // "-J)", tmptrace, \
                            tracename, dims="auto", \
                            intype="auto", outtype="auto", refim="auto", \
                            bwidth=0, btype="nearest", bpixval=0, \
                            rangecheck+, verbose-, exprdb="none")
                    } else {
                        hselect (tracename, "naxis1", yes) | scan (ipeak)
                        imexpr ("a*(" // ipeak // "-I)", tmptrace, \
                            tracename, dims="auto", \
                            intype="auto", outtype="auto", refim="auto", \
                            bwidth=0, btype="nearest", bpixval=0, \
                            rangecheck+, verbose-, exprdb="none")
                    }
                } else {
                    if (1 == dispaxis) {
                        hselect (tracename, "naxis2", yes) | scan (ipeak)
                        gauss (tracename, tmptrace, sigma=ipeak/20, \
                            ratio=0.05, theta=90, nsigma=4, bilinear+, \
                            boundary="constant", constant=0)
                    } else {
                        hselect (tracename, "naxis1", yes) | scan (ipeak)
                        gauss (tracename, tmptrace, sigma=ipeak/20, \
                            ratio=0.05, theta=0, nsigma=4, bilinear+, \
                            boundary="constant", constant=0)
                    }
                    if ("centre" == trproc) {
                        if (l_fl_inter) {
                            printlog ("WARNING - NSFLAT: Trace profile \
                                has an additional smoothed component.", \
                                l_logfile, verbose+) 
                        }
                        tmptrace2 = mktemp ("tmptrace2")
                        imexpr ("a + 2*b", tmptrace2, tracename, \
                            tmptrace, dims="auto", \
                            intype="auto", outtype="auto", refim="auto", \
                            bwidth=0, btype="nearest", bpixval=0, \
                            rangecheck+, verbose-, exprdb="none")
                            imdelete (tmptrace, verify-, >& "dev$null")
                            tmptrace = tmptrace2
                    }
                }
                tracename = tmptrace
            }


            dispaxis_copy = no
            keypar (tracename, l_key_dispaxis, silent+)
            if (no == keypar.found) {
                keypar (tracename//"[0]", l_key_dispaxis, silent+)
                if (keypar.found) {
                    if (debug) print ("copying dispaxis: " // keypar.value)
                    dispaxis_copy = yes
                    gemhedit (tracename, l_key_dispaxis, keypar.value, "", \
                        delete-)
                }
            }

            apall (tracename, output="", apertures="", \
                format="multispec", references="", profile="", \
                inter=l_fl_inter, find=yes, recenter=yes, resize=no, \
                edit=l_fl_inter, trace=yes, fittrace=yes, extract=no, \
                extras=no, review=no, line=apline, nsum=20, \
                lower=aplower, upper=apupper, width=width, \
                radius=2*width, thresho=l_thresh, nfind=1, minsep=5.,
                maxsep=1000., order="increasing",
                t_nsum=l_apsum, t_step=l_tr_step, t_nlost=l_tr_nlost, 
                t_funct=l_tr_function, t_order=l_tr_order, t_sample=sdsample, 
                t_naver=l_tr_naver, t_niter=l_tr_niter, t_low_r=l_tr_lowrej, 
                t_high_r=l_tr_highrej, t_grow=l_tr_grow, background="none", 
                b_sample=bgsample, >> tmpapalllog)

            if (dispaxis_copy) {
                if (debug) print ("erasing dispaxis copy")
                gemhedit (tracename, l_key_dispaxis, "", "", delete+)
            }

            apflatten (tmpsubsec, tmpnorm, apertures="", \
                reference=tracename, interactive=l_fl_inter, find=no, \
                recenter=no, resize=no, edit=no, trace=no, fittrace=no, \
                flatten=yes, fitspec=yes, line=apline, nsum=l_apsum, \
                threshold=apthresh, pfit="fit1d", clean=no, \
                saturation=INDEF, readnoise="0.", gain="1.", lsigma=4., \
                usigma=4., function=l_function, order=florder, \
                sample=apsample, naverage=5, niterate=3, \
                low_reject=4., high_reject=4., grow=0., >> tmpapalllog)

        } else {
            printlog ("ERROR - NSFLAT: Unexpected process: " \
                // l_process, l_logfile, l_verbose)
            goto clean
        }

        # Copy subsection back into appropriate data

        if ("" != fitused) {
            tmpfullnorm = mktemp("tmpfullnorm")
            if (debug) print ("empty image")
            imcalc (dsub, tmpfullnorm, "1", verbose-)
            if (debug) print ("copy across")
            if (debug) imhead (tmpfullnorm)
            if (1 == dispaxis) {
                tmpstr = substr(fitused, 1, stridx(",", fitused)-1)
                tmpsec = tmpstr//",*]"
            } else {
                tmpstr = substr(fitused, stridx(",",fitused)+1, 
                    strlen(fitused))
                tmpsec = "[*,"//tmpstr
            }
            imcopy (tmpnorm, tmpfullnorm // tmpsec, verbose-)
        } else {
            tmpfullnorm = tmpnorm
        }

        # Divide by fit/smoothed data

        if ("trace" == l_process) {
            # Already normalized, but need variance
            # We infer what the scaling was from data to flat by
            # dividing the data by the flat (thanks Frank!)
            if (l_fl_vardq) {
                imexpr ("a * (c/b)**2", ratio, esub, dsub, tmpfullnorm, 
                    verbose-)
                imdelete (esub, verify-, > "dev$null")
                imrename (ratio, esub, verbose-)
            }
            imdelete (dsub, verify-, >& "dev$null")
            dsub = mktemp ("tmpdsub")
            imcopy (tmpfullnorm, dsub, verbose-)
        } else {
            if (debug) print ("normalize")
            imdelete (ratio, verify-, >& "dev$null")
            if (debug) print ("tmpfullnorm: divisor: 3")
            imarith (dsub, "/", tmpfullnorm, ratio, verbose-)
            imdelete (dsub, verify-, > "dev$null")
            dsub = mktemp ("tmpdsub")
            imrename (ratio, dsub, verbose-)
            if (l_fl_vardq) {
                imarith (esub, "/", tmpfullnorm, ratio, verbose-) 
                imdelete (esub, verify-, > "dev$null")
                imrename (ratio, esub, verbose-)
                imarith (esub, "/", tmpfullnorm, ratio, verbose-) 
                imdelete (esub, verify-, > "dev$null")
                imrename (ratio, esub, verbose-)
            }
        }

        # Get normalization and good pixel limits for flat

        if (debug) print ("normzn, stats")
        norm = 0.0
        dqmask = dqsum

        if (l_fl_inter && l_fl_range) {
            # Do the interactive thing
            imhistogram (dsub // statused, z1 = INDEF, z2 = INDEF, \
                binwidth = INDEF, nbins = 512, autoscale+, \
                top_closed-, hist_type = "normal", listout-, \
                plot_type = "line", logy+, device = "stdgraph") 
            print ("Flat histogram") 
            print ("Hit any key to mark lower limit for good pixels") 
            dum = fscan (gcur, lower) 
            print ("Hit any key to mark upper limit for good pixels") 
            dum = fscan (gcur, upper) 
            print (" ") 
            if (lower > upper) {
                temp = upper
                upper = lower
                lower = temp
            }

            if (l_fl_vardq) {
                mimstat (dsub // statused, fields = l_normstat, \
                    imasks = dqmask // statused, omasks = "", \
                    lower = lower, upper = upper, nclip = 0, \
                    lsigma = INDEF, usigma = INDEF, binwidth = 0.1, \
                    format-, cache-) | scan (norm) 
            } else {
                imstat (dsub // statused, fields = l_normstat, \
                    lower = lower, upper = upper, nclip = 0, \
                    lsigma = INDEF, usigma = INDEF, binwidth = 0.1, \
                    format-, cache-) | scan (norm) 
            }

        } else {

            # Use fixed values, clip by 3-sig
            if (debug) imstat (dsub // statused)
            if (debug) imstat (dsub // statused, fields = "mode,stddev")
            if (l_fl_vardq) {
                mimstat (dsub // statused, fields = "mode,stddev", \
                    imasks = dqmask // statused, omasks = "", \
                    lower = INDEF, upper = satsci, nclip = 0, \
                    lsigma = INDEF, usigma = INDEF, binwidth = 0.1, \
                    format-, cache-) | scan (norm, stddev) 
                if (debug)
                    imstat (dsub // statused, lower = (norm-stddev*3),
                        upper = (norm+stddev*3))
                mimstat (dsub // statused, fields = "midpt", \
                    imasks = dqmask // statused, omasks = "", \
                    lower = (norm-stddev*3), upper = (norm+stddev*3), \
                    nclip = 0, lsigma = INDEF, usigma = INDEF, \
                    binwidth = 0.1, format-, cache-) | scan (norm)
            } else {
                imstat (dsub // statused, fields = "mode,stddev", \
                    lower = INDEF, upper = satsci, nclip = 0, \
                    lsigma = INDEF, usigma = INDEF, binwidth = 0.1, \
                    format-, cache-) | scan (norm, stddev) 
                if (debug)
                    imstat (dsub // statused, lower = (norm-stddev*3),
                        upper = (norm+stddev*3))
                imstat (dsub // statused, fields = "midpt", \
                    lower = (norm-stddev*3), upper = (norm+stddev*3), \
                    nclip = 0, lsigma = INDEF, usigma = INDEF, \
                    binwidth = 0.1, format-, cache-) | scan (norm)
            }
            if (norm == INDEF) {
                printlog ("WARNING - NSFLAT: Statistics failed, possibly \
                    due to a bad statsec.", l_logfile, verbose+) 
                do_norm = no
            } else {
                lower = norm*l_thr_flo
                upper = norm*l_thr_fup
                if (l_fl_vardq) {
                    mimstat (dsub // statused, fields = l_normstat, \
                        imasks = dqmask // statused, omasks = "", \
                        lower = lower, upper = upper, nclip = 0, \
                        lsigma = INDEF, usigma = INDEF, binwidth = 0.1, \
                        form-, cache-) | scan (norm) 
                } else {
                    imstat (dsub // statused, fields = l_normstat, \
                        lower = lower, upper = upper, nclip = 0, \
                        lsigma = INDEF, usigma = INDEF, binwidth = 0.1, \
                        form-, cache-) | scan (norm) 
                }
                if (norm == INDEF) {
                    printlog ("WARNING - NSFLAT: Statistics failed, \
                        possibly due to a bad statsec.", l_logfile, verbose+) 
                    do_norm = no
                }
            }
        }

        # Combine masks
        
        if (debug) print ("combine masks")
        if (l_fl_vardq) {
            
            if (do_norm) {
                imexpr ("((a<" // lower // ") || (a>" // upper // ") || \
                    (b>0)) ? 1 : 0", mask[1], dsub, dqsum, outtype = "ushort",
                    verbose-) 
                if (debug)
                    print ("dq after clip " // lower // "/" // upper)
                if (debug) imstat (mask[1])

                imstat (mask[1] // fitused, fields="mean", lower=INDEF, \
                    upper=INDEF, nclip=0, lsigma=INDEF, usigma=INDEF, \
                    binwidth=0.1, format-, cache-) | scan (temp)
                if ((temp-dqbase)/(1.0-dqbase) > dqlimit) {
                    printlog ("WARNING - NSFLAT: DQ for flat is poor \
                        (thresh: " // lower // " to " // upper // ")", \
                        l_logfile, verbose+)
                }

            } else {
                mask[1] = dqsum
                printlog ("WARNING - NSFLAT: No flagging of out-of-range \
                    flat values (bad stats).", l_logfile, verbose+)
            }

            if (debug) print ("mask stats")
            if (debug) imstat (mask[1])

            if (fl_dark) {
                if (debug) imstat (mask[2])
                addmasks (mask[1] // "," // mask[2], tmpbpmfile // ".pl", 
                    "im1 || im2", flags = " ") 
            } else {
                imcopy (mask[1], tmpbpmfile // ".pl", verbose-) 
            }
            if (debug) print ("final dq stats")
            if (debug) imstat (tmpbpmfile // ".pl")
            dqmask = tmpbpmfile // ".pl"
        }

        # Normalize the flat so that it has a midpt of 1, given
        # the clipping etc above
        if (do_norm) {
            tmpdiv = mktemp ("tmpdiv")
            if (debug) print ("final normalizn " // dsub // ", " \
                // norm // ", " // tmpdiv)
            if (debug) dir (dsub // "*")
            # imarith here crashed - don't know why
            imexpr ("a/" // norm, tmpdiv, dsub, verbose-)
            imdelete (dsub, verify-, >& "dev$null")
            dsub = mktemp ("tmpdsub")
            imrename (tmpdiv, dsub)
        } else {
            printlog ("WARNING - NSFLAT: Cannot normalize, due to bad \
                stats", l_logfile, verbose+)
        }

        # Fix bad pixels in flat image
        if (l_fl_fixbad) {
            imcopy (dsub, tmpsci, verbose-) 
            imdelete (dsub, verify-, >& "dev$null") 
            dsub = mktemp ("tmpdsub")
            imexpr ("(a==1) ? " // l_fixvalue // " : b", dsub, \
                tmpbpmfile, tmpsci, outtype = "real", refim="b", verbose-) 
        }

        # Report on final stats - use mimarith if we have DQ masking
        # (nicer to avoid "ifs" with empty mask in other case - tidy!)
        if (debug) print ("imstat 1 " // satsci // ", " // dqmask)
        if (l_fl_vardq) {

            # intermittent crash occurs here (with debug)
            # files involved are
            #  dsub
            #  dqmask (error still occurs if this is not used)
            if (debug) {
                mimstat (dsub // statused, fields = "midpt,stddev", \
                    imasks = dqmask // statused, omasks = "", \
                    lower = INDEF, upper = satsci, nclip = 0, \
                    lsigma = INDEF, usigma = INDEF, binwidth = 0.1, \
                    format-, cache-)
            }
            mimstat (dsub // statused, fields = "midpt,stddev", \
                imasks = dqmask // statused, omasks = "", \
                lower = INDEF, upper = satsci, nclip = 0, lsigma = INDEF, \
                usigma = INDEF, binwidth = 0.1, format-, cache-) \
                | scan (temp, stddev) 
        } else {
            if (debug) 
                imstat (dsub // statused, fields = "midpt,stddev", \
                    lower = INDEF, upper = satsci, nclip = 0, \
                    lsigma = INDEF, usigma = INDEF, binwidth = 0.1, \
                    format-, cache-)
                imstat (dsub // statused, fields = "midpt,stddev", \
                    lower = INDEF, upper = satsci, nclip = 0, lsigma = INDEF, \
                    usigma = INDEF, binwidth = 0.1, format-, cache-) \
                    | scan (temp, stddev) 
        }
        if (debug && isindef (temp)) print ("indef temp")
        if (debug && isindef (stddev)) print ("indef stddev")
        if (debug && no == isindef (temp) && no == isindef (stddev))
            print ("imstat 2 " // temp //", " // stddev)
        if (l_fl_vardq) {
            if (debug) 
                mimstat (dsub // statused, fields = "mean", \
                    imasks = dqmask // statused, omasks = "", \
                    lower = (temp-4*stddev), upper = (temp+4*stddev), \
                    nclip = 0, lsigma = INDEF, usigma = INDEF, \
                    binwidth = 0.1, format-, cache-)
            mimstat (dsub // statused, fields = "mean", \
                imasks = dqmask // statused, omasks = "", \
                lower = (temp-4*stddev), upper = (temp+4*stddev), \
                nclip = 0, lsigma = INDEF, usigma = INDEF, \
                binwidth = 0.1, format-, cache-) | scan (mean) 
        } else {
            if (debug)
                imstat (dsub // statused, fields = "mean", \
                    lower = (temp-4*stddev), upper = (temp+4*stddev), \
                    nclip = 0, lsigma = INDEF, usigma = INDEF, \
                    binwidth = 0.1, format-, cache-)
            imstat (dsub // statused, fields = "mean", \
                lower = (temp-4*stddev), upper = (temp+4*stddev), \
                nclip = 0, lsigma = INDEF, usigma = INDEF, \
                binwidth = 0.1, format-, cache-) | scan (mean) 
        }

        if (l_fl_vardq) {
            if (debug) print ("var stats:")
            if (debug) print (dqmask // statused)
            if (debug) imstat (dqmask // statused)
            if (debug)
                mimstat (esub // statused, \
                    fields = "image,npix,stddev,min,max,midpt,stddev", \
                    imasks = dqmask // statused, omasks = "", \
                    lower = INDEF, upper = satsci, nclip = 0, \
                    lsigma = INDEF, usigma = INDEF, binwidth = 0.1, \
                    format+, cache-)
            mimstat (esub // statused, fields = "midpt,stddev", \
                imasks = dqmask // statused, omasks = "", \
                lower = INDEF, upper = satsci, nclip = 0, lsigma = INDEF, \
                usigma = INDEF, binwidth = 0.1, format-, cache-) \
                | scan (temp, stddev) 
            if (debug)
                mimstat (esub // statused, \
                    fields = "image,npix,stddev,min,max,midpt,stddev,mean", \
                    imasks = dqmask // statused, omasks = "", \
                    lower = (temp-4*stddev), upper = (temp+4*stddev), \
                    nclip = 0, lsigma = INDEF, usigma = INDEF, \
                    binwidth = 0.1, format+, cache-)
            mimstat (esub // statused, fields = "mean", \
                imasks = dqmask // statused, omasks = "", \
                lower = (temp-4*stddev), upper = (temp+4*stddev), \
                nclip = 0, lsigma = INDEF, usigma = INDEF, \
                binwidth = 0.1, format-, cache-) | scan (stddev) 

            if (0.000001 < stddev) {
                s2n = mean / sqrt(stddev)
                printf ("%7.2f %5.1f  \n", mean, s2n) | scan (sline)
            } else {
                printf ("%7.2f %5s  \n", mean, ">1000") | scan (sline)
            }

        } else {

            printf ("%7.2f %5s  \n", mean, "---") | scan (sline)
            printlog ("WARNING - NSFLAT: No statistics for S/N \
                (fl_vardq=no).", l_logfile, verbose+) 

        }

        statline = statline // sline

        printlog (statline, l_logfile, l_verbose)

        # Copy to output MEF file

        if (debug) print ("final mef file")

        gemhedit (dsub, "i_title", l_flattitle, "", delete-)
        imgout = l_flatfile // "[" // l_sci_ext // "," // version // ",append]"
        imcopy (dsub, imgout, verbose-) 

        if (l_fl_vardq) {
            gemhedit (esub, "i_title", l_flattitle, "", delete-)
            imgout = l_flatfile // "[" // l_var_ext // "," // \
                version // ",append]"
            imcopy (esub, imgout, verbose-)

            gemhedit (tmpbpmfile, "i_title", l_bpmtitle, "", delete-)
            imgout = l_flatfile // "[" // l_dq_ext // "," // \
                version // ",append]"
            # Ensure that the output DQ extension has the correct pixel type
            chpixtype (tmpbpmfile // ".pl", imgout, newpixtype="ushort", \
                oldpixtype="all", verbose-)
        }

        # Update headers in final image PHU (call mkheader many times
        # as there's currently a bug that limits the maximum size)

        if (1 == version) {

            if (debug) print ("PHU headers")

            gemdate ()
            phu = l_flatfile // ".fits[0]"
            gemhedit (phu, "GEM-TLM", gemdate.outdate, 
                "UT Last modification with GEMINI", delete-) 

            printf ("%-8s= \'%-18s\' / %-s\n", "NSFLAT", gemdate.outdate, \
                "UT Time stamp for NSFLAT", >> tmphead) 
            if (do_norm)
                printf ("%-8s= %20.5f / %-s\n", "NSFLTNRM", norm, \
                    "Normalization constant for NSFLAT", >> tmphead) 
            printf ("%-8s= \'%-18s\' / %-s\n", "NSFLTNST", l_normstat, \
                "Normalization statistic for NSFLAT", >> tmphead) 
            if ("fit" == l_process || "trace" == l_process) {
                printf ("%-8s= \'%-18s\' / %-s\n", "NSFLTFUN", \
                    l_function, \
                    "Function fit to illumination by NSFLAT", >> tmphead) 
                printf ("%-8s= \'%-18s\' / %-s\n", "NSFLTORD", l_order, \
                    "Order of function fit by NSFLAT", >> tmphead) 
            } else if ("median" == l_process) {
                printf ("%-8s= \'%-18s\' / %-s\n", "NSFLTFUN", "median", \
                    "Function fit to illumination by NSFLAT", >> tmphead) 
            }
            if (do_norm) {
                printf ("%-8s= %20.5f / %-s\n", "NSFLTFLO", lower, \
                    "Lower limit for good pixels (fraction)", >> tmphead) 
                printf ("%-8s= %20.5f / %-s\n", "NSFLTFUP", upper, \
                    "Upper limit for good pixels (fraction)", >> tmphead) 
            }
            printf ("%-8s= \'%-18s\' / %-s\n", "NSFLTSTA", statused, \
                "Statistics section used by NSFLAT", >> tmphead)
            printf ("%-8s= \'%-18s\' / %-s\n", "NSFLTFIT", fitused, \
                "Fit area ", >> tmphead) 

            if (debug) print ("mkheader: " // phu)
            if (debug) type (tmphead)
            mkheader (phu, tmphead, append+, verbose-) 
            delete (tmphead, verify-, >& "dev$null")

            if (fl_dark) {
                dark1 = (int (dark1*100.) ) /100.
                dark2 = (int (dark2*100.) ) /100.
                printf ("%-8s= %20.5f / %-s\n", "NSFLTDLO", dark1, \
                    "Lower limit for good pixels in dark (ADU)", >> tmphead) 
                printf ("%-8s= %20.5f / %-s\n", "NSFLTDUP", dark2, \
                    "Upper limit for good pixels in dark (ADU)", >> tmphead) 
            }
            if (l_fl_fixbad) {
                printf ("%-8s= %20.5f / %-s\n", "NSFLTFIX", l_fixvalue, \
                    "Bad pixels replaced with this value", >> tmphead) 
            }
            printf ("%-8s= \'%-18s\' / %-s\n", "NSFLTCOM", l_combtype, \
                "Combine method used by NSFLAT", >> tmphead) 
            if (l_zero != "none") {
                printf ("%-8s= \'%-18s\' / %-s\n", "NSFLTZER", l_zero, \
                    "Statistic used to compute additive offsets", \
                    >> tmphead) 
            }
            if (l_scale != "none") {
                printf ("%-8s= \'%-18s\' / %-s\n", "NSFLTSCL", l_scale, \
                    "Statistic used to compute scale factors", >> tmphead) 
            }
            if (l_weight != "none") {
                printf ("%-8s= \'%-18s\' / %-s\n", "NSFLTWEI", l_weight, \
                    "Statistic used to compute relative weights", \
                >> tmphead) 
            }

            if (access (tmphead)) {
                if (debug) print ("mkheader: " // phu)
                if (debug) type (tmphead)
                mkheader (phu, tmphead, append+, verbose-) 
                delete (tmphead, verify-, >& "dev$null")
            }

            if (l_rejtype != "none") {
                printf ("%-8s= \'%-18s\' / %-s\n", "NSFLTREJ", l_rejtype, \
                    "Rejection algorithm used by NSFLAT", >> tmphead) 
                printf ("%-8s= \'%-18s\' / %-s\n", "NSFLTLTH", \
                    str(l_lthreshold), "Lower threshold before combining", \
                    >> tmphead) 
                printf ("%-8s= \'%-18s\' / %-s\n", "NSFLTHTH", \
                    str(l_hthreshold), "Upper threshold before combining", \
                    >> tmphead) 
                printf ("%-8s= %20.5f / %-s\n", "NSFLTGRW", l_grow, \
                    "Radius for additional pixel rejection", >> tmphead) 
            }
            if (l_rejtype == "minmax") {
                printf ("%-8s= %20.0f / %-s\n", "NSFLTNLO", l_nlow, \
                    "Low pixels rejected (minmax)", >> tmphead) 
                printf ("%-8s= %20.0f / %-s\n", "NSFLTNHI", l_nhigh, \
                    "High pixels rejected (minmax)", >> tmphead) 
            }
            if (l_rejtype == "sigclip" || l_rejtype == "avsigclip" \
                || l_rejtype == "ccdclip" || l_rejtype == "crreject" \
                || l_rejtype == "pclip") {
                
                printf ("%-8s= %20.5f / %-s\n", "NSFLTLSI", l_lsigma, \
                    "Lower sigma for rejection", >> tmphead) 
                printf ("%-8s= %20.5f / %-s\n", "NSFLTHSI", l_hsigma, \
                    "Upper sigma for rejection", >> tmphead) 
                printf ("%-8s= %20.0f / %-s\n", "NSFLTNKE", l_nkeep, \
                    "Min(max) number of pixels to keep(reject)", >> tmphead) 
                printf ("%-8s= \'%-18b\' / %-s\n", "NSFLTMCL", l_mclip, \
                    "Use median in clipping algorithms", >> tmphead)
            } 

            if (access (tmphead)) {
                if (debug) print ("mkheader: " // phu)
                if (debug) type (tmphead)
                mkheader (phu, tmphead, append+, verbose-) 
                delete (tmphead, verify-, >& "dev$null")
            }

            if (l_rejtype == "sigclip" || l_rejtype == "avsigclip" \
                || l_rejtype == "ccdclip" || l_rejtype == "crreject") {
                
                printf ("%-8s= %20.5f / %-s\n", "NSFLTSSC", l_sigscale, \
                    "Tolerance for sigma clip scaling correction", >> tmphead) 
            }
            if (l_rejtype == "pclip") {
                printf ("%-8s= %20.5f / %-s\n", "NSFLTPCL", l_pclip, \
                    "Percentile clipping factor used by pclip", >> tmphead) 
            }
            if (l_rejtype == "ccdclip") {
                printf ("%-8s= \'%-18s\' / %-s\n", "NSFLTSNO", l_snoise, \
                    "Sensitivity noise (e) used by ccdclip", >> tmphead) 
            }

            if (access (tmphead)) {
                if (debug) print ("mkheader: " // phu)
                if (debug) type (tmphead)
                mkheader (phu, tmphead, append+, verbose-) 
                delete (tmphead, verify-, >& "dev$null")
            }

            i = 1
            scanin1 = infiles[1]
            while (fscan (scanin1, img) != EOF) {
                if (img != l_flatfile) {
                    printf ("%-8s= \'%-18s\' / %-s\n", 
                        "NSFLON" // str (i), img, "Input flat image",
                        >> tmphead) 
                    i += 1
                }
            }

            if (access (tmphead)) {
                if (debug) print ("mkheader: " // phu)
                if (debug) type (tmphead)
                mkheader (phu, tmphead, append+, verbose-) 
                delete (tmphead, verify-, >& "dev$null")
            }

            if (fl_dark) {
                i = 1
                scanin1 = infiles[2]
                while (fscan (scanin1, img) != EOF) {
                    if (img != l_flatfile) {
                        printf ("%-8s= \'%-18s\' / %-s\n", \
                            "NSFLOF" // str (i), img, "Input dark", >> tmphead)
                        i += 1
                    }
                }
            }

            if (access (tmphead)) {
                if (debug) print ("mkheader: " // phu)
                if (debug) type (tmphead)
                mkheader (phu, tmphead, append+, verbose-) 
                delete (tmphead, verify-, >& "dev$null")
            }

            # These are not relevant for spectral flats
            if (debug) print ("deleting header info")
            gemhedit (phu, l_key_ron, "", "", delete+)
            gemhedit (phu, l_key_sat, "", "", delete+)
            gemhedit (phu, l_key_nonlinear, "", "", delete+)
            gemhedit (phu, l_key_gain, "", "", delete+)
        }

        # Write the output BPM file, with complete header

        if (debug) print ("BPM file " // l_bpmfile)

        if (l_fl_vardq) {
            imgdq = l_flatfile // "[" // l_dq_ext // "," // version // "]"
            keypar (imgdq, "naxis1", silent+)
            xsize = int (keypar.value)
            keypar (imgdq, "naxis2", silent+)
            ysize = int (keypar.value)

            if (iscut) {
                keypar (phu, "ORIGXSIZ", silent+) 
                if (keypar.found) {
                    xsize = int (keypar.value) 
                    keypar (phu, "ORIGYSIZ", silent+) 
                    ysize = int (keypar.value) 
                } else {
                    printlog ("WARNING - NSFLAT: Can't determine original \
                        raw image size.", l_logfile, verbose+) 
                    printlog ("                  Using " // xsize // "x" \
                        // ysize // ".", l_logfile, verbose+) 
                }
            }

            if (substr (l_bpmfile, strlen (l_bpmfile) -2, \
                strlen (l_bpmfile) ) != ".pl") {
                
                l_bpmfile = l_bpmfile // ".pl"
            }
            section = "[*,*]"
            if (debug) print ("iscut? " // iscut)
            if (iscut) {
                img = mdftemplate // "[" // l_sci_ext // "," // version // "]"
                keypar (img, l_key_cut_section, silent+)
                if (keypar.found) {
                    section = keypar.value
                } else {
                    printlog ("WARNING - NSFLAT: Can't determine cut \
                        section from " // l_key_cut_section, \
                        l_logfile, verbose+) 
                    printlog ("                  in " // img, 
                        l_logfile, verbose+) 
                }
                if (debug)
                    print ("section: " // section // " " // keypar.found)
            }

            if (debug) {
                print ("copying " // section // " to file of \
                    size " // xsize // "x" // ysize)
                imhead (imgdq)
            }

            if (1 == version) {
                mkimage (l_bpmfile, "make", 0, ndim = 2, \
                    dims = xsize // " " // ysize, pixtype = "short") 
                imcopy (imgdq, l_bpmfile // section, verbose-) 

#               imgdq = l_flatfile // "[" // l_dq_ext // "," \
#                   // version // ",inherit]"
#               if (debug) print ("copying header from: " // imgdq)
#               keypar (imgdq, "DATE", silent+)
#               if (debug) print (keypar.value)
#               if (keypar.found)
#                   gemhedit (l_bpmfile, "DATE", keypar.value, \
#                       "Date FITS file was generated", delete-)
#               if (debug) print (keypar.value)
#               if (keypar.found)
#               gemhedit (l_bpmfile, "DATE", keypar.value, \
#                   "Date FITS file was generated", delete-)
            } else {
                mkimage (tmpdq, "make", 0, ndim = 2, \
                    dims = xsize // " " // ysize, pixtype = "short")
                imcopy (imgdq, tmpdq // section, verbose-)
                imrename (l_bpmfile, tmpbpm)
                imexpr ("a | b", l_bpmfile, tmpbpm, tmpdq, verbose-)
            }
        }
        imdelete (tmpfullnorm, verify-, >& "dev$null") 
    } # next version!


    if (access (l_bpmfile)) {
        # Copy header info intp BPM file
        delete (tmphead, verify-, >& "dev$null")
        imhead (l_flatfile // "[" // l_dq_ext // ",1,inherit]", 
            long+, > tmphead)
        mkheader (l_bpmfile, tmphead, append+, verbose-)
        delete (tmphead, verify-, >& "dev$null")
    }


    # If we get this far, the status should indicate success
    status = 0

clean:
    if (debug) print ("clean")

    if (l_verbose && access (tmpapalllog)) {
        printlog ("NSFLAT: Output from APALL/APFLATTEN:", l_logfile, l_verbose)
        scanin1 = tmpapalllog
        while (fscan (scanin1, sline) != EOF) {
            printlog (sline, l_logfile, l_verbose)
        }
    }

    # note that tmpdarks amy also be infiles[2], so don't alter order here
    if (access (tmpdarks)) {
        if (debug) print ("tmpdarks: ")
        if (debug) type (tmpdarks)
        if (no == l_fl_save_darks)
        imdelete ("@" // tmpdarks, verify-, default+, >& "dev$null")
        delete (tmpdarks, verify-, default+, >& "dev$null")
    }
    if (access (tmpmdfs)) {
        if (debug) print ("tmpmdfs: ")
        if (debug) type (tmpmdfs)
        imdelete ("@" // tmpmdfs, verify-, default+, >& "dev$null")
        delete (tmpmdfs, verify-, default+, >& "dev$null")
    }
    delete (tmpout // "," // tmphead // "," // tmpdarklist, \
        verify-, >& "dev$null") 
    for (i = 1; i <= 2; i += 1) {
        delete (infiles[i], verify-, >& "dev$null") 
        # save combined darks if necessary
        if (i == 1 || "" == l_darkfile)
            imdelete (combout[i], verify-, >& "dev$null") 
    }
    if (access (tmpmdf)) 
        delete (tmpmdf, verify-, default+, >& "dev$null")
    imdelete (dsub // "," // esub // "," // dqsum // "," \
        // tmpbpmfile // "," // mask[1] // "," // mask[2] \
        // "," // mask[3] // "," // tmpsci // "," // tmpdq \
        // "," // tmpbpm // "," // ratio // "," // tmpfullnorm \
        // "," // project // "," // projfit // "," // tmpnorm \
        // "," // tmpsubsec \
        // "," // tmptrace, verify-, >& "dev$null") 
    delete (aprow, verify-, >& "dev$null")
    scanin1 = ""
    scanin2 = ""
    specred.database = databasesave
    specred.logfile = logfilesave

    printlog (" ", l_logfile, l_verbose) 
    if (status == 0) {
        printlog ("NSFLAT exit status: good.", l_logfile, l_verbose) 
    } else {
        printlog ("NSFLAT exit status: error.", l_logfile, l_verbose) 
        if (status == 1) imdelete (l_flatfile, verify-, >& "dev$null") 
    }
    printlog ("---------------------------------------------------------\
        ----------------------", l_logfile, l_verbose) 

end

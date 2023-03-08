# Copyright(c) 2001-2011 Association of Universities for Research in Astronomy, Inc.

procedure nsslitfunction (skyflats, output) 

# makes an illumination corrected GNIRS/NIRI spectroscopic flat field

string  skyflats    {prompt = "Input sky flats"}
string  output      {prompt = "Output normalized flat (MEF)"}
string  flat        {"", prompt = "Input lamp flat"} 
string  flexflat    {"", prompt = "Optional lamp flat at sky position"}
string  dark        {"", prompt = "Input dark image (at sky position, if flexflat used)"}
string  title       {"", prompt = "Title for output SCI plane"}
string  combine     {"median", prompt = "Combination operation"}
string  reject      {"none", prompt = "Rejection algorithm"}
string  masktype    {"goodvalue", enum="none|goodvalue", prompt="Mask type"}
real    maskvalue   {0., prompt="Mask value"}
string  scale       {"median", prompt = "Image scaling"}
string  zero        {"none", prompt = "Image zeropoint offset"}
string  weight      {"none", prompt = "Image weights"}
string  statsec     {"[*,*]", prompt = "Statistics section"}
real    lthreshold  {INDEF, prompt = "Lower threshold"}
real    hthreshold  {INDEF, prompt = "Upper threshold"}
int     nlow        {1, min=1, prompt = "minmax: Number of low pixels to reject"}
int     nhigh       {1, min=1, prompt = "minmax: Number of high pixels to reject"}
int     nkeep       {0, min=0, prompt = "Minimum to keep or maximum to reject"}  
bool    mclip       {yes, prompt = "Use median in sigma clipping algorithms?"}
real    lsigma      {3., prompt = "Lower sigma clipping factor"}
real    hsigma      {3., prompt = "Upper sigma clipping factor"}
real    ron         {0.0, min = 0., prompt = "Readout noise rms in electrons"}
real    gain        {1.0, min = 0.00001, prompt = "Gain in e-/ADU"}
string  snoise      {"0.0", prompt = "ccdclip: Sensitivity noise (electrons)"}
real    sigscale    {0.1, prompt = "Tolerance for sigma clipping scaling correction"}
real    pclip       {-0.5, prompt = "pclip: Percentile clipping parameter"}
real    grow        {0.0, prompt = "Radius (pixels) for neighbor rejection"}
string  function    {"spline3", prompt = "Fitting function"}
int     order       {3, min = 1, prompt = "Order of fitting function"}
bool    fl_vary     {no, prompt = "Allow scaling to vary with order?"}
bool    fl_dqprop   {no, prompt = "Propagate all DQ values?"}
bool    fl_inter    {no, prompt = "Fit interactively?"}
bool    fl_gemarith {yes, prompt = "Use gemarith (faster)?"}
string  logfile     {"", prompt = "Logfile name"}
bool    verbose     {yes, prompt = "Verbose output?"}
int     status      {0, prompt = "Exit status (0=good)"}
struct  *scanfile   {"", prompt = "Internal use only"} 

begin

    string  l_skyflats = ""
    string  l_dark = ""
    string  l_flat = ""
    string  l_flexflat = ""
    string  l_output = ""
    string  l_title = ""
    string  l_combine = ""
    string  l_reject = ""
    string  l_masktype = ""
    real    l_maskvalue
    string  l_scale = ""
    string  l_zero = ""
    string  l_weight = ""
    string  l_statsec = ""
    real    l_lthreshold
    real    l_hthreshold
    int     l_nlow
    int     l_nhigh
    int     l_nkeep
    bool    l_mclip
    real    l_lsigma
    real    l_hsigma
    real    l_ron
    real    l_gain
    string  l_snoise = ""
    real    l_sigscale
    real    l_pclip
    real    l_grow
    string  l_function = ""
    int     l_order
    bool    l_fl_vary
    bool    l_fl_dqprop
    bool    l_fl_inter
    bool    l_fl_gemarith
    string  l_logfile = ""
    bool    l_verbose

    string  l_sci_ext
    string  l_var_ext
    string  l_dq_ext
    string  l_key_ron
    string  l_key_gain
    string  l_key_dispaxis
    string  l_key_exptime

    string  filelist, scilist, dqlist, specsec1, sflat, phulist
    string  img, dqsum, dqsumold, suf, response, illumvar, tmplog, phu
    string  tmpdq, combout, combsig, combdq, flatsub, illum, badhdr
    string  tmplis1, tmplis2, tmplis3, tmplis, tmpsuper, tmpshift
    real    norm, stddev
    int     i, k, n, j, l, idum, nbad, len, junk, nsci, iextn, nextn
    int     iver, nver, nx, ny, dispaxis
    string  keyfound
    struct  sdate, line
    bool    debug, first, ok, havedark, havedqvar
    string  ssrclist, sdestlist, vsrclist, vdestlist, comblist, scratch
    string  dqsrclist, dqdestlist
    string  ssrc, sdest, vsrc, vdest, extn, skyextn, darkextn, flatextn
    string  comblist2, dqsrc, dqdest, sect
    real    mdresp, mddata, ratio, correcn, oldcorrecn, oldmddata
    real    shiftx, shifty, shx, shy, exptime, prevexp, sigma
    real    flx_shx, flx_shy, dx, dy
    char    comboutptr, previmg, expr


    junk = fscan (  skyflats, l_skyflats)
    junk = fscan (  dark, l_dark)
    junk = fscan (  flat, l_flat)
    junk = fscan (  flexflat, l_flexflat)
    junk = fscan (  output, l_output)
    junk = fscan (  title, l_title)
    junk = fscan (  combine, l_combine)
    junk = fscan (  reject, l_reject)
    junk = fscan (  masktype, l_masktype)
    l_maskvalue =   maskvalue
    junk = fscan (  scale, l_scale)
    junk = fscan (  zero, l_zero)
    junk = fscan (  weight, l_weight)
    junk = fscan (  statsec, l_statsec)
    l_lthreshold =  lthreshold
    l_hthreshold =  hthreshold
    l_nlow =    nlow
    l_nhigh =   nhigh
    l_nkeep =   nkeep
    l_mclip =   mclip
    l_lsigma =  lsigma
    l_hsigma =  hsigma
    l_ron =     ron
    l_gain =    gain
    junk = fscan (  snoise, l_snoise)
    l_sigscale =    sigscale
    l_pclip =       pclip
    l_grow =        grow
    junk = fscan (  function, l_function)
    l_order =       order
    l_fl_vary =     fl_vary
    l_fl_dqprop =   fl_dqprop
    l_fl_inter =    fl_inter
    l_fl_gemarith = fl_gemarith
    junk = fscan (  logfile, l_logfile)
    l_verbose =     verbose

    badhdr = ""
    junk = fscan (  nsheaders.sci_ext, l_sci_ext)
    if ("" == l_sci_ext) badhdr = badhdr + " sci_ext"
    junk = fscan (  nsheaders.var_ext, l_var_ext)
    if ("" == l_var_ext) badhdr = badhdr + " var_ext"
    junk = fscan (  nsheaders.dq_ext, l_dq_ext)
    if ("" == l_dq_ext) badhdr = badhdr + " dq_ext"
    junk = fscan (  nsheaders.key_ron, l_key_ron)
    if ("" == l_key_ron) badhdr = badhdr + " key_ron"
    junk = fscan (  nsheaders.key_gain, l_key_gain)
    if ("" == l_key_gain) badhdr = badhdr + " key_gain"
    junk = fscan (  nsheaders.key_dispaxis, l_key_dispaxis)
    if ("" == l_key_dispaxis) badhdr = badhdr + " key_dispaxis"
    junk = fscan (  nsheaders.key_exptime, l_key_exptime)
    if ("" == l_key_exptime) badhdr = badhdr + " key_exptime"

    status = 1
    debug = no

    havedark = no
    havedqvar = yes
    prevexp = INDEF
    shiftx = INDEF; shifty = INDEF
    flx_shx = INDEF; flx_shy = INDEF
    if ("" == l_flexflat) l_flexflat = l_flat

    cache ("keypar", "gnirs", "gemextn", "gemdate")

    filelist = mktemp ("tmpfilelist") 
    phulist = mktemp ("tmpphulist") 
    if (debug) print (phulist)
    scilist = mktemp ("tmpscilist") 
    combout = mktemp ("tmpcombout") 
    flatsub = mktemp ("tmpflatsub") 
    sflat = mktemp ("tmpsflat") 
    response = mktemp ("tmpresponse") 
    ssrclist = mktemp ("tmpssrclist")
    sdestlist = mktemp ("tmpsdestlist")
    vsrclist = mktemp ("tmpvsrclist")
    vdestlist = mktemp ("tmpvdestlist")
    dqsrclist = mktemp ("tmpdqsrclist")
    dqdestlist = mktemp ("tmpdqdestlist")
    comblist = mktemp ("tmpcomblist")
    comblist2 = mktemp ("tmpcomblist2")
    scratch = mktemp ("tmpscratch")
    tmplis = mktemp ("tmplis1")
    tmplis1 = mktemp ("tmplis1")
    tmplis2 = mktemp ("tmplis2")
    tmplis3 = mktemp ("tmplis3")
    tmpsuper = mktemp ("tmpsuper")
    tmpshift = mktemp ("tmpshift")


    # Start logging

    if ("" == l_logfile) {
        junk = fscan (gnirs.logfile, l_logfile) 
        if (l_logfile == "") {
            l_logfile = "gnirs.log"
            printlog ("WARNING - NSSLITFUNCTION: Both nsslitfunc.logfile \
                and gnirs.logfile are empty.", l_logfile, verbose+) 
            printlog ("                   Using default file " \
                // l_logfile // ".", l_logfile, verbose+) 
        }
    }
    date | scan (sdate) 
    printlog ("-----------------------------------------------------\
        -------------------------", l_logfile, verbose = l_verbose) 
    printlog ("NSSLITFUNCTION -- " // sdate, l_logfile, l_verbose) 
    printlog (" ", l_logfile, l_verbose) 
    printlog ("skyflats = " // l_skyflats, l_logfile, l_verbose) 
    printlog ("dark = " // l_dark, l_logfile, l_verbose) 
    printlog ("output = " // l_output, l_logfile, l_verbose) 
    printlog ("title = " // l_title, l_logfile, l_verbose) 
    printlog ("function = " // l_function, l_logfile, l_verbose) 
    printlog ("order = " // l_order, l_logfile, l_verbose) 
    printlog ("sci_ext = " // l_sci_ext, l_logfile, l_verbose) 
    printlog ("var_ext = " // l_var_ext, l_logfile, l_verbose) 
    printlog ("dq_ext = " // l_dq_ext, l_logfile, l_verbose) 
    printlog (" ", l_logfile, l_verbose) 

    if ("" != badhdr) {
        printlog ("ERROR - NSSLITFUNCTION: Parameter(s) missing from \
            nsheaders: " // badhdr, l_logfile, verbose+) 
        goto clean
    }


    # Check/expand sky flats

    gemextn (l_skyflats, check = "exists,mef", process = "none", \
        index = "", extname = "", extversion = "", ikparams = "", \
        omit = "", replace = "", outfile = filelist, logfile="",
        glogpars="", verbose=l_verbose)
    if (gemextn.fail_count != 0 || gemextn.count == 0) {
        printlog ("ERROR - NSSLITFUNCTION: Bad input images.", l_logfile, 
            verbose+) 
        goto clean
    }
    n = gemextn.count
    gemextn (l_skyflats, check = "exists,mef", process = "expand", \
        index = "0", extname = "", extversion = "", ikparams = "", \
        omit = "", replace = "", outfile = phulist, logfile="",
        glogpars="", verbose=l_verbose)
    if (gemextn.fail_count != 0 || gemextn.count != n) {
        printlog ("ERROR - NSSLITFUNCTION: Bad or missing PHUs.", \
            l_logfile, verbose+) 
        goto clean
    }
    head (filelist, nlines=1) | scan (img)
    gemextn (img, check = "", process = "expand", index = "", \
        extname = l_sci_ext, extversion = "1-", ikparams = "", \
        omit = "", replace = "", outfile = "dev$null", \
        logfile="dev$null", glogpars="", verbose-)
    nsci = gemextn.count


    # Check/expand dark
    # Doesn't matter if this doesn't have dq/var

    gemextn (l_dark, check = "exists,mef", process = "none", \
        index = "", extname = "", extversion = "", ikparams = "", \
        omit = "", replace = "", outfile = "dev$null", logfile="",
        glogpars="", verbose=l_verbose)
    if (gemextn.fail_count == 0 && gemextn.count == 0) {
        printlog ("WARNING - NSSLITFUNCTION: No dark supplied.", \
            l_logfile, verbose+)
    } else if (gemextn.fail_count != 0 || gemextn.count != 1) {
        printlog ("ERROR - NSSLITFUNCTION: Bad or missing dark.", \
            l_logfile, verbose+) 
        goto clean
    } else {
        havedark = yes
        gemextn (l_dark, check = "exists,mef", process = "none", \
            index = "", extname = "", extversion = "", ikparams = "", \
            omit = "", replace = "", outfile = "STDOUT", \
            logfile="dev$null", glogpars="", verbose-) | scan (line)
        l_dark = line
        gemextn (l_dark, check = "", process = "expand", index = "", \
            extname = l_sci_ext, extversion = "1-", ikparams = "", \
            omit = "", replace = "", outfile = "dev$null", \
            logfile="dev$null", glogpars="", verbose-)
        if (nsci != gemextn.count) {
            printlog ("ERROR - NSSLITFUNCTION: Skyflats have different no \
                of extensions (" // nsci // ")", l_logfile, verbose+) 
            printlog ("                        to combined dark (" \
                // gemextn.count // ")", l_logfile, verbose+) 
            goto clean
        }
        gemextn (l_dark, check = "exists,mef", process = "expand", \
            index = "0", extname = "", extversion = "", \
            ikparams = "", omit = "", replace = "", outfile = "STDOUT", \
            logfile="dev$null", glogpars="", verbose-) | scan (line)
        hselect (line, "MDF_XSHF", yes) | scan (flx_shx)
        hselect (line, "MDF_YSHF", yes) | scan (flx_shy)
        hselect (line, l_key_exptime, yes) | scan (prevexp)
        previmg = line
        if (isindef (prevexp)) {
            printlog ("WARNING - NSSLITFUNCTION: missing exposure time \
                in " // line, l_logfile, verbose+) 
        }
    }


    # Check/expand flat
    # This can be at a different shift to the dark/flexflat

    gemextn (l_flat, check = "exists,mef", process = "none", \
        index = "", extname = "", extversion = "", ikparams = "", \
        omit = "", replace = "", outfile = "dev$null", logfile="",
        glogpars="", verbose=l_verbose)
    if (gemextn.fail_count != 0 || gemextn.count != 1) {
        printlog ("ERROR - NSSLITFUNCTION: Bad or missing flat.", \
            l_logfile, verbose+) 
        goto clean
    }
    gemextn (l_flat, check = "exists,mef", process = "none", \
        index = "", extname = "", extversion = "", ikparams = "", \
        omit = "", replace = "", outfile = "STDOUT", \
        logfile="dev$null", glogpars="", verbose-) | scan (line)
    l_flat = line
    gemextn (l_flat, check = "", process = "expand", index = "", \
        extname = l_sci_ext, extversion = "1-", ikparams = "", \
        omit = "", replace = "", outfile = "dev$null", \
        logfile="dev$null", glogpars="", verbose-)
    if (nsci != gemextn.count) {
        printlog ("ERROR - NSSLITFUNCTION: Flats have different no \
            of extensions (" // nsci // ")", l_logfile, verbose+) 
        printlog ("                        to input flat (" \
            // gemextn.count // ")", l_logfile, verbose+) 
        goto clean
    }
    gemextn (l_flat, check = "exists,mef", process = "expand", \
        index = "0", extname = "", extversion = "", \
        ikparams = "", omit = "", replace = "", outfile = "STDOUT", \
        logfile="dev$null", glogpars="", verbose-) | scan (line)
    hselect (line, "MDF_XSHF", yes) | scan (shiftx)
    hselect (line, "MDF_YSHF", yes) | scan (shifty)

    if (havedqvar) {
        gemextn (l_flat, proc="expand", check="", index="", \
            extname=l_var_ext, extversion="-", ikparams="", omit="", \
            replace="", outfile="dev$null", glogpars="", verbose-)
        if (no == (nsci == gemextn.count)) havedqvar = no
        gemextn (l_flat, proc="expand", check="", index="", \
            extname=l_dq_ext, extversion="-", ikparams="", omit="", \
            replace="", outfile="dev$null", glogpars="", verbose-)
        if (no == (nsci == gemextn.count)) havedqvar = no
        if (no == havedqvar) {
            printlog ("WARNING - NSSLITFUNCTION: no dq/variance data \
                in " // l_flat, l_logfile, verbose+) 
        }
    }

    # Check/expand flexflat
    # This has to have the same shift as the dark
    # (if no flexflat is given, this is the same as the flat,
    # so we get the right checks "for free" but have to be 
    # careful that error messages make sense)
    # Doesn't matter if this doesn't have dq/var

    gemextn (l_flexflat, check = "exists,mef", process = "none", \
        index = "", extname = "", extversion = "", ikparams = "", \
        omit = "", replace = "", outfile = "dev$null", logfile="",
        glogpars="", verbose=l_verbose)
    if (gemextn.fail_count != 0 || gemextn.count != 1) {
        printlog ("ERROR - NSSLITFUNCTION: Bad or missing flexflat.", \
            l_logfile, verbose+) 
        goto clean
    }
    gemextn (l_flexflat, check = "exists,mef", process = "none", \
        index = "", extname = "", extversion = "", ikparams = "", \
        omit = "", replace = "", outfile = "STDOUT", \
        logfile="dev$null", glogpars="", verbose-) | scan (line)
    l_flexflat = line
    gemextn (l_flexflat, check = "", process = "expand", index = "", \
        extname = l_sci_ext, extversion = "1-", ikparams = "", \
        omit = "", replace = "", outfile = "dev$null", \
        logfile="dev$null", glogpars="", verbose-)
    if (nsci != gemextn.count) {
        printlog ("ERROR - NSSLITFUNCTION: Skyflats have different no \
            of extensions (" // nsci // ")", l_logfile, verbose+) 
        printlog ("                        to flexflat (" \
            // gemextn.count // ")", l_logfile, verbose+) 
        goto clean
    }
    gemextn (l_flexflat, check = "exists,mef", process = "expand", \
        index = "0", extname = "", extversion = "", \
        ikparams = "", omit = "", replace = "", outfile = "STDOUT", \
        logfile="dev$null", glogpars="", verbose-) | scan (line)
    shx = INDEF; shy = INDEF
    hselect (line, "MDF_XSHF", yes) | scan (shx)
    hselect (line, "MDF_YSHF", yes) | scan (shy)
    if (! havedark) {
        flx_shx = shx
        flx_shy = shy
    } else {
        ok = isindef (flx_shx) && isindef (shx)
        if (no == ok && no == isindef (flx_shx) && no == isindef (shx)) {
            ok = abs (flx_shx - shx) < 0.1
        }
        if (no == ok) {
            printlog ("ERROR - NSSLITFUNCTION: X MDF shift has changed \
                for " // l_dark // "/" // l_flexflat // ".", \
                l_logfile, verbose+)
            goto clean
        }
        ok = isindef (flx_shy) && isindef (shy)
        if (no == ok && no == isindef (flx_shy) && no == isindef (shy)) {
            ok = abs (flx_shy - shy) < 0.1
        }
        if (no == ok) {
            printlog ("ERROR - NSSLITFUNCTION: Y MDF shift has changed \
                for " // l_dark // "/" // l_flexflat // ".", \
                l_logfile, verbose+)
            goto clean
        }
    }


    # Check output not already present

    gemextn (l_output, check = "absent", process = "none", \
        index = "", extname = "", extversion = "", ikparams = "", \
        omit = "", replace = "", outfile = "dev$null", logfile="",
        glogpars="", verbose=l_verbose)
    if (gemextn.fail_count != 0 || gemextn.count != 1) {
        printlog ("ERROR - NSSLITFUNCTION: Bad or existing output.", \
            l_logfile, verbose+) 
        goto clean
    }
    gemextn (l_output, check = "absent", process = "none", \
        index = "", extname = "", extversion = "", ikparams = "", \
        omit = "", replace = "", outfile = "STDOUT", \
        logfile="dev$null", glogpars="", verbose-) | scan (line)
    l_output = line


    # Check each skyflat for various preconditions

    scanfile = phulist
    while (fscan (scanfile, phu) != EOF) {

        keyfound = ""
        hselect(phu, "*PREPAR*", yes) | scan (keyfound)
        if (keyfound == "") {
            printlog ("ERROR - NSSLITFUNCTION: Image " // phu \
                // " not PREPAREd.", l_logfile, verbose+) 
            goto clean
        }

        shx = INDEF; shy = INDEF
        hselect (phu, "MDF_XSHF", yes) | scan (shx)
        hselect (phu, "MDF_YSHF", yes) | scan (shy)
        ok = isindef (flx_shx) && isindef (shx)
        if (no == ok && no == isindef (flx_shx) && no == isindef (shx)) {
            ok = abs (flx_shx - shx) < 0.1
        }
        if (no == ok) {
            printlog ("ERROR - NSSLITFUNCTION: X MDF shift has changed \
                for " // phu // "/" // l_flexflat // ".", l_logfile, \
                verbose+)
            goto clean
        }
        ok = isindef (flx_shy) && isindef (shy)
        if (no == ok && no == isindef (flx_shy) && no == isindef (shy)) {
            ok = abs (flx_shy - shy) < 0.1
        }
        if (no == ok) {
            printlog ("ERROR - NSSLITFUNCTION: Y MDF shift has changed \
                for " // phu // "/" // l_flexflat // ".", l_logfile, \
                verbose+)
            goto clean
        }

        exptime = INDEF
        hselect (phu, l_key_exptime, yes) | scan (exptime)
        if (isindef (exptime)) {
            printlog ("WARNING - NSSLITFUNCTION: missing exposure time \
                in " // phu, l_logfile, verbose+) 
        } else if (no == isindef (prevexp)) {
            if (abs (exptime - prevexp) > 0.1 \
                && no == ("exposure" == strlwr (l_scale))) {

                printlog ("WARNING - NSSLITFUNCTION: exposure time \
                    changed", l_logfile, verbose+)
                printf ("                  %-28s: %6.2f\n", \
                    previmg, prevexp) | scan (line)
                printlog (line,  l_logfile, verbose+)
                printf ("                  %-28s: %6.2f\n", \
                    phu, exptime) | scan (line)
                printlog (line,  l_logfile, verbose+)
            }
        }
        previmg = phu; prevexp = exptime

    }


    # Combine sky flats

    if (n == 1) {

        gemextn ("@" // filelist, proc="none", check="exists", \
            index="", extname="", extversion="", ikparams="", omit="", \
            replace="", outfile="STDOUT", logfile="dev$null", \
            glogpars="", verbose=l_verbose) | scan (comboutptr)

        if (havedqvar) {
            gemextn (comboutptr, proc="expand", check="", index="", \
                extname=l_sci_ext, extversion="-", ikparams="", omit="", \
                replace="", outfile="dev$null", glogpars="", verbose-)
            if (0 == gemextn.count || 0 < gemextn.fail_count) {
                printlog ("ERROR - NSSLITFUNCTION: missing science data \
                    in " // comboutptr, l_logfile, verbose+) 
                goto clean
            }
            nextn = gemextn.count
            gemextn (comboutptr, proc="expand", check="", index="", \
                extname=l_var_ext, extversion="-", ikparams="", omit="", \
                replace="", outfile="dev$null", glogpars="", verbose-)
            if (no == (nextn == gemextn.count)) havedqvar = no
            gemextn (comboutptr, proc="expand", check="", index="", \
                extname=l_dq_ext, extversion="-", ikparams="", omit="", \
                replace="", outfile="dev$null", glogpars="", verbose-)
            if (no == (nextn == gemextn.count)) havedqvar = no
            if (no == havedqvar) {
                printlog ("WARNING - NSSLITFUNCTION: no dq/variance data \
                    in " // comboutptr, l_logfile, verbose+) 
            }
        }

    } else {

        printlog ("NSSLITFUNCTION: Combining sky flats.", \
            l_logfile, l_verbose)

        gemcombine ("@" // filelist, combout, title="", \
            combine=l_combine, reject=l_reject, \
            masktype=l_masktype, maskvalue=l_maskvalue, \
            offsets="none", scale=l_scale, zero=l_zero, \
            weight=l_weight, statsec=l_statsec, expname=l_key_exptime, \
            lthreshold=l_lthreshold, hthreshold=l_hthreshold, \
            nlow=l_nlow, \
            nhigh=l_nhigh, nkeep=l_nkeep, mclip=l_mclip, lsigma=l_lsigma, \
            hsigma=l_hsigma, key_ron=l_key_ron, key_gain=l_key_gain, \
            ron=l_ron, gain=l_gain, snoise=l_snoise, sigscale=l_sigscale, \
            pclip=l_pclip, grow=l_grow, bpmfile="", nrejfile="", \
            sci_ext=l_sci_ext, var_ext=l_var_ext, dq_ext="JUNK", \
            fl_vardq+, logfile=l_logfile, fl_dqprop=l_fl_dqprop, \
            verbose-)

        if (gemcombine.status != 0) {
            printlog ("ERROR - NSSLITFUNCTION: Error combining sky \
                flats.", l_logfile, verbose+)
            goto clean
        }

        comboutptr = combout

    }


    printlog ("NSSLITFUNCTION: Subtracting dark and applying \
        (flex)flat.", l_logfile, l_verbose)

    # Note that we don't need to carry dq/var info through here
    # We are going to fit to this data and use the resulting
    # model.  Not only do we have no way of calculating the errors
    # in that model, they are (1) horribly correlated and (2)
    # hopefully magnitudes smaller than the noise in the flat,
    # which comes directly from raw data

    if (l_fl_gemarith) {

        # Subtract dark

        if (havedark) {
            if (debug) print ("subtracting dark")
            gemarith (comboutptr, "-", l_dark, flatsub, \
                sci_ext=l_sci_ext, var_ext=l_var_ext, dq_ext=l_dq_ext, \
                fl_vardq-, logfile=l_logfile, verbose-)
            if (! imaccess (flatsub)) {
                printlog ("ERROR - NSSLITFUNCTION: Error subtracting \
                    dark.", l_logfile, verbose+)
                goto clean
            }
        } else {
            flatsub = comboutptr
        }

        # Divide by flex flat

        if (debug) print ("dividing by flex flat")
        gemarith (flatsub, "/", l_flexflat, sflat, sci_ext=l_sci_ext, \
            var_ext=l_var_ext, dq_ext=l_dq_ext, fl_vardq-, \
            logfile=l_logfile, verbose-)
        if (! imaccess (sflat)) {
            printlog ("ERROR - NSSLITFUNCTION: Error dividing by \
                (flex)flat.", l_logfile, verbose+)
            goto clean
        }

    } else {

        # Use imarith instead

        imcopy (comboutptr // "[0]", sflat // ".fits[0,overwrite]",
            >& "dev$null")
        if (debug) gemextn (sflat)

        delete (tmplis, verify-, >& "dev$null")
        delete (tmplis1, verify-, >& "dev$null")
        delete (tmplis2, verify-, >& "dev$null")
        delete (tmplis3, verify-, >& "dev$null")
        tmplis = mktemp ("tmplis1")
        tmplis1 = mktemp ("tmplis1")
        tmplis2 = mktemp ("tmplis2")
        tmplis3 = mktemp ("tmplis3")

        gemextn (comboutptr, check="exists,mef", process="expand", \
            index="", extname=l_sci_ext, extver="1-", ikparams="", \
            omit="", replace="", outfile=tmplis1, logfile=l_logfile, \
            glogpars="", verbose-)
        if (0 != gemextn.fail_count) {
            printlog ("ERROR - NSSLITFUNCTION: Error expanding " \
                // combout, l_logfile, verbose+)
            goto clean
        }
        if (havedark) {
            gemextn (l_dark, check="exists,mef", process="expand", \
                index="", extname=l_sci_ext, extver="1-", ikparams="", \
                omit="", replace="", outfile=tmplis2, \
                logfile=l_logfile, glogpars="", verbose-)
            if (0 != gemextn.fail_count) {
                printlog ("ERROR - NSSLITFUNCTION: Error expanding " \
                    // l_dark, l_logfile, verbose+)
                goto clean
            }
        } else {
            touch (tmplis2)
        }
        gemextn (l_flexflat, check="exists,mef", process="expand", \
            index="", extname=l_sci_ext, extver="1-", ikparams="", \
            omit="", replace="", outfile=tmplis3, logfile=l_logfile, \
            glogpars="", verbose-)
        if (0 != gemextn.fail_count) {
            printlog ("ERROR - NSSLITFUNCTION: Error expanding " \
                // l_flexflat, l_logfile, verbose+)
            goto clean
        }

        if (debug) print ("list1:")
        if (debug) type (tmplis1)
        if (debug && havedark) print ("list2:")
        if (debug && havedark) type (tmplis2)
        if (debug) print ("list3:")
        if (debug) type (tmplis3)

        joinlines (tmplis1 // "," // tmplis2 // "," // tmplis3, "", \
            output=tmplis, delim=" ", missing="none", maxchars=1000, \
            shortest-, verbose-)
        if (debug) print ("list:")
        if (debug) type (tmplis)

        scanfile = tmplis
        iver = 0
        while (EOF != fscan (scanfile, skyextn, darkextn, flatextn)) {

            iver = iver + 1

            if (debug)
                print ("(" // skyextn // " - " // darkextn // \
                    ") / " // flatextn)
            if (debug)
                print (sflat // "[" // l_sci_ext // "," // iver // ",append]")
            if (debug)
                gemextn (sflat // "[" // l_sci_ext // "," // iver // ",\
                    append]")

            if (havedark) expr = "(a-b)/c"
            else          expr = "a/c"

            imexpr (expr,
                sflat // "[" // l_sci_ext // "," // iver // ",append]", \
                skyextn, darkextn, flatextn, dims="auto", \
                intype="auto", outtype="auto", refim="auto", \
                bwidth=0, btype="nearest", bpixval=0, rangecheck+, \
                verbose-, exprdb="none")

        }
    }


    # Do the illumination thing to each version

    printlog ("NSSLITFUNCTION: Calculating illumination response.", \
        l_logfile, l_verbose)

    gemextn (sflat, check = "exists,mef", process = "expand", \
        index = "", extname = l_sci_ext, extversion = "1-", \
        ikparams = "", omit = "", replace = "", outfile = scilist, \
        logfile=l_logfile, glogpars="", verbose=l_verbose)
    if (gemextn.fail_count != 0 || gemextn.count == 0) {
        printlog ("ERROR - NSSLITFUNCTION: Missing internal data.", \
            l_logfile, verbose+) 
        goto clean
    }
    nver = gemextn.count


    # Output PHU

    imcopy (l_flat // ".fits[0]", l_output // ".fits[0,overwrite]", \
        >& "dev$null")

    # Output MDF

    extn = INDEF
    gemextn (comboutptr, check="exists,mef", process="expand", \
        index="0-", extname="MDF", extver="", ikparams="", \
        omit="name", replace="", outfile="STDOUT", \
        logfile="dev$null", glogpars="", verbose-) | scan (extn)
    if (no == isindef (extn))
        fxinsert (extn, l_output // ".fits[0]", groups="", verbose-,
            >& "dev$null")
    if (debug) gemextn (l_output)


    # List of source frames

    gemextn (l_flat, check="", process="expand", index="", \
        extname=l_sci_ext, extversion="1-", ikparams="", \
        omit="", replace="", outfile=ssrclist, logfile=l_logfile, \
        glogpars="", verbose=l_verbose)
    if (gemextn.fail_count != 0 || gemextn.count != nver) {
        printlog ("ERROR - NSSLITFUNCTION: Sky flats inconsistent with \
            lamp flats (cut before processing?).", l_logfile, verbose+) 
        goto clean
    }
    gemextn (l_flat, check="", process="expand", index="", \
        extname=l_var_ext, extversion="1-", ikparams="", \
        omit="", replace="", outfile=vsrclist, logfile=l_logfile,
        glogpars="", verbose=l_verbose)
    gemextn (l_flat, check="", process="expand", index="", \
        extname=l_dq_ext, extversion="1-", ikparams="", \
        omit="", replace="", outfile=dqsrclist, logfile=l_logfile,
        glogpars="", verbose=l_verbose)


    # List of output frames

    gemextn (l_output, check="", process="append", index="", \
        extname=l_sci_ext, extversion="1-" // nver, ikparams="append", \
        omit="", replace="", outfile=sdestlist, logfile=l_logfile, \
        glogpars="", verbose=l_verbose)
    gemextn (l_output, check="", process="append", index="", \
        extname=l_var_ext, extversion="1-" // nver, ikparams="append", \
        omit="", replace="", outfile=vdestlist, logfile=l_logfile,
        glogpars="", verbose=l_verbose)
    gemextn (l_output, check="", process="append", index="", \
        extname=l_dq_ext, extversion="1-" // nver, ikparams="append", \
        omit="", replace="", outfile=dqdestlist, logfile=l_logfile,
        glogpars="", verbose=l_verbose)

    joinlines (scilist // "," // ssrclist // "," // sdestlist // "," \
        // vsrclist // "," // vdestlist // "," // dqsrclist // "," \
        // dqdestlist, "", output=comblist, delim=" ", \
        missing="none", maxchars=1000, shortest-, verbose-)


    # Calculate shifts for use below

    if (isindef (shiftx)) shiftx = 0
    if (isindef (shifty)) shifty = 0
    if (isindef (flx_shx)) flx_shx = 0
    if (isindef (flx_shy)) flx_shy = 0
    # May be a sign error here
    #        dx = -1.0 * (flx_shx - shiftx)
    #        dy = -1.0 * (flx_shy - shifty)
    dx = flx_shx - shiftx
    dy = flx_shy - shifty

    if (dx == 0 && dy == 0) {
        printlog ("NSSLITFUNCTION: No flexure shift.", \
            l_logfile, l_verbose)
    } else {
        printf ("NSSLITFUNCTION: Flexure shift X=%4.1f Y=%4.1f\n", dx, dy) \
            | scan (line)
        printlog (line, l_logfile, l_verbose)
        dx = (nint ((dx - int (dx)) * 10)) / 10.0
        dy = (nint ((dy - int (dy)) * 10)) / 10.0
        printf ("NSSLITFUNCTION: Residual shift X=%4.1f Y=%4.1f\n", dx, dy) \
            | scan (line)
        printlog (line, l_logfile, l_verbose)
    }


    # Processing and normalisation in two steps to allow for normzn by
    # middle order

    scanfile = comblist
    i = 0
    n = nver / 2 + 1

    while (fscan (scanfile,img,ssrc,sdest,vsrc,vdest,dqsrc,dqdest) != EOF) {
    i = i + 1

        if (debug) {
            print ("img:      " // img)
            print ("ssrc:     " // ssrc)
            print ("sdest:    " // sdest)
        }

        if (debug) print ("illumination")
        if (debug)
            imstat (img, fields="image,npix,mean,stddev,min,max,midpt")
        illumination (img, response, low_rej=2., high_rej=2, \
            niter=3, order=l_order, nbins=5, function=l_function, \
            interpolator="poly3", interactive=l_fl_inter)
        if (debug)
            imstat (response, fields="image,npix,mean,stddev,min,max,midpt")

        # Do shift

        if (abs (dx) > epsilon || abs (dy) > epsilon) {
            if (debug) print ("shifting")
            imdelete (tmpshift, verify-, >& "dev$null")
            tmpshift = mktemp ("tmpshift")
            #               hselect (response, "i_naxis1", yes) | scan (nx)
            #               hselect (response, "i_naxis2", yes) | scan (ny)
            imshift (response, tmpshift, dx, dy, shifts_file="", \
                interp_type="spline3", boundary_type="nearest")
            # swap file names, so delete/keep right file
            line = tmpshift; tmpshift = response; response = line
            if (debug) print (response)
            if (debug)
                imstat (response, fields="image,npix,mean,stddev,min,max,\
                    midpt")
        }

        # fails if dispaxis is set for 1-D image, but code below will
        # handle that peculiar case for any likely value of dispaxis
        dispaxis = -1; nx = -1; ny = -1
        hselect (img, "i_naxis1,i_naxis2,"//l_key_dispaxis, yes) |\
            scan (nx, ny, dispaxis)

        # if dispaxis is missing, assume GNIRS default (along columns)
        # if image extension is not 2-D, does something appropriate
        # if axis is shorter than needed, or algorithm becomes confused,
        #    do nothing (i.e., sect="")
        sect = ""
        if (dispaxis == 1 && nx/2 > 10)
            printf ("[%d:%d,*]\n", (nx/2-10), (nx/2+10)) | scan (sect)
        else if (ny/2 > 10)
            printf ("[*,%d:%d]\n", (ny/2-10), (ny/2+10)) | scan (sect)

        imstat (img//sect, fields="midpt,stddev", lower=0., upper=INDEF, \
            nclip=0, lsigma=3, usigma=3, binwidth=0.1, format-,
            cache-) | scan (mddata,sigma)
        if (debug) {
            print ("1: mddata: " // mddata)
            print ("   sigma:  " // sigma)
        }

        imstat (img//sect, fields="midpt", \
            lower=max(0.3*mddata,(mddata-3*sigma)),upper=(mddata+3*sigma), \
            nclip=0, lsigma=3, usigma=3, binwidth=0.1, format-,
            cache-) | scan (mddata)
        if (debug) {
            print ("2: mddata: " // mddata)
        }

        imstat (response, fields="midpt", lower=INDEF, upper=INDEF, \
            nclip=0, lsigma=3, usigma=3, binwidth=0.1, format-,
            cache-) | scan (mdresp)
        if (debug) print ("mddata, mdresp: " // mddata // ", " // mdresp)

        print (img // " " // ssrc // " " // sdest // " " // response \
            // " " // mddata // " " // mdresp // " " // vsrc // " " \
            // vdest // " " // dqsrc // " " // dqdest, >> comblist2)

        response = mktemp ("tmpresponse")

        if (no == fl_vary && i == n) oldmddata = mddata
    }


    # So now we do the normlzn and assemble the output

    printlog ("NSSLITFUNCTION: Normalizing response.", l_logfile, l_verbose)

    scanfile = comblist2
    vsrc = ""
    vdest = ""
    first = yes

    while (fscan (scanfile, img, ssrc, sdest, response, mddata, mdresp, \
        vsrc, vdest, dqsrc, dqdest) != EOF) {

        if (debug) {
            print ("img:      " // img)
            print ("ssrc:     " // ssrc)
            print ("sdest:    " // sdest)
            print ("response: " // response)
            print ("mddata:   " // mddata)
        }

        # if fl_vary, scaling is just 1/response to get unity
        # otherwise, that's only true for the middle extension
        # and others must reflect the change in mddata
        if (fl_vary) {
            correcn = 1.0 / mdresp
        } else {
            correcn = (1.0 / mdresp) * (mddata / oldmddata)
        }
        if (debug) print ("correction: " // correcn)
        imarith (response, "*", correcn, response, verbose=debug)

        if (debug) print ("scaling sci")
        imarith (ssrc, "*",  response, scratch, verbose-)
        imcopy (scratch, sdest, >& "dev$null")
        imdelete (scratch, verify-, >& "dev$null")
        scratch = mktemp ("tmpscratch")

        if (havedqvar) {
            if (debug) print ("scaling var")
            imarith (vsrc, "*",  response, scratch, verbose-)
            imarith (scratch, "*",  response, scratch, verbose-)
            imcopy (scratch, vdest, >& "dev$null")
            imdelete (scratch, verify-, >& "dev$null")
            scratch = mktemp ("tmpscratch")

            if (debug) print ("copying dq")
            imcopy (dqsrc, dqdest, >& "dev$null")
        }

        first = no
    }


    gemdate () 
    if (debug) print (l_output)
    gemhedit (l_output // "[0]", "NSSLIT", gemdate.outdate, 
        "UT Time stamp for NSSLITFUNCITON", delete-)
    gemhedit (l_output // "[0]", "GEM-TLM", gemdate.outdate, 
        "UT Last modification with GEMINI", delete-)
    if (l_title != "") {
        gemhedit (l_output // "[0]", "i_title", l_title, "", delete-) 
    }

    status = 0


clean:

    if (debug) print ("clean")
    if (access (comblist2)) {
        scanfile = comblist2
        while (fscan (scanfile, img, ssrc, sdest, response) != EOF) {
            imdelete (response, verify-, >& "dev$null")
        }
    }
    delete (filelist, verify-, >& "dev$null")
    delete (phulist, verify-, >& "dev$null")
    delete (scilist, verify-, >& "dev$null")
    delete (ssrclist, verify-, >& "dev$null")
    delete (sdestlist, verify-, >& "dev$null")
    delete (vsrclist, verify-, >& "dev$null")
    delete (vdestlist, verify-, >& "dev$null")
    delete (dqsrclist, verify-, >& "dev$null")
    delete (dqdestlist, verify-, >& "dev$null")
    delete (comblist, verify-, >& "dev$null")
    delete (comblist2, verify-, >& "dev$null")
    delete (tmplis, verify-, >& "dev$null")
    delete (tmplis1, verify-, >& "dev$null")
    delete (tmplis2, verify-, >& "dev$null")
    delete (tmplis3, verify-, >& "dev$null")

    imdelete (combout, verify-, >& "dev$null")
    imdelete (sflat, verify-, >& "dev$null")
    imdelete (response, verify-, >& "dev$null")
    imdelete (scratch, verify-, >& "dev$null")
    imdelete (tmpshift, verify-, >& "dev$null")

    if (havedark)
        imdelete (flatsub, verify-, >& "dev$null")

    printlog (" ", l_logfile, l_verbose) 
    printlog ("NSSLITFUNCTION done", l_logfile, l_verbose) 
    printlog ("------------------------------------------------------\
        --------------------------", l_logfile, l_verbose)
        
end

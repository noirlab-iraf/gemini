# Copyright(c) 2000-2013 Association of Universities for Research in Astronomy, Inc.

procedure nsreduce (inimages) 

# Basic reductions of GNIRS/NIRI spectroscopic images.
# Requires a flat, bad pixel mask, sky and/or dark image to subtract.  
# Sky images can be constructed from the input images at different nod
# positions. 
#
# A log of the reductions is written in the logfile
# The header of the reduced image will contain information about
# the reduction as a series of header entries.
#
# Variance handling:
#   var = (df/da)^2 * siga^2 + (df/db)^2 * sigb^2
#   f=a-k*b   var(f) = var(a) + k*k*var(b)
#   f=a/b     var(f) = [sci(a)^2/sci(b)^4]*var(b) + [1/sci(b)^2]*var(a)
#   f=(a-b)/c var(f) = [(sci(a)^2+sci(b)^2-2sci(a)sci(b))/sci(c)^4]*var(c) +
#                      (1/sci(c)^2)*[var(a)+var(b)]
#
# Data Quality Array handling:
#   DQ = dq(a) OR dq(b)
#
# Read noise header keyword in output image is the quadrature sum of the
# values in the input images and the sky or dark images.  No adjustment
# has been made for division by the flat field, assuming the read noise in
# the flat field is negligible.  The saturation header keyword in the 
# the output image is the input value minus the median value in the sky
# or dark image + the constant sky value added back on after.
#
# Version  Sept 20, 2002 JJ v1.4 release
#          Aug 19, 2003  KL IRAF2.12 - new parameters
#                              hedit: addonly
#                              imstat: nclip,lsigma,usigma,cache
#          Oct 29, 2003  KL moved from niri to gnirs package
#          Nov 19, 2003  AC major rewrite - gemextn, temp files, nssky


# Basic input/output parameters

char    inimages    {prompt = "Input image(s)"}                             # OLDP-1-primary-single-prefix=r
char    outimages   {"", prompt = "Output image(s)"}                        # OLDP-1-output
char    outprefix   {"r", prompt = "Prefix for output image(s)"}            # OLDP-4


# Cut options

bool    fl_cut      {yes, prompt = "Cut the data?"}   # OLDP-2
char    section     {"", prompt = "Alternative section or keyword (blank for MDF)"} # OLDP-2-fl_cut-optional
bool    fl_corner   {yes, prompt = "Zero corners, if specified in MDF?"}    # OLDP-2
bool    fl_process_cut  {yes, prompt = "Do processing to cut data (otherwise, use uncut)?"} # OLDP-2-fl_cut-optional
char    gradimage   {"", prompt = "Image to use for finding slit edges using the gradient method (F2 MOS)"}
char    outslitim   {"slits.fits", prompt = "Output image to be used to determine s-distortion correction (F2 MOS)"}
real    edgediff    {2.0, prompt = "Allowed difference in pixels when associating slit edges to a slit (F2 MOS)"}
char    refimage    {"", prompt = "Reference image for slit positions (F2 MOS)"}
char    database    {"", prompt = "Directory for files containing slit trace (F2 MOS)"}

# Appwave options

bool    fl_nsappwave    {yes, prompt="Call nsappawave?"}                    # OLDP-2
char    nsappwavedb {"gnirs$data/nsappwave.fits", prompt = "nsappwave calibration table"} # OLDP-2-fl_nsappwave-optional
real    crval       {INDEF, prompt = "Central wavelength"}                  # OLDP-3-fl_nsappwave-optional
real    cdelt       {INDEF, prompt = "Resolution in wavelength per pixel"}  # OLDP-3-fl_nsappwave-optional


# Dark options

bool    fl_dark     {no, prompt = "Do dark subtraction?"}                   # OLDP-2
char    darkimage   {"", prompt = "Dark current image to subtract"}         # OLDP-2-input-fl_dark-optional
bool    fl_save_dark    {no, prompt = "Save dark after cutting?"}     # OLDP-3-fl_dark-optional


# Sky options

bool    fl_sky      {yes, prompt = "Do sky subtraction?"}                   # OLDP-2
char    skyimages   {"", prompt = "Sky image(s) from other nod positions"}  # OLDP-2-input-fl_sky-optional
char    skysection  {"", prompt = "Level, sample area, or header keyword for sample area"} # OLDP-2-fl_sky-optional
char    combtype    {"median", enum = "average|median", prompt = "Type of combine operation for sky"} # OLDP-2-fl_sky-optional
char    rejtype     {"avsigclip", enum = "none|minmax|ccdclip|crreject|sigclip|avsigclip|pclip", prompt = "Type of rejection for combining sky"} # OLDP-2-fl_sky-optional
char    masktype    {"goodvalue", enum="none|goodvalue", prompt="Mask type"}    # OLDP-3
real    maskvalue   {0., prompt="Mask value"}                               # OLDP-3
char    scale       {"none", prompt = "Image scaling for combining sky (see imcombine.scale)"} # OLDP-3-fl_sky-optional
char    zero        {"median", prompt = "Image zero-point offset for combining sky (see imcombine.zero)"} # OLDP-3-fl_sky-optional
char    weight      {"none", prompt = "Image weights for combining sky (see imcombine.weight)"} # OLDP-3-fl_sky-optional
char    statsec     {"[*,*]", prompt = "Statistics section"}                # OLDP-3-fl_sky-optional
real    lthreshold  {INDEF, prompt = "Lower threshold"}                     # OLDP-3-fl_sky-optional
real    hthreshold  {INDEF, prompt = "Upper threshold"}                     # OLDP-3-fl_sky-optional
int     nlow        {1, min = 0, prompt = "minmax: Number of low pixels to reject"} # OLDP-3-fl_sky-optional
int     nhigh       {1, min = 0, prompt = "minmax: Number of high pixels to reject"} # OLDP-3-fl_sky-optional
int     nkeep       {0, prompt = "Minimum to keep or maximum to reject"}    # OLDP-3-fl_sky-optional
bool    mclip       {yes, prompt = "Use median in sigma clipping algorithms?"}  # OLDP-3-fl_sky-optional
real    lsigma      {3., min = 0., prompt = "Lower sigma clipping factor"}  # OLDP-3-fl_sky-optional
real    hsigma      {3., min = 0., prompt = "Upper sigma clipping factor"}  # OLDP-3-fl_sky-optional
char    snoise      {"0.0", prompt = "ccdclip: Sensitivity noise (electrons)"}  # OLDP-3-fl_sky-optional
real    sigscale    {0.1, min = 0., prompt = "Tolerance for sigma clipping scaling correction"} # OLDP-3-fl_sky-optional
real    pclip       {-0.5, prompt = "Percentile clipping parameter"}        # OLDP-3-fl_sky-optional
real    grow        {0.0, min = 0., prompt = "Radius (pixels) for neighbor rejection"} # OLDP-3-fl_sky-optional
real    skyrange    {INDEF, prompt = "Time window for includinging sky frame (seconds)"} # OLDP-2-fl_sky-optional
real    nodsize     {3., prompt = "Minimum separation of nod positions in arcsec"} # OLDP-2-fl_sky-optional

# Flat-field options

bool    fl_flat     {yes, prompt = "Do flat-fielding?"}                     # OLDP-2
char    flatimage   {"", prompt = "Spectral flat field image to divide"}    # OLDP-2-input-fl_flat-required
real    flatmin     {0.0, prompt = "Lower limit to flat (avoiding overflows)"}  # OLDP-3


# Output options

bool    fl_vardq    {yes, prompt = "Create variance and data quality frames?"}  # OLDP-2


# Standard log options, status, structs

char    logfile     {"", prompt = "Logfile"}                                # OLDP-1
bool    verbose     {yes, prompt = "Verbose"}                               # OLDP-2
bool    debug       {no, prompt = "Very verbose"}                           # OLDP-3
bool    force       {no, prompt = "Force use with earlier IRAF versions?"}  # OLDP-3
int     status      {0, prompt = "Exit status (0=good)"}                    # OLDP-4

struct *scanin1    {prompt = "Internal use"}                               # OLDP-4
struct *scanin2    {prompt = "Internal use"}                               # OLDP-4

begin

    char    l_inimages = ""
    char    l_outimages = ""
    char    l_outprefix = ""

    bool    l_fl_cut
    char    l_section = ""
    bool    l_fl_corner
    bool    l_fl_process_cut
    char    l_gradimage = ""
    char    l_outslitim = ""
    real    l_edgediff
    char    l_refimage = ""
    char    l_database = ""

    bool    l_fl_nsappwave
    char    l_nsappwavedb = ""
    real    l_crval, l_cdelt

    bool    l_fl_dark
    char    l_darkimage = ""
    bool    l_fl_save_dark

    bool    l_fl_sky
    char    l_skyimages = ""
    char    l_skysection = ""
    char    l_combtype = ""
    char    l_rejtype = ""
    char    l_masktype = ""
    real    l_maskvalue
    char    l_scale = ""
    char    l_zero = ""
    char    l_weight = ""
    char    l_statsec = ""
    real    l_lthreshold, l_hthreshold
    int     l_nlow, l_nhigh, l_nkeep
    bool    l_mclip
    real    l_lsigma, l_hsigma
    char    l_snoise = ""
    real    l_sigscale, l_pclip, l_grow, l_skyrange, l_nodsize

    bool    l_fl_flat
    char    l_flatimage = ""
    real    l_flatmin

    bool    l_fl_vardq
    char    l_sci_ext = ""
    char    l_var_ext = ""
    char    l_dq_ext = ""

    char    l_key_ron = ""
    char    l_key_gain = ""
    char    l_key_date = ""
    char    l_key_time = ""
    char    l_key_camera = ""
    char    l_key_grating = ""
    char    l_key_filter = ""
    char    l_key_decker = ""
    char    l_key_fpmask = ""
    char    l_key_order = ""
    char    l_key_dispaxis = ""
    char    l_key_wave = ""
    char    l_key_waveorder = ""
    char    l_key_sat = ""
    char    l_key_nonlinear = ""
    char    l_key_instrument = ""

    char    l_logfile = ""
    bool    l_verbose
    bool    l_debug
    bool    l_force


    int     junk, count, istart, iend, version, nfiles, nsky, dqcount
    int     imdf, nxt, dispaxis, daxis
    real    shiftx, shifty, shx, shy, distance
    char    tmpstageC, tmpcutin, tmpcutout, tmpstageS, tmpstageD
    char    tmpstageE, tmpstageO, tmpstageF, tmpdq, tmpnsappwavein
    char    tmpnsappwaveout, tmpstageDimg, unused, diststring
    char    tmpdqlist, tmpskyin, tmpskyout, tmpexpand, tmpskyinsky
    char    tmpforsub, tmpsort, tmpdiff, badhdr, mdfdark, tmpdark
    char    inimg, outimg, root, skyimg, line, tmpimg, phu, result
    char    darkimg, flatimg, dqimg, ignore, skylist
    char    flatfilter, flatslit, flatdecker, skyslit, skyfilter, skydecker
    char    expression, subexpr, varexpression, subvarexpr, keyfound
    char    filter, slit, secn, decker
    char    varlist, scilist, dqin, dqexpression
    char    ima, imb, imc, imd, ime, imf, img, imh
    bool    donecut, intdbg, dovardq, dowarn, dosky, dodark
    bool    doscale, doflat, first, first2, dolog, dupsky, ok
    bool    call_cut, call_nsappwave, cut_after_proc
    bool    scaleisvalue, scaleissection, f2
    struct  sdate, sline
    real    scalevalue, temp, stddev, sat

    intdbg = no
    status = 1
    skyfilter = "none"
    skyslit = "none"
    skydecker = "none"
    unused = "(unused)"
    f2 = no
    tmpstageD = mktemp ("tmpstageD") 
    tmpstageDimg = mktemp ("tmpstageDimg") 
    print (tmpstageDimg, >> tmpstageD)
    tmpstageS = mktemp ("tmpstageS") 
    print (tmpstageS, >> tmpstageD)
    tmpstageC = mktemp ("tmpstageC") 
    print (tmpstageC, >> tmpstageD)
    tmpstageE = mktemp ("tmpstageE") 
    print (tmpstageE, >> tmpstageD)
    tmpstageO = mktemp ("tmpstageO") 
    print (tmpstageO, >> tmpstageD)
    tmpstageF = mktemp ("tmpstageF") 
    print (tmpstageF, >> tmpstageD)
    tmpcutin = mktemp ("tmpcutin") 
    print (tmpcutin, >> tmpstageD)
    tmpcutout = mktemp ("tmpcutout") 
    print (tmpcutout, >> tmpstageD)
    tmpnsappwavein = mktemp ("tmpnsappwavein") 
    print (tmpnsappwavein, >> tmpstageD)
    tmpnsappwaveout = mktemp ("tmpnsappwaveout") 
    print (tmpnsappwaveout, >> tmpstageD)
    tmpexpand = mktemp ("tmpexpand") 
    print (tmpexpand, >> tmpstageD)
    tmpsort = mktemp ("tmpsort") 
    print (tmpsort, >> tmpstageD)
    tmpdiff = mktemp ("tmpdiff") 
    print (tmpdiff, >> tmpstageD)
    tmpdq = mktemp ("tmpdq") 
    print (tmpdq, >> tmpstageDimg)
    tmpdqlist = mktemp ("tmpdqlist") 
    print (tmpdqlist, >> tmpstageD)
    tmpskyin = mktemp ("tmpskyin") 
    print (tmpskyin, >> tmpstageD)
    tmpskyout = mktemp ("tmpskyout") 
    print (tmpskyout, >> tmpstageD)
    tmpskyinsky = mktemp ("tmpskyinsky") 
    print (tmpskyinsky, >> tmpstageD)
    tmpforsub = mktemp ("tmpforsub") 
    print (tmpforsub, >> tmpstageD)
    mdfdark = mktemp ("mdfdark")
    print (mdfdark, >> tmpstageDimg)
    tmpdark = mktemp ("tmpdark")
    print (tmpdark, >> tmpstageD)

    junk = fscan (inimages, l_inimages)
    junk = fscan (outimages, l_outimages)
    junk = fscan (outprefix, l_outprefix)

    l_fl_cut = fl_cut
    junk = fscan (section, l_section)
    l_fl_corner = fl_corner
    l_fl_process_cut = fl_process_cut
    junk = fscan (gradimage, l_gradimage)
    junk = fscan (outslitim, l_outslitim)
    l_edgediff = edgediff
    junk = fscan (refimage, l_refimage)
    junk = fscan (database, l_database)

    l_fl_nsappwave = fl_nsappwave
    junk = fscan (nsappwavedb, l_nsappwavedb)
    l_crval = crval
    l_cdelt = cdelt

    l_fl_dark = fl_dark
    junk = fscan (darkimage, l_darkimage)
    l_fl_save_dark = fl_save_dark

    l_fl_sky = fl_sky
    junk = fscan (skyimages, l_skyimages)
    junk = fscan (skysection, l_skysection)
    junk = fscan (combtype, l_combtype)
    junk = fscan (rejtype, l_rejtype)
    junk = fscan (masktype, l_masktype)
    l_maskvalue = maskvalue
    junk = fscan (scale, l_scale)
    junk = fscan (zero, l_zero)
    junk = fscan (weight, l_weight)
    junk = fscan (statsec, l_statsec)
    l_lthreshold = lthreshold
    l_hthreshold = hthreshold
    l_nlow = nlow
    l_nhigh = nhigh
    l_nkeep = nkeep
    l_mclip = mclip
    l_lsigma = lsigma
    l_hsigma = hsigma
    junk = fscan (snoise, l_snoise)
    l_sigscale = sigscale
    l_pclip = pclip
    l_grow = grow
    l_skyrange = skyrange
    l_nodsize = nodsize

    l_fl_flat = fl_flat
    junk = fscan (flatimage, l_flatimage)
    l_flatmin = flatmin

    l_fl_vardq = fl_vardq

    badhdr = ""
    junk = fscan (nsheaders.sci_ext, l_sci_ext)
    if ("" == l_sci_ext) badhdr = badhdr + " sci_ext"
    junk = fscan (nsheaders.var_ext, l_var_ext)
    if ("" == l_var_ext) badhdr = badhdr + " var_ext"
    junk = fscan (nsheaders.dq_ext, l_dq_ext)
    if ("" == l_dq_ext) badhdr = badhdr + " dq_ext"
    junk = fscan (nsheaders.key_camera, l_key_camera)
    if ("" == l_key_camera) badhdr = badhdr + " key_camera"
    junk = fscan (nsheaders.key_grating, l_key_grating)
    if ("" == l_key_grating) badhdr = badhdr + " key_grating"
    junk = fscan (nsheaders.key_filter, l_key_filter)
    if ("" == l_key_filter) badhdr = badhdr + " key_filter"
    junk = fscan (nsheaders.key_decker, l_key_decker)
    if ("" == l_key_decker) badhdr = badhdr + " key_decker"
    junk = fscan (nsheaders.key_fpmask, l_key_fpmask)
    if ("" == l_key_fpmask) badhdr = badhdr + " key_fpmask"
    junk = fscan (nsheaders.key_order, l_key_order)
    if ("" == l_key_order) badhdr = badhdr + " key_order"
    junk = fscan (nsheaders.key_dispaxis, l_key_dispaxis)
    if ("" == l_key_dispaxis) badhdr = badhdr + " key_dispaxis"
    junk = fscan (nsheaders.key_wave, l_key_wave)
    if ("" == l_key_wave) badhdr = badhdr + " key_wave"
    junk = fscan (nsheaders.key_waveorder, l_key_waveorder)
    if ("" == l_key_waveorder) badhdr = badhdr + " key_waveorder"
    junk = fscan (nsheaders.key_sat, l_key_sat)
    if ("" == l_key_sat) badhdr = badhdr + " key_sat"
    junk = fscan (nsheaders.key_nonlinear, l_key_nonlinear)
    if ("" == l_key_nonlinear) badhdr = badhdr + " key_nonlinear"
    junk = fscan (nsheaders.key_ron, l_key_ron)
    if ("" == l_key_ron) badhdr = badhdr + " key_ron"
    junk = fscan (nsheaders.key_gain, l_key_gain)
    if ("" == l_key_gain) badhdr = badhdr + " key_gain"
    junk = fscan (nsheaders.key_date, l_key_date)
    if ("" == l_key_date) badhdr = badhdr + " key_date"
    junk = fscan (nsheaders.key_time, l_key_time)
    if ("" == l_key_time) badhdr = badhdr + " key_time"
    junk = fscan (nsheaders.key_instrument, l_key_instrument)
    if ("" == l_key_instrument) badhdr = badhdr + " key_instrument"

    junk = fscan (logfile, l_logfile)
    l_debug = debug
    l_verbose = verbose || l_debug
    l_force = force

    cache ("niri", "gemextn", "keypar", "nscut", "nsappwave", "nssky")
    cache ("gemdate") 

    if ("" == l_logfile) {
        junk = fscan (gnirs.logfile, l_logfile)
        if ("" == l_logfile) {
            l_logfile = "gnirs.log"
            printlog ("WARNING - NSREDUCE: Both nsreduce.logfile and \
                gnirs.logfile are", l_logfile, verbose+) 
            printlog ("                     undefined. Using " \
                // l_logfile, l_logfile, verbose+)
        }
    }

    date | scan (sdate) 
    printlog ("----------------------------------------------------------\
        ---------------------", l_logfile, l_verbose) 
    printlog ("NSREDUCE -- " // sdate, l_logfile, l_verbose) 
    printlog (" ", l_logfile, l_verbose) 

    if ("" != badhdr) {
        printlog ("ERROR - NSREDUCE: Parameter(s) missing from \
            nsheaders: " // badhdr, l_logfile, verbose+) 
        goto clean
    }

    # Convert nodsize/2. to distance. Check for nodsize=INDEF as
    # INDEF/2. leads to an UNDEF value => causes a whole bunch of 
    # problems!
    if (l_nodsize == INDEF)
        distance = INDEF
    else
        distance = l_nodsize/2.

    # Basic checks before starting processing

    if (intdbg) print ("basic checks")

    # Do a little parsing of skysection before the main loop
    doscale = no
    scaleisvalue = no
    scaleissection = no
    if (l_fl_sky) {
        if ("" != l_skysection) {
            doscale = yes
            junk = fscan (l_skysection, scalevalue)
            scaleisvalue = (0 != junk)
            if (scaleisvalue) {
                l_skysection = "value"
            } else {
                scaleissection = ("[" == substr (l_skysection, 1, 1))
            }
        }
    }

    if (l_fl_sky && l_fl_dark && ! doscale) { 
        printlog ("WARNING - NSREDUCE: Both a sky frame and a dark \
            frame will be subtracted.", l_logfile, verbose+) 
        printlog ("                    Since most sky frames contain \
            dark current, this is ", l_logfile, verbose+) 
        printlog ("                    probably a bad idea! Continuing \
            anyway.", l_logfile, verbose+) 
    }

    if ("" == l_sci_ext) {
        printlog ("ERROR - NSREDUCE: Extension sci_ext is undefined.", \
            l_logfile, verbose+) 
        goto clean
    }
    if (l_fl_vardq && "" == l_var_ext) {
        printlog ("ERROR - NSREDUCE: Extension var_ext is undefined.", \
            l_logfile, verbose+) 
        goto clean
    }
    if (l_fl_vardq && "" == l_dq_ext) {
        printlog ("ERROR - NSREDUCE: Extension dq_ext is undefined.", \
            l_logfile, verbose+) 
        goto clean
    }

    # Stage I - check input/output

    # Generate list of input files - must exist and be mef
    # Strip off extension info and kernel data - will expand 
    # extensions later

    if (intdbg) print ("inimages check")

    gemextn (l_inimages, check="exists,mef", process="none", index="", \
        extname="", extversion="", ikparams="", omit="extension,kernel", \
        replace="", outfile=tmpexpand, logfile="", glogpars="", \
        verbose=l_verbose)
    if (0 < gemextn.fail_count) {
        printlog ("ERROR - NSREDUCE: Bad or missing file in inimages.", \
            l_logfile, verbose+) 
        goto clean
    }

    sort (tmpexpand, column = 0, ignore_white = no, numeric_sort = no, \
        reverse_sort = no, > tmpsort)
    unique (tmpsort, > tmpstageC)
    diff (tmpsort, tmpstageC, > tmpdiff)
    first = yes
    scanin1 = tmpdiff
    while (fscan (scanin1, line) != EOF) {
        if (first) first = no
        else {
            printlog ("WARNING - NSREDUCE: Ignoring duplicate input " \
                // substr (line, 3, strlen (line)), l_logfile, l_verbose)
        }
    }

    gemextn ("@" // tmpstageC, check="", process="none", index="", \
        extname="", extversion="", ikparams="", omit="", \
        replace="", outfile="dev$null", logfile="", glogpars="",
        verbose=l_verbose)
    nfiles = gemextn.count

    if (0 == nfiles) {
        printlog ("ERROR - NSREDUCE: No input files.", \
            l_logfile, verbose+) 
        goto clean
    }

    count = 0
    dispaxis = INDEF
    first = yes # first with dispaxis, not first file
    shiftx = INDEF
    shifty = INDEF
    scanin1 = tmpstageC
    while (fscan (scanin1, inimg) != EOF) {
        phu = inimg // "[0]"
        keyfound=""
        hselect(phu, "*PREPAR*", yes) | scan (keyfound)
        if (keyfound == "") {
            printlog ("WARNING - NSREDUCE: File not prepared " \
                // inimg, l_logfile, l_verbose)
            count = count + 1
        }
        if (isindef(dispaxis)) {
            hselect (phu, l_key_dispaxis, yes) | scan (dispaxis)
            if (isindef (dispaxis)) {
                printlog ("WARNING - NSREDUCE: Dispaxis (" \
                    // l_key_dispaxis // ") undefined in " // inimg \
                    // ".", l_logfile, verbose+)
            }
        } else {
            daxis = INDEF
            hselect (phu, l_key_dispaxis, yes) | scan (daxis)
            if (isindef (daxis)) {
                printlog ("WARNING - NSREDUCE: Dispaxis (" \
                    // l_key_dispaxis // ") undefined in " // inimg \
                    // ".", l_logfile, verbose+)
            } else if (daxis != dispaxis) {
                printlog ("ERROR - NSREDUCE: Dispaxis has changed for " \
                    // inimg // ".", l_logfile, verbose+)
                goto clean
            }
        }
        if (no ==isindef (dispaxis)) {
            if (first) {
                hselect (phu, "MDF_XSHF", yes) | scan (shiftx)
                hselect (phu, "MDF_YSHF", yes) | scan (shifty)
                first = no
            } else {
                shx = INDEF; shy = INDEF
                hselect (phu, "MDF_XSHF", yes) | scan (shx)
                hselect (phu, "MDF_YSHF", yes) | scan (shy)
                ok = isindef (shiftx) && isindef (shx)
                if ((no == ok) && (no == isindef (shiftx)) && \
                    (no == isindef (shx))) {
                    ok = abs (shiftx - shx) < 0.1
                }
                if (no == ok) {
                    printlog ("ERROR - NSREDUCE: X MDF shift has changed \
                        for " // inimg // ".", l_logfile, verbose+)
                    goto clean
                }
                ok = isindef (shifty) && isindef (shy)
                if ((no == ok) && (no == isindef (shifty)) && \
                    (no == isindef (shy))) {
                    ok = abs (shifty - shy) < 0.1
                }
                if (no == ok) {
                    printlog ("ERROR - NSREDUCE: Y MDF shift has changed \
                        for " // inimg // ".", l_logfile, verbose+)
                    goto clean
                }
            }
        }
    }
    if (count > 0) {
        printlog ("ERROR - NSREDUCE: " // count // " input files not \
            prepared ", l_logfile, verbose+)
        goto clean
    }
    
    # Check output names

    if (intdbg) print ("checking output")
    gemextn (l_outimages, proc="none", check="absent", index="", \
        extname="", extver="", ikparams="", omit="", replace="", \
        outfile="dev$null", logfile="", glogpars="", verbose-)
    if (0 < gemextn.fail_count) {
        printlog ("ERROR - NSREDUCE: Output data already exist.", \
            l_logfile, verbose+) 
        goto clean
    }
    if (0 != gemextn.count) {
        if (nfiles != gemextn.count) {
            printlog ("ERROR - NSREDUCE: Incorrect number of output \
                files.", l_logfile, verbose+) 
            goto clean
        }
    } else {
        gemextn ("%^%" // l_outprefix // "%" // "@" // tmpexpand,
            proc="none", check="absent", index="", extname="", \
            extver="", ikparams="", omit="kernel", replace="", \
            outfile="dev$null", logfile="", glogpars="", verbose-)
        if (0 < gemextn.fail_count || nfiles != gemextn.count) {
            printlog ("ERROR - NSREDUCE: Bad output files.", \
                l_logfile, verbose+) 
            goto clean
        }
    }

    # Stage C - handle cutting, if required

    if (intdbg)
        print ("cut: " // l_fl_cut // "; process: " // l_fl_process_cut)

    scanin1 = tmpstageC
    call_cut = no
    call_nsappwave = no
    cut_after_proc = no
    first = yes
    while (fscan (scanin1, inimg) != EOF) {
        phu = inimg // "[0]"

        keypar (phu, l_key_instrument, silent+)
        if (keypar.found)
             if (keypar.value == "F2" || keypar.value == "Flam")
                 f2 = yes

        keypar (phu, "NSCUT", silent+)
        donecut = keypar.found
            
        if (intdbg) print (phu // " donecut: " // donecut)
        if (no == donecut) {
            if (l_fl_process_cut) {
                if (l_fl_cut) {
                    if (intdbg) print ("need to cut " // inimg)
                    call_cut = yes
                    print (inimg, >> tmpcutin)
                    outimg = mktemp (inimg // "-")
                    print (outimg, >> tmpcutout)
                    print (outimg, >> tmpstageDimg)
                    if (l_fl_nsappwave) {
                        if (intdbg) print ("need to call nsappwave for " \
                            // inimg)
                        call_nsappwave = yes
                        outimg = mktemp (inimg // "-")
                        print (outimg, >> tmpnsappwaveout)
                        print (outimg, >> tmpstageDimg)
                    }
                } else {
                    printlog ("ERROR - NSREDUCE: Input file "//inimg//" has \
                        not been previously cut and", l_logfile, verbose+)
                    printlog ("                  will not be cut (since \
                        fl_cut=no), but fl_process_cut=yes", l_logfile, \
                        verbose+)
                    printlog ("                  (did you mean to run \
                        nsreduce with fl_process_cut=no?).", l_logfile, \
                        verbose+)
                    goto clean
                }
            } else {
                # fl_process_cut=no
                if (l_fl_cut) {
                    # Process data without cutting, then cut afterwards
                    if (first)
                        printlog ("NSREDUCE: Processing data before cutting",
                            l_logfile, l_verbose)
                    cut_after_proc = yes
                    outimg = inimg
                } else {
                    # No need to cut data at all
                    outimg = inimg
                }
            }
        } else {
            # Data have been cut previously
            if (l_fl_process_cut) {
                if (l_fl_cut) {
                    printlog ("WARNING - NSREDUCE: Data have already been \
                        cut, but fl_cut=yes, setting fl_cut=no", l_logfile,
                        l_verbose)
                    l_fl_cut = no
                }
                outimg = inimg
                if (l_fl_nsappwave) {
                    call_nsappwave = yes
                    outimg = mktemp (inimg // "-")
                    print (inimg, >> tmpnsappwavein)
                    print (outimg, >> tmpnsappwaveout)
                    print (outimg, >> tmpstageDimg)
                }

            } else {
                printlog ("ERROR - NSREDUCE: Input file "//inimg//" has \
                    already been cut", l_logfile, verbose+)
                printlog ("                  (did you mean to run nsreduce \
                    with fl_process_cut=yes?).", l_logfile, verbose+)
                goto clean
            }
        }
        line = inimg // " " // outimg
        if (intdbg) print ("stage C: " // line)
        print (line, >> tmpstageS)
        first = no
    }

    if (call_cut) {
        if (intdbg) print ("call cut")
        if (f2) {
            f2cut (inimages="@"//tmpcutin, outimages="@"//tmpcutout, \
                outprefix="", fl_vardq=l_fl_vardq, gradimage=l_gradimage, \
                outslitim=l_outslitim, edgediff=l_edgediff, \
                refimage=l_refimage, database=l_database, logfile=l_logfile, \
                verbose=l_verbose)
            if (f2cut.status != 0) {
                printlog ("ERROR - NSREDUCE: Error in f2cut.", l_logfile, \
                    verbose+)
                goto clean
            }
        } else {
            nscut (inimages="@"//tmpcutin, outspectra="@"//tmpcutout, \
                outprefix="", section=l_section, fl_corner=l_fl_corner, \
                logfile=l_logfile, verbose=l_verbose, debug=l_debug)
            if (nscut.status != 0) {
                printlog ("ERROR - NSREDUCE: Error in nscut.", l_logfile, \
                    verbose+)
                goto clean
            }
        }
        if (call_nsappwave)
            tmpnsappwavein = tmpcutout
    }
    if (call_nsappwave) {
        if (intdbg) print ("call nsappwave")
        nsappwave (inimages="@"//tmpnsappwavein, \
            outspectra="@"//tmpnsappwaveout, outprefix="", \
            nsappwavedb=l_nsappwavedb, fl_phu-, crval=l_crval, cdelt=l_cdelt, \
            logfile=l_logfile, verbose=l_verbose, debug=l_debug)
        if (0 != nsappwave.status) {
            printlog ("ERROR - NSREDUCE: Error in nsappwave.", l_logfile, \
                verbose+) 
            goto clean
        }
    }

    # Now repeat the above for sky images, checking and cutting,
    # if required, before calling nssky

    if (intdbg) print ("sky validation")

    if (l_fl_sky) {

        # Stage I - check input (sky)

        delete (tmpexpand, verify-, >& "dev$null")
        delete (tmpsort, verify-, >& "dev$null")
        delete (tmpdiff, verify-, >& "dev$null")
        delete (tmpstageC, verify-, >& "dev$null")

        gemextn (l_skyimages, check="exists,mef", process="none", index="", \
            extname="", extversion="", ikparams="", omit="extension,kernel", \
            replace="", outfile=tmpexpand, logfile="", glogpars="",
            verbose=l_verbose)
        if (0 < gemextn.fail_count) {
            printlog ("ERROR - NSREDUCE: Bad or missing file in skyimages.", \
                l_logfile, verbose+) 
            goto clean
        }

        if (gemextn.count == 0) {
            nsky = 0
        } else {
            # Don't ignore duplicate sky images
            copy (tmpexpand, tmpstageC)

            gemextn ("@" // tmpstageC, check="", process="none", \
                index="", extname="", extversion="", ikparams="", \
                omit="", replace="", outfile="dev$null", logfile="",
                glogpars="", verbose=l_verbose)
            nsky = gemextn.count
        }

        if (0 != nsky) {
            count = 0
            scanin1 = tmpexpand
            while (fscan (scanin1, root) != EOF) {
                phu = root // "[0]"
                keypar (phu, "PREPARE", silent+)
                if (no == keypar.found) {
                    printlog ("WARNING - NSREDUCE: File not prepared " \
                        // root, l_logfile, l_verbose)
                    count = count + 1
                }
            }
            if (count > 0) {
                printlog ("ERROR - NSREDUCE: " // count // " sky \
                    files not prepared ", l_logfile, verbose+)
                goto clean
            }

            # Stage C - handle cutting (sky), if required
            # Don't bother with nsappwave for sky files

            if (intdbg) print ("cut: " // l_fl_cut // "; process: " \
                // l_fl_process_cut)

            delete (tmpcutin, verify-, >& "dev$null")
            delete (tmpcutout, verify-, >& "dev$null")

            scanin1 = tmpstageC
            call_cut = no
            while (fscan (scanin1, inimg) != EOF) {
                phu = inimg // "[0]"
                keypar (phu, "NSCUT", silent+) 
                donecut = keypar.found
                if (intdbg) print (phu // " donecut: " // donecut)
                if (no == donecut) {
                    if (l_fl_process_cut) {
                        if (l_fl_cut) {
                            if (intdbg) print ("need to cut " // inimg)
                            call_cut = yes
                            print (inimg, >> tmpcutin)
                            outimg = mktemp (inimg // "-")
                            print (outimg, >> tmpcutout)
                            print (outimg, >> tmpstageDimg)
                        } else {
                            printlog ("ERROR - NSREDUCE: Sky file "//inimg//" \
                                has not been previously cut and", l_logfile, \
                                verbose+)
                            printlog ("                  will not be cut \
                                (since fl_cut=no), but fl_process_cut=yes", \
                                l_logfile, verbose+)
                            printlog ("                  (did you mean to run \
                                nsreduce with fl_process_cut=no?).", \
                                l_logfile, verbose+)
                            goto clean
                        }
                    } else {
                        # fl_process_cut=no, so no need to cut sky
                        outimg = inimg
                    }
                } else {
                    # Data have been cut previously
                    if (l_fl_process_cut) {
                        if (l_fl_cut) {
                            printlog ("WARNING - NSREDUCE: Sky file \
                                "//inimg//" has already been cut, but \
                                fl_cut=yes, setting fl_cut=no", \
                                l_logfile, l_verbose)
                            l_fl_cut = no
                        }
                        outimg = inimg
                    } else {
                        printlog ("ERROR - NSREDUCE: Sky file "//inimg//" has \
                            already been cut", l_logfile, verbose+)
                        printlog ("                  (did you mean to run \
                            nsreduce with fl_process_cut=yes?).", l_logfile, \
                            verbose+)
                        goto clean
                    }
                }
                print (outimg, >> tmpskyinsky)
            }

            if (call_cut) {
                if (intdbg) print ("call cut")
                if (f2) {
                    f2cut (inimages="@"//tmpcutin, outimages="@"//tmpcutout, \
                        outprefix="", fl_vardq=l_fl_vardq, \
                        gradimage=l_gradimage, outslitim=l_outslitim, \
                        edgediff=l_edgediff, refimage=l_refimage, \
                        database=l_database, logfile=l_logfile, \
                        verbose=l_verbose)
                    if (f2cut.status != 0) {
                        printlog ("ERROR - NSREDUCE: Error in f2cut.", \
                            l_logfile, verbose+)
                        goto clean
                    }
                } else {
                    nscut (inimages="@"//tmpcutin, outspectra="@"//tmpcutout, \
                        outprefix="", section=l_section, \
                        fl_corner=l_fl_corner, logfile=l_logfile, \
                        verbose=l_verbose, debug=l_debug)
                    if (0 != nscut.status) {
                        printlog ("ERROR - NSREDUCE: Error in nscut.", \
                            l_logfile, verbose+)
                        goto clean
                    }
                }
            }

        } else {
            if (nfiles > 1) {
                diststring = str(distance)
                printlog ("NSREDUCE: Generating the sky frame(s) using \
                    the other nod positions", l_logfile, l_verbose) 
                if (isindef (l_skyrange)) {
                    printlog ("          (separations greater than " \
                        // diststring // " arcsec) from neighbouring", \
                        l_logfile, l_verbose) 
                    printlog ("          exposures.", \
                        l_logfile, l_verbose)
                } else {
                    printlog ("          (separations greater than " \
                        // diststring // " arcsec) and taken within", \
                        l_logfile, l_verbose) 
                    printlog ("          " // l_skyrange // " seconds.", \
                        l_logfile, l_verbose)
                }
                printlog ("", l_logfile, l_verbose)
            } else {
                printlog ("WARNING - NSREDUCE: Only one input image and \
                    no sky images specified.", l_logfile, verbose+) 
                printlog ("                    Turning sky subtraction \
                    off.", l_logfile, verbose+) 
                l_fl_sky = no
            }
        }
    }

    # Stage S - process sky

    # Check input files not already sky-subtracted
    if (l_fl_sky) {
        count = 0
        scanin1 = tmpstageS
        while (fscan (scanin1, root, inimg) != EOF) {
            phu = root // "[0]"
            keypar (phu, "SKYIMAGE", silent+)
            if (no == keypar.found) {
                print (inimg, >> tmpskyin)
                count = count + 1
            } else {
                printlog ("WARNING - NSREDUCE: " // root \
                    // " already sky-subtracted.", l_logfile, verbose+) 
            }
        }
        l_fl_sky = (0 != count)
    }

    # Generate sky frames
    if (l_fl_sky) {

        # Avoid non-existent file error if no sky data
        if (no == access (tmpskyinsky)) print ("", > tmpskyinsky)

        nssky (inimages = "@" // tmpskyin, \
            skyimages = "@" // tmpskyinsky, \
            distance = distance, age = l_skyrange, \
            combtype = l_combtype, rejtype = l_rejtype, \
            masktype=l_masktype, maskvalue=l_maskvalue, scale = l_scale, \
            zero = l_zero, weight = l_weight, statsec = l_statsec, \
            lthreshold = l_lthreshold, \
            hthreshold = l_hthreshold, nlow = l_nlow, nhigh = l_nhigh, \
            nkeep = l_nkeep, mclip = l_mclip, lsigma = l_lsigma, \
            hsigma = l_hsigma, snoise = l_snoise, \
            sigscale = l_sigscale, pclip = l_pclip, grow = l_grow, \
            fl_vardq = l_fl_vardq, index = tmpskyout, \
            logfile = l_logfile, verbose = l_verbose, debug = l_debug, \
            force = l_force)

        if (0 != nssky.status) {
            printlog ("ERROR - NSREDUCE: Problem making skies.", \
                l_logfile, verbose+)
            goto clean 
            #printlog ("WARNING - NSREDUCE: Continuing without sky.", \
            #    l_logfile, verbose+)
            #l_fl_sky = no
        }
    }

    if (l_fl_sky) {

        scanin1 = tmpstageS
        scanin2 = tmpskyout
        tmpimg = ""
        skyimg = ""
        skylist = ""
        while (fscan (scanin1, root, inimg) != EOF) {
            if (intdbg) print ("from tmpStageS:" // root // ", " // inimg)
            if ("" == tmpimg) {
                junk = fscan (scanin2, tmpimg, skyimg, ignore, skylist)
                if (intdbg) 
                print ("from tmpskyout: " // tmpimg // ", " // skyimg \
                    // ", " // skylist)
                if (intdbg) print ("access: " // access (skyimg))
            }
            if (inimg == tmpimg) {
                line = root // " " // inimg // " " // skyimg // " " \
                    // skylist
                tmpimg = ""
            } else {
                line = root // " " // inimg // " " // "none none"
            }
            if (intdbg) print ("stage S: " // line)
            print (line, >> tmpstageE)
        }

    } else {

        scanin1 = tmpstageS
        while (fscan (scanin1, root, inimg) != EOF) {
            line = root // " " // inimg // " " // "none none"
            if (intdbg) print ("stage S: " // line)
            print (line, >> tmpstageE)
        }

    }

    # Stage E - expand extensions, check whether var and dq
    # processing required and possible

    if (intdbg) print ("expansion, vardq check")

    scanin1 = tmpstageE
    count = 0
    nfiles = 0
    while (fscan (scanin1, root, inimg, skyimg, skylist) != EOF) {
        nfiles = nfiles + 1

        delete (tmpexpand, verify-, >& "dev$null")
        gemextn (inimg, check="exists", process="expand", index="", \
            extname="SCI", extversion="1-", ikparams="", \
            omit="index,name,params,section", \
            replace="", outfile=tmpexpand, logfile="", glogpars="",
            verbose=l_verbose)

        if ("none" != skyimg) {
            junk = gemextn.count
            gemextn (inimg, check="exists", process="expand", index="", \
                extname="SCI", extversion="1-", ikparams="", \
                omit="index,name,params,section", \
                replace="", outfile="dev$null", logfile="", glogpars="",
                verbose=l_verbose)
            if (junk != gemextn.count) {
                printlog ("ERROR - NSREDUCE: Extensions for image " \
                    // inimg // " (from " // root // ") do not match \
                    sky " // skyimg, l_logfile, verbose+)
                goto clean
            }
        }

        if (intdbg) type (tmpexpand)
        scanin2 = tmpexpand
        dowarn = yes # only warn once per input root
        while (fscan (scanin2, outimg) != EOF) {

            # ugly - need a parse task or change to gemextn output
            if (intdbg) print (outimg)
            istart = stridx ("=", outimg)
            iend = stridx ("]", outimg)
            version = int (substr (outimg, istart+1, iend-1))

            dovardq = l_fl_vardq

            if (l_fl_vardq) {
                if (no ==imaccess (inimg // "[" // l_var_ext // "," \
                    // version // "]")) {
                    
                    if (dowarn) {
                        printlog ("WARNING - NSREDUCE: Cannot compute VAR \
                            planes because", l_logfile, verbose+) 
                        printlog ("                    input image " \
                            //root // " does not have a VAR plane.", \
                            l_logfile, verbose+) 
                        printlog ("                    Resetting \
                            l_vardq=no for some data in this file.", \
                            l_logfile, verbose+) 
                    }
                    dovardq = no
                } else if (no ==imaccess (inimg // "[" // l_dq_ext // "," \
                    // version // "]")) {
                    
                    if (dowarn) {
                        printlog ("WARNING - NSREDUCE: Cannot compute \
                            DQ planes because", l_logfile, verbose+) 
                        printlog ("                    input image " \
                            // root // " does not have a DQ plane.", \
                            l_logfile, verbose+) 
                        printlog ("                    Resetting \
                            fl_vardq=no for some data in this file.", \
                            l_logfile, verbose+) 
                    }
                    dovardq = no
                }
                dowarn = dowarn && dovardq
            }

            line = root // " " // inimg // " " // skyimg // " " \
                // skylist // " " // version // " " // dovardq
            if (intdbg) print ("stage E: " // line)
            print (line, >> tmpstageO)
            count = count + 1
        }
    }
    printlog ("NSREDUCE: Processing " // count // " extension(s) from " \
        // nfiles // " file(s).", l_logfile, l_verbose) 

    # TODO - value for skyfilter, skyslit (nssky below?)

    # Stage O - Check/sort output images

    if (intdbg) print ("outimages check")

    delete (tmpexpand, verify-, >& "dev$null")
    delete (tmpsort, verify-, >& "dev$null")
    gemextn (l_outimages, check="absent", process="none", index="", \
        extname="", extversion="", ikparams="", omit="extension,kernel", \
        replace="", outfile=tmpexpand, logfile="", glogpars="",
        verbose=l_verbose)
    if (0 == gemextn.count) {
        if (intdbg) print ("output from substitution")
        scanin1 = tmpstageS
        while (fscan (scanin1, root, inimg) != EOF)
            print (root, >> tmpforsub)
        delete (tmpexpand, verify-, >>& "dev$null")
        gemextn ("%^%" // l_outprefix // "%" // "@" // tmpforsub, \
            check="absent", process="none", index="", \
            extname="", extversion="", ikparams="", omit="kernel,exten", \
            replace="", outfile=tmpexpand, logfile="", glogpars="",
            verbose=l_verbose)
        if (gemextn.fail_count != 0 || gemextn.count == 0) {
            printlog ("ERROR - NSREDUCE: No or incorrectly formatted \
                output files", l_logfile, verbose+) 
            goto clean
        }
    }
    sort (tmpexpand, column = 0, ignore_white = no, numeric_sort = no, \
        reverse_sort = no, > tmpsort)

    scanin1 = tmpstageO
    scanin2 = tmpsort
    tmpimg = ""
    while (fscan (scanin1, root, inimg, skyimg, skylist, version, \
        dovardq) != EOF) {

        if (tmpimg != root) {
            junk = fscan (scanin2, outimg)
            tmpimg = root
        }
        line = root // " " // inimg // " " // skyimg // " " // skylist \
            // " " // version // " " // dovardq // " " // outimg
        if (intdbg) print ("stage O: " // line)
        print (line, >> tmpstageF)
    }

    # Stage F - Check flat and dark images

    if (intdbg) print ("flat check " // l_fl_flat)

    if (l_fl_flat) {
        delete (tmpexpand, verify-, >>& "dev$null")
        gemextn (l_flatimage, check="exists,mef", process="none", \
            index="", extname="", extversion="", ikparams="", \
            omit="extension,kernel", replace="", outfile=tmpexpand, \
            logfile="", glogpars="", verbose=l_verbose)
        if (0 < gemextn.fail_count) {
            printlog ("ERROR - NSREDUCE: Bad or missing file in \
                flatimage.", l_logfile, verbose+) 
            goto clean
        }
        if (1 != gemextn.count) {
            printlog ("ERROR - NSREDUCE: Missing or multiple \
                flatimage.", l_logfile, verbose+) 
            goto clean
        }
        scanin1 = tmpexpand
        junk = fscan (scanin1, l_flatimage)
        if (intdbg) print ("flat field: " // l_flatimage)

        phu = l_flatimage // "[0]"
        shx = INDEF; shy = INDEF
        hselect (phu, "MDF_XSHF", yes) | scan (shx)
        hselect (phu, "MDF_YSHF", yes) | scan (shy)
        ok = isindef (shiftx) && isindef (shx)
        if (no == ok && no == isindef (shiftx) && no == isindef (shx)) {
            ok = abs (shiftx - shx) < 0.1
        }
        if (no == ok) {
            printf ("ERROR - NSREDUCE: Inconsistent X MDF shift for \
                flat: %6f != %6f\n", shiftx, shx) | scan (sline)
            printlog (sline, l_logfile, verbose+)
            goto clean
        }
        ok = isindef (shifty) && isindef (shy)
        if (no == ok && no == isindef (shifty) && no == isindef (shy)) {
            ok = abs (shifty - shy) < 0.1
        }
        if (no == ok) {
            printf ("ERROR - NSREDUCE: Inconsistent Y MDF shift for \
                flat: %6f != %6f\n", shifty, shy) | scan (sline)
            printlog (sline, l_logfile, verbose+)
            goto clean
        }

        if (l_fl_vardq) {
            gemextn (l_flatimage, check="exists", process="expand", \
                index="", extname=l_var_ext, extversion="1-", \
                ikparams="", omit="", replace="", outfile="dev$null", 
                logfile="dev$null", glogpars="", verbose-)
            if (0 == gemextn.count || 0 != gemext.fail_count) {
                printlog ("WARNING - NSREDUCE: Flat field image does \
                    not contain VAR data.", l_logfile, verbose+) 
                printlog ("                    Setting fl_vardq=no \
                    and proceeding.", l_logfile, verbose+) 
                l_fl_vardq = no
            } else {
                gemextn (l_flatimage, check="exists", process="expand", \
                    index="", extname=l_dq_ext, extversion="1-", \
                    ikparams="", omit="", replace="", outfile="dev$null",
                    logfile="dev$null", glogpars="", verbose-)
                if (0 == gemextn.count || 0 != gemext.fail_count) {
                    printlog ("WARNING - NSREDUCE: Flat field image does \
                        not contain DQ data.", l_logfile, verbose+) 
                    printlog ("                    Setting fl_vardq=no \
                        and proceeding.", l_logfile, verbose+) 
                    l_fl_vardq = no
                }
            }
        }
        keypar (l_flatimage // "[0]", l_key_filter, silent+)
        if (no == keypar.found) {
            if (l_key_filter != unused) {
                printlog ("WARNING - NSREDUCE: Cannot read filter from \
                    flat field image.", l_logfile, verbose+) 
            }
            flatfilter = "none"
        } else {
            flatfilter = keypar.value
        }
        keypar (l_flatimage // "[0]", l_key_fpmask, silent+)
        if (no == keypar.found) {
            if (l_key_fpmask != unused) {
                printlog ("WARNING - NSREDUCE: Cannot read slit from flat \
                    field image.", l_logfile, verbose+) 
            }
            flatslit = "none"
        } else {
            flatslit = keypar.value
        }
        keypar (l_flatimage // "[0]", l_key_decker, silent+)
        if (no == keypar.found) {
            if (l_key_decker != unused) {
                printlog ("WARNING - NSREDUCE: Cannot read decker from \
                    flat field image.", l_logfile, verbose+) 
            }
            flatdecker = "none"
        } else {
            flatdecker = keypar.value
        }
    }

    if (intdbg) print ("dark check " // l_fl_dark)

    if (l_fl_dark) { 
        delete (tmpexpand, verify-, >>& "dev$null")
        gemextn (l_darkimage, check="exists,mef", process="none", \
            index="", extname="", extversion="", ikparams="", \
            omit="extension,kernel", replace="", outfile=tmpexpand, \
            logfile="", glogpars="", verbose=l_verbose)
        if (0 < gemextn.fail_count) {
            printlog ("ERROR - NSREDUCE: Bad or missing file in \
                dark image.", l_logfile, verbose+) 
            goto clean
        }
        if (1 != gemextn.count) {
            printlog ("ERROR - NSREDUCE: Missing or multiple \
                dark image.", l_logfile, verbose+) 
            goto clean
        }
        scanin1 = tmpexpand
        junk = fscan (scanin1, l_darkimage)
        if (intdbg) print ("dark: " // l_darkimage)

        if (l_fl_vardq) {
            gemextn (l_darkimage, check="exists", process="expand", \
                index="", extname=l_var_ext, extversion="1-", \
                ikparams="", omit="", replace="", outfile="dev$null",
                logfile="dev$null", glogpars="", verbose-)
            if (0 == gemextn.count || 0 != gemext.fail_count) {
                printlog ("WARNING - NSREDUCE: Dark image does \
                    not contain VAR data.", l_logfile, verbose+) 
                printlog ("                    Setting fl_vardq=no \
                    and proceeding.", l_logfile, verbose+) 
                l_fl_vardq = no
            } else {
                gemextn (l_darkimage, check="exists", process="expand", \
                    index="", extname=l_dq_ext, extversion="1-", \
                    ikparams="", omit="", replace="", outfile="dev$null",
                    logfile="dev$null", glogpars="", verbose-)
                if (0 == gemextn.count || 0 != gemext.fail_count) {
                    printlog ("WARNING - NSREDUCE: Dark image does \
                        not contain DQ data.", l_logfile, verbose+) 
                    printlog ("                    Setting fl_vardq=no \
                        and proceeding.", l_logfile, verbose+) 
                    l_fl_vardq = no
                }
            }
        }

        keypar (l_darkimage // "[0]", l_key_filter, silent+)
        if (no == keypar.found) {
            printlog ("WARNING - NSREDUCE: Cannot read filter from dark \
                image.", l_logfile, verbose+) 
        } else if ("blank" == keypar.value || "BLANK" == keypar.value) {
            printlog ("WARNING - NSREDUCE: Dark current image not \
                taken with filter wheel blanked.", l_logfile, verbose+) 
        }

        keypar (l_darkimage // "[0]", "NSCUT", silent+)
        if (keypar.found && ! l_fl_process_cut) {

            printlog ("ERROR - NSREDUCE: Dark image has been processed by \
                NSCUT.", l_logfile, verbose+) 
            goto clean

        } else if (no == keypar.found && l_fl_process_cut) {

            # Need to cut for these data

            gemextn (l_darkimage, proc="expand", check="exists,mef", \
                index="", extver="", extname="MDF", ikparams="", omit="", \
                replace="", outfile="dev$null", logfile="dev$null",
                glogpars="", verbose-)
            if (0 == gemextn.count) {

                # Need to copy MDF across

                scanin1 = tmpstageF
                junk = fscan (scanin1, root, inimg, skyimg, skylist, \
                    version, dovardq, outimg)

                # MDF index
                gemextn (root, proc="expand", index="0-", extname="MDF", \
                    extver="", check="exists", omit="kernel", \
                    ikparams="", replace="", outfile="STDOUT", \
                    logfile="dev$null", glogpars="", verbose-) | \
                    scan (sline)
                if (1 != gemextn.count || 0 != gemextn.fail_count) {
                    printlog ("ERROR - NSREDUCE: bad MDF in " // root, \
                        l_logfile, verbose+) 
                    goto clean
                }
                if (intdbg) print (sline)
                if (intdbg)
                    print (substr (sline, stridx ("[", sline) + 1, \
                        stridx ("]", sline) - 1))
                imdf = int (substr (sline, stridx ("[", sline) + 1, \
                    stridx ("]", sline) - 1))

                if (intdbg) print (l_darkimage // " --> " // mdfdark)
                copy (l_darkimage // ".fits", mdfdark // ".fits")
                fxinsert (input = root // ".fits[" // imdf // "]", \
                    output=mdfdark // ".fits[0]", groups="", verbose=intdbg)

            } else {
                mdfdark = l_darkimage
            }

            if (l_fl_save_dark) {
                gemextn ("c@" // l_darkimage, proc="none", check="absent", 
                    omit="", replace="", ikparams="", outfile=tmpdark, 
                    logfile="", glogpars="", verbose=l_verbose)
                if (gemextn.fail_count > 0) {
                    printlog ("ERROR - NSREDUCE: can't save dark - \
                        already exists.", l_logfile, verbose+) 
                    goto clean
                }
            } else {
                tmpimg = mktemp ("tmpdark")
                print (tmpimg, >> tmpdark)
                print (tmpimg, >> tmpstageDimg)
            }
            if (intdbg) print ("call cut")
            if (f2) {
                printlog ("ERROR - NSREDUCE: Cannot cut dark images for \
                    FLAMINGOS-2.", l_logfile, verbose+)
                goto clean
            } else {
                nscut (inimages=mdfdark, outspectra="@"//tmpdark, \
                    outprefix="", section="", fl_corner = l_fl_corner, \
                    logfile=l_logfile, verbose=l_verbose, debug=l_debug)
                if (nscut.status != 0) {
                    printlog ("ERROR - NSREDUCE: problem cutting darks.", \
                        l_logfile, verbose+)
                    goto clean
                }
            }

            type (tmpdark) | scan (l_darkimage)
        }
    }

    # Stage M - Main loop, process data

    if (intdbg) print ("stage M: process data")

    scanin1 = tmpstageF
    tmpimg = ""
    count = 0
    dolog = no
    first2 = yes
    while (fscan (scanin1, root, inimg, skyimg, skylist, version, \
        dovardq, outimg) != EOF) {
        
        ## old delete position

        # catch global vardq failure
        dovardq = dovardq && l_fl_vardq

        # For each new file

        if (tmpimg != root) {
            ## New delete position
            delete (tmpdqlist, verify-, >>& "dev$null")

            # Start a new list
            print (inimg, >> tmpdqlist)

            if (intdbg) print ("new root: " // root)

            tmpimg = root
            expression = "a"
            varexpression = "e"
            phu = inimg // "[0]"
            slit = "none"
            filter = "none"
            count = count + 1
            dolog = yes

            keypar (phu, l_key_filter, silent+)
            if (keypar.found) {
                filter = keypar.value
            } else {
                if (l_key_filter != unused) {
                    printlog ("WARNING - NSREDUCE: Cannot read filter \
                        from image.", l_logfile, verbose+) 
                }
            }

            keypar (phu, l_key_fpmask, silent+)
            if (keypar.found) {
                slit = keypar.value
            } else {
                if (l_key_fpmask != unused) {
                    printlog ("WARNING - NSREDUCE: Cannot read slit \
                        from image.", l_logfile, verbose+) 
                }
            }

            if (intdbg) print ("checking dark")

            # Check for previous dark subtraction
            dodark = l_fl_dark
            darkimg = "none"
            if (dodark) {
                keypar (phu, "DARKIMAG", silent+)
                if (keypar.found) {
                    printlog ("WARNING - NSREDUCE: Image " // root \
                        // " has already been dark-subtracted.", \
                        l_logfile, verbose+) 
                    printlog ("                    Dark-subtraction \
                        NOT performed.", l_logfile, verbose+) 
                    dodark = no
                }
            }

            # Extend expression
            if (dodark) {
                darkimg = l_darkimage
                expression = "(" // expression // "-c)"
                varexpression = "(" // varexpression // "+g)"
                print (darkimg, >> tmpdqlist)
            }

            if (intdbg) print ("checking sky: " // phu)

            # Check for previous sky subtraction
            dosky = l_fl_sky && "none" != skyimg
            if (dosky) {
                keypar (phu, "SKYIMAGE", silent+)
                if (keypar.found) {
                    printlog ("WARNING - NSREDUCE: Image " // root \
                        // " has already been sky-subtracted.", \
                        l_logfile, verbose+) 
                    printlog ("                    Sky-subtraction \
                        NOT performed.", l_logfile, verbose+) 
                    dosky = no
                    skyimg = "none"
                }
            }

            # Check filter, slit, and extend expression
            if (dosky) {

            if (intdbg) print ("checking sky: " // skyimg)

                keypar (skyimg // "[0]", l_key_filter, silent+)
                if (no == keypar.found) {
                    if (l_key_filter != unused) {
                        printlog ("WARNING - NSREDUCE: Cannot read filter \
                            from sky image.", l_logfile, verbose+) 
                    }
                    skyfilter = "none"
                } else {
                    skyfilter = keypar.value
                }

                if (filter != skyfilter) {
                    printlog ("WARNING - NSREDUCE: The grism/filter \
                        for image " // root // " does not match", 
                        l_logfile, verbose+) 
                    printlog ("                    the sky image \
                        grism/filter. Proceeding anyway.", 
                        l_logfile, verbose+) 
                }

                keypar (skyimg // "[0]", l_key_fpmask, silent+)
                if (no == keypar.found) {
                    if (l_key_fpmask != unused) {
                        printlog ("WARNING - NSREDUCE: Cannot read slit \
                            from sky image.", l_logfile, verbose+) 
                    }
                    skyslit = "none"
                } else {
                    skyslit = keypar.value
                }

                if (slit != skyslit) {
                    printlog ("WARNING - NSREDUCE: The slit for image " \
                        // root // " does not match", l_logfile, verbose+) 
                    printlog ("                    the sky image \
                        slit. Proceeding anyway.", l_logfile, verbose+) 
                }

                keypar (skyimg // "[0]", l_key_decker, silent+)
                if (no == keypar.found) {
                    if (l_key_decker != unused) {
                        printlog ("WARNING - NSREDUCE: Cannot read decker \
                            from sky image.", l_logfile, verbose+) 
                    }
                    skydecker = "none"
                } else {
                    skydecker = keypar.value
                }

                keypar (phu, l_key_decker, silent+)
                if (no == keypar.found) {
                    if (l_key_decker != unused) {
                        printlog ("WARNING - NSREDUCE: Cannot read decker \
                            from image.", l_logfile, verbose+) 
                    }
                    decker = "none"
                } else {
                    decker = keypar.value
                }
                if (decker != skydecker) {
                    printlog ("WARNING - NSREDUCE: The decker for image " \
                        // root // " does not match", l_logfile, verbose+) 
                    printlog ("                    the sky image \
                        decker. Proceeding anyway.", l_logfile, verbose+) 
                } 

                if (doscale) {
                    expression = "(" // expression // "-(b*%.5f))"
                    varexpression = "(" // varexpression // "+(f*%.5f))"
                } else {
                    expression = "(" // expression // "-b)"
                    varexpression = "(" // varexpression // "+f)"
                }
                print (skyimg, >> tmpdqlist)
            } 

            if (intdbg) print ("checking flat")

            # Check for previous flat fielding
            doflat = l_fl_flat
            flatimg = "none"
            if (doflat) {
                keypar (phu, "FLATIMAG", silent+)
                if (keypar.found) {
                    printlog ("WARNING - NSREDUCE: Image " // root \
                        // " has already been flat-fielded.", \
                        l_logfile, verbose+) 
                    printlog ("                    Flat-fielding \
                        NOT performed.", l_logfile, verbose+) 
                    doflat = no
                    flatimg = "none"
                }
            }

            # Check filter, slit, and extend expression
            if (doflat) {
                keypar (phu, l_key_filter, silent+)
                if (no == keypar.found) {
                    if (l_key_filter != unused) {
                        printlog ("WARNING - NSREDUCE: Cannot read filter \
                            from image.", l_logfile, verbose+) 
                    }
                } else {
                    filter = keypar.value
                }
                if (filter != flatfilter) {
                    printlog ("WARNING - NSREDUCE: The grism/filter \
                        for image " // root // " does not match", 
                        l_logfile, verbose+) 
                    printlog ("                    the flat image \
                        grism/filter. Proceeding anyway.", 
                        l_logfile, verbose+) 
                }
                keypar (phu, l_key_fpmask, silent+)
                if (no == keypar.found) {
                    if (l_key_fpmask != unused) {
                        printlog ("WARNING - NSREDUCE: Cannot read slit \
                            from image.", l_logfile, verbose+) 
                    }
                } else {
                    slit = keypar.value
                }
                if (slit != flatslit) {
                    printlog ("WARNING - NSREDUCE: The slit for image " \
                        // root // " does not match", l_logfile, verbose+) 
                    printlog ("                    the flat image \
                        slit. Proceeding anyway.", l_logfile, verbose+) 
                } 
                keypar (phu, l_key_decker, silent+)
                if (no == keypar.found) {
                    if (l_key_decker != unused) {
                        printlog ("WARNING - NSREDUCE: Cannot read decker \
                            from image.", l_logfile, verbose+) 
                    }
                    decker = "none"
                } else {
                    decker = keypar.value
                }
                if (decker != flatdecker) {
                    printlog ("WARNING - NSREDUCE: The decker for image " \
                        // root // " does not match", l_logfile, verbose+) 
                    printlog ("                    the flat image \
                        decker. Proceeding anyway.", l_logfile, verbose+) 
                } 

                flatimg = l_flatimage

                expression = expression // "/d"
                varexpression = "("+varexpression+"/(d*d))+"
                varexpression = varexpression // "(h/(d*d*d*d))*("
                varexpression = varexpression // expression // "*" \
                    // expression // ")"
                print (flatimg, >>& tmpdqlist)

            }

            # Saturation value
            if ((doscale && ! scaleisvalue)) {
                if (intdbg) print ("reading saturation value " // phu)
                keypar (phu, l_key_sat, silent+) 
                if (no == keypar.found) {
                    printlog ("ERROR - NSREDUCE: Could not get \
                        saturation level from header", l_logfile, verbose+)
                    printlog ("                  of image " // root, \
                        l_logfile, verbose+)
                    goto clean
                } else {
                    sat = real (keypar.value)
                }
            }

            # Create the output file

            if (intdbg) print ("copying phu " // phu // " to " // outimg)
            imcopy (phu, outimg, verbose-)

            if (no == fl_vardq) {
                gemextn (inimg, check="exists,mef", process="expand", 
                    index="", extname=l_var_ext, extversion="1-", \
                    ikparams="", omit="", replace="", outfile="dev$null", \
                    logfile="dev$null", glogpars="", verbose-)
                if (gemextn.count != 0) {
                    keypar (phu, "NEXTEND", silent+)
                    if (keypar.found) {
                        nxt = int (keypar.value)
                        if (intdbg) print ("nextend: " // nxt)
                        nxt = nxt / 3
                        gemhedit (outimg // "[0]", "NEXTEND", nxt, "", delete-)
                    }
                }
            }

            gemextn (inimg, check="exists,table", process="expand", 
                index="", extname="MDF", extversion="", \
                ikparams="", omit="", replace="", outfile="dev$null", \
                logfile="dev$null", glogpars="", verbose-)
            if (gemextn.count == 1) {
                tcopy (inimg // ".fits[MDF]", outimg // ".fits[MDF]", \
                    verbose=intdbg)
            }

        } # For new file

        # Do the math, finally!

        ima = inimg // "[" // l_sci_ext // "," // version // "]"
        imb = "INDEF"; imc = "INDEF"; imd = "INDEF"
        if (dosky) 
            imb = skyimg // "[" // l_sci_ext // "," // version // "]"
        if (dodark)
            imc = darkimg // "[" // l_sci_ext // "," // version // "]"
        if (doflat)
            imd = flatimg // "[" // l_sci_ext // "," // version // "]"

        result = outimg // "[" // l_sci_ext // "," // version // ",append]"

        # Measure scale section

        if (intdbg)
            print ("scale: isvalue: " // scaleisvalue // \
                "; issection: " // scaleissection)
        if (doscale && dosky) {
            if (no == scaleisvalue) {
                secn = l_skysection
                if (no == scaleissection) {
                    keypar (ima, l_skysection, silent+)
                    if (no == keypar.found) {
                        printlog ("ERROR - NSREDUCE: No value found for \
                            skysection in header (" // l_skysection \
                            // ") for " // root, l_logfile, verbose+)
                        goto clean
                    }
                    secn = keypar.value
                }

                # Statistics (3-sig clipped around midpt)
                # (process taken from nsflat.cl)
                #printlog ("Scale level from " // secn // " for " \
                #    // root // "[" // l_sci_ext // "," // version \
                #    // "]", l_logfile, l_verbose)
                imstat (ima // secn, fields = "midpt,stddev", \
                    lower = INDEF, upper = sat, nclip = 0, \
                    lsigma = INDEF, \
                    usigma = INDEF, binwidth = 0.1, format-, cache-) \
                    | scan (temp, stddev) 
                printf ("Median, sigma %51s: %6.2f %6.2f\n", \
                    ima // secn, temp, stddev) | scan (sline)
                printlog (sline, l_logfile, l_debug)
                imstat (ima // secn, fields = "midpt", \
                    lower = (temp-3*stddev), upper = (temp+3*stddev), \
                    nclip = 0, lsigma = INDEF, usigma = INDEF, \
                    binwidth = 0.1, format-, cache-) | scan (scalevalue)
                printf ("Data level                                    \
                                       : %6.2f\n", scalevalue) | scan (sline)
                printlog (sline, l_logfile, l_debug)

                # If we are subtracting a dark then we need to reduce
                # the scale by the relevant amount
                if (dodark) {
                    imstat (imc // secn, fields = "midpt,stddev", \
                        lower = INDEF, upper = sat, nclip = 0, \
                        lsigma = INDEF, usigma = INDEF, binwidth = 0.1, \
                        format-, cache-) | scan (temp, stddev) 
                    printf ("Median, sigma %51s: %6.2f %6.2f\n", \
                        imc // secn, temp, stddev) | scan (sline)
                    printlog (sline, l_logfile, l_debug)
                    imstat (imc // secn, fields = "midpt", \
                        lower = (temp-3*stddev), \
                        upper = (temp+3*stddev), nclip = 0, \
                        lsigma = INDEF, usigma = INDEF, binwidth = 0.1, \
                        format-, cache-) | scan (temp)
                    scalevalue = scalevalue - temp
                    printf ("Level without dark                        \
                                           : %6.2f %6.2f\n", scalevalue, \
                        temp) | scan (sline)
                    printlog (sline, l_logfile, l_debug)
                }

                # And normalize by the sky level
                imstat (imb // secn, fields = "midpt,stddev", \
                    lower = INDEF, upper = sat, nclip = 0, \
                    lsigma = INDEF, usigma = INDEF, binwidth = 0.1, \
                    format-, cache-) | scan (temp, stddev) 
                printf ("Median, sigma %51s: %6.2f %6.2f\n", \
                    imb // secn, temp, stddev) | scan (sline)
                printlog (sline, l_logfile, l_debug)
                imstat (imb // secn, fields = "midpt", \
                    lower = (temp-3*stddev), \
                    upper = (temp+3*stddev), nclip = 0, \
                    lsigma = INDEF, usigma = INDEF, binwidth = 0.1, \
                    format-, cache-) | scan (temp)
                scalevalue = scalevalue / temp
                printf ("Sky level and final scaling                   \
                                       : %6.2f %6.2f\n", temp, scalevalue) \
                    | scan (sline)
                printlog (sline, l_logfile, l_debug)
            }
        } else {
            scalevalue = 1.0
        }

        # Substitute scale level
        if (intdbg)
            print ("flags: " // dosky // ", " // doscale // ", " // doflat)
        if (intdbg) print ("expression: " // expression)
        if (intdbg && dosky && doscale) print ("scalevalue: " // scalevalue)

        if (dosky && doscale) {
            if (intdbg) print ("evaluating " // expression)
            printf (expression // "\n", scalevalue) | scan (sline)
            subexpr = sline
            if (doflat) {
                if (intdbg) print ("evaluating (var, flat) " // varexpression)
                printf (varexpression // "\n", scalevalue * \
                    scalevalue, scalevalue, scalevalue) | scan (sline)
            } else {
                if (intdbg) print ("evaluating (var) " // varexpression)
                printf (varexpression // "\n", scalevalue * \
                    scalevalue) | scan (sline)
            }
            subvarexpr = sline
        } else {
            subexpr = expression
            subvarexpr = varexpression
        }

        if (intdbg) print ("sci: " // subexpr)
        if (intdbg) {
            print ("a: " // ima)
            print ("b: " // imb)
            print ("c: " // imc)
            print ("d: " // imd)
        }
        imexpr (subexpr, result, ima, imb, imc, imd, \
            outtype = "real", verbose-) 

        if (dovardq) {
            ime = inimg // "[" // l_var_ext // "," // version // "]"
            imf = "INDEF"; img = "INDEF"; imh = "INDEF"
            if (dosky) 
                imf = skyimg // "[" // l_var_ext // "," // version // "]"
            if (dodark)
                img = darkimg // "[" // l_var_ext // "," // version // "]"
            if (doflat)
                imh = flatimg // "[" // l_var_ext // "," // version // "]"

            result = outimg // "[" // l_var_ext // "," // version \
                // ",append]"

            # protect from overflow
            if (doflat && l_flatmin > 0.0) {
                subvarexpr = "d > " // flatmin // " ? " // subvarexpr \
                    // " : 0.0"
            }

            if (intdbg) print ("var: " // subvarexpr)
            if (intdbg) {
                print ("e: " // ime)
                print ("f: " // imf)
                print ("g: " // img)
                print ("h: " // imh)
            }
            imexpr (subvarexpr, result, \
                ima, imb, imc, imd, ime, imf, img, imh, \
                outtype = "real", verbose-) 

            dqin = ""
            dqexpression = ""
            dqcount = 0
            first = yes

            scanin2 = tmpdqlist
            while (fscan (scanin2, dqimg) != EOF) {
                dqcount = dqcount + 1
                if (first) {
                    first = no
                } else {
                    dqin = dqin // ","
                    dqexpression = dqexpression // " || "
                }
                dqin = dqin // dqimg // "[" // l_dq_ext // "," // version \
                    // "]"
                dqexpression = dqexpression // "im" // dqcount
            }

            imdelete (tmpdq, verify-, >>& "dev$null") 
            result = tmpdq // ".fits"

            if (intdbg) {
                print ("dq: " // dqexpression)
                print (dqin)
                print (dqcount // ", " // result)
            }
            if (dqcount > 1) {
                addmasks (dqin, result, dqexpression)
            } else {
                imcopy (dqin, result, verbose-)
            }

            if (intdbg) print (outimg)
            imcopy (result, outimg // "[" // l_dq_ext // "," \
                // version // ",append]", verbose-) 
        }
        if (intdbg) print ("vardq done")

        if (dolog) {
            if (intdbg) print ("dolog")

            # Log the file
            # TODO - should we also log extensions?

            if (first2) {
                printlog ("NSREDUCE: slit: " // slit, l_logfile, l_verbose)
                if (f2)
                    printlog ("NSREDUCE: filter: " // filter, l_logfile, \
                        l_verbose)
                else
                    printlog ("NSREDUCE: grism: " // filter, l_logfile, \
                        l_verbose)
                printlog ("", l_logfile, l_verbose)
                printlog ("   n input                   --> output   " \
                    // "               sky image", l_logfile, l_verbose) 
                printlog (" dark                    \
                    flat                    scale", l_logfile, l_verbose) 
                printlog ("", l_logfile, l_verbose) 
                first2 = no
            }

            # TODO - this is odd (poorly defined, and why "- 6"?)
            if ("" != slit) slit = substr (slit, 1, strlen (slit) - 6)
            if ("none" == skylist) {
                printf ("%4.0d %-23s --> %-23s %-23s\n", \
                    count, root, outimg, skyimg) | scan (sline) 
                printlog (sline, l_logfile, l_verbose) 
            } else {
                if (intdbg) print ("opening " // skylist)
                if (intdbg) type (skylist)
                first = yes
                scanin2 = skylist
                while (fscan (scanin2, skyimg) != EOF) {
                    if (first) {
                        printf ("%4.0d %-23s --> %-23s %-23s\n", count, \
                            root, outimg, skyimg) | scan (sline) 
                        first = no
                    } else {
                        printf ("%4s %-23s --> %-23s %-23s\n", \
                            "", "", "", skyimg) | scan (sline) 
                    }
                    printlog (sline, l_logfile, l_verbose) 
                }
            }

            printf (" %-23s %-23s %6.4f %s\n", \
            darkimg, flatimg, scalevalue, l_skysection) | scan (sline) 
            printlog (sline, l_logfile, l_verbose) 

            dolog = no
        }

        if ("none" == skyimg && l_fl_sky) {
            printlog ("WARNING - NSREDUCE: No sky found for " \
                // root, l_logfile, verbose+) 
        }

        # Update the header
        gemdate ()
        result = outimg // "[0]"
        gemhedit (result, "NSREDUCE", gemdate.outdate, 
            "UT Time stamp for NSREDUCE", delete-)
        gemhedit (result, "GEM-TLM", gemdate.outdate, 
            "UT Last modification with GEMINI", delete-)
        if (dosky)
            gemhedit (result, "SKYIMAGE", skyimg, \
                "Sky image subtracted from raw data", delete-) 
        if (dodark)
            gemhedit (result, "DARKIMAG", l_darkimage, \
                "Dark current image subtracted from raw data", delete-) 
        if (doflat)
            gemhedit (result, "FLATIMAG", l_flatimage, \
                "Flat field image used", delete-) 

        # No longer meaningful
        if (dosky || dodark || doflat) {
            gemhedit (result, l_key_ron, "", "", delete+)
            gemhedit (result, l_key_sat, "", "", delete+)
            gemhedit (result, l_key_nonlinear, "", "", delete+)
            gemhedit (result, l_key_gain, "", "", delete+)
            # For some very strange reason, the GAIN keyword magically 
            # re-appears in the output image (as defined by result above) 
            # after the imexpr call on line 1793 below, but only when 
            # using pyraf, even though the keyword is deleted above. 
            # Deleting it again here seems to solve the problem - EH.
            gemhedit (result, l_key_gain, "", "", delete+)
        }

    } # Next extension

    # Stage C: cut after processing (if l_fl_process_cut=no)

    if (intdbg)
        print ("cut: " // l_fl_cut // "; process: " // l_fl_process_cut)

    delete (tmpcutin, verify-, >& "dev$null")
    delete (tmpcutout, verify-, >& "dev$null")

    scanin1 = tmpstageF
    call_cut = no
    while (fscan (scanin1, root, inimg, skyimg, skylist, version, \
        dovardq, outimg) != EOF) {

        if (cut_after_proc) {

            if (intdbg) print ("need to cut " // outimg)
            call_cut = yes
            tmpimg = mktemp (outimg // "-")
            print (tmpimg, >> tmpcutin)
            print (tmpimg, >> tmpstageDimg)
            imrename (outimg, tmpimg, verbose-, >& "dev$null")
            print (outimg, >> tmpcutout)
        }
    }

    if (call_cut) {

        printlog (" ", l_logfile, l_verbose) 
        printlog ("NSREDUCE: Cutting processed data", l_logfile, l_verbose)
        if (intdbg) print ("call cut")
        if (f2) {
            f2cut (inimages="@"//tmpcutin, outimages="@"//tmpcutout, \
                outprefix="", fl_vardq=l_fl_vardq, gradimage=l_gradimage, \
                outslitim=l_outslitim, edgediff=l_edgediff, \
                refimage=l_refimage, database=l_database, logfile=l_logfile, \
                verbose=l_verbose)
            if (f2cut.status != 0) {
                printlog ("ERROR - NSREDUCE: Error in f2cut.", l_logfile, \
                    verbose+)
                goto clean
            }
        } else {
            nscut (inimages="@"//tmpcutin, outspectra="@"//tmpcutout, \
                outprefix="", section=l_section, fl_corner=l_fl_corner, \
                logfile=l_logfile, verbose=l_verbose, debug=l_debug)
            if (0 != nscut.status) {
                printlog ("ERROR - NSREDUCE: Error in nscut.", \
                    l_logfile, verbose+) 
                goto clean
            }
        }
    }

    scanin1 = tmpstageF
    call_nsappwave = no
    while (fscan (scanin1, root, inimg, skyimg, skylist, version, \
        dovardq, outimg) != EOF) {

        if (cut_after_proc) {

            if (l_fl_nsappwave) {
                if (intdbg) print ("need to call nsappwave for " // outimg)
                call_nsappwave = yes
                tmpimg = mktemp (outimg // "-")
                print (tmpimg, >> tmpnsappwavein)
                print (tmpimg, >> tmpstageDimg)
                imrename (outimg, tmpimg, verbose-, >& "dev$null")
                print (outimg, >> tmpnsappwaveout)
            }
        }
    }

    if (call_nsappwave) {

        if (intdbg) print ("call nsappwave")
        nsappwave (inimages="@"//tmpnsappwavein, \
            outspectra="@"//tmpnsappwaveout, outprefix="", \
            nsappwavedb=l_nsappwavedb, fl_phu-, crval=l_crval, \
            cdelt=l_cdelt, logfile=l_logfile, verbose=l_verbose, \
            debug=l_debug)
        if (0 != nsappwave.status) {
            printlog ("ERROR - NSREDUCE: Error in nsappwave.", \
                l_logfile, verbose+) 
            goto clean
        }
    }

    # We're done, so clear status
    status = 0

clean:

    # Stage D - delete files

    if (intdbg) print ("stageD: delete")

    if (l_fl_sky) {
        if (intdbg) {
            print ("deleting output from nssky")
            type (tmpskyout)
        }
        if (access (tmpskyout)) {
            scanin1 = tmpskyout
            while (fscan (scanin1, tmpimg, skyimg, ignore, \
                skylist, dupsky) != EOF) {
                
                if (no == dupsky) {
                    scanin2 = skylist
                    count = 0
                    while (fscan (scanin2, ignore) != EOF) {
                        count = count + 1
                    }
                    if (count > 1) {
                        if (intdbg) print ("composite sky for " // tmpimg)
                        imdelete (skyimg, verify-, >& "dev$null")
                    }
                    if (intdbg) print ("deleting " // skylist)
                    delete (skylist, verify-, >& "dev$null")
                }
            }
        }
    }
    
    if (access (tmpstageDimg)) {
        imdelete ("@" // tmpstageDimg, verify-, >& "dev$null")
    }

    if (access (tmpstageD)) {
        delete ("@" // tmpstageD, verify-, >& "dev$null")
        delete (tmpstageD, verify-, >& "dev$null")
    }

    printlog (" ", l_logfile, l_verbose) 
    if (status == 0) {
        printlog ("NSREDUCE exit status: good.", l_logfile, l_verbose) 
    }
    printlog ("-------------------------------------------------------\
        -----------------------", l_logfile, l_verbose) 

    scanin1 = ""
    scanin2 = ""

end

# Copyright(c) 2001-2016 Association of Universities for Research in Astronomy, Inc.

procedure nsprepare (inimages) 

# Pre-process GNIRS data
#
# The following steps are taken:
# - information is read from config.fits and array.fits, based on headers
# - saturation level and linear limits are calculated
# - data values are corrected by number of lnrs if this isn't already done
# - the offset ("random bias") level is measured and subtracted
# - the data are corected for non-linearity
# - variance data are generated from science data (photon stats)
# - data quality data are generated from data levels and bpm
# - mdf is attached
# - wcs is checked/corrected
#
# See "help nsprepare" for more details

char    inimages    {prompt = "Input GNIRS image(s)"}                   # OLDP-1-primary-single-prefix=n
char    rawpath     {"", prompt = "Path for input raw images"}          # OLDP-4
char    outimages   {"", prompt = "Output image(s)"}                    # OLDP-1-output
char    outprefix   {"n", prompt = "Prefix for output image(s)"}        # OLDP-4
char    bpm         {"gnirs$data/gnirsn_2012dec05_bpm.fits", prompt = "Bad pixel mask file"}                # OLDP-1

char    logfile     {"", prompt = "Logfile"}                            # OLDP-1

bool    fl_vardq    {yes, prompt = "Create variance and data quality frames?"}  # OLDP-2
bool    fl_cravg    {no, prompt = "Attempt to flag cosmic rays and radiation events?"} # OLDP-2
int     crradius    {0, prompt = "Grow radius for bad pixels"}          # OLDP-2
bool    fl_dark_mdf {no, prompt = "Force the attachment of MDFs to dark frames?"} # OLDP-3
bool    fl_correct  {yes, prompt = "Correct for non-linearity in the data?"}    # OLDP-3
bool    fl_saturated    {yes, prompt = "Flag saturated pixels in DQ?"}  # OLDP-3
bool    fl_nonlinear    {yes, prompt = "Flag non-linear pixels in DQ?"} # OLDP-3
bool    fl_checkwcs {yes, prompt = "Check the WCS?"}                    # OLDP-3
bool    fl_forcewcs {yes, prompt = "Force the WCS to be what we expect?"}    # OLDP-3

char    arraytable  {"gnirs$data/array.fits", prompt = "Table of array settings"} # OLDP-3
char    configtable {"gnirs$data/config.fits", prompt = "Table of configuration settings"} # OLDP-3
char    specsec     {"[*,*]", prompt = "Section or MDF file if config table missing"} # OLDP-3
char    offsetsec   {"none", prompt = "Offset region (none|config|section|value)"} # OLDP-3
real    pixscale    {0.15, prompt = "Pixel scale if config table missing"}   # OLDP-3
char    shiftimage  {"", prompt = "Image to copy shift from"}           # OLDP-3
real    shiftx      {0., prompt = "Shift in X for MDF (subtracted from x_ccd)"}  # OLDP-3
real    shifty      {0., prompt = "Shift in Y for MDF (subtracted from y_ccd)"}  # OLDP-3
char    obstype     {"FLAT", prompt = "Observation type for cross-correlation"} # OLDP-3
bool    fl_inter    {no, prompt = "Fit offset interactively?"}          # OLDP-3
bool    verbose     {yes, prompt = "Verbose"}                           # OLDP-2
int     status      {0, prompt = "Exit status (0=good)"}                # OLDP-4

struct  *scanin1    {"", prompt = "Internal use only"}                  # OLDP-4
struct  *scanin2    {"", prompt = "Internal use only"}                  # OLDP-4

begin
    char    l_inimages = ""
    char    l_rawpath = ""
    char    l_outimages = ""
    char    l_outprefix = ""
    char    l_bpm = ""
    char    l_logfile = ""
    char    l_sci_ext = ""
    char    l_var_ext = ""
    char    l_dq_ext = ""
    bool    l_fl_vardq
    bool    l_fl_cravg
    int     l_crradius
    bool    l_fl_dark_mdf
    bool    l_fl_correct
    bool    l_fl_saturated
    bool    l_fl_nonlinear
    bool    l_fl_checkwcs
    bool    l_fl_forcewcs
    char    l_key_ndavgs = ""
    char    l_key_coadds = ""
    char    l_key_lnrs = ""
    char    l_key_bias = ""
    char    l_key_ron = ""
    char    l_key_gain = ""
    char    l_key_sat = ""
    char    l_key_nonlinear = ""
    char    l_key_arrayid = ""
    char    l_key_filter = ""
    char    l_key_grating = ""
    char    l_key_decker = ""
    char    l_key_camera = ""
    char    l_key_prism = ""
    char    l_arraytable = ""
    char    l_configtable = ""
    char    l_specsec = ""
    char    l_key_section = ""
    char    l_key_exptime = ""
    char    l_key_dispaxis = ""
    char    l_key_dark = ""
    char    l_key_mode = ""
    char    l_val_dark = ""
    char    l_offsetsec = ""
    real    l_pixscale
    char    l_shiftimage = ""
    real    l_shiftx
    real    l_shifty
    char    l_obstype = ""
    bool    l_fl_inter
    bool    l_verbose
    char    l_key_obstype = ""
    char    l_key_instrument = ""

    char    inimg, outimg, rootimg, extn, phu, sections, expr, oldmdf
    char    tmpexpand, tmpfile, tmpfile2, tmproot, tmpsummary, tmpout
    char    tmphead, tmpopen, tmpmdfrow, tmparray, tmpsections
    char    tmpmdffile, tmpoffset, tmpnotset
    bool    already, domdf, havemdf, havedark, correlate
    int     nimages, junk
    struct  sdate
    int     bpmxsize, bpmysize, xsize, ysize, sx1, sx2, sy1, sy2, len
    char    tmpsci, tmpvar, tmpdq, tmpinimage, tmpcoordsin
    int     noutimages, nbad, dispaxis
    bool    bad, oldtimer, isgnirs, issouth, debug, hasdate, firstmdf, first
    real    well, buas, ronref, gain
    real    linearlimit, ron, sat, ndavgs, coadds, lnrs, exptime
    real    nonlinlimit, coeff1, coeff2, coeff3
    real    vdduc, vdet, linlimit
    real    biasvolt
    char    arrayid
    char    filter, filter1, filter2, filter3, pupil, fpmask, shiftfile
    char    varexpression, tempexpression, temp, word, badhdr, dark
    real    pixscl, cval1, cval2, x
    char    camera, prism, grating, decker, offsetsection, specsection
    char    prismsearch
    char    date, keyfound, procmode, NONE
    int     count, year, day, month
    bool    fixshift, fixwcs, wcsok, defshift, found, others
    real    dec, ra, decoff, raoff, pix, pa, cd11, cd12, cd21, cd22
    real    ranod, decnod, xobj, yobj, xpix, ypix, rmax
    real    cd11obs, cd12obs, cd21obs, cd22obs, jacob
    real    pi, cospa, sinpa, cosdec, torad, abserr, tiny, errlim
    real    offset, stddev, midpt
    int     offsettype, OFF_NONE, OFF_CONF, OFF_SECN, OFF_VAL

    struct  l_struct

    status = 1
    debug = no
    pi = 3.1415926
    torad = 2*pi/360.0
    tiny = 1e-10
    errlim = 0.05 # 5% of a pixel error
    OFF_NONE = 1
    OFF_CONF = 2
    OFF_SECN = 3
    OFF_VAL = 4
    NONE = "none"
    offsettype = 0 # none of the above
    offset = 0
    offsetsection = NONE
    xobj = INDEF
    yobj = INDEF
    rmax = 2.0 ** 2

    cache ("gemextn", "keypar", "tinfo", "tabpar", "nsoffset", "gemdate")

    junk = fscan (  inimages, l_inimages)
    junk = fscan (  rawpath, l_rawpath)
    junk = fscan (  outimages, l_outimages)
    junk = fscan (  outprefix, l_outprefix)
    junk = fscan (  bpm, l_bpm)
    junk = fscan (  logfile, l_logfile)
    l_fl_vardq =    fl_vardq
    l_fl_cravg =    fl_cravg
    l_crradius =    crradius
    l_fl_dark_mdf = fl_dark_mdf
    l_fl_correct =  fl_correct
    l_fl_saturated =    fl_saturated
    l_fl_nonlinear =    fl_nonlinear
    l_fl_checkwcs = fl_checkwcs
    l_fl_forcewcs = fl_forcewcs
    junk = fscan (  arraytable, l_arraytable)
    junk = fscan (  configtable, l_configtable)
    junk = fscan (  specsec, l_specsec)
    l_offsetsec =   offsetsec # may contain spaces, text
    l_pixscale =    pixscale
    junk = fscan (  shiftimage, l_shiftimage)
    l_shiftx =      shiftx
    l_shifty =      shifty
    junk = fscan (  obstype, l_obstype)
    l_fl_inter =    fl_inter
    l_verbose =     verbose

    # Section may contain spaces
    if ("" != l_specsec)
        if (substr (l_specsec, 1, 1) == "[")
            l_specsec = specsec

    tmpexpand = mktemp("tmpexpand")
    tmpfile = mktemp ("tmpfile") 
    tmproot = mktemp ("tmproot") 
    tmpout = mktemp ("tmpout") 
    tmpsummary = mktemp ("tmpsummary") 
    tmpfile2 = mktemp ("tmpfile2") 
    tmphead = mktemp ("tmphead") 
    tmpopen = mktemp ("tmpopen") 
    tmpmdfrow = mktemp ("tmpmdfrow") //".fits"
    tmparray = mktemp ("tmparray")//".fits"
    tmpsections = mktemp ("tmpsections")
    tmpmdffile = mktemp ("tmpmdffile")
    tmpcoordsin = mktemp ("tmpcoordsin")
    tmpnotset = mktemp ("tmpnotset")

    tmpvar = ""
    tmpsci = ""
    tmpdq = ""
    tmpinimage = ""
    firstmdf = yes

    # Shared definitions
    badhdr = ""
    junk = fscan (nsheaders.sci_ext, l_sci_ext)
    if ("" == l_sci_ext) badhdr = badhdr + " sci_ext"
    junk = fscan (nsheaders.var_ext, l_var_ext)
    if ("" == l_var_ext) badhdr = badhdr + " var_ext"
    junk = fscan (nsheaders.dq_ext, l_dq_ext)
    if ("" == l_dq_ext) badhdr = badhdr + " dq_ext"
    junk = fscan (nsheaders.key_ndavgs, l_key_ndavgs)
    if ("" == l_key_ndavgs) badhdr = badhdr + " key_ndavgs"
    junk = fscan (nsheaders.key_coadds, l_key_coadds)
    if ("" == l_key_coadds) badhdr = badhdr + " key_coadds"
    junk = fscan (nsheaders.key_lnrs, l_key_lnrs)
    if ("" == l_key_lnrs) badhdr = badhdr + " key_lnrs"
    junk = fscan (nsheaders.key_bias, l_key_bias)
    if ("" == l_key_bias) badhdr = badhdr + " key_bias"
    junk = fscan (nsheaders.key_ron, l_key_ron)
    if ("" == l_key_ron) badhdr = badhdr + " key_ron"
    junk = fscan (nsheaders.key_gain, l_key_gain)
    if ("" == l_key_gain) badhdr = badhdr + " key_gain"
    junk = fscan (nsheaders.key_sat, l_key_sat)
    if ("" == l_key_sat) badhdr = badhdr + " key_sat"
    junk = fscan (nsheaders.key_nonlinear, l_key_nonlinear)
    if ("" == l_key_nonlinear) badhdr = badhdr + " key_nonlinear"
    junk = fscan (nsheaders.key_arrayid, l_key_arrayid)
    if ("" == l_key_arrayid) badhdr = badhdr + " key_arrayid"
    junk = fscan (nsheaders.key_filter, l_key_filter)
    if ("" == l_key_filter) badhdr = badhdr + " key_filter"
    junk = fscan (nsheaders.key_grating, l_key_grating)
    if ("" == l_key_grating) badhdr = badhdr + " key_grating"
    junk = fscan (nsheaders.key_decker, l_key_decker)
    if ("" == l_key_decker) badhdr = badhdr + " key_decker"
    junk = fscan (nsheaders.key_camera, l_key_camera)
    if ("" == l_key_camera) badhdr = badhdr + " key_camera"
    junk = fscan (nsheaders.key_prism, l_key_prism)
    if ("" == l_key_prism) badhdr = badhdr + " key_prism"
    junk = fscan (nsheaders.key_section, l_key_section)
    if ("" == l_key_section) badhdr = badhdr + " key_section"
    junk = fscan (nsheaders.key_exptime, l_key_exptime)
    if ("" == l_key_exptime) badhdr = badhdr + " key_exptime"
    junk = fscan (nsheaders.key_dispaxis, l_key_dispaxis)
    if ("" == l_key_dispaxis) badhdr = badhdr + " key_dispaxis"
    junk = fscan (nsheaders.key_dark, l_key_dark)
    if ("" == l_key_dark) badhdr = badhdr + " key_dark"
    junk = fscan (nsheaders.key_mode, l_key_mode)
    if ("" == l_key_mode) badhdr = badhdr + " key_mode"
    junk = fscan (nsheaders.val_dark, l_val_dark)
    if ("" == l_val_dark) badhdr = badhdr + " val_dark"
    junk = fscan (nsheaders.key_obstype, l_key_obstype)
    if ("" == l_key_obstype) badhdr = badhdr + " key_obstype"
    junk = fscan (nsheaders.key_instrument, l_key_instrument)
    if ("" == l_key_instrument) badhdr = badhdr + " key_instument"

    if (debug) print ("validating input")

    if ("" == l_logfile) {
        junk = fscan (gnirs.logfile, l_logfile)
        if ("" == l_logfile) {
            l_logfile = "gnirs.log"
            printlog ("WARNING - NSPREPARE: Both nsprepare.logfile and \
                gnirs.logfile are", l_logfile, verbose+) 
            printlog ("                     undefined.  Using " \
                // l_logfile, l_logfile, verbose+)
        }
    }
    date | scan (sdate) 
    printlog ("----------------------------------------------------------\
        ---------------------", l_logfile, l_verbose) 
    printlog ("NSPREPARE -- " // sdate, l_logfile, l_verbose) 
    printlog (" ", l_logfile, l_verbose) 

    if ("" != badhdr) {
        printlog ("ERROR - NSPREPARE: Parameter(s) missing from \
            nsheaders: " // badhdr, l_logfile, verbose+) 
        goto clean
    }
    if (l_fl_vardq && ! l_fl_saturated) {
        printlog ("WARNING - NSPREPARE: Saturated pixels not flagged.", \
            l_logfile, verbose+) 
    }
    if (l_fl_vardq && ! l_fl_nonlinear) {
        printlog ("WARNING - NSPREPARE: Non-linear pixels not flagged.", \
            l_logfile, verbose+) 
    }

    if (debug) print ("shift options")

    shiftfile = ""
    gemextn (l_shiftimage, check="exists", process="expand", \
        index="0", extname="", extversion="", ikparams="", omit="", \
        replace="", outfile="STDOUT", logfile="dev$null", glogpars="", \
        verbose=l_verbose) | scan (shiftfile)
    if ("" != shiftfile) {
        printlog ("NSPREPARE: Shift from " // shiftfile, \
            l_logfile, l_verbose)
        hselect (shiftfile, "MDF_XSHF", yes) | scan (l_shiftx)
        hselect (shiftfile, "MDF_YSHF", yes) | scan (l_shifty)
    } else if (no == ("" == l_shiftimage)) {
        printlog ("NSPREPARE: ERROR - Cannot open shiftimage " \
            // l_shiftimage, l_logfile, verbose+)
        goto clean
    }

    correlate = isindef (l_shiftx) && isindef (l_shifty)
    defshift = ! isindef (l_shiftx) && ! isindef (l_shifty)
    if (defshift)
        defshift = l_shiftx == 0 && l_shifty == 0

    if (debug) print ("check configtable")

    domdf = access (l_configtable)
    if (no == domdf) {
        printlog ("WARNING - NSPREPARE: Config table " // l_configtable \
            // " not found.", l_logfile, verbose+) 
    }

    if (debug) print ("check offset section")

    l_offsetsec = strlwr (l_offsetsec)
    if (l_offsetsec == NONE) {
        offsettype = OFF_NONE
    } else if (l_offsetsec == "config") {
        if (no == domdf) {
            printlog ("ERROR - NSPREPARE: Cannot read offset section \
                from missing config table.", l_logfile, verbose+) 
            goto clean
        }
        offsettype = OFF_CONF
    } else {
        # drop surrounding whitespace
        while (l_offsetsec != "" && substr (l_offsetsec, 1, 1) == " ") {
            l_offsetsec = substr (l_offsetsec, 2, strlen (l_offsetsec))
        }
        while (l_offsetsec != "" && \
            substr (l_offsetsec, strlen (l_offsetsec), \
            strlen (l_offsetsec)) == " ") {
            
            l_offsetsec = substr (l_offsetsec, 1, strlen (l_offsetsec) - 1)
        }
        if (debug) print ("trimmed l_offsetsec to " // l_offsetsec)
        if (l_offsetsec == "") {
            printlog ("ERROR - NSPREPARE: offsetsec must be 'none', \
                'config', a section, or a value.", l_logfile, verbose+)
            goto clean
        }
        if (substr (l_offsetsec, 1, 1) == "[" && \
            substr (l_offsetsec, strlen (l_offsetsec), \
            strlen (l_offsetsec)) == "]") {
            
            offsettype = OFF_SECN
            offsetsection = l_offsetsec
        } else {
            offset = INDEF
            print (l_offsetsec) | scan (offset)
            if (isindef (offset) || 0 == nscan ()) {
                printlog ("ERROR - NSPREPARE: offsetsec must be 'none', \
                    'config', a section, or a value.", l_logfile, verbose+)
                goto clean
            }
            offsettype = OFF_VAL
            offsetsection = "value"
        }
    }

    if (offsettype == OFF_NONE) {
        printlog ("NSPREPARE: No offset subtraction.", \
            l_logfile, l_verbose)
    } else if (offsettype == OFF_CONF) {
        printlog ("NSPREPARE: Offset section from configtable.", \
            l_logfile, l_verbose)
    } else if (offsettype == OFF_SECN) {
        printlog ("NSPREPARE: Offset section is " // offsetsection, \
            l_logfile, l_verbose)
    } else if (offsettype == OFF_VAL) {
        printlog ("NSPREPARE: Offset value is " // offset, \
            l_logfile, l_verbose)
    } else {
        printlog ("ERROR - NSPREPARE: Uncertain offset.", \
            l_logfile, verbose+)
        goto clean
    }

    if (debug) print ("process in and output files")

    # Check the rawpath name for a final /
    if ("" != l_rawpath) {
        len = strlen (rawpath)
        if (substr (l_rawpath, len, len) != "/")
            l_rawpath = l_rawpath // "/"
    }

    # Expand list (in case comma separated etc)
    gemextn (l_inimages, check="", process="none", index="", \
        extname="", extversion="", ikparams="", omit="", \
        replace="", outfile=tmpexpand, logfile="", glogpars="",
        verbose=l_verbose)
    if (0 < gemextn.fail_count) {
        printlog ("ERROR - NSPREPARE: Bad syntax in inimages.", \
            l_logfile, verbose+) 
        goto clean
    }
    if (0 == gemextn.count) {
        printlog ("ERROR - NSPREPARE: No input images.", \
            l_logfile, verbose+) 
        goto clean
    }

    # Add directory to start and validate
    gemextn ("@" // tmpexpand, check="exists,mef", process="none", \
        index="", extname="", extversion="", ikparams="", \
        omit="kernel,exten", replace="^%%" // l_rawpath // "%", \
        outfile=tmpfile, logfile="", glogpars="", verbose=l_verbose)
    if (0 < gemextn.fail_count) {
        printlog ("ERROR - NSPREPARE: Input images missing.", \
            l_logfile, verbose+) 
        goto clean
    }
    nimages = gemextn.count
    
    # Confirm that this is GNIRS data.  Users are known to try nsprepare
    # on NIRI data, since everything else starts with "ns"
    scanin1 = tmpfile
    while (fscan (scanin1, inimg) != EOF) {
        keypar (inimg//"[0]", l_key_instrument, silent+) 
        if (no == keypar.found) {
            printlog ("ERROR - NSPREAPRE: Image " // inimg, l_logfile,
                verbose+)
            printlog ("ERROR - NSPREPARE: Could not read " // l_key_instrument \
                // " from header.", l_logfile, verbose+) 
            goto clean
        } else if (keypar.value != "GNIRS") {
            printlog ("ERROR - NSPREAPRE: Image " // inimg, l_logfile,
                verbose+)
            printlog ("ERROR - NSPREPARE: nsprepare is for GNIRS data only. \n" \
                // "                   Use the 'prepare' task in the package "\
                // " for your instrument.",
                l_logfile, verbose+)
            # NIRI is a common mistake.
            if (keypar.value == "NIRI") {
                printlog ("ERROR - NSPREPARE: This is NIRI data, use 'nprepare'",
                l_logfile, verbose+)
            }
            goto clean
        }
    }
    scanin1 = ""

    # Generate output images
    gemextn (l_outimages, check="absent", process="none", index="", \
        extname="", extversion="", ikparams="", omit="kernel,exten", \
        replace="", outfile=tmpout, logfile="", glogpars="",
        verbose=l_verbose)
    if (gemextn.fail_count != 0) {
        printlog ("ERROR - NSPREPARE: Existing or incorrectly formatted \
            output files", l_logfile, verbose+) 
        goto clean
    }

    # If tmpout is empty, the output files names should be 
    # created with the prefix parameter
    if (gemextn.count == 0) {
        gemextn ("@" // tmpfile, check="", process="none", index="", \
            extname="", extversion="", ikparams="", omit="path", \
            replace="", outfile=tmproot, logfile="", glogpars="",
            verbose=l_verbose)
        if (0 < gemextn.fail_count) {
            printlog ("ERROR - NSPREPARE: Unexpected error.", \
                l_logfile, verbose+) 
            goto clean
        }
        gemextn ("%^%" // l_outprefix // "%" // "@" // tmproot, \
            check="absent", process="none", index="", \
            extname="", extversion="", ikparams="", omit="kernel,exten", \
            replace="", outfile=tmpout, logfile="", glogpars="",
            verbose=l_verbose)
        if (gemextn.fail_count != 0 || gemextn.count == 0) {
            printlog ("ERROR - NSPREPARE: No or incorrectly formatted \
                output files", l_logfile, verbose+) 
            goto clean
        }
    }

    # Check for previous processing, matching of input to output,
    # and write to a summary file
    scanin1 = tmpfile
    scanin2 = tmpout
    while (fscan (scanin1, inimg) != EOF) {
        if (fscan (scanin2, outimg) == EOF) {
            printlog ("ERROR - NSPREPARE: Too few output files.", \
                l_logfile, verbose+) 
            goto clean
        }

        keyfound=""
        already=no
        hselect(inimg//"[0]", "*PREPAR*", yes) | scan(keyfound)
        if (keyfound != "") {
            already=yes
            printlog ("WARNING - NSPREPARE: Image " // inimg \
                // " already fixed using *PREPARE.", l_logfile, verbose+) 
            printlog ("                     No futher processing done.",
                l_logfile, verbose+) 
        }

        print (inimg // " " // outimg // " " // already, >> tmpsummary)
    }
    if (fscan (scanin2, outimg) != EOF) {
        printlog ("ERROR - NSPREPARE: Too many output files.", \
            l_logfile, verbose+) 
        goto clean
    }

    printlog ("Processing " // nimages // " files", l_logfile, l_verbose) 


    # Check for existence of bad pixel mask

    # TODO - add check for .pl or .fits, and convert .pl to fits
    # if required by addmasks

    if (debug) print ("checking bpm")

    if (l_fl_vardq && "" != l_bpm && ! imaccess (l_bpm)) {
        printlog ("WARNING - NSPREPARE: Bad pixel mask " // l_bpm \
            // " not found.", l_logfile, verbose+) 
        l_bpm = ""
    } else if (l_fl_vardq && "" == l_bpm) {
        printlog ("WARNING - NSPREPARE: Bad pixel mask not specified.", \
            l_logfile, verbose+) 
    }
    if ("" == l_bpm)
        l_bpm = "none" # For logging

    # Check to make sure BPM and input images are same size

    if ("none" != l_bpm) {
        if (debug) print ("checking sizes")
        bad = no
        keypar (l_bpm, "i_naxis1", silent+) 
        if (keypar.found)
            bpmxsize = int (keypar.value) 
        else
            bad = yes
        keypar (l_bpm, "i_naxis2", silent+) 
        if (keypar.found)
            bpmysize = int (keypar.value) 
        else
            bad = yes
        if (bad) {
            printlog ("WARNING - NSPREPARE: No size information in \
                bad pixel mask header.", l_logfile, verbose+) 
            printlog ("                    Not using BPM to generage \
                data quality planes.", l_logfile, verbose+) 
            l_bpm = "none"
        } else {

            scanin1 = tmpsummary
            while (fscan (scanin1, inimg, outimg, already) != EOF) {
                bad = no
                extn = inimg // "[1]"
                keypar (extn, "i_naxis1")
                if (keypar.found)
                    xsize = int (keypar.value)
                else
                    bad = yes
                keypar (extn, "i_naxis2")
                if (keypar.found)
                    ysize = int (keypar.value)
                else
                    bad = yes
                if (bad) {
                    printlog ("ERROR - NSPREPARE: No size information \
                        in " // extn // ".", l_logfile, verbose+) 
                    goto clean
                }
                if (xsize != bpmxsize || ysize != bpmysize) {
                    printlog ("WARNING - NSPREPARE: Input images and \
                        BPM are not the same size.", l_logfile, verbose+) 
                    printlog ("                    Not using BPM to \
                        generage data quality plane.", l_logfile, verbose+)
                    l_bpm = "none"
                    break
                }
            }
        }
    }

    # If we're shifting, we want the first item in the input list to
    # have an obstype that matches that given (to avoid asymmetric
    # data).

    if (correlate && no == ("" == l_obstype) && \
        no == ("(unused)" == l_key_obstype)) {

        tmpoffset = mktemp ("tmpsummary")
        scanin1 = tmpsummary
        found = no
        others = no
        while (fscan (scanin1, inimg, outimg, already) != EOF) {
            if (no == found && no == already) {
                phu = inimg // "[0]"
                hselect (phu, l_key_obstype, yes) | scan (l_struct)
                if (l_struct == l_obstype) {
                    print (inimg // " " // outimg // " " // already, \
                        > tmpoffset)
                    found = yes
                    printlog ("NSPREPARE: Will measure MDF offset using " \
                        // inimg, l_logfile, l_verbose)
                    printlog ("           (" // l_key_obstype // "=" \
                        // l_obstype // ")", l_logfile, l_verbose)
                } else {
                    print (inimg // " " // outimg // " " // already, \
                        >> tmpnotset)
                    others = yes
                }
            } else {
                print (inimg // " " // outimg // " " // already, \
                    >> tmpnotset)
                others = yes
            }
        }

        if (no == found) {
            printlog ("ERROR - NSPREPARE: No suitable candidate was found \
                for cross-correlation", l_logfile, verbose+)
            printlog ("                   (" // l_key_obstype // "!=" \
                // l_obstype // ")", l_logfile, verbose+)
            goto clean
        }

        delete (tmpsummary, verify-, >& "dev$null")
        tmpsummary = tmpoffset

        if (others) {
            scanin1 = tmpnotset
            while (fscan (scanin1, inimg, outimg, already) != EOF) {
                print (inimg // " " // outimg // " " // already, \
                    >> tmpoffset)
            }
        }
    }

    # Help users with automatic shifting etc
    if (correlate) {
        printlog ("NSPREPARE: Both shiftx and shifty are INDEF.  The \
            shift for MDF will be", l_logfile, l_verbose)
        printlog ("           measured using cross-correlation.", \
            l_logfile, l_verbose)
        printlog ("           You should verify that the measured \
            shift is correct by", l_logfile, l_verbose)
        printlog ("           examining the cut data.", \
            l_logfile, l_verbose)
    } else if (defshift) {
        printlog ("NSPREPARE: Both shiftx and shifty are 0.  This will \
            place the MDF at the", l_logfile, l_verbose)
        printlog ("           default position on the detector, which \
            may not be correct (if,", l_logfile, l_verbose)
        printlog ("           for example, the slit has shifted on the \
            detector because of", l_logfile, l_verbose)
        printlog ("           flexure).  You can verify this by \
            displaying a cut (using nscut,", l_logfile, l_verbose)
        printlog ("           or nsreduce with fl_cut+) image and \
            checking that all data are", l_logfile, l_verbose)
        printlog ("           visible.  Alternatively, by setting shiftx \
            ans shifty to INDEF,", l_logfile, l_verbose)
        printlog ("           the correct shift can be measured \
            automatically.", l_logfile, l_verbose)
    }
    if (correlate || defshift) {
        printlog ("           NOTE that it is important to ensure that \
            the same shift is used", l_logfile, l_verbose)
        printlog ("           for all data processed together, including \
            calibration frames.", l_logfile, l_verbose)
        printlog ("           So, if you use the automatic shift, either \
            process all the data", l_logfile, l_verbose)
        printlog ("           in just one call to nsprepare or use the \
            shiftimage parameter", l_logfile, l_verbose)
        printlog ("           (or shiftx and shifty directly) to copy \
            the shift from a ", l_logfile, l_verbose)
        printlog ("           previously processed image on subsequent \
            calls to nsprepare.", l_logfile, l_verbose)
    }

    filter1 = ""; filter2 = ""; filter3 = ""
    pupil = "none"; fpmask = ""

    first = yes
    scanin1 = tmpsummary
    count = 0
    while (fscan (scanin1, inimg, outimg, already) != EOF) {
        count = count + 1
        if (debug) print (count // ": " // inimg // " -> " // outimg)

        # tmp files used within this loop
        if (no == already) {
            tmpvar = mktemp ("tmpvar")
            tmpsci = mktemp ("tmpsci")
            tmpdq = mktemp ("tmpdq")
            tmpinimage = mktemp ("tmpinimage")
        }

        # Set read noise, gain, and saturation for the data

        phu = inimg // "[0]"
        havedark = no

        dispaxis = INDEF
        hselect (phu, l_key_dispaxis, yes) | scan (dispaxis)
        if (isindef (dispaxis)) {
            dispaxis = 2
        } else {
            printlog ("WARNING - NSPREPARE: " // l_key_dispaxis \
                // " already set to " // dispaxis, \
            l_logfile, verbose+) 
        }

        keypar (phu, l_key_dark, silent+) 
        if (no == keypar.found) {
            printlog ("WARNING - NSPREPARE: Could not read " \
                // l_key_dark // " for " // inimg, l_logfile, verbose+) 
            dark = ""
        } else {
            dark = keypar.value
        }

        if (debug)
            print (l_val_dark//", "//dark//", "//strstr (l_val_dark, dark))
        if (strstr (l_val_dark, dark) > 0) {
            if (no == l_fl_dark_mdf) {
                printlog ("WARNING - NSPREPARE: " // inimg \
                    // " is a dark", l_logfile, verbose+)
                printlog ("                     and so will not have \
                    an MDF added.", l_logfile, verbose+)
                havedark = yes
            } else {
                printlog ("WARNING - NSPREPARE: " // inimg // " is a dark, \
                    but fl_dark_mdf is set", l_logfile, verbose+)
                printlog ("                     so it will have an \
                    MDF added.", l_logfile, verbose+)
                havedark = no
            }
        }

        if (domdf && ! havedark) {

            keypar (phu, l_key_prism, silent+) 
            if (no == keypar.found) {
                printlog ("WARNING - NSPREPARE: Could not read prism \
                    for " // inimg, l_logfile, verbose+) 
                prism = ""
            } else {
                prism = keypar.value
                # Someone had the brilliant idea to put a plus-sign in
                # the value for the prism.  tselect doesn't like that.
                # So create a "search" string pattern that will be use later
                # tselect.
                print(prism) | translit("STDIN","+","?") | scan(prismsearch)
            }

            keypar (phu, l_key_decker, silent+) 
            if (no == keypar.found) {
                printlog ("WARNING - NSPREPARE: Could not read decker \
                    for " // inimg, l_logfile, verbose+) 
                decker = ""
            } else {
                decker = keypar.value
            }

            keypar (phu, l_key_camera, silent+) 
            if (no == keypar.found) {
                printlog ("WARNING - NSPREPARE: Could not read camera \
                    for " // inimg, l_logfile, verbose+) 
                camera = ""
            } else {
                camera = keypar.value
            }

            keypar (phu, l_key_grating, silent+) 
            if (no == keypar.found) {
                printlog ("WARNING - NSPREPARE: Could not read grating \
                    for " // inimg, l_logfile, verbose+) 
                grating = ""
            } else {
                grating = keypar.value
            }
        }

        keypar (phu, l_key_ndavgs, silent+) 
        if (no == keypar.found) {
            printlog ("ERROR - NSPREPARE: Could not read number of \
                digital averages from header", l_logfile, verbose+) 
            printlog ("                  of file " // inimg, 
                l_logfile, verbose+) 
            goto clean
        } else {
            ndavgs = real (keypar.value) 
        }

        keypar (phu, l_key_coadds, silent+) 
        if (no == keypar.found) {
            printlog ("ERROR - NSPREPARE: Could not read number of \
                coadds from header", l_logfile, verbose+) 
            printlog ("                  of file " // inimg, 
                l_logfile, verbose+) 
            goto clean
        } else {
            coadds = real (keypar.value) 
        }

        keypar (phu, l_key_lnrs, silent+) 
        if (no == keypar.found) {
            printlog ("ERROR - NSPREPARE: Could not read number of \
                non-destructive reads from header", l_logfile, verbose+) 
            printlog ("                  of file " // inimg, 
                l_logfile, verbose+) 
            goto clean
        } else {
            lnrs = real (keypar.value) 
        }

        keypar (phu, l_key_bias, silent+) 
        if (no == keypar.found) {
            printlog ("ERROR - NSPREPARE: Could not read " // l_key_bias \
                // " from header.", l_logfile, verbose+) 
            goto clean
        } else {
            biasvolt = real (keypar.value) 
            if (debug) print ("bias: " // biasvolt)
        }
        
        keypar(phu, l_key_arrayid, silent+)
        if (no == keypar.found) {
            printlog ("ERROR - NSPREPARE: Could not read " // l_key_arrayid \
                // " from header.", l_logfile, verbose+) 
            goto clean
        } else {
            arrayid = keypar.value
            if (debug) print ("arrayid: " // arrayid)
        }

        # Read array info from table

        if (debug) print ("reading arraytable")

        if (no == access (l_arraytable) ) {
            printlog ("ERROR - NSPREPARE: Array table " // l_arraytable \
                // " not found.", l_logfile, verbose+) 
            goto clean
        } else {

            tdelete (tmparray, go_ahead+, verify-, >& "dev$null")

            expr = "arrayid=='"// arrayid // \
                "'&&abs( bias - abs ( " // biasvolt // " )) < 0.1"
            if (debug) print (expr)
            if (debug) print (l_arraytable)
            tselect (intable = l_arraytable, outtable = tmparray, expr = expr)

            tinfo (table = tmparray, >& "dev$null")
            if (1 != tinfo.nrows) {
                printlog ("ERROR - NSPREPARE: No array calibration data \
                    for bias voltage of " // biasvolt, l_logfile, verbose+)
                printlog ("                   in " // l_arraytable, \
                    l_logfile, verbose+) 
                goto clean
            }

            tabpar (table = tmparray, column = "readnoise", row = 1, format-)
            if (tabpar.undef) {
                printlog ("ERROR - NSPREPARE: No readnoise \
                    for bias voltage of " // biasvolt, l_logfile, verbose+)
                printlog ("                   in " // l_arraytable, \
                    l_logfile, verbose+) 
                goto clean
            } else ronref = real (tabpar.value)

            ron = (ronref * sqrt (coadds) ) / (sqrt (lnrs) * sqrt (ndavgs) ) 

            tabpar (table = tmparray, column = "gain", row = 1, format-)
            if (tabpar.undef) {
                printlog ("ERROR - NSPREPARE: No gain \
                    for bias voltage of " // biasvolt, l_logfile, verbose+)
                printlog ("                   in " // l_arraytable, \
                    l_logfile, verbose+) 
                goto clean
            } else gain = real (tabpar.value)

            tabpar (table = tmparray, column = "well", row = 1, \
            format-)
            if (tabpar.undef) {
                printlog ("ERROR - NSPREPARE: No well defined in \
                    table " // l_arraytable, l_logfile, verbose+) 
                goto clean
            } else well = real (tabpar.value)

            tabpar (table = tmparray, column = "linearlimit", \
            row = 1, format-)
            if (tabpar.undef) {
                printlog ("ERROR - NSPREPARE: No linearlimit defined in \
                    table " // l_arraytable, l_logfile, verbose+) 
                goto clean
            } else linearlimit = real (tabpar.value)

            tabpar (table = tmparray, column = "nonlinearlimit", \
            row = 1, format-)
            if (tabpar.undef) {
                printlog ("ERROR - NSPREPARE: No nonlinearlimit defined \
                    in table " // l_arraytable, l_logfile, verbose+) 
                goto clean
            } else nonlinlimit = real (tabpar.value)

            if (l_fl_correct) {
                tabpar (table = tmparray, column = "coeff1", row = 1, format-)
                if (tabpar.undef) {
                    printlog ("ERROR - NSPREPARE: No coeff1 defined in \
                        table " // l_arraytable, l_logfile, verbose+) 
                    goto clean
                } else coeff1 = real (tabpar.value)

                tabpar (table = tmparray, column = "coeff2", row = 1, format-)
                if (tabpar.undef) {
                    printlog ("ERROR - NSPREPARE: No coeff2 defined in \
                        table " // l_arraytable, l_logfile, verbose+) 
                    goto clean
                } else coeff2 = real (tabpar.value)

                tabpar (table = tmparray, column = "coeff3", row = 1, format-)
                if (tabpar.undef) {
                    printlog ("ERROR - NSPREPARE: No coeff3 defined in \
                        table " // l_arraytable, l_logfile, verbose+) 
                    goto clean
                } else coeff3 = real (tabpar.value)
            }

        }

        # Read the config table, if available

        specsection = l_specsec
        pixscl = l_pixscale
        havemdf = no
        sections = "none"
        procmode = "LS"

        if (domdf && ! havedark) {
            if (debug) print ("reading mdf from table")
            if ("" != decker && "" != prism && "" != grating \
                && "" != camera && "" != arrayid) {

                tdelete (table = tmpmdfrow, go_ahead+, verify-, >& "dev$null")
                expr = "prism?='" // prismsearch // "'&&decker='" // decker \
                    // "'&&grating='" // grating // "'&&camera='" // camera \
                    // "'&&arrayid='" // arrayid // "'"
                tselect (intable = l_configtable, outtable = tmpmdfrow, \
                    expr = expr)
                if (debug) tinfo (table = tmpmdfrow)
                tinfo (table = tmpmdfrow, >>& tmpopen)
                havemdf = (0 != tinfo.nrows)
                if (havemdf) {
                    tabpar (table = tmpmdfrow, column = "mdf", row = 1, \
                        format-)
                    sections = tabpar.value
                    if (offsettype == OFF_CONF) {
                        tabpar (table = tmpmdfrow, column = "offsetsection",
                            row = 1, format-)
                        if (NONE == strlwr (tabpar.value)) {
                            offsettype = OFF_NONE
                        } else {
                            offsetsection = tabpar.value
                        }
                    }
                    tabpar (table = tmpmdfrow, column = "pixscale", \
                        row = 1, format-)
                    pixscl = real (tabpar.value)
                    tabpar (table = tmpmdfrow, column = "mode", \
                        row = 1, format-)
                    procmode = strupr (tabpar.value)
                }
            }
        }

        if (domdf && ! havedark && ! havemdf) {
            if (offsettype == OFF_CONF) {
                printlog ("ERROR - NSPREPARE: No data found using \
                    table " // l_configtable, l_logfile, verbose+) 
                printlog ("                   for prism " // prism // \
                    " / decker " // decker // " / grating " // grating \
                    // " / camera " // camera, l_logfile, verbose+)
                printlog ("                   but offsetsec='conf'.", \
                    l_logfile, verbose+)
                goto clean
            } else {
                printlog ("WARNING - NSPREPARE: No data found using \
                    table " // l_configtable, l_logfile, verbose+)
                printlog ("                     for prism " // prism // \
                    " / decker " // decker // " / grating " // grating \
                    // " / camera " // camera, l_logfile, verbose+)
            }
        }

        if (havedark && offsettype == OFF_CONF) {
            printlog ("NSPREPARE: Using full frame for dark offset.", \
                l_logfile, l_verbose)
            offsetsection = "[*,*]"
        }

        # specsec can be a section or an MDF file
        if (no == havemdf && no == havedark && "" != l_specsec) {
            if (debug) print ("checking specsec")
            havemdf = ! (substr (l_specsec, 1, 1) == "[")
            if (havemdf)
                sections = l_specsec
        }

        # Removed "/ gain" here - well in ADU
        if (debug) print ("well, coadds: " // well // " " // coadds)
        sat = int (well * coadds)
        if (l_fl_correct)
            linlimit = int (sat * linearlimit) 
        else
            linlimit = int (sat * nonlinlimit) 
        if (debug) print ("sat, linlimit " // sat // " " // linlimit)

        # Read date/instrument to control various fixes

        date = INDEF
        hasdate = no
        hselect (phu, "DATE-OBS", yes) | scan (date)
        if (isindef (date)) {
            printlog ("WARNING - NSPREPARE: Cannot read \
                observation date.  No date-specific fixes applied.", 
                l_logfile, verbose+) 
        } else {
            year = int (substr (date, 1, 4))
            month = int (substr (date, 6, 7))
            day = int (substr (date, 9, 10))
            hasdate = yes
        }

        keypar (phu, "INSTRUME", silent+)
        isgnirs = (keypar.found && strlwr (keypar.value) == "gnirs") 

        issouth = no
        if (isgnirs && hasdate) {issouth = year < 2008}

        extn = inimg // "[1]"

        if (no == already) {

            # Fix the data if number of non-destructive reads > 1 
            # for data before Nov 2001 or gnirs before 17 Mar 2004
            if (hasdate && lnrs != 1) {
                oldtimer = no
                if (isgnirs) {
                    oldtimer = (year < 2004 || \
                        (year == 2004 && \
                        (month < 3 || \
                        (month == 3 && day < 17))))
                } else {
                    oldtimer = (year < 2001 || \
                        (year == 2001 && month < 11))
                }
                if (oldtimer) {
                    if (isgnirs) {
                        printlog ("WARNING - NSPREPARE: Data taken prior \
                            to 17 Mar 2004", l_logfile, l_verbose)
                    } else {
                        printlog ("WARNING - NSPREPARE: Data taken prior \
                            to Nov 2001", l_logfile, l_verbose)
                    }
                    printlog ("                     Divided " // inimg \
                        // " by number of non-dest. reads: " // lnrs, \
                        l_logfile, l_verbose) 
                    imexpr ("a/"//lnrs, tmpinimage, extn, 
                        outtype = "real", >& "dev$null") 
                } else {
                    imcopy (extn, tmpinimage, verbose-) 
                }
            } else {
                imcopy (extn, tmpinimage, verbose-) 
            }

            # Remove offset
            if (offsettype == OFF_CONF || offsettype == OFF_SECN) {
                imstat (tmpinimage // offsetsection, \
                    fields="midpt,stddev", lower=INDEF, upper=INDEF, \
                    nclip=0, lsigma=INDEF, \
                    usigma=INDEF, binwidth=0.1, format-, cache-) \
                    | scan (midpt, stddev) 
                if (debug) print (midpt // ", " // stddev)
                imstat (tmpinimage // offsetsection, fields="midpt", \
                    lower=(midpt-3*stddev), upper=(midpt+3*stddev), \
                    nclip=0, lsigma=INDEF, usigma=INDEF, \
                    binwidth=0.1, format-, cache-) | scan (offset)
            }
            if (offsettype == OFF_CONF || offsettype == OFF_SECN \
                || offsettype == OFF_VAL) {
                
                if (debug) print ("subtracting offset of " // offset)
                imarith (tmpinimage, "-", offset, tmpinimage, title="", \
                    divzero=0, hparams="", pixtype="", calctype="", verbose-)
            }

            # Correct for non-linearity
            if (l_fl_correct) {

                if (debug) print (coeff1 // "," // coeff2 // "," // coeff3)
                if (debug) print (coadds)

                coeff1 = coeff1
                coeff2 = coeff2 / coadds
                coeff3 = coeff3 / (coadds * coadds)

                x = sat / 32767.0
                x = (((coeff3 * x) + coeff2) * x) + coeff1
                sat = sat * x
                x = linlimit / 32767.0
                x = (((coeff3 * x) + coeff2) * x) + coeff1
                linlimit = int (linlimit * x)
                if (debug)
                    print ("nonlin correcn " // x // " sat: " // sat // \
                        " linlimit: " // linlimit)
                if (debug) print (coeff1 // "," // coeff2 // "," // coeff3)

                irlincor (tmpinimage, tmpinimage, coeff1=coeff1, 
                    coeff2=coeff2, coeff3=coeff3)

            } else {
                printlog ("WARNING - NSPREPARE: Not correcting for \
                    non-linearity.", l_logfile, verbose+) 
            }

            # Check header for values
            keypar (phu, l_key_gain, silent+) 
            if (keypar.found) {
                gain = real (keypar.value) 
                printlog ("WARNING - NSPREPARE: Gain already set \
                    in image header.  Using a gain", l_logfile, verbose+) 
                printlog ("                    of " // gain \
                    // " electrons per ADU.", l_logfile, verbose+) 
            }

            keypar (phu, l_key_ron, silent+) 
            if (keypar.found) {
                ron = real (keypar.value) 
                printlog ("WARNING - NSPREPARE: Read noise \
                already set in image header.", l_logfile, verbose+) 
                printlog ("                    Using a read \
                noise of " // ron // " electrons.", l_logfile, verbose+) 
            }

            keypar (phu, l_key_sat, silent+) 
            if (keypar.found) {
                sat = real (keypar.value) 
                printlog ("WARNING - NSPREPARE: Saturation \
                level already set in image header.", l_logfile, verbose+) 
                printlog ("                    Using a \
                saturation level of " // sat // " ADU.", l_logfile, verbose+) 
            }

            if (l_fl_vardq) {

                # Create the variance frame, if it doesn't already exist
                # The variance frame is generated as:
                # var = (read noise/gain)^2 + max(data,0.0)/gain

                if (debug) print ("vardq")
                extn = inimg // "[" // l_var_ext // "]"

                if (no == imaccess (extn) ) {

                    varexpression = "((max(a,0.0))/" // gain \
                        // "+(" // ron // "/" // gain // ")**2)"
                    imexpr (varexpression, tmpvar, tmpinimage, \
                        outtype = "real", verbose-) 

                } else {
                    printlog ("WARNING - NSPREPARE: Variance frame already \
                        exists for " // inimg // ".", l_logfile, verbose+) 
                    printlog ("                     New variance frame \
                        not created.", l_logfile, verbose+) 
                    imcopy (extn, tmpvar, verbose-) 
                }

                # Create the DQ frame, if it doesn't already exist
                # The preliminary DQ frame is constructed by using 
                # the bad pixel mask to set bad pixels to 1, pixels 
                # in the non-linear regime to 2, and saturated pixels to 4.

                extn = inimg // "[" // l_dq_ext // "]"

                if (no == imaccess (extn)) {
                    if (no == l_fl_saturated && no == l_fl_nonlinear) {
                        tempexpression = "0 * a"
                    } else if (l_fl_saturated && l_fl_nonlinear) {
                        tempexpression = "(a>" // sat // ") ? 4 : \
                            ((a>" // linlimit // ") ? 2 : 0)"
                    } else if (l_fl_saturated) {
                        tempexpression = "(a>" // sat // ") ? 4 : 0"
                    } else if (l_fl_nonlinear) {
                        tempexpression = "(a>" // linlimit // ") ? 2 : 0"
                    }
                    imexpr (tempexpression, tmpsci, tmpinimage, \
                        outtype = "short", verbose-) 

                    # If there's no BPM, then just keep the saturated 
                    # pixels
                    if (l_bpm == "none") {
                        addmasks (tmpsci // ".fits", tmpdq // ".fits", "im1") 
                    } else {
                        addmasks (tmpsci // ".fits," // l_bpm, \
                            tmpdq // ".fits", "im1 || im2") 
                    }
                    imdelete (tmpsci, verify-, >& "dev$null") 
                } else {
                    printlog ("WARNING - NSPREPARE: Data quality frame \
                        already exists for " // inimg // ".", \
                        l_logfile, verbose+) 
                    printlog ("                     New DQ frame not \
                        created.", l_logfile, verbose+) 
                    imcopy (extn, tmpdq // ".fits", verbose-) 
                }

                # Pack up the results and clean up

                wmef (tmpinimage // "," // tmpvar // ".fits," \
                    // tmpdq // ".fits", outimg, extnames = l_sci_ext \
                    // "," // l_var_ext // "," // l_dq_ext, \
                    phu = inimg // ".fits[0]", verbose-) 
                # don't hide error messages - this may be a source
                # of errors with old iraf versions and we're not
                # going to know if the error message is hidden

                if (wmef.status != 0) {
                    printlog ("ERROR - NSPREPARE: Could not write final \
                        MEF file (WMEF).", l_logfile, verbose+) 
                    goto clean
                }

                gemhedit (outimg // "[" // l_var_ext // "]", "EXTVER", 1, 
                    "", delete-) 
                gemhedit (outimg // "[" // l_dq_ext // "]", "EXTVER", 1,
                    "", delete-) 
                imdelete (tmpvar // "," // tmpdq, verify-, >& "dev$null") 

            } else {

                # No var/dq

                wmef (tmpinimage, outimg, extnames = l_sci_ext, 
                    phu = inimg // ".fits[0]", verbose-, >& "dev$null") 
                if (wmef.status != 0) {
                    printlog ("ERROR - NSPREPARE: Could not write final \
                        MEF file (WMEF).", l_logfile, verbose+) 
                    goto clean
                }
            }

            imdelete (tmpinimage, verify-, >& "dev$null") 

            # Add the information

            phu = outimg // "[0]"

            if (no == havedark) {

                if (havemdf) {
                    if (debug) print ("adding mdf")
                    delete (tmpmdffile, verify-, >& "dev$null")
                    tmpmdffile = mktemp ("tmpmdffile")
                    gemextn (sections, check="exists,mef,table", \
                        process="expand", index="1", extname="", \
                        extversion="", ikparams="", omit="exten", \
                        replace="", outfile=tmpmdffile, logfile="", \
                        glogpars="", verbose-)
                    if (0 != gemextn.fail_count || 1 != gemextn.count) {
                        printlog ("ERROR - NSPREPARE: MDF file " \
                            // sections // " not found, or incorrect \
                            format.", l_logfile, verbose+) 
                        goto clean
                    }
                    scanin2 = tmpmdffile
                    junk = fscan (scanin2, sections)

                    # Auto-shift
                    if (correlate && firstmdf) {
                        oldmdf = sections
                        imdelete (tmpsections, verify-, >& "dev$null")
                        gemextn (sections, check="", process="none", \
                            index="", extname="", extversion="", \
                            ikparams="", omit="path,exten,kernel,index", \
                            replace="", outfile="STDOUT", \
                            logfile="dev$null", glogpars="", verbose-) | \
                            scan (tmpsections)
                        tmpsections = mktemp (tmpsections // "-")
                        if (debug) print (tmpsections)
                        nsoffset (outimg, sections, outimages=tmpsections, \
                            outprefix="", axis=3-dispaxis, pixscale=pixscl,
                            fl_apply+, fl_all-, 
                            fl_project=(no==("XD"==procmode)), \
                            fl_inter=l_fl_inter, logfile=l_logfile, \
                            verbose=l_verbose)
                        sections = tmpsections
                        if (0 != nsoffset.status) {
                            printlog ("ERROR - NSPREPARE: Auto-shift \
                                failed for MDF.", l_logfile, verbose+) 
                            goto clean
                        }
                        if (2 == dispaxis) {
                            l_shiftx = nsoffset.shift
                                gemhedit (phu, "MDF_XSHF", l_shiftx, \
                                    "MDF X shift (subtracted from x_ccd)", \
                                    delete-)
                        } else {
                            l_shifty = nsoffset.shift
                            gemhedit (phu, "MDF_YSHF", l_shifty, \
                                "MDF Y shift (subtracted from y_ccd)", \
                                delete-)
                        }
                        firstmdf = no
                    } else if (correlate) {
                        if (sections == oldmdf) {
                            sections = tmpsections
                        } else {
                            printlog ("ERROR - NSPREPARE: MDF cannot \
                                change if using auto-shifting.", \
                                l_logfile, verbose+)
                            goto clean
                        }
                        if (2 == dispaxis) {
                            gemhedit (phu, "MDF_XSHF", l_shiftx, \
                                "MDF X shift (subtracted from x_ccd)", \
                                delete-)
                        } else {
                            gemhedit (phu, "MDF_YSHF", l_shifty, \
                                "MDF Y shift (subtracted from y_ccd)", \
                                delete-)
                        }
                    } else {
                        oldmdf = sections
                    }

                    fxinsert (sections, outimg // ".fits[0]", groups="1", \
                        verbose-)
                    thedit (outimg // ".fits[1]", "EXTNAME", "MDF", \
                        delete-, show+, >& "dev$null")

                    if (no == correlate && no == isindef (l_shiftx) && \
                        0.0 != l_shiftx) {

                        tcalc (outimg // ".fits[MDF]", outcol="x_ccd", \
                            equals="x_ccd-" // l_shiftx, datatype="real")
                        printlog ("NSPREPARE: Shifting MDF in X by " \
                            // l_shiftx // " pixels (subtracted from \
                            x_ccd)", l_logfile, l_verbose) 
                        gemhedit (phu, "MDF_XSHF", l_shiftx, \
                            "MDF X shift (subtracted from x_ccd)", \
                            delete-)
                    }
                    if (no == correlate && no == isindef (l_shifty) && \
                        0.0 != l_shifty) {

                        tcalc (outimg // ".fits[MDF]", outcol="y_ccd", \
                            equals="y_ccd-" // l_shifty, datatype="real")
                        printlog ("NSPREPARE: Shifting MDF in Y by " \
                            // l_shifty // " pixels (subtracted from \
                            y_ccd)", l_logfile, l_verbose) 
                        gemhedit (phu, "MDF_YSHF", l_shifty, \
                            "MDF Y shift (subtracted from y_ccd)", \
                            delete-)
                    }

                    gemhedit (phu, "MDF_FILE", oldmdf, "MDF source", delete-)

                } else {
                    printlog ("Using section " // specsection, \
                        l_logfile, l_verbose)
                    gemhedit (phu, l_key_section // "1", specsection,
                        "", delete-) 
                }
            }

        } else {
            phu = inimg // "[0]"
        }

        if (first) {
            # Start output

            printlog (" ", l_logfile, l_verbose) 
            printlog ("  n      input file -->      output file", 
                l_logfile, l_verbose) 
            printlog ("                 filter     focal plane       \
                input BPM   RON  gain     sat", l_logfile, l_verbose) 
            printlog ("                 offset   MDF", \
                l_logfile, l_verbose) 
            printlog (" ", l_logfile, l_verbose) 
            first = no
        }

        # Fix WCS problems

        wcsok = yes

        # !!!!!!!!!!!!!!!!!!!! needs fixing once we know this works!
        # It would really help if anyone knew what this was meant to do!
        fixshift = (isgnirs && (year < 2004 || \
            (year == 2004 && (month < 7 || \
            month == 7 && day < 16))) \
            && ! already)

        # This date should probably be earlier
        fixwcs = (isgnirs && (year < 2004 || \
            (year == 2004 && (month < 3 || \
            month == 3 && day < 12))) \
            && ! already)

        if (fixshift || fixwcs || l_fl_checkwcs || l_fl_forcewcs) {

            keypar (phu, "RA"); ra = real (keypar.value)
            keypar (phu, "DEC"); dec = real (keypar.value)
            keypar (phu, "RAOFFSET"); raoff = real (keypar.value)
            keypar (phu, "DECOFFSE"); decoff = real (keypar.value)
            keypar (phu, "PA"); pa = real (keypar.value)

            if (debug) {
                print "RA: " // ra // "; DEC: " // dec
                print "RAOFFSET: " // raoff // "; DECOFF: " // decoff
                print "PA: " // pa // "; PIXSCALE: " // pixscl
            }

	    decoff = decoff / 3600.0 # arcsec
	    if (isgnirs && (year*100+month)<201403) {
	        decoff = -decoff
	    }
            raoff = raoff / (3600.0 * cos ((dec+decoff)*torad))
            ranod = ra + raoff
            decnod = dec + decoff

            xsize = 0
            ysize = 0

            # First data frame
            extn = INDEF
            gemextn (inimg, proc="expand", check="image", index="1-", \
                extname="", extver="", ikparams="", omit="", replace="", \
                outfile="STDOUT", logfile="dev$null", glogpars="", \
                verbose-) | scan (extn)
            keypar (extn, "i_naxis1")
            if (keypar.found)
                xsize = int (keypar.value)
            keypar (extn, "i_naxis2")
            if (keypar.found)
                ysize = int (keypar.value)

        }

        if (fixwcs || l_fl_checkwcs || l_fl_forcewcs) {

            #pa = pa - 90
            pa = pa * torad
            pix = pixscl / 3600.0

            cospa = cos (pa)
            sinpa = sin (pa)

            cd11 = -pix*sinpa
            cd12 = pix*cospa
            cd21 = -pix*cospa
            cd22 = -pix*sinpa

            if (debug) {
                printf ("%11.5g %11.5g\n", cd11, cd21) | scan (l_struct)
                print (l_struct)
                printf ("%11.5g %11.5g\n", cd12, cd22) | scan (l_struct)
                print (l_struct)
            }
        }

        if (l_fl_checkwcs) {

            if (no == fixshift) {

                # Check offsets

                keypar (phu, "CRVAL1", silent+)
                if (no == keypar.found) {
                    printlog ("WARNING - NSPREPARE: Missing CRVAL1.", \
                        l_logfile, verbose+)
                    cval1 = INDEF
                } else {
                    cval1 = real (keypar.value)
                    abserr = abs (cval1 - ranod) / pix
                    if (abserr > xsize) {
                        printf ("%11.5g %11.5g\n", ranod, cval1) \
                            | scan (l_struct)
                        printlog ("WARNING - NSPREPARE: RA and CRVAL1 \
                            differ: " // l_struct, l_logfile, verbose+)
                        wcsok = no
                    }
                }

                keypar (phu, "CRVAL2", silent+)
                if (no == keypar.found) {
                    printlog ("WARNING - NSPREPARE: Missing CRVAL2.", \
                        l_logfile, verbose+)
                    cval2 = INDEF
                } else {
                    cval2 = real (keypar.value)
                    abserr = abs (cval2 - decnod) / pix
                    if (abserr > ysize) {
                        printf ("%11.5g %11.5g\n", decnod, cval2) \
                            | scan (l_struct)
                        printlog ("WARNING - NSPREPARE: DEC and CRVAL2 \
                            differ: " // l_struct, l_logfile, verbose+)
                        wcsok = no
                    }
                }

                # Check that object always at same place on detector

                # Take ra + dec and invert to find pixel coords

                delete (tmpcoordsin, verify-, >& "dev$null")
                tmpcoordsin = mktemp ("tmpcoordsin")
                printf ("%f %f\n", (ra-raoff), (dec-decoff), >> tmpcoordsin)
                # First data frame with inheritance
                extn = INDEF
                gemextn (inimg, proc="expand", check="image", index="1-", \
                    extname="", extver="", ikparams="inherit", omit="", \
                    replace="", outfile="STDOUT", logfile="dev$null", \
                    glogpars="", verbose-) | scan (extn)
                wcsctran (tmpcoordsin, "STDOUT", extn, \
                    "world", "logical", columns="1 2", units="", \
                    formats="", min_sigdigits=7, verbose-) |
                    scan (xpix, ypix)

                if (no == (isindef (xobj) || isindef (yobj))) {
                    if ((xobj-xpix)**2 + (yobj-ypix)**2 > rmax) {
                        printf ("(%6.1f,%6.1f) -> (%6.1f,%6.1f)\n", \
                            xobj, yobj, xpix, ypix) | scan (l_struct)
                        printlog ("WARNING - NSPREPARE: Object shift: " \
                            // l_struct, l_logfile, verbose+)
                        wcsok = no
                    }
                }
                xobj = xpix
                yobj = ypix
            }

            if (no == fixwcs) {

                # Check matrix values

                keypar (phu, "CD1_1", silent+)
                if (no == keypar.found) {
                    printlog ("WARNING - NSPREPARE: Missing CD1_1.", \
                        l_logfile, verbose+)
                } else {
                    cd11obs = real (keypar.value)
                    abserr = abs (cd11obs - cd11) / pix
                    if (abserr > errlim) {
                        printf ("%11.5g %11.5g\n", cd11, cd11obs) \
                            | scan (l_struct)
                        printlog ("WARNING - NSPREPARE: CD1_1 error \
                            (calc, obs):" // l_struct, l_logfile, verbose+)
                        wcsok = no
                    }
                }

                keypar (phu, "CD1_2", silent+)
                if (no == keypar.found) {
                    printlog ("WARNING - NSPREPARE: Missing CD1_2.", \
                        l_logfile, verbose+)
                } else {
                    cd12obs = real (keypar.value)
                    abserr = abs (cd12obs - cd12) / pix
                    if (abserr > errlim) {
                        printf ("%11.5g %11.5g\n", cd12, cd12obs) \
                            | scan (l_struct)
                        printlog ("WARNING - NSPREPARE: CD1_2 error \
                            (calc, obs):" // l_struct, l_logfile, verbose+)
                        wcsok = no
                    }
                }

                keypar (phu, "CD2_1", silent+)
                if (no == keypar.found) {
                    printlog ("WARNING - NSPREPARE: Missing CD2_1.", \
                        l_logfile, verbose+)
                } else {
                    cd21obs = real (keypar.value)
                    abserr = abs (cd21obs - cd21) / pix
                    if (abserr > errlim) {
                        printf ("%11.5g %11.5g\n", cd21, cd21obs) \
                            | scan (l_struct)
                        printlog ("WARNING - NSPREPARE: CD2_1 error \
                            (calc, obs):" // l_struct, l_logfile, verbose+)
                        wcsok = no
                    }
                }

                keypar (phu, "CD2_2", silent+)
                if (no == keypar.found) {
                    printlog ("WARNING - NSPREPARE: Missing CD2_2.", \
                        l_logfile, verbose+)
                } else {
                    cd22obs = real (keypar.value)
                    abserr = abs (cd22obs - cd22) / pix
                    if (abserr > errlim) {
                        printf ("%11.5g %11.5g\n", cd22, cd22obs) \
                            | scan (l_struct)
                        printlog ("WARNING - NSPREPARE: CD2_2 error \
                            (calc, obs):" // l_struct, l_logfile, verbose+)
                        wcsok = no
                    }
                }

                jacob = sqrt ( abs(cd11obs*cd22obs - \
                    cd12obs*cd21obs)) * 3600
                abserr = abs (pixscl - jacob) / (tiny + abs (pixscl))
                if (abserr > errlim) {
                    printf ("%11.5g %11.5g\n", pixscl, jacob) \
                        | scan (l_struct)
                    printlog ("WARNING - NSPREPARE: Pixel scale error \
                        (calc, obs):" // l_struct, l_logfile, verbose+)
                    wcsok = no
                }
            }
        }

        if (l_fl_checkwcs && ((! fixshift) || (! fixwcs))) {
            if (wcsok) {
                printlog ("NSPREPARE: Header checks OK.", l_logfile, \
                    l_verbose)
            }
        }

        if (fixshift || l_fl_forcewcs) {

            if (l_fl_checkwcs && fixshift) {
                printlog ("WARNING - NSPREPARE: Cannot check pointing.", \
                    l_logfile, l_verbose)
            }
            printlog ("WARNING - NSPREPARE: Adding telescope \
                pointing.", l_logfile, l_verbose)

            gemhedit (phu, "CRVAL1", ranod, "", delete-)
            gemhedit (phu, "CRVAL2", decnod, "", delete-)

        }

        if (fixwcs || l_fl_forcewcs) {

            if (l_fl_checkwcs && fixwcs) {
                printlog ("WARNING - NSPREPARE: Cannot check WCS \
                    matrix.", l_logfile, l_verbose)
            }
            printlog ("WARNING - NSPREPARE: Supplying WCS \
                matrix.", l_logfile, l_verbose)

            gemhedit (phu, "CD1_1", cd11, "", delete-)
            gemhedit (phu, "CD1_2", cd12, "", delete-)
            gemhedit (phu, "CD2_1", cd21, "", delete-)
            gemhedit (phu, "CD2_2", cd22, "", delete-)
        }

        # Fix up the headers

        if (no == already) {

            if (debug) print ("headers")

            if (l_fl_vardq) {
                delete (tmphead, verify-, >& "dev$null") 
                printf ("%-8s=                    %d / %-s\n", \
                    l_key_dispaxis, dispaxis, \
                    "Dispersion axis (along columns)", >> tmphead) 
                printf ("%-8s= %20.5f / %-s\n", "PIXSCALE", \
                    pixscl, "Pixel scale in arcsec/pixel", >> tmphead) 
                mkheader (outimg // "[" // l_var_ext // "]", tmphead, \
                    append+, verbose-) 
                mkheader (outimg // "[" // l_dq_ext // "]", tmphead, \
                    append+, verbose-) 
            }

            gemhedit (outimg // "[" // l_sci_ext // "]", "EXTVER", 1, 
                "", delete-) 

            delete (tmphead, verify-, >& "dev$null") 
            printf ("%-8s=                    %d / %-s\n", \
                l_key_dispaxis, dispaxis, \
                "Dispersion axis (along columns)", >> tmphead) 
            printf ("%-8s= %20.5f / %-s\n", "PIXSCALE", \
                pixscl, "Pixel scale in arcsec/pixel", >> tmphead) 
            mkheader (outimg // "[" // l_sci_ext // "]", tmphead, \
                append+, verbose-) 
            delete (tmphead, verify-, >& "dev$null") 

            if (l_fl_correct) {
                printf ("%7.5f %7.5f %7.5f\n", coeff1, coeff2, coeff3) \
                    | scan (l_struct)
                gemhedit (phu, "NONLINCR", l_struct, 
                    "Non-linear correction applied", delete-)

                if (debug)
                    print (coeff1 // "," // coeff2 // "," // coeff3)
            }

            if (offsettype != OFF_NONE) {
                if (offsettype == OFF_SECN || offsettype == OFF_CONF) {
                    printf ("%-8s= %20s / %-s\n", "OFFSETSE", \
                        offsetsection, "Offset section", >> tmphead) 
                }
                printf ("%-8s= %20.2f / %-s\n", "OFFSET", \
                    offset, "Offset value (subtracted from data)", \
                    >> tmphead) 
                mkheader (phu, tmphead, append+, verbose-) 
                delete (tmphead, verify-, >& "dev$null") 
            }

            gemhedit (phu, "WMEF", "", "", delete+)
            gemhedit (phu, "GEM-TLM", "", "", delete+)

            # Add comment to NEXTEND
            keypar (phu, "NEXTEND", silent+) 
            if (keypar.found) {
                gemhedit (phu, "NEXTEND", "", "", delete+)
                delete (tmphead, verify-, >& "dev$null") 
                printf ("%-8s= %20.0f / %-s\n", "NEXTEND", \
                    real (keypar.value), "Number of extensions", \
                    >> tmphead) 

                if (debug) print (phu)
                if (debug) type (tmphead)
                mkheader (phu, tmphead, append+, verbose-) 
            }

#            delete (tmphead, verify-, >& "dev$null") 
#            printf ("%-8s= %20s / %-s\n", l_key_mode, procmode, \
#                "Reduction mode", >> tmphead) 
#            mkheader (phu, tmphead, append+, verbose-) 

            gemhedit (phu, l_key_mode, procmode, "Reduction mode", delete-)

            delete (tmphead, verify-, >& "dev$null") 

            printf ("%-8s= %20.2f / %-s\n", l_key_ron, ron, \
                "Estimated read noise (electrons)", >> tmphead) 

            if (debug) print (phu)
            if (debug) type (tmphead)
            mkheader (phu, tmphead, append+, verbose-) 
            delete (tmphead, verify-, >& "dev$null") 

            printf ("%-8s= %20.2f / %-s\n", l_key_gain, gain, \
                "Gain (electrons/ADU)", >> tmphead) 

            if (debug) print (phu)
            if (debug) type (tmphead)
            mkheader (phu, tmphead, append+, verbose-) 
            delete (tmphead, verify-, >& "dev$null") 

            printf ("%-8s= %20.3f / %-s\n", "BIASVOLT", biasvolt, \
                "Array bias voltage (V)", >> tmphead) 

            if (debug) print (phu)
            if (debug) type (tmphead)
            mkheader (phu, tmphead, append+, verbose-) 
            delete (tmphead, verify-, >& "dev$null") 

            printf ("%-8s= %20.0f / %-s\n", l_key_sat, sat, \
                "Saturation level in ADU", >> tmphead) 

            if (debug) print (phu)
            if (debug) type (tmphead)
            mkheader (phu, tmphead, append+, verbose-) 
            delete (tmphead, verify-, >& "dev$null") 

            printf ("%-8s= %20.0f / %-s\n", l_key_nonlinear, linlimit, \
                "Non-linear regime in ADU", >> tmphead) 

            if (debug) print (phu)
            if (debug) type (tmphead)
            mkheader (phu, tmphead, append+, verbose-) 
            delete (tmphead, verify-, >& "dev$null") 

            if (coadds != 1) {
                keypar (phu, l_key_exptime, silent+) 
                exptime = real (keypar.value) 
                printf ("%-8s= %20.5f / %-s\n", "COADDEXP", exptime, \
                    "Exposure time (s) for each frame", >> tmphead) 

                if (debug) print (phu)
                if (debug) type (tmphead)
                mkheader (phu, tmphead, append+, verbose-) 
                delete (tmphead, verify-, >& "dev$null") 

                exptime = exptime * coadds
                gemhedit (phu, l_key_exptime, "", "", delete+)
                printf ("%-8s= %20.5f / %-s\n", l_key_exptime, exptime, \
                    "Exposure time (s) for sum of all coadds", >> tmphead) 
                if (debug) {
                    printlog ("Scaling exposure time for " // outimg \
                        // ". Coadded exposure time is " // exptime \
                        // " s.", l_logfile, l_verbose) 
                }

                if (debug) print (phu)
                if (debug) type (tmphead)
                mkheader (phu, tmphead, append+, verbose-) 
                delete (tmphead, verify-, >& "dev$null") 

            }

            # TODO - what should this be for GNIRS?

            # Parse the individual filter entries into a single 
            # header keyword

            keypar (phu, l_key_filter, silent+) 
            if (no == keypar.found) {
                keypar (phu, "FILTER1", silent+) 
                filter1 = keypar.value
                keypar (phu, "FILTER2", silent+) 
                filter2 = keypar.value
                keypar (phu, "FILTER3", silent+) 
                filter3 = keypar.value
                if (filter1 == "0" || filter1 == "OPEN" \
                    || filter1 == "open" || filter1 == "Open" \
                    || filter1 == "INVALID") filter1 = ""
                if (filter2 == "0" || filter2 == "OPEN" \
                    || filter2 == "open" || filter2 == "Open" \
                    || filter2 == "INVALID") filter2 = ""
                if (filter3 == "0" || filter3 == "OPEN" \
                    || filter3 == "open" || filter3 == "Open" \
                    || filter3 == "INVALID") filter3 = ""
                if (substr (filter3, 1, 3) == "pup") {
                    pupil = filter3
                    if (substr (pupil, strlen (pupil) -5, \
                        strlen (pupil) -4) == "_G") {
                        
                        pupil = substr (pupil, 1, (strlen (pupil) -6) )
                    }
                    filter3 = ""
                }

                # Only strip the Gemini filter numbers before making the 
                # combined filter keyword if the data is from GNIRS-S. For 
                # data from GNIRS-N (which have an arrayid), keep the Gemini 
                # filter numbers.

                if (arrayid == "UNKNOWN") {
                    if (substr (filter1, strlen (filter1) -5, \
                        strlen (filter1) -4) == "_G") 
                        filter1 = substr (filter1, 1, (strlen (filter1) -6) ) 

                    if (substr (filter2, strlen (filter2) -5, \
                        strlen (filter2) -4) == "_G") 
                        filter2 = substr (filter2, 1, (strlen (filter2) -6) ) 

                    if (substr (filter3, strlen (filter3) -5, \
                        strlen (filter3) -4) == "_G") 
                        filter3 = substr (filter3, 1, (strlen (filter3) -6) ) 
                }

                if (filter1 == "blank" || filter2 == "blank" \
                    || filter3 == "blank") {
                    
                    print ("blank", > tmpfile2) 
                    pupil = "none"
                } else {
                    printf ("%s\n%s\n%s\n", filter1, filter2, filter3) | \
                        sort ("STDIN", col = 0, ignore+, num-, rev-, \
                        > tmpfile2) 
                }

                filter = ""
                scanin2 = tmpfile2
                while (fscan (scanin2, temp) != EOF)
                    filter = filter + temp

                printf ("%-8s= \'%-18s\' / %-s\n", l_key_filter, \
                    filter, "Filter name combined from all 3 wheels", \
                    >> tmphead) 

                if (debug) print (phu)
                if (debug) type (tmphead)
                mkheader (phu, tmphead, append+, verbose-) 
                delete (tmphead, verify-, >& "dev$null") 

                printf ("%-8s= \'%-18s\' / %-s\n", "PUPILMSK", pupil, \
                    "Name of pupil mask", >> tmphead) 

                if (debug) print (phu)
                if (debug) type (tmphead)
                mkheader (phu, tmphead, append+, verbose-) 
                delete (tmphead, verify-, >& "dev$null") 

                scanin2 = ""; delete (tmpfile2, verify-, >& "dev$null") 

            } else {
                printlog ("WARNING - NSPREPARE: filter keyword " \
                    // l_key_filter // " already exists.", l_logfile, verbose+)
                printlog ("                     GNIRS filters not \
                    parsed.", l_logfile, verbose+) 
                filter = ""
            }

            printf ("%-8s= %20.5f / %-s\n", "PIXSCALE", \
                pixscl, "Pixel scale in arcsec/pixel", >> tmphead) 

            if (debug) print (phu)
            if (debug) type (tmphead)
            mkheader (phu, tmphead, append+, verbose-) 
            delete (tmphead, verify-, >& "dev$null") 

            printf ("%-8s=                    1 / %-s\n", "NSCIEXT", \
                "Number of science extensions", >> tmphead) 

            if (debug) print (phu)
            if (debug) type (tmphead)
            mkheader (phu, tmphead, append+, verbose-) 
            delete (tmphead, verify-, >& "dev$null") 

            printf ("%-8s=                    %d / %-s\n", \
                l_key_dispaxis, dispaxis, \
                "Dispersion axis (along columns)", >> tmphead) 

            if (debug) print (phu)
            if (debug) type (tmphead)
            mkheader (phu, tmphead, append+, verbose-) 
            delete (tmphead, verify-, >& "dev$null") 

            # Time stamps
            gemdate ()
            printf ("%-8s= \'%-18s\' / %-s\n", "GEM-TLM", gemdate.outdate, \
                "UT Last modification with GEMINI", >> tmphead) 

            if (debug) print (phu)
            if (debug) type (tmphead)
            mkheader (phu, tmphead, append+, verbose-) 
            delete (tmphead, verify-, >& "dev$null") 

            printf ("%-8s= \'%-18s\' / %-s\n", "PREPARE", gemdate.outdate, \
                "UT Time stamp for NSPREPARE", >> tmphead) 

            if (debug) print (phu)
            if (debug) type (tmphead)
            mkheader (phu, tmphead, append+, verbose-) 
            delete (tmphead, verify-, >& "dev$null") 

            printf ("%-8s= \'%-18s\' / %-s\n", "BPMFILE", l_bpm, \
                "Input bad pixel mask file", >> tmphead) 

            if (debug) print (phu)
            if (debug) type (tmphead)
            mkheader (phu, tmphead, append+, verbose-) 
            delete (tmphead, verify-, >& "dev$null") 

            # print output to logfile
            keypar (phu, "FPMASK", silent+) 
            if (substr (keypar.value, strlen (keypar.value) -5, \
                strlen (keypar.value) -4) == "_G") {
                
                fpmask = substr (keypar.value, 1, (strlen (keypar.value) -6) ) 
            } else {
                fpmask = keypar.value
            }

            printf ("%3.0d %15s --> %16s \n", count, inimg, outimg) \
                | scan (l_struct) 
            printlog (l_struct, l_logfile, l_verbose) 
            printf ("      %17s %15s %15s %5.1f %5.1f %7.0d \n", \
                filter, fpmask, l_bpm, ron, gain, sat) | scan (l_struct) 
            printlog (l_struct, l_logfile, l_verbose) 
            if (havemdf) {
                printf ("        %-8s=%6.2f %40s\n", offsetsection, \
                    offset, oldmdf) | scan (l_struct) 
            } else {
                printf ("        %-8s=%6.2f %40s\n", offsetsection, \
                    offset, sections) | scan (l_struct) 
            }
            printlog (l_struct, l_logfile, l_verbose) 

        } else {

            printf ("%3.0d %15s already processed\n", count, inimg) \
                | scan (l_struct) 
            printlog (l_struct, l_logfile, l_verbose) 

        }

        # Detect cosmic rays and radiation events
        # Currently experimental, with fixed params
        if (l_fl_cravg) {
            printlog ("NSPREPARE: Detecting cosmic rays and radiation \
                events", l_logfile, l_verbose)
            craverage (outimg // "[" // l_sci_ext // ",1]", "",
                crmask=outimg // "[" // l_dq_ext // ",1]", \
                average="", sigma="", navg=9, nrej=48, nbkg=11, nsig=25, \
                var0=(ron**2), var1=gain, crval=8, lcrsig=5, hcrsig=10, \
                crgrow=1, objval=0, lobjsig=10, hobjsig=10, objgrow=0)
        }

        # Grow bad pixel detections
        if (crradius > 0) {
            printlog ("NSPREPARE: Extending bad pixel mask", \
                l_logfile, l_verbose)
            crgrow (outimg // "[" // l_dq_ext // ",1]", \
                outimg // "[" // l_dq_ext // ",1]", radius=crradius,
                inval=INDEF, outval=INDEF)
        }

    }

    status = 0

clean:
    if (status == 0) {
        printlog (" ", l_logfile, l_verbose) 
        printlog ("NSPREPARE exit status: good.", l_logfile, l_verbose) 
    }
    printlog ("-----------------------------------------------------\
        -----------------------", l_logfile, l_verbose) 

    scanin1 = "" 
    scanin2 = "" 
    delete (tmpexpand, verify-, >& "dev$null")
    delete (tmpfile, verify-, >& "dev$null")
    delete (tmproot, verify-, >& "dev$null")
    delete (tmpout, verify-, >& "dev$null")
    delete (tmpsummary, verify-, >& "dev$null")
    delete (tmpfile2, verify-, >& "dev$null")
    delete (tmphead, verify-, >& "dev$null")
    delete (tmpopen, verify-, >& "dev$null")
    delete (tmpinimage, verify-, >& "dev$null")
    delete (tmpmdffile, verify-, >& "dev$null")
    delete (tmpcoordsin, verify-, >& "dev$null")
    delete (tmpnotset, verify-, >& "dev$null")
    imdelete (tmpsections, verify-, >& "dev$null")
    tdelete (tmpmdfrow, go_ahead+, verify-, >& "dev$null")
    tdelete (tmparray, go_ahead+, verify-, >& "dev$null")
    imdelete (tmpsci // "," // tmpvar // "," // tmpdq, verify-, >& "dev$null") 

end

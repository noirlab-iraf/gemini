# Copyright(c) 2006-2012 Association of Universities for Research in Astronomy, Inc.

# NFPREPARE - Pre-process NIFS data
#
# Original author: Tracy Beck

procedure nfprepare (inimages) 

# The following steps are taken:
# - information is read from config.fits and array.fits, based on headers
# - saturation level and linear limits are calculated
# - data values are corrected by number of lnrs if this isn't already done
# - the reference pixel values are subtracted
# - the data are corected for non-linearity
# - variance data are generated from science data (photon stats)
# - data quality data are generated from data levels and bpm
# - mdf is attached


char    inimages    {prompt = "Input NIFS image(s)"}                # OLDP-1-primary-single-prefix=n
char    rawpath     {"", prompt = "Path for input raw images"}      # OLDP-4
char    outimages   {"", prompt = "Output image(s)"}                # OLDP-1-output
char    outprefix   {"n", prompt = "Prefix for output image(s)"}    # OLDP-4
char    bpm         {"", prompt = "Bad pixel mask file"}            # OLDP-1

bool    fl_vardq    {yes, prompt = "Create variance and data quality frames?"}  # OLDP-2
bool    fl_cravg    {no, prompt = "Attempt to flag cosmic rays and radiation events?"}    # OLDP-2
int     crradius    {0, prompt = "Grow radius for bad pixels"}      # OLDP-2
bool    fl_dark_mdf {no, prompt = "Force the attachment of MDFs to dark frames?"}    # OLDP-3
bool    fl_subtract {yes, prompt="Subtract off reference pixels?"} 
bool    fl_correct  {no, prompt = "Correct for non-linearity in the data?"}    # OLDP-3
bool    fl_saturated    {yes, prompt = "Flag saturated pixels in DQ?"}  # OLDP-3
bool    fl_nonlinear    {no, prompt = "Flag non-linear pixels in DQ?"} # OLDP-3

char    arraytable  {"nifs$data/nifsarray.fits", prompt = "Table of array settings"}    # OLDP-3
char    configtable {"nifs$data/nifsconfig.fits", prompt = "Table of configuration settings"}    # OLDP-3
char    specsec     {"[*,*]", prompt = "Section or MDF file if config table missing"}    # OLDP-3
real    pixscale    {0.043, prompt = "Pixel scale if if config table missing"}  # OLDP-3
char    shiftimage  {"", prompt = "Image to copy shift from"}       # OLDP-3
real    shiftx      {0, prompt = "Shift in X for MDF (subtracted from x_ccd)"}  # OLDP-3
real    shifty      {0, prompt = "Shift in Y for MDF (subtracted from y_ccd)"}  # OLDP-3
char    obstype     {"FLAT", prompt = "Observation type for cross-correlation"} # OLDP-3
bool    fl_inter    {no, prompt = "Fit offset interactively?"}      # OLDP-3

char    logfile     {"", prompt = "Logfile"}                        # OLDP-1
bool    verbose     {yes, prompt = "Verbose"}                       # OLDP-2
int     status      {0, prompt = "Exit status (0=good)"}            # OLDP-4

struct  *scanin1    {"", prompt = "Internal use only"}              # OLDP-4
struct  *scanin2    {"", prompt = "Internal use only"}              # OLDP-4

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
    bool    l_fl_subtract
    bool    l_fl_correct
    bool    l_fl_saturated
    bool    l_fl_nonlinear
    char    l_key_coadds = ""
    char    l_key_lnrs = ""
    char    l_key_bias = ""
    char    l_key_ron = ""
    char    l_key_gain = ""
    char    l_key_sat = ""
    char    l_key_nonlinear = ""
    char    l_key_filter = ""
    char    l_key_grating = ""
    char    l_key_fpmask = ""
    char    l_arraytable = ""
    char    l_configtable = ""
    char    l_specsec = ""
    char    l_key_section = ""
    char    l_key_exptime = ""
    char    l_key_dispaxis = ""
    char    l_key_dark = ""
    char    l_key_mode = ""
    char    l_val_dark = ""
    real    l_pixscale
    char    l_shiftimage = ""
    real    l_shiftx
    real    l_shifty
    char    l_obstype = ""
    bool    l_fl_inter
    bool    l_verbose
    char    l_key_obstype = ""

    char    inimg, outimg, extn, phu, sections, expr, oldmdf
    char    tmpexpand, tmpfile, tmproot, tmpsummary, tmpout
    char    tmphead, tmpopen, tmpmdfrow, tmparray, tmpsections
    char    tmpmdffile, tmpoffset, tmpnotset
    char    l_refsec, refexpr, dimstr, tmp1, tmp2
    bool    already, domdf, havemdf, havedark, correlate
    int     nimages, junk
    struct  sdate
    int     bpmxsize, bpmysize, xsize, ysize, len
    char    tmpsci, tmpvar, tmpdq, tmpinimage, tmpcoordsin
    char    tmpsci2, tmpvar2, tmpdq2, tmptest
    int     dispaxis, nxpix, nypix, k
    bool    bad, oldtimer, debug, hasdate, firstmdf, first
    real    well, ronref, gain
    real    linearlimit, ron, sat, coadds, lnrs, exptime
    real    nonlinlimit, coeff1, coeff2, coeff3
    real    linlimit, ref1, ref2, slope, yint
    real    biasvolt
    char    filter, shiftfile
    char    varexpression, tempexpression, badhdr, dark
    real    pixscl, cval1, cval2, x
    char    grating, fpmask, offsetsection, specsection
    char    date, keyfound, procmode, NONE
    int     count, year, day, month
    bool    fixshift, fixwcs, wcsok, defshift, found, others
    real    dec, ra, decoff, raoff, pix, pa, cd11, cd12, cd21, cd22
    real    ranod, decnod, xobj, yobj, xpix, ypix, rmax
    real    cd11obs, cd12obs, cd21obs, cd22obs, jacob
    real    pi, cospa, sinpa, torad, abserr, tiny, errlim
    real    offset, stddev, midpt

    struct  l_struct

    status  = 1
    debug   = no
    pi      = 3.1415926
    torad   = 2*pi/360.0
    tiny    = 1e-10
    errlim  = 0.05    # 5% of a pixel error
    NONE        = "none"
    xobj    = INDEF
    yobj    = INDEF
    rmax    = 2.0 ** 2        

    cache ("gemextn", "keypar", "tinfo", "tabpar", "nsoffset")

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
    l_fl_subtract = fl_subtract
    l_fl_correct =  fl_correct
    l_fl_saturated =    fl_saturated
    l_fl_nonlinear =    fl_nonlinear
    junk = fscan (  arraytable, l_arraytable)
    junk = fscan (  configtable, l_configtable)
    junk = fscan (  specsec, l_specsec)
    l_pixscale =    pixscale
    junk = fscan (  shiftimage, l_shiftimage)
    l_shiftx =  shiftx
    l_shifty =  shifty
    junk = fscan (  obstype, l_obstype)
    l_fl_inter =    fl_inter
    l_verbose =     verbose

    # Section may contain spaces
    if ("" != l_specsec)
        if (substr (l_specsec, 1, 1) == "[") l_specsec = specsec

    tmpexpand = mktemp("tmpexpand")
    tmpfile = mktemp ("tmpfile") 
    tmproot = mktemp ("tmproot") 
    tmpout = mktemp ("tmpout") 
    tmpsummary = mktemp ("tmpsummary") 
    tmphead = mktemp ("tmphead") 
    tmpopen = mktemp ("tmpopen") 
    tmpmdfrow = mktemp ("tmpmdfrow") //".fits"
    tmparray = mktemp ("tmparray")//".fits"
    tmpsections = mktemp ("tmpsections")
    tmpmdffile = mktemp ("tmpmdffile")
    tmpcoordsin = mktemp ("tmpcoordsin")
    tmpnotset = mktemp ("tmpnotset")
    tmptest = mktemp ("tmptest")

    tmpvar = ""
    tmpsci = ""
    tmpdq = ""
    tmpvar2 = ""
    tmpsci2 = ""
    tmpdq2 = ""
    tmpinimage = ""
    firstmdf = yes
    filter = ""
    fpmask = ""
    grating = ""

    # Shared definitions, define parameters from nsheaders

    badhdr = ""
    junk = fscan (  nsheaders.sci_ext, l_sci_ext)
    if ("" == l_sci_ext) badhdr = badhdr + " sci_ext"
    junk = fscan (  nsheaders.var_ext, l_var_ext)
    if ("" == l_var_ext) badhdr = badhdr + " var_ext"
    junk = fscan (  nsheaders.dq_ext, l_dq_ext)
    if ("" == l_dq_ext) badhdr = badhdr + " dq_ext"
    junk = fscan (  nsheaders.key_coadds, l_key_coadds)
    if ("" == l_key_coadds) badhdr = badhdr + " key_coadds"
    junk = fscan (  nsheaders.key_lnrs, l_key_lnrs)
    if ("" == l_key_lnrs) badhdr = badhdr + " key_lnrs"
    junk = fscan (  nsheaders.key_bias, l_key_bias)
    if ("" == l_key_bias) badhdr = badhdr + " key_bias"
    junk = fscan (  nsheaders.key_ron, l_key_ron)
    if ("" == l_key_ron) badhdr = badhdr + " key_ron"
    junk = fscan (  nsheaders.key_gain, l_key_gain)
    if ("" == l_key_gain) badhdr = badhdr + " key_gain"
    junk = fscan (  nsheaders.key_sat, l_key_sat)
    if ("" == l_key_sat) badhdr = badhdr + " key_sat"
    junk = fscan (  nsheaders.key_nonlinear, l_key_nonlinear)
    if ("" == l_key_nonlinear) badhdr = badhdr + " key_nonlinear"
    junk = fscan (  nsheaders.key_filter, l_key_filter)
    if ("" == l_key_filter) badhdr = badhdr + " key_filter"
    junk = fscan (  nsheaders.key_grating, l_key_grating)
    if ("" == l_key_grating) badhdr = badhdr + " key_grating"
    junk = fscan (  nsheaders.key_fpmask, l_key_fpmask)
    if ("" == l_key_fpmask) badhdr = badhdr + " key_fpmask"
    junk = fscan (  nsheaders.key_section, l_key_section)
    if ("" == l_key_section) badhdr = badhdr + " key_section"
    junk = fscan (  nsheaders.key_exptime, l_key_exptime)
    if ("" == l_key_exptime) badhdr = badhdr + " key_exptime"
    junk = fscan (  nsheaders.key_dispaxis, l_key_dispaxis)
    if ("" == l_key_dispaxis) badhdr = badhdr + " key_dispaxis"
    junk = fscan (  nsheaders.key_dark, l_key_dark)
    if ("" == l_key_dark) badhdr = badhdr + " key_dark"
    junk = fscan (  nsheaders.key_mode, l_key_mode)
    if ("" == l_key_mode) badhdr = badhdr + " key_mode"
    junk = fscan (  nsheaders.val_dark, l_val_dark)
    if ("" == l_val_dark) badhdr = badhdr + " val_dark"
    junk = fscan (  nsheaders.key_obstype, l_key_obstype)
    if ("" == l_key_obstype) badhdr = badhdr + " key_obstype"


    if (debug) print ("validating input")

    if ("" == l_logfile) {
        junk = fscan (nifs.logfile, l_logfile)
        if ("" == l_logfile) {
            l_logfile = "nifs.log"
            printlog ("WARNING - NFPREPARE: Both nfprepare.logfile and \
                nifs.logfile are", l_logfile, verbose+) 
            printlog ("                     undefined.  Using " \
                // l_logfile, l_logfile, verbose+)
        }
    }
    date | scan (sdate) 
    printlog ("----------------------------------------------------------\
        ---------------------", l_logfile, l_verbose) 
    printlog ("NFPREPARE -- " // sdate, l_logfile, l_verbose) 
    printlog (" ", l_logfile, l_verbose) 

    if ("" != badhdr) {
        printlog ("ERROR - NFPREPARE: Parameter(s) missing from \
            nsheaders: " // badhdr, l_logfile, verbose+) 
        goto clean
    }
    if (l_fl_vardq && ! l_fl_saturated) {
        printlog ("WARNING - NFPREPARE: Saturated pixels not flagged.", \
            l_logfile, verbose+) 
    }
    if (l_fl_vardq && ! l_fl_nonlinear) {
        printlog ("WARNING - NFPREPARE: Non-linear pixels not flagged.", \
            l_logfile, verbose+) 
    }

    if (debug) print ("shift options")

    shiftfile = ""
    gemextn (l_shiftimage, check="exists", process="expand", index="0",
        extname="", extversion="", ikparams="", omit="", replace="", 
        outfile="STDOUT", logfile="dev$null", glogpars="",
        verbose=l_verbose) | scan (shiftfile)
    if ("" != shiftfile) {
        printlog ("NFPREPARE: Shift from " // shiftfile, \
            l_logfile, l_verbose)
        hselect (shiftfile, "MDF_XSHF", yes) | scan (l_shiftx)
        hselect (shiftfile, "MDF_YSHF", yes) | scan (l_shifty)
    } else if (no == ("" == l_shiftimage)) {
        printlog ("NFPREPARE: ERROR - Cannot open shiftimage " \
            // l_shiftimage, l_logfile, verbose+)
        goto clean
    }

    correlate = isindef (l_shiftx) && isindef (l_shifty)
    defshift = ! isindef (l_shiftx) && ! isindef (l_shifty)
    if (defshift) defshift = l_shiftx == 0 && l_shifty == 0

    if (debug) print ("check configtable")

    domdf = access (l_configtable)
    if (no == domdf) {
        printlog ("WARNING - NFPREPARE: Config table " // l_configtable \
            // " not found.", l_logfile, verbose+) 
    }

    if (debug) print ("process in and output files")        

    # Check the rawpath name for a final /
    if ("" != l_rawpath) {
        len = strlen (rawpath)
        if (substr (l_rawpath, len, len) != "/")
            l_rawpath = l_rawpath // "/"
    }

    # Expand list (in case comma separated etc)       
    gemextn (l_inimages, check="", process="none", index="", extname="",
        extversion="", ikparams="", omit="", replace="", outfile=tmpexpand,
        logfile="", glogpars="", verbose=l_verbose)
    if (0 < gemextn.fail_count) {
        printlog ("ERROR - NFPREPARE: Bad syntax in inimages.", \
            l_logfile, verbose+) 
        goto clean
    }
    if (0 == gemextn.count) {
        printlog ("ERROR - NFPREPARE: No input images.", l_logfile, verbose+) 
        goto clean
    }

    # Add directory to start and validate
    gemextn ("@" // tmpexpand, check="exists,mef", process="none", index="",
        extname="", extversion="", ikparams="", omit="kernel,exten",
        replace="^%%" // l_rawpath // "%", outfile=tmpfile, logfile="",
        glogpars="", verbose=l_verbose)
    if (0 < gemextn.fail_count) {
        printlog ("ERROR - NFPREPARE: Input images missing.", \
            l_logfile, verbose+) 
        goto clean
    }
    nimages = gemextn.count

    # Generate output images
    gemextn (l_outimages, check="absent", process="none", index="", extname="",
        extversion="", ikparams="", omit="kernel,exten", replace="",
        outfile=tmpout, logfile="", glogpars="", verbose=l_verbose)
    if (gemextn.fail_count != 0) {
        printlog ("ERROR - NFPREPARE: Existing or incorrectly formatted \
            output files", l_logfile, verbose+) 
        goto clean
    }

    # If tmpout is empty, the output files names should be 
    # created with the prefix parameter
    if (gemextn.count == 0) {
        gemextn ("@" // tmpfile, check="", process="none", index="", 
            extname="", extversion="", ikparams="", omit="path", replace="",
            outfile=tmproot, logfile="", glogpars="", verbose=l_verbose)
        if (0 < gemextn.fail_count) {
            printlog ("ERROR - NFPREPARE: Unexpected error.", \
                l_logfile, verbose+) 
            goto clean
        }
        gemextn ("%^%" // l_outprefix // "%" // "@" // tmproot, check="absent",
            process="none", index="", extname="", extversion="", ikparams="",
            omit="kernel,exten", replace="", outfile=tmpout, logfile="",
            glogpars="", verbose=l_verbose)
        if (gemextn.fail_count != 0 || gemextn.count == 0) {
            printlog ("ERROR - NFPREPARE: No or incorrectly formatted \
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
            printlog ("ERROR - NFPREPARE: Too few output files.", \
                l_logfile, verbose+) 
            goto clean
        }

        keyfound = ""
        already = no
        hselect (inimg//"[0]", "*PREPAR*", yes) | scan (keyfound)
        if (keyfound != "") {
            already = yes
            printlog ("WARNING - NFPREPARE: Image " // inimg \
                // " already fixed using *PREPARE.", l_logfile, verbose+) 
            printlog ("                     No futher processing done.",
                l_logfile, verbose+) 
        }

        print (inimg // " " // outimg // " " // already, >> tmpsummary)
    }

    if (fscan (scanin2, outimg) != EOF) {
        printlog ("ERROR - NFPREPARE: Too many output files.", \
            l_logfile, verbose+) 
        goto clean
    }

    printlog ("Processing " // nimages // " files", l_logfile, l_verbose) 


    # Check for existence of bad pixel mask

    # TODO - add check for .pl or .fits, and convert .pl to fits
    # if required by addmasks

    if (debug) print ("checking bpm")

    if (l_fl_vardq && "" != l_bpm && ! imaccess (l_bpm)) {
        printlog ("WARNING - NFPREPARE: Bad pixel mask " // l_bpm \
            // " not found.", l_logfile, verbose+) 
        l_bpm = ""
    } else if (l_fl_vardq && "" == l_bpm) {
        printlog ("NFPREPARE: Bad pixel mask not specified.", \
            l_logfile, verbose+) 
    }
    if ("" == l_bpm) l_bpm = "none"    # For logging


    # Check size of BPM

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
            printlog ("WARNING - NFPREPARE: No size information in \
                bad pixel mask header.", l_logfile, verbose+) 
            printlog ("                    Not using BPM to generate \
                data quality planes.", l_logfile, verbose+) 
            l_bpm = "none"
        } 

#else {

#            scanin1 = tmpsummary
#            while (fscan (scanin1, inimg, outimg, already) != EOF) {
#                bad = no
#                extn = inimg // "[1]"
#                keypar (extn, "i_naxis1")
#                if (keypar.found)
#                    xsize = int (keypar.value)
#                else
#                    bad = yes
#                keypar (extn, "i_naxis2")
#                if (keypar.found)
#                    ysize = int (keypar.value)
#                else
#                    bad = yes
#                if (bad) {
#                    printlog ("ERROR - NFPREPARE: No size information \
#                        in " // extn // ".", l_logfile, verbose+) 
#                    goto clean
#                }
#                if ((xsize != bpmxsize) || (ysize != bpmysize)) {
#                    printlog ("WARNING - NFPREPARE: Input images and \
#                        BPM are not the same size.", l_logfile, verbose+) 
#                    printlog ("                    Not using BPM to \
#                        generage data quality plane.", l_logfile, verbose+)
#                    l_bpm = "none"
#                    break
#                }
#            }
#        }        
    }


    # If we're shifting, we want the first item in the input list to
    # have an obstype that matches that given (to avoid asymmetric
    # data).

    if (correlate && (no == ("" == l_obstype)) && \
        (no == ("(unused)" == l_key_obstype))) {

        tmpoffset = mktemp ("tmpsummary")
        scanin1 = tmpsummary
        found = no
        others = no
        while (fscan (scanin1, inimg, outimg, already) != EOF) {
            if ((no == found) && (no == already)) {
                phu = inimg // "[0]"
                hselect (phu, l_key_obstype, yes) | scan (l_struct)
                if (l_struct == l_obstype) {
                    print (inimg // " " // outimg // " " // already, \
                        > tmpoffset)
                    found = yes
                    printlog ("NFPREPARE: Will measure MDF offset using " \
                        // inimg, l_logfile, l_verbose)
                    printlog ("           (" // l_key_obstype // "=" \
                        // l_obstype // ")", l_logfile, verbose=l_verbose)
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
            printlog ("ERROR - NFPREPARE: No suitable candidate was found \
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
                print (inimg // " " // outimg // " " // already, >> tmpoffset)
            }
        }
    }


    # Help users with automatic shifting etc
    if (correlate) {
        printlog ("NFPREPARE: Both shiftx and shifty are INDEF.  The \
            shift for MDF will be", l_logfile, l_verbose)
        printlog ("           measured using cross-correlation.", \
            l_logfile, l_verbose)
        printlog ("           You should verify that the measured \
            shift is correct by", l_logfile, l_verbose)
        printlog ("           examining the cut data.", \
            l_logfile, l_verbose)
    } else if (defshift) {
        printlog ("NFPREPARE: Both shiftx and shifty are 0.  This will \
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
        printlog ("           in just one call to nfprepare or use the \
            shiftimage parameter", l_logfile, l_verbose)
        printlog ("           (or shiftx and shifty directly) to copy \
            the shift from a ", l_logfile, l_verbose)
        printlog ("           previously processed image on subsequent \
            calls to nfprepare.", l_logfile, l_verbose)
    }


    fpmask = ""

    first = yes
    scanin1 = tmpsummary
    count = 0
    while (fscan (scanin1, inimg, outimg, already) != EOF) {
        count = count + 1
        if (debug) print (count // ": " // inimg // " -> " // outimg)

        # tmp files used within this loop
        if (no == already) {
            tmpvar = mktemp ("tmpvar")
            tmp1 = mktemp ("tmp1")
            tmp2 = mktemp ("tmp2")
            tmpsci = mktemp ("tmpsci")
            tmpdq = mktemp ("tmpdq")
            tmpvar2 = mktemp ("tmpvar2")
            tmpsci2 = mktemp ("tmpsci2")
            tmpdq2 = mktemp ("tmpdq2")
            tmpinimage = mktemp ("tmpinimage")
        }

        # Set read noise, gain, and saturation for the data

        phu = inimg // "[0]"
        havedark = no

        dispaxis = INDEF
        hselect (phu, l_key_dispaxis, yes) | scan (dispaxis)
        if (isindef (dispaxis)) {
            dispaxis = 1
        } else {
            printlog ("WARNING - NFPREPARE: " // l_key_dispaxis \
                // " already set to " // dispaxis, l_logfile, verbose+) 
        }

        keypar (phu, l_key_dark, silent+) 
        if (no == keypar.found) {
            printlog ("WARNING - NFPREPARE: Could not read " \
                // l_key_dark // " for " // inimg, l_logfile, verbose+) 
            dark = ""
        } else {
            dark = keypar.value
        }

        if (debug)
            print (l_val_dark//", "//dark//", "//strstr (l_val_dark, dark))
        if (strstr (l_val_dark, dark) > 0) {
            if (no == l_fl_dark_mdf) {
                printlog ("WARNING - NFPREPARE: " // inimg \
                    // " is a dark", l_logfile, verbose+)
                printlog ("                     and so will not have \
                    an MDF added.", l_logfile, verbose+)
                havedark = yes
            } else {
                printlog ("WARNING - NFPREPARE:"//inimg//" is a dark, but \
                    fl_dark_mdf is set", l_logfile, verbose+)
                printlog ("                     so it will have an MDF added.",
                    l_logfile, verbose+)
                havedark = no
            }
        }

        if (domdf && ! havedark) {

            keypar (phu, l_key_fpmask, silent+) 
            if (no == keypar.found) {
                printlog ("WARNING - NFPREPARE: Could not read fpmask \
                    for " // inimg, l_logfile, verbose+) 
                fpmask = ""
            } else {
                fpmask = keypar.value
            }

            keypar (phu, l_key_filter, silent+) 
            if (no == keypar.found) {
                printlog ("WARNING - NFPREPARE: Could not read filter \
                    for " // inimg, l_logfile, verbose+) 
                filter = ""
            } else {
                filter = keypar.value
            }

            keypar (phu, l_key_grating, silent+) 
            if (no == keypar.found) {
                printlog ("WARNING - NFPREPARE: Could not read grating \
                    for " // inimg, l_logfile, verbose+) 
                grating = ""
            } else {
                grating = keypar.value
            }
        }


        keypar (phu, l_key_coadds, silent+) 
        if (no == keypar.found) {
            printlog ("ERROR - NFPREPARE: Could not read number of \
                coadds from header", l_logfile, verbose+) 
            printlog ("                  of file " // inimg, 
                l_logfile, verbose+) 
            goto clean
        } else {
            coadds = real (keypar.value) 
        }

        keypar (phu, l_key_lnrs, silent+) 
        if (no == keypar.found) {
            printlog ("ERROR - NFPREPARE: Could not read number of \
                non-destructive reads from header", l_logfile, verbose+) 
            printlog ("                  of file " // inimg, 
                l_logfile, verbose+) 
            goto clean
        } else {
            lnrs = real (keypar.value) 
        }

        keypar (phu, l_key_bias, silent+) 
        if (no == keypar.found) {
            printlog ("ERROR - NFPREPARE: Could not read " // l_key_bias \
                // " from header.", l_logfile, verbose+) 
            goto clean
        } else {
            biasvolt = real (keypar.value) 
            if (debug) print ("bias: " // biasvolt)
        }

        # Read array info from table

        if (debug) print ("reading arraytable")

        if (no == access (l_arraytable) ) {
            printlog ("ERROR - NFPREPARE: Array table " // l_arraytable \
                // " not found.", l_logfile, verbose+) 
            goto clean
        } else {

            tdelete (tmparray, go_ahead+, verify-, >& "dev$null")

            expr = "abs( bias - abs ( " // biasvolt // " )) < 0.1"
            if (debug) print (expr)
            tselect (intable = l_arraytable, outtable = tmparray, expr = expr)

            tinfo (table = tmparray, >& "dev$null")
            if (1 != tinfo.nrows) {
                printlog ("ERROR - NFPREPARE: No array calibration data \
                    for bias voltage of " // biasvolt, l_logfile, verbose+)
                printlog ("                   in " // l_arraytable, \
                    l_logfile, verbose+) 
                goto clean
            }

            tabpar (table = tmparray, column = "readnoise", row = 1, format-)
            if (tabpar.undef) {
                printlog ("ERROR - NFPREPARE: No readnoise for bias \
                    voltage of " // biasvolt, l_logfile, verbose+)
                printlog ("                   in " // l_arraytable, \
                    l_logfile, verbose+) 
                goto clean
            } else ronref = real (tabpar.value)

            ron = (ronref * sqrt (coadds) ) / (sqrt (lnrs)) 

            tabpar (table = tmparray, column = "gain", row = 1, format-)
            if (tabpar.undef) {
                printlog ("ERROR - NFPREPARE: No gain for bias voltage \
                    of " // biasvolt, l_logfile, verbose+)
                printlog ("                   in " // l_arraytable, \
                    l_logfile, verbose+) 
                goto clean
            } else gain = real (tabpar.value)

            tabpar (table = tmparray, column = "well", row = 1, format-)
            if (tabpar.undef) {
                printlog ("ERROR - NFPREPARE: No well defined in \
                    table " // l_arraytable, l_logfile, verbose+) 
                goto clean
            } else well = real (tabpar.value)

            tabpar (table = tmparray, column = "linearlimit", row = 1, format-)
            if (tabpar.undef) {
                printlog ("ERROR - NFPREPARE: No linearlimit defined in \
                    table " // l_arraytable, l_logfile, verbose+) 
                goto clean
            } else linearlimit = real (tabpar.value)

            tabpar (table = tmparray, column = "nonlinearlimit", row = 1,
                format-)
            if (tabpar.undef) {
                printlog ("ERROR - NFPREPARE: No nonlinearlimit defined \
                    in table " // l_arraytable, l_logfile, verbose+) 
                goto clean
            } else nonlinlimit = real (tabpar.value)

            if (l_fl_correct) {
                tabpar (table = tmparray, column = "coeff1", row = 1, format-)
                if (tabpar.undef) {
                    printlog ("ERROR - NFPREPARE: No coeff1 defined in \
                        table " // l_arraytable, l_logfile, verbose+) 
                    goto clean
                } else coeff1 = real (tabpar.value)

                tabpar (table = tmparray, column = "coeff2", row = 1, format-)
                if (tabpar.undef) {
                    printlog ("ERROR - NFPREPARE: No coeff2 defined in \
                        table " // l_arraytable, l_logfile, verbose+) 
                    goto clean
                } else coeff2 = real (tabpar.value)

                tabpar (table = tmparray, column = "coeff3", row = 1, format-)
                if (tabpar.undef) {
                    printlog ("ERROR - NFPREPARE: No coeff3 defined in \
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
            if (("" != fpmask) && ("" != grating) && ("" != filter)) {
                tdelete (table = tmpmdfrow, go_ahead+, verify-, >& "dev$null")
                tselect (intable = l_configtable, outtable = tmpmdfrow,
                    expr = "fpmask='" // fpmask // "'&&fpmask='" \
                    // fpmask // "'&&grating='" // grating \
                    // "'&&filter='" // filter // "'")
                if (debug) tinfo (table = tmpmdfrow)
                tinfo (table = tmpmdfrow, >>& tmpopen)
                havemdf = (0 != tinfo.nrows)
                if (havemdf) {
                    tabpar (table = tmpmdfrow, column = "mdf", row = 1, 
                        format-)
                    sections = tabpar.value
                    tabpar (table = tmpmdfrow, column = "pixscale", row = 1,
                        format-)
                    pixscl = real (tabpar.value)
                    tabpar (table = tmpmdfrow, column = "mode", row = 1,
                        format-)
                    procmode = strupr (tabpar.value)
                } else {
                    printlog ("ERROR - NFPREPARE: Cannot find entry in " //
                        l_configtable // " matching", l_logfile, verbose+)
                    printlog ("                     fpmask  = " // fpmask,
                        l_logfile, verbose+)
                    printlog ("                     grating = " // grating,
                        l_logfile, verbose+)
                    printlog ("                     filter  = " // filter, 
                        l_logfile, verbose+)
                    goto clean
                }
            }
        }

        # specsec can be a section or an MDF file
        if (no == havemdf && no == havedark && "" != l_specsec) {
            if (debug) print ("checking specsec")
            havemdf = ! (substr (l_specsec, 1, 1) == "[")
            if (havemdf) sections = l_specsec
        }

        # Removed "/ gain" here - well in ADU
        if (debug) print ("well, coadds: " // well // " " // coadds)
        sat = int (well * coadds)
        if (l_fl_correct) linlimit = int (sat * linearlimit) 
        else              linlimit = int (sat * nonlinlimit) 
        if (debug) print ("sat, linlimit " // sat // " " // linlimit)

        # Read date/instrument to control various fixes

        date = INDEF
        hasdate = no
        hselect (phu, "DATE-OBS", yes) | scan (date)
        if (isindef (date)) {
            printlog ("WARNING - NFPREPARE: Cannot read observation date.  \
                No date-specific fixes applied.", l_logfile, verbose+) 
        } else {
            year = int (substr (date, 1, 4))
            month = int (substr (date, 6, 7))
            day = int (substr (date, 9, 10))
            hasdate = yes
        }

        keypar (phu, "INSTRUME", silent+)

        extn = inimg // "[1]"

        if (no == already) {

            # Substract reference pixels

    printlog ("-----------------------------------------------------\
        -----------------------", l_logfile, l_verbose) 
  if (l_fl_subtract) {
    printlog ("NFPREPARE: Subtracting reference pixels.", l_logfile, l_verbose)
    imgets(inimg//"[1]","i_naxis1", >& "dev$null")
    nxpix=int(imgets.value)
    imgets(inimg//"[1]","i_naxis2", >& "dev$null")
    nypix=int(imgets.value)
    # Loop over amplifiers.
    for (k=1; k<=4; k+=1) {
      l_refsec="["//str(int((k-1)*nxpix/4+1))//":"//str(int(k*nxpix/4))//\
        ",1:4]"
      imstat(extn//l_refsec,fields="mean,stddev",
        lower=0.,upper=sat,format-) | scan(ref1,stddev)
      imstat(extn//l_refsec,fields="mean",
        lower=(ref1-3*stddev),upper=(ref1+3*stddev),format-) | scan(ref1)
      # Ignore erroneous top row of reference pixels.
      l_refsec="["//str(int((k-1)*nxpix/4+1))//":"//str(int(k*nxpix/4))//","//
        str(nypix-3)//":"//str(nypix-1)//"]"
      imstat(extn//l_refsec,fields="mean,stddev",
        lower=0.,upper=sat,format-) | scan(ref2,stddev)
      imstat(extn//l_refsec,fields="mean",
        lower=(ref2-3*stddev),upper=(ref2+3*stddev),format-) | scan(ref2)
      slope=(ref2-ref1)/real(nxpix/4*(nypix-4))
      yint=ref1-slope*real(2*nxpix/4)
      if (k == 1 || k == 3) {
        # Amplifiers 1 & 3 read left-to-right.
        refexpr="( (I > "//str((k-1)*nxpix/4)//") && (I < "//str(k*nxpix/4+1)//
          ") ) ? a+b*real("//str(nxpix/4)//"*(J-1)+I-"//str((k-1)*nxpix/4)//
          ") : 0"
      } else {
        # Amplifiers 2 & 4 read right-to-left.
        refexpr="( (I > "//str((k-1)*nxpix/4)//") && (I < "//str(k*nxpix/4+1)//
          ") ) ? a+b*real("//str(nxpix/4)//"*(J-1)+"//str(nxpix/4)//"-I+1+"//
          str((k-1)*nxpix/4)//") : 0"
      }
      dimstr=str(nxpix)//","//str(nypix)
      if (k == 1) {
        imexpr(refexpr, tmp1, yint, slope, dims=dimstr, intype="real", 
          outtype="real", refim="auto", bwidth=0, btype="nearest", 
          rangecheck+, verbose-)
      } else {
        imexpr(refexpr, tmp2, yint, slope, dims=dimstr, intype="real",
          outtype="real", refim="auto", bwidth=0, btype="nearest", 
          rangecheck+, verbose-)
        imarith(tmp1,"+", tmp2, tmp1, title="", divzero=0.0, hparams="",
          pixtype="", calctyp="", verbose-, noact-)
        imdelete(tmp2, verify-, >& "dev$null")
      }
    }
    imarith(extn,"-",tmp1,tmpinimage,
      title="",divzero=0.0,hparams="",pixtype="",calctyp="",verbose-,noact-)
    imdelete(tmp1,verify-, >& "dev$null")
  } else {
    imcopy(extn,tmpinimage, >& "dev$null")
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
                printlog ("WARNING - NFPREPARE: Not correcting for \
                    non-linearity.", l_logfile, verbose+) 
            }


            # Check header for values
            keypar (phu, l_key_gain, silent+) 
            if (keypar.found) {
                gain = real (keypar.value) 
                printlog ("WARNING - NFPREPARE: Gain already set in image \
                    header.  Using a gain", l_logfile, verbose+) 
                printlog ("                    of " // gain // \
                    " electrons per ADU.", l_logfile, verbose+) 
            }

            keypar (phu, l_key_ron, silent+) 
            if (keypar.found) {
                ron = real (keypar.value) 
                printlog ("WARNING - NFPREPARE: Read noise already set in \
                    image header.", l_logfile, verbose+) 
                printlog ("                    Using a read noise of \
                    " // ron // " electrons.", l_logfile, verbose+) 
            }

            keypar (phu, l_key_sat, silent+) 
            if (keypar.found) {
                sat = real (keypar.value) 
                printlog ("WARNING - NFPREPARE: Saturation level already set \
                    in image header.", l_logfile, verbose+) 
                printlog ("                    Using a saturation level of \
                    " // sat // " ADU.", l_logfile, verbose+) 
            }


            if (l_fl_vardq) {

                # Pad the science frame for overfilled pixels.

                nfpad (tmpinimage, tmpsci2, exttype="SCI", logfile=l_logfile,
                    verbose=l_verbose)

                # Create the variance frame, if it doesn't already exist
                # The variance frame is generated as:
                # var = (read noise/gain)^2 + max(data,0.0)/gain

                if (debug) print ("vardq")
                extn = inimg // "[" // l_var_ext // "]"

                if (no == imaccess (extn) ) {

                    varexpression = "((max(a,0.0))/" // gain // \
                        "+(" // ron // "/" // gain // ")**2)"
                    imexpr (varexpression, tmpvar, tmpinimage, 
                        outtype = "real", verbose-) 

                } else {
                    printlog ("WARNING - NFPREPARE: Variance frame already \
                        exists for " // inimg // ".", l_logfile, verbose+) 
                    printlog ("                     New variance frame not \
                        created.", l_logfile, verbose+) 
                    imcopy (extn, tmpvar, verbose-) 
                }

                nfpad (tmpvar, tmpvar2, exttype="VAR", logfile=l_logfile,
                    verbose=l_verbose)

                # Create the DQ frame, if it doesn't already exist
                # The preliminary DQ frame is constructed by using 
                # the bad pixel mask to set bad pixels to 1, pixels 
                # in the non-linear regime to 2, and saturated pixels to 4.
 
    if ("none" != l_bpm) {

               keypar (tmpsci2, "i_naxis1")
                if (keypar.found)
                    xsize = int (keypar.value)
                else
                    bad = yes
                keypar (tmpsci2, "i_naxis2")
                if (keypar.found)
                    ysize = int (keypar.value)
                else
                    bad = yes
                if (bad) {
                    printlog ("ERROR - NFPREPARE: No size information \
                        in " // extn // ".", l_logfile, verbose+) 
                    goto clean
                }
                if ((xsize != bpmxsize) || (ysize != bpmysize)) {
                   printlog ("WARNING - NFPREPARE: Input images and \
                        BPM are not the same size.", l_logfile, verbose+) 
                    printlog ("                    Not using BPM to \
                        generate data quality plane.", l_logfile, verbose+)
                    l_bpm = "none"
                    break
                }
     }
                
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
                    imexpr (tempexpression, tmpsci, tmpinimage, 
                        outtype="short", verbose-) 
                    nfpad (tmpsci, tmptest, exttype="DQ", logfile=l_logfile,
                        verbose=l_verbose)

                    # If there's no BPM, then just keep the saturated 
                    # pixels
                    if (l_bpm == "none") {
                        addmasks (tmptest // ".fits", tmpdq2 // ".fits", \ 
                            "im1") 
                    } else {
                        addmasks (tmptest // ".fits," // l_bpm, tmpdq2 // \
                            ".fits", "im1 || im2") 
                    }

                    imdelete (tmpsci, verify-) 

                } else {
                    printlog ("WARNING - NFPREPARE: Data quality frame \
                        already exists for " // inimg // ".", l_logfile, 
                        verbose+) 
                    printlog ("                     New DQ frame not created.",
                        l_logfile, verbose+) 
                    imcopy (extn, tmpdq // ".fits", verbose-) 
                    nfpad (tmpdq, tmpdq2, exttype="DQ", logfile=l_logfile,
                        verbose=l_verbose)
                }

                # Pack up the results and clean up

                wmef (tmpsci2 // "," // tmpvar2 // ".fits," // tmpdq2 // \
                    ".fits", outimg, 
                    extnames=l_sci_ext//","//l_var_ext//","// l_dq_ext,
                    phu=inimg//".fits[0]", verbose-) 
                imdelete (tmpsci2 // "," // tmpvar2 // "," // tmpdq2 //"," \
                    // tmptest, verify-, >& "dev$null")

                # don't hide error messages - this may be a source
                # of errors with old iraf versions and we're not
                # going to know if the error message is hidden

                if (wmef.status != 0) {
                    printlog ("ERROR - NFPREPARE: Could not write final \
                        MEF file (WMEF).", l_logfile, verbose+) 
                    goto clean
                }

                gemhedit (outimg // "[" // l_var_ext // "]", "EXTVER", 1, "")
                gemhedit (outimg // "[" // l_dq_ext // "]", "EXTVER", 1, "")
                imdelete (tmpvar // "," // tmpdq, verify-, >& "dev$null") 

            } else {

                # No var/dq

                nfpad (tmpinimage, tmpsci2, exttype="SCI", logfile=l_logfile,
                    verbose=l_verbose)

                wmef (tmpsci2, outimg, extnames = l_sci_ext,
                    phu = inimg // ".fits[0]", verbose-, >& "dev$null") 
                if (wmef.status != 0) {
                    printlog ("ERROR - NFPREPARE: Could not write final \
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
                    gemextn (sections, check="exists,mef,table", 
                        process="expand", index="1", extname="", 
                        extversion="", ikparams="", omit="exten", replace="", 
                        outfile=tmpmdffile, logfile="", glogpars="", verbose-)
                    if ((0 != gemextn.fail_count) || (1 != gemextn.count)) {
                        printlog ("ERROR - NFPREPARE: MDF file \
                            " // sections // " not found, or incorrect \
                            format.", l_logfile, verbose+)
                        goto clean
                    }
                    scanin2 = tmpmdffile
                    junk = fscan (scanin2, sections)

                    # Auto-shift
                    if (correlate && firstmdf) {
                        oldmdf = sections
                        imdelete (tmpsections, verify-, >& "dev$null")
                        gemextn (sections, check="", process="none", index="",
                            extname="", extversion="", ikparams="",
                            omit="path,exten,kernel,index", replace="",
                            outfile="STDOUT", logfile="dev$null", glogpars="",
                            verbose-) | scan (tmpsections)
                        tmpsections = mktemp (tmpsections // "-")
                        if (debug) print (tmpsections)
                        nsoffset (outimg, sections, outimages=tmpsections,
                            outprefix="", axis=3-dispaxis, pixscale=pixscl,
                            fl_apply+, fl_all-, 
                            fl_project=(no==("XD"==procmode)),
                            fl_inter=l_fl_inter, logfile=l_logfile,
                            verbose=l_verbose)
                        sections = tmpsections
                        if (0 != nsoffset.status) {
                            printlog ("ERROR - NFPREPARE: Auto-shift failed \
                                for MDF.", l_logfile, verbose+) 
                            goto clean
                        }
                        if (2 == dispaxis) {
                            l_shiftx = nsoffset.shift
                            gemhedit (phu, "MDF_XSHF", l_shiftx, \
                                "MDF X shift (subtracted from x_ccd)")
                        } else {
                            l_shifty = nsoffset.shift
                            gemhedit (phu, "MDF_YSHF", l_shifty, \
                                "MDF Y shift (subtracted from y_ccd)")
                        }
                        firstmdf = no
                    } else if (correlate) {
                        if (sections == oldmdf) {
                            sections = tmpsections
                        } else {
                            printlog ("ERROR - NFPREPARE: MDF cannot change \
                                if using auto-shifting.", l_logfile, l_verbose)
                            goto clean
                        }
                        if (2 == dispaxis) {
                            gemhedit (phu, "MDF_XSHF", l_shiftx, \
                                "MDF X shift (subtracted from x_ccd)")
                        } else {
                            gemhedit (phu, "MDF_YSHF", l_shifty, \
                                "MDF Y shift (subtracted from y_ccd)")
                        }
                    } else {
                        oldmdf = sections
                    }

                    fxinsert (sections, outimg // ".fits[0]", groups="1", \
                        verbose-)
                    thedit (outimg // ".fits[1]", "EXTNAME", "MDF", delete-, \
                        show+, >& "dev$null")

                    if ((no == correlate) && (no == isindef (l_shiftx)) && \
                        (0.0 != l_shiftx)) {

                        tcalc (outimg // ".fits[MDF]", outcol="x_ccd", \
                            equals="x_ccd-" // l_shiftx, datatype="real")
                        printlog ("NFPREPARE: Shifting MDF in X by \
                            "//l_shiftx// " pixels (subtracted from x_ccd)", \
                            l_logfile, l_verbose) 
                        gemhedit (phu, "MDF_XSHF", l_shiftx, \
                            "MDF X shift (subtracted from x_ccd)")
                    }
                    if ((no == correlate) && (no == isindef (l_shifty)) && \
                        (0.0 != l_shifty)) {

                        tcalc (outimg // ".fits[MDF]", outcol="y_ccd", \
                            equals="y_ccd-" // l_shifty, datatype="real")
                        printlog ("NFPREPARE: Shifting MDF in Y by \
                            "//l_shifty//" pixels (subtracted from y_ccd)", 
                            l_logfile, l_verbose) 
                        gemhedit (phu, "MDF_YSHF", l_shifty, \
                            "MDF Y shift (subtracted from y_ccd)")
                    }

                    gemhedit (phu, "MDF_FILE", oldmdf, "MDF source")

                } else {
                    printlog ("Using section " // specsection, l_logfile, \
                        verbose+)
                    gemhedit (phu, l_key_section // "1", specsection, "")
                }
            }

        } else {
            phu = inimg // "[0]"
        }

        # Fix up the headers

        if (no == already) {

            if (debug) print ("headers")

            if (l_fl_vardq) {
                delete (tmphead, verify-, >& "dev$null") 
                printf ("%-8s=                    %d / %-s\n", l_key_dispaxis,
                    dispaxis, "Dispersion axis (along columns)", >> tmphead) 
                printf ("%-8s= %20.5f / %-s\n", "PIXSCALE", pixscl,
                    "Pixel scale in arcsec/pixel", >> tmphead) 
                mkheader (outimg // "[" // l_var_ext // "]", tmphead, append+,
                    verbose-) 
                mkheader (outimg // "[" // l_dq_ext // "]", tmphead, append+,
                    verbose-) 
            }

            gemhedit (outimg // "[" // l_sci_ext // "]", "EXTVER", 1, "")

            delete (tmphead, verify-, >& "dev$null") 
            printf ("%-8s=                    %d / %-s\n", l_key_dispaxis,
                dispaxis, "Dispersion axis (along columns)", >> tmphead) 
            printf ("%-8s= %20.5f / %-s\n", "PIXSCALE", pixscl,
                "Pixel scale in arcsec/pixel", >> tmphead) 
            mkheader (outimg // "[" // l_sci_ext // "]", tmphead, append+,
                verbose-) 
            delete (tmphead, verify-, >& "dev$null") 

            if (l_fl_correct) {
                printf ("%7.5f %7.5f %7.5f\n", coeff1, coeff2, coeff3) \
                    | scan (l_struct)
                gemhedit (phu, "NONLINCR", l_struct,
                    "Non-linear correction applied")
                if (debug) print (coeff1 // "," // coeff2 // "," // coeff3)
            }


            gemhedit (phu, "WMEF", "", "", delete+)
            gemhedit (phu, "GEM-TLM", "", "", delete+)

            # Add comment to NEXTEND
            keypar (phu, "NEXTEND", silent+) 
            if (keypar.found) {
                gemhedit (phu, "NEXTEND", "", "", delete+)
                delete (tmphead, verify-, >& "dev$null") 
                printf ("%-8s= %20.0f / %-s\n", "NEXTEND", real (keypar.value),
                    "Number of extensions", >> tmphead) 

                if (debug) print (phu)
                if (debug) type (tmphead)
                mkheader (phu, tmphead, append+, verbose-) 
            }

            delete (tmphead, verify-, >& "dev$null") 
            printf ("%-8s= '%-18s' / %-s\n", l_key_mode, procmode,
                "Reduction mode", >> tmphead) 
            mkheader (phu, tmphead, append+, verbose-) 

            delete (tmphead, verify-, >& "dev$null") 

            printf ("%-8s= %20.2f / %-s\n", l_key_ron, ron,
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


            printf ("%-8s= %20.5f / %-s\n", "PIXSCALE", pixscl,
                "Pixel scale in arcsec/pixel", >> tmphead) 

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

            printf ("%-8s=                    %d / %-s\n", l_key_dispaxis,
                dispaxis, "Dispersion axis (along columns)", >> tmphead) 

            if (debug) print (phu)
            if (debug) type (tmphead)
            mkheader (phu, tmphead, append+, verbose-) 
            delete (tmphead, verify-, >& "dev$null") 

            # Time stamps
            gemdate ()
            gemhedit (phu, "PREPARE", gemdate.outdate, 
                "UT Time stamp for NFPREPARE")
            gemhedit (phu, "GEM-TLM", gemdate.outdate, 
                "UT Last modification with GEMINI")
            
            printf ("%-8s= \'%-18s\' / %-s\n", "BPMFILE", l_bpm, \
                "Input bad pixel mask file", >> tmphead) 

            if (debug) print (phu)
            if (debug) type (tmphead)
            mkheader (phu, tmphead, append+, verbose-) 
            delete (tmphead, verify-, >& "dev$null") 

            # print output to logfile
            keypar (phu, "APERTURE", silent+) 
            if (substr (keypar.value, strlen (keypar.value) -5, \
                strlen (keypar.value) -4) == "_G") {
                
                fpmask = substr (keypar.value, 1, (strlen (keypar.value) -6) ) 
            } else {
                fpmask = keypar.value
            }

            printlog (" ", l_logfile, l_verbose) 
            printlog ("  n      input file -->      output file", 
                l_logfile, l_verbose) 
            printf ("%3.0d %15s --> %16s \n", count, inimg, outimg) \
                | scan (l_struct) 
            printlog (l_struct, l_logfile, l_verbose) 
            printlog ("                 filter     focal plane       \
                input BPM    RN  gain     sat", l_logfile, l_verbose) 
            printf ("      %17s %15s %15s %5.1f %5.1f %7.0d \n", \
                filter, fpmask, l_bpm, ron, gain, sat) | scan (l_struct) 
            printlog (l_struct, l_logfile, l_verbose)
             
        } else {

            printf ("%3.0d %15s already processed\n", count, inimg) |\
                 scan (l_struct) 
            printlog (l_struct, l_logfile, l_verbose) 

        }


        # Detect cosmic rays and radiation events
        # Currently experimental, with fixed params
        if (l_fl_cravg) {
            printlog ("NFPREPARE: Detecting cosmic rays and radiation events",
                l_logfile, l_verbose)
            craverage (outimg // "[" // l_sci_ext // ",1]", "",
                crmask=outimg // "[" // l_dq_ext // ",1]", average="",
                sigma="", navg=9, nrej=48, nbkg=11, nsig=25, var0=(ron**2),
                var1=gain, crval=8, lcrsig=5, hcrsig=10, crgrow=1, objval=0,
                lobjsig=10, hobjsig=10, objgrow=0)
        }

        # Grow bad pixel detections
        if (l_crradius > 0) {
            printlog ("NFPREPARE: Extending bad pixel mask",
                l_logfile, l_verbose)
            crgrow (outimg // "[" // l_dq_ext // ",1]",
                outimg // "[" // l_dq_ext // ",1]", radius=l_crradius,
                inval=INDEF, outval=INDEF)
        }

        if (imaccess(outimg//"["//l_dq_ext//",1]")) {
            chpixtype (outimg//"["//l_dq_ext//",1]", \
                outimg//"["//l_dq_ext//",1,overwrite]", newpixtype="ushort", \
                oldpixtype="all", verbose-)
            gemhedit (outimg//"["//l_dq_ext//",1]", "OBJECT", "", "", \
                delete-)
        }
    }

    status = 0


clean:
    if (status == 0) {
        printlog (" ", l_logfile, l_verbose) 
        printlog ("NFPREPARE exit status:  good.", l_logfile, l_verbose) 
    }
    printlog ("-----------------------------------------------------\
        -----------------------", l_logfile, l_verbose) 

    scanin1 = "" 
    scanin2 = "" 
    delete (tmpexpand, verify-, >& "dev$null")  
    delete (tmpfile, verify-, >& "dev$null")
    delete (tmproot, verify-, >& "dev$null")
    delete (tmpout, verify-, >& "dev$null")
    delete (tmpsci, verify-, >& "dev$null")
    delete (tmpvar, verify-, >& "dev$null")
    delete (tmpdq, verify-, >& "dev$null")
    delete (tmpsummary, verify-, >& "dev$null")
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
    imdelete (tmpsci2 // "," // tmpvar2 // "," // tmpdq2, verify-, 
        >& "dev$null") 

end

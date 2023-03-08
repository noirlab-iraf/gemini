# Copyright(c) 2002-2017 Association of Universities for Research in Astronomy, Inc

procedure giflat (inflats, outflat)

# Derive normalize imaging flat for GMOS
#
# Version  Feb 28, 2002  MT,IJ,BM  v1.3 release
#          Sept 20, 2002           v1.4 release
#          Oct 23, 2002  BM        generalized to 6-amp mode, changed
#                                  low/high gain criteria for GMOS-S,
#                                  added fl_bias and fl_trim keywords
#          Nov 20, 2002  IJ   fixed sctype of first image
#          Jan 22, 2003  ML   fixed normsec,statsec to work for any
#                             binning and number of amps
#          Feb 16, 2003  BM   fix bug when inputing gireduced images
#          Mar 1, 2003   ??   added crreject to reject options
#          Aug 26, 2003  KL   IRAF2.12 - new/modified parameters
#                               hedit: addonly
#                               imcombine: headers,bpmasks,expmasks,outlimits
#                                    rejmask->rejmasks, plfile->nrejmasks
#                               imstat: nclip,lsigma,usigma,cache
#          Apr 18, 2008  JH   redo file naming conventions. update
#                             variance creation

char   inflats      {prompt="Input flat field images"}                     # OLDP-1-input-primary-combine-suffix=_flat
char   outflat      {prompt="Output flat field image"}                     # OLDP-1-output
char   normsec      {"default",prompt="Image section to get the normalization."} # OLDP-3
bool   fl_scale     {yes,prompt="Scale the flat images before combining?"}    # OLDP-2
char   sctype       {"mean",prompt="Type of statistics to compute for scaling", enum="mean|mode|midpt"} # OLDP-2
char   statsec      {"default",prompt="Image section for relative intensity scaling"} # OLDP-2
char   key_gain     {"GAIN",prompt="Header keyword for gain [e-/ADU]"}      # OLDP-3
bool   fl_stamp     {no, prompt="Input is stamp image"}

char   sci_ext      {"SCI",prompt="Name of science extension"}                # OLDP-3
char   var_ext      {"VAR",prompt="Name of variance extension"}               # OLDP-3
char   dq_ext       {"DQ",prompt="Name of data quality extension"}            # OLDP-3
bool   fl_vardq     {no,prompt="Create variance and data quality frames?"}    # OLDP-2
char   sat          {"default", prompt="Saturation level in raw images [ADU]"}  # OLDP-3
bool   verbose      {yes,prompt="Verbose output?"}                            # OLDP-4
char   logfile      {"",prompt="Name of logfile"}                             # OLDP-1
int    status       {0,prompt="Exit status (0=good)"}                         # OLDP-4

char   combine      {"average",prompt="Type of combine operation",enum="average|median"} # OLDP-3
char   reject       {"avsigclip",prompt="Type of rejection in flat average",enum="none|minmax|avsigclip|crreject"} # OLDP-3
real   lthreshold   {INDEF,prompt="Lower threshold when combining"}           # OLDP-3
real   hthreshold   {INDEF,prompt="Upper threshold when combining"}           # OLDP-3
int    nlow         {0,min=0,prompt="minmax: Number of low pixels to reject"} # OLDP-3
int    nhigh        {1,min=0,prompt="minmax: Number of high pixels to reject"} # OLDP-3
int    nkeep        {1,prompt="avsigclip: Minimum to keep (pos) or maximum to reject (neg)"} # OLDP-3
bool   mclip        {yes,prompt="avsigclip: Use median in clipping algorithm?"} # OLDP-3
real   lsigma       {3.,prompt="avsigclip: Lower sigma clipping factor"}       # OLDP-3
real   hsigma       {3.,prompt="avsigclip: Upper sigma clipping factor"}       # OLDP-3
real   sigscale     {0.1,prompt="avsigclip: Tolerance for clipping scaling corrections"} # OLDP-3
real   grow         {0.,prompt="minmax or avsigclip: Radius (pixels) for neighbor rejection"} # OLDP-3

char   gp_outpref   {"g",prompt="Gprepare prefix for output images"}           # OLDP-4
char   rawpath      {"",prompt="GPREPARE: Path for input raw images"}          # OLDP-4
char   key_ron      {"RDNOISE",prompt="Header keyword for readout noise"}      # OLDP-3
char   key_datasec  {"DATASEC",prompt="Header keyword for data section"}       # OLDP-3
char   gaindb       {"default",prompt="Database with gain data"} # OLDP-2
char   bpm          {"",prompt="Bad pixel mask"}                             # OLDP-2-input

char   gi_outpref   {"r",prompt="Gireduce prefix for output images"}           # OLDP-4
char   bias         {"",prompt="Bias calibration image"}                       # OLDP-1-input-fl_bias-required
bool   fl_over      {yes,prompt="Subtract overscan level?"}                    # OLDP-2
bool   fl_trim      {yes,prompt="Trim images?"}                                # OLDP-2
bool   fl_bias      {yes,prompt="Bias-subtract images?"}                       # OLDP-2
bool   fl_qecorr    {no, prompt="QE correct the input images?"}
bool   fl_inter     {no,prompt="Interactive overscan fitting?"}               # OLDP-3
char   nbiascontam  {"default", prompt="Number of columns removed from overscan region"}
char   biasrows     {"default", prompt="Rows to use for overscan region"}
char   key_biassec  {"BIASSEC",prompt="Header keyword for overscan image section"} # OLDP-3
bool   median       {no,prompt="Use median instead of average in column bias?"} # OLDP-3
char   function     {"chebyshev",prompt="Overscan fitting function.",enum="legendre|chebyshev|spline1|spline3"} # OLDP-3
char   order        {"default", prompt="Order of overscan fitting function."}           # OLDP-3
real   low_reject   {3.,prompt="Low sigma rejection factor."}                  # OLDP-3
real   high_reject  {3.,prompt="High sigma rejection factor."}                 # OLDP-3
int    niterate     {2,prompt="Number of rejection iterations."}             # OLDP-3
char    qe_data     {"gmosQEfactors.dat", prompt="Data file that contains QE information."}
char    qe_datadir  {"gmos$data/", prompt="Directory containg QE data file."}

struct  *scanfile   {"",prompt="Internal use only"}                            # OLDP-4

###############################################################################
begin

    char    l_inflats, l_outflat, l_bias, l_bpm, l_normsec
    char    l_sci_ext, l_var_ext, l_dq_ext, l_combine, l_sctype
    char    l_statsec, l_reject, l_biasrows
    char    l_key_biassec, l_key_datasec, l_gp_outpref, l_key_ron, l_key_gain
    char    l_gi_outpref, l_function, l_logfile, l_gaindb,l_rawpath,l_sat
    char    l_nbiascontam, l_qe_data, l_qe_datadir, l_order
    bool    l_fl_vardq, l_fl_scale, l_mclip, l_fl_inter, l_median
    bool    l_verbose, l_fl_over, l_fl_bias, l_fl_trim, l_fl_stamp, l_fl_qecorr
    real    l_lthreshold, l_hthreshold, l_lsigma, l_hsigma, l_sigscale, l_grow
    real    l_low_reject, l_high_reject
    int     l_nlow, l_nhigh, l_nkeep, l_niterate

    # Script variables
    int     nbad, nfiles, ii, jj, nccd, sx1, sx2, sy1, sy2
    int     nx1, nx2, ny1, ny2, nsci, satvalue, naxis[2]
    real    l_gains[12], l_gainsl[12], l_gainfl[12], l_gainfh[12], norm, stddev
    real    gain_norm, gain1, gain2, gain3, gain4, gain5, gain6
    real    gain7, gain8, gain9, gain10, gain11, gain12
    real    meanval, firstval, newgain, gainmode_threshold, readmode_slowval
    real    readmode_fastval, tnorm
    char    s_empty, img, trimsec, datasec, ccdsum[12], fccd[12]
    int     stat_ext, norm_ext
    char    keyfound
    char    l_readmode, ll_readmode, l_gainmode, ll_gainmode
    char    tmpfile, tmpinlist, combsig[12], combdq[12]
    char    tmpraw, tmpgprep, tmpgired, tmpvardqin
    char    obsm, masktyp, scaletype, dqlist
    bool    fl_gprepare, ocreate
    bool    ock, bck, fck, tck
    bool    ishamamatsu, ise2vDD
    struct  l_struct
    int     xbin, ybin, ext_to_use
    int     secx1, secx2, secy1, secy2
    int     secx1_bin, secx2_bin, secy1_bin, secy2_bin
    int     sx2_test_value, nx2_test_value, test_width, binned_test_width
    char    sec_to_use, extn
    real    amplifier_factor
    char    gaindata, tmpgtile, stat_file, tmpmask, tmp2mask, dqexp
    int     n_lines, l_maxfiles
    char    tmplist

    struct  l_datasec, pixels_to_use
    ###########################################################################

    cache ("imgets", "gimverify", "keypar", "fparse", "gemdate", "gtile")

    ###########################################################################
    # Set local variables
    ###########################################################################

    l_inflats = inflats
    l_outflat = outflat
    l_bias = bias
    l_bpm = bpm
    l_normsec = normsec
    l_fl_stamp = fl_stamp
    l_sci_ext = sci_ext
    l_var_ext = var_ext
    l_dq_ext = dq_ext
    l_fl_vardq = fl_vardq
    l_combine = combine
    l_fl_scale = fl_scale
    l_sctype = sctype
    l_statsec = statsec
    l_reject = reject
    l_lthreshold = lthreshold
    l_hthreshold = hthreshold
    l_nlow = nlow
    l_nhigh = nhigh
    l_nkeep = nkeep
    l_mclip = mclip
    l_lsigma = lsigma
    l_hsigma = hsigma
    l_sigscale = sigscale
    l_grow = grow
    l_key_biassec = key_biassec
    l_key_datasec = key_datasec
    l_gp_outpref = gp_outpref
    l_key_ron = key_ron
    l_key_gain = key_gain
    l_sat = sat
    l_gi_outpref = gi_outpref
    l_fl_inter = fl_inter
    l_median = median
    l_function = function
    l_nbiascontam = nbiascontam
    l_biasrows = biasrows
    l_order = order
    l_low_reject = low_reject
    l_high_reject = high_reject
    l_niterate = niterate
    l_verbose = verbose
    l_logfile = logfile
    l_gaindb = gaindb
    l_rawpath = rawpath
    l_fl_over = fl_over
    l_fl_bias = fl_bias
    l_fl_trim = fl_trim
    l_fl_qecorr = fl_qecorr
    l_qe_data = qe_data
    l_qe_datadir = qe_datadir

    ###########################################################################
    # Initialize variables
    ###########################################################################

    status = 0
    nccd = 3
    ock = no
    bck = no
    tck = no
    fck = no
    ishamamatsu = no
    ise2vDD = no
    ocreate = no

    # Maximum number of files gireduce can handle at any one time
    l_maxfiles = gireduce.maxfiles

    ###########################################################################
    # Make temporary files
    ###########################################################################

    tmpfile = mktemp ("tmpfile")
    tmpinlist = mktemp ("tmpinlist")
    tmpraw = mktemp ("tmpraw")
    tmpgprep = mktemp ("tmpgprep")
    tmpgired = mktemp ("tmpgired")
    tmplist = mktemp("tmplist")

    if (l_fl_vardq) { #don't' need these to be arrays, they are not in niflat
        for (ii=1; ii<=12; ii+=1) {
            combsig[ii] = mktemp ("tmpsig")
            combdq[ii] = mktemp ("tmpcombdq")
        }
    } else {
        for (ii=1; ii<=12; ii+=1) {
            combsig[ii] = ""
            combdq[ii] = ""
        }
    }

    ###########################################################################
    # Check string parameters
    ###########################################################################

    ###########################################################################
    # The logfile:
    ###########################################################################

    s_empty = ""
    print (l_logfile) | scan (s_empty)
    l_logfile = s_empty
    if (l_logfile == "") {
        l_logfile = gmos.logfile
        if (l_logfile == "") {
            l_logfile = "gmos.log"
            printlog ("WARNING - GIFLAT: both giflat.logfile and gmos.logfile \
                are empty.", l_logfile, verbose+)
            printlog ("                Using default file gmos.log.",
                l_logfile, verbose+)
        }
    }

    # Now start logging: (there are 76 dashes)
    date | scan(l_struct)
    printlog ("-------------------------------------------------------------\
        ---------------", l_logfile, l_verbose)
    printlog ("GIFLAT -- "//l_struct, l_logfile, l_verbose)
    printlog (" ", l_logfile, l_verbose)
    printlog ("Flat field images = "//l_inflats, l_logfile, l_verbose)
    printlog ("Output flat field image = "//l_outflat, l_logfile, l_verbose)
    printlog (" ", l_logfile, l_verbose)

    ###########################################################################
    # Check that output flat field image does not exist.
    ###########################################################################

    gimverify (l_outflat)
    if (gimverify.status != 1) {
        printlog ("ERROR - GIFLAT: Output flat field image ("//l_outflat//") \
            already exists.", l_logfile, verbose+)
        goto crash
    }
    l_outflat = gimverify.outname
    ocreate = yes

    ###########################################################################
    # Check for empty strings:
    ###########################################################################

    # SCI and DQ extensions -

    s_empty = ""
    print (l_sci_ext) | scan (s_empty)
    l_sci_ext = s_empty
    if (l_sci_ext == "") {
        printlog ("ERROR - GIFLAT: Science extension name sci_ext is missing",
            l_logfile, verbose+)
        goto crash
    }
    s_empty = ""
    print (l_dq_ext) | scan (s_empty)
    l_dq_ext = s_empty
    if (l_dq_ext == "") {
        printlog ("ERROR - GIFLAT: Data quality extension name dq_ext is \
            missing", l_logfile, verbose+)
        goto crash
    }

    ###########################################################################
    # Check input and output files:
    ###########################################################################

    # Input files

    if (substr(l_inflats,1,1) == "@") {
        img = substr (l_inflats, 2, strlen(l_inflats))
        if (access(img) == no) {
            printlog ("ERROR - GIFLAT: Input file "//img//" not found.",
                l_logfile, verbose+)
            goto crash
        }
    }

    files (l_inflats, sort-, > tmpfile)

    nbad = 0
    nfiles = 0
    for (ii=1; ii<=12 ; ii+=1) {
        fccd[ii] = ""
    }
    scanfile = tmpfile
    while (fscan (scanfile,img) != EOF) {
        # Check if image has not been gprepare'd - input could be
        # gprepare'd or we could find prefixed prepared file in
        # current dir - either way we use it.
        gimverify (l_gi_outpref//l_gp_outpref//img)
        if (gimverify.status == 0) {
            keypar (l_gi_outpref//l_gp_outpref//img//"[0]", "GIREDUCE",
                silent+)
            if (keypar.found) {
                img = gimverify.outname//".fits"
                print (gimverify.outname, >> tmpinlist)
                l_rawpath = ""
                printlog ("WARNING - GIFLAT: using previously gireduce'd \
                    image: "//img, l_logfile, l_verbose)
                nfiles += 1
                goto nextimage
            } else {
                printlog ("ERROR - GIFLAT: prefixed image exists: "//img//" \
                    but has not been gireduce'd.", l_logfile, l_verbose)
                goto crash
            }
        }
        gimverify (l_gp_outpref//img)
        if (gimverify.status == 0) {
            keyfound = ""
            hselect (l_gp_outpref//img//"[0]", "*PREPAR*", yes) | \
                scan (keyfound)
            if (keyfound != "") {
                img = gimverify.outname//".fits"
                print (gimverify.outname, >> tmpgired)
                l_rawpath = ""
                printlog ("WARNING - GIFLAT: using previously prepare'd \
                    image: "//img, l_logfile, l_verbose)
                nfiles += 1
                goto nextimage
            } else {
                printlog ("ERROR - GIFLAT: prefixed image exists: "//img//" \
                    but has not been gprepare'd.", l_logfile, l_verbose)
                goto crash
            }
        }

        # Check if image is "raw".
        gimverify (l_rawpath//img)
        if (gimverify.status == 0) {
            keypar (l_rawpath//img//"[0]", "GIFLAT", silent+)
            if (keypar.found) {
                printlog("WARNING - GIFLAT: Image "//img//" has been \
                    run through giflat before.", l_logfile, verbose+)
                goto nextimage
            }
            keypar (l_rawpath//img//"[0]", "GIREDUCE", silent+)
            if (keypar.found == no) {
                keyfound = ""
                hselect (l_rawpath//img//"[0]", "*PREPAR*", yes) | \
                    scan (keyfound)
                if (keyfound == "") {
                    print (img, >> tmpgprep)
                    nfiles +=1
                    goto nextimage
                } else {
                    printlog ("WARNING - GIFLAT: using previously prepare'd \
                        image: "//img, l_logfile, l_verbose)
                    print (img, >> tmpgired)
                    nfiles +=1
                    goto nextimage
                }
            } else {
                printlog ("WARNING - GIFLAT: using previously gireduce'd \
                    image: "//img, l_logfile, l_verbose)
                print (img, >> tmpinlist)
                nfiles +=1
                goto nextimage
            }
        } else if (gimverify.status == 1) {
            printlog ("ERROR - GIFLAT: Input file "//img//" not found.",
                l_logfile, verbose+)
            nbad += 1
            goto nextimage
        } else if (gimverify.status > 1) {
            printlog ("ERROR - GIFLAT: Input file "//img//" not a MEF FITS \
                image.", l_logfile, verbose+)
            nbad += 1
            goto nextimage
        }
nextimage:
    }

    # Check for empty file list
    if (nfiles == 0) {
        printlog ("ERROR - GIFLAT: No input images to process.", \
            l_logfile, verbose+)
        goto crash
    } else if ((nfiles==2) && (l_reject != "none") && (l_nlow > 0) \
        && (l_nhigh > 0)) {

        printlog ("WARNING - GIFLAT: Combining only 2 images with rejection.",
            l_logfile, verbose+)
        printlog ("WARNING - GIFLAT: Resetting nlow to 0. and nhigh to 1",
            l_logfile, verbose+)
        printlog ("WARNING - GIFLAT: However, the user should really review \
            the input parameters.", l_logfile, verbose+)
        l_nlow = 0
        # l_nhigh = 1
    }

    if (l_bias == "" || l_bias == " ")
        l_fl_bias=no

    # Check to ensure variance plane can be created
    # Ways variance is created: sigma from combining multiple images
    # if single image: to create variance using photon statistics in gireduce
    #                  need to either be bias subtracting or have previous bias
    #                  subtraction
    if (l_fl_vardq && (nfiles==1) && !l_fl_bias) {
        keypar (l_rawpath//img//"[0]", "BIASIM", silent+)
        if (keypar.found == no) {
            printlog("WARNING - GIFLAT: For single image, variance is created \
                using", l_logfile, verbose+)
            printlog("                  photon statistics requiring bias \
                subtraction.", l_logfile, verbose+)
            printlog("                  Image "//img//" has not been bias \
                subtracted.", l_logfile, verbose+)
            printlog("                  cannot create single variance. \
                Setting fl_vardq=no", l_logfile, verbose+)
            l_fl_vardq = no
        }
    }

    # Exit if problems found with input files
    if (nbad > 0) {
        printlog ("ERROR - GIFLAT: "//nbad//" image(s) either do not exist or \
            are not MEF files.", l_logfile, verbose+)
        goto crash
    } # end if (nbad > 0)

    if (access(tmpgprep)) {
        printlog ("Running gprepare on input files:", l_logfile, l_verbose)
        if (l_verbose)
            type (tmpgprep)
        type (tmpgprep, >> l_logfile)

        # RDNOISE and GAIN keywords -
        s_empty = ""
        print (l_key_ron) | scan (s_empty)
        l_key_ron = s_empty
        if (l_key_ron == "") {
            printlog ("WARNING - GIFLAT: Readout noise keyword parameter \
                key_ron is missing. Using default values.", l_logfile,
                verbose+)
        }
        s_empty = ""
        print (l_key_gain) | scan (s_empty)
        l_key_gain = s_empty
        if (l_key_gain == "") {
            printlog ("WARNING - GIFLAT: Gain keyword parameter key_gain is \
                missing. Using default values.", l_logfile, verbose+)
        }

        # Input BPM file -
        if (l_fl_vardq) {
            s_empty = ""
            print (l_bpm) | scan (s_empty)
            l_bpm = s_empty
            if (l_bpm == "") {
                printlog ("WARNING - GIFLAT: Bad pixel mask is an empty \
                    string", l_logfile, verbose+)
            } else {
                gimverify (l_bpm)
                l_bpm = gimverify.outname

                if (gimverify.status == 1) {
                    printlog ("ERROR - GIFLAT: Bad pixel mask "//l_bpm//" not \
                        found.", l_logfile, verbose+)
                    goto crash
                } else if (gimverify.status > 1) {
                    printlog ("ERROR - GIFLAT: Bad pixel mask "//l_bpm//" not \
                        a MEF FITS image.", l_logfile, verbose+)
                    goto crash
                } else {
                    imgets (l_bpm//"[0]", "GBPM", >& "dev$null")
                    if (imgets.value == "0") {
                        printlog ("WARNING - GIFLAT: Bad pixel mask "//\
                            l_bpm//" was not created with GMOS GBPM task.",
                            l_logfile, verbose+)
                    }
                }
            }
        }

        gprepare ("@"//tmpgprep, rawpath=l_rawpath, outimages="",
            fl_addmdf=no, logfile=l_logfile, outpref=l_gp_outpref,
            sci_ext=l_sci_ext, key_ron=l_key_ron, key_gain=l_key_gain,
            gaindb=l_gaindb, verbose=l_verbose)

        if (gprepare.status != 0)
            nbad += 1
        else {
            sections (l_gp_outpref//"@"//tmpgprep, option="fullname", \
                >> tmpgired)
            l_rawpath = "" #images are now in current directory
        }
    }

    if (access(tmpgired)) {
        if ((l_fl_bias || l_fl_trim || l_fl_over || l_fl_qecorr)) {

            # Input bias image -
            s_empty = ""
            print (l_bias) | scan (s_empty)
            l_bias = s_empty
            gimverify (l_bias)
            l_bias = gimverify.outname
            if (gimverify.status != 0 && l_fl_bias) {
                printlog ("ERROR - GIFLAT: Bias image ("//l_bias//") does not \
                    exist.", l_logfile, verbose+)
                goto crash
            }
            printlog ("Running gireduce on input files:", l_logfile, l_verbose)
            if (l_verbose)
                type (tmpgired)
            type (tmpgired, >> l_logfile)

            count (tmpgired) | scan (n_lines)
            while (n_lines > 0) {
                tail (tmpgired, nlines=n_lines) |
                    head ("STDIN", nlines=l_maxfiles, > tmplist)

                # call gireduce
                gireduce ("@"//tmplist, outpref=l_gi_outpref, outimages="",
                    fl_over=l_fl_over, fl_trim=l_fl_trim, fl_bias=l_fl_bias,
                    fl_dark=no, fl_qecorr=l_fl_qecorr, fl_flat=no, \
                    fl_vardq=l_fl_vardq, fl_addmdf=no,
                    bias=l_bias, dark="", flat1="", flat2="", flat3="",
                    flat4="", \
# MS - Commented out QE parameters removed from GIREDUCE in version 1.79
#                    qe_refim="", fl_keep_qeim=no, qe_corrimage="", \
                    qe_data=l_qe_data, qe_datadir=l_qe_datadir, \
                    key_biassec=l_key_biassec, key_datasec=l_key_datasec,
                    rawpath=l_rawpath, sci_ext=l_sci_ext, var_ext=l_var_ext,
                    dq_ext=l_dq_ext, bpm=l_bpm, gaindb=l_gaindb, sat=l_sat,
                    fl_mult=no, fl_inter=l_fl_inter, median=l_median,
                    function=l_function, nbiascontam=l_nbiascontam,
                    biasrows=l_biasrows, order=l_order, low_reject=l_low_reject,
                    high_reject=l_high_reject, niterate=l_niterate,
                    logfile=l_logfile, verbose=l_verbose)

                delete (tmplist, verify-, >& "dev$null")
                n_lines -= l_maxfiles

                if (gireduce.status != 0) {
                    printlog ("ERROR - GIFLAT: gireduce returned with error.",
                        l_logfile, verbose+)
                    goto crash
                }
            }

            sections (l_gi_outpref//"@"//tmpgired, option="fullname", \
               >> tmpinlist)
            l_rawpath = "" #images are now in current directory

        } else { # Image not gireduced and user still doesn't want to gireduce
            sections ("@"//tmpgired, option="fullname", >> tmpinlist)
        }
    }
    delete (tmpgprep, verify-, >& "dev$null")
    #delete (tmpgired, verify-, >& "dev$null")

    nbad = 0
    head (tmpinlist, nlines=1) | scan(img)

    imgets (l_rawpath//img//"[0]", "NSCIEXT", >& "dev$null")
    nsci = int(imgets.value)
    if (nsci == 0) {
        printlog ("ERROR - GIFLAT: NSCIEXT keyword not found", l_logfile,
            verbose+)
        goto crash
    }

    # Check if data is from Hamamatsu CCDs
    imgets (l_rawpath//img//"[0]", "DETTYPE", >& "dev$null")
    if (imgets.value == "S10892" || imgets.value == "S10892-N") {
        ishamamatsu = yes
        # Set hardcoded values for ll_readmode and ll_gainmode
        # GMOS-S
        gainmode_threshold = 3.
        readmode_slowval = 11880.
        readmode_fastval = 4000.
        # Set saturation level
        satvalue = 65000 #hcode ##M
    } else if (imgets.value == "SDSU II e2v DD CCD42-90") {
        # New e2vDD CCDs
        ise2vDD = yes
        # Set hardcoded values for ll_readmode and ll_gainmode #hcode
        gainmode_threshold = 3.
        readmode_slowval = 5000.
        readmode_fastval = 1000.
        # Set saturation level
        satvalue = 65000 # This is a default value that never gets used
                         ##M This needs to be fixed.
    } else if (imgets.value == "SDSU II CCD") { # Current EEV2 CCDs Jul-11
        # Set hardcoded values for ll_readmode and ll_gainmode #hcode
        gainmode_threshold = 3.
        readmode_slowval = 5000.
        readmode_fastval = 1000.
        # Set saturation level
        satvalue = 65000
    }

    # Set saturation level if not default
    if (l_sat != "default") {
        # Option for user-defined value
        satvalue = int(l_sat)
    }

    for (ii=1; ii<=nsci; ii+=1) {
        imgets (l_rawpath//img//"["//l_sci_ext//","//ii//"]", "CCDSUM",
            >& "dev$null")
        fccd[ii] = imgets.value
        ccdsum[ii] = fccd[ii]
    }
    count (tmpinlist) | scan (nfiles)
    delete (tmpfile, verify-, >& "dev$null")
    system.tail (tmpinlist, nlines=nfiles-1, > tmpfile)
    scanfile = tmpfile
    while (fscan (scanfile, img) != EOF) {

        for (ii=1; ii<=nsci; ii+=1) {
            imgets (l_rawpath//img//"["//l_sci_ext//","//ii//"]", "CCDSUM",
                >& "dev$null")
            ccdsum[ii] = imgets.value
            if (ccdsum[ii] != fccd[ii]) {
                printlog ("ERROR - GIFLAT: input images have different \
                    binning.", l_logfile, verbose+)
                nbad += 1
            }
        }
    } # end while loop
    # Exit if problems found with CCDSUM
    if (nbad > 0)
        goto crash

    printlog ("Using input files:", l_logfile, l_verbose)
    if (l_verbose)
        type (tmpinlist)
    type (tmpinlist, >> l_logfile)
    delete (tmpfile, verify-, >& "dev$null")

    ###########################################################################
    # Check l_statsec #hcode
    s_empty = ""
    print (l_statsec) | scan (s_empty)
    l_statsec = s_empty
    print (l_statsec) | ucase | scan (l_statsec)

    # Check l_normsec #hcode
    s_empty = ""
    print (l_normsec) | scan (s_empty)
    l_normsec = s_empty
    print (l_normsec) | ucase | scan (l_normsec)

    # Read the CCDSUM values into variables. Can use ccdsum[1] as all images
    # must have the same binning to get to this point. - MS
    # Initiate xbin and ybin
    xbin=1
    ybin=1

    # Scan the binning values in
    print (ccdsum[1]) | scanf ("%d %d", xbin, ybin)
    if (xbin != ybin){
        printlog ("WARNING - GIFLAT: Input flat images do not have square \
                binning.", l_logfile, verbose+)
    }

    # Initialise stat_ext and norm_ext
    stat_ext = 2
    norm_ext = 2

    # Set parameters for DEFAULT norm and stat secs
    if (l_fl_stamp) {
        secx1 = 28
        secx2 = 280 #140
        secy1 = 28
        secy2 = 280
        # Had to put this in to get around tests below with naxis1 -MS
        test_width = 300
    } else {
        secx1 = 100
        secx2 = 1800 # 900 (old 1900)
        secy1 = 100
        if (ishamamatsu) {
            # To account for the reduced number of pixels in y direction
            secy2 = 4000
        } else {
            secy2 = 4500
        }
        # Had to put this in to get around tests below with naxis1 - MS
        test_width = 2048
    }

    # Calculate default values to use for stat_ext, l_statsec, norm_ext and
    # l_normsec
    binned_test_width = nint(test_width / xbin)
    ext_to_use = 0 # Due to using gtile with fl_stats_only later on
    secx1_bin = nint(secx1 / xbin)
    secx2_bin = nint(secx2 / xbin)
    secy1_bin = nint(secy1 / ybin)
    secy2_bin = nint(secy2 / ybin)

    # Set up the default section to use
    pixels_to_use = "["//secx1_bin//":"\
        //secx2_bin//","//secy1_bin//":"//secy2_bin//"]"
    sec_to_use = "["//ext_to_use//"]"//pixels_to_use

    if (l_statsec == "DEFAULT") {
        l_statsec = sec_to_use
    } else if (substr(l_statsec,2,4) != l_sci_ext) {
        # Check user defined statsec to check extension type is correct
        printlog ("ERROR - GIFLAT: illegal parameter STATSEC \
            ("//l_statsec//") format is e.g [SCI,2][100:1900,100:4500]", \
            l_logfile, verbose+)
        goto crash
    }

    # Get norm_ext number, x & y range from user defined variable/default value
    # for statsec. Then check to make sure that they are valid - double checks
    # DEFAULT values
    if (l_statsec == sec_to_use) {
        print (l_statsec) | \
            scanf("[%d][%d:%d,%d:%d]", stat_ext, sx1, sx2, \
            sy1, sy2)

        # To sort out extension indexing - MS
        stat_ext = ext_to_use + 2
        if (l_fl_stamp) {
            stat_ext -= 1
        }
    } else {
        print (l_statsec) | \
            scanf("["//l_sci_ext//",%d][%d:%d,%d:%d]", stat_ext, sx1, sx2, \
            sy1, sy2)
    }

    # NAXIS1 and NAXIS2 values for comparison only needs to be done once
    if ((sx1 > sx2) || (sy1 > sy2) || (sx1 < 0) || (sy1 < 0) || \
        (stat_ext > nsci)) {
        printlog ("ERROR - GIFLAT: Illegal parameter STATSEC: "//\
            l_statsec, l_logfile, verbose+)
        goto crash
    }

    if (l_normsec == "DEFAULT") {
        l_normsec = sec_to_use
    } else if (substr(l_normsec,2,4) != l_sci_ext) {
        # Check user defined normsec to check extension type is correct
        printlog ("ERROR - GIFLAT: illegal parameter NORMSEC \
            ("//l_normsec//") format is e.g [SCI,2][100:1900,100:4500]", \
            l_logfile, verbose+)
        goto crash
    }

    # Get norm_ext number, x & y range from user defined variable/default value
    # for statsec. Then check to make sure that they are valid - double checks
    # DEFAULT values
    if (l_normsec == sec_to_use) {
        print (l_normsec) | \
            scanf("[%d][%d:%d,%d:%d]", norm_ext, nx1, nx2, \
            ny1, ny2)
        # To sort out extension indexing
        norm_ext = ext_to_use + 2
        if (l_fl_stamp) {
            norm_ext -= 1
        }
    } else {
        print (l_normsec) | \
            scanf("["//l_sci_ext//",%d][%d:%d,%d:%d]", norm_ext, nx1, nx2, \
            ny1, ny2)
    }

    if ((nx1 > nx2) || (ny1 > ny2) || (nx1 < 0) || (ny1 < 0) || \
        (norm_ext > nsci)) {
        printlog ("ERROR - GIFLAT: Illegal parameter NORMSEC: "//l_normsec,
            l_logfile, verbose+)
        goto crash
    }

    # Print out appropriate sectinos to be used -MS
    if (l_statsec == sec_to_use) {
        printlog ("Using STATSEC: "//pixels_to_use//" from CCD2", \
            l_logfile, l_verbose)
    } else {
        printlog ("Using STATSEC: "//l_statsec, l_logfile, l_verbose)
    }
    if (l_normsec == sec_to_use) {
        printlog ("Using NORMSEC: "//pixels_to_use//" from CCD2", \
            l_logfile, l_verbose)
    } else {
        printlog ("Using NORMSEC: "//l_normsec, l_logfile, l_verbose)
    }

    ###########################################################################
    # Create an empty MEF flat field -

    imcopy (l_rawpath//img//"[0]", l_outflat//".fits", verbose-)

    ###########################################################################
    # The input flats are in tmpinlist
    # The output flat is named: l_outflat
    ###########################################################################
    # Check that for flat field images OBSMODE=IMAGE and MASKTYPE=0
    ###########################################################################

    rename (tmpinlist, tmpfile, field="all")
    scanfile = tmpfile
    l_readmode = "slow"  # default readmode
    l_gainmode = "low"   # default gainmode
    ii = 1
    while (fscan(scanfile ,img) != EOF ) {
        imgets (l_rawpath//img//"[0]", "OBSMODE", >& "dev$null")

        obsm = imgets.value
        imgets (l_rawpath//img//"[0]", "MASKTYP", >& "dev$null")
        masktyp = imgets.value
        if ((obsm == "IMAGE") && (masktyp == "0"))
            print (img, >> tmpinlist)
        else
            printlog ("WARNING - GIFLAT: Input flat image "//img//" is not an \
                imaging flat.", l_logfile, verbose+)

        imgets (l_rawpath//img//"[0]", "AMPINTEG", >& "dev$null") #hcode
        if (real(imgets.value) == readmode_slowval)
            ll_readmode = "slow" #hcode
        if (real(imgets.value) == readmode_fastval)
            ll_readmode = "fast" #hcode
        if (ii == 1) {
            l_readmode = ll_readmode
        } else if (l_readmode != ll_readmode) {
            printlog ("ERROR - GIFLAT: Cannot mix slow and fast readout \
                images", l_logfile, verbose+)
            goto crash
        }
        imgets (l_rawpath//img//"[SCI,1]", "GAIN", >& "dev$null") #hcode
        if (real(imgets.value) > gainmode_threshold) #hcode
            ll_gainmode = "high"
        if (real(imgets.value) < gainmode_threshold) #hcode
            ll_gainmode = "low"
        if (ii == 1) {
            l_gainmode = ll_gainmode
        } else if (l_gainmode != ll_gainmode) {
            printlog ("ERROR - GIFLAT: Cannot mix high and low gain images",
                l_logfile, verbose+)
            goto crash
        }
        ii += 1
    }
    delete (tmpfile, verify-, >& "dev$null")
    if (access(tmpinlist) == no) {
        printlog ("ERROR - GIFLAT: No flat field images available.",
            l_logfile, verbose+)
        goto crash
    }

    # Check the image size: statsec and normsec should be within the image
    # (lower range already checked to be greater than zero)
    head (tmpinlist, nline=1) | scan(img)
    keyfound = ""
    hselect (l_rawpath//img//"["//l_sci_ext//","//str(stat_ext)//"]", \
        "i_naxis1", yes) | scan (keyfound)
    naxis[1] = int(keyfound)
    keyfound = ""
    hselect (l_rawpath//img//"["//l_sci_ext//","//str(stat_ext)//"]", \
        "i_naxis2", yes) | scan (keyfound)
    naxis[2] = int(keyfound)

    if (l_statsec == sec_to_use) {
        sx2_test_value = binned_test_width
    } else {
        sx2_test_value = naxis[1]
    }

    if ((sx2 > sx2_test_value) || (sy2 > naxis[2])) {
        print (sx2//" "//sx2_test_value//" "//sy2)
        printlog ("ERROR - GIFLAT: Statistics section is outside the image",\
            l_logfile, verbose+)
        goto crash
    }

    keyfound = ""
    hselect (l_rawpath//img//"["//l_sci_ext//","//norm_ext//"]", "i_naxis1", \
        yes) | scan (keyfound)
    naxis[1] = int(keyfound)
    keyfound = ""
    hselect (l_rawpath//img//"["//l_sci_ext//","//norm_ext//"]", "i_naxis2", \
        yes) | scan (keyfound)
    naxis[2] = int(keyfound)

    if (l_normsec == sec_to_use) {
        nx2_test_value = binned_test_width
    } else {
        nx2_test_value = naxis[1]
    }

    if ((nx2 > nx2_test_value) || (ny2 > naxis[2])) {
        printlog ("ERROR - GIFLAT: Normalization section is outside the image",
            l_logfile, verbose+)
        goto crash
    }

    ###########################################################################
    # Combine reduced images - they are in tmpinlist
    ###########################################################################

    scaletype = "none"
    count (tmpinlist) | scan (nfiles)
    if (l_fl_scale && (nfiles > 1)) {
        scaletype = "!RELINT"
        head (tmpinlist, nline=1) | scan(img)

        if (l_statsec == sec_to_use) {
            tmpgtile = mktemp("tmpgtile")

            # Tile the requested extension
            gtile (inimages=l_rawpath//img, outimages=tmpgtile, out_ccds="2", \
                fl_stats_only=yes, fl_tile_det=no, sci_ext=l_sci_ext, \
                var_ext="VAR", dq_ext="DQ", mdf_ext="MDF", \
                key_detsec="DETSEC", key_ccdsec="CCDSEC", \
                key_datasec="DATASEC", key_biassec="BIASSEC", \
                key_ccdsum="CCDSUM", fl_verbose=l_verbose, logfile=l_logfile,\
                >& "dev$null")

            if (gtile.status != 0) {
                printlog ("ERROR - GIFLAT: GTILE retunred a non-zero status.",\
                    l_logfile, verbose+)
                goto crash
            }

            stat_file = tmpgtile//l_statsec
        } else {
            stat_file = l_rawpath//img//l_statsec
        }

        firstval = INDEF

        imstatistics (stat_file, field=l_sctype, lower=INDEF,
            upper=INDEF, nclip=0, lsigma=INDEF, usigma=INDEF, binwidth=0.1,
            format-, cache-) | scan (firstval)

        if (isindef(firstval)) {
            printlog ("ERROR - GIFLAT: Cannot determine "//l_sctype//" for "//\
                img, l_logfile, verbose+)
            goto crash
        }

        if (l_statsec == sec_to_use) {
            delete (tmpgtile//"*.fits", verify-, >& "dev$null")
        }

        for (ii=1; ii<=nsci; ii+=1) {
            gemhedit (l_rawpath//img//"["//l_sci_ext//","//ii//"]", "RELINT", \
                1., "", delete-)
        }
        count (tmpinlist) | scan (ii)
        system.tail (tmpinlist, nlines=ii-1, > tmpfile)
        scanfile = tmpfile
        meanval = 1.0
        while (fscan(scanfile,img) != EOF) {

            if (l_statsec == sec_to_use) {
                tmpgtile = mktemp("tmpgtile")

                # Tile the requested extension
                gtile (inimages=l_rawpath//img, outimages=tmpgtile, \
                    out_ccds="2", \
                    fl_stats_only=yes, fl_tile_det=no, sci_ext=l_sci_ext, \
                    var_ext="VAR", dq_ext="DQ", mdf_ext="MDF", \
                    key_detsec="DETSEC", key_ccdsec="CCDSEC", \
                    key_datasec="DATASEC", key_biassec="BIASSEC", \
                    key_ccdsum="CCDSUM", fl_verbose=l_verbose, \
                    logfile=l_logfile, >& "dev$null")
                stat_file = tmpgtile//l_statsec
            } else {
                stat_file = l_rawpath//img//l_statsec
            }

            imstatistics (stat_file, field=l_sctype,
                lower=INDEF, upper=INDEF, nclip=0, lsigma=INDEF, usigma=INDEF,
                binwidth=0.1, format-, cache-) | scan(meanval)
            meanval = firstval / meanval
            for (ii=1; ii<=nsci; ii+=1) {
                gemhedit (l_rawpath//img//"["//l_sci_ext//","//ii//"]",
                    "RELINT", meanval, "", delete-)
            }
            if (l_statsec == sec_to_use) {
                delete (tmpgtile//"*.fits", verify-, >& "dev$null")
            }
        }
        scanfile = ""
        delete (tmpfile, verify-, >& "dev$null")
    }

    # Run imcombine through all nsci science extensions
    for (ii=1; ii<=nsci; ii+=1) {
        l_gains[ii] = 0.0

        # Place all [sci,ii] extensions in file
        delete (tmpfile, verify-, >& "dev$null")
        scanfile = tmpinlist
        while (fscan(scanfile,img) != EOF) {
            print (l_rawpath//img//"["//l_sci_ext//","//ii//"]", >> tmpfile)
            if (l_gains[ii] == 0.0) {
                imgets (l_rawpath//img//"["//l_sci_ext//","//ii//"]",
                    l_key_gain, >& "dev$null")
                l_gains[ii] = real(imgets.value)
            }
        }
        flpr
        ##M This should be a gemcombine call
        imcombine ("@"//tmpfile,
            l_outflat//"["//l_sci_ext//","//ii//",append]", headers="",
            bpmasks="", rejmasks="", nrejmasks=combdq[ii], expmasks="",
            sigmas=combsig[ii], logfile=l_logfile, combine=l_combine,
            reject=l_reject, project-, outtype="real", outlimits="",
            offsets="none", masktype="none", maskvalue=0., blank=0.,
            scale=scaletype, zero="none", weight="none", statsec="",
            lthres=l_lthreshold, hthres=l_hthreshold, nlow=l_nlow,
            nhigh=l_nhigh, nkeep=l_nkeep, mclip=yes, lsigma=l_lsigma,
            hsigma=l_hsigma, sigscale=l_sigscale, grow=l_grow)

       if (imaccess(l_outflat//"["//l_sci_ext//","//ii//"]") == no) {
            printlog ("ERROR - GIFLAT: Unable to combine images for extension \
                ["//l_sci_ext//","//ii//"]", l_logfile, verbose+)
            goto crash
        }
        delete (tmpfile, verify-, >& "dev$null")

        # Change intensities to electrons - will do the var into electrons
        # conversion later.

        imarith (l_outflat//"["//l_sci_ext//","//ii//"]", "*", l_gains[ii],
            l_outflat//"["//l_sci_ext//","//ii//",overwrite]", title="",
            hparams="", pixtype="", calctype="", verbose-, noact-)
        imgets (l_outflat//"["//l_sci_ext//","//ii//"]", "NCOMBINE")
        gemhedit (l_outflat//"["//l_sci_ext//","//ii//"]", "GAIN",
            real(imgets.value), "", delete-)
        gemhedit (l_outflat//"["//l_sci_ext//","//ii//"]", "GAINORIG",
            l_gains[ii], "Orig. Gain [e-/ADU]", delete-)
    }

    # Get the normalization from the center CCD
    keypar (l_outflat//"["//l_sci_ext//","//norm_ext//"]", "GAINORIG", silent+)
    newgain = real(keypar.value)

    keypar (l_outflat//"["//l_sci_ext//","//norm_ext//"]","OVERSCAN", \
        silent+)
    if (keypar.found) {
        satvalue = (satvalue - int(keypar.value)) * newgain
    } else
        satvalue = satvalue * newgain

    if (l_normsec == sec_to_use) {
        tmpgtile = mktemp("tmpgtile")

        # Tile the requested extension
        gtile (inimages=l_outflat, outimages=tmpgtile, out_ccds="2", \
            fl_stats_only=yes, fl_tile_det=no, sci_ext=l_sci_ext, \
            var_ext="VAR", dq_ext="DQ", mdf_ext="MDF", \
            key_detsec="DETSEC", key_ccdsec="CCDSEC", \
            key_datasec="DATASEC", key_biassec="BIASSEC", \
            key_ccdsum="CCDSUM", fl_verbose=l_verbose, >& "dev$null")

        stat_file = tmpgtile//l_normsec
        # For e2vDD CCD mask
        # Parse the extension requested
        extn = substr(l_normsec,1,stridx("]",l_normsec))

    } else {
        stat_file = l_outflat//l_normsec
        # For e2vDD CCD mask
        # Parse the extension requested
        extn = substr(l_normsec,1,stridx("]",l_normsec))
    }

    if (ise2vDD) {
        tmpmask = mktemp("tmpmask")//".fits"
        fxcopy (l_outflat, output=tmpmask, groups="0", verbose=no)

        # Create new saturation mask - e2vDD CCDs ##M Need to use in DQ too
        dqexp = "((a >= b) ? 1 : 0)"
        for (ii=1; ii<=nsci; ii+=1) {
            if (l_sat == "default") {
                # Obtain the saturation for this extension
                gsat (l_outflat, extension="["//l_sci_ext//\
                    ","//ii//"]", gaindb=l_gaindb, bias="static", \
                    pixstat="midpt", statsec="default", \
                    gainval="default", logfile=l_logfile, \
                    verbose=no)
                if (gsat.status != 0) {
                    printlog ("ERROR - GIREDUCE: GSAT returned a "//\
                        "non-zero status. Exiting", \
                        l_logfile, verbose+)
                    status = 1
                    goto crash
                }
                satvalue = -999999.0
                satvalue = real(gsat.saturation)
                if (satvalue == -999999.0 || satvalue == INDEF) {
                    printlog ("ERROR - GIFLAT: GSAT calculated an invalid "//\
                        "saturation value.", l_logfile, verbose+)
                }
            }

            # Create mask
            imexpr (dqexp, tmpmask//"[SCI,"//ii//",append]", \
                l_outflat//"[SCI,"//ii//"]", satvalue,  dims="auto", \
                intype="auto", outtype="short", refim="auto", \
                bwidth=0, btype="nearest", bpixval=0., rangecheck+, \
                verbose=no, exprdb="none", lastout="", mode="ql")
             ##M Add masks together if BPM present?
            if (!access(tmpmask)) {
                printlog ("ERROR - GIFLAT: Can't access "//tmpmask, \
                    l_logfile, verbose+)
                status = 1
                goto crash
            }
        }
        if (l_normsec == sec_to_use) {
            tmp2mask = mktemp("tmp2mask")//".fits"
            # Tile up the sat mask
            gtile (inimages=tmpmask, outimages=tmp2mask, out_ccds="2", \
                fl_stats_only=yes, fl_tile_det=no, sci_ext=l_sci_ext, \
                var_ext="VAR", dq_ext="DQ", mdf_ext="MDF", \
                key_detsec="DETSEC", key_ccdsec="CCDSEC", \
                key_datasec="DATASEC", key_biassec="BIASSEC", \
                key_ccdsum="CCDSUM", fl_verbose=l_verbose, logfile=l_logfile, \
                >& "dev$null")

            if (!access(tmp2mask)) {
                printlog ("ERROR - GIFLAT: Can't access "//tmp2mask, \
                    l_logfile, verbose+)
                status = 1
                goto crash
            }

        } else {
            tmp2mask = tmpmask
        }
        tmp2mask = tmp2mask//extn
    }

    if (ise2vDD) {
        norm = -99999.0
        stddev = -99999.0
        mimstatistics (stat_file, imasks=tmp2mask, field="mode,stddev", \
            lower=0., upper=INDEF, lsigma=INDEF, usigma=INDEF, \
            binwidth=0.1, format-, cache-) | scan (norm, stddev)

        if (norm == -99999.0 || stddev == -99999.0) {
            printlog ("ERROR - GIFLAT: norm or stddev values not "//\
                "determined for "//stat_file, l_logfile, verbose+)
            status = 1
            goto crash
        }

        tnorm = norm
        mimstatistics (stat_file, imasks=tmp2mask, field="midpt", \
            lower=(norm-3*stddev),\
            upper=(norm+3*stddev), lsigma=INDEF, usigma=INDEF, binwidth=0.1, \
            format-, cache-) | scan (norm)

        if (norm == tnorm) {
            printlog ("ERROR - GIFLAT: normalization not determined for "//\
                stat_file, l_logfile, verbose+)
            status = 1
            goto crash
        }

        delete (tmpmask, verify-, >& "dev$null")
        if (l_normsec == sec_to_use) {
            fparse (tmp2mask)
            imdelete (fparse.directory//fparse.root//fparse.extension, \
                verify-, >& "dev$null")
        }
    } else {
        ## This will need updating for Hamamatsu data
        imstatistics (stat_file, field="mode,stddev", lower=0.,
            upper=satvalue, lsigma=INDEF, usigma=INDEF, binwidth=0.1,
            format-, cache-) | scan (norm, stddev)

        imstatistics (stat_file, field="midpt", lower=(norm-3*stddev),\
            upper=(norm+3*stddev), lsigma=INDEF, usigma=INDEF, binwidth=0.1, \
            format-, cache-) | scan (norm)
    }

    printlog ("GIFLAT: Normalization factor "//norm, l_logfile, l_verbose)

    if (l_normsec == sec_to_use) {
        delete (tmpgtile//"*.fits", verify-, >& "dev$null")
    }

    # Normalize the flats, change var into electrons and normalize
    # add dq
    head (tmpinlist, nlines=1) | scan(img)

    for (ii=1; ii<=nsci; ii+=1) {
        imarith (l_outflat//"["//l_sci_ext//","//ii//"]", "/", norm,
            l_outflat//"["//l_sci_ext//","//ii//",overwrite]", title="",
            hparams="", pixtype="", calctype="", verbose-, noact-)
    }


    if (l_fl_vardq && nfiles==1) { # copy poisson variance created in gireduce
        if (access(tmpgired)) { # image was gireduced here
            type (tmpgired) | scan(tmpvardqin)
            tmpvardqin = l_gi_outpref//tmpvardqin
        } else if (access(tmpinlist)) { # used previously gireduced image
            type (tmpinlist) | scan(tmpvardqin)
        }

        for (ii=1;ii<=nsci;ii+=1) {
            if (imaccess(tmpvardqin//"["//l_var_ext//","//ii//"]")) {
            # divide by norm squared, multiply by gain squared,
            # and add tmpvardqin to final flat
                imexpr ("a* (c**2)/(b**2)",
                    l_outflat//"["//l_var_ext//","//ii//",append]", \
                    tmpvardqin//"["//l_var_ext//","//ii//"]", norm, \
                    l_gains[ii], dims="auto", intype="real", \
                    outtype="real", refim="auto", bwidth=0, btype="nearest", \
                    bpixval=0., rangecheck+, verbose-, exprdb="none")
                imcopy (tmpvardqin//"["//l_dq_ext//","//ii//"]", \
                    l_outflat//"["//l_dq_ext//","//ii//",append]", verbose-)
            } else {
                printlog ("WARNING - GIFLAT: Cannot access "//tmpvardqin// \
                    "["//l_var_ext//","//ii//"]",l_logfile, l_verbose)
                printlog ("                  Output will not have var/dq", \
                    l_logfile, l_verbose)
            }
        }
    } else if (l_fl_vardq && nfiles > 1) {
        # using combsig[ii] from imcombine
        # Divide sigma by normalization squared - multiply by gain squared
        # Square the sigma file to get the variance
        for (ii=1; ii<=nsci; ii+=1) {
            imexpr ("(a*c/b)**2",
                l_outflat//"["//l_var_ext//","//ii//",append]", \
                combsig[ii], norm, l_gains[ii], dims="auto", intype="real", \
                outtype="real", refim="auto", bwidth=0, btype="nearest", \
                bpixval=0., rangecheck+, verbose-, exprdb="none")
            # Add bad pixel masks
            ##M this is a quick hack and needs to be rewritten to use
            ##M gemcombine as mentioned above. Also add fl_dqprop flag
            imexpr ("a == "//nfiles//" ? 1 : 0", \
                l_outflat//"["//l_dq_ext//","//ii//"append]", \
                combdq[ii]//".pl", dims="auto", intype="auto", \
                outtype="ushort", refim="auto", bwidth=0, \
                btype="nearest", bpixval=0., rangecheck=yes,\
                verbose=no)
            delete (combdq[ii] // ".pl", verify-, >& "dev$null")
        }
    }

    # Now update the PHU
    gemdate ()
    gemhedit (l_outflat//"[0]", "GIFLAT", gemdate.outdate,
        "UT Time stamp for GIFLAT", delete-)
    gemhedit (l_outflat//"[0]", "GIFLNORM", norm,
        "Normalization factor for GIFLAT", delete-)
    gemhedit (l_outflat//"[0]", "GEM-TLM", gemdate.outdate,
        "UT Last modification with GEMINI", delete-)
    goto clean

    ###########################################################################
    # Exit with error
    ###########################################################################
crash:
    status = 1
    if (ocreate)
        imdelete (l_outflat, verify-, >& "dev$null")
    goto clean

    ###########################################################################
    # Clean up and exit
    ###########################################################################
clean:
    delete (tmpfile//","//tmpinlist, verify-, >& "dev$null")
    delete (tmpgprep, verify-, >& "dev$null")
    delete (tmpgired, verify-, >& "dev$null")

    for (ii=1; ii<=12; ii+=1) {
        delete (combsig[ii]//".fits", verify-, >& "dev$null")
        delete (combdq[ii]//".fits", verify-, >& "dev$null")
        if (access(combdq[ii]//".pl")) {
            delete (combdq[ii] // ".pl", verify-, >& "dev$null")
        }
    }

    if (status == 0) {
        printlog ("", l_logfile, l_verbose)
        printlog ("GIFLAT exit status: good.", l_logfile, l_verbose)
        printlog ("----------------------------------------------------------\
            ------------------", l_logfile, l_verbose)
    } else {
        printlog ("", l_logfile, l_verbose)
        printlog ("GIFLAT exit status: error.", l_logfile, l_verbose)
        printlog ("-----------------------------------------------------------\
            -----------------", l_logfile, l_verbose)
    }

end

# Copyright(c) 2002-2017 Association of Universities for Research in Astronomy, Inc.

procedure gireduce (inimages)

# Basic reductions of GMOS imaging data
#
# Version Feb 28, 2002  MT,IJ,BM  v1.3 release
#         Aug 16, 2002  IJ debug version: "Cannot coerce string `' to real" bug (Linux bug)
#         Aug 26, 2002  IJ parameter encoding
#         Sep 20, 2002  IJ v1.4 release
#         Oct 31, 2002  IJ overvar[6,200] bug fixed
#         Feb 12, 2003  BM Use imgets to read GAIN instead of keypar to avoid problem on Linux
#         Mar  3, 2003  BM User imgets instead of keypar to read TRIMSEC, was the cause of segmentation fault on some Linux systems
#         Jun 11, 2003  BM Generalize biassec modifications
#         Aug 26, 2003  KL IRAF2.12 - new parameters
#                          hedit: addonly;  imstat: nclip,lsigma,usigma,cache
#         Oct 17, 2003  KL Fix invalid array index when fnum==0

char    inimages    {prompt="Input GMOS images"}                # OLDP-1-input-primary-single-prefix=r
char    outpref     {"r", prompt="Prefix for output images"}        # OLDP-4
char    outimages   {"", prompt="Output images"}                    # OLDP-1-output
bool    fl_over     {yes, prompt="Subtract overscan level"}         # OLDP-3
bool    fl_trim     {yes, prompt="Trim off the overscan section"}   # OLDP-2
bool    fl_bias     {yes, prompt="Subtract bias image"}             # OLDP-2
bool    fl_dark     {no, prompt="Subtract (scaled) dark image"}     # OLDP-3
bool    fl_qecorr   {no, prompt="QE correct the input images?"}
bool    fl_flat     {yes, prompt="Do flat field correction?"}       # OLDP-2
bool    fl_vardq    {no, prompt="Create variance and data quality frames"}  # OLDP-3
bool    fl_addmdf   {no, prompt="Add Mask Definition File? (LONGSLIT/MOS/IFU modes)"} # OLDP-3
char    bias        {"", prompt="Bias image name"}                  # OLDP-1-input-fl_bias-required
char    dark        {"", prompt="Dark image name"}                  # OLDP-1-input-fl_dark-required
char    flat1       {"", prompt="Flatfield image 1"}                # OLDP-1-input-fl_flat-required
char    flat2       {"", prompt="Flatfield image 2"}                # OLDP-1-input-fl_flat-optional
char    flat3       {"", prompt="Flatfield image 3"}                # OLDP-1-input-fl_flat-optional
char    flat4       {"", prompt="Flatfield image 4"}              # OLDP-1-input-fl_flat-optional
char    qe_refim    {"", prompt="QE wavelength reference image (spectroscopy only)"}
bool    fl_keep_qeim {no, prompt="Keep QE correction?"}
char    qe_corrpref {"qecorr", prompt="Prefix for QE correction files"}
char    qe_corrimages {"", prompt="Name for QE correction data"}
char    qe_data     {"gmosQEfactors.dat", prompt="Data file that contains QE information"}
char    qe_datadir  {"gmos$data/", prompt="Directory containg QE data file"}
char    key_exptime {"EXPTIME", prompt="Header keyword of exposure time"}   # OLDP-3
char    key_biassec {"BIASSEC", prompt="Header keyword for bias section"}   # OLDP-3
char    key_datasec {"DATASEC", prompt="Header keyword for data section"}   # OLDP-3
char    rawpath     {"", prompt="GPREPARE: Path for input raw images"}      # OLDP-4
char    gp_outpref  {"g", prompt="GPREPARE: Prefix for output images"}      # OLDP-4

char    sci_ext     {"SCI", prompt="Name of science extension"}     # OLDP-3
char    var_ext     {"VAR", prompt="Name of variance extension"}    # OLDP-3
char    dq_ext      {"DQ", prompt="Name of data quality extension"} # OLDP-3
char    key_mdf     {"MASKNAME", prompt="Header keyword for the Mask Definition File"} # OLDP-3
char    mdffile     {"", prompt="MDF file to use if keyword not found"}     # OLDP-3-xinput
char    mdfdir      {"gmos$data/", prompt="MDF database directory"} # OLDP-3

char    bpm         {"", prompt="Bad pixel mask"}                   # OLDP-2-input
char    gaindb      {"default", prompt="Database with gain data"} # OLDP-3
char    sat         {"default", prompt="Saturation level in raw images [ADU]"}  # OLDP-3

char    key_nodcount {"NODCOUNT", prompt="Header keyword with number of nod cycles"}
char    key_nodpix   {"NODPIX", prompt="Header keyword with shuffle distance"}


char    key_filter  {"FILTER2", prompt="Header keyword of filter"}  # OLDP-3
char    key_ron     {"RDNOISE", prompt="Header keyword for readout noise"}  # OLDP-3
char    key_gain    {"GAIN", prompt="Header keyword for gain [e-/ADU]"}     # OLDP-3
real    ron         {3.5, prompt="Readout noise in electrons"}      # OLDP-3
real    gain        {2.2, prompt="Gain in e-/ADU"}                  # OLDP-3
bool    fl_mult     {yes, prompt="Multiply by gains to get output in electrons"}

bool    fl_inter    {no, prompt="Interactive overscan fitting?"}    # OLDP-3
bool    median      {no, prompt="Use median instead of average in column bias?"} # OLDP-3
char    function    {"chebyshev", prompt="Overscan fitting function",enum="spline3|legendre|chebyshev|spline1"}  # OLDP-3
char    nbiascontam {"default", prompt="Number of columns removed from overscan region"}
string  biasrows    {"default", prompt="Rows to use for overscan region"}
char    order       {"default", prompt="Order of overscan fitting function"}# OLDP-3
real    low_reject  {3., prompt="Low sigma rejection factor in overscan fit"}    # OLDP-3
real    high_reject {3., prompt="High sigma rejection factor in overscan fit"}   # OLDP-3
int     niterate    {2, prompt="Number of rejection iterations in overscan fit"} # OLDP-3

int     maxfiles    {100, min=100, max=100, prompt="For external access only. Maximum number of input files allowed"}
char    logfile     {"", prompt="Logfile"}                          # OLDP-1
bool    verbose     {yes, prompt="Verbose?"}                        # OLDP-4
int     status      {0, prompt="Exit status (0=good)"}              # OLDP-4
char    *scanfile   {"", prompt="For internal use only"}            # OLDP-4

###############################################################################

begin

    # Variable declaration
    char    l_inimages, l_outpref, l_outimages, l_sci_ext, l_var_ext, l_dq_ext
    char    l_key_biassec, l_function, l_key_datasec, l_bias, l_dark
    char    l_key_filter,l_key_exptime, l_key_nodpix, l_key_nodcount
    char    l_flat[4], l_gp_outpref, l_bpm, l_key_ron, l_key_gain, l_rawpath
    char    l_key_mdf, l_mdffile, l_mdfdir, l_logfile, l_gaindb, l_biasrows

    bool    l_fl_vardq, l_fl_over, l_median, l_fl_inter, l_fl_mult
    bool    l_fl_trim, l_fl_bias, l_fl_dark, l_fl_flat, l_fl_addmdf, l_verbose
    int     l_niterate
    real    l_low_reject, l_high_reject, l_ron, l_gain, satfrac

    char    l_key_dsec, l_sat, l_nbiascontam, l_order
    char    varexp, dqexp, tmpd[50], tmpdq
    char    sci_tsec, bpm_tsec, sci_dsec, bpm_dsec, tsec
    #char    obs_bpm
    int     next_bpm

    # General variables:
    char    img, s_empty, input[100], output[100]
    char    tmpinlist, tmpoutlist, tmplog, tmpfile
    char    obsmodeflat[4]
    char    l_biassec, l_datasec, tsec_img, tsec_bias[60], tsec_flat[60,4]
    char    tsec_dark[60], l_nodpix, l_nodcount
    char    thisfilter, flatfilter[4]
    char    srms, subflat, subbias, subdark
    char    junk, keyfound, incrpix1, l_bname, l_dname, l_fname[4]
    int     naxis1, ncheck, checkcount, osc_order
    bool    extnum_mismatch

    # For removing 48 rows at bottom of new Hamamatsu CCDs
    char    l_datasec_new, l_ccdsec, l_detsec, l_ccdsec_new, l_detsec_new
    int     datax1, datax2, datay1, datay2
    int     row_removal_value, row_removal_value_binned, row_removal_value_new
    int     xbin, ybin
    int     detsecx1, detsecx2, detsecy1, detsecy2
    int     ccdsecx1, ccdsecx2, ccdsecy1, ccdsecy2

    # QE correction variables
    char    l_qe_refim, l_qe_corrpref, qecorrapp, l_qecorr_data
    char    l_qecorr_datadir, l_qe_corrimages, tmpqecorrimgs
    char    files_refimage[100], filename, existed, ll_corrimg[100]
    bool    ll_corrimg_exists[100], inqecorr_state[100], flatqecorr_state[4]
    bool    l_fl_qecorr, l_fl_keep_qeim, canqecorr
    bool    isimaging, globalcanqecorr
    string  ll_qecorr, ll_varqecorr , ll_qecorr_int, ll_qecorr_int1
    string  ccdname, l_key_qeimfactor, l_key_qecorrim
    string  l_qeimfactor_comment, l_qecorrim_comment, qecorrmssg
    int     loopcounter, counter_refimage, prev_qecorr, num_qecorrimgs

    # local image names and expressions
    char    ll_bias, ll_flat, ll_dark, ll_sciexp, ll_varexp, ll_dqexp
    char    ll_varbias, ll_varflat, ll_vardark, ll_dqbias, ll_dqflat, ll_dqdark
    int     ll_ndq, nbiascontamvalue
    real    l_darkexptime, l_darkscale, l_darknodpix, l_darknodcount, satvalue
    char    s_bias, s_flat, s_dark  # for printlog only
    real    n_gain[60], overvar[60,100], dummy
    char    tmpmdf, tmptrim
    char    dqlist, tmppl

    int     nbad, nfiles, nimages, ii, jj, kk, biasx1, biasx2, biasy1, biasy2
    int     n_extflat[4], n_extbias, n_ext[100], n_extdark
    int     fnum, l_maxfiles, atposition, lower_bias_y
    char    ccd_flat[4,60], ccd_bias[60], ccd_dark[60], ccd[100,60]
    char    flatmssg, biasmssg, darkmssg, pathtest, rphend, testfile, img_root

    bool    l_fl_raw, overstate_cal, overstate_img
    bool    canbias, candark, canflat[4], n_do[100]
    bool    hasmdf[100]
    bool    fl_tsec_img, fl_tsec_bias[60], fl_tsec_flat[60,4], badflag
    bool    fl_tsec_dark[60]
    bool    ishamamatsu = no

    # For updating the DQ headers
    int     hdw_num = 6
    char    hdwords_to_read[6], name_of_bpm, path_to_bpm, hdwords_comments[6]

    struct  l_struct, l_date_struct

    ###########################################################################
    # Get parameters
    l_inimages = inimages
    l_outpref = outpref
    l_outimages = outimages
    l_fl_over = fl_over
    l_fl_mult = fl_mult
    l_logfile = logfile
    l_verbose=verbose
    l_fl_bias = fl_bias
    l_bias = bias
    l_fl_flat = fl_flat
    l_key_filter = key_filter
    l_flat[1] = flat1
    l_flat[2] = flat2
    l_flat[3] = flat3
    l_flat[4] = flat4
    l_fl_vardq = fl_vardq
    l_sci_ext = sci_ext
    l_var_ext = var_ext
    l_dq_ext = dq_ext
    l_key_ron = key_ron
    l_key_gain = key_gain
    l_key_biassec = key_biassec
    l_nbiascontam = nbiascontam
    l_biasrows = biasrows
    l_fl_dark = fl_dark
    l_dark = dark
    l_fl_trim = fl_trim
    l_key_exptime = key_exptime
    l_sat = sat
    l_key_nodpix = key_nodpix
    l_key_nodcount = key_nodcount
    l_bpm = bpm
    l_fl_qecorr = fl_qecorr

    # GQECORR Parameters:
    l_qe_refim = qe_refim
    l_fl_keep_qeim = fl_keep_qeim
    l_qe_corrpref = qe_corrpref
    l_qecorr_data = qe_data
    l_qecorr_datadir = qe_datadir
    l_qe_corrimages = qe_corrimages

    # GPREPARE PARAMETERS:
    l_gp_outpref = gp_outpref
    l_key_mdf = key_mdf
    l_rawpath = rawpath
    l_fl_addmdf = fl_addmdf
    l_mdffile = mdffile
    l_mdfdir = mdfdir
    l_key_datasec = key_datasec
    l_gaindb = gaindb
    l_ron = ron
    l_gain = gain

    # COLBIAS PARAMETERS
    l_median = median
    l_fl_inter = fl_inter
    l_function = function
    l_order = order
    l_niterate = niterate
    l_low_reject = low_reject
    l_high_reject = high_reject

    # Keywords to add to headers of DQ extensions
    # Done here instead of the it's decleration - due to datasec and biassec
    # being able to be changed - MS
    hdwords_to_read[1] = "DETSEC"
    hdwords_to_read[2] = "CCDSEC"
    hdwords_to_read[3] = l_key_datasec
    hdwords_to_read[4] = "TRIMSEC"
    hdwords_to_read[5] = "CCDSUM"
    hdwords_to_read[6] = l_key_biassec

    # Their comments
    hdwords_comments[1] = "Detector section(s)"
    hdwords_comments[2] = "CCD section(s)"
    hdwords_comments[3] = "Data section"
    hdwords_comments[4] = "Trimmed section(s)"
    hdwords_comments[5] = "CCD sum"
    hdwords_comments[6] = "Bias section(s)"

    ###########################################################################
    # Initialize variables
    status = 0
    l_maxfiles = maxfiles
    badflag = no
    qecorrapp = "/" # How to apply the qecorrection
    globalcanqecorr = no # Needs to be set for clean up at end, regardless.
    canqecorr = no # Needs to be set for clean up at end, regardless.
    l_key_qeimfactor = "QEFACTOR"
    l_key_qecorrim = "QECORRIM"
    l_qeimfactor_comment = "Imaging QE correction factor"
    l_qecorrim_comment = "Image used to QE correct data"
    isimaging = no # Flag to say whether data is image data or not; changed if
                   #     needed after call to qecorr, where this is tested.
    l_qecorr_data = l_qecorr_datadir//l_qecorr_data

    for (ii = 1; ii <= l_maxfiles; ii += 1) {
        ll_corrimg[ii] = ""
        ll_corrimg_exists[ii] = no
        inqecorr_state[ii] = no
    }

    cache ("imgets", "keypar", "gimverify", "gemdate")

    # Number of unused rows at bottom of new Hamamatsu CCDs to remove
    row_removal_value = 48

    # For scaling Hamamatsu saturation values [default is for non-Hamamatsu
    # CCDs]
    satfrac = 1.0

    # Make temp files
    tmpfile = mktemp ("tmpfile")
    tmpinlist = mktemp ("tmpinlist")
    tmpoutlist = mktemp ("tmpoutlist")
    tmplog = mktemp ("tmplog")
    tmppl = mktemp ("tmppl")
    tmpqecorrimgs = mktemp("tmpqecorrimgs")

    ###########################################################################
    # Define the name of the logfile
    s_empty = ""
    print (l_logfile) | scan (s_empty)
    l_logfile = s_empty
    if (l_logfile == "") {
        l_logfile = gmos.logfile
        if (l_logfile == "") {
            l_logfile = "gmos.log"
            printlog ("WARNING - GIREDUCE: both gireduce.logfile and \
                gmos.logfile are empty.", l_logfile, l_verbose)
            printlog ("                Using default file gmos.log.",
                l_logfile,l_verbose)
        }
    }

    # Write to the logfile
    date | scan(l_date_struct)
    printlog ("---------------------------------------------------------------\
        -------------", l_logfile, l_verbose)
    printlog ("GIREDUCE -- "//l_date_struct, l_logfile, l_verbose)
    printlog ("", l_logfile, l_verbose)

    ###########################################################################
    # Check that the rawpath has a trailing slash and is a valid entry
    if (l_rawpath != "") {
        rphend = substr(l_rawpath,strlen(l_rawpath),strlen(l_rawpath))
        if (rphend == "$") {
            show (substr(l_rawpath,1,strlen(l_rawpath)-1)) | scan (pathtest)
            rphend = substr(pathtest,strlen(pathtest),strlen(pathtest))
        }
        if (rphend != "/") {
            l_rawpath = l_rawpath//"/"
        }
        if (!access(l_rawpath)) {
            printlog ("ERROR - GIREDUCE: Cannot access rawpath: "//l_rawpath, \
                l_logfile, verbose+)
            goto clean
        }
    }

    # Check the input images and place them in a file.
    nimages = 0
    atposition = stridx("@",l_inimages)
    if (atposition > 0) {

        testfile = substr(l_inimages,atposition+1,strlen(l_inimages))
        if (!access(testfile)) {
            printlog ("ERROR - GIREDUCE: The input list "//testfile//\
                " does not exist.", logfile=l_logfile, verbose+)
            goto clean
        }
    }

    files (l_inimages, sort-, > tmpfile)

    nbad = 0
    nfiles = 0
    prev_qecorr = 0
    extnum_mismatch = no
    scanfile = tmpfile
    while (fscan(scanfile, img) != EOF) {
        l_fl_raw = no
        fparse (img)
        img_root = fparse.root
        gimverify (img)
        img = gimverify.outname
        if ((gimverify.status >= 1) && (l_rawpath != "")) {
            gimverify (l_rawpath//img)
            l_fl_raw = yes
        }
        if (gimverify.status == 1) {
            printlog ("ERROR - GIREDUCE: File "//img//" not found.",
                l_logfile, verbose+)
            nbad += 1
        } else if (gimverify.status>1) {
            printlog ("ERROR - GIREDUCE: File "//img//" not a MEF FITS image.",
                l_logfile, verbose+)
            nbad += 1
        } else {

            # Check for BIAS / DARK subtraction if applying VARDQ or QECORR
            if ((l_fl_vardq || l_fl_qecorr) && !fl_bias) {
                # Check if BIAS subtracted - MS
                if (l_fl_raw)
                    keypar (l_rawpath//img//"[0]", "BIASIM", >& "dev$null")
                else
                    keypar (img//"[0]", "BIASIM", >& "dev$null")
                if (keypar.found == no && !fl_dark) {
                    # Check if DARK subtracted. Assumes bias removed
                    # when the DARK was subtracted (typical for N&S - MS
                    if (l_fl_raw)
                        keypar (l_rawpath//img//"[0]", "DARKIM", >& "dev$null")
                    else
                        keypar (img//"[0]", "DARKIM", >& "dev$null")

                    if (keypar.found == no) {
                        if (l_fl_vardq) {
                            # Cannot do vardq calculations
                            printlog("WARNING - GIREDUCE: Variance is "//\
                                "created using", l_logfile, verbose+)
                            printlog("                    photon statistics"//\
                                " requiring bias or dark subtraction.",
                                l_logfile, verbose+)
                            printlog("                    Image "//img//\
                                " has not been bias or dark subtracted.", \
                                l_logfile, verbose+)
                            printlog("                    cannot create "//\
                                "poisson variance. Setting fl_vardq=no", \
                                 l_logfile, verbose+)
                            # Reset vardq flag to no
                            l_fl_vardq=no
                        }

                        # For QE correction the same rules apply as for VAR/DQ
                        # - MS
                        if (l_fl_qecorr) {
                            printlog ("\nWARNING - GIREDUCE: QE correction "//\
                                "can only be applied if the image is to ", \
                                l_logfile, verbose+)
                            printlog ("                    be or has "//\
                                "previously been BIAS or DARK subtracted.",
                                l_logfile, verbose+)
                            printlog ("                    Reseting "//\
                                "fl_qecorr to no.", \
                                l_logfile, verbose+)
                            # Reset qe correction flag
                            l_fl_qecorr = no
                        }
                    }
                }

            }

            keyfound = ""
            if (l_fl_raw)
                hselect (l_rawpath//img//"[0]", "*PREPAR*", yes) | \
                    scan (keyfound)
            else
                hselect (img//"[0]", "*PREPAR*", yes) | scan (keyfound)

            if (l_fl_raw && (keyfound != "")) {
                printlog ("ERROR - GIREDUCE: Input image "//img//" is in raw \
                    data directory "//l_rawpath//" and is gprepared",
                    l_logfile, verbose+)
                goto crash
            }
            if (keyfound == "") {
                # Use already gprepared image if it exists
                if (access (l_gp_outpref//img_root//".fits")) {
                    hselect (l_gp_outpref//img_root//"[0]", "*PREPAR*", yes) |\
                        scan (keyfound)
                    if (keyfound == "") {
                        printlog ("ERROR - GIREDUCE: gprepare outpref \
                            parameter points to existing un-gprepared file",
                            l_logfile, verbose+)
                        goto crash
                    } else {
                        printlog ("Using already gprepared image "//\
                            l_gp_outpref//img_root, l_logfile, l_verbose)
                    }
                } else {
                    gprepare (img, rawpath=l_rawpath, outimages="",
                        outpref=l_gp_outpref, fl_addmdf=l_fl_addmdf,
                        key_mdf=l_key_mdf, gaindb=l_gaindb, mdffile=l_mdffile,
                        mdfdir=l_mdfdir, logfile=l_logfile, sci_ext=l_sci_ext,
                        key_ron=l_key_ron, key_gain=l_key_gain, ron=l_ron,
                        gain=l_gain, verbose=l_verbose)

                    if (gprepare.status != 0) {
                        printlog ("ERROR - GIREDUCE: gprepare returned with \
                            error", l_logfile, verbose+)
                        goto crash
                    }
                }
                img = l_gp_outpref//img_root
            }
            nfiles += 1
            input[nfiles] = img
            print (img, >> tmpinlist)
            # Do this here to read the number of sci extension - allows
            # checks against the input flats
            keypar (img//"[0]", "NSCIEXT", silent+)
            n_ext[nfiles] = int (keypar.value)
            if (nfiles > 1) {
                if (n_ext[nfiles] != n_ext[nfiles]) {
                    extnum_mismatch = yes
                }
            }
        }
    } # end while loop

    # Loop over input images and test for QE correction
    for (ii = 1; ii <= nfiles; ii += 1) {

        img = input[ii]

        # Perform OBSMODE checks for QE correction
        # Only allow Hamamatsu imaging data to be QE corrected
        if (l_fl_qecorr) {
            keypar (img//"[0]", "OBSMODE", silent+)
            if (keypar.value == "IMAGE") {
                keypar (img//"[0]", "DETTYPE", silent+)
                if (keypar.value != "S10892" && keypar.value != "S10892-N") {
                    printlog ("WARNING - GIREDUCE: Cannot QE "//\
                        "correct imaging data from EEV or e2vDD CCDs"//\
                        "                    Setting fl_qecorr=no", \
                        l_logfile, verbose+)
                    l_fl_qecorr = no
                }
            }
        }

        # Check if image has already been QE corrected
        # Test phu for a QE correction image (spectroscopy)
        keypar (img//"[0]", l_key_qecorrim, silent+)

        if (keypar.found) {
            # Record the QE state of the inout images (default is no)
            inqecorr_state[ii] = yes

            # Increment previously QE corrected counter
            prev_qecorr += 1

            # This allows a list of previously QE corrected images
            # to be shown
            printlog ("WARNING - GIREDUCE: "//img//" has "//\
                "already been QE corrected.", \
                l_logfile, verbose+)
        } else {
            # Test first SCI extension for a QE factor (imaging)
            keypar (img//"["//l_sci_ext//",1]", \
                l_key_qeimfactor, silent+)

            if (keypar.found) {
                # Record the QE state of the inout images (default is no)
                inqecorr_state[ii] = yes

                # Increment previously QE corrected counter
                prev_qecorr += 1

                # This allows a list of previously QE corrected images
                # to be shown
                printlog ("WARNING - GIREDUCE: "//img//" has "//\
                    "already been QE corrected.", \
                    l_logfile, verbose+)
            }
        }
    }

    # Switch off QE correction if any of the files have been QE corrected
    if (l_fl_qecorr && (prev_qecorr > 0)) {
        printlog ("WARNING - GIREDUCE: One or more input images have "//\
            "previously been QE corrected. \n"//\
            "                    Setting fl_qecorr=no for all images.", \
            l_logfile, verbose+)
        l_fl_qecorr = no
    }

    # Is there anything else to do?
    if (!l_fl_over && !l_fl_trim && !l_fl_bias && !l_fl_dark && !l_fl_flat \
        && !l_fl_qecorr) {
        printlog ("GIREDUCE: There is nothing left for GIREDUCE to do; "//\
            "fl_over, fl_trim, fl_bias, fl_dark", \
            l_logfile, verbose+)
        printlog ("          fl_qecorr and fl_flat are set to no. "//\
            "Exiting.", l_logfile, verbose+)
        goto clean
    }

    # Check for empty file list
    if (nfiles == 0) {
        printlog ("ERROR - GIREDUCE:  No input images to process.",
            l_logfile, verbose+)
        nimages = 0
        goto crash
    } else if (nfiles > l_maxfiles) {
        printlog ("ERROR - GIREDUCE:  Maximum number of input files is "//\
            l_maxfiles, l_logfile, verbose+)
        goto crash
    }

    # Exit if problems found with input files
    if (nbad > 0) {
        printlog ("ERROR - GIREDUCE: "//nbad//" image(s) either do not exist \
            or are not MEF files.", l_logfile, verbose+)
        goto crash
    } # end if (nbad > 0)

    printlog ("", l_logfile, l_verbose)
    printlog ("Input files:", l_logfile, l_verbose)
    if (l_verbose)
        type (tmpinlist)
    type (tmpinlist, >> l_logfile)
    printlog ("", l_logfile, l_verbose)
    delete (tmpfile, verify-, >& "dev$null")
    # Input images are in tmpinlist file and an array input[ii]
    nimages = nfiles

    # Test the VAR and DQ extension names, the RDNOISE, GAIN and DATASEC
    # keywords
    if (l_fl_vardq) {
        if ((l_var_ext == "") || (stridx(" ", l_var_ext) > 0) || \
            (l_dq_ext == "") || (stridx(" ",l_dq_ext) > 0) ) {

            printlog ("ERROR - GIREDUCE: Variance and/or data quality \
                extension names var_ext/dq_ext not defined", l_logfile,
                verbose+)
            goto crash
        }

        if ((l_key_ron == "") || (stridx(" ", l_key_ron) > 0)) {
            printlog ("ERROR - GIREDUCE: Readout noise keyword key_ron is not \
                defined", l_logfile, verbose+)
            goto crash
        }

        if ((l_key_gain == "") || (stridx(" ", l_key_gain) > 0)) {
            printlog ("ERROR - GIREDUCE: Gain keyword key_gain is not defined",
                l_logfile, verbose+)
            goto crash
        }

        if ((l_key_datasec == "") || (stridx(" ", l_key_datasec) > 0)) {
            printlog ("ERROR - GIREDUCE: Data section keyword key_datasec is \
                not defined", l_logfile, verbose+)
            goto crash
        }
    }

    ###########################################################################
    # Check the output images. Check that they have been defined, and that they
    # do not exist already. Count the number of output images. If they are not
    # the same as the input, exit with an error (only check if they are given
    # as outimage).

    s_empty = ""
    print (l_outimages) | scan (s_empty)
    l_outimages = s_empty
    ii = 0
    if (l_outimages != "") {
        atposition = stridx("@",l_outimages)
        if (atposition > 0) {
            testfile = substr(l_outimages,atposition+1,strlen(l_outimages))
            if (!access(testfile)) {
                printlog ("ERROR - GIREDUCE: The ouput list "//testfile//\
                    " does not exist.", logfile=l_logfile, verbose+)
                goto clean
            }
        }

        files (l_outimages, sort-, > tmpfile)
        scanfile = tmpfile
        while (fscan(scanfile,img) != EOF) {
            nfiles -= 1
            ii += 1
            gimverify (img)
            img = gimverify.outname
            if (gimverify.status != 1) {
                printlog ("WARNING - GIREDUCE: File "//img//" already exists.",
                    l_logfile, l_verbose)
                n_do[ii] = no
                nfiles += 1     # add 1 back to nfiles
            } else {
                output[ii] = img
                n_do[ii] = yes
                print (img, >> tmpoutlist)
            }
        } # end while loop
        scanfile = ""
        delete (tmpfile, verify-, >& "dev$null")
        ;
    } else {  # l_outimages==""
        scanfile = tmpinlist
        while (fscan(scanfile,img) != EOF) {
            nfiles -= 1
            ii += 1
            output[ii] = l_outpref//img
            gimverify (output[ii])
            output[ii] = gimverify.outname
            if (gimverify.status != 1) {
                printlog ("WARNING - GIREDUCE: File "//output[ii]//\
                    " already exists", l_logfile, l_verbose)
                n_do[ii] = no
                nfiles += 1     # add 1 back to nfiles
            } else {
                n_do[ii] = yes
                print (output[ii], >> tmpoutlist)
            }
        } # end while loop
        scanfile = ""
    }

    if (nfiles != 0) {
        printlog ("ERROR - GIREDUCE: problem with input and output images.", \
            l_logfile, verbose+)
        goto crash
    }

    printlog ("Output files:", l_logfile, l_verbose)
    if (l_verbose)
        type (tmpoutlist)
    type (tmpoutlist, >> l_logfile)
    printlog ("", l_logfile, l_verbose)

    # Output images are in array output[ii]
    # Don't delete the tmpoutlist here - used to pass files to qecorr - MS
    delete (tmpinlist, verify-, >& "dev$null")

    ###########################################################################
    # Check the extension names

    # SCI extension -
    s_empty = ""
    print (l_sci_ext) | scan (s_empty)
    l_sci_ext = s_empty
    if (l_sci_ext == "") {
        printlog ("ERROR - GIREDUCE: Science extension parameter sci_ext is \
            missing", l_logfile, verbose+)
        goto crash
    }

    ###########################################################################
    # If flat-fielding: Check that filter key is defined. Check that flat
    # field images are ok. Get the number of science extensions.
    # The obsmode for the flat fields is the same as that of the input images.
    for (ii=1; ii<=4; ii+=1) {
        canflat[ii] = no
        flatqecorr_state[ii] = no
    }
    if (l_fl_flat) {

        # The header keyword
        s_empty = ""
        print (l_key_filter) | scan (s_empty)
        l_key_filter = s_empty
        if (l_key_filter == "") {
            printlog ("WARNING - GIREDUCE: Cannot flat field correct.",
                l_logfile, l_verbose)
            printlog ("                    Header keyword parameter \
                key_filter is not defined", _logfile, l_verbose)
            l_fl_flat = no
            goto noflats
        }

        # The flat field images -
        for (ii=1; ii<=4; ii+=1) {
            s_empty = ""
            print (l_flat[ii]) | scan (s_empty)
            l_flat[ii] = s_empty
            if (l_flat[ii]!="") {
                gimverify (l_flat[ii])
                l_flat[ii] = gimverify.outname

                # Used to update the output name
                fparse(l_flat[ii])
                l_fname[ii] = fparse.root

                if (gimverify.status == 0) {

                    keypar (l_flat[ii]//"[0]", "NSCIEXT", silent+)
                    n_extflat[ii] = int(keypar.value)

                    canflat[ii] = yes
                    keypar (l_flat[ii]//"[0]", "GIFLAT", silent+)
                    if (!keypar.found) {

                        # Check fo spectroscopic flats that have the same
                        # number of SCI extensions - This is done instead of
                        # updating the flat image on disk in gsreduce - MS
                        keypar (l_flat[ii]//"[0]", "GSFLAT", silent+)
                        if (keypar.found) {
                            if (extnum_mismatch) {
                                ncheck = nimages
                            } else {
                                ncheck = 1
                            }
                            checkcount = 0
                            for (kk = 1; kk <= ncheck; kk += 1) {
                                if (n_extflat[ii] != n_ext[kk]) {
                                    checkcount += 1
                                }
                            }
                            if (checkcount == ncheck) {
                                printlog ("ERROR - GSFLAT: Flat field: "//\
                                    l_flat[ii]//" was created with "//\
                                    "gsflat but the number of "//\
                                    l_sci_ext//" extensions does not "//\
                                    "match the number "//\
                                    "of "//l_sci_ext//" in any of the "//\
                                    " input images.", \
                                    l_logfile, verbose+)
                                    badflag = yes
                            }
                        } else {
                            printlog ("ERROR - GIREDUCE: Flat field: "//\
                                l_flat[ii]//" was not created with giflat \
                                task.",
                                l_logfile, verbose+)
                            badflag = yes
                        }
                    }

                    if (l_fl_vardq) {
                        if ((!imaccess(l_flat[ii]//"["//l_var_ext//"]")) \
                            || (!imaccess(l_flat[ii]//"["//l_dq_ext//"]"))) {
                            printlog ("WARNING - GIREDUCE: Flat field image \
                                does not contain both VAR and DQ planes.", \
                                l_logfile, verbose+)
                            printlog ("                    Setting \
                                fl_vardq=no and proceeding.", l_logfile, \
                                verbose+)
                            l_fl_vardq = no
                        }
                    }

                    # Read the QE correction state of the FLAT image
                    # Default is no
                    # Spectroscopy data
                    keypar (l_flat[ii]//"[0]", l_key_qecorrim, silent+)
                    if (keypar.found) {
                        flatqecorr_state[ii] = yes

                    } else {
                        # Imaging data
                        keypar (l_flat[ii]//"["//l_sci_ext//",1]", \
                            l_key_qeimfactor, silent+)

                        if (keypar.found) {
                            flatqecorr_state[ii] = yes
                        }
                    }

                    # Do a quick check here for QE corrected FLAT but none
                    # QE corrected or not to be asked to QE corrected
                    # More checks are done in the canflat loops
                    if (flatqecorr_state[ii] && (prev_qecorr == 0)) {

                        if (!l_fl_qecorr) {
                            printlog ("WARNING - GIREDUCE: FLAT image: "//\
                                l_flat[ii]//" has been QE corrected but\n"//\
                                "                    the input data have"//\
                                " not been QE corrected and QE correction\n"//\
                                "                    of the inpout data has"//\
                                " not been requested.\n"//\
                                "                    Switching fl_flat to "//\
                                "no.",\
                                l_logfile, verbose+)
                            l_fl_flat = no
                        }
                    } else if (!flatqecorr_state[ii] && \
                        (prev_qecorr == nimages)) {
                        # Flat not QE corrected but all of the inputs are

                        printlog ("ERROR - GIREDUCE: FLAT image: "//\
                            l_flat[ii]//" has not been QE corrected but\n"//\
                            "                    the input data have all"//\
                            "been QE corrected.\n"//\
                            "                    Switching fl_flat to no.", \
                            l_logfile, verbose+)
                        l_fl_flat = no
                    }

                    imgets (l_flat[ii]//"[0]", "OBSMODE", >& "dev$null")
                    obsmodeflat[ii] = imgets.value

                    # Get the name of the filter in array flatfilter[4]
                    keypar (l_flat[ii]//"[0]", l_key_filter, silent+)
                    flatfilter[ii] = keypar.value
                    if (!keypar.found)
                        canflat[ii] = no

                    # Place the sizes of the flats in the tsec_flat[60,4] and
                    # fl_tsec_flat[60,4] arrays.
                    for (jj=1; jj<=n_extflat[ii]; jj+=1) {
                        keypar (l_flat[ii]//"["//l_sci_ext//","//jj//"]",
                            "TRIMSEC", silent+)
                        tsec_flat[jj,ii] = keypar.value
                        fl_tsec_flat[jj,ii] = keypar.found
                        imgets (l_flat[ii]//"["//l_sci_ext//","//jj//"]",
                            "DETSEC", >& "dev$null")
                        ccd_flat[ii,jj] = imgets.value
                        imgets (l_flat[ii]//"["//l_sci_ext//","//jj//"]",
                            "CCDSUM", >& "dev$null")
                        ccd_flat[ii,jj] = ccd_flat[ii,jj]//imgets.value
                     }
                } else {
                    printlog ("ERROR - GIREDUCE: Can not find flat-field \
                        frame: "//l_flat[ii], l_logfile, verbose+)
                    badflag = yes
                }
            }
        } # end of checking flat field images

noflats:

    } # end of   if(l_fl_flat)

    ###########################################################################
    # If overscan subt.: Check that bias section keywords is defined and is ok.
    if (l_fl_over) {
        s_empty = ""
        print (l_key_biassec) | scan(s_empty)
        l_key_biassec = s_empty
        if (l_key_biassec == "") {
            printlog ("ERROR - GIREDUCE: bias section keyword (key_biassec) \
                is an empty string.", l_logfile, verbose+)
            goto crash
        }
    }

    ###########################################################################
    # If bias subtracting: check that the bias image is defined and is ok.
    # If fl_vardq=yes, need to check that the bias image has var and dq planes.
    # Get the number of science extensions.
    canbias = no
    if (l_fl_bias) {
        s_empty = ""
        print (l_bias) | scan (s_empty)
        l_bias = s_empty
        gimverify (l_bias)
        l_bias = gimverify.outname

        # Used in updating the output header
        fparse (l_bias)
        l_bname = fparse.root

        if (gimverify.status == 0) {
            canbias = yes
            keypar (l_bias//"[0]", "GBIAS", silent+)
            if (!keypar.found) {
                printlog ("ERROR - GIREDUCE: Bias frame: "//l_bias//" was \
                    not created with gbias task.", l_logfile, verbose+)
                badflag = yes
            }

            imgets (l_bias//"[0]", "NSCIEXT", >& "dev$null")
            n_extbias = int(imgets.value)
            for (ii=1; ii<=n_extbias; ii+=1) {
                keypar (l_bias//"["//l_sci_ext//","//ii//"]", "TRIMSEC",
                    silent+)
                tsec_bias[ii] = keypar.value
                fl_tsec_bias[ii] = keypar.found
                imgets (l_bias//"["//l_sci_ext//","//ii//"]", "DETSEC",
                    >& "dev$null")
                ccd_bias[ii] = imgets.value
                imgets (l_bias//"["//l_sci_ext//","//ii//"]", "ccdsum",
                    >& "dev$null")
                ccd_bias[ii] = ccd_bias[ii]//imgets.value
            }

        } else {
            printlog ("ERROR - GIREDUCE: Can not find bias frame: "//l_bias,
                l_logfile, verbose+)
            badflag = yes
        }
    }

    ###########################################################################
    # If dark current correcting: check that the dark image is defined and is
    # ok. If fl_vardq=yes, need to check that the dark image has var and dq
    # planes. Get the number of science extensions. Check that the exposure
    # time keyword is defined
    candark = no
    n_extdark = 0
    if (l_fl_dark) {
        s_empty = ""
        print (l_dark) | scan (s_empty)
        l_dark = s_empty
        gimverify (l_dark)
        l_dark = gimverify.outname

        # Used to update output image header
        fparse(l_dark)
        l_dname = fparse.root

        if (gimverify.status == 0) {
            candark = yes
            keypar (l_dark//"[0]", "GNSDARK", silent+)
            if (!keypar.found) {
                printlog ("ERROR - GIREDUCE: Dark frame: "//l_dark//\
                    " was not created with gndark task.", l_logfile, verbose+)
                badflag = yes
            }

            imgets (l_dark//"[0]", "NSCIEXT", >>& "dev$null")
            n_extdark = int(imgets.value)
            for (ii=1; ii<=n_extdark; ii+=1) {
                keypar (l_dark//"["//l_sci_ext//","//ii//"]", "TRIMSEC",
                    silent+)
                tsec_dark[ii] = keypar.value
                fl_tsec_dark[ii] = keypar.found
                imgets (l_dark//"["//l_sci_ext//","//ii//"]", "DETSEC",
                    >& "dev$null")
                ccd_dark[ii] = imgets.value
                imgets (l_dark//"["//l_sci_ext//","//ii//"]", "ccdsum",
                    >& "dev$null")
                ccd_dark[ii] = ccd_dark[ii]//imgets.value
            }

            imgets (l_dark//"[0]", l_key_exptime, >>& "dev$null")
            if (imgets.value == "0") {
                printlog ("ERROR - GIREDUCE: Dark image "//l_dark//" contains \
                    no exposure time", l_logfile, verbose+)
                badflag = yes
            } else {
                l_darkexptime = real(imgets.value)
            }

            imgets (l_dark//"[0]", l_key_nodpix, >>& "dev$null")
            if (imgets.value == "0") {
                printlog ("ERROR - GIREDUCE: Dark image "//l_dark//" contains \
                    no shuffle distance", l_logfile, verbose+)
                badflag = yes
            } else {
                l_darknodpix = real(imgets.value)
            }

            imgets (l_dark//"[0]", l_key_nodcount, >>& "dev$null")
            if (imgets.value == "0") {
                printlog ("ERROR - GIREDUCE: Dark image "//l_dark//" contains \
                    no number of cycles", l_logfile, verbose+)
                badflag = yes
            } else {
                l_darknodcount = real(imgets.value)
            }

        } else {
            printlog ("ERROR - GIREDUCE: Can not find dark frame: "//l_dark,
                l_logfile, verbose+)
            badflag = yes
        }
    }

    # Now check for existence of Bad Pixel Mask
    if (l_fl_vardq && (imaccess(l_bpm) == no) && (l_bpm != "") && \
        (stridx(" ", l_bpm) <= 0)) {

        printlog ("WARNING - GIREDUCE: Bad Pixel Mask "//l_bpm//" not found",
            l_logfile, verbose+)
        printlog ("                    Only saturated pixels will be flagged",
            l_logfile, verbose+)
        l_bpm = ""
        next_bpm = 0
    } else if (l_fl_vardq && ((l_bpm == "") || (stridx(" ",l_bpm) > 0))) {
        printlog ("WARNING - GIREDUCE: Bad Pixel Mask filename is an empty \
            string", l_logfile, verbose+)
        printlog ("                    Only saturated pixels will be flagged",
            l_logfile, verbose+)
        l_bpm = ""
        next_bpm = 0
    }

    # This is to store the BPM name without any paths included when updating
    # headers later on
    fparse (l_bpm)
    name_of_bpm = str(fparse.root)//str(fparse.extension)

    # Check if BPM is a MEF and keep the number of extensions for later use
    if (l_fl_vardq && (l_bpm != "")) {
        gimverify (l_bpm)
        if (gimverify.status == 0) {
            # Must be .fits (gimverify strips the .fits if it was there, so we
            # have to put it back)
            l_bpm = gimverify.outname//".fits"
            printlog ("GIREDUCE: bmp = "//l_bpm, l_logfile, l_verbose)
        } else if (gimverify.status == 1) {
            printlog ("ERROR - GIREDUCE: BPM "//l_bpm//" does not exist",
                l_logfile, verbose+)
            goto crash
        } else {
            printlog ("ERROR - GIREDUCE: BPM "//l_bpm//" is not a MEF file.",
                l_logfile, verbose+)
            goto crash
        }

        # Now get the number of science extensions. We can change this back to
        # imget(NEXTEND) as we are keeping the keyword in the PHU (it is
        # updated at the end of the task).
        # Assumes all extensions are to be used as BPM, IJ

        imgets (l_bpm//"[0]", "NEXTEND", >& "dev$null")
        next_bpm = int(imgets.value)

        # Get the OBSMODE of the BPM
        #imgets (l_bpm//"[0]", "OBSMODE", >& "dev$null")
        #obs_bpm = imgets.value
    }

    ###########################################################################
    # If any of the calibration files is missing, exit the program.
    if (badflag)
        goto crash

    ###########################################################################
    # Redefine the flags:
    if (!canbias)
        l_fl_bias = no
    if (!candark)
        l_fl_dark = no
    if (!canflat[1] && !canflat[2] && !canflat[3] && !canflat[4])
        l_fl_flat = no

    ###########################################################################
    # Do the actual job:
    # The input images name are in an array input[ii]
    # The output images name are in an array: output[ii]
    # If there is a MDF file, place it in the output MEF file.
    # Do the overscan of biassec with colbias.
    # Trim the image to datasec
    # Do the bias subtraction
    # Do the dark subtraction
    # Do the flat field correction
    # Do the job with imexpr and the expressions above.

    for (ii=1; ii<=nimages; ii+=1) {
        # Make the output file - place header, only if the image does not exist
        if (n_do[ii]) {
            imcopy (input[ii]//"[0]", output[ii]//".fits", verbose-)
        } else {
            input[ii] = output[ii]
        }

        # Set default satvalue and nbiascontamvalue depending on detector type
        imgets (input[ii]//"[0]", "DETTYPE", >& "dev$null")
        if (imgets.value == "SDSU II CCD") { # Current CCDs 2011-07-25
            osc_order = 1
            nbiascontamvalue = 4 ##M NBIASCONTAM
        } else if (imgets.value == "SDSU II e2v DD CCD42-90") {
            # New e2vDD CCDs
            osc_order = 1
            nbiascontamvalue = 5 ##M NBIASCONTAM
        } else if (imgets.value == "S10892") { #Hamamatsu CCDs SOUTH
            # Set Hamamatsu flag
            osc_order = 7
            ishamamatsu = yes
            nbiascontamvalue = 4 ##M NBIASCONTAM
            lower_bias_y = row_removal_value
            satfrac = 0.98
        } else if (imgets.value == "S10892-N") { # Hamamatsu NORTH
            # Set Hamamatsu flag
            osc_order = 7
            ishamamatsu = yes
            nbiascontamvalue = 4 ##M NBIASCONTAM
            lower_bias_y = row_removal_value
            satfrac = 0.98
        }

        # Adjust bias contamination if user supplies one
        if (l_nbiascontam != "default") {
            # Option for user-defined value
            nbiascontamvalue = int(l_nbiascontam)
        }

        if (l_order != "default") {
            # OPtion for user-defined value
            osc_order = int(l_order)
        }

        # The number of extensions will be needed
        imgets (input[ii]//"[0]", "NEXTEND", >& "dev$null")
        n_ext[ii] = int(imgets.value)

        # Add the MDF to the output file if appropriate
        hasmdf[ii] = no
        imgets (input[ii]//"[0]", "OBSMODE", >>& "dev$null")
        if (imgets.value != "IMAGE") {
            for (jj = 1; jj<=n_ext[ii]; jj+=1) {
                keypar (input[ii]//".fits["//jj//"]", "EXTNAME", silent=yes,
                    >& "dev$null")
                if (keypar.value == "MDF") {
                    hasmdf[ii] = yes
                }
            }
        }
        # Need number of SCI extensions, not number of total extensions
        # for loops
        hselect (input[ii]//"[0]", "NSCIEXT", yes) | scan (n_ext[ii])

        # Moved the loop to populate the ccs[ii,jj] array with DETSEC and
        # CCDSUM numbers into the if (l_fl_bias || l_fl_flat || l_fl_dark ||
        # l_fl_over) loop. - MS
    }

    ###########################################################################
    # Populate the overvar[60,100] matrix, rms of overscan region
    for (ii=1; ii<=nimages; ii+=1) {
        for (jj=1; jj<=n_ext[ii]; jj+=1)
            overvar[jj,ii] = 0.
    }

    ###########################################################################
    if (l_fl_over) {
        for (ii=1; ii<=nimages; ii+=1) {
            for (jj=1; jj<=n_ext[ii]; jj+=1) {

                # Test if already overscan subtracted
                keypar (input[ii]//"["//l_sci_ext//","//jj//"]", \
                    "OVERSCAN", silent+)
                if (keypar.found) {
                    # Need to copy image if any of the following flags
                    # are true - MS
                    if (l_fl_trim || l_fl_bias || l_fl_flat || l_fl_dark || \
                        l_fl_qecorr) {

                        printlog ("WARNING - GIREDUCE: "//input[ii]//\
                            " already OVERSCAN corrected", l_logfile, \
                            l_verbose)
                        for (jj=1; jj<=n_ext[ii]; jj+=1) {
                            imcopy (input[ii]//"["//l_sci_ext//","//\
                                jj//"]",output[ii]//"["//l_sci_ext//","//jj//\
                                ",append]", verbose-)
                        }
                        # Skip to next image as already overscan subtracted -MS
                        goto next_overscan_image
                    } else {
                        # Nothing else to do go crash - MS
                        printlog ("ERROR - GIREDUCE: "//input[ii]//\
                            " already OVERSCAN corrected and no other"//\
                            " reduction steps requested. Exiting.", l_logfile,\
                            verbose+)
                        goto crash
                    }
                }

                # Access the BIASSEC keyword value
                keypar (input[ii]//"["//l_sci_ext//","//jj//"]",
                    l_key_biassec, silent+)
                l_biassec = keypar.value
                if (!keypar.found) {
                    printlog ("ERROR - GIREDUCE: "//input[ii]//"[\
                        "//l_sci_ext//","//jj//"] does not have biassec", \
                        l_logfile, l_verbose)
                    goto crash
                }

                # Use a subsection of the BIASSEC for calulating the overscan
                # RMS, if nbiascontamvalue > 0.
                if (nbiascontamvalue > 0) {
                    # Parse the BIASSEC value to local variables
                    print (l_biassec) | translit ("STDIN", ":, [, ]", " ") \
                        | fields ("STDIN", 1-5) | scan (biasx1, biasx2, \
                        biasy1, biasy2)

                    # Previously the test of the BIASSEC was testing hardcoded
                    # values and missed certain binning values for the
                    # different amplifier modes. With the introduction of the
                    # new Hamamatsu CCDs with more amplifier configurations
                    # this hardcoded test was changed to test the NAXIS1
                    # header keyword value against the righthand edge of the
                    # the BIASSEC instead. - MS

                    # Access the NAXIS1 header keyword value
                    keypar (input[ii]//"["//l_sci_ext//","//jj//"]",
                        "i_naxis1", silent+)
                    if (keypar.found) {
                        naxis1 = int(keypar.value)
                    } else {
                        # File should have NAXIS1 if it's a FITS file
                        printlog ("ERROR - GIREDUCE: "//input[ii]//"[\
                            "//l_sci_ext//","//jj//"] does not have NAXIS1 \
                            keyword.", l_logfile, l_verbose)
                        goto crash
                    }

                    if (biasx2 == naxis1) {
                        # The overscan section is on the right side
                        biasx1 = biasx1 + nbiascontamvalue
                        biasx2 = biasx2 - 1 # Take out edge pixel
                    } else {
                        # The overscan section is on the left side
                        biasx1 = biasx1 + 1 # Take out edge pixel
                        biasx2 = biasx2 - nbiascontamvalue
                    }

                    # Added by JT: also restrict the rows used
                    if (l_biasrows!="default") {
                        print (l_biasrows) | scanf ("%d:%d", biasy1, biasy2)
                    } else if (ishamamatsu) {

                        keypar (input[ii]//"["//l_sci_ext//","//jj//"]",
                            "CCDSUM", silent+)
                        if (keypar.found) {
                            ybin = int(keypar.value)
                        } else {
                            # File should have CCDSUM if it's a FITS file
                            printlog ("ERROR - GIREDUCE: "//input[ii]//"[\
                                "//l_sci_ext//","//jj//"] does not have CCDSUM \
                                keyword.", l_logfile, l_verbose)
                            goto crash
                        }
                        biasy1 = int(lower_bias_y / ybin) + 1
                    }
                    l_biassec = "["//biasx1//":"//biasx2//","//biasy1//":"//\
                        biasy2//"]"
                }

                # call colbias to identify and subtract overscan value
                print ("yes") | \
                    colbias (input[ii]//"["//l_sci_ext//","//jj//"]", \
                    output[ii]//"["//l_sci_ext//","//jj//",append]", \
                    bias=l_biassec, trim="[]", median=l_median, \
                    interactive=l_fl_inter, function=l_function, \
                    order=osc_order, low_reject=l_low_reject, \
                    high_reject=l_high_reject, niterate=l_niterate, \
                    logfile=tmplog, graphics="stdgraph", cursor="", \
                    >& "dev$null")

                # Update the current extension of the output file with
                # the new keyword OVERSEC which details the region of
                # BIASSEC used to calculate the overscan RMS value
                # - MS
                gemhedit (output[ii]//"["//l_sci_ext//","//jj//"]", \
                    "OVERSEC", l_biassec, \
                    "Section used for overscan calculation", delete-)

                match ("RMS", tmplog, stop-) | scan (junk, junk, junk, srms)
                overvar[jj,ii] = real(srms)
                overvar[jj,ii] = overvar[jj,ii]**2  # overscan sub variance
                gemhedit (output[ii]//"["//l_sci_ext//","//jj//"]",
                    "OVERRMS", srms, "Overscan RMS value from colbias",
                    delete-)

                # colbias log does not provide mean value of overscan
                # calling imstat just to put overscan value in header
                # OVERSCAN key is used by giflat, gdark, and gbias
                ##M This should be updated to use the curve used to subtract
                ##M the overscan itself and possibly update to use median
                ##M (midpt)
                imstat (input[ii]//"["//l_sci_ext//","//jj//"]"//l_biassec,
                    fields="mean", lower=INDEF, upper=INDEF, nclip=0,
                    lsigma=INDEF, usigma=INDEF, binwidth=0.1, format=no,
                    cache=no) | scan (dummy)
                gemhedit (output[ii]//"["//l_sci_ext//","//jj//"]", \
                    "OVERSCAN", dummy, "Overscan mean value", delete-)
                delete (tmplog, verify-, >& "dev$null")
            } # End of loop over extensions jj in overscan subtract section

            printlog ("GIREDUCE: Image "//output[ii]//" overscan subtracted", \
                l_logfile, l_verbose)

# For use when the image has already been overscan subtracted
next_overscan_image:

            flpr #to avoid a nasty seg fault when calling from shell
        } # End of loop over ii images in overscan subtraction section
    } # End of l_fl_over section

    ###########################################################################
    if (l_fl_trim) {
        for (ii=1; ii<=nimages; ii+=1) {
            # Test if already trimmed - MS
            keypar (input[ii]//"[0]", "TRIMMED", silent+)
            if (keypar.found) {

                # Need to copy image if any of the following flags
                # are true - MS
                if (l_fl_bias || l_fl_flat || l_fl_dark || l_fl_qecorr) {
                    printlog ("WARNING - GIREDUCE: "//input[ii]//\
                        " already TRIMMED", l_logfile, l_verbose)
                    # Only need to copy if the l_fl_over flag was false - MS
                    if (!l_fl_over) {
                        for (jj=1; jj<=n_ext[ii]; jj+=1) {
                            imcopy (input[ii]//"["//l_sci_ext//","//jj//"]",
                                output[ii]//"["//l_sci_ext//","//jj//\
                                ",append]", verbose-)
                        }
                    }
                    # Skip to next image as already trimmed subtracted -MS
                    goto next_trim_image
                } else {
                    # Nothing else to do go crash - MS
                    printlog ("ERROR - GIREDUCE: "//input[ii]//\
                        " already trimmed corrected and no futher"//\
                        " reduction steps requested. Exiting.", l_logfile,\
                        verbose+)
                    goto crash
                }
            }

            for (jj=1; jj<=n_ext[ii]; jj+=1) {
                # Inititate data section variable
                l_datasec = ""
                # Access the datasec keyword value for this extension
                keypar (input[ii]//"["//l_sci_ext//","//jj//"]",
                    l_key_datasec, silent+)
                if (keypar.found) {
                    l_datasec = keypar.value
                } else {
                    # Image should have DATASEC - if not crash
                    printlog ("ERROR - GIREDUCE:  "//l_key_datasec//" not \
                        found.", l_logfile, verbose+)
                    goto crash
                }

                if (ishamamatsu) {
                    # For Hamamatsu data: a test is required to determine if
                    # the bottom 48 rows (unbinned value) are to be removed
                    # (full frame) or not (e.g. stamp images) - MS

                    # DETSEC and CCDSUM are required to determine whether
                    # bottom 48 rows (unbinned) need removing -MS

                    # Obtain CCDSUM (binning) value for comparison
                    # calculation - MS
                    keypar (input[ii]//"["//l_sci_ext//","//jj//"]",
                       "CCDSUM", silent+)
                    if (keypar.found){
                        print (keypar.value) | scanf ("%d %d", xbin, ybin)
                    } else {
                        # It should exist - if not crash.
                        printlog ("ERROR - GIREDUCE: CCDSUM not found.", \
                            l_logfile, verbose+)
                        goto crash
                    }

                    # Initialise l_(ccd/det)sec(_new) local variables for
                    # updating headers
                    l_ccdsec = ""
                    l_ccdsec_new = ""
                    l_detsec = ""
                    l_detsec_new = ""

                    # Obtain DETSEC to determine if bottom 48 rows require
                    # removal - MS
                    keypar (input[ii]//"["//l_sci_ext//","//jj//"]",
                        "DETSEC", silent+)
                    if (keypar.found) {
                        l_detsec = keypar.value
                        print (keypar.value) | translit ("STDIN",\
                            ":, [, ]", " ") | fields ("STDIN", 1-5) | \
                            scan (detsecx1, detsecx2, detsecy1, \
                            detsecy2)
                    }

                    # Calculate the binned equivalent of the 48 rows that
                    # should be removed if necessary - MS
                    row_removal_value_binned = nint(row_removal_value / ybin)

                    # Remove lower 48 (unbinned) rows if detsecy1 is
                    # less than or equal to the binned equivalent of
                    # the bottom 48 rows - MS
                    if (detsecy1 <= row_removal_value) {

                        print (l_datasec) | translit ("STDIN", ":, [, ]",\
                             " ") | fields ("STDIN", 1-5) | scan (datax1,\
                              datax2, datay1, datay2)

                        datay1 = datay1 + ((row_removal_value - \
                            (detsecy1 - 1)) / ybin)
                        l_datasec = "["//datax1//":"//datax2//","//\
                            datay1//":"//datay2//"]"

                        # Obtain CCDSEC and update to l_ccdsec_new, so it shows
                        # the region of the CCD that the new l_datasec_new
                        # represents - MS
                        keypar (input[ii]//"["//l_sci_ext//","//jj//"]",
                            "CCDSEC", silent+)
                        if (keypar.found) {
                            l_ccdsec = keypar.value
                            print (keypar.value) | translit ("STDIN",\
                                 ":, [, ]", " ") | fields ("STDIN", 1-5) \
                                 | scan (ccdsecx1, ccdsecx2, ccdsecy1,\
                                  ccdsecy2)

                            ccdsecy1 = ccdsecy1 + row_removal_value
                            l_ccdsec_new = "["//ccdsecx1//":"//ccdsecx2//\
                                ","//ccdsecy1//":"//ccdsecy2//"]"
                        }

                        # Update DETSEC to l_detsec_new, so it shows
                        # the region of the CCD that the new l_datasec_new
                        # represents - MS
                        if (l_detsec != "") {
                            detsecy1 = detsecy1 + row_removal_value
                            l_detsec_new = "["//detsecx1//":"//detsecx2//\
                                ","//detsecy1//":"//detsecy2//"]"
                        }
                    } else {
                        l_detsec_new = l_detsec
                        l_ccdsec_new = l_ccdsec
                    }
                }

                ##M What about CRPIX2 for all CCD types including Hamamatsu?
                ##M Stop the overwriting...

                # Do the trimming
                # If the image has been overscan subtracted in this pass
                # the output image science extensions have already been created
                # and need overwriting, else they will need appending to the
                # output image - MS
                if (l_fl_over) {
                    # Need to copy to a temporary file and sleep briefly
                    # before overwriting to prevent segmentation violation
                    # errors from occuring in cl - EH
                    tmptrim = mktemp ("tmptrim")
                    imcopy (output[ii]//"["//l_sci_ext//","//jj//"]"//\
                        l_datasec, tmptrim, verbose-)
                    sleep(1)
                    imcopy (tmptrim, \
                        output[ii]//"["//l_sci_ext//","//jj//",overwrite]", \
                        verbose-)
                    imdelete (tmptrim, verify-, >& "dev$null")
                } else {
                    imcopy (input[ii]//"["//l_sci_ext//","//jj//"]"//\
                        l_datasec, output[ii]//"["//l_sci_ext//","//jj//\
                        ",append]", verbose-)
                }

                # Update the current extension header with the keyword TRIMSEC
                # and the area (l_datasec) trimmed.
                gemhedit (output[ii]//"["//l_sci_ext//","//jj//"]", \
                    "TRIMSEC", l_datasec, "Trimmed section(s)", delete-)

                if (ishamamatsu) {
                    if (l_ccdsec_new != l_ccdsec) {
                        # Update CCDSEC keyword
                        gemhedit (output[ii]//"["//l_sci_ext//","//jj//"]", \
                            "CCDSEC", l_ccdsec_new, "", delete-)
                    }
                    if (l_detsec_new != l_detsec) {
                        # Update DETSEC keyword
                        gemhedit (output[ii]//"["//l_sci_ext//","//jj//"]", \
                            "DETSEC", l_detsec_new, "", delete-)
                    }
                }
            } # End of loop over extensions jj in trimming section

            # Update datasec for all extensions
            gsetsec (output[ii], key_datsec=l_key_datasec)
            # Add TRIMMED keyword to phu of output image
            gemhedit (output[ii]//"[0]", "TRIMMED", "yes",
                "Overscan section trimmed", delete-)
            printlog ("GIREDUCE: Image "//output[ii]//" trimmed",
                l_logfile, l_verbose)

# For use when the image has already been trimmed
next_trim_image:

        } # End of loop over ii images in trimming section

    } else if (!l_fl_over) { # If no overscan subtraction or trimming
        for (ii=1; ii<=nimages; ii+=1) {
            for (jj=1; jj<=n_ext[ii]; jj+=1) {
                imcopy (input[ii]//"["//l_sci_ext//","//jj//"]",
                    output[ii]//"["//l_sci_ext//","//jj//",append]", verbose-)
            }
        }
    } # End of l_fl_trim loop

    ###########################################################################
    # There are 6 cases:
    # 1) bias subtraction only
    # 2) dark correction only
    # 3) flat field correction only
    # 4) bias + dark correction
    # 5) bias + flat field
    # 6) dark + flat field

    # Merge all 6 cases
    if (l_fl_bias || l_fl_flat || l_fl_dark || l_fl_over || l_fl_qecorr) {

        # QE Correction
        # Call GQECORR - it does all of the required checking and logic
        # In most cases the user will only supply one reference image/data
        # with the same filter - so, for most cases it is quicker to call
        # GQECORR with all of the images being passed to it at once.
        # If spectroscopy - QECORR will write to disk the same number of
        # correction images as there are reference images.
        # For imaging - if all inputs have the same filters, then only one file
        # will be written to disk, else there will be one file for each input.
        # Files to pass to qecorr are in tmpoutlist
        if (l_fl_qecorr) {

            printlog ("\nGIREDUCE: calling GQECORR", l_logfile, l_verbose)

            # Use tmpoutlist to call GQECORR
            gqecorr (inimages="@"//tmpoutlist, outimages="", outpref="qe",\
                refimages=l_qe_refim, fl_correct=no, \
                corrimages=l_qe_corrimages,\
                corrimpref=l_qe_corrpref, fl_vardq=l_fl_vardq, \
                sci_ext=l_sci_ext, var_ext=l_var_ext, dq_ext=l_dq_ext, \
                mdf_ext="MDF", key_detsec="DETSEC", key_ccdsec="CCDSEC", \
                key_datasec=l_key_datasec, key_biassec=l_key_biassec, \
                key_ccdsum="CCDSUM", qecorr_data=l_qecorr_data, \
                logfile=l_logfile, verbose=l_verbose, \
                scanfile2=tmpqecorrimgs)

            # Check exit status of GQECORR and set flags accordingly
            # and that the tmp file containing the list of correction files
            # and whether they existed previous to the call to GQECORR, exist
            if (gqecorr.status == 0 && access(tmpqecorrimgs)) {
                printlog ("GIREDUCE: Call to GQECORR successful", l_logfile,
                    l_verbose)

                scanfile = ""
                filename = ""
                num_qecorrimgs = 1
                scanfile = tmpqecorrimgs

                # Read in the correction data file names and the boolean
                # saying whether the correction data existed on disk before
                # the call to gqecorr
                while (fscanf(scanfile,"%s %s",filename,existed) != EOF) {

                    ll_corrimg[num_qecorrimgs] = filename
                    if (existed == "yes") {
                        ll_corrimg_exists[num_qecorrimgs] = yes
                    }
                    num_qecorrimgs += 1
                }

                # Adjust count accordingly
                num_qecorrimgs -= 1

                # Set can QE correction flags
                canqecorr = yes
                globalcanqecorr = yes

                # Check the obsmode keyword of first output image
                # Will determine the naming convention for the correction data
                # Can be done due to checks performed by GQECORR
                keypar (output[1]//"[0]", "OBSMODE", silent+)
                if (keypar.value == "IMAGE") {
                    isimaging = yes

                    # make sure fl_keep_qeim is set to no as it's imaging data
                    l_fl_keep_qeim = no
                }

            } else {
                # If it fails exit
                printlog ("ERROR - GIREDUCE: Call to GQECORR returned "//\
                    "a non zero status. Exiting.", l_logfile, verbose+)
                goto crash
            }
        } # End of setting up first part to QE correction

        canflat[1] = l_fl_flat
        canbias = l_fl_bias
        candark = l_fl_dark
        printlog ("", l_logfile, l_verbose)
        if (candark)
            printlog ("If dark subtraction: scale = exptime_input/\
                exptime_dark", l_logfile, l_verbose)

        for (ii=1; ii<=nimages; ii+=1) {

            canflat[1] = l_fl_flat
            canbias = l_fl_bias
            candark = l_fl_dark

            # Obtain DETSEC and CCDSUM for all science extensions in all images
            # Moved due to Hamamatsu changes - CCDSEC and DETSEC are
            # updated after trimming for Hamamatsu. For old CCDs the output
            # science extensions have the original DETSEC and CCDSEC values
            # from the input image so it is OK to use output[ii] in both cases.
            # The loop can be placed here as the checks against ccd[ii,jj]
            # are done further down in the
            # 'if (l_fl_bias || l_fl_flat || l_fl_dark || l_fl_over)' loop.
            for (jj = 1; jj<=n_ext[ii]; jj+=1) {
                imgets (output[ii]//"["//l_sci_ext//","//jj//"]", "DETSEC",
                    >& "dev$null")
                ccd[ii,jj] = imgets.value
                imgets (output[ii]//"["//l_sci_ext//","//jj//"]", "CCDSUM",
                    >& "dev$null")
                ccd[ii,jj] = ccd[ii,jj]//imgets.value
            }

            # BIAS: Check the image and bias frame have the same number of
            # extensions and detector section
            if (canbias) {
                if (n_ext[ii] != n_extbias) {
                    printlog ("WARNING - GIREDUCE: No bias subtraction. "//\
                        input[ii]//" and bias image ("//l_bias//") have \
                        different number of extensions", l_logfile, verbose+)
                    canbias = no
                }
                for (jj=1; jj<=n_ext[ii]; jj+=1) {
                    if (ccd[ii,jj] != ccd_bias[jj]) {
                        printlog ("WARNING - GIREDUCE: Mismatch detector \
                            section or binning of bias", l_logfile, verbose+)
                        printlog ("                    ("//l_bias//") and \
                            image ("//input[ii]//")", l_logfile, verbose+)
                        printlog ("                    No bias subtraction "//\
                            "taking place for "//input[ii], \
                            l_logfile, verbose+)
                        canbias = no
                    }
                }
                # Check the input image has not been bias corrected before
                keypar (output[ii]//"[0]", "BIASIM", silent+)
                if (keypar.found)  {
                    printlog ("WARNING - GIREDUCE: Image "//input[ii]//\
                        " already BIAS subtracted", l_logfile, verbose+)
                    printlog ("                    No bias subtraction"//\
                        "taking place for "//input[ii], \
                        l_logfile, verbose+)
                    canbias = no
                }
            }  # end of canbias checks

            # DARK: Check the image and dark frame have the same number of
            # extensions
            l_darkscale = 0.
            if (candark) {
                if (n_ext[ii]!=n_extdark) {
                    printlog ("WARNING - GIREDUCE: No dark count \
                        correction. "//input[ii]//" and dark image \
                        ("//l_dark//") have different number of extensions",\
                        l_logfile, verbose+)
                    candark = no
                }

                for (jj=1; jj<=n_ext[ii]; jj+=1) {
                    if (ccd[ii,jj] != ccd_dark[jj]) {
                        printlog ("WARNING - GIREDUCE: Can not dark count \
                            correct.", l_logfile, verbose+)
                        printlog ("                    Mismatch detector \
                            section or binning of dark", l_logfile, verbose+)
                        printlog ("("//l_dark//") and image \
                            ("//input[ii]//").", l_logfile, verbose+)
                        candark = no
                    }
                }

                # Check the input image has not been dark corrected before
                keypar (output[ii]//"[0]", "DARKIM", silent+)
                if (keypar.found)  {
                    printlog ("WARNING - GIREDUCE: Can not dark count \
                        correct.", l_logfile, verbose+)
                    printlog ("                    Image "//input[ii]//\
                        " already dark corrected.", l_logfile, verbose+)
                    candark = no
                }
                # Get the exposure time from the image
                keypar (output[ii]//"[0]", l_key_exptime, silent+)
                if (!keypar.found) {
                    printlog ("ERROR - GIREDUCE: Can not dark count \
                        correct.", l_logfile, verbose+)
                    printlog ("                    ("//output[ii]//") does \
                        not have an exposure time in header from keyword "//\
                        l_key_exptime, l_logfile, verbose+)
                    candark = no
                    goto crash
                } else
                    l_darkscale = real(keypar.value)/l_darkexptime

                # Get the Shuffle distance from the image
                keypar (output[ii]//"[0]", l_key_nodpix, silent+)
                if (!keypar.found) {
                    printlog ("ERROR - GIREDUCE: Can not dark count \
                        correct.", l_logfile, verbose+)
                    printlog ("                 ("//output[ii]//") does \
                        not have a shuffle distance in ",l_logfile, verbose+)
                    printlog ("                 header from keyword "//\
                        l_key_nodpix, l_logfile, verbose+)

                    candark = no
                    goto crash
                } else
                    l_nodpix = real(keypar.value)

                # Get the number of nod cycles from the image
                keypar (output[ii]//"[0]", l_key_nodcount, silent+)
                if (!keypar.found) {
                    printlog ("                    Can not dark count \
                        correct.", l_logfile, verbose+)
                    printlog ("ERROR - GIREDUCE: ("//output[ii]//") does \
                        not have number of nod cycles in header from keyword \
                        "// l_key_nodcount, l_logfile, verbose+)
                    candark = no
                    goto crash
                } else
                    l_nodcount = real(keypar.value)

            }  # end of candark checks

            if (globalcanqecorr) {
                ll_qecorr = "INDEF"
                ll_qecorr_int = ""

                # All is good at this point in terms of QE correction
                # Already, checked if BIAS/DARK correcting/corrected.
                # Call to GQECORR has been tested and the globalcanqecorr flag
                # set accordingly. - MS

                # Set the name according to the name stored in ll_corrimg[]
                # See note above at call to gqecorr - MS
                if (num_qecorrimgs == 1) {
                    ll_qecorr_int = ll_corrimg[1]
                } else {
                    ll_qecorr_int = ll_corrimg[ii]
                }

                # Doouble check the correction files exist on disk
                if (isimaging){
                    if (access(ll_qecorr_int)) {
                        ll_qecorr = ll_qecorr_int
                    } else {
                        printlog ("ERROR - GIREDUCE: Cannot access "//\
                            ll_qecorr_int//". Cannot QE correct "//\
                            l_output[ii], l_logfile, verbose+)
                        canqecorr = no
                    }
                } else {
                    if (imaccess(ll_qecorr_int)) {
                        ll_qecorr = ll_qecorr_int
                    } else {
                        printlog ("ERROR - GIREDUCE Cannot access "//\
                            ll_qecorr_int//". Cannot QE correct "//\
                            l_output[ii], l_logfile, verbose+)
                        canqecorr = no
                    }
                }

            } # End of qecorr setup and final tests


            # FLAT checks + one last final check for QE correcting
            if (canflat[1]) {
                # Check the input image has not been flat corrected before
                keypar (output[ii]//"[0]", "FLATIM", silent+)
                if (keypar.found)  {
                    printlog ("WARNING - GIREDUCE: Image "//input[ii]//\
                        " already FLAT FIELD corrected.", l_logfile, verbose+)
                    printlog ("                    No flat division"//\
                        "taking place for "//input[ii], \
                        l_logfile, verbose+)
                    canflat[1] = no
                }
                # Get the filter of the input image
                keypar (input[ii]//"[0]", l_key_filter, silent+)
                thisfilter = keypar.value
                if (!keypar.found) {
                    printlog ("WARNING - GIREDUCE: Filter of "//input[ii]//\
                        " is unknown.", l_logfile, verbose+)
                    printlog ("                    No flat division"//\
                        "taking place for "//input[ii], \
                        l_logfile, verbose+)
                    canflat[1] = no
                }

                # Pick the flat that has the same filter as the input image
                fnum = 0
                for (kk=1; kk<=4; kk+=1) {
                    if (thisfilter==flatfilter[kk]) {
                        fnum = kk
                        kk = 5
                    }
                }

                # Check that FLATS and inputs have the same number of
                # extension
                # KL - Do this check only if number of flat, fnum, is
                # at least one; array[0] is not valid as arrays are
                # 1-indexed.
                if (fnum == 0) {
                    printlog ("WARNING - GIREDUCE: Can not apply flat field \
                        correction to image "//input[ii]//" - no match in \
                        filters", l_logfile, l_verbose)
                    canflat[1] = no
                    # Skip the rest of the checks - there are too many that
                    # rely on fnum and if it's set to 0 it breaks everything -
                    # MS
                    goto FLAT_CHECK_END
                } else if (n_ext[ii] != n_extflat[fnum]) {
                    printlog ("WARNING - GIREDUCE: Can not apply flat field \
                        correction to image "//input[ii]//" - images have \
                        different sizes", l_logfile, verbose+)
                    canflat[1] = no
                }

                for (jj=1; jj<=n_ext[ii]; jj+=1) {
                    if (ccd[ii,jj] != ccd_flat[fnum,jj]) {
                        printlog ("WARNING - GIREDUCE: Mismatch detector \
                            section or binning of flat", l_logfile, verbose+)
                        printlog ("("//l_flat[fnum]//") and image \
                            ("//input[ii]//")", l_logfile, verbose+)
                        canflat[1] = no
                        printlog ("                    No flat division"//\
                            "taking place for "//input[ii], \
                             l_logfile, verbose+)
                    }
                }

                # Check obsmode of image and flat
                keypar (input[ii]//"[0]", "OBSMODE", silent+)
                if (obsmodeflat[fnum] != keypar.value) {
                    printlog ("WARNING - GIREDUCE: OBSMODE parameter of "//\
                        input[ii]//" and flat field do not match",
                        l_logfile,l_verbose)
                    printlog ("                    No flat division "//\
                        "taking place for "//input[ii], \
                        l_logfile, verbose+)
                    canflat[1] = no
                }

                # The case of QE corrected FLATS, non-QE corrected inputs and
                # no request to QE correct the inputs has been checked for
                # already. Also, the case of non-QE corrected FLATS and
                # where 'all' of the inputs have already been QE corrected
                # has already been checked for.
                # Check the status of the QE correction of flat and input
                if (flatqecorr_state[fnum] == inqecorr_state[ii]) {

                    if (canqecorr) {
                        # QE correcting the input data (inqecorr_state = no)
                        # Therefore cannot FLAT divide as FLAT not QE corrected
                        printlog ("WARNING GIREDUCE - Cannot flat field "//\
                            output[ii]//". Flat image "//l_flat[fnum]//\
                            " has not been QE corrected.", \
                            l_logfile, verbose+)
                        canflat[1] = no
                    }
                } else {
                    # Flat is QE corrected, input not (yet?) corrected
                    if (flatqecorr_state[fnum]) {
                        # Cannot QE correct input for some reason; though call
                        # to GQECORR was made
                        if (!canqecorr) {
                            # flat is qe corrected inout is not and not
                            printlog ("WARNING GIREDUCE - Cannot flat "//\
                                "field "//output[ii]//". Flat image "//\
                                 l_flat[fnum]//" has been QE corrected and"//\
                                 " cannot QE correct "//output[ii]//".", \
                                 l_logfile, verbose+)

                            canflat[1] = no
                        }
                    } else if (inqecorr_state[ii]) {
                        # Cannot FLAT field as QE image already QE corrected
                        # but FLAT isn't
                        printlog ("WARNING GIREDUCE - Cannot flat "//\
                            "field "//output[ii]//". Flat image "//\
                             l_flat[fnum]//" has not been QE corrected, "//\
                             " whereas, "//output[ii]//"has been.", \
                             l_logfile, verbose+)

                            canflat[1] = no
                    }
                } # End of the Flat QE corr checks

FLAT_CHECK_END:
            }  # end of canflat checks

            if (!canbias && !canflat[1] && !candark && !canqecorr)
                goto nextimage1

            if (canbias)
                s_bias = l_bias
            else
                s_bias = "INDEF"
            if (candark)
                s_dark = l_dark
            else
                s_dark = "INDEF"
            if (canflat[1])
                s_flat = l_flat[fnum]
            else
                s_flat = "INDEF"

            l_struct = ""
            printf ("%-18s %-20s %-18s %-10s %-6s\n","Output image",
                "Bias", "Flat", "Dark", "Scale") | scan (l_struct)
            printlog (l_struct, l_logfile, l_verbose)
            printf ("%-18s %-20s %-18s %-10s %5.2f\n", output[ii], s_bias,
                s_flat, s_dark, l_darkscale) | scan (l_struct)
            printlog (l_struct, l_logfile, l_verbose)

            for (jj=1; jj<=n_ext[ii]; jj+=1) {

                imgets (output[ii]//"["//l_sci_ext//","//jj//"]", l_key_gain,
                    >& "dev$null")
                if (imgets.value == "0") {
                    printlog ("WARNING - GIREDUCE: Cannot find "//\
                        l_key_gain//" in ext "//str(jj), l_logfile, l_verbose)
                    n_gain[jj] = 1
                } else {
                    n_gain[jj] = real(imgets.value)
                }

                # Get the trim section in the input/output image extension

                imgets (output[ii]//"["//l_sci_ext//","//jj//"]", "TRIMSEC",
                    >& "dev$null")
                if (imgets.value == "0") {
                    fl_tsec_img = no
                    tsec_img = imgets.value
                } else {
                    fl_tsec_img = yes
                    tsec_img = imgets.value
                }
                ll_sciexp = "a"

                if (canbias) {
                    # Check the BIAS and inputs have the same dimensions in
                    # each extension
                    subbias = ""

                    if (fl_tsec_img && fl_tsec_bias[jj] && \
                        (tsec_bias[jj] != tsec_img)) {

                        printlog ("ERROR - GIREDUCE: "//input[ii]//" and \
                            bias image ("//l_bias//") have been trimmed to \
                            different sections in the CCD.",
                            l_logfile, l_verbose)
                        goto nextextension1
                    } else if (!fl_tsec_img && fl_tsec_bias[jj]) {
                        printlog ("ERROR - GIREDUCE: Can not bias subtract.",
                            l_logfile, verbose+)
                        printlog ("                Trim flag is not set but \
                            bias image ("//l_bias//") is trimmed to "//\
                            tsec_bias[jj], l_logfile, verbose+)
                        status = 1
                        goto nextextension1
                    } else if (fl_tsec_img && !fl_tsec_bias[jj]) {
                        printlog ("WARNING - GIREDUCE: "//input[ii]//" is \
                            trimmed but bias image ("//l_bias//") is not.",
                            l_logfile, l_verbose)
                        printlog ("                    Using a section of the \
                            bias image", l_logfile, l_verbose)
                        subbias = tsec_img
                    }

                    # Check if the BIAS was created with the raw BIASes being
                    # overscan trimmed. Check overscan state of output image
                    # too. Only print a warning if they do not match. Only
                    # check first SCI extension as keypar will crash due to
                    # the previous extension not having finished writing
                    # to disk.
                    if (jj == 1) {
                        keypar (l_bias//"["//l_sci_ext//","//jj//"]", \
                            "OVERSEC", silent+)
                        if (keypar.found) {
                            overstate_cal = yes
                        } else {
                            keypar (l_bias//"["//l_sci_ext//","//jj//"]", \
                                "OVERSCAN", silent+)
                            if (keypar.found) {
                                overstate_cal = yes
                            } else {
                                overstate_cal = no
                            }
                        }

                        # Output image
                        keypar (output[ii]//"["//l_sci_ext//","//jj//"]", \
                            "OVERSEC", silent+)
                        if (keypar.found) {
                            overstate_img = yes
                        } else {
                            keypar (output[ii]//"["//l_sci_ext//\
                                ","//jj//"]", "OVERSCAN", silent+)

                            if (keypar.found) {
                                overstate_img = yes
                            } else {
                                overstate_img = no
                            }
                        }

                        # Perform the comparison test
                        if ((overstate_img != overstate_cal)) {
                            if (!overstate_cal) {
                                # Image overscan subtracted but BIAS is not
                                printlog ("ERROR - GIREDUCE: image "//\
                                    input[ii]//\
                                    " has been overscan subtracted."//\
                                    "\n                  "//\
                                    "However, BIAS image "//l_bias//\
                                    " appears"//\
                                    "\n                  "//\
                                    "to have been created from "//\
                                    "biases that were not"//\
                                    "\n                  "//\
                                    "overscan subtracted.", \
                                    l_logfile, verbose+)
                            } else {
                                # Image not overscan subtracted but BIAS is
                                printlog ("ERROR - GIREDUCE: image "//\
                                    input[ii]//" has not been overscan "//\
                                    "subtracted."//\
                                    "\n                  "//\
                                    "However, BIAS image "//l_bias//\
                                    " appears"//\
                                    "\n                  "//\
                                    "to have been created from "//\
                                    "biases that were"//\
                                    "\n                  "//\
                                    "overscan subtracted.", \
                                    l_logfile, verbose+)
                            }
                            goto crash
                        }
                    }

                    ll_bias = l_bias//"["//l_sci_ext//","//jj//"]"//subbias
                    ll_sciexp = ll_sciexp//"-b"
                    # DQ frame is generated as using BPM to set 1 bit and
                    # saturation level to set 4 bit
                } else {
                    ll_bias = "INDEF"
                }

                if (candark) {
                    # Check that DARK and inputs have the same dimensions
                    # in each extension

                    # darkimg = l_dark//"["//l_sci_ext//","//jj//"]"

                    subdark = ""

                    if (fl_tsec_img && fl_tsec_dark[jj] && \
                        (tsec_dark[jj] != tsec_img)) {
                        printlog ("ERROR - GIREDUCE: "//input[ii]//" and \
                            dark image ("//l_dark//") have been trimmed to \
                            different sections in the CCD.",
                            l_logfile, l_verbose)
                        goto nextextension1
                    } else if (!fl_tsec_img && fl_tsec_dark[jj]) {
                        printlog ("ERROR - GIREDUCE: Can not dark count \
                            correct.", l_logfile, verbose+)
                        printlog ("                Trim flag is not set but \
                            dark image ("//l_dark//") is trimmed to "//\
                            tsec_dark[jj], l_logfile, l_verbose)
                        status = 1
                        goto nextextension1
                    } else if (fl_tsec_img && !fl_tsec_dark[jj]) {
                        printlog ("WARNING - GIREDUCE: "//input[ii]//" is \
                            trimmed but dark image ("//l_dark//") is not.",
                            l_logfile, l_verbose)
                        printlog ("                    Using a section of the \
                            dark image", l_logfile, l_verbose)
                        subdark = tsec_img
                    }

                    if (int(l_darknodpix) != int(l_nodpix)) {
                        printlog ("ERROR - GIREDUCE: Can not dark count \
                            correct.", l_logfile, verbose+)
                        printlog ("                    "//input[ii]//" and \
                            dark image ("//l_dark//") have different shuffle \
                            distances.",
                            l_logfile, verbose+)
                        status = 1
                        goto nextextension1
                    }

                    if (real(l_darknodcount) != real(l_nodcount)) {
                        printlog ("ERROR - GIREDUCE: Can not dark count \
                            correct.", l_logfile, verbose+)
                        printlog ("                    "//input[ii]//" and \
                            dark image ("//l_dark//") have different shuffle \
                            cycles.",
                            l_logfile, l_verbose)
                        status = 1
                        goto nextextension1
                    }

                    ll_dark = l_dark//"["//l_sci_ext//","//jj//"]"//subdark
                    ll_sciexp = ll_sciexp//"-c*"//str(l_darkscale)

                } else {
                    ll_dark = "INDEF"
                }

                # Variables for QE correction
                if (canqecorr) {
                    # Set sciexp
                    ll_sciexp = "("//ll_sciexp//")"//qecorrapp//"f"
                    # qecorrapp is how to apply the correction - set internally
                    # in gireduce (matched to how gqecorr does it)
                    # f = output from gqecorr (image/number)
                    # If image, no BPM, places in image not to be corrected
                    # set to one in the correction image - spectroscopy only.
                    # - MS

                    if (!isimaging) {
                        # Set extension of spectroscopic image to use
                        ll_qecorr = ll_qecorr_int//\
                            "["//l_sci_ext//","//jj//"]"
                        if (jj == 1) {
                            qecorrmssg  = "GIREDUCE: Used "//ll_qecorr_int//\
                                " to QE correct "//output[ii]

                            # Write QE correction image to phu of output
                            gemhedit (output[ii]//"[0]", \
                                l_key_qecorrim, \
                                ll_qecorr_int, l_qecorrim_comment, \
                                delete-)

                        }
                    } else {
                        # Parse out the QE factor using CCDNAME for imaging
                        keypar (output[ii]//\
                            "["//l_sci_ext//","//jj//"]", \
                            "CCDNAME", silent+)

                        if (keypar.found) {
                            # Due to CCDNAME naming convention the
                            # gmosQEfactors.dat file has the CCDNAME stored
                            # underscores not spaces in the string. So need
                            # to add them
                            print (keypar.value) | translit ("STDIN", \
                                from_string=" ", to_string="_", delete=no, \
                                collapse=no) | \
                                    scan(ccdname)

                            # Read the QE factor
                            match (ccdname, ll_qecorr_int, \
                                metacharacters=no, stop=no, \
                                print_file_names=no) | \
                                fields ("STDIN", 1, lines="1-",\
                                    quit_if_missing=no, \
                                    print_file_names=no) |\
                                        scanf ("%s", ll_qecorr)

                            # Write factor to header of current output
                            # extension
                            gemhedit (output[ii]//"["//l_sci_ext//\
                                ","//jj//"]", l_key_qeimfactor, \
                                ll_qecorr, l_qeimfactor_comment, \
                                delete-)
                        } else {
                            # Crash
                            printlog ("ERROR - GIREDUCE: "//\
                                "CCDNAME keyword not found in "//\
                                output[ii]//\
                                "["//l_sci_ext//","//jj//"]", \
                                l_logfile, verbose+)
                            goto crash
                        }

                        printlog ("GIREDUCE: QE correction factor for "//\
                            "extension ["//l_sci_ext//","//jj//"] is: "//\
                            str(1/real(ll_qecorr)), \
                            l_logfile, l_verbose)
                    }

                } else {
                    # Set up QE correction here
                     ll_qecorr = "INDEF"
                } # End of canqecorr loop


                if (canflat[1]) {
                    # Check the FLATS and inputs have the same dimensions
                    # in each extension

                    #flatimg = l_flat[fnum]//"["//l_sci_ext//","//jj//"]"

                    subflat = ""

                    if (fl_tsec_img && fl_tsec_flat[jj,fnum] && \
                        (tsec_flat[jj,fnum] != tsec_img)) {

                        printlog ("WARNING - GIREDUCE: "//input[ii]//" and \
                            flat image ("//l_flat[fnum]//") have been trimmed \
                            to different sections in the CCD.",
                            l_logfile, l_verbose)
                        goto nextextension1
                    } else if (!fl_tsec_img && fl_tsec_flat[jj,fnum]) {
                        printlog ("ERROR - GIREDUCE: Can not apply flat field \
                            correction.", l_logfile, verbose+)
                        printlog ("                Trim flag is not set but \
                            flat image ("//l_flat[fnum]//") is trimmed to "//\
                            tsec_flat[jj,fnum], l_logfile, l_verbose)
                        status = 1
                        goto nextextension1
                    } else if (fl_tsec_img && !fl_tsec_flat[jj,fnum]) {
                        printlog ("WARNING - GIREDUCE: "//input[ii]//" is \
                            trimmed but flat image ("//l_flat[fnum]//") is \
                            not", l_logfile, l_verbose)
                        printlog ("                    Using subsection of \
                            the flat image", l_logfile, l_verbose)
                        subflat = tsec_img
                    }
                    ll_flat = l_flat[fnum]//"["//l_sci_ext//","//jj//"]"//\
                        subflat
                    ll_sciexp = "(("//ll_sciexp//")/d)"

                } else {
                    ll_flat = "INDEF"
                }

                # The bias image is:
                #      l_bias//"["//l_sci_ext//","//jj//"]"//tsec_img
                # The flat field image is:
                #      l_flat[fnum]//"["//l_sci_ext//","//jj//"]"//tsec_img

                imexpr (ll_sciexp,
                    output[ii]//"["//l_sci_ext//","//jj//",overwrite]",
                    output[ii]//"["//l_sci_ext//","//jj//"]",
                    ll_bias, ll_dark, ll_flat, n_gain[jj], ll_qecorr, \
                    dims="auto",
                    intype="real", outtype="real", refim="auto", bwidth=0,
                    btype="nearest", bpixval=0., rangecheck+, verbose-,
                    exprdb="none")

                # !!!!!!!!!!!!START VAR/DQ LOOP!!!!!!!!!!!!!!
                if (l_fl_vardq) {
                    # SCI should be the same number in the BPM
                    if ((n_ext[ii] != next_bpm) && (l_bpm != "")) {
                        printlog ("ERROR - GIREDUCE: Number of \
                        "//l_sci_ext//" extensions in "//img, l_logfile, \
                            verbose+)
                        printlog ("                  is not the same as in \
                            the Bad Pixel Mask." , l_logfile, verbose+)
                        imdelete (output[ii], verify-, >& "dev$null")
                        goto crash
                    }
                    tmpd[jj] = mktemp ("tmpd")
                    tmpdq = mktemp ("tmpdq")
                    imgets (output[ii]//"["//l_sci_ext//","//jj//"]",
                        l_key_ron, >& "dev$null")
                    if (imgets.value != "0")
                        l_ron = real(imgets.value)
                    else {
                        printlog ("WARNING - GIREDUCE: Read noise not found \
                            in image header.", l_logfile, verbose+)
                        printlog ("                    Using a read noise of \
                            "//l_ron//" electrons.", l_logfile, verbose+)
                    }

                    # Create DQ frame and combine with the BPM if one is given
                    # Identify saturated pixels

                    # The use of the satfrac is variable is for Hamamatsu
                    # data only; this is because the overscan subtraction is
                    # not a straight line (i.e., order > 1). Ideally, one would
                    # directly use the overscan curve instead of using the
                    # OVERSCAN keyword. Basically this is a horrible hack. This
                    # hack was chosen over a complete re-write of gireduce and
                    # gsat - MS

                    ##M Should be changed to '>='
                    ##M This calculation should also include the bias image!
                    dqexp = "((a > b) ? 4 : 0)"

                    # Adjust saturation level if user supplies one
                    if (l_sat != "default") {
                        # Option for user-defined value
                        satvalue = real(l_sat)
                    } else {
                        gsat (output[ii], extension="["//l_sci_ext//\
                            ","//jj//"]", gaindb=l_gaindb, bias="static", \
                            pixstat="midpt", statsec="default", \
                            gainval="default", scale=satfrac, \
                            logfile=l_logfile, verbose=no)
                        if (gsat.status != 0) {
                            printlog ("ERROR - GIREDUCE: GSAT returned a "//\
                                "non-zero status. Exiting", \
                                l_logfile, verbose+)
                            status = 1
                            goto crash
                        }
                        satvalue = real(gsat.saturation)
                    }

                    imexpr (dqexp, tmpdq, output[ii]//"["//l_sci_ext//",\
                        "//jj//"]", satvalue, dims="auto", intype="auto",\
                        outtype="ushort", refim="auto", bwidth=0, \
                        btype="nearest", bpixval=0., rangecheck+, verbose-, \
                        exprdb="none", lastout="", mode="ql")

                    # Check if correspond to the same region of the CCD
                    # (DATASEC and TRIMSEC keywords)
                    # Needs to be modified to handle windowed readout images.

                    if (l_bpm != "") {
                        imgets (output[ii]//"["//l_sci_ext//","//jj//"]",
                            l_key_datasec, >& "dev$null")
                        sci_dsec = imgets.value
                        imgets (l_bpm//"["//l_dq_ext//","//jj//"]", \
                            l_key_datasec, >& "dev$null")
                        bpm_dsec = imgets.value

                        # will exit here if image is in windowed readout mode.
                        if (sci_dsec != bpm_dsec) {
                            printlog ("ERROR - GIREDUCE: Image "//img//" and \
                                bad pixel mask "//l_bpm, l_logfile,verbose+)
                            printlog ("                  contain different \
                                data sections.", l_logfile, verbose+)
                            imdelete (output[ii]//","//tmpdq, verify-,
                                >& "dev$null")
                            goto crash
                        }

                        imgets (output[ii]//"["//l_sci_ext//","//jj//"]",
                            "TRIMSEC", >& "dev$null")
                        sci_tsec = imgets.value

                        imgets (l_bpm//"["//l_dq_ext//","//jj//"]", "TRIMSEC",
                            >& "dev$null")
                        bpm_tsec = imgets.value

                        if ((bpm_tsec == "0") && (sci_tsec != "0")) {
                            tsec = sci_tsec
                        } else if ((sci_tsec == "0") && (bpm_tsec != "0")) {
                            tsec = bpm_tsec
                        } else if ((sci_tsec != "0") &&
                            (sci_tsec != bpm_tsec)) {
                            printlog ("ERROR - GIREDUCE: Image "//img//" and \
                                bad pixel mask "//l_bpm, l_logfile, verbose+)
                            printlog ("                  contain different \
                                trimming sections.", l_logfile, verbose+)
                            imdelete (output[ii]//","//tmpdq, verify-,
                                >& "dev$null")
                            goto crash
                        } else
                            tsec = ""

                        addmasks (tmpdq//tsec//","//l_bpm//"["//l_dq_ext//",\
                            "//jj//"]"//tsec, tmpd[jj]//".fits","im1 || im2")

                        if (tsec != "" ) {
                            gemhedit (tmpd[jj], "TRIMSEC", tsec, "", delete-)
                        }

                        gemhedit (tmpd[jj], "BPMFILE", name_of_bpm//\
                            "["//l_dq_ext//","//jj//"]", \
                             "User defined input BPM", delete-)

                        imdelete (tmpdq, verify-, >& "dev$null")

                    } else
                        imrename (tmpdq, tmpd[jj]//".fits", >& "dev$null")

                    # Define the expressions for imexpr -- To start:
                    # Make initial variance image based on photon statistics
                    # if fl_over+, science has been overscan subtracted already
                    # (max(counts,0.)/gain) counts must only be actual
                    # photons detected (not bias/dark counts)

                    ll_varexp ="((a/b)**2)" # This is in ADU squared
                    # add small overscan error
                    if (l_fl_over)
                        ll_varexp = ll_varexp//"+i"

                    # i is the variance from the overscan subtraction, zero
                    # if not overscan subtracting
                    dqlist = tmpd[jj]//".fits"
                    ll_ndq = 1
                    imgets (output[ii]//"["//l_sci_ext//","//jj//"]",
                        "GAINMULT", >& "dev$null")
                    if ((ll_bias == "INDEF") ||
                        (!imaccess(l_bias//"["//l_var_ext//","//jj//"]"))) {
                        # image was previously bias subtracted.
                        ll_varbias = "INDEF"

                        # We must assume bias variance was already added to
                        # science variance
                        ll_dqbias = "INDEF"
                        # using previously bias subtracted image for photon
                        # counts
                        if (imgets.value == "0") {
                            # All in ADU still
                            ll_varexp = ll_varexp//" + max(c,0.0)/b"
                        } else {
                            #convert back to ADU
                            ll_varexp = ll_varexp//" + max(c,0.0)/(b**3)"
                        }

                    } else if (imaccess(l_bias//"["//l_var_ext//","//jj//"]")){
                        ll_varbias = l_bias//"["//l_var_ext//","//jj//"]"//\
                            subbias
                        # If bias subtracting sci frame already bias subtracted
                        # by this point, so just use c, not c-d. - MS
                        ll_varexp = ll_varexp//" + max(c,0.0)/b + e"
                        dqlist = dqlist//","//l_bias//"["//l_dq_ext//","\
                            //jj//"]"//subbias
                        ll_ndq += 1
                    }

                    if (ll_dark == "INDEF") {
                        ll_vardark = "INDEF"
                        ll_dqdark = "INDEF"
                    } else if (imaccess(l_dark//"["//l_var_ext//","//jj//"]")){
                        ll_vardark = l_dark//"["//l_var_ext//","//jj//"]"//\
                            subdark
                        if (canbias) {
                            # science minus artificial bias level will give
                            # good estimate of photons detected
                            ll_varexp = ll_varexp//"+j*"//str(l_darkscale)//\
                                "*"//str(l_darkscale)
                        } else {
                            # science minus artificial dark level will give
                            # good estimate of photons detected
                            # If dark subtracting sci frame already dark
                            # subtracted by this point, so just use c, not
                            # c-h. - MS
                            ll_varexp = ll_varexp//" + max(c,0.0)/b + j*"//\
                                str(l_darkscale)//"*"//str(l_darkscale)
                            dqlist = dqlist//","//l_dark//"["//l_dq_ext//","\
                                //jj//"]"//subdark
                            ll_ndq += 1
                        }
                    } else {
                        #dark did not have var/dq attached
                        ll_vardark = "INDEF"
                        ll_dqdark = "INDEF"
                    }


                    # Set up QE correction variance input here
                    ll_varqecorr = "INDEF"
                    if (canqecorr) {
                        ll_varexp = "(("//ll_varexp//")"//qecorrapp//"(k**2))"
                        # k is the output from gqecorr
                        # qecorrapp is how to apply the qe correection
                        # They are scalar values, either an image or value
                        ll_varqecorr = ll_qecorr
                    }

                    if (ll_flat == "INDEF") {
                        ll_varflat = "INDEF"
                        ll_dqflat = "INDEF"
                        ll_varexp = ll_varexp
                    } else {
                        ll_varflat = l_flat[fnum]//"["//l_var_ext//","\
                            //jj//"]"//subflat
                        if (ll_dark == "INDEF" && imgets.value=="0")
                            # Science frame has aleady been bias subtracted -MS
#                                (g/(b**2)*((c-d)**2))/(f**4))"
                            ll_varexp = "((("//ll_varexp//")/(f*f)) + \
                                (g/(b**2)*((c)**2))/(f**4))"
                        else if (ll_dark != "INDEF" && imgets.value=="0")
                            # Science frame has aleady been dark subtracted -MS
#                                ((g/(b**2))*((c-h)**2))/(f**4))"
                            ll_varexp = "((("//ll_varexp//")/(f**2)) + \
                                ((g/(b**2))*((c)**2))/(f**4))"
                        dqlist = dqlist//","//l_flat[fnum]//"["//l_dq_ext//\
                            ","//jj//"]"//subflat
                        ll_ndq += 1
                    }

                    # Variance expression
                    # If varflat created with gsflat or giflat it is in
                    # electrons, need to convert back to ADU for now
                    # (((ron/gain)**2 + max(sci-bias,0.0)/gain) + varover +
                    # varbias + vardark * darkscale**2)/(flat**2)) +
                    # ((varflat/(gain**2))*(sci-bias(or dark))**2)/(flat**4))
                    imexpr (ll_varexp, \
                        output[ii]//"["//l_var_ext//","//jj//",append]", \
                        l_ron, n_gain[jj], output[ii]//"["//l_sci_ext//"," \
                        //jj//"]", ll_bias, ll_varbias, ll_flat, ll_varflat, \
                        ll_dark, overvar[jj,ii], ll_vardark, ll_varqecorr,
                        dims="auto", \
                        intype="real", outtype="real", refim="auto", \
                        bwidth=0, btype="nearest", bpixval=0., \
                        rangecheck=yes, verbose=no, exprdb="none")

                    # Add up DQ planes using addmasks
                    ll_dqexp = "im1"
                    for (kk=2; kk<=ll_ndq; kk+=1)
                        ll_dqexp = ll_dqexp//" || im"//str(kk)

                    addmasks (dqlist, tmppl//".pl", ll_dqexp)

                    # This is a to get it to output 16-bit with a
                    # reduced set of header keywords - MS
                    imarith (tmppl//".pl", "+", 0, \
                        output[ii]//"["//l_dq_ext//","\
                        //jj//",append]", pixtype="ushort", verbose-)

                    imdelete (tmppl//".pl", verify-, >& "dev$null")

                    # Add keywprds to header - MS
                    for (kk = 1; kk <= hdw_num; kk += 1) {
                        # Update the section keywords
                        keypar (output[ii]//"["//l_sci_ext//","//jj//"]", \
                            hdwords_to_read[kk], silent+)
                        if (keypar.found) {
                            gemhedit (output[ii]//"["//l_dq_ext//","//jj//"]",\
                                hdwords_to_read[kk], keypar.value, \
                                hdwords_comments[kk], delete-)
                        }
                    }

                    # Add BPM name if present - MS
                    keypar (tmpd[jj]//".fits", "BPMFILE", silent+)
                    if (keypar.found) {
                        gemhedit (output[ii]//"["//l_dq_ext//","//jj//"]",\
                            "BPMFILE", keypar.value, "User defined BPM", \
                            delete-)
                    }

                    imdelete (tmpd[jj], verify-, >& "dev$null")

                } #end vardq loop
            } #end of jj extension loop over

            #if (canflat[ii])
            if (canflat[1])
                flatmssg = "GIREDUCE: Divided image "//output[ii]//" by \
                    flatfield "//l_flat[fnum]

            biasmssg = "GIREDUCE: Subtracted bias "//l_bias//" from \
                "//output[ii]
            darkmssg = "GIREDUCE: Subtracted dark " //l_dark//" from \
                "//output[ii]

            if (strstr("b", ll_sciexp) > 0 ) printlog (biasmssg, l_logfile, \
                l_verbose)
            if (strstr("c", ll_sciexp) > 0 ) printlog (darkmssg, l_logfile, \
                l_verbose)
            if (canqecorr && isimaging) {
                printlog ("GIREDUCE: QE corrected using displayed "//\
                    "coefficients", \
                    l_logfile, l_verbose)
            } else if (canqecorr && !isimaging) {
                printlog (qecorrmssg, l_logfile, l_verbose)
            }
            if (strstr("d", ll_sciexp) > 0 ) printlog (flatmssg, l_logfile, \
                l_verbose)

            if (l_fl_mult) {
                imgets (output[ii]//"["//l_sci_ext//", 1]", "GAINMULT", \
                    >& "dev$null")
                if (imgets.value == "0") {
                    printlog ("GIREDUCE: multiplying image "//output[ii]//" \
                        by gain", l_logfile, l_verbose)
                    flpr # bug is a problem for multiple images with \
                         #fl_bias+, used a magic flpr
                    ggain (output[ii], gaindb=l_gaindb, logfile=l_logfile, \
                        key_gain=l_key_gain, key_ron=l_key_ron, gain=l_gain, \
                        ron=l_ron, fl_update+, fl_mult=l_fl_mult, verbose=no, \
                        sci_ext=l_sci_ext, var_ext=l_var_ext)
                }
                printlog ("GIREDUCE: output counts in electrons", l_logfile,
                    l_verbose)
            } else {
                printlog ("GIREDUCE: output counts in ADU\n", l_logfile,
                    l_verbose)
            }
nextextension1:

            if (canbias)
                gemhedit (output[ii]//"[0]", "BIASIM", l_bname,
                    "Bias image used by GIREDUCE", delete-)
            if (candark)
                gemhedit (output[ii]//"[0]", "DARKIM", l_dname,
                    "Dark image used by GIREDUCE", delete-)
            if (canflat[1])
                gemhedit (output[ii]//"[0]", "FLATIM", l_fname[fnum],
                    "Flat field used by GIREDUCE", delete-)
nextimage1:
        } # End for(ii=1; i<=nimages ...
    }

finish:
    # Copy MDF
    for (ii=1; ii<=nimages; ii+=1) {
        if (hasmdf[ii]) {
            tmpmdf = mktemp ("tmpmdf")
            tcopy (input[ii]//".fits[MDF]", tmpmdf//".fits", verbose-)
            fxinsert (tmpmdf//".fits", output[ii]//".fits["//n_ext[ii]//"]",
                "1", verbose-, >& "dev$null")
            delete (tmpmdf//".fits", verify-, >& "dev$null")
        }
    }

    # Update the headers: GIREDUCE, GEM-TLM
    gemdate ()
    for (ii=1; ii<=nimages; ii+=1) {
        gemhedit (output[ii]//"[0]", "GIREDUCE", gemdate.outdate,
            "UT Time stamp for GIREDUCE", delete-)
        gemhedit (output[ii]//"[0]", "GEM-TLM", gemdate.outdate,
            "UT Last modification with GEMINI", delete-)
    }
    goto clean

    # Exit with error
crash:
    status = 1
    for (ii=1; ii<=nimages; ii+=1) {
        if (n_do[ii]) imdelete (output[ii], verify-, >& "dev$null")
    }

    # Clean up and exit
clean:
    scanfile = ""

    date | scan (l_date_struct)
    printlog ("", l_logfile, l_verbose)
    printlog ("GIREDUCE - Cleaning up -- "//l_date_struct, \
        l_logfile, l_verbose)

    # Clean up QE files created by GQECORR if not wanted and
    # the list of their names
    if (access(tmpqecorrimgs)) {
        delete (tmpqecorrimgs, verify-, >& "dev$null")
    }

    if (globalcanqecorr) {
        if (!l_fl_keep_qeim) {
            for (ii = 1; ii <= num_qecorrimgs; ii +=1) {
                if (!ll_corrimg_exists[ii]) {
                    delete (ll_corrimg[ii], verify-, >& "dev$null")
                }
            }
        }
    }

    delete (tmpfile//","//tmpinlist//","//tmpoutlist//","//tmplog, verify-,
        >& "dev$null")
    printlog ("", l_logfile, l_verbose)

    date | scan (l_date_struct)

    printlog ("GIREDUCE -- "//l_date_struct, l_logfile, l_verbose)

    if (status == 0) {
        printlog ("GIREDUCE exit status:  good.", l_logfile, l_verbose)
        printlog ("----------------------------------------------------------\
            ------------------", l_logfile, l_verbose)
    } else {
        printlog ("GIREDUCE exit status:  error.", l_logfile, l_verbose)
        printlog ("-----------------------------------------------------------\
            -----------------", l_logfile, l_verbose)
    }
end

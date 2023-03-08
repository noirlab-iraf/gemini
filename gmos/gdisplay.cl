# Copyright(c) 2002-2017 Association of Universities for Research in Astronomy, Inc.

procedure gdisplay (image, frame)

# Display GMOS (raw) images
#
# Version Feb 28, 2002   IJ  v1.3 release
#         Apr 26, 2002   IJ  change warning comment, handle [1] images
#         Jul  7, 2002   IJ  handle dir-names+short image names correctly
#         Sept 20, 2002  IJ   v1.4 release
#         Oct 17, 2002   BM  generalized for GMOS-S, 6amps, and all readmodes
#         Oct 22, 2002   BM  fix problem with defining l_gainmode when fl_bias-
#         Jan 31, 2003   BM  Added output file
#         Sep 13, 2006   JH  Add image number feature as input image if current UT
#         Aug 02, 2009   AS  Set saturation limits based on binning
#         Sep 22, 2009   AS  Set saturation limits based on twilight
#

char    image       {prompt="GMOS image to display, can use number if current UT"} # OLDP-1-input-primary-none
int     frame       {1,prompt="Frame to write to"}                                 # OLDP-1
char    output      {"",prompt="Save pasted file to this name if not blank"}       # OLDP-3
bool    fl_paste    {no,prompt="Paste images to one for imexamine"}                # OLDP-2
bool    fl_bias     {no,prompt="Rough bias subtraction"}                           # OLDP-2
char    rawpath     {"adata$",prompt="Path for input image if not included in name"}
char    gap         {"default",prompt="Size of the gap between the CCDs (in pixels)"}     # OLDP-3
real    z1          {0.,prompt="Lower limit if not autoscaling"}                   # OLDP-2
real    z2          {0.,prompt="Upper limit if not autoscaling"}                   # OLDP-2
bool    fl_sat      {no,prompt="Flag saturated pixels"}                            # OLDP-2
char    satvalue    {"default",prompt="Saturation value"}                            # OLDP-2
bool    fl_imexam   {yes,prompt="If possible, run imexam"}
real    signal      {INDEF,prompt="Flag pixels with signal above this limit"}      # OLDP-2
bool    fl_pretty   {no,prompt="Acquitions only - level amplifiers"}
bool    ret_roi     {yes,prompt="Return only the data within specified ROIs"}
int     req_roi     {0,prompt="Requested ROI (0=all)"}
char    bias_type   {"default", enum="default|calc|static", prompt="Bias value type to determine."}
char    sci_ext     {"all",prompt="Name of extension(s) to display"}               # OLDP-3
char    observatory {")_.observatory", prompt="Observatory (gemini-north or gemini-south)"}
char    prefix      {"auto",prompt="File prefix, (N/S)YYYYMMDDS if not auto"}
char    key_detsec  {"DETSEC",prompt="Header keyword for detector section"}        # OLDP-3
char    key_datasec {"DATASEC",prompt="Header keyword for data section"}           # OLDP-3
char    key_biassec {"BIASSEC",prompt="Header keyword for overscan section"}             # OLDP-3
char    key_ccdsec  {"CCDSEC",prompt="Header keyword for CCD section"}  # OLDP-3
char    key_ccdsum  {"CCDSUM",prompt="Header keyword for CCD binning"}  # OLDP-3
char    gaindb      {"default",prompt="Database with gain data"}    # OLDP-2
bool    verbose     {no,prompt="Verbose"}                                          # OLDP-4
char    logfile     {"",prompt="Logfile"}
int     status      {0,prompt="Exit status (0=good)"}
struct  *scanfile   {"",prompt="Internal use only-."}                               # OLDP-4

#####

begin

    char    l_image = ""
    char    l_output = ""
    char    l_rawpath = ""
    char    l_sci_ext = ""
    char    l_observatory = ""
    char    l_prefix = ""
    char    l_key_detsec = ""
    char    l_key_datasec = ""
    char    l_key_biassec = ""
    char    l_key_ccdsec = ""
    char    l_key_ccdsum = ""
    char    l_gaindb = ""
    char    l_bias_type
    bool    l_verbose, l_fl_paste, l_fl_bias, l_fl_sat, l_fl_imexam
    bool    l_ret_roi
    int     l_frame, l_req_roi
    real    l_z1, l_z2, l_signal

    char    tmpin, l_dsec, l_datasec, l_stdimage, junk, tmpout, l_redirect
    char    l_imext, display_mask, l_gap, l_logfile
    char    l_var_ext, l_dq_ext, l_mdf_ext, pr_str, l_satval
    char    extn, satextn, signalextn, sciextn, l_ccdsum, outextn
    char    lastchar, siteprefix, utdate, imstring, fileprefix
    char    l_sat_ext, l_sig_ext, inextn, msk_extn, disp_msk_extn, gem_sci_ext
    char    sat_image, sig_image, chk_mask, extn_root

    string  dbpm, doverlay, l_express

    int     n_ext, Xmax, Ymax, x1, x2, y1, y2, n_k  #j
    int     Xoff, Yoff, Xbin, Ybin, gapvalue, Xlen, xcenval, pr_len
    int     Ymaxtrim, row_removal_value, det_roi_y_start
    int     l_gapvalue, counter, l_nsample, xwidth, ywidth, ndetroi
    int     l_min, l_max, l_max_nsample, l_nbiascontam

    real    l_xcenter, l_ycenter, l_mag, imX, imY, biasval, chk_value
    real    l_stand_sat, l_reduced_sat

    bool    fl_erase, fl_zscale, fl_zrange, l_fl_mag, l_fl_chmag, l_fl_sec
    bool    ishamamatsu, l_fl_pretty, profile, processed, usestatic, istiled
    bool    fl_check_masks, not_said, iscut, ismosaiced, debug
    bool    app_masks, ise2vDD

    struct  sdate

    # Default parameters
    status = 0
    ise2vDD = no

    # For display command - the number of pixels to sample for z values
    l_max_nsample = 5*250*250
    l_nsample = 1000

    # For saturation / sigmal masks
    fl_check_masks = yes
    not_said = yes
    sat_image = ""
    sig_image = ""
    display_mask = ""

    # Used when creating an output image
    l_sat_ext = "DQSAT"
    l_sig_ext = "DQSIGN"

    # For new GMOS-N Hamamatsu CCDs: number of unused rows at the bottom of the
    # detector
    row_removal_value = 48
    ishamamatsu = no

    # GMULTIAMP
    l_var_ext = "VAR"
    l_dq_ext = "DQ"
    l_mdf_ext = "MDF"
    l_nbiascontam = 4  ##M NBIASCONTAM
    processed = no

    # Temporary variable names
    tmpin = mktemp ("tmpin")
    tmpout = mktemp ("tmpout")

    # Initialize magnification
    l_fl_chmag = yes  # check magnification on first round
    l_fl_mag = no     # initialize check for changed magnification

    # Append saturation and signal masks to output image?
    app_masks = no

    l_verbose = verbose
    debug = no
    profile = no
    if (debug) {
        profile = yes
        l_verbose = yes
    }

    # set the redirect for output to screen if not requested
    if (l_verbose) {
        l_redirect = "STDOUT"
    } else {
        l_redirect = "dev$null"
    }
            
    # Read user parameters
    junk = fscan (image,       l_image)
    junk = fscan (output,      l_output)
    junk = fscan (rawpath,     l_rawpath)
    junk = fscan (observatory, l_observatory)
    junk = fscan (prefix,      l_prefix)
    junk = fscan (key_detsec,  l_key_detsec)
    junk = fscan (key_datasec, l_key_datasec)
    junk = fscan (key_biassec, l_key_biassec)
    junk = fscan (key_ccdsec,  l_key_ccdsec)
    junk = fscan (key_ccdsum,  l_key_ccdsum)
    junk = fscan (gaindb,      l_gaindb)
    junk = fscan (sci_ext,     l_sci_ext)
    junk = fscan (gap,         l_gap)
    l_bias_type = bias_type
    l_req_roi = req_roi
    l_ret_roi = ret_roi
    l_satval = satvalue
    l_fl_pretty = fl_pretty
    l_logfile = logfile

    if (profile) {
        date ("+%H:%M:%S.%N%n") | scan (sdate)
        printlog ("____GDISPLAY Started: "//sdate, l_logfile, verbose+)
    }

    if (l_logfile == "") {
        l_logfile = "dev$null"
    } else if (stridx(" ", l_logfile) > 0) {
        l_logfile = gmos.logfile
        if ((l_logfile == "") || (stridx (" ", l_logfile) > 0)) {
            l_logfile = "gmos.log"
            printlog ("WARNING - GDISPLAY: Both gdisplay.logfile and "//\
                "gmos.logfile fields are empty or not set properly.", \
                l_logfile, verbose+)
            printlog ("                   Using default file gmos.log", \
                l_logfile, verbose+)
        }
    }

    # Parse sci_ext
    print (l_sci_ext) | \
        translit ("STDIN", "a-z", "A-Z", delete-, collapse-) | \
        scan (l_sci_ext)

    l_z1 = z1
    l_z2 = z2
    l_frame = frame
    l_fl_paste = fl_paste
    l_fl_bias = fl_bias
    l_fl_sat = fl_sat
    l_fl_imexam = fl_imexam
    l_signal = signal

    # Print start time
    gemdate (zone="local")
    printlog ("\nGDISPLAY - Started: "//gemdate.outdate//"\n", l_logfile, \
        verbose+)

    # Set z-scale / range parameters for display command
    fl_zscale = yes
    fl_zrange = yes
    if ((l_z1 != 0.) || (l_z2 != 0.)) {
        fl_zscale = no
        fl_zrange = no
    }

    cache ("imgets", "gimverify", "fparse", "gemisnumber", "getfakeUT", \
        "gmultiamp")

    if (debug) {
        printlog ("...checking that the paths have a trailing slash...", \
            l_logfile, verbose+)
    }

    lastchar = substr(l_rawpath, strlen(l_rawpath), strlen(l_rawpath))

    if ((l_rawpath != "") && (lastchar != "/") && (lastchar != "$"))
        l_rawpath = l_rawpath//"/"

    #-----------------------------------------------------------------------
    # Construct image name

    # Check if it's a number

    gemisnumber (l_image, ttest="integer", verbose=no)

    if (gemisnumber.fl_istype) {
        printlog ("GDISPLAY: Constructing image name based on today's UT \
            date or given file prefix.", l_logfile, l_verbose)

        # Build file prefix in 'auto' mode
        if (l_prefix == "auto") {

            #### Convert observatory value to lower case
            l_observatory = strlwr (l_observatory)

            #### Determine site prefix
            if (l_observatory == "gemini-north") {
                siteprefix = "N"
            } else if (l_observatory == "gemini-south") {
                siteprefix = "S"
            } else {
                printlog ("ERROR - GDISPLAY: Parameter 'observatory' must be \
                    either 'gemini-north' or 'gemini-south'", \
                    l_logfile, verbose+)
                status = 121
                goto crash
            }

            #### Build and set file prefix
            getfakeUT()
            utdate = getfakeUT.fakeUT

            fileprefix = siteprefix//utdate//"S"
        } else {
            fileprefix = l_prefix
        }

        printlog ("GDISPLAY: prefix = "//fileprefix, l_logfile, l_verbose)

        printf ("%04d\n", int(l_image)) | scan(imstring)

        l_image = fileprefix//imstring//".fits"

        printlog ("GDISPLAY: image = "//l_image, l_logfile, verbose+)
    }
    #-----------------------------------------------------------------------

    # Is there a directory name
    gimverify (l_image)
    fparse (l_image)
    if (fparse.directory == "") {
        gimverify (l_image) # current directory

        if (gimverify.status != 0) {
            gimverify (l_rawpath//l_image) # try rawpath
        }
    }

    if (gimverify.status == 1) {
        printlog ("ERROR - GDISPLAY: Cannot access image "//l_image, \
            l_logfile, verbose+)
        goto crash
    } else if (gimverify.status != 0) {
        printlog ("ERROR - GDISPLAY: Image "//l_image//" is not a MEF file", \
            l_logfile, verbose+)
        goto crash
    }
    l_image = gimverify.outname # make sure .fits is gone

    printlog ("GDISPLAY: Using image "//l_image//"\n", l_logfile, l_verbose)

    # Set default values used for Xoff Yoff, l_<x/y>center calculations
    # and gapvalue
    imgets (l_image//"[0]", "DETTYPE", >& "dev$null")
    if (imgets.value == "S10892" ) { # Hamamatsu CCDs SOUTH
        ishamamatsu = yes
        l_gapvalue = 61 ##M CHIP_GAP
        xcenval = 2048 #hcode
        Xlen = 6144 #hcode
        imX = 6400. #hcode
        imY = 4644. #hcode
        Ymax = 4224 #hcode # Changed to reflect new CCD size
        Ymaxtrim = 4176 #hcode # Set for trimmed Hamamatsu data
    } else if (imgets.value == "S10892-N") { #Hamamatsu NORTH
        ishamamatsu = yes
        l_gapvalue = 67 ##M CHIP_GAP
        xcenval = 2048 #hcode
        Xlen = 6144 #hcode
        imX = 6400. #hcode
        imY = 4644. #hcode
        Ymax = 4224 #hcode # Changed to reflect new CCD size
        Ymaxtrim = 4176 #hcode # Set for trimmed Hamamatsu data
    } else if (imgets.value == "SDSU II CCD") { # Current EEV2 CCDs
        #Default values for GMOS S and old GMOS N
        l_gapvalue = 37
        l_stand_sat = 58000. # Standard saturation limit
        l_reduced_sat = 45000. # Twilight reduced saturation limit
        xcenval = 2048
        Xlen = 6144
        imX = 6400.
        imY = 4644.
        Ymax = 4608
    } else if (imgets.value == "SDSU II e2v DD CCD42-90") { \
        # New e2vDD CCDs
        ise2vDD = yes
        l_gapvalue = 37 ##M
        xcenval = 2048
        Xlen = 6144
        imX = 6400.
        imY = 4644.
        Ymax = 4608
        l_nbiascontam = 5 ##M NBIASCONTAM
    } else {
        keypar (l_image//"[0]", "INSTRUME", silent+)
        if (keypar.found) {
            if (strstr("GMOS",str(keypar.value)) == 0) {
                printlog ("ERROR - GDISPLAY: Image is not a GMOS image.", \
                    l_logfile, verbose+)
                goto crash
            } else {
                printlog ("ERROR - GDISPLAY: Unregognized DETTYPE value in "//\
                    l_image//"[0]", l_logfile, verbose+)
                goto crash
            }
        } else {
            printlog ("ERROR - GDISPLAY: INSTRUME keyword missing from "//\
                l_image//"[0]", l_logfile, verbose+)
            goto crash
        }
    }

    # Set gapvalue
    if (l_gap == "default") {
        gapvalue = l_gapvalue
    } else {
        gapvalue = int(l_gap)
    }

    # Test printing
    if (debug) {
        printlog ("*** Xlen="//Xlen//" imX="//imX//" imY="//imY//\
            " Ymax="//Ymax//" gap="//gapvalue//" xcenval ="//xcenval, \
            l_logfile, verbose+)
    }

    # Read the number of ROIs
    ndetroi = 1
    keypar (l_image//"[0]", "DETNROI", silent+)
    if (keypar.found) {
        ndetroi = int (keypar.value)
    }

    # Print out image information if verbose
    if (l_verbose) {
        printlog ("GDISPLAY: Image information...", l_logfile, verbose+)
        keypar (l_image//"[0]", "OBJECT", silent+)
        if (keypar.found) {
            pr_str = keypar.value
        } else {
            pr_str = "INDEF"
        }

        printlog ("              OBJECT  = "//pr_str, \
            l_logfile, l_verbose)

        keypar (l_image//"[0]", "GRATING", silent+)
        if (keypar.found) {
            pr_str = keypar.value
        } else {
            pr_str = "INDEF"
        }

        printlog ("              GRATING = "//pr_str, \
            l_logfile, l_verbose)
        if (keypar.found) {
            pr_str = keypar.value
        } else {
            pr_str = "INDEF"
        }

        keypar (l_image//"[0]", "FILTER1", silent+)
        if (keypar.found) {
            pr_str = keypar.value
        } else {
            pr_str = "INDEF"
        }
        printlog ("              FILTER1 = "//pr_str, \
            l_logfile, l_verbose)

        keypar (l_image//"[0]", "FILTER2", silent+)
        if (keypar.found) {
            pr_str = keypar.value
        } else {
            pr_str = "INDEF"
        }
        printlog ("              FILTER2 = "//pr_str, \
            l_logfile, l_verbose)

        if (ndetroi > 1 && l_req_roi == 0) {
            printlog ("              DETNROI = "//ndetroi, \
                l_logfile, l_verbose)
        }
        printlog ("", l_logfile, l_verbose)
    }

    # Form the tmporaray output name
    fparse (gimverify.outname)
    l_imext = fparse.root//".fits"
    tmpout = tmpout//l_imext
    fparse (l_image)

    if (fparse.extension != ".fits") {
        l_image = l_image//".fits"
    }

    # Check stdimage is set properly
    if (defvar("stdimage")) {
        show ("stdimage") | scan (l_stdimage)
    } else {
        printlog ("ERROR - GDISPLAY: environmental variable stdimage is not \
            set", l_logfile, verbose+)
        goto crash
    }

    printlog ("GDISPLAY: stdimage = "//l_stdimage//"\n", l_logfile, l_verbose)

    if ((substr(l_stdimage,1,7) == "imtgmos") && \
        (l_stdimage != "imtgmosccd")) {
        Xmax = Xlen + (2. * gapvalue)
        Xoff = nint(imX - Xmax) / 2.
        Yoff = nint(imY - Ymax) / 2.
        Xmax = imX
        Ymax = imY
    } else {
        printlog ("ERROR - GDISPLAY: Please set stdimage to a full display \
            GMOS device", l_logfile, verbose+)
        printlog ("                  imtgmos, imtgmos2, imtgmos4 or imtgmos8",\
            l_logfile, verbose+)
        goto crash
    }

    l_mag = 1.
    if (strlen(l_stdimage) == 8) {
        l_mag = 1./real(substr(l_stdimage,8,8))
    }

    if (l_sci_ext != "ALL") {
        fxhead (l_image) | match ("IMAGE", "STDIN", stop-) | \
            translit ("STDIN", "a-z", "A-Z", delete-, collapse-) | \
            match (l_sci_ext, "STDIN", stop-, > tmpin)
    } else {
        n_k = 0
        fxhead (l_image) | match("IMAGE", "STDIN", stop-) | \
            translit ("STDIN", "a-z", "A-Z", delete-, collapse-) | \
            match ("SCI", "STDIN", stop-) | scan (n_k)
        if (n_k > 0) {
            fxhead (l_image) | match ("IMAGE", "STDIN", stop-) | \
                translit ("STDIN", "a-z", "A-Z", delete-, collapse-) | \
                match ("SCI", "STDIN", stop-, > tmpin)
            printlog ("GDISPLAY: Displaying the SCI extensions", \
                l_logfile, l_verbose)
            l_sci_ext = "SCI"
        } else {
            fxhead (l_image) | match ("IMAGE", "STDIN", stop-, > tmpin)
            # raw data
            l_sci_ext = ""
        }
    }

    count (tmpin) | scan (n_ext) # number of extensions

    delete (tmpin, verify-, >& "dev$null")

    # Set "default" values according to detetctor type
    # Saturation value only for EEV CCDs currently - rest are defined in gsat
    ##M This can be removed once gsat has all of the CORRECT values
    if (l_satval == "default" && (!ise2vDD && !ishamamatsu)) {
        print (l_stand_sat) | scan (l_satval)

        # Adjust if 1x1 binning
        if (l_sci_ext == "") {
            keypar (l_image//"[1]", l_key_ccdsum, silent+)
            if (!keypar.found) {
                keypar (l_image//"[2]", l_key_ccdsum, silent+)
            }
        } else {
            keypar (l_image//"["//l_sci_ext//",1]", l_key_ccdsum, silent+)
        }

        Xbin = 0
        Ybin = 0
        if (keypar.found) {
            print (keypar.value) | scanf ("%d %d", Ybin, Ybin)
            if (Xbin == 1 && Ybin == 1) {
                print (l_reduced_sat) | scan (l_satval)
            }
        }
    }

    # bias_type
    if (l_bias_type == "default") {
        if (ise2vDD || ishamamatsu) {
            l_bias_type = "calc"
        } else {
            l_bias_type = "static"
        }
    }

    # Check for unecessary flags
    if (n_ext == 1 && l_fl_paste) {
        printlog ("WARNING - GDISPLAY: only one science extension found, \
            turning paste off", l_logfile, verbose+)
        l_fl_paste = no
    }

    # Check to switch off flags if image is processed - they are in order - MS
    iscut = no
    ismosaiced = no
    istiled = no

    keypar (l_image//"[0]", "GSCUT", silent+)
    if (keypar.found) {
        iscut = yes
        pr_str = "GSCUT"
        goto NEXT_CHK
    }

    keypar (l_image//"[0]", "GMOSAIC", silent+)
    if (keypar.found) {
        ismosaiced = yes
        pr_str = "GMOSAIC"
        goto NEXT_CHK
    }

    keypar (l_image//"[0]", "GTILE", silent+)
    if (keypar.found) {
        istiled = yes
        pr_str = "GTILE"
        goto NEXT_CHK
    }

NEXT_CHK:

    # Switch certain things off if required (reduces the number of statements
    # printed to screen from GMULTIAMP)
    if (iscut || ismosaiced || istiled) {

        if (l_fl_bias || l_fl_pretty || l_fl_paste || l_ret_roi) {

            printlog ("\nWARNING - GDISPLAY: Input image has been processed"//\
                " by "//pr_str//"."//\
                "\n                    Switching off the following flags:", \
                l_logfile, verbose+)

            if (l_fl_bias) {
                l_fl_bias = no
                printlog ("                          fl_bias=no", \
                    l_logfile, verbose+)
            }
            if (l_fl_pretty) {
                l_fl_pretty = no
                printlog ("                          fl_pretty=no", \
                    l_logfile, verbose+)
            }

            # Can paste gtiled files - it will handle the state appropriately
            if (l_fl_paste && (iscut || ismosaiced)) {
                l_fl_paste = no
                printlog ("                          fl_paste=no", \
                    l_logfile, verbose+)
            }
            if (l_ret_roi) {
                l_ret_roi = no
                printlog ("                          ret_roi=no", \
                    l_logfile, verbose+)
            }
            printlog ("", l_logfile, verbose+)
        }
    }

    # Set the sciextn parameter
    if (l_sci_ext != "") {
        sciextn = l_sci_ext//","
    } else {
        sciextn = ""
    }

    # Now act on the following items
    if (l_fl_paste || l_fl_bias || l_fl_pretty || !isindef(signal) \
        || l_fl_sat || l_ret_roi) {

        printlog ("GDISPLAY: Temporary image name "//tmpout, \
            l_logfile, verbose=yes)

        # Call gmultiamp this will do everything else (bar display)
        gmultiamp (l_image, outimages=tmpout, \
            outprefix="", ret_roi=l_ret_roi, req_roi=l_req_roi, \
            fl_goversub=l_fl_bias, rawpath="", sci_ext=l_sci_ext, \
            var_ext=l_var_ext, dq_ext=l_sig_ext, mdf_ext=l_mdf_ext, \
            key_biassec=l_key_biassec, bias_type=l_bias_type, pixstat="midpt",\
            nbiascontam=l_nbiascontam, statsec="default", \
            fl_gtile=l_fl_paste, out_ccds="all", fl_stats_only=no, \
            fl_tile_det=l_fl_paste, fl_pad=yes, key_detsec=l_key_detsec, \
            key_ccdsec=l_key_ccdsec, key_datasec=l_key_datasec, \
            key_ccdsum=l_key_ccdsum, fl_pretty=l_fl_pretty, \
            fl_sat=l_fl_sat, satval=l_satval, signal=l_signal, \
            satprefix="gmasat", sigprefix="gmasig", logfile=l_logfile, \
            gaindb=l_gaindb, verbose=yes)

        if (gmultiamp.status >= 2) {
            printlog ("ERROR - GDISPLAY: GMULTIAMP returned an ERROR status.",\
                l_logfile, verbose+)
            goto crash
        } else if (gmultiamp.status == 1) {
            # Nothing was done to the image
            tmpout = l_image
        } else {
            processed = yes
        }

        # Reset the number of extensions
        keypar (tmpout//"[0]", "NSCIEXT", silent+)
        if (keypar.found) {
            n_ext = int(keypar.value)
        } else {
            keypar (tmpout//"[0]", "NEXTEND", silent+)
            if (keypar.found) {
                n_ext = int(keypar.value)
            }
        }

        # Set saturation and signal mask names automatically
        sat_image = "gmasat"//tmpout
        sig_image = "gmasig"//tmpout
        l_image = tmpout

        # Check masks made...

        if (l_fl_sat && !processed) {
            printlog ("WARNING - GDISPLAY: Not displaying saturation masks", \
                l_logfile, verbose+)
            l_fl_sat = no
        } else if (!imaccess(sat_image) && l_fl_sat) {
            printlog ("WARNING - GDISPLAY: Not displaying saturation masks", \
                l_logfile, verbose+)
            l_fl_sat = no
        }

        if (!isindef(l_signal) && !processed) {
            printlog ("WARNING - GDISPLAY: Not displaying signal masks", \
                l_logfile, verbose+)
            l_signal = INDEF
        } else if (!imaccess(sig_image) && !isindef(l_signal)) {
            printlog ("WARNING - GDISPLAY: Not displaying signal masks", \
                l_logfile, verbose+)
            l_signal = INDEF
        }

        # Reset sci_ext in certain circumstances
        # Set saturation and ignal extension variables
        if (l_sci_ext == "") {
            if (processed && (fl_pretty || ((l_fl_sat || !isindef(l_signal)) \
                && l_fl_bias))) {
                sciextn = "SCI"
            }
            satextn = "SCI"
            signalextn = "SCI"
        } else {
            satextn = l_sci_ext
            signalextn = l_sci_ext
            sciextn = l_sci_ext
        }

        if (sciextn != "") {
            sciextn = sciextn//","
        }

        if (satextn != "") {
            satextn = satextn//","
        }

        if (signalextn != "") {
            signalextn = signalextn//","
        }

        # Need to create a display mask if ndetroi > 1; can be reset by
        # gmultiamp depending on requested ROI
        keypar (tmpout//"[0]", "DETNROI", silent+)
        ndetroi = int(keypar.value)
    } # End of GMULTIAMP loop

    # Double check if it's now been gtiled - MS
    istiled = no
    keypar (l_image//"[0]", "GTILE", silent+)
    if (keypar.found) {
        istiled = yes
    }

    # Only create a display mask (for z values) if needed
    if (ndetroi > 1 && istiled && !iscut && fl_zscale && fl_zrange) {
        display_mask = mktemp ("tmpdispmask")//".fits"

        gem_sci_ext = l_sci_ext
        if (gem_sci_ext == "") {
            gem_sci_ext = "SCI"
        }
        # Create a mask for display command to tell it which pixels to
        # determine z1 z2 values from; 1 means use those pixels - MS
        gemexpr ("a == 0 ? 0 : 1", display_mask, \
            l_image, var_expr="", dq_expr="", sci_ext=gem_sci_ext, \
            var_ext=l_var_ext, dq_ext=l_dq_ext, mdf_ext=l_mdf_ext, \
            fl_vardq=no, dims="default", intype="default", \
            outtype="ushort", refim="default", rangecheck=yes, \
            logfile=l_logfile, verbose=no, exprdb="none")

        if (gemexpr.status != 0) {
            display_mask = ""
        }

        # Always check masks in this case
        fl_check_masks = yes

        # update the nsample parameter
        l_nsample = int (l_max_nsample / ndetroi)

    }

    # Loop over extnsions and display them
    for (counter = 1; counter <= n_ext; counter += 1) {

        # Set all of the extension parameters
        sciextn = l_sci_ext
        if (sciextn != "") {
            sciextn = sciextn//","
            msk_extn = sciextn
        } else {
            msk_extn = "SCI,"
        }

        # Set the current extesnions display mask
        disp_msk_extn = ""
        if (display_mask != "") {
            disp_msk_extn = display_mask//"["//msk_extn//counter//"]"
        }

        # Set input names
        extn = "["//sciextn//counter//"]"
        inextn = l_image//extn
        fparse (l_image)
        extn_root = fparse.root//fparse.extension//extn

        # Reset saturation and signal mask displaying options
        if (l_fl_sat) {
            dbpm = sat_image//"["//satextn//counter//"]"

            if (dbpm != "" && not_said) {
                printlog ("GDISPLAY: Saturated pixels flagged as red", \
                    l_logfile, verbose+)
                not_said = no
            }

        } else {
            dbpm = ""
        }

        if (!isindef(l_signal)) {
            doverlay =  sig_image//"["//signalextn//counter//"]"

            if (counter == 1) {
                printlog ("GDISPLAY: Pixels with signal >= "//str(l_signal)//\
                    " above bias flagged as green", l_logfile, verbose+)
            }
        } else {
            doverlay = ""
        }

        # Check the masks will not crash display
        if (fl_check_masks && (dbpm != "" || disp_msk_extn != "")) {

            chk_value = 1
            if (dbpm != "" && disp_msk_extn != "") {
                chk_mask = mktemp ("tmpchk_mask")//".fits"
                imexpr ("a == 0 ? 1 : b == 1 ? 1 : 0", chk_mask, disp_msk_extn,
                    dbpm, dims="auto", intype="auto", outtype="ushort", \
                    refim="auto", bwidth=0, btype="nearest", bpixval=0., \
                    rangecheck=yes, verbose=no, exprdb="none")

            } else if (dbpm != "") {
                chk_mask = dbpm
            } else {
                chk_mask = disp_msk_extn
                chk_value = 0 # If all zero there are no pixels to read
            }

            l_min = INDEF
            l_max = INDEF

            if (debug) {
                print (chk_mask)
                imstat (chk_mask, fields="min,max", lower=INDEF, upper=INDEF, \
                    nclip=0, lsigma=3., usigma=3., binwidth=0.1, format=no, \
                    cache=no)
            }

            imstat (chk_mask, fields="min,max", lower=INDEF, upper=INDEF, \
                nclip=0, lsigma=3., usigma=3., binwidth=0.1, format=no, \
                cache=no) | scan (l_min, l_max)

            if (isindef(l_min) || isindef(l_max)) {
                printlog ("WARNING - GDISPLAY: Imstat could not determine \
                    the min and max values for the display mask "//\
                    "\n                    and / or saturation mask"//\
                    "\n                    Switching both off", \
                    l_logfile, verbose+)

                dbpm = ""
                disp_msk_extn = ""

            } else if (l_min == chk_value && (l_min == l_max)) {
                # The extension is completely saturated
                if (dbpm != "") {
                    printlog ("WARNING - GDISPLAY: Extension "//inextn//\
                        " is COMPLETELY saturated. "//\
                        "\n                    Cannot display overlay of \
                        saturated pixels.", l_logfile, verbose+)
                }
                dbpm = ""
                disp_msk_extn = ""
            }

            if (dbpm != "" && disp_msk_extn != "") {
                imdelete (chk_mask, verify-, >& "dev$null")
            }
        }

        # Read DATASEC
        keypar (inextn, l_key_datasec, silent+)
        if (keypar.found) {
            l_datasec = str(keypar.value)

        } else {
            l_datasec = ""
        }

        if (counter == 1) {
            printlog ("GDISPLAY: Displaying image...", l_logfile, verbose+)
            # Set length of printf string used later on
            pr_len = strlen(extn_root) + strlen(str(n_ext)) - 1
        }

        # Read all of the informtion again if not pasted and more than one
        # extension
        if (!l_fl_paste && n_ext > 1) {
            # Read CCDSUM
            keypar (inextn, l_key_ccdsum, silent+)
            if (keypar.found) {
                l_ccdsum = str(keypar.value)
                print (l_ccdsum) | scanf ("%d %d", Xbin, Ybin)
            } else {
                # Assume 1 by 1 binning
                Xbin = 1
                Ybin = 1
            }

            # Read DETSEC
            keypar (inextn, l_key_detsec, silent+)
            if (keypar.found) {
                l_dsec = str(keypar.value)
                print (l_dsec) | scanf ("[%d:%d,%d:%d]", x1, x2, y1, y2)
            } else {
                x1 = 1
                y1 = 1

                keypar (inextn, "i_naxis1", silent+)
                if (!keypar.found) {
                    printlog ("ERROR - GDSIPLAY: Cannot find NAXIS1 keyword \
                        in "//inextn, l_logfile, verbose+)
                    goto crash
                }
                x2 = int(keypar.value) * Xbin

                keypar (inextn, "i_naxis2", silent+)
                if (!keypar.found) {
                    printlog ("ERROR - GDSIPLAY: Cannot find NAXIS2 keyword \
                        in "//inextn, l_logfile, verbose+)
                    goto crash
                }
                y2 = int(keypar.value) * Ybin

            }

            # Set up display options for this extension
            l_xcenter = ((x2 - x1) / 2.) + x1
            if (l_xcenter > xcenval) {
                l_xcenter += gapvalue
            }

            if ( l_xcenter > ((2 * xcenval) + gapvalue)) {
                l_xcenter += gapvalue
            }

            l_xcenter += Xoff

            l_xcenter = l_xcenter / Xmax

            l_ycenter = (((y2 - y1)/ 2.) + y1 + Yoff) / Ymax

            # Reset xcenter,ycenter if image is binned
            l_xcenter = ((l_xcenter - 0.5) / Xbin) + 0.5
            l_ycenter = ((l_ycenter - 0.5) / Ybin) + 0.5

            # Reset magnification if higher resolution can be displayed
            # Never display more than 1pix to 1pix
            if ((min(Xbin,Ybin) > 1.) && (l_mag*min(Xbin,Ybin) <= 1.) && \
                l_fl_chmag) {

                l_mag = l_mag * min(Xbin,Ybin)
                l_fl_chmag = no  # don't check/change l_mag in the next rounds
                l_fl_mag = yes   # update xcenter,ycenter in all rounds
            }

            # Adjust xcenter,ycenter if l_mag has been changed
            if (l_fl_mag) {
                l_xcenter = (l_xcenter-0.5) * min(Xbin,Ybin) + 0.5
                l_ycenter = (l_ycenter-0.5) * min(Xbin,Ybin) + 0.5
            }

            printf ("    %-"//pr_len//"s %4d %4d %4d %4d %5.1f ",
                    extn_root, x1, x2, y1, y2, l_mag, >& l_redirect)

        } else {
            # Either pasted or only one input extension
            # If displaying a previously mosacied image (n_ext==1) reset
            # (x/y)center values to 0.5, as in the l_fl_paste is true (gmosaic)
            # case. Any change in magnification should not affect the centering
            # - MS
            l_xcenter = 0.5
            l_ycenter = 0.5
            l_mag = 1

            printf ("    %s ", inextn, >& l_redirect)

        } # End of setting up display commands for this extension

        # Set th edisplay frame erase flag
        if (counter == 1) {
            fl_erase = yes
        } else {
            fl_erase = no
        }

        # Display the image
        display (inextn//l_datasec, l_frame, \
            erase=fl_erase, xcenter=l_xcenter, ycenter=l_ycenter, z1=l_z1, \
            z2=l_z2, zrange=fl_zrange, zscale=fl_zscale, zmask=disp_msk_extn,\
            nsample=l_nsample, xmag=l_mag, ymag=l_mag, fill=no, bpm=dbpm, \
            bpdisplay="overlay", bpcolors="red", overlay=doverlay, \
            ocolors="green", >& l_redirect)
    }

    # Imexam
    if (n_ext == 1 && l_fl_imexam) {
        print ("\nGDISPLAY: Starting imexamine - quit with q")
        imexam (inextn//l_datasec, frame=l_frame, display="")
    }

    # Copy the image if requested
    if (l_output != "") {

        if (imaccess(l_output)) {
            printlog ("ERROR - GDISPLAY: output "//l_output//\
                " already exists", l_logfile, verbose+)
            goto crash
        }

        copy (l_image, l_output, verbose-)
        # Append any masks
        if (l_fl_sat && app_masks) {
            for (counter = 1; counter <= n_ext; counter += 1) {
                imcopy (sat_image//"["//satextn//counter//"]", \
                    l_output//"["//l_sat_ext//","//counter//",append+]", \
                    verbose-)

                gemhedit (l_output//"["//l_sat_ext//","//counter//"]", \
                    "i_title", "saturation mask", "", delete-)
            }
        }
        if (!isindef(l_signal) && app_masks) {
            for (counter = 1; counter <= n_ext; counter += 1) {
                imcopy (sig_image//"["//satextn//counter//"]", \
                    l_output//"["//l_sig_ext//","//counter//",append+]", \
                    verbose-)

                gemhedit (l_output//"["//l_sig_ext//","//counter//"]", \
                    "i_title", "signal mask", "", delete-)
            }
        }
    }

    goto clean

crash:
    # Exit with error subroutine
    status = 1

clean:

    if (processed) {
        imdelete (l_image, verify-, >& "dev$null")
    }
    if (l_fl_sat) {
        imdelete (sat_image, verify-, >& "dev$null")
    }
    if (!isindef(l_signal)) {
        imdelete (sig_image, verify-, >& "dev$null")
    }

    if (imaccess(display_mask)) {
        imdelete (display_mask, verify-, >& "dev$null")
    }

    # Print start time
    gemdate (zone="local")
    printlog ("\nGDISPLAY - Finished: "//gemdate.outdate//"\n", \
        l_logfile, verbose+)

    if (status != 0) {
        printlog ("GDISPLAY -- exit staus: ERROR", l_logfile, verbose+)
    }

    if (profile) {
        date ("+%H:%M:%S.%N%n") | scan (sdate)
        printlog ("____GDISPLAY Finished: "//sdate, l_logfile, verbose+)
    }

end

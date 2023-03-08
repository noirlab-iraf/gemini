# Copyright(c) 2011-2015 Association of Universities for Research in Astronomy, Inc.

procedure gmultiamp (inimages)

# This is a wrapper script for goversub and gtile, for use only with gemlocal
# tasks like gacq, ghartmann, ggrating etc. This script assumes that all input
# files have the same number of amplifiers per CCD; that is the keyword, NAMPS
# is the same in each input file's PHU.

# goversub is called by default when NAMPS >= 2.
# gtile is called by default

# Updated to create saturation masks and signal masks (formily in gdisplay)

# Updated with fl_pretty option for operations

# Standard user input parameters
char    inimages     {prompt="Input image to be roughly overscan subtracted and or gtiled."}
char    outimages    {"", prompt="Output filenames"}
char    outprefix    {"gma", prompt="Prefix for output images"}

# General ROI parameter
bool    ret_roi      {yes, prompt="Return only the data within specified ROI"}
int     req_roi      {0, prompt="Requested ROI (0=all)"}

# goversub parameters
char    fl_goversub  {"default", enum="default|yes|no", prompt="Call goversub? (default|yes|no)"}
char    rawpath      {"", prompt="Path to raw data set"}
char    sci_ext      {"SCI", prompt="Name of science extension(s)"}
char    var_ext      {"VAR", prompt="Name of variance extension(s)"}
char    dq_ext       {"DQ", prompt="Name of data quality extension(s)"}
char    mdf_ext      {"MDF", prompt="Mask definition file extension name"}
char    key_biassec  {"BIASSEC", prompt="Header keyword for overscan section"}
char    bias_type    {"calc", enum="calc|static", prompt="Bias value type to determine."}
char    pixstat      {"midpt", prompt="Statistic to be calculated by imstatistics to use to determine the rough BIAS level"}
int     nbiascontam  {5, prompt="Number of columns to strip from BIASSEC on the side closest to the DATASEC",min=0}
char    statsec      {"default", prompt="Section of BIASSEC to use. Coordinates are relative to the BIAS section. In the format: [%d:%d,%d:%d]"}

# gtile parameters
bool    fl_gtile     {yes, prompt="Call gtile?"}
char    out_ccds     {"all", enum="all|1|2|3", prompt="Which CCDs to return. (all|1|2|3)"}
bool    fl_stats_only{no, prompt="Return only the science data for use with statistics tasks"}
bool    fl_tile_det  {no, prompt="To tile all three CCDs"}
bool    fl_pad       {no, prompt="Include chip gaps in tiled detector image?"}
bool    fl_app_rois  {no, prompt="Tile and append ROIs separately?"}
char    key_detsec   {"DETSEC", prompt="Header keyword for detector section"}
char    key_ccdsec   {"CCDSEC", prompt="Header keyword for CCD section"}
char    key_datasec  {"DATASEC", prompt="Header keyword for data section"}
char    key_ccdsum   {"CCDSUM", prompt="Header keyword for CCD binning"}

# Pretty / mask parameters
bool    fl_pretty    {no, prompt="For aquistion / gemlocal only"}
bool    fl_sat       {no, prompt="Create mask containing saturated pixles?"}
char    satval       {"default", prompt="saturation value to use"}
real    signal       {INDEF, prompt="Create mask of pixels with a value above signal?"}
char    satprefix    {"default", prompt="Prefix for saturation mask image"}
char    sigprefix    {"default", prompt="Prefix for saturation mask image"}

# joint parameters
char    gaindb       {"default", prompt="Database with gain data"}
char    logfile      {"dev$null", prompt="Logfile for this task"}
bool    verbose      {no, prompt="Verbose, yes or no"}
int     status       {0, prompt="Exit status, 2=ERROR; 0=GOOD, performed as expected; 1=good, no processing done"}
struct  *scanfile    {"", prompt="For internal use only"}
struct  *scanfile2   {"", prompt="For internal use only"}
struct  *scanfile3   {"", prompt="For internal use only"}
struct  *scanfile4   {"", prompt="For internal use only"}

#####
begin

    # Declare local variables used within script not set from user input
    char    tmpinlist, tmpsub, file_to_check, tile_inname, l_crash_list
    char    inname, tmpoutlist, outname, testfile, rphend, pathtest
    char    go_in, go_out, go_path, p_in, p_out, p_path, proc_task
    char    tile_in, tile_out, tile_path, tile_list, edit_image
    char    outimg, nextin, nextout, nextpath, l_gkey, sig_out
    char    tmpoutlist2, debias_info, tmpdb, tile_imtodel, tmpunique
    char    l_inum_str, l_fnum_str, l_char_str, l_key_oscan
    char    sat_key, sig_key, outextn, sig_ext, sat_ext, tmpdet_info
    char    sat_out, tmp_sat_out, l_sat_out, tmp_sig_out, l_sig_out
    char    grroi_in, grroi_out, grroi_path, tmproi, orig_inimage
    char    gexpr_sexpr, gexpr_vexpr, gexpr_dexpr, all_oscan_vals
    char    tmposcan, extn, sciextn, go_del, p_del, l_out_ext
    char    l_gain_key, l_key_gain, l_key_rdnoise, img_root, l_gem_outtype
    char    all_gain_vals, tmpgains, l_gain, l_detsec, l_ampname
    char    l_extn, gain, oscan, l_key_ampname, l_satdb, gexpr_svexpr
    char    sigeqn, inimg, sateqn, tmpheader, infoimg, tmpamps, tmppretty
    char    l_rhsval, l_rhscol, l_lhsval, l_lhscol, l_nwidth

    real    l_repval, l_oscan, l_signal, sig_value, sat_value
    real    sub_value, gtotal, lhs_value, rhs_value

    int     num_files, junk, atlocation, test_fn, dx1, dx2
    int     ii, jj, xbin, num_unique_amps, l_naxis1, l_naxis2, test, tx1, tx2
    int     nlines, ncols_lhs, ncols_rhs, col_limit
    int     oscan_cnt, nsci, gain_cnt, namps, ccdwidth

    bool    gexpr_vardq, go_calc_only, oscan_sub, cangoversub
    bool    fl_debias, isgainmult, call_gretroi, call_gemexpr, isoverscansub
    bool    isprocessed, use_prefix, read_debias, call_2_gretroi
    bool    debug, profile

    struct sdate

    # Declare local variables for user input
    char l_inimages
    char l_outimages
    char l_outprefix
    char l_rawpath
    bool l_ret_roi
    int  l_req_roi
    char l_bias_type
    char l_pixstat
    char l_sci_ext
    char l_var_ext
    char l_dq_ext
    char l_mdf_ext
    char l_fl_goversub
    int  l_nbiascontam
    char l_key_biassec
    char l_statsec
    bool l_fl_gtile
    bool l_fl_app_rois
    char l_key_detsec
    char l_key_ccdsec
    char l_key_datasec
    char l_key_ccdsum
    char l_out_ccds
    char l_sigprefix
    char l_satprefix
    bool l_fl_stats_only
    bool l_fl_tile_det
    bool l_fl_pretty
    char l_logfile
    bool l_verbose
    bool l_fl_sat
    bool l_fl_signal
    char l_satval
    char l_gaindb
    bool l_fl_pad

    # Read user input
    l_inimages = inimages
    l_outprefix = outprefix
    l_outimages = outimages
    l_rawpath = rawpath
    l_ret_roi = ret_roi
    l_req_roi = req_roi
    l_fl_goversub = fl_goversub
    l_bias_type = bias_type
    l_pixstat = pixstat
    l_sci_ext = sci_ext
    l_var_ext = var_ext
    l_dq_ext = dq_ext
    l_mdf_ext = mdf_ext
    l_nbiascontam = nbiascontam
    l_statsec = statsec
    l_key_biassec = key_biassec
    l_key_detsec = key_detsec
    l_key_ccdsec = key_ccdsec
    l_key_datasec = key_datasec
    l_key_ccdsum = key_ccdsum
    l_out_ccds = out_ccds
    l_fl_gtile = fl_gtile
    l_fl_app_rois = fl_app_rois
    l_fl_stats_only = fl_stats_only
    l_fl_tile_det = fl_tile_det
    l_sigprefix = sigprefix
    l_satprefix = satprefix
    l_fl_pretty = fl_pretty
    l_satval = satval
    l_gaindb = gaindb
    l_fl_sat = fl_sat
    l_fl_pad = fl_pad
    l_signal = signal
    l_logfile = logfile
    l_verbose = verbose

    # Printf strings for gemhedit upfile creation
    l_inum_str = "%-8s %d \"%s\"\n"
    l_fnum_str = "%-8s %.5f \"%s\"\n"
    l_char_str = "%-8s \"%s\" \"%s\"\n"

    # Key words names either used as inputs or updates for use with gemexpr
    l_key_ampname = "AMPNAME"
    l_key_rdnoise = "RDNOISE"
    l_key_gain = "GAIN"
    l_key_oscan = "OVERSUB"
    sat_key = "SATURATI"
    sig_key = "SIGNAL"
    l_nwidth = "IMWIDTH"
    l_lhsval = "LHSVAL"
    l_lhscol = "LHSNCOL"
    l_rhsval = "RHSVAL"
    l_rhscol = "RHSNCOL"

    # Create temporary variable names / files
    l_crash_list = mktemp("tmplcrashlist")
    tmpinlist = mktemp ("tmpinlist")
    tmpoutlist = mktemp ("tmpoutlist")
    tmpoutlist2 = mktemp ("tmpoutlistII")
    tmproi = mktemp("tmproi")
    tmpsub = mktemp("tmpsub")
    tmppretty = mktemp("tmppretty")
    tmpheader = mktemp ("tmpheader")
    tmp_sat_out = mktemp("tmpsat_out")//".fits"
    tmp_sig_out = mktemp("tmpsig_out")//".fits"
    tmpunique = mktemp ("tmpunique")
    tmposcan = mktemp ("tmposcan")
    tmpgains = mktemp ("tmpgains")
    tmpamps = mktemp ("tmpamps")
    tile_list = mktemp ("tmptile_list")
    tmpdet_info = mktemp ("tmpdetinfo")

    # Database / lookup files
    debias_info = "gemlocal$data/gmultiamp_info.dat"
    # Open up this to user?
    l_satdb = "gmos$data/gmosFWdepths.dat"

    # Set default parameters
    status = 0
    go_del = ""
    p_del = ""
    ccdwidth = 2048
    debug = no
    profile = no

    if (debug) {
        l_verbose = yes
        profile = yes
    }

    if (!isindef(l_signal)) {
        l_fl_signal = yes
    } else {
        l_fl_signal = no
    }

    # Cache the parameters for tasks used in this script
    cache ("goversub", "gtile", "gemarith", "ggain", "gretroi")

    if (profile) {
        date ("+%H:%M:%S.%N%n") | scan (sdate)
        printlog ("____GMULTIAMP Started: "//sdate, l_logfile, verbose+)
    }

    # Print start time
    gemdate (zone="local")
    printlog ("\nGMULTIAMP - Started: "//gemdate.outdate//"\n", \
        l_logfile, l_verbose)

    ######## Start checks ########

    # Only allow fl_pretty if inside gemini
    if (!deftask("gemlocal") && l_fl_pretty) {
        printlog ("WARNING - GMULTIAMP: fl_pretty is for observation \
            purposes only. Switching fl_pretty off", l_logfile, verbose+)
        l_fl_pretty = no
    }

    # Set the debias
    fl_debias = no
    if (l_fl_pretty) {
        if (access(debias_info)) {
            fl_debias = yes
        } else {
            printlog ("WARNING - GMULTIAMP: Cannot access gmultiamp data \
                file."//\
                "\n          Additional bias levels next to amplifier \
                boarders will not be removed.", l_logfile, verbose+)
        }
    }

    if (l_fl_sat) {
        if ((l_satprefix == "") || (stridx(l_satprefix," ") > 0)) {
            printlog ("ERROR - GMULTIAMP: fl_sat=yes but satprefix is not set \
                correctly.", l_logfile, verbose+)
            status = 2
            goto crash
        } else if (l_satprefix == "default") {
            l_satprefix = l_outprefix//"sat"
        }
    }

    if (l_fl_signal) {
        if ((l_sigprefix == "") || (stridx(l_sigprefix," ") > 0)){
            printlog ("ERROR - GMULTIAMP: fl_sat=yes but satprefix is not set \
                correctly.", l_logfile, verbose+)
            status = 2
            goto crash
        } else if (l_sigprefix == "default") {
            l_sigprefix = l_outprefix//"sig"
        }
    }

    ##M This is may be a slight hack right now - but it will require removing
    ##M logic from GTILE and putting it in GRETROI, which should probably be
    ##M done but it will take too long right now - MS (2012-11-05)
    if (l_out_ccds != "all") {
        if (!l_fl_gtile) {
            printlog ("WARING - GMULTIAMP: out_ccds!=\"all\" setting \
                fl_gtile=yes", \
                l_logfile, verbose+)
            l_fl_gtile = yes
        }
    }

    #### Rawpath ####
    if ((l_rawpath != "") && (stridx(" ", l_rawpath) == 0)) {
        rphend = substr(l_rawpath,strlen(l_rawpath),strlen(l_rawpath))
        if (rphend == "$") {
            show (substr(l_rawpath,1,strlen(l_rawpath)-1)) | scan (pathtest)
            rphend = substr(pathtest,strlen(pathtest),strlen(pathtest))
        }
        if (rphend != "/") {
            l_rawpath = l_rawpath//"/"
        }
        if (!access(l_rawpath)) {
            printlog ("ERROR - GMULTIAMP: Cannot access rawpath: "//\
                l_rawpath, l_logfile, verbose+)
            goto crash
        }
    }

    # Check that the input list exits (if it's a list)
    atlocation = stridx("@",l_inimages)
    if (atlocation > 0) {
        testfile = substr(l_inimages,atlocation+1,strlen(l_inimages))
        if (!access(testfile)){
            printlog ("ERROR - GMULTIAMP: Cannot access file "//testfile//\
                ". Exiting", l_logfile, verbose=yes)
            status = 2
            goto crash
        }
    }

    # Sort out input images
    files (l_inimages, sort-, > tmpinlist)
    print (tmpinlist, > l_crash_list)

    # Set the output name for gtile if outnames is ""; also create the tmp file
    # containing the output names to check to see if the exist.
    use_prefix = no
    if (l_outimages == "") {
        use_prefix = yes
        scanfile = tmpinlist
        while (fscan(scanfile,inname) != EOF) {
            fparse (inname)
            outname = l_outprefix//str(fparse.root)//str(fparse.extension)
            print (outname, >> tmpoutlist)
        }
        scanfile = ""
    } else {
        files (l_outimages, sort-, > tmpoutlist)
    }

    # Re-assign the l_outimages variable to tmpoutlist2 so the gimverify file
    # outnames are written correctly
    l_outimages = tmpoutlist2

    # Check the output files do not already exist
    scanfile = tmpoutlist
    while (fscan(scanfile,outname) != EOF) {
        gimverify(outname)
        if (gimverify.status != 1) {
            printlog ("ERROR - GMULTIAMP: "//outname//" already exists. \
                Exiting.", l_logfile, verbose=yes)
            status = 2
            goto crash
        } else {
            print (gimverify.outname//".fits", >> l_outimages)
        }
    }
    l_outimages = "@"//l_outimages

    file_to_check = ""
    scanfile = tmpinlist
    scanfile2 = tmpoutlist

    if (debug) {
        date ("+%H:%M:%S.%N%n") | scan (sdate)
        printlog ("____Starting loop over images: "//sdate, \
            l_logfile, verbose+)
    }

    # Loop over images and reduce them separately
    while (fscan(scanfile,inimg) != EOF) {

        # Read the output name (even if using prefix)
        test_fn = fscan (scanfile2,outimg)
        if (test_fn == 0 || test_fn == EOF) {
            printlog ("ERROR - GMULTIAMP: Incorrect number of output images", \
                l_logfile, verbose+)
            status = 2
            goto crash
        }

        # Check input exists
        file_to_check = l_rawpath//inimg

        orig_inimage = file_to_check
        if (!imaccess(file_to_check)) {
            printlog ("ERROR - GMULTIAMP: Cannot access "//rawpath//inimg//". \
                Exiting.", l_logfile, verbose+)
            status = 2
            goto crash
        }

        fparse (outimg)
        if (use_prefix) {
            l_sat_out = l_satprefix//substr(fparse.root,strlen(l_outprefix)+1,\
                strlen(fparse.root))//".fits"
            l_sig_out = l_sigprefix//substr(fparse.root,strlen(l_outprefix)+1,\
                strlen(fparse.root))//".fits"
        } else {
            l_sat_out = l_satprefix//fparse.root//".fits"
            l_sig_out = l_sigprefix//fparse.root//".fits"
        }

        if (l_fl_sat && imaccess(l_sat_out)) {
            printlog ("ERROR - GMULTIAMP: Saturation mask file "//l_sat_out//\
                " already exists", l_logfile, verbose+)
            status = 2
            goto crash
        }

        if (l_fl_signal && imaccess(l_sig_out)) {
            printlog ("ERROR - GMULTIAMP: Signal mask file "//l_sig_out//\
                " already exists", l_logfile, verbose+)
            status = 2
            goto crash
        }


        # Check the requested operations against the reduction state of the
        # image

        isprocessed = no
        if ((l_fl_sat && l_satval == "default") || l_fl_pretty || \
            l_ret_roi || l_fl_goversub != "no") {

            keypar (file_to_check//"[0]", "GMOSAIC", silent+)

            if (keypar.found) {
                isprocessed = yes
                proc_task = "GMOSAIC"
                goto SET_FLAGS
            }

            keypar (file_to_check//"[0]", "GTILE", silent+)
            if (keypar.found) {
                isprocessed = yes
                proc_task = "GTILE"
                goto SET_FLAGS
            }

            keypar (file_to_check//"[0]", "GSCUT", silent+)
            if (keypar.found) {
                isprocessed = yes
                proc_task = "GSCUT"
                goto SET_FLAGS
            }

            if (!isprocessed) {
                goto FLAGS_OK
            }

SET_FLAGS:
            printlog ("WARNING - GMULTIAMP: "//inimg//" has already been \
                prcossed with "//proc_task//\
                "\n                     Switching the following flags off:", \
                l_logfile, verbose+)

            if (l_fl_sat && l_satval == "default") {
                l_fl_sat = no
                printlog ("                     fl_sat=no", \
                    l_logfile, verbose+)
            }

            if (l_fl_goversub != "no") {
                l_fl_goversub = "no"
                printlog ("                     fl_goversub=no", \
                    l_logfile, verbose+)
            }

            if (l_fl_pretty) {
                l_fl_pretty = no
                printlog ("                     fl_pretty=no", \
                    l_logfile, verbose+)
            }

            if (l_ret_roi) {
                l_ret_roi = no
                printlog ("                     ret_roi=no", \
                    l_logfile, verbose+)
            }

FLAGS_OK:

        }

        # Set the goversub variable
        cangoversub = no
        isoverscansub = no # Used/reset later on

        if (l_fl_goversub == "default" || l_fl_pretty) {
            # If default check number of amplifiers
            keypar (file_to_check//"[0]", "NAMPS", silent+)
            if (!keypar.found) {
                printlog ("ERROR - GMULTIAMP: NAMPS keyword not found in "//\
                    inimg, l_logfile, verbose+)
                status = 2
                goto crash
            } else {
                namps = int(keypar.value)
            }
        }

        if (l_fl_goversub == "default") {

            # If number of amps per CCD is 2 or more or calling fl_pretty
            # set to yes
            if (namps > 1 || l_fl_pretty) {
                cangoversub = yes
            }

        } else if (l_fl_goversub == "yes") {
            cangoversub = yes
        }

        # Flag to perform the overscan subtraction with gemexpr
        oscan_sub = no
        if (cangoversub && (fl_pretty || l_fl_sat || l_fl_signal)) {
            oscan_sub = yes
        }

        # Set the flag to see if gemexpr will be called
        call_gemexpr = no

        if (l_fl_pretty || l_fl_sat || l_fl_signal) {
            call_gemexpr = yes

            # Reset canoversub to yes for it to be called - but call_gemexpr
            # flag stops it writing to disk
            cangoversub = yes
        }

        # Do we need to call gretroi?
        call_gretroi = no
        call_2_gretroi = no
        if (l_ret_roi) {

            if (cangoversub || call_gemexpr) {
                call_gretroi = yes
            }

            if (!l_fl_gtile) {
                # Need to call gretroi at the end regardless if not gtiling
                call_2_gretroi = yes
            }

            # Check if it's been rand before
            keypar (file_to_check//"[0]", "GRETROI", silent+)
            if (keypar.found) {
                keypar (file_to_check//"[0]", "RETROI", silent+)
                if (keypar.found) {
                    if (int(keypar.value) != l_req_roi && \
                        int(keypar.value) != 0) {
                            printlog ("ERROR - GMULTIAMP: GRETROI already \
                                ran on "//inimg//" and requested ROIs do \
                                match", l_logfile, verbose+)
                            status = 2
                            goto crash
                    } else if (int(keypar.value) == l_req_roi) {
                        call_gretroi = no
                    }
                } else {
                    printlog ("ERROR - GMULTIAMP: GREROI already \
                        ran on "//inimg//" but RETROI not found in PHU", \
                        l_logfile, verbose+)
                    status = 2
                    goto crash
                }
            }
        }

        if (debug) {
            printlog ("____flags: "//call_gretroi//" "//cangoversub//\
                " "//call_gemexpr//" "//l_fl_gtile//" "//call_2_gretroi//" "//\
                l_fl_pretty//" "//l_fl_signal, \
                l_logfile, verbose+)
        }

        # Check if there is anything to do for this image
        if (!cangoversub && !call_gemexpr && !l_fl_pretty && !l_fl_gtile \
            && !l_fl_signal && !call_2_gretroi) {
            printlog ("GMULTIAMP: There is nothing to be done for "//inimg, \
                l_logfile, verbose+)
            status = 1
            goto NEXTIMG
        }

        # Set up the input and output variables for each of the next steps and
        fparse (inimg)
        img_root = fparse.root//fparse.extension
        nextin = img_root
        nextout = outimg
        nextpath = l_rawpath//fparse.directory

        # Check if you can access l_sci_ext extensions - reset appropriately
        gemextn (file_to_check, check="exists,mef", index="", \
            extname=l_sci_ext, extversion="1-", ikparams="", \
            omit="", outfile="dev$null", logfile=l_logfile, glogpars="", \
            verbose=debug)

        if (gemextn.status != 0) {
            printlog ("ERROR - GMULTIAMP: GEMEXTn returned a non-zero \
                status.", l_logfile, verbose+)
            status = 2
            goto crash
        } else if (gemextn.count <= 0) {
            l_sci_ext = ""
        }

        if (call_gretroi) {
            grroi_in = nextin
            grroi_path = nextpath


            if (!cangoversub && !l_fl_pretty && !l_fl_gtile \
                && !call_2_gretroi) {
                grroi_out = nextout
            } else {
                grroi_out = tmproi//img_root
            }
            nextin = grroi_out
            nextpath = ""
        }

        if (cangoversub) {
            # Default values
            go_in = nextin
            go_path = nextpath
            go_calc_only = no
            go_out = nextout
            go_del = go_path//go_in

            # Check if l_fl_pretty to set go_calc_only properly and also the
            # input and output names
            if (call_gemexpr) {
                go_calc_only = yes
                go_del = ""
            } else if (l_fl_gtile || call_2_gretroi) {
                # Reset outputs and path
                go_out = tmpsub//img_root
                nextin = go_out
                nextpath = ""
            }
        }

        if (call_gemexpr) {
            p_path = nextpath
            nextpath = ""

            p_in = nextin
            p_in = p_path//p_in

            p_out = ""

            if (l_fl_gtile || call_2_gretroi) {
                p_out = tmppretty//img_root

                nextin = p_out

            } else if (!call_gretroi) {
                p_out = tmppretty//img_root
            } else {
                p_out = nextout
            }

            if (call_gretroi) {
               edit_image = p_in
            } else {
               edit_image = p_out
               nextin = tmppretty//p_out
            }

            if (p_out != "") {
                p_del = p_path//p_in
            }
        }

        # Need to call gretroi again in some circustances!
        if (l_fl_gtile || call_2_gretroi) {
            tile_in = nextin
            tile_out = nextout
            tile_path = nextpath

            print (tile_in//" "//tile_out//" 0", >> tile_list)
        }

        # Call the tasks

        all_oscan_vals = ""

        if (debug) {
            printlog ("____"//cangoversub//" "//call_gemexpr//" "//fl_debias//\
                " "//l_fl_pretty//" "//oscan_sub//" "//l_fl_gtile, \
                l_logfile, verbose+)
        }

        if (call_gretroi) {
            # Call gretroi - we need any bias sections so don't return just the
            # ROI.
            printlog ("GMULTIAMP: Removing any superfluous extensions...", \
                l_logfile, verbose=l_verbose)

            gretroi (grroi_in, outimage=grroi_out, outfile="", \
                rawpath=grroi_path, req_roi=l_req_roi, roi_only=no, \
                fl_vardq=no, sci_ext=l_sci_ext, var_ext=l_var_ext,
                dq_ext=l_dq_ext, mdf_ext=l_mdf_ext, key_detsec=l_key_detsec, \
                key_ccdsec=l_key_ccdsec, key_datasec=l_key_datasec, \
                key_ccdsum=l_key_ccdsum, logfile=l_logfile, verbose=debug)

            if (gretroi.status != 0) {
                printlog ("ERROR - GMULTIAMP: GRETROI returned an non-zero \
                    status", l_logfile, verbose+)
                status = 2
                goto crash
            }

            # If input sci_ext variable is not set reset it to SCI as this is
            # what happends in gretroi
            if (l_sci_ext == "") {
                l_sci_ext = "SCI"
            }
        }

        # Call goversub
        if (cangoversub) {
            printlog ("GMULTIAMP: Calculating overscan levels...", \
                l_logfile, l_verbose)

            # Perform the overscan - all checking done in here
            flpr; sleep 1  # FR 33692
            goversub (go_in, outimages=go_out, outprefix="", \
                rawpath=go_path, pixstat=l_pixstat, bias_type=l_bias_type, \
                sci_ext=l_sci_ext, var_ext=l_var_ext, dq_ext=l_dq_ext, \
                mdf_ext=l_mdf_ext, nbiascontam=l_nbiascontam, \
                key_biassec=l_key_biassec, statsec=l_statsec, \
                calc_only=go_calc_only, gaindb=l_gaindb, logfile=l_logfile, \
                verbose=debug)

            # Check the output status of goversub
            if (goversub.status >= 2) {
                # It failed
                printlog ("ERROR - GMULTIAMP: GOVERSUB returned a status \
                    value of 2 or more", l_logfile, verbose+)
                status = 2
                goto crash
            } else if (goversub.status == 1) {
                printlog ("GMULTIAMP: GOVERSUB returned a status value of \
                    1 - Image(s) were already overscan subtracted.", \
                    l_logfile, l_verbose)

                # There is no need to scan_sub if fl_pretty in this case
                oscan_sub = no
                isoverscansub = yes
            } else {
                isoverscansub = no
            }

            all_oscan_vals = goversub.oscan_val

            # Count the number of extensions (use the second field to determine
            # the number of words
            oscan_cnt = INDEF
            print (all_oscan_vals) | tokens ("STDIN", \
                ignore_comments=yes, begin_comment=yes, end_comment="eol", \
                newlines=yes, >> tmposcan)

            count (tmposcan) | fields ("STDIN", fields="2",\
                lines="1-", quit_if_missing=no, print_file_name=no) | \
                scan (oscan_cnt)

            if (isindef(oscan_cnt)) {
                printlog ("ERROR - GMULTIAMP: Cannot determine number of \
                    extensions from goverscan output", l_logfile, verbose+)
                status = 2
                goto crash
            }

            if (go_del != "" && go_del != file_to_check) {

                imdelete (go_del, verify-, >& "dev$null")
            }

        } # End of goversub loop

        if (debug) {
            date ("+%H:%M:%S.%N%n") | scan (sdate)
            printlog ("____At call gemexpr loop: "//sdate, l_logfile, verbose+)
        }

        # Call gemexpr
        if (call_gemexpr) {
            if (debug) {
                printlog ("____In call_gemexpr loop", l_logfile, verbose+)
            }

            # Need to obtain information from each extension

            # Need ---

            # May or may not need overscan values

            # Unique amplifier names - number of amps

            # Gain values (if gainmult and not isoverscansub
            # (/debiased/darksub) need to adjust ovserscan value?)

            # Detsecs to uniquely idenify extensions

            # Also, an additional debias values for amp boundries

            # If l_fl_sat need saturation values

            # Need to set up data base file - called tmpdb before calling
            # gemexpr

            # Check if image is prepared & determine number of extensions
            # Can use gain from headers - quicker

            isgainmult = no

            keypar (p_in//"[0]", "GPREPARE", silent+)
            if (keypar.found) {
                keypar (p_in//"[0]", "NSCIEXT", silent+)
                if (keypar.found) {
                    nsci = int(keypar.value)
                    
                    # Science extension names will be set
                    sciextn= l_sci_ext//","
                } else {
                    printlog ("ERROR - GMULTIAMP: Image is prepared but \
                        cannot access NSCIEXT keyword.", l_logfile, verbose+)
                    status = 2
                    goto crash
                }

                if (nsci != oscan_cnt) {
                    printlog ("ERROR - GMULTIAMP: Number of science \
                        counted do not match keyword.", \
                        l_logfile, verbose+)
                    status = 2
                    goto crash
                }

                # Determine which gain keyword to use
                l_gain_key = "GAINORIG"
                keypar (p_in//"["//sciextn//"1]", l_gain_key, silent+)
                if (!keypar.found) {
                    l_gain_key = "GAIN"
                }

                # Check if it's gain multiplied
                keypar (p_in//"["//sciextn//"1]", "GAINMULT", silent+)
                if (keypar.found) {
                    isgainmult = yes
                }

            } else {

                # If not prepared, assume the following:

                #     GAINS are wrong - need to call ggain
                #     It's not gain mulitplied

                l_gain_key = ""
                # Call ggain to get correct ggains and use it's out put to
                # determine number of extensions?
                ggain (p_in, gaindb=l_gaindb, logfile=l_logfile, \
                    key_gain=l_key_gain, key_ron=l_key_rdnoise, gain=2.2, \
                    ron=3.5, fl_update=no, fl_mult=no, \
                    sci_ext="", var_ext=l_var_ext, verbose=debug)

                if (ggain.status != 0) {
                    printlog ("ERROR - GMULTIAMP: GGAIN returned a non-zero \
                        status.", l_logfile, verbose+)
                    status = 2
                    goto crash
                }
                all_gain_vals = ggain.gainout

                # Count the number of extensions (use the second field to
                # determine the number of words
                gain_cnt = INDEF
                print (all_gain_vals) | tokens ("STDIN", \
                    ignore_comments=yes, begin_comment=yes, \
                    end_comment="eol", newlines=yes, >> tmpgains)

                count (tmpgains) | \
                    fields ("STDIN", fields="2", lines="1-", \
                    quit_if_missing=no, print_file_name=no) | \
                    scan (gain_cnt)

                if (isindef(gain_cnt)) {
                    printlog ("ERROR - GMULTIAMP: Cannot determine number of \
                        extensions from goverscan output", l_logfile, verbose+)
                    status = 2
                    goto crash
                }

                if (gain_cnt != oscan_cnt) {
                    printlog ("ERROR - GMULTIAMP: Number of science \
                        counted do not match keyword.", \
                        l_logfile, verbose+)
                    status = 2
                    goto crash
                } else {
                    nsci = gain_cnt
                }
            }

            sciextn = l_sci_ext
            if (sciextn != "") {
                sciextn = sciextn//","
            }

            if (debug) {
                printlog ("____nsci: "//nsci, l_logfile, verbose+)
            }

            if (l_gain_key == "") {
                l_gkey = "GMAGAIN"
            } else {
                l_gkey = l_gain_key
            }

            # Only copy the image if not ran through gretroi
            l_out_ext = l_sci_ext
            if (!call_gretroi) {
                imcopy (p_in//"[0]", p_out, verbose-)
            }

            if (l_out_ext == "") {
                l_out_ext = "SCI"
            }

            # Loop over images and form expression!
            for (ii = 1; ii <= nsci; ii += 1) {


                # Print the default nhedit parameters to tmpheader
                print ("default_pars update+ add+ before=\"\" \
                    after=\"\" delete- ", >> tmpheader)

                # Set current extension
                extn = p_in//"["//sciextn//ii//"]"

                if (debug) {
                    printlog ("____"//extn, l_logfile, verbose+)
                }

                # Get gain
                if (l_gain_key == "") {
                    fields (tmpgains, fields="1", lines=ii, \
                        quit_if_missing=no, print_file_name=no) | \
                        scan (l_gain)
                } else {
                    keypar (extn, l_gain_key, silent+)
                    ##M Put check in here
                    l_gain = keypar.value
                }

                # OSCAN value
                fields (tmposcan, fields="1", lines=ii, \
                    quit_if_missing=no, print_file_name=no) | \
                    scan (l_oscan)

                # Sort out the saturation value if needed
                sat_value = INDEF
                if (l_fl_sat && l_satval == "default") {
                    if (ii == 1) {
                        printlog ("GMULTIAMP: Determining saturation \
                            information...", l_logfile, l_verbose)
                    }

                    if (debug) {
                        date ("+%H:%M:%S.%N%n") | scan (sdate)
                        printlog ("____Calling gsat: "//sdate, \
                            l_logfile, verbose+)
                    }

                    gsat (p_in, "["//sciextn//ii//"]", satdb=l_satdb, \
                        gaindb=l_gaindb, bias=l_oscan, pixstat="midpt", \
                        statsec="default", gainval=l_gain, logfile=l_logfile, \
                        verbose=debug)

                    if (gsat.status != 0) {
                        printlog ("ERROR - GMULTIAMP: GSAT returned a \
                            non-zero status.", l_logfile, verbsoe+)
                        status = 2
                        goto crash
                    }
                    sat_value = gsat.saturation

                    if (debug) {
                        date ("+%H:%M:%S.%N%n") | scan (sdate)
                        printlog ("____Back from gsat: "//sdate, \
                            l_logfile, verbose+)
                    }

                } else if (l_fl_sat && l_satval != "default") {

                    test = fscan(l_satval,sat_value)

                    if (test != 1) {
                        printlog ("ERROR - GMULTIAMP: satval "//l_satval//\
                            " in not a number", l_logfile, verbose+)
                        status = 2
                        goto crash
                    }
                }

                if (!isindef(sat_value)) {

                    printf(l_fnum_str, sat_key, real(sat_value), "Saturation",\
                        >> tmpheader)
                } else {
                    l_fl_sat = no
                } # End of saturation loop

                # Sort out the signal value if needed
                if (l_fl_signal) {
                    if (ii == 1) {
                        printlog ("GMULTIAMP: Determining signal \
                            information...", l_logfile, l_verbose)
                    }

                    sig_value = l_signal
                    if (!isoverscansub) {
                        sig_value += l_oscan
                    }

                    printf(l_fnum_str, sig_key, sig_value, "Signal", \
                        >> tmpheader)
                } # End of signal loop

                if (l_fl_pretty || oscan_sub) {
                    printf (l_fnum_str, l_key_oscan, l_oscan, \
                        "Rough OVSERSAN value", >> tmpheader)
                }

                if (l_fl_pretty) {

                    if (ii == 1) {
                        printlog ("GMULTIAMP: Determining information to \
                            level amplifiers...", l_logfile, l_verbose)
                    }

                    # AMPNAME
                    keypar (extn, l_key_ampname, silent-)

                    ##M Put check in
                    print (keypar.value) | translit ("STDIN", \
                        from_string=" ", to_string="_", delete=no, \
                        collapse=no) | scan (l_ampname)

                    printf (l_fnum_str, l_gkey, real(l_gain), "GAIN", \
                        >> tmpheader)

                    # To determine sum of unique gains
                    print (l_ampname//" "//l_gain, >> tmpamps)

                    # Detemine the additional debias info
                    if (fl_debias) {

                        dx1 = INDEF
                        dx2 = INDEF

                        keypar (extn, l_key_datasec, silent+)
                        if (!keypar.found) {
                            printlog ("ERROR - GMULTIAMP: Cannot access "//\
                                l_key_datasec//" keyword in "//extn, \
                                l_logfile, verbose+)
                            goto crash
                        }
                        print (keypar.value) | scanf ("[%d:%d,%*d:%*d]", \
                            dx1, dx2)

                        if (isindef(dx1) || isindef(dx2)) {
                            printlog ("ERROR - GMULTIAMP: Cannot \
                                determine "//\
                                l_key_datasec//" X dimension values", \
                                l_logfile, verbose+)
                        }


                        # Read the ccdsec and see if at ampboundry
                        keypar (extn, l_key_ccdsec, silent+)
                        if (!keypar.found) {
                            printlog ("ERROR - GMULTIAMP: Cannot "//\
                                l_key_ccdsec//" in "//extn, \
                                l_logfile, verbose+)
                            status = 2
                            goto crash
                        }

                        print (keypar.value) | scanf ("[%d:%d,%*d:%*d]", \
                            tx1, tx2)

                        # xbin = 1 always - the extra debias doesn't is not
                        # dependant on binning
                        xbin = 1

                        read_debias = no

                        for (jj = 1; jj <= namps; jj += 1) {
                            if (tx1 == (1 + ((jj - 1) * (ccdwidth / namps)))) {
                                read_debias = yes
                                break
                            } else if (tx2 == (jj * (ccdwidth / namps))) {
                                read_debias = yes
                                break
                            }
                        }

                        nlines = 0

                        if (read_debias) {
                            if (debug) {
                                printlog ("Matching debias level", \
                                    l_logfile, verbose+)
                            }

                            match (l_ampname, files=debias_info, stop=no, \
                                print_file_name=no, metacharacters=no) | \
                                match ("#", files="STDIN", stop=yes, \
                                print_file_name=no, metacharacters=no, \
                                > tmpdet_info)

                             count (tmpdet_info) | scan (nlines)
                        }

                        if (nlines == 0 || !read_debias) {
                            if (debug) {
                                printlog ("Blank debias", \
                                    l_logfile, verbose+)
                            }

                            print ("INDEF 0 0 0 0", >> tmpdet_info)
                        }

                        ncols_lhs = INDEF
                        lhs_value = INDEF
                        ncols_rhs = INDEF
                        rhs_value = INDEF

                        if (debug) {
                            type (tmpdet_info)
                        }

                        scanfile4 = tmpdet_info
                        while (fscanf(scanfile4, "%s %d %f %d %f",\
                            l_ampname, ncols_lhs, lhs_value, ncols_rhs, \
                            rhs_value) != EOF){

                            if (!isindef(ncols_lhs) && !isindef(lhs_value)) {
                                col_limit = int (ncols_lhs / xbin) + dx1 - 1
                                sub_value = lhs_value
                            } else {
                                col_limit = 0
                                sub_value = 0
                            }

                            printf (l_inum_str, l_lhscol, col_limit, \
                                "Number of columns on LHS to remove extra \
                                bias from", >> tmpheader)
                            printf (l_fnum_str, l_lhsval, sub_value, \
                                "Bias value to subtract from the LHSNCOL", \
                                >> tmpheader)

                            if (!isindef(ncols_rhs) && !isindef(rhs_value)) {
                                col_limit = dx2 - (int(ncols_rhs / xbin) - 1)
                                sub_value = rhs_value
                            } else {
                                col_limit = 0
                                sub_value = 0.
                            }

                            printf (l_inum_str, l_rhscol, col_limit, \
                                "Number of columns on RHS to remove extra \
                                bias from", >> tmpheader)

                            printf (l_fnum_str, l_rhsval, sub_value, \
                                "Bias value to subtract from the RHSNCOL", \
                                >> tmpheader)

                            ncols_lhs = INDEF
                            lhs_value = INDEF
                            ncols_rhs = INDEF
                            rhs_value = INDEF

                        } # End of adding sections for additional bias to be

                        scanfile4 = ""
                        delete (tmpdet_info, verify-, >& "dev$null")
                    }
                } # End of fl_pretty loop debias

                # Create a fake image to store the new tmpheader
                outextn = edit_image//"["//l_out_ext//","//ii//",append+]"
                if (!call_gretroi) {
                    imcopy (extn, outextn, verbose-)
                }

                gemhedit (outextn, "", "", "", delete-, upfile=tmpheader)

                if (debug) {
                    type (tmpheader)
                }

                delete (tmpheader, verify-, >& "dev$null")

            } # End of for loop over extensions (ii)

            if (!call_gretroi) {
                # reset p_in and p_out
                p_in = p_out
                if (l_fl_gtile || call_2_gretroi) {
                    p_out = tmppretty//p_in
                } else {
                    p_out = outimg
                }
            }

            #### Set up the expression for gemexpr

            gexpr_sexpr = ""
            gexpr_vexpr = "a[VAR]"
            gexpr_dexpr = "a[DQ]"
            gexpr_svexpr = "a[VAR]"
            gexpr_vardq = yes

            # If pretty
            if (l_fl_pretty) {

                # Determine the number of unique amplifiers
                sort (tmpamps, column=0, ignore_white=yes, numeric_sort=no, \
                    reverse_sort=yes) | unique ("STDIN", > tmpunique)

                # Calculate the average gain from the unique amplifiers
                gtotal = INDEF

                fields (tmpunique, fields="2", lines="1-", \
                    quit_if_missing=no, print_file_name=no) | \
                    average | scan (gtotal)

                if (isindef(gtotal) || gtotal == 0.) {
                    printlog ("ERROR - GMULTIAMP: Cannot determine the number \
                        of unique amplifiers.", l_logfile, verbose+)
                    status = 2
                    goto crash
                }

                gexpr_sexpr = "(a)"

                if (!isgainmult) {
                    if (!isoverscansub && oscan_sub) {
                        gexpr_sexpr = "("//gexpr_sexpr//" - a."//\
                            l_key_oscan//")"
                    }

                    if (fl_debias) {

                        gexpr_sexpr = "(I <= a."//l_lhscol//" ? ("//\
                            gexpr_sexpr//" - a."//l_lhsval//") : (I >= ("//\
                            "(a."//l_rhscol//\
                            ")) ? ("//\
                            gexpr_sexpr//" - a."//\
                            l_rhsval//") : "//gexpr_sexpr//"))"
                    }

                    gexpr_sexpr = "("//gexpr_sexpr//" * a."//l_gkey//")"
                }

                gexpr_sexpr = "("//gexpr_sexpr//" / "//gtotal//")"
                gexpr_svexpr = "("//gexpr_svexpr//" / ("//gtotal//\
                    " * "//gtotal//"))"

            } else if (oscan_sub) {
                gexpr_sexpr = "a - a."//l_key_oscan
            } else {
                gexpr_sexpr = "a"
            }

            if (debug) {
                printlog ("____sci_expr: "//gexpr_sexpr, \
                    l_logfile, verbose+)

                printlog ("____var_expr: "//gexpr_vexpr, \
                    l_logfile, verbose+)

                printlog ("____dq_expr: "//gexpr_dexpr, \
                    l_logfile, verbose+)

            }

            # Print equations to database file
            if (l_fl_sat) {

                if (debug) {
                    date ("+%H:%M:%S.%N%n") | scan (sdate)
                    printlog ("____At second gemexpr call: "//sdate, \
                        l_logfile, verbose+)
                }

                gexpr_vexpr = "(a >= a."//sat_key//" ? 1 : 0)"

                if (l_fl_gtile || call_2_gretroi) {
                    sat_out = tmp_sat_out
                    print (sat_out//" "//l_sat_out//" 0", >> tile_list)
                } else {
                    sat_out = l_sat_out
                }

                gemexpr (gexpr_vexpr, sat_out, p_in, \
                    var_expr=gexpr_vexpr, dq_expr=gexpr_dexpr, \
                    sci_ext=l_out_ext,
                    var_ext=l_var_ext, dq_ext=l_dq_ext, mdf_ext=l_mdf_ext, \
                    fl_vardq=no, dims="default", intype="default", \
                    outtype="ushort", refim="default", rangecheck=yes, \
                    verbose=debug, exprdb="none", logfile=l_logfile, \
                    glogpars="")

                if (gemexpr.status != 0) {
                    printlog ("ERROR - GMULTIAMP: GEMEXPR returned a non-zero \
                        status.", l_logfile, verbose+)
                    status = 2
                    goto crash
                }

            } else {
                gexpr_vexpr = ""
            }

            if (l_fl_signal) {

                if (debug) {
                    date ("+%H:%M:%S.%N%n") | scan (sdate)
                    printlog ("____At third gemexpr call: "//sdate, \
                        l_logfile, verbose+)
                }

                gexpr_dexpr = "(a >= a."//sig_key//" ? 1 : 0)"

                if (l_fl_gtile || call_2_gretroi) {
                    sig_out = tmp_sig_out
                    print (sig_out//" "//l_sig_out//" 0", >> tile_list)
                } else {
                    sig_out = l_sig_out
                }

                gemexpr (gexpr_dexpr, sig_out, p_in, \
                    var_expr=gexpr_vexpr, dq_expr=gexpr_dexpr, \
                    sci_ext=l_out_ext, \
                    var_ext=l_var_ext, dq_ext=l_dq_ext, mdf_ext=l_mdf_ext, \
                    fl_vardq=no, dims="default", intype="default", \
                    outtype="ushort", refim="default", rangecheck=yes, \
                    verbose=debug, exprdb="none", logfile=l_logfile, \
                    glogpars="")

                if (gemexpr.status != 0) {
                    printlog ("ERROR - GMULTIAMP: GEMEXPR returned a non-zero \
                        status.", l_logfile, verbose+)
                    status = 2
                    goto crash
                }

            } else {
                gexpr_dexpr = ""
            }

            # Set the outype - have to go through gemexpr
            l_gem_outtype = "ref"
            if (l_fl_pretty || oscan_sub) {
                l_gem_outtype = "real"
            }

            if (debug) {
                date ("+%H:%M:%S.%N%n") | scan (sdate)
                printlog ("____At first gemexpr call: "//sdate, \
                    l_logfile, verbose+)
            }

            gemexpr (gexpr_sexpr, p_out, p_in, \
                var_expr=gexpr_svexpr, dq_expr="a[DQ]", \
                sci_ext=l_out_ext, \
                var_ext=l_var_ext, dq_ext=l_dq_ext, mdf_ext=l_mdf_ext, \
                fl_vardq=no, dims="default", intype="default", \
                outtype=l_gem_outtype, refim="default", rangecheck=yes, \
                verbose=debug, exprdb="none", logfile=l_logfile, \
                glogpars="")

            if (gemexpr.status != 0) {
                printlog ("ERROR - GMULTIAMP: GEMEXPR returned a non-zero \
                    status.", l_logfile, verbose+)
                status = 2
                goto crash
            }

            delete (tmpunique//","//tmpamps, verify-, >& "dev$null")

            if (p_del != "" && p_del != file_to_check) {

                imdelete (p_del, verify-, >& "dev$null")
            }

            # Delete extra file if needed
            if (!call_gretroi && p_path//p_in != file_to_check) {

                imdelete (p_in, verify-, >& "dev$null")
            }
        } # End of call_gemexpr loop

        delete (tmposcan//","//tmpgains, verify-, >& "dev$null")

        if (debug) {
            date ("+%H:%M:%S.%N%n") | scan (sdate)
            printlog ("____Finished gemexpr: "//sdate, l_logfile, verbose+)
        }

        # Call gtile (or call _gretroi again to remove biassections etc.)
        if (l_fl_gtile || call_2_gretroi) {

            if (l_fl_gtile && (l_fl_sat || l_fl_signal)) {
                printlog ("GMULTIAMP: Tiling image and mask(s)...", \
                    l_logfile, l_verbose)
            } else if (l_fl_gtile) {
                printlog ("GMULTIAMP: Tiling image...", \
                    l_logfile, l_verbose)
            }

            # Counter for print statements
            count (tile_list) | scan (nlines)
            jj = 1

            scanfile3 = tile_list
            while (fscan(scanfile3,tile_in,tile_out,l_repval) != EOF) {

                # call gtile
                if (l_fl_gtile) {

                    gtile (tile_in, outimages=tile_out, outpref="", \
                        out_ccds=l_out_ccds, ret_roi=l_ret_roi, \
                        req_roi=l_req_roi, fl_stats_only=l_fl_stats_only, \
                        fl_tile_det=l_fl_tile_det, fl_app_rois=l_fl_app_rois, \
                        fl_pad=l_fl_pad, sci_padval=l_repval, var_padval=0., \
                        dq_padval=16., sci_fakeval=l_repval, \
                        var_fakeval=0., dq_fakeval=16., chipgap="default", \
                        sci_ext=l_sci_ext, var_ext=l_var_ext, dq_ext=l_dq_ext,\
                        mdf_ext=l_mdf_ext, key_detsec=l_key_detsec, \
                        key_ccdsec=l_key_ccdsec, key_datasec=l_key_datasec, \
                        key_biassec=l_key_biassec, key_ccdsum=l_key_ccdsum, \
                        rawpath=tile_path, logfile=l_logfile, \
                        fl_verbose=debug)

                    if (gtile.status != 0) {
                        printlog ("ERROR - GMULTIAMP: GTILE returned a \
                            non-zero status", \
                            l_logfile, verbose+)
                        status = 2
                        goto crash
                    } else if (jj == 1) {
                        # Image will always be first
                        printlog ("               Image tiled.", \
                            l_logfile, l_verbose)
                    } else if (jj == nlines && nlines > 1) {
                        # Masks
                        printlog ("               Mask(s) tiled.", \
                            l_logfile, l_verbose)
                    }

                } else if (call_2_gretroi) {

                    # Return only the ROI data
                    gretroi (tile_in, outimage=tile_out, outfile="", \
                        rawpath=tile_path, req_roi=l_req_roi, \
                        roi_only=l_ret_roi, \
                        fl_vardq=no, sci_ext=l_sci_ext, var_ext=l_var_ext,
                        dq_ext=l_dq_ext, mdf_ext=l_mdf_ext, \
                        key_detsec=l_key_detsec, \
                        key_ccdsec=l_key_ccdsec, key_datasec=l_key_datasec, \
                        key_ccdsum=l_key_ccdsum, logfile=l_logfile,
                        verbose=debug)

                    if (gretroi.status != 0) {
                        printlog ("ERROR - GMULTIAMP: GRETROI returned a \
                            non-zero status", \
                            l_logfile, verbose+)
                        status = 2
                        goto crash
                    }
                }

                # Delete input images only if not the input image
                if (tile_path//tile_in != file_to_check) {

                    imdelete (tile_in, verify-, >& "dev$null")
                }

                # Counter for print statements
                jj += 1
            }

            delete (tile_list, verify-, >& "dev$null")
        }

        if (debug) {
            date ("+%H:%M:%S.%N%n") | scan (sdate)
            printlog ("____Finished gtile: "//sdate, l_logfile, verbose+)
        }

NEXTIMG:
    } # End of loop over input images

    goto clean

crash:
    if (access(l_crash_list)) {
        delete ("@"//l_crash_list, verify-, >& "dev$null")
        delete (l_crash_list, verify-,  >& "dev$null")
    }

    if (status != 2) {
        status = 2
    }

clean:

    delete (tmpinlist//","//tile_list, verify-, >& "dev$null")

    if (access(tmpoutlist)) {
        delete (tmpoutlist//", "//l_crash_list, verify-, >& "dev$null")
    }

    if (access(tmpoutlist2)) {
        delete (tmpoutlist2//", "//l_crash_list, verify-, >& "dev$null")
    }

    scanfile = ""
    scanfile2 = ""
    scanfile3 = ""
    scanfile4 = ""

    # Print finish time
    gemdate (zone="local")
    printlog ("\nGMULTIAMP - Finished: "//gemdate.outdate, \
        l_logfile, l_verbose)

    # This needs delaing with
    if (status == 0) {
        printlog ("GMULTIAMP - Exit status: GOOD\n", l_logfile, l_verbose)
    } else if (status == 1) {
        printlog ("GMULTIAMP - Exit status: GOOD - no processing performed \
            on one or more images\n", l_logfile, verbose+)
    } else {
        printlog ("GMULTIAMP - Exit status: ERROR\n", l_logfile, verbose+)
    }

    if (profile) {
        date ("+%H:%M:%S.%N%n") | scan (sdate)
        printlog ("____GMULTIAMP Finished: "//sdate, l_logfile, verbose+)
    }

end

# Copyright(c) 2010-2013 Association of Universities for Research in Astronomy, Inc.

procedure gaprepare (inimages)

##M Set up flag for ODGW if they will always be parked in the corner of the
##M     physical array if in detector or array mode?
##M Check for frozen ODGW - flag here
##M     -- Write ODGW CCDSEC to SCI extension -- only if guiding
##M Make it backwards compatible for commissioning -- timestamps
##M     -- Have it doing a check at moment but not entirely reliable
##M Add CCDSEC and DATASEC parameters...

# INPUTS:
# - GSAOI raw data with 4 extensions

# What this task does:
# - update any header keyword that needs updating
# - do the proper trimming as decided by the instrument team
# - applies the non-linear correction for each array
#   the linearity correction is: ADU_out = a*ADU_in + b*ADU_in**2 + c*ADU_in**3
#   where a, b and c are in the file given in arraysdb. Currently the same
#   coefficients are used for all pixels values.
# - create the basic VAR and DQ planes
# - flags non-linear and saturated pixels in the DQ plane if requested
# - adds a header keyword (METACONF) we think useful for the data reduction

# The input file names can be in the form of:
# - "@list" containing a list of files one per line, ".fits" is not required
#   but can be present. The "@" can be anywhere in the parameter.
# - comma-separated list
# - one or more range of image numbers, in the format allowed by gemlist, in
#   which case current UT date is assumed for the rootname
# - Technically all of the above (including file names) can be in a
#   comma-separated list!

# Outputs:
# - The output frames are by default prefixed with a "g", this can be changed
#   with the outpref parameter.
# - The keyword GAPREPAR is added with the time stamp when the task is last run

# NOTE: that there is NO outimages list

# Version 2010 May 20 V0.0 CW - created
# Version 2010 Jun 10 V1.0 CW - working version, but with TODOs
# Version 2011 Feb 24 V1.5 CW - added DQ plane from static BPM, plus
#                               sat/nonlinear
# Version 2011 May 17 V1.6 CW - added non-linearity correction for testing
#                               (per detector)
# Version 2011 Jun 23 V1.7 CW - final linearity correction, saturation and 95%
#                               linearity levels
# Version 2011 Sep 27 V1.8 CW - handling of ROIs
# Version 2011 Dec 07 V1.9 CW - modify METACONF to handle ROIs
# Version 2012 Jan 08 V2.0 CW - correct filter unknown to clear if eng position
#                               is correct but warn that the other filter may
#                               be wrong. If the eng position is wrong, skip
#                               with a warning. Also flags the active and
#                               parked ODGW in the DQ.
#
# Descriptions of future revisions are recorded in CVS logs

######################################################

char    inimages    {prompt="Raw GSAOI images; image name(s), range or list"}
char    rawpath     {"",prompt="Path to input raw images"}
char    outpref     {"g",prompt="Prefix for output images"}
char    rootname    {"",prompt="Root name if supplying image number(s); blank for today's UT"}
bool    fl_trim     {yes,prompt="Trim the images?"}
bool    fl_nlc      {yes,prompt="Apply non-linear correction to data?"}
bool    fl_vardq    {no,prompt="Create variance and data quality frames?"}
bool    fl_sat      {no,prompt="Include non-linear and saturated pixels in data quality planes?"}
char    badpix      {"gsaoi$data/gsaoibpm_high_full.fits",prompt="Static (MEF) Bad Pixel Mask"}
char    sci_ext     {"SCI",prompt="Name of science extensions"}
char    var_ext     {"VAR",prompt="Name of variance extensions"}
char    dq_ext      {"DQ",prompt="Name of data quality extensions"}
char    arraysdb    {"gsaoi$data/gsaoiAMPS.dat",prompt="Database file containing array information"}
char    non_lcdb    {"gsaoi$data/gsaoiNLC.dat",prompt="Database file containing non-linearity coefficients"}
char    logfile     {"",prompt="Logfile"}
bool    verbose     {yes,prompt="Verbose?"}
int     status      {0,prompt="Exit status (0=good)"}
struct  *scanfile   {"",prompt="Internal use only"}

######################################################

begin

    ########
    # Declare and set local variables; set default values; initiate temporary
    # files

    ####
    # Variable declaration
    char l_inimages, l_rawpath, l_outpref, l_logfile, l_rootname, l_badpix
    char trimsec, inphu, outphu, datasec, filelist, inimg, outimg, utdate
    char nlexpr, varexpr, metakey, obskey, filt1key, filt2key, itimekey
    char coaddkey, gcalkey, lnrskey, ccdsec, roiid, objct, inlist, t_string
    char rphend, pathtest, inimg_root, filter1, filter2, tmplist, l_sci_ext
    char l_var_ext, l_dq_ext, tmpcrash, outname, outimg_orig, inname, readmode
    char gwfs_value, gwfs_cfg, inimg_orig, sciext, sciappend, varextn
    char varappend, dqextn, dqappend, bpmsec, bpmexp, dqexpr, ccdname
    char tmpnotproc, sline, op_toprint, bpm_toprint, detmode, filt_string
    char time_string, object_raw, outgaimchk_list, tmpfile, l_arraysdb
    char high_range, l_non_lcdb, bpmext
    char l_timename, uttime, l_datename, datevalue

    int ii, jj, comma_pos, coadds, lnrs, num_created, satlev, nonlin, satval
    int lcount, low_nlc_val, odgwx[4], odgwy[4], odgws, odgwfd
    int atlocation, ccdx1, ccdx2, ccdy1, ccdy2, naxis1
    int datax1, datax2, datay1, datay2, detx1, detx2, dety1, dety2
    int bpmccdx1, bpmccdx2, bpmccdy1, bpmccdy2, bpmx1, bpmx2, bpmy1, bpmy2
    int odgw_x1, odgw_x2, odgw_y1, odgw_y2, extn, nextend, trimcount, nlccount
    int bpmnext, trimx1, trimx2, trimy1, trimy2, pad, shift, obs_time

    real nla, nlb, nlc, nlla, nllb, nllc, nonln5, nonln2
    real rdnval, gainval
    real coeff1, coeff2, coeff3, xx, lcoeff1, lcoeff2, lcoeff3
    real rdnoise, gain, rdnorig, rdcor, gainext

    bool l_verbose,l_fl_trim, l_fl_nlc, l_fl_vardq, l_fl_sat, fixed_headers
    bool fl_odgw, prepared, istrimmed, alltrimmed, allnlc, commissioning
    bool fl_trim_orig, fl_nlc_orig, fl_vardq_orig, do_low_nlc_sep
    bool debug

    struct comm_struct
    struct comm_test="Number of Fowler samples"

    ####
    # Set temporary files
    filelist = mktemp("tmpfile")
    tmplist = mktemp("tmplist")
    tmpcrash = mktemp("tmpcrash")
    tmpnotproc = mktemp("tmpnotproc")
    outgaimchk_list = mktemp("tmpimchk")

    ####
    # Set local variables
    l_inimages = inimages
    l_rawpath = rawpath
    l_outpref = outpref
    l_rootname = rootname
    l_fl_trim = fl_trim
    l_fl_nlc = fl_nlc
    l_fl_vardq = fl_vardq
    l_fl_sat = fl_sat
    l_badpix = badpix
    l_sci_ext = sci_ext
    l_var_ext = var_ext
    l_dq_ext = dq_ext
    l_arraysdb = arraysdb
    l_non_lcdb = non_lcdb
    l_logfile = logfile
    l_verbose = verbose

    ####
    # Set default values
    debug = no
    status = 0
    comma_pos = 0
    nonlin = 0
    satval = 0
    fl_odgw = no
    prepared = no

    # Set this to yes if the headers are fixed that Folwer read become true
    # LNRs for commissioning data - may be better to this with a time stamp -MS
    fixed_headers = no

    # If it is decided to apply the non-linear correction to <2000 ADU
    # separartely set this to yes; change all occurances of ALL in the non_lcdb
    # file to HIGH and add lines for the <2000 coeffients with LOW in the range
    # column - MS
    do_low_nlc_sep = no
    low_nlc_val = 2000

    ########
    # Information on the values that characterise the instrument - they are
    # stored in files given in arraysdb and non_lcdb.

    # Non-lineraity corection coefficients
    #
    # The non-linear parameters for each array as derived using
    # the GNIRS equivalent method - the correction is a quadratic fit
    # a+b*cts+c*cts*cts - no normalization applied (not using irlincor).
    # Note that the applied correction is actually a third order poly.

    # Non-linearity percentages
    #
    # nonln5 correspond to the % of saturation value where counts depart >5%
    # from the linearity BEFORE non-linearity correction.
    #
    # nonln2 correspond to the % of saturation value where counts depart >2%
    # from the linearity AFTER the non-linearity correction is applied
    # This is indexed by extension

    # Saturation values
    #
    # The saturation values are BEFORE non-linearity correction in [ADU]

    # Readnoise and gain are in the headers - RDNOISE, GAIN (each extension)
    #
    # Readnoise is dependent of read mode (bright, faint, very faint)
    # Bright -> LNRS = 2
    # Faint -> LNRS =  8
    # Very faint -> LNRS = 16
    #
    # IMPORTANT - these values are in [ADU] -- convert to [e-] when using
    # IMPORTANT - these values are for COADD = 1 -- final rdnoise is
    #             rdnorig*sqrt(coadds)
    #
    # Gain is units [e-/ADU]

    ########
    # Here is where the actual work starts

    # Test the log file
    if ((l_logfile == "") || (stridx(" ", l_logfile) > 0)) {

        l_logfile = gsaoi.logfile
        if ((l_logfile == "") || (stridx (" ", l_logfile) > 0)) {
            l_logfile = "gsaoi.log"
            printlog ("WARNING - GAPREPARE: Both gaprepare.logfile and \
                gsaoi.logfile fields are empty", l_logfile, l_verbose)
            printlog ("                     Using default file gsaoi.log", \
                l_logfile, l_verbose)
        }
    }

    # Print start time
    gemdate (zone="local")
    printlog ("\n------------------------------------------------------------",
        l_logfile, l_verbose)
    printlog ("GAPREPARE -- Started: "//gemdate.outdate//"\n", \
        l_logfile, l_verbose)

    # Record all input parameters
    printlog ("GAPREPARE: Input parameters...", l_logfile, l_verbose)
    printlog ("    inimages  =  "//l_inimages, l_logfile, l_verbose)
    printlog ("    rawpath   =  "//l_rawpath, l_logfile, l_verbose)
    printlog ("    outpref   =  "//l_outpref, l_logfile, l_verbose)
    printlog ("    rootname  =  "//l_rootname, l_logfile, l_verbose)
    printlog ("    fl_trim   =  "//l_fl_trim, l_logfile, l_verbose)
    printlog ("    fl_nlc    =  "//l_fl_nlc, l_logfile, l_verbose)
    printlog ("    fl_vardq  =  "//l_fl_vardq, l_logfile, l_verbose)
    printlog ("    fl_sat    =  "//l_fl_sat, l_logfile, l_verbose)
    printlog ("    badpix    =  "//l_badpix, l_logfile, l_verbose)
    printlog ("    sci_ext   =  "//l_sci_ext, l_logfile, l_verbose)
    printlog ("    var_ext   =  "//l_var_ext, l_logfile, l_verbose)
    printlog ("    dq_ext    =  "//l_dq_ext, l_logfile, l_verbose)
    printlog ("    arraysdb  =  "//l_arraysdb, l_logfile, l_verbose)
    printlog ("    non_lcdb  =  "//l_non_lcdb, l_logfile, l_verbose)
    printlog ("    logfile   =  "//l_logfile, l_logfile, l_verbose)
    printlog ("    verbose   =  "//l_verbose, l_logfile, l_verbose)
    printlog ("", l_logfile, l_verbose)

    ####
    # Start user input checks and data validation

    # Check the l_badpix value
    if ((l_badpix == "") || (stridx(" ", l_badpix) > 0)) {
        # This assumes that the bpm covers all of the areas readout in each
        # detector in each image
        l_badpix = "none"
    } else if (!imaccess(l_badpix)) {
        printlog ("ERROR - GAPREPARE: cannot access bpm "//l_badpix//\
            ". Exiting", l_logfile, verbose+)
        goto crash
    } else {
        # Check it hasn't been mosaiced
        keypar (l_badpix//"[0]", "GAMOSAIC", silent+)
        if (keypar.found) {
            printlog ("ERROR - GAPREPARE: "//l_badpix//" has been mosaiced."//\
                " Please supply a static BPM that has not been mosaiced.", \
                l_logfile, verbose)
            goto crash
        }
    }
    printlog ("GAPREPARE: badpix value after checks is "//l_badpix, \
        l_logfile, l_verbose)

    # Check the validity of the arraysdb database file
    if ((l_arraysdb == "") || (stridx(" ", l_arraysdb) > 0)) {
        printlog ("ERROR - GAPREPARE: arraysdb not set."//\
            ". Exiting", l_logfile, verbose+)
        goto crash
    } else if (!access(l_arraysdb)) {
        printlog ("ERROR - GAPREPARE: cannot access "//l_arraysdb//\
            ". Exiting", l_logfile, verbose+)
        goto crash
    }

    # Check the validity of the non_lcdb database file
    if ((l_non_lcdb == "") || (stridx(" ", l_non_lcdb) > 0)) {
        printlog ("ERROR - GAPREPARE: non_lcdb not set."//\
            ". Exiting", l_logfile, verbose+)
        goto crash
    } else if (!access(l_non_lcdb)) {
        printlog ("ERROR - GAPREPARE: cannot access "//l_non_lcdb//\
            ". Exiting", l_logfile, verbose+)
        goto crash
    }

    # Call gaimchk to perform input checks
    gaimchk (inimages=l_inimages, rawpath=l_rawpath, rootname=l_rootname, \
        obstype_allowed="", object_allowed="", \
        key_allowed="", key_forbidden="GAPREPAR", key_exists="", \
        fl_prep_check=yes, gaprep_pref=l_outpref, fl_redux_check=no, \
        garedux_pref="r", fl_fail=yes, fl_out_check=no, fl_vardq_check=no, \
        sci_ext=l_sci_ext, var_ext=l_var_ext, dq_ext=l_dq_ext, \
        outlist=outgaimchk_list, logfile=l_logfile, verbose=no)

    if (gaimchk.status != 0) {
        printlog ("ERROR - GAPREAPRE: GAIMCHK returned a non-zero status. \
            Exiting.", l_logfile, verbose+)
        goto crash
    } else if (!access(outgaimchk_list)) {
        printlog ("ERROR - GAPREPARE: Cannot access output list from \
            GAIMCHK", l_logfile, verbose+)
        goto crash
    }

    # Files that have been processed by GAPREPARE
    tmpfile = ""
    match ("tmpGAPREPAR", outgaimchk_list, stop=no, print_file_name=no, \
        metacharacter=no) | scan (tmpfile)
    if (tmpfile != "") {
        printlog ("WARNING - GAPREPARE: The following files have already been \
            processed by GAPREPARE."//\
            "\n                     Ignoring these files.\n", \
            l_logfile, verbose+)
        fields (tmpfile, fields=1, lines="1-", quit_if_miss=no, \
            print_file_names=no) | \
            tee (l_logfile, out_type="text", append+)
    }

    # Files that need to be prepared
    tmpfile = ""
    match ("tmpNotFound", outgaimchk_list, stop=no, print_file_name=no, \
        metacharacter=no) | scan (tmpfile)
    if (tmpfile != "") {

        fields (tmpfile, fields=1, lines="1-", quit_if_miss=no, \
            print_file_names=no, >> filelist)
    } else {
        printlog ("ERROR - GAPREPARE: There are no files to prepare. \
            Exiting", l_logfile, verbose+)
        goto crash
    }

    # Consistency check for flagging saturated and non-linear pixels
    if (l_fl_sat && !l_fl_vardq) {
        printlog ("WARNING - GAPREPARE: fl_sat is yes but fl_vardq is no. \
            fl_vardq must be set to\n"//\
            "                     yes to flag saturated and non-linear \
            pixels.\n"//\
            "                     Setting fl_sat to no.", \
            l_logfile, verbose+)
        fl_sat = no
    }

    fl_trim_orig = l_fl_trim
    fl_nlc_orig = l_fl_nlc
    fl_vardq_orig = l_fl_vardq

    num_created = 0

    # Loop over input images to process them
    scanfile = filelist
    while (fscan(scanfile, inimg) != EOF) {

        # Reset the flags
        l_fl_trim = fl_trim_orig
        l_fl_nlc = fl_nlc_orig
        l_fl_vardq = fl_vardq_orig

        # Set istrimmed flag to no - this is for an individual extension
        istrimmed = no

        # Set alltrimmed flag to no - this is for the entire image
        alltrimmed = no

        # Set allnlc flag to no - this is for the entire image
        allnlc = no

        # Parse the file name, set input name and set output name
        fparse (inimg, verbose-)
        inimg_root = fparse.root
        inname = inimg_root//fparse.extension
        outname = l_outpref//inimg_root//".fits"

        # Store original input name - inimg gets updated later on depending on
        # operations desired
        inimg_orig = inimg

        printlog ("\n----", l_logfile, l_verbose)
        printlog ("GAPREPARE: Processing image "//inname//" into "//\
            outname//"\n",\
            l_logfile, l_verbose)

        # Check if it's been trimmed already
        if (l_fl_trim) {
            keypar (inimg//"[0]","TRIMMED", silent+)
            if (keypar.found) {
                printlog ("WARNING - GAPREPARE: Image "//inname//\
                    " already trimmed.", \
                    l_logfile, l_verbose)
                alltrimmed = yes
                l_fl_trim = no
            }
        }

        # Check if it's already been non-linear corrected
        if (l_fl_nlc) {
            keypar (inimg//"[0]","LINRCORR", silent+)
            if (keypar.found) {
                printlog ("WARNING - GAPREPARE: Image "//inname//\
                    " already corrected for non-linearity.", \
                    l_logfile, l_verbose)
                allnlc = yes
                l_fl_nlc = no
            }
        }

        # Set up the output names here
        if (l_fl_nlc) {
            # Create temporary filename and store it in outimg_orig
            outimg_orig = mktemp("tmpout")//".fits"
        } else {
            outimg_orig = outname
        }

        # Set outimg to outimg_orig
        outimg = outimg_orig

        # Set PHU to read required information
        outphu = outimg//"[0]"

        # If the image is prepared (shouldn't actually get here) all keywords
        # that need to be updated regardless of trimming, nlc and vardq will
        # have already been done in the previous pass
        if (prepared && !l_fl_trim && !l_fl_nlc && !l_fl_vardq) {
            printlog ("WARNING - GAPREPARE: There is nothing to be done for \
                image "//inname//\
                ".\n                     Moving on to next image.", \
                l_logfile, l_verbose)
            print (inimg, >> tmpnotproc)
            goto NEXTIMAGE
        }

        # Check for the existance of VAR and DQ planes
        gemextn (inimg, check="exists", process="expand", index="", \
            extname=l_var_ext//","//l_dq_ext, extversion="1-", \
            ikparams="", omit="", replace="", outfile="dev$null", \
            logfile=l_logfile, glogpars="", verbose=no)

        if (gemextn.status != 0) {
            printlog ("ERROR- GAPREPARE: GEMEXTN returned a non-zero \
                status. Exiting.", l_logfile, verbose+)
                goto crash
        } else if (gemextn.count > 0 && l_fl_vardq) {
            printlog ("WARNING - GAPREPARE: Image "//inname//\
                " already contains "//l_var_ext //" and "//l_dq_ext//\
                " extensions.\n"//\
                "                     These extensions will be recalculated"//\
                " for this image.", l_logfile, verbose+)
        }

        # Copy PHU to outimg to allow updating of FILTER keywords
        fxcopy (inimg, outimg, groups="0", new_file=yes, verbose-)
        print (outimg, > tmpcrash)

        # Check for the BUNIT keyword (old data has BUNITS)
        keypar (outimg//"[0]", "BUNIT", silent+)
        if (!keypar.found) {
            keypar (outimg//"[0]", "BUNITS", silent+)
            if (keypar.found) {
                # Delete BUNITS keyword
                gemhedit (outimg//"[0]", "BUNITS", "", "", delete+)
            }

            # This assumes only raw data comes through this point - MS
            gemhedit (outimg//"[0]", "BUNIT", "ADU", \
                "Units of quantity read from detector array", delete-)
        }

        # Check FILTER1 and FILTER2 for the UNKNOWN issue.
        # This is a pending issue - sometimes the exposure is taken before the
        # filters are in position. The indication is that usually one of them
        # is Unknown (usually corresponding to clear), but sometimes the other
        # may not be in position either.

        # Read filter1 keyword
        keypar (outphu, "FILTER1", silent+)
        if (!keypar.found) {
            printlog ("ERROR - GAPREPARE: FILTER1 keyword not found in "//\
                inname//"[0]", l_logfile, verbose+)
            goto crash
        }

        filter1 = keypar.value
        if (filter1 == "Unknown") {

            keypar (outphu, "FILT1POS", silent+)
            if (keypar.value == "200") {
                printlog ("WARNING - GAPREPARE: FILTER1 in image "//inname//\
                    " was Unknown but at Clear Eng position.", \
                    l_logfile, l_verbose)
                printlog ("                     Header set to clear, but \
                    be aware that FILTER2 may also be wrong.", \
                    l_logfile, l_verbose)

                gemhedit (outphu, "FILTER1", "Clear", "", delete-)

            } else {
                printlog ("WARNING - GAPREPARE: FILTER1 in image "//inname//\
                    " was Unknown and Eng positions does not", \
                    l_logfile, l_verbose)
                printlog ("                     correspond to Clear. Skipping \
                    this image",l_logfile, l_verbose)
                print (inimg, >> tmpnotproc)
                goto NEXTIMAGE
            }
        }

        # Read FILTER2 keyword
        keypar (outphu, "FILTER2", silent+)
        if (!keypar.found) {
            printlog ("ERROR - GAPREPARE: FILTER2 keyword not found in "//\
                inname//"[0]", l_logfile, verbose+)
            goto crash
        }

        filter2 = keypar.value
        if (filter2 == "Unknown") {
            keypar (outphu, "FILT2POS", silent+,>& "dev$null")
            if (keypar.value == "6933") {
                printlog ("WARNING - GAPREPARE: FILTER2 in image "//inname//\
                    " was Unknown but at Clear Eng position.", \
                    l_logfile, l_verbose)
                printlog ("                     Header set to clear, but be \
                    aware that FILTER1 may also be wrong.", \
                    l_logfile, l_verbose)
                gemhedit (outphu, "FILTER2", "Clear", "", delete-)

            } else {
                printlog ("WARNING - GAPREPARE: FILTER2 in image "//inname//\
                    " was Unknown and Eng positions does not", \
                    l_logfile, l_verbose)
                printlog ("                     correspond to Clear. Skipping \
                    this image",l_logfile, l_verbose)
                print (inimg, >> tmpnotproc)
                goto NEXTIMAGE
            }
        }

        # Read the number of extensions - for NSCIEXT first though
        keypar (inimg//"[0]", "NSCIEXT", silent+)
        if (!keypar.found) {
            keypar (inimg//"[0]", "NEXTEND", silent+)
            if (!keypar.found) {
                printlog ("WARNING - GAPREPARE: NEXTEND and NSCIEXT keywords \
                    not found in "//inname//"[0], assuming 4 extensions \
                    present.", l_logfile, verbose+)
                nextend = 4
            } else {
                nextend = int (keypar.value)
            }
        } else {
            nextend = int (keypar.value)
        }

        ####
        # Determining the ODGW configuration (only if creating / propagting
        # var/dq planes) - All of this information is in the PHU hence why
        # determining the configuration here.

        if (l_fl_vardq) {

            # If GWFS?CFG is ODGWa, where '?' is 1 to 4 and 'a' is 1 to 4, but
            # not necessarily the same. '?' is the GWFS configuration, 'a' is
            # the detector where the ODGW that is actively guiding is located.
            # Need to read GWFS?X, GSWFS?Y as the central pixel of the ODGW and
            # GWFS?SIZ as the size of the box. That goes in the BPM as 1 later
            # on if creating vardq planes. The ones that are not active will
            # be parked but reading, so need to be flagged as well. They will
            # be parked in the the corner furthest away from the centre of the
            # entire array.

            # This flag is used later on in the creation of the var / dq planes
            fl_odgw = no

            # Set default values for ODGW X and Y positions
            for (ii = 1; ii <= nextend; ii += 1) {
                # Can use zero as coordinate system is 1 based
                odgwx[ii] = 0
                odgwy[ii] = 0
            }

            # Need to loop over the number of extensions at least; not sure if
            # they
            # will always be within the first 4? I think this could be 5?
            for (ii = 1; ii <= nextend; ii += 1) {

                # Set default values
                gwfs_value = ""

                # Read the gwfs configuration
                gwfs_cfg = "GWFS"//ii//"CFG"

                keypar (outphu, gwfs_cfg, silent+)
                if (keypar.found) {
                    gwfs_value = str(keypar.value)
                }

                # Test the gwfs_value for the ODGW string
                if (strstr("ODGW", gwfs_value) != 0) {

                    # Determine the size of the ODGW, it is the side length of
                    # the box. All ODGW parked or guiding are teh same size.
                    keypar (outphu, "GWFS"//ii//"SIZ", silent+)
                    if (!keypar.found) {
                        printlog ("ERROR - GAPREPARE: GWFS"//ii//"CFG has a \
                            value of "//gwfs_value//" but GWFS"//ii//\
                            "SIZ is not present. Exiting\n", \
                            l_logfile, verbose+)
                        goto crash
                    }

                    if (!fl_odgw) {
                        # This assumes that all ODGW are the same size and a
                        odgws = int(keypar.value) / 2
                    } else if (fl_odgw && (odgws != int(keypar.value) / 2)) {
                        # This is to check that they are the same - they should
                        # always be the same - to get here means there is a
                        # problem. Only check if fl_odgw is true - MS
                        printlog ("ERROR - GAPREPARE: ODGW sizes do not \
                            match. Exiting.", l_logfile, verbose+)
                        goto crash
                    }

                    # Parse the extension described by this gwfs configuration
                    # To be here gwfs_value will be "ODGW"//a as described
                    # above
                    extn = int(substr(gwfs_value, 5, strlen (gwfs_value)))

                    # Read the X and Y positions of the ODGW
                    keypar (outphu, "GWFS"//ii//"X", silent+)
                    odgwx[extn] = int(keypar.value)

                    keypar (outphu, "GWFS"//ii//"Y", silent+)
                    odgwy[extn] = int(keypar.value)

                    # Set the odgw flag
                    fl_odgw = yes
               }
            } # End of loop over extensions
        } # End of l_fl_vardq loop to determine ODGW configuration

        # Read COADDS and LNRS kewords to be used later on to update keywords
        # and calculate updated readnoises later on
        keypar (outphu, "COADDS", silent+)
        if (!keypar.found) {
            printlog ("ERROR - GAPREPARE: COADDS keyword not found in "//\
                inname//"[0]. Exiting", l_logfile, verbose+)
            goto crash
        }
        coadds = int(keypar.value)

        keypar (outphu, "LNRS", silent+)
        if (!keypar.found) {
            printlog ("ERROR - GAPREPARE: LNRS keyword not found in "//\
                inname//"[0]. Exiting", l_logfile, verbose+)
            goto crash
        }
        lnrs = int(keypar.value)

        # Make it backwards compatible - LNRS is number of low noise reads
        # which is 2 times the number of Fowler reads but in addition the
        # definitions before this change were changed too - MS
        commissioning = no

        if (!fixed_headers) {

            l_timename = "UT"
            # The observation time of the flat
            keypar (outphu, l_timename, silent+)
            if (!keypar.found) {
               printlog ("ERROR - GAPREPARE: "//l_timename//" keyword not \
                   found in "//inname//"[0]", \
                   l_logfile, verbose+)
               goto crash
            }
            uttime = keypar.value

            l_datename = "DATE-OBS"
            keypar (outphu, l_datename, silent+)
            if (!keypar.found) {
               printlog ("ERROR - GAPREPARE: "//l_datename//" keyword not \
                   found in "//inname//"[0]",
                   l_logfile, verbose+)
               goto crash
            }
            datevalue = keypar.value

            cnvtsec (datevalue, uttime) | scan (obs_time)

            # Need to deal with the LNRS keyword differently depending on when
            # the file is from.
            if (obs_time >= 1036195200) {
                # "2012-11-01" "00:00:00"
                # LNRS refers correctly to the number of low read noises for
                # each mode.
                lnrs = lnrs

            } else {

                # LNRS refers to number of Fowler reads (in additon there is
                # the fact that the very faint and very, very modes had the
                # number of reads reduced by a factor of 2 at some point

                if (lnrs == 1 || lnrs == 2) {
                    lnrs = 2
                } else {
                    lnrs = 2 * lnrs
                }

                commissioning = yes
            }
        }

        # Determine the readmode - these values are for data after
        # commissioning - may not neccessarily be backwards compatible. - MS
        if (lnrs == 2) {
            readmode = "bright"
        } else if (lnrs == 8) {
            readmode = "faint"
        } else if (lnrs == 16 || lnrs == 32) {
            readmode = "very_faint"
        } else {
            printlog ("ERROR - GAPREPARE: LNRS value "//lnrs//", not \
                recognised. Exiting", l_logfile, verbose+)
            goto crash
        }

        # Print warning that the information may be wrong for this readmode
        if (commissioning && !fixed_headers) {

            if (lnrs == 16 || lnrs == 32) {
                printlog ("WARNING - GAPREPARE: Image \""//inimg//\
                    " appears to be commissing data with 8 or 16 Folwer \
                    reads."//\
                    "\n                     Gain, readnoise and non-linear "//\
                    "correction coefficients may be incorrect.", \
                    l_logfile, verbose+)
            }

            if (lnrs == 16) {
                printlog ("                     If readmode is known "//\
                    "to be "//\
                    "\"faint\" reset LNRS keyword to 4 to obtain the "//\
                    "latest values and run again.\n", l_logfile, verbose+)
            }
        }

        # Reset counters
        trimcount = 0
        nlccount = 0

        # Loop over extensions then apply things that may need to be done
        for (ii = 1; ii <= nextend; ii +=1) {

            # Set the output image to write to
            outimg = outimg_orig

            # Reset inimg
            inimg = inimg_orig

            # Set PHU to read required information
            outphu = outimg//"[0]"

            # Set the image extension to read
            if (prepared) {
                sciext = "["//l_sci_ext//","//ii//"]"
            } else {
                sciext = "["//ii//"]"
            }

            # Set the image extension to append to
            sciappend = "["//l_sci_ext//","//ii//",append]"

            # Read CCDSEC
            keypar (inimg//sciext, "CCDSEC", silent+)
            if (!keypar.found) {
                printlog ("ERROR - GAPREPARE: CCDSEC not found in "//\
                    inimg//sciext, l_logfile, verbose+)
                goto crash
            }
            ccdsec = keypar.value
            print (ccdsec) | \
                scanf ("[%d:%d,%d:%d]", ccdx1, ccdx2, ccdy1, ccdy2)

            # Read DATASEC
            keypar (inimg//sciext, "DATASEC", silent+)
            if (!keypar.found) {
                printlog ("ERROR - GAPREPARE: DATASEC not found in "//\
                    inimg//sciext, l_logfile, verbose+)
                goto crash
            }
            print (keypar.value) | \
                scanf ("[%d:%d,%d:%d]", datax1, datax2, datay1, datay2)

            # Read the charaecteristic values from the database files
            #
            # Read the CCDNAME of the input extension
            keypar (inimg//sciext, "CCDNAME", silent+)
            if (!keypar.found) {
                printlog ("ERROR - GAPREPARE: CCDNAME keyword not \
                    found in "//outimg//sciext//". Exiting", \
                    l_logfile, verbose+)
                goto crash
            }
            ccdname = keypar.value

            # Read in the GAIN, READNOISE, SATURATION and nonln5 and nonln2
            # values. If there are changes to the arrays in the future add more
            # match statements to filter be keywords that identify the changes.
            # The * in the format string means skip this field
            # This reads everything except the READNOISE.
            match ("info", l_arraysdb, stop-, print_file_names=no, \
                metacharacters=no) | match (ccdname, "STDIN", stop-, \
                print_file_names=no, metacharacters=no) | \
                scanf ("%*s %*s %f %d %f %f\n", gainval, satlev, nonln5, \
                nonln2)

            # This reads the READNOISE.
            match (readmode, l_arraysdb, stop-, print_file_names=no, \
                metacharacters=no) | match (ccdname, "STDIN", stop-, \
                print_file_names=no, metacharacters=no) | \
                scanf ("%*s %*s %f\n", rdnval)

            # See if applying non-linear correction belwo 2000 ADU sepeartely
            if (do_low_nlc_sep) {
                # This reads the low NLC coefficients
                match (ccdname, l_non_lcdb, stop-, print_file_names=no, \
                    metacharacters=no) | match ("#", "STDIN", stop+, \
                    print_file_names=no, metacharacters=no) | \
                    scanf ("%*s %*s %g %g %g\n", nlla, nllb, nllc)

                high_range = "HIGH"
            } else {
                high_range = "ALL"
                nlla = INDEF
                nllb = INDEF
                nllc = INDEF
            }

            # This reads the NLC coefficients
            match (ccdname, l_non_lcdb, stop-, print_file_names=no, \
                metacharacters=no) | match (high_range, "STDIN", stop-, \
                print_file_names=no, metacharacters=no) | \
                match ("#", "STDIN", stop+, print_file_names=no, \
                metacharacters=no) | \
                scanf ("%*s %*s %g %g %g\n", nla, nlb, nlc)

            # Set the ROIID
            if (ii == 1) {

                # Read NAXIS1
                keypar (inimg//sciext, "i_naxis1", silent+)
                if (!keypar.found) {
                    printlog ("ERROR - GAPREPARE: NAXIS1 not found in "//\
                        inname//sciext, l_logfile, verbose+)
                    goto crash
                }
                naxis1 = int(keypar.value)

                if (ccdx1 != 1 && ccdx2 != 2048) {
                    # This is for array mode - ROIs in each extension that do
                    # not form a central square when mosaiced, e.g., 4 seperate
                    # stamp images
                    ##M Should this be a CCDSEC instead of naxis 1? If array
                    ##M mode can be centred in different places?
                    roiid = "A"//naxis1
                    detmode = "ARRAY"

                } else if (ccdx1 == 1 && ccdx2 == 2048) {
                    # Full frame
                    roiid = "FF"
                    detmode = "FULL FRAME"

                } else {
                    # Detector mode - a subregion of an extension that when
                    # mosaiced to forms a square centred on the centre of the
                    # array
                    roiid = "D"//naxis1
                    detmode = "DETECTOR"
                }

                printlog ("GAPREPARE: Image "//inname//" is in "//detmode//\
                    " mode.", l_logfile, l_verbose)

            } # End of setting ROIID value

            # Trim the data if required - calculate trimec (relavtive to
            # datasec) and update CCDSEC
            if (l_fl_trim) {

                # Use the istrimmed flag to determine if things need to be
                # updated in the header / what sentence gets printed to screen
                # It also helps with teh ODGW positions later on

                # This assumes no binning which is the case for GSAOI though
                # there is CCDSUM keyword. Should this change then this loop
                # will need to be updated - MS
                if (ccdx1 <= 4) {
                    datax1 += (5 - ccdx1)
                    ccdx1 = 5
                    istrimmed = yes
                }
                if (ccdx2 > 2044) {
                    datax2 -= (ccdx2 - 2044)
                    ccdx2 = 2044
                    istrimmed = yes
                }
                if (ccdy1 <= 4) {
                    datay1 += (5 - ccdy1)
                    ccdy1 = 5
                    istrimmed = yes
                }
                if (ccdy2 > 2044) {
                    datay2 -= (ccdy2 - 2044)
                    ccdy2 = 2044
                    istrimmed = yes
                }

                # The image doesn't need to be trimmed
                if (!istrimmed) {
                    printlog ("GAPREPARE: Image extension "//inimg//sciext//\
                        " does not require any trimming.", \
                        l_logfile, l_verbose)
                }

                # Update counter for checks later on to update PHU
                trimcount += 1

            } # End of trim loop; setting the trimsec and new ccdsec

            # Reset the CCDSEC
            ccdsec = "["//ccdx1//":"//ccdx2//","//ccdy1//":"//ccdy2//"]"

            # Set trimsec to extract data
            trimsec = "["//datax1//":"//datax2//","//datay1//":"//datay2//"]"

            ##M Should the detsec not be updated?

            # Copy the current input extensions trimsec (if not trimming this
            # will be the entire array
            imcopy (inimg//sciext//trimsec, outimg//sciappend, verbose-)

            # Update sciext here
            sciext = "["//l_sci_ext//","//ii//"]"

            # Update the CCDSEC
            gemhedit (outimg//sciext, "CCDSEC", ccdsec, "", delete-)

            # Update the TRIMSEC
            # Only update it if it doesn't exist and image has been trimmed
            # Thus if fl_trim+ and in array mode no trimsec written to header
            keypar (outimg//sciext, "TRIMSEC", silent+)
            if (!keypar.found && l_fl_trim && istrimmed) {
                gemhedit (outimg//sciext, "TRIMSEC", trimsec, \
                    "Trimmed section of input frame", delete-)
            }

            # Update the DATASEC
            datasec = "[1:"//(datax2 - (datax1 - 1))//",1:"//\
                (datay2 - (datay1 - 1))//"]"
            gemhedit (outimg//sciext, "DATASEC", datasec, "", delete-)

            # Check for the prepsence of NONLINEA and SATRUATI kewords
            # If not present, use numbers defined at the top of this file
            keypar (outimg//sciext, "NONLINEA", silent+)
            if (!keypar.found) {
                nonlin = int(nonln5 * satlev)

                printlog ("\nGAPREPARE: Adding NONLINEA keyword with value "//\
                    "of "//nonlin//" to "//outname//sciext, \
                    l_logfile, verbose-)

                gemhedit (outimg//sciext, "NONLINEA", nonlin, \
                    "5% Non-Linear before lin. correction (ADU)", delete-)
            }

            keypar (outimg//sciext, "SATURATI", silent+)
            if (!keypar.found) {
                printlog ("GAPREPARE: Adding SATURATI keyword with value "//\
                    "of "//satlev//" to "//outname//sciext, \
                    l_logfile, verbose-)

                gemhedit (outimg//sciext, "SATURATI", satlev,\
                    "Saturation before lin. correction (ADU)", delete-)
            }

            # Read and check (update if required) the gain keyword
            gainext = 0.
            keypar (outimg//"["//ii//"]", "GAIN", silent+)
            if (!keypar.found) {
                printlog ("WARNING - GAPREPARE: GAIN keyword not found in "//\
                    " in "//outimg//sciext//".\n", l_logfile, verbose+)
            } else {
                gainext = real(keypar.value)
            }

            if (gainext == 0.) {
                printlog ("WARNING - GAPREPARE: Updating GAIN keyword with "//\
                    "a value of "//gainval//" in "//outimg//sciext, \
                    l_logfile, verbose+)

                gemhedit (outimg//sciext,"GAIN", gainval, \
                    "Array gain (e-/ADU)", delete-)
            }

            # Read and check (update if needed) the read noise keyword
            rdnorig = 0.
            keypar (outimg//sciext, "RDNOISE", silent+)
            if (!keypar.found) {
                printlog ("WARNING - GAPREPARE: RDNOISE keyword not found "//\
                    " in "//outimg//sciext//".\n", l_logfile, verbose+)
            } else {
                rdnorig = real(keypar.value)
            }

            # If RDNOISE not in the headers, or 0. convert the values above to
            # e- using the values also given above for the gains - there is a
            # possible error here, but I'll leave to the user to figure it
            # out.

            if (rdnorig == 0.) {
                rdnorig = rdnval * gainval

                gemhedit (outimg//sciext, "RDNOISE", rdnorig, \
                    "Readout Noise in electrons", delete-)

            } # End of updating the rdnoise accordning to number of reads

            if (l_fl_nlc) {
                # Next is applying the linear correction - there should be one
                # per array again, made it simpler by hardcoding it here, best
                # solution is to read from a file this should use imexpr on
                # each extension I'll assume for the moment that the
                # non-linearity can be approximated by
                # ADU' = ADU*(a + b*ADU + c*ADU**2)
                # I'm not calling irlincor to avoid the normalization factor
                # (which is obsolete?)
                # NOTE - multiple coadds are AVERAGED (not summed) in the
                # output image.

                # Reset some of the variables
                inimg = outimg
                outimg = outname

                # Set PHU to read required information
                outphu = outimg//"[0]"

                if (ii == 1) {
                    # Copy out the phu again
                    fxcopy (inimg, outimg, groups="0", new_file=yes, verbose-)
                    print (outimg, >> tmpcrash)
                }

                # Set up the coeficients for ALL/HIGH ranges
                coeff1 = nla
                coeff2 = nlb
                coeff3 = nlc

                # Set up the expression low_nlc_val
                if (do_low_nlc_sep) {
                    nlexpr = "(a < "//low_nlc_val//") ? \
                       ((e + (f * a) + (g * a* a)) * a) : \
                       ((b + (c * a) + (d * a * a)) * a)"

                    lcoeff1 = nlla
                    lcoeff2 = nllb
                    lcoeff3 = nllc

                } else {
                    nlexpr = "(b + (c * a) + (d * a * a)) * a"

                    lcoeff1 = INDEF
                    lcoeff2 = INDEF
                    lcoeff3 = INDEF

                }

                imexpr (nlexpr, outimg//sciappend, inimg//sciext, \
                    coeff1, coeff2, coeff3, lcoeff1, lcoeff2, lcoeff3, \
                    dims="auto", intype="auto", \
                    outtype="real", refim="auto", bwidth=0, btype="nearest", \
                    rangecheck=yes, verbose=no, exprdb="none")

                # Calculate adjusted saturation value and non-linear level too
                xx = satlev
                satval = int((coeff1 + (coeff2 * xx) + \
                    (coeff3 * xx * xx)) * xx)
                nonlin = int(nonln2 * satval)

                # Update their header values too
                gemhedit (outimg//sciext, "NONLINEA", nonlin, \
                    "2% Non-Linear after lin. correction (ADU)", delete-)
                gemhedit (outimg//sciext, "SATURATI", satval, \
                    "Saturation after lin. correction (ADU)", delete-)

                # Update the nlccounter for check later
                nlccount += 1

            } # End of non-linear correction loop

            # Add the BUNIT keyword
            keypar (outimg//sciext, "BUNIT", silent+)
            if (!keypar.found) {
                # This assumes only raw data comes through this point - MS
                gemhedit (outimg//sciext, "BUNIT", "ADU", \
                    "Units of pixel data", delete-)
            }

            # Add the VAR and DQ planes if desired. This is done after the
            # linearity correction as it may be desirable to include non-linear
            # pixels in the DQ uses imcopy from a static BPM to create the
            # initial DQ, then imexpr to flag the new
            # bad pixels, then combine the static BPM with the floating one
            # bad pixels are set to 1,2,4, good pixels are set to 0
            # but note that processing down the line will lose that info, and
            # set all to 1

            if (l_fl_vardq) {

                varappend = "["//l_var_ext//","//ii//",append]"
                varextn = "["//l_var_ext//","//ii//"]"
                dqappend = "["//l_dq_ext//","//ii//",append]"
                dqextn = "["//l_dq_ext//","//ii//"]"

                # Calculate the variance plane -- Observational Astronomy
                # D. Scott Birney, Guillermo González, David Oesper page 178:
                # var = sigma^2 = gain*ADU + readnoise^2, where the readnoise
                # is given in electrons which is the case for GSAOI
                # a = sci extension, b = read noise, c = gain
                # Data are in ADU here; Variance should be science units
                # squared (ADU**2). so:
                # variance = (pixel_value[ADU] / gain) + (readnoise[e-]/gain)**2

                # Set the variance expression
                varexpr = "((max(a,0.0) / c) + (b / c)**2)"

                # The following two keywords have be previously checked for
                keypar (outimg//sciext, "GAIN", silent+)
                gain = real(keypar.value)

                keypar (outimg//sciext, "RDNOISE", silent+)
                rdnoise = real(keypar.value)

                # Correct READNOISE by coadds rn = rn / sqrt(coadds) as they
                # are averaged. Readnoise is in electrons, converted in VAR
                # expression to ADU
                rdcor = (rdnoise / sqrt(coadds))

                # Create the variance plane
                imexpr (varexpr, outimg//varappend, outimg//sciext, \
                    rdcor, gain, dims="auto", intype="auto", outtype="ref", \
                    refim="a", bwidth=0, btype="nearest", bpixval=0., \
                    rangecheck=yes, verbose=no, exprdb="none")

                # Update the title of this extension
                gemhedit (outimg//varextn, "i_title", "Variance", "", \
                    delete-)

                # Add the BUNIT keyword
                keypar (outimg//varextn, "BUNIT", silent+)
                # This assumes only raw data comes through this point - MS
                gemhedit (outimg//varextn, "BUNIT", "ADU*ADU", \
                    "Units of pixel data", delete-)

                # Update the GAIN value of the VAR extension
                keypar (outimg//varextn, "GAIN", silent+)
                gemhedit (outimg//varextn, "GAIN", \
                    (real(keypar.value) * real(keypar.value)), \
                    "Gain*Gain (e-/ADU)**2", delete-)

                # Calculate the DQ plane
                # First check if the static bpm exists
                # then check if the sat/non-lin are to be flagged
                # saturated = 4, non-linear = 2, bad pixels = 1
                # if both the above are no, DQ is all 0 unless ODGW present
                # if first = no, but second = yes, calculate DQ from sci
                # if first = yes and second = yes, combine the two
                # if first yes and second no, just copy the static in the DQ

                # Remember that some tasks use bad value=1, need to use good
                # value=0, so they flag saturated/non-linear as well.

                # Read / parse all required values for DQ plane
                if (l_fl_sat) {
                    # Read the saturation and non-linear values - they will
                    # exist by this point
                    keypar (outimg//sciext, "NONLINEA", silent+)
                    nonlin = real(keypar.value)

                    keypar (outimg//sciext, "SATURATIO", silent+)
                    satval = real(keypar.value)
                } else {
                    nonlin = INDEF
                    satval = INDEF
                } # End of l_fl_sat DQ loop

                # Check that the static bpm covers the sci extensions
                if (l_badpix != "none") {

                    keypar (l_badpix//"[0]", "NEXTEND", silent+)
                    if (!keypar.found) {
                        printlog ("ERROR - GAPREPARE: NEXTEND keyword not \
                            found in "//l_badpix//"[0]. Exiting", \
                            l_logfile, verbose+)
                        goto crash
                    }
                    bpmnext = int(keypar.value)

                    # Loop over BPM extensions and match ccdname
                    for (jj = 1; jj <= bpmnext; jj += 1) {
                        keypar (l_badpix//"["//jj//"]", "CCDNAME", silent+)

                        if (!keypar.found) {
                            printlog ("ERROR - GAPREPARE: CCDNAME keyword \
                                not found in "//outimg//sciext//". Exiting", \
                                l_logfile, verbose+)
                            goto crash

                        } else if (keypar.value == ccdname) {

                            bpmext = str(jj)

                            keypar (l_badpix//"["//jj//"]", "EXTNAME", silent+)
                            if (keypar.found && keypar.value != "") {
                                bpmext = keypar.value

                                keypar (l_badpix//"["//jj//"]", \
                                    "EXTVER", silent+)
                                if (keypar.found && str(keypar.value) != "") {
                                    bpmext = bpmext//","//str(keypar.value)
                                } else {
                                    bpmext = str(jj)
                                }
                            }

                            break
                        } else if (jj == bpmnext) {
                            # CCDNAME not matched at all
                            printlog ("ERROR - GAPREPARE: could not match \
                                CCDNAME of "//ccdname//" from "//\
                                outimg//sciext//" in "//l_badpix//\
                                ". Exiting.", \
                                l_logfile, verbose+)
                            goto crash
                        }
                    }

                    # Call gadimschk to return the correct section of the BPM
                    gadimschk (outimg//sciext, section="", chkimage=l_badpix//\
                        "["//bpmext//"]", key_check="CCDSEC", \
                        logfile=l_logfile, verbose=no)

                    if (gadimschk.status != 0) {
                        printlog ("\nERROR - GAPREPARE: GADIMSCHK returned \
                            a non-zero status."//\
                            "\n                   Possible ROI mismatch \
                            between static BPM "//l_badpix//"["//bpmext//"]"//\
                            " and "//outimg//sciext, \
                            l_logfile, verbose+)
                        goto crash
                    } else {
                        bpmexp = gadimschk.out_chkimage
                    }

                } else {
                    bpmexp = INDEF
                } # End of BPM checks

                # Parse the coordinates for ODGWs - all information has been
                # read already
                ##M Possible that for ROIs parked ODGW will stay in corner of
                ##M array not ROI
                ##M Also might be possible that ODGW can guide outside science
                ##M ROI!
                if (fl_odgw) {
                    # Test if this extensions guide window was guiding

                    if (debug) {
                        printlog ("\nODGW Size: "//odgws, l_logfile, verbose+)
                        printlog ("odgwx: "//odgwx[ii]//" odgwy: "//odgwy[ii],\
                            l_logfile, verbose+)
                    }

                    if (odgwx[ii] != 0 && odgwy[ii] != 0) {
                        # Guiding
                        # This should account for trimming / ROIs
                        # Assumes ODGW coordinates in header are CCD coordinate
                        # frame

                        ##M May be able to guide with ODGW outside science ROI
                        ##M so this may need to be handled here somehow.
                        odgw_x1 = odgwx[ii] - (odgws - 1) - (ccdx1 - 1)
                        odgw_x2 = odgwx[ii] + odgws - (ccdx1 - 1)
                        odgw_y1 = odgwy[ii] - (odgws - 1) - (ccdy1 - 1)
                        odgw_y2 = odgwy[ii] + odgws - (ccdy1 - 1)

                        # Catch any out of bounds numbers - This should never
                        # happen except if the ODGW has some area in the
                        # trimmed section - MS
                        odgw_x1 = max (1, odgw_x1)
                        odgw_x2 = min ((ccdx2 - ccdx1 + 1), odgw_x2)
                        odgw_y1 = max (1, odgw_y1)
                        odgw_y2 = min ((ccdy2 - ccdy1 + 1), odgw_y2)

                    } else {
                        # Assume not guiding - use corner furthest away from
                        # the centre of the mosaiced array
                        ##M - This may not be true for ROIs. Need to check in
                        ##M the future. -- Could be the case where they are
                        ##M parked in the corner but there is still some
                        ##M overlap

                        # Parse the DETSEC for this extension and determine
                        # which corner it's in.
                        keypar (outimg//sciext, "DETSEC", silent+)
                        if (!keypar.found) {
                            printlog ("ERROR - GAPREPARE: DETSEC keyword not \
                                found in "//outimg//sciext//". Exiting", \
                                l_logfile, verbose+)
                            goto crash
                        }
                        print (keypar.value) | scanf ("[%d:%d,%d:%d]", detx1, \
                            detx2, dety1, dety2)

                        # Parse the datasec
                        keypar (outimg//sciext, "DATASEC", silent+)
                        if (!keypar.found) {
                            printlog ("ERROR - GAPREPARE: DATASEC keyword not \
                                found in "//outimg//sciext//". Exiting", \
                                l_logfile, verbose+)
                            goto crash
                        }
                        print (keypar.value) | scanf ("[%d:%d,%d:%d]", datax1,\
                            datax2, datay1, datay2)

                        # Appears that the ODGW are not quite where first
                        # thought when parked. They are parked in the corners
                        # but should be outside of the trimmed area. However,
                        # it appears that with the exception of extension 2
                        # they overlap by one pixel one or two edges. All
                        # parked ODGW have been shifted by 1 pixel in both x
                        # any y relative to the bottom left corner of the
                        # entire array. Assume this is the case for all modes
                        # - MS

                        # Calculate the trim to add on if not trimmed
                        if (ccdx1 < 5) {
                            trimx1 = 5 - ccdx1
                        } else {
                            trimx1 = 0
                        }

                        if (ccdx2 > 2044) {
                            trimx2 = ccdx2 - 2044
                        } else {
                            trimx2 = 0
                        }

                        if (ccdy1 < 5) {
                            trimy1 = 5 - ccdy1
                        } else {
                            trimy1 = 0
                        }

                        if (ccdy2 > 2044) {
                            trimy2 = ccdy2 - 2044
                        } else {
                            trimy2 = 0
                        }

                        # The extra addition of (shift/pad) is to compensate
                        # for the shifted ODGW positions when parked. If only
                        # masking the ODGW and not including the additional
                        # rows to the for extensions 2 and 3, and below for
                        # extensions 1 and 2, then add 1 (shift) to all
                        # values (see *** below). If including these extra
                        # rows use pad in the if loops below - MS
                        pad = 1
                        if (detx1 > 2048) {
                            odgw_x1 = datax2 - (((2 * odgws) - 1) + trimx2)
                            odgw_x2 = datax2 - trimx2
                        } else {
                            odgw_x1 = datax1 + trimx1 - pad
                            odgw_x2 = datax1 + (((2 * odgws) - 1) + trimx1)
                        }

                        if (dety1 > 2048) {
                            odgw_y1 = datay2 - (((2 * odgws) - 1) + trimy2)
                            odgw_y2 = datay2 - trimy2
                        } else {
                            odgw_y1 = datay1 + trimy1 - pad
                            odgw_y2 = datay1 + (((2 * odgws) - 1) + trimy1)
                        }

                        #***
                        shift = 1
                        odgw_x1 += shift
                        odgw_x2 += shift
                        odgw_y1 += shift
                        odgw_y2 += shift

                        # Catch any values above the limits of the array
                        # Possible due to shift of ODGW - MS
                        odgw_x1 = max (datax1, odgw_x1)
                        odgw_x2 = min (datax2, odgw_x2)
                        odgw_y1 = max (datay1, odgw_y1)
                        odgw_y2 = min (datay2, odgw_y2)

                        if (debug) {
                            printlog ("datax1: "//datax1//", datax2: "//\
                                datax2//", datay1: "//datay1//", datay2: "//\
                                datay2, l_logfile, verbose+)
                        }

                    } # End of determining ODGW positions
                } else {
                    odgw_x1 = INDEF
                    odgw_x2 = INDEF
                    odgw_y1 = INDEF
                    odgw_y2 = INDEF
                } # End of l_fl_odgw loop

                # Set the DQ expression - only need to do this for the first
                # extension - dqexpr will not change after this extension for
                # this image
                if (ii == 1) {
                    # a = output sci extension (outimg//sciext),
                    # b = saturation level (satval?),
                    # c = non-linear level (nonlin?),
                    # d = bpm (bpmexp),
                    # e - h = odgw coordinates

                    # Set up default dqexpr - this is to keep the header
                    # information from the SCI extension always - the int is to
                    # allow bitwise or's
                    dqexpr = "(int(a * 0))"

                    # Use section of BPM if supplied
                    if (l_badpix != "none") {
                        dqexpr = dqexpr // " | d"
                        printlog ("\nGAPREPARE: Adding static BPM to DQ \
                            plane", l_logfile, l_verbose)
                    } else {
                        printlog ("GAPREPARE: No static BPM to be included in \
                            DQ plane", l_logfile, l_verbose)
                    }

                    # Set the saturation / non-linear flagging expression
                    if (l_fl_sat) {
                        printlog ("GAPREPARE: Flagging saturated and \
                            non-linear pixels in DQ plane", \
                            l_logfile, l_verbose)

                        dqexpr = dqexpr // " | ((a > b) ? 4 : (a > c) ? 2 : 0)"

                    } else {
                        printlog ("GAPREPARE: Not flagging saturated and \
                            non-linear pixel in DQ plane", \
                            l_logfile, l_verbose)
                    }

                    # Flag all ODGW (including parked ones) with a 1
                    if (fl_odgw) {
                        printlog ("GAPREPARE: Flagging ODGW areas with a \
                            value of 16 in the DQ plane", l_logfile, l_verbose)

                        dqexpr = dqexpr // \
                            " | ((I >= e && I <= f && J >= g && J <= h) ? \
                            16 : 0)"
                    } else {
                        printlog ("GAPREPARE: No ODGW areas to flag in \
                            DQ plane", l_logfile, l_verbose)
                    }
                } # End of determining DQ expression

                if (debug) {
                    printlog ("\nDQEXPR: \""//dqexpr//"\"", \
                        l_logfile, verbose+)
                    printlog ("Operands:\n", \
                        l_logfile, verbose+)
                    printlog ("    a = "//outimg//sciext, \
                        l_logfile, verbose+)

                    if (isindef(satval)) {
                        op_toprint = "INDEF"
                    } else {
                        op_toprint = str(satval)
                    }
                    printlog ("    b = "//op_toprint, \
                        l_logfile, verbose+)

                    if (isindef(nonlin)) {
                        op_toprint = "INDEF"
                    } else {
                        op_toprint = str(nonlin)
                    }
                    printlog ("    c = "//op_toprint, \
                        l_logfile, verbose+)

                    if (isindef(bpmexp)) {
                        op_toprint = "INDEF"
                    } else {
                        op_toprint = bpmexp
                    }
                    printlog ("    d = "//op_toprint, \
                        l_logfile, verbose+)

                    if (isindef(odgw_x1)) {
                        op_toprint = "INDEF"
                    } else {
                        op_toprint = str(odgw_x1)
                    }
                    printlog ("    e = "//op_toprint, \
                        l_logfile, verbose+)

                    if (isindef(odgw_x2)) {
                        op_toprint = "INDEF"
                    } else {
                        op_toprint = str(odgw_x2)
                    }
                    printlog ("    f = "//op_toprint, \
                        l_logfile, verbose+)

                    if (isindef(odgw_y1)) {
                        op_toprint = "INDEF"
                    } else {
                        op_toprint = str(odgw_y1)
                    }
                    printlog ("    g = "//op_toprint, \
                        l_logfile, verbose+)

                    if (isindef(odgw_y2)) {
                        op_toprint = "INDEF"
                    } else {
                        op_toprint = str(odgw_y2)
                    }
                    printlog ("    h = "//op_toprint//"\n", \
                        l_logfile, verbose+)
                } # End of debug

                # Create output DQ plane
                imexpr (dqexpr, outimg//dqappend, outimg//sciext, satval, \
                    nonlin, bpmexp, odgw_x1, odgw_x2, odgw_y1, odgw_y2,
                    dims="auto", intype="auto", outtype="ushort", \
                    refim="auto", bwidth=0, btype="nearest", bpixval=0., \
                    rangecheck=yes, verbose=no, exprdb="none")

                # Update the title of this extension
                gemhedit (outimg//dqextn, "i_title", "Data Quality", "", \
                    delete-)

                # Delete the BUNIT keyword
                gemhedit (outimg//dqextn, "BUNIT", "", "", delete+)

                # If BPM used add to name and region used to DQ plane
                if (!isindef(bpmexp)) {
                    # Parse the BPM to remove any directory names
                    fparse (bpmexp)
                    bpm_toprint = substr(bpmexp, \
                        (strstr(fparse.directory,bpmexp) + \
                            strlen(fparse.directory)), strlen(bpmexp))
                    gemhedit (outimg//dqextn, "BPMASK", \
                        bpm_toprint, "Static BPM used", delete-)
                }

            } # End of vardq loop

        } # End of loop over extensions

        printlog ("", l_logfile, l_verbose)

        # Update phu if image was trimmed
        if (trimcount == nextend && l_fl_trim) {
            # This means that if in array mode - the PHU will be updated with
            # the TRIMMED keyword if fl_trim+ but the science extensions will
            # not contain a TRIMSEC keyword
            gemhedit (outphu,"TRIMMED", "yes", \
                "Image has been trimmed", delete-)

            printlog ("GAPREPARE: Image "//outimg//" trimmed", \
                l_logfile, l_verbose)

        } else if (trimcount != nextend && l_fl_trim) {
            printlog ("ERROR - GAPREPARE: Only some extensions of \
                Image "//outimg//" have bene trimmed. Exiting", \
                l_logfile, l_verbose)
            goto crash

        } else if (!l_fl_trim && !alltrimmed) {
            printlog ("GAPREPARE: Image "//outimg//" NOT trimmed",\
                l_logfile, l_verbose)
        }

        # Update PHU if non-linear corrected
        if (nlccount == nextend && l_fl_nlc) {
            gemhedit (outphu, "LINRCORR", "yes", \
                "Linearity correction applied", delete-)

            printlog ("GAPREPARE: Image "//outimg//" corrected for \
                non-linearity",l_logfile, l_verbose)

        } else if (nlccount != nextend && l_fl_nlc) {
            printlog ("ERROR - GAPREPARE: Only some extensions of \
                Image "//outimg//" have bene non-linear corrected. Exiting", \
                l_logfile, l_verbose)

            goto crash
        } else if (!l_fl_nlc && !allnlc) {
            printlog ("GAPREPARE: Image "//outimg//" NOT corrected \
                for non-linearity", l_logfile, l_verbose)
        }

        # Check validity of VAR and DQ planes
        if (l_fl_vardq) {
            gemextn (outimg, check="exists", process="expand", index="", \
                extname=l_var_ext//","//l_dq_ext, extversion="1-"//nextend, \
                ikparams="", omit="", replace="", outfile="dev$null", \
                logfile=l_logfile, glogpars="", verbose=no)

            if (gemextn.status != 0) {
                printlog ("ERROR- GAPREPARE: GEMEXTN returned a non-zero \
                    status. Exiting.", l_logfile, verbose+)
                goto crash
            } else {
                if ((gemextn.count != (2 * nextend)) && (gemextn.count > 0)) {
                    printlog ("ERROR - GAPREPARE: Image "//outimg//\
                        " does not contain the correct number of VAR and DQ \
                        extensions. Exiting.", l_logfile, verbose+)
                    goto crash
                } else {
                    printlog ("GAPREPARE: VAR/DQ planes added to image "//\
                        outimg, l_logfile, l_verbose)
                }
            }
        } else {
            printlog ("GAPREPARE: VAR/DQ planes NOT added to image "//outimg, \
                l_logfile, l_verbose)
        }

        printlog ("", l_logfile, l_verbose)

        # Finally, add two keywords: METACONF and GAPREPAR

        # - METACONF, which is a metaconfiguration keyword to identify
        # the full configuration used (filter, plus any other mechanism
        # position / parameter relevant for data reduction and calibration
        # association)
        #     for objects:
        #         FILTER1_FILTER2+EXPTIME_LNRS_NCOADDS+'OBJ'+ROIID(+TRIM+NLC)
        #     for flats:
        #         FILTER1_FILTER2+EXPTIME_NCOADDS_LNRS+TYPE(+GCALSHUT)//\
        #         +ROIID(+TRIM+NLC)
        #     for darks: EXPTIME_NCOADDS_LNRS+ROIID(+TRIM+NLC)
        # Note the use of "+" rather than "_" to be easier to distinguish
        # between metaconf elements

        # Initiate metaconfig keyword
        metakey = ""

        # Read OBSTYPE
        keypar (outphu,"OBSTYPE", silent+)
        if (!keypar.found) {
            printlog ("ERROR - GAPREPARE: OBSTYPE keyword not found in "//\
                outphu, l_logfile, verbose+)
            goto crash
        }
        obskey = keypar.value

        # Read FILTER1
        keypar (outphu, "FILTER1", silent+)
        if (!keypar.found) {
            printlog ("ERROR - GAPREPARE: FILTER1 keyword not found in "//\
                outphu, l_logfile, verbose+)
            goto crash
        }
        filt1key = keypar.value

        # Read FILTER2
        keypar (outphu, "FILTER2", silent+)
        if (!keypar.found) {
            printlog ("ERROR - GAPREPARE: FILTER2 keyword not found in "//\
                outphu, l_logfile, verbose+)
            goto crash
        }
        filt2key = keypar.value

        # Read EXPTIME
        keypar (outphu, "EXPTIME", silent+)
        if (!keypar.found) {
            printlog ("ERROR - GAPREPARE: EXPTIME keyword not found in "//\
                outphu, l_logfile, verbose+)
            goto crash
        }
        itimekey = keypar.value

        # Due to differences bewteen PyRAF anf CL numbers that end in .0 may or
        # may not have that ".0" in the exposure time metaconfiguration section
        if (stridx(".",itimekey) == strlen(itimekey)) {
            itimekey = itimekey//"0"
        }

        # Read COADDS
        keypar (outphu,"COADDS", silent+)
        if (!keypar.found) {
            printlog ("ERROR - GAPREPARE: COADDS keyword not found in "//\
                outphu, l_logfile, verbose+)
            goto crash
        }
        coaddkey = keypar.value

        # Read LNRS
        keypar (outphu, "LNRS", silent+)
        if (!keypar.found) {
            printlog ("ERROR - GAPREPARE: LNRS keyword not found in "//\
                outphu, l_logfile, verbose+)
            goto crash
        }
        lnrskey = keypar.value

        # Read OBJECT
        keypar (outphu, "OBJECT", silent+)
        if (!keypar.found) {
            printlog ("ERROR - GAPREPARE: OBJECT keyword not found in "//\
                outphu, l_logfile, verbose+)
            goto crash
        } else {
            object_raw = keypar.value
            if (object_raw == "Twilight") {
                objct = "TWLT+"
            } else if (strstr("Domeflat",object_raw) > 0) {
                # Updated to always have OPEN in the string due to the
                # possibility of having OFF flats too
                objct = "DOME+"
                # Check for dome OFF flats
                if (strstr("OFF",object_raw) > 0) {
                    objct = objct//"CLOS+"
                } else {
                    objct = objct//"OPEN+"
                }

            } else {
                # Although OBJECT can be many other things than a true object
                # and those listed above, because the METACONF keyword is setup
                # below in the way it is can set objct = to "OBJ+" here - MS
                objct = "OBJ+"
            }
        }

        filt_string = filt1key//"_"//filt2key//"+"
        time_string = itimekey//"_"//lnrskey//"_"//coaddkey//"+"

        # Form the METACONF keyword
        if (obskey == "OBJECT") {
            # Science data, twilights and dome flats

            metakey = filt_string//objct

        } else if (obskey == "FLAT") {
            # Only == FLAT if GCAL

            # Read GCALSHUT (only present if GCAL FLAT)
            keypar (outphu, "GCALSHUT", silent+)
            if (!keypar.found && (objct == "GCALflat")) {
                printlog ("ERROR - GAPREPARE: GCALSHUT keyword not found \
                    in "//outphu, l_logfile, verbose+)
                goto crash
            }
            gcalkey = substr(keypar.value,1,4)

            metakey = filt_string//"GCAL+"//gcalkey//"+"

        } else if (obskey != "DARK") {

            metakey = "UNKNOWN"
        }

        if (metakey != "UNKNOWN") {

            metakey = time_string//metakey//roiid

            keypar (outphu, "TRIMMED", silent+)
            if (keypar.found) {
                metakey = metakey//"+TRIM"
            }

            keypar (outphu, "LINRCORR", silent+)
            if (keypar.found) {
                metakey = metakey//"+NLC"
            }
        }

        printlog ("GAPREPARE: Metaconfiguration for "//outimg//" is:\n"//\
            "           "//metakey, l_logfile, l_verbose)

        printlog ("\nGAPREPARE: The frame "//outimg//" has OBSTYPE="//obskey//\
            " and OBJECT="//object_raw, l_logfile, l_verbose)

        gemhedit (outphu, "METACONF", metakey, \
            "Data Reduction Metadata", delete-)

        # Add number of science extensions to phu
        gemhedit (outphu, "NSCIEXT", nextend, "Number of science extensions", \
            delete-)

        # Update the number of extensions
        gemextn (outimg, check="", process="expand", index="1-", extname="", \
            extver="", ikparam="", omit="", replace="", outfile="dev$null", \
            logfile=l_logfile, glogpars="", verbose=no)

        if (gemextn.status != 0) {
            printlog ("ERROR- GAPREPARE: GEMEXTN returned a non-zero \
                status. Exiting.", l_logfile, verbose+)
            goto crash
        } else {
            gemhedit (outphu, "NEXTEND", gemextn.count, "", delete-)
        }

        # Add GAPREPAR with the UT time stamp to indicate when/if the file has
        # been prepared
        gemdate()
        gemhedit (outphu, "GAPREPAR", gemdate.outdate, \
            "UT Time stamp for GAPREPARE", delete-)

        # Update last time (editted/accessed?) by gemini software
        gemhedit (outphu, "GEM-TLM", gemdate.outdate, \
            "UT Last modification with GEMINI", delete-)

        num_created += 1

NEXTIMAGE:

        printlog ("----", l_logfile, l_verbose)

        # Remove the temporary file
        if (l_fl_nlc && imaccess(outimg_orig)) {
            imdelete (outimg_orig, verify-, >& "dev$null")
        }

        # Remove the tmpcrash list file
        if (access(tmpcrash)) {
            delete (tmpcrash, verify-, >& "dev$null")
        }

    } # End of loop over images

goto clean

#--------------------------------------------------------------------------
crash:
    # Exit with error subroutine
    status = 1

    # Remove the tmpcrash list file and any files it contains
    if (access(tmpcrash)) {
        delete ("@"//tmpcrash, verify-, >& "dev$null")
        delete (tmpcrash, verify-, >& "dev$null")
    }
    goto clean

clean:

    # Remove the file containinig the input list
    delete (filelist, verify-, >& "dev$null")

    scanfile = ""

    if (access(outgaimchk_list)) {
        delete ("@"//outgaimchk_list//","//outgaimchk_list, verify-, \
            >& "dev$null")
    }

    # Check for files that were not processed
    if (access(tmpnotproc)) {
        printlog ("\nWARNING - GAPREPARE: The files listed below were skipped \
            and therefore not processed.\n", l_logfile, verbose+)

        # Print the files to the screen
        scanfile = tmpnotproc
        while (fscan(scanfile, inimg) != EOF) {
            printlog ("    "//inimg, l_logfile, verbose+)
        }

        delete (tmpnotproc, verify-, >& "dev$null")
        scanfile = ""
        if (num_created == 0) {
            printlog ("ERROR - GAPREPARE: No out files were created", \
                l_logfile, verbose+)
            status = 1
        }
    }

    # Print finish time
    gemdate(zone="local")
    printlog ("\nGAPREPARE -- Finished: "//gemdate.outdate, \
        l_logfile, l_verbose)

    if (status == 0) {
        printlog ("\nGAPREPARE -- Exit status: GOOD", \
            l_logfile, verbose=l_verbose)
    } else {
        printlog ("\nGAPREPARE -- Exit status: ERROR", \
            l_logfile, verbose+)
    }

    if (status != 0) {
        printlog ("          -- Please read the logfile \""//l_logfile//\
                "\" for more information.", l_logfile, verbose+)
    }

    printlog ("------------------------------------------------------------\n",
        l_logfile, l_verbose)

end

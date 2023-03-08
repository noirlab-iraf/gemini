# Copyright(c) 2010-2016 Association of Universities for Research in Astronomy, Inc.

procedure gaflat (inimages)

##M ignore_exp - this potentially has a bug when mixing different observing
##M modes and coadds - maybe it should just exclude up to and including the first
##M "_"

# TODO
# - optimize parameters for combine for the various types - default option

# INPUTS:
#     GSAOI raw or gaprepare'd flat data
#     This task is smart and will find all types flats in the given location if
#     *.fits is supplied. See gaimchk for more information.

# This task will:
# - combine all flats for same filter/(if requested exposure time too)
# - combine all OFF flats for those required to have OFF flats same filter/(if
#       requested exposure time too)
# - subtract ON - OFF; if required
# - normalize using the options available from gastat
# - combine twilight/dome flats at this point without previous dark
#   subtraction. There is no warning!

# IT DOES NOT WORK TO CREATE SKY FLATS (it requires either GCALSHUT=OPEN
# or OBJECT=Twilight/Domeflat - use gasky instead and normalize with gemexpr

# The output frame inherit the appended specified prefix from GAPREPARE and an
# _flat is added at the end (name is that of the first image on the stack)
# a keyword GAFLAT is added with the time stamp when the task is last run

# CAVEATS:
# - Flats off are MANDATORY for GCAL
# - default for ignore_exp is only true for twilight flights, other types must
#   have matching exposure times

# Version 2010 May 20 V0.0 CW - created
#         2010 Jun 24 V1.0 CW - GCAL flats only
#         2010 Jul 05 V1.1 CW - Added basic sky flat (rejection only, no
#                               masking)
#         2011 Mar 04 V1.2 CW - removed sky flats until dark subtraction can be
#                               implemented; removed flats off option
#                               (mandatory)
#         2011 Apr 15 V1.3 CW - implemented option to normalize by fit1d
#                               smoothing
#         2011 Apr 21 V1.4 CW - twilight flat option added. Does not require
#                               darks, if dark subtraction require will have \
#                               to call gareduce due to the different exptimes
#         2011 Oct 06 V1.5 CW - handling of ROIs
#         2011 Dec 07 V1.6 CW - modify metaconf to handle ROIs
#         2012 Jan 23 V1.7 CW - include Domeflats as we don't know what the
#                               problem is with twilights vs GCAL.
#
# Descriptions of future revisions are recorded in CVS logs

######################################################
# input, output and selection parameters
char    inimages    {prompt="Input GSAOI images"}
char    rawpath     {"",prompt="Path for raw input images"}
char    outsufx     {"flat",prompt="Suffix for output flat"}
char    rootname    {"",prompt="Root name for images; blank=today UT"}
int     minflat     {5,min=2,prompt="Minimum number of flat images to combine"}
char    ignore_exp  {"default",enum="yes|no|default",prompt="Ignore exposure times when combining"}
#bool    fl_scale    {no,prompt="Scale images before combining"}
char    stattype    {"mean",enum="mean|mode|midpt",prompt="Type of statistics to compute for normalization"}
char    statextn    {"ARRAY",prompt="How to calculate normalization. (DETECTOR|ARRAY|<extnname,extnver>|<index>)"}
char    statsec     {"[*,*]",prompt="Statistics section for scaling when combing (relative to an array)"}
bool    fl_mask     {yes,prompt="Mask non-good pixels when calculating normalization statistics?"}
char    badpix      {"gsaoi$data/gsaoibpm_high_full.fits", prompt="Static Bad Pixel Mask - not mosaic"}
bool    fl_vardq    {no,prompt="Create variance and data quality frames"}
real    maxtime     {INDEF,prompt="Maximum time interval from first image in the list"}
char    datename    {"DATE-OBS",prompt="Date header keyword"}
char    timename    {"UT",prompt="Time stamp header keyword"}
bool    ignore_nlc  {no,prompt="Ignore NLC state of the input files?"}
char    use_offs    {"default",enum="default|yes|no",prompt="Use off flats?"}
# gaprepare parameters
char    gaprep_pref {"g",prompt="Prefix for GAPREPARE output images"}
bool    fl_trim     {yes,prompt="Trim the images?"}
bool    fl_nlc      {yes,prompt="Apply non-linear correction to each array?"}
bool    fl_sat      {yes,prompt="Include non-linear and saturated pixels in data quality planes"}
char    arraysdb    {"gsaoi$data/gsaoiAMPS.dat",prompt="Database file for characteristics of each array"}
char    non_lcdb    {"gsaoi$data/gsaoiNLC.dat",prompt="Database file for non-linearity correction coefficients"}
# gemcombine parameters
char    combine     {"default",enum="default|average|median",prompt="Combination operation"}
char    reject      {"avsigclip",enum="none|minmax|ccdclip|crreject|sigclip|avsigclip|pclip",prompt="Rejection algorithm"}
char    masktype    {"goodvalue",enum="none|goodvalue",prompt="Bad Pixel Mask type"}
real    maskvalue   {0.,prompt="Good pixel value in the BPM"}
char    zero        {"none",prompt="Image zero point offset (none|mode|median|mean|@<file>|!<keyword>)"}
char    weight      {"none",prompt="Image weights (none|mode|median|mean|exposure|@<file>|!<keyword>"}
char    expname     {"EXPTIME",prompt="Exposure time header keyword"}
real    lthreshold  {INDEF,prompt="Lower threshold"}
real    hthreshold  {INDEF,prompt="Upper threshold"}
int     nlow        {1,min=0,prompt="minmax: Number of low pixels to reject"}
int     nhigh       {1,min=0,prompt="minmax: Number of high pixels to reject"}
int     nkeep       {1,prompt="Minimum to keep or maximum to reject"}
bool    mclip       {yes,prompt="Use median in sigma clipping algorithms?"}
real    lsigma      {3.,min=0,prompt="Lower sigma clipping factor"}
real    hsigma      {3.,min=0,prompt="Upper sigma clipping factor"}
char    key_ron     {"RDNOISE",prompt="Keyword for readout noise in e-"}
char    key_gain    {"GAIN",prompt="Keyword for gain in electrons/ADU"}
real    ron         {0.0,min=0.,prompt="Readout noise rms in electrons"}
real    gain        {1.0,min=0.00001,prompt="Gain in e-/ADU"}
char    snoise      {"0.0",prompt="ccdclip: Sensitivity noise (electrons)"}
real    sigscale    {0.1,min=0.,prompt="Tolerance for sigma clipping scaling correction"}
real    pclip       {-0.5,prompt="pclip: Percentile clipping parameter"}
real    grow        {0.0,min=0.,prompt="Radius (pixels) for neighbor rejection"}
bool    fl_dqprop   {no,prompt="Propagate all DQ values?"}
# other stuff
char    sci_ext     {"SCI",prompt="Name of science extensions"}
char    var_ext     {"VAR",prompt="Name of variance extensions"}
char    dq_ext      {"DQ",prompt="Name of data quality extensions"}
char    logfile     {"",prompt="Logfile"}
bool    verbose     {yes,prompt="Verbose?"}
int     status      {0,prompt="Exit status (0=good)"}
struct  *scanfile   {"",prompt="Internal use only"}
struct  *scanfile2  {"",prompt="Internal use only"}
struct  *scanfile3  {"",prompt="Internal use only"}

######################################################

begin

    ########
    # Declare and set local variables; set default values; initiate temporary
    # files

    ####
    # Variable declaration
    char l_inimages, l_rawpath, l_outsufx, l_logfile, l_rootname
    char l_badpix, l_combine, l_reject, l_masktype, metaconf, flatname
    char l_scale, l_zero, l_weight, l_statsec, l_expname, inimg
    char l_key_ron, l_key_gain, l_snoise, l_ignore_exp, filelist, utdate
    char l_datename, l_timename, print_ext, current_type
    char inlist, t_string, image_root, conflist, inphu
    char tmpconflist, l_sci_ext, l_var_ext, l_dq_ext, inimg_root, inname
    char l_gaprep_pref, off_metaconf, off_conflist, tmp_skipped, tmpnotproc
    char tmpgaprep, todelete, tmplist, nextlist, sub_flat, refname, tmpinsert
    char onlist, offlist, oncomb_list, offcomb_list, out_oncomb_list
    char out_offcomb_list, l_title, comblist, tmpflat, norm_list, inflat
    char out_flat, imtmptodel, tmpnorm_list, onflat, off_flat, pr_string
    char on_uttime, on_datevalue, off_uttime, off_datevalue, closest_off
    char stat_file, dq_file, tile_list, dqtile_list, l_stattype, scale_list
    char tmponcomb_list, tmpoffcomb_list, tmpflatin, l_statextn, first_val
    char orig_first_val, next_statvalue, sciext, outgaimchk_list, tmpfile
    char l_redirect, l_obs_allowed, l_obj_allowed, new_metaconf, l_arraysdb
    char l_non_lcdb, tmpfix, proc_state, l_orig_reject, l_orig_combine
    char tmpbpm, bpmgemcomb, area_bpm, next_stddev, l_use_offs

    int comma_pos, lcount, num_created, open_location, comb_num, location_type
    int nmatch, i, j, k, l_test, l_test2, last_char, last2_char, nsci
    int l_minflat, l_nlow, l_nhigh, l_nkeep, nconf, nstat, total_num_created

    real l_maxtime, norm_factor, on_compare, off_compare, time_diff, t_min
    real statval, relint_val, norm_upper, norm_lower, stddev, first_num
    real l_maskvalue, l_lthreshold, l_hthreshold, l_lsigma, l_hsigma
    real l_ron, l_gain, l_sigscale, l_pclip, l_grow, next_num, norm_err

    bool l_fl_igexp, dooffs, processing, skipped, l_fl_dqprop, l_fl_scale
    bool l_fl_mask, skip_conf, finished, off_flats, list_created, l_ignore_nlc
    bool l_fl_vardq, l_mclip, l_verbose, l_fl_trim, l_fl_nlc, l_fl_sat
    bool debug

    # Used in the code to test for flat types
    char types[4]="GCAL","DOME","TWLT","OBJ"
    int num_types=3

    # Used as the default gemcombine settings - indicies are as for types
    # There are default nhigh/low settings dependant on number of files in the
    # l_combine=default loop too. - MS
    char def_combine[3]="average","median","median"
    char def_reject[3]="avsigclip","minmax","minmax"

    # Uncomment the user parameter if scaling is wished to be allowed and then
    # comment this out
    bool fl_scale = no

    int num_off_types = 2     # Number of flat types which require subtracting ON - OFF
    int num_filt_types = 15   # Number of filters which require subtracting ON - OFF
    char off_types[2]         # Types of flat field which require subtracting ON - OFF
    off_types[1] = "GCAL"
    off_types[2] = "DOME"
    char filt_types[2,15]     # [OFF_TYPE, FILTER]
    filt_types[1, 1] = "ALL"  # All GCAL flats must subtract ON - OFF
    filt_types[2, 1] = "H_G1103"
    filt_types[2, 2] = "Kprime_G1104"
    filt_types[2, 3] = "Kshort_G1105"
    filt_types[2, 4] = "K_G1106"
    filt_types[2, 5] = "CH4short_G1109"
    filt_types[2, 6] = "CH4long_G1110"
    filt_types[2, 7] = "Kcntshrt_G1111"
    filt_types[2, 8] = "Kcntlong_G1112"
    filt_types[2, 9] = "FeII_G1118"
    filt_types[2,10] = "H2O_G1119"
    filt_types[2,11] = "HeI-2p2s_G1120"
    filt_types[2,12] = "H2(1-0)_G1121"
    filt_types[2,13] = "BrG_G1122"
    filt_types[2,14] = "H2(2-1)_G1123"
    filt_types[2,15] = "CO2360_G1124"

    ####
    ##M Make sure these are in the same order as the parameters
    # Set local variables
    l_inimages = inimages
    l_rawpath = rawpath
    l_outsufx = outsufx
    l_logfile = logfile
    l_rootname = rootname
    l_badpix = badpix
    l_combine = combine
    l_reject = reject
    l_masktype = masktype
    l_zero = zero
    l_weight = weight
    l_statsec = statsec
    l_expname = expname
    l_key_ron = key_ron
    l_key_gain = key_gain
    l_snoise = snoise
    l_fl_vardq = fl_vardq
    l_mclip = mclip
    l_verbose = verbose
    l_minflat = minflat
    l_nlow = nlow
    l_nhigh = nhigh
    l_nkeep = nkeep
    l_maskvalue = maskvalue
    l_lthreshold = lthreshold
    l_hthreshold = hthreshold
    l_lsigma = lsigma
    l_hsigma = hsigma
    l_ron = ron
    l_gain = gain
    l_sigscale = sigscale
    l_pclip = pclip
    l_grow = grow
    l_fl_trim = fl_trim
    l_fl_sat = fl_sat
    l_fl_nlc = fl_nlc
    l_maxtime = maxtime
    l_datename = datename
    l_timename = timename
    l_ignore_exp = ignore_exp
    l_gaprep_pref = gaprep_pref
    l_fl_dqprop = fl_dqprop
    l_stattype = stattype
    l_fl_scale = fl_scale
    l_fl_mask = fl_mask
    l_statextn = statextn
    l_arraysdb = arraysdb
    l_non_lcdb = non_lcdb
    l_ignore_nlc = ignore_nlc
    l_use_offs = use_offs
    l_sci_ext = sci_ext
    l_var_ext = var_ext
    l_dq_ext = dq_ext

    ####
    # Set temporary files
    filelist = mktemp("tmpfile")
    tmpgaprep = mktemp ("tmpgaprep")
    nextlist = mktemp ("tmpnextlist")
    tmplist = mktemp ("tmplist")
    todelete = mktemp ("tmptodelete")
    imtmptodel = mktemp ("tmpimtodel")
    tmp_skipped = mktemp ("tmp_skipped")
    tmpnotproc = mktemp ("tmpnotproc")
    outgaimchk_list = mktemp ("tmpgaimchk")
    tmpbpm = mktemp ("tmpbpm")//".fits"

    ####
    # Set default values
    debug = no
    status = 0
    nconf = 1
    l_orig_combine = l_combine
    l_orig_reject = l_reject

    if (l_verbose) {
        l_redirect = "STDOUT"
    } else {
        l_redirect = "dev$null"
    }

    ########
    # Here is where the actual work starts

    # Test the log file
    if ((l_logfile == "") || (stridx(" ", l_logfile) > 0)) {
       l_logfile = gsaoi.logfile
       if ((l_logfile == "") || (stridx (" ", l_logfile) > 0)) {
          l_logfile = "gsaoi.log"
          printlog ("WARNING - GAFLAT: Both gaflat.logfile and \
              gsaoi.logfile fields are empty", l_logfile, l_verbose)
          printlog ("                     Using default file gsaoi.log", \
              l_logfile, l_verbose)
       }
    }

    printlog ("", l_logfile, l_verbose)
    printlog ("------------------------------------------------------------", \
        l_logfile, l_verbose)
    gemdate (zone="local")
    printlog ("GAFLAT -- Started: "//gemdate.outdate//"\n", \
        l_logfile, l_verbose)

    # Record all input parameters relavent to this task only - other tasks will
    # print their inputs to log
    printlog ("GAFLAT: Input parameters...", l_logfile, l_verbose)
    printlog ("    inimages    = "//l_inimages, l_logfile, l_verbose)
    printlog ("    rawpath     = "//l_rawpath, l_logfile, l_verbose)
    printlog ("    outsufx     = "//l_outsufx, l_logfile, l_verbose)
    printlog ("    rootname    = "//l_rootname, l_logfile, l_verbose)
    printlog ("    minflat     = "//l_minflat, l_logfile, l_verbose)
    printlog ("    datename    = "//l_datename, l_logfile, l_verbose)
    if (isindef(l_maxtime)) {
        printlog ("    maxtime     = INDEF", l_logfile, l_verbose)
    } else {
        printlog ("    maxtime     = "//l_maxtime, l_logfile, l_verbose)
    }
    printlog ("    timename    = "//l_timename, l_logfile, l_verbose)
    printlog ("    ignore_exp  = "//l_ignore_exp, l_logfile, l_verbose)
    printlog ("    stattype    = "//l_stattype, l_logfile, l_verbose)
    printlog ("    statextn    = "//l_statextn, l_logfile, l_verbose)
    printlog ("    statsec     = "//l_statsec, l_logfile, l_verbose)
    printlog ("    fl_mask     = "//l_fl_mask, l_logfile, l_verbose)
    printlog ("    badpix      = "//l_badpix, l_logfile, l_verbose)
#    printlog ("    fl_scale   = "//l_fl_scale, l_logfile, l_verbose)
    printlog ("    fl_vardq    = "//l_fl_vardq, l_logfile, l_verbose)
    printlog ("    ignore_nlc  = "//l_ignore_nlc, l_logfile, l_verbose)
    printlog ("    gaprep_pref = "//l_gaprep_pref, l_logfile, verbose-)
    printlog ("    fl_trim     = "//l_fl_trim, l_logfile, verbose-)
    printlog ("    fl_nlc      = "//l_fl_nlc, l_logfile, verbose-)
    printlog ("    fl_sat      = "//l_fl_sat, l_logfile, verbose-)
    printlog ("    sci_ext     = "//l_sci_ext, l_logfile, verbose-)
    printlog ("    var_ext     = "//l_var_ext, l_logfile, verbose-)
    printlog ("    dq_ext      = "//l_dq_ext, l_logfile, verbose-)
    printlog ("    arraysdb    = "//l_arraysdb, l_logfile, verbose-)
    printlog ("    non_lcdb    = "//l_non_lcdb, l_logfile, verbose-)
    printlog ("    combine     = "//l_combine, l_logfile, l_verbose)
    printlog ("    reject      = "//l_reject, l_logfile, l_verbose)
    printlog ("    fl_dqprop   = "//l_fl_dqprop, l_logfile, l_verbose)
    printlog ("    logfile     = "//l_logfile, l_logfile, l_verbose)
    printlog ("", l_logfile, l_verbose)

    ####
    # Start user input checks and data validation

    # Check the l_badpix value
    if ((l_badpix == "") || (stridx(" ", l_badpix) > 0)) {
        # This assumes that the bpm covers all of the areas readout in each
        # detector in each image
        l_badpix = "none"
    } else if (!imaccess(l_badpix)) {
        printlog ("ERROR - GAFLAT: cannot access bpm "//l_badpix//\
            ". Exiting", l_logfile, verbose+)
    } else {
        # Check it hasn't been mosaiced
        keypar (l_badpix//"[0]", "GAMOSAIC", silent+)
        if (keypar.found) {
            printlog ("ERROR - GAFLAT: "//l_badpix//" has been mosaiced."//\
                " Please supply a static BPM that has not been mosaiced.", \
                l_logfile, verbose)
            goto crash
        }
    }
    printlog ("GAFLAT: badpix value after checks is \""//l_badpix//"\"", \
        l_logfile, l_verbose)

    # The checks in gaimchk convert these strings to lower case!
    l_obs_allowed = "OBJECT,FLAT"
    l_obj_allowed = "TWILIGHT,DOMEFLAT,DOMEFLAT OFF,GCALFLAT"

    # Call gaimchk to perform inout checks
    gaimchk (inimages=l_inimages, rawpath=l_rawpath, rootname=l_rootname, \
        obstype_allowed=l_obs_allowed, object_allowed=l_obj_allowed, \
        key_allowed="GAPREPAR", key_forbidden="GAFLAT", key_exists="", \
        fl_prep_check=yes, gaprep_pref=l_gaprep_pref, fl_redux_check=no, \
        garedux_pref="r", fl_fail=no, fl_out_check=no, \
        fl_vardq_check=l_fl_vardq, \
        sci_ext=l_sci_ext, var_ext=l_var_ext, dq_ext=l_dq_ext, \
        outlist=outgaimchk_list, \
        logfile=l_logfile, verbose=no)

    if (gaimchk.status != 0) {
        printlog ("ERROR - GAFLAT: GAIMCHK returned a non-zero status. \
            Exiting.", l_logfile, verbose+)
        goto crash
    } else if (!access(outgaimchk_list)) {
        printlog ("ERROR - GAFLAT: Cannot access output list from GAIMCHK", \
            l_logfile, verbose+)
        goto crash
    }

    # Files that have been processed by GALFAT
    tmpfile = ""
    match ("tmpGAFLAT", outgaimchk_list, stop=no, print_file_name=no, \
        metacharacter=no) | scan (tmpfile)
    if (tmpfile != "") {
        printlog ("WARNING - GAFLAT: The following files have already been \
            processed by GAFLAT. Ignoring these files.", l_logfile, verbose+)
        fields (tmpfile, fields=1, lines="1-", quit_if_miss=no, \
            print_file_names=no) | \
            tee (l_logfile, out_type="text", append+)
    }

    # Prepared files
    tmpfile = ""
    match ("tmpGAPREPAR", outgaimchk_list, stop=no, print_file_name=no, \
        metacharacter=no) | scan (tmpfile)
    if (tmpfile != "") {
        type (tmpfile, >> nextlist)
        print (nextlist, >> todelete)
    }

    # Files that need to be prepared
    tmpfile = ""
    match ("tmpNotFound", outgaimchk_list, stop=no, print_file_name=no, \
        metacharacter=no) | scan (tmpfile)
    if (tmpfile != "") {

        fields (tmpfile, fields=1, lines="1-", quit_if_miss=no, \
            print_file_names=no, >> tmpgaprep)

        print (tmpgaprep, >> todelete)

        printlog ("\nGAFLAT: Calling gaprepare to process unprepared input \
            files", l_logfile, l_verbose)

        gaprepare ("@"//tmpgaprep, rawpath="",outpref=l_gaprep_pref, \
            rootname="", fl_trim=l_fl_trim, fl_nlc=l_fl_nlc, \
            fl_vardq=l_fl_vardq, fl_sat=l_fl_sat, badpix=l_badpix, \
            sci_ext=l_sci_ext, var_ext=l_var_ext, dq_ext=l_dq_ext, \
            arraysdb=l_arraysdb, non_lcdb=l_non_lcdb, \
            logfile = l_logfile, verbose=l_verbose)

        if (gaprepare.status != 0) {
            printlog ("ERROR - GAFLAT: GAPREPARE returned an non-zero status. \
                Exiting.", l_logfile, verbose+)
            goto crash
        }

        fields (tmpfile, fields=2, lines="1-", quit_if_miss=no, \
            print_file_names=no, >> tmplist)
        sections (l_gaprep_pref//"//@"//tmplist, \
            option="fullname", >> nextlist)

        delete (tmplist, verify-, >& "dev$null")
        delete (tmpgaprep, verify=no, >& "dev$null")

    } # End of prepare loop

    # Check if there are any images in the next list list
    if (!access(nextlist)) {
        printlog ("ERROR - GAFLAT: No input images can be used. \
            Please try again.", l_logfile, verbose+)
        goto crash
    } else {
        print (nextlist, >> todelete)
    }

    # Loop over nextlist to obtain the METACONFIG keyword. This contains
    # information that will group files such that a) they are in the same
    # configuration and b) that they have been reduced in the same way.
    printlog ("\nGAFLAT: Accessing METACONF keyword values", \
        l_logfile, l_verbose)

    scanfile = nextlist
    while (fscan(scanfile, inimg) != EOF) {

        inphu = inimg//"[0]"

        # Double check that you can access the file...
        if (!imaccess(inimg)) {
            printlog ("WARNING - GAFLAT: Image "//inimg//\
                " does not exist. Skipping this image.", l_logfile, verbose)
            print (inimg, >> tmpnotproc)

        } else {
            # Read the METACONFIG keyword
            keypar (inphu, "METACONF", silent+)
            if (!keypar.found) {
                printlog ("WARNING - GAFLAT: METACONF keyword not found in "//\
                    inimg//" Skipping this image.", l_logfile, verbose+)
                print (inimg, >> tmpnotproc)

            } else {
                metaconf = keypar.value

                # Check it's valid
                if (metaconf == "UNKOWN" || metaconf == "" || \
                    (stridx(" ",metaconf) > 0)) {
                    printlog ("WARNING - GAFLAT: METACONF keyword is "//\
                        "\"UNKNOWN\" or bad in"//inimg//".\n"//\
                        "                  Skipping this image.", \
                        l_logfile, verbose+)
                    print (inimg, >> tmpnotproc)
                } else {


                    # Decide on whether the exposure time can be ignored
                    # Default is to only include exposure time if gcal
                    l_fl_igexp = no
                    if (l_ignore_exp == "yes") {
                        l_fl_igexp = yes

                    } else if (l_ignore_exp == "default") {

                        if (strstr("TWLT",metaconf) > 0) {
                            l_fl_igexp = yes
                        }
                    }

                    # Remove the exposure time if requested to ignore it
                    if (l_fl_igexp) {
                        metaconf = substr(metaconf,stridx("+",metaconf)+1,\
                            strlen(metaconf))
                    }

                    if (l_ignore_nlc) {
                        l_test = strstr("NLC",metaconf)
                        if (l_test > 0) {
                            # NLC is always at the end
                            metaconf = substr(metaconf,1,l_test-2)
                        }
                    }

                    print (inimg//" "//metaconf, >> filelist)
                }
            }
        }

    } # End of loop to read METACONF keyword

    # Check if there are any images in the newly created filelist
    if (!access(filelist)) {
        printlog ("ERROR - GAFLAT: No images with a valid METACONF keyword \
            value. Please try again.", l_logfile, verbose+)
        goto crash
    }

    # Determine the unique METACONF values for the file that made it
    printlog ("GAFLAT: Determining unique METACONFIG values", \
        l_logfile, l_verbose)

    tmpconflist = mktemp ("tmpconflist")

    # Make sure the reverse_sort is true to have closed GCAL flats after open
    # ones - MS
    fields (filelist, fields=2, lines="1-", quit_if_miss=no, \
        print_file_names=no) | sort ("STDIN", column=0, ignore_white=no, \
        numeric_sort=no, reverse_sort=yes) | unique ("STDIN", > tmpconflist)

    count (tmpconflist) | scan (lcount)

    # Check the list was created
    if (access(tmpconflist) && lcount > 0) {
        print (tmpconflist, >> todelete)
    } else {
        printlog ("ERROR - GAFLAT: No unique metaconfigurations values "//\
            "found.\n                Exiting", l_logfile, verbose+)
        goto crash
    }

    ####
    # Create the FLATs...

    total_num_created = 0

    # This becomes a maze of lists! Reader beware.

    # Create separate lists for each unique METACONF keyword
    printlog ("GAFLAT: Creating METACONFIG lists and processing images...", \
        l_logfile, l_verbose)

    tmplist = mktemp("tmplist")

    scanfile = tmpconflist
    while (fscan(scanfile, metaconf) != EOF) {

        printlog ("\n----\nGAFLAT: Processing configuration: "//\
            metaconf//"\n", l_logfile, l_verbose)
        # Check for types that require off flats to be subtracted
        ##M Do we fail if no off lamps? Different for GCALs / DOMEs? (If we
        ##M does this with DOMEs?
        ##M Put check in for CLOS but no open!!

        dooffs = no
        if (l_use_offs == "yes") {
            dooffs = yes
            printlog ("GAFLAT: Will look for OFF flats too...", \
                l_logfile, l_verbose)

        } else if (l_use_offs == "default") {

            for (i = 1; i <= num_off_types; i += 1) {
                if (strstr(off_types[i],metaconf) > 1) {

                    # Dome flats are separated out by FILTER for on off flats
                    for (j = 1; j <= num_filt_types; j += 1) {
                        if ((strstr(filt_types[i,j],metaconf) > 1) || \
                            (filt_types[i,j] == "ALL")) {

                            printlog ("GAFLAT: FLAT type of "//off_types[i]//\
                                " will look for OFF flats too...", \
                                l_logfile, l_verbose)
                            dooffs = yes
                            break
                        }
                    }
                }
            }
        } # Finished setting dooffs flag

        # Determine type for title later
        for (i = 1; i <= num_types; i += 1) {
            if (strstr(types[i],metaconf) > 1) {
                current_type = types[i]
                location_type = i
                break
            }
        }

        # Set up temporary configuration list
        conflist = tmplist//"_on_conf"//nconf//".lis"

        # Extract the files that have the same METACONF keyword value
        match (metaconf, filelist, stop=no, print_file_name=no, \
            metacharacter=no) | fields ("STDIN", fields=1, lines="1-", \
            quit_if_miss=no, print_file_names=no, > conflist)

        # Keep track of files to delete at the end
        print (conflist, >> todelete)

        bpmgemcomb = ""
        # Check the BPM dimensions against the first image in the list
        if (l_badpix != "none" && l_masktype != "none") {

            # Use the first file in the list - the rest should all be the same
            # dimensions - MS
            refname = ""
            head (conflist, nlines=1) | scan (refname)

            keypar (refname//"[0]", "NSCIEXT", silent+)
            nsci = int(keypar.value)
            for (i = 1; i <= nsci; i += 1) {
                gadimschk (refname//"["//l_sci_ext//","//i//"]", \
                    section="", \
                    chkimage=l_badpix//"["//l_dq_ext//","//i//"]", \
                    key_check="CCDSEC", logfile=l_logfile, verbose-)

                if (gadimschk.status != 0) {
                    printlog ("ERROR - GAFLAT: GADIMSCHK returned a non-zero \
                        status. Exiting.", l_logfile, verbose-)
                    goto crash
                } else {
                    # Output from gadimschk has always a sections appended
                    area_bpm = substr(gadimschk.out_chkimage,\
                        strldx("[",gadimschk.out_chkimage),\
                        strlen(gadimschk.out_chkimage)) == "[*,*]"

                    if (area_bpm == "[*,*]") {
                        bpmgemcomb = l_badpix
                    } else if (i == 1) {
                        imcopy (l_badpix//"[0]", tmpbpm, verbose-)
                        imcopy (gadimschk.out_chkimage, tmpbpm//"[append]", \
                            verbose-)
                        bpmgemcomb = tmpbpm
                    } else if (!imaccess(tmpbpm) && i > 1) {
                        printlog ("ERROR - GAFLAT: The array dimensions are \
                            different sizes for input image. Exiting.", \
                            l_logfile, verbose+)
                    } else {
                        imcopy (gadimschk.out_chkimage, tmpbpm//"[append]", \
                            verbose-)
                    }
                }
            } # End of loop over sci extensions
            refname = ""
        } # End of bpm dimensions check

        # Reset the num_created counter
        num_created = 0

        # Check for OFF flats for this configuration
        if (dooffs) {

            # This will not work for twilights - MS
            open_location = strstr("OPEN",metaconf)
            if (open_location > 0) {
                off_metaconf = substr(metaconf,1,strstr("OPEN",metaconf)-1)//\
                    "CLOS"//substr(metaconf,strstr("OPEN",metaconf)+4,\
                    strlen(metaconf))

                # Set up temporary configuration list
                off_conflist = tmplist//"_off_conf"//nconf//".lis"

                # Extract the files that have the same METACONF keyword value
                match (off_metaconf, filelist, stop=no, print_file_name=no, \
                    metacharacter=no) | fields ("STDIN", fields=1, lines="1-",\
                    quit_if_miss=no, print_file_names=no, > off_conflist)

                # Keep track of files to delete at the end
                print (off_conflist, >> todelete)

                # Check there are OFFs
                count (off_conflist) | scan (lcount)
                if (!access(off_conflist) || lcount <= 0) {
                    printlog ("WARNING - GAFLAT: No OFF flats found for \
                        configuration: "//\
                        "\n                  "//metaconf, \
                        l_logfile, verbose+)
                    dooffs = no
                    # Print additional warning if l_use_offs == yes
                    if (l_use_offs == "yes") {
                        printlog ("                  use_flats is yes. \
                            Confirm that the METACONF \n"//\
                            "                  values of the OFF flats \
                            contain \"+CLOS\" immediately after \
                            the flat type", \
                            l_logfile, verbose+)
                    }

                    printlog ("                  Skipping configuration\n", \
                        l_logfile, verbose+)

                    goto NEXT_CONFIGURATION
                } else {
                    printlog ("        Found associated OFF flats:"//\
                        "\n            "//off_metaconf, \
                        l_logfile, verbose+)
                }
            } else {
                dooffs = no

                if (l_use_offs == "yes" || l_use_offs == "default") {
                    printlog ("WARNING - GAFLAT: METACONF doesn't contain \
                       the \"OPEN\" string.\n"//\
                       "                  Will not subtract OFF flats.", \
                       l_logfile, verbose+)
                }
            }
        }

        oncomb_list = mktemp("tmpONCOMB")
        offcomb_list = mktemp("tmpOFFCOMB")

        # Do an initial check of the number of input file in the two lists
        count (conflist) | scan (lcount)
        if (debug) {
            printlog ("___conflist lcount = "//lcount, l_logfile, verbose+)
        }

        if (lcount < l_minflat) {
            printlog ("WARNING - GAFLAT: Number of input files is \
                less than "//l_minflat//" (minflat) for \
                configuration: "//\
                "\n                      "//metaconf//\
                "\n                  Skipping configuration.", \
                l_logfile, verbose+)

            goto NEXT_CONFIGURATION

        } else if (dooffs) {
            count (off_conflist) | scan (lcount)
            if (debug) {
                printlog ("___off_conflist lcount = "//lcount, \
                    l_logfile, verbose+)
            }

            if (lcount < l_minflat) {
                printlog ("WARNING - GAFLAT: Number of input files is \
                    less than "//l_minflat//" (minflat) for \
                    configuration: "//\
                    "\n                      "//off_metaconf//\
                    "\n                  Skipping configuration:"//\
                    "\n                      "//metaconf,
                    l_logfile, verbose+)

                goto NEXT_CONFIGURATION
            }
        }

        if (!isindef(l_maxtime)) {
            # Set flags, temporary files to be used inside the next loop
            processing = yes
            off_flats = no
            finished = no
            list_created = no
            onlist = conflist
            inlist = mktemp("tmpinlist")
            offlist = mktemp("tmpofflist")
            tmponcomb_list = mktemp("tmponcomblist")
            tmpoffcomb_list = mktemp("tmpoffcomblist")

            skipped = no # This says whether the current configuation list has
                         # been processed more than once

            # Keep track of files to delete at the end
            print (inlist, >> todelete)
            print (offlist, >> todelete)

            printlog ("GAFLAT: Sorting OFF flats by time", \
                l_logfile, l_verbose)

            skip_conf = no

            # Create the flats
            while (processing) {

                # Determine the number of files in onlist
                count (onlist) | scan (lcount)

                if (debug) {
                    printlog ("HERE 3 - lcount = "//lcount, \
                        l_logfile, verbose+)
                }

                # Only continue if the number of files in onlist are greater
                # than the minimum requested number of flats required to create
                # a master flat number of flats to combine
                if (lcount < l_minflat) {

                    if (skipped) {
                        tmpinsert = "now "
                    } else {
                        tmpinsert = ""
                    }

                    printlog ("\nWARNING - GAFLAT: Number of input files for \
                        configuration "//\
                        "\n                  "//metaconf//\
                        "\n                  is "//tmpinsert//\
                        "less than the minimum "//\
                        "requested "//minflat//" files.", \
                        l_logfile, verbose+)

                    if (!skipped) {
                        printlog ("                  Skipping this \
                            configuration", l_logfile, verbose+)
                    } else {
                        printlog ("                  Moving on from this \
                            configuration", l_logfile, verbose+)
                    }

                    processing = no
                    finished = yes
                    skip_conf = yes
                    goto SKIP_PROCESSING
                }

                # Read the reference file
                head (onlist, nlines=1) | scan (refname)

                if (debug) {
                    printlog ("Refname -- "//refname, l_logfile, verbose+)
                }

                # Find any matching files by time
                # l_distance is INDEF so only seperating by time
                gemoffsetlist ("@"//onlist, refname, INDEF, l_maxtime, \
                    inlist, offlist, fl_nearer=yes, direction=3, \
                    fl_younger=yes, fl_noref=no, wcs_source="none", \
                    key_xoff="XOFFSET",  key_yoff="YOFFSET", \
                    key_date=l_datename, key_time=l_timename, \
                    logfile=l_logfile, verbose=no, force=no, >& "dev$null")

                # Check the output status
                if (gemoffsetlist.status != 0) {
                    printlog ("ERROR - GAFLAT: GEMOFFSETLIST returned a \
                        non-zero status. Exiting.", l_logfile, verbose+)
                    goto crash
                }

                # Number of matches found -- will always contain the reference
                # file
                nmatch = gemoffsetlist.count

                if (debug) {
                    printlog ("NMATCH -- "//nmatch, l_logfile, verbose+)

                    printlog ("HERE 3a matched list...", l_logfile, verbose+)

                    type (inlist) | tee (l_logfile, out_type="text", append+)

                    if (access(offlist)) {
                        printlog ("HERE 3b offlist...", l_logfile, verbose+)
                        type (offlist) | \
                            tee (l_logfile, out_type="text", append+)
                    }
                }

                # Check the number of matches against the requested minimum
                # number of flats to combine
                if (nmatch < l_minflat) {

                    # If not enough, remove the first and use the next as
                    # reference
                    ##M Edit this warning!
                    printlog ("WARNING - GAFLAT: Number of images in"//\
                        " list within "//l_maxtime//" seconds of", \
                        l_logfile, l_verbose)
                    printlog ("                  reference "//refname//\
                        " less than requested minimum.",l_logfile, l_verbose)

                    # Remove the temporary files
                    delete (offlist//","//inlist, verify-, >& "dev$null")

                    # Remove the first file (reffile) from onlist and reset
                    # onlist
                    tail (onlist, nl=-1, >> offlist)

                    # Rename offlist to be the new onlist
                    delete (onlist, verify-, >& "dev$null")
                    rename (offlist, onlist, field="all")

                    if (debug) {
                        printlog ("HERE 3c new onlist...", l_logfile, verbose+)
                        type (onlist) | \
                            tee (l_logfile, out_type="text", append+)
                    }

                    # Set the skipped flag
                    skipped = yes

                } else {
                    # Enough were selected to combine [there may be files that
                    # were not selected too]

                    if (debug) {
                        printlog ("HERE 4", l_logfile, verbose+)
                    }

                    # Print list to a tmporary file to be printed into the
                    # oncomb_list later
                    if (!off_flats) {
                        # Print list to a tmporary file to be printed into the
                        # oncomb_list later
                        print (inlist, >> tmponcomb_list)
                        print (tmponcomb_list, >> todelete)
                    } else {
                        # Ditto for the off flats
                        print (inlist, >> tmpoffcomb_list)
                        print (tmpoffcomb_list, >> todelete)
                    }

                    list_created = yes
                    inlist = mktemp("tmpinlist")

                    print (inlist, >> todelete)

                    # If there's a list of files that didn't get selected
                    # rename it and set it back to be the next onlist list
                    if (access(offlist)) {
                        if (debug) {
                            printlog ("HERE 6", l_logfile, verbose+)
                        }

                        delete (onlist, verify-, >& "dev$null")
                        rename (offlist, onlist, field="all")
                        delete (inlist, verify-, >& "dev$null")
                        skipped = yes

                    } else {
                        processing = no
                        finished = yes
                    }
                }

SKIP_PROCESSING:

                if (debug) {
                    printlog ("____finished: "//finished//", "//\
                        "skip_conf: "//skip_conf//", "//\
                        "dooffs: "//dooffs//", "//\
                        "\n____off_flats: "//off_flats//", "//\
                        "list_created: "//list_created, \
                        l_logfile, verbose+)
                }

                # Figure out what is to be done
                if (finished && skip_conf) {
                    if (!dooffs && list_created && !off_flats) {
                        if (debug) printlog ("____IN A", l_logfile, verbose+)

                        processing = no
                        type (tmponcomb_list, >> oncomb_list)
                    } else if (dooffs && list_created && !off_flats) {
                        if (debug) printlog ("____IN B", l_logfile, verbose+)

                        finished = no
                        off_flats = yes
                        processing = yes
                        onlist = off_conflist
                        skip_conf = no
                        printlog ("GAFLAT: Sorting OFF flats by time", \
                            l_logfile, l_verbose)
                    } else if (dooffs && list_created && off_flats) {
                        if (debug) printlog ("____IN C", l_logfile, verbose+)

                        processing = no
                        type (tmponcomb_list, >> oncomb_list)
                        type (tmpoffcomb_list, >> offcomb_list)
                    } else if (!list_created && !off_flats) {
                        if (debug) printlog ("____IN D", l_logfile, verbose+)

                        goto NEXT_CONFIGURATION
                    } else if (!list_created && dooffs && off_flats) {
                        if (debug) printlog ("____IN E", l_logfile, verbose+)

                        goto NEXT_CONFIGURATION
                    }
                } else if (dooffs && finished && !off_flats) {
                    if (debug) printlog ("____IN F", l_logfile, verbose+)

                    finished = no
                    off_flats = yes
                    processing = yes
                    onlist = off_conflist
                    skip_conf = no
                    printlog ("GAFLAT: Sorting OFF flats by time", \
                        l_logfile, l_verbose)
                } else if (finished && !dooffs) {
                    if (debug) printlog ("____IN G", l_logfile, verbose+)

                    type (tmponcomb_list, >> oncomb_list)
                } else if (finished && off_flats) {
                    if (debug) printlog ("____IN H", l_logfile, verbose+)

                    type (tmponcomb_list, >> oncomb_list)
                    type (tmpoffcomb_list, >> offcomb_list)
                }

            } # End of while processing loop


            print (tmponcomb_list, >> todelete)
            print (tmpoffcomb_list, >> todelete)

        } else {

            # No time offsets so just dump the ecurrent configuation list inot
            # the oncomb_list to combine
            print (conflist, >> oncomb_list)
            print (oncomb_list, >> todelete)

            if (dooffs) {
               # Ditto for the off flats
               print (off_conflist, >> offcomb_list)
               print (offcomb_list, >> todelete)
            }
        }

        if (debug) {
            printlog ("____oncomblist...", l_logfile, verbose+)
            type (oncomb_list) | tee (l_logfile, out_type="text", append+)
        }

        print (oncomb_list, >> todelete)

        if (debug) printlog ("____dooffs = "//dooffs, l_logfile, verbose+)

        comb_num = 2
        if (!access(offcomb_list) && dooffs) {
            printlog ("ERROR - GAFLAT: Cannot find list of OFF flats to \
                combine.", l_logfile, verbose+)
            goto crash
        } else if (!dooffs) {
            comb_num = 1
        } else {
            print (offcomb_list, >> todelete)
        }

        if (debug) {
            if (comb_num == 2) {
                printlog ("____offcomblist...", l_logfile, verbose+)
                type (offcomb_list) | tee (l_logfile, out_type="text", append+)
            }
        }

        out_oncomb_list = mktemp("tmp_out_oncomblist")
        print (out_oncomb_list, >> todelete)

        out_offcomb_list = mktemp("tmp_out_offcomblist")
        print (out_oncomb_list, >> todelete)

        nstat = 1
        for (i = 1; i <= comb_num; i += 1) {

            if (i == 2) {
                scanfile2 = offcomb_list

                printlog ("GAFLAT: Working on OFF flats...", \
                    l_logfile, l_verbose)
            } else {
                scanfile2 = oncomb_list
            }

            # Determine the closest one in time for OFF flats?

            scale_list = mktemp ("tmpscale_list")
            while (fscan(scanfile2, comblist) != EOF) {

                # Determine the output name
                head (comblist) | scan (tmpflatin)
                fparse (tmpflatin, verbose-)

                if (i == 1) {
                    # Check that the output doesn't already exist - Can do this
                    # here because output name is based on input name always
                    out_flat = fparse.root//"_"//l_outsufx//".fits"

                    if (debug) printlog ("___TEST output name is "//out_flat, \
                        l_logfile, verbose+)

                    if (imaccess(out_flat)) {
                        printlog ("ERROR - GAFLAT: "//out_flat//\
                            " already exists.", l_logfile, verbose+)
                        goto crash
                    }

                    # Set up tmp names for (ON) combined flats
                    tmpflat = fparse.root//"_GAtmpon"//mktemp("tmp")//".fits"

                } else {
                    # Set up tmp names for OFF combined flats
                    tmpflat = fparse.root//"_GAtmpoff"//mktemp("tmp")//".fits"
                }

                count (comblist) | scan (lcount)

                # Scaling is off and not open to the user - MS
                # Need to scale some how if required
                # Only do this is scaling using the whole detector
                # Need somechecks for stats sec etc.
                if (l_fl_scale && lcount > 1) {

                    printlog ("GAFLAT: Determining scaling factors...", \
                        l_logfile, l_verbose)

                    k = 0
                    scanfile3 = comblist
                    while (fscan(scanfile3, img) != EOF) {
                        k += 1

                        # Call gastat
                        gastat (img, stattype=l_stattype, statsec=l_statsec, \
                            statextn=l_statextn, fl_mask=l_fl_mask, \
                            badpix=l_badpix, calc_stddev=no, lower=INDEF, \
                            upper=INDEF, nclip=0, lsigma=INDEF, usigma=INDEF, \
                            sci_ext=l_sci_ext, dq_ext=l_dq_ext, \
                            logfile=l_logfile, verbose-)

                        if (gastat.status != 0) {
                            printlog ("ERROR - GAFLAT: GASTAT returned a \
                                non-zero status. Exiting.", \
                                l_logfile, verbose+)
                            goto crash
                        }

                        next_statvalue = gastat.outstat

                        if (debug) printlog ("____k: "//k//\
                            " next_statvalue: "//\
                            next_statvalue, l_logfile, verbose+)

                        # Now caclculate relative intensity
                        if (k == 1) {
                            orig_first_val = next_statvalue
                        } else {
                            first_val = orig_first_val

                            l_test = 1
                            l_test2 = 1
                            j = 0

                            while (l_test > 0 && l_test2 > 0) {
                                j += 1

                                l_test = stridx (",",first_val)
                                l_test2 = stridx (",",next_statvalue)

                                if (l_test == 0) {
                                    last_char = strlen(first_val)
                                } else {
                                    last_char = l_test - 1
                                }

                                if (l_test2 == 0) {
                                    last2_char = strlen(next_statvalue)
                                } else {
                                    last2_char = l_test - 1
                                }

                                # Parse the current value
                                first_num = real(substr(first_val,1,last_char))
                                next_num = real(substr(next_statvalue,1,\
                                    last2_char))

                                relint_val = first_num / next_num

                                if (k == 2) {
                                    print ("1.0", >> scale_list//"_"//j)
                                    print (scale_list//"_"//j, >> todelete)
                                }

                                if (debug) {

                                    printlog ("____first_num: "//\
                                        first_num//" next_num: "//\
                                        next_num//" relint_val: "//\
                                        relint_val, l_logfile, verbose+)
                                }

                                print (relint_val, >> \
                                    scale_list//"_"//j)

                                # Rest the output statistics for next value
                                first_val = substr(first_val,l_test+1,\
                                    strlen(first_val))

                                next_statvalue = substr(next_statvalue,\
                                    l_test+1, strlen(next_statvalue))
                            }

                        }
                    } # End of loop over input images to be scaled

                    l_scale = scale_list//"_1,"//scale_list//"_2,"//\
                        scale_list//"_3,"//scale_list//"_4"

                    printlog ("        Done.", l_logfile, l_verbose)

                } else {
                    l_scale = ""
                } # End of do scale loop

                if (dooffs && i == 1) {
                    l_title = "Combined "//current_type//" ON FLAT"
                } else if (dooffs && i == 2) {
                    l_title = "Combined "//current_type//" OFF FLAT"
                } else {
                    l_title = "Combined "//current_type//" FLAT"
                }

                print (tmpflat, >> imtmptodel)

                if (lcount > 1) {

                    printlog ("\nGAFLAT: Calling GEMCOMBINE to combine \
                        flats...", l_logfile, l_verbose)

                    # Set the default gemcombine parameters if required
                    # Regardless of flat type this will always be what the
                    # user sets it to for each configuration. - MS
                    if (l_orig_combine == "default") {

                        printlog ("\nGAFLAT: combine=default. Setting default \
                            GEMCOMBINE\n        combine and rejection \
                            parameters.", \
                            l_logfile, l_verbose)

                        # Refine by type - defined by location_type - defaults
                        # set at top of script
                        l_combine = def_combine[location_type]
                        l_reject = def_reject[location_type]

                        l_nhigh = 1
                        l_nlow = 1

                        # This overides the default settings if the number of
                        # inputs is too small
                        if (lcount < 5) {
                            l_combine = "average"
                            l_reject = "minmax"
                            l_nlow -= 1
                            printlog("\nWARNING - GAFLAT: Averaging 4 or \
                                fewer images with no low pixels rejected"//\
                                "\n                  and only 1 high pixel \
                                rejected", l_logfile, verbose+)
                        } else if (lcount >= 8) {
                            l_nhigh += 1
                        }
                    } # End of default section

                    # Inform the user of the combination and rejection methods
                    pr_string = strupr(substr(l_combine,1,1))//\
                        strlwr(substr(l_combine,2,strlen(l_combine)))

                    printlog ("\nGAFLAT: "//pr_string//\
                        " combining "//current_type//" flat frames...", \
                        l_logfile, l_verbose)

                    if (l_reject == "minmax") {
                        pr_string = " with "//l_nlow//" low and "//l_nhigh//\
                        " high values rejected"
                    } else {
                        pr_string = ""
                    }

                    printlog ("        Rejection type is "//l_reject//\
                        pr_string, l_logfile, l_verbose)

                    printlog ("\nGAFLAT: Input files are...", \
                        l_logfile, l_verbose)

                    type (comblist) | tee (l_logfile, out_type="text", \
                        append+, >> l_redirect)

                    printlog ("", l_logfile, l_verbose)

                    # Call gemcombine
                    gemcombine ("@"//comblist, tmpflat, title=l_title, \
                        combine=l_combine, offsets="none", reject=l_reject, \
                        masktype=l_masktype, maskvalue=l_maskvalue, \
                        scale=l_scale, zero=l_zero, weight=l_weight, \
                        statsec=l_statsec, expname=l_expname, \
                        lthreshold=l_lthreshold, hthreshold=l_hthreshold, \
                        nlow=l_nlow, nhigh=l_nhigh, nkeep=l_nkeep, \
                        mclip=l_mclip, lsigma=l_lsigma, hsigma=l_hsigma, \
                        key_ron=l_key_ron, key_gain=l_key_gain, \
                        ron=l_ron, gain=l_gain, snoise=l_snoise, \
                        sigscale=l_sigscale, pclip=l_pclip, grow=l_grow, \
                        bpmfile=bpmgemcomb, nrejfile="", sci_ext=l_sci_ext, \
                        var_ext=l_var_ext, dq_ext=l_dq_ext, \
                        fl_vardq=l_fl_vardq, logfile=l_logfile, \
                        fl_dqprop=l_fl_dqprop, \
                        verbose=l_verbose)

                    if (l_scale != "") {
                        if (nstat == 1) {
                            scale_list = substr(l_scale,2,strlen(l_scale))
                        }
                        delete (l_scale, verify-, >& "dev$null")
                    }

                    type (imtmptodel, >> todelete)

                    # Check the output status of gemcombine
                    if (gemcombine.status != 0) {
                        printlog ("ERROR - GAFLAT: gemcombine returned a \
                            non-zero status. Exiting.", l_logfile, verbose+)
                        goto crash
                    } else {
                        printlog ("GAFLAT: Returned from GEMCOMBINE.\n", \
                            l_logfile, l_verbose)
                    }

                } else {

                    printlog ("WARNING - GAFLAT: Only one file in list to be \
                        combined.\n"//\
                        "                  just copying the image.", \
                        l_logfile, l_verbose)

                    fxcopy (tmpflatin, tmpflat, groups="", new_file=yes, \
                        verbose-)
                }


                if (imaccess(tmpflat)) {
                    if (i == 1) {
                        print (tmpflat, >> out_oncomb_list)
                    } else {
                        print (tmpflat, >> out_offcomb_list)
                    }
                } else {
                    printlog ("ERROR - GAFLAT: Cannot access "//tmpflat//\
                       ". Exiting", l_logfile, verbose+)
                    goto crash
                }
            }

        }

        delete (imtmptodel, verify-, >& "dev$null")

        if (debug) {
            printlog ("____out_oncomb_list:", l_logfile, verbose+)
            type (out_oncomb_list) | tee (l_logfile, out_type="text", append+)
        }

        # On and off flats
        ##M Should this not be done before the combining step for the OFF flats
        if (comb_num == 2) {

            if (debug) {
                printlog ("____out_offcomb_list:", l_logfile, verbose+)
                type (out_offcomb_list) | \
                    tee (l_logfile, out_type="text", append+)
            }

            tmpnorm_list = mktemp("tmpnormlist")

            # Loop over the on combined flats
            scanfile2 = out_oncomb_list
            while (fscan(scanfile2, onflat) != EOF) {

                # The observation time of the flat
                keypar (onflat//"[0]", l_timename, silent+)
                if (!keypar.found) {
                   printlog ("ERROR - GAFLAT: "//l_timename//" keyword not \
                       found in "//onflat//". Exiting", l_logfile, verbose+)
                   goto crash
                }
                on_uttime = keypar.value

                keypar (onflat//"[0]", l_datename, silent+)
                if (!keypar.found) {
                   printlog ("ERROR - GAFLAT: "//l_datename//" keyword not \
                       found in "//onflat//". Exiting", l_logfile, verbose+)
                   goto crash
                }
                on_datevalue = keypar.value

                cnvtsec (on_datevalue, on_uttime) | scan (on_compare)

                # Loop over the OFF combined flats
                j = 0
                scanfile3 = out_offcomb_list
                while (fscan(scanfile3, off_flat) != EOF) {

                    j += 1

                    # The observation time of the flat
                    keypar (off_flat//"[0]", l_timename, silent+)
                    if (!keypar.found) {
                       printlog ("ERROR - GAFLAT: "//l_timename//\
                           " keyword not found in "//off_flat//". Exiting", \
                           l_logfile, verbose+)
                       goto crash
                    }
                    off_uttime = keypar.value

                    keypar (off_flat//"[0]", l_datename, silent+)
                    if (!keypar.found) {
                       printlog ("ERROR - GAFLAT: "//l_datename//\
                           " keyword not found in "//off_flat//". Exiting", \
                           l_logfile, verbose+)
                       goto crash
                    }
                    off_datevalue = keypar.value

                    cnvtsec (off_datevalue, off_uttime) | scan (off_compare)

                    time_diff = abs (on_compare - off_compare)

                    if (j == 1) {
                        t_min = time_diff
                        closest_off = off_flat
                    } else if (t_min > time_diff) {
                        t_min = time_diff
                        closest_off = off_flat
                    }
                }


                if (debug) {
                    printlog ("____OFF FLAT for: "//onflat//\
                        " is: "//closest_off, l_logfile, verbose+)
                }

                # Determine the output name
                fparse (onflat, verbose-)
                sub_flat = substr(fparse.root,1,\
                    strlstr("_GAtmpon",onflat)-1)//\
                    "_GAtmpsub"//mktemp("tmp")//".fits"

                # Store name incase of clean up
                print (sub_flat, >> imtmptodel)

                printlog ("GAFLAT: Subtracting combined OFF flat from \
                    combined ON flat...", l_logfile, l_verbose)

                # Check the DQ planes are getting updated appropriately
                gemarith (onflat, "-", closest_off, sub_flat, \
                    sci_ext=l_sci_ext, var_ext=l_var_ext, dq_ext=l_dq_ext, \
                    mdf="MDF", fl_vardq=l_fl_vardq, dims="default", \
                    intype="default", outtype="ref", refim="operand1", \
                    rangecheck=yes, verbose=no, logfile=l_logfile, \
                    glogpars="glogpars")

                # Print the names of the files that may get left over to a
                # temporary todel file
                if (gemarith.status != 0) {
                    printlog ("ERROR - GAFLAT: GEMARITH returned a non-zero \
                        status. Exiting.", l_logfile, verbose+)

                    type (imtmptodel, >> todelete)
                    type (out_oncomb_list, >> todelete)
                    type (out_offcomb_list, >> todelete)
                    goto crash
                }

                gemhedit (sub_flat//"[0]", "LOFFSUB", "yes", \
                    "Combined lamps off subtracted", delete-, upfile="")

                printlog ("        Done.\n", l_logfile, l_verbose)

                print (sub_flat, >> tmpnorm_list)

            }

            norm_list = tmpnorm_list
            delete ("@"//out_oncomb_list//",@"//out_offcomb_list, \
                verify-, > "dev$null")
            delete (out_oncomb_list//","//out_offcomb_list, \
                verify-, > "dev$null")
            print (norm_list, >> todelete)
            delete (imtmptodel, verify-, >& "dev$null")

        } else {
            # Just ON FLATs
            norm_list = out_oncomb_list
        }

        # reset nstat
        nstat = 1

        # Loop over norm_list
        scanfile2 = norm_list
        while (fscan(scanfile2, inflat) != EOF) {

            printlog ("GAFLAT: Calculating normalization factors...", \
                l_logfile, l_verbose)

            # Determine the output name
            fparse (inflat, verbose-)
            out_flat = substr(fparse.root,1,strlstr("_GAtmp",inflat))//\
                l_outsufx//".fits"

            # Call gastat
            gastat (inflat, stattype=l_stattype, statsec=l_statsec, \
                statextn=l_statextn, fl_mask=l_fl_mask, \
                badpix=bpmgemcomb, calc_stddev=yes, lower=INDEF, upper=INDEF, \
                nclip=1, lsigma=3., usigma=3., sci_ext=l_sci_ext, \
                dq_ext=l_dq_ext, logfile=l_logfile, verbose-)

            if (gastat.status != 0) {
                printlog ("ERROR - GAFLAT: GASTAT returned a \
                    non-zero status. Exiting.", \
                    l_logfile, verbose+)
                goto crash
            }

            next_statvalue = gastat.outstat
            next_stddev = gastat.stddev

            l_test = 1
            l_test2 = 1
            j = 0

            while (l_test > 0) {
                j += 1

                l_test = stridx (",",next_statvalue)
                l_test2 = stridx (",",next_stddev)

                if (l_test == 0) {
                    last_char = strlen(next_statvalue)
                } else {
                    last_char = l_test
                }

                if (l_test2 == 0) {
                    last2_char = strlen(next_stddev)
                } else {
                    last2_char = l_test2
                }

                # Parse the current value
                norm_factor = real(substr(next_statvalue,1,\
                    last_char))

                norm_err = real(substr(next_stddev,1,last2_char))

                if (debug) {

                    printlog ("____norm_factor: "//\
                        norm_factor, l_logfile, verbose+)
                }

                next_statvalue = substr(next_statvalue,\
                    l_test+1, strlen(next_statvalue))

                next_stddev = substr(next_stddev,\
                    l_test2+1, strlen(next_stddev))

                sciext = "["//l_sci_ext//","//j//"]"

                printlog ("        Normalization factor "//\
                    "for "//out_flat//sciext//": "//\
                    norm_factor, l_logfile, l_verbose)

                printlog ("        Normalization factor "//\
                    "error for "//out_flat//sciext//": "//\
                    norm_err, l_logfile, verbose-)

                gemhedit (inflat//sciext, "NORMFACT", norm_factor, \
                    "Normalization value",  delete-)

                gemhedit (inflat//sciext, "NORMFERR", norm_err, \
                    "Normalization value stddev",  delete-)
            }

            # Store name incase of clean up
            print (inflat, >> imtmptodel)

            # Normalise
            gemexpr ("a / a.NORMFACT", out_flat, inflat, \
                var_expr="(a[VAR] + (((a[SCI] * a[SCI].NORMFERR)**2)"//\
                " / ((a[SCI].NORMFACT)**2))) / ((a[SCI].NORMFACT)**2)", \
                dq_expr="a[DQ]", sci_ext=l_sci_ext, var_ext=l_var_ext, \
                dq_ext=l_dq_ext, mdf="MDF", fl_vardq=l_fl_vardq, \
                dims="default", intype="default", outtype="ref", refim="a", \
                rangecheck=yes, verbose=no, logfile=l_logfile)

            if (gemexpr.status != 0) {
                printlog ("ERROR - GAFLAT: GEMEXPR returned a non-zero \
                    status. Exiting.", l_logfile, verbose+)
                type (imtmptodel, >> todelete)
                goto crash
            }

            # Update the metaconf
            new_metaconf = "DUMMY" # Set as a warning
            for (j = 1; j <= num_types; j += 1) {

                l_test = strstr(types[j],metaconf)

                # Find the flat type
                if (l_test > 0) {
                    proc_state = ""
                    if (substr(metaconf,l_test+strlen(types[j]),\
                        l_test+strlen(types[j])) == "_") {
                        #It's been processed somehow by gareduce
                        proc_state = substr(metaconf,\
                            l_test+strlen(types[j]),\
                            strlen(metaconf))
                        proc_state = substr(proc_state,1,\
                            stridx("+",proc_state)-1)
                    }

                    if (l_test == (stridx("+",metaconf) + 1)) {
                        # Exposure time has been ignored!
                        new_metaconf = substr(metaconf,1,\
                            l_test+strlen(types[j])-1)
                    } else {
                        # Remove exposure time
                        new_metaconf = substr(metaconf,\
                            stridx("+",metaconf)+1,\
                            l_test+strlen(types[j])-1)
                    }

                    new_metaconf = new_metaconf//"_FLAT"//proc_state

                    if (strstr("OPEN",metaconf) > 0) {
                        new_metaconf = new_metaconf//substr(metaconf,\
                            strstr("OPEN",metaconf)+4,strlen(metaconf))
                    } else {
                        new_metaconf = new_metaconf//substr(metaconf,\
                            l_test+strlen(proc_state)+strlen(types[j]),\
                            strlen(metaconf))
                    }
                    break
                }
            } ##M BUG here if we allow sky flats?

            # Update the gain and units of the variance planes and remove
            # BUNIT keyword from DQ planes!
            if (l_fl_vardq) {
                keypar (out_flat//"[0]", "NSCIEXT", silent+)
                nsci = int(keypar.value)

                for (i = 1; i <= nsci; i += 1) {
                    # Update the units of the variance plane
                    keypar (out_flat//"["//l_var_ext//","//i//"]", \
                        "BUNIT", silent+)

                    if (stridx("*",str(keypar.value)) == 0) {
                        gemhedit (out_flat//"["//l_var_ext//","//i//"]", \
                            "BUNIT", keypar.value//"*"//keypar.value, "", \
                            delete-)
                    }

                    # Delete the units of the DQ plane
                    gemhedit (out_flat//"["//l_dq_ext//","//i//"]", \
                        "BUNIT", "", "", delete+)

                    # Update the title of the DQ plane
                    gemhedit (out_flat//"["//l_dq_ext//","//i//"]", \
                        "i_title", "Data Quality", "", delete-)

                }
            }

            gemhedit (out_flat//"[0]", "METACONF", new_metaconf, "", delete-)

            gemdate ()
            gemhedit (out_flat//"[0]", "GEM-TLM", gemdate.outdate,\
                "UT Last modification with GEMINI", delete-)
            gemhedit (out_flat//"[0]", "GAFLAT", gemdate.outdate, \
                "UT Time stamp for GAFLAT", delete-)

            printlog ("\nGAFLAT: Created "//out_flat, l_logfile, l_verbose)
        }

        delete (imtmptodel, verify-, >& "dev$null")

        imdelete ("@"//norm_list, verify-, >& "dev$null")

        printlog ("\nGAFLAT: Finished processing configuration "//\
            metaconf//"\n----", l_logfile, l_verbose)

        # Add to the tally of files created
        num_created += 1

NEXT_CONFIGURATION:

        imdelete (tmpbpm, verify-, >& "dev$null")

        # Counter for number of configurations for creating lists
        nconf += 1

        # Keep track if any output flats were created
        if (num_created == 0) {
            printf ("%22s"//metaconf//"\n", " ", >> tmp_skipped)
        }

        # Remove the current configuration from the tmpconflist.
        match (metaconf, tmpconflist, stop=yes, print_file_name=no, \
            metacharacter=no, > tmplist)

        delete (tmpconflist, verify-, >& "dev$null")

        # This need to be set like this for CL - MS
        scanfile = ""

        # Remove off flats if required
        if (dooffs) {
            match (off_metaconf, tmplist, stop=yes, print_file_name=no, \
                metacharacter=no, > tmpconflist)

            delete (tmplist, verify-, >& "dev$null")
        } else {
            rename (tmplist, tmpconflist, field="all")
        }

        scanfile = tmpconflist

        total_num_created += num_created

    } # End of loop over unique configurations

    goto clean

    #--------------------------------------------------------------------------
crash:

    # Exit with error subroutine
    status = 1

clean:

    if (access(imtmptodel)) {
        delete (imtmptodel, verify-, >& "dev$null")
    }

    imdelete (tmpbpm, verify-, >& "dev$null")

    if (access(outgaimchk_list)){
        delete ("@"//outgaimchk_list//","//outgaimchk_list, \
            verify-, >& "dev$null")
    }

    delete (filelist, verify-, >& "dev$null")

    if (access(todelete)) {
        delete ("@"//todelete, verify-, >& "dev$null")
        delete (todelete, verify-, >& "dev$null")
    }

    scanfile = ""
    scanfile2 = ""

    # Check for configurations that were skipped and if any fils were created
    if (status == 0) {
        if (access(tmp_skipped)) {
            printlog ("\nWARNING - GAFLAT: The following configurations \
                were not processed - \n", l_logfile, verbose+)
            type (tmp_skipped) | tee (l_logfile, out_type="text", append+)
            delete (tmp_skipped, verify-, >& "dev$null")

            if (total_num_created == 0) {
                printlog ("\nERROR - GAFLAT: No output files created", \
                    l_logfile, verbose+)
                status = 1
            }
        }
    }

    # Print finish time
    gemdate(zone="local")
    printlog ("\nGAFLAT -- Finished: "//gemdate.outdate, \
        l_logfile, l_verbose)

    if (status == 0) {
        printlog ("\nGAFLAT -- Exit status: GOOD", l_logfile, l_verbose)

    } else {
        printlog ("\nGAFLAT -- Exit status: ERROR", l_logfile, l_verbose)
    }

    if (status != 0) {
        printlog ("       -- Please read the logfile \""//l_logfile//\
                "\" for more information.", l_logfile, l_verbose)
    }

    printlog ("------------------------------------------------------------", \
        l_logfile, l_verbose)

end




# Copyright(c) 2010-2015 Association of Universities for Research in Astronomy, Inc.

procedure gareduce (inimages)

##M Add CCDSEC and DATASEC parameters?
##M Remove creation of temporary database and related variables

# Takes GSAOI raw or gaprepared or reduced data, 4 extensions
# This task will, eventually
# - take a list of science frames already GPREPAREd or run it if the keyword is
#   not there
# - then, on each extension
# - subtract dark if required (check exposure time)
# - subtract the sky using offset+time or running average of dithered positions
#   or a sky file
# - divide by the proper flat (check filter)
# - convert each extension to ADU by multiplying by the gain
# - when it does that, it should add the keyword GAINMULT to the header of the
#   extension and set the value of the GAIN keyword to 1.
# - add to the header GAREDUCE time stamp, flat, dark and sky file names as
#   needed

# Now assuming that a typical GSAOI sequence will be filter - some type of
# dither sequence:
#   filter - dither sequence - filter - etc.
# - gaprepare
# - create the lists for the different filter+exposure times
# - dark subtract first if required - note that the way this is written
#   the dark current is removed either with the dark OR with the sky
#   so if building the sky from the images, dark subtraction is turned off
# - skysubtract - this can be done several different ways
#   - use a pre-generated sky frame
#   - create a sky frame from the sequence using time as parameter (dithered
#     sequences for point or point-like sources) and re-doing the sky frame for
#     each science frame (running average)
#   - create a sky frame based on offset distance (for extended objects),
#     redoing the sky frame for each science frame based on time interval (if
#     possible)
# - then flat field
# - adds back a constant sky level or uses the median value of the sky frame
# - it uses gemoffsetlist to build the list of sky images
# - this tasks does not mosaic and does not stack
# - most of the work is in the sky creation (second option above).
# - dark subtraction is needed for the hot pixels, not as much as for the dark
#   current itself (otherwise a significant number of pixels needs to be masked
#   in the BPM), so it should go away with the sky subtraction

# Version 2010 Oct 14  V0.0 CW - created
#         2010 Nov 19  V0.5 CW - mostly done, not tested
#         2011 Sep 26  V1.0 CW - tested without gacalfind; adding back median
#                                of sky frame as a constant after processing.
#         2011 Oct 03  V1.1 CW - added ROI for all frames with consistency
#                                checks
#         2011 Dec 06  V1.2 CW - changed sky creation to use gasky (copy of
#                                nisky)
#         2012 Feb 29  V1.3 CW - initial implementation of VAR/DQ
#         2012 Apr 10  V1.4 CW - consistency checks for calibrations,
#                                implemented code to allow no-processing (gain
#                                multiplication only).
#         2012 Apr 11  V1.5 CW - activated automatic calibration finding
#
# Descriptions of future revisions are recorded in CVS logs

######################################################

# Inputs
char    inimages      {prompt="Input GSAOI images"}
char    rawpath       {"",prompt="Path for input images"}
char    outpref       {"r",prompt="Prefix for output processed images"}
char    rootname      {"",prompt="Root name for images, blank=today UT"}
# Flags and calibration info
bool    fl_dark       {no,prompt="Subtract dark image?"}
bool    fl_flat       {no,prompt="Do flat field correction?"}
bool    fl_sky        {no,prompt="Subtract sky image?"}
bool    fl_autosky    {yes,prompt="Use median of the sky to add back constant to processed frame"}
bool    fl_mult       {yes,prompt="Multiply by gain to convert to electrons?"}
bool    fl_vardq      {no,prompt="Create variance and data quality frames"}
char    darkimg       {"",prompt="Dark image to be used (\"\"|find|<filename>)"}
char    flatimg       {"",prompt="Flat field image to be used (\"\"|find<{DOME|TWLT|GCAL}>|<filename>)"}
char    skyimg        {"",prompt="Sky image to be used (\"\"|time|distance|both|<filename>)"}
char    calpath       {"./",prompt="Path for calibration images"}
char    caltable      {"gsaoical.fits",prompt="Name of table with calibration information (from GACALFIND)"}
bool    fl_calrun     {no,prompt="When looking for calibrations, first re-create existing table"}
int     cal_maxtime   {INDEF,prompt="Maximum time difference between calibration and observation"}
char    badpix        {"gsaoi$data/gsaoibpm_high_full.fits",prompt="Static bad pixel mask"}
real    lsigma        {3.,min=0,prompt="Lower sigma clipping factor when fl_autosky=yes"}
real    hsigma        {3.,min=0,prompt="Upper sigma clipping factor when fl_autosky=yes"}
# gaprepare
char    gaprep_pref   {"g",prompt="Prefix for gaprepare output images"}
bool    fl_trim       {yes,prompt="Trim the images?"}
bool    fl_nlc        {yes,prompt="Apply non-linear correction to each array?"}
bool    fl_sat        {yes,prompt="Include non-linear and saturated pixels in DQ"}
char    arraysdb      {"gsaoi$data/gsaoiAMPS.dat",prompt="Database file containing array information"}
char    non_lcdb      {"gsaoi$data/gsaoiNLC.dat",prompt="Database file containing non-linearity coefficients"}
# Combine for sky frame - using gasky: values may need tweaking
char    sky_sufx      {"sky",prompt="Output suffix if creating sky frames from inputs"}
real    maxtime       {900.,prompt="For sky=\"time\", combine frames within this time interval"}
real    minoffs       {90.,prompt="For sky=\"distance\", combine frames with offsets larger than this"}
int     minsky        {5,prompt="Minimum number of sky images to combine"}
char    combine       {"default",enum="default|average|median",prompt="Combination operation for sky frames (default|average|media)"}
char    reject        {"minmax",enum="none|minmax|avsigclip",prompt="Rejection algorithm for sky frames (none|minmax|avsigclip)"}
int     nlow          {1,min=0,prompt="minmax: Number of low pixels to reject"}
int     nhigh         {1,min=0,prompt="minmax: Number of high pixels to reject"}
char    statsec       {"[5%]",prompt="Statistics section"}
char    masktype      {"goodvalue",enum="none|goodvalue",prompt="Bad Pixel Mask type (none|goodvalue)"}
real    maskvalue     {0.,prompt="Good pixel value in the BPM"}
char    datename      {"DATE-OBS",prompt="Date header keyword"}
char    timename      {"UT",prompt="Time stamp header keyword"}
char    expname       {"EXPTIME",prompt="Exposure time header keyword"}
char    key_ron       {"RDNOISE",prompt="Keyword for readout noise in e-"}
char    key_gain      {"GAIN",prompt="Keyword for gain in electrons/ADU"}
real    ron           {0.0,min=0.,prompt="Readout noise rms in electrons"}
real    gain          {1.0,min=0.00001,prompt="Gain in e-/ADU"}
real    lthreshold    {3.,min=1.5,prompt="Threshold in sigma for object detection"}
int     ngrow         {3,prompt="Number of iterations to grow objects into the wings"}
real    agrow         {3.,prompt="Area limit for growing objects into the wings"}
int     minpix        {6,prompt="Minimum number of pixels to be identified as an object"}
bool    fl_mask       {yes,prompt="Mask bad pixels during source detection and when fl_autosky=yes"}
bool    qreduce_gasky {yes,prompt="Quickly reduce images, with sky from GAFASTSKY for object masking?"}
real    minflat_val   {INDEF,prompt="Minimum allowed pixel value in flat image"}
bool    fl_keepsky    {no,prompt="Keep the generated sky frames"}
#bool    fl_keepmasks  {no,prompt="Keep object masks for each input image?"}
#char    masksuffix    {"msk",prompt="Mask name suffix"}
bool    fl_dqprop     {yes,prompt="Retain sky DQ information in output sky frame?"}
# Other stuff
char    sci_ext       {"SCI",prompt="Name of science extensions"}
char    var_ext       {"VAR",prompt="Name of variance extensions"}
char    dq_ext        {"DQ",prompt="Name of data quality extensions"}
char    logfile       {"",prompt="Logfile"}
bool    verbose       {yes,prompt="Verbose?"}
int     status        {0,prompt="Exit status (0=good)"}
struct  *scanfile     {"",prompt="Internal use only"}
struct  *scanfile2    {"",prompt="Internal use only"}
struct  *scanfile3    {"",prompt="Internal use only"}

######################################################

begin

    ####
    # Variable declaration
    char l_inimages, l_rawpath, l_outpref, l_rootname, l_datename
    char l_timename, l_calpath, l_darkimg, l_flatimg, l_skyimg
    char l_badpix, l_combine, l_reject, l_masktype, l_scale, l_zero
    char l_weight, l_statsec, l_expname, l_key_ron, l_key_gain, l_snoise
    char l_logfile, l_caltable, l_temp, filelist, utdate, inimg, metaconf
    char tmpkeeplst, darkfile, flatfile, selexpr, skyname, skyfile, skylist
    char skysuf, tmpskylst, tmpdiflst, tmplst, skylevel, varexpr, dqexpr
    char l_gaprep_pref, l_key_allowed, l_key_forbidden, l_object_allowed
    char l_obstype_allowed, tmpfile, tmpconflist, outgaimchk_list, l_sci_ext
    char l_var_ext, l_dq_ext, nextlist, out_metaconf, imgtype, meta_proc
    char testflat1, testdark1, proc_test, tmpflttab, tmpflat_caltab, tmp2file
    char tmpdark_caltab, tmptabflt, tmpgaprep, tmplist, inphu, conflist
    char inname, uttime, val_date, cal_table, tmpproc_test, l_non_lcdb
    char skymeta1, skymeta2, testsky1, l_pwd, img_root, l_masksuffix
    char a_operand, b_operand, c_operand, e_operand, f_operand, tmpnotproc
    char d_operand, gexpr, calcextn[3], tmpsky, tmpdark, tmpflat, l_units
    char darkexpr, tmpexprdb, outphu, l_out_key_gain, l_arraysdb, sum_str
    char skyexpr, flatexpr, outname, outextn, outappend, s_obj_exp
    char skyconst_expr, sciexpr, t_string, l_redirect, d_expr, l_key_froz_odgw
    char l_darkimg_orig, l_flatimg_orig, l_skyimg_orig, inflat_gasky
    char auto_skyval, auto_orig_skyval, sky_comment, skyvarexpr, testfilephu
    char darkgain, skygain, testfile, testoutfile, tmpout, auto_val, gemout
    char l_fcalpath, tmpfcal, tmpdcal, l_dcalpath, chkincal, chkpath
    char chktab, gconstexpr, darkvarexpr, skyflatimg, darkname, flatname
    char tmpgcorr_sky, asky_file, csky_type, skymidpt[4], auto_sky_stat
    char full_skymeta, nscicoadds, nskycoadds, csky_flat, l_sky_sufx
    char scirep_val, varrep_val, dqrep_val, tmpkeepskylst, pr_str
    char tmpskyfile, tmpskyinnames, chkcaltype

    int l_minsky, l_nlow, l_nhigh, l_minpix, l_ngrow, l_cal_maxtime
    int i, total_num_created, num_created, num_counter, obs_time, lcount, t_min
    int tvalmin, npathcols, num_calc, star_pt, neg_pt, pos_skyconst, nsci
    int l_test, ii, jj, ntotal, nconf, exist, nrows, nsky1, nsky2, last_char
    int ndarksci, nskysci, nflatsci, sarray_index, im_type_len

    real l_maxtime, l_minoffs, l_maskvalue, l_lthreshold, l_hthreshold
    real l_lsigma, l_hsigma, l_ron, l_gain, l_sigscale, l_pclip, l_agrow
    real l_minflat_val, l_ingain, sky_exp, obj_exp, exp_tolerance
    real sdistance, stime

    bool l_fl_dark, l_fl_flat, l_fl_sky, l_fl_vardq, l_verbose
    bool l_fl_calrun, l_fl_keepsky, l_fl_autosky, l_fl_dqprop, l_fl_trim
    bool l_fl_nlc, l_fl_sat, l_fl_mult, ignore_cal_nlc, processed, candark
    bool canflat, cansky, cangainmult, isgainmult, skip_if_skysub, isdarksub
    bool isflat, isskysub, doskycomb, l_fl_mask, l_fl_keepmasks, someprepared
    bool l_qreduce, dovardq, dark_isgainmult, sky_isgainmult, isdarkflat
    bool isskyflat, isskydark, isflatdark, debug, caldebug, update_metaconf
    bool first, fl_usersky, fixed_gemexpr, chkcal, skyscaled, determine_sky
    bool ignore_coadds, allow_tolerance, unmatched_coadds, alt_sky
    bool snearer, syounger

    struct l_time_struct

    # There is a hiercahy to this list: upto and including TWLT are added by
    # gaprepapre with the rest being added by other tasks - MS
    char poss_imgtype[7]="OBJ","DOME","GCAL","TWLT","DARK","SKY","FLAT"
    int num_imgtypes=7

    ####
    # Temporary files
    filelist = mktemp("tmpfile")
    skylist = mktemp("skylist")
    outgaimchk_list = mktemp ("tmpgaimchk")
    nextlist = mktemp ("tmpnextlist")
    tmpkeeplst = mktemp("tmpkeep")
    tmpkeepskylst = mktemp("tmpkeepskylst")
    tmpconflist = mktemp ("tmptmpconflist")
    tmpgaprep = mktemp ("tmpgaprep")
    tmpflat_caltab = mktemp("tmpflat_caltab")//".fits"
    tmpdark_caltab = mktemp ("tmpdark_caltab")//".fits"
    tmptabflt = mktemp ("tmptabflt")
    tmplist = mktemp ("tmplist")
    conflist = mktemp("tmpconflist")
    tmpexprdb = mktemp ("tmpexprdb")
    tmpnotproc = mktemp("tmpnotproc")
    tmpdark = mktemp ("tmpdark")
    tmpflat = mktemp ("tmpflat")
    tmpsky = mktemp ("tmpsky")
    tmpout = mktemp ("tmpout")//".fits"
    tmpfcal = mktemp ("tmpfcal")//".fits"
    tmpdcal = mktemp ("tmpdcal")//".fits"
    tmpgcorr_sky = mktemp ("tmpautosky")//".fits"
    tmpskyinnames = mktemp("tmpskyinnames")
    tmpdiflst = mktemp("tmpdiflst")
    tmplst = mktemp("tmplst")

    ####
    # Local variables
    l_inimages = inimages
    l_rawpath = rawpath
    l_outpref = outpref
    l_rootname = rootname
    l_datename = datename
    l_timename = timename
    l_calpath = calpath
    l_darkimg = darkimg
    l_flatimg = flatimg
    l_skyimg = skyimg
    l_badpix = badpix
    l_combine = combine
    l_reject = reject
    l_statsec = statsec
    l_expname = expname
    l_key_ron = key_ron
    l_key_gain = key_gain
    l_logfile = logfile
    l_caltable = caltable
    l_minsky = minsky
    l_nlow = nlow
    l_nhigh = nhigh
    l_minpix = minpix
    l_sky_sufx = sky_sufx
    l_maxtime = maxtime
    l_minoffs = minoffs
    l_lthreshold = lthreshold
    l_lsigma = lsigma
    l_hsigma = hsigma
    l_ron = ron
    l_gain = gain
    l_fl_dark = fl_dark
    l_fl_flat = fl_flat
    l_fl_sky = fl_sky
    l_fl_vardq = fl_vardq
    l_verbose = verbose
    l_fl_calrun = fl_calrun
    l_fl_keepsky = fl_keepsky
    l_fl_autosky = fl_autosky
    l_ngrow = ngrow
    l_agrow = agrow
    l_fl_dqprop=fl_dqprop
    l_fl_trim = fl_trim
    l_fl_nlc = fl_nlc
    l_fl_sat = fl_sat
    l_cal_maxtime = cal_maxtime
    l_arraysdb = arraysdb
    l_non_lcdb = non_lcdb
#    l_fl_keepmasks = fl_keepmasks
#    l_masksuffix = masksuffix
    l_fl_mask = fl_mask
    l_qreduce = qreduce_gasky
    l_fl_mult = fl_mult
    l_gaprep_pref = gaprep_pref
    l_sci_ext = sci_ext
    l_var_ext = var_ext
    l_dq_ext = dq_ext
    l_masktype = masktype
    l_maskvalue = maskvalue

    # Set the redirect when typing and teeing lists
    if (l_verbose) {
        l_redirect = "STDOUT"
    } else {
        l_redirect = "dev$null"
    }

    ####
    # Default values
    debug = no
    caldebug = no
    status = 0
    someprepared = no
    determine_sky = no

    # This used to be a parameter; only used when using science inputs to
    # create sky frames. The masks are likely to change due to the way the
    # files are processed in GASKY. So, don't allow user to keep
    # masks. Information on the location of sources if they affect the output
    # skyframe are stored in the DQ plane if fl_vardq+ and fl_dqprop - MS
    l_fl_keepmasks = no
    l_masksuffix = "msk"

    # Set this whether to add the output processing state of the image to the
    # configuration keywprd - adds '_' then 'd' for dark subtracted, 's' for
    # sky subtracted, then 'f' for flat fielded, and g for gain multiplied
    # after the image type
    update_metaconf = yes

    # Set this to ignore NLC state of calibrations - is parameter to gacalfind
    ignore_cal_nlc = no

    # Set this whether to skip the image if its already been sky subtracted
    skip_if_skysub = no

    # Change this to yes if the gemexpr exprdb bug gets fixed
    fixed_gemexpr = yes

    # Set these two parameters to determine how the VAR and DQ values are
    # replaced if using minflat_val
#    varrep_val = "0"
#    dqrep_val = "256"
    varrep_val = "d[VAR]"
    dqrep_val = "d[DQ]"

    # Keyword in skyfile used to say ODGWs frozen - only it's presence not what
    # it's value is say they are frozen
    l_key_froz_odgw = "SKYFODGW"

    # Set this to ignore the number of coadds in the object frame and the sky
    # frame. NOTE: This will cause the incorrect dark subtraction if using sky
    # to dark subtract and may cause incorrect noise calculations.
    ignore_coadds = yes

    # There could be rounding differences in the exposure times so lets allow a
    # tolerance between sky frame and object frames
    allow_tolerance = yes
    exp_tolerance = 0.5

    # If we want to allow sky frames to be scaled set this appropriately later
    # in the code
    skyscaled = no

    ########
    # Here is where the actual work starts

    cache ("gemexpr", "gastat", "gadimschk", "gasky", "gacalfind")

    # Test the log file
    if ((l_logfile == "") || (stridx(" ", l_logfile) > 0)) {
       l_logfile = gsaoi.logfile
       if ((l_logfile == "") || (stridx (" ", l_logfile) > 0)) {
          l_logfile = "gsaoi.log"
          printlog ("WARNING - GAREDUCE: Both gareduce.logfile and \
              gsaoi.logfile fields are empty", l_logfile, l_verbose)
          printlog ("                     Using default file gsaoi.log", \
              l_logfile, l_verbose)
       }
    }

    # Print start time
    printlog ("", l_logfile, verbose+)
    gemdate (zone="local")
    printlog ("------------------------------------------------------------", \
        l_logfile, l_verbose)
    printlog ("GAREDUCE -- Started: "//gemdate.outdate//"\n", \
        l_logfile, l_verbose)

    # Record all input parameters
    printlog ("GAREDUCE: Input parameters...\n", l_logfile, l_verbose)
    printlog ("    inimages    = "//l_inimages, l_logfile, l_verbose)
    printlog ("    outpref     = "//l_outpref, l_logfile, l_verbose)
    printlog ("    rootname    = "//l_rootname, l_logfile, l_verbose)
    printlog ("    rawpath     = "//l_rawpath, l_logfile, l_verbose)
    printlog ("    datename    = "//l_datename, l_logfile, l_verbose)
    printlog ("    timename    = "//l_timename, l_logfile, l_verbose)
    printlog ("    fl_dark     = "//l_fl_dark, l_logfile, l_verbose)
    printlog ("    fl_flat     = "//l_fl_flat, l_logfile, l_verbose)
    printlog ("    fl_sky      = "//l_fl_sky, l_logfile, l_verbose)
    printlog ("    fl_autosky  = "//l_fl_autosky, l_logfile, l_verbose)
    printlog ("    fl_mult     = "//l_fl_mult, l_logfile, l_verbose)
    printlog ("    fl_vardq    = "//l_fl_vardq, l_logfile, l_verbose)
    printlog ("    darkimg     = "//l_darkimg, l_logfile, l_verbose)
    printlog ("    flatimg     = "//l_flatimg, l_logfile, l_verbose)
    printlog ("    skyimg      = "//l_skyimg, l_logfile, l_verbose)
    printlog ("    calpath     = "//l_calpath, l_logfile, l_verbose)
    printlog ("    caltable    = "//l_caltable, l_logfile, l_verbose)
    printlog ("    fl_calrun   = "//l_fl_calrun, l_logfile, l_verbose)
    if (isindef(l_cal_maxtime)) {
        printlog ("    cal_maxtime = INDEF", l_logfile, l_verbose)
    } else {
        printlog ("    cal_maxtime = "//l_cal_maxtime, l_logfile, l_verbose)
    }
    printlog ("    badpix      = "//l_badpix, l_logfile, l_verbose)
    printlog ("    gaprep_pref = "//l_gaprep_pref, l_logfile, l_verbose)
    printlog ("    fl_trim     = "//l_fl_trim, l_logfile, l_verbose)
    printlog ("    fl_nlc      = "//l_fl_nlc, l_logfile, l_verbose)
    printlog ("    fl_vardq    = "//l_fl_vardq, l_logfile, l_verbose)
    printlog ("    fl_sat      = "//l_fl_sat, l_logfile, l_verbose)
    printlog ("    arraysbd    = "//l_arraysdb, l_logfile, l_verbose)
    printlog ("    non_lcdb    = "//l_non_lcdb, l_logfile, l_verbose)
    printlog ("    sci_ext     = "//l_sci_ext, l_logfile, l_verbose)
    printlog ("    var_ext     = "//l_var_ext, l_logfile, l_verbose)
    printlog ("    dq_ext      = "//l_dq_ext, l_logfile, l_verbose)
    printlog ("    logfile     = "//l_logfile, l_logfile, l_verbose)

    ####
    # Check parameters...

    # Create the expr database used by imexpr
    if (!isindef(minflat_val)) {
        # BUG in IRAF (2.14 min) that means that cannot use str() function in
        # first argument to fscan. Also BUG in PyRAF (v1.10) that means you
        # have to use a variable that is a string or char as the first argument
        # - MS
        t_string = str(minflat_val)
        l_test = fscan (t_string,l_minflat_val)
        if (l_test != 1) {
            printlog ("ERROR - GAREDUCE: minflat_val is not a number. Exiting")
            goto crash
        } else {
            # Create a temporary expr databae for imexpr
            # Basically checking for near zero numbers
            print ("chkval(a)        ((abs(a) < "//l_minflat_val//") ? "//\
                "((a < 0) ? (-"//l_minflat_val//") : "//l_minflat_val//\
                ") : a)", > tmpexprdb)
            printlog ("    minflat_val = "//l_minflat_val, \
                l_logfile, l_verbose)
        }
        scirep_val = str(l_minflat_val)
    } else {
        # Create a temporary expr databae for imexpr
        print ("chkval(a)        (a)", > tmpexprdb)
        l_minflat_val = minflat_val
        scirep_val = "0"
    }
    printlog ("", l_logfile, l_verbose)

    if (debug) type (tmpexprdb)

    # Check outpref
    if (l_outpref == "" || stridx(" ",l_outpref) > 0) {
        printlog ("WARNING - GAREDUCE: outpref is an empty string or contains \
            spaces. Setting to default \"r\"", l_logfile, verbose+)
        l_outpref = "r"
    }

    # Check gaprepare out prefix
    if (l_gaprep_pref == "" || stridx(" ",l_gaprep_pref) > 0) {
        printlog ("WARNING - GAREDUCE: outpref is an empty string or contains \
            spaces. Setting to default \"g\"", l_logfile, verbose+)
        l_gaprep_pref = "g"
    }

    if (l_fl_vardq) {
        if (l_dq_ext == "" || (stridx(" ",l_dq_ext) > 0) || \
            l_var_ext == "" || (stridx(" ",l_var_ext) > 0)) {

            printlog("WARNING - GAREDUCE: var_ext or dq_ext have not been \
                set.\n                     Output image will not have VAR or \
                DQ planes.", l_logfile,verbose+)
            l_fl_vardq=no
        }
    }

    # Check the l_badpix value
    if ((l_badpix == "") || (stridx(" ", l_badpix) > 0)) {
        # This assumes that the bpm covers all of the areas readout in each
        # detector in each image
        l_badpix = "none"
    } else if (!imaccess(l_badpix)) {
        printlog ("ERROR - GAREDUCE: cannot access bpm "//l_badpix//\
            ". Exiting", l_logfile, verbose+)
    } else {
        # Check it hasn't been mosaiced
        keypar (l_badpix//"[0]", "GAMOSAIC", silent+)
        if (keypar.found) {
            printlog ("ERROR - GAREDUCE: "//l_badpix//" has been mosaiced."//\
                " Please supply a static BPM that has not been mosaiced.", \
                l_logfile, verbose)
            goto crash
        }
    }
    printlog ("GAREDUCE: badpix value after checks is "//l_badpix, \
        l_logfile, l_verbose)

    # Check that the caltable has an extension ".fits"
    fparse (l_caltable)
    if (fparse.extension != ".fits") {
        printlog ("WARNING - GAREDUCE: Setting caltable extension to \
            \".fits\"", l_logfile, verbose+)
        l_caltable = fparse.directory//fparse.root//".fits"
    }

    # Check the setting for darkimg and fl_dark
    if (l_fl_dark) {
        if (l_darkimg == "" || stridx(" ",l_darkimg) > 0) {
            printlog ("ERROR - GAREDUCE: fl_dark=yes but darkimg is empty or \
                contains spaces. "//\
                "\n                  Please supply a filename or set \
                darkimg=find. Exiting.", l_logfile, verbose+)
            goto crash
        } else if (l_darkimg != "find") {
            # Check the input file exists
            if (!imaccess(l_darkimg)) {

                if (!imaccess(l_calpath//l_darkimg)) {

                    printlog ("ERROR - GAREDUCE: Cannot access dark "//\
                        " image \""//l_darkimg//"\". Exiting", \
                        l_logfile, verbose+)
                    goto crash
                } else {
                    l_dcalpath = l_calpath
                }
            } else {
                l_dcalpath = "./"
            }

            # Check it's actually a dark
            keypar (l_dcalpath//l_darkimg//"[0]", "GADARK", silent+)
            if (!keypar.found) {
                printlog ("ERROR - GAREDUCE: Dark image \""//l_darkimg//\
                    "\" has not been processed by GADARK. Exiting", \
                    l_logfile, verbose+)
                goto crash
            } else {
                # Check it's not mosaiced
                keypar (l_dcalpath//darkimg//"[0]", "GAMOSAIC", silent+)
                if (keypar.found) {
                    printlog ("ERROR - GAREDUCE: Dark image \""//\
                        l_darkimg//"\" has been mosaiced. Exiting", \
                        l_logfile, verbose+)
                    goto crash
                }
            }

            # Send it off to gacalfind for checking later
            gacalfind (inimages=l_darkimg, calpath=l_dcalpath, \
                caltable=tmpdcal, fl_calrun=yes, fl_find=no, caltype="", \
                sciimg="", maxtime=INDEF, ignore_nlc=ignore_cal_nlc, \
                sci_ext=l_sci_ext, var_ext=l_var_ext, dq_ext=l_dq_ext, \
                logfile=l_logfile, verbose=no)

            if (gacalfind.status != 0) {
                printlog ("ERROR - GAREDUCE: GACALFIND returned a non-zero \
                    status. Exiting.", l_logfile, verbose+)
                goto crash
            } else if (!access(tmpdcal)) {
                printlog ("ERROR - GAREDUCE: Cannot access \""//tmpdcal//\
                    "\". Exiting.", l_logfile, verbose+)
                goto crash
            }
        }
    } else {
        l_darkimg = INDEF
    } # End of fl_dark checks

    # Check the setting for flatimg and fl_flat
    if (l_fl_flat) {
        if (l_flatimg == "" || stridx(" ",l_flatimg) > 0) {
            printlog ("ERROR - GAREDUCE: fl_flat=yes but flatimg is empty or \
                contains spaces. "//\
                "\n                  Please supply a filename or set \
                flatimg=find. Exiting.", l_logfile, verbose+)
            goto crash
        } else if (strstr("find", l_flatimg) == 0) {
            # Check the input file exists
            if (!imaccess(l_flatimg)) {

                # Check the calpath
                if (!imaccess(l_calpath//l_flatimg)) {
                    printlog ("ERROR - GAREDUCE: Cannot access flat "//\
                        "image \""//l_flatimg//"\". Exiting", \
                         l_logfile, verbose+)
                    goto crash
                } else {
                    l_fcalpath = l_calpath
                }
            } else {
                l_fcalpath = "./"
            }

            # Check it's actually a flat
            keypar (l_fcalpath//l_flatimg//"[0]", "GAFLAT", silent+)
            if (!keypar.found) {
                printlog ("ERROR - GAREDUCE: Flat image \""//l_flatimg//\
                    "\" has not been processed by GAFLAT. Exiting", \
                    l_logfile, verbose+)
                goto crash
            } else {
                # Check it's not mosaiced
                keypar (l_fcalpath//l_flatimg//"[0]", "GAMOSAIC", silent+)
                if (keypar.found) {
                    printlog ("ERROR - GAREDUCE: Flat image \""//\
                        l_flatimg//"\" has been mosaiced. Exiting", \
                        l_logfile, verbose+)
                    goto crash
                }
            }

            # Send it off to gacalfind for checking later
            gacalfind (inimages=l_flatimg, calpath=l_fcalpath, \
                caltable=tmpfcal, fl_calrun=yes, fl_find=no, caltype="", \
                sciimg="", maxtime=INDEF, ignore_nlc=ignore_cal_nlc, \
                sci_ext=l_sci_ext, var_ext=l_var_ext, dq_ext=l_dq_ext, \
                logfile=l_logfile, verbose=debug)

            if (gacalfind.status != 0) {
                printlog ("ERROR - GAREDUCE: GACALFIND returned a non-zero \
                    status. Exiting.", l_logfile, verbose+)
                goto crash
            } else if (!access(tmpfcal)) {
                printlog ("ERROR - GAREDUCE: Cannot access \""//tmpfcal//\
                    "\". Exiting.", l_logfile, verbose+)
                goto crash
            }
        }
    } else {
        l_flatimg = ""
    } # End of fl_flat checks

    # Check the setting for skyimg and fl_sky
    if (l_fl_sky) {
        if (l_skyimg == "" || stridx(" ",l_skyimg) > 0) {
            printlog ("ERROR - GAREDUCE: fl_sky=yes but skyimg is empty or \
                contains spaces. "//\
                "\n                  Please supply a filename or set \
                skyimg=time or skyimg=distamnce. Exiting.", \
                l_logfile, verbose+)
            goto crash
        } else if (l_skyimg != "time" && l_skyimg != "distance" \
            && l_skyimg != "both") {

            # Check the input file exists
            if (!imaccess(l_skyimg)) {
                printlog ("ERROR - GAREDUCE: Cannot access sky image \""//\
                    l_skyimg//"\". Exiting", l_logfile, verbose+)
                goto crash
            } else {
                # Check it's been processed by gasky
                keypar (l_skyimg//"[0]", "GASKY", silent+)
                if (!keypar.found) {
                    printlog ("ERROR - GAREDUCE: Sky image \""//l_skyimg//\
                        "\" has not been processed by GASKY. Exiting", \
                        l_logfile, verbose+)
                    goto crash
                } else {
                    # Check it's not mosaiced
                    keypar (l_skyimg//"[0]", "GAMOSAIC", silent+)
                    if (keypar.found) {
                        printlog ("ERROR - GAREDUCE: Sky image \""//\
                            l_skyimg//"\" has been mosaiced. Exiting", \
                            l_logfile, verbose+)
                        goto crash
                    }
                }
            }
        } else {
            determine_sky = yes
        }
    } else {
        l_skyimg = ""
    } # End of fl_sky checks

    # Check input images...
    printlog ("\nGAREDUCE: Calling GAIMCHK to check input files...", \
        l_logfile, l_verbose)

    # In decending processing steps... keep them this way for GAIMCHK - -MS
    l_key_allowed = "GAREDUCE,GAPREPAR"

    # Cannot be mosaiced
    l_key_forbidden = "GAMOSAIC"

    l_obstype_allowed = ""
    l_object_allowed = ""

    # Call gaimchk to perform input checks
    gaimchk (inimages=l_inimages, rawpath=l_rawpath, rootname=l_rootname, \
        obstype_allowed=l_obstype_allowed, object_allowed=l_object_allowed, \
        key_allowed=l_key_allowed, key_forbidden=l_key_forbidden, \
        key_exists="", fl_prep_check=yes, gaprep_pref=l_gaprep_pref, \
        fl_redux_check=yes, garedux_pref=l_outpref, fl_fail=no, \
        fl_out_check=no, fl_vardq_check=no, sci_ext=l_sci_ext, \
        var_ext=l_var_ext, dq_ext=l_dq_ext, outlist=outgaimchk_list, \
        logfile=l_logfile, verbose=no)

    if (gaimchk.status != 0) {
        printlog ("ERROR - GAREDUCE: GAIMCHK returned a non-zero status. \
            Exiting.", l_logfile, verbose+)
        goto crash
    } else if (!access(outgaimchk_list)) {
        printlog ("ERROR - GAREDUCE: Cannot access output list from \
            GAIMCHK", l_logfile, verbose+)
        goto crash
    }

    printlog ("GAREDUCE: Returned from GAIMCHK.\n", l_logfile, l_verbose)

    # Files that have have already been processed by GAREDUCE
    tmpfile = ""
    match ("tmpGAREDUCE", outgaimchk_list, stop=no, print_file_name=no, \
        metacharacter=no) | scan (tmpfile)
    if (tmpfile != "") {
        printlog ("WARNING - GAREDUCE: The following files have already been \
            processed by GAREDUCE.\n", \
            l_logfile, verbose+)
        fields (tmpfile, fields=1, lines="1-", quit_if_miss=no, \
            print_file_names=no) | \
            tee (l_logfile, out_type="text", append+)

        fields (tmpfile, fields=1, lines="1-", quit_if_miss=no, \
            print_file_names=no, >> nextlist)
    }

    # Files that have been processed by GAPREPARE
    someprepared = no
    tmpfile = ""
    match ("tmpGAPREPAR", outgaimchk_list, stop=no, print_file_name=no, \
        metacharacter=no) | scan (tmpfile)
    if (tmpfile != "") {
        someprepared = yes
        tmp2file = tmpfile

        fields (tmpfile, fields=1, lines="1-", quit_if_miss=no, \
            print_file_names=no, >> nextlist)
    }

    # Files that need to be prepared
    tmpfile = ""
    match ("tmpNotFound", outgaimchk_list, stop=no, print_file_name=no, \
        metacharacter=no) | scan (tmpfile)
    if (tmpfile != "") {

        fields (tmpfile, fields=1, lines="1-", quit_if_miss=no, \
            print_file_names=no, >> tmpgaprep)

        # Now check if there were any GAPREPAREd files too (not GAREDUCEd
        # though and print warning if there are.
        if (someprepared) {
            printlog ("WARNING - GAREDUCE: The following files have \
                already been processed by GAPREPARE.\n", \
                l_logfile, verbose+)
            fields (tmp2file, fields=1, lines="1-", quit_if_miss=no, \
                print_file_names=no) | \
                tee (l_logfile, out_type="text", append+)
        }

        printlog ("\nGAREDUCE: Calling gaprepare to process unprepared input \
            files", l_logfile, l_verbose)

        gaprepare ("@"//tmpgaprep, rawpath="",outpref=l_gaprep_pref, \
            rootname="", fl_trim=l_fl_trim, fl_nlc=l_fl_nlc, \
            fl_vardq=l_fl_vardq, fl_sat=l_fl_sat, badpix=l_badpix, \
            sci_ext=l_sci_ext, var_ext=l_var_ext, dq_ext=l_dq_ext, \
            arraysdb=l_arraysdb, non_lcdb=l_non_lcdb, \
            logfile = l_logfile, verbose=l_verbose)

        if (gaprepare.status != 0) {
            printlog ("ERROR - GAREDUCE: GAPREPARE returned an non-zero \
                status. Exiting.", l_logfile, verbose+)
            goto crash
        }
        # The status can be 2 which means a warning - i.e., some files were not
        # processed which would mean sorting out the nextlist ##M This is
        # problem I may want to fix in the future

        fields (tmpfile, fields=2, lines="1-", quit_if_miss=no, \
            print_file_names=no, >> tmplist)
        sections (l_gaprep_pref//"//@"//tmplist, \
            option="fullname", >> nextlist)

        delete (tmplist, verify-, >& "dev$null")
        delete (tmpgaprep, verify=no, >& "dev$null")

    } # End of prepare loop

    # Check if there are any images in the nextlist list
    lcount = 0

    if (access(nextlist)) {
        count (nextlist) | scan (lcount)
    }

    if (lcount == 0) {
        printlog ("ERROR - GAREDUCE: No input images can be used. \
            Please try again.", l_logfile, verbose+)
        goto crash
    }

    # Loop over nextlist to obtain the METACONFIG keyword. This contains
    # information that will group files such that a) they are in the same
    # configuration and b) that they have been reduced in the same way.
    printlog ("\nGAREDUCE: Accessing METACONF keyword values", \
        l_logfile, l_verbose)

    scanfile = nextlist
    while (fscan(scanfile, inimg) != EOF) {

        inphu = inimg//"[0]"

        # Double check that you can access the file...
        if (!imaccess(inimg)) {
            printlog ("WARNING - GAREDUCE: Image "//inimg//\
                " does not exist. Skipping this image.", l_logfile, verbose)
            print (inimg, >> tmpnotproc)

        } else {
            # Read the METACONFIG keyword
            keypar (inphu, "METACONF", silent+)
            if (!keypar.found) {
                printlog ("WARNING - GAREDUCE: METACONF keyword not found \
                    in "//inimg//" Skipping this image.", l_logfile, verbose+)
                print (inimg, >> tmpnotproc)

            } else {
                metaconf = keypar.value

                # Check it's valid
                if (metaconf == "UNKOWN" || metaconf == "" || \
                    (stridx(" ",metaconf) > 0)) {
                    printlog ("WARNING - GAREDUCE: METACONF keyword is "//\
                        "\"UNKNOWN\" or bad in"//inimg//".\n"//\
                        "                  Skipping this image.", \
                        l_logfile, verbose+)
                    print (inimg, >> tmpnotproc)
                } else {

                    print (inimg//" "//metaconf, >> filelist)
                }
            }
        }
    } # End of loop to read METACONF keyword

    # Check if there are any images in the newly created filelist
    lcount = 0

    if (access(filelist)) {
        count (filelist) | scan (lcount)
    }

    if (lcount == 0) {
        printlog ("ERROR - GAREDUCE: No images with a valid METACONF keyword \
            value. Please try again.", l_logfile, verbose+)
        goto crash
    }

    # Determine the unique METACONF values for the file that made it
    printlog ("GAREDUCE: Determining unique METACONFIG values", \
        l_logfile, l_verbose)

    # Find the unique metaconf values
    fields (filelist, fields=2, lines="1-", quit_if_miss=no, \
        print_file_names=no) | sort ("STDIN", column=0, ignore_white=no, \
        numeric_sort=no, reverse_sort=yes) | unique ("STDIN", > tmpconflist)

    # Check if there are any images in the newly created conflist
    lcount = 0

    if (access(tmpconflist)) {
        count (tmpconflist) | scan (lcount)
    }

    if (lcount == 0) {
        printlog ("ERROR - GAREDUCE: No images with a valid METACONF keyword \
            value. Please try again.", l_logfile, verbose+)
        goto crash
    }

    ####
    # Find calibrations... (does not include Sky frames)

    # Call gacalfind if required
    if (l_darkimg == "find" || strstr("find",l_flatimg) == 1) {

        printlog ("\nGAREDUCE: Calling GACALFIND to find processed \
            calibrations...", l_logfile, l_verbose)

        gacalfind (inimages="*.fits", calpath=l_calpath, caltable=l_caltable, \
            fl_calrun=l_fl_calrun, fl_find=no, caltype="", sciimg="", \
            maxtime=l_cal_maxtime, ignore_nlc=ignore_cal_nlc, \
            sci_ext=l_sci_ext, var_ext=l_var_ext, dq_ext=l_dq_ext, \
            logfile=l_logfile, verbose=debug)

        if (gacalfind.status != 0) {
            printlog ("ERROR - GAREDUCE: GACALFIND returned a non-zero \
                status. Exiting.", l_logfile, verbose+)
            goto crash
        } else if (!access(l_caltable)) {
            printlog ("ERROR - GAREDUCE: Cannot access \""//l_caltable//\
                "\". Exiting.", l_logfile, verbose+)
            goto crash
        } else {
            printlog ("GAREDUCE: Returned from GACALFIND\n", \
                l_logfile, l_verbose)
        }
    }

    ####
    # Process the unique configurations...

    total_num_created = 0

    # This becomes a maze of lists! Reader beware.

    # Create separate lists for each unique METACONF keyword
    printlog ("GAREDUCE: Creating METACONFIG lists and processing images...", \
        l_logfile, l_verbose)

    # Set orig variables for calibration settings -  allows tmpfiles to be
    # created to speed things up
    l_darkimg_orig = l_darkimg
    l_flatimg_orig = l_flatimg
    l_skyimg_orig = l_skyimg

    tmplist = mktemp("tmplist")
    nconf = 0

    scanfile = tmpconflist
    while (fscan(scanfile, metaconf) != EOF) {

        nconf += 1

        printlog ("\n----\nGAREDUCE: Processing configuration: "//\
            metaconf, l_logfile, l_verbose)
        # Extract the files that have the same METACONF keyword value
        match (metaconf, filelist, stop=no, print_file_name=no, \
            metacharacter=no) | fields ("STDIN", fields=1, lines="1-", \
            quit_if_miss=no, print_file_names=no, >> conflist)

        # Set a flag that gets used later on to warn uer of mismatched
        # COADDS between sky and image; due to the fact COADDS are averaged.
        unmatched_coadds = no

        # Determine which img type it is so as to parse the metaconf later
        for (jj = 1; jj <= num_imgtypes; jj += 1) {
            if (strstr(poss_imgtype[jj],metaconf) > 0) {
                imgtype = poss_imgtype[jj]
                break
            } else if (jj == num_imgtypes) {
                printlog ("ERROR - GAREDUCE: METACONF "//metaconf//\
                    " does not\n                  contain a valid image type. \
                    Exiting.", \
                    l_logfile, verbose+)
                goto crash
            }
        }

        # List input files for this configuration
        printlog ("GAREDUCE: Input files for this configuration are:\n", \
            l_logfile, l_verbose)

        type (conflist) | tee (l_logfile, out_type="text", append+, \
            >& l_redirect)

        # Form part of the expression to test the sky image later
        testsky1 = substr(metaconf,1,strstr(imgtype,metaconf)-2)
        if (ignore_coadds) {
            nscicoadds = substr(testsky1,1,stridx("+",testsky1)-1)
            nscicoadds = substr(nscicoadds,strldx("_",nscicoadds)+1,\
                strlen(nscicoadds))

            testsky1 = substr(testsky1,1,\
                strstr("_"//nscicoadds//"+",testsky1)-1)//\
                substr(testsky1, \
                   strstr("_"//nscicoadds//"+",testsky1) + \
                   int (1+strlen(nscicoadds)),\
                   strlen(testsky1))

            if (debug) {
                printlog ("____"//testsky1//" {} "//nscicoadds, \
                    l_logfile, verbose+)
            }
        }

        # Strip from metaconf from OBJ
        im_type_len = strlen(imgtype)
#        if (imgtype == "DOME" || imgtype == "GCAL") {
#            im_type_len += 5
#        }
        proc_test = substr(metaconf,\
            strstr(imgtype,metaconf)+im_type_len,\
            strlen(metaconf))

        if (debug) printlog ("___proc_test before: "//proc_test, \
            l_logfile, verbose+)

        proc_test = substr(proc_test,strstr("+",proc_test)+1,\
            strlen(proc_test))

        if (debug) printlog ("___proc_test after: "//proc_test, \
            l_logfile, verbose+)

        ##M Can we ignore NLC ?
        tmpproc_test = proc_test
        if (ignore_cal_nlc) {
            tmpproc_test = substr(tmpproc_test,1,strstr("NLC",tmpproc_test)-2)
        }

        # Reset the num_created counter
        num_created = 0
        num_counter = 0
        first = yes

        # Reset calibrations for the next configuration
        l_darkimg = l_darkimg_orig
        l_flatimg = l_flatimg_orig
        l_skyimg = l_skyimg_orig
        fl_usersky = no

        # Loop over the current configuration list
        scanfile2 = conflist
        while (fscan(scanfile2, inimg) != EOF) {

            # The file in conflist are to be treated the same way...
            num_counter += 1

            if (num_counter > 1) {
                first = no
            }

            fparse (inimg)
            img_root = fparse.root
            inname = img_root//fparse.extension
            inphu = inimg//"[0]"

            outname = l_outpref//img_root//".fits"

            printlog ("\n--\nGAREDUCE: Processing image - "//inname//"...", \
                l_logfile, l_verbose)

            # Reset the processed part of the metaconfiguration
            meta_proc = ""

            # Reset the call to gasky flag to be no - gets set to yes later if
            # needed
            doskycomb = no

            if (debug) {
                date ("+%H:%M:%S.%N%n") | scan (l_time_struct)
                printf ("____TIMESTAMP - %s started: %s\n", \
                    inname, l_time_struct, >> l_logfile)
            }

            keypar (inphu, "NSCIEXT", silent+)
            if (keypar.found) {
                nsci = int(keypar.value)
            } else {
                printlog ("ERROR - GAREDUCE: CAnnot access NSCIEXT keyword \
                    in "//inphu//". Exiting", l_logfile, verbose+)
            }

            keypar (inphu, "GAREDUCE", silent+)
            if (keypar.found) {
                processed = yes
            } else {
                processed = no
            }

            # Determine the UT in seconds of the file
            l_timename = "UT"
            # The observation time of the flat
            keypar (inphu, l_timename, silent+)
            if (!keypar.found) {
               printlog ("ERROR - GAREDUCE: "//timename//" keyword not \
                   found in "//off_flat//". Exiting", \
                   l_logfile, verbose+)
               goto crash
            }
            uttime = keypar.value

            l_datename = "DATE-OBS"
            keypar (inphu, l_datename, silent+)
            if (!keypar.found) {
               printlog ("ERROR - GAREDUCE: "//datename//" keyword not \
                   found in "//off_flat//". Exiting", \
                   l_logfile, verbose+)
               goto crash
            }
            val_date = keypar.value

            cnvtsec (val_date, uttime) | scan (obs_time)

            # Set the vardq flag for each image!
            if (l_fl_vardq) {
                dovardq = yes
            } else {
                dovardq = no
            }

            # Determine the reduction state of the image
            isdarksub = no
            isflat = no
            isskysub = no
            isgainmult = no
            if (processed) {
                keypar (inphu, "DARKIMG", silent+)
                if (keypar.found) {
                    isdarksub = yes
                }

                keypar (inphu, "FLATIMG", silent+)
                if (keypar.found) {
                    isflat = yes
                }

                keypar (inphu, "SKYIMG", silent+)
                if (keypar.found) {
                    isskysub = yes
                }

                keypar (inphu, "BUNIT", silent+)
                if (keypar.found) {
                    if (keypar.value == "electrons")
                    isgainmult = yes
                } else {
                    printlog ("WARNING - GAREDUCE: BUNIT keword not found \
                        in "//inphu//". Assuming units of image are ADU", \
                        l_logfile, verbose+)
                }

                if (isdarksub && isflat && isskysub && isgainmult) {
                    printlog ("WARNING - GARREDUCE: There are no steps \
                        to perform on image \""//inname//"\" it's already \
                        reduced.", l_logfile, verbose+)
                    goto NEXT_IMAGE
                }
            }

            # Check for multiplying by gain twice
            if (l_fl_mult) {
                cangainmult = yes

                if (isgainmult) {
                    cangainmult = no
                }
            } else {
                cangainmult = no
            }

            if (cangainmult) {
                printlog ("GAREDUCE: Will multiply image by GAIN to convert \
                    to electrons.", l_logfile, verbose+)
            }

            # For now run the same things for all images - Can cut this to just
            # the first image later if needed.
            if (l_fl_dark) {

                chkcal = no

                # Check if the current image needs dark subtracting -
                if (isdarksub) {
                    printlog ("WARNING - GAREDUCE: Image \""//inname//\
                        "\" has already been dark subtracted."//\
                        "\n                    Switching dark \
                        subtraction off for this image.", \
                        l_logfile, verbose+)
                    candark = no
                } else if (first) {

                    # Will get set to no further on in the l_fl_flat if
                    # statement if there is a problem. Gets reset for each
                    # configuration.
                    chkcal = yes
                    candark = yes
                    chkcaltype = "DARK"

                    # Let gacalfind find it
                    if (l_darkimg == "find") {

                        chkincal = "*.fits"
                        chkpath = l_calpath
                        chktab = l_caltable

                    } else {

                        # Check user input (again?)
                        chkincal = l_darkimg
                        chkpath = l_dcalpath
                        chktab = tmpdcal
                    }
                }

                if (chkcal) {
                    # Call gacalfind to find the correct / check the flat
                    gacalfind (chkincal, calpath=chkpath, \
                        caltable=chktab, \
                        fl_calrun=no, fl_find=yes, caltype=chkcaltype, \
                        sciimg=inimg, maxtime=INDEF, \
                        ignore_nlc=ignore_cal_nlc, sci_ext=l_sci_ext, \
                        var_ext=l_var_ext, dq_ext=l_dq_ext, \
                        logfile=l_logfile, verbose=debug)

                    if (gacalfind.status != 0) {
                        printlog ("ERROR - GAREDUCE: GACALFIND returned a \
                            non-zero status. Exiting.", l_logfile, verbose+)
                        goto crash
                    } else if (gacalfind.out_calimg == "") {
                        printlog ("WARNING - GAREDUCE: GACALFIND could not \
                            find an appropriate DARK IMAGE for "//inname//\
                            ".\n                    Setting fl_dark = no for \
                            this image.", l_logfile, verbose+)
                        candark = no
                        darkfile = ""
                    } else {
                        darkfile = gacalfind.out_calimg
                    }
                }

                # Double check the presence of the file
                if (candark) {
                    if (!imaccess(darkfile)) {
                        printlog ("ERROR - GAREDUCE: Cannot access dark "//\
                            "image \""//darkfile//"\". Exiting.", \
                            l_logfile, verbose+)
                        goto crash
                    } else {
                        printlog ("GAREDUCE: Using \""//darkfile//"\" \
                            as DARK image.", l_logfile, l_verbose)
                    }
                }

            } else {
                candark = no
                darkfile = ""
            } # End of l_fl_dark check

            if (l_fl_flat) {

                chkcal = no

                # Check if the current image needs flat dividing subtracting -
                if (isflat) {
                    printlog ("WARNING - GAREDUCE: Image \""//inname//\
                        "\" has already been flat divided."//\
                        "\n                    Switching flat \
                        division off for this image.", \
                        l_logfile, verbose+)
                    canflat = no
                } else if (first) {

                    # Will get set to no further on in the l_fl_flat if
                    # statement if there is a problem. Gets reset for each
                    # configuration.
                    canflat = yes
                    chkcal = yes
                    chkcaltype = "FLAT"

                    # Let gacalfind find it
                    if (strstr("find",l_flatimg) == 1) {

                        chkincal = "*.fits"
                        chkpath = l_calpath
                        chktab = l_caltable
                        chkcaltype = chkcaltype//substr(l_flatimg,\
                                         strlen("find")+1,strlen(l_flatimg))

                    } else {
                        # Check user input (again?)
                        chkincal = l_flatimg
                        chkpath = l_fcalpath
                        chktab = tmpfcal
                    }
                }

                if (chkcal) {

                    # Call gacalfind to find the correct / check the flat
                    gacalfind (chkincal, calpath=chkpath, \
                        caltable=chktab, \
                        fl_calrun=no, fl_find=yes, caltype=chkcaltype, \
                        sciimg=inimg, maxtime=INDEF, \
                        ignore_nlc=ignore_cal_nlc, sci_ext=l_sci_ext, \
                        var_ext=l_var_ext, dq_ext=l_dq_ext, \
                        logfile=l_logfile, verbose=debug)

                    if (gacalfind.status != 0) {
                        printlog ("ERROR - GAREDUCE: GACALFIND returned a \
                            non-zero status. Exiting.", l_logfile, verbose+)
                        goto crash
                    } else if (gacalfind.out_calimg == "") {
                        printlog ("WARNING - GAREDUCE: GACALFIND could not \
                            find an appropriate FLAT IMAGE for "//inname//\
                            ".\n                    Setting fl_flat = no for \
                            this image.", l_logfile, verbose+)
                        canflat = no
                        flatfile = ""
                    } else {
                        flatfile = gacalfind.out_calimg
                    }
                }

                # Double check the presence of the file
                if (canflat) {
                    if (!imaccess(flatfile)) {
                        printlog ("ERROR - GAREDUCE: Cannot access flat "//\
                            " image \""//flatfile//"\". Exiting.", \
                            l_logfile, verbose+)
                        goto crash
                    } else {
                        printlog ("GAREDUCE: Using \""//flatfile//"\" \
                            as FLAT field.", l_logfile, l_verbose)
                    }
                }

            } else {
                canflat = no
                flatfile = ""
            } # End of l_fl_flat check

            if (l_fl_sky) {

                # Check if the current image needs dark subtracting -
                if (isskysub) {
                    printlog ("WARNING - GAREDUCE: Image \""//inname//\
                        "\" has already been sky corrected."//\
                        "\n                    Switching sky \
                        correction off for this image.", \
                        l_logfile, verbose+)
                    cansky = no
                    # This is dangerous - it depends on the way the sky was
                    # built for the moment skip - may want to change to a
                    # warning and continue

                    if (skip_if_skysub) {
                        printlog ("                    Processing depends on \
                            sky used, so skipping.", l_logfile, verbose+)
                        goto NEXT_IMAGE
                    }
                } else if (!determine_sky && first) {
                    # now let's see what is to be done in terms of sky

                    # use a pre-generated sky frame (skyimg = <filename>)
                    # need to check if it exists in the calpath, if it matches
                    # the image configuration
                    # for filter only, if it is of the same size as the image,
                    # if it has NOT been dark subtracted and flat fielded, and
                    # NOT  been converted to e-

                    skyname = l_skyimg
                    skyfile = l_skyimg
                    fparse (skyname, verbose-)
                    if (fparse.extension != ".fits") {
                        skyfile = skyfile//".fits"
                    }

                    # exists in the pwd
                    if (!access(skyfile)) {
                        skyfile = l_calpath//skyfile

                        # exists in the calpath
                        if (!access(skyfile)) {
                            printlog ("WARNING - GAREDUCE: Cannot access sky \
                                file "//l_skyimg//" in PWD or in "//\
                                l_calpath//". Switching sky subtraction off", \
                                l_logfile, verbose+)
                            cansky = no
                        }

                    } else {
                        cansky = yes
                    }

                    if (cansky) {
                        # matches filter configuration
                        # ITIME_LNRS_NCOADDS+FILTER1_FILTER2; ROIID+TRIM+NLC

                        # ERROR if metaconf not found
                        keypar (skyfile//"[0]", "METACONF", silent+)
                        if (!keypar.found) {
                            printlog ("ERROR - GAREDUCE: METACONF keyword not \
                                found in "//skyname//". Exiting", \
                                l_logfile, verbose+)
                            goto crash
                        }
                        full_skymeta = keypar.value
                        skymeta1 = substr(full_skymeta,1,\
                            strstr("SKY",full_skymeta)-2)

                        skymeta2 = substr(full_skymeta,\
                            strstr("SKY",full_skymeta)+strlen("SKY")+1,\
                            strlen(full_skymeta))

                        # Allow a tolerance in the exposure time
                        # Note the number of coadds is important to the
                        # determination of the noise. However, the counts in an
                        # image with > 1 coadd is an average.

                        # Ignore coadds is set at the top of the script
                        if (ignore_coadds) {
                            nskycoadds = substr(skymeta1,1,\
                                stridx("+",skymeta1)-1)

                            nskycoadds = substr(nskycoadds,\
                                strldx("_",nskycoadds)+1,\
                                strlen(nskycoadds))

                            skymeta1 = substr(skymeta1,1,\
                                strstr("_"//nskycoadds//"+",skymeta1)-1)//\
                                substr(skymeta1, \
                                   strstr("_"//nskycoadds//"+",skymeta1) + \
                                   int (1+strlen(nskycoadds)),\
                                   strlen(skymeta1))

                            # Print a warning to user if not dark subtracting
                            # or dark subtracted if coadds don't match
                            if (nskycoadds != nscicoadds) {
                                unmatched_coadds = yes
                            }

                            if (debug) {
                                printlog ("____"//skymeta1//" {} "//\
                                    nskycoadds, \
                                    l_logfile, verbose+)
                            }
                        }

                        # This falg is set at the top of the script
                        if (allow_tolerance) {
                            sky_exp = real(substr(skymeta1,1,\
                                stridx("_",skymeta1)-1))
                            s_obj_exp = substr(testsky1,1,\
                                stridx("_",testsky1)-1)
                            obj_exp = real(s_obj_exp)

                            # exp_tolerance is set at the top of the script
                            if (abs(sky_exp - obj_exp) < exp_tolerance) {
                                skymeta1 = s_obj_exp//substr(skymeta1,\
                                    stridx("_",skymeta1),strlen(skymeta1))
                            }
                        }

                        if (debug) printlog ("____"//skymeta1//" "//\
                            skymeta2//" "//testsky1//" "//tmpproc_test, \
                            l_logfile, verbose+)

                        if (skymeta1 != testsky1) {

                            printlog ("WARNING - GAREDUCE: sky image "//\
                                skyfile//" filter combination "//\
                                "\n                    and or exposure \
                                configuration do not match. "//\
                                "\n                    Switching sky \
                                subtraction off.", l_logfile, verbose+)
                            cansky = no
                        } else if (strstr(tmpproc_test,skymeta2) == 0) {

                            printlog ("WARNING - GAREDUCE: sky image "//\
                                skyfile//" has not been"//\
                                "\n                    processed in the same \
                                way as "//inname//" by gaprepare."//\
                                "\n                    Switching \
                                sky subtraction off", l_logfile, verbose+)
                            cansky = no
                        }
                    }

                   # Set a flag to say using sky input not making sky
                   fl_usersky = yes

                } else if (determine_sky) {

                    cansky = yes

                    # Set up the parameters to select the sky
                    sdistance = l_minoffs
                    stime = l_maxtime
                    snearer = no
                    syounger = yes

                    if (l_skyimg == "time") {
                        pr_str = "time"
                        sdistance = INDEF
                    } else if (l_skyimg == "distance") {
                        pr_str = "distance"
                        stime = INDEF
                    } else {
                        pr_str = "both time and distance"
                    }

                    # skyimg = time
                    # create a sky frame from the sequence using time as
                    # parameter (dithered  sequences for point or point-like
                    # sources) and re-doing the sky frame for each
                    # science frame (running average)

                    printlog ("GAREDUCE: Deriving sky frame from data \
                        based on "//pr_str, l_logfile, l_verbose)

                    skyfile = ""

                    # use gemoffsetlist to build the list of sky images

                    # now, I'd like a way to avoid re-creating the same sky for
                    # every image - in other words
                    # only create a sky if none of the already existing work  -
                    # have to compare the combine lists.
                    # Initially, create the list based on the name of the image
                    # being processed, to be sure it does not exist already

                    fparse (inimg, verbose-)
                    skylist = fparse.root//".lst"
                    delete(skylist, verify-, >& "dev$null")

                    # the parameter maxtime defines the largest time difference
                    # within which a candidate is acceptable. There is no
                    # offset restriction - I want to include all pointings, but
                    # there should be a way of avoiding too many of those at
                    # the same position as the reference for the moment, coding
                    # this as a warning with the number of images within
                    # 1arcsec of reference

                    gemoffsetlist ("@"//conflist, inimg, distance=sdistance, \
                        age=stime, targetlist=skylist, \
                        offsetlist="dev$null", \
                        fl_nearer=snearer, direction=3, fl_younger=syounger, \
                        fl_noref=yes, wcs_source="none", key_xoff="XOFFSET", \
                        key_yoff="YOFFSET", key_date=l_datename, \
                        key_time=l_timename, logfile=l_logfile,\
                        verbose=debug, force=no)

                    if (gemoffsetlist.status != 0) {
                        printlog ("ERROR - GAREDUCE: GEMOFFSETLIST returned \
                            a non-zero status. Exiting.", l_logfile, verbose+)
                        goto crash
                    }

                    nsky1 = gemoffsetlist.count

                    # cannot do it if only one image (should use filename
                    # instead)

                    if (nsky1 < 2) {
                        printlog ("WARNING - GAREDUCE: cannot combine less \
                            than 2 images.\n                    will not \
                            sky subtract image "//inimg, l_logfile, verbose+)
                        cansky = no
                        doskycomb = no

                        delete (skylist, verify-, >& "dev$null")

                    } else if (nsky1 < l_minsky) {

                        # have less than the minimum number requested
                        printlog ("WARNING - GAREDUCE: number of sky \
                            images less than minimum number requested.\n"//\
                            "                    will not sky subtract \
                            image "//inimg, l_logfile, verbose+)

                        delete (skylist, verify-, >& "dev$null")

                        cansky = no
                        doskycomb = no

                    } else {

                        # this is the skylist check, silencing the output from
                        # gemoffsetlist
                        if (l_skyimg == "time") {
                            gemoffsetlist ("@"//skylist, inimg, distance=1.0, \
                                age=INDEF, targetlist="dev$null", \
                                offsetlist="dev$null", fl_nearer=yes, \
                                direction=3, fl_younger=yes, \
                                fl_noref=yes, wcs_source="none", \
                                key_xoff="XOFFSET", key_yoff="YOFFSET", \
                                key_date=l_datename, key_time=l_timename, \
                                logfile=l_logfile, verbose-, force=no)

                            if (gemoffsetlist.status != 0) {
                                printlog ("ERROR - GAREDUCE: GEMOFFSETLIST \
                                    returned a non-zero status. Exiting.", \
                                l_logfile, verbose+)
                                goto crash
                            }

                            nsky2 = gemoffsetlist.count
                            pr_str = "within "//l_maxtime//" seconds \
                                     of reference image "//inimg//\
                                     ",\n          of which "//nsky2//\
                                     " are located within 1 arcsec of the \
                                     reference"

                        } else if (l_skyimg == "distance") {
                            pr_str = "greater than "//l_minoffs//" arcseconds \
                                     from reference image "//inimg
                        } else {
                            pr_str = "within "//l_maxtime//" seconds \
                            and greater than "//l_minoffs//" arcseconds away \
                            from of reference image "//inimg
                        }
                        printlog ("GAREDUCE: There are "//nsky1//" images \
                            to combine "//pr_str, l_logfile, l_verbose)

                        # since I go one image at a time, I can compare the
                        # current list with all existing ones then delete all
                        # of them at the end of the config loop. just use diff
                        # to a tempfile and count the lines, if not empty there
                        # is a difference The name of the sky frame is based on
                        # the first image in the list. Note that I can have
                        # two lists starting with the same image, but not
                        # having the same content, so I have to be able to
                        # differentiate the sky generated - use a mktemp as
                        # sufix. at the same time, I need to be able to
                        # associate the list with the sky frame generated, so
                        # here I rename the input list using the same
                        # root+sufix as the sky

                        if (debug)
                            print ("____doskycomb="//doskycomb)

                        skyname = img_root//"_"//l_sky_sufx//".fits"
                        if (access(skyname)) {
                            printlog ("ERROR - GAREDUCE: "//skyname//\
                                " already exists", l_logfile, verbose+)
                            goto crash
                        }

                        # Check stored sky files from this run for their inputs
                        # to see if they match the current list
                        if (access(tmpkeepskylst)) {
                            if (debug)
                                type(tmpkeepskylst)

                            scanfile3 = tmpkeepskylst
                            while (fscan(scanfile3, tmpskyfile) != EOF) {
                                if (debug)
                                    print (tmpskyfile)

                                hselect (tmpskyfile//"[0]", "GASKIM*", yes) |\
                                    token (ignore_comments=yes, \
                                    begin_comment="\#", end_comment="eol", \
                                    newlines=yes, >> tmplst)

                                # remove the final new line due to tokens
                                nrows = 0
                                count (tmplst) | scan (nrows)
                                head (tmplst, nlines=(nrows - 1), \
                                    >> tmpskyinnames)

                                diff (skylist, tmpskyinnames, >> tmpdiflst)
                                if (debug) {
                                    print ("+++++0")
                                    type (skylist)
                                    print ("+++++1")
                                    type (tmpskyinnames)
                                    print ("+++++2")
                                    type (tmpdiflst)
                                }

                                nrows = 0
                                count (tmpdiflst) | scan (nrows)

                                delete (tmpskyinnames//","//tmpdiflst//","//\
                                    tmplst, verify-, >& "dev$null")

                                if (nrows == 0) {
                                    fparse (tmpskyfile, verbose-)
                                    skyfile = fparse.root//".fits"
                                    printlog ("WARNING - GAREDUCE: Sky frame \
                                        for this skylist already exists: "//\
                                        skyfile, l_logfile, verbose+)
                                    break
                                } else {
                                    #Perform skycomb
                                    doskycomb = yes
                                    skyfile = skyname
                                }
                            }
                        } else {
                            # should be the first time around and therefore
                            # have enough input files and no outputsky files to
                            # compare to
                            doskycomb = yes
                            skyfile = skyname
                        } # End of checking for previous versions
                    } # End of nsky checks

                    if (debug) print ("____doskycomb="//doskycomb)
                } # End of determine_sky loop

                if (doskycomb) {

                    # the combine list is skylist, the sky name is skyname
                    # this is for later deletion of the skylists and skyframes
                    # if fl_keepsky = no
                    print (skylist, >> tmpkeeplst)
                    print (skyname, >> tmpkeeplst)
                    print (skyname, >> tmpkeepskylst)

                    printlog ("GAREDUCE: Calling GASKY", l_logfile, l_verbose)

                    if (l_qreduce && !isindef(flatfile)) {
                        inflat_gasky = flatfile
                    } else {
                        inflat_gasky = ""
                    }

                    gasky ("@"//skylist, outimages=skyname, \
                        outsufx=l_sky_sufx, \
                        key_exptime=l_expname, combine=l_combine, \
                        statsec=l_statsec, reject=l_reject, nlow=l_nlow, \
                        nhigh=l_nhigh, masktype=l_masktype, \
                        maskvalue=l_maskvalue, \
                        sci_ext=l_sci_ext, var_ext=l_var_ext, \
                        dq_ext=l_dq_ext, key_ron=l_key_ron, \
                        key_gain=l_key_gain, masksuffix=l_masksuffix, \
                        fl_keepmasks=l_fl_keepmasks, threshold=l_lthreshold, \
                        ngrow=l_ngrow, agrow=l_agrow, minpix=l_minpix, \
                        qreduce=l_qreduce, flatimg=inflat_gasky, \
                        fl_mask=l_fl_mask, fl_vardq=l_fl_vardq, \
                        fl_dqprop=l_fl_dqprop, minflat_val=l_minflat_val, \
                        logfile=l_logfile, verbose=l_verbose)

                    delete (skylist, verify-, >& "dev$null")

                    if (gasky.status != 0) {
                        printlog ("ERROR - GAREDUCE: GASKY returned a \
                            non-zero status. Exiting.", l_logfile, verbose+)
                        goto crash
                    }
                }

                if (cansky) {
                    printlog ("GAREDUCE: Using \""//\
                        skyfile//"\" as sky frame", l_logfile, l_verbose)
                    # Check for frozen ODGWs
                    keypar (skyfile//"[0]", l_key_froz_odgw, silent+)
                    if (keypar.found) {
                        printlog ("WARNING - GAREDUCE: Skyfile \""//\
                            skyfile//"\" may contain areas with no data "//\
                            "due to frozen ODGWs.", l_logfile, verbose+)
                    }
                }
            } else {
                cansky = no
                skyfile = ""
            } # End of l_fl_sky check

            # Check if there is anything to do
            if (!candark && !cansky && !canflat && !cangainmult) {
                printlog ("WARNING - GAREDUCE: There is nothing to be done "//
                    "to "//inname//".", \
                    l_logfile, verbose+)
                if (update_metaconf) {
                    printlog ("                    Going to next \
                        configuration", l_logfile, verbose+)
                    goto NEXT_CONFIGURATION
                } else {
                    printlog ("                  Skipping to next image", \
                        l_logfile, verbose+)
                    goto NEXT_IMAGE
                }
            }

            if (debug) printlog ("____candark: "//candark//" canflat: "//\
                canflat//" cansky: "//cansky, l_logfile, verbose+)
            if (debug) printlog ("____isdarksub: "//isdarksub//" isflat: "//\
                isflat//" isskysub: "//isskysub, l_logfile, verbose+)
            if (debug) printlog ("____darkfile: \""//darkfile//\
                "\"\n____flatfile: \""//flatfile//"\"\n____skyfile: \""//\
                skyfile//"\"", l_logfile, verbose+)

            # Determine the reduction states of the calibration images
            # Dark image
            isdarkflat = no
            dark_isgainmult = no
            # Will fail if sky subtracted

            if (candark) {
                # Is it flat?
                keypar (darkfile//"[0]", "FLATIMG", silent+)
                if (keypar.found) {
                    isdarkflat = yes
                }

                # Is it sky subtracted if so - complain and error!
                keypar (darkfile//"[0]", "SKYIMG", silent+)
                if (keypar.found) {
                    printlog ("ERROR - GAREDUCE: Dark image "//darkfile//\
                        " has been sky subtracted. Exiting.", \
                        l_logfile, verbose+)
                    goto crash
                }

                # Check units of image
                keypar (darkfile//"[0]", "BUNIT", silent+)
                if (keypar.found && keypar.value == "electrons") {
                    dark_isgainmult = yes

                } else if (!keypar.found) {
                    printlog ("ERROR - GAREDUCE: Cannot access BUNIT keyword \
                        in "//darkfile//"[0]. Exiting", \
                        l_logfile, verbose+)
                    goto crash
                }
            } # End of dark checks

            # Sky image
            isskydark = no
            isskyflat = no
            sky_isgainmult = no

            if (cansky) {

                # Is it dark subtracted
                keypar (skyfile//"[0]", "DARKIMG", silent+)
                if (keypar.found) {
                    isskydark = yes
                } else if (candark) {
                    if (determine_sky) {
                        printlog ("GAREDUCE: Switching off dark subtraction \
                            as sky frame determined from science inputs", \
                            l_logfile, verbose+)
                    } else {
                        printlog ("WARNING - GAREDUCE: Sky subtracting a \
                            non-dark subtracted sky frame and also attempting \
                            \n                    to dark subtract. Switching \
                            off dark subtraction.", l_logfile, verbose+)
                    }
                    candark = no
                } else if (unmatched_coadds) {
                    printlog ("WARNING - GAREDUCE: Sky subtracting a non-dark \
                        corrected sky image \
                        \n                    from a non-dark corrected \
                        science image, \
                        \n                    each with a different \
                        number of COADDs.", \
                        l_logfile, verbose+)
                }

                # Is it flat divided if so.
                keypar (skyfile//"[0]", "FLATIMG", silent+)
                if (keypar.found) {
                    isskyflat = yes
                    skyflatimg = keypar.value
                    fparse (skyflatimg)
                    skyflatimg = fparse.root

                    if (isflat) {
                        # Science is flat too
                        # Check that the flat fields match in name at least
                        keypar (inphu, "FLATIMG", silent+)
                        if (keypar.value != skyflatimg) {

                            printlog ("WARNING - GAREDUCE: Both "//inname//\
                                " and "//skyfile//"\n"//\
                                "                    are flat divided but \
                                they flat fielded using different flat \
                                images,"//\
                                "\n                   "//keypar.value//\
                                " and "//\
                                skyflatimg//", respectively", \
                                l_logfile, verbose+)
                        }
                    } else if (canflat) {
                        # Will be flat fielding - check the filenames match
                        fparse (flatfile)
                        if (skyflatimg != fparse.root) {
                            printlog ("ERROR - GAREDUCE: Sky image "//\
                                skyfile//" was flat fielded by "//\
                                "\n                  "//skyflatimg//\
                                " but the image\n                  to flat "//\
                                "field the science input is "//\
                                "\n                  "//flatfile//\
                                ".\n                  "//"Exiting.", \
                                l_logfile, verbose+)
                            goto crash
                        }
                    } else if (!canflat) {
                        # Can't sky subtract
                        printlog ("WARNING - GAREDUCE: Sky image is \
                            flat fielded but science is not and not flat \
                            fielding in this call.\n                    \
                            Switching sky subtraction off.", \
                            l_logfile, verbose+)
                        cansky = no
                    }
                } else if (isflat) {
                    printlog ("WARNING - GAREDUCE: Input image is flat \
                        fielded but the sky frame is not.\n"//\
                        "                    Switching sky subtaction off", \
                        l_logfile, verbose+)
                    cansky = no
                }

                if (cansky) {
                    # Check units of image
                    keypar (skyfile//"[0]", "BUNIT", silent+)
                    if (keypar.found && keypar.value == "electrons") {
                        sky_isgainmult = yes

                    } else if (!keypar.found) {
                        printlog ("ERROR - GAREDUCE: Cannot access BUNIT \
                            keyword in "//skyfile//"[0]. Exiting", \
                            l_logfile, verbose+)
                        goto crash
                    }
                }
            } # End of sky image checks

            # Flat image
            isflatdark = no
            # Flat's are unitless

            if (canflat) {

                # Is it dark subtracted - print warning if not the same as sci
                # Or not dark subtrcating
                keypar (flatfile//"[0]", "DARKIMG", silent+)
                if (keypar.found) {
                    if (!isdarksub && !isskysub && !candark && !cansky) {
                        printlog ("WARNING - GAREDUCE: Flat image "//\
                            flatfile//" has been dark subtracted but "//\
                            "\n                     "//inname//\
                            " is not and is not going to be dark \
                            or skysubtracted.", l_logfile, verbose+)
                    }
                } # I think all of the other the other cases are caught earlier
                  # on due to the switching on and off of things

                # Is it sky subtracted if so - complain and error!
                keypar (flatfile//"[0]", "SKYIMG", silent+)
                if (keypar.found) {
                    printlog ("ERROR - GAREDUCE: Flat image "//flatfile//\
                        " has been sky subtracted. Exiting.", \
                        l_logfile, verbose+)
                    goto crash
                }
            } # End of flat checks

            if (debug) printlog ("____candark: "//candark//" canflat: "//\
                canflat//" cansky: "//cansky, l_logfile, verbose+)
            if (debug) printlog ("____isdarksub: "//isdarksub//" isflat: "//\
                isflat//" isskysub: "//isskysub, l_logfile, verbose+)
            if (debug) printlog ("____darkfile: \""//darkfile//\
                "\"\n____flatfile: \""//flatfile//"\"\n____skyfile: \""//\
                skyfile//"\"", l_logfile, verbose+)

            # Double check if there is anything to do
            if (!candark && !cansky && !canflat && !cangainmult) {
                printlog ("WARNING - GAREDUCE: There is nothing to be done "//
                    "to "//inname//".", \
                    l_logfile, verbose+)
                if (update_metaconf) {
                    printlog ("                    Going to next \
                        configuration",\
                        l_logfile, verbose+)
                    goto NEXT_CONFIGURATION
                } else {
                    printlog ("                    Skipping to next image", \
                        l_logfile, verbose+)
                    goto NEXT_IMAGE
                }
            }

            # Check the GAIN keyword to use for inputs - check only first
            # extension - MS
            if (!isgainmult) {
                testfile = inimg//"["//l_sci_ext//",1]"
                keypar (testfile, "GAINORIG", silent+)
                if (keypar.found) {
                    l_out_key_gain =  "GAINORIG"
                } else {
                    l_out_key_gain =  "GAIN"
                }
            } else {
                l_out_key_gain = "GAINMULT"
            }

            # Now that the can flags have been set check the dimensions of the
            # input calibration frames and check which GAIN keyword to use.
            # Can check only the first due to smart sorting by metaconf - MS
            if (candark && first) {

                darkname = darkfile

                keypar (darkfile//"[0]", "NSCIEXT", silent+)
                ndarksci = int (keypar.value)
                if (ndarksci != nsci) {
                    printlog ("ERROR - GAREDUCE: Number of "//l_sci_ext//\
                        " extensiosn do not match between "//darkfile//\
                        " and "//inimg//". Exiting.", l_logfile, verbose+)
                    goto crash
                }

                # Call gacaltrim to trim the file as needed
                gacaltrim (inimg, darkfile, key_check="CCDSEC", \
                    sci_ext=l_sci_ext, var_ext=l_var_ext, dq_ext=l_dq_ext, \
                    logfile=l_logfile, verbose=debug)

                # If the calibration file doesn't have the same number of
                # extensions it will fail!
                if (gacaltrim.status != 0) {
                    printlog ("ERROR - GAREDUCE: GACALTRIM returned a \
                        non-zero status.", l_logfile, verbose+)
                    goto crash
                } else if (gacaltrim.outimg != darkfile) {
                    # Only update files as required if outname is not equal to
                    # the darkfile name
                    # reset tmpdark so it gets cleaned up
                    tmpdark = gacaltrim.outimg
                    darkfile = tmpdark
                }

                # Check which GAIN keyword to use
                # Assumes all extension have been operated on
                # uniformily - MS
                testfile = darkfile//"["//l_sci_ext//",1]"

                if (dark_isgainmult) {
                    darkgain = "GAINMULT"
                } else {
                    keypar (testfile, "GAINORIG", silent+)
                    if (keypar.found) {
                        darkgain = "GAINORIG"
                    } else {
                        keypar (testfile, "GAIN", silent+)
                        if (!keypar.found) {
                            printlog ("ERROR - GAREDUCE: Cannot \
                                access GAIN keyword in "//\
                                testfile, l_logfile, verbose+)
                        }
                    }
                } # End of check of gain keyword to use

            } # End of checking input dark file (dimensions and vardq)

            # Sky frame - if user inputs sky frame only do it for first image
            if (cansky && first && fl_usersky) {
                keypar (skyfile//"[0]", "NSCIEXT", silent+)
                nskysci = int (keypar.value)
                if (nskysci != nsci) {
                    printlog ("ERROR - GAREDUCE: Number of "//l_sci_ext//\
                        " extensiosn do not match between "//skyfile//\
                        " and "//inimg//". Exiting.", l_logfile, verbose+)
                    goto crash
                }

                # Call gacaltrim to trim the file as needed
                gacaltrim (inimg, skyfile, key_check="CCDSEC", \
                    sci_ext=l_sci_ext, var_ext=l_var_ext, dq_ext=l_dq_ext, \
                    logfile=l_logfile, verbose=debug)

                if (gacaltrim.status != 0) {
                    printlog ("ERROR - GAREDUCE: GACALTRIM returned a \
                        non-zero status.", l_logfile, verbose+)
                    goto crash
                } else if (gacaltrim.outimg != skyfile) {
                    # Only update files as required if outname is not equal to
                    # the skyfile name
                    # reset tmpsky so it gets cleaned up
                    tmpsky = gacaltrim.outimg
                    skyfile = tmpsky
                }

                # Check which GAIN keyword to use
                # Assumes all extension have been operated on
                # uniformily - MS
                testfile = skyfile//"["//l_sci_ext//",1]"

                if (sky_isgainmult) {
                    skygain = "GAINMULT"
                } else {
                    keypar (testfile, "GAINORIG", silent+)
                    if (keypar.found) {
                        skygain = "GAINORIG"
                    } else {
                        keypar (testfile, "GAIN", silent+)
                        if (!keypar.found) {
                            printlog ("ERROR - GAREDUCE: Cannot \
                                access GAIN keyword in "//\
                                testfile, l_logfile, verbose+)
                        }
                    }
                } # End of check of gain keyword to use
                # End of checking input flat file (dimensions and vardq)
            } else if (cansky && !fl_usersky) {
                # Using skies made with inputs (fl_vardq passed to gareduce and
                # gasky) - MS
                skygain = l_out_key_gain
            }

            # Flat image
            if (canflat && first) {
                flatname = flatfile

                keypar (flatfile//"[0]", "NSCIEXT", silent+)
                nflatsci = int (keypar.value)
                if (nflatsci != nsci) {
                    printlog ("ERROR - GAREDUCE: Number of "//l_sci_ext//\
                        " extensiosn do not match between "//flatfile//\
                        " and "//inimg//". Exiting.", l_logfile, verbose+)
                    goto crash
                }

                # Call gacaltrim to trim the file as needed
                gacaltrim (inimg, flatfile, key_check="CCDSEC", \
                    sci_ext=l_sci_ext, var_ext=l_var_ext, dq_ext=l_dq_ext, \
                    logfile=l_logfile, verbose=debug)

                # GACALTRIM will fail if the files have different number of
                # extensions
                if (gacaltrim.status != 0) {
                    printlog ("ERROR - GAREDUCE: GACALTRIM returned a \
                        non-zero status.", l_logfile, verbose+)
                    goto crash
                } else if (gacaltrim.outimg != flatfile) {
                    # Only update files as required if outname is not equal to
                    # the flatfile name
                    # reset tmpflat so it gets cleaned up
                    tmpflat = gacaltrim.outimg
                    flatfile = tmpflat
                }
            } # End of checking input flat file (dimensions and vardq)

            if (debug) printlog ("____candark: "//candark//" canflat: "//\
                canflat//" cansky: "//cansky, l_logfile, verbose+)
            if (debug) printlog ("____isdarksub: "//isdarksub//" isflat: "//\
                isflat//" isskysub: "//isskysub, l_logfile, verbose+)
            if (debug) printlog ("____darkfile: \""//darkfile//\
                "\"\n____flatfile: \""//flatfile//"\"\n____skyfile: \""//\
                skyfile//"\"", l_logfile, verbose+)

            ## Form the expressions
            # Will set expression according to units of input file

            # Expression uses a = img, b = dark, c = sky, d = flat

            # Reset the operands - this is being ultra sure - MS
            a_operand = inimg

            # Dark image
            if (candark) {
                b_operand = darkfile
            } else {
                b_operand = INDEF
            }

            # Sky image
            if (cansky) {
                c_operand = skyfile
            } else {
                c_operand = INDEF
            }

            # Flat frame
            if (canflat) {
                d_operand = flatfile
            } else {
                d_operand = INDEF
            }

            # Only form the expression if first image or not updating the
            # reduction status of the file in the metaconfiguration

#            if (!first && update_metaconf) {
#                goto SKIP_EXPR
#            }

            sciexpr = "a"
            varexpr = "a[VAR]"
            dqexpr = "a[DQ]"

            if (candark) {

                if (isgainmult == dark_isgainmult) {
                    gexpr = ""
                    darkvarexpr = " + b[VAR]"
                } else if (isgainmult) {
                    gexpr = "b."//darkgain//" * "
                    darkvarexpr = " + (b[VAR]."//darkgain//" * b[VAR])"
                } else {
                    gexpr = "(1 / (b."//darkgain//")) * "
                    darkvarexpr = " + ((1 / (b[VAR]."//darkgain//")) * b[VAR])"
                }

                sciexpr = sciexpr//" - ("//gexpr//"b)"
                varexpr = varexpr//darkvarexpr
                dqexpr = dqexpr//" | b[DQ]"
            }

            if (cansky) {

                if (isgainmult == sky_isgainmult) {
                    gexpr = ""
                    skyvarexpr = " + c[VAR]"
                } else if (isgainmult) {
                    gexpr = "c."//skygain//" * "
                    skyvarexpr = " + (c[VAR]."//skygain//" * c[VAR])"
                } else {
                    gexpr = "(1 / (c."//skygain//")) * "
                    skyvarexpr = " + ((1/(c[VAR]."//skygain//")) * c[VAR])"
                }

                skyexpr = " - ("//gexpr//"c)"
                dqexpr = dqexpr//" | c[DQ]"

            } else {
                skyexpr = ""
                skyvarexpr = ""
            }

            if (canflat) {
                d_operand = flatfile

                if (fixed_gemexpr) {
                    d_expr = "chkval(d,"//scirep_val//")"

                } else if (isindef(l_minflat_val)) {
                    d_expr = "(d)"
                } else {
                    d_expr = "((d < "//l_minflat_val//") ? "//l_minflat_val//\
                        " : d)"
                }

                if (cansky && isskyflat) {
                    # Form the var expression first in this case
                    varexpr = "(("//varexpr//" + ((("//sciexpr//") / "//\
                        d_expr//")**2) * chkvval(d,"//scirep_val//\
                        ",d[VAR],"//varrep_val//")) / "//d_expr//"**2)"//\
                        skyvarexpr
                    sciexpr = "((("//sciexpr//") / "//d_expr//")"//skyexpr//")"
                } else {
                    sciexpr = "(("//sciexpr//skyexpr//") / "//d_expr//")"
                    varexpr = "("//varexpr//skyvarexpr//" + ("//sciexpr//\
                        "**2) * chkvval(d,"//scirep_val//\
                        ",d[VAR],"//varrep_val//")) / "//d_expr//"**2"
                }

                dqexpr = dqexpr//" | chkdval(d,"//scirep_val//\
                    ",d[DQ],"//dqrep_val//")"
            } else if (cansky) {
                # Var and DQ expressions already handled
                sciexpr = sciexpr//skyexpr
                varexpr = varexpr//skyvarexpr
            }

            if (cangainmult) {
                sciexpr = "(a."//l_out_key_gain//" * ("//sciexpr//"))"
                varexpr = "(a[VAR]."//l_out_key_gain//" * ("//varexpr//"))"
            }

#SKIP_EXPR:
            # reset the output name if adding sky back in
            if (l_fl_autosky && cansky) {
                gemout = tmpout
            } else {
                gemout = outname
            }

            if (debug) {
                printlog ("____SCIEXPR 2: "//sciexpr, \
                    l_logfile, verbose+)
                printlog ("____VAREXPR 2: "//varexpr, \
                    l_logfile, verbose+)
                printlog ("____DQEXPR: "//dqexpr, \
                    l_logfile, verbose+)

                if (!isindef(a_operand)) {
                    printlog ("____a_operand: "//a_operand, \
                        l_logfile, verbose+)
                    if (caldebug) {
                        imstat (a_operand)
                    }
                }
                if (!isindef(b_operand)) {
                    printlog ("____b_operand: "//b_operand, \
                        l_logfile, verbose+)
                    if (caldebug) {
                        imstat (b_operand)
                    }
                }
                if (!isindef(c_operand)) {
                    printlog ("____c_operand: "//c_operand, \
                        l_logfile, verbose+)
                    if (caldebug) {
                        imstat (c_operand)
                    }
                }
                if (!isindef(d_operand)) {
                    printlog ("____d_operand: "//d_operand, \
                        l_logfile, verbose+)
                    if (caldebug) {
                        imstat (d_operand)
                    }
                }
            }

            # Perform calculation
            gemexpr (sciexpr, gemout, a_operand, b_operand, c_operand, \
                d_operand, var_expr=varexpr, dq_expr=dqexpr, \
                sci_ext=l_sci_ext, var_ext=l_var_ext, dq_ext=l_dq_ext, \
                mdf_ext="MDF", fl_vardq=dovardq, dims="default", \
                intype="default", outtype="ref", refim="a", rangecheck=yes, \
                verbose=debug, exprdb="gsaoi$data/gsaoiEDB.dat", \
                logfile=l_logfile)

            if (gemexpr.status != 0) {
                printlog ("ERROR - GAREDUCE: GEMEXPR returned a non-zero \
                    status. Exiting.", l_logfile, verbose+)
                goto crash
            }

            # Add back in the sky
            if (l_fl_autosky && cansky) {
                # The correct way of adding back in the sky is under discussion
                # The following implimentation will calculate the median sky in
                # electrons using the entire dector statistic and will try and
                # match the reduction state of the science for this pass.
                # The correction value will be adjusted to the correct units
                # when it is applied. An unflattened sky can only be flattened
                # if flat dividing in this pass. The comment for the SKYMIDPT
                # keyword states what state the sky frame was in when it was
                # calulated.

                alt_sky = no
                csky_type = "a"
                sarray_index = 1

                # Set the comment for the SKYMIDPT keyword
                # Value is always stored in electrons
                sky_comment = "Added sky const (electrons; "

                # Always do the following calulation in electrons
                if (!sky_isgainmult) {
                    alt_sky = yes
                    sarray_index = 2
                    csky_type = csky_type//" * a."//skygain
                }

                # We can only flat divide the image if flat dividing input
                # Better to have flat divided your sky first - MS
                ##M  Well we could use the information from above to determine
                ##M  if it can be flat divided.... somehow -- maybe some of teh
                ##M  logic for this should be done above... ??? MS
                ##M Is there any way to flat divide the inputs before sending
                ##M to gasky and storing them?
                # Update sky_comment too
                csky_flat = INDEF

                if (canflat && !isskyflat) {
                    alt_sky = yes
                    csky_type = "("//csky_type//") / b"
                    csky_flat = flatfile
                    sarray_index = 3
                    sky_comment = sky_comment//"flat; "
                } else if (isskyflat) {
                    sky_comment = sky_comment//"flat; "
                } else {
                    sky_comment = sky_comment//"not flat; "
                }

                # Do the calculation on the entire detector - should be the
                # same value that gets added back to all arrays. This could be
                # dependant on how the QE correction is applied. - MS
                auto_sky_stat = "DETECTOR"
                sky_comment = sky_comment//auto_sky_stat//")"

                # Set how the Skymidpt is added back in - adjust for gain
                skyconst_expr = "(a.SKYMIDPT)"

                if (!cangainmult && !isgainmult) {
                    # Output science not in electrons
                    skyconst_expr = "("//skyconst_expr//" / a."//\
                        l_out_key_gain//")"
                }

                if (debug) printlog ("____skyconst_expr: "//skyconst_expr//\
                    ", comment: "//sky_comment//", "//sky_isgainmult//\
                    ", csky_type: "//csky_type, \
                    l_logfile, verbose+)

                # Initialise the skymidpt information
                # There has to be some logic in here about checking if a
                # different sky is being used i.e., it's be scaled

                # If we're doing this calulation in electrons it is the same
                # each time except for if you can flat divide or not; then if
                # first/scaled/!user defined sky
                if (first || skyscaled || !fl_usersky) {

                    # This information is for the various combination so as to
                    # not have to repeat steps when not neccessary.
                    skymidpt[1] = "INDEF"
                    skymidpt[2] = "INDEF"
                    skymidpt[3] = "INDEF"
                }

                # Caluclate skymidpt
                if (skymidpt[sarray_index] == "INDEF") {

                    # Act on sky if required
                    if (alt_sky) {

                        imdelete (tmpgcorr_sky, verify-, >& "dev$null")

                        gemexpr (csky_type, tmpgcorr_sky, skyfile, csky_flat,\
                            var_expr="a[VAR]", dq_expr="a[DQ]", \
                            sci_ext=l_sci_ext, var_ext=l_var_ext, \
                            dq_ext=l_dq_ext, mdf_ext="MDF", fl_vardq=dovardq, \
                            dims="default", intype="default", outtype="ref", \
                            refim="a", rangecheck=yes, verbose=no, \
                            exprdb="none", logfile=l_logfile)

                        if (gemexpr.status != 0) {
                            printlog ("ERROR - GAREDUCE: GEMEXPR returned \
                                a non-zero status. Exiting.", \
                                l_logfile, verbose+)

                            goto crash
                        } else {
                            asky_file = tmpgcorr_sky
                        }
                    } else {
                        asky_file = skyfile
                    }

                    gastat (asky_file, stattype="midpt", statsec=l_statsec, \
                        statextn=auto_sky_stat, fl_mask=l_fl_mask, \
                        badpix=l_badpix, \
                        calc_stddev=no, lower=INDEF, upper=INDEF, nclip=1, \
                        lsigma=3., usigma=3., sci_ext=l_sci_ext, \
                        dq_ext=l_dq_ext, logfile=l_logfile, verbose=no)

                    if (gastat.status != 0) {
                        printlog ("\nERROR - GAREDUCE: GASTAT \
                            returned a non-zero status. Exiting", \
                            l_logfile, verbose+)
                        goto crash
                    } else {
                        skymidpt[sarray_index] = gastat.outstat
                    }
                }
                auto_skyval = skymidpt[sarray_index]

                # Write the SKYMIDPT to the science header
                # Number of extensions checked earlier on.
                for (ii = 1; ii <= nsci; ii += 1) {

                    outextn = gemout//"["//l_sci_ext//","//ii//"]"

                    last_char = stridx(",",auto_skyval) - 1

                    if (last_char <= 0 ) {
                       last_char = strlen (auto_skyval)
                    }

                    auto_val = substr(auto_skyval,1,last_char)

                    if (debug) printlog (auto_val//" "//auto_skyval, \
                        l_logfile, verbose+)

                    if (last_char != strlen(auto_skyval)) {
                        auto_skyval = substr(auto_skyval,last_char+2,\
                                          strlen(auto_skyval))
                    }

                    # Write to header
                    gemhedit (outextn, "SKYMIDPT", auto_val, \
                       sky_comment, delete-)

                } # End of loop over extensions

                # Add back in the sky - no need to worry about VAR/DQ only
                # adding a constant - MS

                if (debug) printlog ("outname: "//outname//" gemout: "//\
                    gemout, l_logfile, verbose+)

                gemexpr ("a + "//skyconst_expr, outname, gemout, \
                    var_expr="a[VAR]", dq_expr="a[DQ]", \
                    sci_ext=l_sci_ext, var_ext=l_var_ext, dq_ext=l_dq_ext, \
                    mdf_ext="MDF", fl_vardq=dovardq, \
                    dims="default", intype="default", outtype="ref", \
                    refim="a", rangecheck=yes, verbose=no, exprdb="none", \
                    logfile=l_logfile)

                if (gemexpr.status != 0) {
                    printlog ("ERROR - GAREDUCE: GEMEXPR returned a non-zero \
                        status. Exiting.", l_logfile, verbose+)
                    imdelete (outname, verify-, >& "dev$null")
                    goto crash
                } else {
                    imdelete (gemout, verify-, >& "dev$null")
                }

            } # End of auto sky

            # Add to the tally of files created
            num_created += 1

            # Update the PHU of the outfile
            outphu = outname//"[0]"
            if (candark) {
                fparse (darkname)
                gemhedit(outphu, "DARKIMG", fparse.root, \
                    "Dark image used by GAREDUCE", delete-)
                meta_proc = "d"
            } else if (isdarksub) {
                meta_proc = "d"
            }

            if (cansky) {
                fparse (skyfile)
                gemhedit (outphu, "SKYIMG", fparse.root, \
                    "Sky image used by GAREDUCE", delete-)
                meta_proc = meta_proc//"s"
            } else if (isskysub) {
                meta_proc = meta_proc//"s"
            }

            if (canflat) {
                fparse (flatname)
                gemhedit (outphu, "FLATIMG", fparse.root, \
                    "Flat image used by GAREDUCE", delete-)
                meta_proc = meta_proc//"f"
            } else if (isflat) {
                meta_proc = meta_proc//"f"
            }

            if (cangainmult) {
                gemhedit (outphu, "BUNIT", "electrons", \
                    "Units of pixel data arrays", delete-)
                meta_proc = meta_proc//"g"


                num_calc = 1
                if (dovardq) {
                   num_calc = 2
                }

                calcextn[1] = "["//l_sci_ext
                calcextn[2] = "["//l_var_ext

                l_units = "electrons"

                for (ii = 1; ii <= num_calc; ii += 1) {

                    if (ii == 2) {
                        l_units = l_units//"*"//l_units
                    }

                    for (jj = 1; jj <= nsci; jj += 1) {
                        outextn = outname//calcextn[ii]//","//jj//"]"

                        gemhedit (outextn, "BUNIT", l_units, "", \
                            delete-)

                        keypar (outextn, l_out_key_gain, silent+)
                        l_ingain = real(keypar.value)

                        if (l_out_key_gain == "GAINORIG") {
                            keypar (outextn, "GAIN", silent+)
                            gemhedit (outextn, "GAIN", \
                                (real(keypar.value) / l_ingain), "", delete-)
                        } else {
                            gemhedit (outextn, "GAIN", 1., "", delete-)
                        }

                        gemhedit (outextn, "GAINMULT", l_ingain, \
                            "GAIN value used to convert to electrons", delete-)
                    }
                }
            } else if (isgainmult) {
                meta_proc = meta_proc//"g"
            }

            # Update metaconf
            if (update_metaconf) {
                out_metaconf = ""
                out_metaconf = substr(metaconf,1,strstr(imgtype,metaconf)-1)
                out_metaconf = out_metaconf//imgtype//"_"//meta_proc//"+"//\
                    proc_test

                gemhedit (outphu, "METACONF", out_metaconf, "", delete-)
                printlog ("GAREDUCE: Updated METACONF keyword value is:"//\
                      "\n          "//out_metaconf, l_logfile, l_verbose)
            }

            gemdate (zone="UT")
            gemhedit (outphu, "GAREDUCE", gemdate.outdate, \
                "UT Time stamp for GAREDUCE", delete-)
            gemhedit (outphu, "GEM-TLM", gemdate.outdate, \
                "UT Last modification with GEMINI", delete-)

            printlog ("GAREDUCE: Created "//outname//"\n--", \
                l_logfile, verbose+)

NEXT_IMAGE:

            if (debug) {
                date ("+%H:%M:%S.%N%n") | scan (l_time_struct)
                printf ("____TIMESTAMP - %s finished: %s\n", \
                    inname, l_time_struct, >> l_logfile)
            }

        } # End of while of current configuration input list

        # fl_keepsky allows keeping the sky frames at the end of each config
        # loop. if set to no, remove not only the skyframe but also the
        # skylists, so the task will not think the sky already exists, when it
        # has been deleted.
        if (access(tmpkeeplst)) {
            if (!l_fl_keepsky) {
                delete ("@"//tmpkeeplst, verify-, >& "dev$null")
            }
            delete (tmpkeeplst,verify-, >& "dev$null")
        }

        printlog ("\nGAREDUCE: Finished processing configuration "//\
            metaconf//"\n----", l_logfile, l_verbose)

NEXT_CONFIGURATION:

        if (imaccess(tmpdark)) {
            imdelete (tmpdark, verify-, >& "dev$null")
        }

        if (imaccess(tmpflat)) {
            imdelete (tmpflat, verify-, >& "dev$null")
        }

        if (imaccess(tmpsky)) {
            imdelete (tmpsky, verify-, >& "dev$null")
        }

        # Counter for number of configurations for creating lists
        nconf += 1

        total_num_created += num_created

        if (debug) printlog ("____deleting tmp configuration files", \
            l_logfile, verbose+)

        delete (conflist, verify-, >& "dev$null")

    } # End of loop over unique configurations

    goto clean

    #--------------------------------------------------------------------------

crash:
    # Exit with error subroutine
    status = 1

clean:

    delete (filelist, verify-, >& "dev$null")

    delete (nextlist//","//tmpconflist//","//tmpflat_caltab//","//\
        tmpexprdb//","//tmpkeeplst//","//conflist//","//\
        tmpfcal//","//tmpdcal//","//tmpkeepskylst//","//skylist//","//\
        tmplst, verify-, \
        >& "dev$null")

    imdelete (tmpdark//","//tmpflat//","//tmpsky//","//tmpout//","//\
        tmpgcorr_sky, verify-, >& "dev$null")

    if (access(outgaimchk_list)) {
        delete ("@"//outgaimchk_list//","//outgaimchk_list, verify-, \
            >& "dev$null")
    }

    scanfile = ""
    scanfile2 = ""
    scanfile3 = ""

    # Print finish time
    gemdate(zone="local")
    printlog ("\nGAREDUCE -- Finished: "//gemdate.outdate, \
        l_logfile, l_verbose)

    if (status == 0) {
        printlog ("\nGAREDUCE -- Exit status:  GOOD", l_logfile, l_verbose)
    } else {
        printlog ("\nGAREDUCE -- Exit status:  ERROR", l_logfile, l_verbose)

    }

    if (status != 0) {
        printlog ("         -- Please read the logfile \""//l_logfile//\
                "\" for more information.", l_logfile, l_verbose)
    }

    printlog ("------------------------------------------------------------", \
        l_logfile, l_verbose)

end

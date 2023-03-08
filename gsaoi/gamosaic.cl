# Copyright(c) 2010-2015 Association of Universities for Research in Astronomy, Inc.

procedure gamosaic (inimages)

##M TODO -- FIX DATASEC and CCDSEC / DETSEC?

# THIS TASK IS DESIGNED FOR GSAOI COMMISSIONING WORK
# The algorithm provides a correct reduction, but may not be optimal AND error
# handling in the task is not up to package standards

# Takes GSAOI raw or GAPREPAREd data, 4 extensions -
#     - FULL FRAME ONLY if fl_paste-

# This task currently does:
# - mosaic the GSAOI images:
#   The fl_paste+ simply piles them up as
#   [3]  [4]
#   [2]  [1]
#   with a predefined gap between each quadrant
#   The fl_paste- applies shift and rotation based on WCS calculation
# - create the mosaiced DQ plane if VAR/DQ requested
# - gap is defined as integer - no fractional pixels!
# - ROIs only work with fl_paste+

# AFTER COMMISSIONING:
# correct propagation of variance plane in the gap region? Is VAR=0 correct?

# the output files are appended a "m" prefix
# a keyword GAMOSAIC is added with the time stamp when the task is last run
# should add a SATURATIO keyword in the header with the lower value from the
# four detectors?
# Need to add proper values of GAIN and RDNOISE to the PHU (take those of the
# reference detector for consistency with gmosaic).

# the input can be in the form of
# - @list containing a list of files one per line, no .fits required, but
#   should not complain about it either
# - comma-separated list of images
# - one or more range of image numbers, in the format allowed by gemlist, in
# which case current UT date is assumed for the rootname, or can be given as
# a parameter

# Note that there is NO outimages list - it just keeps the original filename
# Extension [2] is the reference WCS

# Version 2010 Jun 09 V0.0 CW - created
#         2010 Jun 11 V1.0 CW - added bias subtraction (based on users manual)
#         2011 Mar 10 V1.5 CW - removed normalization and bias subtraction. for
#                               quick visualization, gadisplay should be used.
#         2011 Oct 06 V1.6 CW - ROI handling, also fixed the gaps to be 131pix
#                               when not trimmed
#         2011 Dec 10 V1.7 CW - added shift and rotation based on WCS; moved
#                               gaprepare flags to parameters; fixed crash when
#                               output image already existed.
#         2012 Jan 06 V1.8 CW - modified the fl_paste- image to use the new WCS
#                               for up-looking port
#         2012 Feb 24 V1.9 CW - modified the fl_paste- to handle non-uniform
#                               gaps; fixed DQ transformation and included gaps
#                               -- fl_paste- is NOT correct for ARRAY and DET
#                               ROIs!!
#         2012 Mar 19 V2.0 CW - added graceful exit for ROIs with fl_paste-
#         2012 Apr 03 V2.1 CW - added a fl_fluxcon option for geotran. Defaults
#                               to NO - there is something here weird with the
#                               way geotran works...
#
# Future revisions recorded in CVS logs

######################################################

# Inputs
char    inimages    {prompt="Input GSAOI images or list"}
char    rawpath     {"",prompt="Path for input raw images"}
char    outpref     {"m",prompt="Prefix for output images"}
char    rootname    {"",prompt="Root name for images; blank=today UT"}
bool    fl_paste    {no,prompt="Paste images instead of mosaic"}
bool    fl_vardq    {no,prompt="Create variance and data quality frames"}
#bool    fl_distcor  {no,prompt="Correct array distortion (vs tilt/rotation)?"}
char    gapfile     {"",prompt ="Bad Pixel file with gap information"}
int     colgap      {153,prompt="Average vertical gap between arrays (trimmed, tiling only)"}
int     linegap     {148,prompt="Average horizontal gap between arrays (trimmed, tiling only)"}
#bool    fl_fixpix   {no,prompt="Interpolate across chip gaps - NOT recommended"}
char    geointer    {"linear",enum="linear|nearest|poly3|poly5|spline3|sinc",prompt="Interpolation method to use with geotran"}
bool    fl_fluxcon  {no,prompt="Preserve image flux in geotran?"}
char    key_datasec {"DATASEC",prompt="Header keyword for data section"}
char    key_ccdsec  {"CCDSEC",prompt="Header keyword for array section"}
# gaprepare
char    gaprep_pref {"g",prompt="Prefix for gaprepare output images"}
bool    fl_trim     {yes,prompt="Trim the images?"}
bool    fl_nlc      {yes,prompt="Apply non-linear correction to each array?"}
bool    fl_sat      {no,prompt="Include non-linear and saturated pixels in BPM"}
char    badpix      {"gsaoi$data/gsaoibpm_high_full.fits",prompt="Static Bad Pixel Mask - not mosaiced"}
char    arraysdb    {"gsaoi$data/gsaoiAMPS.dat",prompt="Database file containing array information"}
char    non_lcdb    {"gsaoi$data/gsaoiNLC.dat",prompt="Database file containing non-linearity coefficients"}
# Other
char    sci_ext     {"SCI",prompt="Name of science extensions"}
char    var_ext     {"VAR",prompt="Name of variance extensions"}
char    dq_ext      {"DQ",prompt="Name of data quality extensions"}
char    logfile     {"",prompt="Logfile"}
bool    verbose     {yes,prompt="Verbose?"}
int     status      {0,prompt="Exit status (0=good)"}
struct  *scanfile   {"",prompt="Internal use only"}

######################################################

begin

    ####
    # Variable declaration
    char l_inimages, l_rawpath, l_outpref, l_rootname
    char l_key_datasec, l_logfile, filelist, l_gapfile
    char utdate, img, inimg, outimg, sciphu
    char normsec,scitile, vartile, dqtile, datasec
    char tmpout, tmpbpm, tmphead, l_dir, l_badpix
    char junk1, junk2, l_geointer, l_key_ccdsec
    char l_temp, tmpscpy, tmpvcpy, tmpdcpy, geodbfile
    char geodbtran[4]
    char tmposci[4], tmpovar[4], tmpodq[4]
    char tmptsci[4], tmptvar[4], tmptdq[4]
    char tmpjnk
    char outgaimchk_list, l_sci_ext, l_var_ext, l_dq_ext, l_arraysdb
    char l_non_lcdb, l_gaprep_pref, nextlist, tmpfile, tmpgaprep, tmp2file
    char l_key_allowed, l_key_forbidden, l_obstype_allowed, l_object_allowed
    char tmplist, img_root, inname, inphu, tmpnotproc, metaconf, roiid

    int ii, l_test, l_colgap, l_linegap
    int nctile, nltile, nc1, nc2, nl1, nl2
    int gc1, gc2, gl1, gl2, gc4, gl4, jnk1, jnk2
    int jnk3, jnk4, jnk5, jnk6
    int newxgap, newygap, newnx, newny, addg
    int pos_type, lcount

    real gxshift[4],gyshift[4],rotation[4],gmag[4]

    bool l_fl_vardq, l_fl_fixpix, l_verbose, l_fl_paste, l_fl_orig_paste
    bool l_fl_trim, l_fl_nlc, l_fl_sat, l_fl_distcor, fixed_roiid
    bool fakevar, tmpflp, l_fl_fluxcon, someprepared, debug

    struct pr_struct

    # There is a hiercahy to this list: upto and including TWLT are added by
    # gaprepapre with the rest being added by other tasks - MS
    char types[7]="OBJ","DOME","GCAL","TWLT","DARK","SKY","FLAT"
    int num_types=7

    ####
    # local variables
    l_inimages = inimages
    l_rawpath = rawpath
    l_outpref = outpref
    l_rootname = rootname
    l_key_datasec = key_datasec
    l_logfile = logfile
    l_fl_vardq = fl_vardq
    #l_fl_fixpix = fl_fixpix
    l_gapfile = gapfile
    l_badpix = badpix
    l_verbose = verbose
    l_colgap = colgap
    l_linegap = linegap
    l_fl_orig_paste = fl_paste
    l_geointer = geointer
    l_key_ccdsec = key_ccdsec
    l_fl_trim = fl_trim
    l_fl_nlc = fl_nlc
    l_fl_sat = fl_sat
    l_fl_fluxcon = fl_fluxcon
    l_gaprep_pref = gaprep_pref
    l_sci_ext = sci_ext
    l_var_ext = var_ext
    l_dq_ext = dq_ext
    l_arraysdb = arraysdb
    l_non_lcdb = non_lcdb

    ####
    # Temporary files
    filelist = mktemp("tmpfile")
    nextlist = mktemp ("tmpnextlist")
    outgaimchk_list = mktemp ("tmpgaimchk")
    tmpgaprep = mktemp ("tmpgaprep")
    tmpout = mktemp("tmpout")
    tmpbpm = mktemp("tmpbpm")
    tmphead = mktemp("tmphead")
    tmpscpy = mktemp("tmpcpys")
    tmpvcpy = mktemp("tmpcpyv")
    tmpdcpy = mktemp("tmpcpyd")
    tmplist = mktemp ("tmplist")
    tmpnotproc = mktemp("tmpnotproc")
    tmpjnk = ""

    ####
    # Default values
    debug = no
    status = 0
    l_test = 0
    nctile = 2
    nltile = 2
    fakevar = no
    someprepared = no
    fixed_roiid = no

    # Set like this for a reason ##M Is this needed?
    l_fl_fixpix = no

    # Set this to yes if we can tranform ROIs - MS
    fixed_roiid = no

    # Transformation is relative to Array [2] (lower left corner)

    # Original calculation for the SIDE PORT (5)
    #rotation[1] = 0.1037
    #rotation[2] = 0.0   # reference!
    #rotation[3] = -0.3362
    #rotation[4] = 0.3066

    #gxshift[1] = -2.8905
    #gxshift[2] = 0.0   # reference!
    #gxshift[3] = 10.3448
    #gxshift[4] = 4.3011

    #gyshift[1] = -10.7355
    #gyshift[2] = 0.0   # reference!
    #gyshift[3] = -7.4838
    #gyshift[4] = -9.7460

    # Instrument in the UPLOOKING PORT, values as of 2012 Feb 24
    rotation[1] = 0.033606
    rotation[2] = 0.0  # reference
    rotation[3] = -0.582767
    rotation[4] = -0.769542

    gxshift[1] = 43.604018
    gxshift[2] = 0.0 # reference
    gxshift[3] = 0.029961
    gxshift[4] = 43.420524

    gyshift[1] = -1.248663
    gyshift[2] = 0.0  # reference
    gyshift[3] = 41.102256
    gyshift[4] = 41.722921

    gmag[1] = 1.0  # 1.0013
    gmag[2] = 1.0  # reference
    gmag[3] = 1.0  # 1.0052
    gmag[4] = 1.0  # 1.0159


    geodbtran[1] = "gsaoisol_det1"
    geodbtran[2] = "gsaoisol_det2"
    geodbtran[3] = "gsaoisol_det3"
    geodbtran[4] = "gsaoisol_det4"

    ########
    # Here is where the actual work starts

#    cache ()

    # Test the log file
    if ((l_logfile == "") || (stridx(" ", l_logfile) > 0)) {
       l_logfile = gsaoi.logfile
       if ((l_logfile == "") || (stridx (" ", l_logfile) > 0)) {
          l_logfile = "gsaoi.log"
          printlog ("WARNING - GAMOSAIC: Both gaprepare.logfile and \
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
    printlog ("GAMOSAIC -- Started: "//gemdate.outdate//"\n", \
        l_logfile, l_verbose)

    # Record all input parameters
    printlog ("GAMOSAIC: Input parameters...\n", l_logfile, l_verbose)
    printlog ("    inimages    = "//l_inimages, l_logfile, l_verbose)
    printlog ("    rawpath     = "//l_rawpath, l_logfile, l_verbose)
    printlog ("    outpref     = "//l_outpref, l_logfile, l_verbose)
    printlog ("    rootname    = "//l_rootname, l_logfile, l_verbose)
    printlog ("    fl_paste    = "//l_fl_orig_paste, l_logfile, l_verbose)
    printlog ("    fl_vardq    = "//l_fl_vardq, l_logfile, l_verbose)
#    printlog ("    fl_distcor  = "//l_fl_distcor, l_logfile, l_verbose)
    printlog ("    gapfile     = "//l_gapfile, l_logfile, l_verbose)
    printlog ("    colgap      = "//l_colgap, l_logfile, l_verbose)
    printlog ("    linegap     = "//l_linegap, l_logfile, l_verbose)
#    printlog ("    fl_fixpix   = "//l_fl_fixpix, l_logfile, l_verbose)
    printlog ("    geointer    = "//l_geointer, l_logfile, l_verbose)
    printlog ("    fl_fluxcon  = "//l_fl_fluxcon, l_logfile, l_verbose)
    printlog ("    badpix      = "//l_badpix, l_logfile, l_verbose)
    printlog ("    gaprep_pref = "//l_gaprep_pref, l_logfile, l_verbose)
    printlog ("    fl_trim     = "//l_fl_trim, l_logfile, l_verbose)
    printlog ("    fl_nlc      = "//l_fl_nlc, l_logfile, l_verbose)
    printlog ("    fl_vardq    = "//l_fl_vardq, l_logfile, l_verbose)
    printlog ("    fl_sat      = "//l_fl_sat, l_logfile, l_verbose)
    printlog ("    arraysbd    = "//l_arraysdb, l_logfile, l_verbose)
    printlog ("    non_lcdb    = "//l_non_lcdb, l_logfile, l_verbose)
    printlog ("    key_datasec = "//l_key_datasec, l_logfile, l_verbose)
    printlog ("    key_ccdsec  = "//l_key_ccdsec, l_logfile, l_verbose)
    printlog ("    sci_ext     = "//l_sci_ext, l_logfile, l_verbose)
    printlog ("    var_ext     = "//l_var_ext, l_logfile, l_verbose)
    printlog ("    dq_ext      = "//l_dq_ext, l_logfile, l_verbose)
    printlog ("    logfile     = "//l_logfile, l_logfile, l_verbose)
    printlog ("", l_logfile, l_verbose)

    ####
    # Check parameters...
    # Check outpref
    if (l_outpref == "" || stridx(" ",l_outpref) > 0) {
        printlog ("WARNING - GAMOSAIC: outpref is an empty string or contains \
            spaces. Setting to default \"m\"", l_logfile, verbose+)
        l_outpref = "m"
    }

    # Check gaprepare out prefix
    if (l_gaprep_pref == "" || stridx(" ",l_gaprep_pref) > 0) {
        printlog ("WARNING - GAMOSAIC: outpref is an empty string or contains \
            spaces. Setting to default \"g\"", l_logfile, verbose+)
        l_gaprep_pref = "g"
    }

    if (l_fl_vardq) {
        if (l_dq_ext == "" || (stridx(" ",l_dq_ext) > 0) || \
            l_var_ext == "" || (stridx(" ",l_var_ext) > 0)) {

            printlog("WARNING - GAMOSAIC: var_ext or dq_ext have not been \
                set.\n                     Output image will not have "//\
                l_var_ext//" or "//l_dq_ext//" planes.", l_logfile,verbose+)
            l_fl_vardq=no
        }
    }

    # Check the l_badpix value
    if ((l_badpix == "") || (stridx(" ", l_badpix) > 0)) {
        # This assumes that the bpm covers all of the areas readout in each
        # detector in each image
        l_badpix = "none"
    } else if (!imaccess(l_badpix)) {
        printlog ("ERROR - GAMOSAIC: cannot access bpm "//l_badpix//\
            ". Exiting", l_logfile, verbose+)
        goto crash
    } else {
        # Check it hasn't been mosaiced
        keypar (l_badpix//"[0]", "GAMOSAIC", silent+)
        if (keypar.found) {
            printlog ("ERROR - GAMOSAIC: "//l_badpix//" has been mosaiced."//\
                " Please supply a static BPM that has not been mosaiced.", \
                l_logfile, verbose)
            goto crash
        }
    }
    printlog ("GAMOSAIC: badpix value after checks is "//l_badpix, \
        l_logfile, l_verbose)

    # Check input images...
    printlog ("\nGAMOSAIC: Calling GAIMCHK to check input files...", \
        l_logfile, l_verbose)

    # In decending processing steps... keep them this way for GAIMCHK - -MS
    l_key_allowed = "GAPREPAR"

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
        fl_out_check=no, fl_vardq_check=l_fl_vardq, sci_ext=l_sci_ext, \
        var_ext=l_var_ext, dq_ext=l_dq_ext, outlist=outgaimchk_list, \
        logfile=l_logfile, verbose=no)

    if (gaimchk.status != 0) {
        printlog ("ERROR - GAMOSAIC: GAIMCHK returned a non-zero status. \
            Exiting.", l_logfile, verbose+)
        goto crash
    } else if (!access(outgaimchk_list)) {
        printlog ("ERROR - GAMOSAIC: Cannot access output list from \
            GAIMCHK", l_logfile, verbose+)
        goto crash
    }

    printlog ("GAMOSAIC: Returned from GAIMCHK.\n", l_logfile, l_verbose)

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

        # Now check if there were any GAPREPAREd files too
        if (someprepared) {
            printlog ("WARNING - GAMOSAIC: The following files have \
                already been processed by GAPREPARE.\n", \
                l_logfile, verbose+)
            fields (tmp2file, fields=1, lines="1-", quit_if_miss=no, \
                print_file_names=no) | \
                tee (l_logfile, out_type="text", append+)
        }

        printlog ("\nGAMOSAIC: Calling GAPREPARE to process unprepared input \
            files", l_logfile, l_verbose)

        gaprepare ("@"//tmpgaprep, rawpath="",outpref=l_gaprep_pref, \
            rootname="", fl_trim=l_fl_trim, fl_nlc=l_fl_nlc, \
            fl_vardq=l_fl_vardq, fl_sat=l_fl_sat, badpix=l_badpix, \
            sci_ext=l_sci_ext, var_ext=l_var_ext, dq_ext=l_dq_ext, \
            arraysdb=l_arraysdb, non_lcdb=l_non_lcdb, \
            logfile = l_logfile, verbose=l_verbose)

        if (gaprepare.status != 0) {
            printlog ("ERROR - GAMOSAIC: GAPREPARE returned an non-zero \
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
        printlog ("ERROR - GAMOSAIC: No input images can be used. \
            Please try again.", l_logfile, verbose+)
        goto crash
    }


    # Loop over the input images
    scanfile = ""
    scanfile = nextlist
    while (fscan(scanfile, inimg) != EOF) {

        fparse (inimg)
        img_root = fparse.root
        inname = img_root//fparse.extension
        inphu = inimg//"[0]"
        outimg = l_outpref//img_root//".fits"

        # Double check that you can access the file...
        if (!imaccess(inimg)) {
            printlog ("WARNING - GAMOSAIC: Image "//inimg//\
                " does not exist. Skipping this image.", l_logfile, verbose)
            print (inimg, >> tmpnotproc)
            goto NEXT_IMAGE
        } else if (imaccess(outimg)) {
            printlog ("ERROR - GAMOSAIC: Image "//outimg//\
                " already exists. Exiting.", l_logfile, verbose)
            goto crash
        }

        printlog ("----\nGAMOSAIC: Processing image - "//inname//"...", \
            l_logfile, l_verbose)

        # Read the METACONFIG keyword
        keypar (inphu, "METACONF", silent+)
        if (!keypar.found) {
            printlog ("WARNING - GAMOSAIC: METACONF keyword not found \
                in "//inimg//" Skipping this image.", l_logfile, verbose+)
            print (inimg, >> tmpnotproc)
            goto NEXT_IMAGE
        } else {
            metaconf = keypar.value

            # Check it's valid
            if (metaconf == "UNKOWN" || metaconf == "" || \
                (stridx(" ",metaconf) > 0)) {
                printlog ("WARNING - GAMOSAIC: METACONF keyword is "//\
                    "\"UNKNOWN\" or bad in"//inimg//".\n"//\
                    "                  Skipping this image.", \
                    l_logfile, verbose+)
                print (inimg, >> tmpnotproc)
                goto NEXT_IMAGE
            }
        }

        if (debug) printlog ("____pre_roiid: "//metaconf, l_logfile, verbose+)

        # Parse the ROIID - and determine what kind of image it is
        for (ii = 1; ii <= num_types; ii += 1) {

            pos_type = strstr(types[ii],metaconf)

            if (debug) printlog ("____ii: "//ii//", pos_type: "//pos_type, \
                l_logfile, verbose+)

            if (pos_type > 0) {
                roiid = substr(metaconf,pos_type,strlen(metaconf))
                roiid = substr(roiid,stridx("+",roiid)+1,strlen(roiid))
                if (types[ii] == "GCAL" || types[ii] == "DOME") {
                    roiid = substr(roiid,stridx("+",roiid)+1,strlen(roiid))
                }
                roiid = substr(roiid,1,stridx("+",roiid)-1)
                break
            } else if (ii == num_types) {

                ##M May need to fix this
                # It's likely to be a dark by this point
                roiid = metaconf
                roiid = substr(roiid,stridx("+",roiid)+1,strlen(roiid))
                roiid = substr(roiid,1,stridx("+",roiid)-1)

                if (debug) printlog ("____roiid(2): "//roiid, l_logfile, \
                    verbose+)

                # If it is a DARK the roiid should begin with either F,D,A.
                if (stridx("FDA",substr(roiid,1,1)) == 0) {
                    printlog ("ERROR - GAMOSAIC: Unable to determine type and \
                        hence ROIID. Exiting", l_logfile, verbose+)
                    goto crash
                }
            }
        }

        if (debug) printlog ("____roiid: "//roiid, l_logfile, verbose+)

        # Check if ROI and fl_paste- in which case turn it off (does not work)
        ##M Can we not fix this? -- Rotation relative to centre of array for FF
        ##M need to adjust centre for ROIs I think - MS
        l_fl_paste = l_fl_orig_paste

        if (strstr("FF",roiid) == 0 && !l_fl_paste && !fixed_roiid) {
            printlog ("WARNING - GAMOSAIC: Correcting tilt/rotation only \
                works with"//\
                "\n                    full-frame images. Setting \
                fl_paste=yes.", l_logfile, l_verbose)
            l_fl_paste = yes
        }

        if (l_fl_paste) {

            printlog ("GAMOSAIC: Not correcting for tilt/rotation.", \
                l_logfile, l_verbose)

            # full frame and ROIs - re-calculate the image size and the gaps
            # for fl_paste+
            # if full frame or DET ROI (trimmed) gaps are the original
            # if full frame or DET ROI (untrimmed) gaps are the original-8 (the
            # inner four non-illuminated
            # pixels in each side are part of the measured gap)
            # if ARRAY ROI (never trimmed), ncols/nlines = 2040+size in x or
            # y+original gap, gap = 2040-size+original gap

            keypar (inimg//"["//l_sci_ext//",1]",l_key_ccdsec, silent+)
            junk1 = substr(keypar.value,1,3)
            if (junk1 == "[1:" || junk1 == "[5:") {

               # this is either full frame or DET ROI
                keypar (inimg//"[0]", "TRIMMED", silent+)
                if (keypar.found) {
                    # image has been trimmed
                    addg = 0
                } else {
                    addg = 8
                }

                newnx = INDEF
                newny = INDEF
                newxgap = l_colgap - addg
                newygap = l_linegap - addg

            } else {
                # this is ARRAY ROI, and I need the actual size
                hselect (inimg//"["//l_sci_ext//",1]", "naxis[12]", \
                    expr="yes") | scan (jnk1,jnk2)
                newnx = 2040 + jnk1 + l_colgap
                newny = 2040 + jnk2 + l_linegap
                newxgap = 2040 - jnk1 + l_colgap
                newygap = 2040 - jnk2 + l_linegap
            }

            # make the mosaic
            scitile = inimg//"["//l_sci_ext//",2],"//\
                inimg//"["//l_sci_ext//",1],"//inimg//\
                "["//l_sci_ext//",3],"//inimg//"["//l_sci_ext//",4]"
            vartile = inimg//"["//l_var_ext//",2],"//\
                inimg//"["//l_var_ext//",1],"//inimg//\
                "["//l_var_ext//",3],"//inimg//"["//l_var_ext//",4]"
            dqtile =  inimg//"["//l_dq_ext//",2],"//inimg//\
                "["//l_dq_ext//",1],"//inimg//\
                "["//l_dq_ext//",3],"//inimg//"["//l_dq_ext//",4]"

            # imtile the science
            printlog ("GAMOSAIC: Creating the output image "//outimg, \
                l_logfile, l_verbose)

            fxcopy (inimg, outimg, groups="0", new_file=yes, verbose-)

            printlog ("GAMOSAIC: Tiling the "//l_sci_ext//" extensions", \
                l_logfile, l_verbose)

            imtile (scitile, outimg//"[1,append]", nctile, nltile, \
               ncoverlap=-newxgap, nloverlap=-newygap, \
               ncols = newnx, nlines = newny, trim_section="", \
               missing_input="", start_tile="ll", \
               row_order=yes, raster_order=no, median_section="", \
               opixtype="r", ovalue=0.0, verbose-)

            if (l_fl_vardq) {

                # imtile the variance
                # imtile the DQ

                printlog ("GAMOSIAC: Tiling the "//l_var_ext//" and "//\
                    l_dq_ext//" extensions", l_logfile, l_verbose)

                imtile (vartile, outimg//"[2,append]", nctile, nltile, \
                     ncoverlap=-newxgap, nloverlap=-newygap, \
                     ncols = newnx, nlines = newny, trim_section="", \
                     missing_input="", start_tile="ll", \
                     row_order=yes, raster_order=no, median_section="", \
                     opixtype="r", ovalue=0.0, verbose-)

                imtile (dqtile, outimg//"[3,append]", nctile, nltile, \
                     ncoverlap=-newxgap, nloverlap=-newygap, \
                     ncols = newnx, nlines = newny, trim_section="", \
                     missing_input="", start_tile="ll", \
                     row_order=yes, raster_order=no, median_section="", \
                     opixtype="r", ovalue=0.0, verbose-)
            }

            # Store the header of EXTENSION 2 so the WCS is correct at the end
            imhead (inimg//"["//l_sci_ext//",2]", long+) | match ("CTYPE1", \
                "STDIN", stop-, >> tmphead)
            imhead (inimg//"["//l_sci_ext//",2]", long+) | match ("CRPIX1", \
                "STDIN", stop-, >> tmphead)
            imhead (inimg//"["//l_sci_ext//",2]", long+) | match ("CRVAL1", \
                "STDIN", stop-, >> tmphead)
            imhead (inimg//"["//l_sci_ext//",2]", long+) | match ("CTYPE2", \
                "STDIN", stop-, >> tmphead)
            imhead (inimg//"["//l_sci_ext//",2]", long+) | match ("CRPIX2", \
                "STDIN", stop-, >> tmphead)
            imhead (inimg//"["//l_sci_ext//",2]", long+) | match ("CRVAL2", \
                "STDIN", stop-, >> tmphead)
            imhead (inimg//"["//l_sci_ext//",2]", long+) | match ("CD1_1", \
                "STDIN", stop-, >> tmphead)
            imhead (inimg//"["//l_sci_ext//",2]", long+) | match ("CD1_2", \
                "STDIN", stop-, >> tmphead)
            imhead (inimg//"["//l_sci_ext//",2]", long+) | match ("CD2_1", \
                "STDIN", stop-, >> tmphead)
            imhead (inimg//"["//l_sci_ext//",2]", long+) | match ("CD2_2", \
                "STDIN", stop-, >> tmphead)
            imhead (inimg//"["//l_sci_ext//",2]", long+) | match ("GAIN1", \
                "STDIN", stop-, >> tmphead)
            imhead (inimg//"["//l_sci_ext//",2]", long+) | match ("RDNOISE1",\
                "STDIN", stop-, >> tmphead)
            imhead (inimg//"["//l_sci_ext//",2]", long+) | match ("RADECSYS",\
                "STDIN", stop-, >> tmphead)
            imhead (inimg//"[0]", long+) | match ("EQUINOX", \
                "STDIN", stop-, >> tmphead)
            imhead (inimg//"[0]", long+) | match ("MJD-OBS", \
                "STDIN", stop-, >> tmphead)

        } else {
            # If I'm going to implement proper correction of tilt/rotation,
            # here is where it goes
            # with an extra if that outputs to another temporary as geotrans
            # cannot write in place
            # optionally, this entire block of code can be colapsed into 4
            # calls to geotran using
            # the example 5 in the help - I may do it once I've figured how it
            # works

            printlog ("GAMOSAIC: Transforming each array for tilt / \
                rotation", l_logfile, l_verbose)

            # transform each array - with these parameters in geotran,
            # the size of the transformed image is the same as the
            # size of the input image

            for (ii = 1; ii <= 4; ii += 1) {
                tmposci[ii] = inimg//"["//l_sci_ext//","//ii//"]"
                tmptsci[ii] = mktemp("tmptsci")

                printlog  ("GAMOSAIC: Transforming EXTVERSION: "//ii, \
                    l_logfile, l_verbose)

                geotran (tmposci[ii], tmptsci[ii], "", "", xin=INDEF, \
                    yin=INDEF, xshift=gxshift[ii], \
                    yshift=gyshift[ii], xout=INDEF, yout=INDEF, xmag=gmag[ii],\
                    ymag=gmag[ii], xrotation=rotation[ii], \
                    yrotation=rotation[ii], xmin=INDEF, xmax=INDEF, \
                    ymin=INDEF, ymax=INDEF, xscale=1., yscale=1., ncols=INDEF,\
                    nlines=INDEF, xsample=1., ysample=1., \
                    interpolant=l_geointer, boundary="constant", constant=0, \
                    fluxcon=l_fl_fluxcon, nxblock=2048, nyblock=2048, verbose-)

                # Don't need to transform var/dq if those are the fake ones.

                if (l_fl_vardq) {
                    tmpovar[ii] = inimg//"["//l_var_ext//","//ii//"]"
                    tmptvar[ii] = mktemp("tmptvar")//".fits"
    #                tmpodq[ii] = inimg//"["//l_dq_ext//","//ii//"]"
                    tmptdq[ii] = mktemp("tmptdq")//".fits"
                    if (!fakevar) {
                        tmpjnk = mktemp("tmpjnk")//".fits"
                        delete (tmpjnk, verify-, >& "dev$null")

                        geotran (tmpovar[ii], tmptvar[ii], "", "", xin=1, \
                            yin=1, xshift=gxshift[ii], \
                            yshift=gyshift[ii], xout=1, yout=1, xmag=gmag[ii],\
                            ymag=gmag[ii], \
                            xrotation=rotation[ii], yrotation=rotation[ii], \
                            xmin=INDEF, xmax=INDEF, \
                            ymin=INDEF, ymax=INDEF, xscale=1., yscale=1., \
                            ncols=INDEF, nlines=INDEF, \
                            xsample=1., ysample=1., interpolant=l_geointer, \
                            boundary="constant", constant=0, \
                            fluxcon=l_fl_fluxcon, nxblock=2048, nyblock=2048, \
                            verbose-)

                        # this is needed to actually get 1 and 0 values in the
                        # DQ - unfortunatelly it loses the sat and non-linear
                        # flags
                        imcopy (inimg//"["//l_dq_ext//","//ii//"]", tmpjnk, \
                            verbose-)

                        imarith (tmpjnk, "*", 1000, tmpjnk, title="", \
                            divzero=0., hparams="", pixtype="", \
                            calctype="", verbose-)

                        # now transform the DQ

                        geotran (tmpjnk, tmptdq[ii], "", "", \
                            xshift=gxshift[ii], yshift=gyshift[ii], \
                            xmag=gmag[ii], \
                            ymag=gmag[ii], xrotation=rotation[ii], \
                            yrotation=rotation[ii], xmin=INDEF, xmax=INDEF, \
                            ymin=INDEF, ymax=INDEF, ncols=INDEF, nlines=INDEF,\
                            interpolant=l_geointer, \
                            boundary="constant", constant=0, \
                            fluxcon=l_fl_fluxcon, nxblock=2048, nyblock=2048, \
                            verbose-)

                        # and change it back to actually contain 0 and 1. The
                        # 500 value is sort of arbitrary.
                        delete (tmpjnk, verify-, >& "dev$null")
                        imexpr ("((a < 500) ? 0 : 1)", tmpjnk, tmptdq[ii], \
                           dims="auto", intype="auto", outtyp="ref", \
                           refim="a", bwidth=0, btype="nearest", \
                           bpixval=0., rangecheck=yes, verbose=no, \
                           exprdb="none")

                        delete (tmptdq[ii], verify-, >& "dev$null")
                        tmptdq[ii] = tmpjnk
                    } else {
                        tmptvar[ii] = tmpovar[ii]
                        tmptdq[ii] = inimg//"["//l_dq_ext//","//ii//"]"
                    }
                }
            }

            # The gaps used here are (TRIMMED):
            # outer corner: 1-2 = 159; 2-3 = 147
            # inner corner: 3-4 = 148; 4-1 = 150px
            # but this includes the blank part in each image
            # with the current transformation, the actual gaps for image
            # size are 116 and 112 - I'll make it square and mask
            # the extra pixels

            # untrimmed gaps need to subtract 8pix

            # full frame
            # x = size*2+116-addg
            # y = size*2+112-addg
            # This DOES NOT work for DET and ARRAY

            # I need to calculate the actual size of the final image - get the
            # size of
            # each array from [1]
            # Optionally this can be gathered from the METACONF

            hselect (inimg//"["//l_sci_ext//",1]", "naxis[12]", expr="yes") \
                | scan (jnk1,jnk2)

            # figure out the ROI and define the size accordingly
            keypar (inimg//"["//l_sci_ext//",1]",l_key_ccdsec, silent+)
            junk1 = substr(keypar.value,1,3)

            if (junk1 == "[1:" || junk1 == "[5:") {
                # this is either full frame, and I need to check for trimming
                keypar (inimg//"[0]", "TRIMMED", silent+)
                if (keypar.found) {
                    # image has been trimmed
                    addg = 0
                } else {
                    addg = 8
                }

                # so the size is
                newnx = 2*jnk1 - addg + 116
                newny = 2*jnk2 - addg + 116
            }

    #else {
    # this is ARRAY ROI, and its never trimmed
    #            addg=0
    # the size of ARRAY ROI is
    ### NOT RIGHT AS THE GAP IS not uniform
    #            newnx = 2040 + jnk1 + 116
    #            newny = 2040 + jnk2 + 116
    #        }

            # create the PHU

            printlog ("\nGAMOSAIC: Creating the output image - "//outimg, \
                l_logfile, l_verbose)

            fxcopy (inimg, outimg, groups="0", new_file=yes, verbose-)

            # create the sci,var,ext extensions with header from [2]
            imexpr ("repl(0,"//newnx//")", output=tmpscpy, \
                dims=newnx//","//newny, \
                intype="auto", outtype="real", refim="auto", bwidth=0, \
                btype="nearest", bpixval=0., rangecheck=yes, verbose=no, \
                exprdb="none")

            if (l_fl_vardq) {

                imexpr ("repl(0,"//newnx//")", output=tmpvcpy, \
                    dims=newnx//","//newny, \
                    intype="auto", outtype="real", refim="auto", bwidth=0, \
                    btype="nearest", bpixval=0., rangecheck=yes, verbose=no, \
                    exprdb="none")

                imexpr ("repl(0,"//newnx//")", output=tmpdcpy, \
                    dims=newnx//","//newny, \
                    intype="auto", outtype="long", refim="auto", bwidth=0, \
                    btype="nearest", bpixval=0., rangecheck=yes, verbose=no, \
                    exprdb="none")
            }

            # now start positioning the pieces where I want them
            # [2] is the reference, lower left corner goes in (1,1)

            printlog ("GAMOSAIC: Creating the mosaicked extensions", \
                l_logfile, l_verbose)

            hselect (tmptsci[2], "naxis[12]", expr="yes") | scan (jnk1,jnk2)
            if (debug) printlog ("____tmptsci[2]: "//jnk1//","//jnk2, \
                l_logfile, verbose+)

            imcopy (tmptsci[2], tmpscpy//"[1:"//jnk1//",1:"//jnk2//"]", \
                verbose-)

            if (l_fl_vardq) {
                imcopy (tmptvar[2], tmpvcpy//"[1:"//jnk1//",1:"//jnk2//"]", \
                    verbose-)
                imcopy (tmptdq[2], tmpdcpy//"[1:"//jnk1//",1:"//jnk2//"]", \
                    verbose-)
            }

            # then put [1] in place: for FF, (1,1) goes to 2040+116-addg, 1

            hselect (tmptsci[1], "naxis[12]", expr="yes") | scan (jnk1,jnk2)

            if (debug) printlog ("____tmptsci[1]: "//jnk1//","//jnk2, \
                l_logfile, verbose+)

            jnk3 = jnk1+116-addg
            jnk4 = jnk1+jnk3
            imcopy (tmptsci[1], \
                tmpscpy//"["//jnk3//":"//jnk4//",1:"//jnk2//"]", verbose-)

            if (l_fl_vardq) {
                imcopy (tmptvar[1], \
                    tmpvcpy//"["//jnk3//":"//jnk4//",1:"//jnk2//"]", verbose-)
                imcopy (tmptdq[1], \
                    tmpdcpy//"["//jnk3//":"//jnk4//",1:"//jnk2//"]", verbose-)
            }

            # now [3]: for FF (1,1) goes to 1, 2040+97-addg  112

            hselect (tmptsci[3], "naxis[12]", expr="yes") | scan (jnk1,jnk2)
            jnk3 = jnk2+97-addg
            jnk4 = jnk2+jnk3-1
            imcopy (tmptsci[3], \
                tmpscpy//"[1:"//jnk1//","//jnk3//":"//jnk4//"]", verbose-)

            if (l_fl_vardq) {
                imcopy (tmptvar[3], \
                    tmpvcpy//"[1:"//jnk1//","//jnk3//":"//jnk4//"]", verbose-)
                imcopy (tmptdq[3], \
                    tmpdcpy//"[1:"//jnk1//","//jnk3//":"//jnk4//"]", verbose-)
            }

            # and [4]: for FF (1,1) goes to 2040+110-addg, 2040+95+addg  77 106

            hselect (tmptsci[4], "naxis[12]", expr="yes") | scan (jnk1,jnk2)
            jnk3 = jnk1+110-addg
            jnk4 = jnk1+jnk3-1
            jnk5 = jnk2+95-addg
            jnk6 = jnk2+jnk5-1
            imcopy (tmptsci[4], \
                tmpscpy//"["//jnk3//":"//jnk4//","//jnk5//":"//jnk6//"]", \
                verbose-)

            if (l_fl_vardq) {
                imcopy (tmptvar[4], \
                    tmpvcpy//"["//jnk3//":"//jnk4//","//jnk5//":"//jnk6//"]", \
                    verbose-)
                imcopy (tmptdq[4], \
                    tmpdcpy//"["//jnk3//":"//jnk4//","//jnk5//":"//jnk6//"]", \
                    verbose-)
            }

            # put it in place
            printlog ("GAMOSAIC: Creating final MEF - "//outimg, \
                l_logfile, l_verbose)

            imcopy (tmpscpy, outimg//"[1,append]", verbose-)
            imdelete (tmpscpy, verify-, >& "dev$null")

            if (l_fl_vardq) {
                imcopy (tmpvcpy, outimg//"[2,append]", verbose-)
                imcopy (tmpdcpy, outimg//"[3,append]", verbose-)
                imdelete (tmpvcpy//","//tmpdcpy, verify-, >& "dev$null")
            }

            # and these are to fix the header later

            imhead (tmptsci[2], long+) | match ("CTYPE1","STDIN", \
                stop-, >> tmphead)
            imhead (tmptsci[2], long+) | match ("CRPIX1","STDIN", \
                stop-, >> tmphead)
            imhead (tmptsci[2], long+) | match ("CRVAL1","STDIN", \
                stop-, >> tmphead)
            imhead (tmptsci[2], long+) | match ("CTYPE2","STDIN", \
                stop-, >> tmphead)
            imhead (tmptsci[2], long+) | match ("CRPIX2","STDIN", \
                stop-, >> tmphead)
            imhead (tmptsci[2], long+) | match ("CRVAL2","STDIN", \
                stop-, >> tmphead)
            imhead (tmptsci[2], long+) | match ("CD1_1","STDIN", \
                stop-, >> tmphead)
            imhead (tmptsci[2], long+) | match ("CD1_2","STDIN", \
                stop-, >> tmphead)
            imhead (tmptsci[2], long+) | match ("CD2_1","STDIN", \
                stop-, >> tmphead)
            imhead (tmptsci[2], long+) | match ("CD2_2","STDIN", \
                stop-, >> tmphead)
            imhead (tmptsci[2], long+) | match ("GAIN","STDIN", \
                stop-, >> tmphead)
            imhead (tmptsci[2], long+) | match ("RDNOISE","STDIN", \
                stop-, >> tmphead)
            imhead (tmptsci[2], long+) | match ("RADECSYS","STDIN", \
                stop-, >> tmphead)
            imhead (inimg//"[0]", long+) | match ("EQUINOX","STDIN", \
                stop-, >> tmphead)
            imhead (inimg//"[0]", long+) | match ("MJD-OBS","STDIN", \
                stop-, >> tmphead)
        }

        # ensure extensions properly named SCI, "//l_var_ext//" and DQ, and
        # numbered all 1

        gemhedit (outimg//"[1]", "EXTNAME", l_sci_ext, "", delete-)
        gemhedit (outimg//"[1]", "EXTVER", 1, "", delete-)

        if (l_fl_vardq) {
            # VAR
            gemhedit (outimg//"[2]", "EXTNAME", l_var_ext, "", delete-)
            gemhedit (outimg//"[2]", "EXTVER", 1, "", delete-)

            # DQ
            gemhedit (outimg//"[3]", "EXTNAME", l_dq_ext, "", delete-)
            gemhedit (outimg//"[3]", "EXTVER", 1, "", delete-)
        }

        # fix headers - need to adjust the WCS, taking one of the extensions as
        # reference
        # note that taking [2] as reference is the easier as it allows to keep
        # CRPIX1 = 1

        ##M TODO -- FIX DATASEC and CCDSEC / DETSEC?

        sciphu = outimg//"[0]"

        printlog ("GAMOSAIC: Updating the header of the output image", \
           l_logfile, l_verbose)

        mkheader (sciphu,tmphead,append+, verbose-)
        mkheader (outimg//"["//l_sci_ext//",1]", tmphead, append+, verbose-)

        if (l_fl_vardq) {
           mkheader (outimg//"["//l_var_ext//",1]", tmphead, append+, verbose-)
           mkheader (outimg//"["//l_dq_ext//",1]", tmphead, append+, verbose-)
        }

        delete (tmphead, verify-, >& "dev$null")

        # flag the gaps if VARDQ or FIXPIX

        if ((l_fl_vardq) || (l_fl_fixpix)) {

            # check the gapfile - if empty, use the defaults - gap position is
            # assumed to be for unbinned after gamosaic
            # the gap positions will  be in the tmpbpm file

            if ((l_gapfile == "") || (stridx(" ", l_gapfile) > 0)) {
                printlog ("WARNING - GAMOSAIC: File with array gap position \
                    does not exist.", l_logfile, l_verbose)
                printlog ("                    Calculating from gap \
                    parameters and image size.\n", l_logfile, l_verbose)

                print ("#GSAOI default gap position, \
                    array unbinned after mosaic", >> tmpbpm)

                # need to read the datasec for extension [2] and calculate the
                # gap position then need to read the size of the final image

                keypar (inimg//"["//l_sci_ext//",2]", l_key_datasec, silent+)
                datasec = keypar.value
                jnk1 = stridx(":",datasec)
                nc1 = int(substr(datasec,2,jnk1-1))
                jnk2 = stridx(",",datasec)
                nc2 = int(substr(datasec,jnk1+1,jnk2-1))
                jnk1 = strldx(":",datasec)
                nl1 = int(substr(datasec,jnk2+1,jnk1-1))
                jnk2 = strlen(datasec)
                nl2 = int(substr(datasec,jnk1+1,jnk2-1))

                keypar (outimg//"["//l_sci_ext//",1]", "i_naxis1", silent+)
                gc4 = int(keypar.value)
                keypar (outimg//"["//l_sci_ext//",1]", "i_naxis2", silent+)
                gl4 = int(keypar.value)

                # Now, the gaps are different if correcting the distortion
                # the fl_paste+ case is easy

                if (l_fl_paste) {

                    gc1 = int(nc2 + 1)
                    gc2 = int(nc2 + newxgap)
                    gl1 = int(nl2 + 1)
                    gl2 = int(nl2 + newygap)

                    print (gc1//"    "//gc2//"    1    "//gl4, >> tmpbpm)
                    print ("1    "//gc4//"    "//gl1//"    "//gl2, >> tmpbpm)
                    printlog ("GAMOSAIC: Using the following for the gaps \
                        positions:\n", l_logfile, l_verbose)
                    printf ("%15d%8d%8d%8d\n", gc1, gc2, 1, gl4) | \
                        scan (pr_struct)
                    printlog (pr_struct, l_logfile, l_verbose)
                    printf ("%15d%8d%8d%8d\n", 1, gc4, gl1, gl2) | \
                        scan (pr_struct)
                    printlog (pr_struct, l_logfile, l_verbose)

                } else {

                    # if correcting the distortion, the gaps are wedged, so I
                    # need to oversize them a bit
                    printlog ("GAMOSAIC: Using the following for the gaps \
                        positions:\n",l_logfile, l_verbose)

                    # gap between [2] and [3]
                    gc1 = int(nc1)
                    gc2 = int(nc2-11)
                    gl1 = int(nl2+1)
                    gl2 = int(nl2+148)
                    print (gc1//"    "//gc2//"    "//gl1//"    "//gl2, \
                        >> tmpbpm)

                    printf ("%15d%8d%8d%8d\n", gc1, gc2, gl1, gl2) | \
                        scan (pr_struct)
                    printlog (pr_struct, l_logfile, l_verbose)

                    # gap between [2] and [1]
                    gc1 = int(nc2+1)
                    gc2 = int(nc2+161)
                    gl1 = int(nl1)
                    gl2 = int(nl2)
                    print (gc1//"    "//gc2//"    "//gl1//"    "//gl2, \
                        >> tmpbpm)

                    printf ("%15d%8d%8d%8d\n", gc1, gc2, gl1, gl2) | \
                        scan (pr_struct)
                    printlog (pr_struct, l_logfile, l_verbose)

                    # gap between [1] and [4]
                    gc1 = int(nc2+162)
                    gc2 = gc4
                    gl1 = int(nl2-1)
                    gl2 = int(nl2+150)
                    print (gc1//"    "//gc2//"    "//gl1//"    "//gl2, \
                        >> tmpbpm)

                    printf ("%15d%8d%8d%8d\n", gc1, gc2, gl1, gl2) | \
                        scan (pr_struct)
                    printlog (pr_struct, l_logfile, l_verbose)

                    # gap between [3] and [4]
                    gc1 = int(nc2-10)
                    gc2 = int(nc2+166)
                    gl1 = int(nl2+151)
                    gl2 = gl4
                    print (gc1//"    "//gc2//"    "//gl1//"    "//gl2, \
                        >> tmpbpm)

                    printf ("%15d%8d%8d%8d\n", gc1, gc2, gl1, gl2) | \
                        scan (pr_struct)
                    printlog (pr_struct, l_logfile, l_verbose)

                    # and the squarish area in the middle
                    gc1 = int(nc2-10)
                    gc2 = int(nc2+161)
                    gl1 = int(nl2+1)
                    gl2 = int(nl2+150)

                    print (gc1//"    "//gc2//"    "//gl1//"     "//gl2, \
                        >> tmpbpm)

                    # also, with the *current transformation* the following
                    # edges need to be masked for full frame only
                    print ("1    18    2188    "//gl4, >> tmpbpm)
                    print ("19    "//gc1-1//"    4177    "//gl4, >> tmpbpm)
                    print (""//gc2+1//"    "//gc4//"    4175    "//gl4, \
                        >> tmpbpm)
                    print ("4190    "//gc4//"    "//gl2+1//"    "//gl4, \
                        >> tmpbpm)
                    print (""//gc4-2//"   "//gc4//"    1    "//gl1-1, \
                        >> tmpbpm)
                }
            } else {
                printlog ("GAMOSAIC: Using "//l_gapfile//" for the gaps \
                    positions:", l_logfile, l_verbose)
                copy (l_gapfile, tmpbpm, verbose-)
            }
        }

        # add the gaps to the DQ - gap position is in the tmpbpm file
        # VAR in the gaps remain as 0 - should it be something else?

        if (l_fl_vardq) {
            printlog ("\nGAMOSAIC: Updating the "//l_dq_ext//" plane with \
                the gaps", l_logfile, l_verbose)

            badpiximage (tmpbpm, outimg//"["//l_sci_ext//",1]", tmpbpm//".pl",\
                goodvalue=0,badvalue=1)

            imexpr ("((b > 0) ? c : a)", outimg//"[3,overwrite]", \
                outimg//"[3]", tmpbpm//".pl", 1, dims="auto", intype="auto", \
                outtyp="ref", refim="a", bwidth=0, btype="nearest", \
                bpixval=0., rangecheck=yes, verbose=no, exprdb="none")
        }

        # fix chip gaps -  use badpiximage and fixpix
        # the variance of the gaps remain as 0.

        # this is not a good idea, so I hardcoded this one to l_fl_fixpix=no

        #    if (l_fl_fixpix) {
        #        badpiximage (tmpbpm,outimg//"["//l_sci_ext//",1]", \
        #            tmpbpm//".pl", \
        #            goodvalue=0,badvalue=1)
        #        proto.fixpix (outimg//"["//l_sci_ext//",1]", tmpbpm//".pl", \
        #            linter="INDEF", cinterp="INDEF", verbose-, pixels-)
        #        if (l_fl_vardq) {
        #            imexpr ("((b > 0) ? c : a)", outimg//"[3,overwrite]", \
        #                outimg//"[3]", tmpbpm//".pl", 1, verbose-)
        #        }
        #    }


        # TODO
        # delete header keywords not needed and fix what need to be fixed

        # then update remaining headers as gem-tlm and gamosaic
        # clean up and go to the next image

        gemhedit (sciphu, "NSCIEXT", 1, "", delete-)

        gemextn (outimg, check="", process="expand", index="1-", extname="", \
            extver="", ikparam="", omit="", replace="", outfile="dev$null", \
            logfile=l_logfile, glogpars="", verbose=no)

        if (gemextn.status != 0) {
            printlog ("ERROR - GAMOSAIC: GEMEXTN returned a non-zero status. \
                Exiting.", l_logfile, verbose+)
            imdelete (outimg, verify-, >& "dev$null")
            goto crash
        }

        gemhedit (sciphu, "NEXTEND", gemextn.count, "", delete-)

        gemdate()
        gemhedit (sciphu, "GEM-TLM", gemdate.outdate, \
            "UT Last modification with GEMINI", delete-)
        gemhedit (sciphu, "GAMOSAIC", gemdate.outdate,\
            "UT Time stamp for GAMOSAIC", delete-)

        imdelete (tmpjnk, verify-, >& "dev$null")

        for (ii=1;ii<=4;ii+=1) {
            imdelete (tmptsci[ii], verify-, >& "dev$null")
            imdelete (tmptvar[ii], verify-, >& "dev$null")
            imdelete (tmptdq[ii], verify-, >& "dev$null")
        }

        delete (tmpbpm//".pl,"//tmpbpm,verify-, >& "dev$null")

        printlog ("\nGAMOSAIC: Created \""//outimg//"\".\n---", \
            l_logfile, l_verbose)

# move on to the next image
NEXT_IMAGE:
    }

    goto clean

#--------------------------------------------------------------------------
crash:
    # Exit with error subroutine
    status = 1
    goto clean

clean:

    if (access(outgaimchk_list)) {
        delete ("@"//outgaimchk_list, verify-, >& "dev$null")
        delete (outgaimchk_list, verify-, >& "dev$null")
    }

    delete (filelist//","//nextlist, verify-, >& "dev$null")

    scanfile = ""

    if (status == 0) {
        printlog ("\nGAMOSAIC -- Exit status:  GOOD", l_logfile, l_verbose)
    } else {
        printlog ("\nGAMOSAIC -- Exit status:  ERROR", l_logfile, l_verbose)
    }

    if (status != 0) {
        printlog ("         -- Please read the logfile \""//l_logfile//\
                "\" for more information.", l_logfile, l_verbose)
    }

    printlog ("------------------------------------------------------------", \
        l_logfile, l_verbose)

end

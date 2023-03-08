# Copyright(c) 2010-2013 Association of Universities for Research in Astronomy, Inc.

procedure gadisplay (inimages, frame)

##M Add gastat call for normalization -- include statextn in parameters

# Takes GSAOI raw or GAPREPAREd data, 4 extensions, full frame and ROI

# Quick display of full frame or ROI raw GSAOI images
# If displaying non gaprepared images, no trimming is applied
# so since the arrays are 2048 x 2048 the final image will
# be 2048*2 + gap

# This task is meant for quick visualization only
# - mosaic the GSAOI images using imtile
# this means the WCS of the mosaic does not take into account any rotation or
# shift between the individual arrays. It simply tile them up as
# [3]  [4]
# [2]  [1]
# with a predefined gap between each quadrant
# current gap values as measured from WCS with trimming (includes the blank
# pixels)
# - gap is defined as integer - no fractional pixels!
# - allows for normalization of each array
# - overlay of static bad pixel mask or DQ planes OR flags non-lin/saturated
# - to combine BPM with non-lin/sat, run gaprepare first.
# - by default, it enters imexamine after display, but that can be turned off
# - however it moves to the next image in the list without query,
#   and each output image is deleted
# - non-linearity and saturation values are for before nlc and very
# conservative

# the input can be in the form of
# - @list containing a list of files one per line, no .fits required, but
# should not complain about it either
# - comma-separated list of images
# - one or more range of image numbers, in the format allowed by gemlist, in
# which case current UT date is assumed for the rootname, or a rootname
# can be given as parameter.

# AFTER COMMISSIONING:
# set correct values for the gaps
# convert from ADU to e- to take into account different gains between arrays?
#   the fl_norm should account for that, as it normalizes each array.

# Version 2010 Jun 10 V0.0 CW - created based on gamosaic V1.0
#         2010 Jul 29 V1.0 CW - working without BPM overlay
#         2011 Mar 09 V1.1 CW - full BPM overlay implemented
#         2011 Oct 04 V1.2 CW - Handling of ROI
#         2012 Feb 23 V1.3 CW - added the option to turn off all dq
#                               overlay with one flag only; updated gap values
#         2012 Mar 19 V1.4 CW - fixed non-linearity and saturation values
#                               if the image has been processed (in e-)

######################################################

char    inimages    {prompt="GSAOI images to display"}
int     frame       {1, prompt="Frame to write to"}
char    rawpath     {"",prompt="Path for input raw images"}
char    rootname    {"",prompt="Root name for images; blank=today UT"}
char    gapfile     {"",prompt="Bad Pixel file with gap information"}
char    badpix      {"gsaoi$data/gsaoibpm_high_full.fits",prompt="Static Bad Pixel Mask - not mosaic"}
int     colgap      {153,prompt="Vertical gap between arrays (in full frame image)"}
int     linegap     {148,prompt="Horizontal gap between arrays (in full frame image)"}
bool    fl_fixpix   {no,prompt="Interpolate across chip gaps?"}
bool    fl_norm     {no,prompt="Normalize each array if not flat fielded?"}
bool    fl_nodqo    {no,prompt="Turn off all DQ/BPM overlay?"}
bool    fl_bpmover  {no,prompt="Overlay static Bad Pixel Mask?"}
bool    fl_satnl    {yes,prompt="Flag non-linear/saturated pixels?"}
bool    fl_usedq    {no,prompt="Use DQ planes instead of static BPM?"}
bool    fl_imexam   {no,prompt="Enter imexamine when displaying image?"}
int     time_gap    {0,min=0,prompt="Time in seconds to leave before displaying images, when fl_imexam=no"}
real    z1          {0.,prompt="Lower limit if not auto scaling"}
real    z2          {0.,prompt= "Upper limit if not auto scaling"}
char    ztrans      {"linear",enum="linear|log|none|user",prompt="Grey level transformation (linear|log|none|user)"}
int     nonlin      {32000,prompt="Non-linearity level to use if keyword not found."}
int     satura      {50000,prompt="Saturation level to use if keyword not found."}
char    gaprep_pref {"g",prompt="Prefix for gaprepared images"}
char    key_datasec {"DATASEC",prompt="Header keyword for data section (size of the data area)"}
char    key_ccdsec  {"CCDSEC",prompt="Header keyword for array section"}
char    key_nlc     {"NONLINEA",prompt="Header keyword for non-linearity level."}
char    key_sat     {"SATURATI",prompt="Header keyword for saturation level."}
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
    char l_inimages, l_rawpath, l_rootname, l_bpmask, l_stdimage
    char l_key_datasec, filelist, l_gapfile
    char utdate, img, inimg, inname, inphu, outimg
    char normsec, scitile, datasec, datakey
    char tmpout, l_dir, l_ztrans, tmpsci, tmpbpm, tmpnorm
    char tmphead, tmpgap, ext1, l_key_nlc, l_key_sat
    char l_key_ccdsec, trimsec[4], tmpdq, dqtile, bpmexpr, tmpfile
    char junk1, junk2, tmpbadpix, l_logfile, l_sci_ext, l_var_ext, l_dq_ext
    char outgaimchk_list, l_gaprep_pref

    int i, ii, l_test, l_colgap, l_linegap, l_frame
    int nctile, nltile, nc1, nc2, nl1, nl2
    int newxgap, newygap, newnx, newny
    int gc1, gc2, gl1, gl2, gc4, gl4, jnk1, jnk2
    int l_nonlin, l_satura, temp, nsci, nfiles, l_sleep

    real meansky, devsky, l_z1, l_z2, skyplus, skyminus
    real avgsky, gainmult

    bool l_fl_fixpix, l_fl_norm, l_verbose, l_fl_usedq
    bool fl_zscale, fl_zrange, l_fl_bpmover, l_fl_satnl
    bool fl_notinh, l_fl_imexam, l_fl_nodqo, prep_check
    bool debug

    struct pr_struct

    # For checking the correct stdimage is set
    char imt_type[3]="imt7","imt8192","imtgmos"
    int n_imt_types = 3

    ####
    # Temporary files
    filelist = mktemp("tmpfile")
    tmpout = mktemp("tmpout")
    tmpsci = mktemp("tmpsci")
    tmpnorm = mktemp("tmpnorm")//".fits"
    tmpgap = mktemp("tmpgap")
    tmpbadpix = mktemp("tmpbadpix")//".fits"
    tmphead = mktemp("tmphead")
    outgaimchk_list = mktemp ("tmpoutgaimchk_list")
    scitile = mktemp ("tmpscitile")

    ####
    # Local variables
    l_inimages = inimages
    l_frame = frame
    l_rawpath = rawpath
    l_rootname = rootname
    l_key_datasec = key_datasec
    l_fl_fixpix = fl_fixpix
    l_gapfile = gapfile
    l_fl_norm = fl_norm
    l_verbose = verbose
    l_colgap = colgap
    l_linegap = linegap
    l_ztrans = ztrans
    l_bpmask = badpix
    l_fl_usedq = fl_usedq
    l_fl_bpmover = fl_bpmover
    l_fl_satnl = fl_satnl
    l_key_nlc = key_nlc
    l_key_sat = key_sat
    l_nonlin = nonlin
    l_satura = satura
    l_key_ccdsec = key_ccdsec
    l_fl_imexam = fl_imexam
    l_fl_nodqo = fl_nodqo
    l_logfile = logfile
    l_sleep = time_gap

    ####
    # Defalut values
    debug = no
    status = 0
    l_test = 0
    nctile = 2
    nltile = 2
    avgsky = 1.0
    gainmult = 1.0

    l_sci_ext = sci_ext
    l_var_ext = var_ext
    l_dq_ext = dq_ext

    fl_notinh = no
    fl_zscale = yes
    fl_zrange = yes
    l_z1 = z1
    l_z2 = z2

    trimsec[1] = ""; trimsec[2] = ""; trimsec[3] = ""; trimsec[4] = ""

    # Used in the case where rootname is set to this looks for gaprepared
    # version - MS
    l_gaprep_pref = gaprep_pref

    if ((l_z1 != 0.) || (l_z2 != 0.) ) {
        fl_zscale = no
        fl_zrange = no
    }

    # Check stdimage is set properly
    if (defvar("stdimage")) {
        show ("stdimage") | scan (l_stdimage)
    } else {
        printlog ("ERROR - GADISPLAY: environmental variable stdimage is not \
            set", l_logfile, verbose+)
        goto crash
    }

    for (i = 1; i <= n_imt_types; i += 1) {
        if (l_stdimage == imt_type[i]) {
            break
        } else if (i == n_imt_types) {
            printlog ("ERROR - GADISPLAY: Please set stdimage to an \
                appropriate value for GSAOI images.\n"//\
                "                   imt7 or imt8192 or imtgmos", \
                l_logfile, verbose+)
            goto crash
        }
    }

    ####
    # Here is where the actual work starts

    # Test the log file - unless specifically set by the user this will send
    # everything to dev$null
    if ((l_logfile == "") || (stridx(" ", l_logfile) > 0)) {
       l_logfile = "dev$null"
    }

    # Print start time
    gemdate (zone="local")
    printlog ("\n------------------------------------------------------------",
        l_logfile, l_verbose)
    printlog ("GADISPLAY -- Started: "//gemdate.outdate//"\n", \
        l_logfile, l_verbose)

    # Record all input parameters
    printlog ("GADISPLAY: Input parameters...", l_logfile, l_verbose)
    printlog ("    inimages     = "//l_inimages, l_logfile, l_verbose)
    printlog ("    rawpath      = "//l_rawpath, l_logfile, l_verbose)
    printlog ("    rootname     = "//l_rootname, l_logfile, l_verbose)
    printlog ("    gapfile      = "//l_gapfile, l_logfile, l_verbose)
    printlog ("    badpix       = "//l_bpmask, l_logfile, l_verbose)
    printlog ("    colgap       = "//l_colgap, l_logfile, l_verbose)
    printlog ("    linegap      = "//l_linegap, l_logfile, l_verbose)
    printlog ("    fl_fixpix    = "//l_fl_fixpix, l_logfile, l_verbose)
    printlog ("    fl_norm      = "//l_fl_norm, l_logfile, l_verbose)
    printlog ("    fl_nodqo     = "//l_fl_nodqo, l_logfile, l_verbose)
    printlog ("    fl_bpmover   = "//l_fl_bpmover, l_logfile, l_verbose)
    printlog ("    fl_satnl     = "//l_fl_satnl, l_logfile, l_verbose)
    printlog ("    fl_usedq     = "//l_fl_usedq, l_logfile, l_verbose)
    printlog ("    fl_imexam    = "//l_fl_imexam, l_logfile, l_verbose)
    printlog ("    time_gap     = "//l_sleep, l_logfile, l_verbose)
    printlog ("    z1           = "//l_z1, l_logfile, l_verbose)
    printlog ("    z2           = "//l_z2, l_logfile, l_verbose)
    printlog ("    ztrans       = "//l_ztrans, l_logfile, l_verbose)
    printlog ("    nonlin       = "//l_nonlin, l_logfile, l_verbose)
    printlog ("    satura       = "//l_satura, l_logfile, l_verbose)
    printlog ("    gaprep_pref  = "//l_gaprep_pref, l_logfile, l_verbose)
    printlog ("    key_datasec  = "//l_key_datasec, l_logfile, l_verbose)
    printlog ("    key_ccdsec   = "//l_key_ccdsec, l_logfile, l_verbose)
    printlog ("    key_nlc      = "//l_key_nlc, l_logfile, l_verbose)
    printlog ("    key_sat      = "//l_key_sat, l_logfile, l_verbose)
    printlog ("    sci_ext      = "//l_sci_ext, l_logfile, l_verbose)
    printlog ("    var_ext      = "//l_var_ext, l_logfile, l_verbose)
    printlog ("    dq_ext       = "//l_dq_ext, l_logfile, l_verbose)
    printlog ("", l_logfile, l_verbose)

    if (l_rootname == l_gaprep_pref) {
       prep_check = yes
    } else {
       prep_check = no
    }

    # Call gaimchk to perform input checks
    gaimchk (inimages=l_inimages, rawpath=l_rawpath, rootname=l_rootname, \
        obstype_allowed="", object_allowed="", \
        key_allowed="", key_forbidden="", key_exists="", \
        fl_prep_check=prep_check, gaprep_pref=l_gaprep_pref, \
        fl_redux_check=no, \
        garedux_pref="r", fl_fail=no, fl_out_check=no, fl_vardq_check=no, \
        sci_ext=l_sci_ext, var_ext=l_var_ext, dq_ext=l_dq_ext, \
        outlist=outgaimchk_list, logfile=l_logfile, verbose=no)

    if (gaimchk.status != 0) {
        printlog ("ERROR - GADISPLAY: GAIMCHK returned a non-zero status. \
            Exiting.", l_logfile, verbose+)
        goto crash
    } else if (!access(outgaimchk_list)) {
        printlog ("ERROR - GADISPLAY: Cannot access output list from \
            GAIMCHK", l_logfile, verbose+)
        goto crash
    }

    # Input files
    tmpfile = ""
    match ("tmpNotFound", outgaimchk_list, stop=no, print_file_name=no, \
        metacharacter=no) | scan (tmpfile)
    if (tmpfile != "") {

        fields (tmpfile, fields=1, lines="1-", quit_if_miss=no, \
            print_file_names=no, >> filelist)

        nfiles = 0
        count (filelist) | scan (nfiles)
        if (nfiles == 0) {
            printlog ("ERROR - GADISPLAY: There are no files to display. \
                Exiting", l_logfile, verbose+)
            goto crash
        }
    } else {
        printlog ("ERROR - GADISPLAY: There are no files to display. \
            Exiting", l_logfile, verbose+)
        goto crash
    }

    # Loop over input images
    scanfile = filelist

    while (fscan(scanfile, img) != EOF) {

        # Make sure scitile doesn't already exist
        delete (scitile, verify-, >& "dev$null")

        fparse (img, verbose-)
        inname = fparse.root//".fits"
        l_dir = fparse.directory
        outimg = tmpout//"_"//inname
        inimg = l_dir//inname
        inphu = inimg//"[0]"

        print ("\n---\nGADISPLAY: Processing File "//inname//".")

        # Find out the number of extensions - most of this code assumes 4!
        keypar (inphu, "NSCIEXT", silent+)
        if (!keypar.found) {
            keypar (inimg//"[0]", "NEXTEND", silent+)
            if (!keypar.found) {
                nsci = 4
            } else {
                nsci = int (keypar.value)
            }
        } else {
            nsci = int (keypar.value)
        }

        # Check if file already gamosaiced - in which case, just display
        # will only overlay BPM if the DQ plane is present.

        keypar (inphu,"GAMOSAIC", silent+)
        if (keypar.found) {

            # Check if displaying overlay has been requested along with either
            # static BPM, saturated/non-linear pixels or use DQ plane - In any
            # case just use pixels from DQ plane.
            if (!l_fl_nodqo && (l_fl_bpmover || l_fl_satnl || l_fl_usedq)) {
               print ("GADISPLAY: Mosaicked image - bad pixels from "//\
                   l_dq_ext//" planes only")

               gemextn (inimg, check="ext=exists", process="expand", \
                   index="", extname=l_dq_ext, extver="1-"//nsci, ikparam="", \
                   omit="", replace="", outfile="dev$null", logfile=l_logfile,\
                   glogpars="", verbose=no)

               if (gemextn.count == 0) {
                   print ("GADISPLAY: No "//l_dq_ext//" plane present - will \
                       not overlay bad pixels")
                   tmpbpm=""
               } else if (gemextn.count != nsci) {
                   print ("GADISPLAY: The number of "//l_dq_ext//" planes \
                       present does \
                       not match the number of science extensions - will not \
                       overlay bad pixels")
                   tmpbpm=""
               } else {
                   tmpbpm = inimg//"["//l_dq_ext//",1]"
               }
            } else {
               tmpbpm=""
            }

            display (inimg//"["//l_sci_ext//",1]", l_frame, bpmask=tmpbpm, \
                bpdisplay="overlay", ocolor="red", erase=yes, \
                fill=no, zscale=fl_zscale, zrange=fl_zrange, z1=l_z1, \
                z2=l_z2, ztrans=l_ztrans)

            if (l_fl_imexam) {
                imexamine (inimg//"["//l_sci_ext//",1]", frame=l_frame)
            } else if (nfiles > 1) {
                sleep (l_sleep)
            }
            goto NEXT_IMAGE
        }

        # this is to allow display of non-gaprepared images

        keypar (inphu, "GAPREPAR", silent+)
        if (!keypar.found) {
            ext1 = "["
        } else {
            ext1 = "["//l_sci_ext//","
        }

        fl_notinh = no

        if (l_fl_satnl) {
            keypar (inimg//ext1//"1]", "NONLINEA", silent+)
            if (!keypar.found) {
                print ("GADISPLAY: non-linearity level keyword not \
                    found, using a value of "//l_nonlin)
                fl_notinh = yes
            }

            keypar (inimg//ext1//"1]", "SATURATI", silent+)
            if (!keypar.found) {
                print ("GADISPLAY: saturation level keyword not found, \
                    using a value of "//l_satura)
                fl_notinh = yes
            }
        }

        ##M Issue with pathlengths! May have to write to a file

        # build the science tiling pattern
        print (inimg//ext1//"2]\n"//\
            inimg//ext1//"1]\n"//\
            inimg//ext1//"3]\n"//\
            inimg//ext1//"4]\n", >> scitile)

        if (debug) {
            printlog ("___scitile(1): ", l_logfile, verbose+)
            type (scitile, >> l_logfile)
            type (scitile)
        }

        # if normalization is requested check if flatfielded
        # this is done checking for a keyword FLATIMG, which GAREDUCE should \
        # add when applying the flat field.
        # note that it multiplies back by the average value to preserve counts
        # when flagging saturated/non-linear

        if (l_fl_norm) {
            keypar (inphu,"FLATIMG", silent+)
            if (keypar.found) {
                print ("GADISPLAY: File already flatfielded. NOT normalizing")
            } else {
                avgsky = 0.
                fxcopy (inimg, tmpnorm, groups="0", new_file=yes, verbose-)
                for (ii = 1; ii <= 4; ii += 1) {
                    keypar (inimg//ext1//ii//"]", "DATASEC", silent+)

                    normsec = keypar.value
                    imstat (inimg//ext1//ii//"]"//normsec, \
                        fields="mean,stddev",lower=INDEF,upper=INDEF, \
                        nclip=0, lsigma=INDEF, usigma=INDEF, binwidth=0.1, \
                        format-,cache-) | scan (meansky,devsky)

                    skyplus = meansky + 3*devsky
                    skyminus = meansky - 3*devsky
                    imstat (inimg//ext1//ii//"]"//normsec, \
                        fields="mean", lower=skyminus, upper=skyplus, \
                        nclip=0, lsigma=INDEF, usigma=INDEF, binwidth=0.1, \
                        format-,cache-) | scan (meansky)

                    avgsky = avgsky+meansky
                    imexpr ("a/b", tmpnorm//ext1//ii//",append]", \
                        inimg//ext1//ii//"]",meansky, verbose-)
                }
                avgsky = avgsky / 4.0
                gemexpr ("a*b", tmpsci, tmpnorm, avgsky, verbose-)

                # Make sure scitile doesn't already exist
                delete (scitile, verify-, >& "dev$null")

                # build the science tiling pattern
                print (tmpsci//ext1//"2]\n"//\
                    tmpsci//ext1//"1]\n"//\
                    tmpsci//ext1//"3]\n"//\
                    tmpsci//ext1//"4]\n", >> scitile)

                if (debug) {
                    printlog ("___scitile(2): ", l_logfile, verbose+)
                    type (scitile, >> l_logfile)
                    type (scitile)
                }

                print ("GADISPLAY: Sky normalization applied.")
            }
        }

        # imtile the science - here is where I handle the ROIs
        # Need to obtain it explicitly so the task can display raw images
        # without having to call gaprepare to create the metaconf

        newnx = INDEF
        newny = INDEF
        newxgap = l_colgap
        newygap = l_linegap

        # if full frame or DET ROI (trimmed) gaps are the original
        # if full frame or DET ROI (untrimmed) gaps are the original-8 (the
        # inner four non-illuminated
        # pixels in each side are part of the measured gap)
        # if ARRAY ROI (never trimmed), ncols/nlines = 2040+size in x or
        # y+original gap, gap = 2040-size+original gap

        keypar (inimg//ext1//"1]","CCDSEC", silent+)
        junk1 = substr(keypar.value,1,3)
        if (junk1 == "[1:" || junk1 == "[5:") {
            # this is either full frame or DET ROI
            keypar (inphu, "TRIMMED", silent+)
            if (keypar.found) {
                # image has been trimmed
                newxgap = l_colgap
                newygap = l_linegap
            } else {
                newxgap = l_colgap - 8
                newygap = l_linegap - 8
            }

            newnx = INDEF
            newny = INDEF

        } else {
            # this is ARRAY ROI, and I need the actual size
            hselect (inimg//ext1//"1]", "naxis[12]", expr="yes") | \
                scan (jnk1,jnk2)
            newnx = 2040 + jnk1 + l_colgap
            newny = 2040 + jnk2 + l_linegap
            newxgap = 2040 - jnk1 + l_colgap
            newygap = 2040 - jnk2 + l_linegap
        }

        fxcopy (inimg, outimg, groups="0", new_file=yes, verbose-)
        imtile ("@"//scitile, outimg//ext1//"1,append]", nctile, nltile, \
            ncoverlap=-newxgap, nloverlap=-newygap, ncols=newnx, nlines=newny,\
            trim_section="", missing_input="", start_tile="ll", \
            row_order=yes, raster_order=no, median_section="", opixtype="r", \
            ovalue=0.0, verbose=no)

        # Bad pixel overlay - use the same trick as above for the ROIs?
        # can use the external BPM, the DQ planes, OR flag the sat/nonlinear

        if (l_fl_usedq) {
            gemextn (inimg, check="ext=exists", process="expand", \
                index="1-",extname=l_dq_ext, extver="", ikparam="", omit="", \
                replace="", outfile="dev$null", logfile=l_logfile, \
                glogpars="", verbose=no)
            if (gemextn.count == 0) {
                print ("GADISPLAY: "//l_dq_ext//" overlay requested, but \
                    no "//l_dq_ext//" plane \
                    present - will not overlay bad pixels")
                tmpbpm=""
            } else {
                if (l_fl_bpmover || l_fl_satnl) {
                    print ("GADISPLAY: "//l_dq_ext//" overlay AND Static \
                        BPM/Saturated \
                        pixels requested - will use only "//l_dq_ext//" plane")
                }
                dqtile = inimg//"["//l_dq_ext//",2],"//\
                    inimg//"["//l_dq_ext//",1],"//\
                    inimg//"["//l_dq_ext//",3],"//\
                    inimg//"["//l_dq_ext//",4]"

                imtile (dqtile, outimg//"["//l_dq_ext//",1,append]", \
                    nctile, nltile, \
                    ncoverlap=-newxgap, nloverlap=-newygap, ncols=newnx, \
                    nlines=newny, trim_section="", missing_input="", \
                    start_tile="ll", row_order=yes, raster_order=no, \
                    median_section="", opixtype="s", ovalue=0.0, verbose=no)

                tmpbpm = outimg//"["//l_dq_ext//",1]"
            }
        } else if (l_fl_bpmover) {
            if (!access(l_bpmask)){
                print ("GADISPLAY: Static bad pixel mask "//l_bpmask//\
                    " not found -  will not overlay bad pixels")
                tmpbpm=""
            } else {
                gemextn (l_bpmask, check="", process="expand", index="1-", \
                    extname="", extver="", ikparam="", omit="", replace="", \
                    outfile="dev$null", logfile=l_logfile, glogpars="", \
                    verbose=no)

                temp = gemextn.count
                if (temp != 4) {
                    print ("GADISPLAY: Static bad pixel mask "//l_bpmask//\
                        "has only "//temp//" extensions (need 4)")
                    print ("           will not overlay bad pixels.")
                    tmpbpm = ""
                } else {
                    if (l_fl_satnl) {
                        print ("GADISPLAY: Static BPM AND non-linear/\
                            saturated pixels requested - will use only \
                            static BPM")
                    }
                    # This needs some work - assumes the static BPM is full
                    # frame, not trimmed (2048x2048)

                    # the copy is to avoid getting a dqtile string that is too
                    # long for imtile to handle

                    copy (l_bpmask, tmpbadpix, verbose-)

                    hselect (l_bpmask//"[1]", "naxis[12]", expr="yes") | \
                        scan (jnk1,jnk2)

                    if (jnk1 != 2048 || jnk2 != 2048) {
                        print ("GADISPLAY: Static BPM is not full frame \
                            untrimmed - will not overlay bad pixels")
                        tmpbpm = ""
                    } else {
                       # if the BPM is 2048, then CCDSEC is the region I need
                       # in the science image, trimmed or not (gaprepare
                       # updates CCDSEC after trimming - the problem is that
                       # each extension is a different one...

                        for (ii = 1; ii<=4; ii+=1) {
                            keypar (inimg//ext1//ii//"]", l_key_ccdsec, \
                                silent+)
                            trimsec[ii] = keypar.value
                        }
                        dqtile = tmpbadpix//"[2]"//trimsec[2]//","//\
                            tmpbadpix//"[1]"//trimsec[1]//","
                        dqtile = dqtile//tmpbadpix//"[3]"//trimsec[3]//","//\
                            tmpbadpix//"[4]"//trimsec[4]

                        imtile (dqtile, outimg//"["//l_dq_ext//",1,append]", \
                            nctile, nltile, ncoverlap=-newxgap, \
                            nloverlap=-newygap, ncols=newnx, nlines=newny, \
                            trim_section="", missing_input="", \
                            start_tile="ll", row_order=yes, raster_order=no, \
                            median_section="", opixtype="s", ovalue=0.0, \
                            verbose-)
                        tmpbpm = outimg//"["//l_dq_ext//",1]"
                    }
                }
            }
        } else if (l_fl_satnl) {
            print ("GADISPLAY: Flagging non-linear/saturated pixels in \
                image "//inimg )

            ### TODO - values for nonlinearity and saturation need to be fixed
            #   if image has been processed
            ## (multiply saturation and nonlinearity limits by the gain)

            bpmexpr = "((a > b) ? 4 : (a > c) ? 2 : 0 )"
            if (fl_notinh) {
                imexpr (bpmexpr, outimg//"["//l_dq_ext//",1,append]", \
                    outimg//ext1//"1]", l_satura, l_nonlin, outtype="ushort", \
                    verbose-)
                tmpbpm = outimg//"["//l_dq_ext//",1]"
            } else {
                tmpdq = mktemp("tmpdq")
                fxcopy (inimg//"[0]", tmpdq, groups="0", new_file=yes, \
                    verbose-)

                for (ii = 1; ii<=4; ii+=1) {

                    # check if the image has been processed (in e-) or not (in
                    # ADU)

                    keypar (inphu, "GAINMULT", silent+)
                    if (keypar.found) {
                        gainmult = real(keypar.value)
                    } else {
                        gainmult = 1.0
                    }

                    keypar (inimg//ext1//ii//"]", l_key_nlc, silent+)
                    l_nonlin = int(keypar.value)*gainmult
                    keypar (inimg//ext1//ii//"]", l_key_sat, silent+)
                    l_satura = int(keypar.value)*gainmult
                    imexpr (bpmexpr, tmpdq//"["//l_dq_ext//","//ii//\
                        ",append]", \
                        inimg//ext1//ii//"]", l_satura, l_nonlin, \
                        outtype="ushort",verbose-)
                }


                dqtile = tmpdq//"["//l_dq_ext//",2],"//tmpdq//"["//l_dq_ext//\
                    ",1],"//tmpdq//\
                    "["//l_dq_ext//",3],"//tmpdq//"["//l_dq_ext//",4]"

                imtile (dqtile, outimg//"["//l_dq_ext//",1,append]", nctile, \
                    nltile, \
                    ncoverlap=-newxgap, nloverlap=-newygap, ncols=newnx, \
                    nlines=newny, trim_section="", missing_input="", \
                    start_tile="ll", row_order=yes, raster_order=no, \
                    median_section="", opixtype="s", ovalue=0.0, verbose-)

                tmpbpm = outimg//"["//l_dq_ext//",1]"
                delete (tmpdq//","//tmpdq//".fits", verify-, >& "dev$null")
            }
        } else {
            tmpbpm=""
        }
        del (tmpbadpix, verify-, >& "dev$null")

         # fix chip gaps -  use badpiximage and fixpix

        if (l_fl_fixpix) {

            # check the gapfile - if empty, use the defaults - gap position is
            # assumed to be for unbinned after gamosaic the gap positions
            # should be in the tmpbpm file

            if ((l_gapfile == "") || (stridx(" ", l_gapfile) > 0)) {
                print ("GADISPLAY: file with array gap position does not \
                    exist.")
                print ("           Calculating from gap parameters and \
                    image size.")
                print ("#GSAOI default gap position, array unbinned after \
                    mosaic", >> tmpbpm)

                # need to read the datasec for extension [2], and the size of
                # the full image, then calculate the gap size/position taking
                # into account the ROIs

                keypar (inimg//ext1//"2]", l_key_datasec, silent+)

                datasec = keypar.value
                jnk1 = stridx(":",datasec)
                jnk2 = stridx(",",datasec)
                nc2 = int(substr(datasec,jnk1+1,jnk2-1))
                jnk1 = strldx(":",datasec)
                jnk2 = strlen(datasec)
                nl2 = int(substr(datasec,jnk1+1,jnk2-1))
                gc1 = int(nc2 + 1)
                gc2 = int(nc2 + newxgap)
                gl1 = int(nl2 + 1)
                gl2 = int(nl2 + newygap)
                imgets (outimg//ext1//"1]","i_naxis1")
                gc4 = int(imgets.value)
                imgets (outimg//ext1//"1]","i_naxis2")
                gl4 = int(imgets.value)
                print (gc1//"    "//gc2//"    1    "//gl4, >> tmpgap)
                print ("1    "//gc4//"    "//gl1//"    "//gl2, >> tmpgap)
                print ("GADISPLAY: using the following for the gaps \
                    positions:")
                printf ("%15d%8d%8d%8d\n", gc1, gc2, 1, gl4) | scan (pr_struct)
                printlog (pr_struct, l_logfile, l_verbose)
                printf ("%15d%8d%8d%8d\n", 1, gc4, gl1, gl2) | scan (pr_struct)
                printlog (pr_struct, l_logfile, l_verbose)

            } else {
                print ("GADISPLAY: using "//l_gapfile//" for the gaps \
                    positions:")
                copy (l_gapfile, tmpgap, verbose-)
            }

            badpiximage (tmpgap, outimg//ext1//"1]", tmpgap//".pl", \
                goodvalue=0, badvalue=1)
            proto.fixpix (outimg//ext1//"1]", tmpgap//".pl", linter="INDEF", \
                cinterp="INDEF",verbose-, pixels-)
            delete (tmpgap//".pl,"//tmpgap, verify-, >& "dev$null")
        }


        # fix headers - need to adjust the WCS, taking one of the extensions as
        # reference note that taking [2] as reference is the easier as it
        # allows to keep CRPIX1 = 1 however, it propagates the coordinate
        # system to the full frame, ignoring possible differences in pixel
        # scale, rotation, etc. between the arrays.

        imhead (inimg//ext1//"2]", long+) | \
            match ("CTYPE","STDIN",stop-, >> tmphead)

        imhead (inimg//ext1//"2]", long+) | \
            match ("CRPIX","STDIN",stop-, >> tmphead)

        imhead (inimg//ext1//"2]", long+) | \
            match ("CRVAL","STDIN",stop-, >> tmphead)

        imhead (inimg//ext1//"2]", long+) | \
            match ("CD","STDIN",stop-, >> tmphead)

        mkheader (outimg//"[0]",tmphead,append+, verbose-)
        mkheader (outimg//ext1//"1]",tmphead,append+, verbose-)
        delete (tmphead, verify-, >& "dev$null")

        # do the display and enter imexamine - since the WCS is in the
        # extension as well as in the PHU ds9/ximtool should be happy.
        # adding the image name as i_title to allow the user to recognize which
        # image it is displaying as the file name is a tmpout<something>

        # if no DQ/BPM overlay is required, just ignore all the bpmask
        # calculated above

        if (l_fl_nodqo) {
            tmpbpm = ""
        }

        gemhedit (outimg//ext1//"1]", "i_title", img, "", delete-)

        display (outimg//ext1//"1]", l_frame, bpmask=tmpbpm, \
            bpdisplay="overlay", ocolor="red", erase=yes, fill=no, \
            zscale=fl_zscale, zrange=fl_zrange, z1=l_z1, z2=l_z2, \
            ztrans=l_ztrans)

        if (l_fl_imexam) {
            imexamine (outimg//ext1//"1]", frame=l_frame)
        } else if (nfiles > 1) {
            sleep (l_sleep)
        }

        imdelete (outimg//","//tmpsci//","//tmpnorm, verify-, >& "dev$null")

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

    delete (filelist//","//scitile, verify-, >& "dev$null")
    scanfile = ""

    gemdate (zone="local")
    printlog ("\n---\nGADISPLAY -- Finished: "//gemdate.outdate, \
        l_logfile, l_verbose)

    if (status == 0) {
        printlog ("\nGADISPLAY -- Exit status:  GOOD", l_logfile, l_verbose)
    } else {
        printlog ("\nGADISPLAY -- Exit status:  ERROR", l_logfile, verbose+)

    }

#    if (status != 0) {
#        printlog ("          -- Please read the logfile \""//l_logfile//\
#                "\" for more information.", l_logfile, l_verbose)
#    }

    printlog ("--------------------------------------------------------\
        ---\n", l_logfile, l_verbose)

end

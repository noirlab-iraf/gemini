# Copyright(c) 2003-2017 Association of Universities for Research in Astronomy, Inc.

procedure gifringe (inimages, outimage)

####
# Create fringe frame
#
# TODO Record the filters that rquire fring correction for the various CCD
# TODO types
#
# TODO More section checking?
#
####

####
# Parameters

char    inimages    {prompt="Input GMOS images"}
char    outimage    {prompt="Output fringe frame"}

# Sky Scaling paramters
char    typezero    {"mean",enum="mean|midpt|skyfile|keyword",prompt="Operation to determine the sky level or zero point"}
char    skysec      {"default", prompt="Zero point statistics section"}
char    skyfile     {"", prompt="File with zero point values for each input image"}
char    key_zero    {"OFFINT", prompt="Keyword for zero level"}

# Mask parameters
real    msigma      {4.,prompt="Sigma threshold above sky for mask"}
char    bpm         {"",prompt="Name of bad pixel mask file or image"}
bool    fl_mask     {no,prompt="Use DQ planes to mask bad pixels during calculation of skylevel"}

# GEMCOMBINE parameters
char    combine     {"median",enum="average|median",prompt="Combination operation"}
char    reject      {"avsigclip",enum="none|minmax|sigclip|avsigclip",prompt="Rejection algorithm"}
char    scale       {"none",enum="none|exposure",prompt="Image scaling"}
char    weight      {"none",enum="none|exposure",prompt="Image Weights"}
char    statsec     {"[*,*]",prompt="Statistics section for image scaling"}
char    expname     {"EXPTIME",prompt="Exposure time header keyword for image scaling"}
int     nlow        {1,min=0,prompt="minmax: Number of low pixels to reject"}
int     nhigh       {1,min=0,prompt="minmax: Number of high pixels to reject"}
int     nkeep       {0,min=0,prompt="Minimum to keep or maximum to reject"}
bool    mclip       {yes,prompt="Use median in sigma clipping algorithms?"}
real    lsigma      {3.,prompt="Lower sigma clipping factor"}
real    hsigma      {3.,prompt="Upper sigma clipping factor"}
real    sigscale    {0.1,prompt="Tolerance for sigma clipping scaling correction"}

# General parameters
char    sci_ext     {"SCI", prompt="Name of science extension"}
char    var_ext     {"VAR", prompt="Name of variance extension"}
char    dq_ext      {"DQ", prompt="Name of data quality extension"}
bool    fl_vardq    {no,prompt="Make variance and data quality planes?"}
char    logfile     {"", prompt="Name of the logfile"}
pset    glogpars    {"", prompt="Logging preferences"}
bool    verbose     {yes, prompt="Verbose output"}
int     status      {0, prompt="Exit status (0=good)"}
struct  *scanfile   {"", prompt="Internal use only"}

####
# Script

begin
    char    l_inimages
    char    l_outimage = ""
    char    l_logfile = ""
    char    l_sci_ext = ""
    char    l_var_ext = ""
    char    l_dq_ext = ""
    char    l_key_zero = ""
    char    l_skyfile = ""
    char    l_skysec = ""
    char    l_typezero = ""
    char    l_combine = ""
    char    l_reject = ""
    char    l_scale = ""
    char    l_weight = ""
    char    l_statsec = ""
    char    l_expname = ""
    char    l_bpm = ""

    bool    l_verbose, l_fl_vardq

    # Gemcombine parameters
    real    l_lsigma, l_hsigma
    real    l_sigscale
    int     l_nlow, l_nhigh, l_nkeep
    bool    l_mclip

    # Object mask parameter
    real    l_msigma

    char    paramstr, errmsg, outputstr, tmpfile, tmpexpand, comblist
    char    tmpimg, tmpfringe, tmpsub, tmpmask, filelist, tmpdel, tmpobjlog
    char    l_sec[3,12]
    char    img, filename, ccdsum, stat_ext, tmp_skysec, mask_file
    char    tmpmerged, diffzero, dqexists, tmpadd, l_out_ccds, l_bpm_file
    char    tmp_mask_file, npix, tmpoutmaskimg, tmpgtoutfile, tmpdqtile
    char    tmpgtile, stat_file, pixels_to_use
    int     imarea, naxis1, naxis2, obj_mask_flag
    int     nimages, maximages, noutimage, nzerovalues
    int     junk, nsci, ii, ccdint, min_x, min_y, max_x, max_y
    real    levelzero, min_area_frac
    bool    l_tile_det, do_tile, do_del, l_stats_only, ishamamatsu, l_fl_mask

    ####
    # Assign input parameters
    l_inimages = inimages
    junk = fscan (outimage, l_outimage)

    junk = fscan (key_zero, l_key_zero)
    junk = fscan (skyfile, l_skyfile)
    junk = fscan (skysec, l_skysec)
    junk = fscan (typezero, l_typezero)
    l_msigma = msigma
    junk = fscan (bpm, l_bpm)
    l_fl_mask = fl_mask

    # gemcombine
    junk = fscan (combine, l_combine)
    junk = fscan (reject, l_reject)
    junk = fscan (scale, l_scale)
    junk = fscan (weight, l_weight)
    junk = fscan (statsec, l_statsec)
    junk = fscan (expname, l_expname)
    l_lsigma = lsigma
    l_hsigma = hsigma
    l_sigscale = sigscale
    l_nhigh = nhigh
    l_nlow = nlow
    l_nkeep = nkeep
    l_mclip = mclip

    junk = fscan (sci_ext, l_sci_ext)
    junk = fscan (var_ext, l_var_ext)
    junk = fscan (dq_ext, l_dq_ext)
    l_fl_vardq = fl_vardq
    junk = fscan (logfile, l_logfile)
    l_verbose = verbose

    # Create temp lists
    comblist = mktemp ("tmpcomblist")
    filelist = mktemp ("tmpfilelist")
    tmpfile = mktemp("tmpfile")
    tmpexpand = mktemp("tmpexpand")

    # Initialize
    status = 1 # Assume task will fail
    nimages = 0
    nzerovalues = 0
    maximages = 100
    diffzero = "no"
    dqexists = "no"
    min_area_frac = 0.5
    obj_mask_flag = 256
    ishamamatsu = no

    cache ("gemextn", "gemdate", "gemcombine", "imcombine", "gtile")

    # Create the list of parameter/value pairs. One pair per line.
    # All lines combined into one string.  Line delimiter is '\n'.
    paramstr =  "inimages       = "//inimages.p_value//"\n"
    paramstr += "outimage       = "//outimage.p_value//"\n"
    paramstr += "typezero       = "//typezero.p_value//"\n"
    paramstr += "key_zero       = "//key_zero.p_value//"\n"
    paramstr += "skyfile        = "//skyfile.p_value//"\n"
    paramstr += "skysec         = "//skysec.p_value//"\n"
    paramstr += "msigma         = "//msigma.p_value//"\n"
    paramstr += "combine        = "//combine.p_value//"\n"
    paramstr += "bpm            = "//bpm.p_value//"\n"
    paramstr += "sci_ext        = "//sci_ext.p_value//"\n"
    paramstr += "var_ext        = "//var_ext.p_value//"\n"
    paramstr += "dq_ext         = "//dq_ext.p_value//"\n"
    paramstr += "fl_vardq       = "//fl_vardq.p_value//"\n"
    paramstr += "logfile        = "//logfile.p_value//"\n"
    paramstr += "glogpars       = "//glogpars.p_value//"\n"
    paramstr += "verbose        = "//verbose.p_value//"\n"

    # Assign a logfile name if not given.  Open logfile and start log.
    # Write parameter/value pairs ("paramstr") to log.
    gloginit (l_logfile, "gifringe", "gmos", paramstr, fl_append+,
        verbose=l_verbose)
    if (gloginit.status != 0) {
        status = gloginit.status
        goto exitnow
    }
    l_logfile = gloginit.logfile.p_value

    # Generate list of input files - must exist and be MEF
    gemextn (l_inimages, check="", process="none", index="", extname="",
        extversion="", ikparams="", omit="extension", replace="",
        outfile=tmpexpand, logfile=l_logfile, verbose=l_verbose)
    gemextn ("@"//tmpexpand, check="exist,mef", process="none", index="",
        extname="", extversion="", ikparams="", omit="", replace="",
        outfile=tmpfile, logfile=l_logfile, verbose=l_verbose)
    nimages = gemextn.count

    if ((gemextn.fail_count > 0) || (nimages == 0) || \
        (nimages > maximages)) {

        if (gemextn.fail_count > 0) {
            errmsg = gemextn.fail_count//" images were not found."
            status = 101
        } else if (nimages == 0) {
            errmsg = "No input images defined."
            status = 121
        } else if (nimages > maximages) {
            errmsg = "Maximum number of input images ["//str(maximages)//"] \
                has been exceeded."
            status = 121
        }
        glogprint (l_logfile, "gifringe", "status", type="error",
            errno=status, str=errmsg, verbose+)
        goto clean
    }
    # Get number of science extensions
    nsci = 0
    head (tmpexpand, nlines=1) | scan(img)
    gemextn (img, check = "", process = "expand", index = "", \
        extname = l_sci_ext, extversion = "1-", ikparams = "", \
        omit = "", replace = "", outfile = "dev$null", \
        logfile="dev$null", glogpars="", verbose=no)
    nsci = gemextn.count

    # Check that output image does not exist
    gimverify (l_outimage)
    if (gimverify.status != 1) {
        status = 102
        glogprint (l_logfile, "gifringe", "status", type="error",
            errno=status, str="Output fringe image \
            ("//l_outimage//") already exists.", verbose+)
        goto clean
    }
    l_outimage = gimverify.outname

    # Check that bpm exists
    if (l_bpm != "" && !imaccess(l_bpm)) {
        status = 101
        glogprint (l_logfile, "gifringe", "status", type="error", \
            errno=status, str="ERROR - GIFRINGE: bad pixel \
            mask "//l_bpm//" does not exist", verbose=yes)
            goto clean
    }

    # Check that skyfile exists and nzerovalues = nimages
    if (typezero == "skyfile") {
        if (l_skyfile == "" || l_skyfile == " ") {
            glogprint (l_logfile, "gifringe", "status", type="error",
                errno=121, str="Must input a skyfile name", verbose+)
            status = 101
            goto clean
        }
        if (access (l_skyfile)) {
            count (l_skyfile) | scan(nzerovalues)
        } else {
            errmsg = "Cannot access skyfile image ("//l_skyfile//")."
            status = 101
            glogprint (l_logfile, "gifringe", "status", type="error",
                errno=status, str=errmsg, verbose+)
            goto clean
        }
        if (nzerovalues != nimages) {
            errmsg = "number of entries in "//l_skyfile//" not equal \
                to number of input images"
            status = 121
        } else if (nzerovalues == 0) {
            errmsg = l_skyfile//" does not contain any zero values"
            status = 121
        }
        if (status > 100) {
            glogprint (l_logfile, "gifringe", "status", type="error",
                errno=status, str=errmsg, verbose+)
            goto clean
        }
    }
    # Find zero offsets for each image
    # Three different ways:
    # 1) User supplied file
    # 2) User supplied keyword
    # 3) Imstat on image section in middle CCD
    scanfile = tmpfile

    while (fscan(scanfile, img) != EOF) {
        tmpimg = mktemp ("tmp"//img)//".fits"
        if (l_fl_mask) {
            copy (img//".fits", tmpimg, verbose=no)
        } else {
            imcopy (img//"[0]", tmpimg, verbose=no, >& "dev$null")
            gemhedit (tmpimg, "NEXTEND", nsci, "", delete-)
            # only copy science extensions over for now.
            for (ii = 1; ii <= nsci; ii+=1) {
                imcopy (img//".fits["//l_sci_ext//","//ii//"]",
                        tmpimg//"[append]", \
                        verbose-)
            }
        }
        print (tmpimg, >> comblist)
        print (tmpimg//"  "//img, >> filelist)
    }

    if (l_typezero == "skyfile") {
        tmpmerged = mktemp("tmpmerge")
        joinlines (filelist, skyfile, output=tmpmerged)
        scanfile = tmpmerged
        while (fscan(scanfile, tmpimg, img, levelzero) != EOF) {
            glogprint (l_logfile, "gifringe", "status", type="string", \
                str="using zero level "//levelzero//" for image "//img)
            # Add each zero value to the headers
            if (l_key_zero == "") {
                   # keyname doesn't matter, it gets deleted later
                   l_key_zero = "OFFINT"
            }
            for (ii = 1; ii <= nsci; ii += 1) {
                gemhedit (tmpimg//"["//l_sci_ext//","//ii//"]", l_key_zero, \
                    levelzero, "", delete-)
            }
        }
        scanfile = ""
        delete (tmpmerged, verify=no, >& "dev$null")

    } else if (l_typezero == "keyword") {
        scanfile = filelist
        while (fscan(scanfile, tmpimg, img) != EOF) {
            for (i = 1; i <= nsci-1; i += 1) {
                keypar (img//"["//l_sci_ext//","//i//"]", l_key_zero, silent+)
                if (keypar.found) {
                    levelzero = real(keypar.value)
                    hdiff (img//"["//l_sci_ext//","//i//"]", \
                        img//"["//l_sci_ext//","//i+1//"]", \
                        keywords=l_key_zero) | scan(diffzero)
                    if (diffzero != "no") {
                        status = 132
                        glogprint (l_logfile, "gifringe", "status", \
                            type="error", errno=status, str= "Keyword \
                            "//l_key_zero//" needs to be identical for \
                            all "//l_sci_ext//" extensions", verbose+)
                        glogprint (l_logfile, "gifringe", "status", \
                            type="error", str="keyword "//l_key_zero//" \
                            differs for "//img//"["//l_sci_ext//","//i//"] \
                            and "//img//"["//l_sci_ext//","//i+1//"]", \
                            verbose+)
                        goto clean
                    }
                } else {
                    status = 131
                    glogprint (l_logfile, "gifringe", "status", type="error",
                        errno=status, str= "Keyword "//l_key_zero//" not \
                        found in "//img//"["//l_sci_ext//","//i//"]", verbose+)
                    goto clean
                }
            }
            glogprint (l_logfile, "gifringe", "status", type="string", \
                str="using zero level "//levelzero//" for image "//img)
        }
        scanfile = ""

    } else {
        # mean or midpoint
        scanfile = filelist
        while (fscan(scanfile, tmpimg, img) != EOF) {
            do_tile = yes
            if (l_skysec == "default") {

                # Set hamamatsu flag
                imgets (img//"[0]", "DETTYPE", >& "dev$null")
                if (imgets.value == "S10892" || imgets.value == "S10892-N") {
                    ishamamatsu = yes
                }

                keypar (img//"["//l_sci_ext//"]", "CCDSUM", silent=yes)
                
                # Set binning
                if (keypar.found) {
                    print(keypar.value) | scan (ccdsum)
                    ccdint = int(ccdsum)

                    min_x = 100
                    min_y = 100
                    max_x = 1800
                    max_y = 4500

                    if (ishamamatsu) {
                        max_y = 4000
                    }
                    pixels_to_use = "[0]["//str(nint(min_x / ccdint))//":"//\
                        str(nint(max_x / ccdint))//","\
                        //str(nint(min_y / ccdint))//\
                        ":"//str(nint(max_y / ccdint))//"]"

                } else {
                    status = 132
                    glogprint (l_logfile, "gifringe", "status", type="error", \
                        errno=status, str="ERROR - GIFRINGE: ccd binning not \
                        found in image "//img, verbose=yes)
                    goto clean
                }
                l_out_ccds = "2"
                l_tile_det = no
            } else {
                # Parse and use user input skysec

                # Assumes extension name is of length 3
                stat_ext = substr(l_skysec,2,4)
                pixels_to_use = substr(l_skysec, \
                    stridx("]", l_skysec)+1, strlen(l_skysec))

                # Allowed values are SCi, DET, CCD
                if (stat_ext != l_sci_ext && stat_ext != "DET" \
                    && stat_ext != "CCD") {
                    status = 122
                    glogprint (l_logfile, "gifringe", "status", type="error", \
                        errno=status, str="ERROR - GIFRINGE: illegal parameter \
                        SKYSEC ("//l_skysec//") format is e.g \
                        [SCI,2][100:1900,100:4500] or \
                        [DET][100:1900,100:4500]", verbose=yes)
                    goto clean

                } else if (stat_ext == "DET") {
                    # Use entire detecor
                    l_out_ccds = "all"
                    l_tile_det = yes

                } else if (stat_ext == "CCD") {
                    # Use one CCD
                    l_out_ccds = substr(l_skysec,6,stridx("]",l_skysec-1))
                    l_tile_det = no
                    # Check that the requested CCD is 1-3.
                    if (int(l_out_ccds) > 3 || int(l_out_ccds) < 1) {
                        status = 122
                        glogptint (l_logfile, "gifringe", "status", \
                            type="error", errno=status, \
                            str="ERROR - GIFRINGE: CCD value can only be \
                                1,2,3", \
                            verbose=yes)
                        got clean
                    }
                } else {
                    # SCI requested no tiling to be done
                    do_tile = no
                }
            }

            if (do_tile) {
                # Calling gtile
                tmpgtoutfile = mktemp("tmpgtileout")//".fits"

                if (l_fl_mask) {
                    l_stats_only = no
                    tmpgtile = mktemp("tmpgtile")
                    tmpdqtile = mktemp("tmpdqtile")
                } else {
                    l_stats_only = yes
                    tmpgtile = tmpgtoutfile
                }

                gtile (inimages=tmpimg, outimages=tmpgtoutfile, \
                    out_ccds=l_out_ccds, \
                    fl_stats_only=l_stats_only, fl_tile_det=l_tile_det,
                    sci_ext=l_sci_ext, \
                    var_ext=l_var_ext, dq_ext=l_dq_ext, mdf_ext="MDF", \
                    key_detsec="DETSEC", key_ccdsec="CCDSEC", \
                    key_datasec="DATASEC", key_biassec="BIASSEC", \
                    key_ccdsum="CCDSUM", fl_verbose=no, logfile=l_logfile)

                if (gtile.status != 0){
                    status = 129
                    glogprint (l_logfile, "gifringe", "status", \
                        type="error", errno=status,
                        str="ERROR - GIFRINGE: GTILE returned "//\
                        "no zero status", verbose=yes)
                    goto clean
                }

                if (l_fl_mask) {
                    # Copy out the appropriate file into simple FITS
                    imcopy (tmpgtoutfile//"["//l_sci_ext//"]", tmpgtile, \
                                verbose=no)
                    imcopy (tmpgtoutfile//"["//l_dq_ext//"]", tmpdqtile, \
                                verbose=no)
                    imdelete (tmpgtoutfile, verify-, >& "dev$null")
                }
            }

            glogprint (l_logfile, "gifringe", "status", type="string", \
                str="Using image section "//l_skysec//" for zero level \
                statistics", verbose=l_verbose)

            # Change mask file according to requested mode
            mask_file = ""
            if (do_tile) {
                stat_file = tmpgtile//pixels_to_use
                if (l_fl_mask) {
                    mask_file = tmpdqtile//pixels_to_use
                }
            } else {
                stat_file = tmpimg//l_skysec
                if (l_fl_mask) {
                    mask_pixels = "["//l_dq_ext//\
                        substr(l_skysec,5,strlen(l_skysec))
                    mask_file = tmpimg//mask_pixels
                }
            }

            print ("GIFRINGE: Calculating statistics for "//tmpimg//l_skysec)
            if (mask_file != "") {
                print ("              Using maskfile: "//mask_file)
            }

            # Calculate statistics
            mimstatistics (stat_file, imasks=mask_file, omasks="", \
                field=l_typezero, lower=INDEF, \
                upper=INDEF, nclip=0, lsigma=INDEF, usigma=INDEF, \
                binwidth=0.1, format=no, cache=no) | scan (levelzero)

            if (do_tile) {
                imdelete (tmpgtile, verify-, >& "dev$null")
                if (l_fl_mask) {
                    imdelete (tmpdqtile, verify-, >& "dev$null")
                }
            }

            glogprint (l_logfile, "gifringe", "status", type="string", \
                str="using zero level "//levelzero//" for image "//img)

            # Add zerolevel to headers
            for (ii = 1; ii <= nsci; ii += 1) {
                gemhedit (tmpimg//"["//l_sci_ext//","//ii//"]", l_key_zero, \
                    -levelzero, "", delete-)
            }
        }
    }
    delete (filelist, verify=no, >& "dev$null")

    # At this point we should have a levelzero
    # make a first rough fringe frame on images with zerokey word in headers
    l_key_zero = "!"//l_key_zero
    tmpfringe = mktemp ("tmpfringe")
    glogprint( l_logfile, "gifringe", "status", type="fork",
        fork="forward", child="gemcombine", verbose=l_verbose )

    gemcombine ("@"//comblist, tmpfringe, title="", logfile=l_logfile, \
        combine=l_combine, reject=l_reject, offsets="none", \
        masktype="none", maskvalue=0., scale=l_scale, \
        zero=l_key_zero, weight=l_weight, statsec=l_statsec, \
        expname=l_expname, lthreshold=INDEF, hthreshold=INDEF, \
        nlow=l_nlow, nhigh=l_nhigh, nkeep=l_nkeep, mclip=l_mclip, \
        lsigma=l_lsigma, hsigma=l_hsigma, key_ron="RDNOISE", \
        key_gain="GAIN", ron=0.0, gain=1., snoise="0.0", \
        sigscale=l_sigscale, pclip=-0.5, sci_ext=l_sci_ext, \
        var_ext=l_var_ext, dq_ext=l_dq_ext, nrejfile="", \
        fl_vardq=l_fl_mask, fl_dqprop=no, verbose=l_verbose)

    if (gemcombine.status != 0) {
        status = 129
        glogprint (l_logfile, "gifringe", "status", \
            type="error", errno=status,
            str="ERROR - GIFRINGE: GEMCOMBINE (first call) returned "//\
            "no zero status", verbose=yes)
        goto clean
    }

    glogprint( l_logfile, "gifringe", "status", type="fork",
        fork="backward", child="gemcombine", verbose=l_verbose)

    # Subtract rough fringe frame off of each image
    scanfile = comblist
    while (fscan(scanfile, tmpimg) != EOF) {
        # Subtract fringe frame off of images that have OFFINT headerkeys
        tmpsub = mktemp ("tmpsub")
        glogprint( l_logfile, "gifringe", "status", type="fork",
            fork="forward", child="gemarith", verbose=l_verbose )

        gemarith (tmpimg, "-", tmpfringe, tmpsub, sci_ext=l_sci_ext, \
            var_ext=l_var_ext, dq_ext=l_dq_ext, mdf_ext="MDF", \
            fl_vardq=yes, verbose=l_verbose, lastout="", logfile=l_logfile)

        glogprint( l_logfile, "gifringe", "status", type="fork",
            fork="backward", child="gemarith", verbose=l_verbose)

        # Create an object mask for each sci ext using subtracted images
        for (ii = 1; ii <= nsci; ii += 1) {
            tmpmask = mktemp ("tmpmask")
            glogprint( l_logfile, "gifringe", "status", type="fork",
                fork="forward", child="objmasks", verbose=l_verbose )

            # Incorporate DQ / BPM as requested
            do_del = no
            if (l_fl_mask && \
                    imaccess(tmpsub//".fits["//l_dq_ext//","//ii//"]")) {
                mask_file = tmpsub//".fits["//l_dq_ext//","//ii//"]"
            } else {
                mask_file = ""
            }

            if (l_bpm != "") {
                l_bpm_file = l_bpm//"["//l_dq_ext//","//ii//"]"
                if (mask_file != "") {
                    tmp_mask_file = mktemp("tmpmaskfile")
                    # TODO expand the parameters
                    imexpr ("a | b", tmp_mask_file, mask_file, l_bpm_file, \
                        verbose-)
                    do_del = yes
                } else {
                    tmp_mask_file = l_bpm_file
                }
            } else if (mask_file != "") {
                tmp_mask_file = mask_file
            } else {
                tmp_mask_file = ""
            }
            
            # Before we run objmasks we need to check if there any good pixels
            if (tmp_mask_file != "") {
                mimstatistics (tmp_mask_file, imasks=tmp_mask_file, \
                    omasks="", fields="npix", format=no) | scan (npix)
                keypar (tmp_mask_file, "i_naxis1", silent+)
                naxis1 = int(keypar.value)
                keypar (tmp_mask_file, "i_naxis2", silent+)
                naxis2 = int(keypar.value)
                imarea = naxis1 * naxis2

                if (real(npix) < (min_area_frac * imarea)) {
                    print ("Skipping "//\
                           tmpsub//".fits["//l_sci_ext//","//ii//"] "//\
                           "not enough good pixels")
                    goto SKIP_OBJMASKS
                }
            }
            tmpobjlog = mktemp("tmpobjmasks")//".log"
            objmasks (tmpsub//".fits["//l_sci_ext//","//ii//"]", \
                tmpmask//".fits", \
                masks=tmp_mask_file, \
                hsigma=l_msigma, hdetect=yes, ldetect=no, \
                omtype="numbers", >& tmpobjlog)

            if (!imaccess(tmpmask//".fits")) {
                status = 129
                glogprint (l_logfile, "gifringe", "status", \
                    type="error", errno=status,
                    str="ERROR - GIFRINGE: No output from OBJMASKS",
                    verbose=yes)
                goto clean
            }
                
            glogprint (l_logfile, "gifringe", "science", type="file",
                str=tmpobjlog)
            glogprint( l_logfile, "gifringe", "status", type="fork",
                fork="backward", child="objmasks", verbose=l_verbose )
            delete (tmpobjlog)

            # Copy object mask to subtracted images as dq plane
            if (imaccess(tmpimg//"["//l_dq_ext//","//ii//"]")) {
                tmpoutmaskimg = mktemp("tmpoutmaskimg")
                imexpr ("((a > 0) ? "//str(obj_mask_flag)//" : 0) | b", \
                    tmpoutmaskimg, \
                    tmpmask//".fits[pl]", \
                    tmpimg//"["//l_dq_ext//","//ii//"]", \
                    verbose=no)
                # An aweful work around due to 2.16 issues with
                # FITS kernel names e.g., overwrite...
                imcopy (tmpoutmaskimg,
                        tmpimg//"["//l_dq_ext//","//ii//"][*,*]", verbose=no)
                imdelete (tmpoutmaskimg, verify-, >& "dev$null")
            } else {
                imcopy (tmpmask//".fits[pl]", \
                    tmpimg//"["//l_dq_ext//","//ii//",append+]", verbose=no)
            }
SKIP_OBJMASKS:
            if (do_del) {
                imdelete (tmp_mask_file, verify-, >& "dev$null")
            }
            delete (tmpmask//".fits", verify=no, >& "dev$null")
        }
        imdelete (tmpsub//".fits", verify=no, >& "dev$null")
    }

    imdelete (tmpfringe, verify=no, >& "dev$null")

    # Use object mask while combining input list to create final fringe frame
    glogprint( l_logfile, "gifringe", "status", type="fork",
        fork="forward", child="gemcombine", verbose=l_verbose )

    gemcombine ("@"//comblist, l_outimage, title="", logfile=l_logfile,
        combine=l_combine, reject=l_reject, offsets="none",
        masktype="goodvalue", maskvalue=0., scale=l_scale,
        zero=l_key_zero, weight=l_weight, statsec=l_statsec, expname=l_expname,
        lthreshold=INDEF, hthreshold=INDEF, nlow=l_nlow,
        nhigh=l_nhigh, nkeep=l_nkeep, mclip=l_mclip, lsigma=l_lsigma,
        hsigma=l_hsigma, key_ron="RDNOISE", key_gain="GAIN", ron=0.0,
        gain=1.0, snoise="0.0", sigscale=l_sigscale, pclip=-0.5,
        sci_ext=l_sci_ext, var_ext=l_var_ext, dq_ext=l_dq_ext, nrejfile="",
        fl_vardq=l_fl_vardq, fl_dqprop=no, verbose=l_verbose)

    if (gemcombine.status != 0) {
        status = 129
        glogprint (l_logfile, "gifringe", "status", \
            type="error", errno=status,
            str="ERROR - GIFRINGE: GEMCOMBINE (second call) returned "//\
            "no zero status", verbose=yes)
        goto clean
    }

    glogprint(l_logfile, "gifringe", "status", type="fork",
        fork="backward", child="gemcombine", verbose=l_verbose)

    # Delete zero science ext header keyword from output fringe frame
    l_key_zero = substr(l_key_zero, 2, strlen(l_key_zero))
    for (ii = 1; ii <= nsci; ii += 1) {
        gemhedit (l_outimage//"["//l_sci_ext//","//ii//"]", l_key_zero, \
            -levelzero, "", delete=yes)
    }

    # Set time stamps
    gemdate ()
    gemhedit (l_outimage//"[0]", "GEM-TLM", gemdate.outdate,
        "UT Last modification with GEMINI", delete-)
    gemhedit (l_outimage//"[0]", "GIFRINGE", gemdate.outdate,
        "UT Time stamp for GIFRINGE", delete-)

    # If we get here all is good
    status = 0

clean:
    scanfile = ""
    delete (tmpfile, verify=no, >& "dev$null")
    delete (filelist, verify=no, >& "dev$null")
    delete (tmpexpand, verify=no, >& "dev$null")

    if (access(comblist)) {
        imdelete ("@"//comblist, verify=no, >& "dev$null")
        delete (comblist, verify=no, >& "dev$null")
    }

    if (status == 0) {
        glogclose (l_logfile, "gifringe", fl_success=yes, verbose=l_verbose)
    } else {
        glogclose (l_logfile, "gifringe", fl_success=no, verbose=yes)
    }

exitnow:
    ;

end



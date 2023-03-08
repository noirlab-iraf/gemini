# Copyright(c) 2000-2012 Association of Universities for Research in Astronomy, Inc.

procedure qsky(inimages, outimage)

# Make an average or median sky image for QUIRC images, flagging
# the objects.
# 
# Version  Apr 23, 2001 JJ,IJ release v1.2
#          Aug 20, 2003 KL  IRAF2.12 - new/modified parameters
#                             hedit: addonly
#                             imcombine: headers,bpmasks,expmasks,outlimits
#                                 rejmask->rejmasks, plfile->nrejmasks
#                             imstat: nclip,lsigma,usigma,cache
#                             apphot.daofind: wcsout,cache

char    inimages  {prompt="Raw QUIRC images to combine"}
char    outimage  {prompt="Output sky image"}
char    outtitle  {"default", prompt="Title for output image"}
char    combtype  {"default", prompt="Type of combine operation",
                    enum="default|median|average"}
char    rejtype   {"avsigclip", prompt="Type of rejection", enum="none|avsigclip|minmax"}
char    logfile   {"", prompt="Name of log file"}
int     nlow      {0, min=0, prompt="minmax: Number of low pixels to reject"}
int     nhigh     {1, min=0, prompt="minmax: Number of high pixels to reject"}
real    lsigma    {3., min=0., prompt="avsigclip: Lower sigma clipping factor"}
real    hsigma    {3., min=0., prompt="avsigclip: Upper sigma clipping factor"}
real    threshold {4.5, min=1.5, prompt="Threshold in sigma for object detection"}
real    fwhmpsf   {4.0, min=0.5, prompt="Estimate of PSF FWHM in pixels"}
real    rebin     {7.0, min=4.0, prompt="Minimum PSF FWHM (pix) for rebinning"}
real    datamax   {50000., prompt="Saturation level for detector"}
char    key_ron   {"RON", prompt="Header keyword for read noise in electrons"}
char    key_gain  {"GAIN", prompt="Header keyword for gain in electrons/ADU"}
real    ron       {15., min=0., prompt="Default read noise if header keyword absent"}
real    gain      {1.85, min=0.01, prompt="Default gain to use if header keyword absent"}
char    key_exptime {"EXPTIME", prompt="Header keyword for exposure time"}
char    key_filter {"FILTER", prompt="Header keyword for filter id"}
char    key_airmass {"AIRMASS", prompt="Header keyword for airmass"}
char    masksuffix {"msk", prompt="Mask name suffix"}
real    maskfactor {1., min=0.01, prompt="Scaling factor to change the size of mask holes"}
bool    fl_keepmasks {no, prompt="Keep object masks for each input image?"}
bool    fl_gemseeing {no, prompt="Use gemseeing to refine the PSF FWHM?"}
bool    verbose   {no, prompt="Verbose actions"}
int     status    {0, prompt="Exit status (0=good)"}
struct* scanfile  {"", prompt="Internal use only"}

begin

    # declare local variables
    char    l_inimages, l_outimage, l_combtype, l_rejtype, l_logfile
    char    l_outtitle
    char    l_key_ron, l_key_gain, l_key_filter, l_key_airmass, l_masksuffix
    char    l_key_exptime
    real    l_fwhmpsf, l_threshold, l_datamax, l_lsigma, l_hsigma
    bool    l_verbose, s_verbose, l_fl_gemseeing, l_fl_keepmasks, l_fl_objects
    int     l_nlow, l_nhigh, l_xsize, l_ysize
    int     i, nimages, l_binfactor, l_nempty, l_nbadpix
    real    skyconst, l_ron, l_gain, l_maskfactor, l_psffactor
    char    img, mask, l_origbpm, l_temp
    char    tmpfile1, tmpfile2, tmpcoo, tmpexp, tmpreg
    char    tmpsky, tmpflat, tmpmsk, tmpimg
    struct  l_struct
    real    l_skylevel, l_skysigma, l_x, l_y, l_r, l_rebin, l_mean, l_sig
    real    l_fwhmmin

    status=0
    scanfile=""
    cache ("imgets", "tinfo", "gemdate")

    # Factor by which to multiply the PSF FWHM prior to sending to daofind.
    # This helps guarantee that faint fuzzy galaxies are identified and masked.
    l_psffactor = 2.

    nimages=0
    skyconst=0.

    # Set local variables to input parameters
    l_inimages=inimages ; l_outimage=outimage
    l_outtitle=outtitle
    l_combtype=combtype ; l_rejtype=rejtype
    l_fwhmpsf=fwhmpsf ; l_threshold=threshold ; l_datamax=datamax
    l_key_ron=key_ron ; l_ron=ron ; l_key_gain=key_gain ; l_gain=gain
    l_key_exptime=key_exptime
    l_key_filter=key_filter ; l_key_airmass=key_airmass
    l_masksuffix=masksuffix ; l_maskfactor=maskfactor 
    l_fl_keepmasks=fl_keepmasks
    l_logfile=logfile ; l_verbose=verbose ; l_fl_gemseeing=fl_gemseeing
    l_nlow=nlow ; l_nhigh=nhigh
    l_lsigma=lsigma ; l_hsigma=hsigma ; l_rebin=rebin

    # Make temporary files
    tmpfile1 = mktemp("tmpfl")
    tmpfile2 = mktemp("tmpfl")
    tmpcoo = mktemp("tmpcoo")
    tmpexp = mktemp("tmpexp")
    tmpreg = mktemp("tmpreg")
    tmpsky = mktemp("tmpsky")
    tmpflat= mktemp("tmpflat")
    tmpmsk = mktemp("tmpmsk")
    tmpimg= mktemp("tmpimg")

    # Check for package log file or user-defined log file
    cache("quirc")
    if ((l_logfile=="") || (l_logfile==" ")) {
        l_logfile=quirc.logfile
        if ((l_logfile=="") || (l_logfile==" ")) {
            l_logfile="quirc.log"
            printlog("WARNING - QSKY:  Both qsky.logfile and quirc.logfile \
                are empty.", logfile=l_logfile, verbose+)
            printlog("                 Using default file quirc.log.",
                logfile=l_logfile, verbose+)
        }
    }

    # Open log file
    date | scan(l_struct)
    printlog("--------------------------------------------------------------\
        --------------", logfile=l_logfile, verbose=l_verbose)
    printlog("QSKY -- "//l_struct, logfile=l_logfile, verbose=l_verbose)
    printlog(" ", logfile=l_logfile, verbose=l_verbose)

    if (l_inimages=="" || l_inimages==" ") {
        printlog("ERROR - QSKY: No input images defined",
            logfile=l_logfile, verbose+)
        status=1
        goto clean
    }
    if (l_outimage=="" || l_outimage==" ") {
        printlog("ERROR - QSKY: No output image defined",
            logfile=l_logfile, verbose+)
        status=1
        goto clean
    }

    # Check mask suffix
    if ((stridx(" ", l_masksuffix) > 0) || \
        (stridx(",", l_masksuffix) > 0) || (l_masksuffix == "")) {
        printlog("WARNING - QSKY: Illigal masksuffix "//l_masksuffix//\
            " changed to msk", logfile=l_logfile, verbose+)
        l_masksuffix="msk"
    }

    # Warn user about PSF FWHM too small
    if (l_fwhmpsf < 1.5) {
        printlog("WARNING - QSKY:  qsky is not designed to work on \
            undersampled", logfile=l_logfile, verbose+)
        printlog("                 data where the PSF FWHM < 1.5 pixels.",
            logfile=l_logfile, verbose+)
    }

    #----------------------------------------------------------------------
    #check to see if output image already exists; strip .fits if present
    if (substr(l_outimage, strlen(l_outimage) - 4, \
        strlen(l_outimage)) == ".fits")
        l_outimage=substr(l_outimage, 1, (strlen(l_outimage)-5))

    if (imaccess(l_outimage//".fits") || \
        imaccess(l_outimage//l_masksuffix//".pl")) {
        printlog("ERROR - QSKY: Output image "//l_outimage//".fits and/or "//\
            l_outimage//l_masksuffix//".pl already exists",
            logfile=l_logfile, verbose+)
        status=1
        goto clean
    }
    #----------------------------------------------------------------------
    # Put all images in a temporary image list file: tmpfile1
    if (substr(l_inimages, 1, 1)=="@") 
        type(substr(l_inimages, 2, strlen(l_inimages)), > tmpfile1)
    else if (stridx("*", l_inimages)>0)
        files(l_inimages, sort-, > tmpfile1)
    else 
        files(l_inimages, sort-, > tmpfile1)

    printlog("Using input files:", logfile=l_logfile, verbose=l_verbose)
    if (l_verbose) type(tmpfile1)
    type(tmpfile1, >> l_logfile)
    printlog("Output image: "//l_outimage, logfile=l_logfile, \
        verbose=l_verbose)

    # Verify that input images actually exist. (Redundant for * though)
    # and remove .fits if present
    scanfile = tmpfile1
    while (fscan(scanfile, img) != EOF) {
        if (imaccess(img)) {
            if (substr(img, strlen(img)-4, strlen(img)) == ".fits")
                img=substr(img, 1, (strlen(img)-5))
            print(img, >> tmpfile2)
        } else {
            printlog("ERROR - QSKY: image "//img//" does not exist.",
                logfile=l_logfile, verbose+)
            status=1
            goto clean
        }
    }
    scanfile=""

    #-------------------------------------------------------------------------
    # Get the number of images, check to make sure there are at least 2

    count(tmpfile2) | scan(nimages)

    if (nimages == 1) {
        printlog("ERROR - QSKY: Cannot combine a single image",
            l_logfile, verbose+)
        delete(tmpfile1//","//tmpfile2, verify-, >>& "dev$null")
        status=1
        goto clean
    } else if (nimages == 0) {
        printlog("ERROR - QSKY: No images to combine", l_logfile, verbose+)
        delete(tmpfile1//","//tmpfile2, verify-, >>& "dev$null")
        status=1
        goto clean
    }

    #--------------------------------------------------------------------------
    # Make a quick temporary sky image using qfastsky.cl and a flat 
    # from the sky
    printlog("QSKY calling QFASTSKY", l_logfile, l_verbose)
    qfastsky("@"//tmpfile2, tmpsky, outtitle="default", combtype="default",
        verbose=l_verbose, key_exptime=l_key_exptime, logfile=l_logfile)
    printlog("Returning to QSKY", l_logfile, l_verbose)

    imstat(tmpsky, fields="midpt,stddev", lower=INDEF, upper=INDEF, nclip=0,
        lsigma=INDEF, usigma=INDEF, binwidth=0.01, format-, cache-) | \
        scan(l_mean, l_sig)
    imstat(tmpsky, fields="mean", lower=(l_mean-4*l_sig), \
        upper=(l_mean+4*l_sig), nclip=0, lsigma=INDEF, usigma=INDEF, \
        binwidth=0.01, format-, cache-) | scan(l_mean)

    if (l_mean <=0.) l_mean=0.0001
    if (sqrt(1./(l_mean*nimages*l_gain))<=0.02) {
        imcalc(tmpsky, tmpflat, "if (im1 > "//(l_mean+4.*l_sig)//") || \
            (im1 < "//(l_mean-4.*l_sig)//") then 1. else (im1/"//l_mean//")",
            verbose-)
    } else
        imcalc(tmpsky, tmpflat, "1.", verbose-)

    #
    # create bad pixel mask in case none is specified in the header
    imcalc(tmpflat, tmpmsk, "if((im1>1.3)||(im1<0.5)) then 1. else 0.", \
        verbose-)
    #
    # set the image size variables
    imgets(tmpflat, "i_naxis1", >>& "dev$null")
    l_xsize = int(imgets.value)
    imgets(tmpflat, "i_naxis2", >>& "dev$null")
    l_ysize = int(imgets.value)
    #--------------------------------------------------------------------------
    # Reduce the images, identify objects, and make a mask for each input file
    #
    # if the object mask file already exists then that file is used and task
    # continues

    cache("tinfo")
    scanfile=tmpfile2
    while (fscan(scanfile, img) != EOF) {

        l_origbpm="none"
        if (!access(img//l_masksuffix//".pl")) {

            # get the input bad pixel mask, if it is specified
            imgets(img, "BPM", >>& "dev$null")
            if ((imgets.value != "0") && (!access(imgets.value)))
                printlog("WARNING - QSKY: Could not find mask " // \
                    imgets.value // " specified in the header", l_logfile, \
                    verbose+)

            if ((imgets.value != "0") && (access(imgets.value))) {
                l_origbpm=imgets.value
                gemhedit(img, "ORIGBPM", imgets.value, "Original BPM", delete-)
            }

            # Reduce sky images using the temporary sky and flat images, fixing
            # bad pixels to help gemseeing and daofind work better
            imstat(img, fields="midpt", lower=INDEF, upper=INDEF, nclip=0,
                lsigma=INDEF, usigma=INDEF, binwidth=0.1, format-, cache-) | \
                scan(skyconst)
            if (l_origbpm != "none")
                imcalc(img//","//tmpsky//","//tmpflat//","//l_origbpm, tmpimg,
                    "if (im4>0) then "//skyconst//" else ((im1-im2)/im3)+"\
                    //skyconst, verbose-)
            else
                imcalc(img//","//tmpsky//","//tmpflat//","//tmpmsk, tmpimg,
                    "if (im4>0) then "//skyconst//" else ((im1-im2)/im3)+"\
                    //skyconst, verbose-)

            # get the read noise and gain from the header for later
            imgets(tmpimg, l_key_ron, >>& "dev$null")

            if (imgets.value == "0")
                printlog("WARNING - QSKY: Adopting a default read noise of "//\
                    l_ron//" electrons for image "//img, l_logfile, verbose+)
            else
                l_ron=real(imgets.value)
            
            imgets(tmpimg, l_key_gain, >>& "dev$null")
            if (imgets.value == "0") {
                printlog("WARNING - QSKY: Adopting a default gain of "//l_gain,
                    l_logfile, verbose+)
                printlog("                electrons/ADU for image"//img,
                    l_logfile, verbose+)
            } else
                l_gain=real(imgets.value)

            # Rebin, if FWHM > (rebin) pix, to speed up the daofind part
            if (l_fwhmpsf > l_rebin) {
                l_binfactor = nint(l_fwhmpsf/3.)
                printlog("Binning by a factor of "//l_binfactor,
                    l_logfile, verbose=l_verbose)
                blkavg(tmpimg, tmpimg, l_binfactor, l_binfactor,
                    option="average")
                l_fwhmpsf = l_fwhmpsf/l_binfactor
                l_ron = l_ron*l_binfactor         # ron=ron*sqrt(nsamples)
                l_gain = l_gain*l_binfactor**2    # gain=gain*nsamples
            } else
                l_binfactor=1

            # Determine the seeing
            if (l_fl_gemseeing) {
                l_fwhmmin=max(1.5, l_fwhmpsf/4.)

                gemseeing(tmpimg, fl_keep=no, fl_update=yes, psffile="default",
                    fl_useall-, maxnobj=100, fl_overwrite+, 
                    sigthres=l_threshold, threshold=INDEF, datamax=l_datamax,
                    ron=l_ron, gain=l_gain, key_ron="", key_gain="",
                    fwhmin=l_fwhmpsf, width=4., mdaofaint=0.,
                    fwhmmin=l_fwhmmin, fwhmmax=0., verbose-, fl_database=no, 
                    fl_strehl-, fl_cleed+, instrument="", camera="",
                    key_inst="INSTRUME", key_camera="CAMERA",
                    key_filter=l_key_filter, sharplow=0.3, sharphigh=1,
                    sharpsig=3, dmagsig=2, fwhmsig=2.5, pixscale=1.,
                    fl_inter=no, logfile="", >>& "dev$null")

                imgets(tmpimg, "FWHMPSF", >>& "dev$null")
                if (real(imgets.value) == 0.) {
                    printlog("WARNING - QSKY:  gemseeing failed to find any \
                        objects.", logfile=l_logfile, verbose+)
                    printlog("                 Staying with previous PSF FWHM \
                        value.", logfile=l_logfile, verbose+)
                } else if (real(imgets.value) > l_fwhmmin) {
                    l_fwhmpsf = real(imgets.value)
                    printf("Adopting new PSF FWHM value from gemseeing: %7.2f \
                        pix for image %s\n", (l_fwhmpsf*l_binfactor), img) | \
                        scan(l_struct)
                    printlog(l_struct, l_logfile, l_verbose)
                } else {
                    printlog("WARNING - QSKY:  gemseeing returned a PSF FWHM \
                        value that", logfile=l_logfile, verbose+)
                    printlog("                 was apparently too small.\
                        Staying with", logfile=l_logfile, verbose+)
                    printlog("                 the previous PSF FWHM value.",
                        logfile=l_logfile, verbose+)
                }
            }

            # Identify stars in reduced sky images using daofind

            # soften the PSF estimate to help find and remove faint fuzzies
            l_fwhmpsf = l_fwhmpsf * l_psffactor

            # Get a reasonable estimate of the noise in the background,
            # just use statistics on the whole image
            imstat(tmpimg, fields="stddev,midpt", lower=INDEF, upper=INDEF,
                nclip=0, lsigma=INDEF, usigma=INDEF, binwidth=0.02, format-,
                cache-) | scan(l_skysigma, l_skylevel)
            imstat(tmpimg, fields="stddev,midpt", \
                lower=(l_skylevel-8*l_skysigma), \
                upper=(l_skylevel+8*l_skysigma), nclip=0, lsigma=INDEF, \
                usigma=INDEF, binwidth=0.02, format-, cache-) | \
                scan(l_skysigma, l_skylevel)

            l_fl_objects=yes
            daofind (tmpimg, output=tmpimg//".coo", starmap="", skymap="",
                datapars="", findpars="", boundary="nearest", constant=0.,
                interactive=no, icommands="", gcommands="", wcsout="logical",
                cache=no, verify=no, update=no, verbose-, graphics="stdgraph",
                display="stdimage",
                datapars.scale=1, datapars.fwhmpsf=l_fwhmpsf,
                datapars.emission=yes, datapars.sigma=l_skysigma,
                datapars.datamin=INDEF, datapars.datamax=INDEF,
                datapars.noise="poisson", datapars.ccdread="", 
                datapars.gain="", datapars.readnoise=l_ron,
                datapars.epadu=l_gain, datapars.airmass=l_key_airmass,
                datapars.filter=l_key_filter, findpars.threshold=l_threshold,
                findpars.nsigma=1.5, findpars.ratio=1., findpars.theta=0.,
                findpars.sharplo=0., findpars.sharphi=100.,
                findpars.roundlo=-100., findpars.roundhi=100.,
                findpars.mkdetections=no)

            # put the PSF FWHM back 
            l_fwhmpsf = l_fwhmpsf / l_psffactor

            # classify these "objects" into objects and cosmic-ray-hits
            # and assign radii
            # CRs    : radius=1        sharpness outside 0.2-1.0 
            #                       OR round outside -1 -> 1.
            # objects: radius ....
            #   Based on an empirical formula developed by IJ and including 
            # factor of 1.5 added by JJ.  Additional scaling applied via 
            # maskfactor.
            #
            # Make the coordinate file an STSDAS table
            pconvert(tmpimg//".coo", tmpcoo//".tab", "*", expr=yes, append-)
            tinfo(tmpcoo, >>& "dev$null")

            if (tinfo.nrows < 1) l_fl_objects=no

            if (l_fl_objects) {
                tcalc(tmpcoo, "peak",
                    "10**(-mag/2.5)*"//str(l_threshold*l_skysigma),
                    colfmt="f6.0", colunit="")
                tcalc(tmpcoo, "peak", "if sharpness/0.7>1 then "//\
                    "max(sharpness/0.7,0.3)*peak else peak")
                tcalc(tmpcoo, "radius", "(1./(2**0.25-1.)) * ((peak/"\
                    //str(l_skysigma)//")**0.25-1.)", colfmt="f6.3", \
                    colunit="")
                tcalc(tmpcoo, "radius",
                    "if radius>0 then sqrt(radius)*"//l_fwhmpsf//"*1.40/2. "//\
                    "else 0.5")
                tcalc(tmpcoo, "radius",
                    "if radius>0 then "//l_maskfactor//"*1.5*radius else "//\
                    "radius")

                # objects with peak above l_datamax and low sharp are mostly
                # saturated stars, make radius for these larger
                tcalc(tmpcoo, "radius",
                    "if sharpness<0.6 && peak>="//l_datamax//" then "//\
                    "radius*1.8 else radius")
                # set radii of CRs
                tcalc(tmpcoo, "radius",
                    "if sharpness<0.2 || sharpness>1 then 1 else radius",
                    colfmt="f6.3", colunit="")
                tcalc(tmpcoo, "radius", "if ground<-1 || ground>1 then 1 "//\
                    "else radius", colfmt="f6.3", colunit="")

                # Get back into the input coordinate frame
                tcalc(tmpcoo, "radius", "radius*"//l_binfactor)
                tcalc(tmpcoo, "xcenter", "xcenter*"//l_binfactor)
                tcalc(tmpcoo, "ycenter", "ycenter*"//l_binfactor)
                l_ron = l_ron/l_binfactor
                l_gain = l_gain/(l_binfactor**2)
                l_fwhmpsf=l_fwhmpsf*l_binfactor

                tprint(tmpcoo, col="xcenter,ycenter,radius", showrow-, \
                    showhdr-, showunits-, > tmpexp)

                print ("circle") | joinlines("STDIN", tmpexp, \
                    output="STDOUT", delim=" ", missing="circle", \
                    maxchars=161, shortest=no, verbose=no, > tmpreg)

            # end if(l_fl_objects)

            } else # make an empty mask
                print(" ", > tmpreg)

            mskregions (tmpreg, masks=img//l_masksuffix//".pl", refimages="", \
                dims=l_xsize//","//l_ysize, depth=0, regnumber="constant", \
                regval=1, exprdb="none", append=no, verbose=l_verbose)

            # Include the original BPM if specified in the header
            if (l_origbpm != "none") {
                imarith(l_origbpm, "max", img//l_masksuffix//".pl",
                    img//l_masksuffix//".pl")
                printlog("Combining the object mask with the badpixel mask as",
                    logfile=l_logfile, verbose=l_verbose)
                printlog("specified in the header: "//l_origbpm,
                    logfile=l_logfile, verbose=l_verbose)
            }

        } else
            printlog("WARNING - QSKY: Using existing mask file "//\
                img//l_masksuffix//".pl", logfile=l_logfile, verbose+)

        # edit header to include mask file as BPM
        gemhedit (img, "BPM", img//l_masksuffix//".pl", "", delete-)
        if (l_fl_keepmasks) {
            gemhedit(img, "QSKYMASK", img//l_masksuffix//".pl",
                "Bad pixel mask from qsky", delete-)
        }

        # get rid of the temporary reduced image and the coordinate files
        delete(tmpexp // "," // tmpreg, verify-, >>& "dev$null")
        delete(tmpcoo//".tab", verify-, >>& "dev$null")
        delete(tmpimg//".coo", verify-, >>& "dev$null")
        delete(tmpimg//"_coo", verify-, >>& "dev$null")
        delete(tmpimg//"_see", verify-, >>& "dev$null")
        delete(tmpimg//"_see.tab", verify-, >>& "dev$null")
        imdelete(tmpimg, verify-, >>& "dev$null")

    } # end while loop
    scanfile=""

    #---------------------------------------------------------------------
    # imcombine raw sky frames with star masks as bad pixel masks

    # save the user's parameters for imcombine and set them to our values
    delete("uparm$imhimcome.par.org", verify-, >>& "dev$null")
    if (access("uparm$imhimcome.par"))
        copy("uparm$imhimcome.par", "uparm$imhimcome.par.org", \
            verbose=l_verbose)

    cache("imcombine")
    imcombine.headers=""
    imcombine.bpmasks=""
    imcombine.rejmask=""
    imcombine.nrejmasks=l_outimage//l_masksuffix//".pl"
    imcombine.expmasks=""
    imcombine.sigmas=""
    imcombine.logfile = "STDOUT" # so that it can be tee'd later
    # imcombine.combine defined later
    # imcombine.reject defined later
    imcombine.project = no
    imcombine.outtype = "real"
    imcombine.outlimits=""
    imcombine.offsets = "none"
    imcombine.masktype = "goodvalue"
    imcombine.maskvalue = 0.
    imcombine.blank = 0.
    imcombine.scale = "none"
    imcombine.zero = "median"
    imcombine.weight = "none"
    imcombine.statsec = ""
    imcombine.expname = ""
    imcombine.lthreshold = INDEF
    imcombine.hthreshold = INDEF
    # imcombine.nlow defined later
    # imcombine.nhigh defined later
    imcombine.nkeep = 1
    # imcombine.lsigma defined later
    # imcombine.hsigma defined later
    imcombine.grow = 0

    if (l_combtype=="default") {
        imcombine.combine = "average"
        if (nimages ==2) {
            imcombine.reject = "minmax"
            imcombine.nlow=0
            imcombine.nhigh=1
            printlog("WARNING - QSKY: Combining two images by taking the \
                minimum.", logfile=l_logfile, verbose+)
        } else {
            imcombine.reject = "avsigclip"
            imcombine.lsigma = 3.
            imcombine.hsigma = 3.
        }
    } else {
        imcombine.combine = l_combtype
        imcombine.reject= l_rejtype
        imcombine.nlow = l_nlow
        imcombine.nhigh = l_nhigh
        imcombine.lsigma = l_lsigma
        imcombine.hsigma = l_hsigma
        if (nimages < 5) {
            printlog("WARNING - QSKY: Combining 4 or fewer images using "//\
                l_combtype, logfile=l_logfile, verbose+)
            if (l_rejtype == "minmax")
                printlog("                with "//l_nlow//" low and "//\
                    l_nhigh//" high pixels rejected.", logfile=l_logfile, \
                    verbose+)
            else if (l_rejtype == "avsigclip")
                printlog("                with "//l_lsigma//"=lower sigma \
                    and "//l_hsigma//"=upper sigma.", logfile=l_logfile, \
                    verbose+)
            else
                printlog("                   with no pixels rejected.",
                    logfile=l_logfile, verbose+)
        } 
        if ((nimages <= (l_nlow+l_nhigh)) && (l_rejtype=="minmax")) {
            printlog("ERROR - QSKY: Cannot reject more pixels than the \
                number of images.", logfile=l_logfile, verbose+)
            status=1
            goto clean
        }

    } # end not-default section

    # Do the combining (suppress output from imcombine)
    imcombine("@"//tmpfile2, l_outimage, >>& "dev$null")

    printlog(" ", logfile=l_logfile, verbose=l_verbose)
    printlog("Combining "//str(nimages)//" images, using "//imcombine.combine,
        logfile=l_logfile, verbose=l_verbose)
    if (imcombine.reject=="minmax")
        printlog("Rejection is minmax, with "//imcombine.nlow//" low and "//\
            imcombine.nhigh//" high values rejected.",
            logfile=l_logfile, verbose=l_verbose)

    if (imcombine.reject=="avsigclip")
        printlog("Rejection is avsigclip, with lsigma="//imcombine.lsigma//\
            " and hsigma="//imcombine.hsigma,
            logfile=l_logfile, verbose=l_verbose)

    # The new BPM should now be used as it will have any empty regions 
    # that were removed by qsky included.  Warn the user.
    # bad pixels have val=n
    imcalc(l_outimage//l_masksuffix//".pl", l_outimage//l_masksuffix//".pl",
        "if (im1=="//nimages//") then 1 else 0", verbose-)
    imstat(l_outimage//l_masksuffix//".pl", fields="npix", lower=0.5, \
        upper=INDEF, nclip=0, lsigma=INDEF, usigma=INDEF, binwidth=0.1,
        format-, cache-) | scan(l_nempty)
    if (l_origbpm != "none") {
        imstat(l_origbpm, fields="npix", lower=0.5, upper=INDEF, nclip=0,
            lsigma=INDEF, usigma=INDEF, binwidth=0.1, format-, cache-) | \
            scan(l_nbadpix)
        l_nempty = l_nempty - l_nbadpix
        if (l_nempty > 0) {
            printlog("WARNING - QSKY: The sky image probably has "//l_nempty,
                logfile=l_logfile, verbose+)
            printlog("                empty pixels.  Please use "//\
                l_outimage//l_masksuffix//".pl", logfile=l_logfile, verbose+)
            printlog("                as the bad pixel mask from now on when \
                using ", logfile=l_logfile, verbose+)
            printlog(l_outimage//".fits", logfile=l_logfile, verbose+)
        }
    } else {
        printlog("There are "//l_nempty//" empty pixels in the output sky \
            image.", logfile=l_logfile, verbose=l_verbose)
        if (l_nempty > 0) {
            printlog("Consult "//l_outimage//l_masksuffix//".pl to see \
                which pixels are empty.", logfile=l_logfile, verbose=l_verbose)
            printlog("If this number is larger than the number of bad pixels \
                in the bad pixel masks,", logfile=l_logfile, verbose=l_verbose)
            printlog("specified in the headers of the input images, or larger \
                than 0 if no input", l_logfile, l_verbose)
            printlog("masks were specified in the headers, then there are \
                regions where", logfile=l_logfile, verbose=l_verbose)
            printlog("all "//nimages//" input images were masked.",
                logfile=l_logfile, verbose=l_verbose)
        }
        printlog(l_outimage//l_masksuffix//".pl should be used as the BPM \
            with this sky frame.", logfile=l_logfile, verbose=l_verbose)
    }

    # Time stamp output, and put critical parameters in the header
    gemdate ()
    gemhedit(l_outimage//".fits", "GEM-TLM", gemdate.outdate,
        "UT Last modification with GEMINI", delete-)
    gemhedit(l_outimage//".fits", "QSKY", gemdate.outdate, \
        "UT Time stamp for qsky", delete-)
    gemhedit(l_outimage//l_masksuffix//".pl", "GEM-TLM", gemdate.outdate,
        "UT Last modification with GEMINI", delete-)
    gemhedit(l_outimage//l_masksuffix//".pl", "QSKY", gemdate.outdate,
            "UT Time stamp for qsky", delete-)

    gemhedit(l_outimage//".fits", "QSKYCOMB", imcombine.combine,
        "Type of combine operation, qsky", delete-)
    gemhedit(l_outimage//".fits", "QSKYREJE", imcombine.reject,
        "Type of rejection, qsky", delete-)
    if (imcombine.reject=="minmax") {
        gemhedit(l_outimage//".fits", "QSKYNLOW", imcombine.nlow,
            "Number of low pixels to reject, qsky", delete-)
        gemhedit(l_outimage//".fits", "QSKYNHIG", imcombine.nhigh,
            "Number of high pixels to reject, qsky", delete-)
    }
    if (imcombine.reject=="avsigclip") {
        gemhedit(l_outimage//".fits", "QSKYLSIG", imcombine.lsigma,
            "Lower sigma clipping factor, qsky", delete-)
        gemhedit(l_outimage//".fits", "QSKYHSIG", imcombine.hsigma,
            "Upper sigma clipping factor, qsky", delete-)
    }

    # Reset BPM keyword in headers of raw images, put names of raw images
    # in header of output image
    i=1
    scanfile=tmpfile2
    while (fscan(scanfile, img) != EOF) {
        imgets(img, "ORIGBPM", >>& "dev$null")
        if (imgets.value != "0") {
            gemhedit (img, "BPM", imgets.value, "", delete-)
            gemhedit (img, "ORIGBPM", "", "", delete=yes, >>& "dev$null")
        } else 
            gemhedit (img, "BPM", "", "", delete=yes, >>& "dev$null")

        gemhedit(l_outimage//".fits", "QSKYIM"//str(i), img,
            "Input image for qsky", delete-)
        i+=1
    }
    scanfile=""

    # update BPM related keywords in output image
    gemhedit(l_outimage//".fits", "QSKYMASK", l_outimage//l_masksuffix//".pl",
        "Bad pixel mask from qsky", delete-)
    gemhedit (l_outimage//".fits", "ORIGBPM", "", "", delete=yes)
    gemhedit (l_outimage//".fits", "BPM", "", "", delete=yes)
    gemhedit (l_outimage//l_masksuffix//".pl", "BPM", "", "", delete=yes)

    # fix the title
    if (l_outtitle=="default" || l_outtitle=="" || l_outtitle==" ") {
        gemhedit (l_outimage//".fits", "i_title", 
            "SKY IMAGE from gemini.quirc.qsky", "", delete-)
        gemhedit (l_outimage//l_masksuffix//".pl", "i_title",
            "BPM for SKY IMAGE from gemini.quirc.qsky", "", delete-)
    } else {
        gemhedit (l_outimage//".fits", "i_title", l_outtitle, "", delete-)
        gemhedit (l_outimage//l_masksuffix//".pl", "i_title", 
            "BPM for "//l_outtitle, "", delete-)
    }

    #------------------------------------------------------------------------

    # Clean up

clean:
    if (status==0)
        printlog("QSKY exit status:  good", logfile=l_logfile, \
            verbose=l_verbose)

    printlog("--------------------------------------------------------------\
        --------------", logfile=l_logfile, verbose=l_verbose)

    if (access(tmpfile2)) {
        scanfile=tmpfile2
        while (fscan(scanfile, img) != EOF) {
            if (!l_fl_keepmasks) {
            imdelete(img//l_masksuffix//".pl", verify-, >>& "dev$null")
            gemhedit (img, "QSKYMASK", "", "", delete=yes, >>& "dev$null")
            }
        }
        scanfile=""
    }

    delete(tmpfile1//","//tmpfile2, verify-, >>& "dev$null")
    imdelete(tmpsky//","//tmpflat//","//tmpmsk, verify-, >>& "dev$null")

    # return to default parameters for imcombine
    unlearn("imcombine")
    # restore the user's parameters for imcombine
    if (access("uparm$imhimcome.par.org"))
        rename("uparm$imhimcome.par.org", "uparm$imhimcome.par", field="all")

end

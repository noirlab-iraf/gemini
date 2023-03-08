# Copyright(c) 2002-2012 Association of Universities for Research in Astronomy, Inc.

procedure gscrrej (inimage, outimage)

# Replace cosmic rays in raw GMOS data by fitted values.
#
# Date     Sept 20, 2002 JT v1.4 release

string  inimage     {prompt="Input image"}
string  outimage    {prompt="Output image with cosmic rays removed"}
real    datares     {4.0,min=1.0, prompt="Instrumental FWHM in x-direction"}
real    fnsigma     {8.0,min=1.0, prompt="Sigma clipping threshold for fitting"}
int     niter       {5,min=1, prompt="Number of fitting iterations"}
real    tnsigma     {10.0,min=1.0, prompt="Sigma clipping threshold for mask"}
bool    fl_inter    {no,prompt="Examine spline fit interactively?"}
string  logfile     {"",prompt="Logfile"}
bool    verbose     {yes,prompt="Verbose?"}
int     status      {0,prompt="Exit status (0=good)"}

begin

    string  l_inimage, l_outimage, l_logfile
    bool    l_fl_inter, l_verbose
    real    l_fnsigma, l_datares, l_tnsigma
    int     l_niter

    int     n, n_ext, nknots, xpix, ypix, bboxsize
    string  timage, timage2, timage3, timage4, tmedimage, tmedimage2
    string  tnimage, tbpm, tbpm2, tout
    string  inext, outext
    char    keyfound
    bool    ismef
    struct  tstruct

    # Keep parameters from changing by outside world:
    cache ("imgets", "gimverify", "gemlogname", "gemdate")

    # Read the parameters:
    l_inimage = inimage
    l_outimage = outimage
    l_datares = datares
    l_fnsigma = fnsigma
    l_niter = niter
    l_tnsigma = tnsigma
    l_fl_inter= fl_inter
    l_verbose = verbose
    l_logfile = logfile

    # Default exit status:
    status = 1

    # Request temp image names:
    tout = mktemp("tmpout")//".fits"

    # Determine the logfile name:
    gemlogname (logpar=l_logfile, package="gmos")
    l_logfile = gemlogname.logname

    # Start log entry:
    date | scan(tstruct)
    printlog ("------------------------------------------------------------\
        ----------------", l_logfile, l_verbose)
    printlog ("GSCRREJ -- "//tstruct//"\n", l_logfile, l_verbose)

    # Warn if logfile name had to revert to default:
    if (gemlogname.status == 1) {
        printlog ("WARNING - GSCRREJ: both gscrrej.logfile and gmos.logfile \
            are empty;", l_logfile, verbose+)
        printlog ("                   using "//l_logfile//"\n", l_logfile, \
            verbose+)
    } else if (gemlogname.status == 2)
        printlog ("WARNING - GSCRREJ: bad logfile name, "//\
            gemlogname.logpar//" - using "//l_logfile//"\n", l_logfile, \
                verbose+)

    # Log the parameter set:
    printlog ("Input GMOS image                       inimage  = "//l_inimage,
        l_logfile, l_verbose)
    printlog ("Output image with cosmic rays removed  outimage = "//l_outimage,
        l_logfile, l_verbose)
    printlog ("Instrumental FWHM in x-direction       datares  = "//l_datares,
        l_logfile, l_verbose)
    printlog ("Sigma clipping threshold for fitting   fnsigma  = "//l_fnsigma,
        l_logfile, l_verbose)
    printlog ("Number of fitting iterations           niter    = "//l_niter,
        l_logfile, l_verbose)
    printlog ("Sigma clipping threshold for mask      tnsigma  = " \
        //l_tnsigma//"\n", l_logfile, l_verbose)

    # Check that input file exists with .fits extension and whether it is MEF:
    gimverify (l_inimage)
    l_inimage = gimverify.outname//".fits"
    gimverify (l_inimage) # (in case the extension is not .fits)

    # Set flag and perform error checking according to whether input is MEF:
    if (gimverify.status == 1) {
        printlog ("ERROR - GSCRREJ: "//l_inimage//" does not exist",
            l_logfile, verbose+)
        goto error
    } else if (gimverify.status == 0) {
        ismef = yes
        keyfound = ""
        hselect (l_inimage//"[0]", "*PREPAR*", yes) | scan (keyfound)
        if (keyfound == "") {
            printlog ("ERROR - GSCRREJ: need to run gprepare on "//l_inimage,
                l_logfile, verbose+)
            goto error
        }
        imgets (l_inimage//"[0]","NSCIEXT")
        n_ext = int(imgets.value)
        if (n_ext < 1) {
            printlog ("ERROR - GSCRREJ: bad NSCIEXT value in " \
                //l_inimage//"[0]", l_logfile, verbose+)
            goto error
        }
    } else {
        ismef = no
        n_ext = 1
    }

    # Check that the output filename is valid:
    if (l_outimage == "" || (stridx(" ", l_outimage) > 0)) {
        printlog ("ERROR - GSCRREJ: invalid output filename", l_logfile, \
            verbose+)
        goto error
    }

    # Make sure output file has a .fits extension:
    gimverify (l_outimage)
    l_outimage = gimverify.outname//".fits"

    # Check that output file does not exist, before wasting time processing:
    if (gimverify.status != 1) {
        printlog ("ERROR - GSCRREJ: "//l_outimage//" already exists",
            l_logfile, verbose+)
        goto error
    }

    # Warn user if task run previously on these data (since may waste time):
    imgets (l_inimage//"[0]","GSCRREJ", >>& "dev$null")
    if (imgets.value != "0")
        printlog ("WARNING - GSCRREJ: "//l_inimage//" already cleaned\n",
            l_logfile, verbose+)

    # Copy input image to output file, so as to keep all non-sci extensions:
    # (use temporary filename, so no output in case of error/crash)
    if (ismef) fxcopy (l_inimage, tout, groups="", new_file+, verbose-)

    # Determine & log size of box for background estimates, based on FWHM:
    bboxsize = nint(max(2.5*l_datares,7))
    printlog ("Using ~ "//bboxsize//"x"//bboxsize//" pixel box for noise \
        estimates\n", l_logfile, l_verbose)

    # Loop over the extensions:
    for (n=1; n <= n_ext; n=n+1) {
        # Create tmp FITS file name used only within this loop
        timage = mktemp("tmpimage")
        timage2 = mktemp("tmpimage2")
        timage3 = mktemp("tmpimage3")
        timage4 = mktemp("tmpimage4")
        tmedimage = mktemp("tmpmedimage")
        tmedimage2 = mktemp("tmpmedimage2")
        tnimage = mktemp("tmpnimage")
        tbpm = mktemp("tmpbpm")
        tbpm2 = mktemp("tmpbpm2")

        if (ismef) {
            inext = l_inimage//"[SCI,"//n//"]"
            outext = tout//"[SCI,"//n//",overwrite]"
        } else {
            inext = l_inimage
            outext = tout
        }

        #printlog("\nProcessing "//inext//":", l_logfile, l_verbose)
        #printlog("------", l_logfile, l_verbose)

        # Determine image size and set number of knots to match resolution:
        imhead (inext,longheader-) | scanf(inext//"[%d,%d][", xpix, ypix)
        nknots = nint(2 * real(xpix) / l_datares)

        printlog ("Processing "//inext//"; cubic spline order "//nknots,
            l_logfile, l_verbose)

        #printlog("Fitting cubic spline, order "//nknots//", "//l_niter//
        #     " iterations, "//l_fnsigma//" sigma rej.", l_logfile, l_verbose)

        if (l_fl_inter)
            printlog ("(interactive mode: fit parameters subject to change)",
                l_logfile, no)

        # Fit the spectra at the correct bandwidth:
        fit1d (inext, timage, type="fit", axis=1, interactive=l_fl_inter,
            sample="*", naverage=1, function="spline3", order=nknots, 
            low_reject=l_fnsigma, high_reject=l_fnsigma, niterate=l_niter,
            grow=0,graphics="stdgraph", cursor="") 

        # Generate an image of the fitting residuals, for making a BPM:
        imarith (inext, "-", timage, timage2, title="", divzero=0, hparams="",
            pixtype="", calctype="", verbose-, noact-)
        imdelete (timage, go_ahead+, verify-, >& "dev$null")

        #printlog("Thresholding residuals at +/- "//l_tnsigma//" sigma... {",
        #    l_logfile, l_verbose)

        # Generate a mask by thresholding the residual image:
        gscrmask (timage2, tbpm, nsigma=l_tnsigma, boxsize=bboxsize,
            varimage=tnimage, verbose=l_verbose, logfile=l_logfile)
        imdelete (timage2, go_ahead+, verify-, >& "dev$null")

        #printlog("}", l_logfile, l_verbose)
        #printlog("Replacing bad pixels with ~"//bboxsize//"x"//bboxsize//
        #           " median", l_logfile, l_verbose)

        # Replace the CR's with the local median, generating a clean image 
        # which can be used to derive a better cubic spline fit:
        imcopy (inext, tmedimage, verbose-)
        median (tmedimage, tmedimage, bboxsize, bboxsize, zloreject=INDEF,
            zhireject=INDEF, boundary="constant", constant=0.0,verbose-)
        imexpr ("c < 0.0001 ? a : b", tmedimage2, inext, tmedimage, tbpm,
            dims="auto",intype="auto", outtype="auto", refim="auto", bwidth=0,
            btype="constant", bpixval=0.0, rangecheck+, verbose-,exprdb="none")
        imdelete (tmedimage,go_ahead+,verify-, >& "dev$null")

        #printlog("Re-fitting cubic spline, order "//nknots//", "//l_niter//
        #    " iterations, "//l_fnsigma//" sigma rej.", l_logfile, l_verbose)
        #if (l_fl_inter)
        #    printlog("(interactive mode: fit parameters subject to change)",
        #        l_logfile, no)

        # Fit the spectra again at the correct bandwidth:
        fit1d (tmedimage2, timage3, type="fit", axis=1, interactive=l_fl_inter,
            sample="*", naverage=1, function="spline3", order=nknots, 
            low_reject=l_fnsigma, high_reject=l_fnsigma, niterate=l_niter,
            grow=0,graphics="stdgraph", cursor="") 
        imdelete (tmedimage2,go_ahead+,verify-, >& "dev$null") #??

        # Generate an image of the data-fit residuals, for making a new BPM:
        imarith (inext, "-", timage3, timage4, title="", divzero=0, hparams="",
            pixtype="", calctype="", verbose-, noact-)

        #printlog("Thresholding new residuals at +/- "//l_tnsigma//" \
        #    sigma... {", l_logfile, l_verbose)

        # Generate a new mask by thresholding the residual image:
        # (the 2nd iteration handles multi-pixel hits significantly better)
        imdelete (tbpm, go_ahead+, verify-, >& "dev$null")
        gscrmask (timage4, tbpm2, nsigma=l_tnsigma, boxsize=bboxsize,
            varimage=tnimage, verbose=l_verbose, logfile=l_logfile)
        imdelete (timage4, go_ahead+, verify-, >& "dev$null")

        #printlog("}", l_logfile, l_verbose)
        #printlog("Replacing bad pixels with cubic spline values:",
        #    l_logfile, l_verbose)

        # Replace the CRs with the fitted cubic spline values:
        imexpr ("c < 0.0001 ? a : b", outext, inext, timage3, tbpm2,
            dims="auto", intype="auto", outtype="auto", refim="auto", bwidth=0,
            btype="constant", bpixval=0.0, rangecheck+, verbose-, 
            exprdb="none")
        # NB. the bad pixel mask should be saved in the DQ extension in future

        # Clear up temp files for this extension:
        imdelete (timage3, go_ahead+, verify-, >& "dev$null")
        imdelete (tbpm2, go_ahead+, verify-, >& "dev$null")
        imdelete (tnimage, go_ahead+, verify-, >& "dev$null")

    } # end for n <= n_ext

    # Update the GMOS header keywords:
    gemdate ()
    gemhedit (tout//"[0]", "GSCRREJ", gemdate.outdate,
        "UT Time stamp for GSCRREJ", delete-)
    gemhedit (tout//"[0]", "GEM-TLM", gemdate.outdate, 
        "UT Last modification with GEMINI", delete-)
    gemhedit (tout//"[0]", "GSCRRES", l_datares,
        "GSCRREJ equivalent instrumental FWHM", delete-)
    gemhedit(tout//"[0]", "GSCRTHRE", l_tnsigma,
        "GSCRREJ sigma rejection threshold", delete-)

    # Rename output file on successful completion:
    imrename (tout, l_outimage, verbose-)

    # Set exit status good:
    status = 0

error:
    # Jump here directly in the event of an error:

    if (status == 0)
        printlog ("\nGSCRREJ exit status: good", l_logfile, l_verbose)
    else
        printlog ("\nGSCRREJ exit status: error", l_logfile, l_verbose)
    printlog ("------------------------------------------------------------\
        ----------------", l_logfile, l_verbose)

    # Clear up
    if (n_ext >= 1) {
        if (imaccess(timage))   imdelete (timage, go_ahead+, verify-)
        if (imaccess(timage2))  imdelete (timage2, go_ahead+, verify-)
        if (imaccess(timage3))  imdelete (timage3, go_ahead+, verify-)
        if (imaccess(timage4))  imdelete (timage4, go_ahead+, verify-)
        if (imaccess(tmedimage))    imdelete (tmedimage, go_ahead+, verify-)
        if (imaccess(tmedimage2))   imdelete (tmedimage2, go_ahead+, verify-)
        if (imaccess(tnimage))  imdelete (tnimage, go_ahead+, verify-)
        if (imaccess(tbpm))     imdelete (tbpm, go_ahead+, verify-)
        if (imaccess(tbpm2))    imdelete (tbpm2, go_ahead+, verify-)
    }

end

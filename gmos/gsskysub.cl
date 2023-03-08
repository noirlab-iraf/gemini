# Copyright(c) 2002-2017 Association of Universities for Research in Astronomy, Inc.

procedure gsskysub (input, fl_answer)

# Skysubtract (rectified) GMOS LONGSLIT or MOS 2D-spectra
#
# Version  Feb 28, 2002  RC,IJ  v1.3 release
#          Sept 20, 2002    v1.4 release
#          Mar 20, 2003  BM pixel scales for both instruments
#          Jun 3, 2003   IJ GMOS=GMOS-N support for old data
#          Aug 26, 2003  KL IRAF2.12 - new parameter, addonly, in hedit

string  input       {prompt="Input GMOS spectra"}
string  output      {"",prompt="Output spectra"}
string  outpref     {"s",prompt="Output prefix"}
string  sci_ext     {"SCI",prompt="Name of science extension"}
string  var_ext     {"VAR",prompt="Name of variance extension"}
string  dq_ext      {"DQ",prompt="Name of data quality extension"}
bool    fl_vardq    {no,prompt="Propagate VAR/DQ planes"}
bool    fl_oversize {yes,prompt="Use 1.05x slit length to accommodate distortion?"}

# Parameters for fit1d
string  long_sample {"*",prompt="Sky sample for LONGSLIT"}
real    mos_sample  {0.9,min=0.01,max=1,prompt="MOS: Maximum fraction of slit length to use as sky sample"}
real    mosobjsize  {1.,min=0,prompt="MOS: Size of object aperture in arcsec"}
int     naverage    {1,prompt="Number of points in sample averaging"}
string  function    {"chebyshev",prompt="Function to fit",enum="spline3|legendre|chebyshev|spline1"}
int     order       {1,min=0,prompt="Order for fit"}
real    low_reject  {2.5,min=0,prompt="Low rejection in sigma of fit"}
real    high_reject {2.5,min=0,prompt="High rejection in sigma of fit"}
int     niterate    {2,min=0,prompt="Number of rejection iterations"}
real    grow        {0.,min=0,prompt="Rejection growing radius in pixels"}
bool    fl_inter    {no,prompt="Fit interactively"}
bool    fl_answer   {prompt="Continue with interactive fitting"}
string  logfile     {"",prompt="Logfile name"}
bool    verbose     {yes,prompt="Verbose?"}
int     status      {0,prompt="Exit status (0=good)"}
struct  *scanfile   {"",prompt="For internal use only"}

begin

    # local variable definitions

    string  l_input, l_output, l_outprefix, l_sci_ext
    string  l_var_ext, l_dq_ext, l_func, l_logfile, l_lsample
    int     l_navrg, l_order, l_niter
    real    l_msample, l_low_rej, l_high_rej, l_grow, l_mosobjsize
    bool    l_fl_vardq, l_fl_inter, l_verbose, l_fl_answer, l_fl_oversize

    # other variable definitions 
    file    tmpin,tmpout,joinlst,tmpfit,tmpsci,tmptab,tmpvar,tmpdq,mdffile
    string  inlst, outlst, inimg, outimg, suf, slittype, l_bsample
    string  obsmode[200]
    int     nerror, ninpimg, noutimg, j, i, k
    int     numext, msktype, nxpix, nypix, inst, iccd
    int     mdfrow, y1, y2, y3, y4, nextnd, xbin, ybin
    int     mdfpos[200], nsciext[200]
    real    asecmm, xscale, yscale, xcen, ycen, sl_size, s_pos, slenfact
    real    rmsval, tmpval, spos_my1, spos_my3, yccd1, yccd3, ll_mosobjsize
    real    lsample, s_cen, slitwidth, pixscale[2,3]
    bool    pref, mdf, ll_fl_inter
    struct  sdate

    # Set local variables

    l_input=input; l_output=output; l_outprefix=outpref
    l_sci_ext=sci_ext; l_dq_ext=dq_ext; l_var_ext=var_ext
    l_fl_vardq=fl_vardq; l_fl_oversize=fl_oversize; l_lsample=long_sample;
    l_msample=mos_sample; l_mosobjsize=mosobjsize
    l_navrg=naverage; l_func=function
    l_order=order; l_low_rej=low_reject
    l_high_rej=high_reject; l_niter=niterate; l_grow=grow
    l_fl_inter=fl_inter; l_logfile=logfile
    l_verbose=verbose

    # Temporary files and some variable definitions

    tmpin = mktemp("tmpin")
    tmpout = mktemp("tmpout")
    joinlst = mktemp("tmpjoin")
    tmptab = mktemp("tmptable")

    status = 0
    asecmm = 1.611444
    # Pixel scales [inst,iccd]
    # inst: 1 - GMOS-N, 2 - GMOS-S
    # iccd: 1 - EEV2 CCDs, 2 - e2vDD CCDs, 3 - Hamamatsu
    pixscale[1,1] = 0.0727
    pixscale[1,2] = 0.07288
    pixscale[1,3] = 0.0807
    pixscale[2,1] = 0.073
    pixscale[2,3] = 0.0800 ##M PIXEL_SCALE

    cache ("imgets", "tabpar", "gimverify")
    cache ("tprint", "gextverify", "gemdate")

    #
    # First of all: check logfile
    #

    if (l_logfile == "STDOUT") {
        l_logfile = ""
    } else if (l_logfile == "" || stridx(" ",l_logfile)>0) {
        l_logfile = gmos.logfile
        if (l_logfile == "" || stridx(" ",l_logfile)>0) {
            l_logfile = "gmos.log"
            printlog ("WARNING - GSSKYSUB: both gsskysub.logfile and \
                gmos.logfile are empty.", l_logfile, l_verbose)
            printlog ("                    Using default file gmos.log.",
                l_logfile, l_verbose)
        }
    }

    # Set up date

    date | scan(sdate)

    # Logfile : what will be done
    printlog ("------------------------------------------------------------\
        ----------------", l_logfile, l_verbose)
    printlog ("GSSKYSUB -- "//sdate, l_logfile, l_verbose)
    printlog ("",l_logfile,l_verbose)
    printlog ("input        = "//l_input, l_logfile, l_verbose)
    printlog ("output       = "//l_output, l_logfile, l_verbose)
    printlog ("outpref      = "//l_outprefix, l_logfile, l_verbose)
    printlog ("fl_vardq     = "//l_fl_vardq, l_logfile, l_verbose)
    printlog ("fl_oversize  = "//l_fl_oversize, l_logfile, l_verbose)
    printlog ("sci_ext      = "//l_sci_ext, l_logfile, l_verbose)
    if (l_fl_vardq) {
        printlog ("var_ext      = "//l_var_ext, l_logfile, l_verbose)
        printlog ("dq_ext       = "//l_dq_ext, l_logfile, l_verbose)
    }
    printlog ("long_sample  = "//l_lsample, l_logfile, l_verbose)
    printlog ("mos_sample   = "//l_msample, l_logfile, l_verbose)
    printlog ("mosobjsize   = "//l_mosobjsize, l_logfile, l_verbose)
    printlog ("naverage     = "//l_navrg, l_logfile, l_verbose)
    printlog ("function     = "//l_func, l_logfile, l_verbose)
    printlog ("order        = "//l_order, l_logfile, l_verbose)
    printlog ("low_reject   = "//l_low_rej, l_logfile, l_verbose)
    printlog ("high_reject  = "//l_high_rej, l_logfile, l_verbose)
    printlog ("niterate     = "//l_niter, l_logfile, l_verbose)
    printlog ("grow         = "//l_grow, l_logfile, l_verbose)
    printlog ("fl_inter     = "//l_fl_inter, l_logfile, l_verbose)
    printlog ("", l_logfile, l_verbose)

    # Now, we start with a lot of verifications.

    nerror = 0

    # Check if the input file is not empty string

    if ((l_input=="") || (l_input==" ")) {
        printlog ("ERROR - GSSKYSUB: input image(s) or list are not specified",
            l_logfile, yes)
        nerror = nerror+1
    }

    # Check existence of input list 

    if (substr(l_input,1,1) == "@") {
        inlst = substr (l_input, 2, strlen(l_input))
        if (!access(inlst)) {
            printlog ("ERROR - GSSKYSUB: Input list "//inlst//" not found",
                l_logfile, yes)
            nerror = nerror+1
        }
    }

    if (nerror > 0) goto outerror

    # Check if the output file is not an empty string. If so, check outprefix.

    if ((l_output=="") || (l_output==" ") )
        pref = yes
    else if ((l_output!="") || (l_output!=" "))
        pref = no

    if (pref) {
        if ((l_outprefix=="") || (l_outprefix==" ")) {
            printlog ("ERROR - GSSKYSUB: outprefix is not specified",
                l_logfile, yes)
            nerror = nerror+1
        } else if (substr(l_output,1,1)=="@") {
            outlst = substr(l_output,2,strlen(l_output))
            if (!access(outlst)) {
                printlog ("ERROR - GSSKYSUB: Output list "//outlst//\
                    " not found", l_logfile, yes)
                nerror = nerror+1
            }
        }
    }

    # Check if input images exists

    if (strstr("@",l_input) != 0) {
#        inlst = substr (l_input, 2, strlen(l_input))
#        sections ("@"//inlst, > tmpin)
        sections (l_input, > tmpin)
    } else
        files (l_input, sort-, > tmpin)

    count (tmpin) | scan (ninpimg)
    scanfile = tmpin
    while (fscan(scanfile,inimg) != EOF) {
        gimverify (inimg)
        if (gimverify.status>0) {
            printlog ("ERROR - GSSKYSUB: Input image "//inimg//" does not \
                exist", l_logfile, yes)
            nerror = nerror+1
        } 
    }

    scanfile = ""

    # Check if outimages exist and if the the number of input images are 
    # equal to the number of output images

    if (!pref) {
        if (substr(l_output,1,1)=="@") {
            outlst = substr (l_output, 2, strlen(l_output))
            sections ("@"//outlst, > tmpout)
        } else
            files (l_output, sort-, > tmpout)
            
        count(tmpout) | scan(noutimg)
        scanfile = tmpout
        if (ninpimg != noutimg) {
            printlog ("ERROR - GSSKYSUB: Number of input and output images \
                does not match", l_logfile, yes)
            nerror = nerror+1
        }
        while (fscan(scanfile,outimg) !=EOF) {
            if (imaccess(outimg)) {
                printlog ("ERROR - GSSKYSUB: Output image "//outimg//" exists",
                    l_logfile, yes)
                nerror = nerror+1
            }
        }
    } else {
        sections (l_outprefix//"@"//tmpin, > tmpout)
        scanfile = tmpout
        while (fscan(scanfile,outimg) !=EOF) {
            if (imaccess(outimg)) {
                printlog ("ERROR - GSSKYSUB: Output image "//outimg//" exist",
                    l_logfile,yes)
                nerror = nerror+1
            }
        }
    }

    scanfile = ""

    gextverify (l_sci_ext)
    l_sci_ext = gextverify.outext
    if (gextverify.status==1) {
        printlog ("ERROR - GSSKYSUB: sci_ext is an empty string",
            l_logfile, yes)
        nerror = nerror+1
    }
    if (l_fl_vardq) {
        gextverify (l_var_ext)
        l_var_ext = gextverify.outext
        if (gextverify.status==1) {
            printlog ("ERROR - GSSKYSUB: var_ext is an empty string",
                l_logfile, yes)
            nerror = nerror+1
        }
        gextverify (l_dq_ext)
        l_dq_ext = gextverify.outext
        if (gextverify.status==1) {
            printlog ("ERROR - GSSKYSUB: dq_ext is an empty string",
                l_logfile, yes)
            nerror = nerror+1
        }
    }

    # If nerror > 0, go out, else continue

    if (nerror > 0)
        goto outerror

    nerror = 0

    # The slit length used is 1.05 times bigger than the slit width defined
    # in the MDF file by default, to take into account the approximations in
    # the distortion mapping
    if (l_fl_oversize)
        slenfact = 1.05
    else
        slenfact = 1.0

    # Check input images

    scanfile = tmpin

    i = 0

    while (fscan(scanfile,inimg) != EOF) {
        i+=1
        suf = substr(inimg,strlen(inimg)-3,strlen(inimg))
        if (suf!="fits")
            inimg = inimg//".fits"

        # Check if the images are MEF
        imgets (inimg//"[0]", "EXTEND", >& "dev$null")
        if (imgets.value=="F" || imgets.value=="0") {
            printlog ("ERROR - GSSKYSUB: image "//inimg//" is not a MEF file",
                l_logfile,yes)
            nerror = nerror+1
        }
        # Check image type
        imgets (inimg//"[0]", "MASKTYP", >& "dev$null")
        msktype = int(imgets.value)
        if (msktype != 1) {
            printlog ("ERROR - GSSKYSUB: "//inimg//" has MASKTYP other than \
                LONGSLIT or MOS mode", l_logfile, yes)
            nerror = nerror+1
        } 
        # Check the MDF file.  Check where is it.
        imgets (inimg//"[0]", "NEXTEND", >& "dev$null")
        numext = int(imgets.value)
        for (k=1; k<=numext; k+=1) {
            keypar (inimg//"["//k//"]", keyword="EXTNAME", silent=yes)
            if (keypar.value == "MDF" && keypar.found) {
                mdf = yes
                mdfpos[i] = k
            }
        }
        if (!mdf) {
            printlog ("ERROR - GSSKYSUB: Input file "//inimg//" does not \
                have an attached MDF table", l_logfile, yes)
            nerror = nerror+1
        }
        # Get obsmode
        imgets (inimg//"[0]", "OBSMODE", >& "dev$null")
        obsmode[i] = imgets.value
        if (imgets.value == "" || imgets.value == " " || imgets.value == "0") {
            printlog ("ERROR - GSSKYSUB: could not find the OBSMODE keyword \
                parameter in "//inimg, l_logfile, yes)
            nerror = nerror+1
        }
        # Check how many SCI extension we have
        imgets (inimg//"[0]", "NSCIEXT", >& "dev$null")
        nsciext[i] = int(imgets.value)
        if (nsciext[i] == 0) {
            printlog ("ERROR - GSSKYSUB: Number of SCI extensions unknown \
                for image "//inimg, l_logfile, yes)
            nerror = nerror+1
        } else {
            for (j=1; j<=nsciext[i]; j+=1) {
                if (!imaccess(inimg//"["//l_sci_ext//","//str(j)//"]")) {
                    printlog ("ERROR - GSSKYSUB: Could not access "//inimg//\
                        "["//l_sci_ext//","//str(j)//"]", l_logfile, yes)
                    nerror = nerror+1
                } else {
                    # Check if the MDFROW keyword exist inside each extension
                    # (MOS only)
                    if (obsmode[i] == "MOS") {
                        imgets (inimg//"["//l_sci_ext//","//str(j)//"]",
                            "MDFROW", >& "dev$null")
                        if ((imgets.value == "") || (imgets.value == " ") || \
                            (imgets.value == "0")) {
                            printlog ("ERROR - GSSKYSUB: could not find the \
                                MDFROW keyword parameter in "//inimg//\
                                "["//l_sci_ext//","//str(j)//"]",
                                l_logfile,yes)
                            nerror = nerror+1
                        }
                    }
                    # If fl_vardq=yes, check if there are access the images
                    if (l_fl_vardq) {
                        if (!imaccess(inimg//"["//l_var_ext//","//str(j)//\
                            "]")) {
                            printlog ("ERROR - GSSKYSUB: Could not access "//\
                                inimg//"["//l_var_ext//","//str(j)//"]",
                                l_logfile, yes)
                            nerror = nerror+1
                        } 
                        if (!imaccess(inimg//"["//l_dq_ext//","//str(j)//\
                            "]")) {
                            printlog ("ERROR - GSSKYSUB: Could not access "//\
                                inimg//"["//l_dq_ext//","//j//"]",
                                l_logfile, yes)
                            nerror = nerror+1
                        }
                    }
                }
            }
        }
        # Check if the images where GSTRANSFORMED. If not, out!
        imgets (inimg//"[0]", "GSTRANSF", >& "dev$null")
        if (imgets.value == "" || imgets.value == " " || imgets.value=="0") {
            printlog ("WARNING - GSSKYSUB: Image "//inimg//" has not been \
                transformed", l_logfile, yes)
            printlog ("                    It is recommended to run \
                GSTRANSFORM first", l_logfile, yes)
        }
    } #end of while loop over input images to check the parameters.

    # If no error found, continue, else stop

    if (nerror > 0)
        goto outerror

    scanfile = ""

    joinlines (tmpin//","//tmpout, output=joinlst, delim=" ", \
        missing="Missing", maxchar=161, shortest+, verbose-)
    delete (tmpin//","//tmpout, verify-, >& "dev$null")

    printlog (" ", l_logfile, l_verbose)
    printlog ("GSSKYSUB: Sky-subtraction started", l_logfile, l_verbose)
    printlog (" ", l_logfile, l_verbose)

    scanfile = joinlst

    k = 0

    while (fscan(scanfile,inimg,outimg) != EOF) {
        # Create tmp FITS file names used within this loop
        mdffile = mktemp("tmpmdf")

        k+=1
        suf = substr (outimg, strlen(outimg)-3, strlen(outimg))
        if (suf!="fits")
            outimg = outimg//".fits"
        suf = substr (inimg, strlen(inimg)-3, strlen(inimg))
        if (suf != "fits")
            inimg = inimg//".fits"
        #
        # copy MDF to outfile
        #
        tcopy (inimg//"["//mdfpos[k]//"]", mdffile//".fits", verbose-)
        wmef (mdffile//".fits", outimg, extnames="MDF", verbose-, phu=inimg,
            >& "dev$null")
        # Define binning and scale
        # Which instrument?
        imgets (inimg//"[0]", "INSTRUME", >& "dev$null")
        if (imgets.value == "0") {
            printlog ("ERROR - GSSKYSUB: Instrument keyword not found.",
                l_logfile, verbose+)
            goto outerror
        }
        if (imgets.value == "GMOS-S")
            inst = 2
        else
            inst = 1

        # Type of detector?
        imgets (inimg//"[0]", "DETTYPE", >& "dev$null")
        if (imgets.value == "0") {
            printlog ("ERROR - GSEXTRACT: DETTYPE keyword not found.",
                l_logfile, verbose+)
            goto outerror
        }
        if (imgets.value == "SDSU II CCD") { # EEV CCds
            iccd = 1
        } else if (imgets.value == "SDSU II e2v DD CCD42-90") {# New e2VDD CCDs
            iccd = 2
        } else if (imgets.value == "S10892" || imgets.value == "S10892-N") { # Hamamatsu CCDs
            iccd = 3
        }

        imgets (inimg//"[0]", "CCDSUM")
        xbin = int(substr(imgets.value,1,1))
        ybin = int(substr(imgets.value,3,3))
        xscale = pixscale[inst,iccd] * real(xbin)
        yscale = pixscale[inst,iccd] * real(ybin)
        ll_mosobjsize = l_mosobjsize/yscale
        ll_fl_inter = l_fl_inter
        printlog ("  Input image     = "//inimg, l_logfile, l_verbose)
        printlog ("  Output image    = "//outimg, l_logfile, l_verbose)
        printlog ("  Observing mode  = "//obsmode[k], l_logfile, l_verbose)
        printlog ("  Number of slits = "//nsciext[k],
            l_logfile, verbose=l_verbose)
        # Loops over science extensions.
        for (j=1; j<=nsciext[k]; j+=1) {
            # Create tmp FITS file names used within this loop
            tmpfit = mktemp("tmpfit")
            tmpsci = mktemp("tmpsci")
            if (l_fl_vardq) {
                tmpvar = mktemp("tmpvar")
                tmpdq = mktemp("tmpdq")
            }

            # X and Y sizes, in pixel, of each slit
            imgets (inimg//"["//l_sci_ext//","//j//"]", "i_naxis1")
            nxpix = int(imgets.value)
            imgets (inimg//"["//l_sci_ext//","//j//"]", "i_naxis2")
            nypix = int(imgets.value)
            # Center of each slit
            xcen = real(nxpix)/2.
            ycen = real(nypix)/2.
            if (obsmode[k] == "MOS") {
                # Read the MDFROW keyword associated with the row in the MDF
                # Obtaining MDFROW value from current science extension
                imgets (inimg//"["//l_sci_ext//","//j//"]", "MDFROW",
                    >& "dev$null")
                mdfrow = int(imgets.value)
                # Read slit type
                tabpar (mdffile//".fits", "slittype", row=mdfrow, \
                    >& "dev$null")
                slittype = tabpar.value
                tabpar (mdffile//".fits","slitsize_my", row=mdfrow,
                    >& "dev$null")
                # Transform slitsize from mm to pixel coordinates
                sl_size = real(tabpar.value) * asecmm/yscale
                
                # Calc slit center
                s_cen = sl_size * slenfact/2.
                
                # read the spec. pos. inside the slit. Transform to pixels 
                # (from arcsec).
                # This is in fact the slit position relative to the object
                
                tabpar (mdffile//".fits", "slitpos_y", row=mdfrow,
                    >& "dev$null")
                s_pos = s_cen - real(tabpar.value)/yscale
                # if slittype=RECTANGLE do the fit else do nothing
                if (slittype == "rectangle") {

                    # Determine the sky sample
                    slitwidth = slenfact*sl_size
                    y1 = int(max(1.,0.5*((1.-l_msample)*slitwidth))+0.5)
                    y4 = int(slitwidth-y1+1.5)
                    y2 = int(max(y1,s_pos-ll_mosobjsize/2.+0.5))
                    y3 = int(min(y4,s_pos+ll_mosobjsize/2.+0.5))

                    # Set background sample
                    l_bsample = y1//":"//y2//","//y3//":"//y4
                    # Sky subtracting using fit1d
                    if (l_fl_inter && j==1) {
                        print ("")
                        print ("Using fit1d for fitting the sky level")
                        print ("  Select the a column to fit")
                        print ("  Exit fitting of that column with q")
                        print ("  To fit the rest of the columns \
                            on-interactively, select q as column number")
                        print ("")
                    }
                    printlog ("  Slit #"//str(j)//"; slit type = "//\
                        slittype//" ; background sample = ["//l_bsample//"]",
                        l_logfile, l_verbose)
                    if (j>1 && ll_fl_inter) {
                        l_fl_answer = fl_answer
                        ll_fl_inter = l_fl_answer
                    } else
                        l_fl_answer = ll_fl_inter
                        
                    fit1d (inimg//"["//l_sci_ext//","//str(j)//"]", tmpfit,
                        "fit", axis=2, interactive=l_fl_answer,
                        sample=l_bsample, naverage=l_navrg, function=l_func,
                        order=l_order, low_reject=l_low_rej,
                        high_reject=l_high_rej, niterate=l_niter, grow=l_grow,
                        graphics="stdgraph", cursor="")
                        
                    # Final sky subtracted image
                    imarith (inimg//"["//l_sci_ext//","//j//"]", "-", tmpfit,
                        tmpsci, verbose-)
                        
                    # If fl_vardq, propagate the variance and DQ plane.
                    # New variance frame is var = var+(rms of the fit)**2
                    # From fit1d it is impossible to obtain the rms of the fit.
                    # Then, we use gfit1d to obtain it
                    if (l_fl_vardq) {
                        # Copy variation plane to a temporary image
                        imcopy (inimg//"["//l_var_ext//","//j//"]", tmpvar,
                            verbose-)
                        # Calculate rms of the fit using gfit1d.
                        gfit1d (inimg//"["//l_sci_ext//","//j//"]", \
                            tmptab//".tab",
                            function=l_func, order=l_order, xmin=1., 
                            xmax=nxpix, ps=yes, errorpars.errcolum="",
                            errorpars.errtype="uniform", errorpars.resample=no,
                            samplepars.axis=2, samplepars.sample=l_bsample,
                            samplepars.naverage=l_navrg,
                            samplepars.low_reject=l_low_rej,
                            samplepars.high_reject=l_high_rej,
                            samplepars.niterate=l_niter, 
                            samplepars.grow=l_grow, interactive=l_fl_answer, 
                            device="stdgraph", cursor="", >& "dev$null")
                        #Get rms value from tmptab using tabpar
                        tabpar (tmptab//".tab", "rms", 1)
                        rmsval = real(tabpar.value)
                        # New var image
                        tmpval = rmsval*rmsval
                        printlog ("         RMS of the fit (for VAR plane \
                            propagation) = "//rmsval, l_logfile, l_verbose)
                        imarith (tmpvar, "+", tmpval, tmpvar, verbose=no)
                        # DQ extension - the same input dq plane
                        imcopy (inimg//"["//l_dq_ext//","//j//"]", tmpdq,
                            verbose-)
                        tdelete (tmptab, verify-, >& "dev$null")
                    }
                } else {
                    printlog ("Slit #"//j//" ; type = "//slittype//" --- not \
                        sky-subtracted", l_logfile, l_verbose)
                    imcopy (inimg//"["//l_sci_ext//","//j//"]", tmpsci, \
                        verbose-)
                    if (l_fl_vardq) {
                        imcopy (inimg//"["//l_var_ext//","//j//"]", tmpvar, \
                            verbose-)
                        imcopy (inimg//"["//l_dq_ext//","//j//"]", tmpdq, \
                            verbose-)
                    }
                }
            } 
            if (obsmode[k] == "LONGSLIT") {
                # The background sample must of the format y1:y2,y3:y4 in 
                # long_sample. Defined by the user.
                
                l_bsample = l_lsample
                printlog ("  Slit #"//str(j)//" ; background sample = \
                    ["//l_bsample//"]", l_logfile, l_verbose)
                fit1d (inimg//"["//l_sci_ext//","//j//"]", tmpfit, "fit",
                    axis=2, interactive=l_fl_inter, sample=l_bsample,
                    naverage=l_navrg, function=l_func, order=l_order,
                    low_reject=l_low_rej, high_reject=l_high_rej,
                    niterate=l_niter, grow=l_grow, graphics="stdgraph",
                    cursor="")
                    
                # Final sky subtracted image
                imarith (inimg//"["//l_sci_ext//"]", "-", tmpfit, tmpsci,
                    verbose-)
                    
                # If fl_vardq, propagate the variance and DQ plane.
                # New variance frame is var = var+(rms of the fit)**2
                # From fit1d it is impossible to obtain the rms of the fit. 
                # Then, we use gfit1d to obtain it.
                
                if (l_fl_vardq) {
                    # Copy variation plane to a temporary image
                    imcopy (inimg//"["//l_var_ext//","//j//"]", tmpvar,
                        verbose-)
                    # Calculate rms of the fit using gfit1d.
                    gfit1d (inimg//"["//l_sci_ext//","//j//"]", tmptab//".tab",
                        function=l_func, order=l_order, xmin=1., xmax=nxpix,
                        ps=yes, errorpars.errcolum="",
                        errorpars.errtype="uniform",
                        errorpars.resample=no, samplepars.axis=2,
                        samplepars.sample=l_bsample,
                        samplepars.naverage=l_navrg,
                        samplepars.low_reject=l_low_rej,
                        samplepars.high_reject=l_high_rej,
                        samplepars.niterate=l_niter, samplepars.grow=l_grow,
                        interactive=l_fl_inter, device="stdgraph", cursor="",
                        >& "dev$null")
                    #Get rms value from tmptab using tabpar
                    tabpar (tmptab//".tab", "rms", 1)
                    rmsval = real(tabpar.value)
                    # New var image
                    tmpval = rmsval*rmsval
                    printlog ("         RMS of the fit (for VAR plane \
                        propagation) = "//rmsval, l_logfile, l_verbose)
                    imarith (tmpvar, "+", tmpval, tmpvar, verbose=no)
                    # DQ extension - the same input dq plane
                    imcopy (inimg//"["//l_dq_ext//","//j//"]", tmpdq, verbose-)
                    tdelete (tmptab, verify-, >& "dev$null")
                }
            }
            imgets (outimg//"[0]", "NEXTEND")
            nextnd = int(imgets.value)
            if (l_fl_vardq) {
                fxinsert (tmpsci//".fits,"//tmpvar//".fits,"//tmpdq//".fits",
                    outimg//"["//nextnd//"]", "0", verbose-, >& "dev$null")
                gemhedit (outimg//"[0]", "NEXTEND", (nextnd+3), "", delete-)
                gemhedit (outimg//"["//(nextnd+1)//"]", "EXTNAME", l_sci_ext, 
                    "", delete-)
                gemhedit (outimg//"["//(nextnd+1)//"]", "EXTVER", j, 
                    "", delete-)
                gemhedit (outimg//"["//(nextnd+2)//"]", "EXTNAME", l_var_ext, 
                    "", delete-)
                gemhedit (outimg//"["//(nextnd+2)//"]", "EXTVER", j, 
                    "", delete-)
                gemhedit (outimg//"["//(nextnd+3)//"]", "EXTNAME", l_dq_ext, 
                    "", delete-)
                gemhedit (outimg//"["//(nextnd+3)//"]", "EXTVER", j, 
                    "", delete-) 
                imdelete (tmpfit//".fits,"//tmpsci//".fits,"//tmpvar//\
                    ".fits,"//tmpdq//".fits", verify-, >& "dev$null")
            } else {
                fxinsert (tmpsci//".fits",outimg//"["//nextnd//"]", "0",
                    verbose-, >& "dev$null")
                gemhedit (outimg//"[0]", "NEXTEND", (nextnd+1), "", delete-)
                gemhedit (outimg//"["//(nextnd+1)//"]", "EXTNAME", l_sci_ext, 
                    "", delete-)
                gemhedit (outimg//"["//(nextnd+1)//"]", "EXTVER", j, 
                    "", delete-)
                imdelete (tmpfit//".fits,"//tmpsci//".fits", verify-,
                    >& "dev$null")
            }
        } # end of for-loop over science extensions
        gemdate ()
        gemhedit (outimg//"[0]", "GSSKYSUB", gemdate.outdate,
            "UT Time stamp for GSSKYSUB")
        gemhedit (outimg//"[0]", "GEM-TLM", gemdate.outdate,
            "UT Last modification with GEMINI")
        gemhedit (outimg//"[0]", "NSCIEXT", nsciext[k], "", delete-)
        imdelete (mdffile//".fits", verify-, >& "dev$null")
        printlog ("", l_logfile, l_verbose)
    } # End of big while loop

    status = 0
    goto clean

outerror:
    # Exit with error
    status = 1
    goto clean

clean:
    # clean up
    scanfile = ""
    date | scan (sdate)
    if (status != 0) {
        delete (tmpin//","//tmpout, verify-, >& "dev$null")
        printlog ("ERROR - GSSKYSUB: "//nerror//" error(s) found. Exiting.",
            l_logfile, yes)
        printlog ("", l_logfile, l_verbose)
        printlog ("GSSKYSUB done. Exit status bad. -- "//sdate,
            l_logfile, l_verbose)
    } else {
        delete (joinlst, verify-, >& "dev$null")
        printlog ("GSSKYSUB done. Exit status good -- "//sdate,
            l_logfile, l_verbose)
    }
    printlog ("------------------------------------------------------------\
        -------------------", l_logfile, l_verbose)
    printlog ("", l_logfile, l_verbose)
    fl_answer = INDEF

end

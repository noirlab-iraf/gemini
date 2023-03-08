# Copyright(c) 2005-2011 Association of Universities for Research in Astronomy, Inc.

procedure msabsflux (inimages, cal, stdname)

# This routine corrects the science target spectra for telluric features.  Both the
# target spectrum and the calibrator for the correction must have been extracted
# using "nsextract".  The standard needs to be one of the known calibration stars.
# 
# Version:   October 26, 2005    KV created this task, mostly a copy of "mstelluric"
#            February 3, 2006    KV bug fix, path to file in my account not changed to midir$data
#            December 4, 2006    KV added small changes for use with "msslice", added fl_plots parameter

    
char    inimages    {prompt="Input T-ReCS or Michelle spectra"}     # OLDP-1-input-primary-single-prefix=s
char    cal         {prompt = "Telluric calibration spectrum"}
char    stdname     {prompt="Name of the calibration standard star, or 'list'"}
char    outimages   {"",prompt="Output spectra"}                    # OLDP-1-output
char    outpref     {"a",prompt="Prefix for output images"}         # OLDP-4
char    outtype     {"fnu",prompt="Type of output spectrum: fnu|flambda|lambda*flambda"}
bool    fl_bbody    {no,prompt="use a blackbody shape for the calibration object"}
real    bbody       {1.0, prompt="Temperature of calibrator for B-body fit"}
bool    xcorr       {no, prompt="Cross correlate for shift (telluric)"}
int     lag         {10, prompt="Cross correlation lag (in pixels) (telluric)"}
real    shift       {0., prompt="Initial shift of calibration spectrum (pixels) (telluric)"}
real    scale       {1., prompt="Initial scale factor multiplying airmass ratio (telluric)"}
real    dshift      {0.5, prompt="Initial shift search step (telluric)"}
real    dscale      {0.1, prompt="Initial scale factor search step (telluric)"}
real    threshold   {0.01, prompt="Threshold value for telluric calibration (telluric)"}
bool    fl_inter    {no,prompt="Interactive tweaking of correction (telluric)"}
bool    fl_plots    {yes,prompt="Show plots during execution"}
char    logfile     {"",prompt="Log file name"}                     # OLDP-1
bool    verbose     {yes,prompt="Verbose"}                          # OLDP-4
int     status      {0,prompt="Exit status: (0=good, >0=bad)"}      # OLDP-4
struct  *scanfile   {"",prompt="Internal use only"}                 # OLDP-4

begin

    char    l_inimages = ""
    char    l_cal = ""
    char    l_stdname = ""
    char    l_outimages = ""
    char    l_outpref = ""
    char    l_outtype = ""
    char    l_logfile = ""
    int     l_lag
    real    l_bbody, l_shift, l_scale, l_dshift, l_dscale, l_threshold
    bool    l_fl_bbody, l_xcorr, l_fl_inter, l_fl_plots, l_verbose
    
    int     junk
    
    char    airmass1, airmass2, filename, keyfound
    char    in[100], out[100], instrument, temp, errmsg
    char    tmpfile, tmpinimg, tmpoutimg, outputstr
    char    tmpim1, tmpim2, tmpim3, tmpim4, tmpim5, tmpim6, tmpim7
    char    tmpfile1, tmpfile2, tmpfile3
    char    paramstr, tmplog
    char    val1,val2
    int     i, j, k
    int     nimages, maximages, noutimages, nbad, ngood, ist, ifin
    int     wheremin
    real    wave[320], bbflux[320], value[320], c1, c2, diff, min
    real    wavel[461],fnu[461],wl[320]
    real    crval1, cdelt1, meanval, bbvalue, ymin, ymax
    real    cutlevel1, cutlevel2, wavelnorm

    char    standardname, stdfile
    bool    there
    int     istd, itype
    real    xmag

    cache ("gemdate", "gemextn")

    tmpfile = mktemp("tmpfile")
    tmpinimg = mktemp("tmpinimg")
    tmpoutimg = mktemp("tmpoutimg")

    junk = fscan (inimages, l_inimages)
    junk = fscan (cal, l_cal)
    junk = fscan (stdname, l_stdname)
    junk = fscan (outimages, l_outimages)
    junk = fscan (outpref, l_outpref)
    junk = fscan (outtype, l_outtype)
    l_fl_bbody  = fl_bbody
    l_bbody     = bbody
    l_xcorr     = xcorr
    l_lag       = lag
    l_shift     = shift
    l_scale     = scale
    l_dshift    = dshift
    l_dscale    = dscale
    l_threshold = threshold
    l_fl_inter  = fl_inter
    l_fl_plots  = fl_plots
    junk = fscan (logfile, l_logfile)
    l_verbose   = verbose
    
    # Initialize
    status=0
    ngood=0
    maximages=100
    nimages=0

    # Create the list of parameter/value pairs.  One pair per line.
    # All lines combined into one string.  Line delimiter is '\n'.
    paramstr =  "inimages       = "//inimages.p_value//"\n"
    paramstr += "cal            = "//cal.p_value//"\n"
    paramstr += "stdname        = "//stdname.p_value//"\n"
    paramstr += "outimages      = "//outimages.p_value//"\n"
    paramstr += "outpref        = "//outpref.p_value//"\n"
    paramstr += "outtype          = "//outtype.p_value//"\n"
    paramstr += "fl_bbody   = "//fl_bbody.p_value//"\n"
    paramstr += "bbody          = "//bbody.p_value//"\n"
    paramstr += "xcorr          = "//xcorr.p_value//"\n"
    paramstr += "lag            = "//lag.p_value//"\n"
    paramstr += "shift          = "//shift.p_value//"\n"
    paramstr += "scale          = "//scale.p_value//"\n"
    paramstr += "dshift         = "//dshift.p_value//"\n"
    paramstr += "dscale         = "//dscale.p_value//"\n"
    paramstr += "threshold      = "//threshold.p_value//"\n"
    paramstr += "fl_inter       = "//fl_inter.p_value//"\n"
    paramstr += "fl_plots       = "//fl_plots.p_value//"\n"
    paramstr += "logfile        = "//logfile.p_value//"\n"
    paramstr += "verbose        = "//verbose.p_value

    # Assign a logfile name if not given.  Open logfile and start log.
    # Write parameter/value pairs ("paramstr") to log.
    gloginit (l_logfile, "msabsflux", "midir", paramstr, fl_append+,
        verbose=l_verbose)
    if (gloginit.status != 0) {
        status = gloginit.status
        goto exitnow
    }
    l_logfile = gloginit.logfile.p_value

    l_outtype = strlwr(l_outtype)
    itype = 0
    if (l_outtype == "fnu") itype=1
    if (l_outtype == "flambda") itype=2
    if (l_outtype == "lambda*flambda") itype=3

    if (itype == 0) {
        glogprint (l_logfile,"msabsflux","status",type="error",errno=121,
            str="ERROR - MSABSFLUX: unrecognised outtype option "/\
            l_outtype//".",verbose+)
        status = 121
        goto clean
    }

    if (l_fl_bbody) {
        if (l_bbody == 1.0) {
            glogprint (l_logfile, "mstelluric", "status", type="error",
                errno=121, str="ERROR - MSTELLURIC: Blackbody Temperature is \
                not set!", verbose+)
            status = 1
            goto clean
        }
    } else {

        # Read in the standard star flux density values

        # check for a blank name string

        if (l_stdname == "") {
            glogprint (l_logfile, "msabsflux", "status", type="error",
                errno=121, str="ERROR - MSABSFLUX: no stadard star name \
                entered.", verbose+)
            status = 1
            goto clean
        }

        stdfile = "midir$data/standards.list"

        scanfile = stdfile

        there = access(stdfile)

        if (!there) {
            glogprint (l_logfile, "msabsflux", "status", type="error",
                errno=121, str="ERROR - Cannot find the standard star fluxes \
                file.", verbose+)
            status = 121
            goto clean
        }

        i = 1
        istd = 0

        # Send back the list of standard star names 
        # if it is entered as "list"

        if (strstr(l_stdname,"list") > 0) {
            while(i > 0) {
                i = fscan (scanfile, standardname, stdfile)
                if (i == 2) {
                    print (standardname)
                }
            }
            goto clean  
        }

        while ((i > 0) && (istd == 0)) {
            i = fscan (scanfile, standardname, stdfile)
            if (i == 2) {
                if (strstr (l_stdname, standardname) > 0) {
                    istd = 1
                }
            }
        }

        if (istd == 0) {
            glogprint (l_logfile,"msabsflux","status",type="error",errno=121,
                str="ERROR - Unrecognised standard star name "//l_stdname//" .",
                verbose+)
            status = 1
            goto clean
        }

        stdfile = "midir$data/"//substr(stdfile,3,strlen(stdfile))

        scanfile = stdfile

        there = access (stdfile)

        if (!there) {
            glogprint (l_logfile, "msabsflux", "status", type="error",
                errno=121, str="ERROR - Cannot find the standard star fluxes \
                file.", verbose+)
            status = 1
            goto clean
        }

        i = fscan (scanfile, standardname, stdfile, junk)
        i = fscan (scanfile, standardname, stdfile, junk)
        i = fscan (scanfile, standardname, stdfile, junk)
        i = fscan (scanfile, standardname, stdfile, junk)

        i = 1
        j = 1

        while (i > 0) {
            i = fscan (scanfile, standardname, stdfile, junk)
            if (i == 3) {
                wavel[j] = real(standardname)
                xmag = real(stdfile)
                fnu[j] = 3631. / 10.**(xmag*0.4)
                if (itype == 2) {
                    fnu[j] = fnu[j] * 1.e-23 * 299792458. * 1.e+06 / \
                        (wavel[j]*wavel[j])
                }
                if (itype == 3) {
                    fnu[j] = fnu[j] * 1.e-23 * 299792458. * 1.e+06 / wavel[j]
                }
                j = j+1
            }
        }

        if (j != 462) {
            glogprint (l_logfile, "msabsflux", "status", type="error",
                errno=121, str="ERROR - There was a problem reading the \
                standard star fluxes.", verbose+)
            status = 1
            goto clean
        }
    }

    # Load up the array of input file names
    gemextn (l_inimages, check="", process="none", index="", extname="",
        extversion="", ikparams="", omit="extension", replace="", 
        outfile=tmpfile, logfile=l_logfile, glogpars="", verbose=l_verbose)
    gemextn ("@"//tmpfile, check="exists,mef", process="none", index="",
        extname="", extversion="", ikparams="", omit="", replace="",
        outfile=tmpinimg, logfile=l_logfile, glogpars="", verbose=l_verbose)
    nimages = gemextn.count
    delete (tmpfile, ver-, >& "dev$null")

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

        glogprint (l_logfile, "msabsflux", "status", type="error", errno=status,
            str=errmsg, verbose+)
        goto clean

    } else {
        scanfile = tmpinimg
        i = 0
        while (fscan(scanfile,filename) != EOF) {
            i += 1
            in[i] = filename
        }
        scanfile = ""
        if (i != nimages) {
            glogprint (l_logfile, "msabsflux", "status", type="error",
                errno=99, str="Error while counting the input images.",
                verbose+)
            status = 1
            goto clean
        }
    }

    # Load up the array of output file names
    if (l_outimages != "")
        outputstr = l_outimages
    else if (l_outpref != "") {
        gemextn ("@"//tmpinimg, check="", process="none", index="", extname="",
            extversion="", ikparams="", omit="path", replace="",
            outfile=tmpoutimg, logfile=l_logfile, glogpars="",
            verbose=l_verbose)
        outputstr = l_outpref//"@"//tmpoutimg
    } else {
        glogprint (l_logfile, "msabsflux", "status", type="error", errno=121,
            str="Neither output image name nor output prefix is defined.",
            verbose+)
        status = 1
        goto clean
    }

    gemextn (outputstr, check="", process="none", index="", extname="",
        extversion="", ikparams="", omit="extension", replace="",
        outfile=tmpfile, logfile=l_logfile, glogpars="", verbose=l_verbose)
    delete (tmpoutimg, ver-, >& "dev$null")
    gemextn ("@"//tmpfile, check="absent", process="none", index="",
        extname="", extversion="", ikparams="", omit="", replace="",
        outfile=tmpoutimg, logfile=l_logfile, glogpars="", verbose=l_verbose)
    noutimages = gemextn.count
    delete (tmpfile, ver-, >& "dev$null")

    if ((gemextn.fail_count > 0) || (noutimages == 0) || \
        (noutimages != nimages)) {

        if (gemextn.fail_count > 0) {
            errmsg = gemextn.fail_count//" image(s) already exist(s)."
            status = 102
        } else if (noutimages == 0) {
            errmsg = "Maximum number of output images exceeded:"//str(maximages)
            status = 121
        } else if (noutimages != nimages) {
            errmsg = "Different number of input images ("//nimages//") and \
                output images ("//noutimages//")."
            status = 121
        }

        glogprint (l_logfile, "msabsflux", "status", type="error", 
            errno=status, str=errmsg, verbose+)
        goto clean

    } else {
        scanfile = tmpoutimg
        i = 0
        while (fscan(scanfile,filename) != EOF) {
            i += 1
            out[i] = filename//".fits"
        }
        scanfile = ""
        if (i != noutimages) {
            glogprint (l_logfile, "msabsflux", "status", type="erro",
                errno=99, str="Error while counting the output images.",
                verbose+)
            status = 1
            goto clean
        }
    }

    # Make sure the calibration spectrum exists
    gemextn (l_cal, check="exists,mef", process="none", index="", extname="",
        extversion="", ikparams="", omit="", replace="", outfile="",
        logfile=l_logfile, glogpars="", verbose=l_verbose)
    if ((gemextn.fail_count > 0) || (gemextn.count == 0)) {
        if (gemextn.fail_count > 0) {
            errmsg = "Calibration file not found."
            status = 101
        } else if (gemextn.count == 0) {
            errmsg = "ERROR - MSABSFLUX: No calibration file defined."
            status = 121
        }

        glogprint (l_logfile, "msabsflux", "status", type="error",
            errno=status, str=errmsg, verbose+)
        goto clean
    }

    hselect (l_cal//"[SCI,1]", "CRVAL1", yes) | scan (standardname)
    crval1 = real(standardname)
    stdfile = ""
    hselect (l_cal//"[SCI,1]", "CD1_1", yes) | scan(stdfile)
    cdelt1 = real(stdfile)

    if (!l_fl_bbody) {

        # get the standard star spectrum values

        standardname = ""
        crval1 = crval1 / 10000.
        cdelt1 = cdelt1 / 10000.

        if ((crval1 < 5.) || (crval1 > 25.) || (cdelt1 <= 0.)) {
            glogprint (l_logfile, "msabsflux", "status", type="error",
                errno=121, str="ERROR - bad wavelength keywords CRVAL1/CD1_1 \
                in the spectrum header.", verbose+)
            status = 1
            goto clean
        }

        tmpfile1 = mktemp ("tmpfile1")
        for (j=1; j <= 320; j=j+1) {
            wl[j] = crval1 + (j-1)*cdelt1
            k = int(20. * (wl[j]-7.0)) 
            if ((k < 1) || (k > 460)) {
                glogprint (l_logfile, "msabsflux", "status", type="string",
                    errno=0, str="WARNING: Spectrum point outside standard \
                    star template wavelength", verbose+)
                glogprint (l_logfile, "msabsflux", "status", type="string",
                    errno=0, str="WARNING: range; setting to zero.", verbose+)
                value[j] = 0.
                print (wl[j], value[j], >> tmpfile1)
            } else {
                value[j] = fnu[k] + \
                    (fnu[k+1]-fnu[k]) * (wl[j]-wavel[k])/(wavel[k+1]-wavel[k])
                print (wl[j], value[j], >> tmpfile1)
            }
        }

        tmpim7 = mktemp ("tmpim7")
        rspectext (tmpfile1, tmpim7)
        delete (tmpfile1, >& "dev$null")
    }

    # Do the work
    i = 1
    while (i <= nimages) {
        glogprint (l_logfile, "msabsflux", "status", type="string",
            str="Processing "//in[i], verbose=l_verbose)

        keyfound = ""
        hselect (in[i]//"[0]", "MSABSFLUX", yes) | scan(keyfound)
        if (keyfound != "") {
            glogprint (l_logfile, "msabsflux", "status", type="warning",
                errno=0, str="File "//in[i]//" has already been telluric \
                corrected.", verbose=l_verbose)
            goto nextimage
        }

        glogprint (l_logfile, "msabsflux", "task", type="string",
            str="  "//in[i]//" --> "//out[i], verbose=l_verbose)

        # check the primary FITS header
        hselect (in[i]//"[0]", "INSTRUME", yes) | scan(instrument)
        glogprint (l_logfile, "msabsflux", "engineering", type="string",
            str="Instrument is: "//instrument, verbose=l_verbose)

        meanval = 1.0
        tmpim1 = mktemp ("tmpim1")
        tmpim2 = mktemp ("tmpim2")
        tmpim3 = mktemp ("tmpim3")
        tmpim4 = mktemp ("tmpim4")
        tmpim5 = mktemp ("tmpim5")
        tmpfile1 = mktemp ("tmpfile1")

        if (l_fl_bbody) {
            c1 = 3.74185E-5
            c2 = 1.43883

            bbvalue = 0
            for (j=1; j<=320; j+=1){
                # bb calculations done in cm
                wave[j] = (crval1 + (cdelt1*(j-1))) / 1.0e8
                value[j] = c2 / (wave[j]*l_bbody)
                bbflux[j] = c1 / (wave[j]**5 * exp(value[j]-1.0))
                bbvalue = bbvalue + (bbflux[j]/320.0)
                # convert wavelength back to microns
                wave[j] = (crval1 + (cdelt1*(j-1))) / 1.0e4
            }

            # Normalize and write to file
            delete (tmpfile1, ver-, >& "dev$null")
            for (j=1; j<=320; j+=1){
                bbflux[j] = wave[j] * (bbflux[j]/bbvalue)
                print (wave[j], bbflux[j], >> tmpfile1)
            }

            glogprint (l_logfile, "mstelluric", "science", type="string",
                str="Removing the blackbody shape from the standard",
                verbose=l_verbose)

            rspectext (tmpfile1, tmpim2)

            imarith (l_cal//"[sci,1]","/",tmpim2//".fits", tmpim1, ver-)

            # Set the normalization wavelength for low-resolution N-band,
            # low-resolution Q-band or higher resolution N-band according to 
            # the range of wavelengths in the spectrum.

            glogprint (l_logfile, "msabsflux", "science", type="string",
                str="Wavelength range: "//str(wave[1])//" "//str(wave[320]),
                verbose=l_verbose)
            wavelnorm = wave[320] - wave[1]
            if ((wavelnorm > 3.0) && (wave[1] < 10.)) wavelnorm = 11.7
            if ((wavelnorm > 3.0) && (wave[1] > 15.)) wavelnorm = 20.5
            if (wavelnorm < 1.5) wavelnorm = wave[160]
            glogprint (l_logfile, "mstelluric", "science", type="string",
                str="Normalizing the standard to a value of 1.0 at "//\
                str(wavelnorm)//" microns", verbose=l_verbose)
            glogprint (l_logfile, "mstelluric", "visual", type="visual",
                vistype="empty", verbose=l_verbose)

            min = 5000.0
            wheremin = 0

            for (j=1; j<=320; j+=1){
                diff = abs(wave[j]-wavelnorm)
                if (diff < min) {
                    wheremin = j 
                    min = diff 
                }
            }

            imstatistics (tmpim1//".fits["//(wheremin-2)//":"//(wheremin+2)//"]",
                field="mean", lower=INDEF, upper=INDEF, nclip=0, lsigma=INDEF,
                usigma=INDEF, format-, cache-) | scan (meanval)

            imarith (tmpim1//".fits", "/", meanval, tmpim3, ver-)
            imcopy (in[i]//"[sci,1]", tmpim5, ver-, >& "dev$null")

        } else {

            imcopy (l_cal//"[SCI,1]", tmpim1, ver-, >& "dev$null")

            glogprint (l_logfile, "msabsflux", "visual", type="visual",
                vistype="empty", verbose=l_verbose)

            tmpfile2 = mktemp ("tempfile2")
            tmpfile3 = mktemp ("tempfile3")
            wspectext (l_cal//"[SCI,1]", tmpfile2, header=no, wformat="")

            scanfile = tmpfile2

            for (j=1; j <= 320; j=j+1) {
                value[j] = 0.
                k = fscan (scanfile, val1, val2)
                if (k == 2) {
                    value[j] = real(val2)
                    wave[j] = real(val1)
                }
            }

            imstatistics (l_cal//"[sci,1]", field="min,max", lower=INDEF,
                upper=INDEF, nclip=0, lsigma=INDEF, usigma=INDEF, format-,
                cache-) | scan (ymin, ymax)
            if (ymin < 0.) {
                cutlevel2 = 2.*ymin
            } else {
                cutlevel2 = 0.
            }
            cutlevel1 = 0.01*ymax

            glogprint (l_logfile, "msabsflux", "science", type="string",
                str="Cut levels: "//str(cutlevel2)//" "//str(cutlevel1),
                verbose=l_verbose)
            for (j=1; j <= 320; j=j+1) {
                if ((value[j] < cutlevel1) && (value[j] > cutlevel2))
                    value[j] = 0.
            }

            cutlevel1 = cutlevel1 / 2.
            cutlevel2 = ymax * 2.
            glogprint (l_logfile, "msabsflux", "science", type="string",
                str="Cut levels: "//str(cutlevel1)//" "//str(cutlevel2),
                verbose=l_verbose)
            for (j=1; j <= 320; j=j+1) {
                if ((value[j] > cutlevel1) && (value[j] < cutlevel2))
                    value[j] = 1.
            }

            for (j=1; j <= 320; j=j+1) {
                print (wave[j], value[j], >> tmpfile3)
            }
            rspectext (tmpfile3, tmpim2)

            imarith (l_cal//"[1]", "*", tmpim2, tmpim3)
            if (l_fl_plots) {
                print ("     The windowed CAL spectrum:")
                splot (tmpim3)
            }

            imarith (in[i]//"[sci,1]", "*", tmpim2, tmpim5)
            if (l_fl_plots) {
                print ("     The windowed OBJECT spectrum:")
                splot (tmpim5)
            }

            imstatistics (l_cal//"[sci,1]", field="mean", lower=INDEF,
                upper=INDEF, nclip=0, lsigma=INDEF, usigma=INDEF, format-,
                cache-) | scan (meanval)
            glogprint (l_logfile, "msabsflux", "science", type="string",
                str="Cal spectrum mean value = "//str(meanval),
                verbose=l_verbose)

            imstatistics (tmpim3, field="mean", lower=INDEF, upper=INDEF,
                nclip=0, lsigma=INDEF, usigma=INDEF, format-, cache-) | \
                scan (meanval)
            glogprint (l_logfile, "msabsflux", "science", type="string",
                str="Windowed Cal spectrum mean value = "//str(meanval),
                verbose=l_verbose)

        }

        #telluric correct the spectra

        # first, add airmass keyword to 1st header extension so telluric doesn't
        # ask for it interactively.

        hselect (in[i]//"[0]", "AIRMASS", yes) | scan (airmass1)
        hselect (l_cal//"[0]", "AIRMASS", yes) | scan (airmass2)
        gemhedit (tmpim5//".fits", "AIRMASS", airmass1, "Airmass", delete-)
        gemhedit (tmpim3//".fits", "AIRMASS", airmass2, "Airmass", delete-)

        tmplog = mktemp ("tmptelluriclog")
        ymax = real(airmass2) / real(airmass1)
        if (fl_inter) {
            print (" ")
            print ("     Entering telluric.")
            print ("     To have no airmass correction, put scale = "//str(ymax))
            print (" ")
        }
        telluric (tmpim5//".fits", tmpim4, tmpim3//".fits", xcorr=l_xcorr,
            lag=l_lag, shift=l_shift, scale=l_scale, dshift=l_dshift,
            offset=1., smooth=1., cursor="", dscale=l_dscale,
            interac=l_fl_inter, thresho=l_threshold, >& tmplog)
        glogprint (l_logfile, "msabsflux", "engineering", type="file",
            str=tmplog, verbose=l_verbose)
        delete (tmplog, ver-, >& "dev$null")
        glogprint (l_logfile, "msabsflux", "visual", type="visual", 
            vistype="empty", verbose=l_verbose)

        imdelete (tmpim3, ver-, >& "dev$null")

        if (!l_fl_bbody) {
            if (l_fl_plots) {
                # now do the scaling with the SED function
                print ("     The standard star spectrum")
                splot (tmpim7)
            }

            tmpim6 = mktemp ("tmpim6")

            filename = "ratio_"//in[i]
            junk = 1
            while ((junk > 0) && (junk < 20)) {
                if (no == imaccess(filename)) {
                    glogprint (l_logfile, "msabsflux", "task", type="string",
                        str="Saving the TELLURIC output spectrum as "//filename,
                        verbose=l_verbose)
                    imcopy (tmpim4, filename, verbose-)
                    junk = 0
                } else {
                    filename = "x"//filename
                    if (junk == 20) {
                        glogprint (l_logfile, "msabsflux", "status",
                            type="warning", errno=0, str="Failed to assign a \
                            file for the TELLURIC ratio spectrum.",
                            verbose=l_verbose)
                    }
                }
            }
            imarith (tmpim4, "*", tmpim7, tmpim6)
            imarith (tmpim6, "/", meanval, tmpim6)
            hselect (tmpim4, "CRVAL1", yes) | scan (airmass1)
            glogprint (l_logfile, "msabsflux", "engineering", type="string",
                str="CRPIX1 = "//airmass1, verbose=l_verbose)
            hselect (tmpim4, "CD1_1", yes) | scan (airmass2)
            glogprint (l_logfile, "msabsflux", "engineering", type="string",
                str="CD1_1 = "//airmass2, verbose=l_verbose)
            gemhedit (tmpim7//".fits", "CRVAL1", airmass1, "", delete-)
            gemhedit (tmpim7//".fits", "CDELT1", airmass2, "", delete-)
            
            if (l_fl_plots) {
                filename = l_stdname
                junk = 1
                while ((junk > 0) && (junk < 20)) {
                    if (no == imaccess(filename)) {
                        glogprint (l_logfile, "msabsflux", "task", type="string",
                            str="Saving the standard star spectrum as "//filename,
                            verbose=l_verbose)
                        imcopy (tmpim7, filename, verbose-)
                        junk = 0
                    } else {
                        filename = "n"//filename
                        if (junk == 20) {
                            glogprint (l_logfile, "msabsflux", "status",
                                type="warning", errno=0, str="Failed to assign a \
                                file for the standard star spectrum.",
                                verbose=l_verbose)
                        }
                    }
                }
            }
        }

        if (no == l_fl_bbody && l_fl_plots) {
            print ("     The output spectrum:")
            splot (tmpim6)
        }
        tmplog = mktemp("tmpfxlog")
        fxcopy (in[i]//"[0]", out[i], >& tmplog)
        if (l_fl_bbody) {
            fxinsert (tmpim4//".fits", out[i]//"[0]", groups="", ver-)
        } else {
            fxinsert(tmpim6//".fits",out[i]//"[0]", groups="", ver-)
        }
        glogprint (l_logfile, "msabsflux", "task", type="file", str=tmplog,
            verbose=l_verbose)
        delete (tmplog, ver-, >& "dev$null")

        gemhedit (out[i]//"[1]", "EXTNAME", "SCI", "Extension name", delete-)
        gemhedit (out[i]//"[1]", "EXTVER", 1, "Extension version", delete-)

        gemdate ()
        gemhedit (out[i]//"[0]", "MSABSFLUX", gemdate.outdate,
            "UT Time stamp for MSABSFLUX", delete-)
        gemhedit (out[i]//"[0]", "GEM-TLM", gemdate.outdate,
            "UT Last modification with GEMINI", delete-)

        # jump to here if there is a problem

nextimage:
        i += 1
        ngood += 1
        imdelete (tmpim1, ver-, >& "dev$null")
        imdelete (tmpim2, ver-, >& "dev$null")
        imdelete (tmpim4, ver-, >& "dev$null")
        imdelete (tmpim5, ver-, >& "dev$null")
        if (!l_fl_bbody) {
            imdelete (tmpim6, ver-, >& "dev$null")
            imdelete (tmpim7, ver-, >& "dev$null")
            delete (tmpfile2, ver-, >& "dev$null")
            delete (tmpfile3, ver-, >& "dev$null")
        }
        delete (tmpfile1, ver-, >& "dev$null")

        glogprint (l_logfile, "msabsflux", "visual", type="visual",
            vistype="shortdash", verbose=l_verbose)

    }    

clean:
    scanfile = "" 
    delete (tmpfile, ver-, >& "dev$null")
    delete (tmpinimg, ver-, >& "dev$null")
    delete (tmpoutimg, ver-, >& "dev$null")

    if (status == 0) {
        glogprint (l_logfile, "msabsflux", "status", type="string",
            str="All "//str(nimages)//" images successfully processed.",
            verbose=l_verbose)
        glogclose (l_logfile, "msabsflux", fl_success+, verbose=l_verbose)
    } else {
        glogprint (l_logfile, "msabsflux", "status", type="string",
            str=str(ngood)//" out of "//str(nimages)//" successfully \
            processed.",verbose=l_verbose)
        glogclose (l_logfile, "msabsflux", fl_success-, verbose=l_verbose)
    }

exitnow:
    ;

end

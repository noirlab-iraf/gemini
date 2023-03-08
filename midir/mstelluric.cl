# Copyright(c) 2004-2011 Association of Universities for Research in Astronomy, Inc.

procedure mstelluric (inimages, cal, bbody)

# This routine corrects the science target spectra for telluric features.  Both the
# target spectrum and the calibrator for the correction must have been extracted
# using "nsextract".  To correct the calibrator for continuum shape, the blackbody
# temperature of the star is a necessary input value.
# 

# Version:   Jun  25, 2004   TB first write
#            Jul  26, 2004   TB finish the task, incorporate call to "telluric"
#            Oct  13, 2004   TB fix unit inconsistency in bbody division.
#            Nov  11, 2004   KV fix small bug: assumption that 11.7 microns will
#                               always be in the range of the spectrum
    
char    inimages    {prompt="Input T-ReCS or Michelle spectra"}  # OLDP-1-input-primary-single-prefix=s
char    cal         {prompt = "Telluric calibration spectrum"}
real    bbody       {prompt="Temperature of calibrator for B-body fit"}
char    outimages   {"",prompt="Output spectra"}                    # OLDP-1-output
char    outpref     {"a",prompt="Prefix for out image(s)"}          # OLDP-4
bool    xcorr       {no, prompt="Cross correlate for shift? (telluric)"}
int     lag         {10, prompt="Cross correlation lag (in pixels)? (telluric)"}
real    shift       {0., prompt="Initial shift of calibration spectrum (pixels)? (telluric)"}
real    scale       {1., prompt="Initial scale factor multiplying airmass ratio? (telluric)"}
real    dshift      {0.5, prompt="Initial shift search step? (telluric)"}
real    dscale      {0.1, prompt="Initial scale factor search step? (telluric)"}
real    threshold   {0.01, prompt="Threshold value for telluric calibration"}
bool    fl_inter    {no, prompt="Interactive tweaking of correction?"}
char    logfile     {"",prompt="Log file name"}                     # OLDP-1
bool    verbose     {yes,prompt="Verbose"}                          # OLDP-4
int     status      {0,prompt="Exit status: (0=good, >0=bad)"}      # OLDP-4
struct* scanfile    {"",prompt="Internal use only"}                 # OLDP-4

begin

    char    l_inimages, l_cal, l_logfile
    char    l_outimages = ""
    char    l_outpref = ""
    int     l_lag
    real    l_bbody, l_shift, l_scale, l_dshift, l_dscale, l_thresh
    bool    l_verbose, l_xcorr, l_fl_inter

    char    airmass1, airmass2, filename, keyfound
    char    in[100], out[100], instrument, temp, errmsg
    char    tmpfile, tmpinimg, tmpoutimg, outputstr
    char    tmpim1, tmpim2, tmpim3, tmpim4, tmpfile1
    char    paramstr, tmplog
    int     i, j, k, junk
    int     nimages, maximages, noutimages, nbad, ngood
    int     wheremin
    real    wave[320], bbflux[320], value[320], c1, c2, diff, min
    real    crval1, cdelt1, meanval, bbvalue
    real    wavelnorm

    tmpfile = mktemp("tmpfile")
    tmpinimg = mktemp("tmpinimg")
    tmpoutimg = mktemp("tmpoutimg")

    l_inimages=inimages
    l_cal=cal
    l_bbody=bbody
    junk = fscan (outimages, l_outimages)	#ensures no extra spaces
    junk = fscan (outpref, l_outpref)
    l_xcorr=xcorr
    l_lag=lag
    l_shift=shift
    l_scale=scale
    l_dshift=dshift
    l_dscale=dscale
    l_fl_inter=fl_inter
    l_logfile=logfile
    l_thresh=threshold
    l_verbose=verbose 

    # Initialize
    status=0
    ngood=0
    maximages=100

    cache ("gemextn", "gemdate")

    # Create the list of parameter/value pairs.  One pair per line.
    # All lines combined into one string.  Line delimiter is '\n'.
    paramstr =  "inimages       = "//inimages.p_value//"\n"
    paramstr += "cal            = "//cal.p_value//"\n"
    paramstr += "bbody          = "//bbody.p_value//"\n"
    paramstr += "outimages      = "//outimages.p_value//"\n"
    paramstr += "outpref        = "//outpref.p_value//"\n"
    paramstr += "xcorr          = "//xcorr.p_value//"\n"
    paramstr += "lag            = "//lag.p_value//"\n"
    paramstr += "shift          = "//shift.p_value//"\n"
    paramstr += "scale          = "//scale.p_value//"\n"
    paramstr += "dshift         = "//dshift.p_value//"\n"
    paramstr += "dscale         = "//dscale.p_value//"\n"
    paramstr += "threshold      = "//threshold.p_value//"\n"
    paramstr += "fl_inter       = "//fl_inter.p_value//"\n"
    paramstr += "logfile        = "//logfile.p_value//"\n"
    paramstr += "verbose        = "//verbose.p_value

    # Assign a logfile name if not given.  Open logfile and start log.
    # Write parameter/value pairs ("paramstr") to log.
    gloginit (l_logfile, "mstelluric", "midir", paramstr, fl_append+,
        verbose=l_verbose)
    if (gloginit.status != 0) {
        status = gloginit.status
        goto exitnow
    }
    l_logfile = gloginit.logfile.p_value

    if (l_bbody == 1.0) {
        glogprint (l_logfile, "mstelluric", "status", type="error", errno=121,
            str="ERROR - MSTELLURIC: Blackbody Temperature is not set!",
            verbose+)
        status=1
        goto clean
    }

    # Load up the array of input file names
    gemextn(l_inimages, check="", process="none", index="", extname="",
        extversion="", ikparams="", omit="extension", replace="", 
        outfile=tmpfile, logfile=l_logfile, glogpars="", verbose=l_verbose)
    gemextn("@"//tmpfile, check="exists,mef", process="none", index="",
        extname="", extversion="", ikparams="", omit="", replace="",
        outfile=tmpinimg, logfile=l_logfile, glogpars="", verbose=l_verbose)
    nimages = gemextn.count
    delete(tmpfile, ver-, >& "dev$null")

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

        glogprint (l_logfile, "mstelluric", "status", type="error", 
            errno=status, str=errmsg, verbose+)
        goto clean

    } else {
        scanfile=tmpinimg
        i=0
        while (fscan(scanfile,filename) != EOF) {
            i+=1
            in[i] = filename
        }
        scanfile=""
        if (i != nimages) {
            glogprint (l_logfile, "mstelluric", "status", type="error",
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
        gemextn("@"//tmpinimg, check="", process="none", index="", extname="",
            extversion="", ikparams="", omit="path", replace="",
            outfile=tmpoutimg, logfile=l_logfile, glogpars="",
            verbose=l_verbose)
        outputstr = l_outpref//"@"//tmpoutimg
    } else {
        glogprint (l_logfile, "mstelluric", "status", type="error", errno=121,
            str="Neither output image name nor output prefix is defined.",
            verbose+)
        status=1
        goto clean
    }

    gemextn(outputstr, check="", process="none", index="", extname="",
        extversion="", ikparams="", omit="extension", replace="",
        outfile=tmpfile, logfile=l_logfile, glogpars="", verbose=l_verbose)
    delete(tmpoutimg, ver-, >& "dev$null")
    gemextn("@"//tmpfile, check="absent", process="none", index="",
        extname="", extversion="", ikparams="", omit="", replace="",
        outfile=tmpoutimg, logfile=l_logfile, glogpars="", verbose=l_verbose)
    noutimages = gemextn.count
    delete(tmpfile, ver-, >& "dev$null")

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

        glogprint (l_logfile, "mstelluric", "status", type="error", 
            errno=status, str=errmsg, verbose+)
        goto clean

    } else {
        scanfile=tmpoutimg
        i=0
        while (fscan(scanfile,filename) != EOF) {
            i+=1
            out[i] = filename//".fits"
        }
        scanfile=""
        if (i != noutimages) {
            glogprint (l_logfile, "mstelluric", "status", type="error",
                errno=99, str="Error while counting the output images.",
                verbose+)
            status = 1
            goto clean
        }
    }

    # Make sure the calibration spectrum exists
    gemextn(l_cal, check="exists,mef", process="none", index="", extname="",
        extversion="", ikparams="", omit="", replace="", outfile="",
        logfile=l_logfile, glogpars="", verbose=l_verbose)
    if ((gemextn.fail_count > 0) || (gemextn.count == 0)) {
        if (gemextn.fail_count > 0) {
            errmsg = "Calibration file not found."
            status = 101
        } else if (gemextn.count == 0) {
            errmsg = "ERROR - MSTELLURIC: No calibration file defined."
            status = 121
        }

        glogprint (l_logfile, "mstelluric", "status", type="error",
            errno=status, str=errmsg, verbose+)
        goto clean
    }

    # Do the work
    i = 1
    while (i <= nimages) {
        glogprint (l_logfile, "mstelluric", "status", type="string",
            str="Processing "//in[i], verbose=l_verbose)

        keyfound = ""
        hselect (in[i]//"[0]", "MSTELLUR", yes) | scan(keyfound)
        if (keyfound != "") {
            glogprint (l_logfile, "mstelluric", "status", type="warning",
                errno=0, str="File "//in[i]//" has already been telluric \
                corrected.",verbose=l_verbose)
            goto nextimage
        }

        glogprint (l_logfile, "mstelluric", "task", type="string",
            str="  "//in[i]//" --> "//out[i],verbose=l_verbose)

        # check the primary FITS header
        hselect (in[i]//"[0]", "INSTRUME", yes) | scan(instrument)
        glogprint (l_logfile, "mstelluric", "engineering", type="string",
            str="Instrument is: "//instrument,verbose=l_verbose)

        keyfound = ""
        hselect (in[i]//"[0]", "NSEXTRAC", yes) | scan(keyfound)
        if (keyfound == "") {
            glogprint (l_logfile, "mstelluric", "status", type="error",
                errno=121, str="Image "//in[i]//" not EXTRACTed.",verbose+)
            status=1
            goto clean
        } else
            glogprint (l_logfile, "mstelluric", "engineering", type="string",
                str="Image "//in[i]//" has been EXTRACTed.",verbose=l_verbose)

        keyfound = ""
        hselect (l_cal//"[0]", "NSEXTRAC", yes) | scan(keyfound)
        if (keyfound == "") {
            glogprint (l_logfile, "mstelluric", "status", type="error",
                errno=121, str="Image "//l_cal//" not EXTRACTed.",verbose+)
            status=1
            goto clean
        } else
            glogprint (l_logfile, "mstelluric", "engineering", type="string",
                str="Image "//l_cal//" has been EXTRACTed.",verbose+)

        meanval=1.0
        tmpim1=mktemp("tmpim1")
        tmpim2=mktemp("tmpim2")
        tmpim3=mktemp("tmpim3")
        tmpim4=mktemp("tmpim4")
        tmpfile1=mktemp("tmpfile1")

        keyfound = ""
        hselect (in[i]//"[sci,1]", "CRVAL1", yes) | scan(keyfound)
        crval1 = real(keyfound)
        keyfound = ""
        hselect (in[i]//"[sci,1]", "CD1_1", yes) | scan(keyfound)
        cdelt1 = real(keyfound)

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
        delete (tmpfile1,ver-, >& "dev$null")
        for (j=1; j<=320; j+=1){
            bbflux[j] = wave[j]*(bbflux[j]/bbvalue)
            print(wave[j], bbflux[j], >> tmpfile1)
        }

        glogprint (l_logfile, "mstelluric", "science", type="string",
            str="Removing the blackbody shape from the standard",
            verbose=l_verbose)

        rspectext(tmpfile1,tmpim2)

        imarith(l_cal//"[sci,1]","/",tmpim2//".fits",tmpim1,ver-)

        # Set the normalization wavelength for low-resolution N-band, low-resolution 
        # Q-band or higher resolution N-band according to the range of wavelengths in 
        # the spectrum.
        #
        print("wavelength range: "//str(wave[1])//" "//str(wave[320]))
        wavelnorm=wave[320]-wave[1]
        if (wavelnorm > 3.0 && wave[1] < 10.) wavelnorm=11.7
        if (wavelnorm > 3.0 && wave[1] > 15.) wavelnorm=20.5
        if (wavelnorm < 1.5) wavelnorm=wave[160]
        glogprint (l_logfile, "mstelluric", "science", type="string",
            str="Normalizing the standard to a value of 1.0 at "//str(wavelnorm)//" microns",
            verbose=l_verbose)
        glogprint (l_logfile, "mstelluric", "visual", type="visual",
            vistype="empty", verbose=l_verbose)

        min=5000.0
        wheremin=0

        for (j=1; j<=320; j+=1){
            diff = abs(wave[j]-wavelnorm)
            if (diff < min) {
                wheremin=j 
                min=diff 
            }
        }

        imstatistics (tmpim1//".fits["//(wheremin-2)//":"//(wheremin+2)//"]",
            field="mean",lower=INDEF,upper=INDEF,nclip=0,lsigma=INDEF,
            usigma=INDEF,format-,cache-) | scan(meanval)

        imarith(tmpim1//".fits","/",meanval,tmpim3,ver-)

        #telluric correct the spectra

        # first, add airmass keyword to 1st header extension so telluric doesn't
        # ask for it interactively.

        hselect (in[i]//"[0]", "AIRMASS", yes) | scan(airmass1)
        hselect (l_cal//"[0]", "AIRMASS", yes) | scan(airmass2)
        gemhedit(in[i]//"[sci,1]","AIRMASS",airmass1,
            "Average AIRMASS for observation")
        gemhedit (tmpim3//".fits", "AIRMASS", airmass2, "", delete-)

        tmplog = mktemp("tmptelluriclog")
        telluric(in[i]//"[sci,1].fits",tmpim4,tmpim3//".fits",xcorr=l_xcorr,
            lag=l_lag,shift=l_shift,scale=l_scale,dshift=l_dshift,
            dscale=l_dscale,interac=l_fl_inter,thresho=l_thresh, >& tmplog)
        glogprint (l_logfile, "mstelluric", "science", type="file", str=tmplog, 
            verbose=l_verbose)
        glogprint (l_logfile, "mstelluric", "visual", type="visual", 
            vistype="empty", verbose=l_verbose)
        delete (tmplog, ver-, >& "dev$null")

        tmplog = mktemp("tmpfxlog")
        fxcopy(in[i]//"[0]",out[i], >& tmplog)
        fxinsert(tmpim4//".fits",out[i]//"[1]", groups="",ver-)
        glogprint (l_logfile, "mstelluric", "task", type="file", str=tmplog,
            verbose=l_verbose)
        delete (tmplog, ver-, >& "dev$null")

        gemhedit (out[i]//"[1]", "EXTNAME", "SCI", "Extension name", delete-)
        gemhedit (out[i]//"[1]", "EXTVER", 1, "Extension version", delete-)
            
        gemdate ()
        gemhedit (out[i]//"[0]", "MSTELLUR", gemdate.outdate,
            "UT Time stamp for MSTELLURIC", delete-)
        gemhedit (out[i]//"[0]", "GEM-TLM", gemdate.outdate,
            "UT Last modification with GEMINI", delete-)

        # jump to here if there is a problem

nextimage:
        i += 1
        ngood += 1
        imdelete(tmpim1,ver-,>& "dev$null")
        imdelete(tmpim2,ver-,>& "dev$null")
        imdelete(tmpim3,ver-,>& "dev$null")
        imdelete(tmpim4,ver-,>& "dev$null")
        delete(tmpfile1,ver-,>& "dev$null")

        glogprint (l_logfile, "mstelluric", "visual", type="visual",
            vistype="shortdash", verbose=l_verbose)

    }    

clean:
    scanfile="" 
    delete(tmpfile, ver-, >& "dev$null")
    delete(tmpinimg, ver-, >& "dev$null")
    delete(tmpoutimg, ver-, >& "dev$null")

    if (status==0) {
        glogprint (l_logfile, "mstelluric", "status", type="string",
            str="All "//str(nimages)//" images successfully processed.",
            verbose=l_verbose)
        glogclose (l_logfile, "mstelluric", fl_success+, verbose=l_verbose)
    } else {
        glogprint (l_logfile, "mstelluric", "status", type="string",
            str=str(ngood)//" out of "//str(nimages)//" successfully \
            processed.",verbose=l_verbose)
        glogclose (l_logfile, "mstelluric", fl_success-, verbose=l_verbose)
    }

exitnow:
    ;

end



#Copyright (c) 2006-2011 Association of Universities for Research in Astronomy, Inc.

procedure msflatcor (inimage, flat, bias)

char    inimage         {prompt="Input Michelle spectrum (raw data file)"}
char    flat            {prompt="Input flat-field image (raw data file)"}
char    bias            {prompt="Input bias image (raw data file)"}
char    outimage        {"", prompt="Output image name"}
char    outpref         {"f", prompt="Output prefix"}
char    rawpath         {"", prompt="path to raw data files"}
bool    fl_bias         {yes, prompt="Do bias subtraction?"}
bool    fl_writeflat    {no, prompt="Write out the normalized flat-field image?"}
char    normflat        {"", prompt="Output normalized flat-field image name?"}
char    logfile         {"", prompt="Logfile name"}
bool    verbose         {yes, prompt="Verbose"}
int     status          {0, prompt="Exit status: (0=good, >0=bad)"}
struct  *scanfile       {"", prompt="Internal use only"}

begin

    char    l_inimage = ""
    char    l_flat = ""
    char    l_bias = ""
    char    l_outimage = ""
    char    l_outpref = ""
    char    l_rawpath = ""
    char    l_normflat = ""
    char    l_logfile = ""
    bool    l_fl_bias, l_fl_writeflat, l_verbose
    
    char    in, out, inflat, inbias, outnorm
    char    paramstr, lastchar, errmsg, outputstr
    char    tmpinimg, tmpoutimg, tmpfile, tmpflat, tmpimage
    char    tmp1, tmp2, tmp3
    int     nimages, noutimages, ncalib
    int     junk, iext, l, lmin, lmax, k, iaxis, naxis[4]
    real    ave, a1, ratio, vmin, vmax
    
    junk = fscan (inimage, l_inimage)
    junk = fscan (flat, l_flat)
    junk = fscan (bias, l_bias)
    junk = fscan (outimage, l_outimage)
    junk = fscan (outpref, l_outpref)
    junk = fscan (rawpath, l_rawpath)
    l_fl_bias = fl_bias
    l_fl_writeflat = fl_writeflat
    junk = fscan (normflat, l_normflat)
    junk = fscan (logfile, l_logfile)
    l_verbose = verbose

    status = 0

    cache ("gemdate", "gloginit", "fparse")

    # Temporary file names
    tmpinimg = mktemp ("tmpinimg")
    tmpoutimg = mktemp ("tmpoutimg")
    tmpfile = mktemp ("tmpfile")
    tmpflat = mktemp ("tmpflatA")

    # Create the list of parameter/value pairs.  One pair per line.
    # All lines combined into one string.  Line delimiter is '\n'.
    paramstr  = "inimage        = "//inimage.p_value//"\n"
    paramstr += "flat           = "//flat.p_value//"\n"
    paramstr += "bias           = "//bias.p_value//"\n"
    paramstr += "outimage       = "//outimage.p_value//"\n"
    paramstr += "outpref        = "//outpref.p_value//"\n"
    paramstr += "rawpath        = "//rawpath.p_value//"\n"
    paramstr += "fl_bias        = "//fl_bias.p_value//"\n"
    paramstr += "fl_writeflat   = "//fl_writeflat.p_value//"\n"
    paramstr += "normflat       = "//normflat.p_value//"\n"
    paramstr += "logfile        = "//logfile.p_value//"\n"
    paramstr += "verbose        = "//verbose.p_value

    # Assign a logfile name if not given.  Open logfile and start log.
    # Write parameter/value pairs ("paramstr") to log.
    gloginit (l_logfile, "msflatcor", "midir", paramstr, fl_append+,
        verbose=l_verbose)
    if (gloginit.status != 0) {
        status = gloginit.status
        goto exitnow
    }
    l_logfile = gloginit.logfile.p_value


    # Check the rawpath name for a final / (or $)
    lastchar = substr (l_rawpath, strlen(l_rawpath), strlen(l_rawpath))
    if ((l_rawpath != "") && (lastchar != "$") && (lastchar != "/"))
        l_rawpath = l_rawpath//"/"
    
    #--------------------------------------------------
    # Check input image
    gemextn (l_inimage, check="exist,mef", process="none", index="",
        extname="", extversion="", ikparam="", omit="kernel,exten",
        replace="^%%"//l_rawpath//"%", outfile=tmpinimg, logfile=l_logfile,
        verbose=l_verbose)
    nimages = gemextn.count
    
    # Only one input allowed
    if ((gemextn.fail_count > 0) || (nimages == 0) || (nimages > 1)) {
        
        if (gemextn.fail_count > 0) {
            errmsg = gemextn.fail_count//" image was not found."
            status = 101
        } else if (nimages == 0) {
            errmsg = "No input image defined."
            status = 121
        } else if (nimages > 1) {
            errmsg = "Only one input image allowed."
            status = 121
        }
        
        glogprint (l_logfile, "msflatcor", "status", type="error",
            errno=status, str=errmsg, verbose+)
        goto exit
    } else {
        scanfile = tmpinimg
        junk = fscan (scanfile, in)
        scanfile = ""
    }
    delete (tmpinimg, ver-, >& "dev$null")

    #-------------------------------------------------------
    # Check input flat and bias image
    gemextn (l_flat, check="exist,mef", process="none", index="", extname="",
        extversion="", ikparam="", omit="kernel,exten",
        replace="^%%"//l_rawpath//"%", outfile=tmpinimg, logfile=l_logfile,
        verbose=l_verbose)
    ncalib = gemextn.count
    
    if ((gemextn.fail_count > 0) || (ncalib == 0) || (ncalib > 1)) {
        if (gemextn.fail_count > 0) {
            errmsg = gemextn.fail_count//" flat was not found."
            status = 101
        } else if (nimages == 0) {
            errmsg = "No input flat defined."
            status = 121
        } else if (nimages > 1) {
            errmsg = "Only one flat image allowed."
            status = 121
        }
        
        glogprint (l_logfile, "msflatcor", "status", type="error",
            errno=status, str=errmsg, verbose+)
        goto exit
    } else {
        scanfile = tmpinimg
        junk = fscan (scanfile, inflat)
        scanfile = ""
    }
    delete (tmpinimg, ver-, >& "dev$null")

    gemextn (l_bias, check="exist,mef", process="none", index="", extname="",
        extversion="", ikparam="", omit="kernel,exten",
        replace="^%%"//l_rawpath//"%", outfile=tmpinimg, logfile=l_logfile,
        verbose=l_verbose)
    ncalib = gemextn.count
    
    if ((gemextn.fail_count > 0) || (ncalib == 0) || (ncalib > 1)) {
        if (gemextn.fail_count > 0) {
            errmsg = gemextn.fail_count//" bias was not found."
            status = 101
        } else if (nimages == 0) {
            errmsg = "No input bias defined."
            status = 121
        } else if (nimages > 1) {
            errmsg = "Only one bias image allowed."
            status = 121
        }
        
        glogprint (l_logfile, "msflatcor", "status", type="error",
            errno=status, str=errmsg, verbose+)
        goto exit
    } else {
        scanfile = tmpinimg
        junk = fscan (scanfile, inbias)
        scanfile = ""
    }
    delete (tmpinimg, ver-, >& "dev$null")
    

    #-------------------------------------------------------
    # Check output image
    if (l_outimage != "")
        outputstr = l_outimage
    else if (l_outpref != "") {
        fparse (in, ver-)
        outputstr = l_outpref//fparse.root//fparse.extension
    } else {
        status = 121
        glogprint (l_logfile, "msflatcor", "status", type="error",
            errno=status, str="Neither output image name nor output prefix \
            defined.", verbose+)
        goto exit
    }
    gemextn (outputstr, check="", process="none", index="", extname="",
        extversion="", ikparam="", omit="extension", replace="",
        outfile=tmpfile, logfile=l_logfile, verbose=l_verbose)
    gemextn ("@"//tmpfile, check="absent", process="none", index="",
        extname="", extversion="", ikparam="", omit="", replace="",
        outfile=tmpoutimg, logfile=l_logfile, verbose=l_verbose)
    noutimages = gemextn.count
    delete (tmpfile, ver-, >& "dev$null")

    if ((gemextn.fail_count > 0) || (noutimages == 0) || \
        (noutimages != nimages)) {
        
        if (gemextn.fail_count > 0) {
            errmsg = gemextn.fail_count//" output image(s) already exist(s)."
            status = 102
        } else if (noutimages == 0) {
            errmsg = "No output images defined."
            status = 121
        } else if (noutimages != nimages) {
            errmsg = "There should be exactly one output image defined."
            status = 121
        }
        
        glogprint (l_logfile, "msflatcor", "status", type="error",
            errno=status, str=errmsg, verbose+)
        goto exit
        
    } else {
        scanfile = tmpoutimg
        junk = fscan (scanfile, out)
        scanfile = ""
    }
    delete (tmpoutimg, ver-, >& "dev$null")
    

    # If fl_writeflat is yes, make sure the output flat does not exist already.
    if (l_fl_writeflat) {
        if (l_normflat == "")   {
            outnorm = "normflat"
            glogprint (l_logfile, "msflatcor", "status", type="warning",
                errno=120, str="Setting normflat='normflat'", verbose=l_verbose)
        } else
            outnorm = l_normflat
            
        gemextn (outnorm, check="", process="none", index="", extname="",
            extversion="", ikparam="", omit="extension", replace="",
            outfile=tmpfile, logfile=l_logfile, verbose=l_verbose)
        gemextn ("@"//tmpfile, check="absent", process="none", index="",
            extname="", extversion="", ikparam="", omit="", replace="",
            outfile=tmpoutimg, logfile=l_logfile, verbose=l_verbose)
        ncalib = gemextn.count
        delete (tmpfile, ver-, >& "dev$null")
        
        if ((gemextn.fail_count > 0) || (ncalib > 1)) {
        
            if (gemextn.fail_count > 0) {
                errmsg = gemextn.fail_count//" 'normflat' already exist(s)"
                status = 102
            } else if (ncalib > 1) {
                errmsg = "Only one output normalized flat allowed."
                status = 121
            }
            
            glogprint (l_logfile, "msflatcor", "status", type="error",
                errno=status, str=errmsg, verbose+)
            goto exit
        } else {
            scanfile = tmpoutimg
            junk = fscan (scanfile, outnorm)
            scanfile = ""
        }
        delete (tmpoutimg, ver-, >& "dev$null")
    }


    # Verify NUMEXT in input image
    iext = 0
    hselect (in//"[0]", "NUMEXT", yes) | scanf ("%d", iext)
    if (iext <= 0) {
        errmsg = "The input file "//in//" has a missing or bad NUMEXT value \
            in the header ("//str(iext)//")."
        status = 132
        glogprint (l_logfile, "msflatcor", "status", type="error",
            errno=status, str=errmsg, verbose+)
        goto exit
    }


    # Get 'ave', average flux of flat field
    if (l_fl_bias) {
        imarith (inflat//"[1]", "-", inbias//"[1]", tmpflat, title="",
            divzero=0., hparams="", pixtype="", calctype="", verbose-,
            noact-, >& "dev$null")
    } else {
        imcopy (inflat//"[1]", tmpflat, ver-, >& "dev$null")
    }
    imstat (tmpflat, field="mean", format-, lower=INDEF, upper=INDEF) | \
        scanf ("%f", ave)


    lmin = 0
    lmax = 0
    for (l=1; l <= 240; l=l+1) {    # 240 is the length of the y-axis
        imstat (tmpflat//"[*,"//str(l)//":"//str(l)//"]", field="mean",
            format-, lower=INDEF, upper=INDEF) | scanf ("%f", a1)
        ratio = a1 / ave
        if ((ratio >= 0.5) && (lmin == 0)) lmin = l
        if ((ratio < 0.5) && (lmax == 0) && (lmin > 0)) lmax = l-1
    }

    if ((lmin > 0) && (lmax > 0)) {
        imstat (tmpflat//"[*,"//str(lmin)//":"//str(lmax)//"]", field="mean",
            format-, lower=INDEF, upper=INDEF) | scanf ("%f", ave) 
        imreplace (tmpflat//"[*,1:"//str(lmin)//"]", ave, lower=INDEF,
            upper=INDEF)
        imreplace (tmpflat//"[*,"//str(lmax)//":240]", ave, lower=INDEF,
            upper=INDEF)

        glogprint (l_logfile, "msflatcor", "engineering", type="string",
            str="Rows 1 to "//str(lmin)//" and "//str(lmax)//\
            " to 240 in "//l_flat//" replaced by average of inner regions \
            ("//str(ave)//").", verbose=l_verbose)
    } else {
        glogprint (l_logfile, "msflatcor", "status", type="warning",
            errno=0, str="The task could not locate any unilluminated regions \
            of the array; using all rows.", verbose=l_verbose)
    }

    imstat (tmpflat, fields="min,max", format-, lower=INDEF, upper=INDEF) | \
        scanf ("%f %f", vmin, vmax)
    if (vmin == 0) {
        # This is data with the dead channel, exclude channel in imstat
        imstat (tmpflat//"[41:320,*]", fields="min,max", format-, lower=INDEF,
            upper=INDEF) | scanf ("%f %f", vmin, vmax)
    }
    ratio = vmax / vmin

    # If the range is large, fit a surface and normalize.
    # Otherwise just divide by the mean without fitting.
    # The threshold value has to be better determined later.

    if ((ratio > 2.) || (ratio < -2.)) {
        tmpimage = mktemp ("tmpimage")
        imcopy (tmpflat, tmpimage, >& "dev$null")
        imdelete (tmpflat, yes, verify-, >& "dev$null")
        
        # new tmpflat name to avoid FITS kernel cache problem
        tmpflat = mktemp ("tmpflatB")   
        if (lmax == 0) lmax = 320   # 320 is length of x-axis
        imsurfit (tmpimage, tmpflat, xorder=4, yorder=4, type_output="response",
            function="legendre", cross_terms+, xmedian=1, ymedian=1,
            median_perce=50., lower=0., upper=0., ngrow=0, niter=0,
            regions="rows", rows=str(lmin+1)//"-"//str(lmax-1),
            columns="*", border="50", sections="", circle="", div_min=INDEF)
        imdelete (tmpimage, yes, verify-, >& "dev$null")
        imreplace (tmpflat//"[*,1:"//str(lmin)//"]", 1., lower=INDEF,
            upper=INDEF)
        imreplace (tmpflat//"[*,"//str(lmax)//":240]", 1., lower=INDEF,
            upper=INDEF)
    } else {
        if (ave <= 0.) {
            errmsg = "The flat field "//l_flat//" has a zero or negative mean."
            status = 121
            glogprint (l_logfile, "msflatcor", "status", type="error",
                errno=status, str=errmsg, verbose+)
            goto exit
        } else {
            imarith (tmpflat, "/", ave, tmpflat, ver-, >& "dev$null")
            glogprint (l_logfile, "msflatcor", "science", type="string",
                str="Normalizing the flat field frame by "//str(ave), ver-)
        }
    }

    # display(tmpflat//"[*,*,1,1]",1,zs-,zr+)
    if (l_fl_writeflat) {
        wmef (tmpflat, outnorm, extname="SCI", phu=inflat, >& "dev$null")
    }

    for (l=1; l <= iext; l=l+1) {
        tmpimage = mktemp ("tmpimage")
        imcopy (in//"["//l//"]", tmpimage, ver-, >& "dev$null")

        hselect (tmpimage, "i_naxis", yes) | scanf ("%d", iaxis)
        naxis[1] = 0
        naxis[2] = 0
        naxis[3] = 0
        naxis[4] = 0

        for (k=1; k <= iaxis; k+=1) {
            hselect (tmpimage, "naxis"//str(k), yes) | scan (naxis[k])
            # note: scanf won't work for arrays
        }

        if ((naxis[1] != 320) || (naxis[2] != 240)) {
            errmsg = "Bad input image dimensions. Exiting."
            status = 123
            glogprint (l_logfile, "msflatcor", "status", type="error",
                errno=status, str=errmsg, verbose+)
            goto exit
        }

        # Bias subtraction is not used in chop mode, where the image dimension
        # is [320,240,3,1].  For stare mode spectra the dimensions are 
        # [320,240,1,1].
        #
        # This should also work for prepared or stacked Michelle files.  The 
        # prepared files just have naxis[4]=0, and in stare mode the image is 
        # not re-arranged.  Stacked Michelle files have naxis[3] of zero and 
        # so are not bias subtracted either.

        if (l_fl_bias) {
            if (naxis[3] == 1) {
                imarith (tmpimage, "-", inbias//"[1]", tmpimage, title="",
                    divzero=0., hparams="", pixtype="", calctype="",
                    verbose-, noact-, >& "dev$null")
            }
        }
        if (naxis[3] == 1) {
            imarith (tmpimage, "/", tmpflat, tmpimage, title="", divzero=0.,
                hparams="", pixtype="", calctype="", verbose-, noact-,
                >& "dev$null")
        }
        if (naxis[3] == 3) {
            # tmp files required because sections not allowed in
            # imarith output files.
            
            tmp1 = mktemp ("tmp1")
            imarith (tmpimage//"[*,*,1,1]", "/", tmpflat, tmp1)
            tmp2 = mktemp ("tmp2")
            imarith (tmpimage//"[*,*,2,1]", "/", tmpflat, tmp2)
            tmp3 = mktemp ("tmp3")
            imarith (tmpimage//"[*,*,3,1]", "/", tmpflat, tmp3)
            imcopy (tmp1, tmpimage//"[*,*,1,1]", ver-, >& "dev$null")
            imcopy (tmp2, tmpimage//"[*,*,2,1]", ver-, >& "dev$null")
            imcopy (tmp3, tmpimage//"[*,*,3,1]", ver-, >& "dev$null")
            imdelete (tmp1, ver-, >& "dev$null")
            imdelete (tmp2, ver-, >& "dev$null")
            imdelete (tmp3, ver-, >& "dev$null")
        }
        if (l == 1) {
            wmef (tmpimage, out, extname="SCI", phu=in, >& "dev$null")
        } else {
            k = l-1
            fxinsert (tmpimage//".fits", out//".fits["//k//"]",
                groups="", ver-)
        }
        imdelete (tmpimage, ver-, >& "dev$null")
    }

    gemdate()
    gemhedit (out//"[0]", "MSFLATCO", gemdate.outdate,
        "UT Time stamp for MSFLATCOR", delete-)
    gemhedit (out//"[0]", "GEM-TLM", gemdate.outdate,
        "UT Last modification with GEMINI", delete-)

exit:
    scanfile = ""
    delete (tmpinimg, ver-, >& "dev$null")
    delete (tmpoutimg, ver-, >& "dev$null")
    delete (tmpfile, ver-, >& "dev$null")
    imdelete (tmpflat, ver-, >& "dev$null")
    imdelete (tmpimage, ver-, >& "dev$null")

    if (status == 0)
        glogclose (l_logfile, "msflatcor", fl_success+, verbose=l_verbose)
    else
        glogclose (l_logfile, "msflatcor", fl_success-, verbose=l_verbose)

exitnow:
    ;

end

# Copyright(c) 2007-2011 Association of Universities for Research in Astronomy, Inc.

procedure mipql (inimages)

# This task takes an output image from MIPSTACK, with the images for 
# the four waveplate positions, and produce the Stokes I, U, and Q images 
# (as does MIPTRANS) plus four more extenstions containting the unpolarized
# flux, the polarized flux, the % polarization, and the polarization angle 
# in degrees.  Hence the output file has 7 extensions.
#

char    inimages    {prompt="Input stacked Michelle polarimetry image(s)"}    # OLDP-1-input-primary-single-prefix=s
char    rawpath     {"", prompt="Path for in raw images"}    # OLDP-4
char    outimages   {"", prompt="Output image(s)"}    # OLDP-1-output
char    outpref     {"p", prompt="Prefix for out image(s)"}    # OLDP-4
bool    fl_mask     {yes, prompt="Mask out intensity values near zero?"}    # OLDP-4
real    threshold   {0.5, min=0.5, prompt="Percent polarization value for masking intensities"}    # OLDP-4
real    noise       {0., prompt="Masking value in ADU (used if > 0.)"}    # OLDP-4
char    blankarea   {"[*,*]", prompt="Image area used to estimate the noise level"}    # OLDP-4
bool    fl_negative {yes, prompt="mask negative pixels"}
bool    fl_register {no, prompt="Register images?"}    # OLDP-4
char    regions     {"[*,*]", prompt="Reference image region used for registration (xregister)"}    # OLDP-4
char    logfile     {"", prompt="Logfile"}    # OLDP-1
bool    verbose     {yes, prompt="Verbose?"}    # OLDP-4
int     status      {0, prompt="Exit status: (0=good, >0=bad)"}    # OLDP-4
struct  *scanfile   {"", prompt="Internal use only"}    # OLDP-4

begin
    char    l_inimages = ""
    char    l_rawpath = ""
    char    l_outimages = ""
    char    l_outpref = ""
    char    l_blankarea = ""
    char    l_regions = ""
    char    l_logfile = ""
    real    l_threshold, l_noise
    bool    l_fl_mask, l_fl_negative, l_fl_register, l_verbose

    char    paramstr, lastchar, errmsg, filename, outputstr, keyfound
    char    tmpinimg, tmpoutimg, tmpfile, tmplog
    char    tmpim1, tmpim2, tmpim3, tmpim4, tmpim5, tmpim6, tmpim7
    char    tmpmask, tmpimage
    char    in[100], out[100]
    int     junk, i, j, nimages, noutimages, maxfiles
    real    sigma, floor, maxval

    cache ("gemdate", "gloginit", "gemextn")

    junk = fscan (inimages, l_inimages)
    junk = fscan (rawpath, l_rawpath)
    junk = fscan (outimages, l_outimages)
    junk = fscan (outpref, l_outpref)
    l_fl_mask = fl_mask
    l_threshold = threshold
    l_noise = noise
    junk = fscan (blankarea, l_blankarea)
    l_fl_negative = fl_negative
    l_fl_register = fl_register
    junk = fscan (regions, l_regions)
    junk = fscan (logfile, l_logfile)
    l_verbose = verbose

    # Initialize
    maxfiles = 100
    status = 0

    # Create temp file names
    tmpinimg = mktemp ("tmpinimg")
    tmpoutimg = mktemp ("tmpoutimg")
    tmpfile = mktemp ("tmpfile")

    # Create the list of parameter/value pairs.  One pair per line.
    # All lines combined into one string.  Line delimiter is '\n'.
    paramstr =  "inimages       = "//inimages.p_value//"\n"
    paramstr += "rawpath        = "//rawpath.p_value//"\n"
    paramstr += "outimages      = "//outimages.p_value//"\n"
    paramstr += "outpref        = "//outpref.p_value//"\n"
    paramstr += "fl_mask        = "//fl_mask.p_value//"\n"
    paramstr += "noise          = "//noise.p_value//"\n"
    paramstr += "threshold      = "//threshold.p_value//"\n"
    paramstr += "blankarea      = "//blankarea.p_value//"\n"
    paramstr += "fl_negative    = "//fl_negative.p_value//"\n"
    paramstr += "fl_register    = "//fl_register.p_value//"\n"
    paramstr += "regions        = "//regions.p_value//"\n"
    paramstr += "logfile        = "//logfile.p_value//"\n"
    paramstr += "verbose        = "//verbose.p_value

    # Assign a logfile name if not given.  Open logfile and start log.
    # Write parameter/value pairs ("paramstr") to log.
    gloginit (l_logfile, "mipql", "midir", paramstr, fl_append+,
        verbose=l_verbose )
    if (gloginit.status != 0) {
        status = gloginit.status
        goto exitnow
    }
    l_logfile = gloginit.logfile.p_value

    # Add the trailing slash to rawpath, if missing
    lastchar = substr (l_rawpath, strlen(l_rawpath), strlen(l_rawpath))
    if ((l_rawpath != "") && (lastchar != "$") && (lastchar != "/"))
        l_rawpath = l_rawpath//"/"


    # Load up the array of input file names    
    gemextn (l_inimages, check="exists,mef", process="none", index="",
        extname="", extversion="", ikparam="", omit="kernel,exten",
        replace="^%%"//l_rawpath//"%", outfile=tmpinimg, logfile=l_logfile,
        verbose=l_verbose)
    nimages = gemextn.count
    
    if ((gemextn.fail_count > 0) || (nimages == 0) || \
        (nimages > maxfiles)) {
        
        if (gemextn.fail_count > 0) {
            errmsg = gemextn.fail_count//" datasets were not found."
            status = 101
        } else if (nimages == 0) {
            errmsg = "No input datasets defined."
            status = 121
        } else if (nimages > maxfiles) {
            errmsg = "Maximum number of input datasets ("//str(maxfiles)//") \
                has been exceeded."
            status = 121
        }
        
    } else {
        scanfile = tmpinimg
        i = 0
        while (fscan (scanfile, filename) != EOF) {
            i += 1
            in[i] = filename
        }
        scanfile = ""
        if (i != nimages) {
            status = 99
            errmsg = "Error while counting the input images."
        }
    }
    if (status != 0) {
        glogprint (l_logfile, "mipql", "status", type="error", errno=status,
            str=errmsg, verbose+)
        goto clean
    }
    # Keep tmpinimg.  Might need it for output names with l_outpref.


    # Load up the array of output file names
    if (l_outimages != "")
        outputstr = l_outimages
    else if (l_outpref != "") {
        gemextn ("@"//tmpinimg, check="", process="none", index="", extname="",
            extversion="", ikparams="", omit="path", replace="",
            outfile=tmpoutimg, logfile=l_logfile, verbose=l_verbose)
        outputstr = l_outpref//"@"//tmpoutimg
    } else {
        status = 121
        errmsg = "Neither output image name nor output prefix is defined."
        glogprint (l_logfile, "mipql", "status", type="error", errno=status,
            str=errmsg, verbose=yes)
        goto clean
    }
    
    gemextn (outputstr, check="", process="none", index="", extname="",
        extversion="", ikparam="", omit="extension", replace="",
        outfile=tmpfile, logfile=l_logfile, verbose=l_verbose)
    delete (tmpoutimg, verify-, >& "dev$null")
    gemextn ("@"//tmpfile, check="absent", process="none", index="",
        extname="", extversion="", ikparams="", omit="", replace="",
        outfile=tmpoutimg, logfile=l_logfile, verbose=l_verbose)
    noutimages = gemextn.count
    delete (tmpfile, verify-, >& "dev$null")
    
    if ((gemextn.fail_count > 0) || (noutimages == 0) || \
        (noutimages != nimages)) {
        
        if (gemextn.fail_count > 0) {
            errmsg = gemextn.fail_count//" dataset(s) already exist(s)."
            status = 102
        } else if (noutimages == 0) {
            errmsg = "No output dataset names defined."
            status = 121
        } else if (noutimages != nimages) {
            errmsg = "Different number of input names ("//ninmages//") and \
                output names ("//noutimages//")."
            status = 121
        }
        
    } else {
        scanfile = tmpoutimg
        i = 0
        while (fscan(scanfile, filename) != EOF) {
            i += 1
            out[i] = filename//".fits"
        }
        scanfile = ""
        if (i != noutimages) {
            errmsg = "Error while counting the output dataset names."
            status = 99
        }
    }
    if (status != 0) {
        glogprint (l_logfile, "mipql", "status", type="error", errno=status,
            str=errmsg, verbose=yes)
        goto clean
    }
    delete (tmpinimg, verify-, >& "dev$null")
    delete (tmpoutimg, verify-, >& "dev$null")

    i = 1
    while (i <= nimages) {
        keyfound = ""
        hselect (in[i]//"[0]", "MIPSTACK", yes) | scan (keyfound)
        if (keyfound == "") {
            glogprint (l_logfile, "mipql", "status", type="warning",
                errno=123, str="Dataset "//in[i]//" has NOT been stacked \
                using MIPSTACK.", verbose=l_verbose)
            goto nextimage
        }

        tmpim1 = mktemp ("tmpim1")
        tmpim2 = mktemp ("tmpim2")
        tmpim3 = mktemp ("tmpim3")
        tmpim4 = mktemp ("tmpim4")
        tmpim5 = mktemp ("tmpim5")
        tmpim6 = mktemp ("tmpim6")
        tmpim7 = mktemp ("tmpim7")
        
        if (l_fl_register == no) {
            
            # Calculate the Q stokes parameter
            #   tmpim1 goes in the 3rd extension
            
            imarith (in[i]//"[2]", "-", in[i]//"[4]", tmpim1,
                title="", divzero=0., hparams="", pixtype="", calctype="",
                verbose=no, noact=no)
            
            # Calculate the U stokes parameter
            #   tmpim2 goes in the 2nd extension
            
            imarith (in[i]//"[1]", "-", in[i]//"[3]", tmpim2,
                title="", divzero=0., hparams="", pixtype="", calctype="",
                verbose=no, noact=no)
            
            # Calculate the I stokes parameter
            #   tmpim3 goes in the first extension
  
            tmplog = mktemp("tmplog")
            imcombine (in[i]//"[1],"//in[i]//"[2],"//in[i]//"[3],"//\
                in[i]//"[4]", tmpim3, combine="sum", 
                headers="", bpmasks="", rejmasks="", nrejmasks="",
                expmasks="", sigmas="", logfile=tmplog, reject="none",
                project=no, outtype="real", outlimits="", offsets="none",
                masktype="none", maskvalue=0., blank=0., scale="none",
                zero="none", weight="none", statsec="", expname="",
                lthreshold=INDEF, hthreshold=INDEF, nlow=1, nhigh=1, nkeep=1,
                mclip=yes, lsigma=3., hsigma=3., rdnoise="0.", gain="1.",
                snoise="0.", sigscale=0.1, pclip=-0.5, grow=0.)
            delete (tmplog, verify-, >& "dev$null")                
            imarith (tmpim3, "/", 2.0, tmpim3,
                title="", divzero=0., hparams="", pixtype="", calctype="",
                verbose=no, noact=no)
            
        } else {
            images.immatch.xregister (in[i]//"[2]", in[i]//"[1]", 
                regions=l_regions, shifts=tmpfile, output=tmpim4, 
                background="none", loreject=INDEF, hireject=INDEF, apodize=0., 
                filter="none", append+, records="", correlation="discrete", 
                xwindow=11, ywindow=11, xcbox=11, ycbox=11, 
                function="centroid", interp_type="poly5", interact-, xlag=0, 
                ylag=0, dxlag=0, dylag=0)
            images.immatch.xregister (in[i]//"[3]", in[i]//"[1]", 
                regions=l_regions, shifts=tmpfile, output=tmpim5, 
                background="none", loreject=INDEF, hireject=INDEF, apodize=0., 
                filter="none", append+, records="", correlation="discrete", 
                xwindow=11, ywindow=11, xcbox=11, ycbox=11, 
                function="centroid", interp_type="poly5", interact-, xlag=0, 
                ylag=0, dxlag=0, dylag=0)
            images.immatch.xregister (in[i]//"[4]", in[i]//"[1]", 
                regions=l_regions, shifts=tmpfile, output=tmpim6, 
                background="none", loreject=INDEF, hireject=INDEF, apodize=0., 
                filter="none", append+, records="", correlation="discrete", 
                xwindow=11, ywindow=11, xcbox=11, ycbox=11, 
                function="centroid", interp_type="poly5", interact-, xlag=0, 
                ylag=0, dxlag=0, dylag=0)
            delete (tmpfile, verify-, >& "dev$null")
                
            # Calculate the Q stokes parameter
            #   tmpim1 goes in the 3rd extension
                       
            imarith (tmpim4, "-", tmpim6, tmpim1,
                title="", divzero=0., hparams="", pixtype="", calctype="",
                verbose=no, noact=no)
            
            # Calculate the U stokes parameter
            #   tmpim2 goes in the 2nd extension

            imarith (in[i]//"[1]", "-", tmpim5, tmpim2,
                title="", divzero=0., hparams="", pixtype="", calctype="",
                verbose=no, noact=no)
            
            # Calculate the I stokes parameter
            #   tmpim3 goes in the first extension
            
            tmplog = mktemp("tmplog")
            imcombine (in[i]//"[1],"//tmpim4//","//tmpim5//","//tmpim6,
                tmpim3, combine="sum",
                headers="", bpmasks="", rejmasks="", nrejmasks="",
                expmasks="", sigmas="", logfile=tmplog, reject="none",
                project=no, outtype="real", outlimits="", offsets="none",
                masktype="none", maskvalue=0., blank=0., scale="none",
                zero="none", weight="none", statsec="", expname="",
                lthreshold=INDEF, hthreshold=INDEF, nlow=1, nhigh=1, nkeep=1,
                mclip=yes, lsigma=3., hsigma=3., rdnoise="0.", gain="1.",
                snoise="0.", sigscale=0.1, pclip=-0.5, grow=0.)
            delete (tmplog, verify-, >& "dev$null")                
            imarith (tmpim3, "/", 2.0, tmpim3,
                title="", divzero=0., hparams="", pixtype="", calctype="",
                verbose=no, noact=no)

            # Clean up of the tmp files used
            imdelete (tmpim4, verify-, >& "dev$null")
            imdelete (tmpim5, verify-, >& "dev$null")
            imdelete (tmpim6, verify-, >& "dev$null")
            tmpim4 = mktemp ("tmpim4")
            tmpim5 = mktemp ("tmpim5")
            tmpim6 = mktemp ("tmpim6")
        }

        # Write the first three extensions to the output file        
        wmef (tmpim3, out[i], extname="SCI", phu=in[i], verbose=l_verbose, 
            >& "dev$null")
        fxinsert (tmpim2, out[i]//"[1]", "", verbose=l_verbose, >& "dev$null")
        fxinsert (tmpim1, out[i]//"[2]", "", verbose=l_verbose, >& "dev$null")

        # Calculate the mask, if requested
        if (l_fl_mask) {
            if (l_noise <= 0.) {
                imstat (tmpim3//l_blankarea, fields="stddev", lower=INDEF,
                    upper=INDEF, nclip=0, lsigma=3., usigma=3., binwidth=0.1,
                    format=no, cache=no) | scan (sigma)
                print (" ")
                print ("Noise value (sigma) is calculated to be "//str(sigma))
                floor = sigma * 70.5 / (l_threshold*l_threshold)
                
                # for this purpose, mask off all negative values and positive 
                # values less than the floor value; 0.5% polarization accuracy 
                # requires S/N of 282.0.

                print ("Floor value is calculated to be "//str(floor))
                imstat (tmpim3, fields="max", lower=INDEF, upper=INDEF,
                    nclip=0, lsigma=3., usigma=3., binwidth=0.1,
                    format=no, cache=no) | scan (maxval)
                if (maxval < floor) {
                    floor = sigma / 2.
                    print ("Masking value is larger than the peak value \
                        ("//str(maxval)//")\n correcting the value to "//\
                        str(floor))
                }
            } else {
                floor = l_noise
                print ("Floor value is set to be "//str(floor))
            }
            tmpmask = mktemp ("tmpmask")
            imcopy (tmpim3, tmpmask, verbose-, >& "dev$null")
            if (l_fl_negative) {
                imreplace (tmpmask, 0., imaginary=0., lower=INDEF, upper=floor,
                    radius=0.)
            } else {
                imreplace (tmpmask, 0., imaginary=0., lower=-floor,
                    upper=floor, radius=0.)
            }
            tmpimage = mktemp("tmpimage")
            imexpr ("5*int(a == 0.)", tmpimage, tmpmask, outtype="int",
                dims="auto", intype="auto", refim="auto", bwidth=0,
                btype="nearest", bpixval=0., rangecheck=yes, verbose=no,
                exprdb="none", >& "dev$null")
            imdelete (tmpmask, verify-, >& "dev$null")
            imrename (tmpimage, tmpmask, verbose-, >& "dev$null")

            imexpr ("(b > 0) ? 0 : a", tmpfile, tmpim1, tmpmask, dims="auto", 
                intype="auto", outtype="auto", refim="auto", bwidth=0, 
                btype="nearest", bpixval=0., rangecheck+, verbose+, 
                exprdb="none", >& "dev$null")
            imdelete (tmpim1, verify-, >& "dev$null")
            imcopy (tmpfile, tmpim1, verbose-, >& "dev$null")
            imdelete (tmpfile, verify-, >& "dev$null")

            imexpr ("(b > 0) ? 0 : a", tmpfile, tmpim2, tmpmask, dims="auto", 
                intype="auto", outtype="auto", refim="auto", bwidth=0, 
                btype="nearest", bpixval=0., rangecheck+, verbose+, 
                exprdb="none", >& "dev$null")
            imdelete (tmpim2, verify-, >& "dev$null")
            imcopy (tmpfile, tmpim2, verbose-, >& "dev$null")
            imdelete (tmpfile, verify-, >& "dev$null")

        }

        # Calculate the Ip component
        #   tmpim4 goes in the 5th extension
        tmpimage = mktemp("tmpimage")
        imexpr ("a**2 + b**2", tmpimage, tmpim2, tmpim1, dims="auto",
            intype="auto", outtype="auto", refim="auto", bwidth=0,
            btype="nearest", bpixval=0., rangecheck=yes, verbose=no,
            exprdb="none")
        imfunction (tmpimage, tmpim4, "sqrt", verbose=no)
        imdelete (tmpimage, verify-, >& "dev$null")
        if (l_fl_mask) {
            imexpr ("(b > 0) ? 0 : a", tmpfile, tmpim4, tmpmask, dims="auto", 
                intype="auto", outtype="auto", refim="auto", bwidth=0, 
                btype="nearest", bpixval=0., rangecheck+, verbose+, 
                exprdb="none", >& "dev$null")
            imdelete (tmpim4, verify-, >& "dev$null")
            imcopy (tmpfile, tmpim4, verbose-, >& "dev$null")
            imdelete (tmpfile, verify-, >& "dev$null")
        }

        # Calculate the Polarization Angle in degrees
        #   tmpim5 goes in the 7th extension
        
        imexpr ("deg(atan2(a,b))/2.", tmpim5, tmpim2, tmpim1, dims="auto",
            intype="auto", outtype="auto", refim="auto", bwidth=0, rangecheck+,
            verbose-, exprdb="none")
        if (l_fl_mask) {
            imexpr ("(b > 0) ? 0 : a", tmpfile, tmpim5, tmpmask, dims="auto", 
                intype="auto", outtype="auto", refim="auto", bwidth=0, 
                btype="nearest", bpixval=0., rangecheck+, verbose+, 
                exprdb="none", >& "dev$null")
            imdelete (tmpim5, verify-, >& "dev$null")
            imcopy (tmpfile, tmpim5, verbose-, >& "dev$null")
            imdelete (tmpfile, verify-, >& "dev$null")
        }

        # Calculate the Unpolarized Intensity
        #   tmpim6 goes in the 4th extension
        
        imarith (tmpim3, "-", tmpim4, tmpim6,
            title="", divzero=0., hparams="", pixtype="", calctype="",
            verbose=no, noact=no)
        if (l_fl_mask) {
            imexpr ("(b > 0) ? 0 : a", tmpfile, tmpim6, tmpmask, dims="auto", 
                intype="auto", outtype="auto", refim="auto", bwidth=0, 
                btype="nearest", bpixval=0., rangecheck+, verbose+, 
                exprdb="none", >& "dev$null")
            imdelete (tmpim6, verify-, >& "dev$null")
            imcopy (tmpfile, tmpim6, verbose-, >& "dev$null")
            imdelete (tmpfile, verify-, >& "dev$null")
        }
        
        
        # Calculate the Fractional Polarization Value in %
        #   tmpim7 goes in the 6th extension
        
        tmpimage = mktemp ("tmpimage")
        imdivide (tmpim4, tmpim3, tmpimage, constant=0., rescale="norescale",
            title="*", mean="1", verbose=no, >& "dev$null")
        imarith (tmpimage, "*", 100.0, tmpim7,
            title="", divzero=0., hparams="", pixtype="", calctype="",
            verbose=no, noact=no)
        imdelete (tmpimage, verify-, >& "dev$null")
        imreplace (tmpim7, 0., lower=100.001, upper=INDEF, radius=0.)
        imreplace (tmpim7, 0., upper=-100.001, lower=INDEF, radius=0.)
        
        
        # Write the last 4 extensions to the output file
        fxinsert (tmpim6, out[i]//"[3]", "", verbose=l_verbose, >& "dev$null")
        fxinsert (tmpim4, out[i]//"[4]", "", verbose=l_verbose, >& "dev$null")
        fxinsert (tmpim7, out[i]//"[5]", "", verbose=l_verbose, >& "dev$null")
        fxinsert (tmpim5, out[i]//"[6]", "", verbose=l_verbose, >& "dev$null")
        
        # We're done with this image, delete all the tmp files
        imdelete (tmpim1//","//tmpim2//","//tmpim3//","//tmpim4//","//\
            tmpim5//","//tmpim6//","//tmpim7, verify-, >& "dev$null")
        imdelete (tmpmask, verify-, >& "dev$null")

        # Identify the extensions
        gemhedit (out[i]//"[1]", "POLRES", "I Stokes Parameter",
            "Polarization result", delete-)
        gemhedit (out[i]//"[2]", "POLRES", "U Stokes Parameter",
            "Polarization result", delete-)
        gemhedit (out[i]//"[3]", "POLRES", "Q Stokes Parameter",
            "Polarization result", delete-)
        gemhedit (out[i]//"[4]", "POLRES", "Unpolarized Intensity",
            "Polarization result", delete-)
        gemhedit (out[i]//"[5]", "POLRES", "Ip Component",
            "Polarization result", delete-)
        gemhedit (out[i]//"[6]", "POLRES", "Fractional Pol. Value in %",
            "Polarization result", delete-)
        gemhedit (out[i]//"[7]", "POLRES", "Pol. Angle in degrees",
            "Polarization result", delete-)
        
        # EXTNAME and EXTVER
        for (j=1; j<=7; j=j+1) {
            gemhedit (out[i]//"["//j//"]", "EXTNAME", "SCI", "Extension name",
                delete-)
            gemhedit (out[i]//"["//j//"]", "EXTVER", j, "Extension version",
                delete-)
        }
        
        # Time stamps
        gemdate()
        gemhedit (out[i]//"[0]", "MIPQL", gemdate.outdate,
            "UT Time stamp for MIPQL", delete-)
        gemhedit (out[i]//"[0]", "GEM-TLM", gemdate.outdate,
            "UT Last modification with GEMINI", delete-)
            
nextimage:
        i = i+1
    }

clean:
    delete (tmpinimg, verify-, >& "dev$null")
    delete (tmpoutimg, verify-, >& "dev$null")
    
    if (status == 0)
        glogclose   (l_logfile, "mipql", fl_success+, verbose=l_verbose )
    else
        glogclose (l_logfile, "mipql", fl_success-, verbose=l_verbose )

exitnow:
    ;
    
end

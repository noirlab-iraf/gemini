# Copyright (c) 2006-2011 Association of Universities for Research in Astronomy, Inc.

procedure msdefringe (inspec)

char    inspec          {prompt="Input Michelle/T-ReCS spectrum name"}
char    outspec         {"", prompt="Output defringed Michelle/T-ReCS spectrum name"}
char    outpref         {"d", prompt="Prefix for output image"}
bool    fl_lowres       {yes, prompt="Low resolution spectrum"}
bool    fl_zerocut      {yes, prompt="Mask out negative points"}
bool    fl_interpolate  {yes, prompt="Interpolate across the masked region"}
int     fmin            {18, prompt="Start of filtering region (pixels)"}
int     fmax            {32, prompt="End of filtering region (pixels)"}
bool    fl_mef          {yes, prompt="MEF file"}
char    logfile         {"", prompt="Log file name"}
bool    verbose         {yes, prompt="Verbose?"}
int     status          {0, prompt="output exit status [0=good,1=bad]"}
struct  *scanfile       {"", prompt="internal use only"}

begin
    
    # Local variables	
    char    l_inspec = ""
    char    l_outspec = ""
    char    l_outpref =""
    char    l_logfile = ""
    bool    l_fl_lowres, l_fl_zerocut, l_fl_mef, l_fl_interpolate, l_verbose
    int     l_fmin, l_fmax

    # Other variables
    char    paramstr, errmsg, in, out, outputstr, keyfound
    char    tmpinimg, tmpoutimg, tmpfile
    int     nimages, noutimages
    int     junk
    
    char    tmp1, tmp2, tmp3, tmp4, tmp5
    char    tmpscaled1, tmpscaled2
    real    value1, value2, slope, a
    int     j, k, l, m, n, iscale, nf1, nl_fmin, nl_fmax, l1, g1, g2
    real    pos1, pos2, x1, x2, y1, y2, val, d1, d2
    char    chin, line

    # Read in parameter values  (fscan removes blank spaces)
    junk = fscan (inspec, l_inspec)
    junk = fscan (outspec, l_outspec)
    junk = fscan (outpref, l_outpref)
    l_fl_lowres = fl_lowres
    l_fl_zerocut = fl_zerocut
    l_fl_interpolate = fl_interpolate 
    l_fmin = fmin 
    l_fmax = fmax
    l_fl_mef = fl_mef
    junk = fscan (logfile, l_logfile)
    l_verbose = verbose

    cache ("gloginit", "gemdate", "fparse")

    # Initialize exit status
    status = 0

    # Temporary file names
    tmpinimg = mktemp ("tmpinimg")
    tmpoutimg = mktemp ("tmpoutimg")
    tmpfile = mktemp ("tmpfile")

    value1 = 0.0

    # Create the list of parameter/value pairs.  One pair per line.
    # All lines combined into one string.  Line delimiter is '\n'.
    paramstr =  "inspec         = "//inspec.p_value//"\n"
    paramstr += "outspec        = "//outspec.p_value//"\n"
    paramstr += "outpref        = "//outpref.p_value//"\n"
    paramstr += "fl_lowres      = "//fl_lowres.p_value//"\n"
    paramstr += "fl_zerocut     = "//fl_zerocut.p_value//"\n"
    paramstr += "fl_interpolate = "//fl_interpolate.p_value//"\n"
    paramstr += "fmin           = "//fmin.p_value//"\n"
    paramstr += "fmax           = "//fmax.p_value//"\n"
    paramstr += "fl_mef         = "//fl_mef.p_value//"\n"
    paramstr += "logfile        = "//logfile.p_value//"\n"
    paramstr += "verbose        = "//verbose.p_value

    # Assign a logfile name if not given.  Open logfile and start log.
    # Write parameter/value pairs ("paramstr") to log.
    gloginit (l_logfile, "msdefringe", "midir", paramstr, fl_append+,
        verbose=yes)
    if (gloginit.status != 0) {
        status = gloginit.status
        goto exitnow
    }
    l_logfile = gloginit.logfile.p_value

    #-----------------------------------------------
    # Check input image
    if (l_fl_mef) {
        gemextn (l_inspec, check="exist,mef", process="none", index="",
            extname="", extversion="", ikparam="", omit="kernel,exten",
            replace="", outfile=tmpinimg, logfile=l_logfile, verbose=l_verbose)
    } else {
        gemextn (l_inspec, check="exist", process="none", index="",
            extname="", extversion="", ikparam="", omit="kernel,exten",
            replace="", outfile=tmpinimg, logfile=l_logfile, verbose=l_verbose)
    }
    nimages = gemextn.count
    
    # Only one input allowed
    if ((gemextn.fail_count > 0) || (nimages == 0) || (nimages > 1)) {
    
        if (gemextn.fail_count > 0) {
            errmsg = gemextn.fail_count//" spectrum was not found."
            status = 101
        } else if (nimages == 0) {
            errmsg = "No input spectrum defined."
            status = 121
        } else if (nimages > 1) {
            errmsg = "Only one input spectrum allowed."
            status = 121
        }
        
        glogprint (l_logfile, "msdefringe", "status", type="error",
            errno=status, str=errmsg, verbose+)
        goto clean
    } else {
        scanfile = tmpinimg
        junk = fscan (scanfile, in)
        scanfile = ""
    }
    delete (tmpinimg, ver-, >& "dev$null")
    
    #-------------------------------------------------
    # Check output image
    if (l_outspec != "")
        outputstr = l_outspec
    else if (l_outpref != "") {
        fparse (in, ver-)
        outputstr = l_outpref//fparse.root//fparse.extension
    } else {
        errmsg = "Neither output spectrum name nor output prefix defined."
        status = 121
        glogprint (l_logfile, "msdefringe", "status", type="error",
            errno=status, str=errmsg, verbose+)
        goto clean
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
            errmsg = gemextn.fail_count//" output spectrum already exists."
            status = 102
        } else if (noutimages == 0) {
            errmsg = "No output spectrum defined."
            status = 121
        } else if (noutimages != nimages) {
            errmsg = "There should be exactly one output spectrum defined."
            status = 121
        }
        
        glogprint (l_logfile, "msdefringe", "status", type="error",
            errno=status, str=errmsg, verbose+)
        goto clean
    } else {
        scanfile = tmpoutimg
        junk = fscan (scanfile, out)
        scanfile = ""
    }
    delete (tmpoutimg, ver-, >& "dev$null")    


    l1 = 0
    while (l1 == 0) {

        tmp1 = mktemp ("tmpf")
        tmp4 = mktemp ("tmpcursor")
        tmp5 = mktemp ("tmpcursor")

        if (l_fl_mef) {
            imcopy (in//"[1]", tmp1, ver-,  >& "dev$null")
        } else {
            imcopy (in, tmp1, ver-, >& "dev$null")
        }

        tmp2 = mktemp ("tmpf")
        tmp3 = tmp2//"new"

        imreplace (tmp1, 0., upper=0., lower=INDEF, radius=0., imaginary=0.)

        forward (tmp1, tmp2, inreal+, inimag-, outreal+, outimag+, coord-,
            center-, inmemory+, verbose=l_verbose, 
            ftpairs="fourier$ftpairs.dat", len_blk=256, >& "dev$null")
        imdelete (tmp1, ver-, >& "dev$null")

        tmpscaled1 = mktemp("tmps")
        tmpscaled2 = mktemp("tmps")
        imcopy (tmp2//"r", tmpscaled1, ver-,  >& "dev$null")
        imcopy (tmp2//"i", tmpscaled2, ver-, >& "dev$null")
        gemhedit (tmpscaled1, "CD1_1", 1., "", delete-)
        gemhedit (tmpscaled2, "CD1_1", 1., "", delete-)

        if (l_fl_lowres) {
            l_fmin = 158
            l_fmax = 164
        }

        if (l_fl_interpolate) {
            imstat (tmp2//"r("//str(l_fmin)//":"//str(l_fmin)//")", 
                fields="mean", lower=INDEF, upper=INDEF, format-, nclip=0, 
                lsigma=3., usigma=3., binwidth=0.1, cache-) | scanf("%f", d1)
            imstat (tmp2//"r("//str(l_fmax)//":"//str(l_fmax)//")", 
                fields="mean", lower=INDEF, upper=INDEF, format-, nclip=0, 
                lsigma=3., usigma=3., binwidth=0.1, cache-) | scanf("%f", d2)
        } else {
            d1 = 0.
            d2 = 0.
        }

        printf ("0 0 1 o\n0 0 1 g\n%d %f 1 x\n%d %f 1 x\n0 0 1 q\n", l_fmin, 
            d1, l_fmax, d2, >> tmp4)
        print ("0 0 1 b\n0 0 1 o\n0 0 1 g\n0 0 1 q\n", >> tmp5)

        n = 0
        while (n == 0) {

            splot (tmpscaled2, options="xydraw", next_image=tmpscaled1,
                cursor=tmp4, band=1, units="", xmin=INDEF, xmax=INDEF,
                ymin=INDEF, ymax=INDEF, save_file="splot.log", nerrsample=0,
                sigma0=INDEF, invgain=INDEF, function="spline3", order=1,
                low_reject=2., high_reject=4., niterate=10, grow=1., markrej+,
                fnuzero=3.68e-20, mode="q")

            chin = " "
            iscale=0
            nl_fmin=0
            nl_fmax=0

            while (chin != "q") {
                line = gcur 
                l = fscanf (line, "%f %f %d %s", pos1, pos2, m, chin)
                if (chin == "l") {
                    nl_fmin = int(pos1)
                    if (nl_fmin > 161) {
                        nl_fmin = 322 - nl_fmin
                    }
                    if (nl_fmin < 2) nl_fmin = 2
                    print ("New lower limit = "//str(nl_fmin))
                }
                if (chin == "u") {
                    nl_fmax = int(pos1)
                    if (nl_fmin > 161) {
                        nl_fmax = 322 - nl_fmax
                    }
                    print ("New upper limit = "//str(nl_fmax))
                }
                if (chin == "p") {
                    print ("X = "//str(pos1)//" Y = "//str(pos2))
                }
                if ((chin == "h") || (chin == "?")) {
                    printf ("Key commands:\n  a:  autoscale plot\n  e:  set \
                        scaling\n")
                    printf ("  l:  set lower filtering frequency\n  u:  set \
                        upper filtering fdrequency\n")
                    printf ("  p:  print cursor position\n  q:  exit plot\n  \
                        h/?:  help\n")
                }
                if (chin == "a") {
                    splot (tmpscaled2, next_image=tmpscaled1, cursor=tmp4,
                        option="xydraw", band=1, units="", xmin=INDEF,
                        xmax=INDEF, ymin=INDEF, ymax=INDEF, 
                        save_file="splot.log", nerrsample=0, sigma0=INDEF, 
                        invgain=INDEF, function="spline3", order=1,
                        low_reject=2., high_reject=4., niterate=10, grow=1.,
                        markrej+, fnuzero=3.68e-20, mode="q")
                }
                if (chin == "e") {
                    if (iscale == 0) {
                        x1 = pos1
                        y1 = pos2
                        iscale = 1
                    } else {
                        iscale = 0
                        x2 = pos1
                        y2 = pos2
                        if (x1 > x2) {
                            val = x2
                            x2 = x1
                            x1 = val
                        }
                        if (y1 > y2) {
                            val = y2
                            y2 = y1
                            y1 = val
                        }
                        splot (tmpscaled2, next_image=tmpscaled1, cursor=tmp4,
                            xmin=x1, xmax=x2, ymin=y1, ymax=y2, option="xydraw",
                            band=1, units="", save_file="splot.log",
                            nerrsample=0, sigma0=INDEF, invgain=INDEF,
                            function="spline3", order=1, low_reject=2.,
                            high_reject=4., niterate=10, grow=1., markrej+,
                            fnuzero=3.68e-20, mode="q")
                    }
                }
            }
            n = 1
        }
        delete (tmp4, ver-, >& "dev$null")


        if ((nl_fmin != 0) || (nl_fmax != 0)) {
            if ((nl_fmin > 0) && (nl_fmax > 0)) {
                if ((nl_fmin > nl_fmax) && !l_fl_lowres) {
                    errmsg = "New lower bound is above the new upper bound.  \
                        The change will not be applied."
                    glogprint (l_logfile, "msdefringe", "status", type="error",
                        errno=121, str=errmsg, verbose+)
                } else {
                    l_fmin = nl_fmin
                    l_fmax = nl_fmax
                }
            } else {
                if (nl_fmin == 0)   l_fmax = nl_fmax
                else                l_fmin = nl_fmin
            }
        }

        if (l_fmin == 1)    l_fmin = 2

        if (l_fl_lowres)    l_fmax = 322 - l_fmin

        
        glogprint (l_logfile, "msdefringe", "science", type="string",
            str="Filtering range: "//str(l_fmin)//":"//str(l_fmax),
            verbose=l_verbose)
        if (l_fl_zerocut) {
            imreplace (tmp2//"r[*]", 0., upper=0., lower=INDEF, radius=0., 
                imaginary=0.)
            imreplace (tmp2//"i[*]", 0., upper=0., lower=INDEF, radius=0., 
                imaginary=0.)
        }

        if (l_fl_lowres) {
            if (l_fl_interpolate) {
                g1 = l_fmin - 1
                imstat (tmp2//"r["//str(g1)//":"//str(g1)//"]", fields="mean",
                    format-, lower=INDEF, 
                    upper=INDEF, nclip=0, lsigma=3., usigma=3., binwidth=0.1, 
                    cache-) | scanf("%f", value1)
                g2 = l_fmax + 1
                imstat (tmp2//"i["//str(g2)//":"//str(g2)//"]", fields="mean",
                    format-, lower=INDEF, 
                    upper=INDEF, nclip=0, lsigma=3., usigma=3., binwidth=0.1, 
                    cache-) | scanf("%f", value2)
            } else {
                value1 = 0.
                value2 = 0.
            }
            imreplace (tmp2//"r["//str(l_fmin)//":"//str(l_fmax)//"]", value1, 
                upper=INDEF, lower=INDEF, radius=0., imaginary=0.)
            imreplace (tmp2//"i["//str(l_fmin)//":"//str(l_fmax)//"]", value2, 
                upper=INDEF, lower=INDEF, radius=0., imaginary=0.)
        } else {
            if (l_fl_interpolate) {
                g1 = l_fmin - 1
                imstat (tmp2//"r["//str(g1)//":"//str(g1)//"]", fields="mean",
                    format-, lower=INDEF, upper=INDEF, nclip=0, lsigma=3.,
                    usigma=3., binwidth=0.1, cache-) | scanf("%f", value1)
                g2 = l_fmax + 1
                imstat (tmp2//"r["//str(g2)//":"//str(g2)//"]", fields="mean",
                    format-, lower=INDEF, upper=INDEF, nclip=0, lsigma=3.,
                    usigma=3., binwidth=0.1, cache-) | scanf("%f", value2)
                slope = (value2 - value1) / (l_fmax - l_fmin + 1)
            } else {
                value1 = 0.
                slope = 0.
            }
            for (j=l_fmin; j <= l_fmax ; j=j+1) {
                k = 322 - j
                a = value1 + (slope * (j - l_fmin + 1))
                imreplace (tmp2//"r["//str(j)//":"//str(j)//"]", a, upper=INDEF,
                    lower=INDEF, radius=0., imaginary=0.)
                imreplace (tmp2//"r["//str(k)//":"//str(k)//"]", a, upper=INDEF,
                    lower=INDEF, radius=0., imaginary=0.)
            }
            if (l_fl_interpolate) {
                imstat (tmp2//"i["//str(l_fmin-1)//":"//str(l_fmin-1)//"]",
                    fields="mean", format-, lower=INDEF, upper=INDEF, nclip=0,
                    lsigma=3., usigma=3., binwidth=0.1, cache-) | \
                    scanf("%f", value1)
                imstat (tmp2//"i["//str(l_fmax+1)//":"//str(l_fmax+1)//"]",
                    fields="mean", format-, lower=INDEF, upper=INDEF, nclip=0,
                    lsigma=3., usigma=3., binwidth=0.1, cache-) | \
                    scanf("%f", value2)
                slope = (value2 - value1) / (l_fmax - l_fmin + 1)
            } else {
                slope = 0.
                value1 = 0.
            }
            for (j=l_fmin; j <= l_fmax ; j=j+1) {
                k = 322 - j
                a = value1 + (slope * (j - l_fmin + 1))
                imreplace (tmp2//"i["//str(j)//":"//str(j)//"]", a, upper=INDEF,
                    lower=INDEF, radius=0., imaginary=0.)
                imreplace (tmp2//"i["//str(k)//":"//str(k)//"]", a, upper=INDEF,
                    lower=INDEF, radius=0., imaginary=0.)
            }
        }

        inverse (tmp2, tmp3, inreal+, inimag+, outreal+, outimag+, coord-,
            decenter+, inmemory+, verbose=l_verbose,
            ftpairs="fourier$ftpairs.dat", len_blk=256,  >& "dev$null")

        imdelete (tmp2, ver-, >& "dev$null")
        imdelete (tmp2//"r,"//tmp2//"i", ver-, >& "dev$null")
        imdelete (tmpscaled1, ver-, >& "dev$null")
        imdelete (tmpscaled2, ver-, >& "dev$null")
        
        if (l_fl_mef)   tmpscaled1 = in//"[1]"
        else            tmpscaled1 = in

        tmpscaled2 = tmp3//"r.fits"

        n = 0
        l1 = 1
        while (n == 0) {

            glogprint (l_logfile, "msdefringe", "task", type="string",
                str="images: "//tmpscaled1//" "//tmpscaled2, verbose=l_verbose)

            splot (tmpscaled1, cursor=tmp5, next_image=tmpscaled2,  band=1,
                units="", xmin=INDEF, xmax=INDEF, ymin=INDEF, ymax=INDEF, 
                save_file="splot.log", nerrsample=0, sigma0=INDEF,
                invgain=INDEF, function="spline3", order=1, low_reject=2.,
                high_reject=4., niterate=10, grow=1., markrej+,
                fnuzero=3.68e-20, mode="q")

            chin = " "
            iscale = 0
            while (chin != "q") {
                line = gcur 
                l = fscanf (line, "%f %f %d %s", pos1, pos2, m, chin)
                if (chin == "p") {
                    print ("X = "//str(pos1)//" Y = "//str(pos2))
                }
                if ((chin == "h") || (chin == "?")) {
                    printf ("Key commands:\n  a:  autoscale plot\n  e:  set \
                        scaling\n")
                    printf ("  p:  print cursor position\n  r:  reset and \
                        execute the filtering again\n")
                    printf ("  i:  exit without writing the output file\n  q:  \
                        exit plot\n  h/?:  help\n")
                }
                if (chin == "a") {
                    splot (tmpscaled1, next_image=tmpscaled2, cursor=tmp5,
                        band=1, units="", xmin=INDEF, xmax=INDEF, ymin=INDEF,
                        ymax=INDEF, save_file="splot.log", nerrsample=0,
                        sigma0=INDEF, invgain=INDEF, function="spline3",
                        order=1, low_reject=2., high_reject=4., niterate=10,
                        grow=1., markrej+, fnuzero=3.68e-20, mode="q")
                }
                if (chin == "r") {
                    print ("Trying again...")
                    chin = "q"
                    l1 = 0
                }
                if (chin == "i") {
                    glogprint (l_logfile, "msdefringe", "status", type="string",
                        str="Exiting without writing the output file",
                        verbose=l_verbose)
                    imdelete (out, verify-,  >& "dev$null")
                    goto clean      
                }
                if (chin == "e") {
                    if (iscale == 0) {
                        x1 = pos1
                        y1 = pos2
                        iscale = 1
                    } else {
                        iscale = 0
                        x2 = pos1
                        y2 = pos2
                        if (x1 > x2) {
                            val = x2
                            x2 = x1
                            x1 = val
                        }
                        if (y1 > y2) {
                            val = y2
                            y2 = y1
                            y1 = val
                        }
                        splot (tmpscaled1, next_image=tmpscaled2, cursor=tmp5,
                            xmin=x1, xmax=x2, ymin=y1, ymax=y2, band=1,
                            units="", save_file="splot.log", nerrsample=0,
                            sigma0=INDEF, invgain=INDEF, function="spline3",
                            order=1, low_reject=2., high_reject=4.,
                            niterate=10, grow=1., markrej+, fnuzero=3.68e-20,
                            mode="q")
                    }
                }
            }
            n = 1
        }
        delete (tmp5, ver-, >& "dev$null")

        if (l1 != 0) {
            keyfound = ""
            if (l_fl_mef) {
                hselect (in//"[0]", "OBJECT", yes) | scan (keyfound)
            } else {
                hselect (in, "OBJECT", yes) | scan (keyfound)
            }
            if (keyfound != "") {
                gemhedit (tmp3//"r.fits", "OBJECT", keyfound, "", delete-)
                gemhedit (tmp3//"i.fits", "OBJECT", keyfound, "", delete-)
            }
            if (l_fl_mef) {
                wmef (tmp3//"r.fits", out, extname="SCI", phu=in, >& "dev$null")
            } else {
                imcopy (tmp3//"r.fits", out//".fits", ver-, >& "dev$null")
            }
        }
        
        imdelete (tmp3, ver-, >& "dev$null")
        imdelete (tmp3//"r,"//tmp3//"i", ver-, >& "dev$null")

    }

    gemdate ()
    if (l_fl_mef) {
        gemhedit (out//"[0]",  "MSDEFRINGE",  gemdate.outdate, 
            "UT Time stamp for MSDEFRINGE",  delete-)
        gemhedit (out//"[0]",  "GEM-TLM",  gemdate.outdate, 
            "UT Last modification with GEMINI",  delete-)
    }
    else {
        gemhedit (out,  "MSDEFRINGE",  gemdate.outdate, 
            "UT Time stamp for MSDEFRINGE",  delete-)
        gemhedit (out,  "GEM-TLM",  gemdate.outdate, 
            "UT Last modification with GEMINI",  delete-)
    }

    status = 0
    
clean:
    scanfile = ""
    delete (tmpinimg, ver-, >& "dev$null")
    delete (tmpoutimg, ver-, >& "dev$null")
    delete (tmpfile, ver-, >& "dev$null")
    
    delete (tmp5, ver-, >& "dev$null")
    imdelete (tmp3, ver-, >& "dev$null")
    imdelete (tmp3//"r,"//tmp3//"i", ver-, >& "dev$null")

    if (status == 0)
        glogclose (l_logfile, "msdefringe", fl_success+, verbose=l_verbose)
    else
        glogclose (l_logfile, "msdefringe", fl_success-, verbose=l_verbose)   

exitnow:
    ;

end

# Copyright(c) 2006-2011 Association of Universities for Research in Astronomy, Inc.

procedure mipsstk (inimage)

# This task stacks up the output files from "mipstokes", with or without 
# registration.  The resulting output file has the same structure as the 
# output files from "mipql" when "fl_variance" is "no".  Otherwise there 
# are three extra extensions with the variances for I, U, and Q respectively.
#
# The registration option is crude at the momemt.
#

char    inimage     {prompt="Input Michelle Stokes polarimetry images"} # OLDP-1-input-primary-single-prefix=s
char    outimage    {"", prompt="Output images"}                        # OLDP-1-output
char    outpref     {"a", prompt="Prefix for output images"}            # OLDP-4
char    rawpath     {"", prompt="Path for input raw images"}            # OLDP-4
bool    fl_register {no,prompt="Register images when combining"}        # OLDP-2
char    regions     {"[*,*]",prompt="Reference image regions used for registration (xregister)"}    # OLDP-4
bool    fl_variance {no,prompt="Output variance images"}                # OLDP-4
bool    fl_stair    {yes, prompt="Correct channel offsets"}             # OLDP-4
bool    fl_mask     {yes,prompt="Mask out intensity values near zero?"} # OLDP-4
real    noise       {0.0,prompt="Masking value in ADU (used if > 0.0)"} # OLDP-4
real    threshold   {0.5,prompt="Percent polarization value for masking intensities"}
char    blankarea   {"[*,*]",prompt="Image area used to estimate the noise level"}    # OLDP-4
char    logfile     {"", prompt="Logfile"}                              # OLDP-1
bool    verbose     {yes, prompt="Verbose"}                             # OLDP-4
int     status      {0, prompt="Exit status: (0=good, >0=bad)"}         # OLDP-4
struct  *scanfile   {"", prompt="Internal use only"}                    # OLDP-4

begin

    char    l_inimage = ""
    char    l_outimage = ""
    char    l_outpref = ""
    char    l_rawpath = ""
    char    l_regions = ""
    char    l_blankarea = ""
    char    l_logfile = ""
    real    l_noise, l_threshold
    bool    l_fl_register, l_fl_variance, l_fl_stair, l_fl_mask, l_verbose
    
    char    paramstr, lastchar
    int     junk    
    char    imagename
    char    imagelist, refimage, tmpimage, tmpwork[3], sigin[3], newimagelist
    char    combimagelist
    char    tmpshift, tmpfinal, tmpshift1, tmpmask, tmplog, tmpfile
    char    tmpshifted, instr, instr1, outstr
    char    tmp1, tmp2, tmp3, tmp4, midval
    int     nexts, i, j, m, nval
    real    chlevel, floor, sigma

    cache ("gemdate")

    junk = fscan (inimage, l_inimage)
    junk = fscan (outimage, l_outimage)
    junk = fscan (outpref, l_outpref)
    junk = fscan (rawpath, l_rawpath)
    l_fl_register = fl_register
    junk = fscan (regions, l_regions)
    l_fl_variance = fl_variance
    l_fl_stair = fl_stair
    l_fl_mask = fl_mask
    l_noise = noise
    l_threshold = threshold
    junk = fscan (blankarea, l_blankarea)
    junk = fscan (logfile, l_logfile)
    l_verbose = verbose
    
    # Initialize
    status = 0

    # Create the list of parameter/value pairs.  One pair per line.
    # All lines combined into one string.  Line delimiter is '\n'.
    paramstr =  "inimage        = "//inimage.p_value//"\n"
    paramstr += "outimage       = "//outimage.p_value//"\n"
    paramstr += "outpref        = "//outpref.p_value//"\n"
    paramstr += "rawpath        = "//rawpath.p_value//"\n"
    paramstr += "fl_register    = "//fl_register.p_value//"\n"
    paramstr += "regions        = "//regions.p_value//"\n"
    paramstr += "fl_variance    = "//fl_variance.p_value//"\n"
    paramstr += "fl_stair       = "//fl_stair.p_value//"\n"
    paramstr += "fl_mask        = "//fl_mask.p_value//"\n"
    paramstr += "noise          = "//noise.p_value//"\n"
    paramstr += "threshold      = "//threshold.p_value//"\n"
    paramstr += "blankarea      = "//blankarea.p_value//"\n"
    paramstr += "logfile        = "//logfile.p_value//"\n"
    paramstr += "verbose        = "//verbose.p_value

    # Assign a logfile name if not given.  Open logfile and start log.
    # Write parameter/value pairs ("paramstr") to log.
    gloginit (l_logfile, "mipsstk", "midir", paramstr, fl_append+,
        verbose=l_verbose )
    if (gloginit.status != 0) {
        status = gloginit.status
        goto exitnow
    }
    l_logfile = gloginit.logfile.p_value

    #----------------------------------------------------------------------
    # Check the rawpath name for a final /
    #
    # Is rawpath really required?  Aren't the inputs supposed to have
    # been mprepare'd, and should therefore no longer be 'raw' and shouldn't
    # they be in the current directory?
    # -- KL

    lastchar = substr (l_rawpath, strlen (l_rawpath), strlen(l_rawpath))
    if ( (l_rawpath != "") && (lastchar != "$") && (lastchar != "/") )
        l_rawpath = l_rawpath//"/"
    

    if (substr(l_inimage,strlen(l_inimage)-4,strlen(l_inimage)) == ".fits") {
        l_inimage = substr (l_inimage, 1, strlen(l_inimage)-5)
    }

    if (l_outpref == "") l_outpref = "s"

    if ((l_outpref != "") && ((l_outimage == "") || (l_outimage == " "))) {
        l_outimage = l_outpref//l_inimage
    }

    if (substr(l_outimage,strlen(l_outimage)-4,strlen(l_outimage)) == ".fits"){
        l_inimage = substr (l_outimage, 1, strlen(l_outimage)-5)
    }

    if (imaccess(l_outimage)) {
        glogprint (l_logfile, "mipsstk", "status", type="error", errno=121,
            str="The output image already exists.")
        status = 1
        goto exit
    }

    imagename = l_inimage
    l_inimage = l_rawpath//l_inimage

    if (no == imaccess(l_inimage)) {
        glogprint (l_logfile, "mipsstk", "status", type="error", errno=121,
            str="The input image does not exist.")
        status = 1
        goto exit
    }

    i = 1
    while (imaccess(l_inimage//"["//i//"]")) {
        i = i+1
    }
    nexts = i-1

    if (nexts < 14) {
        glogprint (l_logfile, "mipsstk", "status", type="error", errno=121,
            str="Error: number of extensions is too small: "//nexts//".")
        status = 1
        goto exit
    }

    tmpfinal = mktemp ("tmpwork")
    tmpwork[1] = mktemp ("tmpwork")
    tmpwork[2] = mktemp ("tmpwork")
    tmpwork[3] = mktemp ("tmpwork")
    tmplog = mktemp ("tmplog")
    tmpfile = mktemp ("tmpfile")
    sigin[1] = mktemp ("tmpsigma")
    sigin[2] = mktemp ("tmpsigma")
    sigin[3] = mktemp ("tmpsigma")

    imagelist = mktemp ("tmpimagelist")
    newimagelist = mktemp ("tmpimagelist")

    for (i=1; i <= nexts; i=i+7) {
        if (l_fl_register) {
            if (i == 1) {
                tmpshift = mktemp("tmpshift")
                tmpshift1 = mktemp("tmpshift")
                refimage = mktemp("tmprefimage")
                #imcopy (l_inimage//"[1][*,*]", refimage, verbose-, 
                #    >& "dev$null")
                imcopy (l_inimage//"[1]", refimage, verbose-, >& "dev$null")
                print (refimage, > imagelist)
            }
            tmpimage = mktemp ("tmpimage")
            imcopy (l_inimage//"["//i//"]", tmpwork[1], verbose-, 
                >& "dev$null")
            images.immatch.xregister (tmpwork[1], refimage, regions=l_regions, 
                shifts=tmpshift, output=tmpimage, background="none", 
                loreject=INDEF, hireject=INDEF, apodize=0., filter="none", 
                append+, records="", correlation="discrete", xwindow=11,
                ywindow=11, function="centroid", interp_type="poly5", 
                interact-, xlag=0, ylag=0, dxlag=0, dylag=0, xcbox=5, ycbox=5, 
                >& tmplog)
            glogprint (l_logfile, "mipsstk", "engineering", type="file",
                str=tmplog, verbose=l_verbose)
            delete (tmplog, verify-, >& "dev$null")
            print (tmpimage, >> imagelist)
            imdelete (tmpwork[1], yes, verify-, >& "dev$null")
        } else {
            tmpimage = mktemp("tmpimage")
            imcopy (l_inimage//"["//i//"]", tmpimage, verbose-, >& "dev$null")
            print (tmpimage, >> imagelist)      
        }
    }

    imcombine ("@"//imagelist, tmpwork[1], combine="average", headers="",
        bpmasks="", rejmasks="", expmask="", sigmas=sigin[1], logfile="STDOUT",
        reject="none", project-, outtype="double", outlimits="", 
        offsets="none", masktype="none", maskvalue=0., blank=0., scale="none",
        zero="none", weight="none", statsec="", expname="", lthreshold=INDEF,
        hthreshold=INDEF, nlow=1, nhigh=1, mclip+, lsigma=3., hsigma=3.,
        rdnoise="0.", gain="1.", snoise="0.", sigscale=0.1, pclip= -0.5,
        grow=0., >& "dev$null")

    for (j=2; j <= 3; j=j+1) {
        delete (imagelist, verify-, >& "dev$null")
        delete (newimagelist, verify-, >& "dev$null")

        for (i=1; i <= nexts; i=i+7) {
            print (l_inimage//"["//str(i+j-1)//"]", >> imagelist)
            if (l_fl_register) {
                tmpshifted = mktemp("tmpshifted")
                print (tmpshifted, >> newimagelist)
            }
        }

        combimagelist = imagelist

        if (l_fl_register) {
            scanfile = tmpshift
            instr = ""
            instr1 = ""
            outstr = ""
            while (fscan(scanfile, instr, instr1) != EOF) {
                if (substr(instr,1,6) == "xshift") {
                    outstr = instr1
                }
                if ((substr(instr,1,6) == "yshift") && (strlen(outstr) > 1)) {
                    outstr = outstr//" "//instr1//" x"
                    outstr = substr (outstr, 1, strldx("x",outstr)-1)
                    print (outstr, >> tmpshift1)
                    outstr = ""
                }
            }
            imshift ("@"//imagelist, "@"//newimagelist, 0., 0.,
                shifts_file=tmpshift1, interp_type="linear",
                boundary_type="constant", constant=0., >& "dev$null")
            combimagelist = newimagelist
            delete (tmpshift1, verify-, >& "dev$null")
        }

        imcombine ("@"//combimagelist, tmpwork[j], combine="average",
            headers="", bpmasks="", rejmasks="", expmask="", sigmas=sigin[j],
            logfile="STDOUT", reject="none", project-, outtype="double",
            outlimits="", offsets="none", masktype="none", maskvalue=0.,
            blank=0., scale="none", zero="none", weight="none", statsec="",
            expname="", lthreshold=INDEF, hthreshold=INDEF, nlow=1, nhigh=1,
            mclip+, lsigma=3., hsigma=3., rdnoise="0.", gain="1.", snoise="0.",
            sigscale=0.1, pclip= -0.5, grow=0., >& "dev$null")
    }


    tmp1 = mktemp ("tmpwork")
    if (l_fl_stair) {
        # In the following, check for midpt => INDEF, which happens when all 
        # the values are the same...although that should mean midpt is the same
        # as any one of the values.  This can happen if the entire region is 
        # blanked off by the threshold.
        #
        for (m=1; m <= 301; m=m+20) {
            chlevel = 0.
            imstatistics (tmpwork[1]//"["//str(m)//":"//str(m+19)//",1:240]",
                fields="midpt", lower=INDEF, upper=INDEF, format-, nclip=0,
                lsigma=3., usigma=3., binwidth=0.1, cache-) | scanf("%s", 
                midval)
            if (stridx(midval,"INDEF") == 0) {
                nval = fscanf (midval, "%f", chlevel)
                imarith (tmpwork[1]//"["//str(m)//":"//str(m+19)//",1:240]",
                    "-", chlevel, tmp1)
                imcopy (tmp1,
                    tmpwork[1]//"["//str(m)//":"//str(m+19)//",1:240]", 
                    verbose-, >& "dev$null")
                imdelete (tmp1, verify-, >& "dev$null")
            }
        }
        for (m=1; m <= 301; m=m+20) {
            chlevel = 0.
            imstatistics (tmpwork[2]//"["//str(m)//":"//str(m+19)//",1:240]",
                fields="midpt", lower=INDEF, upper=INDEF, format-, nclip=0,
                lsigma=3., usigma=3., binwidth=0.1, cache-) | scanf("%s", 
                midval)
            if (stridx(midval,"INDEF") == 0) {
                nval = fscanf (midval, "%f", chlevel)
                imarith (tmpwork[2]//"["//str(m)//":"//str(m+19)//",1:240]",
                    "-", chlevel, tmp1)
                imcopy (tmp1,
                    tmpwork[2]//"["//str(m)//":"//str(m+19)//",1:240]", 
                    verbose-, >& "dev$null")
                imdelete (tmp1, verify-, >& "dev$null")
            }
        }
        for (m=1; m <= 301; m=m+20) {
            chlevel = 0.
            imstatistics (tmpwork[3]//"["//str(m)//":"//str(m+19)//",1:240]",
                fields="midpt", lower=INDEF, upper=INDEF, format-, nclip=0,
                lsigma=3., usigma=3., binwidth=0.1, cache-) | scanf("%s", 
                midval)
            if (stridx(midval,"INDEF") == 0) {
                nval = fscanf (midval, "%f", chlevel)
                imarith (tmpwork[3]//"["//str(m)//":"//str(m+19)//",1:240]",
                    "-", chlevel, tmp1)
                imcopy (tmp1, 
                    tmpwork[3]//"["//str(m)//":"//str(m+19)//",1:240]", 
                    verbose-, >& "dev$null")
                imdelete (tmp1, verify-, >& "dev$null")
            }
        }
    }

    wmef (tmpwork[1], l_outimage, extname="SCI", phu=l_inimage//"[0]",
        verbose-, >& "dev$null")

    if (l_fl_mask) {
        if (l_noise <= 0.) {
            imstat (tmpwork[1]//l_blankarea, fields="stddev", form-, nclip=0,
                upper=INDEF, lower=INDEF) | scanf ("%f",sigma)
            glogprint (l_logfile, "mipsstk", "visual", type="visual",
                vistype="shortdash", verbose=l_verbose)
            glogprint (l_logfile, "mipsstk", "visual", type="visual",
                vistype="empty", verbose=l_verbose)
            glogprint (l_logfile, "mipsstk", "engineering", type="string",
                str="Noise value (sigma) is calculated to be "//str(sigma),
                verbose=l_verbose)
            floor = sigma * 282.0 / (4.*l_threshold*l_threshold)
            
            # S/N required is 282 for 0.5% accuracy in the final image, goes as
            # the inverse square of the accuracy required.
            
            glogprint (l_logfile, "mipsstk", "engineering", type="string",
                str="Masking value is calculated as "//str(floor),
                verbose=l_verbose)
            imstat (tmpwork[1], fields="max", form-, nclip=0, upper=INDEF,
                lower=INDEF) | scanf ("%f",sigma)
            if (sigma < floor) {
                floor = sigma/2.
                glogprint (l_logfile, "mipsstk", "engineering", type="string",
                    str="Masking value is larger than the peak value \
                    ("//str(sigma)//") correcting the value to "//str(floor),
                    verbose=l_verbose)
            }
            glogprint (l_logfile, "mipsstk", "visual", type="visual",
                vistype="empty", verbose=l_verbose)
        } else {
            glogprint (l_logfile, "mipsstk", "visual", type="visual",
                vistype="shortdash", verbose=l_verbose)
            glogprint (l_logfile, "mipsstk", "visual", type="visual",
                vistype="empty", verbose=l_verbose)
            glogprint (l_logfile, "mipsstk", "engineering", type="string",
                str="Masking value is set to be "//str(l_noise),
                verbose=l_verbose)
            glogprint (l_logfile, "mipsstk", "visual", type="visual",
                vistype="empty", verbose=l_verbose)
            floor = l_noise
        }
        
        # for this purpose, mask off all pixels with values less than the 
        # floor value in the Stokes I image, then use this as a mask for the 
        # Stokes U and Q images.  This masks out the negative beams.
        
        tmpmask = mktemp ("tmpmask")
        imcopy (tmpwork[1], tmpmask, verbose-, >& "dev$null")
        imreplace (tmpmask, 0., radius=0., upper=floor, lower=INDEF,
            >& "dev$null")
        imexpr ("5*int(a == 0.)", tmp1, tmpmask, outtype="int", >& "dev$null")
        imdelete (tmpmask, verify-, >& "dev$null")
        imrename (tmp1, tmpmask, verbose-, >& "dev$null")

        imexpr ("(b > 0) ? 0 : a", tmpfile, tmpwork[2], tmpmask, dims="auto", 
            intype="auto", outtype="auto", refim="auto", bwidth=0, 
            btype="nearest", bpixval=0., rangecheck+, verbose+, exprdb="none", 
            >& "dev$null")
        imdelete (tmpwork[2], verify-, >& "dev$null")
        imcopy (tmpfile, tmpwork[2], verbose-, >& "dev$null")
        imdelete (tmpfile, verify-, >& "dev$null")

        imexpr ("(b > 0) ? 0 : a", tmpfile, tmpwork[3], tmpmask, dims="auto", 
            intype="auto", outtype="auto", refim="auto", bwidth=0, 
            btype="nearest", bpixval=0., rangecheck+, verbose+, exprdb="none", 
            >& "dev$null")
        imdelete (tmpwork[3], verify-, >& "dev$null")
        imcopy (tmpfile, tmpwork[3], verbose-, >& "dev$null")
        imdelete (tmpfile, verify-, >& "dev$null")
        imdelete (tmpmask, verify-, >& "dev$null")
    }

    fxinsert (tmpwork[2], l_outimage//".fits[1]", "", >& "dev$null")
    fxinsert (tmpwork[3], l_outimage//".fits[2]", "", >& "dev$null")
    
    tmp2 = mktemp ("tmpwork")
    tmp3 = mktemp ("tmpwork")
    tmp4 = mktemp ("tmpwork")

    imcopy (tmpwork[2], tmp1, >& "dev$null")
    imcopy (tmpwork[3], tmp2, >& "dev$null")
    imfunction (tmp1, tmp1, "square", >& "dev$null")
    imfunction (tmp2, tmp2, "square", >& "dev$null")
    imarith (tmp1, "+", tmp2, tmp1, >& "dev$null")
    imfunction (tmp1, tmp1, "sqrt", >& "dev$null")
    # tmp1 is the Ip component
    imdelete (tmp2, verify-, >& "dev$null")
    tmp2 = mktemp ("tmp2")
    imexpr ("deg(atan2(a,b))/2.", tmp2, tmpwork[2], tmpwork[3], dims="auto",
        intype="auto", outtype="auto", refim="auto", bwidth=0, rangecheck+,
        verbose-, exprdb="none")
    # tmp2 is the polarization angle in degrees
    imarith (tmpwork[1], "-", tmp1, tmp3)
    # tmp3 is the unpolarized intensity
    imdivide (tmp1, tmpwork[1], tmp4, constant=0., rescale="norescale",
        >& "dev$null")
    imarith (tmp4, "*", 100.0, tmp4)
    # mask out values over 100% or under -100% as obviously due to noise
    imreplace (tmp4,0., lower=100.001, upper=INDEF, radius=0.)
    imreplace (tmp4,0., upper=-100.001, lower=INDEF, radius=0.)
    # tmp4 is the fractional polarization value in %

    fxinsert (tmp3, l_outimage//".fits[3]", "", >& "dev$null")
    fxinsert (tmp1, l_outimage//".fits[4]", "", >& "dev$null")
    fxinsert (tmp4, l_outimage//".fits[5]", "", >& "dev$null")
    fxinsert (tmp2, l_outimage//".fits[6]", "", >& "dev$null")

    for (i=1; i<=7; i+=1) {
        gemhedit (l_outimage//".fits["//str(i)//"]", "EXTNAME", "SCI",
            "Extension name", delete-)
        gemhedit (l_outimage//".fits["//str(i)//"]", "EXTVER", i,
            "Extension version", delete-)
    }
    
    if (l_fl_variance) {
        imfunction (sigin[1], sigin[1], "square", verbose-, >& "dev$null")
        fxinsert (sigin[1], l_outimage//".fits[7]", "", >& "dev$null")
        imfunction (sigin[2], sigin[2], "square", verbose-, >& "dev$null")
        fxinsert (sigin[2], l_outimage//".fits[8]", "", >& "dev$null")
        imfunction (sigin[3], sigin[3], "square", verbose-, >& "dev$null")
        fxinsert (sigin[3], l_outimage//".fits[9]", "", >& "dev$null")
        for (i=1; i<=3; i+=1) {
            gemhedit (l_outimage//".fits["//str(i+7)//"]", "EXTNAME", "VAR",
                "Extension name", delete-)
            gemhedit (l_outimage//".fits["//str(i+7)//"]", "EXTVER", i,
                "Extension version", delete-)
        }
        
    }
    
    gemdate()
    gemhedit (l_outimage//".fits[0]", "MIPSSTK", gemdate.outdate,
        "UT Time stamp for MIPSSTK", delete-)
    gemhedit (l_outimage//".fits[0]", "GEM-TLM", gemdate.outdate,
        "UT Last modification with GEMINI", delete-)

    imdelete (tmpwork[1], verify-, >& "dev$null")
    imdelete (tmpwork[2], verify-, >& "dev$null")
    imdelete (tmpwork[3], verify-, >& "dev$null")
    imdelete (tmp1, verify-, >& "dev$null")
    imdelete (tmp2, verify-, >& "dev$null")
    imdelete (tmp3, verify-, >& "dev$null")
    imdelete (tmp4, verify-, >& "dev$null")

    delete ("tmpshift*", verify-, >& "dev$null")
    delete ("tmpimage*", verify-, >& "dev$null")
    delete ("tmpsigma*", verify-, >& "dev$null")
    delete ("tmprefimage*", verify-, >& "dev$null")
    delete ("tmpimagelist*", verify-, >& "dev$null")

exit:

    if (status==0)
        glogclose (l_logfile, "mipsstk", fl_success+, verbose=l_verbose)
    else
        glogclose (l_logfile, "mipsstk", fl_success-, verbose=l_verbose)

exitnow:
    ;

end

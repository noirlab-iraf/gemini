# Copyright(c) 2005-2011 Association of Universities for Research in Astronomy, Inc.

procedure nvnoise(inimages)

# Subtract a constant from each of the 32 amplifier outputs
# to remove the vertical stripe noise that sometimes occurs
# in NIRI and GNIRS data.
#
# No correction for bad pixels is done; adding a BPM flag would
# be a good addition.
#
# This script runs on RAW pre-*prepared images only.  This way 
# the VAR plane can be generated properly after noise-correction.
#
# Currently only tested and useful for GNIRS data where most of
# the array is unilluminated.
# 
# Version  Feb. 14, 2005 JJ, first working version for GNIRS
#          Feb. 15, 2005 JJ, file index bug fixed
#          Feb. 17, 2005 JJ, before or after *PREPARE, variance fixed
#          Apr. 27, 2005 BR, added descriptive line to log output
#          May 2, 2005   JJ, added statsec and bad pixel fixing
#          Jun 2, 2005   JJ, try to match across quadrant boundaries better 
#                           (for NIRI)
#          Jun 30, 2006  JJ, fixed no BPM bug, fixed index problem


char    inimages    {prompt="Input NIRI or GNIRS image(s)"}             # OLDP-1-input-primary-single-prefix=v
char    outimages   {"", prompt="Output image(s)"}                      # OLDP-1-output
char    outprefix   {"v", prompt="Prefix for output image(s)"}          # OLDP-4
char    statsec     {"[20:1000,20:1000]",prompt="Statistics section"}   # OLDP-2
char    stattype    {"mode", prompt="Statistics type (mean,midpt,mode)"}    # OLDP-2
real    sigma       {3., prompt="Rejection limit in sigma for IMSTAT"}  # OLDP-2 
int     nclip       {5, prompt="Number of clipping iterations for IMSTAT"}  # OLDP-2
bool    fl_boundary {no, prompt="Try to match quadrant boundaries (NIRI)"}  # OLDP-2
char    logfile     {"", prompt="Logfile"}                              # OLDP-1
bool    verbose     {yes, prompt="Verbose"}                             # OLDP-4
int     status      {0, prompt="Exit status (0=good)"}                  # OLDP-4
struct  *scanfile   {prompt="Internal use"}                             # OLDP-4

begin

    char    l_inimages = ""
    char    l_outimages = ""
    char    l_outprefix = ""
    char    l_logfile = ""
    char    l_statsec = ""
    char    l_stattype = ""
    real    l_sigma
    int     l_nclip
    bool    l_verbose
    bool    l_fl_boundary

    char    l_sci_ext = ""
    char    l_var_ext = ""
    char    l_dq_ext = ""
    char    l_key_gain = ""
    char    l_key_ron = ""
    
    char    l_temp, badhdr
    char    in[1000], out[1000], keyfound
    char    tmpout, tmpfile
    int     junk
    int     i, ii, j, k, nimages, noutimages, maxfiles, nbad, nnew
    int     xsize, ysize, halfx, halfy
    int     x1, x2, y1, y2
    char    section[32], l_expression
    real    tempmean, mean[32], quadmean[4], l_gain, l_ron
    real    boundary[8], stddev, lower, upper
    struct  l_struct

    status = 0
    nimages = 0
    maxfiles = 1000
    tempmean = 0.
    tmpfile = mktemp("tmpin")
    tmpout = ""

    cache("gnirs", "niri", "imgets", "gemextn", "gemdate")

    # set the local variables
    
    junk = fscan (inimages, l_inimages)
    junk = fscan (outimages, l_outimages)
    junk = fscan (outprefix, l_outprefix)
    junk = fscan (logfile, l_logfile)
    junk = fscan (statsec, l_statsec)
    junk = fscan (stattype, l_stattype)
    l_sigma = sigma
    l_nclip = nclip
    l_verbose = verbose
    l_fl_boundary = fl_boundary
    
    badhdr = ""
    junk = fscan (nsheaders.sci_ext, l_sci_ext)
    if ("" == l_sci_ext) badhdr = badhdr + " sci_ext"
    junk = fscan (nsheaders.var_ext, l_var_ext)
    if ("" == l_var_ext) badhdr = badhdr + " var_ext"
    junk = fscan (nsheaders.dq_ext, l_dq_ext)
    if ("" == l_dq_ext) badhdr = badhdr + " dq_ext"
    junk = fscan (nsheaders.key_gain, l_key_gain)
    if ("" == l_key_gain) badhdr = badhdr + " key_gain"
    junk = fscan (nsheaders.key_ron, l_key_ron)
    if ("" == l_key_ron) badhdr = badhdr + " key_ron"

    #------------------------------------------------------------------------
    # Check for package log file or user-defined log file
    if (l_logfile == "") {
        junk = fscan (gnirs.logfile, l_logfile)
        if (l_logfile == "") {
            l_logfile="gnirs.log"
            printlog ("WARNING - NVNOISE: Both nvnoise.logfile and \
                gnirs.logfile are empty.", l_logfile, verbose+)
            printlog ("                     Using default file gnirs.log.",
                l_logfile, verbose+)
        }
    }
    # Open log file
    date | scan (l_struct)
    printlog ("-------------------------------------------------------------\
        ---------------", l_logfile, l_verbose)
    printlog ("NVNOISE -- "//l_struct, l_logfile, l_verbose)
    printlog(" ",l_logfile, l_verbose)

    # Check for header errors
    if (badhdr != "") {
        printlog ("ERROR - NVNOISE: Parameter(s) missing from \
            nsheaders: "//badhdr, l_logfile, verbose+)
        status = 1
        goto clean
    }

    #-----------------------------------------------------------------------
    # Load up arrays of input name lists

    if (substr (l_inimages, 1, 1) == "@") {
        l_temp = substr(l_inimages, 2, strlen(l_inimages))
        if (access(l_temp) == no) {
            printlog ("ERROR - NVNOISE:  Input file "//l_temp//" not found.",
                l_logfile,verbose+)
            status = 1
            goto clean
        }
    }

    nimages = 0
    nbad = 0
    files (l_inimages, sort-, > tmpfile) 
    scanfile = tmpfile

    while (fscan (scanfile, l_temp) != EOF) {
        gimverify (l_temp)
        if (gimverify.status == 1) {
            printlog ("ERROR - NVNOISE: File "//l_temp//" not found.",
                l_logfile,verbose+)
            nbad += 1
        } else if (gimverify.status > 1) {
            printlog ("ERROR - NVNOISE: File "//l_temp//" not a MEF FITS \
                image.",l_logfile,verbose+)
            nbad += 1
        } else {
            # strip .fits if present
            if (substr (l_temp, strlen(l_temp)-4, strlen(l_temp)) == ".fits")
            l_temp = substr (l_temp, 1, (strlen(l_temp)-5))
            nimages = nimages + 1
            in[nimages] = l_temp 
            if (nimages > 1) {
                for (j=1; j<nimages; j+=1) {
                    if (in[nimages] == in[j]) {
                        printlog ("WARNING - NVNOISE: Input image name "//\
                            in[nimages]//" repeated.", l_logfile, verbose+)
                        printlog("                   Not including it again.",
                            l_logfile, verbose+)
                        nimages = nimages - 1
                    }
                }
            }
        }
    } # end while

    # Exit if problems found
    if (nbad > 0) {
        printlog ("ERROR - NVNOISE: "//nbad//" image(s) either do not exist \
            or are not MEF files,", l_logfile, verbose+)
        status = 1
        goto clean
    }
    if (nimages > maxfiles) {
        printlog ("ERROR - NVNOISE: Maximum number of input images exceeded \
            ("//str(maxfiles)//")", l_logfile, verbose+ )
        status = 1
        goto clean
    }
    if (nimages == 0) {
        printlog ("ERROR - NVNOISE: No valid input images.",
            l_logfile, verbose+ )
        status = 1
        goto clean
    }

    printlog ("Processing "//nimages//" files.", l_logfile, l_verbose)
    printlog (" ", l_logfile, l_verbose)
#print headline output to log
    if (nimages>0) {
        printf ("%s\n",
	  "Repeating 8-column bias correction by quadrant, clockwise from upper left",>> l_logfile)
        if (l_verbose) {
            printf ("%s\n","Repeating 8-column bias correction by quadrant, clockwise from upper left:")
        }
    }
    scanfile = ""
    delete (tmpfile,ver-, >& "dev$null")

    #---------------------------------------------------------------------------
    # Now for the output images
    # outimages could contain legal * if it is of a form like %st%stX%*.imh

    noutimages = 0
    if (l_outimages != "") {
        if (substr (l_outimages, 1, 1) == "@") 
            scanfile = substr (l_outimages, 2, strlen(l_outimages))
        else if (stridx ("*", l_outimages) > 0)  {
            files (l_outimages, sort-) | \
                match (".hhd", stop+, print-, metach-, > tmpfile)
            scanfile = tmpfile
        } else {
            files (l_outimages, sort-, > tmpfile)
            scanfile = tmpfile
        }

        while (fscan (scanfile, l_temp) != EOF) {
            noutimages = noutimages + 1
            if (noutimages > maxfiles) {
                printlog ("ERROR - NVNOISE: Maximum number of output images \
                    exceeded ("//str(maxfiles)//")",l_logfile,verbose+)
                status = 1
                goto clean
            }
            out[noutimages] = l_temp 
            if (imaccess(out[noutimages])) {
                printlog ("ERROR - NVNOISE: Output image "//\
                    out[noutimages]//" already exists", l_logfile, verbose+)
                status = 1
            }
        }
        if (status != 0) goto clean

        # if there are too many or too few output images exit with error
        if (nimages != noutimages) {
            printlog ("ERROR - NVNOISE: Number of input and output images \
                are unequal.", l_logfile, verbose+)
            status = 1
            goto clean
        }
    } else { #use prefix instead
        if (l_outprefix == "") {
            printlog ("ERROR - NVNOISE: Neither output image name nor \
                output prefix is defined.", l_logfile, verbose+)
            status = 1
            goto clean
        }
        i = 1
        nnew = 1
        while (i <= nimages) {
            out[i] = l_outprefix//in[i]
            for (j=1; j<i; j+=1) {
                if (out[i] == out[j]) {
                    printlog ("WARNING - NVNOISE: Output image name "//\
                        out[i]//" repeated.", l_logfile, verbose+)
                    printlog("                   Appending _"//nnew//" to \
                        the output file name.", l_logfile, verbose+)
                    out[i] = out[i]//"_"//nnew
                    nnew = nnew + 1
                }
            }

            if (imaccess(out[i])) {
                printlog("ERROR - NVNOISE: Output image "//out[i]//\
                    " already exists.", l_logfile, verbose+)
                status = 1
            }
            i = i + 1
        }
        if (status != 0) goto clean
    }
    scanfile = ""
    delete (tmpfile, ver-, >& "dev$null")

    #--------------------------------------------------------------------------
    # MAIN LOOP

    i = 1
    while (i <= nimages) {

        # See if input file has been *PREPARED:
        keyfound = ""
        hselect (in[i]//"[0]", "*PREPAR*", yes) | scan (keyfound)
        if (keyfound == "") {
            l_sci_ext = "1"
            l_var_ext = ""
            l_dq_ext = ""
        }

        # Parse statsec
        # Changed 1-19-2007 to be PyRAF-friendly - RIJ

        l_temp = substr(l_statsec, (stridx("[",l_statsec)+1), strlen(l_statsec))
        x1 = int(substr(l_temp, 1, stridx(':', l_temp) - 1))
        l_temp = substr(l_temp, (stridx(":",l_temp)+1), strlen(l_temp))
        x2 = int(substr(l_temp, 1, stridx(',', l_temp) - 1))
        l_temp = substr(l_temp, (stridx(",",l_temp)+1), strlen(l_temp))
        y1 = int(substr(l_temp, 1, stridx(':', l_temp) - 1))
        l_temp = substr(l_temp, (stridx(":",l_temp)+1), strlen(l_temp))
        y2 = int(substr(l_temp, 1, stridx(']', l_temp) - 1))
#        print ("["//x1//":"//x2//","//y1//":"//y2//"]")
        # Get image size.  Note that for GNIRS, two top rows are "missing"
        # We'll hardwire in a small correction... (otherwise assume the
        # quadrants are the same size)
        hselect (in[i]//"["//l_sci_ext//"]", "naxis1", yes) | scan (xsize)
        hselect (in[i]//"["//l_sci_ext//"]", "naxis2", yes) | scan (ysize)
        halfx = int(xsize/2)
        halfy = int(ysize/2)
        hselect (in[i]//"[0]", "INSTRUME", yes) | scan (l_temp)
        if (l_temp == "GNIRS") halfy = halfy + 1

        # check for errors in statsec parameters
        if ( (x1>=(halfx-16)) || (y1>=(halfy-16)) || (x2<=(halfx+16)) || \
            (y2<=(halfy+16))) {
            printlog ("ERROR - NVNOISE: statsec must include all quadrants",
                l_logfile, verbose+)
            status = 1
        }
        if (x2>xsize) x2 = xsize
        if (y2>ysize) y2 = ysize
        if (x1<1)     x1 = 1
        if(y1<1)      y1 = 1

        # Fix bad pixels before statistics
        tmpout = mktemp("tmpout")
        if (imaccess(in[i]//"["//l_dq_ext//"]")) {
            imstat (in[i]//"["//l_sci_ext//"]"//l_statsec,fi="mode",
                nclip=l_nclip, lsigma=l_sigma, usigma=l_sigma, binwidth=0.1,
                lower=INDEF, upper=INDEF, format-) | scan(mean[1])
            l_expression = "(a>0) ? "//mean[1]//" : b"
            imexpr (l_expression, tmpout,in[i]//"["//l_dq_ext//"]",
                in[i]//"["//l_sci_ext//"]", >& "dev$null")
        } else {
            imcopy (in[i]//"["//l_sci_ext//"]", tmpout, >& "dev$null")
        }

        # define statistics sections by amplifier:
        for (j=0; j<=halfx; j+=8) {
            if ((halfx-j) > x1) junk = (halfx-j)
        }
        x1 = junk - 7
        for (j=1; j<=8; j+=1) {
            printf ("[%d:%d:8,%d:%d]\n", (x1+j-1), halfx, y1, halfy) | \
                scan(section[j])
            printf ("[%d:%d:8,%d:%d]\n", (halfx+j), x2, y1, halfy) | \
                scan(section[j+8])
            printf ("[%d:%d:8,%d:%d]\n", (x1+j-1), halfx, (halfy+1), y2) | \
                scan(section[j+16])
            printf ("[%d:%d:8,%d:%d]\n", (halfx+j), x2, (halfy+1), y2) | \
                scan(section[j+24])
        }

        # Get statistics on input image with bad pixels fixed
        for (j=1; j<=32; j+=1) {
            imstat (tmpout//section[j], fi=l_stattype,
                nclip=l_nclip, lsigma=l_sigma, usigma=l_sigma, binwidth=0.1,
                lower=INDEF, upper=INDEF, format-) | scan(mean[j])
        }
        imdelete (tmpout, ver-, >& "dev$null")    

        # Compute deltas; note that we do not force columns to zero, 
        # but total net change=0
        tempmean = 0.
        for (j=1; j<=32; j+=1) {
            tempmean = tempmean + mean[j]
        }
        tempmean = tempmean / 32.
        for (j=1; j<=32; j+=1) {
            mean[j] = mean[j] - tempmean
        }

        # define output sections:
        # Note that this assumes CENTERED subarrays, if subarrays are used
        for (j=1; j<=8; j+=1) {
            printf ("[%d:%d:8,1:%d]\n", j, halfx, halfy) | \
                scan(section[j])
            printf ("[%d:%d:8,1:%d]\n", (halfx+j), xsize, halfy) | \
                scan(section[j+8])
            printf ("[%d:%d:8,%d:%d]\n", j, halfx, (halfy+1), ysize) | \
                scan(section[j+16])
            printf ("[%d:%d:8,%d:%d]\n", (halfx+j), xsize, (halfy+1), ysize) | \
                scan(section[j+24])
        }

        # Make output image
        if (keyfound != "") {
            gemextn (in[i], check="exists,mef,table", process="append",
                index="", extname="MDF", extversion="", ikparams="", omit="",
                replace="", outfile="dev$null", logfile="", glogpars="",
                verbose-, >& "dev$null")
            if (gemextn.fail_count != 0 || gemextn.count == 0) {
                if (imaccess(in[i]//"["//l_var_ext//"]"))
                    fxcopy (in[i], out[i]//".fits", "0-3", new_file+,
                        >& "dev$null")
                else
                    fxcopy (in[i], out[i]//".fits", "0-1", new_file+,
                        >& "dev$null")
            } else {
                if (imaccess(in[i]//"["//l_var_ext//"]"))
                    fxcopy (in[i], out[i]//".fits", "0-4", new_file+,
                        >& "dev$null")
                else
                    fxcopy (in[i], out[i]//".fits", "0-2", new_file+,
                        >& "dev$null")
            }
        } else   
            fxcopy (in[i], out[i]//".fits", "0-1", new_file+, >& "dev$null")

        for (j=1; j<=32; j+=1) {
            tmpout = mktemp("tmpout")
            l_expression = "(a-"//mean[j]//")"
            imexpr (l_expression, tmpout,
                in[i]//"["//l_sci_ext//"]"//section[j], verbose-, >& "dev$null")
            imcopy (tmpout, out[i]//"["//l_sci_ext//",overwrite]"//section[j],
                >& "dev$null")
            imdelete (tmpout, ver-, >& "dev$null")    
        }

        # Try to match quadrant boundaries to take out quadrant offsets better
        if (l_fl_boundary) {
            tempmean = 0.
            imstat (out[i]//"["//l_sci_ext//"]"//l_statsec, fi="midpt,stddev",
                binwidth=0.1, lower=INDEF, upper=INDEF, format-) | \
                scan (tempmean, stddev)
            imstat (out[i]//"["//l_sci_ext//"]"//l_statsec, fi="mode,stddev",
                lower=(tempmean-3*stddev), upper=(tempmean+3*stddev), 
                binwidth=0.1, format-) | scan(tempmean,stddev)
            imstat (out[i]//"["//l_sci_ext//"]"//l_statsec, fi="mode,stddev",
                lower=(tempmean-3*stddev), upper=(tempmean+3*stddev), 
                binwidth=0.1, format-) | scan(tempmean,stddev)
            lower = tempmean - 3*stddev
            upper = tempmean + 3*stddev

            section[1] = "["//(halfx-50)//":"//(halfx-10)//","//\
                              (halfy-200)//":"//(halfy-10)//"]"
            section[2] = "["//(halfx+10)//":"//(halfx+50)//","//\
                              (halfy-200)//":"//(halfy-10)//"]"
            section[3] = "["//(halfx-50)//":"//(halfx-10)//","//\
                              (halfy+10)//":"//(halfy+200)//"]"
            section[4] = "["//(halfx+10)//":"//(halfx+50)//","//\
                              (halfy+10)//":"//(halfy+200)//"]"
            section[5] = "["//(halfx-200)//":"//(halfx-10)//","//\
                              (halfy-50)//":"//(halfy-10)//"]"
            section[6] = "["//(halfx-200)//":"//(halfx-10)//","//\
                              (halfy+10)//":"//(halfy+50)//"]"
            section[7] = "["//(halfx+10)//":"//(halfx+200)//","//\
                              (halfy-50)//":"//(halfy-10)//"]"
            section[8] = "["//(halfx+10)//":"//(halfx+200)//","//
                              (halfy+10)//":"//(halfy+50)//"]"
            for (j=1; j<=8; j+=1) {
                imstat (out[i]//"["//l_sci_ext//"]"//section[j], fi="mode",
                    binwidth=0.1, lower=lower, upper=upper, format-) | \
                    scan (boundary[j])
            } 
            quadmean[1] = boundary[1]-boundary[2] + boundary[7]-boundary[8]
            quadmean[1] = (quadmean[1] + boundary[5]-boundary[6] + \
                           boundary[3]-boundary[4]) / 2.
            quadmean[2] = boundary[7] - boundary[8]
            quadmean[3] = boundary[3] - boundary[4]
            quadmean[4] =  0.
            tempmean = (quadmean[1]+quadmean[2]+quadmean[3]+quadmean[4]) / 4.
            quadmean[1] = quadmean[1] - tempmean
            quadmean[2] = quadmean[2] - tempmean
            quadmean[3] = quadmean[3] - tempmean
            quadmean[4] = quadmean[4] - tempmean

            section[1] = "[1:"//halfx//",1:"//halfy//"]"
            section[2] = "["//(halfx+1)//":"//xsize//",1:"//halfy//"]"
            section[3] = "[1:"//halfx//","//(halfy+1)//":"//ysize//"]"
            section[4] = "["//(halfx+1)//":"//xsize//","//\
                              (halfy+1)//":"//ysize//"]"
            for (j=1; j<=4; j+=1) {
                tmpout = mktemp("tmpout")
                l_expression = "(a-"//quadmean[j]//")"
                imexpr (l_expression, tmpout,
                    out[i]//"["//l_sci_ext//"]"//section[j], verbose-,
                    >& "dev$null")
                imcopy (tmpout,
                    out[i]//"["//l_sci_ext//",overwrite]"//section[j],
                    >& "dev$null")
                imdelete (tmpout, ver-, >& "dev$null")                    
            }

            # adjust the values for printed output
            for (j=1; j<=32; j+=1) {
                if (j<9)        mean[j] = mean[j] + quadmean[1]
                else if (j<17)  mean[j] = mean[j] + quadmean[2]
                else if (j<25)  mean[j] = mean[j] + quadmean[3]
                else            mean[j] = mean[j] + quadmean[4]
            }
        }

        # Regenerate the variance plane from the original data
        if ((keyfound != "") && imaccess(in[i]//"["//l_var_ext//"]")) {
            hselect (in[i]//"[0]", l_key_gain, yes) | scan (l_gain)
            hselect(in[i]//"[0]", l_key_ron, yes) | scan (l_ron)
            l_ron = (l_ron / l_gain)**2
            l_expression = "max(a/"//l_gain//",0.) + "//l_ron
            tmpout = mktemp("tmpout")
            imexpr (l_expression, tmpout, out[i]//"["//l_sci_ext//"]",
                outtype="real", verbose-, >& "dev$null")
            imcopy (tmpout, out[i]//"["//l_var_ext//",overwrite]",
                >& "dev$null")
            gemhedit (out[i]//"["//l_var_ext//"]", field="EXTVER", value=1,
                comment="", delete-)
            imdelete(tmpout, ver-, >& "dev$null")    
        }

        # update the header
        gemdate ()
        gemhedit (out[i]//"[0]", "NVNOISE", gemdate.outdate, 
            "UT Time stamp for NVNOISE")
        gemhedit (out[i]//"[0]","GEM-TLM", gemdate.outdate, 
            "UT Last modification with GEMINI", delete-)
        gemhedit (out[i]//"[0]", "NVNSTAT", l_stattype,
            "Statistics type for NVNOISE")
        gemhedit (out[i]//"[0]", "NVNSIGMA", l_sigma,
            "Sigma clipping for NVNOISE stats")
        gemhedit (out[i]//"[0]", "NVNNCLIP", l_nclip,
            "Sigma clipping iterations for NVNOISE")

        l_temp = ""
        
	#print bias corrections by quadrant, clockwise from upper left
        if (l_verbose) {
            printf ("%19s %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f\n",
                out[i], mean[17], mean[18], mean[19], mean[20], mean[21], 
                mean[22], mean[23], mean[24])
            printf ("%19s %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f\n",
                l_temp, mean[25], mean[26], mean[27], mean[28], mean[29],
                mean[30], mean[31], mean[32])
            printf ("%19s %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f\n",
                l_temp, mean[9], mean[10], mean[11], mean[12], mean[13], 
                mean[14], mean[15], mean[16])
            printf ("%19s %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f\n",
                l_temp, mean[1], mean[2], mean[3], mean[4], mean[5],
                mean[6], mean[7], mean[8])
        }
        printf ("%19s %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f\n",
            out[i], mean[17], mean[18], mean[19], mean[20], mean[21], mean[22],
            mean[23], mean[24], >> l_logfile)
        printf ("%19s %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f\n",
            l_temp, mean[25], mean[26], mean[27], mean[28], mean[29], mean[30],
            mean[31], mean[32], >> l_logfile)
        printf ("%19s %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f\n",
            l_temp, mean[9], mean[10], mean[11], mean[12], mean[13], mean[14],
            mean[15], mean[16], >> l_logfile)
        printf ("%19s %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f\n",
            l_temp, mean[1], mean[2], mean[3], mean[4], mean[5], mean[6],
            mean[7], mean[8], >> l_logfile)
        printlog (" ", l_logfile, l_verbose)

        i += 1
    }
    # end the main loop

    #---------------------------------------------------------------------------
    # Clean up
clean:
    if (status == 0)
        printlog("NVNOISE exit status:  good.", l_logfile, l_verbose)
    else
        printlog("NVNOISE exit status:  failed.", l_logfile, l_verbose)
    printlog("----------------------------------------------------------------------------",l_logfile,l_verbose)

    scanfile = ""
    delete (tmpfile, ver-, >& "dev$null")
end



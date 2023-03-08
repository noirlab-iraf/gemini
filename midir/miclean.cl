# Copyright(c) 2009-2011 Association of Universities for Research in Astronomy, Inc.

procedure miclean (inimages)

# This task is used to remove certain types of pattern noise from 
# Michelle or T-ReCS images.  It would normally be used with stacked 
# files output from "mireduce" but it can be used with raw data files 
# or output files from "mprepare" if desired.
#

char    inimages    {prompt="Input Michelle/T-ReCS image(s)"}
char    rawpath     {"", prompt="Path for input raw images"}
char    outimages   {"", prompt="Output image(s)"}
char    outpref     {"c", prompt="Prefix for output image name(s)"}
bool    fl_stair    {yes, prompt="Do staircase cleaning"}
bool    fl_columns  {no, prompt="Clean columns"}
bool    fl_rows     {no, prompt="Clean rows"}
bool    fl_simage   {no, prompt="Make staircase image(s) yes/no?"}
char    stairpref   {"s", prompt="Optional prefix for output staircase image(s)"}
bool    fl_nimage   {no, prompt="Make row/column noise image(s) yes/no?"}
char    noisepref   {"n", prompt="Optional prefix for output noise image(s)"}
real    threshold   {3.0, prompt="Number of sigmas for masking bright regions"}
char    logfile     {"", prompt="Logfile name"}
bool    verbose     {yes, prompt="Verbose logging yes/no?"}
int     status      {0, prompt="Exit error status: (0=good, >0=bad)"}
struct  *scanfile   {"", prompt="Internal use only"}

begin
    char    l_inimages = ""
    char    l_rawpath = ""
    char    l_outimages = ""
    char    l_outpref = ""
    char    l_stairpref = ""
    char    l_noisepref = ""
    char    l_logfile = ""
    real    l_threshold
    bool    l_fl_stair, l_fl_columns, l_fl_rows, l_fl_simage, l_fl_nimage
    bool    l_verbose

    char    paramstr, lastchar, errmsg, filename, outputstr, keyfound
    char    tmpinimg, tmpoutimg, tmpfile
    char    tmpworkimg, tmpstatimg, tmpgauss
    char    tmpsubimg, tmpwithnoise, tmpstair, tmpnoise
    char    in[100], out[100], stairimage, noiseimage
    char    inextension, outextension, stairextension, noiseextension
    char    section, chsection, outsection
    char    instrument
    int     junk, i, l, m, row, col, maxfiles, nimages, noutimages, tmpstatus
    int     iext, nextensions
    int     iaxis, naxes, naxis[4]
    real    average, sigma, lothresh, upthresh, chlevel, rawoffset

    cache ("gemdate", "gloginit", "gemextn")

    junk = fscan (inimages, l_inimages)
    junk = fscan (rawpath, l_rawpath)
    junk = fscan (outimages, l_outimages)
    junk = fscan (outpref, l_outpref)
    l_fl_stair = fl_stair
    l_fl_columns = fl_columns
    l_fl_rows = fl_rows
    l_fl_simage = fl_simage
    junk = fscan (stairpref, l_stairpref)
    l_fl_nimage = fl_nimage
    junk = fscan (noisepref, l_noisepref)
    l_threshold = threshold
    junk = fscan (logfile, l_logfile)
    l_verbose = verbose

    # Initialize
    status = 0
    tmpstatus = 0
    maxfiles = 100

    # Create tmp file names
    tmpinimg = mktemp ("tmpinimg")
    tmpoutimg = mktemp ("tmpoutimg")
    tmpfile = mktemp ("tmpfile")

    # Create the list of parameter/value pairs.  One pair per line.
    # All lines combined into one string.  Line delimiter is '\n'.
    paramstr  = "inimages       = "//inimages.p_value//"\n"
    paramstr += "rawpath        = "//rawpath.p_value//"\n"
    paramstr += "outimages      = "//outimages.p_value//"\n"
    paramstr += "outpref        = "//outpref.p_value//"\n"
    paramstr += "fl_stair       = "//fl_stair.p_value//"\n"
    paramstr += "fl_columns     = "//fl_columns.p_value//"\n"
    paramstr += "fl_rows        = "//fl_rows.p_value//"\n"
    paramstr += "fl_simage      = "//fl_simage.p_value//"\n"
    paramstr += "stairpref      = "//stairpref.p_value//"\n"
    paramstr += "fl_nimage      = "//fl_nimage.p_value//"\n"
    paramstr += "noisepref      = "//noisepref.p_value//"\n"
    paramstr += "threshold      = "//threshold.p_value//"\n"
    paramstr += "logfile        = "//logfile.p_value//"\n"
    paramstr += "verbose        = "//verbose.p_value
    
    # Assign a logfile name, if not given.  Open logfile and start log.
    # Write parameter/value pairs ("paramstr") to log.
    gloginit (l_logfile, "miclean", "midir", paramstr, fl_append+,
        verbose=l_verbose)
    if (gloginit.status != 0) {
        status = gloginit.status
        goto exitnow
    }
    l_logfile = gloginit.logfile.p_value

    # Add the trailing slash to rawpath, if missing
    lastchar = substr (l_rawpath, strlen(l_rawpath), strlen(l_rawpath))
    if ((l_rawpath != "") && (lastchar != "$") && (lastchar != "/"))
        l_rawpath = l_rawpath//"/"


    # Validate user input parameter values
    if (l_fl_rows && l_fl_columns) {
        errmsg = "Please, select fl_rows OR fl_columns, not both."
        status = 121
        glogprint (l_logfile, "miclean", "status", type="error", errno=status,
            str=errmsg, verbose=yes)
    }    
    if (fl_simage && (stairpref == "")) {
        errmsg = "Stair image(s) requested but stairpref not defined."
        status = 121
        glogprint (l_logfile, "miclean", "status", type="error", errno=status,
            str=errmsg, verbose=yes)
    }
    if (fl_nimage && (noisepref == "")) {
        errmsg = "Noise image(s) requested but stairpref not defined."
        status = 121
        glogprint (l_logfile, "miclean", "status", type="error", errno=status,
            str=errmsg, verbose=yes)
    }    
    if (status != 0)    goto clean
    

    # Adjustment based on user input parameter values
    if (l_threshold < 2.0) {
        l_threshold = 0.0
        glogprint (l_logfile, "miclean", "task", type="warning",
            str="Threshold less than 2.0.  Resetting to 0.0",
            verbose=l_verbose)
    }

    
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
            errmsg = "Error while counting the input images."
            status = 99
        }
    }
    if (status != 0) {
        glogprint (l_logfile, "miclean", "status", type="error", errno=status,
            str=errmsg, verbose+)
        goto clean
    }
    # Keep tmpinimg. Might need it for output names with l_outpref.


    # Load up the array of output file names
    if (l_outimages != "")
        outputstr = l_outimages
    else if (l_outpref != "") {
        gemextn ("@"//tmpinimg, check="", process="none", index="", extname="",
            extversion="", ikparams="", omit="path", replace="",
            outfile=tmpoutimg, logfile=l_logfile, verbose=l_verbose)
        outputstr = l_outpref//"@"//tmpoutimg
    } else {
        errmsg = "Neither output image name nor output prefix is defined."
        status = 121
        glogprint (l_logfile, "miclean", "status", type="error", errno=status,
            str=errmsg, verbose=yes)
        goto clean
    }
    
    gemextn (outputstr, check="", process="none", index="", extname="",
        extversion="", ikparam="", omit="extension", replace="",
        outfile=tmpfile, logfile=l_logfile, verbose=l_verbose)
    delete (tmpoutimg, ver-, >& "dev$null")
    gemextn ("@"//tmpfile, check="absent", process="none", index="",
        extname="", extversion="", ikparams="", omit="", replace="",
        outfile=tmpoutimg, logfile=l_logfile, verbose=l_verbose)
    noutimages = gemextn.count
    delete (tmpfile, ver-, >& "dev$null")
    
    if ((gemextn.fail_count > 0) || (noutimages == 0) || \
        (noutimages != nimages)) {
        
        if (gemextn.fail_count > 0) {
            errmsg = gemextn.fail_count//" dataset(s) already exist(s)."
            status = 102
        } else if (noutimages == 0) {
            errmsg = "No output dataset names defined."
            status = 121
        } else if (noutimages != nimages) {
            errmsg = "Different number of input names ("//ninimages//") and \
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
        glogprint (l_logfile, "miclean", "status", type="error", errno=status,
            str=errmsg, verbose=yes)
        goto clean
    }
    delete (tmpinimg, ver-, >& "dev$null")


    # Check if output stair image names are valid (if fl_simage)
    if (l_fl_simage) {
        gemextn (l_stairpref//"@"//tmpoutimg, check="absent", process="none",
            index="", extname="", extversion="", ikparams="", omit="",
            replace="", outfile="dev$null", logfile=l_logfile,
            verbose=l_verbose)
        if (gemextn.fail_count > 0) {
            errmsg = gemextn.fail_count//" stair image name(s) \
                already exist(s)."
            status = 102
            glogprint (l_logfile, "miclean", "status", type="error",
                errno=status, str=errmsg, verbose=yes)
            goto clean
        }
    }
    
    # Check if output noise image names are valid (if fl_nimage)
    if (l_fl_nimage) {
        gemextn (l_noisepref//"@"//tmpoutimg, check="absent", process="none",
            index="", extname="", extversion="", ikparams="", omit="",
            replace="", outfile="dev$null", logfile=l_logfile,
            verbose=l_verbose)
        if (gemextn.fail_count > 0) {
            errmsg = gemextn.fail_count//" noise image name(s) \
                already exist(s)."
            status = 102
            glogprint (l_logfile, "miclean", "status", type="error",
                errno=status, str=errmsg, verbose=yes)
            goto clean
        }
    }
    delete (tmpoutimg, ver-, >& "dev$null")


    # Okay, all the user inputs have been checked.
    # Start processing.
    
    glogprint (l_logfile, "miclean", "status", type="string",
        str="Processing "//str(nimages)//" image(s).", verbose=l_verbose)


    # Start of main loops.  There are three possible things that are done 
    # depending on the l_fl_stair, l_fl_rows, and l_fl_columns flags.  The 
    # staircasing fix is done first if requested.
    #

    for (i=1; i<=nimages; i=i+1) {
	imgets(in[i]//"[0]","INSTRUME")
	instrument=imgets.value

        # Get the number of extensions
        fxhead (in[i], format_file="", long=no, count=no, > tmpfile)
        system.tail (tmpfile, nlines=1) | \
            fields ("STDIN", "1", lines="1", quit=no, print=no) | \
            scan (nextensions)
        delete (tmpfile, ver-, >& "dev$null")

        # Prepare output images, copy PHU.
        fxcopy (in[i]//".fits", out[i], "0", new_file=yes, verbose=l_verbose)
        if (l_fl_simage) {
            stairimage = l_stairpref//out[i]
            fxcopy (in[i]//".fits", stairimage, "0", new_file=yes,
                verbose=l_verbose)
        }
        if (l_fl_nimage) {
            noiseimage = l_noisepref//out[i]
            fxcopy (in[i]//".fits", noiseimage, "0", new_file=yes,
                verbose=l_verbose)
        }

        for (iext=1; iext<=nextensions; iext=iext+1) {
        
            inextension = in[i]//"["//str(iext)//"]"
            outextension = out[i]//"["//str(iext)//"]"
            if (l_fl_simage)
                stairextension = stairimage//"["//str(iext)//"]"
            if (l_fl_nimage)
                noiseextension = noiseimage//"["//str(iext)//"]"
            
            # Add extension to output files.  This will allow the use
            # of imcopy with section.  This extension will be overwritten
            # with the images slices calculated below.
            
            fxinsert (inextension, out[i]//"["//str(iext-1)//"]", groups="",
                verbose=no)
            if (l_fl_simage)
                fxinsert (inextension, stairimage//"["//str(iext-1)//"]",
                    groups="", verbose=no)
            if (l_fl_nimage) {
                fxinsert (inextension, noiseimage//"["//str(iext-1)//"]",
                    groups="", verbose=no)
            }
        
            # Determine the dimensions of the extension.
            # It can be 2, 3, or 4 depending on whether it is a file
            # output from "mireduce", one from "mprepare" or "tprepare",
            # or a raw data file.

            keyfound = ""
            hselect (inextension, "i_naxis", yes) | scan (keyfound)
            if (keyfound == "") {
                errmsg = "No NAXIS keyword in "//inextension//\
                    ". Skipping image."
                tmpstatus = 132
            } else {
                naxes = int(keyfound)
                if ((naxes < 2) || (naxes > 4)) {
                    errmsg = "Unexpected number of axes ("//naxes//"). \
                        Skipping image."
                    tmpstatus = 132
                } 
            }
            if (tmpstatus != 0) {
                glogprint (l_logfile, "miclean", "status", type="error",
                    errno=tmpstatus, str=errmsg, verbose=yes)
                status = tmpstatus
                tmpstatus = 0
                imdelete (out[i], ver-, >& "dev$null")
                if (l_fl_simage)    imdelete (stairimage, ver-, >& "dev$null")
                if (l_fl_nimage)    imdelete (noiseimage, ver-, >& "dev$null")
                goto nextimage
            }
            
            for (iaxis=1; iaxis<=4; iaxis=iaxis+1) {
                naxis[iaxis] = 1
            }
            for (iaxis=1; iaxis<=naxes; iaxis=iaxis+1) {
                keyfound = ""
                hselect (inextension, "NAXIS"//str(iaxis), yes) | \
                    scan(keyfound)
                if (keyfound == "") {
                    errmsg = "NAXIS"//str(iaxis)//" not found.  Skipping \
                        image."
                    status = 132
                    glogprint (l_logfile, "miclean", "status", type="error",
                        errno=status, str=errmsg, verbose=yes)
                    imdelete (out[i], ver-, >& "dev$null")
                    if (l_fl_simage)    imdelete (stairimage, ver-, >& "dev$null")
                    if (l_fl_nimage)    imdelete (noiseimage, ver-, >& "dev$null")
                    goto nextimage
                } else {
                    naxis[iaxis] = int(keyfound)
                }
            }
            
            if ((naxis[1] != 320) || (naxis[2] != 240)) {
                errmsg = "Extension "//inextension//" has bad dimensions ["//\
                    str(naxis[1])//","//str(naxis[2])//"].  Skipping."
                status = 123
                glogprint (l_logfile, "miclean", "status", type="error",
                    errno=status, str=errmsg, verbose=l_verbose)
                imdelete (out[i], ver-, >& "dev$null")
                if (l_fl_simage)    imdelete (stairimage, ver-, >& "dev$null")
                if (l_fl_nimage)    imdelete (noiseimage, ver-, >& "dev$null")
                goto nextimage
            }
            

            # Now we know the data structure.  Let's get to work.
                        
                
            # Loop over slices (nods, etc.) in the extensions
            #
            # l => index for 3rd dimension
            # m => index for 4th dimension

            for (l=1; l<=naxis[3]; l=l+1) {
                for (m=1; m<=naxis[4]; m=m+1) {
                    section = "[1:320,1:240"
                    if (naxes == 2) {
                        section = section//"]"
                        outsection = section
                    } else if (naxes == 3) {
                        section = section//","//str(l)//"]"
                        outsection = section
                    } else if (naxes == 4) {
                        section = section//","//str(l)//","//str(m)//"]"
                        outsection = section
                    }

                    # Create a copy of input that we will contain the
                    # incremental corrections.

                    tmpworkimg = mktemp("tmpworkimg")
                    imcopy (inextension//section, tmpworkimg//".fits",
                        verbose=no, >& "dev$null")

                    # Create a copy that will be used to calculate
                    # statistics.  This image might or might not be masked.
                    
                    tmpstatimg = mktemp("tmpstatimg")
                    imcopy (inextension//section, tmpstatimg//".fits",
                        verbose=no, >& "dev$null")

                    if ((naxes == 3 && l < 3) ||
                        (naxes == 4 && instrument == "michelle" && l != 1)) {
                            imstatistics (tmpstatimg, fields="mean",
                                lower=INDEF, upper=INDEF, format=no, nclip=0.,
                                lsigma=3., usigma=3., binwidth=0.1,
                                cache=no) | scan (rawoffset)
                    }
                    else {
                        rawoffset=0.0
                    }

                    # Apply threshold masking, if requested
                    # This way of thresholding by setting the pixel values
                    # to zero works only because the average background is 
                    # already very close to zero.

                    # This should not be used with a raw image as the 
                    # mean is only near zero for a difference image.
                    
                    if (l_threshold >= 2.0) {
                        imstatistics (tmpstatimg, fields="mean,stddev",
                            lower=INDEF, upper=INDEF, format=no, nclip=0.,
                            lsigma=3., usigma=3., binwidth=0.1,
                            cache=no) | scan (average, sigma)
                        if ((isindef(average) == no) && 
                            (isindef(sigma) == no)) {

                            upthresh = average + sigma * l_threshold
                            lothresh = average - sigma * l_threshold
                            imreplace (tmpstatimg, 0., lower=INDEF,
                                upper=lothresh)
                            imreplace (tmpstatimg, 0., lower=upthresh,
                                upper=INDEF)
                        }
                    }


                    # Apply staircase correction, if requested
                    if (l_fl_stair) {
                        # It is not clear how well the cleaning will work on 
                        # raw data files or the output files from 
                        # mprepare/tprepare. It will usually be better to 
                        # stack up the images first and then run this routine.

                        # The "raw offset" value is used to preserve the 
                        # average over the whole image after staircasing in 
                        # the cases where a raw data image is being cleaned.
                        
                        # Smooth the input image with a gaussian to displace
                        # the median from any of the actual pixel values.

                        tmpgauss = mktemp("tmpgauss")
                        gauss (tmpstatimg, tmpgauss, 1.0, nsigma=10., ratio=1.,
                            theta=0., bilinear=yes, boundary="nearest",
                            constant=0.)

                        # Find the median of each channel and substract it off
                        for (col=1; col<=naxis[1]; col=col+20) {

                            chsection = "["//str(col)//":"//str(col+19)//\
                                ",1:"//str(naxis[2])//"]"
                            imstatistics (tmpgauss//chsection, fields="midpt",
                                lower=INDEF, upper=INDEF, format=no, nclip=0.,
                                lsigma=3., usigma=3., binwidth=0.1,
                                cache=no) | scan (chlevel)
                            chlevel=chlevel-rawoffset
                            if (isindef(chlevel) == no) {
                                tmpsubimg = mktemp ("tmpsubimg")
                                imarith (tmpworkimg//chsection, "-", chlevel,
                                    tmpsubimg)
                                imcopy (tmpsubimg, tmpworkimg//chsection,
                                    verbose=no)
                                imdelete (tmpsubimg, ver-, >& "dev$null")
                            }
                        }
                        imdelete (tmpgauss, ver-, >& "dev$null")
                        
                        if (l_fl_simage) {
                            tmpstair = mktemp ("tmpstair")
                            imarith (inextension//section, "-", tmpworkimg,
                                tmpstair, verbose=no)   

                            imcopy (tmpstair//outsection,
                                stairextension//section, verbose=no)  
 
                            imdelete (tmpstair, ver-, >& "dev$null")
                        }

                        # Now "tmpworkimg" is the stair corrected slice

                        # Create a new masked tmpstatimg from the staircase
                        # corrected image.  It improves the final result a
                        # bit if you have somewhat large positive or negative
                        # 'sources'.
                        
                        if (l_threshold >= 2.0) {
                            imdelete (tmpstatimg, ver-, >& "dev$null")
                            imcopy (tmpworkimg, tmpstatimg, verbose-)
                            imstatistics (tmpstatimg, fields="mean,stddev",
                                lower=INDEF, upper=INDEF, format=no, nclip=0.,
                                lsigma=3., usigma=3., binwidth=0.1,
                                cache=no) | scan (average, sigma)
                            if ((isindef(average) == no) && 
                                (isindef(sigma) == no)) {

                                upthresh = average + sigma * l_threshold
                                lothresh = average - sigma * l_threshold
                                imreplace (tmpstatimg, 0., lower=INDEF,
                                    upper=lothresh)
                                imreplace (tmpstatimg, 0., lower=upthresh,
                                    upper=INDEF)
                            }
                        }
                    }
                    
                    # Keep a copy of this slice if l_fl_nimage = yes
                    if (l_fl_nimage && (l_fl_rows || l_fl_columns)) {
                        tmpwithnoise = mktemp ("tmpwithnoise")
                        imcopy (tmpworkimg, tmpwithnoise, verbose=no)
                    }
                                        

                    # Fix row or column noise pattern
                    if (l_fl_rows) {
                        for (row=1; row<=naxis[2]; row=row+1) {
                            chsection = "[1:320,"//str(row)//":"//str(row)//"]"
                            imstatistics (tmpstatimg//chsection,
                                fields="midpt", lower=INDEF, upper=INDEF,
                                format=no, nclip=0., lsigma=3., usigma=3.,
                                binwidth=0.1, cache=no) | scan (chlevel)
                            chlevel=chlevel-rawoffset
                            if (isindef(chlevel) == no) {
                                tmpsubimg = mktemp ("tmpsubimg")
                                imarith (tmpworkimg//chsection, "-", chlevel,
                                    tmpsubimg)
                                imcopy (tmpsubimg, tmpworkimg//chsection,
                                    verbose=no)
                                imdelete (tmpsubimg, ver-, >& "dev$null")
                            }
                        }
   
                        # Now "tmpworkimg" is the row corrected slice
                        
                    } else if (l_fl_columns) {
                        for (col=1; col<=naxis[1]; col=col+1) {
                            chsection = "["//str(col)//":"//str(col)//",1:240]"
                            imstatistics (tmpstatimg//chsection,
                                fields="midpt", lower=INDEF, upper=INDEF,
                                format=no, nclip=0., lsigma=3., usigma=3.,
                                binwidth=0.1, cache=no) | scan (chlevel)
                            chlevel=chlevel-rawoffset
                            if (isindef(chlevel) == no) {
                                tmpsubimg = mktemp ("tmpsubimg")
                                imarith (tmpworkimg//chsection, "-", chlevel,
                                    tmpsubimg)
                                imcopy (tmpsubimg, tmpworkimg//chsection,
                                    verbose=no)
                                imdelete (tmpsubimg, ver-, >& "dev$null")
                            }
                        }
                        
                        # Now "tmpworkimg" is the column corrected slice
                    }
 
                    # We're done with tmpstatimg, delete it
                    imdelete (tmpstatimg, ver-, >& "dev$null")
                    
                    # Create the "noise" image, if requested
                    if (l_fl_nimage && (l_fl_rows || l_fl_columns)) {
                        tmpnoise = mktemp ("tmpnoise")
                        imarith (tmpwithnoise, "-", tmpworkimg, tmpnoise,
                            verbose=no)
                        imcopy (tmpnoise//outsection,
                            noiseextension//section, verbose=no)
                        imdelete (tmpnoise, ver-, >& "dev$null")
                        imdelete (tmpwithnoise, ver-, >& "dev$null")
                    }

                    # We're done correcting this slice, write it to disk
                    imcopy (tmpworkimg//"[1:320,1:240]", outextension//section,
                        verbose=no)

                    # We don't need tmpworkimg anymore
                    imdelete (tmpworkimg, ver-, >& "dev$null")
                    
                }   # end for loop over 'm'
            }   # end for loop over 'l'
        }   # end for loop over extensions


        # Time stamp
        gemdate()
        gemhedit (out[i]//"[0]", "MICLEAN", gemdate.outdate,
            "UT Time stamp for MICLEAN", delete-)
        gemhedit (out[i]//"[0]", "GEM-TLM", gemdate.outdate,
            "UT Last modification with GEMINI", delete-)

nextimage:
        ;

    }   # end for loop over 'i' (all images)

clean:
    delete (tmpinimg, ver-, >& "dev$null")
    delete (tmpoutimg, ver-, >& "dev$null")

    if (status == 0)
        glogclose (l_logfile, "miclean", fl_success+, verbose=l_verbose)
    else
        glogclose (l_logfile, "miclean", fl_success-, verbose=l_verbose)

exitnow:
    ;

end

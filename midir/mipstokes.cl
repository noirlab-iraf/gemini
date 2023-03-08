# Copyright(c) 2006-2011 Association of Universities for Research in Astronomy, Inc.

procedure mipstokes (inimages)

# This routine produces stokes parameter maps from the individual frames of a "prepared" 
# Michelle polarimetry file.  There is one output map per NOD position.
#
# Version:  September 13, 2005  KV wrote original script 
#
#           Output file has 4 extensions, these are each 320 by 240 pixel images
#           Waveplate angles are 0, 22.5, 45, and 67.5 degrees for extensions 
#           1 to 4 respectively.
#
#           There is a header keyword POLANGLE for each extension that gives the 
#           values of the polarimetry plate angle in degrees.
#
#           Jan  10, 2006   KV added the WCS parameters to the image extension from the 
#                              primary header.
#
char    inimages    {prompt="Input Michelle polarimetry image(s)"}  # OLDP-1-input-primary-single-prefix=s
char    outimages   {"",prompt="Output image(s)"}                   # OLDP-1-output
char    outpref     {"z",prompt="Prefix for out image(s)"}          # OLDP-4
char    rawpath     {"",prompt="Path for in raw images"}            # OLDP-4
char    frametype   {"dif", enum="src|ref|dif", prompt="Type of frame to combine (src, ref, dif)"}   # OLDP-2
char    combine     {"average",enum="average|sum",prompt="Combining images by average|sum"}    # OLDP-2
bool    fl_mask     {yes,prompt="Mask out intensity values near zero?"}     # OLDP-4
real    noise       {0.0,prompt="Masking value in ADU (used if > 0.0)"}     # OLDP-4
real    threshold   {0.5,min=0.5,prompt="Percent polarization value for masking intensities"}
char    blankarea   {"[*,*]",prompt="Image area used to estimate the noise level"}    # OLDP-4
bool    fl_register {no,prompt="Register images when combining"}    # OLDP-2
char    regions     {"[*,*]",prompt="Reference image regions used for registration (xregister)"}    # OLDP-4
bool    fl_stair    {yes, prompt="Correct channel offsets"}         # OLDP-4
char    logfile     {"",prompt="Logfile"}                           # OLDP-1
bool    verbose     {yes,prompt="Verbose"}                          # OLDP-4
int     status      {0,prompt="Exit status: (0=good, >0=bad)"}      # OLDP-4
struct  *scanfile   {"",prompt="Internal use only"}                 # OLDP-4

begin

    char    l_inimages = ""
    char    l_outimages = ""
    char    l_outpref = ""
    char    l_rawpath = ""
    char    l_frametype = ""
    char    l_combine = ""
    char    l_blankarea = ""
    char    l_regions = ""
    char    l_logfile = ""
    real    l_noise, l_threshold
    bool    l_fl_mask, l_fl_register, l_fl_stair, l_verbose
    
    char    lastchar, paramstr, filename, errmsg, outputstr, keyfound
    char    in[100], out[100]
    char    tmpinimg, tmpoutimg, tmpfile, tmpshifts
    char    tmpreg1, tmpreg2, tmpworkimg, tmpimg1, tmpimg2, tmpimg3, tmpimg4
    char    tmprefer1, tmprefer2, tmprefer3, tmprefer4
    char    tmpims1a, tmpims1b, tmpims2a, tmpims2b, tmpims3a, tmpims3b
    char    tmpims4a, tmpims4b
    char    tmpUstokes, tmpQstokes, tmpIstokes
    char    instrument, imgsec
    int     i, k, ext, nod, coord, junk, maximages
    int     outextn, firstextn, lastextn
    int     nimages, noutimages, nskipped, nextensions, nnods
    int     framevalue, tmpstatus, modeflag, flag[4]
    real    chlevel, sigma, floor, max

    cache ("gemdate")
    
    junk = fscan (inimages, l_inimages)
    junk = fscan (outimages, l_outimages)
    junk = fscan (outpref, l_outpref)
    junk = fscan (rawpath, l_rawpath)
    junk = fscan (frametype, l_frametype)
    junk = fscan (combine, l_combine)
    l_fl_mask     = fl_mask
    l_noise       = noise
    l_threshold   = threshold
    junk = fscan (blankarea, l_blankarea)
    l_fl_register = fl_register
    junk = fscan (regions, l_regions)
    l_fl_stair    = fl_stair
    junk = fscan (logfile, l_logfile)
    l_verbose = verbose

    # Initialize
    maximages = 100
    status = 0

    # Temporary file names
    tmpinimg = mktemp ("tmpinimg")
    tmpoutimg = mktemp ("tmpoutimg")
    tmpfile = mktemp ("tmpfile")
    if (l_fl_register) {
        tmpshifts = mktemp ("tmpshifts")
    }

    # Create the list of parameter/value pairs.  One pair per line.
    # All lines combined into one string.  Line delimiter is '\n'.
    paramstr =  "inimages       = "//inimages.p_value//"\n"
    paramstr += "outimages      = "//outimages.p_value//"\n"
    paramstr += "outpref        = "//outpref.p_value//"\n"
    paramstr += "rawpath        = "//rawpath.p_value//"\n"
    paramstr += "frametype      = "//frametype.p_value//"\n"
    paramstr += "combine        = "//combine.p_value//"\n"
    paramstr += "fl_mask        = "//fl_mask.p_value//"\n"
    paramstr += "noise          = "//noise.p_value//"\n"
    paramstr += "threshold      = "//threshold.p_value//"\n"
    paramstr += "blankarea      = "//blankarea.p_value//"\n"
    paramstr += "fl_register    = "//fl_register.p_value//"\n"
    paramstr += "regions        = "//regions.p_value//"\n"
    paramstr += "fl_stair       = "//fl_stair.p_value//"\n"
    paramstr += "logfile        = "//logfile.p_value//"\n"
    paramstr += "verbose        = "//verbose.p_value

    # Assign a logfile name if not given.  Open logfile and start log.
    # Write parameter/value pairs ("paramstr") to log.
    gloginit (l_logfile, "mipstokes", "midir", paramstr, fl_append+,
        verbose=l_verbose)
    if (gloginit.status != 0) {
        status = gloginit.status
        goto exitnow
    }
    l_logfile = gloginit.logfile.p_value
    
    if (l_frametype == "dif")
        framevalue = 3  
    else if (l_frametype == "src")
        framevalue = 1  
    else if (l_frametype == "ref")
        framevalue = 2 

    # Not sure resetting the user's outpref value is a good thing...
    # I am adding a warning to the user here.
    # -- KL    
    if (l_outpref != "") {
        if (framevalue == 1) {
            l_outpref = "c"
            status = 1
        } else if (framevalue == 2) {
            l_outpref = "a"
            status = 1
        }
        if (status != 0) {
            glogprint (l_logfile, "mipstokes", "status", type="warning",
                errno=0, str="Resetting 'outpref' to '"//l_outpref//"' \
                because 'frametype' not 'dif'", verbose=l_verbose)
        }
    }

    #----------------------------------------------------------------------
    # Check the rawpath name for a final /
    #
    # Is rawpath really required?  Aren't the inputs supposed to have
    # been mprepare'd, and should therefore no longer be 'raw' and shouldn't
    # they be in the current directory?
    # -- KL
    
    lastchar = substr (l_rawpath, strlen (l_rawpath), strlen(l_rawpath))
    if ((l_rawpath != "") && (lastchar != "$") && (lastchar != "/"))
        l_rawpath = l_rawpath//"/"
    
    #-----------------------------------------------------------------------
    # Load up arrays of input name lists

    gemextn (l_inimages, check="exists,mef", process="none", index="",
        extname="", extversion="", ikparam="", omit="kernel,exten",
        replace="^%%"//l_rawpath//"%", outfile=tmpinimg, logfile=l_logfile,
        verbose=l_verbose)
    nimages = gemextn.count
    
    if ((gemextn.fail_count > 0) || (nimages == 0) || \
        (nimages > maximages)) {
        
        if (gemextn.fail_count > 0) {
            errmsg = gemextn.fail_count//" images were not found."
            status = 101
        } else if (nimages == 0) {
            errmsg = "No input images defined."
            status = 121
        } else if (nimages > maximages) {
            errmsg = "Maximum number of input images ("//str(maximages)//") \
                has been exceeded."
            status = 121
        }
        
        glogprint (l_logfile, "mipstokes", "status", type="error", 
            errno=status, str=errmsg, verbose+)
        goto clean
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
            glogprint (l_logfile, "mipstokes", "status", type="error",
                errno=status, str=errmsg, verbose+)
        }
    }

    # Now, do the same counting for the out file
    if (l_outimages != "")
        outputstr = l_outimages
    else if (l_outpref != "") {
        gemextn ("@"//tmpinimg, check="", process="none", index="", extname="",
            extversion="", ikparams="", omit="path", replace="",
            outfile=tmpoutimg, logfile=l_logfile, verbose=l_verbose)
        outputstr = l_outpref//"@"//tmpoutimg
    } else {
        status = 121
        glogprint (l_logfile, "mipstokes", "status", type="error",
            errno=status, str="Neither output image name nor output prefix \
            defined.", verbose+)
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
            errmsg = gemextn.fail_count//" images(s) already exist(s)."
            status = 102
        } else if (noutimages == 0) {
            errmsg = "No output images defined."
            status = 121
        } else if (noutimages != nimages) {
            errmsg = "Different number of input images ("//nimages//") and \
                output images ("//noutimages//")."
            status = 121
        }
        
        glogprint (l_logfile, "mipstokes", "status", type="error", 
            errno=status, str=errmsg, verbose+)
        goto clean
    } else {
        scanfile = tmpoutimg
        i = 0
        while (fscan (scanfile, filename) != EOF) {
            i += 1
            out[i] = filename//".fits"
        }
        scanfile = ""
        if (i != noutimages) {
            status = 99
            glogprint (l_logfile, "mipstokes", "status", type="error",
                errno=status, str=errmsg, verbose+)
            goto clean 
        }
    }
    delete (tmpinimg, verify-, >& "dev$null")
    delete (tmpoutimg, verify-, >& "dev$null")

    ########
    # Start the work.
    
    nskipped = 0
    for (i=1; i<=nimages; i=i+1) {
        
        glogprint (l_logfile, "mipstokes", "status", type="string",
            str="Processing file "//in[i]//" ...", verbose=l_verbose)
        
        # A few more checks are required
        #
        # First, input must not be already stacked
        
        tmpstatus = 0
        keyfound = ""
        hselect (in[i]//"[0]", "MIPSTOKES", yes) | scan (keyfound)
        if (keyfound != "") tmpstatus = 123
        keyfound = ""
        hselect (in[i]//"[0]", "MISTACK", yes) | scan (keyfound)
        if (keyfound != "") tmpstatus = 123
        # Shouldn't it check for MIPSTACK instead? -- KL
        
        if (tmpstatus != 0) {
            glogprint (l_logfile, "mipstokes", "status", type="warning",
                errno=status, str="Input file "//in[i]//" has already been \
                stacked.", verbose=l_verbose)
            nskipped += 1
            status = tmpstatus
            goto nextimage
        }
        
        # Which instrument?
        instrument = ""
        hselect (in[i]//"[0]", "INSTRUME", yes) | scan (instrument)
        instrument = strlwr(instrument)
        if (instrument != "michelle") {
            status = 123
            errmsg = "File "//in[i]//". Unknown or unsupported instrument ("//\
                instrument//")."
            glogprint (l_logfile, "mipstokes", "status", type="error",
                errno=status, str=errmsg, verbose+)
            nskipped += 1
            goto nextimage
        }
        
        # Input must have been PREPARE'd
        
        keyfound = ""
        hselect (in[i]//"[0]", "PREPARE", yes) | scan (keyfound)
        if (keyfound == "") {
            status = 121
            errmsg = "Image "//in[i]//" not "
            if (instrument == "michelle")   errmsg = errmsg//"MPREPAREd."
            glogprint (l_logfile, "mipstokes", "status", type="error", 
                errno=status, str=errmsg, verbose+)
            nskipped += 1
            goto nextimage
        }
        
        # Find the observation mode
        if (instrument == "michelle") {
            keyfound = ""
            hselect (in[i]//"[0]", "MODE", yes) | scan (keyfound)
            if (keyfound == "") {
                status = 131
                errmsg = "Header keyword MODE not found in PHU of "//in[i]//"."
                glogprint (l_logfile, "mipstokes", "status", type="error",
                    errno=status, str=errmsg, verbose+)
                nskipped += 1
                goto nextimage
            }
            
            # Change these according to the Michelle "MODE" keywords.
            # I am not sure whether there are other "non-destructive" (nd) 
            # modes, or whether these keywords are all correct.
            # -- KVolk
            
            if (keyfound == "chop-nod")     modeflag = 1
            else if (keyfound == "ndchop")  modeflag = 1
            else if (keyfound == "chop")    modeflag = 2
            else if (keyfound == "nod")     modeflag = 3
            else if (keyfound == "ndstare") modeflag = 4
            else if (keyfound == "stare")   modeflag = 4
            else {
                status = 132
                errmsg = "Unrecognized MODE ("//keyfound//") in PHU of "//\
                    in[i]//"."
                glogprint (l_logfile, "mipstokes", "status", type="error",
                    errno=status, str=errmsg, verbose+)
                nskipped += 1
                goto nextimage
            }
            
        } else {
            status = 99
            errmsg = "Unsupported instrument ("//instrument//").  Should have \
                been caught earlier."
            glogprint (l_logfile, "mipstokes", "status", type="error",
                errno=status, str=errmsg, verbose+)
            goto clean
        }
        
        # Get number of extensions and check dimension
        nextensions = 0
        fxhead (in[i], format="", long-, count-, >& tmpfile)
        system.tail (tmpfile, nlines=1) | scan (nextensions)
        delete (tmpfile, verify-, >& "dev$null")
        for (ext=1; ext<=nextensions; ext=ext+1) {
            tmpstatus = 0
            keyfound = ""
            hselect (in[i]//"["//str(ext)//"]", "i_naxis", yes) | \
                scan (keyfound)
            if (((modeflag == 1) || (modeflag == 2)) && (int(keyfound) != 3)) {
                tmpstatus = 132
                errmsg = "Extension "//str(ext)//" has "//keyfound//\
                    " dimensions.   It should be 3."                
            } else if (((modeflag == 3) || (modeflag == 4)) && \
                (int(keyfound) != 2)) {
                
                tmpstatus = 132
                errmsg = "Extension "//str(ext)//" has "//keyfound//\
                    " dimensions.   It should be 2."
            }
            
            if (tmpstatus != 0) {
                status = tmpstatus
                glogprint (l_logfile, "mipstokes", "status", type="error",
                    errno=status, str=errmsg, verbose+)
                nskipped += 1
                goto nextimage
            }
        }
        if (nextensions == 0) {
            status = 123
            errmsg = "No data extensions in file "//in[i]//"."
            glogprint (l_logfile, "mipstokes", "status", type="error",
                errno=status, str=errmsg, verbose+)
            nskipped += 1
            goto nextimage
        } else
            glogprint (l_logfile, "mipstokes", "engineering", type="string",
                str="Number of extensions is "//nextensions, verbose=l_verbose)
        
        # Check file structure (number of extensions and nods)
        if (modeflag == 1) {
            if (mod (nextensions, 8) != 0) {
                status = 123
                errmsg = "Number of extensions does not correspond to \
                    complete cycles of the polarimetry mode."
                glogprint (l_logfile, "mipstokes", "status", type="error",
                    errno=status, str=errmsg, verbose+)
                nskipped += 1
                goto nextimage
            }
            
            nnods = int (nextensions / 8)
            if (nnods == 1) {
                modeflag = 2
                
                status = 123
                errmsg = "Based on file structure, resetting observing mode."
                glogprint (l_logfile, "mipstokes", "status", type="warning",
                    errno=status, str=errmsg, verbose=l_verbose)
            } else if (mod (nnods, 2) != 0) {
                status = 123
                errmsg = "Number of extensions does not correspond to \
                    complete nodsets.  Removing last unmatched nod position."
                glogprint (l_logfile, "mipstokes", "status", type="warning",
                    errno=status, str=errmsg, verbose=l_verbose)
                nextensions -= 8
            }
        }
        
        if (modeflag == 1) {
            # in this initial routine I will NOT check for "BADNOD" 
            # flags....everything gets stacked up
            # -- KVolk
            
            # Make Stokes parameters for each pair of AB nod positions
            flag[1] = 0
            flag[2] = 0
            flag[3] = 0
            flag[4] = 0
            tmprefer1 = mktemp ("tmprefer1")
            tmprefer2 = mktemp ("tmprefer2")
            tmprefer3 = mktemp ("tmprefer3")
            tmprefer4 = mktemp ("tmprefer4")
            
            # Looping through the AB nod sets
            for (nod=1; nod<=nnods; nod=nod+2) {
                tmpims1a = mktemp ("tmpims1a")
                tmpims1b = mktemp ("tmpims1b")
                tmpims2a = mktemp ("tmpims2a")
                tmpims2b = mktemp ("tmpims2b")
                tmpims3a = mktemp ("tmpims3a")
                tmpims3b = mktemp ("tmpims3b")
                tmpims4a = mktemp ("tmpims4a")
                tmpims4b = mktemp ("tmpims4b")
                
                glogprint (l_logfile, "mipstokes", "engineering", 
                    type="string", str="Working on nods "//str(nod)//" and "//\
                    str(nod+1), verbose=l_verbose)
                
                for (ext=8*(nod-1)+1; ext<=8*nod; ext=ext+1) {
                    tmpreg1 = mktemp ("tmpreg1")
                    tmpreg2 = mktemp ("tmpreg2")
                    tmpworkimg = mktemp ("tmpworkimg")
                    
                    imgsec = "[*,*,"//framevalue//"]"
                    imcopy (in[i]//"["//ext//"]"//imgsec//".fits",
                        tmpreg1, verbose-)
                    imcopy (in[i]//"["//ext+8//"]"//imgsec//".fits",
                        tmpreg2, verbose-)

                    if (l_fl_stair) {
                        for (coord=1; coord<=301; coord=coord+20) {
                            imgsec = "["//str(coord)//":"//str(coord+19)//\
                                ",1:240]"

                            # tmpreg1
                            tmpimg1 = mktemp ("tmpimg1")
                            imstat (tmpreg1//imgsec, fields="midpt",
                                lower=INDEF, upper=INDEF, format-, nclip=0,
                                lsigma=3., usigma=3., binwidth=0.1,
                                cache-) | scan (chlevel)
                            imarith (tmpreg1//imgsec, "-", chlevel, tmpimg1,
                                title="", divzero=0., hparams="", pixtype="",
                                calctype="", verbose-, noact-)
                            imcopy (tmpimg1, tmpreg1//imgsec, verbose-)
                            imdelete (tmpimg1, verify-, >& "dev$null")

                            # tmpreg2
                            tmpimg1 = mktemp ("tmpimg1")
                            imstat (tmpreg2//imgsec, field="midpt",
                                lower=INDEF, upper=INDEF, format-, nclip=0,
                                lsigma=3., usigma=3., binwidth=0.1,
                                cache-) | scan (chlevel)
                            imarith (tmpreg2//imgsec, "-", chlevel, tmpimg1,
                                title="", divzero=0., hparams="", pixtype="",
                                calctype="", verbose-, noact-)
                            imcopy (tmpimg1, tmpreg2//imgsec, verbose-)
                            imdelete (tmpimg1, verify-, >& "dev$null")
                        }                            
                    }
                    
                    if (l_fl_register) {
                        images.immatch.xregister (tmpreg2, tmpreg1, 
                            regions=l_regions, shifts=tmpshifts,
                            output=tmpworkimg, append+, records="", coords="",
                            xlag=0, ylag=0, dxlag=0, dylag=0,
                            background="none", border=INDEF, loreject=INDEF, 
                            hireject=INDEF, apodize=0., filter="none",
                            correlation="discrete", xwindow=11, ywindow=11,
                            function="centroid", xcbox=11, ycbox=11,
                            interp_type="poly5", boundary_type="nearest", 
                            constant=0., interact-, verbose+)
                        imarith (tmpreg1, "+", tmpworkimg, tmpworkimg, 
                            title="", divzero=0., hparams="", pixtype="", 
                            calctype="", verbose-, noact-)
                    } else {
                        imarith (tmpreg1, "+", tmpreg2, tmpworkimg, title="",
                            divzero=0., hparams="", pixtype="", calctype="",
                            verbose-, noact-)
                    }
                    
                    imdelete (tmpreg1, verify-, >& "dev$null")
                    imdelete (tmpreg2, verify-, >& "dev$null")
                    
                    if (l_combine == "average") {
                        imarith (tmpworkimg, "/", "2.0", tmpworkimg, title="",
                            divzero=0., hparams="", pixtype="", calctype="",
                            verbose-, noact-)
                    }
                                        
                    k = ext - 8*int(ext/8)
                    if (k == 0) k = 8
                    
                    if ((k == 1) || (k == 4)) {     # pair 1
                        if (l_fl_register) {
                            if (flag[1] == 0) {
                                flag[1] = 1
                                imcopy (tmpworkimg, tmprefer1, verbose-)
                                imcopy (tmpworkimg, tmpims1a, verbose-)
                            } else {
                                tmpreg1 = mktemp ("tmpreg1")
                                images.immatch.xregister (tmpworkimg, 
                                    tmprefer1, regions=l_regions,
                                    shifts=tmpshifts, output=tmpreg1, append+,
                                    records="", coords="", xlag=0, ylag=0, 
                                    dxlag=0, dylag=0, background="none", 
                                    border=INDEF, loreject=INDEF,
                                    hireject=INDEF, apodize=0., filter="none",
                                    correlation="discrete", xwindow=11,
                                    ywindow=11, function="centroid", xcbox=11,
                                    ycbox=11, interp_type="poly5", 
                                    boundary_typ="nearest", constant=0., 
                                    interact-, verbose+) 
                                if (k == 1)
                                    imcopy (tmpreg1, tmpims1a, verbose-)
                                else if (k == 4)
                                    imcopy (tmpreg1, tmpims1b, verbose-)
                                imdelete (tmpreg1, verify-, >& "dev$null")
                            }
                        } else {
                            if (k == 1)
                                imcopy (tmpworkimg, tmpims1a, verbose-)
                            else if (k == 4)
                                imcopy (tmpworkimg, tmpims1b, verbose-)
                        }
                    } else if ((k == 2) || (k == 3)) {  # pair 2
                        if (l_fl_register) {
                            if (flag[2] == 0) {
                                flag[2] = 1
                                imcopy (tmpworkimg, tmprefer2, verbose-)
                                imcopy (tmpworkimg, tmpims2a, verbose-)
                            } else {
                                tmpreg1 = mktemp ("tmpreg1")
                                images.immatch.xregister (tmpworkimg, 
                                    tmprefer2, regions=l_regions, 
                                    shifts=tmpshifts, output=tmpreg1, append+, 
                                    records="", coords="", xlag=0, ylag=0, 
                                    dxlag=0, dylag=0, background="none", 
                                    border=INDEF, loreject=INDEF, 
                                    hireject=INDEF, apodize=0., filter="none", 
                                    correlation="discrete", xwindow=11, 
                                    ywindow=11, function="centroid", xcbox=11, 
                                    ycbox=11, interp_type="poly5", 
                                    boundary_typ="nearest", constant=0., 
                                    interact-, verbose+)
                                if (k == 2)
                                    imcopy (tmpreg1, tmpims2a, verbose-)
                                else if (k == 3)
                                    imcopy (tmpreg1, tmpims2b, verbose-)
                                imdelete (tmpreg1, verify-, >& "dev$null")
                            }
                        } else {
                            if (k == 2)
                                imcopy (tmpworkimg, tmpims2a, verbose-)
                            else if (k == 3)
                                imcopy (tmpworkimg, tmpims2b, verbose-)
                        }
                    } else if ((k == 5) || (k == 8)) {  # pair 3
                        if (l_fl_register) {
                            if (flag[3] == 0) {
                                flag[3] = 1
                                imcopy (tmpworkimg, tmprefer3, verbose-)
                                imcopy (tmpworkimg, tmpims3a, verbose-)
                            } else {
                                tmpreg1 = mktemp ("tmpreg1")
                                images.immatch.xregister (tmpworkimg, 
                                    tmprefer3, regions=l_regions, 
                                    shifts=tmpshifts, output=tmpreg1, append+, 
                                    records="", coords="", xlag=0, ylag=0, 
                                    dxlag=0, dylag=0, background="none", 
                                    border=INDEF, loreject=INDEF, 
                                    hireject=INDEF, apodize=0., filter="none", 
                                    correlation="discrete", xwindow=11, 
                                    ywindow=11, function="centroid", xcbox=11, 
                                    ycbox=11, interp_type="poly5", 
                                    boundary_typ="nearest", constant=0., 
                                    interact-, verbose+)
                                if (k == 5)
                                    imcopy (tmpreg1, tmpims3a, verbose-)
                                else if (k == 8)
                                    imcopy (tmpreg1, tmpims3b, verbose-)
                                imdelete (tmpreg1, verify-, >& "dev$null")
                            }
                        } else {
                            if (k == 5)
                                imcopy (tmpworkimg, tmpims3a, verbose-)
                            else if (k == 8)
                                imcopy (tmpworkimg, tmpims3b, verbose-)
                        }
                    } else if ((k == 6) || (k == 7)) {  # pair 4
                        if (l_fl_register) {
                            if (flag[4] == 0) {
                                flag[4] = 1
                                imcopy (tmpworkimg, tmprefer4, verbose-)
                                imcopy (tmpworkimg, tmpims4a, verbose-)
                            } else {
                                tmpreg1 = mktemp ("tmpreg1")
                                images.immatch.xregister (tmpworkimg, 
                                    tmprefer4, regions=l_regions, 
                                    shifts=tmpshifts, output=tmpreg1, append+, 
                                    records="", coords="", xlag=0, ylag=0, 
                                    dxlag=0, dylag=0, background="none", 
                                    border=INDEF, loreject=INDEF, 
                                    hireject=INDEF, apodize=0., filter="none", 
                                    correlation="discrete", xwindow=11, 
                                    ywindow=11, function="centroid", xcbox=11, 
                                    ycbox=11, interp_type="poly5", 
                                    boundary_typ="nearest", constant=0., 
                                    interact-, verbose+)
                                if (k == 6)
                                    imcopy (tmpreg1, tmpims4a, verbose-)
                                else if (k == 7)
                                    imcopy (tmpreg1, tmpims4b, verbose-)
                                imdelete (tmpreg1, verify-, >& "dev$null")
                            }
                        } else {
                            if (k == 6)
                                imcopy (tmpworkimg, tmpims4a, verbose-)
                            else if (k == 7)
                                imcopy (tmpworkimg, tmpims4b, verbose-)
                        }
                    } else {
                        status = 99
                        errmsg = "In AB loop, k="//k//", ext="//ext
                        glogprint (l_logfile, "mipstokes", "status",
                            type="error", errno=status, str=errmsg, verbose+)
                        goto clean
                    }                    

                    imdelete (tmpworkimg, verify-, >& "dev$null")
                    
                } # end for-loop extensions
                
                #-------------
                # Transform to Stokes parameters
                
                tmpUstokes = mktemp ("tmpUstokes")  # Stokes U map, extension 2
                tmpQstokes = mktemp ("tmpQstokes")  # Stokes Q map, extension 3
                tmpIstokes = mktemp ("tmpIstokes")  # Stokes I map, extension 1
                
                tmpimg1 = mktemp ("tmpimg1")
                tmpimg2 = mktemp ("tmpimg2")
                tmpimg3 = mktemp ("tmpimg3")
                tmpimg4 = mktemp ("tmpimg4")

                imcombine (tmpims1a//","//tmpims1b, tmpimg1, headers="",
                    bpmasks="", rejmasks="", nrejmasks="", expmasks="",
                    sigmas="", logfile="STDOUT", combine="average",
                    reject="none", project-, outtype="real", outlimits="",
                    offsets="none", masktype="none", maskvalue=0., blank=0.,
                    scale="none", zero="none", weight="none", statsec="",
                    expname="", lthreshold=INDEF, hthreshold=INDEF, nlow=1,
                    nhigh=1, nkeep=1, mclip=yes, lsigma=3., hsigma=3.,
                    rdnoise="0.", gain="1.", snoise="0.1", sigscale=0.1,
                    pclip=-0.5, grow=0., >& "dev$null")
                imcombine (tmpims2a//","//tmpims2b, tmpimg2, headers="",
                    bpmasks="", rejmasks="", nrejmasks="", expmasks="",
                    sigmas="", logfile="STDOUT", combine="average",
                    reject="none", project-, outtype="real", outlimits="",
                    offsets="none", masktype="none", maskvalue=0., blank=0.,
                    scale="none", zero="none", weight="none", statsec="",
                    expname="", lthreshold=INDEF, hthreshold=INDEF, nlow=1,
                    nhigh=1, nkeep=1, mclip=yes, lsigma=3., hsigma=3.,
                    rdnoise="0.", gain="1.", snoise="0.1", sigscale=0.1,
                    pclip=-0.5, grow=0., >& "dev$null")
                imcombine (tmpims3a//","//tmpims3b, tmpimg3, headers="",
                    bpmasks="", rejmasks="", nrejmasks="", expmasks="",
                    sigmas="", logfile="STDOUT", combine="average",
                    reject="none", project-, outtype="real", outlimits="",
                    offsets="none", masktype="none", maskvalue=0., blank=0.,
                    scale="none", zero="none", weight="none", statsec="",
                    expname="", lthreshold=INDEF, hthreshold=INDEF, nlow=1,
                    nhigh=1, nkeep=1, mclip=yes, lsigma=3., hsigma=3.,
                    rdnoise="0.", gain="1.", snoise="0.1", sigscale=0.1,
                    pclip=-0.5, grow=0., >& "dev$null")
                imcombine (tmpims4a//","//tmpims4b, tmpimg4, headers="",
                    bpmasks="", rejmasks="", nrejmasks="", expmasks="",
                    sigmas="", logfile="STDOUT", combine="average",
                    reject="none", project-, outtype="real", outlimits="",
                    offsets="none", masktype="none", maskvalue=0., blank=0.,
                    scale="none", zero="none", weight="none", statsec="",
                    expname="", lthreshold=INDEF, hthreshold=INDEF, nlow=1,
                    nhigh=1, nkeep=1, mclip=yes, lsigma=3., hsigma=3.,
                    rdnoise="0.", gain="1.", snoise="0.1", sigscale=0.1,
                    pclip=-0.5, grow=0., >& "dev$null")
                    
                imdelete (tmpims1a, verify-, >& "dev$null")
                imdelete (tmpims1b, verify-, >& "dev$null")
                imdelete (tmpims2a, verify-, >& "dev$null")
                imdelete (tmpims2b, verify-, >& "dev$null")
                imdelete (tmpims3a, verify-, >& "dev$null")
                imdelete (tmpims3b, verify-, >& "dev$null")
                imdelete (tmpims4a, verify-, >& "dev$null")
                imdelete (tmpims4b, verify-, >& "dev$null")
               
                if (l_fl_register) {
                    tmpreg1 = mktemp ("tmpreg1")
                    images.immatch.xregister (tmpimg2, tmpimg1, 
                        regions=l_regions, shifts=tmpshifts, output=tmpreg1,
                        append+, records="", coords="", xlag=0, ylag=0, 
                        dxlag=0, dylag=0, background="none", border=INDEF, 
                        loreject=INDEF, hireject=INDEF, apodize=0., 
                        filter="none", correlation="discrete", xwindow=11, 
                        ywindow=11, function="centroid", xcbox=11, ycbox=11,
                        interp_type="poly5", boundary_typ="nearest", 
                        constant=0., interact-, verbose+)
                    imdelete (tmpimg2, verify-, >& "dev$null")
                    imcopy (tmpreg1, tmpimg2, verbose-)
                    imdelete (tmpreg1, verify-, >& "dev$null")
                    
                    tmpreg1 = mktemp ("tmpreg1")
                    images.immatch.xregister (tmpimg3, tmpimg1, 
                        regions=l_regions, shifts=tmpshifts, output=tmpreg1, 
                        append+, records="", coords="", xlag=0, ylag=0, 
                        dxlag=0, dylag=0, background="none", border=INDEF, 
                        loreject=INDEF, hireject=INDEF, apodize=0., 
                        filter="none", correlation="discrete", xwindow=11,
                        ywindow=11, function="centroid", xcbox=11, ycbox=11,
                        interp_type="poly5", boundary_typ="nearest", 
                        constant=0., interact-, verbose+)
                    imdelete (tmpimg3, verify-, >& "dev$null")
                    imcopy (tmpreg1, tmpimg3, verbose-)
                    imdelete (tmpreg1, verify-, >& "dev$null")
                    
                    tmpreg1 = mktemp ("tmpreg1")
                    images.immatch.xregister (tmpimg4, tmpimg1, 
                        regions=l_regions, shifts=tmpshifts, output=tmpreg1, 
                        append+, records="", coords="", xlag=0, ylag=0, 
                        dxlag=0, dylag=0, background="none", border=INDEF, 
                        loreject=INDEF, hireject=INDEF, apodize=0., 
                        filter="none", correlation="discrete", xwindow=11, 
                        ywindow=11, function="centroid", xcbox=11, ycbox=11, 
                        interp_type="poly5", boundary_typ="nearest", 
                        constant=0., interact-, verbose+)
                    imdelete (tmpimg4, verify-, >& "dev$null")
                    imcopy (tmpreg1, tmpimg4, verbose-)
                    imdelete (tmpreg1, verify-, >& "dev$null")
                }
                
                imarith (tmpimg3, "-", tmpimg4, tmpUstokes, title="",
                    divzero=0., hparams="", pixtype="", calctype="", verbose-,
                    noact-)
                imarith (tmpimg1, "-", tmpimg2, tmpQstokes, title="",
                    divzero=0., hparams="", pixtype="", calctype="", verbose-,
                    noact-)
                imcombine (tmpimg1//","//tmpimg2//","//tmpimg3//","//tmpimg4,
                    tmpIstokes, headers="", bpmasks="", rejmasks="",
                    nrejmasks="", expmasks="", sigmas="", logfile="STDOUT",
                    combine="sum", reject="none", project-, outtype="real",
                    outlimits="", offsets="none", masktype="none", 
                    maskvalue=0., blank=0., scale="none", zero="none",
                    weight="none", statsec="", expname="", lthreshold=INDEF,
                    hthreshold=INDEF, nlow=1, nhigh=1, nkeep=1, mclip=yes,
                    lsigma=3., hsigma=3., rdnoise="0.", gain="1.",
                    snoise="0.1", sigscale=0.1, pclip=-0.5, grow=0., 
                    >& "dev$null") 
                imarith (tmpIstokes, "/", 2.0, tmpIstokes, title="", 
                    divzero=0., hparams="", pixtype="", calctype="", verbose-,
                    noact-)
                
                imdelete (tmpimg1, verify-, >& "dev$null")
                imdelete (tmpimg2, verify-, >& "dev$null")
                imdelete (tmpimg3, verify-, >& "dev$null")
                imdelete (tmpimg4, verify-, >& "dev$null")

                if (l_fl_mask) {
                    if (l_noise <= 0.) {
                        imstat (tmpIstokes//l_blankarea, fields="stddev",
                            lower=INDEF, upper=INDEF, nclip=0, lsigma=3.,
                            usigma=3., binwidth=0.1, format-, cache-) | \
                            scan (sigma)
                        glogprint (l_logfile, "mipstokes", "engineering",
                            type="string", str="Noise value (sigma) is \
                            calculated to be "//str(sigma), verbose=l_verbose)
                        floor = sigma * 282.0 / \
                            (4. * sqrt(nnods/2)*(l_threshold*l_threshold))
                            
                        # S/N required is 282 for 0.5% accuracy in the 
                        # final image, goes as the inverse square of the 
                        # accuracy required.  Per consistuant image the 
                        # S/N level needed is lower by the square-root of 
                        # the number of images being combined.
                        
                        glogprint (l_logfile, "mipstokes", "engineering",
                            type="string", str="Masking value is calculated \
                            as "//str(floor), verbose=l_verbose)
                        imstat (tmpIstokes, fields="max", format-, lower=INDEF,
                            upper=INDEF, nclip=0, lsigma=3., usigma=3.,
                            binwidth=0.1,  cache-) | scan (max)
                        if (max < floor) {
                            floor = max / 2.
                            glogprint (l_logfile, "mipstokes", "engineering",
                                type="warning", errno=121, str="Masking value \
                                is larger than the peak value ("//str(max)//\
                                ").  Correcting the value to "//str(floor),
                                verbose=l_verbose)
                        }
                    } else {
                        glogprint (l_logfile, "mipstokes", "engineering",
                            type="string", str="Masking value is set to be "//\
                            str(l_noise), verbose=l_verbose)
                        floor = l_noise
                    }
                    
                    # for this purpose, mask off all pixels with values 
                    # less than the floor value in the Stokes I image, then 
                    # use this as a mask for the Stokes U and Q images.  
                    # This masks out the negative beams.
                    
                    tmpimg1 = mktemp ("tmpimg1")
                    imreplace (tmpIstokes, 0., imaginary=0., lower=INDEF,
                        upper=floor, radius=0.)
                    imexpr ("5*int(a == 0.)", tmpimg1, tmpIstokes, dims="auto",
                        intype="auto", outtype="int", refim="auto", bwidth=0,
                        btype="nearest", bpixval=0., rangecheck+, verbose+,
                        exprdb="none", >& "dev$null")

                    imexpr ("(b > 0) ? 0 : a", tmpfile, tmpQstokes, tmpimg1,
                        dims="auto", intype="auto", outtype="auto", 
                        refim="auto", bwidth=0, btype="nearest", bpixval=0., 
                        rangecheck+, verbose+, exprdb="none", >& "dev$null")
                    imdelete (tmpQstokes, verify-, >& "dev$null")
                    imcopy (tmpfile, tmpQstokes, verbose-, >& "dev$null")

                    imdelete (tmpfile, verify-, >& "dev$null")
                    imexpr ("(b > 0) ? 0 : a", tmpfile, tmpUstokes, tmpimg1,
                        dims="auto", intype="auto", outtype="auto", 
                        refim="auto", bwidth=0, btype="nearest", bpixval=0., 
                        rangecheck+, verbose+, exprdb="none", >& "dev$null")
                    imdelete (tmpUstokes, verify-, >& "dev$null")
                    imcopy (tmpfile, tmpUstokes, verbose-, >& "dev$null")
                    imdelete (tmpfile, verify-, >& "dev$null")

                    imdelete (tmpimg1, verify-, >& "dev$null")
                }
                                
                # Calculate tmpimg1, the Ip component
                tmpimg1 = mktemp ("tmpimg1")
                tmpimg2 = mktemp ("tmpimg2")
                imcopy (tmpUstokes, tmpimg1, verbose-)
                imcopy (tmpQstokes, tmpimg2, verbose-)
                imfunction (tmpimg1, tmpimg1, "square", verbose-)
                imfunction (tmpimg2, tmpimg2, "square", verbose-)
                imarith (tmpimg1, "+", tmpimg2, tmpimg1, title="", divzero=0.,
                    hparams="", pixtype="", calctype="", verbose-, noact-)
                imfunction (tmpimg1, tmpimg1, "sqrt", verbose-)
                imdelete (tmpimg2, verify-, >& "dev$null")
                
                # Calculate tmpimg2, the polarization angle in degrees
                tmpimg2 = mktemp ("tmpimg2")
                imexpr ("deg(atan2(a,b))/2.", tmpimg2, tmpUstokes, tmpQstokes,
                    dims="auto", intype="auto", outtype="auto", refim="auto",
                    bwidth=0, btype="nearest", bpixval=0., rangecheck+,
                    verbose-, exprdb="none")
                
                # Calculate tmpimg3, the unpolarized intensity
                tmpimg3 = mktemp ("tmpimg3")
                imarith (tmpIstokes, "-", tmpimg1, tmpimg3, title="",
                    divzero=0., hparams="", pixtype="", calctype="", verbose-,
                    noact-)
                
                # Calculate tmpimg4, the fractional polarization value in %
                tmpimg4 = mktemp ("tmpimg4")
                imdivide (tmpimg1, tmpIstokes, tmpimg4, title="", constant=0.,
                    rescale="norescale", mean="1", verbose-)
                imarith (tmpimg4, "*", 100.0, tmpimg4, title="", divzero=0.,
                    hparams="", pixtype="", calctype="", verbose-, noact-)
                # mask out values over 100% or under -100% as obviously 
                # due to noise
                imreplace (tmpimg4, 0., imaginary=0., lower=100.001,
                    upper=INDEF, radius=0.)
                imreplace (tmpimg4, 0., imaginary=0., lower=INDEF,
                    upper=-100.001, radius=0.)
                
                # At this point:
                #   tmpimg1 is the Ip component
                #   tmpimg2 is the polarization angle in degrees
                #   tmpimg3 is the unpolarized intensity
                #   tmpimg4 is the fractional polarization value in %
                
                if (nod == 1) {
                    wmef (tmpIstokes, out[i], extname="SCI", phu=in[i]//"[0]",
                        verbose-, >& "dev$null")
                    outextn = 1
                } else {
                    fxinsert (tmpIstokes, out[i]//"["//str(ext)//"]", "", 
                        verbose-)
                    outextn += 1
                }
                firstextn = outextn
              
                fxinsert (tmpUstokes, out[i]//"["//str(outextn)//"]", "", 
                    verbose-)
                outextn += 1
                fxinsert (tmpQstokes, out[i]//"["//str(outextn)//"]", "", 
                    verbose-)
                outextn += 1
                fxinsert (tmpimg3, out[i]//"["//str(outextn)//"]", "", 
                    verbose-)
                outextn += 1
                fxinsert (tmpimg1, out[i]//"["//str(outextn)//"]", "", 
                    verbose-)
                outextn += 1
                fxinsert (tmpimg4, out[i]//"["//str(outextn)//"]", "", 
                    verbose-)
                outextn += 1
                fxinsert (tmpimg2, out[i]//"["//str(outextn)//"]", "", 
                    verbose-)
                outextn += 1
                lastextn = outextn                
                
                imdelete (tmpimg1, verify-, >& "dev$null")
                imdelete (tmpimg2, verify-, >& "dev$null")
                imdelete (tmpimg3, verify-, >& "dev$null")
                imdelete (tmpimg4, verify-, >& "dev$null")
                imdelete (tmpQstokes, verify-, >& "dev$null")
                imdelete (tmpUstokes, verify-, >& "dev$null")
                imdelete (tmpIstokes, verify-, >& "dev$null")
                
                # Propagate WCS information to the stacked images
                for (ext=firstextn; ext<=lastextn; ext=ext+1) {
                    gemhedit (out[i]//"["//str(ext)//"]", "WAXMAP01", "", "",
                        delete=yes)
                    
                    keyfound = ""
                    hselect (in[i]//"[0]", "CTYPE1", yes) | scan (keyfound)
                    if (keyfound != "") {
                        gemhedit (out[i]//"["//str(ext)//"]", "WCSAXES", 2,
                            "Number of WCS axes in the image", delete-)
                        gemhedit (out[i]//"["//str(ext)//"]", "CTYPE1",
                            keyfound, "R.A. in tangent plane projection",
                            delete-)
                    }
                    
                    keyfound = ""
                    hselect (in[i]//"[0]", "CRPIX1", yes) | scan (keyfound)
                    if (keyfound != "")
                        gemhedit (out[i]//"["//str(ext)//"]", "CRPIX1",
                            keyfound, "Ref pix of axis 1", delete-)
                    
                    keyfound = ""
                    hselect (in[i]//"[0]", "CRVAL1", yes) | scan (keyfound)
                    if (keyfound != "")
                        gemhedit (out[i]//"["//str(ext)//"]", "CRVAL1",
                            keyfound, "RA at Ref pix in decimal degrees",
                            delete-)
                    
                    keyfound = ""
                    hselect (in[i]//"[0]", "CTYPE2", yes) | scan (keyfound)
                    if (keyfound != "")
                        gemhedit (out[i]//"["//str(ext)//"]", "CTYPE2",
                            keyfound, "DEC. in tangent plane projection",
                            delete-)
                    
                    keyfound = ""
                    hselect (in[i]//"[0]", "CRPIX2", yes) | scan (keyfound)
                    if (keyfound != "")
                        gemhedit (out[i]//"["//str(ext)//"]", "CRPIX2",
                            keyfound, "Ref pix of axis 2", delete-)
                    
                    keyfound = ""
                    hselect (in[i]//"[0]", "CRVAL2", yes) | scan (keyfound)
                    if (keyfound != "")
                        gemhedit (out[i]//"["//str(ext)//"]", "CRVAL2",
                            keyfound, "DEC at Ref pix in decimal degrees",
                            delete-)
                    
                    keyfound = ""
                    hselect (in[i]//"[0]", "CD1_1", yes) | scan (keyfound)
                    if (keyfound != "")
                        gemhedit (out[i]//"["//str(ext)//"]", "CD1_1",
                            keyfound, "WCS matrix element 1 1", delete-)
                    
                    keyfound = ""
                    hselect (in[i]//"[0]", "CD1_2", yes) | scan (keyfound)
                    if (keyfound != "")
                        gemhedit (out[i]//"["//str(ext)//"]", "CD1_2",
                            keyfound, "WCS matrix element 1 2", delete-)
                    
                    keyfound = ""
                    hselect (in[i]//"[0]", "CD2_1", yes) | scan (keyfound)
                    if (keyfound != "")
                        gemhedit (out[i]//"["//str(ext)//"]", "CD2_1",
                            keyfound, "WCS matrix element 2 1", delete-)
                    
                    keyfound = ""
                    hselect (in[i]//"[0]", "CD2_2", yes) | scan (keyfound)
                    if (keyfound != "")
                        gemhedit (out[i]//"["//str(ext)//"]", "CD2_2",
                            keyfound, "WCS matrix element 2 2", delete-)
                    
                    keyfound = ""
                    hselect (in[i]//"[0]", "RADECSYS", yes) | scan (keyfound)
                    if (keyfound != "")
                        gemhedit (out[i]//"["//str(ext)//"]", "RADECSYS",
                            keyfound, "R.A./DEC. coordinate system reference",
                            delete-)
                    
                    gemhedit (out[i]//"["//str(ext)//"]", "EXTNAME", "SCI",
                        "Extension name", delete-)
                    gemhedit (out[i]//"["//str(ext)//"]", "EXTVER", ext,
                        "Extension number", delete-)
                
                }
                                
            } # end for-loop AB nod
            
        
            imdelete (tmprefer1, verify-, >& "dev$null")
            imdelete (tmprefer2, verify-, >& "dev$null")
            imdelete (tmprefer3, verify-, >& "dev$null")
            imdelete (tmprefer4, verify-, >& "dev$null")
        
        } else {
            status = 122
            errmsg = "Observing mode not yet supported."
            glogprint (l_logfile, "mipstokes", "status", type="error",
                errno=status, str=errmsg, verbose+)
            nskipped += 1
            goto nextimage
        }
        
        # Time stamps
        gemdate ()
        gemhedit (out[i]//"[0]", "MIPSTOKE", gemdate.outdate,
            "UT Time stamp for MIPSTOKES", delete-)
        gemhedit (out[i]//"[0]", "GEM-TLM", gemdate.outdate,
            "UT Last modification with GEMINI", delete-)
        
        
        # Jump to here if there is a problem
nextimage:
        ;        
        
    }

clean:
    scanfile = ""
    delete (tmpinimg, verify-, >& "dev$null")
    delete (tmpoutimg, verify-, >& "dev$null")
    if (l_fl_register) {
        delete (tmpshifts, verify-, >& "dev$null")
    }
    
    if (status == 0)
        glogclose (l_logfile, "mipstokes", fl_success+, verbose=l_verbose)
    else
        glogclose (l_logfile, "mipstokes", fl_success-, verbose=l_verbose)

exitnow:
    ;

end

# Copyright(c) 2003-2011 Association of Universities for Research in Astronomy, Inc.

procedure mprepare (inimages)

# Procedure to take raw Michelle data, with n MEF extensions (where n is the
# number of nodsets), and reformat the data so that it is consistent with
# the higher level data reduction image format.

# Data format is flipped from difference frame (ext 1), nod A (ext 2), nod
# B (ext 3) to nod A (ext 1), nod B (ext2), difference frame (ext 3).

# Version Sep  3, 2003  TB, first write
#         Oct 29, 2003  KL  IRAF 2.12 - new parameters
#                             hedit: addonly
#         Oct 20, 2005  KV added polarimetry and "rescue" support
#         Dec  5, 2005  KV changed polarimetry from a flag to checking the 
#                         Michelle filter name for "GRID".

char    inimages    {prompt="Input Michelle image(s)"}          # OLDP-1-input-primary-single-prefix=n
char    outimages   {"",prompt="Output image(s)"}               # OLDP-1-output
char    outprefix   {"m",prompt="Prefix for output image(s)"}   # OLDP-4
char    rawpath     {"",prompt="Path for input raw images"}     # OLDP-4
bool    fl_rescue   {no,prompt="Try to rescue an incomplete observation"}   # OLDP-1
char    logfile     {"",prompt="Logfile"}                       # OLDP-1
bool    verbose     {yes,prompt="Verbose"}                      # OLDP-4
int     status      {0,prompt="Exit status (0=good)"}           # OLDP-4
struct  *scanfile   {"",prompt="Internal use only"}             # OLDP-4

begin
    char    l_inimages = ""
    char    l_rawpath = ""
    char    l_outimages = ""
    char    l_outprefix = ""
    char    l_logfile = ""
    bool    l_fl_rescue, l_verbose
    
    char    errmsg, filename, outputstr
    char    tmpinimg, tmpoutimg, tmpfile, tmplog
    char    tmpim1, tmpim2, tmpim3, tmpim4
    char    in[500], out[500]
    char    paramstr, keyfound
    int     junk, nerror
    int     i, j, k, nimages, noutimages, maxfiles, n_j, modeflag
    int     next_inp, next_sci, count, next_tot
    bool    alreadyfixed[500], polarimetry

    status = 0
    maxfiles = 500
    nerror = 0
    
    cache ("imgets", "gemdate")

    # set the local variables
    junk = fscan (inimages, l_inimages)
    junk = fscan (rawpath, l_rawpath)
    junk = fscan (outimages, l_outimages)
    junk = fscan (outprefix, l_outprefix)
    l_fl_rescue = fl_rescue
    junk = fscan (logfile, l_logfile)
    l_verbose = verbose

    # open temp files
    tmpinimg = mktemp ("tmpinimg")
    tmpoutimg = mktemp ("tmpoutimg")
    tmplog = mktemp ("tmplog")    
    tmpfile = mktemp ("tmpfile")

    # Create the list of parameter/value pairs.  One pair per line.
    # All lines combined into one string.  Line delimiter is '\n'.
    paramstr =  "inimages       = "//inimages.p_value//"\n"
    paramstr += "outimages      = "//outimages.p_value//"\n"
    paramstr += "outprefix      = "//outprefix.p_value//"\n"
    paramstr += "rawpath        = "//rawpath.p_value//"\n"
    paramstr += "fl_rescue      = "//fl_rescue.p_value//"\n"
    paramstr += "logfile        = "//logfile.p_value//"\n"
    paramstr += "verbose        = "//verbose.p_value
    
    # Assign a logfile name if not given. Open logfile and start log.
    # Write parameter/value pairs ("paramstr") to log.
    gloginit (l_logfile, "mprepare", "midir", paramstr, fl_append+,
        verbose=l_verbose)
    if (gloginit.status != 0) {
        status = gloginit.status
        goto exitnow
    }
    l_logfile = gloginit.logfile.p_value

    #----------------------------------------------------------------------
    # Check the rawpath name for a final /
    if (l_rawpath != "") {
        if (substr(l_rawpath,(strlen(l_rawpath)),(strlen(l_rawpath))) != "/")
            l_rawpath = l_rawpath//"/"
    }

    #-----------------------------------------------------------------------
    # Load up arrays of input name lists

    gemextn (l_inimages, check="exists,mef", process="none", index="",
        extname="", extversion="", ikparam="", omit="kernel,exten",
        replace="^%%" // l_rawpath // "%", outfile=tmpinimg, logfile=l_logfile,
        verbose=l_verbose)
    nimages = gemextn.count
    
    if ((gemextn.fail_count > 0) || (nimages == 0) || \
        (nimages > maxfiles)) {
        
        if (gemextn.fail_count > 0) {
            errmsg = gemextn.fail_count//" images were not found."
            status = 101
        } else if (nimages == 0) {
            errmsg = "No input images defined."
            status = 121
        } else if (nimages > maxfiles) {
            errmsg = "Maximum number of input images ("//str(maxfiles)//") \
                has been exceeded."
            status = 121
        }
        
        glogprint (l_logfile, "mprepare", "status", type="error", errno=status,
            str=errmsg, verbose+)
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
            glogprint (l_logfile, "mprepare", "status", type="error",
                errno=status, str=errmsg, verbose+)
        }
    }
    
    # Check if some of the inputs have already been MPREPAREd
    for (i=1; i<=nimages; i+=1) {
        keyfound = ""
        hselect (in[i]//"[0]", "MPREPARE", yes) | scan (keyfound)
        if (keyfound == "") {
            alreadyfixed[i] = no
        } else {
            alreadyfixed[i] = yes
            glogprint (l_logfile, "mprepare", "status", type="warning",
                errno=0, str="Image "//in[i]//" already prepared with MPREPARE",
                verbose=l_verbose)
            glogprint (l_logfile, "mprepare", "status", type="string",
                str="              Data scaling, correction for number of \
                reads, and header", verbose=l_verbose)
            glogprint (l_logfile, "mprepare", "status", type="string",
                str="              updating not performed.", verbose=l_verbose)
        }
    }

    ##############
    # Now for the output images

    # Note: For simplicity, the already MPREPAREd images names are modified
    #       like for the others (it's simpler when ensuring that the number of
    #       output names is equal to the number of inputs).  Since those files
    #       will simply be copied over in the loop (alreadyfixed[i]=yes), 
    #       it works out fine.
    # 
    if (l_outimages != "")
        outputstr = l_outimages
    else if (l_outprefix != "") {
        gemextn ("@"//tmpinimg, check="", process="none", index="", extname="",
            extversion="", ikparams="", omit="path", replace="",
            outfile=tmpoutimg, logfile=l_logfile, verbose=l_verbose)
        outputstr = l_outprefix//"@"//tmpoutimg
    } else {
        status = 121
        glogprint (l_logfile, "mprepare", "status", type="error",
            errno=status, str="Neither output image name nor output prefix \
            defined.", verbose+)
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
        
        glogprint (l_logfile, "mprepare", "status", type="error",
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
            glogprint (l_logfile, "mprepare", "status", type="error",
                errno=status, str=errmsg, verbose+)
            goto clean
        }
    }
    delete (tmpinimg, ver-, >& "dev$null")
    delete (tmpoutimg, ver-, >& "dev$null")

    ############################3
    # Starting the real work

    glogprint (l_logfile, "mprepare", "status", type="string",
        str="Processing "//nimages//" files",verbose=l_verbose)

    scanfile = ""
    i = 1
    while (i <= nimages) {

        in[i] = in[i]//".fits"

        if (alreadyfixed[i] == yes) {
            fparse (in[i], ver-)
            if (no == access (fparse.root//fparse.extension)) {
                glogprint (l_logfile, "mprepare", "task", type="string",
                    str="Copying already prepared file "//in[i],
                    verbose=l_verbose)
                glogprint (l_logfile, "mprepare", "task", type="string",
                    str="to current directory.", verbose=l_verbose)
                copy (in[i], ".", verbose=l_verbose, >& tmplog)
                glogprint (l_logfile, "mprepare", "task", type="file",
                    str=tmplog, verbose=l_verbose)
                delete (tmplog, ver-, >& "dev$null")
            } else {
                glogprint (l_logfile, "mprepare", "task", type="string",
                    str="Prepared file "//in[i]//" already in current \
                    directory.", verbose=l_verbose)
            }
            goto nextimage
        }
                   
        #--------------------------------
        # Determine the number of extensions/nod sets.

        # Is this polarimetry?
        j = 0
        k = 0
        keyfound = ""
        hselect (in[i]//"[0]", "FILTERA", yes) | scan (keyfound)
        if (substr(keyfound, 1, 4) == "Grid")
            j = 1
        
        keyfound = ""
        hselect (in[i]//"[0]", "FILTERB", yes) | scan (keyfound)
        if (substr(keyfound, 1, 4) == "Grid")
            k = 1
            
        if ((j > 0) || (k > 0))     polarimetry = yes
        else                        polarimetry = no

        # Determine mode
        keyfound = ""
        hselect (in[i]//"[0]", "MODE", yes) | scan (keyfound)
        if (keyfound == "") {
            status = 131
            errmsg = "Could not find the MODE from the primary header."
            glogprint (l_logfile, "mprepare", "status", type="error",
                errno=status, str=errmsg, verbose+)
                
            nerror = nerror+1
            goto nextimage
        }

        if      (keyfound == "chop-nod")    modeflag = 1
        else if (keyfound == "ndchop")      modeflag = 1
        else if (keyfound == "chop")        modeflag = 2
        else if (keyfound == "nod")         modeflag = 3
        else if (keyfound == "ndstare")     modeflag = 4
        else if (keyfound == "stare")       modeflag = 4
        else {
            status = 132
            errmsg =  "Unrecognised MODE ("//imgets.value//") in the primary \
                header."
            glogprint (l_logfile, "mprepare", "status", type="error",
                errno=status, str=errmsg, verbose+)
                
            nerror = nerror+1
            goto nextimage
        }

        # Determine the number of extensions
        keyfound = ""
        hselect (in[i]//"[0]", "NUMEXT", yes) | scan (keyfound)
        
        if (l_fl_rescue) {
            if (keyfound == "") {
                glogprint (l_logfile, "mprepare", "status", type="warning",
                    errno=123, str="NUMEXT keyword missing from the \
                    primary header; recovery will be attempted.", verbose+)
            }
            gemextn (in[i], check="exists", process="expand", index="1-",
                extname="", extversion="", ikparams="", omit="", replace="",
                outfile="dev$null", logfile="", verbose=no)
            next_inp = gemextn.count
            
            glogprint (l_logfile, "mprepare", "engineering", type="string",
                str="              Number of extensions found is "//\
                str(next_inp), verbose=l_verbose) 
            
            if (modeflag == 1) {
                next_inp = 4*int(next_inp/4)
                glogprint (l_logfile, "mprepare", "engineering", type="string",
                    str="   Chop/Nod mode: number of extension \
                    to be used is "//str(next_inp), verbose=l_verbose)
            }
                     
        } else if (keyfound == "") {
            status = 123
            errmsg = "NUMEXT keyword missing from the primary header"
            glogprint (l_logfile, "mprepare", "status", type="error",
                errno=status, str=errmsg, verbose+)
                
            nerror = nerror+1
            goto nextimage
            
        } else
            next_inp = int(keyfound)
        

        if (polarimetry) {
            glogprint (l_logfile, "mprepare", "engineering", type="string",
                str="   Chop/Nod mode: number of extensions in \
                raw data file is  "//str(next_inp), verbose=l_verbose)
            glogprint (l_logfile, "mprepare", "engineering", type="string",
                str="   Polarimetry mode: number of extensions \
                to be created is "//str(8*next_inp), verbose=l_verbose)
        }
                      
        if ((modeflag == 4) && (next_inp != 1)) {
            status = 123
            errmsg = "Number of extensions in the primary header not \
                consistent with the MODE (stare)."
            glogprint (l_logfile, "mprepare", "status", type="error",
                errno=status, str=errmsg, verbose+)
            nerror = nerror+1
            goto nextimage
            
        } else if ((modeflag == 2) && (next_inp != 1)) {
            status = 123
            errmsg = "Number of extensions in the primary header not consistent \
                with the MODE (chop)."
            glogprint (l_logfile, "mprepare", "status", type="error",
                errno=status, str=errmsg, verbose+)
            nerror = nerror+1
            goto nextimage
            
        }

        #----------------------------------------------------------------------
        # start output
        glogprint (l_logfile, "mprepare", "visual", type="visual", 
            vistype="shortdash", verbose=l_verbose )
        glogprint (l_logfile, "mprepare", "task", type="string",
            str="    input file ->    output file", verbose=l_verbose )

        #Start on the main mprepare steps - rearrange the order of images 
        #within the extensions, and convert the 4D format to 3D format.


        #Copy PHU into new file
        fxcopy (in[i], out[i], "0", new_file+, verbose+, >& tmplog)
        glogprint (l_logfile, "mprepare", "task", type="file", str=tmplog,
            verbose=l_verbose )
        delete (tmplog, ver-, >& "dev$null")

        n_j = 1
        next_sci = 0

        while (n_j <= next_inp) {
            # Create tmp FITS file name used within this loop
            tmpim1 = mktemp("tmpim1")
            tmpim2 = mktemp("tmpim2")
            tmpim3 = mktemp("tmpim3")
            tmpim4 = mktemp("tmpim4")

            if (modeflag == 4) {
                imcopy (in[i]//"["//n_j//"]"//"[*,*,1,1]", tmpim1, ver-)      
                fxinsert (tmpim1//".fits", out[i]//"["//next_sci//"]",
                    groups="", ver-)
                glogprint (l_logfile, "mprepare", "engineering", type="string",
                    str="Prepared image for stare "//i, verbose=l_verbose )
                n_j = n_j+1
		
            } else if (modeflag == 3) {
                imcopy (in[i]//"["//n_j//"]"//"[*,*,1,1]", tmpim1, ver-)      
                fxinsert (tmpim1//".fits", out[i]//"["//next_sci//"]",
                    groups="", ver-)
                glogprint (l_logfile, "mprepare", "engineering", type="string",
                    str="Prepared image for nod "//n_j, verbose=l_verbose)
                n_j = n_j+1
                imarith (in[i]//"["//n_j//"]"//"[*,*,1,1]", "*", "-1",
                    tmpim2, ver-)
                fxinsert (tmpim2//".fits", out[i]//"["//next_sci//"]",
                    groups="", ver-)
                glogprint (l_logfile, "mprepare", "engineering", type="string",
                    str="Prepared image for nod "//n_j, verbose=l_verbose)
                n_j = n_j+1
                if  (n_j < next_inp) {
                    imarith (in[i]//"["//n_j//"]"//"[*,*,1,1]", "*",
                        "-1", tmpim3, ver-)
                    fxinsert (tmpim3//".fits", out[i]//"["//next_sci//"]",
                        groups="", ver-)
                    glogprint (l_logfile, "mprepare", "engineering",
                        type="string", str="Prepared image for nod "//n_j,
                        verbose=l_verbose)
                    n_j = n_j+1
                    imcopy (in[i]//"["//n_j//"]"//"[*,*,1,1]", 
                        tmpim4, ver-)      
                    fxinsert (tmpim4//".fits", out[i]//"["//next_sci//"]",
                        groups="", ver-)
                    glogprint (l_logfile, "mprepare", "engineering",
                        type="string", str="Prepared image for nod "//n_j,
                        verbose=l_verbose)
                    n_j = n_j+1
                }

            } else if (((modeflag == 2) || (modeflag ==1)) && (next_inp == 1)) {
                #Re-arrange A nod frame
                imcopy (in[i]//"["//n_j//"]"//"[*,*,*,1]", tmpim1, ver-)      
                imcopy (in[i]//"["//n_j//"]"//"[*,*,2,1]",
                    tmpim1//"[*,*,1]", ver-)      
                imcopy (in[i]//"["//n_j//"]"//"[*,*,3,1]",
                    tmpim1//"[*,*,2]", ver-)      
                imcopy (in[i]//"["//n_j//"]"//"[*,*,1,1]",
                    tmpim1//"[*,*,3]", ver-)      
                fxinsert (tmpim1//".fits", out[i]//"["//next_sci//"]",
                    groups="", ver-)
                n_j = n_j+1
                glogprint (l_logfile, "mprepare", "engineering", type="string",
                    str="Rearranged image for chop frame "//i,
                    verbose=l_verbose)
		
            } else if ((modeflag == 1) && (next_inp != 1)) {
                if (polarimetry) {
                    count = 8
                } else {
                    count = 1
                }
                for (k=1; k <= count; k=k+1) {
                    #Re-arrange A nod frame
                    imcopy (in[i]//"["//n_j//"]"//"[*,*,*,"//str(k)//"]",
                        tmpim1, ver-)      
                    imcopy (in[i]//"["//n_j//"]"//"[*,*,2,"//str(k)//"]",
                        tmpim1//"[*,*,1]", ver-)      
                    imcopy (in[i]//"["//n_j//"]"//"[*,*,3,"//str(k)//"]",
                        tmpim1//"[*,*,2]", ver-)      
                    imcopy (in[i]//"["//n_j//"]"//"[*,*,1,"//str(k)//"]",
                        tmpim1//"[*,*,3]", ver-)      
                    fxinsert (tmpim1//".fits", out[i]//"["//next_sci//"]",
                        groups="", ver-)
                    if (polarimetry) {
                        glogprint (l_logfile, "mprepare", "engineering",
                            type="string", str="Rearranged image for nod "//\
                            n_j//" component "//str(k), verbose=l_verbose)
                    } else {
                        glogprint (l_logfile, "mprepare", "engineering",
                            type="string", str="Rearranged image for nod "//\
                            n_j, verbose=l_verbose)
                    } 
                    imdelete (tmpim1, ver-, >& "dev$null")
                    tmpim1 = mktemp ("tmpim1")
                    next_sci = next_sci+1
                }
                n_j = n_j+1
                for (k=1; k <=count; k=k+1) {
                    imdelete (tmpim2, ver-, >& "dev$null")
                    tmpim2 = mktemp("tmpim2")
                    #Re-arrange and invert beam for B nod frame
                    imcopy (in[i]//"["//n_j//"]"//"[*,*,*,"//str(k)//"]",
                        tmpim1, ver-)      
                    imcopy (in[i]//"["//n_j//"]"//"[*,*,2,"//str(k)//"]",
                        tmpim1//"[*,*,2]", ver-)      
                    imcopy (in[i]//"["//n_j//"]"//"[*,*,3,"//str(k)//"]",
                        tmpim1//"[*,*,1]", ver-)      
                    imarith (in[i]//"["//n_j//"]"//"[*,*,1,"//str(k)//"]",
                        "*", "-1", tmpim2, ver-)
                    imcopy (tmpim2, tmpim1//"[*,*,3]", ver-)      
                    fxinsert (tmpim1//".fits", out[i]//"["//next_sci//"]",
                        groups="", ver-)
                    if (polarimetry) {
                        glogprint (l_logfile, "mprepare", "engineering",
                            type="string", str="Rearranged image for nod "//\
                            n_j//" component "//str(k), verbose=l_verbose)
                    } else {
                        glogprint (l_logfile, "mprepare", "engineering",
                            type="string", str="Rearranged image for nod "//\
                            n_j, verbose=l_verbose)
                    }
                    imdelete (tmpim1, ver-, >& "dev$null")
                    tmpim1 = mktemp("tmpim1")
                    next_sci = next_sci+1
                }
                n_j = n_j+1

                # Take into account the fact that there might be a 
                # non-divisible by 4 number of nod pairs
                
                if (n_j < next_inp) {
                    for (k=1; k <= count; k=k+1) {
                        imdelete (tmpim1, ver-, >& "dev$null")
                        imdelete (tmpim2//".fits", ver-, >& "dev$null")
                        tmpim1 = mktemp ("tmpim1")
                        tmpim2 = mktemp ("tmpim2")

                        #Re-arrange and invert beam for B nod frame
                        imcopy (in[i]//"["//n_j//"]"//"[*,*,*,"//str(k)//"]",
                            tmpim1, ver-)      
                        imcopy (in[i]//"["//n_j//"]"//"[*,*,2,"//str(k)//"]",
                            tmpim1//"[*,*,2]", ver-)      
                        imcopy (in[i]//"["//n_j//"]"//"[*,*,3,"//str(k)//"]",
                            tmpim1//"[*,*,1]", ver-)      
                        imarith (in[i]//"["//n_j//"]"//"[*,*,1,"//str(k)//"]",
                            "*", "-1", tmpim2, ver-)
                        imcopy (tmpim2, tmpim1//"[*,*,3]", ver-)      
                        fxinsert (tmpim1//".fits",
                            out[i]//"["//next_sci//"]", groups="", ver-)
                        if (polarimetry) {
                            glogprint (l_logfile, "mprepare", "engineering",
                                type="string", str="Rearranged image for \
                                nod "//n_j//" component "//str(k),
                                verbose=l_verbose)
                        } else {
                            glogprint (l_logfile, "mprepare", "engineering",
                                type="string", str="Rearranged image for \
                                nod "//n_j, verbose=l_verbose)
                        }
                        next_sci = next_sci+1
                    }

                    #Re-arrange A nod frame
                    n_j = n_j+1
                    imdelete (tmpim1//".fits", ver-, >& "dev$null")
                    tmpim1 = mktemp("tmpim1")
                    for (k=1; k <= count; k=k+1) {
                        imcopy (in[i]//"["//n_j//"]"//"[*,*,*,"//str(k)//"]",
                            tmpim1, ver-)      
                        imcopy (in[i]//"["//n_j//"]"//"[*,*,2,"//str(k)//"]",
                            tmpim1//"[*,*,1]", ver-)
                        imcopy (in[i]//"["//n_j//"]"//"[*,*,3,"//str(k)//"]",
                            tmpim1//"[*,*,2]", ver-)      
                        imcopy (in[i]//"["//n_j//"]"//"[*,*,1,"//str(k)//"]",
                            tmpim1//"[*,*,3]", ver-)      
                        fxinsert (tmpim1//".fits",
                            out[i]//"["//next_sci//"]", groups="", ver-)
                        if (polarimetry) {
                            glogprint (l_logfile, "mprepare", "engineering",
                                type="string", str="Rearranged image for \
                                nod "//n_j//" component "//str(k),
                                verbose=l_verbose)
                        } else {
                            glogprint (l_logfile, "mprepare", "engineering",
                                type="string", str="Rearranged image for \
                                nod "//n_j, verbose=l_verbose)
                        }
                        imdelete (tmpim1//".fits", ver-, >& "dev$null")
                        tmpim1 = mktemp("tmpim1")
                        next_sci= next_sci+1
                    } 
                    n_j = n_j+1
                }
            } #end if-else-if on modeflag
	    
            # Delete tmp images
            imdelete (tmpim1, ver-, >& "dev$null")
            if ((modeflag == 2) || (modeflag == 1)) {
                imdelete (tmpim2, ver-, >& "dev$null")
            }
            if (modeflag == 3) {
                imdelete (tmpim2, ver-, >& "dev$null")
                imdelete (tmpim3, ver-, >& "dev$null")
                imdelete (tmpim4, ver-, >& "dev$null")
            }
        }

        #-------------
        # Pack up the results and clean up	

        # Time stamps
        gemdate ()
        gemhedit (out[i]//"[0]", "GEM-TLM",gemdate.outdate,
            "UT Last modification with GEMINI IRAF")
        gemhedit (out[i]//"[0]", "MPREPARE", gemdate.outdate,
            "UT Time stamp for MPREPARE")
        gemhedit (out[i]//"[0]", "PREPARE", gemdate.outdate, 
            "UT Time stamp - Data have been PREPAREd")
        gemhedit (out[i]//"[0]", "DISPAXIS", 1, "Dispersion direction")

        if (polarimetry)
            next_tot = 8 * next_inp
        else
            next_tot = next_inp
            
        for (n_j=1; n_j<=next_tot; n_j=n_j+1) {
            gemhedit (out[i]//"["//n_j//"]", "AXISLAB4", "", "", delete=yes)
            gemhedit (out[i]//"["//n_j//"]", "EXTVER", n_j, "Extension version")
            gemhedit (out[i]//"["//n_j//"]", "EXTNAME", "SCI", "Extension name")
        }
        
        glogprint (l_logfile, "mprepare", "visual", type="visual",
            vistype="empty", verbose=l_verbose )
        glogprint (l_logfile, "mprepare", "engineering", type="string",
            str="Finished with image number "//i//" at ["//gemdate.outdate//"]",
            verbose=l_verbose )
        glogprint (l_logfile, "mprepare", "visual", type="visual",
            vistype="shortdash", verbose=l_verbose)

nextimage:
        i = i+1

    }

clean:
    if ((status == 0) && (nerror != 0))
        status = nerror
        
    scanfile = "" 
    delete (tmpfile, ver-, >& "dev$null")
    delete (tmpinimg//","//tmpoutimg, ver-, >& "dev$null")
    
    if (status == 0)
        glogclose (l_logfile, "mprepare", fl_success+, verbose=l_verbose)
    else
        glogclose (l_logfile, "mprepare", fl_success-, verbose=l_verbose)
	

# do not put an statement after exitnow.  KL
exitnow:
    ;

end

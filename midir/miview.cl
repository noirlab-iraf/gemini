# Copyright(c) 2001-2011 Association of Universities for Research in Astronomy, Inc.

procedure miview(inimages)

#View frames of raw Michelle data: src, ref or diff (default=dif);
# interactive mode: step through frames, eventually will be able to
# select frames for removal from stack
# non-interactive: just display frames with specified delay
#
# Images are displayed in ximtool; nothing is written to logfile 
#
#NOTES: in current form, need to have ximtool already running;
#       works on single image only
#
# Version 06 Feb 2003 - BR
#
# History:
# 15 Nov 2001,BR- modeled after IDL UF f6movie task; 
#                 started work 22 Oct01 
# 25 Jul 2002, BR - copied from oview, modified to handle MEF
# 05 Feb 2003, BR - major rewrite for Michelle data (will work for T-ReCS too, after 'tprepare' is written)
# 21 Jun 2003, TB - edit to work with real michelle data - works.
# 05 Oct 2003  TB - edit to work with michelle data only.
# 29 Oct 2003  KL - IRAF 2.12 - new parameters
#                     imstat: nclip,lsigma,usigma,cache (no changes required)
#                     hedit: addonly
#
# 10 Dec 2003  KV - changed 'delete("tmp*")' to specific delete of tmp files
#                   used in the script

char    inimages        {prompt="Images to display"}
char    outimages       {"",prompt="Output images(s)"}               
char    outprefix       {"v",prompt="Prefix for out images(s)"} 
char    rawpath         {"",prompt="Path for input raw images"}
char    type            {"dif",prompt="src|ref|dif|sig"}
real    delay           {0.,prompt="update delay in seconds"}
bool    fl_inter        {no,prompt="Run interactively?"}
bool    fl_disp_fill    {no,prompt="Fill display?"}
bool    fl_verbose      {yes,prompt="Verbose?"}
char    logfile         {"",prompt="Logfile"}                        
real    z1              {0.,prompt="Minimum level to be displayed"}
real    z2              {0., prompt="Maximum level to be displayed"}
bool    zscale          {yes,prompt="Auto set grayscale display range"}
bool    zrange          {yes,prompt="Auto set image intensity range"}
char    ztrans          {"linear",prompt="Greyscale transformation"}
int     status          {0,prompt="Exit status (0=good)"}
struct* scanfile        {"",prompt="Internal use only"} 

begin

    char    l_inimages, l_type, l_logfile, l_outimages, l_rawpath
    char    l_prefix, l_temp, lztrans, check
    char    l_filename, in[100], out[100], tmpfile
    char    paramstr, cursinput, header, instrument
    real    l_delay, lz1, lz2

    bool    l_verbose, l_dispfill, l_manual, lzscale, lzrange, testmode
    int     n_nodsets, i, nextns, naxis3, n_a, n_b, nbad, nimages, noutimages
    int     n_i, n_j, ext, ndim, wcs, nbadsets, npars, maximages, modeflag
    real    xpos, ypos

    testmode = no

    tmpfile = mktemp("tmpfile")

    l_inimages=inimages ; l_type=type ; l_delay=delay ;
    l_verbose=fl_verbose ; l_dispfill=fl_disp_fill ; 
    l_manual = fl_inter
    cursinput = ""
    l_logfile = logfile
    l_outimages = outimages
    l_rawpath = rawpath
    l_prefix = outprefix
    lz1=z1 ; lz2=z2
    lztrans=ztrans ; lzscale=zscale ; lzrange=zrange

    if (testmode) time

    cache ("gemdate")

    # Initialize
    nimages = 0
    status = 0
    maximages = 100

    # Create the list of parameter/value pairs.  One pair per line.
    # All lines combined into one string.  Line delimiter is '\n'.
    paramstr =  "inimages       = "//inimages.p_value//"\n"
    paramstr += "outimages      = "//outimages.p_value//"\n"
    paramstr += "outprefix      = "//outprefix.p_value//"\n"
    paramstr += "rawpath        = "//rawpath.p_value//"\n"
    paramstr += "type           = "//type.p_value//"\n"
    paramstr += "delay          = "//delay.p_value//"\n"
    paramstr += "fl_inter       = "//fl_inter.p_value//"\n"
    paramstr += "fl_disp_fill   = "//fl_disp_fill.p_value//"\n"
    paramstr += "fl_verbose     = "//fl_verbose.p_value//"\n"
    paramstr += "logfile        = "//logfile.p_value//"\n"
    paramstr += "z1             = "//z1.p_value//"\n"
    paramstr += "z2             = "//z2.p_value//"\n"
    paramstr += "zscale         = "//zscale.p_value//"\n"
    paramstr += "ztrans         = "//ztrans.p_value

    # Assign a logfile name if not given.  Open logfile and start log.
    # Write parameter/value pairs ("paramstr") to log.
    gloginit (l_logfile, "miview", "midir", paramstr, fl_append+,
        verbose=l_verbose)
    if (gloginit.status != 0) {
        status = gloginit.status
        goto exitnow
    }
    l_logfile = gloginit.logfile.p_value

    # Load up arrays of input name lists
    # This version handles both *s, / and commas in l_inimages

    # Check the rawpath name for a final /
    if (substr(l_rawpath,(strlen(l_rawpath)),(strlen(l_rawpath))) != "/")
        l_rawpath=l_rawpath//"/"
    if (l_rawpath=="/" || l_rawpath==" ")
        l_rawpath=""

    # check that list file exists
    if (substr(l_inimages,1,1)=="@") {
        l_temp=substr(l_inimages,2,strlen(l_inimages))
        if (!access(l_temp) && !access(l_rawpath//l_temp)) {
            glogprint(l_logfile, "miview", "status", type="error", errno=101,
                str="Input file "//l_temp//" not found.",verbose+)
            status=1
            goto clean
        }
    }

    # Count the number of in images
    # First, generate the file list if needed
    if (stridx("*",l_inimages) > 0) {
        files(l_inimages, > tmpfile)
        l_inimages="@"//tmpfile
    }

    if (substr(l_inimages,1,1)=="@")
        scanfile=substr(l_inimages,2,strlen(l_inimages))
    else {
        files(l_inimages,sort-, > tmpfile)
        scanfile=tmpfile
    }

    i = 0
    nimages = 0
    nbad = 0
    noutimages = 0
    npars = 0
    while (fscan(scanfile,l_filename) != EOF && i <= 100) {
        i=i+1

        if (substr(l_filename,strlen(l_filename)-4,strlen(l_filename)) == ".fits")
            l_filename=substr(l_filename,1,strlen(l_filename)-5)

        l_temp = substr(l_inimages,2,strlen(l_inimages))
        if (!access(l_filename//".fits") && !access(l_rawpath//l_filename//".fits")) {
            glogprint(l_logfile, "miview", "status", type="error", errno=101,
                str="Input image "//l_filename//" was not found.",verbose+)
            status=1
            goto clean
        } else {
            nimages=nimages+1
            if (nimages > maximages) {
                glogprint(l_logfile, "miview", "status", type="error",
                    errno=121, str="Maximum number of input images \
                    exceeded:"//maximages,verbose+)
                status=1 
                goto clean
            }
        }
        if (l_rawpath=="" || l_rawpath==" ")
            in[nimages]=l_filename
        else
            in[nimages]=l_rawpath//l_filename
    }

    scanfile = ""
    delete (tmpfile,verify-,>& "dev$null")

    if (nimages == 0) {
        glogprint(l_logfile, "miview", "status", type="error", errno=121,
            str="No input images defined.",verbose+)
        status=1
        goto clean
    }

    # If prefix is to be used instead of filename
    if (l_manual) {
        if (l_outimages=="" || l_outimages==" ") {
            print(l_prefix) | scan(l_prefix)
            if (l_prefix=="" || l_prefix==" ") {
                glogprint(l_logfile, "miview", "status", type="error",
                    errno=121, str="Neither output images name nor \
                    output prefix is defined.",verbose+)
                status=1
                goto clean
            }

            if (l_outimages=="" || l_outimages==" ")
                noutimages = nimages

            i = 1
            while (i <= nimages) {
                # Append prefix to image name
                fparse(in[i])
                out[i] = l_prefix//fparse.root//".fits"

                if (imaccess(out[i])) {
                    glogprint(l_logfile, "miview", "status", type="error",
                        errno=102, str="Output image "//out[i]//" already \
                        exists.",verbose+)
                    nbad+=1
                }
                i+=1
            }
        } else {
            # Now, do the same counting for the out file
            tmpfile=mktemp("tmpfile")

            if (substr(l_outimages,1,1) == "@")
                scanfile=substr(l_outimages,2,strlen(l_outimages))
            else {
                if (stridx("*",l_outimages) > 0) {
                    files(l_outimages,sort-) | match(".hhd",stop+,print-,
                        metach-, > tmpfile)
                    scanfile=tmpfile
                } else {
                    files(l_outimages,sort-, > tmpfile)
                    scanfile=tmpfile
                }
            }

            while (fscan(scanfile,l_filename) != EOF) {
                noutimages=noutimages+1
                if (noutimages > maximages) {
                    glogprint(l_logfile, "miview", "status", type="error",
                        errno=121, str="Maximum number of output images \
                        exceeded:"//maximages,verbose+)
                    status=1
                    goto clean
                }
                if (substr(l_filename,strlen(l_filename)-4,strlen(l_filename)) != ".fits")
                    out[noutimages]=l_filename//".fits"
                else
                    out[noutimages]=l_filename

                if (imaccess(out[noutimages])) {
                    glogprint(l_logfile, "miview", "status", type="error",
                        errno=102, str="Output image "//l_filename//" already \
                        exists.",verbose=l_verbose)
                    nbad+=1
                }
            }

            if (noutimages != nimages) {
                glogprint(l_logfile, "miview", "status", type="error",
                    errno=121, str="Different number of in images \
                    ("//nimages//") and out images ("//noutimages//")",
                    verbose+)
                status =1
                goto clean
            }

            scanfile=""
            delete(tmpfile,verify-, >& "dev$null")
        }
        if (nbad > 0) {
            glogprint(l_logfile, "miview", "status", type="error",
                errno=102, str=nbad//" image(s) already exist.",verbose+)
            status=1
            goto clean
        }
    }

    i=1
    while(i<=nimages) {
        #check instrument
        header=in[i]//"[0]"

        imgets(header,"INSTRUMENT")
        instrument = imgets.value
        glogprint(l_logfile, "miview", "science", type="string",
            str="Instrument is:"//instrument,verbose=l_verbose)

        imgets(header,"MISTACK", >& "dev$null")
        check = imgets.value

        if (check == "0") {
            if (instrument == "michelle") {
                imgets(header,"MPREPARE", >& "dev$null")
                if (imgets.value == "0") {
                    glogprint(l_logfile, "miview", "status", type="error",
                        errno=123, str="Image "//in[i]//" not MPREPAREd.  \
                        Use MVIEW to display.",verbose+)
                    status=1
                    goto clean
                } else {
                    glogprint (l_logfile, "miview", "status", type="string",
                        str="Displaying michelle images.",verbose=l_verbose)
                    imgets(header,"MODE", >& "dev$null")
                    if (imgets.value == "0") {
                        glogprint(l_logfile, "miview", "status", type="error",
                            errno=131, str="Could not find the MODE from \
                            the primary header.",verbose+)
                        status=status+1
                        goto clean
                    }
                    modeflag=0
                    if (imgets.value == "chop-nod") modeflag=1
                    if (imgets.value == "ndchop") modeflag=1
                    if (imgets.value == "chop") modeflag=2
                    if (imgets.value == "nod") modeflag=3
                    if (imgets.value == "ndstare") modeflag=4
                    if (imgets.value == "stare") modeflag=4
                    if (modeflag == 0) {
                        glogprint(l_logfile, "miview", "status", type="error",
                            errno=132, str="Unrecognised MODE \
                            ("//imgets.value//") in the primary header.",
                            verbose+)
                        status=status+1
                        goto clean
                    }
                }
            }

            if (instrument == "TReCS") {
                imgets(header,"TPREPARE", >& "dev$null")
                if (imgets.value == "0") {
                    glogprint(l_logfile, "miview", "status", type="error",
                        errno=123, str="Image "//in[i]//" not TPREPAREd.  \
                        Use TVIEW to display.",verbose+)
                    status=1
                    goto clean
                } else {
                    glogprint(l_logfile, "miview", "status", type="string",
                        str="Displaying TReCS images.", verbose=l_verbose)
                    imgets(header,"OBSMODE", >& "dev$null")
                    if (imgets.value == "0") {
                        glogprint(l_logfile, "miview", "status", type="error",
                            errno=131, str="Could not find the OBSMODE from \
                            the primary header.",verbose+)
                        status=status+1
                        goto nextimage
                    }
                    modeflag=0
                    if (imgets.value == "chop-nod") modeflag=1
                    if (imgets.value == "chop") modeflag=2
                    if (imgets.value == "nod") modeflag=3
                    if (imgets.value == "stare") modeflag=4
                    if (modeflag == 0) {
                        glogprint(l_logfile, "miview", "status", type="error",
                            errno=132, str="Unrecognised MODE \
                            ("//imgets.value//") in the primary header.",
                            verbose+)
                        status=status+1
                        goto clean
                    }
                }
            }


            ##############

            if (modeflag == 1 || modeflag==2) {
                #check type
                if ((l_type!="sig") && (l_type!="ref") && (l_type!="dif") && \
                    (l_type!="src")) {
                    glogprint(l_logfile, "miview", "status", type="error",
                        errno=123, str="Image type keyword invalid.  lpar \
                        miview for valid values.",verbose+)    
                    status=1
                    goto nextimage
                }
            }

            if (modeflag == 1 || modeflag==2) {
                #set filename for display type
                if (l_type == "src")
                    ndim=1
                else if (l_type == "ref")
                    ndim=2
                else if (l_type == "dif")
                    ndim=3
                else if (l_type == "sig") {
                    glogprint(l_logfile, "miview", "status", type="warning",
                        errno=99, str="'sig' type not implemented, using \
                        'dif'", verbose+)
                    ndim=3
                } 
            }

            if ((l_type!="sig") && (l_type!="ref") && (l_type!="dif") && \
                (l_type!="src")) {
                glogprint(l_logfile, "miview", "status", type="error",
                    errno=123, str="Images type invalid.  lpar miview for \
                    valid values.",verbose+)    
                status=1
                goto nextimage
            }

            if (modeflag ==3 || modeflag == 4)
                ndim=1

            cache("imgets","display")

            #get axes
            imgets(in[i]//"[1]","i_naxis3") ; naxis3=int(imgets.value)
            if (naxis3!=3) {
                glogprint(l_logfile, "miview", "status", type="error",
                    errno=123, str="Images "//in[i]//" is not the correct \
                    format.",verbose+)
                glogprint(l_logfile, "miview", "engineering", type="string",
                    str="                     n_choppos= "//naxis3, 
                    verbose=l_verbose)
                status=1
                goto nextimage
            }

            #check extensions
            #check number of extensions
            if (instrument == "michelle") {
                imgets(header,"NUMEXT", >& "dev$null")
                nextns=int(imgets.value)
            }
            if (instrument == "TReCS") {
                imgets(header,"NNODS", >& "dev$null")
                nextns=int(imgets.value)
                imgets(header,"NNODSETS", >& "dev$null")
                j=int(imgets.value)
                nextns=nextns*j
            }

            if (nextns == 0) {
                glogprint(l_logfile, "miview", "status", type="error",
                    errno=123, str="Images "//in[i]//" has zero extensions.",
                    verbose+)
                status=1
                goto nextimage
            }

            nbadsets=0

            #display images in non-interactive mode:
            n_i=1
            if (!l_manual) {
                for (n_i=1;n_i<=nextns;n_i+=1) {
                    glogprint(l_logfile, "miview", "engineering", 
                        type="string", str="Displaying "//l_type//": Nod "\
                        //str(n_i), verbose=l_verbose)
                    display (in[i]//"["//str(n_i)//"][*,*,"//str(ndim)//",1]",
                        1, erase-, zscale=lzscale, z1=lz1, z2=lz2,
                        ztrans=lztrans, zrange=lzrange, fill=l_dispfill, 
                        >& "dev$null")
                }
            }

            #now, make the output file that can have the header edited:
            n_i = 1
            if (l_manual) {
                glogprint(l_logfile, "miview", "task", type="string",
                    str="Copying image "//i//", "//in[i]//",", 
                    verbose=l_verbose)
                glogprint(l_logfile, "miview", "task", type="string",
                    str=" to a new file, "//out[i]//".", verbose=l_verbose)
                fxcopy(in[i]//".fits[0]",out[i])
                glogprint(l_logfile, "miview", "engineering", type="string",
                    str=l_type//": Nod "//str(n_i), verbose=l_verbose)
                for (n_i=1;n_i<=nextns;n_i+=1) {
                    fxinsert (in[i]//".fits["//n_i//"]", 
                        out[i]//"["//n_i-1//"]",  groups="", verbose-)
                }
            }

            #Starting the header editing for interactive mode:
            if (l_manual) {
                for (n_i=1;n_i<=nextns;n_i+=1) {
                    if (npars != 2)
                        npars=0
                    glogprint (l_logfile, "miview", "engineering",
                        type="string", str="Displaying "//l_type//": Nod "\
                        //str(n_i), verbose=l_verbose)
                    display (out[i]//"["//str(n_i)//"][*,*,"//str(ndim)//"]",1,
                        erase-,zscale=lzscale,z1=lz1,z2=lz2,ztrans=lztrans,
                        zrange=lzrange,fill=l_dispfill, >& "dev$null")

                    sleep(l_delay)

                    while (npars == 0) {
                        npars=1
                        glogprint(l_logfile, "miview", "status", type="string",
                            str="Starting interactive cursor input \
                            mode for nod:"//n_i//".", verbose+)
                        printf("MIVIEW Press h for help.\n")
                        if (fscan(imcur,xpos,ypos,wcs,cursinput) != EOF) {

                            if (cursinput == "q" || cursinput == "Q") {
                                glogprint(l_logfile, "miview", "status",
                                    type="string", str="Exiting interactive \
                                    mode for nod: "//n_i//".", 
                                    verbose=l_verbose)
                                npars=1
                            }

                            if (cursinput == "x" || cursinput == "X") {
                                glogprint(l_logfile, "miview", "status", 
                                    type="string", str="Exiting interactive \
                                    mode.", verbose=l_verbose)
                                npars=2
                            }

                            #
                            # For the "h" command, loop to get another input 
                            # value.  Any value aside from the defined ones 
                            # causes the next images to be displayed.
                            #
                            if (cursinput == "h" || cursinput == "H") {
                                printf("\n")
                                printf("--------MIVIEW: INTERACTIVE HELP------\
                                    -- \n ")
                                printf("Available key commands:\n (h) print \
                                    this help \n "//
                                    "(b) mark as a bad frame \n "//
                                    "(u) unmark as a bad frame \n "// 
                                    "(i) run imexamine on this image.  \n "//
                                    "(q) stop interactive mode and move on to \
                                    next nod \n "//
                                    "(s) get images statistics  \n "// 
                                    "(x) exit interactive mode immediately.\n")
                                printf("Key commands can be entered in upper \
                                    or lower case.  Any undefined keystroke \
                                    \n will automatically advance the display \
                                    to the next nod image.\n")
                                npars=0
                                printf("--------------------------------------\
                                    \n ")
                                printf("\n")
                            }
                            if (cursinput == "s" || cursinput == "S") {
                                imstat(out[i]//"["//str(n_i)//"][*,*,"\
                                    //str(ndim)//",1,1]")
                                npars=0
                            }
                            #
                            # For the moment, disable marking bad frames....
                            #

                            if (cursinput == "b" || cursinput == "B") {
                                glogprint(l_logfile, "miview", "science",
                                    type="string", str="Nod "//n_i//" marked \
                                    as bad.",verbose=l_verbose)
                                gemhedit (out[i]//"["//str(n_i)//"]", "BADNOD",
                                    "1", "", delete-)
                                nbadsets=nbadsets+1
                                npars=0
                            }
                            if (cursinput == "u" || cursinput == "U") {
                                glogprint(l_logfile, "miview", "science",
                                    type="string", str="Nod "//n_i//" \
                                    unmarked as bad.",verbose=l_verbose)
                                gemhedit (out[i]//"["//str(n_i)//"]", "BADNOD",
                                    "", "", delete=yes)
                                nbadsets=nbadsets-1
                                npars=0
                            }
                            if (cursinput == "i" || cursinput == "I") {
                                if (l_verbose)
                                    printf("MIVIEW Entering imexam.\n")
                                imexamine()
                                if (l_verbose)
                                    printf("MIVIEW Exiting imexam.\n")
                                npars=0
                            }
                            if (npars == 1) {
                                if (l_verbose) {
                                    printf("MIVIEW Going to next frame.\n")
                                }
                            }
                        }
                    } #end while-loop 
                } #end for-loop

                sleep(l_delay) 

            } #end if (l_manual)
        } #end if (check == "0")

        if (check != "0") {
            imgets(header,"NUMEXT", >& "dev$null")
            nextns=int(imgets.value)

            glogprint(l_logfile, "miview", "task", type="string",
                str="Output image: "//out[i], verbose=l_verbose)
            glogprint(l_logfile, "miview", "status", type="warning", errno=0,
                str="Data has been already been stacked with MISTACK, \
                displaying coadded frame.", verbose=l_verbose)

            #display images in non-interactive mode:
            n_i=1
            if (!l_manual) {
                glogprint(l_logfile, "miview", "engineering", type="string",
                    str="Displaying "//in[i], verbose=l_verbose)
                display (in[i]//"["//str(n_i)//"][*,*]", 1, erase-, 
                    zscale=lzscale, z1=lz1, z2=lz2, ztrans=lztrans,
                    zrange=lzrange, fill=l_dispfill, >& "dev$null")
            }

            if (l_manual) {
                nbadsets=0
                glogprint(l_logfile, "miview", "task", type="string",
                    str="Copying image "//i//", "//in[i]//",", 
                    verbose=l_verbose)
                glogprint(l_logfile, "miview", "task", type="string",
                    str="to a new file, "//out[i]//".", verbose=l_verbose)
                fxcopy(in[i]//".fits[0]",out[i])
                tmpfile1=mktemp("tmpfile1")
                imcopy(in[i]//"["//n_i//"]",tmpfile1,verbose=no)
                fxinsert(tmpfile1//".fits",out[i]//"["//n_i-1//"]",
                    groups="",verbose-)
                imdelete(tmpfile1,verify-,>& "dev$null")
            }

            #Starting the header editing for interactive mode:
            if (l_manual) {
                if (npars != 2)
                    npars=0
                glogprint (l_logfile, "miview", "engineering", type="string",
                    str="Displaying "//out[i], verbose=l_verbose)
                display (out[i]//"["//str(n_i)//"]", 1, erase-, zscale=lzscale,
                    z1=lz1, z2=lz2, ztrans=lztrans, zrange=lzrange,
                    fill=l_dispfill, >& "dev$null")

                sleep(l_delay)

                while (npars == 0) {
                    npars=1
                    glogprint(l_logfile, "miview", "status", type="string",
                        str="Starting interactive cursor input mode for \
                        nod:"//n_i//".", verbose+)
                    printf("     MIVIEW Press h for help.\n")
                    if (fscan(imcur,xpos,ypos,wcs,cursinput) != EOF) {
                        if (cursinput == "q" || cursinput == "Q") {
                            glogprint(l_logfile, "miview", "status",
                                type="string", str="Exiting interactive mode \
                                for  nod: "//n_i//".", verbose=l_verbose)
                            npars=1
                        }
                        if (cursinput == "x" || cursinput == "X") {
                            glogprint(l_logfile, "miview", "status",
                                type="string", str="Exiting interactive mode.",
                                verbose=l_verbose)
                            npars=2
                        }

                        #
                        # For the "h" command, loop to get another input value.
                        # Any value aside from the defined ones causes the 
                        # next images to be displayed.
                        #
                        if (cursinput == "h" || cursinput == "H") {
                            printf("\n")
                            printf("--------MIVIEW: INTERACTIVE HELP--------\
                                \n")
                            printf("Available key commands:\n (h) print this \
                                help \n "//
                                "(b) mark as a bad frame \n "//
                                "(u) unmark as a bad frame \n "// 
                                "(i) run imexamine on this image.  \n "//
                                "(q) stop interactive mode and move on to \
                                next nod \n "//
                                "(s) get images statistics  \n "// 
                                "(x) exit interactive mode immediately for \
                                all images.  \n")
                            printf("Key commands can be entered in upper or \
                                lower case.  Any undefined keystroke \n will \
                                automatically advance the display to the next \
                                nod image.\n")
                            npars=0
                            printf("--------------------------------------\n")
                            printf("\n")
                        }
                        if (cursinput == "s" || cursinput == "S") {
                            imstat(out[i]//"["//str(n_i)//"][*,*,"\
                                //str(ndim)//",1,1]")
                            npars=0
                        }
                        #
                        # For the moment, disable marking bad frames....
                        #

                        if (cursinput == "b" || cursinput == "B") {
                            glogprint(l_logfile, "miview", "science",
                                type="string", str="Nod "//n_i//" marked as \
                                bad.",verbose=l_verbose)
                            gemhedit (out[i]//"["//str(n_i)//"]", "BADNOD", 
                                "1", "", delete-)
                            nbadsets=nbadsets+1
                            npars=0
                        }
                        if (cursinput == "u" || cursinput == "U") {
                            glogprint(l_logfile, "miview", "science",
                                type="string", str="Nod "//n_i//" unmarked as \
                                bad.", verbose=l_verbose)
                            gemhedit (out[i]//"["//str(n_i)//"]", "BADNOD", "",
                                "", delete=yes)
                            nbadsets=nbadsets-1
                            npars=0
                        }
                        if (cursinput == "i" || cursinput == "I") {
                            if (l_verbose)
                                printf("MIVIEW Entering imexam.\n")
                            imexamine()
                            if (l_verbose)
                                printf("MIVIEW Exiting imexam.\n")
                            npars=0
                        }

                        if (npars == 1) {
                            if (l_verbose)
                                printf("\n\n")
                        }
                    }
                } #end of while-loop
            } #end of if (l_manual)

            sleep(l_delay) 

        } #end of if (check != "0")

        if (!l_manual) 
            delete(out[i],verify-, >& "dev$null")

        if (l_manual) {
            if (nbadsets == 0) {
                delete(out[i],verify-, >& "dev$null")
                glogprint(l_logfile, "miview", "task", type="string",
                    str="No bad nodsets identified.  Header not changed, ",
                    verbose=l_verbose)
                glogprint(l_logfile, "miview", "task", type="string",
                    str="so no output image has been written to disk.",
                    verbose=l_verbose)
            } else {
                glogprint(l_logfile, "miview", "science", type="string",
                    str="The number of nodsets marked as bad is:"//\
                    nbadsets//".", verbose=l_verbose)
                gemdate ()
                gemhedit (out[i]//"[0]", "GEM-TLM", gemdate.outdate,
                    "UT Last modification with GEMINI IRAF", delete-)
                gemhedit (out[i]//"[0]", "MIVIEW", gemdate.outdate,
                    "UT Time stamp for MIVIEW", delete-)
            }
        }
        glogprint(l_logfile, "miview", "status", type="string",
            str="Now done displaying image:"//i//".", verbose=l_verbose)

        if (check == "0")
            npars=0

        glogprint(l_logfile, "miview", "visual", type="visual", 
            vistype="longdash", verbose=l_verbose)
nextimage:
        i=i+1
    }

    if (testmode) time

    glogprint(l_logfile, "miview", "visual", type="visual", vistype="longdash",
        verbose=l_verbose)

clean:
    #cleanup

    if (status==0)
        glogclose(l_logfile, "miview", fl_success+, verbose=l_verbose)
    if (status!=0)
        glogclose(l_logfile, "miview", fl_success-, verbose=l_verbose)

    scanfile="" 
    delete("tmpfile*",verify-, >& "dev$null")

exitnow:
    ;

end




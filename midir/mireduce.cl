# Copyright(c) 2003-2011 Association of Universities for Research in Astronomy, Inc.

procedure mireduce(inimages)

# Routine "mireduce"
#
# This script carries out the basic mid-ir pipeline on either (1) raw T-RECS or
# Michelle data files or (2) output files from TPREPARE or MPREPARE.
# It is supposed to apply all the steps needed to reduce an initial raw image 
# to a two-dimensional output image.
#
# This task uses most of the other tasks in the "midir" package 
#
# T-ReCS   : instrument=1
# Michelle : instrument=2
# CanariCam: instrument=3
#
# Version  Oct  17, 2003   KV made changes so it works with Michelle data
#          Oct  10, 2003   KV rewrote internal logic to clean up execution
#          Oct  05, 2003   KV cleaned up warning messages
#          Oct  03, 2003   KV syntax changes in parameter names, added bpm
#          Sept 30, 2003   KV added functionality for MIPBM, MIFLAT, MIREGISTER
#          Sept 11, 2003   KV writes first version
#          Jan  28, 2004   TB updated to include src|ref|dif stacking
#          Aug  19, 2004   KV added "combine" option--average or sum
#          Aug  24, 2004   KV added the "fl_check" flag
#          Oct  21, 2005   KV added the "fl_rescue" flag, for mprepare
#          Oct  21, 2005   KV added the "fl_polarimetry" flag, signals polarimetry mode
#          Dec   7, 2005   KV took out "fl_polarimetry" since mprepare now
#                            checks the header to flag polarimetry mode
#          Feb  18, 2006   KV added the fl_variance flag, passed to mistack and miregister; added "region" as well

char    inimages        {prompt="Input T-ReCS or Michelle image(s)"} # OLDP-1-input-primary-single-prefix=r
char    outimages       {"",prompt="Output image(s)"}                   # OLDP-1-output
char    outpref         {"r",prompt="Prefix for output image name(s)"}  # OLDP-4
char    rawpath         {"",prompt="Path for input raw images"}         # OLDP-4
bool    fl_background   {no,prompt="Apply T/M/MIBACKGROUND to the image"}  # OLDP-2
bool    fl_view         {no,prompt="Apply T/M/MIVIEW to the image"}     # OLDP-2
bool    fl_mask         {no,prompt="Apply a pixel mask to the image"}   # OLDP-2
char    bpm             {"",prompt="bad pixel mask file name"}          # OLDP-1-input-fl_mask-required
bool    fl_flat         {no,prompt="Apply flatfield frame to the image?"}   # OLDP-2
char    flatfieldfile   {"",prompt="Flat field file name(s)"}           # OLDP-1-input-fl_flat-required
char    stackoption     {"stack",prompt="Image combining option: stack|register"}	# OLDP-2
char    frametype       {"dif",prompt="Type of frame to stack (src, ref, dif)"}	# OLDP-2
char    combine         {"average",prompt="Combining images by average|sum"}	# OLDP-2
bool    fl_display      {no,prompt="Display each final image"}          # OLDP-2
int     frame           {1, prompt="Frame to display data"}
bool    fl_check        {yes,prompt="Check images during processing"}   # OLDP-2
bool    fl_rescue       {no,prompt="Try to rescue files (Michelle only)"}    # OLDP-2
char    region          {"[*,*]",prompt="Region to be used for registration (miregister)"}    # OLDP-4
bool    fl_variance     {no,prompt="Output variance frame (mistack, miregister)"}    # OLDP-2
char    logfile         {"",prompt="Log file name"}                     # OLDP-1
bool    verbose         {yes,prompt="Verbose logging yes/no?"}          # OLDP-4
int     status          {0,prompt="Exit status: (0=good, >0=bad)"}      # OLDP-4
struct  *scanfile       {"",prompt="Internal use only"}                 # OLDP-4

begin

    char    l_inputimages, l_prefix, l_image, l_rawpath, l_outputimages
    char    l_flatfieldfile, l_stackoption, l_logfile, l_filename, l_mask
    char    l_frametype, l_combine, l_region
    bool    l_verbose, l_doflat, l_doprepare, l_dobackground, l_domask, l_doview
    bool    l_display, l_check, l_rescue, l_fl_variance
	int     l_frame

    char    in[200], out[200], flat[200]
    char    in1, out1
    char    tmpfile, phu, interfile, tempprefix, interfile1, header, extname
    char    tmpstring
    char    paramstr
    int     nimages, modeflag, maximages, noutputimages, nflats
    int     i, j, k, l, m, n, m1, l1
    int     instrument, ngood, nwarnings
    int     naxis[4], iaxis

    l_inputimages = inimages
    l_outputimages = outimages
    l_rawpath = rawpath
    l_flatfieldfile = flatfieldfile
    l_prefix = outpref
    l_stackoption = stackoption
    l_verbose = verbose
    l_doflat = fl_flat
    l_doview = fl_view
    l_dobackground = fl_background
    l_domask = fl_mask
    l_logfile = logfile
    l_display = fl_display
	l_frame   = frame
    l_mask = bpm
    l_frametype = frametype
    l_combine = combine
    l_check = fl_check
    l_rescue=fl_rescue
    l_fl_variance=fl_variance
    l_region=region

    cache ("gemdate")

    # Initialize
    status = 0
    ngood = 0
    nwarnings = 0

    nimages = 0
    nflats = 0
    maximages = 200

    tmpfile = mktemp("tmplist")

    # Create the list of parameter/value pairs.  One pair per line.
    # All lines combined into one string.  Line delimiter is '\n'.
    paramstr =  "inimages       = "//inimages.p_value//"\n"
    paramstr += "outimages      = "//outimages.p_value//"\n"
    paramstr += "outpref        = "//outpref.p_value//"\n"
    paramstr += "rawpath        = "//rawpath.p_value//"\n"
    paramstr += "fl_background  = "//fl_background.p_value//"\n"
    paramstr += "fl_view        = "//fl_view.p_value//"\n"
    paramstr += "fl_mask        = "//fl_mask.p_value//"\n"
    paramstr += "bpm            = "//bpm.p_value//"\n"
    paramstr += "fl_flat        = "//fl_flat.p_value//"\n"
    paramstr += "flatfieldfile  = "//flatfieldfile.p_value//"\n"
    paramstr += "stackoption    = "//stackoption.p_value//"\n"
    paramstr += "frametype      = "//frametype.p_value//"\n"
    paramstr += "combine        = "//combine.p_value//"\n"
    paramstr += "fl_display     = "//fl_display.p_value//"\n"
	paramstr += "fl_frame       = "//frame.p_value//"\n"
    paramstr += "fl_check       = "//fl_check.p_value//"\n"
    paramstr += "fl_rescue      = "//fl_rescue.p_value//"\n"
    paramstr += "region         = "//region.p_value//"\n"
    paramstr += "fl_variance    = "//fl_variance.p_value//"\n"
    paramstr += "logfile        = "//logfile.p_value//"\n"
    paramstr += "verbose        = "//verbose.p_value
    
    # Assign a logfile name if not given.  Open logfile and start log.
    # Write parameter/value pairs ("paramstr") to log.
    gloginit( l_logfile, "mireduce", "midir", paramstr, fl_append+,
        verbose=l_verbose )
    if (gloginit.status != 0) {
        status = gloginit.status
        goto exitnow
    }
    l_logfile = gloginit.logfile.p_value

    # If fl_check is "no", override the fl_background, fl_mask,
    # and fl_flat flags (until we find out how these tasks 
    # interact with "rescued" files...)
    if (!l_check) {
        l_domask=no
        l_dobackground=no
        l_doflat=no
    }

    # Check stackoption and frametype for consistency.  Only "stack" 
    # allows frametypes of "src" or "ref".
    if (l_stackoption == "register" && l_frametype != "dif") {
        l_stackoption="stack"
        glogprint (l_logfile, "mireduce", "status", type="warning", errno=121,
            str="stackoption = register and frametype = "//l_frametype//\
            " are inconsistent.  Resetting stackoption to stack.",
            verbose=l_verbose)
        glogprint (l_logfile, "mireduce", "visual", type="visual",
            vistype="empty", verbose=l_verbose)
    }
    
    # Edit the output prefix for the calibration src| ref frames
    if ((l_frametype == "src") && (l_prefix != ""))
        l_prefix="c"
    if ((l_frametype == "ref") && (l_prefix != "")) 
        l_prefix="a"

    # Check the rawpath name for a final /
    if (substr(l_rawpath,(strlen(l_rawpath)),(strlen(l_rawpath))) != "/")
        l_rawpath = l_rawpath//"/"
    if (l_rawpath=="/" || l_rawpath==" ")
        l_rawpath = ""

    if (l_combine != "average" && l_combine != "sum") {
        glogprint (l_logfile, "mireduce", "status", type="error", errno=121,
            str="Unrecognized combine parameter ("//l_combine//")", verbose+)
        status = 1
        goto exit
    }

    # Count the number of input images 
    # First, generate the file list if needed
    
    if (stridx("*",l_inputimages) > 0) {
        files(l_rawpath//l_inputimages, > tmpfile)
        l_inputimages="@"//tmpfile
    }

    if (substr(l_inputimages,1,1)=="@")
        scanfile=substr(l_inputimages,2,strlen(l_inputimages))
    else {
        if (stridx(",",l_inputimages)==0) 
            files(l_inputimages, > tmpfile)
        else {
            j=9999
            while (j!=0) {
                j = stridx(",",l_inputimages)
                if (j>0)
                    files(substr(l_inputimages,1,j-1), >> tmpfile)
                else
                    files(l_inputimages, >> tmpfile)
                    
                l_inputimages=substr(l_inputimages,j+1,strlen(l_inputimages))
            }
        }
        scanfile=tmpfile
    }

    i=0
    while (fscan(scanfile,l_filename) != EOF) {
        i=i+1

        if (substr(l_filename,strlen(l_filename)-4,strlen(l_filename)) == ".fits")
            l_filename=substr(l_filename,1,strlen(l_filename)-5)

        j=0
        if (stridx("/",l_rawpath) > 0 && stridx("/",l_filename) > 0) {
            j=stridx("/",l_filename)
            if (j > 0) {
                for (k=1; k < 100 && j > 0; k+=1) {
                    l_filename=substr(l_filename,j+1,strlen(l_filename))
                    j=stridx("/",l_filename)
                }
            }
        }

        if (!imaccess(l_rawpath//l_filename)) {
            glogprint( l_logfile, "mireduce", "status", type="error", errno=101,
                str="Input image "//l_rawpath//l_filename//" not found.",
                verbose+)
            status = status + 1
        } else {
            nimages=nimages+1
            if (nimages > maximages) {
                glogprint( l_logfile, "mireduce", "status", type="error",
                    errno=121, str="Maximum number of input images \
                    ["//str(maximages)//"] has been exceeded.", verbose+)
                status=1
                goto exit
            }
            in[nimages]=l_filename
            out[nimages]=l_filename
            j=stridx("/",out[nimages])
            if (j > 0) {
                for (k=1; k < 100 && j > 0; k+=1) {
                    out[nimages]=substr(out[nimages],j+1,strlen(out[nimages]))
                    j=stridx("/",out[nimages])
                }
            }
        }
    }

    scanfile=""
    delete(tmpfile,ver-,>& "dev$null")

    if (nimages == 0) {
        glogprint( l_logfile, "mireduce", "status", type="error", errno=121,
            str="No input images were defined.", verbose+ )
        status=1
        goto exit
    } else {
        glogprint( l_logfile, "mireduce", "status", type="string",
            str="Processing "//str(nimages)//" images(s).", verbose=l_verbose )
    }

    tmpfile=mktemp("tmplist")

    if (l_doflat) {

        # Count the number of flat field images
        # First, generate the file list if needed

        if (stridx("*",l_flatfieldfile) > 0) {
            files(l_flatfieldfile, > tmpfile)
            l_flatfieldfile="@"//tmpfile
        }

        if (substr(l_flatfieldfile,1,1)=="@")
            scanfile=substr(l_flatfieldfile,2,strlen(l_flatfieldfile))
        else {
            if (stridx(",",l_flatfieldfile)==0) 
                files(l_flatfieldfile, > tmpfile)
            else {
                j=9999
                while (j!=0) {
                    j=stridx(",",l_flatfieldfile)
                    if (j>0)
                        files(substr(l_flatfieldfile,1,j-1), >> tmpfile)
                    else
                        files(l_flatfieldfile, >> tmpfile)

                    l_flatfieldfile = substr(l_flatfieldfile,j+1,
                        strlen(l_flatfieldfile))
                }
            }
            scanfile=tmpfile
        }

        i=0
        while (fscan(scanfile,l_filename) != EOF) {
            i=i+1

            if (substr(l_filename,strlen(l_filename)-4,strlen(l_filename)) == ".fits")
                l_filename=substr(l_filename,1,strlen(l_filename)-5)

            j=0
            if (stridx("/",l_rawpath) > 0 && stridx("/",l_filename) > 0) {
                j=stridx("/",l_filename)
                if (j > 0) {
                    for (k=1; k < 100 && j > 0; k+=1) {
                        l_filename=substr(l_filename,j+1,strlen(l_filename))
                        j=stridx("/",l_filename)
                    }
                }
            }

            if (!imaccess(l_filename)) {
                glogprint( l_logfile, "mireduce", "status", type="error",
                    errno=101, str="Flat field image "//l_filename//
                    " not found.", verbose+ )
                status=status+1
            } else {
                nflats=nflats+1
                if (nflats > maximages) {
                    glogprint( l_logfile, "mireduce", "status", type="error",
                        errno=121, str="Maximum number of flat field images \
                        ["//str(maximages)//"] has been exceeded.",
                        verbose+ )
                    status=1
                    goto exit
                }
                flat[nflats]=l_filename
            }
        }

        scanfile=""
        delete(tmpfile,ver-,>& "dev$null")

        if (nflats == 0) {
            glogprint( l_logfile, "mireduce", "status", type="error", errno=121,
                str="No flat field images were defined.", verbose+ )
            status=1
            goto exit
        }

        if (nflats == 1 && nimages > 1) {
            for (j=2; j <= nimages; j+=1) {
                flat[j]=flat[1]
            }
            nflats=nimages
        }

        if (nflats != nimages) {
            glogprint( l_logfile, "mireduce", "status", type="error", errno=121,
                str="The number flat field images ("//str(nflats)//") does not \
                match the number of input images ("//str(nimages)//").",
                verbose+ )
            status=1
            goto exit
        }
    }  #End if (l_doflat)

    # Now, do the same counting for the output file

    tmpfile=mktemp("tmplist")

    noutputimages=0
    if (l_outputimages != "" && l_outputimages != " ") {
        if (substr(l_outputimages,1,1) == "@")
            scanfile=substr(l_outputimages,2,strlen(l_outputimages))
        else {
            if (stridx("*",l_outputimages) > 0) {
                files(l_outputimages,sort-) | match(".hhd",stop+,print-,
                    metach-, > tmpfile)
                scanfile=tmpfile
            } else {
                files(l_outputimages,sort-, > tmpfile)
                scanfile=tmpfile
            }
        }
  
        while (fscan(scanfile,l_filename) != EOF) {
            if (substr(l_filename,strlen(l_filename)-4,strlen(l_filename)) == ".fits")
                l_filename=substr(l_filename,1,strlen(l_filename)-5)

            noutputimages=noutputimages+1
            if (noutputimages > maximages) {
                glogprint( l_logfile, "mireduce", "status", type="error",
                    errno=121, str="Maximum number of output images "//\
                    str(maximages)//" exceeded.",verbose+ )
                status=1
                goto exit
            }
            out[noutputimages]=l_filename
            if (imaccess(out[noutputimages])) {
                glogprint( l_logfile, "mireduce", "status", type="error",
                    errno=102, str="Output image "//l_filename//\
                    " already exists.", verbose+ )
                status=1
                goto exit
            }
        }
        if (noutputimages != nimages) {
            glogprint( l_logfile, "mireduce", "status", type="error", errno=121,
                str="Different number of input ("//str(nimages)//") and output \
                ("//str(noutputimages)//") image names have been specified.",
                verbose+ )
            status=1
            goto exit
        }

        scanfile=""
        delete(tmpfile,ver-, >& "dev$null")
    } else {
        if (l_prefix == "" || l_prefix == " ")
            l_prefix="r"

        for (i=1; i <= nimages; i+=1) {
            out[i]=l_prefix//out[i]
            if (imaccess(out[i])) {
                glogprint( l_logfile, "mireduce", "status", type="error",
                    errno=102, str="Output image "//out[i]//" already exists.",
                    verbose+ )
                status=1
                goto exit
            }
        }
    }

    l=1
    while (l <= nimages) {

        j=0
        instrument=0
        l_image=in[l]
        if (!imaccess(l_rawpath//l_image)) {
            glogprint( l_logfile, "mireduce", "status", type="error", errno=101,
                str="Input image "//l_rawpath//l_image//\
                " could not be accessed.", verbose+ )
            status=status+1
            goto nextimage
        }
        phu=l_rawpath//l_image//"[0]"
        imgets(phu,"INSTRUME",>& "dev$null")
        if ((imgets.value == "TReCS") || (imgets.value == "CanariCam")) {
            if (imgets.value == "TReCS")
                instrument = 1
            else {   # CanariCam
                instrument = 3
                
                ## stackoption = "register" not supported for CanariCam.
                if (l_stackoption == "register") {
                    l_stackoption = "stack"
                    glogprint (l_logfile, "tprepare", "status", type="warning",
                        errno=121, str="stackoption='register' is not \
                        supported for CanariCam. RESETING to STACK.", 
                        verbose=l_verbose)
                    nwarnings += 1
                }
            }
            imgets(phu,"TPREPARE",>& "dev$null")
            if (imgets.value == "0")
                l_doprepare=yes
            else
                l_doprepare=no
                
            # the following added because mibackground DO NOT EXIST!!!
            if ((l_dobackground==yes) && (l_doprepare==no)) {
                l_dobackground=no
                glogprint (l_logfile, "mireduce", "status", type="warning",
                    errno=121, str="Image already TPREPARE'd. Cannot \
                    process background. Resetting 'fl_background' to 'no'",
                    verbose=l_verbose)
            }
        } else if (imgets.value == "michelle") {
            instrument=2
            imgets(phu,"MPREPARE",>& "dev$null")
            if (imgets.value == "0")
                l_doprepare=yes
            else
                l_doprepare=no
            # the following added because mbackground and mibackground
            # DO NOT EXIST !!!
            if (l_dobackground == yes) {
                l_dobackground=no
                glogprint (l_logfile, "mireduce", "status", type="warning",
                    errno=121, str="No 'background' tasks available for \
                    Michelle.  Resetting 'fl_background' to 'no'.",
                    verbose=l_verbose)
            }
        } else {
            glogprint( l_logfile, "mireduce", "status", type="error", errno=123,
                str="Input image "//in[l]//" is not a T-ReCS, Michelle, or \
                CanariCam image file.", verbose+ )
            status=status+1
            goto nextimage
        }

        interfile=in[l]

        if (l_dobackground) {

            interfile1=mktemp("reducetmpfile")
            if (l_doprepare == no) {
                glogprint( l_logfile, "mireduce", "status", type="fork",
                    fork="forward", child="mibackground", verbose=l_verbose )
                mibackground(interfile,outimages=interfile1,rawpath=l_rawpath)
                glogprint( l_logfile, "mireduce", "status", type="fork",
                    fork="backward", child="mibackground", verbose=l_verbose )
                if (mibackground.status == 0)
                    interfile=interfile1
                else {
                    #Once mibackground has been gemlog'ed, remove this glogprint
                    glogprint( l_logfile, "mireduce", "status", type="error",
                        errno=mibackground.status, str="Task MIBACKGROUND \
                        failed for image "//in[l]//".", verbose+ )
                    status=status+1
                    goto nextimage
                }
            } else {
                if ((instrument == 1) || (instrument == 3)) {
                    glogprint( l_logfile, "mireduce", "status", type="fork",
                        fork="forward", child="tbackground", verbose=l_verbose )
                    tbackground(interfile,outimages=interfile1,
                        rawpath=l_rawpath)
                    glogprint( l_logfile, "mireduce", "status", type="fork",
                        fork="backward", child="tbackground",verbose=l_verbose )
                    if (tbackground.status == 0)
                        interfile=interfile1
                    else {
                        #Once tbackground has been gemlog'ed, remove this
                        #glogprint
                        glogprint( l_logfile, "mireduce", "status", 
                            type="error", errno=tbackground.status, 
                            str="Task TBACKGROUND failed for image "//\
                            in[l]//".", verbose+ )
                        status=status+1
                        goto nextimage      
                    }
                } else {
                    glogprint( l_logfile, "mireduce", "status", type="fork",
                        fork="forward", child="mbackground", verbose=l_verbose )
                    mbackground(interfile,outimages=interfile1,
                        rawpath=l_rawpath)
                    glogprint( l_logfile, "mireduce", "status", type="fork",
                        fork="backward", child="mbackground",verbose=l_verbose )
                    if (mbackground.status == 0)
                        interfile=interfile1
                    else {
                        #Once mbackground has been gemlog'ed, remove this
                        #glogprint
                        glogprint( l_logfile, "mireduce", "status",
                            type="error", errno=mbackground.status,
                            str="Task MBACKGROUND failed for image "//\
                            in[l]//".", verbose+ )
                        status=status+1
                        goto nextimage
                    }
                }
            }
        }  # end if (l_dobackground)

        if (l_doview) {

            interfile1=mktemp("reducetmpfile")
            if (l_doprepare == no) {
                glogprint( l_logfile, "mireduce", "status", type="fork",
                    fork="forward", child="miview", verbose=l_verbose )
		    
                if (interfile == l_image)
                    miview(interfile,outimages=interfile1,rawpath=l_rawpath)
                else
                    miview(interfile,outimages=interfile1)
		    
                glogprint( l_logfile, "mireduce", "status", type="fork",
                    fork="backward", child="miview", verbose=l_verbose )
		    
                if (miview.status == 0) {
                    if (imaccess(interfile1))
                        interfile=interfile1
                } else {
                    #Once miview has been gemlog'ed, remove this glogprint
                    glogprint( l_logfile, "mireduce", "status", type="error",
                        errno=miview.status, str="Task MIVIEW failed for \
                        image "//in[l]//".", verbose+ )
                    status=status+1
                    goto nextimage
                }
            } else {
                if ((instrument == 1) || (instrument == 3)) {
                    glogprint( l_logfile, "mireduce", "status", type="fork",
                        fork="forward", child="tview", verbose=l_verbose )

                    if (interfile == l_image)
                        tview(interfile,outimages=interfile1,rawpath=l_rawpath)
                    else
                        tview(interfile,outimages=interfile1)

                    glogprint( l_logfile, "mireduce", "status", type="fork",
                        fork="backward", child="tview", verbose=l_verbose )

                    if (tview.status == 0) {
                        if (imaccess(interfile1))
                            interfile=interfile1
                    } else {
                        #Once tview has been gemlog'ed, remove this glogprint
                        glogprint( l_logfile, "mireduce", "status", 
                            type="error", errno=tview.status, str="Task TVIEW \
                            failed for image "//in[l]//".", verbose+ )
                        status=status+1
                        goto nextimage
                    }
                } else {
                    glogprint( l_logfile, "mireduce", "status", type="fork",
                        fork="forward", child="mview", verbose=l_verbose )

                    if (interfile == l_image)
                        mview(interfile,outimages=interfile1,rawpath=l_rawpath)
                    else
                        mview(interfile,outimages=interfile1)

                    glogprint( l_logfile, "mireduce", "status", type="fork",
                        fork="backward", child="mview", verbose=l_verbose )

                    if (mview.status == 0) {
                        if (imaccess(interfile1))
                            interfile=interfile1
                    } else {
                        #Once mview has been gemlog'ed, remove this glogprint
                        glogprint( l_logfile, "mireduce", "status",
                            type="error", errno=mview.status, str="Task MVIEW \
                            failed for image "//in[l]//".", verbose+ )
                        status=status+1
                        goto nextimage
                    }
                }
            }
        }  #end if (l_doview)

        if (l_stackoption == "stack") {
            if ((instrument == 1) || (instrument == 3)) {
                imgets(phu,"TPREPARE",>& "dev$null")
                if (imgets.value == "0") {
                    tempprefix="t"
                    interfile1="t"//l_image
                    if (interfile1 == out[l] || imaccess(interfile1)) {
                        k=0
                        while (k < 20) {
                            tempprefix="t"//tempprefix
                            interfile1=tempprefix//l_image
                            if (!imaccess(interfile1)) k=50
                            k=k+1
                        }
                        if (k < 40) {
                            glogprint( l_logfile, "mireduce", "status",
                                type="error", errno=102, str="Could not use \
                                the prefix 't' with file "//in[l]//\
                                " because the output file already exists.",
                                verbose+ )
                            status=status+1
                            goto nextimage
                        }
                    }
                    glogprint( l_logfile, "mireduce", "status", type="fork",
                        fork="forward", child="tprepare", verbose=l_verbose )
                    if (interfile == l_image) {
                        tprepare(interfile,rawpath=l_rawpath,
                            outimage=interfile1,stackoption=l_stackoption,
                            outpref="",combine=l_combine,fl_check=l_check,
                            verbose=l_verbose,logfile=l_logfile)
                    } else {
                        tprepare(interfile,rawpath="",outimage=interfile1,
                            stackoption=l_stackoption,outpref="",
                            combine=l_combine, fl_check=l_check,
                            verbose=l_verbose,logfile=l_logfile)
                    }
                    glogprint( l_logfile, "mireduce", "status", type="fork",
                        fork="backward", child="tprepare", verbose=l_verbose )
                    j=tprepare.status
                    if (j != 0) {
                        #Once tprepare has been gemlog'ed, remove this
                        #glogprint
                        glogprint( l_logfile, "mireduce", "status", 
                            type="error", errno=tprepare.status, 
                            str="Task TPREPARE failed for image "//in[l]//".",
                            verbose+ )
                        status=status+1
                        goto nextimage
                    }
                } else {
                    glogprint( l_logfile, "mireduce", "task", type="warning",
                        errno=0, str="File "//l_image//" has already been \
                        'prepared'",verbose=l_verbose )
                    interfile1=interfile
                }
            }
            if (instrument == 2) {
                imgets(phu,"MPREPARE",>& "dev$null")
                if (imgets.value == "0") {
                    tempprefix="m"
                    interfile1="m"//l_image
                    if (interfile1 == out[l] || imaccess(interfile1)) {
                        k=0
                        while (k < 20) {
                            tempprefix="m"//tempprefix
                            interfile1=tempprefix//l_image
                            if (!imaccess(interfile1)) k=50
                            k=k+1
                        }
                        if (k < 40) {
                            glogprint( l_logfile, "mireduce", "status",
                                type="error", errno=102, str="Could not use \
                                the prefix 'm' with file "//in[l]//\
                                " because the output file already exists.",
                                verbose+)
                            status=status+1
                            goto nextimage
                        }
                    }
                    glogprint( l_logfile, "mireduce", "status", type="fork",
                        fork="forward", child="mprepare", verbose=l_verbose )
                    if (interfile == l_image)
                        mprepare(interfile,rawpath=l_rawpath,
                            outimage=interfile1,outpref="",verbose=l_verbose,
                            fl_rescue=l_rescue,logfile=l_logfile)
                    else
                        mprepare(interfile,rawpath="",outimage=interfile1,
                            fl_rescue=l_rescue,
                            outpref="",verbose=l_verbose,logfile=l_logfile)
                    glogprint( l_logfile, "mireduce", "status", type="fork",
                        fork="backward", child="mprepare", verbose=l_verbose )
			
                    j=mprepare.status
                    if (j != 0) {
                        #Once mprepare has been gemlog'ed, remove this
                        #glogprint
                        glogprint( l_logfile, "mireduce", "status", 
                            type="error", errno=mprepare.status, 
                            str="Task MPREPARE failed for image "//in[l]//".",
                            verbose+ )
                        status=status+1
                        goto nextimage
                    }
                } else {
                    glogprint( l_logfile, "mireduce", "task",
                        str="File "//l_image//" has already been 'prepared'.",
                        verbose=l_verbose )
                    interfile1=interfile
                }
            }
            if (imaccess(interfile1))
                interfile=interfile1

            imgets(phu,"MISTACK", >& "dev$null")
            if (imgets.value == "0") {
                interfile1=mktemp("reducetmpfile")

                glogprint( l_logfile, "mireduce", "status", type="fork",
                    fork="forward", child="mistack", verbose=l_verbose )
                if (interfile == l_image) {
                    mistack(interfile,rawpath=l_rawpath,outimages=interfile1,
                        outpref=l_prefix,frametype=l_frametype,combine=l_combine,
                        fl_variance=l_fl_variance,logfile=l_logfile,verbose=l_verbose)
                } else {
                    mistack(interfile,outimages=interfile1,outpref=l_prefix,
                        frametype=l_frametype,logfile=l_logfile,combine=l_combine,
                        fl_variance=l_fl_variance,verbose=l_verbose)
                }
                glogprint( l_logfile, "mireduce", "status", type="fork",
                    fork="backward", child="mistack", verbose=l_verbose )

                if (mistack.status > 0) {
                    #Once mistack has been gemlog'ed, remove this glogprint
                    glogprint( l_logfile, "mireduce", "status", type="error",
                        errno=mistack.status, str="Task MISTACK failed for \
                        image "//in[l]//".", verbose+ )
                    status=status+1
                    goto nextimage
                }
            }
            if (imaccess(interfile1)) {
                interfile=interfile1
            } else {
                glogprint( l_logfile, "mireduce", "task", type="string",
                    str="File "//l_image//" has already been 'stacked'.",
                    verbose=l_verbose )
                interfile1=interfile
            }
            if (imaccess(interfile1))
                interfile=interfile1
        } #end if (l_stackoption == "stack")
        if (l_domask) {
            if (l_mask == "") {
                if (instrument == 1)
                    l_mask="midir$data/trecs.badpixels"
                else if (instrument == 2)
                    l_mask="midir$data/michelle.badpixels"
                else    # instrument == 3
                    l_mask = "midir$data/canaricam.badpixels"
            }
            if (!access(l_mask)) {
                glogprint( l_logfile, "mireduce", "status", type="error",
                    errno=101, str="Bad pixel mask file "//l_mask//\
                    " not found.", verbose=l_verbose )
                status=status+1
                goto nextimage
            }
            glogprint (l_logfile, "mireduce", "task", type="string",
                str="Using bad pixel mask '"//l_mask//"'", verbose=l_verbose)
    
            tmpfile=mktemp("tmplist")
            interfile1=mktemp("reducetmpfile")

            if (interfile == l_image)
                in1=l_rawpath//interfile
            else
                in1=interfile

            fparse(in1)
            in1 = fparse.directory//fparse.root//".fits"
            fparse(interfile1)
            out1 = fparse.directory//fparse.root//".fits"
    
            copy(in1,out1,verbose-)
            glogprint( l_logfile, "mireduce", "engineering", type="string",
                str="Copying image "//interfile//" to "//interfile1,
                verbose=l_verbose )
            header=interfile1//"[0]"
            if (!imaccess(header)) {
                glogprint( l_logfile, "mireduce", "status", type="error",
                    errno=101, str="Output file "//interfile1//" was not \
                    created properly.",verbose+ )
                status=1
                goto nextimage
            }

            l1=0
            for (m=1; m < 100000 && l1 == 0; m+=1) {
                extname=interfile1//"["//str(m)//"]"
                if (imaccess(extname)) {
                    imgets(extname,"i_naxis")
                    iaxis=int(imgets.value)
                    if (iaxis < 2 || iaxis > 4) {
                        glogprint( l_logfile, "mireduce", "status", 
                            type="error", errno=123, str="Input image "//\
                            l_image//" has a bad number of dimensions \
                            ("//str(iaxis)//") for extension "//str(m)//".",
                            verbose+ )
                        status=status+1
                        goto nextimage
                    }
                    naxis[1]=0
                    naxis[2]=0
                    naxis[3]=0
                    naxis[4]=0

                    for (n=1; n <= iaxis; n+=1) {
                        tmpstring="naxis"//str(n)
                        imgets(extname,tmpstring)
                        naxis[n]=int(imgets.value)
                        if (naxis[n] <= 0) {
                            glogprint( l_logfile, "mireduce", "status",
                                type="error", errno=123, str="Image "//\
                                l_image//" has zero or negative image \
                                dimension(s) in extension "//str(m)//".",
                                verbose+ )
                            status=status+1
                            goto nextimage
                        }
                    }

                    if (naxis[1] != 320 || naxis[2] != 240) {
                        glogprint( l_logfile, "mireduce", "status", 
                            type="error", errno=123, str="Image "//l_image//\
                            " does not have the correct image dimensions in \
                            extension "//str(m)//".", verbose+ )
                        status=status+1
                        goto nextimage
                    }
                    if (naxis[3] == 0) naxis[3]=1
                    if (naxis[4] == 0) naxis[4]=1

                    for (m1=1; m1 <= naxis[3]; m1+=1) {
                        for (n=1; n <= naxis[4]; n+=1) {
                            tmpfile=mktemp("tmplist")
                            imedit(extname,tmpfile,cursor=l_mask,display-)
                            if (iaxis == 2)
                                imcopy(tmpfile,extname//"[1:320,1:240]",
                                    verbose-)
                            else if (iaxis == 3)
                                imcopy(tmpfile,extname//"[*,*,"//str(m1)//"]",
                                    verbose-)
                            else if (iaxis == 4)
                                imcopy(tmpfile,extname//"[*,*,"//str(m1)//\
                                    ","//str(n)//"]",verbose-)
                            delete(tmpfile,ver-,>& "dev$null")
                        }
                    }
                } else
                l1=1
            } #end for-loop

            if (imaccess(interfile1))
                interfile=interfile1

        } # end if (l_domask)

        if (l_doflat) {
            if (!imaccess(flat[l])) {
                glogprint( l_logfile, "mireduce", "status", type="error",
                    errno=121, str="The flat field file "//flat[l]//\
                    " was not found.", verbose+ )
                status=status+1
                goto nextimage
            }
    
            tmpfile=mktemp("tmplist")
            interfile1=mktemp("reducetmpfile")

            if (interfile == l_image)
                in1=l_rawpath//interfile
            else
                in1=interfile

            m=stridx(".",in1)
            n=stridx(".",interfile1)
            if (m == 0)
                in1=in1//".fits"
            if (n == 0)
                out1=interfile1//".fits"
            else
                out1=interfile1

            copy(in1,out1,verbose-)
            glogprint( l_logfile, "mireduce", "engineering", str="Copying \
                image "//interfile//" to "//interfile1, verbose=l_verbose )
            header=interfile1//"[0]"
            if (!imaccess(header)) {
                glogprint( l_logfile, "mireduce", "status", type="error",
                    errno=101, str="Output file "//interfile1//" was not \
                    created properly.", verbose+ )
                status=1
                goto nextimage
            }

            l1=0
            for (m=1; m < 100000 && l1 == 0; m+=1) {
                extname=interfile1//"["//str(m)//"]"
                if (imaccess(extname)) {
                    imgets(extname,"i_naxis")
                    iaxis=int(imgets.value)
                    if (iaxis < 2 || iaxis > 4) {
                        glogprint( l_logfile, "mireduce", "status",
                            type="error", errno=123, str="Input image "//\
                            l_image//" has a bad number of dimensions \
                            ("//str(iaxis)//") for extension "//str(m)//".",
                            verbose+ )
                        status=status+1
                        goto nextimage
                    }
                    naxis[1]=0
                    naxis[2]=0
                    naxis[3]=0
                    naxis[4]=0

                    for (n=1; n <= iaxis; n+=1) {
                        tmpstring="naxis"//str(n)
                        imgets(extname,tmpstring)
                        naxis[n]=int(imgets.value)
                        if (naxis[n] <= 0) {
                            glogprint( l_logfile, "mireduce", "status",
                                type="error", errno=123, str="Image "//\
                                l_image//" has zero or negative image \
                                dimension(s) in extension "//str(m)//".",
                                verbose+ )
                            status=status+1
                            goto nextimage
                        }
                    }

                    if (naxis[1] != 320 || naxis[2] != 240) {
                        glogprint( l_logfile, "mireduce", "status",
                            type="error", errno=123, str="Image "//l_image//\
                            "does not have the correct image dimensions in \
                            extension "//str(m)//".", verbose+ )
                        status=status+1
                        goto nextimage
                    }
                    if (naxis[3] == 0) naxis[3]=1
                    if (naxis[4] == 0) naxis[4]=1

                    for (m1=1; m1 <= naxis[3]; m1+=1) {
                        for (n=1; n <= naxis[4]; n+=1) {
                            tmpfile=mktemp("tmplist")
                            imarith(extname,"*",flat[l]//".fits[1]",tmpfile)
                            if (iaxis == 2)
                                imcopy(tmpfile,extname//"[1:320,1:240]",
                                    verbose-)
                            else if (iaxis == 3)
                                imcopy(tmpfile,extname//"[*,*,"//str(m1)//"]",
                                    verbose-)
                            else if (iaxis == 4)
                                imcopy(tmpfile,extname//"[*,*,"//str(m1)//\
                                    ","//str(n)//"]",verbose-)
                                    
                            delete(tmpfile,ver-,>& "dev$null")
                        }
                    }
                } else
                    l1=1
            } #end for-loop

            if (imaccess(interfile1))
                interfile=interfile1

        } #end if (l_doflat)

        if (l_stackoption == "register") {
            if ((instrument == 1) || (instrument == 3)) {
                imgets(phu,"TPREPARE",>& "dev$null")
                if (imgets.value == "0") {
                    tempprefix="t"
                    interfile1="t"//l_image
                    if (interfile1 == out[l] || imaccess(interfile1)) {
                        k=0
                        while (k < 20) {
                            tempprefix="t"//tempprefix
                            interfile1=tempprefix//l_image
                            if (!imaccess(interfile1)) k=50
                            k=k+1
                        }
                        if (k < 40) {
                            glogprint( l_logfile, "mireduce", "status",
                                type="error", errno=102, str="Could not use \
                                the prefix 't' with file "//in[l]//" because \
                                the output file already exists.", verbose+ )
                            status=status+1
                            goto nextimage
                        }
                    }
                    glogprint( l_logfile, "mireduce", "status", type="fork",
                        fork="forward", child="tprepare", verbose=l_verbose )
                    if (interfile == l_image)
                        tprepare(interfile,rawpath=l_rawpath,
                            outimage=interfile1,stackoption=l_stackoption,
                            combine=l_combine,outpref="",fl_check=l_check,
                            verbose=l_verbose,logfile=l_logfile)
                    else
                        tprepare(interfile,rawpath="",outimage=interfile1,
                            stackoption=l_stackoption,outpref="",
                            combine=l_combine, fl_check=l_check,
                            verbose=l_verbose,logfile=l_logfile)
			    
                    glogprint( l_logfile, "mireduce", "status", type="fork",
                        fork="backward", child="tprepare", verbose=l_verbose )

                    j=tprepare.status
                    if (j != 0) {
                        #Once tprepare has been gemlog'ed, remove this
                        #glogprint
                        glogprint( l_logfile, "mireduce", "status", 
                            type="error", errno=tprepare.status, str="Task \
                            TPREPARE failed for image "//in[l]//".",
                            verbose+ )
                        status=status+1
                        goto nextimage
                    }
                } else {
                    glogprint( l_logfile, "mireduce", "task", type="string",
                        str="File "//l_image//" has already been 'prepared'",
                        verbose=l_verbose )
                    interfile1=interfile
                }
                
            } else if (instrument == 2) {
                imgets(phu,"MPREPARE",>& "dev$null")
                if (imgets.value == "0") {
                    tempprefix="m"
                    interfile1="m"//l_image
                    if (interfile1 == out[l] || imaccess(interfile1)) {
                        k=0
                        while (k < 20) {
                            tempprefix="m"//tempprefix
                            interfile1=tempprefix//l_image
                            if (!imaccess(interfile1)) k=50
                            k=k+1
                        }
                        if (k < 40) {
                            glogprint( l_logfile, "mireduce", "status",
                                type="error", errno=102, str="Could not use the \
                                prefix "//"'m'"//" with file "//in[l]//\
                                " because the output file already exists.",
                                verbose+ )
                            status=status+1
                            goto nextimage
                        }
                    }
                    glogprint( l_logfile, "mireduce", "status", type="fork",
                        fork="forward", child="mprepare", verbose=l_verbose )
                    if (interfile == l_image)
                        mprepare(interfile,rawpath=l_rawpath,
                            outimage=interfile1,outpref="",verbose=l_verbose,
                            fl_rescue=l_rescue,logfile=l_logfile)
                    else
                        mprepare(interfile,rawpath="",outimage=interfile1,
                            fl_rescue=l_rescue,
                            outpref="",verbose=l_verbose,logfile=l_logfile)
                    glogprint( l_logfile, "mireduce", "status", type="fork",
                        fork="backward", child="mprepare", verbose=l_verbose )

                    j=mprepare.status
                    if (j != 0) {
                        #Once mprepare has been gemlog'ed, remove this
                        #glogprint
                        glogprint( l_logfile, "mireduce", "status", 
                            type="error", errno=mprepare.status, str="Task \
                            MPREPARE failed for image "//in[l]//".",
                            verbose+ )
                        status=status+1
                        goto nextimage
                    }
                } else {
                    glogprint( l_logfile, "mireduce", "task", type="string",
                        str="File "//l_image//" has already been 'prepared'.",
                        verbose=l_verbose )
                    interfile1=interfile
                }
            }
            if (imaccess(interfile1))
                interfile=interfile1

            imgets(phu,"MIREGIST", >& "dev$null")
            if (imgets.value == "0") {
                interfile1=mktemp("reducetmpfile")
                glogprint( l_logfile, "mireduce", "status", type="fork",
                    fork="forward", child="miregister", verbose=l_verbose )
                if (interfile == l_image)
                    miregister(interfile,rawpath=l_rawpath,outimages=interfile1,
                    combine=l_combine,outpref=l_prefix,logfile=l_logfile,
                    region=l_region,fl_variance=l_fl_variance,verbose=l_verbose)
                else
                    miregister(interfile,outimages=interfile1,outpref=l_prefix,
                        fl_variance=l_fl_variance,combine=l_combine,
                        region=l_region,logfile=l_logfile,verbose=l_verbose)
                        
                glogprint( l_logfile, "mireduce", "status", type="fork",
                    fork="backward", child="miregister", verbose=l_verbose )

                if (miregister.status > 0) {
                    #Once miregister has been gemlog'ed, remove this glogprint
                    glogprint( l_logfile, "mireduce", "status", type="error",
                        errno=miregister.status, str="Task MIREGISTER failed \
                        for image "//in[l]//".", verbose+ )
                    status=status+1
                    goto nextimage
                }
                if (imaccess(interfile1))
                    interfile=interfile1
            } else {
                glogprint( l_logfile, "mireduce", "task", str="File "//\
                    l_image//" has already been 'stacked'.", verbose=l_verbose )
                interfile1=interfile
            }
        } #end if (l_stackoption == "register")

        k=stridx(".",interfile1)
        j=stridx(".",out[l])
        if (k == 0)
            in1=interfile1//".fits"
        else
            in1=interfile1
        if (j == 0)
            out1=out[l]//".fits"
        else
            out1=out[l]

        copy(in1,out1,verbose-)
        glogprint( l_logfile, "mireduce", "engineering", type="string",
            str="Renaming image "//interfile1//" to "//out[l], 
            verbose=l_verbose )

        if (l_display)
            display(out[l]//"[1]", frame=l_frame, bpdisplay="none", \
                bpcolors="red", bpm="BPM", ocolors="green", erase=yes, \
				fill=no, zscale=yes, zrange=yes, xcenter=0.5, ycenter=0.5, \
                xmag=1., ymag=1., z1=0., z2=0)

        # Time stamp the primary header
        #
        gemdate ()
        gemhedit (out[l]//"[0]", "MIREDUCE", gemdate.outdate,
            "UT Time stamp for MIREDUCE")
        gemhedit (out[l]//"[0]", "GEM-TLM", gemdate.outdate,
            "UT Last modification with GEMINI", delete-)
        gemhedit(out[l]//"[1]","EXTVER",1,"Extension version")
        gemhedit(out[l]//"[1]","EXTNAME","SCI","Extension name")

        ngood=ngood+1
        if ((l_frametype == "src") || (l_frametype == "ref"))
            gemhedit (out[l]//"[0]", "NSCUT", gemdate.outdate, 
                "UT Time stamp - arc frame")


nextimage:
        delete("reducetmpfile*.fits",ver-,>& "dev$null")
        delete("tmplist*",ver-,>& "dev$null")
        l=l+1

    } #end while-loop over the images

exit:
    if (nwarnings > 0) {
        glogprint (l_logfile, "tprepare", "visual", type="visual", 
            vistype="empty", verbose=l_verbose)
        glogprint (l_logfile, "tprepare", "status", type="warning", errno=0,
            str="There were "//str(nwarnings)//" warning(s).", verbose=l_verbose)
        glogprint (l_logfile, "tprepare", "status", type="warning", errno=0,
            str="Please review the logs.", verbose=l_verbose)
        glogprint (l_logfile, "tprepare", "visual", type="visual", 
            vistype="empty", verbose=l_verbose)
    }
    if (status == 0) {
        glogprint( l_logfile, "mireduce", "status", type="string",
            str="All "//str(nimages)//" images successfully reduced.",
            verbose=l_verbose )
        glogclose( l_logfile, "mireduce", fl_success+, verbose=l_verbose )
    } else {
        glogprint( l_logfile, "mireduce", "status", type="string",
            str=str(ngood)//" of "//str(nimages)//" successfully processed.",
            verbose+ )
        glogclose( l_logfile, "mireduce", fl_success-, verbose=l_verbose )
    }

exitnow:
    ;

end

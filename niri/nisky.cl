# Copyright(c) 2001-2015 Association of Universities for Research in Astronomy, Inc.

procedure nisky(inimages)

# Make an average or median sky image for NIRI images, removing
# stars first.  Input images (inimages) are the images to combine.
# 
# Version  June 19, 2001  JJ, beta-release
#          Oct 15, 2001  JJ, v1.2 release
#          Nov 16, 2001  JJ, automatic output filename
#          Jan 3, 2002   JJ, fixed bug with fl_vardq+ but DQ doesn't exist
#                            and you want VAR and DQ in the output (as can 
#                            be the case with skyflats)
#          Jan 15, 2002  JJ, added checks for INDEF in object list; adjusted
#                            minimum radius a bit
#          Feb 28, 2002  JJ, v1.3 release
#          Mar 25, 2002  JJ, fixed bug in fixing nonlinear and saturated pix
#                            that kept daofind from finding bright stars
#          May 29, 2002  JJ, removed WMEF calls,modified scilist
#          May 30, 2002  JJ, headers written all at once, debug
#          Jun 07, 2002  JJ, delete output if exiting with ERROR
#                            added logic to avoid unnecessary nifastsky call
#                            fixed bug for repeated input file names
#          Jun 12, 2002  JJ  fixed bug when input file doesn't exist
#          Sept 10, 2002 IJ  parameter encoding
#          Sept 20, 2002 JJ  v1.4 release
#          Aug 18, 2003  KL  IRAF2.12 - new/modified parameters
#                              hedit: addonly
#                              imstat: nclip, lsigma, usigma, cache
#                              imcombine: headers,bpmasks,expmasks,outlimits
#                                         rejmask->rejmasks, plfile->nrejmasks
#                              apphot.daofind: wcsout, cache
#          Sep 15, 2004  JJ  Major revision: changed daofind (etc.) to objmask
#          May 27, 2005  JJ  Write NISKYMSK always, even if fl_keepmasks-
#          Mar 09, 2006  JH  Converted all hedit and gemhedit calls to nhedit
#          Aug 29, 2007  JH  Changed statsec parameter to "default" and added
#                            support for subarrays
#

char    inimages    {prompt="Raw NIRI image list to combine"}                   # OLDP-1-input-primary-combine-suffix=_sky
char    outimage    {"",prompt="Output sky image"}                              # OLDP-1-output
char    outtitle    {"default",prompt="Title for output image"}                 # OLDP-3
char    combtype    {"default",prompt="Type of combine operation",enum="default|median|average"} # OLDP-2
char    rejtype     {"avsigclip",prompt="Type of rejection",enum="none|avsigclip|minmax"} # OLDP-2
char    logfile     {"",prompt="Name of log file"}                              # OLDP-1
char    statsec     {"default",prompt="Statistics section"}           # OLDP-2
int     nlow        {0,min=0,prompt="minmax: Number of low pixels to reject"}   # OLDP-2
int     nhigh       {1,min=0,prompt="minmax: Number of high pixels to reject"}  # OLDP-2
real    lsigma      {3.,min=0.,prompt="avsigclip: Lower sigma clipping factor"} # OLDP-2
real    hsigma      {3.,min=0.,prompt="avsigclip: Upper sigma clipping factor"} # OLDP-2
real    threshold   {3.,min=1.5,prompt="Threshold in sigma for object detection"} # OLDP-2
int     ngrow       {3,prompt="Number of iterations to grow objects into the wings"} # OLDP-2
real    agrow       {3., prompt="Area limit for growing objects into the wings"} # OLDP-2
int     minpix      {6,prompt="Minimum number of pixels to be identified as an object"} # OLDP-3
char    key_exptime {"EXPTIME",prompt="Keyword for exposure time"}            # OLDP-3
char    key_ron     {"RDNOISE",prompt="Keyword for readout noise in e-"}      # OLDP-3
char    key_gain    {"GAIN",prompt="Keyword for gain in electrons/ADU"}       # OLDP-3
char    masksuffix  {"msk",prompt="Mask name suffix"}                         # OLDP-3
bool    fl_keepmasks {no,prompt="Keep object masks for each input image?"}    # OLDP-2
bool    fl_nifastsky {yes, prompt="Reduce inputs before masking objects"}     # OLDP-3
char    sci_ext     {"SCI",prompt="Name or number of science extension"}      # OLDP-3
char    var_ext     {"VAR",prompt="Name or number of variance extension"}     # OLDP-3
char    dq_ext      {"DQ",prompt="Name or number of data quality extension"}  # OLDP-3
bool    fl_vardq    {yes,prompt="Create variance and data quality frames in output?"} # OLDP-3
bool    fl_dqprop   {no,prompt="Retain input DQ information in output?"}      # OLDP-3
bool    verbose     {yes,prompt="Verbose actions"}                            # OLDP-4
int     status      {0,prompt="Exit status (0=good)"}                         # OLDP-4
struct  *scanfile   {prompt="Internal use only"}                              # OLDP-4

begin

    #-------------------------------------------------------------------------
    # declare local variables
    char    l_inimages, l_outimage, l_combtype, l_rejtype, l_logfile
    char    l_key_ron, l_key_gain, l_masksuffix
    char    l_key_exptime, l_outtitle, l_statsec
    real    l_threshold, l_agrow
    int     l_ngrow, l_minpix
    real    l_lsigma, l_hsigma
    bool    l_verbose, l_fl_keepmasks, l_fl_nifastsky
    int     l_nlow, l_nhigh, nfiles
    int     i, n, nbad, nx, ny
    real    skyconst, l_gain
    char    img, imgsci, imgdq
    char    tmpfile1, tmpfile2, tmpimg
    char    tmpsky, tmpskymef, tmpflat, l_temp
    struct  l_struct
    real    l_mean, l_sig
    real    l_expone, l_ronone, l_gainone
    char    l_sci_ext, l_var_ext, l_dq_ext
    char    scilist, dqlist, combsig, combdq
    char    tmphead, tmpdq, tmpdqcomb, dqsumold, dqsum, suf, keyfound
    real    l_skymean
    bool    l_fl_vardq, l_fl_dqprop, needtoreduce

    #------------------------------------------------------------------------
    # Set local variables
    scanfile=""
    status=0
    skyconst=0.
    n=0

    l_inimages=inimages ; l_outimage=outimage
    l_outtitle=outtitle
    l_combtype=combtype ; l_rejtype=rejtype
    l_threshold=threshold ; l_agrow=agrow
    l_ngrow=ngrow ; l_minpix=minpix
    l_key_ron=key_ron ; l_key_gain=key_gain
    l_key_exptime=key_exptime ; l_masksuffix=masksuffix 
    if ( stridx(" ",l_masksuffix)>0 || stridx(",",l_masksuffix) > 0 ) {
        printlog("WARNING - NISKY: Illegal masksuffix "//l_masksuffix//\
            " changed to msk", l_logfile,verbose+)
        l_masksuffix="msk"
    }
    l_logfile=logfile ; l_verbose=verbose
    l_nlow=nlow ; l_nhigh=nhigh ; l_fl_keepmasks=fl_keepmasks
    l_fl_nifastsky = fl_nifastsky
    l_lsigma=lsigma ; l_hsigma=hsigma
    l_sci_ext=sci_ext ; l_dq_ext=dq_ext ; l_var_ext=var_ext
    l_fl_vardq=fl_vardq ; l_fl_dqprop=fl_dqprop 
    l_statsec=statsec

    #------------------------------------------------------------------------
    # Check for package log file or user-defined log file
    cache("niri", "gemdate")

    print(l_logfile) | scan(l_logfile)
    if (l_logfile=="" || l_logfile==" ") {
        l_logfile=niri.logfile
        print(l_logfile) | scan(l_logfile)
        if (l_logfile=="" || l_logfile==" ") {
            l_logfile="niri.log"
            printlog("WARNING - NISKY:  Both nisky.logfile and niri.logfile \
                are empty.", l_logfile, verbose+)
            printlog("                  Using default file niri.log.",
                l_logfile, verbose+)
        }
    }
    # Open log file
    date | scan(l_struct)
    printlog("--------------------------------------------------------------\
        --------------", l_logfile, l_verbose)
    printlog("NISKY -- "//l_struct, l_logfile, l_verbose)
    printlog(" ",l_logfile, l_verbose)

    # check statsec
    if (l_statsec != "default") {
        if (substr(l_statsec,1,1)!="[" || \
            substr(l_statsec,strlen(l_statsec),strlen(l_statsec))!="]" || \
            stridx(",",l_statsec)==0) {
            printlog("WARNING - NISKY:  statsec has wrong format.  Reverting \
                to [*,*]", l_logfile, verbose+) 
            l_statsec="[*,*]"
        }
    }
    # Keep imgets parameters from changing by outside world
    cache("imgets")

    # Make temporary file names
    tmpfile1  = mktemp("tmpfl")
    tmpfile2  = mktemp("tmpfl")
    tmpsky    = mktemp("tmpsky")
    tmpskymef = mktemp("tmpskymef")
    tmpflat   = mktemp("tmpflat")
    scilist   = mktemp("tmpscilist")
    dqlist    = mktemp("tmpdqlist")
    dqsum     = mktemp("tmpdqsum")
    tmpdq     = mktemp("tmpdq")
    combsig   = mktemp("tmpcombsig") 
    combdq    = mktemp("tmpcombdq")
    tmpdqcomb = mktemp("tmpdqcomb")
    tmphead   = mktemp("tmphead")

    #----------------------------------------------------------------------
    # Start checking for valid input

    # Verify that the extension names are not empty, otherwise exit gracefully
    print(l_sci_ext) | scan(l_sci_ext)
    if (l_sci_ext=="" || l_sci_ext==" ") {
        printlog("ERROR - NISKY: extension name sci_ext is missing.",
            l_logfile, verbose+)
        status=1
        goto clean
    }

    if (l_fl_vardq) {
        print(l_dq_ext) | scan(l_dq_ext)
        print(l_var_ext) | scan(l_var_ext)
        if (l_dq_ext=="" || l_dq_ext==" " || l_var_ext=="" || l_var_ext==" ") {
            printlog("WARNING - NISKY: extension name var_ext or dq_ext is \
                missing.", l_logfile, verbose+)
            printlog("                 Output image will not have VAR or DQ \
                planes.", l_logfile, verbose+)
            l_fl_vardq=no
            l_fl_dqprop=no
        } 
    }
    #----------------------------------------------------------------------
    # Put all input images in a temporary file: tmpfile1

    if (substr(l_inimages,1,1)=="@") {
        l_temp=substr(l_inimages,2,strlen(l_inimages))
        if (!access(l_temp)) {
            printlog("ERROR - NISKY: Input file "//l_temp//" not found.",
                l_logfile, verbose+)
            status=1
            goto clean
        }
    }

    files(l_inimages,sort-) | unique("STDIN", > tmpfile1)

    # Verify that input images actually exist. (Redundant for * though)
    # at the same time check the exposure times and MEF format;
    # strip .fits if present
    nbad=0
    nfiles=0
    scanfile = tmpfile1
    while (fscan(scanfile,img) != EOF) {
        # check output image
        if (nfiles==0) {
            print(l_outimage) | scan(l_outimage)
            if (l_outimage=="" || l_outimage==" ") {
                suf = substr(img,strlen(img)-4,strlen(img))
                if (suf==".fits")
                    l_outimage=substr(img,1,strlen(img)-5)//"_sky" 
                else
                    l_outimage=img//"_sky"
            }
            if (imaccess(l_outimage)) {
                printlog("ERROR - NISKY: Output file "//l_outimage//" already \
                    exists.", l_logfile, verbose+)
                status=2
                goto clean
            }
        }

        # check input images
        gimverify(img)
        if (gimverify.status==1) {
            printlog("ERROR - NISKY: File "//img//" not found.",
                l_logfile, verbose+)
            nbad+=1
        } else if (gimverify.status>1) {
            printlog("ERROR - NISKY: File "//img//" not a MEF FITS image.",
                l_logfile, verbose+)
            nbad+=1
        } else {
            keyfound=""
            hselect(img//"[0]","*PREPAR*",yes) | scan(keyfound)
            if (keyfound == "") {
                printlog("ERROR - NISKY: Image "//img//" not *PREPAREd.",
                    l_logfile, verbose+)
                nbad+=1
            }
        }

        # If the files don't exist cannot do the rest of this loop without a
        # nasty error. Skip the rest of the checks for this file and move on
        # to the next one. It will error out after the while loop. - MS
        if (nbad > 0) {
            status = 1
            goto SKIP_CHECKS
        }
        
        # get the dimensions for statistics section
        if (l_statsec == "default") {
            imgets (img//"["//l_sci_ext//","//1//"]","i_naxis1",
                        >& "dev$null")
            nx = int(imgets.value)
            imgets (img//"["//l_sci_ext//","//1//"]","i_naxis2",
                        >& "dev$null")
            ny = int(imgets.value)
            if (nx == 1024 && ny == 1024) { # full image size
                l_statsec = "[100:924,100:924]" # remove 100 pix border
                printlog ("using section "//l_statsec//" for image \
                    statistics", l_logfile, verbose=l_verbose)
            } else { # using subarray, already has border removed 
                l_statsec = "[*,*]"
                printlog("Using entire image for statistics", l_logfile, \
                    verbose=l_verbose)
            }
       }

        # get the header of the first image and prepare for output image
        if (nfiles==0) imcopy(img//"[0]",l_outimage, verbose-)

        if (l_fl_dqprop && !imaccess(img//"["//l_dq_ext//"]") ) {
            printlog("WARNING - NISKY: Cannot propagate input DQ planes \
                because", l_logfile, verbose+)
            printlog("                 input image "//img//" does not have \
                a DQ plane.", l_logfile, verbose+)
            l_fl_dqprop=no
        }
        nfiles+=1

SKIP_CHECKS:
    } # end while loop

    # check for empty file list
    if (nfiles==0 && nbad == 0) {
        printlog("ERROR - NISKY: No input images meet wildcard criteria.",
            l_logfile, verbose+)
        status=1
        goto clean
    }

    # Exit if problems found with input files
    if (nbad > 0) {
        printlog("ERROR - NISKY: "//nbad//" image(s) either do not exist, \
            are not MEF files, or", l_logfile, verbose+)
        printlog("                 have not been run through *PREPARE.",
            l_logfile, verbose+)
        status=1
        goto clean
    } # end if (nbad > 0)

    printlog("Using input files:", l_logfile, l_verbose)
    if (l_verbose) type(tmpfile1)
    type(tmpfile1, >> l_logfile)
    printlog("Output image: "//l_outimage, l_logfile, l_verbose)

    #------------------------------------------------------------------------
    # strip out science and dq planes
    l_expone=0.0 ; l_ronone=0.0 ; l_gainone=0.0
    n=0
    scanfile = tmpfile1
    while (fscan(scanfile,img) != EOF) {
        n=n+1

        # strip suffix if present
        suf = substr(img,strlen(img)-3,strlen(img))
        if (substr(img,strlen(img)-4,strlen(img)) == "["//l_sci_ext//"]")
            img=substr(img,1,(strlen(img)-5))
        if (substr(img,strlen(img)-2,strlen(img)) == "[1]")
            img=substr(img,1,(strlen(img)-3))
        if (substr(img,strlen(img)-4,strlen(img)) == ".fits")
            img=substr(img,1,(strlen(img)-5))

        # check exposure time
        imgets(img//"[0]",l_key_exptime, >& "dev$null")
        if (imgets.value == "0") {
            printlog("ERROR - NISKY: Image header parameter not found \
                ("//l_key_exptime//")", l_logfile, verbose+)
            status=1
            goto clean
        }
        if (n == 1)
            l_expone=real(imgets.value)
        else
            if (abs(real(imgets.value)-l_expone) > 0.1) {
                printlog("WARNING - NISKY: Exposure times are significantly \
                    different.  Continuing.", l_logfile, verbose+)
            }
 
        # check for gain, readnoise, and saturation values
        imgets(img//"[0]",l_key_ron, >& "dev$null")
        if ((imgets.value == "0") || (imgets.value == "")) {
            printlog("ERROR - NISKY: Could not get read noise from header of \
                image "//img, l_logfile, verbose+)
            status=1
            goto clean
        }
        if (n == 1)
            l_ronone=real(imgets.value)
        else
            if (abs(real(imgets.value)-l_ronone) > 1.) {
                printlog("WARNING - NISKY: read noise values are different.  \
                    Continuing, but the", l_logfile, verbose+)
                printlog("                 read noise in the output header \
                    will be wrong.", l_logfile, verbose+)
            }

        imgets(img//"[0]",l_key_gain, >& "dev$null")
        if ((imgets.value == "0") || (imgets.value == "")) {
            printlog("ERROR - NISKY: Could not get gain from header of \
                image"//img, l_logfile, verbose+)
            status=1
            goto clean
        } else 
            l_gain=real(imgets.value)

        if (n == 1)
            l_gainone=real(imgets.value)
        else
            if (abs(real(imgets.value)-l_gainone) > 0.5) {
                printlog("WARNING - NISKY: gain values are different.  \
                    Continuing, but the gain in", l_logfile, verbose+)
                printlog("                 the output header will be wrong.",
                    l_logfile, verbose+)
            }

        print(img, >>tmpfile2)

        # science extension
        print(img//"["//l_sci_ext//"]", >> scilist)

        # DQ extension
        if (l_fl_vardq) {
            if (!imaccess(img//"["//l_dq_ext//"]")) {
                printlog("WARNING - NISKY: No DQ plane for "//img//", \
                    creating empty one.", l_logfile, verbose+)
                imarith(img//"["//l_sci_ext//"]", "*", "0.0", \
                    tmpdq//"_"//n//".pl", verbose-)
            } else
                imcopy(img//"["//l_dq_ext//"]",tmpdq//"_"//n//".pl", verbose-)
            if (l_fl_dqprop) 
                print(tmpdq//"_"//n//".pl", >> dqlist) 
            gemhedit (img//"["//l_sci_ext//"]", "BPM", tmpdq//"_"//n//".pl", \
                "", delete-)
         }
    } # end of while loop
    scanfile=""

    l_ronone = l_ronone/sqrt(n)
    l_ronone = real(int(l_ronone*10.))/10.
    l_gainone = l_gainone*n

    if (n == 1) {
        printlog("ERROR - NISKY: Cannot combine a single image.", 
            l_logfile, verbose+)
        status=1
        goto clean
    } else if (n == 0) {
        printlog("ERROR - NISKY: No images to combine.", l_logfile, verbose+)
        status=1
        goto clean
    }

    #--------------------------------------------------------------------------
    # Check for mask images that already exist

    if (l_fl_nifastsky) {
        needtoreduce=no
        i=1
        scanfile=tmpfile2
        while (fscan(scanfile,img) != EOF) { 
            if (!access(img//l_masksuffix//".pl")) 
                needtoreduce=yes
        }
        scanfile=""
    } else
        needtoreduce=no

    #--------------------------------------------------------------------------
    # Make a quick temporary sky image using nifastsky.cl and flat field

    if (needtoreduce) {
        printlog("NISKY calling NIFASTSKY", l_logfile, l_verbose)
        nifastsky("@"//tmpfile2, outimage=tmpskymef, outtitle="default",
            combtype="default", sci_ext=l_sci_ext, var_ext=l_var_ext,
            dq_ext=l_dq_ext, fl_vardq=no, fl_dqprop=no,
            key_exptime=l_key_exptime, verbose=no, logfile=l_logfile)
        if (nifastsky.status != 0) {
            printlog("ERROR - NISKY: Error in NIFASTSKY prevents NISKY \
                completion.", l_logfile, verbose+)
            status=1
            goto clean
        }
        printlog("Returning to NISKY",l_logfile, l_verbose)

        imcopy(tmpskymef//"["//l_sci_ext//"]",tmpsky, verbose-)
        imstat(tmpsky//l_statsec,fi="midpt,stddev",lower=INDEF,upper=INDEF,
            nclip=0,lsigma=INDEF,usigma=INDEF,binwidth=0.05,format-,cache-) | \
            scan(l_mean,l_sig)
        imstat(tmpsky//l_statsec,fi="mean",lower=(l_mean-4*l_sig),
            upper=(l_mean+4*l_sig),nclip=0,lsigma=INDEF,usigma=INDEF,
            binwidth=0.05,format-,cache-) | scan(l_mean)
        if (l_mean==INDEF) {
            printlog("ERROR - NISKY: Statistics failed, possibly due to a \
                bad statsec.", l_logfile, verbose+)
            status=1
            goto clean
        }

        l_skymean=l_mean
        if (l_mean <= 0.) l_mean=0.0001
        if (sqrt(1./(l_mean*n*l_gain))<=0.02) {
            l_temp="((a>"//(l_mean+4.*l_sig)//") || \
                (a<"//(l_mean-4.*l_sig)//")) ? 1. : (a/"//l_mean//")"
            imexpr(l_temp,tmpflat,tmpsky,outtype="real", >& "dev$null")
            imreplace(tmpflat, 0.0001, upper=0.0001) # added by KAO
        } else
            imexpr("a*0.+1.",tmpflat,tmpsky,outtype="real", >& "dev$null")

    } # end if(needtoreduce)

    #--------------------------------------------------------------------------
    # Reduce the images, identify objects, and make a mask for each input file
    # (if the object mask file already exists then that file is used and task
    # continues)

    scanfile=tmpfile2

    # Start mask-making loop here
    i=1
    while (fscan(scanfile,img) != EOF) { 
        # Create tmp FITS file names re-used within this loop
        tmpimg = mktemp("tmpim")

        imgsci=img//"["//l_sci_ext//"]"
        if (l_fl_vardq) imgdq=tmpdq//"_"//i//".pl"

        # Check for existing masks
        if (!access(img//l_masksuffix//".pl")) { 
            printlog("Making mask for image "//img, l_logfile, l_verbose)

            # Reduce the input image using the fast sky and flat constructed
            # above
            if (needtoreduce) {
                imstat (imgsci//l_statsec, fields="midpt,stddev", lower=INDEF,
                    upper=INDEF, nclip=0, lsigma=INDEF, usigma=INDEF,
                    binwidth=0.05, format-, cache-) | scan (l_mean, l_sig)
                imstat (imgsci//l_statsec, fields="midpt",
                    lower=(l_mean-4.*l_sig), upper=(l_mean+4.*l_sig), nclip=0,
                    lsigma=INDEF, usigma=INDEF, binwidth=0.05, format-,
                    cache-) | scan(skyconst)
                if (skyconst==INDEF) {
                    printlog("ERROR - NISKY: Statistics failed, possibly due \
                        to a bad statsec.", l_logfile, verbose+)
                    status=1
                    goto clean
                }
                imexpr ("((a-(b*"//(skyconst/l_skymean)//"))/c)+"//skyconst,
                    tmpimg, imgsci, tmpsky, tmpflat, rangecheck=yes,
                    outtype="real", >& "dev$null")
            } else
                imcopy (imgsci, tmpimg, verbose-, >& "dev$null")
                
            if (l_fl_vardq)
                proto.fixpix(tmpimg,imgdq,linterp="INDEF",cinterp="INDEF",
                    >& "dev$null")
    
            # Identify stars in reduced sky images using OBJMASKS
            objmasks(tmpimg, img//l_masksuffix//".pl", omtype="boolean",
                hsigma=l_threshold, minpix=l_minpix, ngrow=l_ngrow, 
                agrow=l_agrow, >& "dev$null")

    ### KL
    ### A floating point error sometimes occurs after this point in the code.
    ### It looks like a classic process buffer corruption cured by 'flpr'.
    ### I do not know where the corruption takes place, but it manifests itself
    ### through the imarith call below.
    ###
    ###  Adding a 'flpr' here.  Hopefully a temporary solution.

            flpr

    ### Commented out since the source of this problem will hopefully be fixed
    ### by the imreplace on line 471 - EH
    ###
    ### Uncommented since it appears to be required for niflat / nisky to 
    ### run without floating point errors - EH

            # Include the original DQ plane
            if (l_fl_vardq)
                imarith(imgdq,"max",img//l_masksuffix//".pl",
                    img//l_masksuffix//".pl",pixtype="short", verbose-)

            # end if(!access(img//l_masksuffix//".pl"))
        } else
            printlog("WARNING - NISKY: Using existing mask file "//\
                img//l_masksuffix//".pl", l_logfile, verbose+)

        # edit header to include mask file as BPM
        gemhedit (img//"["//l_sci_ext//"]", "BPM", img//l_masksuffix//".pl", \
            "", delete-)
        gemhedit (img//"[0]","NISKYMSK", img//l_masksuffix//".pl",
            "Object mask generated by NISKY", delete-)

        # get rid of the temporary reduced image and the coordinate files
        delete(tmpimg//".coo", verify-, >& "dev$null")
        delete(tmpimg//"_coo", verify-, >& "dev$null")
        delete(tmpimg//"_see", verify-, >& "dev$null")
        delete(tmpimg//"_see.tab", verify-, >& "dev$null")
        imdelete(tmpimg, verify-, >& "dev$null")

        i+=1
    } # end while loop
    scanfile=""

    #---------------------------------------------------------------------
    # imcombine raw sky frames with star masks as bad pixel masks

    # save the user's parameters for imcombine
    delete("uparm$imhimcome.par.org", verify-, >& "dev$null")
    if (access("uparm$imhimcome.par"))
        copy("uparm$imhimcome.par","uparm$imhimcome.par.org", verbose-)

    cache("imcombine")

    # Set imcombine parameters
    imcombine.headers=""
    imcombine.bpmasks=""
    imcombine.rejmasks=""
    imcombine.nrejmasks = combdq
    imcombine.expmasks = ""
    imcombine.sigmas = combsig
    imcombine.logfile = l_logfile
    # imcombine.combine defined below
    # imcombine.reject defined below
    imcombine.project = no
    imcombine.outtype = "real"
    imcombine.outlimits = ""
    imcombine.offsets = "none"
    imcombine.masktype = "goodvalue"
    imcombine.maskvalue = 0
    imcombine.blank = 0.
    imcombine.scale = "median"
    imcombine.zero = "none"
    imcombine.weight = "none"
    imcombine.statsec = l_statsec
    imcombine.expname = ""
    imcombine.lthreshold = INDEF
    imcombine.hthreshold = INDEF
    # imcombine.nlow defined below
    # imcombine.nhigh defined below
    imcombine.nkeep = 1
    imcombine.mclip=yes
    # imcombine.lsigma defined below
    # imcombine.hsigma defined below
    imcombine.sigscale=0.1
    imcombine.grow = 0

    if (l_combtype=="default") {
        imcombine.combine = "average"
        if (n == 2) {
            imcombine.reject = "minmax"
            imcombine.nlow=0
            imcombine.nhigh=1
            printlog("WARNING - NISKY: Combining two images by taking the \
                minimum.", l_logfile, verbose+)
        } else {
            imcombine.reject = "avsigclip"
            imcombine.lsigma = 3.
            imcombine.hsigma = 3.
        }
    } else { 
        imcombine.combine = l_combtype
        imcombine.reject= l_rejtype
        imcombine.nlow = l_nlow
        imcombine.nhigh = l_nhigh
        imcombine.lsigma = l_lsigma
        imcombine.hsigma = l_hsigma
        if (n < 5) {
            printlog("WARNING - NISKY: Combining 4 or fewer images using "//\
                l_combtype, l_logfile, verbose+)
            if (l_rejtype == "minmax")
                printlog("                 with "//l_nlow//" low and "//\
                    l_nhigh//" high pixels rejected.", l_logfile, verbose+)
            else if (l_rejtype == "avsigclip")
                printlog("                 with "//l_lsigma//"=lower sigma \
                    and "//l_hsigma//"=upper sigma.", l_logfile, verbose+)
            else
                printlog("                 with no pixels rejected.",
                    l_logfile, verbose+)
        }
        if ((n <= (l_nlow+l_nhigh)) && (l_rejtype=="minmax")) {
            printlog("ERROR - NISKY: Cannot reject more pixels than the \
                number of images.", l_logfile, verbose+)
            status=1
            goto clean
        }
    } #end not-default section

    # Do the combining (supress the output from imcombine)
    imcombine("@"//scilist,l_outimage//"["//l_sci_ext//",append]")

    printlog("Combining "//str(n)//" images, using "//imcombine.combine,
        l_logfile, l_verbose)
    printlog("Rejection type is "//imcombine.reject,l_logfile, l_verbose)
    if (imcombine.reject=="minmax")
        printlog("with "//imcombine.nlow//" low and "//imcombine.nhigh//\
            " high values rejected.", l_logfile, l_verbose)
    if (imcombine.reject=="avsigclip")
        printlog("with lsigma="//imcombine.lsigma//" and hsigma="//\
            imcombine.hsigma, l_logfile, l_verbose)

    if (l_fl_vardq) {
        # make variance image by squaring combsig
        imarith(combsig,"*",combsig,l_outimage//"["//l_var_ext//",append]",
            pixtype="real", verbose-)
        # bad pixels have val=n
        imexpr("(a=="//n//") ? 1 : 0", tmpdqcomb, combdq, outtype="ushort", \
            >& "dev$null")
        if (!l_fl_dqprop) {
            chpixtype (tmpdqcomb, l_outimage//"["//l_dq_ext//",append]", \
                newpixtype="ushort", oldpixtype="all", verbose-)
        } else {
            scanfile=dqlist
            for (j=1; j<=n; j+=1) {
                # img is the input DQ extension (added to dqlist above)
                i=fscan(scanfile, img)

                # Create tmp file names used within this loop
                dqsumold = mktemp("tmpdqsumold")

                if (j==1)
                    imcopy(tmpdqcomb, dqsumold, verbose-)
                else
                    imrename(dqsum, dqsumold, verbose-)

                addmasks(dqsumold//","//img, dqsum//".pl", "im1 || im2", \
                    flags=" ")
                imdelete(dqsumold, verify-, >& "dev$null")
            }
            scanfile=""

            # Append the summed DQ extension (i.e., the output mask from
            # imcombine plus all the input DQ extensions) to the output image.
            chpixtype(dqsum, l_outimage//"["//l_dq_ext//",append]", \
                newpixtype="ushort", oldpixtype="all", verbose-)
        }

        # update headers
        gemhedit (l_outimage//"["//l_sci_ext//"]", "BPM", "", "", delete=yes)
        gemhedit (l_outimage//"["//l_var_ext//"]", "BPM", "", "", delete=yes)
        gemhedit (l_outimage//"["//l_dq_ext//"]", "BPM", "", "", delete=yes)
        gemhedit(l_outimage//"["//l_var_ext//"]", "EXTVER", 1, \
            "Extension version", delete-)
        gemhedit(l_outimage//"["//l_dq_ext//"]", "EXTVER", 1, \
            "Extension version", delete-)
    } # end if(l_fl_vardq)

    #--------------------------------------------------------------------------
    # update headers in final image PHU

    gemdate ()
    gemhedit (l_outimage//"[0]", "GEM-TLM", gemdate.outdate, 
        "UT Last modification with GEMINI", delete-)
    printf("%-8s= \'%-18s\' / %-s\n","NISKY",gemdate.outdate,"UT Time \
    stamp for NISKY", >> tmphead)

    if (l_outtitle=="default" || l_outtitle=="" || l_outtitle==" ")
        gemhedit (l_outimage//"[0]", "i_title", \
            "SKY IMAGE from gemini.niri.nisky", "Image title", delete-)
    else
        gemhedit (l_outimage//"[0]", "i_title", l_outtitle, "Image title", \
            delete-)

    l_ronone=(int(l_ronone*10000.))/10000.
    l_gainone=(int(l_gainone*10.))/10.
    gemhedit (l_outimage//"[0]", l_key_ron, l_ronone, "", delete-)
    gemhedit (l_outimage//"[0]", l_key_gain, l_gainone, "", delete-)

    printf("%-8s= \'%-18s\' / %-s\n","NISKYSTA",l_statsec,
        "Statistics region used by NISKY", >> tmphead)

    printf("%-8s= %20.2f / %-s\n","NISKYTHR",l_threshold,
        "Object detection threshold used by NISKY", >> tmphead)
    printf("%-8s= %20.0f / %-s\n","NISKY",l_minpix,
        "Minimum pixels per object used by NISKY", >> tmphead)
    printf("%-8s= %20.0f / %-s\n","NISKY",l_ngrow,
        "Iterations to extend objects used by NISKY", >> tmphead)
    printf("%-8s= %20.2f / %-s\n","NISKY",l_agrow,
        "Area limit for extending objects used by NISKY", >> tmphead)
    printf("%-8s= \'%-18s\' / %-s\n","NISKYCOM",imcombine.combine,
        "Type of combine used by NISKY", >> tmphead)
    printf("%-8s= \'%-18s\' / %-s\n","NISKYREJ",imcombine.reject,
        "Type of rejection used by NISKY", >> tmphead)
    if (imcombine.reject=="minmax") {
        printf("%-8s= %20.0f / %-s\n","NISKYNLO",imcombine.nlow,
            "Low pixels rejected (minmax)", >> tmphead)
        printf("%-8s= %20.0f / %-s\n","NISKYNHI",imcombine.nhigh,
            "High pixels rejected (minmax)", >> tmphead)
    }
    if (imcombine.reject=="avsigclip") {
        printf("%-8s= %20.2f / %-s\n","NISKYLSI",imcombine.lsigma,
            "Lower sigma for rejection", >> tmphead)
        printf("%-8s= %20.2f / %-s\n","NISKYHSI",imcombine.hsigma,
            "Upper sigma for rejection", >> tmphead)
        printf("%-8s= %20.0f / %-s\n","NISKYNKE",imcombine.nkeep,
            "Min(max) number of pixels to keep(reject)", >> tmphead)
        printf("%-8s= \'%-18b\' / %-s\n","NISKYMCL",imcombine.mclip,
            "Use median in clipping algorithms", >> tmphead)
        printf("%-8s= %20.2f / %-s\n","NISKYSSC",imcombine.sigscale,
            "Tolerance for sigma clip scaling correction", >> tmphead)
    }
    printf("%-8s= \'%-18b\' / %-s\n","NISKYDQP",l_fl_dqprop,
        "Retain input data quality planes?", >> tmphead)

    # Put input image names in header
    i=1
    scanfile=tmpfile1
    while (fscan(scanfile,img) != EOF) {
        printf("%-8s= \'%-18s\' / %-s\n","NISKIM"//str(i),img,
            "Input image combined with NISKY", >> tmphead)
        i=i+1
    }

    # put all the new stuff in the header
    mkheader(l_outimage//"[0]",tmphead,append+, verbose-)

    scanfile=""

    #-------------------------------------------------------------------------
clean:
    if(status==0)
        printlog("NISKY exit status: good.",l_logfile, l_verbose)
    else if (status!=2)
        imdelete(l_outimage, verify-, >& "dev$null")
    
    printlog("------------------------------------------------------------\
        ----------------", l_logfile, l_verbose)

    if (access(tmpfile2)) {
        scanfile=tmpfile2
        while (fscan(scanfile,img) != EOF) {
            gemhedit (img//"["//l_sci_ext//"]", "BPM", "", "", delete=yes)
            if (!l_fl_keepmasks)
                imdelete(img//l_masksuffix//".pl", verify-, >& "dev$null")
        }
    }

    delete(tmpfile1, verify-, >& "dev$null")
    delete(tmpfile2, verify-, >& "dev$null")
    delete(scilist//","//dqlist//","//tmphead, verify-, >& "dev$null")
    imdelete(combsig//","//combdq, verify-, >& "dev$null")
    imdelete(dqsum, verify-, >& "dev$null")
    imdelete(tmpdqcomb, verify-, >& "dev$null")
    imdelete(tmpdq//"_*.pl", verify-, >& "dev$null")
    imdelete(tmpsky, verify-, >& "dev$null")
    imdelete(tmpskymef, verify-, >& "dev$null")
    imdelete(tmpflat, verify-, >& "dev$null")
    unlearn("imcombine")
    # restore the user's parameters for imcombine
    if (access("uparm$imhimcome.par.org"))
        rename("uparm$imhimcome.par.org","uparm$imhimcome.par",field="all")

end

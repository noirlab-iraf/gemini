# Copyright(c) 2001-2012 Association of Universities for Research in Astronomy, Inc.

procedure nifastsky(inimages)

# Make a quick-and-dirty median sky image for NIRI images.
# NOTE: combtype=default overrides the user's settings of nlow and nhigh
#
# Version  June 19, 2001  JJ, beta-release
#          Oct 12, 2001   JJ  v1.2 release
#          Nov 16, 2001   JJ  automatic output filename
#          Feb 28, 2002   JJ  v1.3 release
#          Jun 7, 2002    JJ  direct to MEF, fixed scilist, multiple identical
#                             input files checked, headers written at once
#          Jun 12, 2002   JJ  fixed bug when input file doesn't exist
#          Sep 10, 2002   IJ parameter encoding
#          Sept 20, 2002  JJ v1.4 release
#          Aug 18, 2003   KL  IRAF2.12 - new/modified parameters
#                               hedit: addonly
#                               imcombine: headers,bpmasks,expmasks,outlimits.
#                                          rejmask->rejmasks, plfile->nrejmasks

char    inimages     {prompt="Raw NIRI images to combine"}                      # OLDP-1-input-primary-combine-suffix=_sky
char    outimage     {"",prompt="Output sky image"}                             # OLDP-1-output
char    outtitle     {"default",prompt="Title for output image"}                # OLDP-3
char    key_exptime  {"EXPTIME",prompt="Keyword for exposure time"}             # OLDP-3
char    combtype     {"default",prompt="Type of combine operation",enum="default|median|average"} # OLDP-3
char    statsec      {"[100:924,100:924]",prompt="Statistics section"}          # OLDP-2
char    rejtype      {"minmax",prompt="Type of rejection", enum="none|minmax"}  # OLDP-3
char    logfile      {"",prompt="Name of log file"}                             # OLDP-1
int     nlow         {0,min=0,prompt="minmax: Number of low pixels to reject"}  # OLDP-3
int     nhigh        {1,min=0,prompt="minmax: Number of high pixels to reject"} # OLDP-3
char    sci_ext      {"SCI",prompt="Name or number of science extension"}       # OLDP-3
char    var_ext      {"VAR",prompt="Name or number of variance extension"}      # OLDP-3
char    dq_ext       {"DQ",prompt="Name or number of data quality extension"}   # OLDP-3
char    key_ron      {"RDNOISE",prompt="Header keyword for read noise (e-)"}    # OLDP-3
char    key_gain     {"GAIN",prompt="Header keyword for gain (e-/ADU)"}         # OLDP-3
bool    fl_vardq     {no,prompt="Create variance and data quality frames in output?"}  # OLDP-3
bool    fl_dqprop    {no,prompt="Retain input DQ information in output?"}       # OLDP-3
bool    verbose      {yes,prompt="Verbose actions"}                             # OLDP-4
int     status       {0,prompt="Exit status (0=good)"}                          # OLDP-4
struct  *scanfile    {"",prompt="Internal use only"}                            # OLDP-4

begin

    char    l_inimages, l_outimage, l_combtype, l_rejtype, l_logfile
    char    l_key_exptime, l_outtitle, l_key_filter, l_statsec
    int     l_nlow, l_nhigh, i, n, nbad, len, nfiles
    bool    l_verbose, l_fl_dqprop, l_fl_vardq
    char    temp, img, firstimage, tmpfile1
    struct  l_struct
    real    l_expone, l_ron, l_gain
    char    l_sci_ext, l_var_ext, l_dq_ext, l_key_ron, l_key_gain
    char    scilist, dqlist, combsig, combdq
    char    tmpsci, tmpdq, tmpdqcomb, dqsumold, dqsum, suf, tmphead, keyfound
    char    sformat, fformat

    # Set local variables
    l_inimages = inimages ; l_outimage = outimage ; l_key_exptime = key_exptime
    l_outtitle = outtitle
    l_combtype = combtype ; l_rejtype = rejtype
    l_logfile = logfile ; l_verbose=verbose
    l_nlow = nlow ; l_nhigh = nhigh
    l_sci_ext = sci_ext ; l_dq_ext = dq_ext ; l_var_ext=var_ext
    l_fl_dqprop=fl_dqprop ; l_fl_vardq=fl_vardq
    l_key_ron=key_ron ; l_key_gain=key_gain
    l_statsec=statsec

    status=0

    #------------------------------------------------------------------------
    # Check for package log file or user-defined log file
    cache("niri", "gemdate")
    print(l_logfile) | scan(l_logfile)
    if (l_logfile=="" || l_logfile==" ") {
        l_logfile=niri.logfile
        print(l_logfile) | scan(l_logfile)
        if (l_logfile=="" || l_logfile==" ") {
            l_logfile="niri.log"
            printlog("WARNING - NIFASTSKY: Both nifastsky.logfile and \
                niri.logfile are empty.", l_logfile, verbose+)
            printlog("                     Using default file niri.log.",
                l_logfile, verbose+)
        }
    }

    # Open log file
    date | scan(l_struct)
    printlog("------------------------------------------------------------\
        ----------------", l_logfile, l_verbose)
    printlog("NIFASTSKY -- "//l_struct, l_logfile, l_verbose)
    printlog(" ", l_logfile, l_verbose)

    # Make temporary files
    tmpfile1  = mktemp("tmpfl")
    scilist   = mktemp("tmpscilist")
    dqlist    = mktemp("tmpdqlist")
    dqsum     = mktemp("tmpdqsum")
    tmpsci    = mktemp("tmpsci")
    tmpdq     = mktemp("tmpdq")
    combsig   = mktemp("tmpcombsig") 
    combdq    = mktemp("tmpcombdq")
    tmpdqcomb = mktemp("tmpdqcomb")
    tmphead   = mktemp("tmphead")

    # Keep imgets parameters from changing by outside world
    cache ("imgets")

    # Verify that the extension names are not empty, otherwise exit gracefully
    print(l_sci_ext) | scan(l_sci_ext)
    if (l_sci_ext=="" || l_sci_ext==" ") {
        printlog("ERROR - NIFASTSKY: extension name sci_ext is missing.",
            l_logfile, verbose+)
        status=1
        goto clean
    }
    if (l_fl_vardq) {
        print(l_dq_ext) | scan(l_dq_ext)
        print(l_var_ext) | scan(l_var_ext)
        if (l_dq_ext=="" || l_dq_ext==" " || l_var_ext=="" || l_var_ext==" ") {
            printlog("WARNING - NIFASTSKY: extension name var_ext or dq_ext \
                is missing.", l_logfile, verbose+)
            printlog("                     Output image will not have VAR or \
                DQ planes.", l_logfile, verbose+)
            l_fl_vardq=no
            l_fl_dqprop=no
        } 
    }

    # check existence of list file
    if (substr(l_inimages,1,1) == "@") {
        len=strlen(l_inimages)
        if (!access(substr(l_inimages,2,len))) {
            printlog("ERROR - NIFASTSKY: Input file "//\
                substr(l_inimages,2,len)//" does not exist.", 
                l_logfile, verbose+)
            status=1
            goto clean
        }
    }
    # Put all images in a temporary file: tmpfile1
    files(l_inimages,sort-) | unique("STDIN", > tmpfile1)

    # Verify that input images actually exist. (Redundant for * though)
    # at the same time check the exposure times and MEF format
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
                printlog("ERROR - NIFASTSKY: Output file "//l_outimage//\
                    " already exists.", l_logfile, verbose+)
                status=2
                goto clean
            }
        }

        # check input images
        gimverify(img)
        if (gimverify.status==1) {
            printlog("ERROR - NIFASTSKY: File "//img//" not found.",
                l_logfile, verbose+)
            nbad+=1
        } else if(gimverify.status>1) {
            printlog("ERROR - NIFASTSKY: File "//img//" not a MEF FITS image.",
                l_logfile, verbose+)
            nbad+=1
        } else {
            keyfound=""
            hselect(img//"[0]","*PREPAR*",yes) | scan(keyfound)
            if (keyfound == "") {
                printlog("ERROR - NIFASTSKY: Image "//img//" not *PREPAREd.",
                    l_logfile, verbose+)
                nbad+=1
            }
        }

        # get the header of the first image and prepare for output image
        if (nfiles==0) imcopy(img//"[0]", l_outimage, verbose-)

        if (l_fl_vardq && l_fl_dqprop && !imaccess(img//"["//l_dq_ext//"]") ) {
            printlog("WARNING - NIFASTSKY: Cannot propagate input DQ planes \
                because", l_logfile, verbose+)
            printlog("                     input image "//img//" does not \
                have a DQ plane.", l_logfile, verbose+)
            l_fl_dqprop=no
        }
        nfiles+=1
    } #end of while loop
    scanfile=""

    # check for empty file list
    if (nfiles==0) {
        printlog("ERROR - NIFASTSKY: No input images meet wildcard criteria.",
            l_logfile, verbose+)
        status=1
        goto clean
    }

    # Exit if problems found
    if (nbad > 0) {
        printlog("ERROR - NIFASTSKY: "//nbad//" image(s) either do not \
            exist, are not MEF files, or", l_logfile, verbose+)
        printlog("                   have not been run through *PREPARE.",
            l_logfile, verbose+)
        status=1
        goto clean
    }

    printlog("Using input files:", l_logfile, l_verbose)
    if (l_verbose) type(tmpfile1)
    type(tmpfile1, >> l_logfile)
    printlog("Output image: "//l_outimage, l_logfile, l_verbose)

    #-------------------------------------------------------------------------
    # Files OK, so continue
    l_expone=0.0
    n=0
    scanfile = tmpfile1
    while (fscan(scanfile,img) != EOF) {
        n=n+1
        suf = substr(img,strlen(img)-4,strlen(img))

        # Strip the suffix of the input file:
        if (suf!=".fits" && imaccess(img//".fits"))
            img = img//".fits"
        if (n==1) 
            firstimage=img

        # check exposure time
        imgets(img//"[0]",l_key_exptime, >& "dev$null")
        if (imgets.value == "0" || imgets.value == " ") {
            printlog("ERROR - NIFASTSKY: Image header parameter not found \
                (key_exptime)", l_logfile, verbose+)
            status=1
            goto clean
        }
        if (n == 1)
            l_expone=real(imgets.value)
        else {
            if (abs(real(imgets.value)-l_expone) > 0.1)
                printlog("WARNING - NIFASTSKY: exposure times are \
                    significantly different.  Continuing.", l_logfile, \
                    verbose+)
        }

        # check read noise
        imgets(img//"[0]",l_key_ron, >& "dev$null")
        if (imgets.value == "0" || imgets.value == " ") {
            printlog("ERROR - NIFASTSKY: Image header parameter not found \
                (key_ron)", l_logfile, verbose+)
            status=1
            goto clean
        }
        if (n == 1)
            l_ron=real(imgets.value)
        else {
            if (abs(real(imgets.value)-l_ron) > 1.) {
                printlog("WARNING - NIFASTSKY: read noise values are \
                    different.  Continuing, but the", l_logfile, verbose+)
                printlog("                     read noise in the output \
                    header will be wrong.", l_logfile, verbose+)
            }
        }

        # check gain
        imgets(img//"[0]",l_key_gain, >& "dev$null")
        if (imgets.value == "0" || imgets.value == " ") {
            printlog("ERROR - NIFASTSKY: Image header parameter not found \
                (key_gain)", l_logfile, verbose+)
            status=1
            goto clean
        }
        if (n == 1)
            l_gain=real(imgets.value)
        else {
            if (abs(real(imgets.value)-l_gain) > 0.5) {
                printlog("WARNING - NIFASTSKY: gain values are different.  \
                    Continuing, but the gain in", l_logfile, verbose+)
                printlog("                     the output header will be \
                    wrong.", l_logfile, verbose+)
            }
        }

        # science extension
        #imcopy(img//"["//l_sci_ext//"]",tmpsci//"_"//n//".fits", verbose-)
        #print(tmpsci//"_"//n//".fits", >> scilist)
        print(img//"["//l_sci_ext//"]", >> scilist)


        # DQ extension
        if (imaccess(img//"["//l_dq_ext//"]") && l_fl_dqprop) {
            imcopy(img//"["//l_dq_ext//"]",tmpdq//"_"//n//".pl", verbose-)
            print(tmpdq//"_"//n//".pl", >> dqlist) 
            #hedit(tmpsci//"_"//n//".fits","BPM",tmpdq//"_"//n//".pl",add+,
            #    addonly-,delete-,verify-,show-,update+)
        }
    } # end of while loop
    scanfile=""

    l_ron = l_ron/sqrt(n)
    l_ron = real(int(l_ron*10.))/10.
    l_gain = l_gain*n

    if (n == 1) {
        printlog("ERROR - NIFASTSKY: Cannot combine a single image.",
            l_logfile, verbose+)
        status=1
        goto clean
    } else if (n == 0) {
        printlog("ERROR - NIFASTSKY: No images to combine.", l_logfile, \
            verbose+)
        status=1
        goto clean
    }

    # save the user's parameters for imcombine
    delete("uparm$imhimcome.par.org", verify-, >& "dev$null")
    if (access("uparm$imhimcome.par"))
        copy("uparm$imhimcome.par","uparm$imhimcome.par.org", verbose-)

    cache("imcombine")

    # Set imcombine parameters
    imcombine.headers=""
    imcombine.bpmasks=""
    imcombine.rejmasks=""
    imcombine.nrejmasks=combdq
    imcombine.expmasks=""
    imcombine.sigmas=combsig
    imcombine.logfile = "STDOUT" # so that it can be tee'd later
    imcombine.project = no
    # imcombine.combine defined below
    # imcombine.reject defined below
    imcombine.outtype = "real"
    imcombine.outlimits=""
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
    imcombine.nkeep = 0
    imcombine.grow = 0

    if (l_combtype=="default") {
        if (n < 5) {
            imcombine.combine = "average"
            imcombine.reject = "minmax"
            imcombine.nlow=0
            imcombine.nhigh=1
            printlog("WARNING - NIFASTSKY: Averaging 4 or fewer images \
                with 1 high pixel rejected", l_logfile, verbose+)
        } else if (n < 8) {
            imcombine.combine = "median"
            imcombine.reject = "minmax"
            imcombine.nlow=1
            imcombine.nhigh=1
        } else {
            imcombine.combine = "median"
            imcombine.reject = "minmax"
            imcombine.nlow=1
            imcombine.nhigh=2
        } 
    } else {
        imcombine.combine = l_combtype
        imcombine.reject= l_rejtype
        imcombine.nlow = l_nlow
        imcombine.nhigh = l_nhigh
        if (n < 5) {
            printlog("WARNING - NIFASTSKY: Combining 4 or fewer images \
                using "//l_combtype, l_logfile, verbose+)
            if (l_rejtype != "none") {
                printlog("                      with "//l_nlow//\
                    " low and "//l_nhigh//" high pixels rejected.",
                    l_logfile, verbose+)
            } else {
                printlog("                      with no pixels rejected.",
                    l_logfile, verbose+)
            }
        } # end if(n<5)
        if ((n <= (l_nlow+l_nhigh)) && (l_rejtype=="minmax")) {
            printlog("ERROR - NIFASTSKY: Cannot reject more pixels than \
                the number of images.", l_logfile, verbose+)
            status=1
            goto clean
        }
    } # end not-default section

    # Do the combine, quietly
    imcombine("@"//scilist,l_outimage//"["//l_sci_ext//",append]",
        >& "dev$null")

    printlog("Combining "//str(n)//" images, using "//imcombine.combine,
        l_logfile, l_verbose)
    printlog("Rejection type is "//imcombine.reject, l_logfile, l_verbose)
    if (imcombine.reject=="minmax") {
        printlog("with "//imcombine.nlow//" low and "//imcombine.nhigh//\
            " high values rejected.", l_logfile, l_verbose)
    } 

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

                # Create tmp file name used only within this loop
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

    #-----------------------------------------------------------------------
    # Fix up header
    gemdate()

    # Set printf formats. This is the syntax for 'nhedit' command file 
    sformat = "%-8s \'%-18s\' \'%-s\'\n"
    fformat = "%-8s %20.0f \'%-s\'\n"

    # Setup default parameters for nhedit command file 'tmphead'.
    # Notice that this particular printf to write 'default_pars'
    # requires single quote. The double quotes are require in 
    # 'after' and 'before' parameter values.
    printf ('default_pars add+,after="",before="",addonly-,\
        delete-,verify-,show-,update+\n',>>tmphead)

    printf (sformat,"GEM-TLM",gemdate.outdate,
        "UT Last modification with GEMINI", >> tmphead)
    printf (sformat,"NIFSTSKY",gemdate.outdate,
        "UT Time stamp for NIFASTSKY", >> tmphead)

    l_ron=(int(l_ron*10000.))/10000.
    l_gain=(int(l_gain*10.))/10.
    printf (fformat,l_key_ron,l_ron,"", >> tmphead)
    printf (fformat,l_key_gain,l_gain,"", >> tmphead)

    if (l_outtitle=="default") {
        printf (sformat,"i_title", "SKY IMAGE from gemini.niri.nifastsky",
            "Image title", >> tmphead)
    } else {
        printf (sformat,"i_title",l_outtitle,"Image title", >> tmphead)
    }
    printf (sformat,"NIFSTSTA",l_statsec,
        "Statistics region used by NIFASTSKY", >> tmphead)
    printf (sformat,"NIFSTCOM",imcombine.combine,
        "Type of combine used by NIFASTSKY", >> tmphead)
    printf (sformat,"NIFSTREJ",imcombine.reject,
        "Type of rejection used by NIFASTSKY", >> tmphead)
    if (imcombine.reject=="minmax") {
        printf (fformat,"NIFSTNLO",imcombine.nlow,
            "Low pixels rejected (minmax)", >> tmphead)
        printf (fformat,"NIFSTNHI",imcombine.nhigh,
            "High pixels rejected (minmax)", >> tmphead)
    }
    # Put input image names in header
    i=1
    scanfile = tmpfile1
    while (fscan(scanfile,img) != EOF) {
        printf (sformat,"NIFSIM"//str(i),img,
            "Input image combined with NIFASTSKY", >> tmphead)
        i=i+1
    }

    # put all the new stuff in the header
    # KL (Jan2011): this call to nhedit will work with both the fitsutil
    #       and the imutil version of nhedit.
    nhedit (l_outimage//"[0]",comfile=tmphead)

    scanfile=""
    #-----------------------------------------------------------------------

clean:
    if (status==0)
        printlog("NIFASTSKY exit status: good.", l_logfile,verbose=l_verbose)
    else if (status != 2)
        imdelete(l_outimage, verify-, >& "dev$null")
        
    printlog("-----------------------------------------------------------\
        -----------------", logfile=l_logfile, verbose=l_verbose)

    delete(tmpfile1, verify-, >& "dev$null")
    delete(scilist//","//dqlist//","//tmphead, verify-, >& "dev$null")
    imdelete(combsig//","//combdq, verify-, >& "dev$null")
    imdelete(dqsum, verify-, >& "dev$null")
    imdelete(tmpsci//"*.fits,"//tmpdq//"_*.pl", verify-, >& "dev$null")
    imdelete(tmpdqcomb, verify-, >& "dev$null")
    # return to default parameters for imcombine
    unlearn("imcombine")
    # restore the user's parameters for imcombine
    if (access("uparm$imhimcome.par.org"))
        rename("uparm$imhimcome.par.org","uparm$imhimcome.par",field="all")

end


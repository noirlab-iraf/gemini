# Copyright(c) 2002-2012 Association of Universities for Research in Astronomy, Inc.

procedure gnscombine(inimages,offsets,outimage)

# Task to combine gnsskysub reduced images
#
# Version Jun 03, 2003   RA,KG,IJ First version tested w/ v1.4, not ready for
#                        release the task has known bugs
# Version Aug 02, 2003   RA Known bugs fixed and task sent to Inger for
#                        approval

char    inimages    {prompt="Input GMOS images or list"}
char    offsets     {prompt="Text file with X,Y offsets"}
char    outimage    {prompt="Output image"}
char    outcheckim  {"",prompt="Name of output check (non-crrej) image"}
char    outmedsky   {"",prompt="Name of output median (for sky) image"}
char    sci_ext     {"SCI", prompt="Name of science extension"}
char    var_ext     {"VAR", prompt="Name of variance extension"}
char    dq_ext      {"DQ", prompt="Name of data quality extension"}
char    mdf_ext     {"MDF", prompt="Mask definition file extension name"}
bool    fl_vardq    {no, prompt="Create variance and data quality frames"}
char    logfile     {"",prompt="Logfile"}
bool    verbose     {yes,prompt="Verbose"}
int     status      {0, prompt="Exit status (0=good)"}
struct  *scanfile   {"",prompt="Internal use only"}

begin

    char    l_inimages, l_outimage, l_outcheckim, l_offsets, l_dir, l_outmedsky
    char    l_logfile, l_sci_ext, l_var_ext, l_dq_ext, l_mdf_ext, extn[3]
    char    appendextn, inextn, l_outputs[3], img, filelist, templist
    char    medframe, semifinal, tmpmedframe, inimg[200], temp[200], l_temp
    char    tempWithSky[200]
    int     maxfiles, shuffle, xbin, ybin, i, j, l_test
    int     ninp, noff, nsci, n_ext, num_extn, nextend
    real    x, y, dx[200], dy[200]
    bool    l_verbose, l_fl_vardq, mdfexist[200]
    string  suf
    struct  l_struct

    # Localize global variables
    l_inimages=inimages
    l_offsets=offsets 
    l_outimage=outimage
    l_outcheckim = outcheckim
    l_outmedsky = outmedsky
    l_verbose=verbose
    l_logfile = logfile
    l_sci_ext = sci_ext
    l_var_ext = var_ext
    l_dq_ext = dq_ext
    l_mdf_ext = mdf_ext
    l_fl_vardq = fl_vardq
    status = 0
    maxfiles = 200

    # If adding anymore variables, do not start them with 'out', CL doesn't
    # like it - there are too many variables with 'out' at the start - MS

    # Create any needed temporary file names (more defined later)
    filelist = mktemp("tmpfile")
    templist = ""
    medframe = ""
    semifinal = ""

    #-------------------------------------------------------------------
    # Test the logfile:

    if (l_logfile == "" || stridx(" ",l_logfile)>0) {
        l_logfile = gmos.logfile
        if (l_logfile == "" || stridx(" ",l_logfile)>0) {
            l_logfile = "gmos.log"
            printlog ("WARNING - GNSCOMBINE: Both gnscombine.logfile and \
                gmos.logfile fields are empty", l_logfile,verbose+)
            printlog("                    Using default file gmos.log",
                l_logfile,verbose+)
        }
    }

    #-------------------------------------------------------------------
    # Start logging

    date | scan(l_struct)
    printlog ("-----------------------------------------------------------\
        -----------------", l_logfile, l_verbose)
    printlog ("GNSCOMBINE -- "//l_struct, l_logfile, l_verbose)
    printlog ("", l_logfile, l_verbose)
    printlog ("Input list    = "//l_inimages, l_logfile, l_verbose)
    printlog ("Output list   = "//l_outimage, l_logfile, l_verbose)
    printlog ("", l_logfile, l_verbose)


    #-------------------------------------------------------------------
    # Last chance to do fiddly stuff before things start seriously happening!
    cache ("imgets")

    #-------------------------------------------------------------------
    # Load up input name list: @list, * and ?, comma separated
    if (l_inimages == "") {
        printlog ("ERROR - GNSCOMBINE: Input file not specified",
            l_logfile, verbose+)
        goto crash
    }

    # Test for @filelist
    if ((substr(l_inimages,1,1) == "@") \
        && !access(substr(l_inimages,2,strlen(l_inimages))) ) {

        printlog ("ERROR - GNSCOMBINE: Input list "//\
            substr(l_inimages,2,strlen(l_inimages))//" does not exist",
            l_logfile,verbose+)
        goto crash
    }

    # parse wildcard and comma-separated lists
    if (substr(l_inimages,1,1)=="@") {
        scanfile = substr(l_inimages,2,strlen(l_inimages))
        while (fscan(scanfile,l_temp) != EOF) {
            files (l_temp, >> filelist)
        }
    } else {
        if (stridx(",",l_inimages)==0)
            files (l_inimages, > filelist)
        else {
            l_test = 9999
            while (l_test!=0) {
                l_test = stridx(",",l_inimages)
                if (l_test>0)
                    files(substr(l_inimages,1,l_test-1), >> filelist)
                else
                    files(l_inimages, >> filelist)
                l_inimages = substr(l_inimages,l_test+1,strlen(l_inimages))
            }
        }
    }
    scanfile = ""

    #--------------------------------------------------------------------
    # Check for the presence of of the output files

    l_outputs[1] = l_outimage
    l_outputs[2] = l_outmedsky
    l_outputs[3] = l_outcheckim

    for (i = 1; i <= 3; i += 1) {
        if (l_outputs[i] != "" && (stridx(" ",l_outputs[i]) == 0)) {
            gimverify (l_outputs[i])
            if (gimverify.status == 1) {
                l_outputs[i] = gimverify.outname//".fits"
                gimverify (l_outputs[i])
                if (gimverify.status != 1) {
                    printlog ("ERROR - GNSCOMBINE: "//l_outputs[i]//\
                        " already exits", l_logfile, verbose+)
                    goto crash
                }
            } else {
                printlog ("ERROR - GNSCOMBINE: "//l_outputs[i]//\
                    " already exits", l_logfile, verbose+)
                goto crash
            }
        } else if (i == 1) {
            # All other files are optional outputs
            printlog ("ERROR - GNSCOMBINE: outimage not defined or is not \
                defined properly", l_logfile, verbose+)
        } else {
            l_outputs[i] = ""
        }
    }

    l_outimage = l_outputs[1]
    l_outmedsky = l_outputs[2]
    l_outcheckim = l_outputs[3]

    #--------------------------------------------------------------------
    # Define the input filename list

    ninp = 0
    scanfile = filelist
    while (fscan(scanfile, img) != EOF) {

        # split off directory path
        fparse (img, verbose-)

        img = fparse.root
        l_dir = fparse.directory

        ninp = ninp+1
        if (ninp > maxfiles) {
            printlog ("ERROR - GNSCOMBINE: Maximum number of input images \
                exceeded", l_logfile, verbose+)
            goto crash
        }
        inimg[ninp] = l_dir//img//".fits" # now has path if relevant

    } # end while
    scanfile = ""

    # At the end of all this inimg[ninp] now contains the input images incl.
    # the directory path
    #--------------------------------------------------------------------------
    # Define the offsets lists

    if (access (l_offsets) == no) {
        printlog ("ERROR - GNSCOMBINE: Offset file '"//l_offsets//"' \
            not found.", l_logfile, verbose+)
        goto crash
    }

    scanfile = l_offsets
    noff = 0
    while (fscan(scanfile,x,y) != EOF){
        noff = noff + 1
        dx[noff] = x
        dy[noff] = y
    }

    if (ninp > maxfiles) {
        printlog ("ERROR - GNSCOMBINE: Number of offsets does not correspond \
            to number of images", l_logfile, verbose+)
        goto crash
    }

    #--------------------------------------------------------------------------
    # Run a quick check to see if the input images have been gprepared
    # KNOWN BUG: l_logfile not defined, IJ
    i = 1
    while (i<=ninp) {
        imgets (inimg[i]//"[0]","NSCIEXT", >>& "dev$null")
        nsci = int(imgets.value)
        if (nsci==0) {
            printlog ("ERROR - GNSCOMBINE: Keyword NSCIEXT not found in \
                image.", l_logfile, verbose+)
            printlog ("ERROR - GNSCOMBINE: File "//inimg[i]//" is not \
                gprepared.", l_logfile, verbose+)
            printlog ("ERROR - GNSCOMBINE: Please run gprepare on all input \
                images.", l_logfile, verbose+)
            goto crash
        }
        i = i+1
    }

    #--------------------------------------------------------------------------
    # Making the output image. This is where all the action really is!
    # Step 0. Copy files for templates
    printlog ("GNSCOMBINE\n", l_logfile, verbose+)
    printlog ("Creating temporary files\n", l_logfile, verbose+)

    # Define extensions that could possible be worked on
    extn[1] = l_sci_ext
    extn[2] = l_var_ext
    extn[3] = l_dq_ext

    # Step 1. Shift images
    printlog ("Shifting images\n", l_logfile, verbose+)
    i = 0
    while (i<ninp) {
        i = i + 1

        temp[i] = mktemp("tmpshift")//".fits"
        imcopy (inimg[i]//"[0]", temp[i], verbose-)
        nextend = 0

        # Check for an MDF and attach it
        gemextn (inimg[i], process="expand", check="exists", \
            extver="", index="1-", extname=l_mdf_ext, \
            outfile="dev$null", \
            logfile=l_logfile, verbose=no)

        if (gemextn.count == 1) {
            mdfexist[i] = yes
            tcopy (inimg[i]//"["//l_mdf_ext//"]", \
                temp[i]//"["//l_mdf_ext//"]", verbose=no)
            nextend += 1
        } else {
            mdfexist[i] = no
        }

        num_extn = 1
        if (l_fl_vardq) {
            # Check for VAR and DQ planes and attach them
            gemextn (inimg[i], process="expand", check="exists,mef", \
                extver="1-", index="", extname=l_var_ext//","//l_dq_ext, \
                outfile="dev$null", omit="index", \
                logfile=l_logfile, verbose=l_verbose)

            if ((gemextn.count / 2) == nsci) {
                num_extn = 3
            }
        }

        # Loop over extension names
        for (j = 1; j <= num_extn; j += 1) {

            # Loop over extension versions
            for (n_ext=1; n_ext <= nsci; n_ext += 1) {

                inextn = "["//str(extn[j])//","//str(n_ext)//"]"
                appendextn = "["//str(extn[j])//","//str(n_ext)//",append]"

                # Only imshift if required - copy instead
                if (dx[i] != 0. || dy[i] != 0.) {

                    imshift (input=inimg[i]//inextn, \
                        output=temp[i]//appendextn, \
                        xshift=dx[i], yshift=dy[i], interp="nearest", \
                        shifts_file="", boundary_type="nearest", constant=0.)
                } else {
                    imcopy (input=inimg[i]//inextn, \ 
                        output=temp[i]//appendextn, verbose=no)
                }

                # Make sure DQ OBJECT keyword stay the same!
                if (extn[j] == l_dq_ext) {
                    keypar (inimg[i]//inextn, "OBJECT", silent+)
                    gemhedit (temp[i]//inextn, "OBJECT", keypar.value, \
                        "", delete-)
                }
                nextend += 1
            }
        }
        gemhedit (temp[i]//"[0]", "NEXTEND", nextend, "", delete-)
    }

    # Step 2. Get median sky
    printlog ("Getting median sky\n", l_logfile, verbose+)

    tmpmedframe = mktemp("tmpmedframe")//".fits"

    if (l_outmedsky != "") {
        medframe = l_outmedsky
    } else {
        medframe = mktemp("med")
    }

    # Have to write extension to variable otherwise the copy command will not
    # work
    suf = substr (medframe, strlen(medframe)-4, strlen(medframe))
    if (suf != ".fits") {
        medframe = medframe//".fits"
    }

    templist = mktemp("tmplist")
    i = 0
    while (i<ninp) {
        i = i+1
        print (temp[i], >>templist)
    }

    if (ninp == 1){
        printlog ("WARNING - GNSCOMBINE: Only one image to combine, \
            skipping gemcombine step", l_logfile, verbose+)
        printlog ("                      Input image copied to median sky \
            image", l_logfile, verbose+)
        copy (temp[1], tmpmedframe, verbose-, >& "dev$null")
    } else {
        gemcombine (input="@"//templist, output=tmpmedframe, title="",
            combine="median", reject="none", offsets="none",
            masktype="goodvalue", maskvalue=0., scale="none", zero="none",
            weight="none", statsec="[*,*]", expname="EXPTIME",
            lthreshold=INDEF, hthreshold=INDEF, nlow=1, nhigh=1, nkeep=1,
            mclip=yes, lsigma=3., hsigma=3., key_ron="RDNOISE",
            key_gain="GAIN", ron=0., gain=1., snoise="0.0", sigscale=0.1,
            pclip=-0.5, grow=0., bpmfile="", nrejfile="", sci_ext=l_sci_ext,
            var_ext=l_var_ext, dq_ext=l_dq_ext, fl_vardq=l_fl_vardq, \
            logfile=l_logfile,
            fl_dqprop=no, verbose=l_verbose)

        if (gemcombine.status != 0) {
            printlog ("ERROR - GNSCOMBINE: GEMCOMBINE returned a non-zero \
                status.", l_logfile, verbose+)
            goto crash
        }
    }

    printlog ("", l_logfile, verbose+)

    # Shift the sky frame - this just makes it look pretty with gnsskysub
    # doing A-B
    # The combined image may still have blotches in the +sky or -sky areas
    hselect (tmpmedframe//"[0]", "NODPIX", yes) | scan(shuffle)
    printlog ("Making shifted sky for shuffle = "//shuffle, l_logfile, \
        verbose+)

#    PyRAF doesn't like hselect/translit, prefer keypar
#    hselect (medframe//"["//l_sci_ext//",1]", "CCDSUM", yes) | \
#        translit ("STDIN", '"', "", delete+, collapse-) | scan (xbin, ybin)
    keypar(tmpmedframe//"["//l_sci_ext//",1]", "CCDSUM")
    print(keypar.value) | scan(xbin,ybin)

    # Correct shuffle for bining in Y
    shuffle = shuffle/ybin

    # Shuffle median image
    imcopy (tmpmedframe//"[0]", medframe, verbose-)

    # Keep track of number of extensions
    nextend = 0

    # Check for an MDF and attach it
    gemextn (tmpmedframe, process="expand", check="exists", \
        extver="", index="1-", extname=l_mdf_ext, \
        outfile="dev$null", \
        logfile=l_logfile, verbose=no)

    if (gemextn.count == 1) {
        tcopy (tmpmedframe//"["//l_mdf_ext//"]", \
            medframe//"["//l_mdf_ext//"]", verbose=no)
        nextend += 1
    }

    # Loop over extenversions
    for (n_ext=1; n_ext<=nsci; n_ext+=1) {
        # Shift SCI extensions
        imshift (input=tmpmedframe//"["//l_sci_ext//","//str(n_ext)//"]",
            output=medframe//"["//l_sci_ext//","//str(n_ext)//",append]", \
            xshift=0, yshift=shuffle, interp="nearest", shifts_file="",
            boundary_type="nearest", constant=0.)

        nextend += 1

        # Shift var and dq extensions too
        if (imaccess(tmpmedframe//"["//l_var_ext//","//str(n_ext)//"]")) {
            imshift (input=tmpmedframe//"["//l_var_ext//","//str(n_ext)//"]",
                output=medframe//"["//l_var_ext//","//str(n_ext)//",append]",
                xshift=0, yshift=shuffle, interp="nearest", shifts_file="",
                boundary_type="nearest", constant=0.)

            imshift (input=tmpmedframe//"["//l_dq_ext//","//str(n_ext)//"]",
                output=medframe//"["//l_dq_ext//","//str(n_ext)//",append]",
                xshift=0, yshift=shuffle, interp="nearest", shifts_file="",
                boundary_type="nearest", constant=0.)
            nextend += 2
        }
    }

    gemhedit (medframe//"[0]", "NEXTEND", nextend, "", delete-)

    imdelete (tmpmedframe, verify-, >& "dev$null")

    # Step 3. Slice and dice
    printlog ("Calling gnsskysub", l_logfile, verbose+)

    gnsskysub (inimages="@"//templist, outimages="", outpref="n", fl_fixnc+,
        fl_vardq=l_fl_vardq, logfile=l_logfile, verbose=l_verbose)

    if (gnsskysub.status !=0) {
        printlog ("ERROR - GNSCOMBINE: GNSSKYSUB returned a non-zero status", \
            l_logfile, verbose+)
        got crash
    }

    printlog ("", l_logfile, verbose+)

    # Step 4. Add back the constant sky (so errors are right)
    printlog ("Adding back constant sky so errors are computed correctly \
        in the next step", l_logfile, verbose+)

    delete (templist, verify-, >& "dev$null")

    i = 1
    while (i<=ninp) {
        tempWithSky[i] = mktemp("tmpwithsky")//".fits"

        printlog ("Adding median sky to "//"n"//temp[i]//\
            " -> "//tempWithSky[i], l_logfile, verbose+)

        # Don't use gemarith as it'll add noise to the variance plane
        # Not sure if this is entirely correct - MS
        imcopy ("n"//temp[i]//"[0]", tempWithSky[i], verbose-)

        if (mdfexist[i]) {
            tcopy ("n"//temp[i]//"["//l_mdf_ext//"]", tempWithSky[i], \
                verbose=no)
        }

        for (j = 1; j <= nsci; j += 1) {
            imarith ("n"//temp[i]//"["//l_sci_ext//","//j//"]", "+", \
                medframe//"["//l_sci_ext//","//j//"]",
                tempWithSky[i]//"["//l_sci_ext//","//j//",append]", verbose-)
            if (imaccess("n"//temp[i]//"["//l_var_ext//","//j//"]")) {
                imcopy ("n"//temp[i]//"["//l_var_ext//","//j//"]", \
                    tempWithSky[i]//"["//l_var_ext//","//j//",append]", \
                    verbose-)
                imcopy ("n"//temp[i]//"["//l_dq_ext//","//j//"]", \
                    tempWithSky[i]//"["//l_dq_ext//","//j//",append]", \
                    verbose-)
            }
        }

        print (tempWithSky[i], >>templist)
        i = i+1
    }
    printlog ("", l_logfile, verbose+)

    # Step 5. Now do the imcombine
    printlog ("Running gemcombine", l_logfile, verbose+)
    printlog ("", l_logfile, verbose+)
    semifinal = mktemp("tmpsemi")//".fits"

    if (ninp == 1){
        printlog ("WARNING - GNSCOMBINE: Only one image to combine, \
            skipping gemcombine step", l_logfile, verbose+)
        printlog ("                      Input image copied to median sky \
            image", l_logfile,verbose+)
        copy (tempWithSky[1], semifinal, verbose-, >& "dev$null")
    } else {
        gemcombine (input="@"//templist, output=semifinal, title="",
            combine="average", reject="ccdclip", offsets="none",
            masktype="goodvalue", maskvalue=0., scale="none", zero="none",
            weight="none", statsec="[*,*]", expname="EXPTIME",
            lthreshold=INDEF, hthreshold=INDEF, nlow=1, nhigh=1, nkeep=1,
            mclip=yes, lsigma=10, hsigma=10, key_ron="RDNOISE",
            key_gain="GAIN", ron=0., gain=1., snoise="0.0", sigscale=0.1,
            pclip=-0.5, grow=2., bpmfile="", nrejfile="", sci_ext=l_sci_ext,
            var_ext=l_var_ext, dq_ext=l_dq_ext, fl_vardq=l_fl_vardq, \
            logfile=l_logfile,
            fl_dqprop=no, verbose=l_verbose)

        if (gemcombine.status != 0) {
            printlog ("ERROR - GNSCOMBINE: GEMCOMBINE returned a non-zero \
                status.", l_logfile, verbose+)
            goto crash
        }
    }
    printlog ("", l_logfile, verbose+)

    # Step 6. Now remove the constant sky again
    printlog ("Subtracting median sky", l_logfile, verbose+)

    # Keep track of number of extensions
    nextend = 0

    # Create final output image
    imcopy (semifinal//"[0]", l_outimage, verbose-)

    # Check for an MDF and attach it
    gemextn (semifinal, process="expand", check="exists", \
        extver="", index="1-", extname=l_mdf_ext, \
        outfile="dev$null", \
        logfile=l_logfile, verbose=no)

    if (gemextn.count == 1) {
        tcopy (semifinal//"["//l_mdf_ext//"]", \
            l_outimage//"["//l_mdf_ext//"]", verbose=no)
        nextend += 1
    }

    # Don't use gemarith as it'll add noise to the variance planes
    # Not sure if this is entirely correct - MS
    # Subract medframe from semifinal image
    for (j = 1; j <= nsci; j += 1) {
        # Only do this for SCI extensions
        imarith (semifinal//"["//l_sci_ext//","//j//"]", "-", \
            medframe//"["//l_sci_ext//","//j//"]",
            l_outimage//"["//l_sci_ext//","//j//",append]", verbose-)

        nextend += 1
        # Copy VAR and DQ if present
        if (imaccess(semifinal//"["//l_var_ext//","//j//"]")) {
            imcopy (semifinal//"["//l_var_ext//","//j//"]", \
                l_outimage//"["//l_var_ext//","//j//",append]", verbose-)
            imcopy (semifinal//"["//l_dq_ext//","//j//"]", \
                l_outimage//"["//l_dq_ext//","//j//",append]", verbose-)
            nextend += 2
        }
    }

    gemhedit (l_outimage//"[0]", "NEXTEND", nextend, "", delete-)

    gemdate ()
    gemhedit (l_outimage//"[0]", "GNSCOMB", gemdate.outdate,
        "UT Time stamp for GNSCOMBINE", delete-)
    gemhedit (l_outimage//"[0]", "GEM-TLM", gemdate.outdate,
        "UT Last modification with GEMINI", delete-)
    printlog ("Created cr-rej combined image "//l_outimage, l_logfile, \
        verbose+)

    if (l_outmedsky != "") {
        gemhedit (l_outmedsky//"[0]", "GNSCOMB", gemdate.outdate,
            "UT Time stamp for GNSCOMBINE", delete-)
        gemhedit (l_outmedsky//"[0]", "GEM-TLM", gemdate.outdate,
            "UT Last modification with GEMINI", delete-)
        printlog ("Created median frame for sky "//l_outmedsky, l_logfile, \
            verbose+)
    }

    if (l_outcheckim == "") goto clean

    # Step7 Make check image (no rejection)
    printlog ("Making check [non-crrej] image "//l_outcheckim, l_logfile, \
        verbose+)
    delete (templist, verify-, >>& "dev$null")
    i = 1
    while (i<=ninp) {
        print ("n"//temp[i], >>templist)
        i = i+1
    }
    printlog ("", l_logfile, verbose+)

    if (ninp == 1){
        printlog ("WARNING - GNSCOMBINE: Only one image to combine, \
            skipping gemcombine step", l_logfile, verbose+)
        printlog ("                      Input image copied to median sky \
            image", l_logfile,verbose+)

        # check .fits
        suf = substr (l_outcheckim, strlen(l_outcheckim)-4,
            strlen(l_outcheckim))
        if (suf != ".fits")
            l_outcheckim = l_outcheckim//".fits"

        copy ("n"//temp[1], l_outcheckim, verbose-, >& "dev$null")
    } else {

        gemcombine (input="@"//templist, output=l_outcheckim, title="",
            combine="average", reject="none", offsets="none",
            masktype="goodvalue", maskvalue=0., scale="none", zero="none",
            weight="none", statsec="[*,*]", expname="EXPTIME",
            lthreshold=INDEF, hthreshold=INDEF, nlow=1, nhigh=1, nkeep=1,
            mclip=yes, lsigma=3., hsigma=3., key_ron="RDNOISE",
            key_gain="GAIN", ron=0., gain=1., snoise="0.0", sigscale=0.1,
            pclip=-0.5, grow=0., bpmfile="", nrejfile="", sci_ext=l_sci_ext,
            var_ext=l_var_ext, dq_ext=l_dq_ext, fl_vardq=l_fl_vardq, \
            logfile=l_logfile,
            fl_dqprop=no, verbose=l_verbose)

        if (gemcombine.status != 0) {
            printlog ("ERROR - GNSCOMBINE: GEMCOMBINE returned a non-zero \
                status.", l_logfile, verbose+)
            goto crash
        }
    }
    printlog ("", l_logfile, verbose+)

    gemdate ()
    gemhedit (l_outcheckim//"[0]", "GNSCOMB", gemdate.outdate,
        "UT Time stamp for GNSCOMBINE", delete-)
    gemhedit (l_outcheckim//"[0]", "GEM-TLM", gemdate.outdate,
        "UT Last modification with GEMINI", delete-)
    printlog ("Created non-crrej combined image "//l_outcheckim, l_logfile, \
        verbose+)

    goto clean

crash:
    # Exit with error subroutine
    status = 1
    printlog ("ERROR - GNSCOMBINE: Exit with error", l_logfile, verbose+)
    goto clean

clean:
    # clean up
    delete (filelist, verify-, >& "dev$null")
    delete (templist, verify-, >& "dev$null")
    if (l_outmedsky == "")
        imdelete (medframe,verify-, >& "dev$null")
    imdelete (semifinal, verify-, >& "dev$null")
    scanfile = ""
    i = 1
    while (i<=ninp) {
        imdelete (temp[i], verify-, >& "dev$null")
        imdelete ("n"//temp[i], verify-, >& "dev$null")
        imdelete (tempWithSky[i], verify-, >& "dev$null")
        i = i+1
    }

    if (status == 0) {
        printlog ("GNSCOMBINE: Exit status GOOD", l_logfile, verbose+)
    } else {
        printlog ("GNSCOMBINE: Exit status ERROR", l_logfile, verbose+)
    }

end

# Copyright(c) 2000-2011 Association of Universities for Research in Astronomy, Inc.

procedure nireduce(inimages)

# Basic reductions of NIRI images.
# Requires a flat, bad pixel mask, sky and/or dark image to subtract.  Adding
# a constant sky back on can be handled explicitly or using fl_autosky.
#
# A log of the reductions is written in the logfile
# The header of the reduced image will contain information about
# the reduction as a series of header entries.
#
# Variance handling:
#   var = (df/da)^2 * siga^2 + (df/db)^2 * sigb^2
#   f=a-k*b   var(f) = var(a) + k*k*var(b)
#   f=a/b     var(f) = [sci(a)^2/sci(b)^4]*var(b) + [1/sci(b)^2]*var(a)
#   f=(a-b)/c var(f) = [(sci(a)^2+sci(b)^2-2sci(a)sci(b))/sci(c)^4]*var(c) +
#                      (1/sci(c)^2)*[var(a)+var(b)]
#
# Data Quality Array handling:
#   DQ = dq(a) OR dq(b)
#
# Read noise header keyword in output image is the quadrature sum of the
# values in the input images and the sky or dark images.  No adjustment
# has been made for division by the flat field, assuming the read noise in
# the flat field is negligible.  The saturation header keyword in the 
# the output image is the input value minus the median value in the sky
# or dark image + the constant sky value added back on after.
#
# Version  June 19, 2001  JJ, beta-release
#          Oct 15, 2001 JJ, v1.2 release
#          16 Nov 2001  JJ, remove excess filter and flatimage parameters
#                           add scaling of skyconst if fl_scalesky+
#                           improved filter checking for sky, flat
#                           improved logic checking on sky subtraction
#          Feb 28, 2002 JJ, v1.3 release
#          11 Apr 2002  JJ, fixed bug in how nonlinear value is computed
#                           fixed bug for repeated file names in a list
#          24 May 2002  JJ, fixed ERROR messages and write direct to MEF
#          Sept 10, 2002 IJ parameter encoding
#          Sept 20, 2002 JJ v1.4 release
#          Aug 18, 2003  KL IRAF2.12 - new parameters
#                             hedit: addonly
#                             imstat: nclip,lsigma,usigma,cache

char    inimages    {prompt="Input NIRI image(s)"}                          # OLDP-1-input-primary-single-prefix=r
char    outimages   {"",prompt="Output image(s)"}                           # OLDP-1-output
char    outprefix   {"r",prompt="Prefix for output image(s)"}               # OLDP-4
char    logfile     {"",prompt="Logfile"}                                   # OLDP-1
bool    fl_sky      {yes,prompt="Do sky subtraction?"}                      # OLDP-2
char    skyimage    {"",prompt="Sky image to subtract"}                     # OLDP-1-input-fl_sky-required
real    skylevel    {0.0,prompt="Constant sky level to add"}                # OLDP-3
bool    fl_autosky  {yes,prompt="Add median of the sky frame as a constant?"}   # OLDP-3
bool    fl_scalesky {yes,prompt="Scale the sky before subtracting?"}        # OLDP-3
bool    fl_dark     {no,prompt="Do explicit dark subtraction?"}             # OLDP-3
char    darkimage   {"",prompt="Dark current image to subtract"}            # OLDP-1-input-fl_dark-required
bool    fl_flat     {yes,prompt="Do flat-fielding?"}                        # OLDP-2
char    flatimage   {"",prompt="Flat field image to divide"}                # OLDP-1-input-fl_flat-required
char    statsec     {"[100:924,100:924]",prompt="Statistics section"}       # OLDP-2
char    sci_ext     {"SCI",prompt="Name or number of science extension"}    # OLDP-3
char    var_ext     {"VAR",prompt="Name or number of variance extension"}   # OLDP-3
char    dq_ext      {"DQ",prompt="Name or number of data quality extension"}    # OLDP-3
char    key_filter  {"FILTER",prompt="Keyword for filter id"}               # OLDP-3
char    key_ron     {"RDNOISE",prompt="Header keyword for read noise (e-)"} # OLDP-3
char    key_sat     {"SATURATI",prompt="Keyword for saturation (ADU)"}      # OLDP-3
char    key_nonlinear   {"NONLINEA",prompt="Header keyword for non-linear regime (ADU)"} # OLDP-3
bool    fl_vardq    {yes,prompt="Create variance and data quality frames?"} # OLDP-3
bool    verbose     {yes,prompt="Verbose"}                                  # OLDP-4
int     status      {0,prompt="Exit status (0=good)"}                       # OLDP-4
struct  *scanfile   {prompt="Internal use"}                                 # OLDP-4

begin

    char    l_inimages, l_outimages, l_skyimage, l_darkimage
    char    l_flatimage, l_filter
    char    l_expression, l_varexpression, l_prefix, l_logfile, l_temp
    char    in[1000], out[1000], l_statsec, l_flatfilter, l_flatim, l_skyfilter
    char    l_key_filter, l_sci_ext, l_var_ext, l_dq_ext
    char    tmpdq, tmpfile, keyfound
    char    tmpdqsci, tmpdqsky, tmpdqflat, tmpdqdark 
    char    l_dqlist, l_varlist, l_scilist
    int     i, nimages, noutimages, maxfiles, ndq
    bool    l_fl_sky, l_fl_dark, l_fl_flat, l_verbose
    bool    l_fl_skytemp, l_fl_flattemp
    bool    l_fl_first, l_fl_darktemp, l_fl_autosky, l_fl_vardq, bad
    bool    l_fl_indvardq[1000], l_fl_scalesky
    real    l_skylevel, l_midsky, l_middark, l_mid, l_sig, l_ron, l_nonlinear
    real    l_skyscale, l_middata
    char    l_key_ron, l_key_sat, l_key_nonlinear
    int     l_sat, nbad, nnew
    char    ima, imb, imc, imd, ime, imf, img, imh

    struct l_struct

    status = 0
    nimages = 0
    maxfiles = 1000
    tmpfile = mktemp("tmpin")
    tmpdq = mktemp("tmpdq")

    cache ("imgets", "gemdate")

    # set the local variables
    l_inimages=inimages ; l_outimages=outimages ; l_skyimage=skyimage
    l_skylevel=skylevel ; l_darkimage=darkimage
    l_fl_sky=fl_sky ; l_fl_dark=fl_dark ; l_fl_flat=fl_flat
    l_verbose=verbose ; l_prefix=outprefix
    l_logfile=logfile ; l_fl_autosky=fl_autosky
    l_key_filter=key_filter ; l_flatimage=flatimage
    l_sci_ext=sci_ext ; l_var_ext=var_ext ; l_dq_ext=dq_ext
    l_key_ron=key_ron ; l_key_sat=key_sat
    l_fl_vardq=fl_vardq ; l_statsec=statsec ; l_key_nonlinear=key_nonlinear
    l_fl_scalesky=fl_scalesky

    #------------------------------------------------------------------------
    # Check for package log file or user-defined log file
    cache ("niri")
    print (l_logfile) | scan (l_logfile)
    if (l_logfile=="" || l_logfile==" ") {
        l_logfile = niri.logfile
        print (l_logfile) | scan (l_logfile)
        if (l_logfile=="" || l_logfile==" ") {
            l_logfile = "niri.log"
            printlog ("WARNING - NIREDUCE: Both nireduce.logfile and \
                niri.logfile are empty.", l_logfile, verbose+)
            printlog ("                    Using default file niri.log.",
                l_logfile, verbose+)
        }
    }
    # Open log file
    date | scan(l_struct)
    printlog ("-------------------------------------------------------------\
        ---------------", l_logfile, l_verbose)
    printlog ("NIREDUCE -- "//l_struct, l_logfile, l_verbose)
    printlog (" ", l_logfile, l_verbose)

    #----------------------------------------------------------------------
    # Check for consistent sky/sky level logic
    if (l_fl_autosky) {
        if (l_skylevel != 0.0) {
            printlog ("ERROR - NIREDUCE: You have specified both a sky \
                constant to add to the", l_logfile, verbose+)
            printlog ("                  final image AND set the fl_autosky \
                flag to determine the", l_logfile, verbose+)
            printlog ("                  sky constant from the sky image.  \
                These options are", l_logfile, verbose+)
            printlog ("                  incompatible.", l_logfile, verbose+)
            status = 1
            goto clean
        } else if (!l_fl_sky) {
            printlog ("WARNING - NIREDUCE: You have set the fl_autosky flag \
                to determine the sky", l_logfile, verbose+)
            printlog ("                    constant from the sky image, but \
                have the fl_sky sky", l_logfile, verbose+)
            printlog ("                    subtraction flag off.  Setting \
                fl_autosky=no.", l_logfile, verbose+)
            l_fl_autosky = no
        } 
    } else {
        if ((l_skylevel == 0.0) && (l_fl_sky)) {
            printlog ("WARNING - NIREDUCE: You have set the sky constant to \
                0.0 and have the", l_logfile, verbose+)
            printlog ("                    fl_autosky flag off.  No constant \
                will be added", l_logfile, verbose+)
            printlog ("                    after the sky image is \
                subtracted.", l_logfile, verbose+)
        } else if ((l_skylevel != 0.0) && (!l_fl_sky)) {
            printlog ("WARNING - NIREDUCE: You have specified a sky constant \
                to add on, but", l_logfile, verbose+)
            printlog ("                    the sky image subtraction flag \
                fl_sky is off.  The", l_logfile, verbose+)
            printlog ("                    sky constant will be reset to 0.", \
                l_logfile, verbose+)
            l_skylevel = 0.
        }
    }

    # Check to see if both dark subtraction and sky subtraction are indicated
    if (l_fl_sky && l_fl_dark) { 
        printlog ("WARNING - NIREDUCE: Both a sky frame and a dark frame will \
            be subtracted.", l_logfile, verbose+)
        printlog ("                    Since most sky frames contain dark \
            current, this is ", l_logfile, verbose+)
        printlog ("                    probably a bad idea!  Continuing \
            anyway.", l_logfile, verbose+)
    }

    # Check to see if sci_ext is defined
    print (l_sci_ext) | scan (l_sci_ext)
    if (l_sci_ext=="" || l_sci_ext==" ") {
        printlog ("ERROR - NIREDUCE:  Extension sci_ext is undefined.",
            l_logfile, verbose+)
        status = 1
        goto clean
    }

    #-----------------------------------------------------------------------
    # Load up arrays of input name lists

    if (substr(l_inimages,1,1)=="@") {
        l_temp = substr(l_inimages,2,strlen(l_inimages))
        if (!access(l_temp)) {
            printlog ("ERROR - NIREDUCE:  Input file "//l_temp//" not found.",
                l_logfile, verbose+)
            status = 1
            goto clean
        }
    }

    nimages = 0
    nbad = 0
    files (l_inimages, sort-, > tmpfile) 
    scanfile = tmpfile

    while(fscan(scanfile,l_temp) != EOF) {
        if (l_fl_vardq) l_fl_indvardq[nimages+1] = yes
        gimverify (l_temp)
        if (gimverify.status==1) {
            printlog ("ERROR - NIREDUCE: File "//l_temp//" not found.",
                l_logfile, verbose+)
            nbad+=1
        } else if(gimverify.status>1) {
            printlog ("ERROR - NIREDUCE: File "//l_temp//" not a MEF FITS \
                image.", l_logfile, verbose+)
            nbad+=1
        } else {
            keyfound = ""
            hselect (l_temp//"[0]", "*PREPAR*", yes) | scan (keyfound)
            if (keyfound == "") {
                printlog ("ERROR - NIREDUCE: Image "//l_temp//" not \
                    *PREPAREd.", l_logfile, verbose+)
                nbad+=1
            } else {
                if (l_fl_vardq && !imaccess(l_temp//"["//l_var_ext//"]") ) {
                    printlog ("WARNING - NIREDUCE: Cannot compute VAR planes \
                        because", l_logfile, verbose+)
                    printlog ("                    input image "//l_temp//" \
                        does not have a VAR plane.", l_logfile, verbose+)
                    printlog ("                    Resetting fl_vardq=no for \
                        this file", l_logfile, verbose+)
                    l_fl_indvardq[nimages+1] = no
                }
                if (l_fl_vardq && !imaccess(l_temp//"["//l_dq_ext//"]") ) {
                    printlog ("WARNING - NIREDUCE: Cannot compute DQ planes \
                        because", l_logfile, verbose+)
                    printlog ("                    input image "//l_temp//" \
                        does not have a DQ plane.", l_logfile, verbose+)
                    printlog ("                    Resetting fl_vardq=no for \
                        this file", l_logfile, verbose+)
                    l_fl_indvardq[nimages+1] = no
                }
            }
            # strip .fits if present
            if (substr(l_temp,strlen(l_temp)-4,strlen(l_temp)) == ".fits")
                l_temp=substr(l_temp,1,(strlen(l_temp)-5))
            nimages = nimages+1
            in[nimages] = l_temp 
            if (nimages>1) {
                for (j=1; j<nimages; j+=1) {
                    if (in[nimages]==in[j]) {
                        printlog ("WARNING - NIREDUCE: Input image name \
                            "//in[nimages]//" repeated.", l_logfile, verbose+)
                        printlog ("                    Not including it \
                            again.", l_logfile, verbose+)
                        nimages = nimages-1
                    }
                }
            }
        }
    } # end while

    # Exit if problems found
    if (nbad > 0) {
        printlog ("ERROR - NIREDUCE: "//nbad//" image(s) either do not exist, \
            are not MEF files, or", l_logfile, verbose+)
        printlog ("                   have not been run through *PREPARE.", \
            l_logfile, verbose+)
        status = 1
        goto clean
    }
    if (nimages > maxfiles) {
        printlog ("ERROR - NIREDUCE: Maximum number of input images \
            exceeded ("//str(maxfiles)//")", l_logfile, verbose+ )
        status = 1
        goto clean
    }
    if (nimages == 0) {
        printlog ("ERROR - NIREDUCE: No valid input images.",
            l_logfile, verbose+ )
        status = 1
        goto clean
    }

    printlog ("Processing "//nimages//" file(s).", l_logfile, l_verbose)
    scanfile = ""
    delete (tmpfile, verify-, >& "dev$null")

    #--------------------------------------------------------------------------
    # Now for the output images
    # outimages could contain legal * if it is of a form like %st%stX%*.imh

    noutimages = 0
    if (l_outimages!="" && l_outimages!=" ") {
        if (substr(l_outimages,1,1)=="@") 
            scanfile = substr(l_outimages,2,strlen(l_outimages))
        else if (stridx("*",l_outimages)>0)  {
            files (l_outimages,sort-) | \
                match (".hhd", stop+, print-, metach-, > tmpfile)
            scanfile = tmpfile
        } else {
            files (l_outimages, sort-, > tmpfile)
            scanfile = tmpfile
        }

        while (fscan(scanfile,l_temp) != EOF) {
            noutimages = noutimages+1
            if (noutimages > maxfiles) {
                printlog ("ERROR - NIREDUCE: Maximum number of output images \
                    exceeded ("//str(maxfiles)//")", l_logfile, verbose+)
                status = 1
                goto clean
            }
            out[noutimages] = l_temp 
            if (imaccess(out[noutimages])) {
                printlog ("ERROR - NIREDUCE: Output image "//\
                    out[noutimages]//" already exists", l_logfile, verbose+)
                status = 1
            }
        }
        if (status != 0) goto clean

        # if there are too many or too few output images exit with error
        if (nimages != noutimages) {
            printlog ("ERROR - NIREDUCE: Number of input and output images \
                are unequal.", l_logfile, verbose+)
            status = 1
            goto clean
        }
    } else { #use prefix instead
        print (l_prefix) | scan (l_prefix)
        if (l_prefix=="" || l_prefix==" ") {
            printlog ("ERROR - NIREDUCE: Neither output image name nor \
                output prefix is defined.", l_logfile, verbose+)
            status = 1
            goto clean
        }
        i = 1
        nnew = 1
        while (i<=nimages) {
            out[i] = l_prefix//in[i]
            for (j=1; j<i; j+=1) {
                if (out[i]==out[j]) {
                    printlog ("WARNING - NIREDUCE: Output image name \
                        "//out[i]//" repeated.", l_logfile, verbose+)
                    printlog ("                    Appending _"//nnew//" to \
                        the output file name.", l_logfile, verbose+)
                    out[i] = out[i]//"_"//nnew
                    nnew = nnew+1
                }
            }

            if (imaccess(out[i])) {
                printlog ("ERROR - NIREDUCE: Output image "//out[i]//\
                    " already exists.", l_logfile, verbose+)
                status = 1
            }
            i = i+1
        }
        if (status != 0) goto clean
    }
    scanfile = ""
    delete (tmpfile, verify-, >& "dev$null")

    #-------------------------------------------------------------------------
    # Check for existence of sky, flat, or dark images, if needed, and 
    # check levels

    # SKY
    if (l_fl_sky) {
        if (!imaccess(l_skyimage) && (l_skyimage != "") \
            && stridx(" ",l_skyimage)<=0) {
            
            printlog ("ERROR - NIREDUCE: Sky image "//l_skyimage//" not \
                found.", l_logfile, verbose+)
            status = 1
            goto clean
        } else if (l_skyimage=="" || stridx(" ",l_skyimage)>0 ) {
            printlog ("ERROR - NIREDUCE: Sky image is either an empty string \
                or contains spaces.", l_logfile, verbose+)
            status = 1
            goto clean
        }
        if (l_fl_vardq) {
            if ((!imaccess(l_skyimage//"["//l_var_ext//"]")) \
                || (!imaccess(l_skyimage//"["//l_dq_ext//"]"))) {
                
                printlog ("WARNING - NIREDUCE: Sky image does not contain \
                    both VAR and DQ planes.", l_logfile, verbose+)
                printlog ("                    Setting fl_vardq=no and \
                    proceeding.", l_logfile, verbose+)
                l_fl_vardq = no
            }
        }
        # remove .fits if present
        if (substr(l_skyimage,strlen(l_skyimage)-4,strlen(l_skyimage)) == \
            ".fits")
            l_skyimage=substr(l_skyimage,1,(strlen(l_skyimage)-5))

        imgets (l_skyimage//"[0]", l_key_filter, >& "dev$null")
        if (imgets.value=="0") {
            printlog ("WARNING - NIREDUCE: Cannot read filter from sky image.",
                l_logfile, verbose+)
            l_skyfilter = "none"
        } else {
            l_skyfilter = imgets.value
        }   

        # get median sky level
        imstat (l_skyimage//"["//l_sci_ext//"]"//l_statsec,
            fields="midpt,stddev", lower=INDEF, upper=INDEF, nclip=0,
            lsigma=INDEF, usigma=INDEF, binwidth=0.01, format-, 
            cache-) | scan (l_mid, l_sig)
        imstat (l_skyimage//"["//l_sci_ext//"]"//l_statsec, fields="midpt", \
            upper=(l_mid+3*l_sig), lower=(l_mid-3*l_sig), nclip=0, \
            lsigma=INDEF, usigma=INDEF, binwidth=0.01, format-, cache-) | \
            scan (l_mid)
        if (l_mid==INDEF) {
            printlog ("ERROR - NIREDUCE: Statistics failed, possibly due to \
                a bad statsec.", l_logfile, verbose+)
            status = 1
            goto clean
        }
        l_midsky = l_mid
    }

    # FLAT
    if (l_fl_flat) {
        if (!imaccess(l_flatimage) && (l_flatimage != "") \
            && stridx(" ",l_flatimage)<=0) {
            
            printlog ("ERROR - NIREDUCE: Flat field image "//l_flatimage//\
                " not found.", l_logfile, verbose+)
            nbad+=1
            status = 1
            goto clean
        } else if (l_flatimage=="" || stridx(" ",l_flatimage)>0 ) {
            printlog ("ERROR - NIREDUCE: Flat field image is either an empty \
                string or", l_logfile, verbose+)
            printlog ("                  contains spaces.",l_logfile,verbose+)
            status = 1
            goto clean
        }
        if (l_fl_vardq) {
            if ((!imaccess(l_flatimage//"["//l_var_ext//"]")) \
                || (!imaccess(l_flatimage//"["//l_dq_ext//"]"))) {
                
                printlog ("WARNING - NIREDUCE: Flat field image does not \
                    contain both VAR and DQ planes.", l_logfile, verbose+)
                printlog ("                    Setting fl_vardq=no and \
                    proceeding.", l_logfile, verbose+)
                l_fl_vardq = no
            }
        }
        # remove .fits if present
        if (substr(l_flatimage,strlen(l_flatimage)-4,strlen(l_flatimage)) == \
            ".fits")
            l_flatimage=substr(l_flatimage,1,(strlen(l_flatimage)-5))

        imgets (l_flatimage//"[0]", l_key_filter, >& "dev$null")
        if (imgets.value=="0") {
            printlog ("WARNING - NIREDUCE: Cannot read filter from flat \
                field image.", l_logfile, verbose+)
            l_flatfilter = "none"
        } else {
            l_flatfilter = imgets.value
        }   
    }

    # DARK
    if (l_fl_dark) {
        if (!imaccess(l_darkimage)&& (l_darkimage != "") \
            && stridx(" ",l_darkimage)<=0) {
            
            printlog ("ERROR - NIREDUCE: Dark current image "//l_darkimage//\
                " not found.", l_logfile, verbose+)
            nbad+=1
            status = 1
            goto clean
        } else if (l_darkimage=="" || stridx(" ",l_darkimage)>0 ) {
            printlog ("ERROR - NIREDUCE: Dark current image is either an \
                empty string or", l_logfile, verbose+)
            printlog ("                   contains spaces.", l_logfile, \
                verbose+)
            status = 1
            goto clean
        }
        if (l_fl_vardq) {
            if ((!imaccess(l_darkimage//"["//l_var_ext//"]")) \
                || (!imaccess(l_darkimage//"["//l_dq_ext//"]"))) {
                
                printlog ("WARNING - NIREDUCE: Dark image does not contain \
                    both VAR and DQ planes.", l_logfile, verbose+)
                printlog ("                    Setting fl_vardq=no and \
                    proceeding.", l_logfile, verbose+)
                l_fl_vardq = no
            }
        }    
        # remove .fits if present
        if (substr(l_darkimage,strlen(l_darkimage)-4,strlen(l_darkimage)) == \
            ".fits")
            l_darkimage = substr(l_darkimage,1,(strlen(l_darkimage)-5))

        imgets (l_darkimage//"[0]", l_key_filter, >& "dev$null")
        if (imgets.value=="0") {
            printlog ("WARNING - NIREDUCE: Cannot read filter from dark \
                current image.", l_logfile, verbose+)
        } else {
            # Add the keyword value "Dark" to this list so that this warning 
            # is not printed when reducing FLAMINGOS-II data
            if ((imgets.value != "blank") && (imgets.value != "BLANK") && \
                (imgets.value != "Dark")) {
                printlog("WARNING - NIREDUCE: dark current image not \
                    taken with filter wheel blanked.", l_logfile, verbose+)
            }
        }

        # get median dark level
        imstat (l_darkimage//"["//l_sci_ext//"]"//l_statsec,
            fields="midpt,stddev", lower=INDEF, upper=INDEF, nclip=0,
            lsigma=INDEF, usigma=INDEF, binwidth=0.01, format-, cache-) |\
            scan (l_mid, l_sig)
        imstat (l_darkimage//"["//l_sci_ext//"]"//l_statsec, fields="midpt",
            upper=(l_mid+2*l_sig), lower=(l_mid-2*l_sig), nclip=0,
            lsigma=INDEF, usigma=INDEF, binwidth=0.01, format-, cache-) |\
            scan (l_mid)
        if (l_mid==INDEF) {
            printlog ("ERROR - NIREDUCE: Statistics failed, possibly due to \
                a bad statsec.", l_logfile, verbose+)
            status = 1
            goto clean
        }
        l_middark = l_mid
    }

    #--------------------------------------------------------------------------
    # The math and bookkeeping:  (MAIN LOOP)

    printlog ("  n      input file -->      output file", l_logfile, l_verbose)
    printlog ("          sky image      dark image      flat image          \
        filter     sky   sky scale", l_logfile, l_verbose)
    printlog (" ", l_logfile, l_verbose)

    l_fl_first = yes  # flag for some warnings
    i = 1
    while (i<=nimages) {

        l_expression = "a"
        l_varexpression = "e"
        l_dqlist = ""
        ndq = 0

        #----------------SKY----------------
        # check for previous sky subtraction and turn it off if necessary

        l_skyscale = 1.0
        if (l_fl_sky) {
            imgets (in[i]//"[0]", "SKYIMAGE", >& "dev$null")
            if (imgets.value != "0") {
                l_fl_skytemp = no
                l_skyimage = "none"
                printlog ("WARNING - NIREDUCE: Image "//in[i]//" has already \
                    been sky-subtracted", l_logfile, verbose+)
                printlog ("                    by nireduce.  Sky-subtraction \
                    NOT performed.", l_logfile, verbose+)
            }
            else {
                imgets (in[i]//"[0]", l_key_filter, >& "dev$null")
                l_filter = imgets.value
                if (l_filter!=l_skyfilter) {
                    printlog ("WARNING - NIREDUCE: The filter for image \
                        "//in[i]//" does not match", l_logfile, verbose+)
                    printlog ("                    the sky image filter.  \
                        Proceeding anyway.", l_logfile, verbose+)
                }
                if (l_fl_scalesky) {
                    imstat (in[i]//"["//l_sci_ext//"]"//l_statsec,
                        fields="midpt,stddev", lower=INDEF, upper=INDEF,
                        nclip=0, lsigma=INDEF, usigma=INDEF, binwidth=0.05,
                        format-, cache-) | scan (l_mid, l_sig)
                    imstat (in[i]//"["//l_sci_ext//"]"//l_statsec,
                        fields="midpt", lower=(l_mid-4.*l_sig),
                        upper=(l_mid+4.*l_sig), nclip=0, lsigma=INDEF,
                        usigma=INDEF, binwidth=0.05, format-, cache-) |\
                        scan (l_middata)
                    if (l_middata==INDEF) {
                        printlog ("ERROR - NIREDUCE: Statistics failed, \
                            possibly due to a bad statsec.",l_logfile, \
                            verbose+)
                        status = 1
                        goto clean
                    }
                    l_skyscale = l_middata/l_midsky
                }
                l_fl_skytemp = yes
                l_expression = "("//l_expression//"-b*"//l_skyscale//")"
                l_varexpression = "("//l_varexpression//"+f*"//(l_skyscale*\
                    l_skyscale)//")"
                l_dqlist = l_dqlist+l_skyimage//".fits["//l_dq_ext//"],"
                ndq = ndq+1
                if (l_fl_autosky) l_skylevel = l_midsky*l_skyscale
            }
            
            # end if(l_fl_sky)
        } else {
            l_skyimage = "none"
            l_fl_skytemp = no
        }

        #----------------DARK-----------------
        # check for previous dark subtraction and turn it off if necessary
        if (l_fl_dark) {
            imgets (in[i]//"[0]", "DARKIMAG", >& "dev$null")
            if (imgets.value != "0") {
                l_fl_darktemp = no
                l_darkimage = "none"
                printlog ("WARNING - NIREDUCE: Image "//in[i]//" has already \
                    been dark-subtracted", l_logfile, verbose+)
                printlog ("                    by nireduce.  Dark-subtraction \
                    NOT performed.", l_logfile, verbose+)
            } else {
                l_fl_darktemp = yes
                l_expression = "("//l_expression//"-c)"
                l_varexpression = "("//l_varexpression//"+g)"
                l_dqlist = l_dqlist+l_darkimage//".fits["//l_dq_ext//"],"
                ndq = ndq+1
            }
        } else {
            l_darkimage = "none"
            l_fl_darktemp = no
        }

        #----------------FLAT------------------
        # check for previous flat fielding, and turn it off if necessary
        l_filter = "none"
        l_flatim = "none"

        if (l_fl_flat) {
            imgets (in[i]//"["//l_sci_ext//"]", "FLATIMAG", >& "dev$null")
            if (imgets.value != "0") {
                l_fl_flattemp = no
                l_flatim = "none"
                printlog ("WARNING - NIREDUCE: Image "//in[i]//" has already \
                    been flat-fielded.", l_logfile, verbose+)
                printlog ("                    by nireduce.  Flat-fielding \
                    NOT performed.", l_logfile, verbose+)
                l_filter = "none"  # set l_filter to something, for logging
            } else {
                l_fl_flattemp = yes
                imgets (in[i]//"[0]", l_key_filter, >& "dev$null")
                l_filter = imgets.value
                if (l_filter != l_flatfilter) {
                    printlog ("WARNING - NIREDUCE: The filter for image \
                        "//in[i]//" does not match", l_logfile, verbose+)
                    printlog ("                    the flat field filter.  \
                        Proceeding anyway.", l_logfile, verbose+)
                }
                l_varexpression = "("+l_varexpression+"/(d*d))+"
                l_varexpression = l_varexpression//"(h/(d*d*d*d))*("
                l_varexpression = l_varexpression//l_expression//"*"//\
                    l_expression//")"
                l_expression = l_expression//"/d"
                l_dqlist = l_dqlist+l_flatimage//".fits["//l_dq_ext//"],"
                ndq = ndq+1
                l_flatim = l_flatimage
            }
            # end if(l_fl_flat)
        } else {   # NO flat fielding
            l_fl_flattemp = no
            l_flatim = "none"
        }

        #-------------SKY CONSTANT---------------
        # add sky constant to expression
        l_expression = l_expression//"+"//str(l_skylevel)

        #-------------DO THE MATH----------------

        printf ("%3.0d %15s --> %16s \n", i, in[i], out[i]) | scan (l_struct)
        printlog (l_struct, l_logfile, l_verbose)
        printf ("    %15s %15s %15s %15s %7.1f %9.3f \n", l_skyimage, \
            l_darkimage, l_flatim, l_filter, l_skylevel, l_skyscale) | \
            scan (l_struct)
        printlog (l_struct, l_logfile, l_verbose)

        # Do the math, finally!  Pack up the results in MEF.

        l_scilist = in[i]//"["//l_sci_ext//"] "
        if (l_skyimage=="none")
            l_scilist = l_scilist//"INDEF "
        else
            l_scilist = l_scilist+l_skyimage//"["//l_sci_ext//"] "

        if (l_darkimage=="none")
            l_scilist = l_scilist//"INDEF "
        else
            l_scilist = l_scilist+l_darkimage//"["//l_sci_ext//"] "
        if (l_flatim=="none")
            l_scilist = l_scilist//"INDEF"
        else
            l_scilist = l_scilist+l_flatimage//"["//l_sci_ext//"]"

        imcopy (in[i]//"[0]", out[i], verbose-)
        ima = ""; imb = ""; imc = ""; imd = ""
        ime = ""; imf = ""; img = ""; imh = ""
        print (l_scilist) | scan (ima, imb, imc, imd)
        imexpr (l_expression,out[i]//"["//l_sci_ext//",append]",
            ima, imb, imc, imd, outtype="real", verbose-)

        if (l_fl_vardq && l_fl_indvardq[i]) {

            l_dqlist = l_dqlist+in[i]//".fits["//l_dq_ext//"]"
            ndq = ndq+1

            l_varlist = in[i]//"["//l_sci_ext//"] "
            if (l_skyimage=="none")
                l_varlist = l_varlist//"INDEF "
            else
                l_varlist = l_varlist+l_skyimage//"["//l_sci_ext//"] "

            if (l_darkimage=="none")
                l_varlist = l_varlist//"INDEF "
            else
                l_varlist = l_varlist+l_darkimage//"["//l_sci_ext//"] "

            if (l_flatim=="none")
                l_varlist = l_varlist//"INDEF "
            else
                l_varlist = l_varlist+l_flatimage//"["//l_sci_ext//"] "

            l_varlist = l_varlist+in[i]//"["//l_var_ext//"] "
            if (l_skyimage=="none")
                l_varlist = l_varlist//"INDEF "
            else
                l_varlist = l_varlist+l_skyimage//"["//l_var_ext//"] "

            if (l_darkimage=="none")
                l_varlist = l_varlist//"INDEF "
            else
                l_varlist = l_varlist+l_darkimage//"["//l_var_ext//"] "

            if (l_flatim=="none")
                l_varlist = l_varlist//"INDEF"
            else
                l_varlist = l_varlist+l_flatimage//"["//l_var_ext//"]"


            ima = ""; imb = ""; imc = ""; imd = ""
            ime = ""; imf = ""; img = ""; imh = ""
            print (l_varlist) | scan (ima, imb, imc, imd, ime, imf, img, imh)
            imexpr (l_varexpression, out[i]//"["//l_var_ext//",append]", \
                ima, imb, imc, imd, ime, imf, img, imh, outtype="real", \
                verbose-)

            if (ndq==4) {
                addmasks (l_dqlist, tmpdq//".pl","im1 || im2 || im3 || im4")
            } else if (ndq==3) {
                addmasks (l_dqlist, tmpdq//".pl","im1 || im2 || im3")
            } else if (ndq==2) {
                addmasks (l_dqlist, tmpdq//".pl","im1 || im2")
            } else {
                imcopy (l_dqlist, tmpdq//".pl", verbose-)
            }
            # Get OBJECT keyword from input DQ frame and update output header
            imgets(in[i]//".fits["//l_dq_ext//"]","OBJECT", >>& "dev$null")
            if (imgets.value=="0")
                gemhedit (tmpdq//".pl", "OBJECT", "", \
                    "Name of the object observed")
            else
                gemhedit (tmpdq//".pl", "OBJECT", imgets.value, \
                    "Name of the object observed")

            gemhedit (tmpdq//".pl", "EXTVER", "1", "Extension version")
            imarith(tmpdq//".pl", "*", 1, out[i]//"["//l_dq_ext//",append]", \
                pixtype="short", verbose-)

        } # end if(l_fl_vardq)

        # update the header
        gemdate ()
        gemhedit (out[i]//"[0]", "NIREDUCE", gemdate.outdate,
            "UT Time stamp for NIREDUCE")
        gemhedit (out[i]//"[0]", "GEM-TLM", gemdate.outdate, 
            comment="UT Last modification with GEMINI")
        if (l_fl_skytemp)
            gemhedit (out[i]//"[0]", "SKYIMAGE", l_skyimage,
                "Sky image subtracted from raw data")
        if (l_fl_scalesky)
            gemhedit (out[i]//"[0]", "SKYSCALE", l_skyscale,
                "Scale factor by which sky was scaled")
        if (l_fl_darktemp)
            gemhedit (out[i]//"[0]", "DARKIMAG", l_darkimage,
                "Dark current image subtracted from raw data")
        if (l_fl_flattemp)
            gemhedit (out[i]//"[0]", "FLATIMAG", l_flatimage,
                "Flat field image used")
        if (l_skylevel!=0.0)
            gemhedit (out[i]//"[0]", "SKYCONST", l_skylevel,
                "Constant added after sky image subtraction")

        # fix saturation and non-linear levels in header 
        imgets (out[i]//"[0]", l_key_sat, >& "dev$null")
        if (imgets.value != "0") {
            l_sat = int(imgets.value)
            if (l_fl_sky) l_sat = l_sat - l_midsky + l_skylevel
            if (l_fl_dark) l_sat = l_sat - l_middark
            gemhedit (out[i]//"[0]", l_key_sat, l_sat, 
                comment="Saturation level in ADU")
        }
        imgets (out[i]//"[0]", l_key_nonlinear, >& "dev$null")
        if (imgets.value != "0") {
            l_nonlinear = int(real(imgets.value))
            if (l_fl_sky) l_nonlinear = l_nonlinear - l_midsky + l_skylevel
            if (l_fl_dark) l_nonlinear = l_nonlinear - l_middark
            gemhedit (out[i]//"[0]", l_key_nonlinear, l_nonlinear, 
                comment="Non-linear regime in ADU")
        }
        # fix read noise level in header (for sky/dark subtraction ONLY)
        imgets (out[i]//"[0]", l_key_ron, >& "dev$null")
        if (imgets.value != "0") {
            l_ron = real(imgets.value)
            if (l_fl_sky) {
                imgets (l_skyimage//"[0]", l_key_ron, >& "dev$null")
                if (imgets.value != "0") {
                    l_mid = real(imgets.value)
                    l_ron = sqrt(l_ron*l_ron + l_mid*l_mid)
                    l_ron = real(int(l_ron * 10.0))/10.
                }
            }
            if (l_fl_dark) {
                imgets (l_darkimage//"[0]", l_key_ron, >& "dev$null")
                if (imgets.value != "0") {
                    l_mid = real(imgets.value)
                    l_ron = sqrt(l_ron*l_ron + l_mid*l_mid)
                    l_ron = real(int(l_ron * 10.0))/10.
                }
            }
            gemhedit (out[i]//"[0]", l_key_ron, l_ron, 
                comment="Estimated read noise (electrons)")
        }

        i = i+1
        imdelete (tmpdq, verify-, >& "dev$null")
    }
    # end the main loop

    #--------------------------------------------------------------------------
    # Clean up
clean:
    printlog (" ", l_logfile, l_verbose)
    if (status==0) {
        printlog ("NIREDUCE exit status:  good.", l_logfile, l_verbose)
    }
    printlog ("------------------------------------------------------------\
        ----------------", l_logfile, l_verbose)

    scanfile = ""
    delete (tmpfile, verify-, >& "dev$null")
    imdelete (tmpdq, verify-, >& "dev$null")

end



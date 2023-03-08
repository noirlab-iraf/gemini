# Copyright(c) 2000-2012 Association of Universities for Research in Astronomy, Inc.

procedure nresidual(inimages)

# Propagate the DQ plane into subsequent images, so that any
# non-linear and saturated pixels flagged as bad will be
# recognized in subsequent images.  Residual images are most
# likely to be associated with non-linear or saturated pixels.
# 
# A log of the reductions is written in the logfile
# The header of the image will contain information about
# which other images' DQ planes were included with its DQ plane.
#
# Data Quality Array handling:
#   dq(b) = dq(a) OR dq(b)
#
# Version  Sept. 14, 2004 JJ, first version
#          Sept. 15, 2004 JJ, fixed file not found bug
#          Sept. 16, 2004 KL, fixed gemoffsetlist and l_proptime
#          Sept. 17, 2004 JJ, added key_time and key_date parameters
#          Sept. 17, 2004 JJ, changed logic to make all propagated pixels=1 (bad)

char  inimages     {prompt="Input NIRI image(s)"}                           # OLDP-1-input-primary-single-prefix=b
char  outimages    {"",prompt="Output image(s)"}                            # OLDP-1-output
char  outprefix    {"b",prompt="Prefix for output image(s)"}                # OLDP-4
real  proptime     {2.,prompt="Residual image propagation time (min)"}      # OLDP-2
char  logfile      {"",prompt="Logfile"}                                    # OLDP-1
char  sci_ext      {"SCI",prompt="Name or number of science extension"}     # OLDP-3
char  var_ext      {"VAR",prompt="Name or number of variance extension"}    # OLDP-3
char  dq_ext       {"DQ",prompt="Name or number of data quality extension"} # OLDP-3
char  key_time     {"UT",prompt="Header keyword for the time"}              # OLDP-3
char  key_date     {"DATE-OBS",prompt="Header keyword for the date"}        # OLDP-3
bool  verbose      {yes,prompt="Verbose"}               # OLDP-4
int   status       {0,prompt="Exit status (0=good)"}    # OLDP-4
struct *scanfile   {prompt="Internal use"}              # OLDP-4

begin

    char l_inimages, l_outimages
    char l_prefix, l_logfile, l_temp
    char in[1000], out[1000]
    char l_sci_ext, l_var_ext, l_dq_ext
    char l_key_time, l_key_date
    char tmplist, tmplist2, tmpout, tmpdq, tmpfile, keyfound
    int  i, ii, j, nimages, noutimages, maxfiles, nbad, nnew
    real l_proptime, obstime
    char obsdate, l_expression
    bool l_verbose, before
    char ima,imb,imc,imd,ime,imf,img,imh
    struct l_struct

    status = 0
    nimages = 0
    maxfiles = 1000
    tmpfile = mktemp("tmpin")
    tmpdq = mktemp("tmpdq")
    tmplist = mktemp("tmplist")
    tmplist2 = mktemp("tmplist2")
    tmpout = mktemp("tmpout")

    cache("imgets", "gemdate")

    # set the local variables
    l_inimages = inimages ; l_outimages = outimages
    l_verbose = verbose ; l_prefix = outprefix
    l_logfile = logfile
    l_sci_ext = sci_ext ; l_var_ext = var_ext ; l_dq_ext = dq_ext
    l_proptime = proptime*60. # gemoffsetlist expects seconds, not minutes
    l_key_time = key_time ; l_key_date = key_date

    #------------------------------------------------------------------------
    # Check for package log file or user-defined log file
    cache("niri")
    print(l_logfile) | scan(l_logfile)
    if (l_logfile == "" || l_logfile == " ") {
        l_logfile = niri.logfile
        print(l_logfile) | scan(l_logfile)
        if (l_logfile == "" || l_logfile == " ") {
            l_logfile = "niri.log"
            printlog("WARNING - NRESIDUAL: Both nresidual.logfile and \
                niri.logfile are", l_logfile, verbose+)
            printlog("                     Undefined.  Using " \
                // l_logfile, l_logfile, verbose+)
        }
    }
    # Open log file
    date | scan(l_struct)
    printlog("--------------------------------------------------------------\
        --------------",l_logfile, l_verbose)
    printlog("NRESIDUAL -- "//l_struct, l_logfile, l_verbose)
    printlog(" ",l_logfile, l_verbose)
    
    # Check to see if sci_ext and dq_ext are defined
    print(l_sci_ext) | scan(l_sci_ext)
    if (l_sci_ext == "" || l_sci_ext == " ") {
        printlog("ERROR - NRESIDUAL: Extension sci_ext is undefined.", \
            l_logfile, verbose+)
        status = 1
        goto clean
    }

    print(l_dq_ext) | scan(l_dq_ext)
    if (l_dq_ext == "" || l_dq_ext == " ") {
        printlog("ERROR - NRESIDUAL: Extension dq_ext is undefined.", \
            l_logfile, verbose+)
        status = 1
        goto clean
    }

    #-----------------------------------------------------------------------
    # Load up arrays of input name lists

    if (substr(l_inimages,1,1) == "@") {
        l_temp = substr(l_inimages,2,strlen(l_inimages))
        if (!access(l_temp)) {
            printlog("ERROR - NRESIDUAL: Input file "//l_temp//" not found.", \
                l_logfile, verbose+)
            status = 1
            goto clean
        }
    }

    nimages = 0; nbad = 0
    files(l_inimages,sort-, > tmpfile) 
    scanfile = tmpfile

    while (fscan(scanfile,l_temp) != EOF) {
        gimverify(l_temp)
        if (gimverify.status == 1) {
            printlog("ERROR - NRESIDUAL: File "//l_temp//" not found.", \
                l_logfile, verbose+)
            nbad+=1
        } else if (gimverify.status>1) {
            printlog("ERROR - NRESIDUAL: File "//l_temp//" not a MEF FITS \
                image.", l_logfile, verbose+)
            nbad+=1
        } else {
            keyfound = ""
            hselect(l_temp//"[0]","*PREPAR*",yes) | scan(keyfound)
            if (keyfound == "") {
                printlog("ERROR - NRESIDUAL: Image "//l_temp//" not \
                    *PREPAREd.", l_logfile, verbose+)
                nbad+=1
            } else {
                if (!imaccess(l_temp//"["//l_dq_ext//"]") ) {
                    printlog("ERROR - NRESIDUAL: input image "//l_temp//" \
                        does not have a DQ plane.", l_logfile, verbose+)
                    nbad+=1
                }
            }

            # strip .fits if present
            if (substr(l_temp,strlen(l_temp)-4,strlen(l_temp)) == ".fits")
                l_temp = substr(l_temp,1,(strlen(l_temp)-5))

            nimages = nimages+1
            in[nimages] = l_temp 
            if (nimages>1) {
                for (j=1; j<nimages; j+=1) {
                    if (in[nimages] == in[j]) {
                        printlog("WARNING - NRESIDUAL: Input image name \
                            "//in[nimages]//" repeated.", l_logfile, verbose+)
                        printlog("                     Not including it \
                            again.", l_logfile, verbose+)
                        nimages = nimages - 1
                    }
                }
            }
        }
    } # end while
    
    # Exit if problems found
    if (nbad > 0) {
        printlog("ERROR - NRESIDUAL: "//nbad//" image(s) either do not exist, \
            are not MEF files,", l_logfile, verbose+)
        printlog("                   do not have DQ planes, or have not been \
            run through *PREPARE.", l_logfile, verbose+)
        status = 1
        goto clean
    }
    if (nimages > maxfiles) {
        printlog("ERROR - NRESIDUAL: Maximum number of input images exceeded \
            ("//str(maxfiles)//")", l_logfile, verbose+ )
        status = 1
        goto clean
    }
    if (nimages == 0) {
        printlog("ERROR - NRESIDUAL: No valid input images.", l_logfile, \
            verbose+)
        status = 1
        goto clean
    }
    if (nimages == 1) {
        printlog("ERROR - NRESIDUAL: Only one input image.", l_logfile, \
            verbose+)
        status = 1
        goto clean
    }

    printlog("Processing "//nimages//" files", l_logfile, l_verbose)
    scanfile = "" ; delete(tmpfile, verify-, >& "dev$null")

    #--------------------------------------------------------------------------
    # Now for the output images
    # outimages could contain legal * if it is of a form like %st%stX%*.imh

    noutimages = 0
    if (l_outimages!="" && l_outimages!=" ") {
        if (substr(l_outimages,1,1) == "@") 
            scanfile = substr(l_outimages,2,strlen(l_outimages))
        else if (stridx("*",l_outimages)>0) {
            files(l_outimages,sort-) | 
                match(".hhd",stop+,print-,metach-, > tmpfile)
            scanfile = tmpfile
        } else {
            files(l_outimages,sort-, > tmpfile)
            scanfile = tmpfile
        }

        while (fscan(scanfile,l_temp) != EOF) {
            noutimages = noutimages+1
            if (noutimages > maxfiles) {
                printlog("ERROR - NRESIDUAL: Maximum number of output images \
                    exceeded ("//str(maxfiles)//")", l_logfile, verbose+)
                status = 1
                goto clean
            }

            out[noutimages] = l_temp 
            if (imaccess(out[noutimages])) {
                printlog("ERROR - NRESIDUAL: Output image \
                    "//out[noutimages]//" already exists", l_logfile, verbose+)
                status = 1
            }
        }

        if (status != 0) goto clean

        # if there are too many or too few output images exit with error
        if (nimages != noutimages) {
            printlog("ERROR - NRESIDUAL: Number of input and output images \
                are unequal.", l_logfile, verbose+)
            status = 1
            goto clean
        }
    } else { # use prefix instead
        print(l_prefix) | scan(l_prefix)
        if (l_prefix == "" || l_prefix == " ") {
            printlog("ERROR - NRESIDUAL: Neither output image name nor \
                output prefix is defined.", l_logfile, verbose+)
            status = 1
            goto clean
        }
        i = 1
        nnew = 1
        while (i <= nimages) {
            out[i] = l_prefix//in[i]
            for (j=1; j<i; j+=1) {
                if (out[i] == out[j]) {
                    printlog("WARNING - NRESIDUAL: Output image name \
                        "//out[i]//" repeated.", l_logfile, verbose+)
                    printlog("                     Appending _"//nnew//" to \
                        the output file name.", l_logfile, verbose+)
                    out[i] = out[i]//"_"//nnew
                    nnew = nnew+1
                }
            }

            if (imaccess(out[i])) {
                printlog("ERROR - NRESIDUAL: Output image "//out[i]//" \
                    already exists.", l_logfile, verbose+)
                status = 1
            }
            i = i+1
        }
        if (status != 0) goto clean
    }

    scanfile = "" ; delete(tmpfile, verify-, >& "dev$null")

    #--------------------------------------------------------------------------
    # The mask bookkeeping: (MAIN LOOP)
    
    imcopy(in[1]//"[0]", out[1], verbose-)
    imcopy(in[1]//"["//l_sci_ext//"]", out[1]//"["//l_sci_ext//",append]", \
        verbose-)
    imcopy(in[1]//"["//l_var_ext//"]", out[1]//"["//l_var_ext//",append]", \
        verbose-)
    imcopy(in[1]//"["//l_dq_ext//"]", out[1]//"["//l_dq_ext//",append]", \
        verbose-)
    
    # Print out the names of the images being processed
    i = 1
    while (i <= nimages) {
        printlog("  "//in[i], l_logfile, l_verbose)
        i+=1
    }

    i = 2
    while (i <= nimages) {

        l_expression = "if(im1"
        imcopy(in[i]//"[0]", out[i], verbose-)
        imcopy(in[i]//"["//l_sci_ext//"]", \
            out[i]//"["//l_sci_ext//",append]", verbose-)
        imcopy(in[i]//"["//l_var_ext//"]", \
            out[i]//"["//l_var_ext//",append]", verbose-)
        imcopy(in[i]//"["//l_dq_ext//"]", out[i]//"["//l_dq_ext//",append]", \
            verbose-)

        # get only the images that are BEFORE the one being worked on
        imgets(in[i]//"[0]",l_key_time, >& "dev$null")
        if (imgets.value != "0") {
            obstime = real(imgets.value)
        } else {
            printlog("ERROR - NRESIDUAL: Can't read time in image "//in[i], \
                l_logfile, verbose+)
            status = 1
            goto clean
        }
        imgets(in[i]//"[0]",l_key_date, >& "dev$null")
        if (imgets.value != "0") {
            obsdate = imgets.value
        } else {
            printlog("ERROR - NRESIDUAL: Can't read date in image "//in[i], \
                l_logfile, verbose+)
            status = 1
            goto clean
        }

        for (ii=1; ii<=nimages; ii+=1) {
            imgets(in[ii]//"[0]",l_key_date, >& "dev$null")
            if (imgets.value != "0") {
                if (imgets.value < obsdate) before = yes
                else if (imgets.value > obsdate) before = no
                else {
                    imgets(in[ii]//"[0]",l_key_time, >& "dev$null")
                    if (imgets.value != "0") {
                        if (real(imgets.value) < obstime) before = yes
                        else before = no 
                    } else {
                        printlog("ERROR - NRESIDUAL: Can't read time in \
                            image "//in[ii], l_logfile, verbose+)
                        status = 1
                        goto clean
                    }
                }
            }
            if (before) print(in[ii], >> tmplist)
            before = no
        }

        # keep only the data taken within proptime minutes of the target image
        gemoffsetlist("@"//tmplist, in[i], distance=INDEF, age=l_proptime, \
            fl_younger+, targetlist=tmpout, offsetlist=tmplist2, \
            key_date=l_key_date, key_time=l_key_time, logfile=l_logfile, \
            verbose-)

        # Combine the masks using OR
        if (access(tmpout)) {
            scanfile = tmpout
            ii = 2
            while (fscan(scanfile,l_temp)!=EOF) {
                print(l_temp//"["//l_dq_ext//"]", >> tmpdq)
                l_expression = l_expression//" || im0"+ii
                ii+=1
            }
            l_expression = l_expression//") then 1"
            print(in[i]//"["//l_dq_ext//"]", >> tmpdq)
            addmasks("@"//tmpdq, out[i]//"["//l_dq_ext//",overwrite]", \
                expr=l_expression)

            # put the original mask values back in
            addmasks(in[i]//"["//l_dq_ext//"],"//out[i]//"["//l_dq_ext//"]", \
                out[i]//"["//l_dq_ext//",overwrite]", \
                expr="if(im1>im2) then im1 else im2")
        }

        # update the header
        gemdate ()
        gemhedit(out[i]//"[0]", "NRESIDUA", gemdate.outdate, \
            "UT Time stamp for NRESIDUAL", delete-)
        gemhedit(out[i]//"[0]", "GEM-TLM", gemdate.outdate, \
            "UT Last modification with GEMINI", delete-)
        if (access(tmpout)) {
            scanfile = tmpout
            ii = 1
            while (fscan(scanfile,l_temp)!=EOF) {
                gemhedit(out[i]//"[0]", "INPUTDQ"//str(ii), \
                    l_temp//"["//l_dq_ext//"]", \
                    "Additional included DQ plane", delete-)
                ii+=1
            }
        }

        i = i+1
        delete(tmpdq//","//tmplist//","//tmplist2//","//tmpout, verify-, \
            >& "dev$null")
    } # end the main loop

#---------------------------------------------------------------------------
# Clean up
clean:
    printlog(" ", l_logfile, l_verbose)
    if (status == 0) {
        printlog("NRESIDUAL exit status: good.", l_logfile, l_verbose)
    }
    printlog("---------------------------------------------------------------\
        -------------", l_logfile, l_verbose)

    scanfile = ""
    delete(tmpfile, verify-, >& "dev$null")

end



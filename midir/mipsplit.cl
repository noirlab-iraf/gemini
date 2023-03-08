# Copyright(c) 2005-2011 Association of Universities for Research in Astronomy, Inc.

procedure mipsplit (inimages, outimages)

# This routine extracts the individual polarimetry images and 
# writes them to simple FITS files.
#
# Version:  October 18, 2005  KV wrote original script 
#
#           There is a header keyword POLANGLE for each extension that gives the 
#           values of the polarimetry plate angle in degrees.
#
#           Jan  10, 2006   KV added the WCS parameters to the image extension from the 
#                              primary header.
#

char    inimages    {prompt="Input Michelle polarimetry image(s)"}  # OLDP-1-input-primary-single-prefix=s
char    outimages   {prompt="Output image(s) base name"}            # OLDP-1-output
char    rawpath     {"", prompt="Path for in input images"}               # OLDP-4
bool    fl_register {no,prompt="Register images when combining"}        # OLDP-2
char    regions     {"[*,*]",prompt="Reference image regions used for registration (xregister)"}    # OLDP-4
bool    fl_single   {no,prompt="Register to a single image"}            # OLDP-4
bool    fl_stair    {yes, prompt="Correct channel offsets"}             # OLDP-4
char    logfile     {"",prompt="Logfile"}                           # OLDP-1
bool    verbose     {yes,prompt="Verbose"}                          # OLDP-4
int     status      {0, prompt="Exit status: (0=good, >0=bad)"}         # OLDP-4
struct  *scanfile   {"", prompt="Internal use only"}                    # OLDP-4

begin

    char    l_inputimages, l_outputimages, l_filename,l_rawpath,l_logfile
    char    in[100], out[100], header, exheader, instrument, l_temp, l_regions
    char    tmpfile, filename, paramstr, tmpwork, tmpwork1
    char    tmpregister1, tmpregister2,tmpregister3,tmpregister4,tmpshifts
    int     i, j, k, l, m, n, nnods, nwp[4], jwp, nbad
    int     nimages, maximages, noutimages, l_extensions
    int     modeflag, inod
    real    chlevel
    bool    l_verbose,l_fl_register,l_fl_stair,l_fl_single

    l_inputimages = inimages
    l_outputimages = outimages
    l_rawpath = rawpath
    l_fl_register=fl_register
    l_regions=regions
    l_fl_single=fl_single
    l_fl_stair=fl_stair
    l_verbose=verbose
    l_logfile=logfile

    tmpfile = mktemp ("tmpin")
    if (l_fl_stair || l_fl_register) {
      tmpwork = mktemp ("tmpwork")
      tmpwork1 = mktemp ("tmpwork1")
    }
    if (l_fl_register) {
      tmpshifts= mktemp ("tmpshifts")
      tmpregister1= mktemp ("tmpregister")
      tmpregister2= mktemp ("tmpregister")
      tmpregister3 = mktemp ("tmpregister")
      tmpregister4 = mktemp ("tmpregister")
    }

    cache ("gemdate")

    nimages = 0
    maximages = 100
    status = 0
    nwp[1] = 0
    nwp[2] = 0
    nwp[3] = 0
    nwp[4] = 0

    # Create the list of parameter/value pairs.  One pair per line.
    # All lines combined into one string.  Line delimiter is '\n'.
    paramstr =  "inimages       = "//inimages.p_value//"\n"
    paramstr =  "outimages      = "//outimages.p_value//"\n"
    paramstr += "rawpath        = "//rawpath.p_value//"\n"
    paramstr += "fl_register    = "//fl_register.p_value//"\n"
    paramstr += "regions        = "//regions.p_value//"\n"
    paramstr += "fl_single      = "//fl_single.p_value//"\n"
    paramstr += "fl_stair       = "//fl_stair.p_value//"\n"
    paramstr += "logfile        = "//logfile.p_value//"\n"
    paramstr += "verbose        = "//verbose.p_value

    # Assign a logfile name if not given.  Open logfile and start log.
    # Write parameter/value pairs ("paramstr") to log.
    gloginit( l_logfile, "mipsplit", "midir", paramstr, fl_append+,
        verbose=l_verbose )
    if (gloginit.status != 0) {
        status = gloginit.status
        goto exitnow
    }
    l_logfile = gloginit.logfile.p_value

    #----------------------------------------------------------------------
    # Check the rawpath name for a final /
    if (substr(l_rawpath,(strlen(l_rawpath)),(strlen(l_rawpath))) != "/")
        l_rawpath = l_rawpath//"/"
    if (l_rawpath=="/" || l_rawpath==" ")
        l_rawpath = ""

    #-----------------------------------------------------------------------
    # Load up arrays of input name lists
    # This version handles both *s, / and commas in l_inputimages

    # check that list file exists
    if (substr(l_inputimages,1,1) == "@") {
        l_temp = substr (l_inputimages, 2, strlen(l_inputimages))
        if (!access(l_temp) && !access(l_rawpath//l_temp)) {
            glogprint( l_logfile, "mipsplit", "status", type="error", 
                errno=101, str="Input file "//l_temp//" not found.",verbose+)
            status = 1
            goto clean
        }
    }

    # Count the number of in images
    # First, generate the file list if needed

    if (stridx("*",l_inputimages) > 0) {
        files (l_inputimages, > tmpfile)
        l_inputimages = "@"//tmpfile
    }

    if (substr(l_inputimages,1,1) == "@")
        scanfile = substr (l_inputimages, 2, strlen(l_inputimages))
    else {
        files (l_inputimages, sort-, > tmpfile)
        scanfile = tmpfile
    }

    i = 0
    while ((fscan(scanfile,l_filename) != EOF) && (i <= 100)) {

        i = i+1

        if (substr(l_filename,strlen(l_filename)-4,strlen(l_filename)) == ".fits")
            l_filename = substr (l_filename, 1, strlen(l_filename)-5)

        if (!imaccess(l_filename) && !imaccess(l_rawpath//l_filename)) {
            glogprint( l_logfile, "mipsplit", "status", type="error", 
                errno=101, str="Input image"//l_filename//" was not found.", 
                verbose+)
            status = 1
            goto clean
        } else {
            nimages = nimages + 1
            if (nimages > maximages) {
                glogprint( l_logfile, "mipsplit", "status", type="error",
                    errno=121, str="Maximum number of input images \
                    exceeded:"//maximages,verbose+)
                status = 1 
                goto clean
            }

            if ((l_rawpath == "") || (l_rawpath == " "))
                in[nimages] = l_filename
            else
                in[nimages] = l_rawpath//l_filename
        }
    }

    scanfile = ""
    delete (tmpfile, verify-, >& "dev$null")

    if (nimages == 0) {
        glogprint( l_logfile, "mipsplit", "status", type="error", errno=121,
            str="No input images defined.",verbose+)
        status = 1
        goto clean
    }

    # Now, do the same counting for the out file

    nbad = 0
    noutimages = 0
    if ((l_outputimages != "") && (l_outputimages != " ")) {
        if (substr(l_outputimages,1,1) == "@")
            scanfile = substr (l_outputimages, 2, strlen(l_outputimages))
        else {
            if (stridx("*",l_outputimages) > 0) {
                files (l_outputimages, sort-) | \
                    match (".hhd", stop+, print-, metach-, > tmpfile)
                scanfile = tmpfile
            } else {
                files (l_outputimages, sort-, > tmpfile)
                scanfile = tmpfile
            }
        }

        while (fscan(scanfile,l_filename) != EOF) {
            noutimages = noutimages + 1
            if (noutimages > maximages) {
                glogprint( l_logfile, "mipsplit", "status", type="error",
                    errno=121, str="Maximum number of output images \
                    exceeded:"//maximages,verbose+)
                status = 1
                goto clean
            }
            if (substr(l_filename,strlen(l_filename)-4,strlen(l_filename)) != ".fits") 
                l_filename = l_filename//".fits"
            # test for one of the output files, exit if it is there....
            out[noutimages] = \
                substr (l_filename, 1, strlen(l_filename)-5)//"_wp1im1.fits"
            if (imaccess(out[noutimages])) {
                glogprint( l_logfile, "mipsplit", "status", type="error",
                    errno=102, str="Output image "//l_filename//" already \
                    exists.",verbose+)
                nbad += 1
            }
            # The output file name is now just the first part of the name 
            # from here on.
            out[noutimages] = substr (l_filename, 1, strlen(l_filename)-5)
        }
        if (noutimages != nimages) {
            glogprint( l_logfile, "mipsplit", "status", type="error", 
                errno=121, str="Different number of in images ("//nimages//") \
                and out images ("//noutimages//")",verbose+)
            status = 1
            goto clean
        }

        scanfile = ""
        delete (tmpfile, verify-, >& "dev$null")

    } else {
       glogprint( l_logfile, "mipsplit", "status", type="error", errno=121,
                str="No output images were defined.",verbose+)
        status = 1
        goto clean
    }

    if (nbad > 0) {
        glogprint( l_logfile, "mipsplit", "status", type="error", errno=102,
            str=nbad//" image(s) already exist.",verbose+)
        status = 1
        goto clean
    }

    nbad = 0
    i = 1
    while (i <= nimages) {
        # check the primary header
        header = in[i]//"[0]"
        imgets (header, "INSTRUMENT")
        instrument = imgets.value
        if (instrument == "michelle") {
            imgets (header, "MPREPARE", >& "dev$null")
            if (imgets.value == "0") {
                glogprint( l_logfile, "mipsplit", "status", type="error",
                    errno=123, str="Image "//in[i]//" not MPREPAREd.",
                    verbose+)
                status = 1
                goto clean
            }
        }

        # find the observation mode
        #
        if (instrument == "michelle") {
            imgets (header, "MODE", >& "dev$null")
            if (imgets.value == "0") {
                glogprint (l_logfile, "mipsplit", "status", type="error",
                    errno=131, str="Could not find the MODE from the primary \
                    header.", verbose+)
                status = status+1
                goto nextimage
            }
            # Change these according to the Michelle "MODE" keywords.
            # I am not sure whether there are other "non-destructive" (nd) 
            # modes, or whether these keywords are all correct.
            modeflag = 0
            if (imgets.value == "chop-nod") modeflag=1
            if (imgets.value == "ndchop") modeflag=1
            if (imgets.value == "chop") modeflag=2
            if (imgets.value == "nod") modeflag=3
            if (imgets.value == "ndstare") modeflag=4
            if (imgets.value == "stare") modeflag=4
            if (modeflag == 0) {
                glogprint (l_logfile, "mipsplit", "status", type="error",
                    errno=132, str="Unrecognized MODE ("//imgets.value//") \
                    in the primary header.",verbose+)
                status = status+1
                goto nextimage
            }
        } else {
            glogprint (l_logfile, "mipsplit", "status", type="error",
                errno=138, str="Error: for image "//in[i]//" the instrument \
                used is not MICHELLE", verbose+)
            status = status+1
            goto nextimage
        }

        imgets(in[i]//"[0]","FILTERA", >& "dev$null")
        j=0
        k=0
        if (substr(imgets.value,1,4) == "Grid") j=1
        imgets(in[i]//"[0]","FILTERB", >& "dev$null")
        if (substr(imgets.value,1,4) == "Grid") k=1
        if (j == 0 && k == 0) {
            glogprint (l_logfile, "mipsplit", "status", type="error",
                errno=142, str="Error: input file "//in[i]//" is not a \
                polarimetry file", verbose+)
            status = status+1
            goto nextimage
        }

        # Count the number of extensions

        l_extensions = 1
        while (imaccess(in[i]//"["//l_extensions//"]")) {
            imgets (in[i]//"["//l_extensions//"]", "i_naxis")
            if ((modeflag == 1) || (modeflag == 2)) {
                if ((imgets.value != "3") && (l_extensions > 0)) {
                    status = status+1
                    goto nextimage
                }
            }
            if ((modeflag == 3) || (modeflag == 4)) {
                if ((imgets.value != "2") && (l_extensions > 0)) {
                    status = status+1
                    goto nextimage
                }
            }
            l_extensions = l_extensions+1
        }

        j = l_extensions - 1
        glogprint( l_logfile, "mipsplit", "engineering", type="string",
            str="Number of extensions is "//j, verbose=l_verbose)

        if (j < 1) {
            glogprint (l_logfile, "mipsplit", "status", type="error", 
                errno=123, str="No data extensions in file "//in[i]//".", 
                verbose+)
            status = status + 1
            goto nextimage
        }

        # WARNING: modeflag can change within the if (modeflat==?) blocks
        #        therefore, the if-else-if structure SHOULD NOT be used.

        if (modeflag == 1) {
            if (8*int(j/8) != j) {
                    glogprint( l_logfile, "mipsplit", "status", type="warning",
                        errno=123, str="Number of extensions for input \
                        file "//in[i]//" does not correspond to complete \
                        cycles of polarimetry.  Skipping the file.",
                        verbose+)
                goto nextimage
            }
            nnods = j/8
            if (2*int(nnods/2) != nnods) {
                if (nnods != 1) {
                    glogprint( l_logfile, "mipsplit", "status", type="warning",
                        errno=123, str="Number of extensions for input \
                        file "//in[i]//" does not correspond to complete \
                        nodsets.  Removing last unmatched nod position.",
                        verbose+)
                    l_extensions = l_extensions - 8
                    if (l_extensions == 0) {
                        print ("No useable data extensions in file \
                            "//in[i]//".")
                        goto nextimage
                    }
                } else
                    modeflag = 2
            }
            
            if (nnods != 1) {
                for (m=1; m <= nnods; m=m+1) {
                    for (j=8*(m-1)+1; j <= 8*m ; j=j+1) {
                        jwp = 0
                        k = j - 8*int(j/8)
                        if (k == 0) k=8
                        if (k == 1 || k == 4) jwp=1
                        if (k == 5 || k == 8) jwp=2
                        if (k == 2 || k == 3) jwp=3
                        if (k == 6 || k == 7) jwp=4
                        # These are flags for positions 0, 22.5, 45.0 and 
                        # 67.5 respectively
                        if ((jwp > 0) && (jwp < 5)) {
                            nwp[jwp] = nwp[jwp] + 1
                            filename = out[i]//"_wp"//jwp//"im"//nwp[jwp]

                            if (!imaccess(in[i]//"["//j//"]")) {
                              goto nextimage
                            }

                            imcopy (in[i]//"["//j//"][*,*,3]", filename, 
                                verbose-, >& "dev$null")
                            imcopy (in[i]//"["//j//"][*,*,3]", tmpwork, 
                                verbose-, >& "dev$null")

                            if (l_fl_stair) {
                              for (m=1; m <= 301; m=m+20) {
                                imstatistics(tmpwork//"["//str(m)//":"//str(m+19)//",1:240]",
                                  fields="midpt",lower=INDEF,upper=INDEF,format-,nclip=0,lsigma=3.,
                                  usigma=3.,binwidth=0.1,cache-) | scanf("%f",chlevel)
                                imarith(tmpwork//"["//str(m)//":"//str(m+19)//",1:240]","-",chlevel,
                                  tmpwork1)
                                imcopy(tmpwork1,tmpwork//"["//str(m)//":"//str(m+19)//",1:240]",verbose-, >& "dev$null")
                                imdelete(tmpwork1,verify-, >& "dev$null")
                              }
                              imdelete(filename,verify-, >& "dev$null")
                              imcopy (tmpwork,filename,verbose-, >& "dev$null")
                            }

                            if (l_fl_register) {
                              if (nwp[jwp] == 1) {
                                if (l_fl_single) {
                                  if (jwp == 1) imcopy(tmpwork,tmpregister1)
                                  if (jwp == 2) imcopy(tmpwork,tmpregister2)
                                  if (jwp == 3) imcopy(tmpwork,tmpregister3)
                                  if (jwp == 4) imcopy(tmpwork,tmpregister4)
                                }
                                else {
                                  if (jwp == 1) imcopy(tmpwork,tmpregister1)
                                  if (jwp == 2) imcopy(tmpwork,tmpregister2)
                                  if (jwp == 3) imcopy(tmpwork,tmpregister3)
                                  if (jwp == 4) imcopy(tmpwork,tmpregister4)
                                }
                                imcopy(tmpwork, tmpwork1, verbose-, 
                                    >& "dev$null")
                              }
                              else {
                                if (jwp == 1) {
                                  images.immatch.xregister(tmpwork, 
                                    tmpregister1, regions=l_regions, 
                                    shifts=tmpshifts, output=tmpwork1, 
                                    background="none", loreject=INDEF, 
                                    hireject=INDEF, apodize=0., filter="none", 
                                    append+, records="", 
                                    correlation="discrete", xwindow=11, 
                                    ywindow=11, xcbox=11, ycbox=11, 
                                    function="centroid", interp_type="poly5", 
                                    interact-, xlag=0, ylag=0, dxlag=0, 
                                    dylag=0)
                                }
                                if (jwp == 2) {
                                  images.immatch.xregister(tmpwork, 
                                    tmpregister2, regions=l_regions, 
                                    shifts=tmpshifts, output=tmpwork1, 
                                    background="none", loreject=INDEF, 
                                    hireject=INDEF, apodize=0., filter="none", 
                                    append+, records="", 
                                    correlation="discrete", xwindow=11, 
                                    ywindow=11, xcbox=11, ycbox=11, 
                                    function="centroid", interp_type="poly5", 
                                    interact-, xlag=0, ylag=0, dxlag=0, 
                                    dylag=0)
                                }
                                if (jwp == 3) {
                                  images.immatch.xregister(tmpwork, 
                                    tmpregister3, regions=l_regions, 
                                    shifts=tmpshifts, output=tmpwork1, 
                                    background="none", loreject=INDEF, 
                                    hireject=INDEF, apodize=0., filter="none", 
                                    append+, records="", 
                                    correlation="discrete", xwindow=11, 
                                    ywindow=11, xcbox=11, ycbox=11, 
                                    function="centroid", interp_type="poly5", 
                                    interact-, xlag=0, ylag=0, dxlag=0, 
                                    dylag=0)
                                }
                                if (jwp == 4) {
                                  images.immatch.xregister(tmpwork, 
                                    tmpregister4, regions=l_regions, 
                                    shifts=tmpshifts, output=tmpwork1, 
                                    background="none", loreject=INDEF, 
                                    hireject=INDEF, apodize=0., filter="none", 
                                    append+, records="", 
                                    correlation="discrete", xwindow=11, 
                                    ywindow=11, xcbox=11, ycbox=11, 
                                    function="centroid", interp_type="poly5", 
                                    interact-, xlag=0, ylag=0, dxlag=0, 
                                    dylag=0)
                                }
                              }
                              imdelete(filename,verify-, >& "dev$null")
                              imcopy (tmpwork1, filename, verbose-, 
                                  >& "dev$null")
                            }

                            imdelete(tmpwork,verify-, >& "dev$null")
                            imdelete(tmpwork1,verify-, >& "dev$null")

                            glogprint (l_logfile, "mipsplit", "task", 
                                type="string", str="  "//in[i]//"["//j//"]\
                                [*,*,3] --> "//filename, verbose=l_verbose)
                            if (jwp == 1) {
                                gemhedit (filename, "WPLATE", 0.0, "", delete-)
                                gemhedit (filename, "POLANGLE", 0.0, "", 
                                    delete-)
                            }
                            if (jwp == 2) {
                                gemhedit (filename, "WPLATE", 22.5, "", 
                                    delete-)
                                gemhedit (filename, "POLANGLE", 45.0, "",
                                    delete-)
                            }
                            if (jwp == 3) {
                                gemhedit (filename, "WPLATE", 45.0, "", 
                                    delete-)
                                gemhedit (filename, "POLANGLE", 90.0, "", 
                                    delete-)
                            }
                            if (jwp == 4) {
                                gemhedit (filename, "WPLATE", 67.5, "", 
                                    delete-)
                                gemhedit (filename, "POLANGLE", 135.0, "", 
                                    delete-)
                            }

                            # add static parameters (to be determined on sky..)

                            gemhedit (filename, "ANGROT", 0.0, "", delete-)
                            gemhedit (filename, "T", 1.0, "", delete-)
                            gemhedit (filename, "EPS", 1.0, "", delete-)
                            
                            # Copy WCS information to the stacked image from 
                            # the primary header
                            imgets(header,"CTYPE1")
                            if (imgets.value != "") {
                              gemhedit(filename,"WCSAXES",2,"Number of WCS \
                                  axes in the image")
                              imgets(header,"CTYPE1")
                              gemhedit(filename,"CTYPE1",imgets.value,"R.A. \
                                  in tangent plane projection")
                            }
                            imgets(header,"CRPIX1")
                            if (imgets.value != "")
                                gemhedit(filename, "CRPIX1", imgets.value, 
                                    "Ref pix of axis 1")
                            imgets(header,"CRVAL1")
                            if (imgets.value != "")
                                gemhedit(filename,"CRVAL1", imgets.value,
                                "RA at Ref pix in decimal degrees")
                            imgets(header,"CTYPE2")
                            if (imgets.value != "")
                                gemhedit(filename,"CTYPE2", imgets.value,
                                "DEC. in tangent plane projection")
                            imgets(header,"CRPIX2")
                            if (imgets.value != "")
                                gemhedit(filename,"CRPIX2", imgets.value,
                                "Ref pix of axis 2")
                            imgets(header,"CRVAL2")
                            if (imgets.value != "")
                                gemhedit(filename,"CRVAL2", imgets.value,
                                "DEC at Ref pix in decimal degrees")
                            imgets(header,"CD1_1")
                            if (imgets.value != "")
                                gemhedit(filename,"CD1_1", imgets.value,
                                "WCS matrix element 1 1")
                            imgets(header,"CD1_2")
                            if (imgets.value != "")
                                gemhedit(filename,"CD1_2", imgets.value,
                                "WCS matrix element 1 2")
                            imgets(header,"CD2_1")
                            if (imgets.value != "")
                                gemhedit(filename,"CD2_1", imgets.value,
                                "WCS matrix element 2 1")
                            imgets(header,"CD2_2")
                            if (imgets.value != "")
                                gemhedit(filename,"CD2_2", imgets.value,
                                "WCS matrix element 2 2")
                            imgets(header,"RADECSYS")
                            if (imgets.value != "")
                                gemhedit(filename,"RADECSYS", imgets.value,
                                "R.A./DEC. coordinate system reference")

                            gemdate ()
                            gemhedit (filename, "GEM-TLM", gemdate.outdate,
                             "UT Last modification with GEMINI", delete-)
                            gemhedit (filename, "MIPSPLIT", gemdate.outdate,
                                "UT Time stamp for MISPLIT", delete-)
                        }
                    }
                }
            }

        }
        # jump to here if there is a problem

nextimage:
        i = i+1
    }

clean:
    scanfile = ""
    delete (tmpfile, verify-, >& "dev$null")
    delete ("tmpin*", verify-, >& "dev$null")
    delete ("tmpshift*", verify-, >& "dev$null")
    delete ("tmpregister*", verify-, >& "dev$null")

    if (status==0)
        glogclose( l_logfile, "mipsplit", fl_success+, verbose=l_verbose )
    else
        glogclose( l_logfile, "mipsplit`", fl_success-, verbose=l_verbose )

exitnow:
    ;

end

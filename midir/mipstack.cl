# Copyright(c) 2006-2011 Association of Universities for Research in Astronomy, Inc.

procedure mipstack (inimages)

# This routine stacks the individual frames of a "prepared" 
# Michelle polarimetry file.
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
char    inimages    {prompt="Input Michelle polarimetry image(s)"}      # OLDP-1-input-primary-single-prefix=s
char    outimages   {"", prompt="Output image(s)"}                      # OLDP-1-output
char    outpref     {"s", prompt="Prefix for output image(s)"}          # OLDP-4
char    rawpath     {"", prompt="Path for input raw images"}            # OLDP-4
char    frametype   {"dif", prompt="Type of frame to combine (src, ref, dif)"}  # OLDP-2
char    combine     {"average", prompt="Combining images by average|sum"}   # OLDP-2
bool    fl_register {no,prompt="Register images when combining"}        # OLDP-2
char    regions     {"[*,*]",prompt="Reference image regions used for registration (xregister)"}    # OLDP-4
bool    fl_stair    {yes, prompt="Correct channel offsets"}             # OLDP-4
char    logfile     {"", prompt="Logfile"}                              # OLDP-1
bool    verbose     {yes, prompt="Verbose"}                             # OLDP-4
int     status      {0, prompt="Exit status: (0=good, >0=bad)"}         # OLDP-4
struct  *scanfile   {"", prompt="Internal use only"}                    # OLDP-4

begin

    char    l_inputimages, l_outputimages, l_filename, l_prefix, l_logfile
    char    l_rawpath, l_frametype, l_combine, l_regions
    char    in[100], out[100], header, exheader, instrument, l_temp
    char    tmpon, tmpfile, tmpfinal, tmpfile1, tmpfile2, tmphead, tmpwork
    char    tmpcomp1, tmpcomp2, tmpcomp3, tmpcomp4, tmpregister1, tmpregister2
    char    paramstr, tmplog, tmpshifts
    char    refim1,refim2,refim3,refim4
    char    keyfound
    int     flag1,flag2,flag3,flag4
    int     i, j, k, l, l_nodset, l_saveset, itotal, n_sig, n_ref, nnods, m
    int     nimages, maximages, noutimages, l_frames, l_extensions
    int     source, reference, nbadsets, badsetnumber[100], badflag, framevalue
    int     aframe, nbad
    real    exptime, norm, ave1, ave2, diff1, chlevel
    real    wplate[4], polangle[4]
    int     modeflag, inod
    bool    l_verbose,l_fl_register,l_fl_stair
    struct  l_struct

    tmpfile = mktemp("tmpin")
    tmphead = mktemp("tmphead")
    tmplog = mktemp("tmplog")
    tmpon = mktemp("tmpon")
    tmpcomp1 = mktemp("tmpcomp1")
    tmpcomp2 = mktemp("tmpcomp2")
    tmpcomp3 = mktemp("tmpcomp3")
    tmpcomp4 = mktemp("tmpcomp4")
    tmpshifts = mktemp("tmpshifts")

    l_verbose = verbose 
    l_inputimages = inimages
    l_outputimages = outimages
    l_logfile = logfile
    l_prefix = outpref
    l_rawpath = rawpath
    l_frametype = frametype
    l_combine = combine
    l_fl_register=fl_register
    l_regions=regions
    l_fl_stair=fl_stair

    cache ("gloginit", "gemdate")

    nimages = 0
    maximages = 100
    status = 0
    flag1=0
    flag2=0
    flag3=0
    flag4=0
    
    wplate[1] = 0.0
    wplate[2] = 22.5
    wplate[3] = 45.0
    wplate[4] = 67.5
    polangle[1] = 0.0
    polangle[2] = 45.0
    polangle[3] = 90.0
    polangle[4] = 135.0
    

    if (l_frametype == "dif")
        framevalue = 3  
    else if (l_frametype == "src")
        framevalue = 1  
    else if (l_frametype == "ref")
        framevalue = 2  

    if ((framevalue == 1) && (l_prefix != ""))
        l_prefix = "c"
    else if ((framevalue == 2) && (l_prefix != ""))
        l_prefix = "a"

    # Create the list of parameter/value pairs.  One pair per line.
    # All lines combined into one string.  Line delimiter is '\n'.
    paramstr =  "inimages       = "//inimages.p_value//"\n"
    paramstr += "outimages      = "//outimages.p_value//"\n"
    paramstr += "outpref        = "//outpref.p_value//"\n"
    paramstr += "rawpath        = "//rawpath.p_value//"\n"
    paramstr += "frametype      = "//frametype.p_value//"\n"
    paramstr += "combine        = "//combine.p_value//"\n"
    paramstr += "fl_register    = "//fl_register.p_value//"\n"
    paramstr += "regions         = "//regions.p_value//"\n"
    paramstr += "fl_stair       = "//fl_stair.p_value//"\n"
    paramstr += "logfile        = "//logfile.p_value//"\n"
    paramstr += "verbose        = "//verbose.p_value

    # Assign a logfile name if not given.  Open logfile and start log.
    # Write parameter/value pairs ("paramstr") to log.
    gloginit (l_logfile, "mipstack", "midir", paramstr, fl_append+,
        verbose=l_verbose)
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
        l_rawpath=""

    if (l_combine != "average" && l_combine != "sum") {
        glogprint (l_logfile, "mipstack", "status", type="error", errno=121,
            str="Bad combine parameter ("//l_combine//")", verbose+)
        status = 1
        goto clean
    }

    #-----------------------------------------------------------------------
    # Load up arrays of input name lists
    # This version handles both *s, / and commas in l_inputimages

    # check that list file exists
    if (substr(l_inputimages,1,1)=="@") {
        l_temp=substr(l_inputimages,2,strlen(l_inputimages))
        if (!access(l_temp) && !access(l_rawpath//l_temp)) {
            glogprint (l_logfile, "mipstack", "status", type="error",
                errno=101, str="Input file "//l_temp//" not found.", verbose+)
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
        scanfile = substr(l_inputimages,2,strlen(l_inputimages))
    else {
        files (l_inputimages, sort-, > tmpfile)
        scanfile = tmpfile
    }

    i = 0

    while ((fscan(scanfile,l_filename) != EOF) && (i <= 10)) {

        i = i+1

        if (substr(l_filename,strlen(l_filename)-4,strlen(l_filename)) == ".fits")
            l_filename = substr (l_filename, 1, strlen(l_filename)-5)

        if (!imaccess(l_filename) && !imaccess(l_rawpath//l_filename)) {
            glogprint (l_logfile, "mipstack", "status", type="error",
                errno=101, str="Input image"//l_filename//" was not found.",
                verbose+)
            status = 1
            goto clean
        } else {
            nimages = nimages + 1
            if (nimages > maximages) {
                glogprint (l_logfile, "mipstack", "status", type="error",
                    errno=121, str="Maximum number of input images \
                    exceeded:"//maximages, verbose+)
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
        glogprint (l_logfile, "mipstack", "status", type="error", errno=121,
            str="No input images defined.",verbose+)
        status = 1
        goto clean
    }

    # Now, do the same counting for the out file

    nbad = 0
    noutimages = 0
    if ((l_outputimages != "") && (l_outputimages != " ")) {
        if (substr(l_outputimages,1,1) == "@")
            scanfile = substr(l_outputimages,2,strlen(l_outputimages))
        else {
            if (stridx("*",l_outputimages) > 0) {
                files (l_outputimages, sort-) | \
                    match (".hhd", stop+, print-, metach-, > tmpfile)
                scanfile = tmpfile
            } else {
                files (l_outputimages, sort-, > tmpfile)
                scanfile=tmpfile
            }
        }

        while (fscan(scanfile,l_filename) != EOF) {
            noutimages = noutimages + 1
            if (noutimages > maximages) {
                glogprint (l_logfile, "mipstack", "status", type="error",
                    errno=121, str="Maximum number of output images \
                    exceeded:"//maximages, verbose+)
                status = 1
                goto clean
            }
            if (substr(l_filename,strlen(l_filename)-4,strlen(l_filename)) != ".fits") 
                l_filename=l_filename//".fits"

            out[noutimages]=l_filename

            if (imaccess(out[noutimages])) {
                glogprint (l_logfile, "mipstack", "status", type="error",
                    errno=102, str="Output image "//l_filename//" already \
                    exists.",verbose+)
                nbad += 1
            }
        }
        if (noutimages != nimages) {
            glogprint   (l_logfile, "mipstack", "status", type="error",
                errno=121, str="Different number of in images \
                ("//nimages//") and out images ("//noutimages//")", verbose+)
            status = 1
            goto clean
        }

        scanfile = ""
        delete (tmpfile, verify-, >& "dev$null")

    } else {    # If prefix is to be used instead of filename

        print (l_prefix) | scan (l_prefix)
        if ((l_prefix == "") || (l_prefix == " ")) {
            glogprint (l_logfile, "mipstack", "status", type="error",
                errno=121, str="Neither output image name nor output prefix \
                is defined.", verbose+)
            status = 1
            goto clean
        }

        i = 1
        while (i <= nimages) {
            fparse (in[i])
            out[i] = l_prefix//fparse.root//".fits"

            if (imaccess(out[i])) {
                glogprint (l_logfile, "mipstack", "status", type="error", 
                    errno=102, str="Output image "//out[i]//" already exists.",
                    verbose+)
                nbad += 1
            }
            i = i+1
        }
    }

    if (nbad > 0) {
        glogprint (l_logfile, "mipstack", "status", type="error", errno=102,
            str=nbad//" image(s) already exist.", verbose+)
        status = 1
        goto clean
    }

    nbad = 0
    i = 1
    while (i <= nimages) {
        # tmp images used within this loop
        tmpfinal = mktemp ("tmpfinal")
        tmpwork = mktemp ("tmpwork")

        imgets (in[i]//"[0]", "MIPSTACK", >& "dev$null")
        if (imgets.value != "0") {
            glogprint (l_logfile, "mipstack", "status", type="warning", 
                errno=123, str="File "//in[i]//" has already been stacked.",
                verbose=l_verbose)
            goto nextimage
        }
        imgets (in[i]//"[0]", "MISTACK", >& "dev$null")
        if (imgets.value != "0") {
            glogprint (l_logfile, "mipstack", "status", type="warning", 
                errno=123, str="File "//in[i]//" has already been stacked.",
                verbose=l_verbose)
            goto nextimage
        }

        glogprint (l_logfile, "mipstack", "task", type="string",
            str="  "//in[i]//" --> "//out[i], verbose=l_verbose)

        # check the primary FITS header
        header = in[i]//"[0]"
        imgets (header, "INSTRUMENT")
        instrument = imgets.value
        glogprint (l_logfile, "mipstack", "science", type="string",
            str="Instrument is:"//instrument, verbose=l_verbose)

        if (instrument == "michelle") {
            imgets (header, "MPREPARE", >& "dev$null")
            if (imgets.value == "0") {
                glogprint (l_logfile, "mipstack", "status", type="error",
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
                glogprint (l_logfile, "mipstack", "status", type="error",
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
                glogprint (l_logfile, "mipstack", "status", type="error",
                    errno=132, str="Unrecognized MODE ("//imgets.value//") \
                    in the primary header.",verbose+)
                status = status+1
                goto nextimage
            }
        } else {
            glogprint (l_logfile, "mipstack", "status", type="error",
                errno=131, str="The instrument used is not MICHELLE",
                verbose+)
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
            glogprint (l_logfile, "mipstack", "status", type="error",
                errno=142, str="Error: input file "//in[i]//" is not a \
                polarimetry file",verbose+)
            status = status+1
            goto nextimage
        }

        # Count the number of extensions

        l_extensions = 1
        while (imaccess(in[i]//"["//l_extensions//"]")) {
            imgets (in[i]//"["//l_extensions//"]", "i_naxis")
            if (modeflag == 1 || modeflag == 2) {
                if ((imgets.value != "3") && (l_extensions > 0)) {
                    glogprint (l_logfile, "mipstack", "status", type="error",
                        errno=123, str="Extension "//l_extension//" has "//\
                        imgets.value//" dimensions.  It should be 3.", 
                        verbose+)
                    status = status+1
                    goto nextimage
                }
            }
            if (modeflag == 3 || modeflag == 4) {
                if ((imgets.value != "2") && (l_extensions > 0)) {
                    glogprint (l_logfile, "mipstack", "status", type="error",
                        errno=123, str="Extension "//l_extension//" has "//\
                        imgets.value//" dimensions.  It should be 2.", 
                        verbose+)
                    status = status+1
                    goto nextimage
                }
            }
            l_extensions = l_extensions+1
        }

        j = l_extensions - 1
        glogprint (l_logfile, "mipstack", "engineering", type="string",
            str="Number of extensions is "//j, verbose=l_verbose)

        if (j < 1) {
            glogprint (l_logfile, "mipstack", "status", type="error",
                errno=123, str="No data extensions in file "//in[i]//".",
                verbose+)
            status = status + 1
            goto nextimage
        }

# Only modeflag = 1 is handled at the moment.  If needed the other 
# modes will be added later.

        if (modeflag == 1) {
            if (8*int(j/8) != j) {
                glogprint (l_logfile, "mipstack", "status", type="error",
                    errno=123, str="Number of extensions for input \
                    file "//in[i]//" does not correspond to complete \
                    cycles of the polarimetry mode.  Skipping file.",
                    verbose+)
                goto nextimage
            }
            nnods = j/8
            if (2*int(nnods/2) != nnods) {
                if (nnods != 1) {
                    glogprint (l_logfile, "mipstack", "status",
                        type="warning", errno=123, str="Number of extensions \
                        for input file "//in[i]//" does not correspond to \
                        complete nodsets.  Removing last unmatched nod \
                        position.", verbose+)
                    l_extensions = l_extensions - 8
                    if (l_extensions == 0) {
                        glogprint (l_logfile, "mipstack", "status",
                            type="error", errno=123, str="No useable data \
                            extensions in file "//in[i]//".", verbose+)
                        status = status + 1
                        goto nextimage
                    }
                } else
                    modeflag = 2
            }

            # in this initial routine I will NOT check for "BADNOD" 
            # flags....everything get stacked up
            if (nnods != 1) {
                for (m=1; m <= nnods; m=m+2) {
                    for (j=8*(m-1)+1; j <= 8*m ; j=j+1) {
                        tmpfile2 = mktemp("tmpfile2")
                        tmpfile1 = mktemp("tmpfile1")
                        if (l_fl_register) {
                          tmpregister1=mktemp("tmpregister1")
                          tmpregister2=mktemp("tmpregister2")
                          imcopy(in[i]//"["//j//"]"//"[*,*,"//framevalue//"\
                            ].fits",tmpregister1,verbose-, >& "dev$null")
                          imcopy(in[i]//"["//j+8//"]"//"[*,*,"//framevalue//"\
                            ].fits",tmpregister2,verbose-, >& "dev$null")
                          images.immatch.xregister(tmpregister2, tmpregister1, 
                            regions=l_regions, shifts=tmpshifts, 
                            output=tmpfile2, background="none", loreject=INDEF,
                            hireject=INDEF, apodize=0., filter="none", append+,
                            records="", correlation="discrete", xwindow=11, 
                            ywindow=11, xcbox=11, ycbox=11, 
                            function="centroid", interp_type="poly5", 
                            interact-, xlag=0, ylag=0, dxlag=0, dylag=0)
                          imarith(tmpregister1,"+",tmpfile2,tmpfile2,title="", 
                            divzero=0., hparams="", pixtype="", calctype="", 
                            verbose-, noact-)
                          imdelete(tmpregister1//","//tmpregister2,verify-, 
                            >& "dev$null")
                        }
                        else {
                          imarith(in[i]//"["//j//"]"//"[*,*,"//framevalue//"\
                              ].fits", "+", in[i]//"["//j+8//"]"//"[*,*,\
                              "//framevalue//"].fits", tmpfile2, title="", 
                              divzero=0., hparams="", pixtype="", calctype="", 
                              verbose-, noact-)
                        }
                        if (l_combine == "average") {
                            imarith (tmpfile2,"/", "2.0", tmpfile1, title="",
                                divzero=0., hparams="", pixtype="", 
                                calctype="", verbose-, noact-)
                        } else {
                            imcopy (tmpfile2, tmpfile1, verbose-, 
                                >& "dev$null")
                        }
                        k = j - 8*int(j/8)
                        if (k == 0) k=8
                        if (k == 1 || k == 4) {
                          if (l_fl_register) {
                            if (flag1 == 0) {
                              flag1=1
                              refim1=tmpfile1
                              print(tmpfile1, >> tmpcomp1)
                            }
                            else {
                              tmpregister1=mktemp("tmpregister1")
                              images.immatch.xregister(tmpfile1, refim1, 
                                regions=l_regions, shifts=tmpshifts, 
                                output=tmpregister1, background="none", 
                                loreject=INDEF, hireject=INDEF, apodize=0., 
                                filter="none", append+, records="", 
                                correlation="discrete", xwindow=11, ywindow=11,
                                xcbox=11, ycbox=11, function="centroid", 
                                interp_type="poly5", interact-, xlag=0, ylag=0,
                                dxlag=0, dylag=0)
                              print(tmpregister1, >> tmpcomp1)
                            }
                          }
                          else {
                            print(tmpfile1, >> tmpcomp1)
                          }
                        }
                        if (k == 2 || k == 3) {
                          if (l_fl_register) {
                            if (flag2 == 0) {
                              flag2=1
                              refim2=tmpfile1
                              print(tmpfile1, >> tmpcomp2)
                            }
                            else {
                              tmpregister1=mktemp("tmpregister1")
                              images.immatch.xregister(tmpfile1, refim2, 
                                regions=l_regions, shifts=tmpshifts, 
                                output=tmpregister1, background="none", 
                                loreject=INDEF, hireject=INDEF, apodize=0., 
                                filter="none", append+, records="", 
                                correlation="discrete", xwindow=11, ywindow=11,
                                xcbox=11, ycbox=11, function="centroid", 
                                interp_type="poly5", interact-, xlag=0, ylag=0,
                                dxlag=0, dylag=0)
                              print(tmpregister1, >> tmpcomp2)
                            }
                          }
                          else {
                            print(tmpfile1, >> tmpcomp2)
                          }
                        }
                        if (k == 5 || k == 8) {
                          if (l_fl_register) {
                            if (flag3 == 0) {
                              flag3=1
                              refim3=tmpfile1
                              print(tmpfile1, >> tmpcomp3)
                            }
                            else {
                              tmpregister1=mktemp("tmpregister1")
                              images.immatch.xregister(tmpfile1, refim3, 
                                regions=l_regions, shifts=tmpshifts, 
                                output=tmpregister1, background="none", 
                                loreject=INDEF, hireject=INDEF, apodize=0., 
                                filter="none", append+, records="", 
                                correlation="discrete", xwindow=11, ywindow=11,
                                xcbox=11, ycbox=11, function="centroid", 
                                interp_type="poly5", interact-, xlag=0, ylag=0,
                                dxlag=0, dylag=0)
                              print(tmpregister1, >> tmpcomp3)
                            }
                          }
                          else {
                            print(tmpfile1, >> tmpcomp3)
                          }
                        }
                        if (k == 6 || k == 7) {
                          if (l_fl_register) {
                            if (flag4 == 0) {
                              flag4=1
                              refim4=tmpfile1
                              print(tmpfile1, >> tmpcomp4)
                            }
                            else {
                              tmpregister1=mktemp("tmpregister1")
                              images.immatch.xregister(tmpfile1, refim4, 
                                regions=l_regions, shifts=tmpshifts, 
                                output=tmpregister1, background="none", 
                                loreject=INDEF, hireject=INDEF, apodize=0., 
                                filter="none", append+, records="", 
                                correlation="discrete", xwindow=11, ywindow=11,
                                xcbox=11, ycbox=11, function="centroid", 
                                interp_type="poly5", interact-, xlag=0, ylag=0,
                                dxlag=0, dylag=0)
                              print(tmpregister1, >> tmpcomp4)
                            }
                          }
                          else {
                            print(tmpfile1, >> tmpcomp4)
                          }
                        }
                        imdelete(tmpfile2, verify-, >& "dev$null")
                      }
                    }
              
                delete (tmplog, verify-, >& "dev$null" )
                imcombine ("@"//tmpcomp1, tmpfinal, combine=l_combine,
                    headers="", bpmasks="", rejmasks="", expmask="", sigmas="",
                    logfile="STDOUT", reject="none", project-, 
                    outtype="double", outlimits="", offsets="none",
                    masktype="none", maskvalue=0., blank=0., scale="none",
                    zero="none", weight="none", statsec="", expname="",
                    lthreshold=INDEF, hthreshold=INDEF, nlow=1, nhigh=1,
                    mclip+, lsigma=3., hsigma=3., rdnoise="0.", gain="1.",
                    snoise="0.", sigscale=0.1, pclip= -0.5, grow=0., >& tmplog)
                glogprint (l_logfile, "mipstack", "science", type="file",
                    str=tmplog, verbose=l_verbose )
                delete (tmplog, verify-, >& "dev$null" )

                if (l_fl_stair) {
                  for (m=1; m <= 301; m=m+20) {
                    imstatistics(tmpfinal//"["//str(m)//":"//str(m+19)//"\
                        ,1:240]", fields="midpt", lower=INDEF, upper=INDEF, 
                        format-,nclip=0,lsigma=3., usigma=3., binwidth=0.1, 
                        cache-) | scanf("%f",chlevel)
                    imarith(tmpfinal//"["//str(m)//":"//str(m+19)//",1:240]", 
                        "-", chlevel, tmpwork)
                    imcopy(tmpwork, tmpfinal//"["//str(m)//":"//str(m+19)//"\
                        ,1:240]",verbose-, >& "dev$null")
                    imdelete(tmpwork,verify-, >& "dev$null")
                  }
                }

                wmef (tmpfinal, out[i], extname="SCI", phu=header, verbose-,
                    >& "dev$null")
                delete (tmpfinal//".fits", verify-, >& "dev$null")

                tmpfinal = mktemp ("tmpfinal")
                delete (tmplog, verify-, >& "dev$null" )
                imcombine ("@"//tmpcomp3, tmpfinal, combine=l_combine,
                    headers="", bpmasks="", rejmasks="", expmask="", sigmas="",
                    logfile="STDOUT", reject="none", project-, 
                    outtype="double", outlimits="", offsets="none",
                    masktype="none", maskvalue=0., blank=0., scale="none", 
                    zero="none", weight="none", statsec="", expname="",
                    lthreshold=INDEF, hthreshold=INDEF, nlow=1, nhigh=1, 
                    mclip+, lsigma=3., hsigma=3., rdnoise="0.", gain="1.", 
                    snoise="0.", sigscale=0.1, pclip= -0.5, grow=0., >& tmplog)
                glogprint (l_logfile, "mipstack", "science", type="file",
                    str=tmplog, verbose=l_verbose )
                delete (tmplog, verify-, >& "dev$null" )

                if (l_fl_stair) {
                  for (m=1; m <= 301; m=m+20) {
                    imstatistics(tmpfinal//"["//str(m)//":"//str(m+19)//"\
                        ,1:240]", fields="midpt", lower=INDEF, upper=INDEF, 
                        format-, nclip=0, lsigma=3., usigma=3., binwidth=0.1, 
                        cache-) | scanf("%f", chlevel)
                    imarith(tmpfinal//"["//str(m)//":"//str(m+19)//",1:240]", 
                        "-", chlevel, tmpwork, >& "dev$null")
                    imcopy(tmpwork,tmpfinal//"["//str(m)//":"//str(m+19)//"\
                        ,1:240]",verbose-, >& "dev$null")
                    imdelete(tmpwork,verify-, >& "dev$null")
                  }
                }

                fxinsert (tmpfinal, out[i]//"[1]", "", verbose=l_verbose)
                delete (tmpfinal//".fits", verify-, >& "dev$null")

                tmpfinal = mktemp ("tmpfinal")
                delete (tmplog, verify-, >& "dev$null" )
                imcombine ("@"//tmpcomp2, tmpfinal, combine=l_combine,
                    headers="", bpmasks="", rejmasks="", expmask="", sigmas="",
                    logfile="STDOUT", reject="none", project-, 
                    outtype="double", outlimits="", offsets="none",
                    masktype="none", maskvalue=0., blank=0., scale="none",
                    zero="none", weight="none", statsec="", expname="",
                    lthreshold=INDEF, hthreshold=INDEF,
                    nlow=1,nhigh=1,mclip+,lsigma=3., hsigma=3., rdnoise="0.",
                    gain="1.", snoise="0.", sigscale=0.1, pclip= -0.5, grow=0.,
                     >& tmplog)
                glogprint (l_logfile, "mipstack", "science", type="file",
                    str=tmplog, verbose=l_verbose )
                delete (tmplog, verify-, >& "dev$null" )

                if (l_fl_stair) {
                  for (m=1; m <= 301; m=m+20) {
                    imstatistics(tmpfinal//"["//str(m)//":"//str(m+19)//"\
                        ,1:240]", fields="midpt", lower=INDEF, upper=INDEF, 
                        format-, nclip=0, lsigma=3., usigma=3., binwidth=0.1, 
                        cache-) | scanf("%f", chlevel)
                    imarith(tmpfinal//"["//str(m)//":"//str(m+19)//",1:240]", 
                        "-", chlevel, tmpwork, verbose-, >& "dev$null")
                    imcopy(tmpwork,tmpfinal//"["//str(m)//":"//str(m+19)//"\
                        ,1:240]",verbose-, >& "dev$null")
                    imdelete(tmpwork,verify-, >& "dev$null")
                  }
                }

                fxinsert (tmpfinal, out[i]//"[2]", "", verbose=l_verbose)
                delete (tmpfinal//".fits", verify-, >& "dev$null")

                tmpfinal = mktemp ("tmpfinal")
                delete (tmplog, verify-, >& "dev$null" )
                imcombine ("@"//tmpcomp4, tmpfinal, combine=l_combine,
                    headers="", bpmasks="", rejmasks="", expmask="", sigmas="",
                    logfile="STDOUT", reject="none", project-, 
                    outtype="double", outlimits="", offsets="none",
                    masktype="none", maskvalue=0., blank=0., scale="none",
                    zero="none", weight="none", statsec="", expname="",
                    lthreshold=INDEF, hthreshold=INDEF, nlow=1, nhigh=1,
                    mclip+, lsigma=3., hsigma=3., rdnoise="0.", gain="1.",
                    snoise="0.", sigscale=0.1, pclip= -0.5, grow=0., >& tmplog)
                glogprint (l_logfile, "mipstack", "science", type="file",
                    str=tmplog, verbose=l_verbose )
                delete (tmplog, verify-, >& "dev$null" )

                if (l_fl_stair) {
                  for (m=1; m <= 301; m=m+20) {
                    imstatistics(tmpfinal//"["//str(m)//":"//str(m+19)//"\
                        ,1:240]", fields="midpt", lower=INDEF, upper=INDEF,
                        format-, nclip=0, lsigma=3., usigma=3., binwidth=0.1,
                        cache-) | scanf("%f", chlevel)
                    imarith(tmpfinal//"["//str(m)//":"//str(m+19)//",1:240]", 
                        "-", chlevel, tmpwork, verbose-, >& "dev$null")
                    imcopy(tmpwork,tmpfinal//"["//str(m)//":"//str(m+19)//"\
                        ,1:240]",verbose-, >& "dev$null")
                    imdelete(tmpwork,verify-, >& "dev$null")
                  }
                }

                fxinsert(tmpfinal, out[i]//"[3]", "", verbose=l_verbose)

                #-----------------------------------
                ###### Update header keywords ######
                #-----------------------------------
                
                # Remove axis labels no longer used.
                for (k=1; k<=4; k=k+1) {
                    gemhedit (out[i]//"["//k//"]", "AXISLAB3", "", "",
                        delete=yes)
                    gemhedit (out[i]//"["//k//"]", "AXISLAB4", "", "",
                        delete=yes)
                }
                    
                # The following are the waveplate angle and the polarization 
                # selection angle
                
                for (k=1; k<=4; k=k+1) {
                    gemhedit (out[i]//"["//k//"]", "WPLATE", wplate[k], 
                        comment="Waveplate angle", delete-)
                    gemhedit (out[i]//"["//k//"]", "POLANGLE", polangle[k],
                        comment="Polarization selection angle", delete-)
                }
                
                # Static parameters of the instrument, to be determined later
                # angrot is the polarization angle offset, t is the 
                # waveplate transmission, e is the waveplate efficiency of
                # detecting polarization.
                
                for (k=1; k<=4; k=k+1) {
                    gemhedit (out[i]//"["//k//"]", "ANGROT", 0.0,
                        comment="Polarization angle offset", delete-)
                    gemhedit (out[i]//"["//k//"]", "T", 1.0,
                        comment="Waveplate transmission", delete-)
                    gemhedit (out[i]//"["//k//"]", "EPS", 1.0, 
                        comment="Waveplate efficiency of detecting \
                        polarization", delete-)
                }

                # Copy WCS information to the stacked image from the primary
                # header
                
                keyfound = ""
                hselect (header, "CTYPE1", yes) | scan (keyfound)
                if (keyfound != "") {
                    for (k=1; k <= 4; k=k+1) {
                        gemhedit (out[i]//"["//str(k)//"]", "WCSAXES", 2, 
                            comment="Number of WCS axes in the image",
                            delete-)
                        gemhedit (out[i]//"["//str(k)//"]", "CTYPE1", keyfound,
                            comment="R.A. in tangent plane projection", 
                            delete-)
                    }
                }
                keyfound = ""
                hselect (header, "CRPIX1", yes) | scan (keyfound)
                if (keyfound != "") {
                    for (k=1; k<=4; k=k+1) {
                        gemhedit (out[i]//"["//str(k)//"]", "CRPIX1", keyfound,
                            comment="Ref pix of axis 1", delete-)
                    }
                }
                keyfound = ""
                hselect (header, "CRVAL1", yes) | scan (keyfound)
                if (keyfound != "") {
                    for (k=1; k<=4; k=k+1) {
                        gemhedit (out[i]//"["//str(k)//"]", "CRVAL1", keyfound,
                            comment="RA at Ref pix in decimal degrees",
                            delete-)
                    }
                }
                keyfound = ""
                hselect (header, "CTYPE2", yes) | scan (keyfound)
                if (keyfound != "") {
                    for (k=1; k<=4; k=k+1) {
                        gemhedit (out[i]//"["//str(k)//"]", "CTYPE2", keyfound,
                            comment="DEC. in tangent plane projection",
                            delete-)
                    }
                }
                keyfound = ""
                hselect (header, "CRPIX2", yes) | scan (keyfound)
                if (keyfound != "") {
                    for (k=1; k<=4; k=k+1) {
                        gemhedit (out[i]//"["//str(k)//"]", "CRPIX2", keyfound,
                            comment="Ref pix of axis 2", delete-)
                    }
                }
                keyfound = ""
                hselect (header, "CRVAL2", yes) | scan (keyfound)
                if (keyfound != "") {
                    for (k=1; k<=4; k=k+1) {
                        gemhedit (out[i]//"["//str(k)//"]", "CRVAL2", keyfound,
                            comment="DEC at Ref pix in decimal degrees",
                            delete-)
                    }
                }
                keyfound = ""
                hselect (header, "CD1_1", yes) | scan (keyfound)
                if (keyfound != "") {
                    for (k=1; k<=4; k=k+1) {
                        gemhedit (out[i]//"["//str(k)//"]", "CD1_1", keyfound,
                            comment="WCS matrix element 1 1", delete-)
                    }
                }
                keyfound = ""
                hselect (header, "CD1_2", yes) | scan (keyfound)
                if (keyfound != "") {
                    for (k=1; k<=4; k=k+1) {
                        gemhedit (out[i]//"["//str(k)//"]", "CD1_2", keyfound,
                            comment="WCS matrix element 1 2", delete-)
                    }
                }
                keyfound = ""
                hselect (header, "CD2_1", yes) | scan (keyfound)
                if (keyfound != "") {
                    for (k=1; k<=4; k=k+1) {
                        gemhedit (out[i]//"["//str(k)//"]", "CD2_1", keyfound,
                            comment="WCS matrix element 2 1", delete-)
                    }
                }
                keyfound = ""
                hselect (header, "CD2_2", yes) | scan (keyfound)
                if (keyfound != "") {
                    for (k=1; k<=4; k=k+1) {
                        gemhedit (out[i]//"["//str(k)//"]", "CD2_2", keyfound,
                            comment="WCS matrix element 2 2", delete-)
                    }
                }
                keyfound = ""
                hselect (header, "RADECSYS", yes) | scan (keyfound)
                if (keyfound != "") {
                    for (k=1; k<=4; k=k+1) {
                        gemhedit (out[i]//"["//str(k)//"]", "RADECSYS", 
                            keyfound, comment="R.A./DEC. coordinate system \
                            reference", delete-)
                    }
                }
                
                # EXTNAME and EXTVER
                for (k=1; k<=4; k=k+1) {
                    gemhedit (out[i]//"["//str(k)//"]", "EXTNAME", "SCI", 
                        comment="Extension name", delete-)
                    gemhedit (out[i]//"["//str(k)//"]", "EXTVER", k,
                        comment="Extension version", delete-)
                }
                                                       
                # Time stamps
                gemdate ()
                gemhedit (out[i]//"[0]", "MIPSTACK", gemdate.outdate,
                    "UT Time stamp for MIPSTACK", delete-)
                gemhedit (out[i]//"[0]", "GEM-TLM", gemdate.outdate,
                    "UT Last modification with GEMINI", delete-)
            }
        }
        # jump to here if there is a problem

nextimage:
        i = i+1
        delete (tmpfinal//".fits", verify-, >& "dev$null")
    }

clean:
    scanfile = ""
    imdelete ("@"//tmpcomp4, verify-, >& "dev$null")
    imdelete ("@"//tmpcomp3, verify-, >& "dev$null")
    imdelete ("@"//tmpcomp2, verify-, >& "dev$null")
    imdelete ("@"//tmpcomp1, verify-, >& "dev$null")
    delete (tmpfile//","//tmphead, verify-, >& "dev$null")
    delete ("tmpin*", verify-, >& "dev$null")
    delete ("tmpfinal*", verify-, >& "dev$null")
    delete ("tmpfile*", verify-, >& "dev$null")
    delete ("tmpcomp*", verify-, >& "dev$null")
    delete("tmpregister*",verify-, >& "dev$null")
    delete("tmpshifts*",verify-, >& "dev$null")
    if (modeflag == 3) {
      imdelete("@"//tmpon, verify-, >& "dev$null")
      delete ("tmpon*", verify-, >& "dev$null")
    }

    if (status == 0)
        glogclose (l_logfile, "mipstack", fl_success+, verbose=l_verbose )
    else
        glogclose (l_logfile, "mipstack", fl_success-, verbose=l_verbose )

exitnow:
    ;

end

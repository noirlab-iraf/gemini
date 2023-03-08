# Copyright(c) 2003-2011 Association of Universities for Research in Astronomy, Inc.

procedure miregister(inimages)

# Combine all the frames in a "prepared" file using registration to the
# first good frame.
#
# Version:  Sept 29, 2003   KV changed "mistack.cl" to "miregister.cl" 
#                              by adding registration
#           Oct   1, 2003   KV added code for all Michelle observing modes
#           Oct   3, 2003   KV made small syntax corrections
#           Oct   5, 2003   KV made changes to allowed MODE FITS keywords for Michelle
#           Oct  29, 2003   KL IRAF 2.12 - new/modified parameters
#                                imcombine: rejmask->rejmasks, plfile->nrejmasks
#                                           headers,bpmasks,expmasks,outlimits
#                              hedit: addonly
#           Dec  10, 2003   KV clean up tmp files without using delete("tmp*")
#           Aug  19, 2004   KV added "combine" option--average or sum, 
#                              and explicitly set imcombine options 
#                              weight, scale, offsets to "none"
#           Aug  23, 2004   KV changed the script to handle incomplete observations
#                              if possible.
#           Jan  10, 2006   KV added the WCS parameters to the image extension from the 
#                              primary header.  Added a variance plane as well.
#           Apr  14, 2006   KV bug fix of WCS code.
#

char    inimages    {prompt="Input T-ReCS or Michelle image(s)"} # OLDP-1-input-primary-single-prefix=x
char    outimages   {"",prompt="Output image(s)"}                   # OLDP-1-output
char    outpref     {"x",prompt="Prefix for output image(s)"}       # OLDP-4
char    rawpath     {"",prompt="Path for input raw images"}         # OLDP-4
char    combine     {"average",prompt="Combining images by average|sum"}    # OLDP-2
bool    fl_variance {no,prompt="Output variance frame"}             # OLDP-4
char    regions     {"[*,*]",prompt="Regions to be used for registration"}
char    logfile     {"",prompt="Logfile"}                           # OLDP-1
bool    verbose     {yes,prompt="Verbose"}                          # OLDP-4
int     status      {0,prompt="Exit error status: (0=good, >0=bad)"}# OLDP-4
struct* scanfile    {"",prompt="Internal use only"}                 # OLDP-4

begin

    char    l_inputimages,l_outputimages,l_filename,l_prefix,l_logfile
    char    l_rawpath,l_combine
    char    in[100],out[100],header,exheader,instrument, l_temp
    char    tmpon,tmpfile,tmpfinal,tmpfile1,tmpfile2,tmpvar
    char    tmpreva,tmprevb,tmprevc
    char    tmprefimage,tmpshift1,tmpshift2,tmpregister,tmplog
    char    paramstr, frag, l_regions
    int     i,j,k,l,l_nodset,l_saveset,itotal
    int     n_sig,n_ref,nimages,maximages,noutimages,l_frames,l_extensions
    int     source,reference,nbadsets,badsetnumber[100],badflag
    int     aframe, nbad
    int     refflag
    real    exptime,norm,ave1,ave2,diff1
    int     modeflag
    bool    l_verbose,l_fl_variance

    tmpfile=mktemp("tmpin")
    tmpon=mktemp("tmpon")
    tmplog=mktemp("tmplog")
    tmpregister=mktemp("tmpregister")

    l_verbose=verbose 
    l_inputimages=inimages
    l_outputimages=outimages
    l_logfile=logfile
    l_prefix=outpref
    l_rawpath=rawpath
    l_combine=combine
    l_regions=regions
    l_fl_variance=fl_variance

    cache ("gemdate")

    status=0
    nimages=0
    maximages=100
    status=0

    # Create the list of parameter/value pairs.  One pair per line.
    # All lines combined into one string. Line delimiter is '\n'.
    paramstr =  "inimages       = "//inimages.p_value//"\n"
    paramstr += "outimages      = "//outimages.p_value//"\n"
    paramstr += "outpref        = "//outpref.p_value//"\n"
    paramstr += "rawpath        = "//rawpath.p_value//"\n"
    paramstr += "combine        = "//combine.p_value//"\n"
    paramstr += "fl_variance    = "//fl_variance.p_value//"\n"
    paramstr += "regions        = "//regions.p_value//"\n"
    paramstr += "logfile        = "//logfile.p_value//"\n"
    paramstr += "verbose        = "//verbose.p_value

    # Assign a logfile name if not given.  Open logfile and start log.
    # Write parameter/value pairs ("paramstr") to log.
    gloginit( l_logfile, "miregister", "midir", paramstr, fl_append+,
        verbose=l_verbose )
    if (gloginit.status != 0) {
        status = gloginit.status
        goto exitnow
    }
    l_logfile = gloginit.logfile.p_value

    #----------------------------------------------------------------------
    # Check the rawpath name for a final /
    if (substr(l_rawpath,(strlen(l_rawpath)),(strlen(l_rawpath))) != "/")
        l_rawpath=l_rawpath//"/"
    if (l_rawpath=="/" || l_rawpath==" ")
        l_rawpath=""

    if (l_combine != "average" && l_combine != "sum") {
        glogprint (l_logfile, "miregister", "status", type="error", errno=121,
            str="Unrecognized combine parameter ("//l_combine//")", verbose+)
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
            glogprint( l_logfile, "miregister", "status", type="error",
                errno=101, str="Input file "//l_temp//" not found.",verbose+)
            status=1
            goto clean
        }
    }


    # Count the number of in images
    # First, generate the file list if needed

    if (stridx("*",l_inputimages) > 0) {
        files(l_inputimages, > tmpfile)
        l_inputimages="@"//tmpfile
    }

    if (substr(l_inputimages,1,1)=="@")
        scanfile=substr(l_inputimages,2,strlen(l_inputimages))
    else {
        files(l_inputimages,sort-, > tmpfile)
        scanfile=tmpfile
    }

    i=0
    while (fscan(scanfile,l_filename) != EOF && i <= 100) {
        i=i+1

        if (substr(l_filename,strlen(l_filename)-4,strlen(l_filename)) == ".fits")
            l_filename=substr(l_filename,1,strlen(l_filename)-5)

        if (!imaccess(l_filename) && !imaccess(l_rawpath//l_filename)) {
            glogprint( l_logfile, "miregister", "status", type="error",
                errno=101, str="Input image"//l_filename//" was not found.",
                verbose+)
            status=1
            goto clean
        } else {
            nimages=nimages+1
            if (nimages > maximages) {
                glogprint( l_logfile, "miregister", "status", type="error",
                    errno=121, str="Maximum number of input images \
                    exceeded:"//maximages,verbose+)
                status=1 
                goto clean
            }
            if (l_rawpath=="" || l_rawpath==" ")
                in[nimages]=l_filename
            else
                in[nimages]=l_rawpath//l_filename
        }
    }

    scanfile=""

    delete(tmpfile,verify-,>& "dev$null")

    if (nimages == 0) {
        glogprint( l_logfile, "miregister", "status", type="error", errno=121,
            str="No input images defined.",verbose+)
        status=1
        goto clean
    }

    # Now, do the same counting for the out file

    nbad=0
    noutimages=0
    if ((l_outputimages != "") && (l_outputimages != " ")) {
        if (substr(l_outputimages,1,1) == "@")
            scanfile=substr(l_outputimages,2,strlen(l_outputimages))
        else {
            if (stridx("*",l_outputimages) > 0) {
                files(l_outputimages,sort-) | \
                    match(".hhd",stop+,print-,metach-, > tmpfile)
                scanfile=tmpfile
            } else {
                files(l_outputimages,sort-, > tmpfile)
                scanfile=tmpfile
            }
        }

        while (fscan(scanfile,l_filename) != EOF) {
            noutimages=noutimages+1
            if (noutimages > maximages) {
                glogprint( l_logfile, "miregister", "status", type="error",
                    errno=121, str="Maximum number of output images \
                    exceeded:"//maximages, verbose+)
                status=1
                goto clean
            }
            out[noutimages]=l_filename
            if (imaccess(out[noutimages])) {
                glogprint( l_logfile, "miregister", "status", type="error",
                    errno=102, str="Output image "//l_filename//" already \
                    exists.",verbose+)
                nbad+=1
            }
        }
        if (noutimages != nimages) {
            glogprint( l_logfile, "miregister", "status", type="error",
                errno=121, str="Different number of in images ("//nimages//") \
                and out images ("//noutimages//")", verbose+)
            status =1
            goto clean
        }

        scanfile=""
        delete(tmpfile,verify-, >& "dev$null")

    } else {    # If prefix is to be used instead of filename
    
        print(l_prefix) | scan(l_prefix)
        if (l_prefix=="" || l_prefix==" ") {
            glogprint( l_logfile, "miregister", "status", type="error",
                errno=121, str="Neither output image name nor output prefix \
                is defined.", verbose+)
            status=1
            goto clean
        }
        i=1
        while (i<=nimages) {
            fparse(in[i])
            out[i]=l_prefix//fparse.root//".fits"
            if (imaccess(out[i])) {
                glogprint( l_logfile, "miregister", "status", type="error",
                    errno=102, str="Output image "//out[i]//" already exists.",
                    verbose+)
                nbad+=1
            }
            i=i+1
        }
    }
    if (nbad > 0) {
        glogprint( l_logfile, "miregister", "status", type="error", errno=102,
            str=nbad//" image(s) already exist.", verbose+)
        status=1
        goto clean
    }


    #
    # The main loop: accumulate source and reference images
    #

    refflag=0
    nbad=0

    i=1
    while (i <= nimages) {
        # Create tmp FITS file names used within this loop
        tmpfinal=mktemp("tmpfinal")
        tmpvar=mktemp("tmpval")

        imgets(in[i]//"[0]","MIREGIST",>& "dev$null")
        if (imgets.value != "0") {
            glogprint( l_logfile, "miregister", "status", type="warning",
                errno=123, str="File "//in[i]//" has already been stacked.",
                verbose=l_verbose)
            goto nextimage
        }

        glogprint( l_logfile, "miregister", "task", type="string",
            str="  "//in[i]//" --> "//out[i],verbose=l_verbose)

        # check the primary FITS header
        header=in[i]//"[0]"

        imgets(header,"INSTRUMENT", >& "dev$null")
        instrument=imgets.value
        glogprint( l_logfile, "miregister", "science", type="string",
            str="Instrument is: "//instrument, verbose=l_verbose)

        if (instrument == "michelle") {
            imgets(header,"MPREPARE", >& "dev$null")
            if (imgets.value == "0") {
                glogprint( l_logfile, "miregister", "status", type="error",
                    errno=123, str="Image "//in[i]//" not MPREPAREd.",verbose+)
                status=1
                goto clean
            }
        } else if (instrument == "TReCS") {
            imgets(header,"TPREPARE", >& "dev$null")
            if (imgets.value == "0") {
                glogprint( l_logfile, "miregister", "status", type="error",
                    errno=123, str="Image "//in[i]//" not TPREPAREd.",verbose+)
                status=1
                goto clean
            }
        } else if (instrument == "CanariCam") {
            imgets(header,"TPREPARE", >& "dev$null")
            if (imgets.value == "0") {
                glogprint( l_logfile, "miregister", "status", type="error",
                    errno=123, str="Image "//in[i]//" not TPREPAREd.",verbose+)
                status=1
                goto clean
            }
        }

        # find the observation mode
        #
        if ((instrument == "TReCS") || (instrument == "CanariCam")) {
            imgets(header,"OBSMODE", >& "dev$null")
            if (imgets.value == "0") {
                glogprint( l_logfile, "miregister", "status", type="error",
                    errno=131, str="Could not find the OBSMODE from the \
                    primary header.",verbose+)
                status=status+1
                goto nextimage
            }
            modeflag=0
            if (imgets.value == "chop-nod") modeflag=1
            if (imgets.value == "chop") modeflag=2
            if (imgets.value == "nod") modeflag=3
            if (imgets.value == "stare") modeflag=4
            if (modeflag == 0) {
                glogprint( l_logfile, "miregister", "status", type="error",
                    errno=132, str="Unrecognized OBSMODE ("//imgets.value//") \
                    in the primary header.", verbose+ )
                status=status+1
                goto nextimage
            }
        } else {
            imgets(header,"MODE", >& "dev$null")
            if (imgets.value == "0") {
                glogprint( l_logfile, "miregister", "status", type="error",
                    errno=131, str="Could not find the MODE from the primary \
                    header.", verbose+ )
                status=status+1
                goto nextimage
            }
            # Change these according to the Michelle "MODE" keywords.
            # I am not sure whether there are other "non-destructive" (nd) 
            # modes, or whether these keywords are all correct.
            modeflag=0
            if (imgets.value == "chop-nod") modeflag=1
            if (imgets.value == "ndchop") modeflag=1
            if (imgets.value == "chop") modeflag=2
            if (imgets.value == "nod") modeflag=3
            if (imgets.value == "ndstare") modeflag=4
            if (imgets.value == "stare") modeflag=4
            if (modeflag == 0) {
                glogprint( l_logfile, "miregister", "status", type="error",
                    errno=132, str="Unrecognized MODE ("//imgets.value//") in \
                    the primary header.", verbose+ )
                status=status+1
                goto nextimage
            }
        }

        # Count the number of extensions

        l_extensions=1
        while (imaccess(in[i]//"["//l_extensions//"]")) {
            imgets(in[i]//"["//l_extensions//"]","i_naxis", >& "dev$null")
            if (modeflag == 1 || modeflag == 2) {
                if (imgets.value != "3" && l_extensions > 0) {
                    glogprint( l_logfile, "miregister", "status", type="error",
                        errno=123, str="Extension "//l_extension//" has "//\
                    imgets.value//" dimensions.  It should be 3.",verbose+)
                    status=status+1
                    goto nextimage
                }
            }
            if (modeflag == 3 || modeflag == 4) {
                if (imgets.value != "2" && l_extensions > 0) {
                    glogprint( l_logfile, "miregister", "status", type="error",
                        errno=123, str="Extension "//l_extension//" has "//\
                    imgets.value//" dimensions.  It should be 2.", verbose+)
                    status=status+1
                    goto nextimage
                }
            }
            l_extensions=l_extensions+1
        }

        j=l_extensions-1
        glogprint( l_logfile, "miregister", "engineering", type="string",
            str="Number of extensions is "//j, verbose=l_verbose)

        if (j < 1) {
            glogprint (l_logfile, "miregister", "status", type="error",
                errno=123, str="No data extensions in file "//in[i]//".",
                verbose+)
            status = status + 1
            goto nextimage
        }

        if (modeflag == 1) {
            if (2*(j/2) != j) {
                glogprint( l_logfile, "miregister", "status", type="warning",
                    errno=123, str="Number of extensions for input file "//\
                    in[i]//" does not correspond to complete nodsets.  \
                    Removing last unmatched nod position.", verbose+)
                l_extensions = l_extensions-1
                if (l_extensions == 0) {
                    glogprint (l_logfile, "miregister", "status", type="error",
                        errno=123, str="No useable data extensions in file "//\
                        in[i]//".", verbose+)
                    status = status + 1
                    goto nextimage
                }
            }

            for (j=1; j < l_extensions; j=j+2) {
                # Create tmp FITS file names created within this loop
                # ** DO NOT DELETE tmpfile1 at the bottom of the loop: 
                # ** the images are used later on.
                tmpfile1=mktemp("tmpfile1")
                tmpfile2=mktemp("tmpfile2")

                imgets(in[i]//"["//j//"]","BADNOD", >& "dev$null")
                if (imgets.value != "0") {
                    imgets(in[i]//"["//j+1//"]","BADNOD", >& "dev$null")
                    if (imgets.value != "0")
                        glogprint( l_logfile, "miregister", "engineering",
                            type="warning", errno=0, str="Nods "//j//" and "//\
                            j+1//" are both bad.)",verbose=l_verbose)
                    else
                        glogprint( l_logfile, "miregister", "engineering",
                            type="warning", errno=0, str="Nod "//j//\
                            " is bad.  (Omitting both Nod "//j//" and "//\
                            j+1//")",verbose=l_verbose)
                } else {
                    imgets(in[i]//"["//j+1//"]","BADNOD", >& "dev$null")
                    if (imgets.value != "0")
                        glogprint( l_logfile, "miregister", "engineering",
                            type="warning", errno=0, str="Nod "//j+1//" is \
                            bad.  (Omitting both Nod "//j//" and "//j+1//")",
                            verbose=l_verbose)
                    else {
                        if (refflag == 0) {
                            tmprefimage=mktemp("tmprefimage")
                            imcopy(in[i]//"["//j//"]"//"[*,*,3]",tmprefimage,
                                >& "dev$null")
                            refflag=1
                        }
                        tmpshift1=mktemp("tmpshift1")
                        tmpshift2=mktemp("tmpshift2")
                        images.immatch.xregister(in[i]//"["//j//"]"//"[*,*,3]",
                            tmprefimage, l_regions, tmpregister, 
                            output=tmpshift1, interac-, verbose-)
                        delete(tmpregister,verify-, >& "dev$null")
                        images.immatch.xregister(in[i]//"["//j+1//"]"//\
                            "[*,*,3]", tmprefimage, l_regions, tmpregister, 
                            output=tmpshift2, interac-, verbose-)
                        delete(tmpregister,verify-, >& "dev$null")
                        imarith(tmpshift1,"+",tmpshift2,tmpfile2,verbose-)
                        if (l_combine == "average")
                            imarith(tmpfile2,"/","2.0",tmpfile1,verbose-)
                        else
                            imcopy(tmpfile2,tmpfile1,verbose-, >& "dev$null")

                        print(tmpfile1, >> tmpon)

                        # Delete tmp files not needed any longer
                        # ** DO NOT DELETE tmpfile1, tmprefimage
                        imdelete(tmpfile2, verify-, >& "dev$null")
                        imdelete(tmpshift1//","//tmpshift2, verify-, 
                            >& "dev$null")
                    }
                }
            }

            delete( tmplog, verify-, >& "dev$null" )
            imcombine("@"//tmpon,tmpfinal,headers="",bpmasks="",rejmasks="",
                nrejmasks="",expmasks="",sigmas=tmpvar,logfile="STDOUT",
                combine=l_combine,reject="none",project=no,outtype="double",
                outlimits="", weight="none", offsets="none", scale="none",
                >& tmplog)
            glogprint( l_logfile, "miregister", "science", type="file",
                str=tmplog, verbose=l_verbose )
            delete( tmplog, verify-, >& "dev$null" )

            wmef(tmpfinal, out[i], extname="SCI", phu=header, verbose-, 
                >& "dev$null")
            if (l_fl_variance) {
              imfunction(tmpvar,tmpvar,"square",verbose-, >& "dev$null")
              frag=substr(out[i],strlen(out[i])-4,strlen(out[i]))
              if (frag != ".fits") {
                fxinsert(tmpvar,out[i]//".fits[1]","", verbose-)
              }
              else {
                fxinsert(tmpvar,out[i]//"[1]","",verbose-)
              }
            }

            # Time stamps
            gemdate ()
            gemhedit (out[i]//"[0]", "MIREGIST", gemdate.outdate,
                "UT Time stamp for MIREGISTER", delete-)            
            gemhedit (out[i]//"[0]", "GEM-TLM", gemdate.outdate,
                "UT Last modification with GEMINI", delete-)

            gemhedit (out[i]//"[1]", "AXISLAB4", "", "", delete=yes)
            gemhedit (out[i]//"[1]", "AXISLAB3", "", "", delete=yes)
        }

        if (modeflag == 2) {
            if (j > 1) {
                glogprint( l_logfile, "miregister", "status", type="error",
                    errno=123, str="More than one extension ("//j//") in CHOP \
                    mode.", verbose+ )
                status=status+1
                goto nextimage
            }

            imgets(in[i]//"[1]","BADNOD", >& "dev$null")

            if (imgets.value != "0") {
                glogprint( l_logfile, "miregister", "status", type="warning",
                    errno=0, str="The only NOD is marked as BAD.", verbose+ )
                status=status+1
                goto nextimage
            }

            tmpfile1=mktemp("tmpfile1")    # DO NOT delete tmpfile1
            imcopy(in[i]//"[1][*,*,3].fits",tmpfile1,verbose-)
            print(tmpfile1, >> tmpon)

            delete( tmplog, verify-, >& "dev$null" )
            imcombine("@"//tmpon,tmpfinal,headers="",bpmasks="",rejmasks="",
                nrejmasks="",expmasks="",sigmas=tmpvar,logfile="STDOUT",
                combine=l_combine,reject="none",project=no,outtype="double",
                outlimits="", weight="none", offsets="none", scale="none",
                >& tmplog)
            glogprint( l_logfile, "miregister", "science", type="file",
                str=tmplog, verbose=l_verbose )

            wmef(tmpfinal, out[i], extname="SCI", phu=header, verbose-, 
                >& "dev$null")
            if (l_fl_variance) {
              imfunction(tmpvar,tmpvar,"square",verbose-, >& "dev$null")
              fxinsert(tmpvar,out[i]//"[1]","",verbose-)
            }

            # Time stamps
            gemdate ()
            gemhedit (out[i]//"[0]", "MIREGIST", gemdate.outdate,
                "UT Time stamp for MIREGISTER", delete-)            
            gemhedit (out[i]//"[0]", "GEM-TLM", gemdate.outdate,
                "UT Last modification with GEMINI", delete-)

            gemhedit (out[i]//"[1]", "AXISLAB4", "", "", delete=yes)
            gemhedit (out[i]//"[1]", "AXISLAB3", "", "", delete=yes)
        }

        if (modeflag == 3) {
            if (2*(j/2) != j) {
                glogprint( l_logfile, "miregister", "status", type="warning",
                    errno=123, str="Number of extensions for input file "//\
                    in[i]//" does not correspond to complete nodsets.  \
                    Removing last unmatched nod position.", verbose+)
                l_extensions=l_extensions-1
                if (l_extensions == 0) {
                    glogprint (l_logfile, "miregister", "status", type="error",
                        errno=123, str="no useable data extensions in file "//\
                        in[i]//".", verbose+)
                    status=status+1
                    goto nextimage
                }
            }

            for (j=1; j < l_extensions; j=j+2) {
                # Create tmp FITS file names created within this loop
                # ** DO NOT DELETE tmpfile1 at the bottom of the loop: 
                # ** the images are used later on.
                tmpfile1=mktemp("tmpfile1")
                tmpfile2=mktemp("tmpfile2")

                imgets(in[i]//"["//j//"]","BADNOD", >& "dev$null")
                if (imgets.value != "0") {
                    imgets(in[i]//"["//j+1//"]","BADNOD", >& "dev$null")
                    if (imgets.value != "0")
                        glogprint( l_logfile, "miregister", "engineering",
                            type="warning", errno=0, str="Nods "//j//" and "//\
                            j+1//" are both bad.)", verbose=l_verbose)
                    else
                        glogprint( l_logfile, "miregister", "engineering",
                            type="warning", errno=0, str="Nod "//j//" is bad. \
                            (Omitting both Nod "//j//" and "//j+1//")",
                            verbose=l_verbose)
                } else {
                    imgets(in[i]//"["//j+1//"]","BADNOD", >& "dev$null")
                    if (imgets.value != "0")
                        glogprint( l_logfile, "miregister", "engineering",
                            type="warning", errno=0, str="Nod "//j+1//" is \
                            bad.  (Omitting both Nod "//j//" and "//j+1//")",
                            verbose=l_verbose)
                    else {
                        if (refflag == 0) {
                            tmprefimage=mktemp("tmprefimage")
                            imcopy(in[i]//"["//j//"]"//"[*,*,3]",tmprefimage,
                                >& "dev$null")
                            refflag=1
                        }
                        tmpshift1=mktemp("tmpshift1")
                        tmpshift2=mktemp("tmpshift2")
                        images.immatch.xregister(in[i]//"["//j//"]"//"[*,*,3]",
                            tmprefimage, l_regions, tmpregister, 
                            output=tmpshift1, interac-, verbose-)
                        delete(tmpregister,verify-, >& "dev$null")
                        images.immatch.xregister(in[i]//"["//j+1//"]"//\
                            "[*,*,3]", tmprefimage, l_regions, tmpregister, 
                            output=tmpshift2, interac-, verbose-)
                        delete(tmpregister,verify-, >& "dev$null")
                        imarith(tmpshift1,"+",tmpshift2,tmpfile2,verbose-)
                        if (l_combine == "average")
                            imarith(tmpfile2,"/","2.0",tmpfile1,verbose-)
                        else
                            imcopy(tmpfile2,tmpfile1,verbose-, >& "dev$null")

                        print(tmpfile1, >> tmpon)

                        # Delete tmp files not needed any longer
                        # ** DO NOT DELETE tmpfile1, tmprefimage
                        imdelete(tmpfile2, verify-, >& "dev$null")
                        imdelete(tmpshift1//","//tmpshift2, verify-, 
                            >& "dev$null")
                    }
                }
            }

            delete( tmplog, verify-, >& "dev$null" )
            imcombine("@"//tmpon,tmpfinal,headers="",bpmasks="",rejmasks="",
                nrejmasks="",expmasks="",sigmas=tmpvar,logfile="STDOUT",
                combine=l_combine,reject="none",project=no,outtype="double",
                outlimits="", weight="none", offsets="none", scale="none",
                >& tmplog)
            glogprint( l_logfile, "miregister", "science", type="file",
                str=tmplog, verbose=l_verbose )

            wmef(tmpfinal, out[i], extname="SCI", phu=header, verbose-, 
                >& "dev$null")
            if (l_fl_variance) {
              imfunction(tmpvar,tmpvar,"square",verbose-, >& "dev$null")
              fxinsert(tmpvar,out[i]//"[1]","",verbose-)
            }


            # Time stamps
            gemdate ()
            gemhedit (out[i]//"[0]", "MIREGIST", gemdate.outdate,
                "UT Time stamp for MIREGISTER", delete-)            
            gemhedit (out[i]//"[0]", "GEM-TLM", gemdate.outdate,
                "UT Last modification with GEMINI", delete-)

            gemhedit (out[i]//"[1]", "AXISLAB4", "", "", delete=yes)
            gemhedit (out[i]//"[1]", "AXISLAB3", "", "", delete=yes)
        }

        if (modeflag == 4) {
            if (j > 1) {
                glogprint( l_logfile, "miregister", "status", type="error",
                    errno=123, str="More than one extension ("//j//") in \
                    STARE mode.", verbose+ )
                status=status+1
                goto nextimage
            }

            imgets(in[i]//"[1]","BADNOD", >& "dev$null")
            if (imgets.value != "0") {
                glogprint( l_logfile, "miregister", "status", type="warning",
                errno=0, str="The only NOD is marked as BAD.", verbose+ )
                status=status+1
                goto nextimage
            }

            tmpfile1=mktemp("tmpfile1")    # DO NOT DELETE tmpfile1
            imcopy(in[i]//"[1][*,*].fits",tmpfile1,verbose-)
            print(tmpfile1, >> tmpon)

            delete( tmplog, verify-, >& "dev$null" )
            imcombine("@"//tmpon,tmpfinal,headers="",bpmasks="",rejmasks="",
                nrejmasks="",expmasks="",sigmas=tmpvar,logfile="STDOUT",
                combine=l_combine,reject="none",project=no,outtype="double",
                outlimits="", weight="none", offsets="none", scale="none",
                >& tmplog)
            glogprint( l_logfile, "miregister", "science", type="file",
                str=tmplog, verbose=l_verbose )

            wmef(tmpfinal, out[i], extname="SCI", phu=header, verbose-,
                >& "dev$null")
            if (l_fl_variance) {
              imfunction(tmpvar,tmpvar,"square",verbose-, >& "dev$null")
              fxinsert(tmpvar,out[i]//"[1]","",verbose-)
            }


            # Time stamps
            gemdate ()
            gemhedit (out[i]//"[0]", "MIREGIST", gemdate.outdate,
                "UT Time stamp for MIREGISTER", delete-)            
            gemhedit (out[i]//"[0]", "GEM-TLM", gemdate.outdate,
                "UT Last modification with GEMINI", delete-)

            gemhedit (out[i]//"[1]", "AXISLAB4", "", "", delete=yes)
            gemhedit (out[i]//"[1]", "AXISLAB3", "", "", delete=yes)
        }
        # Add extension number and name to variance plane
        if (l_fl_variance) {
            gemhedit (out[i]//"[2]", "EXTVER", 2, "Extension version")
            gemhedit (out[i]//"[2]", "EXTNAME", "VAR", "Extension name")
        }

        # Copy WCS information to the stacked image from the primary header
        imgets(header,"CTYPE1")
        if (imgets.value != "") {
          gemhedit(out[i]//"[1]","WCSAXES", 2, 
              "Number of WCS axes in the image")
          gemhedit(out[i]//"[1]", "CTYPE1", imgets.value,
              "R.A. in tangent plane projection")
        }
        imgets(header, "CRPIX1")
        if (imgets.value != "") 
            gemhedit(out[i]//"[1]", "CRPIX1", imgets.value, 
                "Ref pix of axis 1")
        imgets(header, "CRVAL1")
        if (imgets.value != "") 
            gemhedit(out[i]//"[1]", "CRVAL1", imgets.value, 
                "RA at Ref pix in decimal degrees")
        imgets(header, "CTYPE2")
        if (imgets.value != "") 
            gemhedit(out[i]//"[1]", "CTYPE2", imgets.value, 
                "DEC. in tangent plane projection")
        imgets(header, "CRPIX2")
        if (imgets.value != "") 
            gemhedit(out[i]//"[1]", "CRPIX2", imgets.value, 
                "Ref pix of axis 2")
        imgets(header, "CRVAL2")
        if (imgets.value != "") 
            gemhedit(out[i]//"[1]", "CRVAL2", imgets.value, 
                "DEC at Ref pix in decimal degrees")
        imgets(header, "CD1_1")
        if (imgets.value != "") 
            gemhedit(out[i]//"[1]", "CD1_1", imgets.value, 
                "WCS matrix element 1 1")
        imgets(header, "CD1_2")
        if (imgets.value != "") 
            gemhedit(out[i]//"[1]", "CD1_2", imgets.value, 
                "WCS matrix element 1 2")
        imgets(header, "CD2_1")
        if (imgets.value != "") 
            gemhedit(out[i]//"[1]", "CD2_1", imgets.value, 
                "WCS matrix element 2 1")
        imgets(header, "CD2_2")
        if (imgets.value != "") 
            gemhedit(out[i]//"[1]", "CD2_2", imgets.value, 
                "WCS matrix element 2 2")
        imgets(header, "RADECSYS")
        if (imgets.value != "") 
            gemhedit(out[i]//"[1]", "RADECSYS", imgets.value, 
                "R.A./DEC. coordinate system reference")

        # jump to here if there is a problem
nextimage:

        i=i+1
        delete(tmpfinal//".fits",verify-,>& "dev$null")
        delete(tmpvar//".fits",verify-,>& "dev$null")
        delete("tmpon*",verify-, >& "dev$null")
    }

clean:
    scanfile="" 
    delete("tmpin*",verify-, >& "dev$null")
    delete("tmpfinal*",verify-, >& "dev$null")
    delete("tmpfile1*",verify-, >& "dev$null")
    delete("tmprefimage*",verify-, >& "dev$null")
    delete("tmplog*",verify-, >& "dev$null")

    if (status==0)
        glogclose( l_logfile, "miregister", fl_success+, verbose=l_verbose )
    else 
        glogclose( l_logfile, "miregister", fl_success-, verbose=l_verbose )

exitnow:
    ;

end

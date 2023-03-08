# Copyright(c) 2003-2011 Association of Universities for Research in Astronomy, Inc.

procedure mistack(inimages)

# This routine stacks the individual frames of a "prepared" T-ReCS or
# Michelle file.
#
# Version:  June 23, 2003  KV wrote original script from BRs prototype 
#           Sept  5, 2003  TB made update for 3D file structure, added outprefix 
#                             option, omitted all TReCS "saveset" info.
#           Oct   1, 2003  KV changed the routine to handle all observing modes
#           Oct   3, 2003  KV made small changes to the parameter syntax
#           Oct   5, 2003  KV made changes to allowed MODE FITS keywords for Michelle
#           Oct  29, 2003  KL IRAF 2.12 - new/modified parameters
#                              imcombine: rejmask->rejmasks, plfile->nrejmasks
#                                         headers,bpmasks,expmasks,outlimits
#			       hedit: addonly
#           Dec   3, 2003  KV removed a 'delete("tmp*")' statement, went back 
#                              to specific deletes of tmp working files.
#           Jan  28, 2004   TB updated to include src|ref|dif stacking
#           Aug  19, 2004   KV added "combine" option--average or sum, 
#                              and explicitly set imcombine options 
#                              weight, scale, offsets to "none"
#           Aug  22, 2004   KV added check for no data extensions
#           Feb  17, 2005   KV fixed a bug caused by someone changing my original
#                              four "if" blocks for "modeflag" to if/else if 
#                              structure.  There was a good reason for my original 
#                              coding structure.
#           Jan  10, 2006   KV added the WCS parameters to the image extension from the 
#                              primary header.
#           Apr  14, 2006   KV bug fix in the WCS code.
#
char    inimages    {prompt="Input T-ReCS or Michelle image(s)"}    # OLDP-1-input-primary-single-prefix=s
char    outimages   {"",prompt="Output image(s)"}                   # OLDP-1-output
char    outpref     {"s",prompt="Prefix for out image(s)"}          # OLDP-4
char    rawpath     {"",prompt="Path for in raw images"}            # OLDP-4
char    frametype   {"dif",prompt="Type of frame to combine (src, ref, dif)"}   # OLDP-2
char    combine     {"average",prompt="Combining images by average|sum"}    # OLDP-2
bool    fl_variance {no,prompt="Output variance frame"}             # OLDP-4
char    logfile     {"",prompt="Logfile"}                           # OLDP-1
bool    verbose     {yes,prompt="Verbose"}                          # OLDP-4
int     status      {0,prompt="Exit status: (0=good, >0=bad)"}      # OLDP-4
struct* scanfile    {"",prompt="Internal use only"}                 # OLDP-4

begin

    char    l_inputimages, l_outputimages, l_filename, l_prefix, l_logfile
    char    l_rawpath, l_frametype, l_combine
    char    in[100],out[100],header,exheader,instrument, l_temp
    char    tmpon,tmpfile,tmpfinal,tmpfile1,tmpfile2,tmplist,tmplist2,tmpvar
    char    tmpinimg
    char    paramstr, tmplog, frag, errmsg, filename
    int     i,j,k,l,l_nodset,l_saveset,itotal,n_sig,n_ref
    int     nimages,maximages,noutimages,l_frames,l_extensions
    int     source,reference,nbadsets,badsetnumber[100],badflag, framevalue
    int     aframe, nbad
    real    exptime,norm,ave1,ave2,diff1
    int     modeflag,inod
    bool    l_verbose, l_fl_variance

    tmpfile=mktemp("tmpin")
    tmplist=mktemp("tmpfile")
    tmplist2=mktemp("tmpfile2")
    tmpon=mktemp("tmpon")
    tmplog=mktemp("tmplog")
    tmpinimg = mktemp ("tmpinimg")

    l_verbose=verbose 
    l_inputimages=inimages
    l_outputimages=outimages
    l_logfile=logfile
    l_prefix=outpref
    l_rawpath=rawpath
    l_frametype=frametype
    l_combine=combine
    l_fl_variance=fl_variance

    cache ("gemdate")

    nimages=0
    maximages=100
    status=0

    if (l_frametype=="dif")
        framevalue=3  
    else if (l_frametype=="src")
        framevalue=1  
    else if (l_frametype=="ref")
        framevalue=2  

    if ((framevalue == 1) && (l_prefix != ""))
        l_prefix="c"
    else if ((framevalue == 2) && (l_prefix != ""))
        l_prefix="a"

    # Create the list of parameter/value pairs.  One pair per line.
    # All lines combined into one string.  Line delimiter is '\n'.
    paramstr =  "inimages       = "//inimages.p_value//"\n"
    paramstr += "outimages      = "//outimages.p_value//"\n"
    paramstr += "outpref        = "//outpref.p_value//"\n"
    paramstr += "rawpath        = "//rawpath.p_value//"\n"
    paramstr += "frametype      = "//frametype.p_value//"\n"
    paramstr += "combine        = "//combine.p_value//"\n"
    paramstr += "fl_variance    = "//fl_variance.p_value//"\n"
    paramstr += "logfile        = "//logfile.p_value//"\n"
    paramstr += "verbose        = "//verbose.p_value

    # Assign a logfile name if not given.  Open logfile and start log.
    # Write parameter/value pairs ("paramstr") to log.
    gloginit( l_logfile, "mistack", "midir", paramstr, fl_append+,
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
        glogprint (l_logfile, "mistack", "status", type="error", errno=121,
            str="Bad combine parameter ("//l_combine//")", verbose+)
        status = 1
        goto clean
    }

    #-----------------------------------------------------------------------
    # Load up arrays of input name lists

    gemextn (l_inputimages, check="exists,mef", process="none", index="",
        extname="", extversion="", ikparam="", omit="kernel,exten",
        replace="^%%"//l_rawpath//"%", outfile=tmpinimg, logfile=l_logfile,
        verbose=l_verbose)
    nimages = gemextn.count
    
    if ((gemextn.fail_count > 0) || (nimages == 0) || \
        (nimages > maximages)) {
        
        if (gemextn.fail_count > 0) {
            errmsg = gemextn.fail_count//" images were not found."
            status = 101
        } else if (nimages == 0) {
            errmsg = "No input images defined."
            status = 121
        } else if (nimages > maximages) {
            errmsg = "Maximum number of input images ("//str(maximages)//") \
                has been exceeded."
            status = 121
        }
        
        glogprint (l_logfile, "mistack", "status", type="error", errno=status,
            str=errmsg, verbose+)
        goto clean
    } else {
        scanfile = tmpinimg
        i = 0
        while (fscan (scanfile, filename) != EOF) {
            i += 1
            in[i] = filename
        }
        scanfile = ""
        if (i != nimages) {
            status 99
            errmsg = "Error while counting the input images."
            glogprint (l_logfile, "mistack", "status", type="error",
                errno=status, str=errmsg, verbose+)
        }
    }

    # Now, do the same counting for the out file

    nbad=0
    noutimages=0
    if ((l_outputimages != "") && (l_outputimages != " ")) {
	if (substr(l_outputimages,1,1)=="@")
	    scanfile=substr(l_outputimages,2,strlen(l_outputimages))
	else { 
	    files(l_outputimages,sort-) | unique("STDIN", > tmpfile)
	    scanfile=tmpfile
	}

        while (fscan(scanfile,l_filename) != EOF) {
            noutimages=noutimages+1
            if (noutimages > maximages) {
                glogprint( l_logfile, "mistack", "status", type="error",
                    errno=121, str="Maximum number of output images \
                    exceeded:"//maximages,verbose+)
                status=1
                goto clean
            }
            out[noutimages]=l_filename
            if (imaccess(out[noutimages])) {
                glogprint( l_logfile, "mistack", "status", type="error",
                    errno=102, str="Output image "//l_filename//" already \
                    exists.",verbose+)
                nbad+=1
            }
        }
        if (noutimages != nimages) {
            glogprint( l_logfile, "mistack", "status", type="error", errno=121,
                str="Different number of in images ("//nimages//") and out \
                images ("//noutimages//")",verbose+)
            status =1
            goto clean
        }

        scanfile=""
        delete(tmpfile,ver-, >& "dev$null")

    } else { 	# If prefix is to be used instead of filename

        print(l_prefix) | scan(l_prefix)
        if (l_prefix=="" || l_prefix==" ") {
            glogprint( l_logfile, "mistack", "status", type="error", errno=121,
                str="Neither output image name nor output prefix is defined.",
                verbose+)
            status=1
            goto clean
        }

        i=1
        while (i<=nimages) {
            fparse(in[i])
            out[i]=l_prefix//fparse.root//".fits"

            if (imaccess(out[i])) {
                glogprint( l_logfile, "mistack", "status", type="error", 
                    errno=102, str="Output image "//out[i]//" already exists.",
                    verbose+)
                nbad+=1
            }
            i=i+1
        }
    }

    if (nbad > 0) {
        glogprint( l_logfile, "mistack", "status", type="error", errno=102,
            str=nbad//" image(s) already exist.",verbose+)
        status=1
        goto clean
    }

    nbad=0
    i=1
    while (i <= nimages) {

        # tmp images used within this loop
        tmpvar=mktemp("tmpvar")
        tmpfinal=mktemp("tmpfinal")

        imgets(in[i]//"[0]","MISTACK", >& "dev$null")
        if (imgets.value != "0") {
            glogprint( l_logfile, "mistack", "status", type="warning", 
                errno=123, str="File "//in[i]//" has already been stacked.",
                verbose=l_verbose)
            goto nextimage
        }

        glogprint( l_logfile, "mistack", "task", type="string",
            str="  "//in[i]//" --> "//out[i],verbose=l_verbose)

        # check the primary FITS header
        header=in[i]//"[0]"
        imgets(header,"INSTRUMENT")
        instrument=imgets.value
        glogprint( l_logfile, "mistack", "science", type="string",
            str="Instrument is:"//instrument,verbose=l_verbose)

        if (instrument == "michelle") {
            imgets(header,"MPREPARE", >& "dev$null")
            if (imgets.value == "0") {
                glogprint( l_logfile, "mistack", "status", type="error",
                    errno=123, str="Image "//in[i]//" not MPREPAREd.",
                    verbose+)
                status=1
                goto clean
            }
        } else if (instrument == "TReCS") {
            imgets(header,"TPREPARE", >& "dev$null")
            if (imgets.value == "0") {
                glogprint( l_logfile, "mistack", "status", type="error",
                    errno=123, str="Image "//in[i]//" not TPREPAREd.",
                    verbose+)
                status=1
                goto clean
            }
        } else if (instrument == "CanariCam") {
            imgets(header,"TPREPARE", >& "dev$null")
            if (imgets.value == "0") {
                glogprint( l_logfile, "mistack", "status", type="error",
                    errno=123, str="Image "//in[i]//" not TPREPAREd.",
                    verbose+)
                status=1
                goto clean
            }
        }

        # find the observation mode
        #
        if ((instrument == "TReCS") || (instrument == "CanariCam")) {
            imgets(header,"OBSMODE", >& "dev$null")
            if (imgets.value == "0") {
                glogprint( l_logfile, "mistack", "status", type="error",
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
                glogprint( l_logfile, "mistack", "status", type="error",
                    errno=132, str="Unrecognized OBSMODE ("//imgets.value//") \
                    in the primary header.",verbose+)
                status=status+1
                goto nextimage
            }
        } else {
            imgets(header,"MODE", >& "dev$null")
            if (imgets.value == "0") {
                glogprint( l_logfile, "mistack", "status", type="error",
                    errno=131, str="Could not find the MODE from the primary \
                    header.",verbose+)
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
                glogprint( l_logfile, "mistack", "status", type="error",
                    errno=132, str="Unrecognized MODE ("//imgets.value//") \
                    in the primary header.",verbose+)
                status=status+1
                goto nextimage
            }
        }

        # Count the number of extensions

        l_extensions=1
        while (imaccess(in[i]//"["//l_extensions//"]")) {
            imgets(in[i]//"["//l_extensions//"]","i_naxis")
            if (modeflag == 1 || modeflag == 2) {
                if (imgets.value != "3" && l_extensions > 0) {
                    glogprint( l_logfile, "mistack", "status", type="error",
                        errno=123, str="Extension "//l_extension//" has "//\
                        imgets.value//" dimensions.  It should be 3.", verbose+)
                    status=status+1
                    goto nextimage
                }
            }
            if (modeflag == 3 || modeflag == 4) {
                if (imgets.value != "2" && l_extensions > 0) {
                    glogprint( l_logfile, "mistack", "status", type="error",
                        errno=123, str="Extension "//l_extension//" has "//\
                        imgets.value//" dimensions.  It should be 2.", verbose+)
                    status=status+1
                    goto nextimage
                }
            }
            l_extensions=l_extensions+1
        }

        j=l_extensions-1
        glogprint( l_logfile, "mistack", "engineering", type="string",
            str="Number of extensions is "//j, verbose=l_verbose)

        if (j < 1) {
            glogprint (l_logfile, "mistack", "status", type="error", errno=123,
                str="No data extensions in file "//in[i]//".", verbose+)
            status = status + 1
            goto nextimage
        }

        # WARNING: modeflag can change within the if (modeflat==?) blocks
        #        therefore, the if-else-if structure SHOULD NOT be used.

        if (modeflag == 1) {
            if (2*int((j/2)) != j) {
                if (j != 1) {
                    glogprint( l_logfile, "mistack", "status", type="warning",
                        errno=123, str="Number of extensions for input \
                        file "//in[i]//" does not correspond to complete \
                        nodsets.  Removing last unmatched nod position.",
                        verbose+)
                    l_extensions=l_extensions-1
                    if (l_extensions == 0) {
                        glogprint (l_logfile, "mistack", "status", type="error",
                            errno=123, str="No useable data extensions in \
                            file "//in[i]//".", verbose+)
                        status=status + 1
                        goto nextimage
                    }
                } else
                    modeflag = 2
            }

            if (j != 1) {
                for (j=1; j < l_extensions; j=j+2) {
                    # tmp image used only within this loop
                    tmpfile2=mktemp("tmpfile2")

                    imgets(in[i]//"["//j//"]","BADNOD",>&"dev$null")
                    if (imgets.value != "0") {
                        imgets(in[i]//"["//j+1//"]","BADNOD",>& "dev$null")
                        if (imgets.value != "0")
                            glogprint( l_logfile, "mistack", "engineering",
                                type="warning", errno=0, str="Nods "//j//\
                                " and "//j+1//" are both bad.)",
                                verbose=l_verbose)
                        else
                            glogprint( l_logfile, "mistack", "engineering",
                                type="warning", errno=0, str="Nod "//j//\
                                " is bad.  (Omitting both Nod "//j//\
                                " and "//j+1//")",verbose=l_verbose)
                    } else {
                        imgets(in[i]//"["//j+1//"]","BADNOD",>& "dev$null")
                        if (imgets.value != "0") {
                            glogprint( l_logfile, "mistack", "engineering",
                                type="warning", errno=0, str="Nod "//j+1//\
                                " is bad.  (Omitting both Nod "//j//" and "//\
                                j+1//")",verbose=l_verbose)
                        } else {
                            # tmpfile1 is used later on through 'tmpon'.
                            # Do not delete!
                            tmpfile1=mktemp("tmpfile1")
                            imarith(in[i]//"["//j//"]"//"[*,*,"//framevalue//"].fits",
                                "+",in[i]//"["//j+1//"]"//"[*,*,"//framevalue//"].fits",
                                tmpfile2,ver-)
                            if (l_combine == "average") {
                                imarith(tmpfile2,"/","2.0",tmpfile1,ver-)
                            }
                            else {
                                imcopy(tmpfile2,tmpfile1,ver-, >& "dev$null")
                            }
                            print(tmpfile1, >> tmpon)
                        }
                    }

                    # delete tmp image
                    imdelete(tmpfile2, ver-, >& "dev$null")
                    
                } #end for-loop over extensions

                imcombine("@"//tmpon,tmpfinal,headers="",bpmasks="",rejmasks="",
                    nrejmasks="",expmasks="",sigmas=tmpvar,logfile="STDOUT",
                    combine=l_combine,reject="none",project=no,outtype="double",
                    outlimits="", weight="none", offsets="none", scale="none",
                    >& tmplog)
                glogprint( l_logfile, "mistack", "science", type="file",
                    str=tmplog, verbose=l_verbose )
                delete( tmplog, verify-, >& "dev$null" )

                wmef(tmpfinal,out[i],extname="SCI",phu=header,verb-,
                    >& "dev$null")
                if (l_fl_variance) {
                  imfunction(tmpvar,tmpvar,"square",ver-, >& "dev$null")
                  frag=substr(out[i],strlen(out[i])-4,strlen(out[i]))
                  if (frag != ".fits") {
                    fxinsert(tmpvar,out[i]//".fits[1]","",ver-)
                  }
                  else {
                    fxinsert(tmpvar//".fits",out[i]//"[1]","",ver-)
                  }
                }

                # Time stamps
                gemdate ()
                gemhedit (out[i]//"[0]", "MISTACK", gemdate.outdate,
                    "UT Time stamp for MISTACK", delete-)
                gemhedit (out[i]//"[0]", "GEM-TLM", gemdate.outdate,
                    "UT Last modification with GEMINI", delete-)
                gemhedit (out[i]//"[1]", "AXISLAB4", "", "", delete=yes)
                gemhedit (out[i]//"[1]", "AXISLAB3", "", "", delete=yes)
            }

        } # see warning above
        if (modeflag == 2) {
            if (j > 1) {
                glogprint( l_logfile, "mistack", "status", type="error",
                    errno=123, str="More than one extension ("//j//") in \
                    CHOP mode.",verbose+)
                status=status+1
                goto nextimage
            }

            imgets(in[i]//"[1]","BADNOD",>& "dev$null")
            if (imgets.value != "0") {
                glogprint( l_logfile, "mistack", "status", type="warning",
                    errno=0, str="The only NOD is marked as BAD.",verbose+)
                    status=status+1
                goto nextimage
            }

            # tmpfile1 is used later on through 'tmpon'
            # Do not delete!
            tmpfile1=mktemp("tmpfile1")
            imcopy(in[i]//"[1][*,*,"//framevalue//"].fits",tmpfile1,ver-)
            print(tmpfile1, >> tmpon)

            imcombine("@"//tmpon,tmpfinal,headers="",bpmasks="",rejmasks="",
                nrejmasks="",expmasks="",sigmas=tmpvar,logfile="STDOUT",
                combine=l_combine,reject="none",project=no,outtype="double",
                outlimits="", weight="none", offsets="none", scale="none",
                >& tmplog)
            glogprint( l_logfile, "mistack", "science", type="file",
                str=tmplog, verbose=l_verbose )
            delete (tmplog, verify-, >& "dev$null")

            wmef(tmpfinal,out[i],extname="SCI",phu=header,verb-,>& "dev$null")
            if (l_fl_variance) {
              imfunction(tmpvar,tmpvar,"square",ver-, >& "dev$null")
              fxinsert(tmpvar,out[i]//"[1]","",ver-)
            }

            # Time stamps
            gemdate ()
            gemhedit (out[i]//"[0]", "MISTACK", gemdate.outdate,
                "UT Time stamp for MISTACK", delete-)
            gemhedit (out[i]//"[0]", "GEM-TLM", gemdate.outdate,
                "UT Last modification with GEMINI", delete-)
            gemhedit (out[i]//"[1]", "AXISLAB4", "", "", delete=yes)
            gemhedit (out[i]//"[1]", "AXISLAB3", "", "", delete=yes)

        } # see warning above
        if (modeflag == 3) {
            inod=0

            if (2*int((j/2)) != j) {
                glogprint( l_logfile, "mistack", "status", type="warning",
                    errno=123, str="Number of extensions for input file "//\
                    in[i]//" does not correspond to complete nodsets.  \
                    Removing last unmatched nod position.", verbose+ )
                l_extensions=l_extensions-1
                if (l_extensions == 0) {
                    glogprint (l_logfile, "mistack", "status", type="error",
                        errno=123, str="no useable data extensions in file "//\
                        in[i]//".", verbose+)
                    status=status+1
                    goto nextimage
                }
            }

            for (j=1; j < l_extensions; j=j+2) {
                # tmp image used only within this loop
                tmpfile2=mktemp("tmpfile2")
                inod=inod+1

                imgets(in[i]//"["//j//"]","BADNOD",>& "dev$null")
                if (imgets.value != "0") {
                    imgets(in[i]//"["//j+1//"]","BADNOD",>& "dev$null")
                    if (imgets.value != "0")
                        glogprint( l_logfile, "mistack", "engineering",
                            type="warning", errno=0, str="Nods "//j//" and "//\
                            j+1//" are both bad.)",verbose=l_verbose)
                    else
                        glogprint( l_logfile, "mistack", "engineering",
                            type="warning", errno=0, str="Nod "//j//" is bad.  \
                            (Omitting both Nod "//j//" and "//j+1//")",
                            verbose=l_verbose)
                } else {
                    imgets(in[i]//"["//j+1//"]","BADNOD",>& "dev$null")
                    if (imgets.value != "0")
                        glogprint( l_logfile, "mistack", "engineering",
                            type="warning", errno=0, str="Nod "//j+1//" is \
                            bad.  (Omitting both Nod "//j//" and "//j+1//")",
                            verbose=l_verbose)
                    else {
                        # tmpfile1 is used later on through 'tmpon'
                        # Do not delete!
                        tmpfile1=mktemp("tmpfile1")
                        # In NOD mode the nod B positions are negative, so taking the straight sum 
                        # of the extension pairs (in ABBA nodding, which is assumed) actually returns 
                        # the difference image.  The "src" and "ref" options are the same in this 
                        # case, and are found by taking nod A - nod B pairs.
                        if (l_frametype == "dif") {
                            imarith(in[i]//"["//j//"]"//"[*,*].fits","+",
                                in[i]//"["//j+1//"]"//"[*,*].fits",tmpfile2,ver-)
                        } else {
                            if (inod == 2*(inod/2)) {
                                imarith(in[i]//"["//j+1//"]"//"[*,*].fits","-",
                                    in[i]//"["//j//"]"//"[*,*].fits",
                                    tmpfile2,ver-)
                            } else {
                                imarith(in[i]//"["//j//"]"//"[*,*].fits","-",
                                    in[i]//"["//j+1//"]"//"[*,*].fits",
                                    tmpfile2,ver-)
                            }
                        }
                        if (l_combine == "average")
                            imarith(tmpfile2,"/","2.0",tmpfile1,ver-)
                        else
                            imcopy(tmpfile2,tmpfile1,ver-, >& "dev$null")

                        print(tmpfile1, >> tmpon)
                    }
                }

                # Delete tmp image
                imdelete(tmpfile2, ver-, >& "dev$null")
            }

            imcombine("@"//tmpon,tmpfinal,headers="",bpmasks="",rejmasks="",
                nrejmasks="",expmasks="",sigmas=tmpvar,logfile="STDOUT",
                combine="average",reject="none",project=no,outtype="double",
                outlimits="", >& tmplog)
            glogprint( l_logfile, "mistack", "science", type="file",
                str=tmplog, verbose=l_verbose )
            delete (tmplog, verify-, >& "dev$null")

            wmef(tmpfinal,out[i],extname="SCI",phu=header,verb-,>& "dev$null")
            if (l_fl_variance) {
              imfunction(tmpvar,tmpvar,"square",ver-, >& "dev$null")
              fxinsert(tmpvar,out[i]//"[1]","",ver-)
            }

            # Time stamps
            gemdate ()
            gemhedit (out[i]//"[0]", "MISTACK", gemdate.outdate,
                "UT Time stamp for MISTACK", delete-)
            gemhedit (out[i]//"[0]", "GEM-TLM", gemdate.outdate,
                "UT Last modification with GEMINI", delete-)
            gemhedit (out[i]//"[1]", "AXISLAB4", "", "", delete=yes)
            gemhedit (out[i]//"[1]", "AXISLAB3", "", "", delete=yes)

        } # see warning above
        if (modeflag == 4) {
            if (j > 1) {
                glogprint( l_logfile, "mistack", "status", type="error",
                    errno=123, str="More than one extension ("//j//") in \
                    STARE mode.",verbose+)
                status=status+1
                goto nextimage
            }

            imgets(in[i]//"[1]","BADNOD",>& "dev$null")
            if (imgets.value != "0") {
                glogprint( l_logfile, "mistack", "status", type="warning",
                    errno=0, str="The only NOD is marked as BAD.",verbose+)
                status=status+1
                goto nextimage
            }

            # tmpfile1 is used later on through 'tmpon'
            # Do not delete!
            tmpfile1=mktemp("tmpfile1")
            imcopy(in[i]//"[1][*,*].fits",tmpfile1,ver-)
            print(tmpfile1, >> tmpon)

            imcombine("@"//tmpon,tmpfinal,headers="",bpmasks="",rejmasks="",
                nrejmasks="",expmasks="",sigmas=tmpvar,logfile="STDOUT",
                combine="average",reject="none",project=no,outtype="double",
                outlimits="", >& tmplog)
            glogprint( l_logfile, "mistack", "science", type="file",
                str=tmplog, verbose=l_verbose )
            delete (tmplog, verify-, >& "dev$null")

            wmef(tmpfinal,out[i],extname="SCI",phu=header,verb-,>& "dev$null")
            if (l_fl_variance) {
              imfunction(tmpvar,tmpvar,"square",ver-, >& "dev$null")
              fxinsert(tmpvar,out[i]//"[1]","",ver-)
            }

            # Time stamps
            gemdate ()
            gemhedit (out[i]//"[0]", "MISTACK", gemdate.outdate,
                "UT Time stamp for MISTACK", delete-)
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
          gemhedit(out[i]//"[1]","WCSAXES",2,"Number of WCS axes in the image")
          imgets(header,"CTYPE1")
# Bug fix, the value is "0" after the first "gemhedit" call.
          gemhedit(out[i]//"[1]","CTYPE1",imgets.value,"R.A. in tangent plane projection")
        }
        imgets(header,"CRPIX1")
        if (imgets.value != "") gemhedit(out[i]//"[1]","CRPIX1",imgets.value,"Ref pix of axis 1")
        imgets(header,"CRVAL1")
        if (imgets.value != "") gemhedit(out[i]//"[1]","CRVAL1",imgets.value,"RA at Ref pix in decimal degrees")
        imgets(header,"CTYPE2")
        if (imgets.value != "") gemhedit(out[i]//"[1]","CTYPE2",imgets.value,"DEC. in tangent plane projection")
        imgets(header,"CRPIX2")
        if (imgets.value != "") gemhedit(out[i]//"[1]","CRPIX2",imgets.value,"Ref pix of axis 2")
        imgets(header,"CRVAL2")
        if (imgets.value != "") gemhedit(out[i]//"[1]","CRVAL2",imgets.value,"DEC at Ref pix in decimal degrees")
        imgets(header,"CD1_1")
        if (imgets.value != "") gemhedit(out[i]//"[1]","CD1_1",imgets.value,"WCS matrix element 1 1")
        imgets(header,"CD1_2")
        if (imgets.value != "") gemhedit(out[i]//"[1]","CD1_2",imgets.value,"WCS matrix element 1 2")
        imgets(header,"CD2_1")
        if (imgets.value != "") gemhedit(out[i]//"[1]","CD2_1",imgets.value,"WCS matrix element 2 1")
        imgets(header,"CD2_2")
        if (imgets.value != "") gemhedit(out[i]//"[1]","CD2_2",imgets.value,"WCS matrix element 2 2")
        imgets(header,"RADECSYS")
        if (imgets.value != "") gemhedit(out[i]//"[1]","RADECSYS",imgets.value,"R.A./DEC. coordinate system reference")

        # jump to here if there is a problem

nextimage:
        i=i+1
        delete(tmpfinal//".fits",ver-,>& "dev$null")
        delete(tmpvar//".fits",ver-,>& "dev$null")
        imdelete("@"//tmpon, ver-, >& "dev$null")
        delete(tmpon,ver-, >& "dev$null")
    }

clean:
    scanfile=""
    delete(tmpfile,ver-, >& "dev$null")
    delete("tmpin*",ver-, >& "dev$null")
    delete("tmpfinal*",ver-, >& "dev$null")
    delete("tmpfile*",ver-, >& "dev$null")
    delete(tmplog, ver-, >& "dev$null")

    if (status==0)
        glogclose( l_logfile, "mistack", fl_success+, verbose=l_verbose )
    else
        glogclose( l_logfile, "mistack", fl_success-, verbose=l_verbose )

exitnow:
    ;

end

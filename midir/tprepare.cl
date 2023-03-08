# Copyright(c) 2003-2011 Association of Universities for Research in Astronomy, Inc.

procedure tprepare(inimages)

#
# This routine converts a raw T-ReCS file to the standard "midir" format, 
# combining savesets for each nod.
#
# Version:  Sept  5, 2003  KV started writing the script from "trgather.cl"
#           Sept  9, 2003  KV added "outpref" based on TBs example
#           Sept 26, 2003  KV register option added, logging changed a bit
#           Oct   3, 2003  KV syntax changes to the parameters, changed to 
#                             using imstack
#           Oct  29, 2003  KL IRAF 2.12 - new/modified parameters
#                             imcombine: headers,bpmasks,expmask,outlimits
#                               rejmask->rejmasks,plfile->nrejmasks
#                               hedit: addonly
#           Nov  12, 2003  KV delete "AXISLAB4" from extension headers
#                             add WCS values to extension headers
#           Dec  02, 2003  KV took out the WCS values since these crash 
#                             imcombine for subsequent steps (i.e. mistack)
#           Aug  19, 2004  KV added "combine" option--average or sum, 
#                             and explicitly set imcombine options 
#                             weight, scale, offsets to "none"
#           Aug  23, 2004  KV made changes to allow "rescuing" of some 
#                             incomplete observations...fl_check added 

char    inimages    {prompt="Input T-ReCS image(s)"}                 # OLDP-1-input-primary-single-prefix=t
char    outimages   {"",prompt="Output image(s)"}                       # OLDP-1-output
char    outpref     {"t",prompt="Prefix for output image name(s)"}      # OLDP-4
char    rawpath     {"",prompt="Path for input raw images"}             # OLDP-4
char    stackoption {"stack",prompt="Image combining option: stack|register"}   # OLDP-2
char    combine     {"average",prompt="Combining images by average|sum"}    # OLDP-2
bool    fl_check    {yes,prompt="Check observation during processing?"} # OLDP-2
char    logfile     {"",prompt="Logfile name"}                          # OLDP-1
bool    verbose     {yes,prompt="Verbose logging yes/no?"}              # OLDP-4
int     status      {0,prompt="Exit error status: (0=good, >0=bad)"}    # OLDP-4
struct  *scanfile   {"",prompt="Internal use only"}                     # OLDP-4

begin

    char    l_inputimages,l_outputimages,l_filename,l_logfile
    char    l_rawpath,l_prefix,l_stack, l_combine
    char    in[100],out[100],header,exheader
    char    tmpin,tmpout,tmpon,tmpoff,tmpfile,tmpfinal,tmpwork
    char    tmplist,tmpshift,tmpregister,tmplog
    char    refimage
    char    paramstr, instrument
    int     i,j,k,l,i_nodset,i_saveset,itotal,n_sig,lext
    int     n_ref,nimages,maximages,noutputimages,i_frames,i_extensions
    int     source, reference, nodcount, nbadsets, badsetnumber[100], badflag
    int     modeflag, nwarnings
    int     filetype,iregister
    bool    l_verbose
    real    exptime,norm,ave1,ave2,diff1,airmass
    bool    l_check

    tmpfile=mktemp("tmpin")
    tmpin=mktemp("tmpin")
    tmplist=mktemp("tmplist")
    tmpshift=mktemp("tmpshift")
    tmpregister=mktemp("tmpregister")
    tmplog=mktemp("tmplog")

    l_inputimages=inimages
    l_outputimages=outimages
    l_prefix=outpref
    l_rawpath=rawpath
    l_logfile=logfile
    l_verbose=verbose
    l_stack=stackoption
    l_combine=combine
    l_check=fl_check

    status=0
    nwarnings = 0
    
    cache ("gemdate")

    # Create the list of parameter/value pairs.  One pair per line.
    # All lines combined into one string. Line delimiter is '\n'.
    paramstr =  "inimages       = "//inimages.p_value//"\n"
    paramstr += "outimages      = "//outimages.p_value//"\n"
    paramstr += "outpref        = "//outpref.p_value//"\n"
    paramstr += "rawpath        = "//rawpath.p_value//"\n"
    paramstr += "stackoption    = "//stackoption.p_value//"\n"
    paramstr += "combine        = "//combine.p_value//"\n"
    paramstr += "fl_check       = "//fl_check.p_value//"\n"
    paramstr += "logfile        = "//logfile.p_value//"\n"
    paramstr += "verbose        = "//verbose.p_value

    # Assign a logfile name if not given.  Open logfile and start log.
    # Write parameter/value pairs ("paramstr") to log.
    gloginit ( l_logfile, "tprepare", "midir", paramstr, fl_append+,
        verbose=l_verbose )
    if (gloginit.status != 0) {
        status = gloginit.status
        goto exitnow
    }
    l_logfile = gloginit.logfile.p_value

    # Check the rawpath name for a final /
    if (substr(l_rawpath,(strlen(l_rawpath)),(strlen(l_rawpath))) != "/")
        l_rawpath=l_rawpath//"/"
    if (l_rawpath=="/" || l_rawpath==" ")
        l_rawpath=""

    if (l_combine != "average" && l_combine != "sum") {
        glogprint (l_logfile, "tprepare", "status", type="error", errno=121,
            str="Bad combine parameter ("//l_combine//")", verbose+)
        status = 1
        goto exit
    }

    nimages=0
    maximages=100

    # Count the number of input images
    # First, generate the file list if needed

    if (stridx("*",l_inputimages) > 0) {
        files(l_rawpath//l_inputimages, > tmpfile)
        l_inputimages="@"//tmpfile
    }

    # The following is copied from the "gprepare.cl" script...but it does 
    # not work properly on aguila.  It complains about not being able to 
    # change the type from string.  
    #
    # Parse wildcard, list file, and comma-separated lists of input images
    #
    if (substr(l_inputimages,1,1)=="@")
        scanfile=substr(l_inputimages,2,strlen(l_inputimages))
    else {
        if (stridx(",",l_inputimages)==0) 
            files(l_inputimages, > tmpfile)
        else {
            j=9999
            while (j!=0) {
                j=stridx(",",l_inputimages)
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

        if (!imaccess(l_rawpath//l_filename))
            glogprint( l_logfile, "tprepare", "status", type="error", 
                errno=101, str="Input image "//l_rawpath//l_filename//" was \
                not found.", verbose+)
        else {
            nimages=nimages+1
            if (nimages > maximages) {
                glogprint( l_logfile, "tprepare", "status", type="error",
                    errno=121, str="Maximum number of input images \
                    ["//str(maximages)//"] has been exceeded.",verbose+)
                status=1
                goto exit
            }
            in[nimages]=l_rawpath//l_filename
            out[nimages]=l_filename
            j=stridx("/",out[nimages])
            if (j > 0) {
                for (k=1; k < 100 && j > 0; k+=1) {
                    out[nimages]=substr(out[nimages],j+1,strlen(out[nimages]))
                    j=stridx("/",out[nimages])
                }
            }
        }
    } #end while-loop

    scanfile=""
    delete(tmpfile//","//tmpin,verify-,>& "dev$null")

    if (nimages == 0) {
        glogprint( l_logfile, "tprepare", "status", type="error", errno=121,
            str="No input images were defined.",verbose+)
        status=1
        goto exit
    } else
        glogprint( l_logfile, "tprepare", "status", type="string",
            str="Processing "//str(nimages)//" image(s).",verbose=l_verbose)

    # Now, do the same counting for the output file
    tmpfile=mktemp("tmpfile")

    noutputimages=0
    if (l_outputimages != "" && l_outputimages != " ") {
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
            if (substr(l_filename,strlen(l_filename)-4,strlen(l_filename)) == ".fits")
                l_filename=substr(l_filename,1,strlen(l_filename)-5)

            noutputimages=noutputimages+1
            if (noutputimages > maximages) {
                glogprint( l_logfile, "tprepare", "status", type="error",
                    errno=121, str="Maximum number of output images "//\
                    str(maximages)//" exceeded.",verbose+)
                status=1
                goto exit
            }
            out[noutputimages]=l_filename
            if (imaccess(out[noutputimages])) {
                glogprint( l_logfile, "tprepare", "status", type="error",
                    errno=102, str="Output image "//l_filename//" already \
                    exists.",verbose+)
                status=1
                goto exit
            }
        }
        if (noutputimages != nimages) {
            glogprint( l_logfile, "tprepare", "status", type="error", 
                errno=121, str="Different number of input ("//str(nimages)//"\
                ) and output ("//str(noutputimages)//" image names have been \
                specified.", verbose+)
            status=1
            goto exit
        }

        scanfile=""
        delete(tmpfile,verify-, >& "dev$null")
    } else {
        if (l_prefix == "" || l_prefix == " ") 
            l_prefix="t"
        for (i=1; i <= nimages; i+=1) {
            out[i]=l_prefix//out[i]
            if (imaccess(out[i])) {
                glogprint( l_logfile, "tprepare", "status", type="error",
                    errno=102, str="Output image "//out[i]//" already exists.",
                    verbose+)
                status=1
                goto exit
            }
        }
    }

    #
    # The main loop: accumulate source and reference images
    #

    i=1
    while (i <= nimages) {
        glogprint( l_logfile, "tprepare", "visual", type="visual",
            vistype="empty", verbose=l_verbose )
        imgets(in[i]//"[0]","TPREPARE", >& "dev$null")
        if (imgets.value != "0") {
            glogprint( l_logfile, "tprepare", "status", type="warning",
                errno=0, str="File "//in[i]//" has already been prepared.",
                verbose=l_verbose)
            goto nextimage
        }
        
        ## for CanariCam, 'stackoption="register"' is not supported
        ## reset to "stack" and issue warning.
        if (l_stack == "register") {
            instrument = ""
            hselect (in[i]//"[0]", "INSTRUME", yes) | scan (instrument)
            if (instrument == "CanariCam") {
                l_stack = "stack"
                glogprint (l_logfile, "tprepare", "status", type="warning",
                    errno=121, str="stackoption='register' is not supported \
                    for CanariCam. RESETING to STACK.", verbose=l_verbose)
                nwarnings += 1
            }
        }

        glogprint( l_logfile, "tprepare", "task", type="string",
            str="Input image "//in[i]//" output image "//out[i]//".",
            verbose=l_verbose)

        # check the file structure

        # The following allows the routine to work with older raw data files.  
        # It should give an error if the value is not "complete" and these 
        # should be handled by another routine once I figure out what need to 
        # be done.
        #
        # For handling older files, I will not check this value.
        #

        imgets(in[i]//"[0]","COMPSTAT", >& "dev$null")
        if (imgets.value == "0") {
            if (l_check) {
                glogprint (l_logfile, "tprepare", "status", type="error",
                    errno=131, str="Image "//in[i]//" does not have the \
                    COMPSTAT header value.", verbose+)
                status = 1
                goto nextimage
            } else
                glogprint( l_logfile, "tprepare", "status", type="warning", 
                    errno=131, str="Image "//in[i]//" does not have the \
                    COMPSTAT header value.",verbose=l_verbose)
        } else {
            if (imgets.value != "COMPLETE") {
                if (l_check) {
                    glogprint (l_logfile, "tprepare", "status", type="error",
                        errno=132, str="Image "//in[i]//" has COMPSTAT = "//\
                        imgets.value//".", verbose+)
                    status=1
                    goto nextimage
                } else
                    glogprint( l_logfile, "tprepare", "status", type="warning",
                        errno=132, str="Image "//in[i]//" has COMPSTAT = "//\
                        imgets.value//".",verbose=l_verbose)
            }
        }

        glogprint( l_logfile, "tprepare", "status", type="fork", 
            fork="forward", child="tcheckstructure", verbose=l_verbose )

        tcheckstructure(in[i],logfile=l_logfile,verbose=l_verbose)

        glogprint( l_logfile, "tprepare", "status", type="fork",
            fork="backward", child="tcheckstructure", verbose=l_verbose )
        j=tcheckstructure.status
        if (l_check && j != 0) {
            glogprint( l_logfile, "tprepare", "status", type="error", 
                errno=123, str="Image "//in[i]//" does not have the expected \
                structure.", verbose+)
            status=1
            goto nextimage
        }

        modeflag=tcheckstructure.modeflag
        if (modeflag > 4) {
            modeflag=modeflag-10
            filetype=2
        } else
            filetype=1

        if (filetype != 1) {
            glogprint( l_logfile, "tprepare", "status", type="error", 
                errno=121, str="Image "//in[i]//" appears to be a 'prepared' \
                T-ReCS file.", verbose+)
            status=1
            goto nextimage
        }

        # check the primary FITS header
        header=in[i]//"[0]"

        if (modeflag == 1 || modeflag == 3) {
            imgets(header,"NNODS", >& "dev$null")
            if (imgets.value != "2") {
                glogprint( l_logfile, "tprepare", "status", type="error",
                    errno=123, str="Image "//in[i]//" has "//imgets.value//\
                    " nod phases (must be 2).",verbose+)
                status=1
                goto nextimage
            }
            imgets(header,"NNODSETS", >& "dev$null")
            i_nodset=int(imgets.value)
            imgets(header,"SAVESETS", >& "dev$null")
            i_saveset=int(imgets.value)
            if (modeflag == 1)
                i_frames=i_saveset*4*i_nodset
            else
                i_frames=i_saveset*2*i_nodset
        }
        if (modeflag == 4 || modeflag == 2) {
            imgets(header,"SAVESETS", >& "dev$null")
            i_saveset=int(imgets.value)
            if (modeflag == 2)
                i_frames=i_saveset*2
            else
                i_frames=i_saveset*2
        }

        if (i_frames == 0) {
            glogprint( l_logfile, "tprepare", "status", type="error", 
                errno=123, str="Either zero nod sets or zero save sets in \
                image "//in[i]//".",verbose+)
            status=1
            goto nextimage
        }

        imgets(header,"NEXTEND", >& "dev$null")
        if (imgets.value != "0")
            i_extensions=int(imgets.value)
        else {
            if (modeflag == 1 || modeflag == 3)
                i_extensions=i_nodset*2
            else
                i_extensions=1
        }

        # Actually count the number of extensions--

        lext=0  
        for (j=1; j <= i_extensions; j=j+1) {
            if (imaccess(in[i]//"["//j//"]"))
                lext=lext+1
        }

        if (l_check) {
            if (lext != i_extensions) {
                glogprint (l_logfile, "tprepare", "status", type="error",
                    errno=123, str="Incorrect number of extensions. Found "//\
                    str(lext)//" in image "//in[i]//": should be "\
                    //str(i_extensions)//".", verbose+)
                status = 1
                goto nextimage
            }
        }

        for (j=1; j <= lext; j=j+1) {
            # Create tmp FITS file names used within this loop
            tmpout=mktemp("tmpout")
            tmpon=mktemp("tmpon")
            tmpoff=mktemp("tmpoff")
            tmpfinal=mktemp("tmpfinal")

            exheader=in[i]//"["//j//"]"
            if (modeflag == 1 || modeflag == 3) {
                imgets(exheader,"NOD", >& "dev$null")
                if (imgets.value == "A") {
                    source=1
                    reference=2
                } else {
                    if (imgets.value == "B") {
                        source=2
                        reference=1
                    } else {
                        glogprint( l_logfile, "tprepare", "status",
                            type="error", errno=123, str="Bad nod position \
                            ("//imgets.value//") in image "//in[i]//".",
                            verbose+)
                        status=1
                        goto nextimage
                    }
                }
                imgets(exheader,"NODSET", >& "dev$null")
                nodcount=int(imgets.value)
                if (imgets.value == "0") {
                    glogprint( l_logfile, "tprepare", "status", type="error",
                        errno=121, str="Nod number recorded as "//\
                        imgets.value//".",verbose+)
                    status=1
                    goto nextimage
                } else {
                    i_nodset=int(imgets.value)
                    if (i_nodset < 1) {
                        gloginit( l_logfile, "tprepare", "status", 
                            type="error", errno=121, str="Nod number \
                            recorded as "//imgets.value//".",verbose+)
                        status=1
                        goto nextimage
                    }
                }
            } else {
                # For chop or stare observations, assume it is in nod 
                # position A unless the header says otherwise.
                #
                source=1
                reference=2
                imgets(exheader,"NOD", >& "dev$null")
                if (imgets.value == "B") {
                    source=2
                    reference=1
                }
            }
            #
            # Check for bad savesets in the nodset...
            #
            imgets(exheader,"NBADSET", >& "dev$null")
            nbadsets=int(imgets.value)
            if (nbadsets < 0 || nbadsets >= i_saveset) {
                glogprint( l_logfile, "tprepare", "status", type="error",
                    errno=132, str="Header records "//str(nbadsets)//\
                    " (should be less than "//str(i_saveset)//").",verbose+)
                goto nextimage
            }
            if (nbadsets > 0) {
                for (k=1; k <= nbadsets; k=k+1) {
                    if (k < 10)
                        imgets(exheader,"BADSET0"//str(k), >& "dev$null")
                    else
                        imgets(exheader,"BADSET"//str(k), >& "dev$null")
                    badsetnumber[k]=int(imgets.value)
                    if (badsetnumber[k] < 1 || badsetnumber[k] > i_saveset) {
                        glogprint( l_logfile, "tprepare", "status",
                            type="error", errno=132, str="Header records bad \
                            save set "//imgets.value//" (range 1 to "//\
                            str(i_savesets)//").",verbose+)
                        goto nextimage
                    }
                }
            }

            iregister=0
            for (k=1; k <= i_saveset; k=k+1) {
                badflag=0
                for (l=1; l <= nbadsets; l=l+1) {
                    if (badsetnumber[l] == k) {
                        badflag=1
                        glogprint( l_logfile, "tprepare", "engineering",
                            type="string", str="Header records bad saveset "//\
                            str(k)//" in nod set; will be omitted from the \
                            combined frame.",verbose=l_verbose)
                    }
                }
                if (badflag == 0) {
                    print(exheader//"[*,*,"//str(source)//","//str(k)//"]",
                        >> tmpin)
                    print(exheader//"[*,*,"//str(reference)//","//str(k)//"]",
                        >> tmpfile)
                    unlearn("imexpr")
                    if (modeflag != 4) {
                        # Create tmp FITS file names used within this loop
                        # DO NOT DELETE tmpwork at the end of this for loop
                        # as the image are used throught @tmplist
                        tmpwork=mktemp("tmpwork")

                        if (iregister == 0) {
                            refimage=tmpwork
                            iregister=1
                        }
                        print(tmpwork, >> tmplist)
                        imexpr("(a-b)",tmpwork,exheader//"[*,*,"//\
                            str(source)//","//str(k)//"]",exheader//\
                            "[*,*,"//str(reference)//","//str(k)//"]",verbose-)
                        print("tmp"//out[i]//"shift"//str(j)//str(k),
                            >> tmpshift)
                    }
                }
            }

            delete( tmplog, verify-, >& "dev$null" )
            imcombine("@"//tmpin,tmpon,headers="",bpmasks="",rejmasks="",
                nrejmasks="",expmasks="",sigmas="",logfile="STDOUT",
                combine=l_combine,reject="none",project=no,outtype="double",
                outlimits="", weight="none", offsets="none", scale="none",
                >& tmplog)
            glogprint( l_logfile, "tprepare", "science", type="file",
                str=tmplog, verbose=l_verbose )

            if (modeflag != 4) {
                delete( tmplog, verify-, >& "dev$null" )
                imcombine("@"//tmpfile, tmpoff, headers="", bpmasks="",
                    rejmasks="", nrejmasks="", expmasks="", sigmas="",
                    logfile="STDOUT", combine=l_combine, reject="none",
                    project=no,outtype="double", outlimits="", weight="none",
                    offsets="none", scale="none", >& tmplog)
                glogprint( l_logfile, "tprepare", "science", type="file",
                    str=tmplog, verbose=l_verbose )

                if (l_stack == "register" && iregister != 0) {
                    images.immatch.xregister("@"//tmplist, refimage, "[*,*]",
                        tmpregister, output="@"//tmpshift, interactive=no, 
                        verbose=no)

                    delete( tmplog, verify-, >& "dev$null" )
                    imcombine("@"//tmpshift,tmpout,headers="",bpmasks="",
                        rejmasks="",nrejmasks="",expmasks="",sigmas="",
                        logfile="STDOUT",combine=l_combine,reject="sigclip",
                        project=no,outtype="double",outlimits="",
                        weight="none", offsets="none", scale="none", >& tmplog)
                    glogprint( l_logfile, "tprepare", "science", type="file",
                        str=tmplog, verbose=l_verbose )
                } else {
                    delete( tmplog, verify-, >& "dev$null" )
                    imcombine("@"//tmplist, tmpout, headers="", bpmasks="",
                        rejmasks="", nrejmasks="", expmasks="", sigmas="",
                        logfile="STDOUT", combine=l_combine, reject="none",
                        project=no, outtype="double", outlimits="", 
                        weight="none", offsets="none", scale="none", >& tmplog)
                    glogprint( l_logfile, "tprepare", "science", type="file",
                        str=tmplog, verbose=l_verbose )
                }
            }
            count(tmpin) | scan(itotal)
            delete(tmpin//","//tmpfile//","//tmplist,verify-, >& "dev$null")
            if (l_stack == "register") {
                delete(tmpshift//","//tmpregister,verify-, >& "dev$null")
                delete("tmp*shift*.fits",verify-, >& "dev$null")
            }
            delete("tmpshift*",verify-, >& "dev$null")
            itotal=itotal/2

            # Write out the combined image(s)
            if (modeflag == 4)
                imcopy(tmpon,tmpfinal,verbose-)
            else
                imstack(tmpon//","//tmpoff//","//tmpout,tmpfinal)

            # If j is one, write the primary header to the output file.  
            # Also write "NEXTEND" if it is not there already...it is useful 
            # for later tasks.
            #
            # Write the image extensions with "fxinsert".
            #
            if (j == 1) {
                imcopy(header,out[i], verbose-)
                gemhedit (out[i]//".fits", "NEXTEND", i_extensions, "",
                    delete-)
            }
            l=j-1
            fxinsert(tmpfinal//".fits", out[i]//".fits["//str(l)//"]", "",
                verbose-)

            if (modeflag == 1 || modeflag == 3) {
                if (source == 1)
                    gemhedit (out[i]//".fits["//str(j)//"]", "NOD", "A", "",
                        delete-)
                else
                    gemhedit (out[i]//".fits["//str(j)//"]", "NOD", "B", "", 
                        delete-)
                gemhedit (out[i]//".fits["//str(j)//"]", "NODSET", nodcount, 
                    "", delete-)
            }
            gemhedit (out[i]//".fits["//str(j)//"]", "SAVESETS", i_saveset,
                "", delete-)
            imgets(exheader,"UTSTART", >& "dev$null")
            if (imgets.value != "0")
                gemhedit (out[i]//".fits["//str(j)//"]", "UTSTART",
                    imgets.value, "", delete-)
            imgets(exheader,"UTEND", >& "dev$null")
            if (imgets.value != "0")
                gemhedit (out[i]//".fits["//str(j)//"]", "UTEND",
                    imgets.value, "", delete-)
            imgets(exheader,"AMSTART", >& "dev$null")
            if (imgets.value != "0") {
                airmass=real(imgets.value)
                gemhedit (out[i]//".fits["//str(j)//"]", "AMSTART",
                    airmass, "", delete-)
            }

            gemhedit (out[i]//".fits["//str(j)//"]", "AXISLAB4", "", "",
                delete=yes)
            # Write extension number and name
            gemhedit (out[i]//"["//j//"]", "EXTVER", j, "Extension version")
            gemhedit (out[i]//"["//j//"]", "EXTNAME", "SCI", "Extension name")

            if (modeflag == 2 || modeflag == 4)
                glogprint( l_logfile, "tprepare", "engineering", type="string",
                    str="Image written to output file extension "//str(j)//".",
                    verbose=l_verbose )
            else {
                if (source == 1)
                    glogprint( l_logfile, "tprepare", "engineering", 
                        type="string", str="Nod A image for set "//\
                        str(nodcount)//" written to output file extension "//\
                        str(j)//".",verbose=l_verbose)
                else
                    glogprint( l_logfile, "tprepare", "engineering",
                        type="string", str="Nod B image for set "//\
                        str(nodcount)//" written to output file extension "//\
                        str(j)//".",verbose=l_verbose)
            }

            # Delete temporary files and images
            imdelete(tmpon//","//tmpoff//","//tmpout//","//tmpfinal//","//\
            tmpshift,verify-, >& "dev$null")
            delete("tmpwork*.fits",verify-, >& "dev$null")

        } #end for-loop on extensions

        # The following is a patch to add spectroscopy keywords that Tracy 
        # needs for her routines

        imgets(in[i]//"[0]","GRATING", >& "dev$null")
        if (imgets.value == "HiRes-10") {
            gemhedit(out[i]//"[0]","DISPAXIS",1,"Dispersion direction")
        }
        if (imgets.value == "LowRes-10") {
            gemhedit(out[i]//"[0]","GRATWAVE",10.5,
                "Grating central position (microns).")
            gemhedit(out[i]//"[0]","GRATDISP",0.022,
                "Approximate grating dispersion (microns).")
            gemhedit(out[i]//"[0]","DISPAXIS",1,"Dispersion direction")
        }
        if (imgets.value == "LowRes-20") {
            gemhedit(out[i]//"[0]","GRATWAVE",20.0,
                "Grating central wavelength (microns).")
            gemhedit(out[i]//"[0]","GRATDISP",0.033,
                "Approximate grating dispersion (microns).")
            gemhedit(out[i]//"[0]","DISPAXIS",1,"Dispersion direction")
        }


        # Time stamp the primary header
        #
        gemdate ()
        gemhedit (out[i]//"[0]", "GEM-TLM", gemdate.outdate,
            "UT Last modification with GEMINI IRAF", delete-)
        gemhedit (out[i]//"[0]", "TPREPARE", gemdate.outdate,
            "UT Time stamp for TPREPARE")
        if (l_check)
            gemhedit(out[i]//"[0]","FILECHK","TRUE ",
                "TPREPARE was run with file checks.")
        else
            gemhedit(out[i]//"[0]","FILECHK","FALSE",
                "TPREPARE was run with NO file checks.")
        gemhedit (out[i]//"[0]", "PREPARE", gemdate.outdate,
            "UT Time stamp for PREPARE")

nextimage:
        # jump to here if there is a problem
        i=i+1
    }

exit:
    delete (tmplog, verify-, >& "dev$null")
    if (nwarnings > 0) {
        glogprint (l_logfile, "tprepare", "visual", type="visual", 
            vistype="empty", verbose=l_verbose)
        glogprint (l_logfile, "tprepare", "status", type="warning", errno=0,
            str="There were "//str(nwarnings)//" warning(s).", 
            verbose=l_verbose)
        glogprint (l_logfile, "tprepare", "status", type="warning", errno=0,
            str="Please review the logs.", verbose=l_verbose)
        glogprint (l_logfile, "tprepare", "visual", type="visual", 
            vistype="empty", verbose=l_verbose)
    }
    if (status == 0)
        glogclose (l_logfile, "tprepare", fl_success+, verbose=l_verbose)
    else
        glogclose (l_logfile, "tprepare", fl_success-, verbose=l_verbose)

exitnow:
    ;

end

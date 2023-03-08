# Copyright(c) 2003-2011 Association of Universities for Research in Astronomy, Inc.

procedure tview (inimages)

#
# Procedure for interactively examining T-ReCS frames
#
# Version:  June  8, 2003  KV routine operational as "trrawview.cl"
#           Sept  9, 2003  KV added logfie entries, copying of input file,
#                             changes to BADSET keyword logic
#           Sept 12, 2003  KV code to handle all modes of T-ReCS
#           Sept 26, 2003  KV syntax changes, add time stamp
#           Oct   3, 2003  KV small syntax changes
#           Oct  29, 2003  KL IRAF 2.12 - new parameter
#                             hedit: addonly
#

char    inimages        {prompt="Input T-ReCS image(s)"}
char    outimages       {"",prompt="Output image(s)"}
char    outpref         {"v",prompt="Prefix for output image name(s)"}
char    rawpath         {"",prompt="Path for input raw images"}
char    type            {"dif",prompt="Type of frame (src|ref|dif)"}
real    delay           {0.,prompt="Update delay in seconds"}
bool    fl_inter        {yes,prompt="Interactive screening of frames?"}
bool    fl_disp_fill    {no,prompt="Fill display"}
bool    fl_label_display    {yes,prompt="Label display"}
char    colour          {"black",prompt="Label colour: black|white|red|green|blue|yellow|cyan|magenta|orange"}
bool    fl_use_imexam   {no,prompt="Use imexam to look at images before screening?"}
bool    fl_sh_change    {no,prompt="Show changes to image headers"}
bool    fl_delete       {yes,prompt="Delete unchanged copy images"}
real    z1              {0.,prompt="Minimum level to be displayed"}
real    z2              {0., prompt="Maximum level to be displayed"}
bool    zscale          {yes,prompt="Display graylevels near the median"}
bool    zrange          {yes,prompt="Display full image intensity range"}
char    ztrans          {"linear",prompt="Greylevel transformation (linear|log|none|user)"}
char    logfile         {"",prompt="Logfile name"}
bool    verbose         {yes,prompt="Verbose"}
int     status          {0,prompt="Exit status: (0=good, >0=bad)"}
struct  *scanfile       {"",prompt="Internal use only"}

begin

    char    l_image,l_colour,inputstring,l_inputimages,l_prefix
    char    l_outputimages,l_rawpath,l_logfile,l_filename,l_type
    real    l_delay
    bool    l_verbose,l_dispfill,l_label,l_manual,l_exam,l_sh_change,l_delete

    char    l_ztrans
    bool    l_zscale
    bool    l_zrange
    real    l_z1,l_z2

    char    fname,tref1,tsig1,tdif1,tmpfile,cursinput,phu
    char    paramstr
    int     nsavesets,nnodsets,nextns,na,nnods,changefile
    int     nbadsets,badsetnumber[100],nbadsetsorig,badoff,badflag,badtoggle
    int     beamval
    int     ncol,npars,maxch,wcs
    int     maxbad
    real    xpos,ypos

    char    in[100],out[100]
    char    in1,out1
    int     nimages,modeflag,maximages,noutputimages
    int     i,j,k,l

    cache ("gemdate")
    
    l_inputimages=inimages
    l_outputimages=outimages
    l_rawpath=rawpath
    l_prefix=outpref
    l_logfile=logfile
    l_delay=delay
    l_colour=colour
    l_verbose=verbose
    l_dispfill=fl_disp_fill
    l_label=fl_label_display
    l_manual=fl_inter
    l_exam=fl_use_imexam
    l_sh_change=fl_sh_change
    l_delete=fl_delete
    l_z1=z1
    l_z2=z2
    l_zscale=zscale
    l_ztrans=ztrans
    l_zrange=zrange
    l_type=type
    cursinput=""
    
    status=0
    maxbad=100

    # Create the list of parameter/value pairs.  One pair per line.
    # All lines combined into one string.  Line delimiter is '\n'.
    paramstr =  "inimages         = "//inimages.p_value//"\n"
    paramstr += "outimages        = "//outimages.p_value//"\n"
    paramstr += "outpref          = "//outpref.p_value//"\n"
    paramstr += "rawpath          = "//rawpath.p_value//"\n"
    paramstr += "type             = "//type.p_value//"\n"
    paramstr += "delay            = "//delay.p_value//"\n"
    paramstr += "fl_inter         = "//fl_inter.p_value//"\n"
    paramstr += "fl_disp_fill     = "//fl_disp_fill.p_value//"\n"
    paramstr += "fl_label_display = "//fl_label_display.p_value//"\n"
    paramstr += "colour           = "//colour.p_value//"\n"
    paramstr += "fl_use_imexam    = "//fl_use_imexam.p_value//"\n"
    paramstr += "fl_sh_change     = "//fl_sh_change.p_value//"\n"
    paramstr += "fl_delete        = "//fl_delete.p_value//"\n"
    paramstr += "z1               = "//z1.p_value//"\n"
    paramstr += "z2               = "//z2.p_value//"\n"
    paramstr += "zscale           = "//zscale.p_value//"\n"
    paramstr += "zrange           = "//zrange.p_value//"\n"
    paramstr += "ztrans           = "//ztrans.p_value//"\n"
    paramstr += "logfile          = "//logfile.p_value//"\n"
    paramstr += "verbose          = "//verbose.p_value

    # Assign a logfile name if not given.  Open logfile and start log.
    # Write parameter/value pairs ("paramstr") to log.
    gloginit( l_logfile, "tview", "midir", paramstr, fl_append+,
    verbose=l_verbose )
    if (gloginit.status != 0) {
        status = gloginit.status
        goto exitnow
    }
    l_logfile = gloginit.logfile.p_value

    # Check the rawpath name for a final /
    if (substr(l_rawpath,(strlen(l_rawpath)),(strlen(l_rawpath))) != "/")
        l_rawpath=l_rawpath//"/"
    if(l_rawpath=="/" || l_rawpath==" ")
        l_rawpath=""

    nimages=0
    maximages=100

    tmpfile=mktemp("tmplist")

    # Count the number of input images
    # First, generate the file list if needed

    if (stridx("*",l_inputimages) > 0) {
        files(l_rawpath//l_inputimages, > tmpfile)
        l_inputimages="@"//tmpfile
    }

    if (substr(l_inputimages,1,1)=="@")
        scanfile=substr(l_inputimages,2,strlen(l_inputimages))
    else {
        if(stridx(",",l_inputimages)==0) 
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
            glogprint( l_logfile, "tview", "status", type="error", errno=101,
                str="Input image "//l_rawpath//l_filename//" was not found.",
                verbose+)
        else {
            nimages=nimages+1
            if (nimages > maximages) {
                glogprint( l_logfile, "tview", "status", type="error",
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
    } # end while-loop

    scanfile=""
    delete(tmpfile,verify-,>& "dev$null")

    if (nimages == 0) {
        glogprint( l_logfile, "tview", "status", type="error", errno=121,
            str="No input images were defined.",verbose+)
        status=1
        goto exit
    } else
        glogprint( l_logfile, "tview", "status", type="string",
            str="Processing "//str(nimages)//" image(s).",verbose=l_verbose)

    # Now, do the same counting for the output file
    tmpfile=mktemp("tmplist")
    noutputimages=0
    if (l_outputimages != "" && l_outputimages != " ") {
        if (substr(l_outputimages,1,1) == "@")
            scanfile=substr(l_outputimages,2,strlen(l_outputimages))
        else {
            if (stridx("*",l_outputimages) > 0) {
                files(l_outputimages,sort-) | match(".hhd",stop+,print-,
                    metach-, > tmpfile)
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
                glogprint( l_logfile, "tview", "status", type="error",
                    errno=121, str="Maximum number of output images "//\
                    str(maximages)//" exceeded.",verbose+)
                status=1
                goto exit
            }
            out[noutputimages]=l_filename
            if (imaccess(out[noutputimages])) {
                glogprint( l_logfile, "tview", "status", type="error",
                    errno=102, str="Output image "//l_filename//" already \
                    exists.",verbose+)
                status=1
                goto exit
            }
        }
        if (noutputimages != nimages) {
            glogprint( l_logfile, "tview", "status", type="error", errno=121,
                str="Different number of input ("//str(nimages)//") and \
                output ("//str(noutputimages)//" image names have been \
                specified.",verbose+)
            status=1
            goto exit
        }

        scanfile=""
        delete(tmpfile,verify-, >& "dev$null")
    } else {
        if (l_prefix == "" || l_prefix == " ")
            l_prefix="v"
        for (i=1; i <= nimages; i+=1) {
            out[i]=l_prefix//out[i]
            if (imaccess(out[i])) {
                glogprint( l_logfile, "tview", "status", type="error",
                    errno=102, str="Output image "//out[i]//" already exists.",
                    verbose+)
                status=1
                goto exit
            }
        }
    }

    #if a label is requested, set the colour value
    ncol=202
    if (l_label) {
        if (l_colour == "black") ncol=202
        else if (l_colour == "white") ncol=203
        else if (l_colour == "red") ncol=204
        else if (l_colour == "green") ncol=205
        else if (l_colour == "blue") ncol=206
        else if (l_colour == "yellow") ncol=207
        else if (l_colour == "cyan") ncol=208
        else if (l_colour == "magenta") ncol=209
        else if (l_colour == "orange") ncol=212
        else ncol=202
    }

    l=1
    while (l <= nimages) {
        changefile=0
        for (i=1; i <= maxbad; i+=1) {
            badsetnumber[i]=0
        }

        fparse(in[l])
        in1=fparse.directory//fparse.root//".fits"
        fparse(out[l])
        out1=fparse.directory//fparse.root//".fits"

        copy(in1,out1,verbose-)
        glogprint( l_logfile, "tview", "task", type="string",
            str="Copying image "//in[l]//" to "//out[l],verbose=l_verbose)
        l_image=out[l]

        #check if image exists
        if (!imaccess(l_image)) {
            glogprint( l_logfile, "tview", "status", type="error", errno=101,
                str="Copied image "//l_image//" is not found",verbose+)
            goto nextimage
        }

        glogprint( l_logfile, "tview", "status", type="string",
            str="Checking image "//l_image,verbose=l_verbose)
        glogprint( l_logfile, "tview", "status", type="fork", fork="forward",
            child="tcheckstructure", verbose=l_verbose )

        tcheckstructure(l_image,logfile=l_logfile,verbose=l_verbose)

        glogprint( l_logfile, "tview", "status", type="fork", fork="backward",
            child="tcheckstructure", verbose=l_verbose )

        j=tcheckstructure.status
        modeflag=tcheckstructure.modeflag

        if (j != 0) {
            glogprint( l_logfile, "tview", "status", type="error", errno=123,
                str="Image "//l_image//" does not have the expected \
                structure.", verbose+)
            status=1
            goto nextimage
        }

        phu=l_image//"[0]"

        if (modeflag == 1 || modeflag == 3) {
            # CHOP-NOD or NOD modes
            imgets(phu,"NNODS",>& "dev$null") ; nnods=int(imgets.value)
            imgets(phu,"NNODSETS",>& "dev$null") ; nnodsets=int(imgets.value)
            imgets(phu,"SAVESETS",>& "dev$null") ; nsavesets=int(imgets.value)
            imgets(phu,"NEXTEND",>& "dev$null") ; nextns=int(imgets.value)
            if (nextns == 0)
                nextns=nnods*nnodsets
        } else {
            # CHOP or STARE mode
            imgets(phu,"SAVESETS",>& "dev$null") ; nsavesets=int(imgets.value)
            nextns=1
        }

        #create temporary images for each type (time-consuming step)
        glogprint( l_logfile, "tview", "visual", type="visual",
            vistype="empty", verbose=l_verbose )
        glogprint( l_logfile, "tview", "status", type="string",
            str="Extracting frames...", verbose=l_verbose)
        glogprint( l_logfile, "tview", "visual", type="visual",
            vistype="empty", verbose=l_verbose )
        na=0

        for (i=1;i<=nextns;i+=1) {
            #
            # Read in any existing "bad set" keywords
            #
            imgets(l_image//"["//str(i)//"]","NBADSET",>& "dev$null")
            nbadsets=int(imgets.value)
            badoff=0
            if (nbadsets < 0 || nbadsets > maxbad || nbadsets > nsavesets) {
                glogprint( l_logfile, "tview", "status", type="warning",
                    errno=0, str="The number of bad sets was read in as "//\
                    str(nbadsets)//" for extension "//str(i),verbose=l_verbose)
                glogprint( l_logfile, "tview", "status", type="warning",
                    errno=0, str="The marking of bad frames is disabled for \
                    this frame.", verbose=l_verbose)
                badoff=1
            }
            nbadsetsorig=nbadsets
            if (nbadsets > 0) {
                for (k=1; k <= nbadsets; k+=1) {
                    if (k < 10)
                        imgets(l_image//"["//str(i)//"]","BADSET0"//str(k))
                    else
                        imgets(l_image//"["//str(i)//"]","BADSET"//str(k))
                    na=int(imgets.value)
                    if (na < 1 || na > nsavesets) {
                        glogprint(l_logfile, "tview", "status", type="warning",
                            errno=0, str="The number of bad sets was read in \
                            as "//str(nbadsets)//" for extension "//str(i),
                            verbose=l_verbose)
                        glogprint(l_logfile, "tview", "status", type="warning",
                            errno=0, str="The marking of bad frames is \
                            disabled for this frame.",verbose=l_verbose)
                        badoff=1
                    }
                    badsetnumber[k]=na
                }
            }

            na=(i+1)/2
            imgets(l_image//"["//str(i)//"]","nod")
            for (k=1;k<=nsavesets;k+=1) {
                # Create tmp FITS file names used within this loop
                tref1=mktemp("tmpref1")
                tsig1=mktemp("tmpsig1")

                badoff=0
                badflag=0
                badtoggle=0
                if (nbadsetsorig > 0) {
                    for (j=1; j <= nbadsets; j+=1) {
                        if (badsetnumber[j] == k) {
                            j=nbadsets+10
                            badoff=1
                            badflag=1
                        }
                    }
                }
                #
                # NOTE: the following assumes that the extenstions are for 
                # nod positions A then B then A then B and so on....it will 
                # fail if this is not the case.  The first extenstion also 
                # must be for nod position A for it to work.
                #
                if (modeflag == 1 || modeflag == 3) {
                    if (imgets.value == "A") {
                        beamval=1
                        imcopy(l_image//"["//str(i)//"][*,*,1,"//str(k)//"]",
                            tsig1, >& "dev$null")
                        imcopy(l_image//"["//str(i)//"][*,*,2,"//str(k)//"]",
                            tref1, >& "dev$null")
                    } else if (imgets.value == "B") {
                        beamval= -1
                        imcopy(l_image//"["//str(i)//"][*,*,2,"//str(k)//"]",
                            tsig1, >& "dev$null")
                        imcopy(l_image//"["//str(i)//"][*,*,1,"//str(k)//"]",
                            tref1, >& "dev$null")
                    } else {
                        glogprint( l_logfile, "tview", "status", type="error",
                            errno=132, str="Unrecognized nod: "//imgets.value,
                            verbose=l_verbose )
                        goto clean
                    }
                    tdif1=mktemp("tmpdif1")
                    imarith(tsig1,"-",tref1,tdif1, >& "dev$null")
                } else {
                    if (modeflag == 4) {
                        tdif1=mktemp("tmpdif1")
                        imcopy(l_image//"["//str(i)//"][*,*,1,"//str(k)//"]",
                            tdif1, >& "dev$null")
                    }
                    if (modeflag == 2) {
                        imcopy(l_image//"["//str(i)//"][*,*,1,"//str(k)//"]",
                            tsig1, >& "dev$null")
                        imcopy(l_image//"["//str(i)//"][*,*,2,"//str(k)//"]",
                            tref1, >& "dev$null")
                        tdif1=mktemp("tmpdif1")
                        imarith(tsig1,"-",tref1,tdif1, >& "dev$null")
                    }
                }
                if (i == 1 && k == 1) {
                    if (l_type == "dif" || modeflag == 4)
                        display(tdif1,1,erase+,fill=l_dispfill,zscale=l_zscale,
                            zrange=l_zrange,z1=l_z1,z2=l_z2,select_frame=yes,
                            >& "dev$null")
                    if (l_type == "sig" && modeflag != 4)
                        display(tsig1,1,erase+,fill=l_dispfill,zscale=l_zscale,
                            zrange=l_zrange,z1=l_z1,z2=l_z2,select_frame=yes,
                            >& "dev$null")
                    if (l_type == "ref" && modeflag != 4)
                        display(tref1,1,erase+,fill=l_dispfill,zscale=l_zscale,
                            zrange=l_zrange,z1=l_z1,z2=l_z2,select_frame=yes,
                            >& "dev$null")
                    if (l_manual)
                        printf("\n")
                    printf("Available key commands:\n  \
                        (h) help\n  (b) mark as a bad frame\n  \
                        (u) unmark as a bad frame\n  \
                        (i) run imexamine on this image\n  \
                        (n) go to next extension\n  \
                        (q) stop interactive mode\n  \
                        (s) get image statistics\n  \
                        (x) exit immediately\nKey commands can be entered \
                        in upper or lower case.\nOther keystrokes mean go \
                        to the next frame.\n\n")
                } else {
                    if (l_type == "dif" || modeflag == 4)
                        display(tdif1,1,erase-,fill=l_dispfill,zscale=l_zscale,
                            zrange=l_zrange,z1=l_z1,z2=l_z2,select_frame=yes,
                            >& "dev$null")
                    if (l_type == "sig" && modeflag != 4)
                        display(tsig1,1,erase-,fill=l_dispfill,zscale=l_zscale,
                            zrange=l_zrange,z1=l_z1,z2=l_z2,select_frame=yes,
                            >& "dev$null")
                    if (l_type == "ref" && modeflag != 4)
                        display(tref1,1,erase-,fill=l_dispfill,zscale=l_zscale,
                            zrange=l_zrange,z1=l_z1,z2=l_z2,select_frame=yes,
                            >& "dev$null")
                }
                if (l_label) {
                    xpos=10
                    ypos=10
                    tmpfile=mktemp("tmplabel")
                    if (badflag == 0) {
                        if (modeflag == 1 || modeflag == 3)
                            print(str(xpos)//" "//str(ypos)//" 'Nod "//\
                                str(na)//" pos. "//imgets.value//" chop "//\
                                str(k)//"'", >> tmpfile)
                        else
                            print(str(xpos)//" "//str(ypos)//" 'Chop "//\
                                str(k)//"'", >> tmpfile)
                    }
                    if (badflag == 1) {
                        if (modeflag == 1 || modeflag == 3)
                            print(str(xpos)//" "//str(ypos)//" 'Nod "//\
                                str(na)//" pos. "//imgets.value//" chop "//\
                                str(k)//" BAD FRAME'", >> tmpfile)
                        else
                            print(str(xpos)//" "//str(ypos)//" 'Chop "//\
                                str(k)//" BAD FRAME'", >> tmpfile)
                    }
                    tvmark(1,coords=tmpfile,label=yes,interactive=no,
                        txsize=2.0,color=ncol)
                    delete(tmpfile,verify-)
                }
                #
                # Cursor input "b" or "B" for bad frames, etc
                #
                if (l_delay > 0.)
                    sleep(l_delay)
                npars=0
                if (l_exam) {
                    if (l_verbose) printf("     TVIEW Entering imexam.\n")
                    imexamine()
                    if (l_verbose) {
                        printf("     TVIEW Exiting imexam.\n")
                        printf("     TVIEW Ready for image screening.\n")
                    }
                }
                if (l_manual) {
                    while (npars == 0) {
                        npars=1
                        if (fscan(imcur,xpos,ypos,wcs,cursinput) != EOF) {
                            if (cursinput == "q" || cursinput == "Q") {
                                glogprint( l_logfile, "tview", "status",
                                    type="string", str="Exiting manual mode.",
                                    verbose=l_verbose)
                                l_manual=no
                                npars=2
                            }
                            if (cursinput == "x" || cursinput == "X") {
                                glogprint( l_logfile, "tview", "status",
                                    type="string", str="Exiting the screening \
                                    loop.",verbose=l_verbose)
                                if (nbadsets != nbadsetsorig && nbadsets > 0) {
                                    gemhedit (l_image//"["//str(i)//"]",
                                        "NBADSET", nbadsets, "", delete-)
                                    for (j=1; j <= nbadsets; j+=1) {
                                        if (j < 10)
                                            gemhedit (l_image//"["//str(i)//\
                                                "]", "BADSET0"//str(j),
                                                badsetnumber[j], "", delete-)
                                        else
                                            gemhedit (l_image//"["//str(i)//\
                                                "]", "BADSET"//str(j),
                                                badsetnumber[j], "", delete-)
                                    }
                                    if (modeflag == 1 || modeflag == 3) {
                                        gemhedit (l_image//"["//str(i+beamval)\
                                            //"]", "NBADSET", nbadsets, "", 
                                            delete-)
                                        for (j=1; j <= nbadsets; j+=1) {
                                            if (j > 9)
                                                gemhedit (l_image//"["\
                                                    //str(i+beamval)//"]",
                                                    "BADSET"//str(j),
                                                    badsetnumber[j], "", 
                                                    delete-)
                                            else
                                                gemhedit (l_image//"["\
                                                    //str(i+beamval)//"]",
                                                    "BADSET0"//str(j),
                                                    badsetnumber[j], "",
                                                    delete-)
                                        }
                                    }
                                    changefile=changefile+1
                                }
                                goto checkchanges
                            } #end if cursinput "x"
                            #
                            # For the "h" command, loop to get another input 
                            # value.  Any value aside from the defined ones 
                            # causes the next image to be displayed.
                            #
                            if (cursinput == "h" || cursinput == "H") {
                                printf("\n")
                                printf("Available key commands:\n  \
                                    (h) help\n  \
                                    (b) mark as a bad frame\n  \
                                    (u) unmark as a bad frame\n  \
                                    (i) run imexamine on this image\n  \
                                    (n) go to next extension\n  \
                                    (q) stop interactive mode\n  \
                                    (s) get image statistics\n  \
                                    (x) exit immediately\nKey commands can be \
                                    entered in upper or lower case.\n\
                                    Other keystrokes mean go to the next \
                                    frame.\n")
                                printf("\n")
                                npars=0
                            }

                            if (cursinput == "N" || cursinput == "n")
                                npars=3

                            if (cursinput == "s" || cursinput == "S") {
                                imstat(tdif1)
                                npars=0
                            }

                            if (cursinput == "b" || cursinput == "B") {
                                if (badoff == 0) {
                                    glogprint( l_logfile, "tview", "science",
                                        type="string", str="Frame "//str(k)//\
                                        " of nod "//str(i)//" marked as bad.",
                                        verbose=l_verbose)
                                    nbadsets=nbadsets+1
                                    badsetnumber[nbadsets]=k
                                    badoff=2
                                    badtoggle=1
                                }
                                npars=0
                            }

                            if (cursinput == "u" || cursinput == "U") {
                                if (badtoggle == 1) {
                                    glogprint( l_logfile, "tview", "science",
                                        type="string", str="Frame "//str(k)//\
                                        " of nod "//str(i)//" unmarked as \
                                        bad.", verbose=l_verbose)
                                    badsetnumber[nbadsets]=0
                                    nbadsets=nbadsets-1
                                    if (badoff == 2) badoff=0
                                    badtoggle=0
                                }
                                npars=0
                            }
                            if (cursinput == "i" || cursinput == "I") {
                                if (l_verbose)
                                    printf("     TVIEW Entering imexam.\n")
                                imexamine()
                                if (l_verbose) {
                                    printf("     TVIEW Exiting imexam.\n")
                                    printf("     TVIEW Ready for image \
                                        screening.\n")
                                }
                                npars=0
                            }
                            if (npars > 0) {
                                if (npars == 1)
                                    glogprint( l_logfile, "tview",
                                        "engineering", type="string",
                                        str="Going to next frame",
                                        verbose=l_verbose)
                                if (npars == 3) {
                                    glogprint( l_logfile, "tview",
                                        "engineering", type="string",
                                        str="Going to next extension",
                                        verbose=l_verbose)
                                    imdelete(tdif1,verify-, >& "dev$null")
                                    imdelete(tsig1,verify-, >& "dev$null")
                                    imdelete(tref1,verify-, >& "dev$null")
                                    goto nextexten  #getting out of the loop
                                } 
                            }
                        } #end if fscan imcur
                    } #end while (npars==0)
                } #end if (l_manual)

nextframe:
                imdelete(tdif1,verify-, >& "dev$null")
                imdelete(tsig1,verify-, >& "dev$null")
                imdelete(tref1,verify-, >& "dev$null")
            } #end for-loop over savesets

nextexten:
            npars=0
            if (nbadsets != nbadsetsorig && nbadsets > 0) {
                gemhedit (l_image//"["//str(i)//"]", "NBADSET", nbadsets, "",
                    delete-)
                for (j=1; j <= nbadsets; j+=1) {
                    if (j < 10)
                        gemhedit (l_image//"["//str(i)//"]", "BADSET0"//str(j),
                            badsetnumber[j], "", delete-)
                    else
                        gemhedit (l_image//"["//str(i)//"]", "BADSET"//str(j),
                            badsetnumber[j], "", delete-)
                }
                if (modeflag == 1 || modeflag == 3) {
                    gemhedit (l_image//"["//str(i+beamval)//"]", "NBADSET",
                        nbadsets, "", delete-)
                    for (j=1; j <= nbadsets; j+=1) {
                        if (j > 9)
                            gemhedit (l_image//"["//str(i+beamval)//"]",
                                "BADSET"//str(j), badsetnumber[j], "", delete-)
                        else
                            gemhedit (l_image//"["//str(i+beamval)//"]",
                                "BADSET0"//str(j), badsetnumber[j], "",
                                delete-)
                    }
                }
                changefile=changefile+1
            }
            # KL - below, why is the value specified if the keyword is to be 
            # deleted?  (just asking in case deleting is not the desired 
            # behaviour :)
            if (nbadsets < nbadsetsorig) {
                if (nbadsets == 0)
                    gemhedit (l_image//"["//str(i)//"]", "NBADSET", nbadsets,
                        "", delete=yes)
                for (j=nbadsets+1; j < nbadsetsorig; j+=1) {
                    if (j > 9)
                        gemhedit (l_image//"["//str(i)//"]", "BADSET"//str(j),
                            badsetnumber[j], "", delete=yes)
                    else
                        gemhedit (l_image//"["//str(i)//"]", "BADSET0"//str(j),
                            badsetnumber[j], "", delete=yes)
                }
                if (modeflag == 1 || modeflag == 3) {
                    if (nbadsets == 0)
                        gemhedit (l_image//"["//str(i+beamval)//"]", "NBADSET",
                            nbadsets, "", delete=yes)
                    for (j=nbadsets+1; j < nbadsetsorig; j+=1) {
                        if (j > 9)
                            gemhedit (l_image//"["//str(i+beamval)//"]",
                                "BADSET"//str(j), badsetnumber[j], "",
                                delete=yes)
                        else
                            gemhedit (l_image//"["//str(i+beamval)//"]",
                                "BADSET0"//str(j), badsetnumber[j], "",
                                delete=yes)
                    }
                }
            }
        } # end for-loop over extensions

checkchanges:
        if (changefile > 0) {
            glogprint( l_logfile, "tview", "task", type="string",
                str="Changes were made to file "//out[l],verbose=l_verbose)
            # Time stamp the primary header
            #
            gemdate ()
            gemhedit (out[l]//"[0]", "TVIEW", gemdate.outdate,
                "UT Time stamp for TVIEW")
            gemhedit (out[l]//"[0]", "GEM-TLM", gemdate.outdate,
                "UT Last modification with GEMINI IRAF")
        } else {
            if (l_delete) {
                imdelete(out[l],verify-)
                glogprint( l_logfile, "tview", "task", type="string",
                    str="No changes were made to file "//out[l]//"; deleting \
                    copy.",verbose=l_verbose)
            } else
                glogprint( l_logfile, "tview", "task", type="string",
                    str="No changes were made to file "//out[l]//".",
                    verbose=l_verbose)
        }

clean:
        if (imaccess(tref1)) imdelete(tref1,verify-, >& "dev$null")
        if (imaccess(tsig1)) imdelete(tsig1,verify-, >& "dev$null")
        if (imaccess(tdif1)) imdelete(tdif1,verify-, >& "dev$null")

nextimage:
        l=l+1
    }

exit:
    if (status == 0)
        glogclose( l_logfile, "tview", fl_success+, verbose=l_verbose )
    else
        glogclose( l_logfile, "tview", fl_success-, verbose=l_verbose )

exitnow:
    ;

end

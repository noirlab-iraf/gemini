# Copyright(c) 2003-2011 Association of Universities for Research in Astronomy, Inc.

procedure tcheckstructure (image)

#
# This is a task that tests MEF file integrity for 
# TReCS files.  It checks the number of extensions and 
# the size of the extensions against the header keywords.
#
# It should work for the output files from "tprepare" as well 
# as the raw T-ReCS files.
#
# Version:   Oct 5,  2003  KV routine changed to not write to the terminal
#                            unless an error oocurs
#            Sept 3, 2003  KV original routine created

char    image       {prompt="Image to check"}
char    logfile     {"",prompt="Logfile name"}
bool    verbose     {no,prompt="Verbose?"}
int     modeflag    {0,prompt="Exit observation mode flag"}
int     status      {0,prompt="Exit error status: (0=good, >0=bad)"}

begin

    char    l_image, l_logfile, extname, tmpstring, paramstr
    struct  l_struct
    bool    l_verbose
    int     nnods, nnodsets, nextends, savesets
    int     iaxis, naxis[4]
    int     i, j, jpre

    l_image=image
    l_logfile=logfile
    l_verbose=verbose

    # Initialize
    status = 0

    # Create the list of parameter/value pairs.  One pair per line.
    # All lines combined into one string.  Line delimiter is '\n'.
    paramstr =  "image          = "//image.p_value//"\n"
    paramstr += "modeflag       = "//modeflag.p_value//"\n"
    paramstr += "logfile        = "//logfile.p_value//"\n"
    paramstr += "verbose        = "//verbose.p_value

    # Assign a logfile name, if not given.  Open logfile and start log.
    # Write parameter/value pairs ("paramstr") to log.
    gloginit (l_logfile, "tcheckstructure", "midir", paramstr, fl_append+,
        verbose=l_verbose)
    if (gloginit.status != 0) {
        status = gloginit.status
        goto exitnow
    }
    l_logfile = gloginit.logfile.p_value

    #check if image exists
    if (!imaccess(l_image)) {
        glogprint (l_logfile, "tcheckstructure", "status", type="error",
            errno=101, str="Image "//l_image//" is not found.",verbose+)
        status = 1
        goto exit
    }


    # Now, read the main header for relevant values

    modeflag = 0
    extname = l_image//"[0]"
    if (!imaccess(extname)) {
        glogprint (l_logfile, "tcheckstructure", "status", type="error",
            errno=101, str="Primary header unit for file "//l_image//\
            " is not found",verbose+)
        status = 1
        goto exit
    }

    imgets(extname,"INSTRUME",>& "dev$null")
    if ((imgets.value != "TReCS") && (imgets.value != "CanariCam")){
        glogprint (l_logfile, "tcheckstructure", "status", type="error",
            errno=121, str="Image "//l_image//" has 'INSTRUME' keyword \
            set to "//imgets.value//" rather than 'TReCS' or 'CanariCam.",
            verbose+)
        status = 1
        goto exit
    }

    imgets(extname,"MISTACK",>& "dev$null")
    if (imgets.value != "0") {
        glogprint (l_logfile, "tcheckstructure", "status", type="error",
            errno=121, str="This is a file from MISTACK.",verbose+)
        status = 1
        goto exit
    }
    imgets(extname,"MIREGISTER",>& "dev$null")
    if (imgets.value != "0") {
        glogprint (l_logfile, "tcheckstructure", "status", type="error",
            errno=121, str="This is a file from MIREGISTER.",verbose+)
        status = 1
        goto exit
    }

    imgets(extname,"OBSMODE",>& "dev$null")
    if (imgets.value == "chop-nod")
        modeflag=1
    else if (imgets.value == "chop")
        modeflag=2
    else if (imgets.value == "nod")
        modeflag=3
    else if (imgets.value == "stare")
        modeflag=4
    else {
        glogprint (l_logfile, "tcheckstructure", "status", type="error",
            errno=132, str="Image "//l_image//" has unrecognised OBSMODE \
            keyword.",verbose+)
        status = 1
        goto exit
    }

    imgets(extname,"TPREPARE",>& "dev$null")
    if (imgets.value == "0")
        jpre=1
    else 
        jpre=2

    imgets(extname,"NNODS",>& "dev$null")
    nnods = int(imgets.value)
    imgets(extname,"NNODSETS",>& "dev$null")
    nnodsets = int(imgets.value)
    if (jpre == 1)
        nextends=nnodsets*nnods
    else
        nextends=1
    imgets(extname,"SAVESETS",>& "dev$null")
    savesets = int(imgets.value)
    if (nnods < 1 || nnods > 2 || nnodsets < 1 || savesets < 1) {
        tmpstring = "Primary header error: NNODS = "//str(nnods)//\
            " NNODSETS = "//str(nnodsets)//" SAVESETS = "//str(savesets)
        glogprint (l_logfile, "tcheckstructure", "status", type="error",
            errno=132, str=tmpstring, verbose+)
        status = 1
        goto exit
    }

    # Check the extensions and the dimension values

    j= -10
    for (i=1; i < 10000 ; i+=1) {
        extname = l_image//"["//i//"]"
        if (!imaccess(extname)) {
            j=i-1
            i=10001
        }
    }

    if (j != nextends) {
        tmpstring = "The file has "//str(j)//" extensions.\nThe expected \
            number of extensions is "//str(nextends)//"."
        glogprint (l_logfile, "tcheckstructure", "status", type="error",
            errno=123, str=tmpstring, verbose+)
        status = 1
        goto exit
    }

    for (i=1; i <= nextends; i+=1) {
        extname = l_image//"["//i//"]"
        if (imaccess(extname)) {
            imgets(extname,"i_naxis")
            iaxis = int(imgets.value)
            if (iaxis < 2 || iaxis > 4) {
                tmpstring = "Image extension "//str(i)//" has NAXIS = "//\
                    str(iaxis)//".\nThe expected number is 2, 3 or 4"
                glogprint (l_logfile, "tcheckstructure", "status", \
                    type="error", errno=123, str=tmpstring, verbose+)
                status = status + 1
            } else {
                naxis[1] = 0
                naxis[2] = 0
                naxis[3] = 0
                naxis[4] = 0

                # add 10 to the modeflag if the file looks like an output 
                # file from "tprepare".

                if (iaxis == 2 || iaxis == 3 && modeflag < 10) 
                    modeflag = modeflag+10

                for (j=1; j <= iaxis; j+=1) {
                    tmpstring = "naxis"//str(j)
                    imgets(extname,tmpstring)
                    naxis[j] = int(imgets.value)
                }

                # The following assumes that NOD mode observations have 
                # structure [320,240,1,nsavesets]: no chopping so only one 
                # chop position gives third dimension 1.
                #
                # This is to be verified.

                if (naxis[1] != 320 || naxis[2] != 240) {
                    if (iaxis == 2)
                        tmpstring = "Image extension "//str(iaxis)//" has bad \
                            image dimensions: "//str(naxis[1])//" "//\
                            str(naxis[2])//".\nThe expected numbers are \
                            [320,240]"
                    else if (iaxis == 3)
                        tmpstring = "Image extension "//str(iaxis)//" has bad \
                            image dimensions: "//str(naxis[1])//" "//\
                            str(naxis[2])//" "//str(naxis[3])//".\nThe \
                            expected numbers are [320,240,3]"
                    else if (iaxis == 4) {
                        if (modeflag == 1 || modeflag == 2)
                            tmpstring = "Image extension "//str(iaxis)//" has \
                                bad image dimensions: "//str(naxis[1])//" "//\
                                str(naxis[2])//" "//str(naxis[3])//" "//\
                                str(naxis[4])//".\nThe expected numbers are \
                                [320,240,2,"//str(savesets)//"]"
                        else if (modeflag == 3 || modeflag == 4)
                            tmpstring = "Image extension "//str(iaxis)//" has \
                                bad image dimensions: "//str(naxis[1])//" "//\
                                str(naxis[2])//" "//str(naxis[3])//" "//\
                                str(naxis[4])//".\nThe expected numbers are \
                                [320,240,1,"//str(savesets)//"]"
                    }
                    glogprint (l_logfile, "tcheckstructure", "status", 
                        type="error", errno=123, str=tmpstring, verbose+)
                    status = status + 1
                }
                if (iaxis == 2 && modeflag != 4) {
                    tmpstring = "Image extension "//str(iaxis)//" has bad \
                        image dimensions: "//str(naxis[1])//" "//\
                        str(naxis[2])//".\nThe OBSMODE is not 'stare'."
                    glogprint (l_logfile, "tcheckstructure", "status",
                        type="error", errno=123, str=tmpstring, verbose+)
                    status = status + 1
                }
                if (iaxis == 3 && naxis[3] != 3) {
                    tmpstring = "Image extension "//str(iaxis)//" has bad \
                        image dimensions: "//str(naxis[1])//" "//\
                        str(naxis[2])//" "//str(naxis[3])//".\nThe expected \
                        numbers are [320,240,3]"
                    glogprint (l_logfile, "tcheckstructure", "status",
                        type="error", errno=123, str=tmpstring, verbose+)
                    status = status + 1
                }
                if (iaxis == 4) {
                    if (modeflag == 2 || modeflag == 1) {
                        if (naxis[3] != 2 || naxis[4] != savesets) {
                            tmpstring = "Image extension "//str(iaxis)//" has \
                                bad image dimensions: "//str(naxis[1])//" "//\
                                str(naxis[2])//" "//str(naxis[3])//" "//\
                                str(naxis[4])//".\nThe expected numbers are \
                                [320,240,2,"//savesets//"]"
                            glogprint (l_logfile, "tcheckstructure", "status",
                                type="error", errno=123, str=tmpstring,
                                verbose+)
                            status = status + 1
                        }
                    }
                    if (modeflag == 3 || modeflag == 4) {
                        if (naxis[3] != 1 || naxis[4] != savesets) {
                            tmpstring = "Image extension "//str(iaxis)//" has \
                                bad image dimensions: "//str(naxis[1])//" "//\
                                str(naxis[2])//" "//str(naxis[3])//" "//\
                                str(naxis[4])//".\nThe expected numbers are \
                                [320,240,1,"//savesets//"]"
                            glogprint (l_logfile, "tcheckstructure", "status",
                                type="error", errno=123, str=tmpstring,
                                verbose+)
                            status = status + 1
                        }
                    }
                }
            }
        } else {
            tmpstring = "Image extension "//str(i)//" was not opened \
                successfully.\nThe expected number of extensions is "//\
                str(nextends)//"."
            glogprint (l_logfile, "tcheckstruture", "status", type="error",
                errno=101, str=tmpstring, verbose+)
            status = 1
            goto exit
        }
    }


exit:
    if (status == 0) {
        glogprint (l_logfile, "tcheckstructure", "status", type="string",
            str="File "//l_image//" appears to be correct.",verbose=l_verbose)
        glogclose (l_logfile, "tcheckstructure", fl_success+, \
            verbose=l_verbose)
    } else {
        glogprint (l_logfile, "tcheckstructure", "status", type="string",
            str="File "//l_image//" has some structural errors.",
            verbose=l_verbose)
        glogclose (l_logfile, "tcheckstructure", fl_success-, \
            verbose=l_verbose)
    }

exitnow:
    ;

end

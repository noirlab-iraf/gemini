# Copyright(c) 2013-2015 Association of Universities for Research in Astronomy, Inc.

procedure gemsecchk (inimage, section)

# Given an image (single extenion) and a section the section is checked and
# return in a way IRAF can understand.
#
# Allowed section inputs are:
#
#     [*]
#     [*,*]
#     [N:M,R:S]
#     [N%]
#     [N%,M%]
#     [N%:M%,R%]
#     [N%:M%,R%:S%]
#     [-N,-M]
#     [-N:-M,-R]
#     [-N:-M,-R:-S]
#     [N%,-M]
#     [N%:M%,-R]
#     [N%:M%,-R:-S]
#     [-N,M%]
#     [-N:-M,R%]
#     [-N:-M,R%:S%]
#     Then you can replace any ?% or -? with ? to represent a specific pixel

char    inimage     {prompt="Input image NOT including section"}
char    section     {prompt="Input Section"}
char    out_imgsect {"",prompt="Output image name with section"}
char    logfile     {"",prompt="Logfile"}
bool    verbose     {yes,prompt="Verbose?"}
int     status      {0,prompt="Exit status (0=good)"}

####

begin

    # Local variables
    char l_inimage, l_section, l_out_imgsect, l_imgsect, l_logfile, tstring
    char imgsec, x_sec, y_sec, x1info, x2info, y1info, y2info, Xinfo, Yinfo

    int  secx1, secx2, secy1, secy2, boarder, ftest
    int  commapos1, commapos2, Xcolonpos, Ycolonpos
    int  l_status, ndim, naxis1, naxis2

    bool sxcheck, sycheck, Xidentical, Yidentical, l_verbose

    # Read input parameters
    l_inimage = inimage
    l_section = section
    l_out_imgsect = ""
    l_logfile = logfile
    l_verbose = verbose
    l_status = 1

    # Set default values
    sxcheck = yes
    sycheck = yes
    Xidentical = no
    Yidentical = no

    secx1 = 0
    secx2 = 0
    secy1 = 0
    secy2 = 0

    cache ("keypar", "gemdate")

    # Test the log file
    if ((l_logfile == "") || (stridx(" ", l_logfile) > 0)) {
       l_logfile = gemtools.logfile
       if ((l_logfile == "") || (stridx (" ", l_logfile) > 0)) {
          l_logfile = "gemtools.log"
          printlog ("WARNING - GEMSECCHK: Both gastat.logfile and \
              gemtools.logfile fields are empty", l_logfile, l_verbose)
          printlog ("                     Using default file gemtools.log", \
              l_logfile, l_verbose)
       }
    }

    # Start logging
    printlog ("", l_logfile, l_verbose)
    printlog ("------------------------------------------------------------", \
        l_logfile, l_verbose)
    gemdate (zone="local")
    printlog ("GEMSECCHK -- Started: "//gemdate.outdate//"\n", \
        l_logfile, l_verbose)

    printlog ("GEMSECCHK: Input parameters...", l_logfile, l_verbose)
    printlog ("    inimage     = "//l_inimage, l_logfile, l_verbose)
    printlog ("    section     = "//l_section, l_logfile, l_verbose)
    printlog ("", l_logfile, l_verbose)

    ####
    # Start input checks
    if (!imaccess(l_inimage)) {
        printlog ("ERROR GEMSECCHK: Cannot access "//l_inimage, \
            l_logfile, verbose+)
        goto crash
    }

    fparse (l_inimage)

    if (fparse.section != "") {

        # More than one section supplied
        printlog ("ERROR - GEMSECCHK: section is "//l_section//" and "//\
            "a section was supplied as part of file name."//\
            "\n                   Please only supply section in the section \
            parameter", l_logfile, verbose+)
        goto crash
    }

    if (l_section == "") {
        imgsec = "*"

    } else {

        # Strip the "[" and "]"
        imgsec = l_section

        if (substr(imgsec,1,1) == "[") {
            imgsec = substr(imgsec,2,strlen(imgsec))
        }
        if (substr(imgsec,strlen(imgsec),strlen(imgsec)) == "]") {
            imgsec = substr(imgsec,1,strlen(imgsec)-1)
        }
    }

    # Determine the dimensionality of the image before assuming below that
    # the section should be 2-D:
    keypar (l_inimage, "i_naxis", silent+)
    if (!keypar.found) {
        # I don't think this ever happens but it's harmless (also below):
        printlog ("ERROR - GEMSECCHK: Keyword NAXIS not found in "//img,
            l_logfile, verbose+)
        goto crash
    }
    ndim = int(keypar.value)

    # Require size of image to help with imgsec parsing and setting
    keypar (l_inimage, "i_naxis1", silent+)
    if (!keypar.found) {
        printlog ("ERROR - GEMSECCHK: Keyword "//\
            "NAXIS1 not found in "//img,
            l_logfile, verbose+)
        goto crash
    }
    naxis1 = int(keypar.value)

    keypar (l_inimage, "i_naxis2", silent+)
    if (!keypar.found) {
        printlog ("ERROR - GEMSECCHK: Keyword "//\
            "NAXIS2 not found in "//img,
            l_logfile, verbose+)
        goto crash
    }
    naxis2 = int(keypar.value)

    if (ndim != 2 || naxis2 < 2) {
        printlog ("WARNING - GEMSECCHK: non-standard syntax is unsupported \
                  for images that aren't \n                     2-D; \
                  passing through the section unchecked & unmodified.", \
                  l_logfile, l_verbose)
        l_out_imgsect = l_section
        l_status = 0
        goto clean
    }

    if (sxcheck && (secx2 > naxis1)) {
        printlog ("ERROR - GEMSECCHK: X2 position of imgsec is \
            bigger than length of NAXIS1", l_logfile, verbose+)
        goto crash
    }

    # Parse the section out for checking later
    # Determine if there is a comma
    commapos1 = stridx (",",imgsec)
    commapos2 = strldx (",",imgsec)

    if (commapos1 != commapos2) {
        printlog ("ERROR - GEMSECCHK: Badly formatted image section \
            supplied"//imgsec, l_logfile, verbose+)
        goto crash
    }

    # Parse X and Y info
    if (commapos1 > 0) {
        # Both X and Y information supplied
        Xinfo = substr(imgsec,1,commapos1-1)
        Yinfo = substr(imgsec,commapos1+1,strlen(imgsec))

    } else {
        # Only one set of information supplied and to be applied to both X & Y
        Xinfo = substr(imgsec,1,strlen(imgsec))
        Yinfo = Xinfo

    } # End of parsing initial X and Y info

    # Parse out the x1info, x2info, y1info and y2info
    Xcolonpos = 0
    Ycolonpos = 0

    Xcolonpos = stridx (":",Xinfo)
    Ycolonpos = stridx (":",Yinfo)

    if (Xcolonpos > 1) {
        x1info = substr(Xinfo,1,Xcolonpos-1)
        x2info = substr(Xinfo,Xcolonpos+1,strlen(Xinfo))
    } else {
        x1info = Xinfo
        x2info = x1info
        Xidentical = yes
    }

    if (Ycolonpos > 1) {
        y1info = substr(Yinfo,1,Ycolonpos-1)
        y2info = substr(Yinfo,Ycolonpos+1,strlen(Yinfo))
    } else {
        y1info = Yinfo
        y2info = y1info
        Yidentical = yes
    }

    # Now parse and set if required x1,x2,y2,y2
    # If the 1 and 2 indicies are identical only do the 1 parsing and use that
    # to set second index.

    # Unfortunately there is a lot of code duplication here

    # secx1:
    if (x1info == "*") {
        sxcheck = no

        if (!Xidentical) {
            printlog ("ERROR - GEMSECCHK: Cannot supply *as part of a #:# \
                range", l_logfile, verbose+)
            goto crash
        }

    } else if (stridx("%",x1info) > 0) {
        #Percentages

        # Needs to be at the end
        if (stridx("%",x1info) != strlen(x1info)) {
            ftest = 0
        } else {
            # It's a percentage
            tstring = substr(x1info,1,stridx("%",x1info)-1)
            ftest = fscan (tstring,boarder)
        }

        if (ftest == 0) {
            printlog ("ERROR - GEMSECCHK: "//x1info//" is not a number", \
                l_logfile, verbose+)
            goto crash
        } else if (boarder < 0) {
            # Test for negative numbers
            printlog ("ERROR - GEMSECCHK: supplied X1 percentage is \
                negative", l_logfile, verbose+)
            goto crash
        }

        secx1 = int ((boarder/100.) * naxis2) + 1
        if (Xidentical) {
            secx2 = naxis2 - (secx1 - 1)
        }

    } else if (stridx("-",x1info) > 0) {
        # Negative number

        # Needs to be at the beginning
        if (stridx("-",x1info) != 1) {
            ftest = 0
        } else {
            # It's a negative number
            tstring  = substr(x1info,2,strlen(x1info))
            ftest = fscan (tstring,boarder)
        }

        if (ftest == 0) {
            printlog ("ERROR - GEMSECCHK: "//x1info//" is not a number", \
                l_logfile, verbose+)
            goto crash
        }
        # No need to check for "%"

        secx1 = int (boarder + 1)
        if (Xidentical) {
            secx2 = naxis2 - (secx1 - 1)
        }

    } else if (Xidentical) {
        printlog ("ERROR - GEMSECCHK: Only one number supplied for X range"//\
            "\n                   Cannot supply only one number unless that"//\
            "\n                   number if negative or a percentage", \
            l_logfile, verbose+)
        goto crash
    } else if (!Xidentical) {
        secx1 = int (x1info)
    }

    # x2:
    if (!Xidentical) {
        if (stridx("%",x2info) > 0) {
            #Percentages

            # Needs to be at the end
            if (stridx("%",x2info) != strlen(x2info)) {
                ftest = 0
            } else {
                # It's a percentage
                tstring = substr(x2info,1,stridx("%",x2info)-1)
                ftest = fscan (tstring,boarder)
            }

            if (ftest == 0) {
                printlog ("ERROR - GEMSECCHK: "//x2info//" is not a number", \
                    l_logfile, verbose+)
                goto crash
            } else if (boarder < 0) {
                # Test for negative numbers
                printlog ("ERROR - GEMSECCHK: supplied X2 percentage is \
                    negative", l_logfile, verbose+)
                goto crash
            }

            secx2 = naxis2 - int ((boarder/100.) * naxis2)

        } else if (stridx("-",x2info) > 0) {
            # Negative number

            # Needs to be at the beginning
            if (stridx("-",x2info) != 1) {
                ftest = 0
            } else {
                # It's a negative number
                tstring = substr(x2info,2,strlen(x2info))
                ftest = fscan (tstring,boarder)
            }

            if (ftest == 0) {
                printlog ("ERROR - GEMSECCHK: "//x2info//" is not a number", \
                    l_logfile, verbose+)
                goto crash
            }
            # No need to check for "%"

            secx2 = naxis2 - int (boarder)

        } else {
            secx2 = int (x2info)
        }
    }

    # y1:
    if (y1info == "*") {
        sycheck = no

        if (!Yidentical) {
            printlog ("ERROR - GEMSECCHK: Cannot supply *as part of a #:# \
                range", l_logfile, verbose+)
            goto crash
        }

    } else if (stridx("%",y1info) > 0) {
        #Percentages

        # Needs to be at the end
        if (stridx("%",y1info) != strlen(y1info)) {
            ftest = 0
        } else {
            # It's a percentage
            tstring = substr(y1info,1,stridx("%",y1info)-1)
            ftest = fscan (tstring,boarder)
        }

        if (ftest == 0) {
            printlog ("ERROR - GEMSECCHK: "//y1info//" is not a number", \
                l_logfile, verbose+)
            goto crash
        } else if (boarder < 0) {
            # Test for negative numbers
            printlog ("ERROR - GEMSECCHK: supplied Y1 percentage is \
                negative", l_logfile, verbose+)
            goto crash
        }

        secy1 = int ((boarder/100.) * naxis2) + 1
        if (Yidentical) {
            secy2 = naxis2 - (secy1 - 1)
        }

    } else if (stridx("-",y1info) > 0) {
        # Negative number

        # Needs to be at the beginning
        if (stridx("-",y1info) != 1) {
            ftest = 0
        } else {
            # It's a negative number
            tstring = substr(y1info,2,strlen(y1info))
            ftest = fscan (tstring,boarder)
        }

        if (ftest == 0) {
            printlog ("ERROR - GEMSECCHK: "//y1info//" is not a number", \
                l_logfile, verbose+)
            goto crash
        }
        # No need to check for "%"

        secy1 = int (boarder + 1)
        if (Yidentical) {
            secy2 = naxis2 - (secy1 - 1)
        }

    } else if (Yidentical) {
        printlog ("ERROR - GEMSECCHK: Only one number supplied for Y range"//\
            "\n                   Cannot supply only one number unless that"//\
            "\n                   number if negative or a percentage", \
            l_logfile, verbose+)
        goto crash
    } else {
        secy1 = int (y1info)
    }

    # y2:
    if (!Yidentical) {
        if (stridx("%",y2info) > 0) {
            #Percentages

            # Needs to be at the end
            if (stridx("%",y2info) != strlen(y2info)) {
                ftest = 0
            } else {
                # It's a percentage
                tstring = substr(y2info,1,stridx("%",y2info)-1)
                ftest = fscan (tstring,boarder)
            }

            if (ftest == 0) {
                printlog ("ERROR - GEMSECCHK: "//y2info//" is not a number", \
                    l_logfile, verbose+)
                goto crash
            } else if (boarder < 0) {
                # Test for negative numbers
                printlog ("ERROR - GEMSECCHK: supplied Y2 percentage is \
                    negative", l_logfile, verbose+)
                goto crash
            }

            secy2 = naxis2 - int ((boarder/100.) * naxis2)

        } else if (stridx("-",y2info) > 0) {
            # Negative number

            # Needs to be at the beginning
            if (stridx("-",y2info) != 1) {
                ftest = 0
            } else {
                # It's a negative number
                tstring = substr(y2info,2,strlen(y2info))
                ftest = fscan (tstring,boarder)
            }

            if (ftest == 0) {
                printlog ("ERROR - GEMSECCHK: "//y2info//" is not a number", \
                    l_logfile, verbose+)
                goto crash
            }
            # No need to test for "%"

            secy2 = naxis2 - int (boarder)

        } else {
            secy2 = int (y2info)
        }
    }

    # Check number for consistency
    if (sxcheck) {
        # Check secy1 and secy2 for consistency
        if (secx1 < 1  || secx1 >= secx2) {
            printlog ("ERROR - GEMSECCHK: Bad X range given in \
                inimage section", l_logfile, verbose+)
            goto crash
        } else if (secx2 > naxis1) {
            printlog ("ERROR - GEMSECCHK: X2 position of imgsec is \
                bigger than length of NAXIS2", l_logfile, verbose+)
            goto crash
        }
        x_sec = "["//secx1//":"//secx2
    } else {
        x_sec = "[*"
    }

    if (sycheck) {
        # Check secy1 and secy2 for consistency
        if (secy1 < 1  || secy1 >= secy2) {
            printlog ("ERROR - GEMSECCHK: Bad Y range given in \
                inimage section", l_logfile, verbose+)
            goto crash
        } else if (secy2 > naxis1) {
            printlog ("ERROR - GEMSECCHK: Y2 position of imgsec is \
                bigger than length of NAXIS2", l_logfile, verbose+)
            goto crash
        }

        y_sec = secy1//":"//secy2//"]"
    } else {
        y_sec = "*]"
    }

    l_out_imgsect = x_sec//","//y_sec

    # All OK
    printlog ("GEMSECCHK: Output image section -- "//l_out_imgsect, \
        l_logfile, l_verbose)
    l_status = 0
    goto clean

crash:

    l_out_imgsect = ""

clean:

    status = l_status
    out_imgsect = l_out_imgsect

    # Print finish time
    gemdate(zone="local")
    printlog ("\nGEMSECCHK -- Finished: "//gemdate.outdate, \
        l_logfile, l_verbose)

    if (status == 0) {
        printlog ("\nGEMSECCHK -- Exit status: GOOD", l_logfile, l_verbose)

    } else {
        printlog ("\nGEMSECCHK -- Exit status: ERROR", l_logfile, verbose+)
    }

    if (status != 0) {
        printlog ("          -- Please read the logfile \""//l_logfile//\
                "\" for more information.", l_logfile, l_verbose)
    }

    printlog ("------------------------------------------------------------", \
        l_logfile, l_verbose)

end

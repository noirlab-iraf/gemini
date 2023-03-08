# Copyright(c) 2012-2013 Association of Universities for Research in Astronomy, Inc.

procedure gadimschk (inimage)

# This check the dimensions and sections (based on requested keyword) and
# returns the section of chkimage that matches the input section of inimage. If
# section of inimage is not wholly contained within chkimage scripts exits with
# an error.

# If no chkimage is supplied, a image and a statsec can be checked on it's own

# Input must be a single extension

char    inimage     {prompt="Input GSAOI image"}
char    section     {"",prompt="Input image section"}
char    chkimage    {"",prompt="Input GSAOI image to check inimage against"}
char    key_check   {"CCDSEC",enum="CCDSEC|DETSEC",prompt="Keyword to check dimensions against"}
char    out_chkimage{"",prompt="chkimage name including corrected section"}
char    logfile     {"",prompt="Logfile"}
bool    verbose     {yes,prompt="Verbose?"}
int     status      {0,prompt="Exit status (0=good)"}


begin

    char l_inimage, l_chkimage, l_key_check, l_logfile, ccdsec, chkccdsec
    char checkimg, imgsec, chk_section, tmpout, l_section

    int secx1, secx2, secy1, secy2, ccdx1, ccdx2, ccdy1, ccdy2
    int chkccdx1, chkccdx2, chkccdy1, chkccdy2
    int chkdatax1, chkdatax2, chkdatay1, chkdatay2
    int spos1, spos2, commapos

    bool l_verbose, debug, sxcheck, sycheck

    l_inimage = inimage
    l_section = section
    l_chkimage = chkimage
    l_key_check = key_check
    l_logfile = logfile
    l_verbose = verbose

    # Temporary files
    tmpout = mktemp ("tmpout")

    # Default values
    debug = no
    status = 0
    out_chkimage = ""
    sxcheck = yes
    sycheck = yes

    secx1 = 0
    secx2 = 0
    secy1 = 0
    secy2 = 0

    # Test the log file
    if ((l_logfile == "") || (stridx(" ", l_logfile) > 0)) {
       l_logfile = gsaoi.logfile
       if ((l_logfile == "") || (stridx (" ", l_logfile) > 0)) {
          l_logfile = "gsaoi.log"
          printlog ("WARNING - GADIMSCHECK: Both gastat.logfile and \
              gsaoi.logfile fields are empty", l_logfile, l_verbose)
          printlog ("                     Using default file gsaoi.log", \
              l_logfile, l_verbose)
       }
    }

    printlog ("", l_logfile, l_verbose)
    printlog ("------------------------------------------------------------", \
        l_logfile, l_verbose)
    gemdate (zone="local")
    printlog ("GADIMSCHK -- Started: "//gemdate.outdate//"\n", \
        l_logfile, l_verbose)

    printlog ("GADIMSCHK: Input parameters...", l_logfile, l_verbose)
    printlog ("    inimage     = "//l_inimage, l_logfile, l_verbose)
    printlog ("    section     = "//l_section, l_logfile, l_verbose)
    printlog ("    chkimage    = "//l_chkimage, l_logfile, l_verbose)
    printlog ("    key_check   = "//l_key_check, l_logfile, l_verbose)
    printlog ("", l_logfile, l_verbose)

    # Check the input images exist
    gemextn (l_inimage, check="exists", process="expand", index="", \
        extname="", extversion="", ikparams="", omit="kernel,section", \
        replace="", outfile=tmpout, logfile=l_logfile, glogpars="", verbose=no)

    if (gemextn.status != 0) {
        printlog ("ERROR- GADIMSCHECK: GEMEXTN returned a non-zero status"//\
            ". Exiting.", l_logfile, verbose+)
        goto crash
    } else {
        if (gemextn.fail_count > 0 || gemextn.count == 0) {
            # This is to catch an input error to gemextn
            printlog ("ERROR - GADIMSCHECK: GEMEXTN returned a count of 0 \
                and a fail_count > 0 for "//l_inimage//". Exiting", \
                l_logfile, verbose+)
            goto crash
        }
    }

    # Parse image section
    gemsecchk (l_inimage, l_section, logfile=l_logfile, verbose=debug)

    if (gemsecchk.status != 0) {
        printlog ("ERROR - GADIMSCHK: GEMSECCHK reurned a non-zero status", \
            l_logfile, verbose+)
        goto crash
    } else {
        imgsec = gemsecchk.out_imgsect
    }

    printlog ("GADIMSCHK: Input image: \""//l_inimage//"\" section: \""//\
        imgsec//\
        "\"", l_logfile, l_verbose)

    if (l_chkimage != "" && stridx (" ", l_chkimage) == 0) {

        delete (tmpout, verify-, >& "dev$null")

        gemextn (l_chkimage, check="exists", process="expand", index="", \
            extname="", extversion="", ikparams="", omit="kernel,section", \
            replace="", outfile=tmpout, logfile=l_logfile, glogpars="", \
            verbose=no)

        if (gemextn.status != 0) {
            printlog ("ERROR - GADIMSCHECK: GEMEXTN returned a non-zero \
                status. Exiting.", l_logfile, verbose+)
            goto crash
        } else {
            if (gemextn.fail_count > 0 || gemextn.count == 0) {
                # This is to catch an input error to gemextn
                printlog ("ERROR - GADIMSCHECK: GEMEXTN returned a count of 0 \
                    and a fail_count of >0 for "//l_chkimage//" . Exiting", \
                    l_logfile, verbose+)
                goto crash
            }
        }

        # Read in the input name
        head (tmpout, nlines=1) | scan (checkimg)

        printlog ("GADIMSCHK: Check Image: "//checkimg, l_logfile, l_verbose)


    } else {
        checkimg = ""
    }

    # Chcek if there is anything to do
    if (imgsec == "[*,*]" && checkimg == "") {
        printlog ("GADIMSCHK: There is nothing to do. Exiting", \
            l_logfile, l_verbose)
        goto clean
    }

    # Set the sxcheck and sycheck variables by parsing the section again
    spos1 = stridx ("*",imgsec)
    spos2 = strldx ("*",imgsec)

    if (spos1 != spos2) {
        # Two *s supplied
        sxcheck = no
        sycheck = no
    } if (spos1 == spos2 && spos1 > 0) {
        # One * supplied
        commapos = stridx (",",imgsec)
        if (spos1 < commapos) {
            # X range is *
            sxcheck = no
            sycheck = yes

            print (substr(imgsec,strldx("[",imgsec),\
                strlen(imgsec))) | scanf ("[*,%d:%d]", secy1, secy2)
        } else {
            # Y range is *
            sxcheck = yes
            sycheck = no

            print (substr(imgsec,strldx("[",imgsec),\
                strlen(imgsec))) | scanf ("[%d:%d,*]", secx1, secx2)
        }
    } else {
        # There is a section to check
        print (substr(imgsec,strldx("[",imgsec),\
            strlen(imgsec))) | \
            scanf ("[%d:%d,%d:%d]", secx1, secx2, secy1, secy2)
    }

    # Determine the region of teh check image to use.
    if (checkimg != "") {

        # Read the two CCDSECs
        keypar (l_inimage, l_key_check, silent+)
        if (!keypar.found) {
            printlog ("ERROR - GADIMSCHK: Keyword "//\
                "CCDSEC not found in "//l_inimage,
                l_logfile, verbose+)
            goto crash
        }
        ccdsec = keypar.value

        keypar (checkimg, l_key_check, silent+)
        if (!keypar.found) {
            printlog ("ERROR - GADIMSCHK: Keyword "//\
                "CCDSEC not found in "//checkimg,
                l_logfile, verbose+)
            goto crash
        }
        chkccdsec = keypar.value

        # Check that they are the same
        if (ccdsec == chkccdsec) {
            chk_section = imgsec

            if (debug) {
                printlog ("____CCDSECs match: "//ccdsec, \
                    l_logfile, verbose+)
            }
        } else {
            # Make sure that the img ccdsec is contained within CHK
            # CCDSEC
            print (ccdsec) | scanf ("[%d:%d,%d:%d]", ccdx1, ccdx2,
                ccdy1, ccdy2)
            print (chkccdsec) | scanf ("[%d:%d,%d:%d]", chkccdx1,
                chkccdx2, chkccdy1, chkccdy2)

            if (ccdx1 < chkccdx1 || ccdx2 > chkccdx2 || \
                ccdy1 < chkccdx1 || ccdy2 > chkccdy2) {
                printlog ("ERROR - GADIMSCHK: CCDSEC for "//checkimg//\
                    " does not contain entirity of CCDSEC for "//l_inimage, \
                    l_logfile, verbose+)
                goto crash
            } else {

                keypar (checkimg, "DATASEC", silent+)
                if (!keypar.found) {
                    printlog ("ERROR - GADIMSCHK: Keyword "//\
                        "DATASEC not found in "//checkimg,
                        l_logfile, verbose+)
                    goto crash
                }
                print (keypar.value) | scanf ("[%d:%d,%d:%d]", \
                    chkdatax1, chkdatax2, chkdatay1, chkdatay2)

                chkdatax1 += ccdx1 - chkccdx1
                chkdatax2 += ccdx2 - chkccdx2
                chkdatay1 += ccdy1 - chkccdy1
                chkdatay2 += ccdy2 - chkccdy2

                if (debug) {
                    printlog ("____chkdata -     "//\
                        "x1: "//chkdatax1//" "//\
                        "    x2: "//chkdatax2//" "//\
                        "    y1: "//chkdatay1//" "//\
                        "    y2: "//chkdatay2, \
                        l_logfile, verbose+)
                }

                # Determine new section for chk file
                if (sxcheck) {
                    # X range given in section
                    chkdatax2 = chkdatax1 + secx2
                    chkdatax1 += secx1
                }

                if (sycheck) {
                    # X range given in section
                    chkdatay2 = chkdatax1 + secy2
                    chkdatay1 += secy1
                }

                if (debug) {
                    printlog ("____New chkdata - "//\
                        "x1: "//chkdatax1//" "//\
                        "    x2: "//chkdatax2//" "//\
                        "    y1: "//chkdatay1//" "//\
                        "    y2: "//chkdatay2//"\n", \
                        l_logfile, verbose+)
                }

                chk_section = "["//chkdatax1//":"//chkdatax2//","//\
                    chkdatay1//":"//chkdatay2//"]"
            }

        } # End of CCDSEC checks
        out_chkimage = l_chkimage//chk_section

        printlog ("GADIMSCHK: Output Check Image: "//out_chkimage, \
            l_logfile, l_verbose)
    }

    goto clean

    #--------------------------------------------------------------------------

crash:

    # Exit with error subroutine
    status = 1

clean:

    delete (tmpout, verify-, >& "dev$null")

    # Print finish time
    gemdate(zone="local")
    printlog ("\nGADIMSCHK -- Finished: "//gemdate.outdate, \
        l_logfile, l_verbose)

    if (status == 0) {
        printlog ("\nGADIMSCHK -- Exit status: GOOD", l_logfile, l_verbose)

    } else {
        printlog ("\nGADIMSCHK -- Exit status: ERROR", l_logfile, l_verbose)
    }

    if (status != 0) {
        printlog ("          -- Please read the logfile \""//l_logfile//\
                "\" for more information.", l_logfile, l_verbose)
    }

    printlog ("------------------------------------------------------------", \
        l_logfile, l_verbose)

end

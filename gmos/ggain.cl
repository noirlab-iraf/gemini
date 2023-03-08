# Copyright(c) 2002-2017 Association of Universities for Research in Astronomy, Inc.

procedure ggain (images)

# Multiply extensions by their gains to help remove chip-chip variations
#
# Version   Feb 28, 2002 BM  v1.3 release
#           Aug 19, 2002 IJ  bugfix when gain not found
#           Aug 27, 2002 IJ  attempt to fool proof for the hedit/Linux bug
#           Sept 20, 2002 IJ v1.4 release
#           Oct 10, 2002 BM  change high gain criteria for GMOS-S
#           Aug 25, 2003 KL  IRAF2.12 - new parameter, addonly, in hedit

string  images      {"",prompt="Input images"}
string  gaindb      {"default",prompt="Database with gain data"}
string  logfile     {"",prompt="Logfile"}
string  key_gain    {"GAIN",prompt="Header keyword for gain (e-/ADU)"}
string  key_ron     {"RDNOISE",prompt="Header keyword for readout noise"}
real    gain        {2.2,prompt="Gain value to use if gain not found in db"}
real    ron         {3.5,prompt="Readout noise value to use if ron not found in db"}
bool    fl_update   {yes,prompt="Update headers?"}
bool    fl_mult     {no,prompt="Multiply by the gains?"}
string  gainout     {"",prompt="Gains found"}
string  ronout      {"",prompt="Read noise found"}
string  sci_ext     {"SCI",prompt="Name of science extension"}
string  var_ext     {"VAR",prompt="Name of variance extension"}
bool    verbose     {yes,prompt="Verbose?"}
int     status      {0,prompt="Exit status (0=good)"}
struct  *scanfile   {"",prompt="Internal use only"}

begin

    string  l_images, l_gaindb, l_logfile
    string  l_key_gain, l_key_ron
    string  l_sci_ext, l_var_ext, pjustify
    string  l_readmode, l_gainmode, dum, readmode_test_val
    string  img, inimages[200], filelist, gaindbname, sciext, ampname
    real    l_gain, l_rdnoise, bias, gainmode_test_val
    int     i, j, n, nbad, nsci, amplen
    bool    l_verbose, l_fl_mult, l_fl_update
    bool    intdbg
    struct  amp, pstruct

    l_images = images
    l_logfile = logfile
    l_sci_ext = sci_ext
    l_var_ext = var_ext
    l_gaindb = gaindb
    l_key_gain = key_gain
    l_key_ron = key_ron
    l_verbose = verbose
    l_fl_mult = fl_mult
    l_fl_update = fl_update
    if (l_fl_mult) {
        l_fl_update = yes
    }

    cache ("imgets", "gimverify", "keypar", "gemdate")
    status = 0
    intdbg = no

    # Test the logfile:
    if (l_logfile == "" || stridx(" ",l_logfile)>0) {
        l_logfile = gmos.logfile
        if (l_logfile == "" || stridx(" ",l_logfile)>0) {
            l_logfile = "gmos.log"
            printlog ("WARNING - GGAIN: Both ggain.logfile and gmos.logfile \
                fields are empty", l_logfile, yes)
            printlog ("                    Using default file gmos.log",
                l_logfile, yes)
        }
    }

    gemdate ()
    printlog ("-------------------------------------------------------------\
        -------------------", l_logfile, verbose=l_verbose)
    printlog ("GGAIN -- Started: "//gemdate.outdate, \
        l_logfile, verbose=l_verbose)

    # Test input images
    if (l_images == "" || stridx(" ",l_images)>0) {
        printlog ("ERROR - GGAIN: No input images given", l_logfile, yes)
        goto error
    }

    # Check if user defined gaindb exists
    if (l_gaindb != "default") {    # user defined gaindb
        gaindbname = osfn(l_gaindb)
        if (access(gaindbname) == no) {
            printlog ("ERROR - GGAIN: Gain database not found", l_logfile, yes)
            goto error
        }
    }

    # Test for @filelist
    if ((substr(l_images,1,1) == "@") \
        && !access(substr(l_images,2,strlen(l_images))) ) {

        printlog ("ERROR - GGAIN: Input list "//\
            substr(l_images,2,strlen(l_images))//" does not exist",
            l_logfile, yes)
        goto error
    }

    filelist = mktemp("tmpfiles")

    files (l_images, sort-, > filelist)
    scanfile = filelist
    n = 0
    nbad = 0
    while (fscan(scanfile, img) != EOF) {
        n = n+1
        gimverify (img)
        if (gimverify.status != 0) {
            nbad = nbad+1
        } else {
            inimages[n] = gimverify.outname//".fits"
        }
        keypar (img//"[0]","GGAIN",silent+)
        if (keypar.found==yes) {
            printlog ("WARNING - GGAIN: "//img//" already gain corrected",
                l_logfile, l_verbose)
            nbad = nbad+1
        }
    }
    scanfile = ""

    if (nbad>0) {
        printlog ("ERROR - GGAIN: "//nbad//" errors or images already gain \
            corrected.", l_logfile, yes)
        goto error
    }


    for (i=1; i<=n; i+=1) {

        # Set test values for gainmode and readmode for each type of detector
        keypar (inimages[i]//"[0]", "DETTYPE", silent+)
        if (keypar.value == "SDSU II CCD") { # Current EEV CCDs
            readmode_test_val = "1000" # Fast readout
            gainmode_test_val = 3.0
        } else if (keypar.value == "SDSU II e2v DD CCD42-90") {
            # New e2vDD CCDs
            readmode_test_val = "1000" ##M # Fast readout
            gainmode_test_val = 3.0 ##M
        } else if (keypar.value == "S10892" || keypar.value == "S10892-N") { #Hamamatsu CCDs
            readmode_test_val = "4000" ##M # Fast readout
            gainmode_test_val = 3.0 ##M
        }

        # Select gaindb if set to 'default'
        if (l_gaindb == "default") {
            ggdbhelper (inimages[i]//"[0]", logfile=l_logfile)
            if (ggdbhelper.status != 0) goto error
            gaindbname = osfn(ggdbhelper.gaindb)
            if (intdbg) print ("DEBUG: gaindbname="//gaindbname)
            if (access(gaindbname) == no) {
                printlog ("ERROR - GGAIN: Gain database not found",
                    l_logfile, yes)
                goto error
            }
        }

        imgets (inimages[i]//"[0]", "NSCIEXT", >& "dev$null")
        nsci = int(imgets.value)
        if (nsci==0) {

            gemextn (inimages[i], check="exists,mef", process="expand", \
                index="1-", extname="", extversion="", ikparams="", omit="", \
                replace="", outfile="dev$null", logfile=l_logfile, \
                glogpars="", verbose=yes)

            if (gemextn.status != 0) {
                printlog ("ERROR - GGAIN: GEMEXTN returned a non-zero \
                    status", l_logfile, verbose+)
                goto error
            } else if (gemextn.count <= 0) {
                printlog ("ERROR - GGAIN: GEMEXTN could not determine the \
                    number of extensions", l_logfile, verbose+)
                goto error
            } else {
                nsci = gemextn.count
            }
        }

        # Get the readmode and set the gains
        l_readmode = "slow"
        l_gainmode = "low"
        keypar (inimages[i]//"[0]", "AMPINTEG", silent+)
        if (keypar.found==no) {
            printlog ("WARNING - GGAIN: Cannot find CCD readmode for image "//\
                inimages[i], l_logfile, verbose+)
            printlog ("                    Assuming slow read",
                l_logfile, verbose+)
        } else if (keypar.value == readmode_test_val) {
                l_readmode = "fast"
        }

        if (l_sci_ext != "") {
            sciext = l_sci_ext//","
        } else {
            sciext = ""
        }

        imgets (inimages[i]//"["//sciext//"1]", "GAIN", >& "dev$null")
        if(real(imgets.value) > gainmode_test_val) {
            l_gainmode = "high"
        } else if(real(imgets.value) < gainmode_test_val) {
            l_gainmode = "low"
        }

        for (j=1; j<=nsci; j+=1) {
            keypar (inimages[i]//"["//sciext//j//"]", "AMPNAME", silent+)

            if (!keypar.found) {
                printlog (inimages[i]//"["//sciext//j//"]", \
                    l_logfile, verbose+)
                goto error
            }
            ampname = keypar.value

            l_gain = 0.0
            l_rdnoise = 0.0
            # Added " " to gainmode match as slow low is matched by both
            # readmode and gainmode calls
            match (l_readmode, gaindbname, stop-) | \
                match (" "//l_gainmode, "STDIN", stop-) | \
                match (ampname, "STDIN", stop-) | \
                scan (dum, dum, l_gain, l_rdnoise, bias, amp)

            amplen = int(strlen(amp) + 3)
            if (strstr("left",amp) > 0) {
                pjustify = "  "
                amplen += 1
            } else {
                pjustify = " "
            }

            if (j==1) {
                printf("%-"//\
                    str(int(strlen(inimages[i])+strlen(l_sci_ext)+5))//\
                    "s %-"//str(amplen)//"s %-6s %s\n", \
                    "File[extension]:", "AMPNAME", "GAIN", "READNOISE") | \
                    scan (pstruct)
                 printlog (pstruct, l_logfile, l_verbose)
            }

            printf ("%s[%s%d]: \"%s\""//pjustify//"%6.3f %5.2f\n", \
                inimages[i], sciext, j, amp, l_gain, l_rdnoise) | \
                scan (pstruct)
            printlog (pstruct, l_logfile, l_verbose)

            if (l_gain==0.0) {
                printlog ("WARNING - GGAIN: gain and readnoise not found "//\
                    " in "//gaindbname//" for this mode", l_logfile, l_verbose)
                l_gain = gain
                l_rdnoise = ron
                printlog ("Using values from parameter", l_logfile, l_verbose)
            }
            if (j==1) {
                gainout = str(l_gain)
                ronout = str(l_rdnoise)
            } else {
                gainout = gainout//" "//str(l_gain)
                ronout = ronout//" "//str(l_rdnoise)
            }
            # update header
            # Until v1.12, this included a workaround for an "hedit Linux bug"
            # where keypar was used to check that the update worked correctly
            # and repeat the hedit if not. However, that keypar step started
            # crashing occasionally on fast machines due to an outdated FITS
            # cache that IRAF thought was still valid because the time stamp
            # resolution is too coarse. So now we use gemhedit/nhedit instead and
            # remove that keypar step.
            if (l_fl_update) {
                gemhedit (inimages[i]//"["//sciext//j//"]", l_key_gain,
                    l_gain, "Amplifier gain", delete-)
                gemhedit (inimages[i]//"["//sciext//j//"]", l_key_ron,
                    l_rdnoise, "Readout noise", delete-)
            }
            # multiply science extension by the gain
            if (l_fl_mult) {
                imarith (inimages[i]//"["//sciext//j//"]", "*", l_gain,
                    inimages[i]//"["//sciext//j//",overwrite+]", verb-)
                gemhedit (inimages[i]//"["//sciext//j//"]", l_key_gain,
                    1.0, "", delete-)
                gemhedit (inimages[i]//"["//sciext//j//"]", "GAINMULT",
                    l_gain, "Gain multiplication")

                gemdate ()
                gemhedit (inimages[i]//"[0]", "GGAIN", gemdate.outdate,
                    "UT Time stamp for GGAIN", delete-)
                gemhedit (inimages[i]//"[0]", "GEM-TLM", gemdate.outdate,
                    "UT Last modification with GEMINI", delete-)
            }
            # If a variance plane, multiply by gain squared
            if (imaccess(inimages[i]//"["//l_var_ext//","//j//"]")) {
                if (l_fl_update) {
                    gemhedit (inimages[i]//"["//l_var_ext//","//j//"]",
                        l_key_gain, l_gain, "", delete-)
                    gemhedit (inimages[i]//"["//l_var_ext//","//j//"]",
                        l_key_ron, l_rdnoise, "", delete-)
                }
                if (l_fl_mult) {
                    imarith (inimages[i]//"["//l_var_ext//","//j//"]", "*",
                        l_gain**2,
                        inimages[i]//"["//l_var_ext//","//j//",overwrite+]",
                        verb-)
                    gemhedit (inimages[i]//"["//l_var_ext//","//j//"]",
                        l_key_gain, 1.0, "", delete-)
                    gemhedit (inimages[i]//"["//l_var_ext//","//j//"]",
                        "GAINMULT", l_gain, "Gain multiplication ")
                }
            }
        }
    }

    goto clean

error:
    status=1

clean:
    scanfile = ""

    delete (filelist, verify-, >& "dev$null")
    # close log file

    gemdate ()
    printlog ("\nGGAIN -- Finished: "//gemdate.outdate, \
        l_logfile, verbose=l_verbose)

    if (status == 0) {
        printlog ("\nGGAIN -- Exit staus: GOOD\n", l_logfile, l_verbose)
    } else {
        printlog ("\nGGAIN -- Exit staus: ERROR\n", l_logfile, l_verbose)
    }
    printlog ("---------------------------------------------------------------\
        -----------------", l_logfile, l_verbose )
end

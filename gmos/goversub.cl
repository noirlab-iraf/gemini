# Copyright(c) 2011-2017 Association of Universities for Research in Astronomy, Inc.

procedure goversub (inimages)

# This is task to do a quick and very rough overscan subtraction on the fly
# with no fitting done. Only imstat is called to determine the overscan. This
# is not available to the public!

# Inputs can be raw or gprepare'd. The task does checks for the existance of a
# BIASSEC and whether it is has aleardy been overscan subtracted (by this or
# gireduce), this includes BIAS and DARK subtraction

# Output statuses are:
#    0 = goversub performed the requested task
#    1 = goversub could not overscan correct due it already being overscan
#        subtracted
#    2 = Error within script; or not overscan subtracted but already trimmed

##M TODO:
##M  Adidtional checks for BIASIM / DARKIM?
##M Update nclip parameter

char    inimages     {prompt="Input images to be roughly overscan subtracted"}
char    outimages    {"", prompt="Output filenames"}
char    outprefix    {"gos", prompt="Prefix for output images"}
char    rawpath      {"", prompt="Path to raw data set"}
char    pixstat      {"midpt", prompt="Statistic to be calculated by imstatistics to use to determine the rough BIAS level"}
char    bias_type    {"calc", enum="calc|static|default", prompt="Bias value type to determine."}
char    sci_ext      {"SCI", prompt="Name of science extension(s)"}
char    var_ext      {"VAR", prompt="Name of variance extension(s)"}
char    dq_ext       {"DQ", prompt="Name of data quality extension(s)"}
char    mdf_ext      {"MDF", prompt="Mask definition file extension name"}
int     nbiascontam  {5, prompt="Number of columns to strip from BIASSEC on the side closest to the DATASEC",min=0}
char    key_biassec  {"BIASSEC", prompt="Header keyword for overscan section"}
char    statsec      {"default", prompt="Section of BIASSEC to use. Coordinates are relative to the BIAS section. In the format: [%d:%d,%d:%d]"}
bool    calc_only    {no, prompt="Only calculate the overscan value"}
char    oscan_val    {"", prompt="Task output; calculated overscan value"}
char    gaindb       {"default",prompt="Database with gain data"}
char    logfile      {"", prompt="Logfile for this task"}
bool    verbose      {yes, prompt="Verbose, yes or no"}
int     status       {0, prompt="Exit status, 0=GOOD, performed as expected; 1=good, input previously overscan subtracted; 2=ERROR"}
struct  *scanfile    {prompt="For internal use only"}
struct  *scanfile2   {prompt="For internal use only"}

begin

    int  i, numext, junk, atlocation, num_inimages
    int  bx1, bx2, by1, by2, ubx1, ubx2, uby1, uby2, file_counter
    int  ox1, ox2, oy1, oy2, xbin, ybin
    int  ROWS_TO_REMOVE
    char outname, tmpinlist, currentimage, extn, listfile, extn_list, filename
    char inname, l_sciext, l_output, imstatsec, inlist, outlist, outextn
    char outimage_files[200], todel_list, rphend, pathtest, tmpextn
    char l_oscan_val, tmplist, tmp2inlist, l_key_ccdsec, l_key_ccdsum
    char l_readmode, readmode_chk, gaindbname, l_gainmode, l_ampname, dummy
    real bias, gainmode_check
    bool mdf_exists, fl_useprefix, istrimmed, profile, ishamamatsu
    struct sdate

    # Declare local variables for user input
    char l_inimages
    char l_outimages
    char l_outprefix
    char l_rawpath
    char l_pixstat
    char l_bias_type
    char l_sci_ext
    char l_var_ext
    char l_dq_ext
    char l_mdf_ext
    int  l_nbiascontam
    char l_key_biassec
    char l_statsec
    char l_gaindb
    bool l_calc_only
    char l_logfile
    bool l_verbose

    # Read user input
    l_inimages = inimages
    l_outprefix = outprefix
    l_outimages = outimages
    l_rawpath = rawpath
    l_pixstat = pixstat
    l_bias_type = bias_type
    l_sci_ext = sci_ext
    l_var_ext = var_ext
    l_dq_ext = dq_ext
    l_mdf_ext = mdf_ext
    l_nbiascontam = nbiascontam
    l_key_biassec = key_biassec
    l_logfile = logfile
    l_statsec = statsec
    l_gaindb = gaindb
    l_calc_only = calc_only
    l_verbose = verbose

    oscan_val = ""
    l_oscan_val = oscan_val

    # Initiate / define default values for local variables
    status = 0
    numext = 0
    profile = no
    ishamamatsu = no
    ROWS_TO_REMOVE = 48
    l_key_ccdsec = "CCDSEC"
    l_key_ccdsum = "CCDSUM"

    # Test the logfile here to start logging:
    if ((l_logfile == "") || (stridx(" ", l_logfile) > 0)) {
        l_logfile = gemlocal.logfile
        if ((l_logfile == "") || (stridx (" ", l_logfile) > 0)) {
            l_logfile = "gmos.log"
            printlog ("GOVERSUB WARNING - Both goversub.logfile and "//\
                "gmos.logfile fields are empty.", l_logfile, verbose+)
            printlog ("                   Using default file gmos.log", \
                l_logfile, verbose+)
        }
    }

    if (profile) {
        date ("+%H:%M:%S.%N%n") | scan (sdate)
        printlog ("____GOVERSUB Started: "//sdate, l_logfile, verbose+)
    }

    # Print start time
    gemdate (zone="local")
    printlog ("\nGOVERSUB - Started: "//gemdate.outdate//"\n", \
        l_logfile, l_verbose)

    # Create temporary file to store lists to be deleted at the end
    todel_list = mktemp("tmptodel_list")

    # Create temporary file for input files
    tmpinlist = mktemp("tmptmpinlist")
    tmp2inlist = mktemp("tmp2tmpinlist")
    print (tmpinlist, > todel_list)
    print (tmp2inlist, >> todel_list)
    inlist = mktemp("tmpinlist")
    print (inlist, >> todel_list)

    #### Inlist ####
    # Check the inimages input
    # Is it a list?
    atlocation = stridx("@", l_inimages)
    if (atlocation > 0) {
        listfile = substr(l_inimages,(atlocation + 1),strlen(l_inimages))
        # Check to se if file exists
        if (!access(listfile)) {
            printlog ("GOVERSUB ERROR - Cannot access file: \""//listfile//\
                "\". Exiting",\
                l_logfile, verbose+)
            status = 2
            goto clean
        }
    }

    # Expand input files into a file to verify files and one to use with
    # output
    sections (l_inimages, > inlist)
    sections (l_inimages, > tmp2inlist)

    #### Rawpath ####
    if ((l_rawpath != "") && (stridx(" ", l_rawpath) == 0)) {
        rphend = substr(l_rawpath,strlen(l_rawpath),strlen(l_rawpath))
        if (rphend == "$") {
            show (substr(l_rawpath,1,strlen(l_rawpath)-1)) | scan (pathtest)
            rphend = substr(pathtest,strlen(pathtest),strlen(pathtest))
        }
        if (rphend != "/") {
            l_rawpath = l_rawpath//"/"
        }
        if (!access(l_rawpath)) {
            printlog ("ERROR - GOVERSUB: Cannot access rawpath: "//l_rawpath, \
                l_logfile, verbose+)
            goto crash
        }
    }

    # Set up parameters for reading input files
    scanfile = tmp2inlist
    filename  = ""
    num_inimages = 0

    # Loop over input list of images
    while (fscan(scanfile,filename) != EOF) {
        num_inimages += 1

        # Verify image
        gimverify (l_rawpath//filename)
        if (gimverify.status > 0) {
             # Image bad/doesn't exist goto to clean
             printlog ("", l_logfile, l_verbose)
             printlog ("GOVERSUB ERROR - File: "//filename//\
                 " does not exist or is bad", l_logfile, verbose+)
             # Change status to bad
             status = 2
             goto clean
        } else {
           print (gimverify.outname//".fits", >> tmpinlist)
        }
    } # End of loop verifying the input files
    delete (tmp2inlist, verify-, >& "dev$null")

    # Check the calc_only flag
    if (num_inimages > 1 && l_calc_only) {
        printlog ("WARNING - GOVERSUB: Cannot only calculate overscan when \
            more than one image is input.\n"//\
            "                    Resetting calc_only to no", \
            l_logfile, verbose+)
        l_calc_only = no
    }

    #### Outlist ####
    # Set default flag for use of the prefix
    fl_useprefix = yes
    outlist = mktemp("tmpoutlist")
    print (outlist, >> todel_list)

    # Check for an output list
    if (l_outimages != "" && (stridx(" ", l_outimages) == 0)) {
        atlocation = stridx("@", l_outimages)
        if (atlocation > 0) {
            listfile = substr(l_outimages,(atlocation + 1),strlen(l_outimages))
            # Check to se if file exists
            if (!access(listfile)) {
                printlog ("ERROR - GOVERSUB: Cannot access file: \""//\
                    listfile//"\". Exiting",\
                    l_logfile, verbose+)
                status = 2
                goto clean
            }
        }

        if (atlocation > 0) {
            # Output is a list
            sections (l_outimages, > outlist)
        } else {
            # Images not in an @list
            files (l_outimages, sort=no, > outlist)
        }

        # Check the number of output files is equal to the number of inputs
        sections ("@"//outlist, option="nolist")
        if (num_inimages == int(sections.nimages)) {
            fl_useprefix = no
        } else {
            delete (outlist, verify-, >& "dev$null")
        } # End of if statement checking number of in and out files
    } # End of output list checks

    # Check if prefix is to be used
    if (fl_useprefix) {
        # Remove any directories in the name
        scanfile = tmpinlist
        while (fscan(scanfile,inname) != EOF) {
            fparse (inname)
            print (l_outprefix//str(fparse.root)//str(fparse.extension), \
                >> outlist)
        }
        scanfile = ""
    }

    # Set up parameters for reading output file list
    scanfile = outlist
    filename = ""
    i =  1

    # Loop over outlist of images to populate outimage_files and check they
    # don't exist
    while (fscan(scanfile,filename) != EOF) {
        # Verify image
        gimverify (filename)
        if (gimverify.status != 1) {
             # If image exits
             printlog ("", l_logfile, verbose+)
             printlog ("GOVERSUB - ERROR: File: "//filename//\
                 " already exists", l_logfile, verbose+)
             # Change status to bad
             status = 2
             goto clean
        }
        outimage_files[i] = gimverify.outname//".fits"
        i += 1
    } # End of while loop verifying and populating outimage_files[]

    # Check l_statec is not an empty string
    if (l_statsec == "" || (stridx (" ",l_statsec) > 0)) {
        printlog ("ERROR - GOVERSUB: satsec is an empty string.", \
            l_logfile, verbose+)
        status = 2
        goto clean
    }

    #### Hard work ####

    # Loop over the in list again to start work on it
    scanfile = tmpinlist
    filename  = ""
    file_counter = 0
    while (fscan(scanfile,filename) != EOF) {
        file_counter += 1

        filename = filename

        # inname
        fparse (filename)
        inname = fparse.root//fparse.extension

        l_output = outimage_files[file_counter]

        extn_list = mktemp("tmpextn_list")

        numext = INDEF

        print (extn_list, >> todel_list)

        ishamamatsu = no
        # Set check values for readmode and gainmode and ishamamatsu flag
        keypar (filename//"[0]", "DETTYPE", silent+)
        if (keypar.found) {
            #KL check if this is valid for gmos N.
            if (keypar.value == "S10892" || keypar.value == "S10892-N") {
                # Hamamatsu CCDs

                # Fast mode
                readmode_chk = "4000"

                # Gain check
                gainmode_check = 3.

            } else if (keypar.value == "SDSU II CCD") {
                # Current EEV2 CCDs

                # Fast mode
                readmode_chk = "1000"

                # Gain check
                gainmode_check = 3.

            } else if (keypar.value == "SDSU II e2v DD CCD42-90") { \
                # e2vDD CCDs

                # Fast mode
                readmode_chk = "1000"

                # Gain check
                gainmode_check = 3.

            } else {
                # Check images...
                keypar (l_image//"[0]", "INSTRUME", silent+)
                if (keypar.found) {
                    if (strstr("GMOS",str(keypar.value) == 0)) {
                        printlog ("ERROR - GOVERSUB: "//filename//\
                            "is not a GMOS image.", l_logfile, verbose+)
                        goto crash
                    } else {
                        printlog ("ERROR - GOVERSUB: Unregognized DETTYPE \
                            value in "//l_image//"[0]", \
                            l_logfile, verbose+)
                        goto crash
                    }
                } else {
                    printlog ("ERROR - GOVERSUB: INSTRUME keyword missing \
                        from "//l_image//"[0]", l_logfile, verbose+)
                    goto crash
                }
            }
        } else {
            printlog ("ERROR - GOVERSUB: DETTYPE keyword missing in "//\
                l_image//"[0]", l_logfile, verbose+)
            goto crash
        }

        # Read the number of SCI extentions
        keypar (filename//"[0]", "NSCIEXT", silent+)
        if (!keypar.found) {

            keypar (filename//"[0]", "NEXTEND", silent+)
            if (!keypar.found) {
                gemextn (filename, check="exists,mef", process="expand", \
                    index="1-", extname="", extversion="", ikparams="", \
                    omit="section", replace="", outfile=extn_list, \
                    logfile=l_logfile, glogpars="", verbose=no)

                if (gemextn.status != 0) {
                    printlog ("ERROR - GOVERSUB: GEMEXTN returned an non-zero \
                        status.", l_logfile, verbose+)
                    status = 2
                    goto crash
                } else if (gemextn.count > 0) {

                    numext = INDEF
                    # Only want extensions with no extname / version
                    match ("][", extn_list, stop=yes, print_file_name=no, \
                        metacharacters=no, > tmp2inlist)

                    count (extn_list) | scan (numext)

                    if (isindef(numext)) {
                        printlog ("ERROR - GOVERSUB: Cannot determine number \
                            of input extensions.", l_logfile, verbose+)
                        status = 2
                        goto crash
                    } else {
                        delete (extn_list, verify-, >& "dev$null")
                    }

                    l_sci_ext = ""
                }
            } else {
                numext = int (keypar.value)
            }

        } else {
            numext = int (keypar.value)
        } # End of checking extnsion numbers etc.

        if (isindef(numext)) {
            printlog ("ERROR - GOVERSUB: Cannot determine number of \
                science extensions",  l_logfile, verbose+)
            status = 2
            goto crash
        }

        # Find any MDFs
        # Check for VARDQ planes exist
        gemextn (filename, check="exists,mef", process="expand", \
            index="", extname=l_mdf_ext, extversion="1-", ikparams="", \
            omit="section,index", replace="", outfile=extn_list, \
            logfile=l_logfile, glogpars="", verbose=no)

        if (gemextn.status != 0) {
            printlog ("ERROR - GOVERSUB: GEMEXTN returned an non-zero \
                status.", l_logfile, verbose+)
            status = 2
            goto crash
        } else if (gemextn.count > 0) {
            printlog ("GOVERSUB: MDF found", \
                l_logfile, l_verbose)
            printlog ("", l_logfile, l_verbose)
        }

        # Check for VARDQ planes exist
        gemextn (filename, check="exists,mef", process="expand", \
            index="", extname=l_var_ext//","//l_dq_ext, extversion="1-", \
            ikparams="", omit="section,index", replace="", outfile=extn_list, \
            logfile=l_logfile, glogpars="", verbose=no)

        if (gemextn.status != 0) {
            printlog ("ERROR - GOVERSUB: GEMEXTN returned an non-zero \
                status.", l_logfile, verbose+)
            status = 2
            goto crash
        } else if (gemextn.count > 0) {
            printlog ("GOVERSUB: Variance and DQ planes found", \
                l_logfile, l_verbose)
            printlog ("", l_logfile, l_verbose)
        }

        # End of checks

        l_sciext = l_sci_ext
        if (l_sciext != "") {
            l_sciext = l_sciext//","
        }

        # Check the files haven't been gmosaic'ed or gtile'd
        keypar (filename//"[0]", "GMOSAIC", silent+)
        if (keypar.found) {
            printlog ("ERROR - GOVERSUB: Cannot overscan subtract: "//inname//\
                " it has been gmosaiced.", l_logfile, verbose+)
            status = 2
            goto crash
        }
        keypar (filename//"[0]", "GTILE", silent+)
        if (keypar.found) {
            printlog ("ERROR - GOVERSUB: Cannot overscan subtract: "//inname//\
                " it has been gtiled.", l_logfile, verbose+)
            status = 2
            goto crash
        }

        keypar (filename//"[0]", "GSCUT", silent+)
        if (keypar.found) {
            printlog ("ERROR - GOVERSUB: Cannot overscan subtract: "//inname//\
                " it has been gscut.", l_logfile, verbose+)
            status = 2
            goto crash
        }

##M What about BIASIM/DARKIM?

        istrimmed = no
        # Check if it has been trimmed too
        keypar (filename//"[0]", "TRIMMED", silent+)
        if (keypar.found) {
            istrimmed = yes
        }

        # Determine what type of bias to return (if trimmed may need to reset
        # to static
        if (l_bias_type == "static" || istrimmed) {

            if (l_gaindb == "default") {
                ggdbhelper (filename//"[0]", logfile=l_logfile)
                if (ggdbhelper.status == 0) {
                    gaindbname = ggdbhelper.gaindb
                } else {
                    printlog ("ERROR - GOVERSUB: GGDBHELPER returned a \
                        non-zero status.", l_logfile, verbose+)
                    goto crash
                }
            } else {
                gaindbname = l_gaindb
            }

            if (!access(gaindbname)) {
                printlog ("ERROR - GOVERSUB: Cannot access gain database \
                    file "//gaindbname, l_logfile, verbose+)
                goto crash
            }

            # default ampread mode is slow
            keypar (filename//"[0]", "AMPINTEG", silent+)
            if (!keypar.found) {
                printlog ("ERROR - GOVERSUB: AMPINTEG keyword not found in "//\
                    filename//"[0].", l_logfile, verbose+)
                status = 2
                goto crash
            } else if (keypar.value == readmode_chk) {
                l_readmode = "fast"
            } else {
                l_readmode = "slow"
            }
        }

        # Copy the phu to the outfile
        if (!l_calc_only) {
            fxcopy (filename//"[0]", l_output, groups="0", verbose=no)
        }

        # Loop over the sci extensions
        for (i = 1; i <= numext; i += 1) {
            # Set up the current image extn and file name
            extn = "["//l_sciext//str(i)//"]"
            currentimage = filename//extn
            outextn = "["//l_sciext//str(i)//",append,dupname+]"

            bias = INDEF

            # Check if we can perform the action
            keypar (currentimage, "OVERSCAN", silent+)
            if (keypar.found) {
                printlog ("WARNING - GOVERSUB: Cannot overscan subtract: "//\
                    inname//" it has been overscan corrected previously."//\
                    "\n                    "//\
                    "No output files will be written to disk.", \
                    l_logfile, verbose+)

                status = 1
                if (!l_calc_only) {
                    goto clean
                } else {
                    bias = real(keypar.value)
                    goto UPDATE
                }
            }

            keypar (currentimage, "OVERSUB", silent+)
            if (keypar.found) {
                printlog ("WARNING - GOVERSUB: Cannot overscan subtract: "//\
                    inname//\
                    "\n                    "//\
                    "It has been overscan corrected"// \
                    " previously by GOVERSUB."//\
                    "\n                    "//\
                    "No output files will be written to disk.", \
                    l_logfile, verbose+)
                status = 1
                if (!l_calc_only) {
                    goto clean
                } else {
                    bias = real(keypar.value)
                    goto UPDATE
                }
            }

            # If the other keywords don't exist and it's trimmed change to
            # static bias_type
            if (isindef(bias) && istrimmed && i == 1) {
                printlog ("WARNING - GOVERSUB: "//filename//" is not "//\
                    "overscan subtracted but is trimmed."//\
                    "\n                    Setting bias_type=static", \
                    l_logfile, verbose+)
                l_bias_type = "static"
            }

            # To get to here all is good and can perform the overscan
            # calculation and subtract

            # Static bias
            if (l_bias_type == "static") {
                # Read the gainmode and then the static bias level

                keypar (currentimage, "GAINORIG", silent+)
                if (!keypar.found) {
                    keypar (currentimage, "GAIN", silent+)
                    if (!keypar.found) {
                        printlog ("ERROR - GOVERSUB: Cannot read GAIN or \
                            GAINORIG keywords in "//currentimage, \
                            l_logfile, verbose+)
                        goto crash
                    }
                }

                if (real(keypar.value) > gainmode_check) {
                    l_gainmode = "high"
                } else {
                    l_gainmode = "low"
                }

                keypar (currentimage, "AMPNAME", silent+)
                if (!keypar.found) {
                    printlog ("ERROR - GOVERSUB: AMPNAME keyword not found \
                        in "//currentimage, l_logfile, verbose+)
                    goto crash
                }
                l_ampname = keypar.value

                match (l_readmode, gaindbname, stop-) | \
                    match (l_gainmode, "STDIN", stop-) | \
                    match (l_ampname, "STDIN", stop-) | \
                    scan (dummy, dummy, dummy, dummy, bias, dummy)

                # Go to the update section
                goto UPDATE
            }

            # Dynamic bias value

            # Read the bias section
            keypar (currentimage, l_key_biassec, silent+)
            if (!keypar.found) {
                printlog ("ERROR - GOVERSUB: "//l_key_biassec//\
                    " not found in: "//inname//extn, \
                    l_logfile, verbose+)
                status = 2
                goto crash
            }

            print (keypar.value) | scanf ("[%d:%d,%d:%d]", bx1, bx2, by1, by2)

            if (l_statsec == "default") {
                # Biassec is on the left
                if (bx1 == 1) {
                    bx2 -= l_nbiascontam
                    bx1 += 1
                } else {
                    # Biassec is on the right
                    bx1 += l_nbiascontam
                    bx2 -= 1
                }

                # Need to remove bottom 48 rows... if needed - MS
                # KL Mar 2017.  This is never run since ishamamatsu is
                # never set to yes.  What the heck?  Don't want to break
                # things right now, so leaving it for now.
                if (ishamamatsu) {
                    keypar (currentimage, l_key_ccdsec, silent+)
                    if (keypar.found) {
                        # Only adjust sections if CCDSEC is found
                        # This can also be done using DETSEC
                        print (keypar.value) | scanf ("[%d:%d,%d:%d]", \
                                                       ox1, ox2, oy1, oy2)
                        keypar(currentimage, l_key_ccdsum, silent+)
                        if (keypar.found) {
                            # Now adjust number of rows to remove by binning
                            print (keypar.value) | scanf ("%d %d", xbin, ybin)

                            # Update the lower limit on the bias section
                            if (oy1 <= ROWS_TO_REMOVE) {
                                by1 += ((ROWS_TO_REMOVE - (oy1 - 1)) / ybin)
                            }
                        }
                    }
                }
            } else {

                print (l_statsec) | scanf ("[%d:%d,%d:%d]", \
                    ubx1, ubx2, uby1, uby2)

                bx2 = bx1 + ubx2
                bx1 += ubx1
                by2 = by1 + uby2
                by1 += uby1

            }

            imstatsec = "["//str(bx1)//":"//str(bx2)//","//str(by1)//":"//\
                str(by2)//"]"

            printlog ("GOVERSUB: Imstatistics section used of input file: "//\
                imstatsec, l_logfile, l_verbose)

            bias = INDEF

            # Run imstatistics on the 'new' biassec.
            imstatistics (currentimage//imstatsec, fields=l_pixstat, \
                lower=INDEF, upper=INDEF, nclip=0, lsigma=3.0, usigma=3.0, \
                binwidth=0.1, format=no, cache=yes) | scan (bias)

            # Check the output was sensible
            if (isindef(bias)) {
                 printlog ("ERROR - GOVERSUB: imstatistics didn't return a \
                     value. Exiting", l_logfile, verbose+)
                 status = 2
                 goto crash
            }

UPDATE:
            printlog ("GOVERSUB: Bias level for "//inname//extn//" is: "//\
                bias, l_logfile, l_verbose)

            if (i == 1) {
                l_oscan_val = str(bias)
            } else {
                l_oscan_val = l_oscan_val//" "//str(bias)
            }

            if (!l_calc_only) {
                # Subtract the bias
                imexpr ("a - b", l_output//outextn, currentimage, bias, \
                    dims="auto", intype="auto", outtype="real", \
                    refim="auto", bwidth=0, btype="nearest", bpixval=0., \
                    rangecheck=yes, verbose=no, exprdb="none")

                gemhedit (l_output//outextn, "OVERSUB", bias, \
                    "Calculated overscan value", delete-)
            }
        }

        ## Copy over any other extnsions present

        if (access(extn_list) && !l_calc_only && status != 1) {
            scanfile2 = extn_list
            while (fscan(scanfile2,tmpextn) != EOF) {
                if (strstr(l_mdf_ext,tmpextn) > 0) {
                    tcopy (tmpextn, l_output, verbose-)
                } else {
                    imcopy (tmpextn, l_output//"[append]", verbose+)
                }
            }
            delete (extn_list, verify-, >& "dev$null")
        }

        if (!l_calc_only) {
            gemdate (zone="UT")
            gemhedit (l_output//"[0]", "GOVERSUB", gemdate.outdate, \
                "UT Time stamp for GOVERSUB", delete-)
            gemhedit (l_output//"[0]", "GEM-TLM", gemdate.outdate, \
                "UT Last modification with GEMINI", delete-)
        }

    } # End of while loop over filenames

    oscan_val = l_oscan_val

    # All good to get to here, so skip crash clean up
    goto clean

crash:

    if (status == 0) {
        status = 2
    }

    oscan_val = ""

    # Clean up any temporary files that may exist
    for (i = 1; i <= file_counter; i += 1) {
        if (access(outimage_files[i])) {
            delete (outimage_files[i], verify-, >& "dev$null")
        }
    }

clean:

    # Delete any temporay files that definitely exist
    if (access(todel_list)) {
        scanfile = todel_list
        while (fscan(scanfile,filename) != EOF) {
            if (access(filename)) {
                delete (filename, verify-, >& "dev$null")
            }
        }

        delete (todel_list, verify-, >& "dev$null")
    }

    scanfile = ""

    # Only delete the output files if status == 1 (headers are written to that
    # file
    if (status == 1) {
        for (i = 1; i <= file_counter; i += 1) {
            delete (outimage_files[i], verify-, >& "dev$null")
        }
    }

    # Print finish time
    gemdate (zone="local")
    printlog ("\nGOVERSUB - Finished: "//gemdate.outdate, \
        l_logfile, l_verbose)

    if (status == 0) {
        printlog ("GOVERSUB - Exit status: GOOD\n", l_logfile, l_verbose)
    } else if (status == 1) {
        printlog ("GOVERSUB - Exit status: GOOD - Though files already \
            overscan subtracted\n", l_logfile, verbose+)
    } else {
        printlog ("GOVERSUB - Exit status: ERROR\n", l_logfile, l_verbose)
    }

    if (profile) {
        date ("+%H:%M:%S.%N%n") | scan (sdate)
        printlog ("____GOVERSUB Finished: "//sdate, l_logfile, verbose+)
    }

end

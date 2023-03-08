# Copyright(c) 2005-2015 Association of Universities for Research in Astronomy, Inc.
#
# CKINPUT -- Check input and output and return a list.
#
# Version: Sept 2, 2005 FV
#          May 2,  2008, NZ   Copy from nifs$nfcheck to this new name

procedure ckinput ()

string  input       {"", prompt = "Input files to check"}
string  output      {"", prompt = "Output files to check"}
string  prefix      {"", prompt = "Output prefix"}
file    outlist     {"", prompt = "Return list of files"}
string  name        {"CKINPUT", prompt = "Task name"}
string  dependflag  {"*PREPAR*", prompt = "Dependent processing flag"}
string  procflag    {"", prompt = "Previous processing flag"}
string  sci_ext     {"SCI", prompt = "Science extension"}
file    logfile     {"STDOUT", prompt = "Logfile"}
bool    verbose     {yes, prompt = "Verbose?"}
bool    vverbose    {no, prompt = "Very verbose?"}

# Output parameters
int     nfiles      {prompt = "Number of files"}
int     nver        {prompt = "Number of extension versions"}
int     status      {prompt = "Return status"}

struct  *scanfile1  {"", prompt="Internal use only"}
struct  *scanfile2  {"", prompt="Internal use only"}

begin
    bool    runbefore
    int     nbad, naxis, junk
    file    img, phu
    file    tmpfile, tmpphu, tmpsci, tmpin, tmpout
    struct  line

    # Local parameter variables
    string l_input
    file l_output
    string l_prefix
    file l_outlist
    string l_name
    string l_dependflag
    string l_procflag
    string l_sci_ext
    string l_logfile
    bool l_verbose
    bool l_vverbose
    int l_nver
    int l_nfiles

    # Set local input parameter variables
    l_input = input
    l_output = output
    l_prefix = prefix
    l_outlist = outlist
    l_name = name
    l_dependflag = dependflag
    l_procflag = procflag
    l_sci_ext = sci_ext
    l_logfile = logfile
    l_verbose = verbose
    l_vverbose = vverbose

    # Define temporary files
    tmpfile = mktemp ("tmpfile")
    tmpphu = mktemp ("tmpphu")
    tmpsci = mktemp ("tmpsci")
    tmpin = mktemp ("tmpin")
    tmpout = mktemp ("tmpout")

    # Cache parameter files
    cache ("gemextn")

    # Initialize
    status = 1
    l_nver = INDEF

    # Check for output list.
    if (access (l_outlist)) {
        printf ("ERROR - %s: Output list %s exists.\n", l_name, l_outlist) | \
            scan (line)
        printlog (line, l_logfile, verbose+)
        goto clean
    }

    # Expand and verify input
    if (l_vverbose) print ("expansion and verification of input")
    gemextn (l_input, check="mef,exists", process="none", index="",
        extname="", extversion="", ikparams="", omit="extension, kernel",
        replace="", outfile=tmpfile, logfile="", glogpars="",
        verbose=l_verbose)
    if (gemextn.count == 0 || gemextn.fail_count != 0) {
        printf ("ERROR - %s: Problem with input.\n", l_name) | scan (line)
        printlog (line, l_logfile, verbose+)
        goto clean
    }
    l_nfiles = gemextn.count

    # Check PHUs
    gemextn ("@" // tmpfile, check="mef,exists", process="append",
        index="0", extname="", extversion="", ikparams="",
        omit="extension, kernel", replace="", outfile=tmpphu,
        logfile="", glogpars="", verbose=l_verbose)
    if (gemextn.count != l_nfiles || gemextn.fail_count != 0) {
        printf ("ERROR - %s: Problem with PHUs.\n", l_name) | scan (line)
        printlog (line, l_logfile, verbose+)
        goto clean
    }

    # Make input list of images not previously run and check keywords.
    scanfile1 = tmpphu
    scanfile2 = tmpfile
    l_nfiles = 0
    nbad = 0
    while (fscan (scanfile1, phu) != EOF) {
        junk = fscan (scanfile2, img)
        if (l_vverbose) print ("file: " // phu)
        runbefore = no

        # Check previous processing flag.
        if (l_procflag != "") {
            line = ""
            hselect (phu, l_procflag, yes) | scan (line)
            if (line != "") {
                printf ("WARNING - %s: Image %s has been run\n",
                    l_name, img) | scan (line)
                printlog (line, l_logfile, verbose+)
                printf ("                       through %s before.\n",
                    l_name) | scan (line)
                printlog (line, l_logfile, verbose+)
                runbefore = yes
            }
        }

        # Check dependent flag.
        if (l_dependflag != "") {
            line = ""
            hselect(phu, l_dependflag, yes) | scan(line)
            if (line == "") {
                printf ("WARNING - %s: Image %s has not been prepared (%s).\n",
                    l_name, img, l_dependflag) | scan (line)
                printlog (line, l_logfile, verbose+)
                nbad += 1
            }
        }

        # Add to list of not previously run input files.
        if (!runbefore) {
            print (img, >> tmpin)
            l_nfiles += 1
        }
    }

    # Check science extensions.
    gemextn ("@" // tmpfile, check="mef,exists", process="expand",
        index="", extname=l_sci_ext, extversion="1-", ikparams="",
        omit="extension,index", replace="", outfile=tmpsci,
        logfile="", glogpars="", verbose=l_verbose)
    if (gemextn.count == 0 || gemextn.fail_count != 0) {
        printf ("ERROR - %s:  Problems with science data.\n", l_name) | \
            scan (line)
        printlog (line, l_logfile, verbose+)
        goto clean
    }

    # Check dimensions.
    scanfile1 = tmpsci
    while (fscan (scanfile1, img) != EOF) {
        hselect (img, "NAXIS", yes) | scan (naxis)
        if (nscan() != 1) {
            printf ("WARNING - %s: Data in %s missing NAXIS.\n",
                l_name, img) | scan (line)
            printlog (line, l_logfile, verbose+)
            nbad += 1
        } else if (naxis != 2) {
            printf ("WARNING - %s: Data in %s are not 2-D.\n", l_name, img) | \
                scan (line)
            printlog (line, l_logfile, verbose+)
            nbad += 1
        }
    }

    # Check for empty file list
    if (l_nfiles == 0) {
        printf ("ERROR - %s:  No input images to process.\n", l_name) | \
            scan (line)
        printlog (line, l_logfile, verbose+)
        goto clean
    }

    # Exit if problems found with input files
    if (nbad > 0) {
        printf ("ERROR - %s: %d bad image(s)\n", l_name, nbad) | scan (line)
        printlog (line, l_logfile, verbose+)
        goto clean
    }

    # Check that all images have the same number of SCI extensions.
    if (l_vverbose) print ("checking version numbers")
    scanfile1 = tmpin
    while (fscan (scanfile1, img) != EOF) {
        if (isindef(l_nver)) {
            gemextn (img, check="exists", process="expand", index="",
                extname=l_sci_ext, extversion="1-",
                ikparams="", omit="", replace="", outfile="dev$null",
                logfile="", glogpars="", verbose-)
            l_nver = gemextn.count
        } else
            gemextn (img, check="exists", process="expand", index="",
                extname=l_sci_ext, extversion="1-" // l_nver,
                ikparams="", omit="", replace="", outfile="dev$null",
                logfile="", glogpars="", verbose-)
        if (l_nver == 0 || gemextn.count != l_nver || gemextn.fail_count != 0) {
            printf ("ERROR - %s: Bad or missing science data in %s.\n",
                l_name, img) | scan (line)
            printlog (line, l_logfile, verbose+)
            goto clean
        }
    }

    # Log input files.
    printlog ("Using input files:", l_logfile, l_verbose)
    if (l_verbose) type (tmpin)
    type (tmpin, >> l_logfile)

    # Expand and verify output
    if (l_output != "" || l_prefix != "") {
        if (l_vverbose) print ("expansion and verification of output")
        gemextn (l_output, check="absent", process="none", index="",
            extname="", extversion="", ikparams="",
            omit="extension,kernel", replace="", outfile=tmpout,
            logfile="", glogpars="", verbose=l_verbose)
        if (gemextn.fail_count != 0) {
            printf ("ERROR - %s: Problems with output.\n", l_name) | scan (line)
            printlog (line, l_logfile, verbose+)
            goto clean
        }
        if (gemextn.count == 0) {
            printf ("%%^%%%s%%@%s\n", l_prefix, tmpin) | scan (line)
            gemextn (line, check="absent", process="none", index="",
                extname="", extversion="", ikparams="",
                omit="kernel,exten", replace="", outfile=tmpout,
                logfile="", glogpars="", verbose=l_verbose)
            if (gemextn.count == 0 || gemextn.fail_count != 0) {
                printf ("ERROR - %s: No or incorrectly formatted output \
                    files\n", l_name) | scan (line)
                printlog (line, l_logfile, verbose+)
                goto clean
            }
        }

        # Check number of input and output images
        if (l_nfiles != gemextn.count) {
            printf ("ERROR - %s: Different number of input and output files\n",
                l_name) | scan (line)
            printlog (line, l_logfile, verbose+)
            goto clean
        }

        # Log output files.
        printlog ("Using output files:", l_logfile, l_verbose)
        if (l_verbose) type (tmpout)
        type (tmpout, >> l_logfile)

        if (l_vverbose) print ("constructing file list")
        joinlines (tmpin, tmpout, output=l_outlist, delim=" ", \
            missing="Missing", maxchar=161, shortest=yes, verbose=no)
    } else
        rename (tmpin, l_outlist)

    # Completed successfully
    status = 0
    nver = l_nver
    nfiles = l_nfiles

clean:
    scanfile1 = ""
    scanfile2 = ""
    delete (tmpfile, verify-, >& "dev$null")
    delete (tmpphu, verify-, >& "dev$null")
    delete (tmpsci, verify-, >& "dev$null")
    delete (tmpin, verify-, >& "dev$null")
    delete (tmpout, verify-, >& "dev$null")
end

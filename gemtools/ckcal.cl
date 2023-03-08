# Copyright(c) 2001-2012 Association of Universities for Research in Astronomy, Inc.
#
# CKCAL -- Check lamp and S-dist files and return an updated list.
#
# Version: Sept 19, 2005 FV
#          May 2,   2008 NZ   Copy from nifs$nsckcal with new name

procedure ckcal () 

string    lampfiles = ""        {prompt = "Lamp files to check"}
string    sdistfiles = ""       {prompt = "S-distortion files to check"}
file      outlist = ""          {prompt = "List to update"}
int       nfiles = 0            {prompt = "Number of files to calibrate"}
int       nver = 0              {prompt = "Number of extension versions"}
string    sci_ext = "SCI"       {prompt = "Science extension"}
file      database = ""         {prompt = "Database"}
string    name = "CKCAL"        {prompt = "Task name"}
file      logfile = "STDOUT"    {prompt = "Logfile"}
bool      verbose = yes         {prompt = "Verbose?"}
bool      vverbose = no         {prompt = "Very verbose?"}

int       nlamps                {prompt = "Number of lamp files"}
int       nsdist                {prompt = "Number of S-distortion files"}
int       status                {prompt = "Return status"}

struct    *scanfile

begin
    int       nbad, version, junk
    file      img, idname, tmplamps, tmpsdist, tmpout, lamp, sdist
    string    sec, l_database = ""
    struct    line

    # Define temporary files
    tmplamps = mktemp ("tmplamps") 
    tmpsdist = mktemp ("tmpsdist") 
    tmpout = mktemp ("tmpout") 

    # Cache parameter files
    cache ("gemextn")

    # Initialize
    nlamps = 0
    nsdist = 0
    status = 1

    # Set the local variables
    junk = fscan (database, l_database)

    # Check for output list.
    if (!access (outlist)) {
        printf ("ERROR - %s: Output list %s does not exist.\n",
            name, outlist) | scan (line)
        printlog (line, logfile, verbose+)
        goto clean
    }

    if (l_database == "") {
        junk = fscan (gnirs.database, l_database)
        if (l_database == "") {
            l_database = "database"
            printlog ("WARNING - NSTELLURIC: Both ckcal.database \
                and gnirs.database are", logfile, verbose+)
            printlog ("                      undefined.  Using " \
                // l_database, logfile, verbose+)
        }
    }

    # Remove any trailing slashes from the database name
    if (strlen (l_database) > 1 && \
        substr (l_database, strlen (l_database), strlen (l_database)) \
        == "/") {
        l_database = substr (l_database, 1, strlen (l_database) - 1)
    }

    # Check lamp calibration files.
    if (lampfiles != "") {
        if (vverbose) print ("checking lamp files")
        gemextn (lampfiles, check="", process="none", index="",
            extname="", extversion="", ikparams="",
            omit="extension, kernel", replace="", outfile=tmplamps,
            logfile="", glogpars="", verbose=verbose)
        if (gemextn.fail_count != 0) {
            printf ("ERROR - %s: Problem with lamp file(s).\n",
                name) | scan (line)
            printlog (line, logfile, verbose+)
            goto clean
        }
        nlamps = gemextn.count
        if (nlamps == 0) {
            printf ("WARNING - %s: lampfile(s) specified but \
                not found (%s)\n", name, lampfile) | scan (line)
            printlog (line, logfile, verbose+)
        }
    }

    if (nlamps>0) {
        if (vverbose) print ("checking number of lamp files")
            if (nlamps != 1 && nlamps != nfiles) {
                printf ("ERROR - %s: Number of wavelength calibration files \
                    is %d.  It must be either 1\n", name, nlamps) | scan (line)
                printlog (line, logfile, verbose+)
                printf ("                     or equal to the number of \
                    input spectra (%d).\n", nfiles) | scan (line)
                printlog (line, logfile, verbose+)
                goto clean
            }
            
            scanfile = tmplamps
            nbad = 0
            while (fscan (scanfile, img) != EOF) {
                for (version = 1; version <= nver; version = version + 1) {
                    sec = "_" // sci_ext // "_" // version // "_"
                    idname = l_database // "/id" // img // sec
                    if (!access (idname) ) {
                        printf ("WARNING - %s: wavelength calibration file \n",
                            name) | scan (line)
                        printlog (line, logfile, verbose+)
                        printf ("                       %s does not exist.\n",
                            idname) | scan (line)
                        printlog (line, logfile, verbose+)
                        nbad = nbad + 1
                    }
                }
            }
            if (nbad == nver) {
                printf ("ERROR - %s: no wavelength information.\n", name) |
                    scan (line)
                printlog (line, logfile, verbose+)
                goto clean
            }
        }
        
        # Check sdist calibration files.
        nsdist = 0
        if (sdistfiles != "") {
            if (vverbose) print ("checking sdist file")
            gemextn (sdistfiles, check="", process="none", index="",
                extname="", extversion="", ikparams="",
                omit="extension, kernel", replace="", outfile=tmpsdist,
                logfile="", glogpars="", verbose=verbose)
            if (gemextn.fail_count != 0) {
                printf ("ERROR - %s: Problem with sdist file(s).\n", name) |
                    scan (line)
                printlog (line, logfile, verbose+)
                goto clean
            }
            nsdist = gemextn.count
            if (nsdist == 0) {
                printf ("WARNING - %s: sdistfile(s) specified but \
                    not found (%s).\n", name, sdisfiles) | scan (line)
                printlog (line, logfile, verbose+)
            }
        }

        if (nsdist>0) {
            if (vverbose) print ("checking number of sdist files")
            if (nsdist != 1 && nsdist != nfiles) {
                printf ("ERROR - %s: Number of S-distortion calibration files \
                    is %d.  It must be either 1\n", name, nsdist) | scan (line)
                printlog (line, logfile, verbose+)
                printf ("                     or equal to the number of \
                    input spectra (%d).\n", nfiles) | scan (line)
                printlog (line, logfile, verbose+)
                goto clean
            }
            
            scanfile = tmpsdist
            nbad = 0
            while (fscan (scanfile, img) != EOF) {
                for (version = 1; version <= nver; version = version + 1) {
                    sec = "_" // sci_ext // "_" // version // "_"
                    idname = l_database // "/id" // img // sec
                    if (!access (idname) ) {
                    printf ("WARNING - %s: S-distortion calibration file \n",
                        name) | scan (line)
                    printlog (line, logfile, verbose+)
                    printf ("                       %s does not exist.\n",
                        idname) | scan (line)
                    printlog (line, logfile, verbose+)
                    nbad = nbad + 1
                }
            }
        }
        if (nbad == nver) {
            printf ("ERROR - %s: no S-distortion information.\n", name) |
                scan (line)
            printlog (line, logfile, verbose+)
            goto clean
        }
    }

    # Check if there is any calibration data.
    if (nlamps == 0 && nsdist == 0) {
        printf ("ERROR - %s: No wavelength or S-distortion \
            transformations defined\n", name) | scan (line)
        printlog (line, logfile, verbose+)
        goto clean
    }

    # At this point tmplamps contains the wavelength calibration files and
    # tmpsdist contains the s-distortion calibration files.  Add input
    # lamps and sdist to the output list, matched with the correct images.

    if (vverbose) print ("constructing main file list")
    if (nlamps > 1)
        joinlines (outlist, tmplamps, output=tmpout, delim=" ",
            missing="Missing", maxchar=161, shortest+, verbose-) 
    else {
        if (nlamps == 0)
            lamp = "none"
        else
            head (tmplamps, nlines=1) | scan (lamp)
        scanfile = outlist
        while (fscan (scanfile, line) != EOF)
            print (line // " " // lamp, >> tmpout) 
    }
    delete (outlist, verify-) 
    rename (tmpout, outlist, field = "all") 

    if (nsdist > 1)
        joinlines (outlist, tmpsdist, output=tmpout, delim=" ",
        missing="Missing", maxchar=161, shortest+, verbose-) 
    else {
        if (nsdist == 0)
            sdist = "none"
        else
            head (tmpsdist, nlines=1) | scan (sdist)
        scanfile = outlist
        while (fscan (scanfile, line) != EOF)
        print (line // " " // sdist, >> tmpout) 
    }
    delete (outlist, verify-) 
    rename (tmpout, outlist, field = "all") 
    
    # Completed successfully
    status = 0

clean:
    scanfile = ""
    delete (tmplamps, verify-, >>& "dev$null") 
    delete (tmpsdist, verify-, >>& "dev$null") 
    delete (tmpout, verify-, >& "dev$null") 
end

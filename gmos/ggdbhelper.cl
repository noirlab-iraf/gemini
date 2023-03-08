# Copyright(c) 2006-2017 Association of Universities for Research in Astronomy, Inc.

procedure ggdbhelper (inimage)

# GGDBHELPER - Using date info in image, select appropriate gaindb file.

char    inimage {prompt="Select gaindb for this image"}
char    gaindb  {"", prompt="Output. Name of gain db for input image"}
char    logfile {"", prompt="Logfile for task"}
int     status  {0, prompt="Exit status (0=good)"}

begin

    char    l_inimage = ""
    char    l_logfile = ""

    char    keyfound, keyfound1, keyfound2, tmpimg
    int     obsdate
    int     junk

    bool    dateobs_fixed=no  # GMOS-N Ham CCD commissioning.

    junk = fscan (inimage, l_inimage)
    l_logfile = logfile

    status = 0

    if (l_logfile == "" || (stridx(" ",l_logfile) > 0)) {
        l_logfile = gmos.logfile
        if (l_logfile == "" || (stridx(" ",l_logfile) > 0)) {
            l_logfile = "gmos.log"
            printlog ("WARNING - GDBHELPER: Both gsreduce.logfile and \
                gmos.logfile are empty", l_logfile, verbose+)
            printlog ("                    Using default; gmos.log",
                l_logfile, verbose+)
        }
    }

    keyfound1 = ""
    keyfound2 = ""

    fparse (l_inimage)
    if (fparse.cl_index == -1 && fparse.ksection == "") {
        tmpimg = fparse.directory//fparse.root//fparse.extension//"[0]"
        printlog ("WARNING - GGDBHELPR: Using PHU for "//l_inimage//": "//\
            tmpimg, l_logfile, verbose+)
        l_inimage = tmpimg
    }

    hselect (l_inimage, "DATE-OBS", yes) | scan (keyfound1)
    ##M Hack for early Hamamatsu commisioning data on GMOS-N: to be removed!!!
    ##M DATE-OBS is 1970-... for early commisioning data use DATE instead
    if (keyfound1 != "" && !dateobs_fixed) {
        keypar (l_inimage, "DETTYPE", silent+)
        if (keypar.value == "S10892-N"){
            keyfound1 = ""  # force use of DATE
        }
    }
    hselect (l_inimage, "DATE", yes) | scan (keyfound2)
    if ((keyfound1 == "") && (keyfound2 == "")) {
        status = 131
        printlog ("GGDBHELPER: DATE-OBS keyword is missing from "//l_inimage, \
            l_logfile, verbose=yes)
        printlog ("GGDBHELPER: DATE keyword is missing from "//l_inimage, \
            l_logfile, verbose=yes)
        printlog ("GGDBHELPER: Unable to find date of observation.", \
            l_logfile, verbose=yes)
    } else if (keyfound1 != "") {
        keyfound = keyfound1
    } else if (keyfound2 != "") {
        keyfound = keyfound2
    }
    
    if (keyfound == "") {
        status = 1
        printlog ("GGDBHELPER: Something went really wrong.", \
            l_logfile, verbose=yes)
    } else {

        # Convert DATE-OBS to second for easier before-after date comparison
        cnvtsec (keyfound, "00:00:00") | scan (obsdate)

        # Calculate new gain data boundary dates with cnvtsec like so:
        #   cnvtsec ("2006-08-31", "00:00:00")
        #
        #   For GMOS-S: UT August 31, 2006 => 841449600 seconds

        if (obsdate < 841449600) {   # Up to UT August 31, 2006
            gaindb = "gmos$data/gmosamps-20060831.dat"
        } else if (obsdate < 1125014400) {   # Before 26 Aug, 2015
            gaindb = "gmos$data/gmosamps-20150826.dat"
        } else if (obsdate < 1172361600) {   # Before 24 Feb, 2017
            gaindb = "gmos$data/gmosamps-20170224.dat"
        } else {                    # Current default
            gaindb = "gmos$data/gmosamps.dat"
        }
    }

    printlog ("GGDBHELPER: gain database selected - "//gaindb, \
            l_logfile, verbose=no)

end

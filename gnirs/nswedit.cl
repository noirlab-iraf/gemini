# Copyright(c) 2004-2012 Association of Universities for Research in Astronomy, Inc.
#
# This code brought to you by acooke, 22 Nov 2003.

# Read/write wavelength calibration information from a binary table
#
# To keep the table compact and flexible, the search:
# - does not try to match empty values (anything in the table will match)
# - tests to see if a key is in the table and, if not, tests again
#   with the part number only (values tend to be of the form text_Gnnnn)
# - if the value is not in the table, it is not used for matching
# - if less than minmatch values match, the task fails
#   (if minmatch=0 then all must match)
#
# Note that the above only occurs when reading the table; when writing
# the values are used as given.

procedure nswedit (table)

char    table          {prompt="Table to read/write"}

char    camera         {"", prompt="Camera used"}
char    grating        {"", prompt="Grating used"}
char    filter         {"", prompt="Filter used"}
char    prism          {"", prompt="Prism used"}
char    mask           {"", prompt="Mask (slit) used"}
char    arrayid        {"", prompt="Array used"}
int     order          {-1, prompt="Spectra order"}
real    centre         {INDEF, prompt="Variable for range test (wavelength, grating angle, etc)"}
real    range          {INDEF, prompt="Range of variable values"}

char    description    {"", prompt="Additional text for table"}
int     minmatch       {0, prompt="Minimum number of keys to match (0=all)"}

bool    append         {no, prompt="Extend the table?"}
bool    overwrite      {no, prompt="Overwrite exiting table entries?"}
bool    create         {no, prompt="Create a new table?"}

char    logfile        {"", prompt = "Logfile"}
bool    verbose        {yes, prompt = "Verbose output?"}

real    lambda         {INDEF, prompt="I/O: Central wavelength (A)"}
real    delta          {INDEF, prompt="I/O: Dispersion (A/pix)"}
real    resoln         {INDEF, prompt="I/O: Resolution (units?)"}
real    cradius        {INDEF, prompt="I/O: Centering radius (pixels)"}

int     status         {0, prompt="Exit status (0 = good)"}

struct* scanin         {prompt="Internal use"}

begin
    char    l_table = ""
    char    l_camera = ""
    char    l_grating = ""
    char    l_filter = ""
    char    l_prism = ""
    char    l_mask = ""
    char    l_arrayid = ""
    int     l_order
    real    l_centre
    real    l_range
    char    l_description = ""
    int     l_minmatch
    bool    l_append, l_overwrite, l_create
    char    l_logfile = ""
    bool    l_verbose
    real    l_lambda
    real    l_delta
    real    l_resoln
    real    l_cradius

    char    tmpcdfile, tmpdatafile, tmprowtable, tmpoldtable, tmpopen
    char    tmpexists
    char    word, expr, sp[8], spaces, expr2, chlambda, devnull, chcradius
    char    lastword, in_tbltype
    int     nkeys = 7
    char    title[7], value[7]    # can't use nkeys as dimen...
    bool    exists, debug, text[7]
    struct  sdate, testopen
    int     junk, nrows, i, j, count, present, idx

    junk = fscan(table, l_table)
    junk = fscan(camera, l_camera)
    junk = fscan(grating, l_grating)
    junk = fscan(filter, l_filter)
    junk = fscan(prism, l_prism)
    junk = fscan(mask, l_mask)
    junk = fscan(arrayid, l_arrayid)
    l_order = order
    l_centre = centre
    l_range = range
    l_description = description    # may contain spaces
    l_minmatch = minmatch
    l_append = append
    l_overwrite = overwrite
    l_create = create
    junk = fscan(logfile, l_logfile)
    l_verbose = verbose
    l_lambda = lambda
    l_delta = delta
    l_resoln = resoln
    l_cradius = cradius

    title[1] = "camera"
    title[2] = "grating"
    title[3] = "filter"
    title[4] = "prism"
    title[5] = "mask"
    title[6] = "order"
    title[7] = "arrayid"
    text[1] = yes
    text[2] = yes
    text[3] = yes
    text[4] = yes
    text[5] = yes
    text[6] = no
    text[7] = yes
    value[1] = l_camera
    value[2] = l_grating
    value[3] = l_filter
    value[4] = l_prism
    value[5] = l_mask
    if (INDEF == l_order) value[6] = "INDEF"
    else value[6] = str (l_order)
    value[7] = l_arrayid

    count = 0
    status = 1
    debug = no
    if (debug) devnull = "STDOUT"
    else       devnull = "dev$null"

    tmpcdfile = mktemp ("tmpcdfile")
    tmpdatafile = mktemp ("tmpdatafile")
    tmprowtable = mktemp ("tmprowtable")//".fits"
    tmpoldtable = mktemp ("tmpoldtable")//".fits"
    tmpopen = mktemp ("tmpopen")
    tmpexists = mktemp ("tmpexists")//".fits"

    cache ("tinfo", "tabpar")

    # Start logging

    if ("" == l_logfile) {
        junk = fscan (gnirs.logfile, l_logfile)
        if (l_logfile == "") {
            l_logfile = "gnirs.log"
            printlog ("WARNING - NSWEDIT: Both gnswedit.logfile and \
                gnirs.logfile are empty.", l_logfile, verbose+)
            printlog ("                   Using default file " \
                // l_logfile // ".", l_logfile, verbose+)
        }
    }
    date | scan (sdate)
    printlog ("----------------------------------------------------\
        --------------------------", l_logfile, l_verbose)
    printlog ("NSWEDIT -- " // sdate, l_logfile, l_verbose)
    printlog (" ", l_logfile, l_verbose)

    # Validate input

    if ("" == l_table) {
        printlog ("ERROR - NSWEDIT: No table.", l_logfile, l_verbose)
        goto clean
    }
    if (debug) {
        print (l_table)
        for (i = 1; i <= nkeys; i = i + 1) {
            print (title[i] // ":" // value[i])
        }
    }

    # Create table description

    if (l_create || l_append || l_overwrite) {

        if (debug) print ("writing description")

        for (i = 1; i <= nkeys; i = i + 1) {
            if (text[i]) {
                print (title[i] // "\tCH*20\ts10\t\"\"", >> tmpcdfile)
            } else {
                print (title[i] // "\tI\ti2\t\"\"", >> tmpcdfile)
            }
        }
        print ("lo           R      f8.2    \"\"", >> tmpcdfile)
        print ("hi           R      f8.2    \"\"", >> tmpcdfile)
        print ("lambda       R      f8.2    \"\"", >> tmpcdfile)
        print ("delta        R      f8.2    \"\"", >> tmpcdfile)
        print ("resoln       R      f8.2    \"\"", >> tmpcdfile)
        print ("cradius      R      f5.2    \"\"", >> tmpcdfile)
        print ("description  CH*40  s-40    \"\"", >> tmpcdfile)
    }

    # Create table if required

    tinfo (table = l_table, >>& tmpopen)
    in_tbltype = tinfo.tbltype
    scanin = tmpopen; junk = fscan (scanin, word); exists = "can't" != word

    if (debug) print ("table exists: " // exists)
    if (no == exists) {
        if (create) {
            printlog ("WARNING - NSWEDIT: Creating new table " \
                // l_table, l_logfile, l_verbose)

            fparse (l_table)
            if (fparse.extension != ".fits") {
                printlog ("WARNING - NSWEDIT: Setting new table extension to \
                    \".fits\"", l_logfile, l_verbose)
                l_table = fparse.directory//fparse.root//".fits"//\
                    fparse.ksection
            }

            delete (files = tmpdatafile, verify-, >& devnull)
            print ("", >> tmpdatafile)
            tcreate (table = l_table, cdfile = tmpcdfile, \
                datafile = tmpdatafile, uparfile = "", nskip = 0, \
                nlines = 0, hist = yes, extrapar = 10, \
                tbltype = "default", extracol = 10)
            if (debug) {
                print ("created:")
                tprint (l_table)
            }
        } else {
            printlog ("ERROR - NSWEDIT: Table " // l_table \
                // " does not exist (and create=no)", \
                l_logfile, l_verbose)
            goto clean
        }
    } else if (in_tbltype != "fits") {

        # Check it's a fits file!
        printlog ("ERROR - NSWEDIT: Input table type must be fits", \
            l_logfile, verbose+)
        goto clean
    }

    # Check whether data aready exist

    expr = ""
    present = 0
    for (i = 1; i <= nkeys; i = i + 1) {
        if ("" != value[i] && "INDEF" != value[i]) {
            present = present + 1
            if ("" != expr) expr = expr // "&&"
            expr = expr // title[i] // "=="
            if (text[i]) expr = expr // "'"
            expr = expr // value[i]
            if (text[i]) expr = expr // "'"
        }
    }
    if (INDEF != l_centre) {
        if ("" != expr) expr = expr // "&&"
        expr = expr // "lo<=" // l_centre // "&&hi>" // l_centre
        present = present + 1
    }

    if ("" == expr) {
        printlog ("ERROR - NSWEDIT: No selection criteria.", \
        l_logfile, l_verbose)
        goto clean
    }
    if (debug) print (expr)

    if (exists) {
        tdelete (table = tmprowtable, go_ahead+, verify-, >& devnull)
        if (debug) print ("inital selection")
        tselect (intable = l_table, outtable = tmprowtable, expr = expr)

        if (debug) tprint (tmprowtable)

        tinfo (table = tmprowtable, >& devnull)
        exists = tinfo.nrows > 0

        if (debug) print ("data exist: " // exists)
    }

    # Modify table if necessary (this leaves new data in tmprowtable)

    if (l_append || l_overwrite) {

        if (exists && ! l_overwrite) {
            printlog ("ERROR - NSWEDIT: Data already exist \
                (and overwrite=no).", l_logfile, l_verbose)
            goto clean
        }
        if (no == exists && no == l_append) {
            printlog ("ERROR - NSWEDIT: Data do not exist \
                (and append=no).", l_logfile, l_verbose)
            goto clean
        }
        if (INDEF == l_delta || INDEF == l_resoln) {
            printlog ("ERROR - NSWEDIT: Dispersion or \
                resolution undefined.", l_logfile, l_verbose)
            goto clean
        }
        if (INDEF == l_lambda) {
            printlog ("WARNING - NSWEDIT: Wavelength undefined (header \
                values or default used in nsappwave).", \
                l_logfile, l_verbose)
        }
        if (INDEF == l_cradius) {
            printlog ("WARNING - NSWEDIT: Cradius undefined (parameter \
                value used in nswavelength).", \
                l_logfile, l_verbose)
        }

        # Create a new table without the data

        if (debug) print ("selecting old data")

        tinfo (table = l_table, >& devnull)
        nrows = tinfo.nrows
        tdelete (table = tmpoldtable, go_ahead+, verify-, >& devnull)
        tselect (intable = l_table, outtable = tmpoldtable, \
        expr = "! (" // expr // ")")

        # And a new table with the data

        if (debug) print ("writing new data")

        tdelete (table = tmpdatafile, go_ahead+, verify-, >& devnull)
        for (i = 1; i <= nkeys; i = i + 1) {
            if (text[i]) {
                print ("'" // value[i] // "'", >> tmpdatafile)
            } else {
                print (value[i], >> tmpdatafile)
            }
        }
        if (INDEF != l_centre && INDEF != l_range) {
            print ((l_centre - 0.5 * l_range), >> tmpdatafile)
            print ((l_centre + 0.5 * l_range), >> tmpdatafile)
        } else {
            print (-9e7, >> tmpdatafile)
            print (9e7, >> tmpdatafile)
        }
        print (l_lambda, >> tmpdatafile)
        print (l_delta, >> tmpdatafile)
        print (l_resoln, >> tmpdatafile)
        print (l_cradius, >> tmpdatafile)
        print ("'" // l_description // "'", >> tmpdatafile)
        tdelete (table = tmprowtable, go_ahead+, verify-, >& devnull)
        tcreate (table = tmprowtable, cdfile = tmpcdfile, \
            datafile = tmpdatafile, uparfile = "", nskip = 0, \
            nlines = 0, hist = yes, extrapar = 10, \
            tbltype = "default", extracol = 10)

        if (debug) {
            print ("merging old and new data")
            tprint (l_table)
            tprint (tmpoldtable)
            tprint (tmprowtable)
        }

        tdelete (table = l_table, go_ahead+, verify-, >& devnull)
        tmerge (intable = tmpoldtable // "," // tmprowtable, \
            outtable = l_table, option = "append", allcols+, \
            tbltype = "default", allrows = nrows + 10, \
            extracol = 10)

    } else {

        # Do full search, with possible truncated values etc

        tdelete (table = tmprowtable, go_ahead+, verify-, >& devnull)
        expr = ""

        for (i = 1; i <= nkeys; i = i + 1) {
            if ("" != value[i] && "INDEF" != value[i]) {
                tdelete (table = tmpexists, go_ahead+, verify-, \
                    >& devnull)
                if (text[i]) expr2 = title[i] // "=='" // value[i] // "'"
                else expr2 = title[i] // "==" // value[i]
                tselect (intable = l_table, outtable = tmpexists, \
                    expr = expr2)
                tinfo (table = tmpexists, >& devnull)
                if (tinfo.nrows == 0) {
                    if (text[i]) {
                        if (debug) print (expr2 // " failed, so abrevn")

                        # this is a little dense, i'm afraid
                        # all we're doing is finding the last
                        # word in a _ separated list
                        idx = strldx("_", value[i])
                        lastword = substr(value[i], idx+1, strlen(value[i]))
                        #print (value[i]) \
                        #| translit ("STDIN", "_", " ", delete-) \
                        #    | scan (sp[8], sp[7], sp[6], \
                        #    sp[5], sp[4], sp[3], sp[2], sp[1])
                        #idx = nscan ()
                        if (debug) print (lastword)
                        tdelete (table = tmpexists, go_ahead+, verify-, \
                            >& devnull)
                        tselect (intable = l_table, outtable = tmpexists, \
                            expr = title[i] // "=='" // lastword // "'")
                        tinfo (table = tmpexists, >& devnull)
                        if (tinfo.nrows != 0) {
                            if (debug) print ("abbrevn matched")
                            if (count > 0) expr = expr // "&&"
                            expr = title[i] // "=='" // lastword // "'"
                            count = count + 1
                        } else {
                            if (debug) print ("abbrevn failed")
                        }
                    }
                } else {
                    if (count > 0) expr = expr // "&&"
                    expr = expr // title[i] // "=="
                    if (text[i]) expr = expr // "'"
                    expr = expr // value[i]
                    if (text[i]) expr = expr // "'"
                    count = count + 1
                }
            }
        }

        if (INDEF != l_centre) {
            if (count > 0) expr = expr // "&&"
            expr = expr // "lo<=" // l_centre // "&&hi>" // l_centre
            count = count + 1
        }

        if (debug) print ("final expr: " // expr)
        #printlog ("NSWEDIT: Search term: " // expr, l_logfile, l_verbose)

        if (count < l_minmatch || (0 == l_minmatch && count != present)) {
            printlog ("NSWEDIT: Too few keys are present in \
                the table.", l_logfile, l_verbose)
            goto clean
        }

        tdelete (table = tmprowtable, go_ahead+, verify-, >& devnull)
        tselect (intable = l_table, outtable = tmprowtable, expr = expr)
    }

    # Read values

    if (debug) tprint (table = tmprowtable)

    tinfo (table = tmprowtable, >& devnull)
    if (1 < tinfo.nrows) {
        printlog ("NSWEDIT: Ambiguous lookup.", l_logfile, l_verbose)
            if (no == l_append && no == l_overwrite)
            printlog ("NSWEDIT: Search term: " // expr, \
                l_logfile, l_verbose)
        goto clean
    } else if (1 > tinfo.nrows) {
        printlog ("NSWEDIT: No matches.", l_logfile, l_verbose)
            if (no == l_append && no == l_overwrite)
            printlog ("NSWEDIT: Search term: " // expr, \
                l_logfile, l_verbose)
        goto clean
    }

    tabpar (table = tmprowtable, column = "lambda", row = 1, format-)
    lambda = INDEF
    if ("INDEF" != tabpar.value) lambda = real (tabpar.value)
    tabpar (table = tmprowtable, column = "delta", row = 1, format-)
    delta = real (tabpar.value)
    tabpar (table = tmprowtable, column = "resoln", row = 1, format-)
    resoln = real (tabpar.value)
    tabpar (table = tmprowtable, column = "cradius", row = 1, format-)
    cradius = INDEF

    if ("INDEF" != tabpar.value) cradius = real (tabpar.value)
    if (isindef (lambda)) chlambda = "INDEF"
    else chlambda = str (lambda)
    if (isindef (cradius)) chcradius = "INDEF"
    else chcradius = str (cradius)
    printlog ("NSWEDIT: " // count // " keys --> lambda: " // chlambda \
        // "; delta: " // delta // "; resoln: " // resoln \
        // "; cradius: " // chcradius, l_logfile, l_verbose)

    status = 0    # Success

clean:

    delete (tmpcdfile, verify-, >& devnull)
    delete (tmpdatafile, verify-, >& devnull)
    tdelete (tmprowtable, go_ahead+, verify-, >& devnull)
    tdelete (tmpoldtable, go_ahead+, verify-, >& devnull)
    delete (tmpopen, verify-, >& devnull)
    tdelete (tmpexists, go_ahead+, verify-, >& devnull)

    scanin = ""

    if (0 == status) {
        printlog (" ", l_logfile, l_verbose)
        printlog ("NSWEDIT exit status:  good.", l_logfile, l_verbose)
    }

end

# Copyright(c) 2004-2015 Association of Universities for Research in Astronomy, Inc.

procedure gemlist (root, range)

char    root     {prompt = "Root name (ie including S)"}
char    range    {prompt = "Integer range"}

begin

    char l_range, l_root, tmprange, tmpfile, tmpout, l_redirect
    int  last, i
    bool debug

    l_root = root
    l_range = range
    l_redirect = "dev$null"

    debug = no
    tmprange = mktemp ("tmprange")
    tmpfile = mktemp ("tmpfile")
    tmpout = mktemp ("tmptab")//".fits"

    cache ("tstat", "tcreate")

    # get the largest number in the range
    print (l_range) \
        | tokens ("STDIN", newlines-) \
        | match (",", "STDIN", stop+) \
        | translit ("STDIN", "-x", " ", collapse=no) \
        | fields ("STDIN", "1", lines="1-", quit_if_miss=no, \
        print_file_name=no, > tmpfile)

    print ("c1 i %d\n") | tcreate (table=tmpout, cdfile="STDIN", \
        datafile=tmpfile, uparfile="", nskip=0, nlines=0, nrows=0, \
        hist=yes, extrapar=5, tbltype="default", extracol=0, >>& l_redirect)

    tstat (tmpout, "c1", outtable="STDOUT", lowlim=INDEF, highlim=INDEF,
        rows="-", >>& l_redirect)

    if (debug) print (l_range // ", " // tstat.vmax // ", " // tmprange)
    if (tstat.nrows == 0) goto empty

    for (i = 1; i <= int (tstat.vmax+0.5); i += 1) {
        printf ("%s%04d\n", l_root, i, >> tmprange)
    }

    fields (tmprange, "1", lines = l_range, quit-, print-)

    delete (tmprange//","//tmpfile//","//tmpout, verify-, >>& l_redirect)

empty:
    last = 0 # keep parser happy
end

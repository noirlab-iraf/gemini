# Copyright(c) 2001-2011 Association of Universities for Research in Astronomy, Inc.
procedure wmef (input, output)

# WMEF creates a MEF file from input FITS or PLIO (.pl) images.
# Input can also be extensions in already existing MEF files.
# There is no attempt to handle other image types.
#
# Version: Aug 17, 2001 MT,BM,IJ v1.2 release
#          Jan 30, 2002 IJ       removed gsetsec, calling task needs to take care of this
#          Feb 28, 2002 IJ       v1.3 release
#          Sept 20,2002 IJ       v1.4 release
#          Jan 10, 2003 IJ       fixed ad+ -> add+ bug
#          Aug 15, 2003 KL       Port to IRAF2.12.  hedit: addonly

char    input       {prompt = "Input images in the order they should appear"}
char    output      {prompt = "Output multi-extension fits file"}
char    extnames    {"SCI,VAR,DQ", prompt = "Extension names (optional)"}
char    phu         {"dummy", prompt = "Source of PHU info"}
bool    verbose     {no, prompt = "Verbose"}
int     status      {0, prompt = "Exit status (0=good)"}
struct  *scanfile   {"", prompt = "Internal use only"}

begin

    char    l_input = ""
    char    l_output = ""
    char    l_extnames = ""
    char    l_phu = ""
    bool    l_verbose

    int     junk,i
    struct  line, tmpstr
    int     nex, n_images, n_match, n_extnames
    char    delim, extname, extn, imgin, imglist
    char    tmpin, tmpinput, tmpext, tmpfits
    bool    havephu, debug

    tmpin = mktemp ("tmpin")
    tmpinput = mktemp ("tmpfil")
    tmpext = mktemp ("tmpext")
    delim = tmpin // "," // tmpinput // "," // tmpext

    debug = no
    havephu = no
    status = 1

    junk = fscan (	input, l_input)
    if ("" != l_input) l_input = input	# can include spaces
    junk = fscan (	output, l_output)
    junk = fscan (	extnames, l_extnames)
    if ("" != l_extnames) l_extnames = extnames	# can include spaces
    junk = fscan (	phu, l_phu)
    l_verbose =	verbose


    cache ("imgets", "tinfo", "gemextn", "gemdate")


	# Check access to phu if not "dummy" or "none"

    if (l_phu != "dummy" && l_phu != "none" && l_phu != "") {
        gemextn (l_phu, check="exists,mef", process="expand", index="0", 
            extname="", extver="", ikparams="", omit="kernel,section,ext", 
            replace="", outfile="STDOUT", logfile="", glogpars="",
            verbose-) | scan (tmpstr)
        gemextn (tmpstr, check="exists", process="none", index="", extname="",
            extver="", ikparams="", omit="", replace="%%.fits%\[",
            outfile="STDOUT", logfile="", glogpars="", verbose-) | scan (line)
        if (1 != gemextn.count || 0 != gemextn.fail_count) {
            print ("ERROR - WMEF: PHU source " // l_phu \
                // " does not exist")
            goto clean
        } else {
            l_phu = line
            if (no == imaccess (l_phu)) {
                print ("ERROR - WMEF: PHU not at " // l_phu \
                    // " (code error?  PHU not FITS?)")
                goto clean
            }
            havephu = yes
        }
    }


    # Check for (absent) output, forcing fits extension

    gemextn (l_output, check="absent", process="none", index="", extname="",
        extver="", ikparams="", omit="kernel,section,ext", replace="", 
        outfile="STDOUT", logfile="", glogpars="", 
        verbose-) | scan (line)
    if (1 == gemextn.count && 0 == gemextn.fail_count) {
        l_output = line // ".fits"
        gemextn (l_output, check="absent", process="none", index="", \
            extname="", extver="", ikparams="", omit="", replace="", 
            outfile="dev$null", logfile="", glogpars="", 
            verbose-)
        if (0 != gemextn.fail_count) {
            print ("ERROR - WMEF: Output " // l_output // " exists \
                (possible code error)")
            goto clean
        }
    } else {
        print ("ERROR - WMEF: Output " // l_output // " exists")
        goto clean
    }	


    # Cheap and cheerful, but not quite reliable check for image
    # sections

    print (l_input) | match (":", "STDIN", stop-, print-, meta+) | \
        count ("STDIN") | scan (n_match)
    if (n_match > 0) {
        print ("ERROR - WMEF: Input files " // l_input \
            // " contain image sections (colons)")
        goto clean
    }


    # Expand input (removes @-files and checks for existence)

    if (debug) print ("input: " // l_input)
    gemextn (l_input, check="exists", process="none", index="", extname="",
        extver="", ikparams="", omit="", replace="", outfile=tmpinput, 
        logfile="", glogpars="", verbose-)
    if (0 != gemextn.fail_count) {
        print ("ERROR - WMEF: Problem with input files.")
        goto clean
    } else {
        n_images = gemextn.count
    }


    # For each input file
    # - if MEF, use index only
    # - force extension
    # - if .pl convert to fits

    scanfile = tmpinput
    i=0
    while (EOF != fscan (scanfile, imgin)) {
        i=i+1
	    
        if (debug) print ("preparing " // imgin)

        gemextn (imgin, check="mef", process="expand", index="0-", \
            extname="", extver="", ikparams="", omit="ext", replace="", 
            outfile="STDOUT", logfile="dev$null", glogpars="", verbose-) | \
            scan (line)
        if (1 == gemextn.count && 0 == gemextn.fail_count) {
            if (debug) print ("filtered to " // line)
            gemextn (line, check="", process="none", index="", extname="", 
                extver="", ikparams="", omit="", replace="%%.fits%\[", 
                outfile="STDOUT", logfile="", glogpars="",
                verbose-) | scan (line)
            imgin = line
            if (debug) print ("filtered to " // imgin)
        } else if (gemextn.count + gemextn.fail_count > 1) {
            # The file is a MEF and there are no extensions.  For backward,
            # compatibility, if the file is a FITS table MDF force the
            # extension to be [1].  For any other cases, the absence of
            # extension is an error.
		
            gemextn(imgin, check="table", process="expand",index="1-", 
                extname="", extver="", ikparams="", 
                omit="extension,kernel,section", replace="", 
                outfile="STDOUT", logfile="dev$null", glogpars="",
                verbose-) | scan (line)
            if ( 0 == gemextn.fail_count ) {
                gemextn (line, check="",process="none", index="", \
                    extname="", extver="", ikparams="", omit="", \
                    replace="%%.fits%\[", outfile="STDOUT", logfile="",
                    glogpars="", verbose-) | scan (line)
                    imgin = line
                    if (debug) print ("filtered to " // imgin)
            } else {
                print ("ERROR - WMEF: Some input files not completely \
                    specified.")
                goto clean
            } 
        }

        gemextn (imgin, check="", process="none", index="0", extname="", 
            extver="", ikparams="", omit="ext", replace="", 
            outfile="STDOUT", logfile="dev$null", glogpars="", 
            verbose-) | scan (line)
        if (line == imgin) {
            if (access (imgin // ".fits")) {
                imgin = imgin // ".fits"
            } else if (access (imgin // ".pl")) {
                imgin = imgin // ".pl"
            } else {
                print ("ERROR - WMEF: File " // imgin \
                    // " has no file extension (missing .fits?)")
                goto clean
            }
        }

        if (substr (imgin, strlen(imgin)-2, strlen(imgin)) == ".pl") {

            tmpfits = mktemp ("tmpfits") // ".fits"
            delim = delim // "," // tmpfits

            imcopy (imgin, tmpfits, verbose-, >& "dev$null")
            imgin = tmpfits
        }

        #print (imgin, >> tmpin)
        if (i == 1)
            imglist = imgin
        else
            imglist = imglist//","//imgin
    }
	

    # Check extension names

    files (l_extnames, > tmpext)
    n_extnames = 0
    if (access (tmpext)) count (tmpext) | scan (n_extnames)
	
    if (n_extnames != 0 && n_extnames != n_images) {
        n_extnames = 0
        delete (tmpext, ver-, >& "dev$null")
        print ("WARNING - WMEF: Number of input images and extensions \
            names differ.")
        print ("                Extensions are not named.")
    }
	


    # Get or make the PHU

    if (no == havephu) {
        fxdummyh (l_output, hdr_file = "")
    } else {
        if (debug) print ("wmef: calling fxcopy")
        fxcopy (l_phu, l_output, groups="", new_file+, verbose=l_verbose)
        if (debug) print ("wmef: fxcopy completed")
    }


	# Write the whole thing
    if (debug) print("imglist = "//imglist)

    if (l_verbose) {
        if (debug) print ("wmef: calling fxinsert")
#	    fxinsert ("@" // tmpin, l_output // "[0]", group="", \
        fxinsert (imglist, l_output // "[0]", group="", \
            verbose=l_verbose)
        if (debug) print ("wmef: fxinsert completed")
    } else {
#	    fxinsert ("@" // tmpin, l_output // "[0]", group="", \
        fxinsert (imglist, l_output // "[0]", group="", \
            verbose=l_verbose, >& "dev$null")
    }

    # Put number of extensions in the header

    if (debug) print ("adding NEXTEND")
    hedit (l_output // "[0]", "NEXTEND", n_images, add+, addonly-, \
        del-, ver-, show-, up+)


    # Name the extensions 

    if (n_extnames > 0) {
        scanfile = tmpext
        nex = 0
        while (EOF != fscan (scanfile, extname)) {
            nex = nex + 1
            extn = l_output // "[" // nex // "]"

            # What's this for?

            tinfo (extn, ttout-, >>& "dev$null")
            if (tinfo.tbltype == "fits") {
                parkey (extname, extn, "EXTNAME", add+)
                parkey ("1", extn, "EXTVER", add+)
            } else {
                hedit (extn, "EXTNAME", extname, add+, addonly-, \
                    delete-, verify-, show-, update+)
                hedit (extn, "EXTVER", "1", add+, addonly-, \
                    delete-, verify-, show-, update+)
            }
        }
    }


    # Gemini default header keywords in PHU

    gemdate ()
    gemhedit (l_output // "[0]", "WMEF", gemdate.outdate, 
        "UT Time stamp for WMEF", delete-)
    gemhedit (l_output // "[0]", "GEM-TLM", gemdate.outdate,
        "UT Time of last modification with Gemini", delete-)


    # set exit status flag
    status = 0


clean:

    # Clean up images and files

    delete (delim, ver-, >>& "dev$null")
    scanfile = ""

end

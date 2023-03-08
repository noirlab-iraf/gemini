# Copyright(c) 2004-2009 Association of Universities for Research in Astronomy, Inc.

procedure gemwcscopy (images, refimages)

char	images		{"", prompt = "List of input images"}
char	refimages	{"", prompt="List of reference images"}
bool	verbose		{no, prompt = "Verbose output?"}
char	logfile		{"", prompt = "Logfile"}

int	status		{0, prompt = "Exit status (0=good)"}
struct*	scanin1		{"", prompt = "Internal use only"}
struct*	scanin2		{"", prompt = "Internal use only"}

begin

	char	l_images = ""
	char	l_refimages = ""
	bool	l_verbose
	char	l_logfile = ""
        char    gn_logfile = ""

	int	junk, nfiles
	bool	debug, ok, manyref, first
	char	tmplog, tmpin, tmpref, tmprefphu, tmpphu, tmpheader
	char	phu, phufrom, phuto, root, reffile, refphu
	struct	line

	junk = fscan (	images, l_images)
	junk = fscan (	refimages, l_refimages)
	l_verbose = 	verbose
	junk = fscan (	logfile, l_logfile)
	
	status = 1
	debug = no

	tmplog = mktemp ("tmplog")
	tmpin = mktemp ("tmpin")
	tmpref = mktemp ("tmpref")
	tmprefphu = mktemp ("tmprefphu") // ".fits"
	tmpphu = mktemp ("tmpphu") // ".fits"
	tmpheader = mktemp ("tmpheader")

	cache ("gemextn")

        if (l_logfile == "") {
            gn_logfile = gnirs.logfile
            junk = fscan (gn_logfile, l_logfile) 
            if (l_logfile == "") {
                l_logfile = "gnirs.log"
                printlog ("WARNING - "//l_logname//": Both nscombine.logfile \
		    and gnirs.logfile are empty.", l_logfile, verbose+) 
                printlog ("                     Using default file " \
		    // l_logfile, l_logfile, verbose+) 
            }
        }


	# Check number of target and reference files

	gemextn (images, check="exists", process="none", index="", \
	    extname="", extversion="", ikparams="", omit="kernel,exten", \
	    replace="", outfile=tmpin, logfile="", glogpars="",
	    verbose=l_verbose)
	if (gemextn.fail_count != 0 || gemextn.count == 0) {
	    printlog ("ERROR - GEMWCSCOPY: Bad images.", l_logfile, verbose+)
	    goto clean
	}
	nfiles = gemextn.count

	gemextn (refimages, check="exists", process="none", index="", \
	    extname="", extversion="", ikparams="", omit="kernel,exten", \
	    replace="", outfile=tmpref, logfile="", glogpars="",
	    verbose=l_verbose)
	if (gemextn.fail_count != 0 \
	    || (gemextn.count != 1 && gemextn.count != nfiles)) {
	    printlog ("ERROR - GEMWCSCOPY: Bad reference.", \
		l_logfile, verbose+)
	    goto clean
	}
	manyref = gemextn.count > 1

	scanin1 = tmpin
	scanin2 = tmpref
	first = yes

	while (fscan (scanin1, root) != EOF) {

	    if (manyref || first) {

		if (debug) print ("making simple file copy of reference phu")

		junk = fscan (scanin2, reffile)
		
		gemextn (reffile, check="", process="expand", index="0", \
		    extname="", extversion="", ikparams="", omit="", \
		    replace="", outfile="STDOUT", logfile="dev$null",
		    glogpars="", verbose-) \
		    | scan (refphu)

		if (gemextn.fail_count != 0 || gemextn.count != 1) {
		    printlog ("ERROR - GEMWCSCOPY: Missing PHU in " \
			// reffile // ".", l_logfile, verbose+)
		    goto clean
		}

		imdelete (tmprefphu, verify-, >& "dev$null")
		tmprefphu = mktemp ("tmprefphu") // ".fits"
		delete (tmpheader, verify-, >& "dev$null")
		tmpheader = mktemp ("tmpheader")

		imheader (refphu, imlist="*.fits", longheader+, userfields+,
		    > tmpheader)
		mknoise (tmprefphu, output="", title="", ncols=2, nlines=2, \
		    header=tmpheader, background=0, gain=1, rdnoise=0, \
		    poisson-, seed=1, cosrays="", ncosrays=0, energy=30000, \
		    radius=0.5, ar=1, pa=0, comments-)

		if (debug) imhead (tmprefphu)
	    }


	    if (debug) print ("making simple file copy of phu")

	    gemextn (root, check="", process="expand", index="0", \
		extname="", extversion="", ikparams="", omit="", \
		replace="", outfile="STDOUT", logfile="dev$null", glogpars="",
		verbose-) \
		| scan (phu)

	    if (gemextn.fail_count != 0 || gemextn.count != 1) {
		printlog ("ERROR - GEMWCSCOPY: Missing PHU in " // root \
		    // ".", l_logfile, verbose+)
		goto clean
	    }

	    imdelete (tmpphu, verify-, >& "dev$null")
	    tmpphu = mktemp ("tmpphu") // ".fits"
	    delete (tmpheader, verify-, >& "dev$null")
	    tmpheader = mktemp ("tmpheader")

	    imheader (phu, imlist="*.fits", longheader+, userfields+,
		> tmpheader)
	    mknoise (tmpphu, output="", title="", ncols=2, nlines=2, \
		header=tmpheader, background=0, gain=1, rdnoise=0, \
		poisson-, seed=1, cosrays="", ncosrays=0, energy=30000, \
		radius=0.5, ar=1, pa=0, comments-)

	    if (debug) imhead (tmpphu)

	    if (debug) print ("updating phu as simple file")
	    
	    wcscopy (tmpphu, tmprefphu, verbose=l_verbose)


	    if (debug) print ("replacing phu")

	    gemextn (tmpphu, check="", process="expand", index="0", \
		extname="", extversion="", ikparams="", omit="", \
		replace="", outfile="STDOUT", logfile="dev$null", glogpars="",
		verbose-) \
		| scan (phufrom)

	    delete (tmpheader, verify-, >& "dev$null")
	    tmpheader = mktemp ("tmpheader")

	    imheader (phufrom, imlist="*.fits", longheader+, userfields+,
		> tmpheader)

	    gemextn (root, check="", process="expand", index="0", \
		extname="", extversion="", ikparams="overwrite", omit="", \
		replace="", outfile="STDOUT", logfile="dev$null", glogpars="",
		verbose-) \
		| scan (phuto)

	    mkheader (phuto, header=tmpheader, append-, verbose=debug)

	}

	status = 0


clean:

	scanin1 = ""
	scanin2 = ""
	delete (tmplog, verify-, >& "dev$null")
	delete (tmpin, verify-, >& "dev$null")
	delete (tmpref, verify-, >& "dev$null")
	delete (tmprefphu, verify-, >& "dev$null")
	delete (tmpphu, verify-, >& "dev$null")
	delete (tmpheader, verify-, >& "dev$null")

end

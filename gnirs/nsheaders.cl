# Copyright(c) 2003-2011 Association of Universities for Research in Astronomy, Inc.

# This code brought to you by acooke, 6 Dec 2003.

# Group all header keywords (and a few other common parameters) in one
# place.

procedure nsheaders (instrument)

char	instrument	{"gnirs", prompt = "The instrument file to read"}
char	directory	{"gnirs$data",prompt = "The directory containing instrument-specific scripts"}

char	logfile		{"", prompt = "Logfile"}
bool	verbose		{no, prompt = "Verbose output?"}

char	sci_ext		{prompt = "Name of science extension"}				# OLDP-3
char	var_ext		{prompt = "Name of variance extension"}				# OLDP-3
char	dq_ext		{prompt = "Name of data quality extension"}			# OLDP-3

char    key_instrument  {prompt = "Header keyword for instrument name"} # OLDP-3
char	key_bias	{prompt = "Header keyword for detector bias (V)"}		# OLDP-3
char	key_ron		{prompt = "Header keyword for read noise (e-)"}			# OLDP-3
char	key_gain	{prompt = "Header keyword for gain (e-/ADU)"}			# OLDP-3
char	key_sat		{prompt = "Header keyword for saturation (ADU)"}		# OLDP-3
char	key_nonlinear	{prompt = "Header keyword for non-linear regime (ADU)"}		# OLDP-3
char    key_arrayid {prompt = "Header keyword for array ID"}            # OLDP-3
char	key_filter	{prompt = "Header keyword for filter"}				# OLDP-3
char	key_decker	{prompt = "Header keyword for decker"}				# OLDP-3
char	key_slit	{prompt = "Header keyword for slit"}				# OLDP-3
char	key_prism	{prompt = "Header keyword for prism"}				# OLDP-3
char	key_order	{prompt = "Header keyword for storing spectral order"}		# OLDP-3
char	key_ndavgs	{prompt = "Header keyword for digital averages"}		# OLDP-3
char	key_coadds	{prompt = "Header keyword for number of co-adds"}		# OLDP-3
char	key_lnrs	{prompt = "Header keyword for number of non-destructive reads"}	# OLDP-3

char	key_camera	{prompt = "Header keyword for camera"}				# OLDP-3
char	key_grating	{prompt = "Header keyword for grating"}				# OLDP-3
char	key_fpmask	{prompt = "Header keyword for focal plane mask (slit)"}		# OLDP-3
char	key_dispaxis	{prompt = "Header keyword for dispersion axis"}			# OLDP-3
char	key_wave	{prompt = "Header keyword for central grating wavelength"}	# OLDP-3
char	key_delta	{prompt = "Header keyword for linear dispersion"}		# OLDP-3
char	key_waveorder	{prompt = "Header keyword for grating wavelength order"}	# OLDP-3
char	key_cradius	{prompt = "Header keyword for wavelength cal search radius"}	# OLDP-3
char	key_exptime 	{prompt = "Header keyword for exposure time"}			# OLDP-3
char	key_section 	{prompt = "Header keyword for image section(s) to cut"}		# OLDP-3
char	key_cut_section	{prompt = "Header keyword for image section that was cut"}	# OLDP-3
char	key_pixscale	{prompt = "Header keyword for pixel scale"}			# ODP-3
char	key_wavevar	{prompt = "Header keyword for nsappwave variable value"}	# ODP-3

char	key_dark	{prompt = "Header keyword to identify darks"}			# ODP-3
char	val_dark	{prompt = "Substring to match against dark header value"}	# ODP-3

char	key_xoff	{prompt = "Header keyword for storing gemoffsetlist x offset (in arcsec)"}									# OLDP-3
char	key_yoff	{prompt = "Header keyword for storing gemoffsetlist y offset (in arcsec)"} 									# OLDP-3
char    key_date        {prompt = "Header keyword for date of observation"}
char	key_time	{prompt = "Header keyword for time of observation"}		# OLDP-3
char	key_airmass	{prompt = "Header keyword for observation airmass"}		# OLDP-3
char    key_mode        {prompt = "Header keyword for reduction mode"}                  # OLDP-3
char    key_obstype     {prompt = "Header keyword for observation type"}                # OLDP-3

int	status		{0, prompt = "Exit status (0=good)"}

begin
	char	l_instrument = ""
	char	l_directory = ""
	char	l_logfile = ""
	bool	l_verbose

	char	filename
	int	junk
	struct	sdate


	status = 1
	junk = fscan (	instrument, l_instrument)
	junk = fscan (	directory, l_directory)
	junk = fscan (	logfile, l_logfile)
	l_verbose =	verbose
	
	date | scan (sdate)
	printf ("---------------------------------------------------------\
	    ---------------------\n")
	printf ("NSHEADERS -- " // sdate //"\n\n")

	filename = l_directory // "/" // l_instrument // ".dat"	    
	if (access (filename)) {
	    cl (< filename)
	    status = 0
	}

	# The default logfile should have been set in the instrument
	# .dat file, and saved to gnirs.logfile.
	if ("" == l_logfile) {
	    junk = fscan (gnirs.logfile, l_logfile) 
	    if (l_logfile == "") {
		l_logfile = "gnirs.log"
		printlog ("WARNING - NSHEADERS: Both nsheaders.logfile and \
		    gnirs.logfile are empty.", l_logfile, verbose+) 
		printlog ("                   Using default file " \
		    // l_logfile // ".", l_logfile, verbose+) 
	    }
	}

	printlog ("------------------------------------------------------------------------------", 
	    l_logfile, verbose-)
	printlog ("NSHEADERS -- " // sdate, l_logfile, verbose-) 
	printlog (" ", l_logfile, verbose-) 
	
	if (0 == status) {
	    printlog ("Set header values from " // filename, \
		l_logfile, l_verbose)
	    printlog ("                   for " // l_instrument // ".",
		l_logfile, l_verbose)
	} else {
	    printlog ("ERROR - NSHEADERS: Cannot load file " // filename, \
		l_logfile, verbose+)
	    printlog ("                   for " // l_instrument // ".",
		l_logfile, verbose+)
	    goto clean
	}
	
clean:

	if (0 == status) {
	    printlog (" ", l_logfile, l_verbose) 
	    printlog ("NSHEADERS exit status: good.", l_logfile, l_verbose) 
	} else {
	    printlog (" ", l_logfile, l_verbose)
	    printlog ("NSHEADERS exit status: failed.", l_logfile, l_verbose)
	}

end

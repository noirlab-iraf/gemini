# Copyright(c) 2006-2011 Association of Universities for Research in Astronomy, Inc.

procedure msreduce (inspec)

# This is a general wrapper task for doing Michelle or T-ReCS spectroscopic
# reductions...it calls various routines in the midir and gnirs packages.
#
# This can do several things:
#
# 1) extract a single spectrum for an object, resulting in a wavelength
#    calibrated spectrum in terms of raw counts
# 2) extract the spectra of an object and a telluric standard, and form
#    an extracted telluric corrected spectrum assuming the standard 
#    follows a blackbody intensity distribution
# 3) extract the spectra of an object and a spectrophotometric standad 
#    (also a telluric reference), and form an extracted absolutely 
#    calibrated spectrum with telluric effects removed [in principle...]
#

# Version:   October 26, 2005   KV creates the task 
#            december 4, 2006   KV modifies the calls to msabsflux (small change)
#            May 2,      2008   NZ split nstransform into nsfitcoords and
#                                  nstransform (with no fitcoords)

char    inspec          {prompt="Name of source spectrum raw data file(s)"}
char    outtype         {"fnu",prompt="Type of output spectrum: fnu|flambda|lambda*flambda"}
char    rawpath         {"",prompt="Path for raw images"}
int     line            {160,prompt="Line to use to identify the aperture"}
bool    fl_std          {no,prompt="Extract both source and standard spectra, do telluric correction"}
char    std             {"",prompt="Name of standard spectrum raw data file(s)"}
char    stdname         {"",prompt="Name of standard star"}
bool    fl_bbody        {no,prompt="Use blackbody option rather than spectrophotometry"}
real    bbody           {10000.0, prompt="Temperature of calibrator for black-body fit"}
bool    fl_flat         {no,prompt="Apply flat/bias to the raw data file(s) before processing"}
char    flat            {"",prompt="Input flat-field image(s) (raw data file(s))"}
char    bias            {"",prompt="Input bias image(s) (raw data file(s))"}
bool    fl_extract      {no,prompt="Run nsextract interactively"}
bool    fl_telluric     {no,prompt="Run msabsflux interactively"}
bool    fl_wavelength   {no,prompt="Run nswavelength interactively"}
bool    fl_fitcoords    {no,prompt="Run nsfitcoords interactively"}
bool    fl_retrace      {no,prompt="Retrace spectrum, if appropriate"}
bool    fl_process      {yes,prompt="Do initial processing with mprepare"}
bool    fl_clear        {yes,prompt="Clear database sub-directory"}
bool    fl_negative     {no,prompt="Extract from a negative spectrum (NOD mode)"}
bool    fl_defringe     {yes,prompt="Defringe the spectra (an intrinsicly interactive step)"}
bool    fl_lowres       {yes,prompt="Low resolution spectrum [msdefringe]"}
bool    fl_skybiassub   {no,prompt="Subtract bias from sky frame"}
int     fmin            {18,prompt="Start of filtering region (pixels) [msdefringe]"}
int     fmax            {32,prompt="End of filtering region (pixels) [msdefringe]"}
bool    fl_dfinterp     {yes,prompt="Interpolate across masked region [msdefringe]"}
bool    fl_zerocut      {yes,prompt="Mask out negative points [msdefringe]"}
bool    fl_reextract    {no,prompt="Just re-extract the spectra"}
char    linelist        {"",prompt="Line list file name"}
char    logfile         {"",prompt="Log file name"}
bool    verbose         {yes,prompt="verbose logging?"}
int     status          {0,prompt="Exit error status: (0=good, >0=bad)"}
struct  *scanfile       {"",prompt="Internal use only"}

begin
    char    l_inspec = ""
    char    l_std = ""
    char    l_outtype = ""
    char    l_rawpath = ""
    char    l_stdname = ""
    char    l_flat = ""
    char    l_bias = ""
    char    l_linelist = ""
    char    l_logfile = ""
    int     l_line, l_fmin, l_fmax
    real    l_bbody
    bool    l_fl_std, l_fl_bbody, l_fl_flat, l_fl_extract, l_fl_telluric
    bool    l_fl_wavelength, l_fl_fitcoords, l_fl_retrace, l_fl_process
    bool    l_fl_clear, l_fl_negative, l_fl_defringe, l_fl_lowres
    bool    l_fl_dfinterp, l_fl_zerocut, l_fl_reextract
    bool    l_verbose, l_fl_plots, l_fl_skybiassub

    char    paramstr, lastchar, errmsg, filename, keyfound
    char    source, reference, flatfield, biasframe
    char    instrument, sinstrument, curpath, prefix, sprefix
    char    tmpinimg, tmpfile, tmpneg, tmpwork, tmpbias
    char    in[20], ref[20], flatim[20], biasim[20]
    int     junk, maxfiles, nsources, nstandards, nflats, nbiases
    int     i, tmpstatus, nskipped, ngood

    status = 0
    ngood = 0
    nskipped = 0
    
    cache ("gemextn", "nsheaders")
    
    junk = fscan (inspec, l_inspec)
    junk = fscan (std, l_std)
    junk = fscan (outtype, l_outtype)
    junk = fscan (rawpath, l_rawpath)
    l_line = line
    l_fl_std = fl_std
    junk = fscan (stdname, l_stdname)
    l_fl_bbody = fl_bbody
    l_bbody = bbody
    l_fl_flat = fl_flat
    junk = fscan (flat, l_flat)
    junk = fscan (bias, l_bias)
    l_fl_extract = fl_extract
    l_fl_telluric = fl_telluric
    l_fl_wavelength = fl_wavelength
    l_fl_fitcoords = fl_fitcoords
    l_fl_retrace = fl_retrace
    l_fl_process = fl_process
    l_fl_clear = fl_clear
    l_fl_negative = fl_negative
    l_fl_defringe = fl_defringe
    l_fl_lowres = fl_lowres
    l_fl_skybiassub = fl_skybiassub
    l_fmin = fmin
    l_fmax = fmax
    l_fl_dfinterp = fl_dfinterp
    l_fl_zerocut = fl_zerocut
    l_fl_reextract = fl_reextract
    junk = fscan (linelist, l_linelist)
    junk = fscan (logfile, l_logfile)
    l_verbose = verbose

    # Initialize
    maxfiles = 20

    # Create temp files
    tmpinimg = mktemp ("tmpinimg")
    tmpfile = mktemp ("tmpfile")

    # Create the list of parameter/value pairs.  One pair per line.
    # All lines combined into one string.  Line delimiter is '\n'.
    paramstr =  "inspec         = "//inspec.p_value//"\n"
    paramstr += "std            = "//std.p_value//"\n"
    paramstr += "outtype        = "//outtype.p_value//"\n"
    paramstr += "rawpath        = "//rawpath.p_value//"\n"
    paramstr += "fl_std         = "//fl_std.p_value//"\n"
    paramstr += "stdname        = "//stdname.p_value//"\n"
    paramstr += "fl_bbody       = "//fl_bbody.p_value//"\n"
    paramstr += "bbody          = "//bbody.p_value//"\n"
    paramstr += "fl_flat        = "//fl_flat.p_value//"\n"
    paramstr += "flat           = "//flat.p_value//"\n"
    paramstr += "bias           = "//bias.p_value//"\n"
    paramstr += "fl_extract     = "//fl_extract.p_value//"\n"
    paramstr += "fl_telluric    = "//fl_telluric.p_value//"\n"
    paramstr += "fl_wavelength  = "//fl_wavelength.p_value//"\n"
    paramstr += "fl_fitcoords   = "//fl_fitcoords.p_value//"\n"
    paramstr += "fl_retrace     = "//fl_retrace.p_value//"\n"
    paramstr += "fl_process     = "//fl_process.p_value//"\n"
    paramstr += "fl_clear       = "//fl_clear.p_value//"\n"
    paramstr += "fl_negative    = "//fl_negative.p_value//"\n"
    paramstr += "fl_defringe    = "//fl_defringe.p_value//"\n"
    paramstr += "fl_lowres      = "//fl_lowres.p_value//"\n"
    paramstr += "fmin           = "//fmin.p_value//"\n"
    paramstr += "fmax           = "//fmax.p_value//"\n"
    paramstr += "fl_dfinterp    = "//fl_dfinterp.p_value//"\n"
    paramstr += "fl_zerocut     = "//fl_zerocut.p_value//"\n"
    paramstr += "fl_lowres      = "//fl_lowres.p_value//"\n"
    paramstr += "fl_reextract   = "//fl_reextract.p_value//"\n"
    paramstr += "logfile        = "//logfile.p_value//"\n"
    paramstr += "verbose        = "//verbose.p_value

    # Assign a logfile name if not given.  Open logfile and start log.
    # Write parameter/value pairs ("paramstr") to log.
    gloginit (l_logfile, "msreduce", "midir", paramstr, fl_append+,
        verbose=l_verbose)
    if (gloginit.status != 0) {
        status = gloginit.status
        goto exitnow
    }
    l_logfile = gloginit.logfile.p_value

    # Check that nsheaders has been run (not full proof, but it should
    # catch some mistakes)
#    if (nsheaders.sci_ext == "") {
#        errmsg = "Package configuration has not been done.  Run NSHEADERS."
#        status = 121
#        glogprint (l_logfile, "msreduce", "status", type="error", 
#            errno=status, str=errmsg, verbose+)
#        goto clean
#    }
# This does not work with the current implementation of nsheaders
#
#    if ((nsheaders.instrument != "michelle") && \
#        (nsheaders.instrument != "trecs")) {
#        
#        errmsg = "Package configuration has not been done.  Run NSHEADERS."
#        status = 121
#        glogprint (l_logfile, "msreduce", "status", type="error", 
#            errno=status, str=errmsg, verbose+)
#        goto clean
#    }

    # Check for the linelist file, if it is defined
    if (l_linelist == "")
        l_linelist="gnirs$data/sky.dat"

    if (no == access(l_linelist)) {
        glogprint (l_logfile, "msreduce", "status", type="error", errno=121,
            str="ERROR - MSREDUCE: line list file "//l_linelist//\
            " was not found.", verbose+)
        status = 1
        goto clean
    }

    glogprint (l_logfile, "msreduce", "status", type="string", 
        str="Line list file name = "//l_linelist, verbose+)

    # Add the trailing slash to rawpath, if missing
    lastchar = substr (l_rawpath, strlen(l_rawpath), strlen(l_rawpath))
    if ((l_rawpath != "") && (lastchar != "$") && (lastchar != "/"))
        l_rawpath = l_rawpath//"/"
    

    # Count the number of source images
    # First, generate the file list if needed

    gemextn (l_inspec, check="exists,mef", process="none", index="",
        extname="", extversion="", ikparam="", omit="kernel,exten",
        replace="^%%"//l_rawpath//"%", outfile=tmpinimg, logfile=l_logfile,
        verbose=l_verbose)
    nsources = gemextn.count
    
    if ((gemextn.fail_count > 0) || (nsources == 0) || \
        (nsources > maxfiles)) {
        
        if (gemextn.fail_count > 0) {
            errmsg = gemextn.fail_count//" source spectra were not found."
            status = 101
        } else if (nsources == 0) {
            errmsg = "No input source spectra defined."
            status = 121
        } else if (nsources > maxfiles) {
            errmsg = "Maximum number of input source spectra \
                ("//str(maxfiles)//") has been exceeded."
            status = 121
        }
        
        glogprint (l_logfile, "msreduce", "status", type="error", errno=status,
            str=errmsg, verbose+)
        goto clean
    } else {
        scanfile = tmpinimg
        i = 0
        while (fscan (scanfile, filename) != EOF) {
            i += 1
            in[i] = filename
        }
        scanfile = ""
        if (i != nsources) {
            status = 99
            errmsg = "Error while counting the input images."
            glogprint (l_logfile, "msreduce", "status", type="error",
                errno=status, str=errmsg, verbose+)
        }
    }
    delete (tmpinimg, ver-, >& "dev$null")

    # Check input standard/telluric spectra
    if (l_fl_std) {
        gemextn (l_std, check="exists,mef", process="none", index="",
            extname="", extversion="", ikparam="", omit="kernel,exten",
            replace="^%%"//l_rawpath//"%", outfile=tmpfile, logfile=l_logfile,
            verbose=l_verbose)
        nstandards = gemextn.count
        
        if ((gemextn.fail_count > 0) || (nstandards == 0) || \
            (nstandards > maxfiles) || (nstandards != nsources)) {
            
            if (gemextn.fail_count > 0) {
                errmsg = gemextn.fail_count//" standard/telluric spectra were \
                    not found."
                status = 101
            } else if (nstandards == 0) {
                errmsg = "No input standard/telluric spectra defined."
                status = 121
            } else if (nstandards > maxfiles) {
                errmsg = "Maximum number of input standard/telluric spectra \
                    ("//str(maxfiles)//") has been exceeded."
                status = 121
            } else if (nstandards != nsources) {
                errmsg = "Different number of input sources ("//nsources//") \
                    and standard/telluric ("//nstandards//") spectra."
                status = 121
            }
            
            glogprint (l_logfile, "msreduce", "status", type="error",
                errno=status, str=errmsg, verbose+)
            goto clean
        } else {
            scanfile = tmpfile
            i = 0
            while (fscan (scanfile, filename) != EOF) {
                i += 1
                ref[i] = filename
            }
            scanfile = ""
            if (i != nstandards) {
                status = 99
                errmsg = "Error while counting the input telluric/standard \
                    spectra."
                glogprint (l_logfile, "msreduce", "status", type="error",
                    errno=status, str=errmsg, verbose+)
            }
        }
        delete (tmpfile, ver-, >& "dev$null")        
    }

    # Check input raw flats and biases
    if (l_fl_flat) {
    
        # Start with the flats
        gemextn (l_flat, check="exists,mef", process="none", index="",
            extname="", extversion="", ikparam="", omit="kernel,exten",
            replace="^%%"//l_rawpath//"%", outfile=tmpfile, logfile=l_logfile,
            verbose=l_verbose)
        nflats = gemextn.count
        
        if ((gemextn.fail_count > 0) || (nflats == 0) || \
            (nflats > maxfiles) || (nflats != nsources)) {
            
            if (gemextn.fail_count > 0) {
                errmsg = gemextn.fail_count//" flats were not found."
                status = 101
            } else if (nflats == 0) {
                errmsg = "No input flats defined."
                status = 121
            } else if (nflats > maxfiles) {
                errmsg = "Maximum number of flats ("//str(maxfiles)//") has \
                    been exceeded."
                status = 121
            } else if (nflats != nsources) {
                errmsg = "Different number of input sources ("//nsources//") \
                    and flats ("//nflats//")."
                status = 121
            }
            
            glogprint (l_logfile, "msreduce", "status", type="error",
                errno=status, str=errmsg, verbose+)
            goto clean
        } else {
            scanfile = tmpfile
            i = 0
            while (fscan (scanfile, filename) != EOF) {
                i += 1
                flatim[i] = filename
            }
            scanfile = ""
            if (i != nflats) {
                status = 99
                errmsg = "Error while counting the input flats."
                glogprint (l_logfile, "msreduce", "status", type="error",
                    errno=status, str=errmsg, verbose+)
            }
        }
        delete (tmpfile, ver-, >& "dev$null")
        
        # Now with the biases
        gemextn (l_bias, check="exists,mef", process="none", index="",
            extname="", extversion="", ikparam="", omit="kernel,exten",
            replace="^%%"//l_rawpath//"%", outfile=tmpfile, logfile=l_logfile,
            verbose=l_verbose)
        nbiases = gemextn.count
        
        if ((gemextn.fail_count > 0) || (nbiases == 0) || \
            (nbiases > maxfiles) || (nbiases != nsources)) {
            
            if (gemextn.fail_count > 0) {
                errmsg = gemextn.fail_count//" biases were not found."
                status = 101
            } else if (nbiases == 0) {
                errmsg = "No input biases defined."
                status = 121
            } else if (nbiases > maxfiles) {
                errmsg = "Maximum number of biases ("//str(maxfiles)//") has \
                    been exceeded."
                status = 121
            } else if (nbiases != nsources) {
                errmsg = "Different number of input sources ("//nsources//") \
                    and biases ("//nbiases//")."
                status = 121
            }
            
            glogprint (l_logfile, "msreduce", "status", type="error",
                errno=status, str=errmsg, verbose+)
            goto clean
        } else {
            scanfile = tmpfile
            i = 0
            while (fscan (scanfile, filename) != EOF) {
                i += 1
                biasim[i] = filename
            }
            scanfile = ""
            if (i != nbiases) {
                status = 99
                errmsg = "Error while counting the input biases."
                glogprint (l_logfile, "msreduce", "status", type="error",
                    errno=status, str=errmsg, verbose+)
            }
        }
        delete (tmpfile, ver-, >& "dev$null")
    } else if (l_fl_skybiassub) {
        delete (tmpfile, ver-, >& "dev$null")

        gemextn (l_bias, check="exists,mef", process="none", index="",
            extname="", extversion="", ikparam="", omit="kernel,exten",
            replace="^%%"//l_rawpath//"%", outfile=tmpfile, logfile=l_logfile,
            verbose=l_verbose)
        nbiases = gemextn.count
         
        if ((gemextn.fail_count > 0) || (nbiases == 0) || \
            (nbiases > maxfiles) || (nbiases != nsources)) {
            
            if (gemextn.fail_count > 0) {
                errmsg = gemextn.fail_count//" biases were not found."
                status = 101
            } else if (nbiases == 0) {
                errmsg = "No input biases defined."
                status = 121
            } else if (nbiases > 1) {
                errmsg = "Only one bias image allowed."
                status = 121
            }
            
            glogprint (l_logfile, "msreduce", "status", type="error",
                errno=status, str=errmsg, verbose+)
            goto clean
        } else {
            scanfile = tmpfile
            i = 0
            while (fscan (scanfile, filename) != EOF) {
                i += 1
                biasim[i] = filename
            }
            scanfile = ""
            if (i != nbiases) {
                status = 99
                errmsg = "Error while counting the input biases."
                glogprint (l_logfile, "msreduce", "status", type="error",
                    errno=status, str=errmsg, verbose+)
            }
        }
        delete (tmpfile, ver-, >& "dev$null")
    }
    ######
    #  Start the work
    ######

    glogprint (l_logfile, "msreduce", "status", type="string", 
        str="MSREDUCE: Processing "//str(nsources)//" image(s).", verbose+)

    ngood = 0
    nskipped = 0
    for (i=1; i<=nsources; i=i+1) {
    
        # Check instrument against nsheader configuration
        instrument = ""
        hselect (in[i]//"[0]", "INSTRUME", yes) | scan (instrument)
        instrument = strlwr(instrument)
        if ((instrument != "michelle") && (instrument != "trecs")) {
            status = 123
            errmsg = "Source file "//in[i]//" not Michelle or T-ReCs data.  \
                Skipping."
            glogprint (l_logfile, "msreduce", "status", type="warning",
                errno=status, str=errmsg, verbose=l_verbose)
            nskipped += 1
            goto nextimage
        }
        else {
            if (i == 1) {
                if (instrument == "michelle") {
#                    printf("Instrument = michelle; running nsheaders\n")
                    nsheaders("michelle")
                }
                if (instrument == "trecs") {
                    nsheaders("trecs")
#                    printf("Instrument = michelle; running nsheaders\n")
                }
            }
        }
# This does not work with the current implementation of nsheaders
#
#        if (instrument != nsheaders.instrument) {
#            status = 123
#            errmsg = "Source "//in[i]//" not "//nsheaders.instrument//" \
#                data.  Skipping."
#            glogprint (l_logfile, "msreduce", "status", type="warning",
#                errno=status, str=errmsg, verbose=l_verbose)
#            nskipped += 1
#            goto nextimage
#        }

        # Ensure input data are spectroscopic observations
        if (instrument == "michelle") {
            keyfound = ""
            hselect (in[i]//"[0]", "CAMERA", yes) | scan (keyfound)
            if (keyfound != "spectroscopy")
                status = 121
        } else if (instrument == "trecs") {
            keyfound = ""
            hselect (in[i]//"[0]", "GRATING", yes) | scan (keyfound)
            if (keyfound == "Mirror")
                status = 121
        }
        if (status != 0) {
            errmsg = "Source file "//in[i]//" not a spectrum.  Skipping."
            glogprint (l_logfile, "msreduce", "status", type="warning",
                errno=status, str=errmsg, verbose=l_verbose)
            nskipped += 1
            goto nextimage
        }

        if ((instrument == "trecs") && (l_fl_flat)) {
            l_fl_flat = no
            glogprint (l_logfile, "msreduce", "status", type="warning",
                errno=121, str="Flat fielding not available for T-ReCS data.  \
                Resetting fl_flat to 'no'.", verbose=l_verbose)
        }

        tmpstatus = 0
        if (l_fl_flat) {
        
            # Check instrument.  Must match input.
            
            sinstrument = ""
            hselect (flatim[i]//"[0]", "INSTRUME", yes) | scan (sinstrument)
            sinstrument = strlwr(sinstrument)
            if (sinstrument != instrument) {
                tmpstatus = 123
                errmsg = "Flat file "//flatim[i]//" not "//instrument//" \
                    data.  Skipping."
                glogprint (l_logfile, "msreduce", "status", type="warning",
                    errno=tmpstatus, str=errmsg, verbose=l_verbose)
            }
#            if (instrument != nsheaders.instrument) {
#                tmpstatus = 123
#                errmsg = "Flat "//flatim[i]//" not "//nsheaders.instrument//"\
#                    data."
#                glogprint (l_logfile, "msreduce", "status", type="warning",
#                    errno=tmpstatus, str=errmsg, verbose=l_verbose)
#            }
            sinstrument = ""
            hselect (biasim[i]//"[0]", "INSTRUME", yes) | scan (sinstrument)
            if (sinstrument != instrument) {
                tmpstatus = 123
                errmsg = "Bias file "//biasim[i]//" not "//instrument//" \
                    data.  Skipping."
                glogprint (l_logfile, "msreduce", "status", type="warning",
                    errno=tmpstatus, str=errmsg, verbose=l_verbose)
            }
#            if (instrument != nsheaders.instrument) {
#                tmpstatus = 123
#                errmsg = "Bias"//biasim[i]//" not "//nsheaders.instrument//" \
#                    data."
#                glogprint (l_logfile, "msreduce", "status", type="warning",
#                    errno=tmpstatus, str=errmsg, verbose=l_verbose)
#            }

            # Ensure spectroscopic observations
            errmsg = ""
            if (sinstrument == "michelle") {
                keyfound = ""
                hselect (flatim[i]//"[0]", "CAMERA", yes) | scan (keyfound)
                if (keyfound != "spectroscopy")
                    tmpstatus = 121
                keyfound = ""
                hselect (biasim[i]//"[0]", "CAMERA", yes) | scan (keyfound)
                if (keyfound != "spectroscopy")
                    tmpstatus = 121
            } else if (sinstrument == "trecs") {
                keyfound = ""
                hselect (flatim[i]//"[0]", "GRATING", yes) | scan (keyfound)
                if (keyfound == "Mirror")
                    tpmstatus = 121
                keyfound = ""
                hselect (biasim[i]//"[0]", "GRATING", yes) | scan (keyfound)
                if (keyfound == "Mirror")
                    tpmstatus = 121
            }
            if (tmpstatus != 0) {
                errmsg = "Flat "//flatim[i]//" and/or bias "//biasim[i]//" \
                    not a spectrum."
                glogprint (l_logfile, "msreduce", "status", type="warning",
                    errno=status, str=errmsg, verbose=l_verbose)
            }
                      
        }
        if (l_fl_std) {
            sinstrument = ""
            hselect (ref[i]//"[0]", "INSTRUME", yes) | scan (sinstrument)
            sinstrument = strlwr(sinstrument)
            if (sinstrument != instrument) {
                tmpstatus = 123
                errmsg = "Standard/telluric "//ref[i]//" not "//instrument//\
                    " data."
            }
#            if (instrument != nsheaders.instrument) {
#                tmpstatus = 123
#                errmsg = "Standard/telluric "//ref[i]//" not "//\
#                    nsheaders.instrument//" data."
#                glogprint (l_logfile, "msreduce", "status", type="warning",
#                    errno=tmpstatus, str=errmsg, verbose=l_verbose)
#            }      
      
            # Ensure spectroscopic observations
            errmsg = ""
            if (sinstrument == "michelle") {
                keyfound = ""
                hselect (ref[i]//"[0]", "CAMERA", yes) | scan (keyfound)
                if (keyfound != "spectroscopy") {
                    tmpstatus = 121
                    errmsg = "Standard/telluric "//ref[i]//" not a spectrum."
                }
            } else if (sinstrument == "trecs") {
                keyfound = ""
                hselect (ref[i]//"[0]", "GRATING", yes) | scan (keyfound)
                if (keyfound == "Mirror") {
                    tpmstatus = 121
                    errmsg = "Standard/telluric "//ref[i]//" not a spectrum."
                }
            }
            if (tmpstatus != 0) {
                glogprint (l_logfile, "msreduce", "status", type="warning",
                    errno=status, str=errmsg, verbose=l_verbose)
            }
        }
        
        if (tmpstatus != 0) {
            status = 121
            errmsg = "Inappropriate calibration.  Skipping source "//in[i]
            glogprint (l_logfile, "msreduce", "status", type="warning",
                errno=status, str=errmsg, verbose=l_verbose)
            nskipped += 1
            goto nextimage
        }
        
        # Parse file names, remove path and extension
        # in[i], ref[i], flatim[i], biasim[i] include l_rawpath
        fparse (in[i], verbose-)
        source = fparse.root
        if (l_fl_std) {
            fparse (ref[i], verbose-)
            reference = fparse.root
        }
        if (l_fl_flat) {
            fparse (flatim[i], verbose-)
            flatfield = fparse.root
            fparse (biasim[i], verbose-)
            biasframe = fpars.root
        }  
        curpath = l_rawpath      # update this to "" after first processing
        
        # If extracting for the first time (fl_reextract == no)        
        if (no == l_fl_reextract) {
        
            # Reset the content of the database for this source
            if (l_fl_clear && access("database")) {
                delete ("database/*"//source//"*", verify-, >& "dev$null")
                if (l_fl_std)
                  delete ("database/*"//reference//"*", verify-, >& "dev$null")
            }
            
            # Do flat correction
            if (l_fl_flat) {
                imdelete ("f"//source, ver-, >& "dev$null")
                imdelete ("normflat_"//source, ver-, >& "dev$null")
                msflatcor (source, flatfield, biasframe, outimage="f"//source,
                    outpref="", rawpath=curpath, fl_bias+, fl_writeflat+,
                    normflat="normflat_"//source, logfile=l_logfile,
                    verbose=l_verbose)
                if (msflatcor.status != 0) {
                    nskipped += 1
                    goto nextimage
                } else
                    source = "f"//source
                
                if (l_fl_std) {
                    imdelete ("f"//reference, ver-, >& "dev$null")
                    imdelete ("normflat_"//reference, ver-, >& "dev$null")
                    msflatcor (reference, flatfield, biasframe, 
                        outimage="f"//reference, outpref="", rawpath=curpath,
                        fl_bias+, fl_writeflat+,
                        normflat="normflat_"//reference, logfile=l_logfile,
                        verbose=l_verbose)
                    if (msflatcor.status != 0) {
                        nskipped += 1
                        goto nextimage
                    } else
                        reference = "f"//reference
                }
                curpath = ""
            }
            
            if (l_fl_process) {
                # Create combined difference (spectrum) and reference
                # (sky lines) images for the target.
                
                imdelete ("r"//source, ver-, >& "dev$null")
                imdelete ("m"//source, ver-, >& "dev$null")
                imdelete ("t"//source, ver-, >& "dev$null")
                mireduce (source, outimages="r"//source, outpref="",
                    rawpath=curpath, fl_background-, fl_view-, fl_mask-, 
                    bpm="", fl_flat-, flatfieldfile="", stackoption="stack",
                    frametype="dif", combine="average", fl_display=no,
                    fl_check=yes, fl_rescue=no, region="[*,*]", fl_variance=no,
                    logfile=l_logfile, verbose=l_verbose)
                
                imdelete ("a"//source, ver-, >& "dev$null")
                imdelete ("m"//source, ver-, >& "dev$null")
                imdelete ("t"//source, ver-, >& "dev$null")
                mireduce (source, outimages="a"//source, outpref="",
                    rawpath=curpath, fl_background=no, fl_view=no, fl_mask=no,
                    bpm="", fl_flat=no, flatfieldfile="", stackoption="stack",
                    frametype="ref", combine="average", fl_display=no,
                    fl_check=yes, fl_rescue=no, region="[*,*]",
                    fl_variance=no, logfile=l_logfile, verbose=l_verbose)
                
                if (l_fl_std) {
                    # Create combined difference (spectrum) and reference
                    # (sky lines) images for the reference star.
                    
                    imdelete ("r"//reference, verify-, >& "dev$null")
                    imdelete ("m"//reference, verify-, >& "dev$null")
                    imdelete ("t"//reference, verify-, >& "dev$null")
                    mireduce (reference, outimages="r"//reference, outpref="",
                        rawpath=curpath, fl_background=no, fl_view=no,
                        fl_mask=no, bpm="", fl_flat=no, flatfieldfile="",
                        stackoption="stack", frametype="dif", 
                        combine="average", fl_display=no, fl_check=yes,
                        fl_rescue=no, region="[*,*]", fl_variance=no,
                        logfile=l_logfile, verbose=l_verbose)
                    
                    imdelete ("a"//reference, verify-, >& "dev$null")
                    imdelete ("m"//reference, verify-, >& "dev$null")
                    imdelete ("t"//reference, verify-, >& "dev$null")
                    mireduce (reference, outimages="a"//reference, outpref="",
                        rawpath=curpath, fl_background=no, fl_view=no,
                        fl_mask=no, bpm="", fl_flat=no, flatfieldfile="",
                        stackoption="stack", frametype="ref", 
                        combine="average", fl_display=no, fl_check=yes,
                        fl_rescue=no, region="[*,*]", fl_variance=no,
                        logfile=l_logfile, verbose=l_verbose)
                }
                    
            } else {
                if (no == access ("a"//source//".fits")) {
                    mistack ("t"//source, outimages="a"//source, outpref="",
                        rawpath="", frametype="ref", combine="average",
                        fl_variance=no, logfile=l_logfile, verbose=l_verbose)
                }
                if (l_fl_std) {
                  if (no == access ("a"//reference//".fits")) {
                    mistack ("t"//reference, outimages="a"//reference,
                        outpref="", rawpath="", frametype="ref",
                        combine="average", fl_variance=no, logfile=l_logfile,
                        verbose=l_verbose)
                  }
                }
            }
            
            # Identify sky lines for the target and standard, every tenth
            # line in the longslit spectrum of the sky reference images.
            if (l_fl_skybiassub) {
                tmpbias = mktemp ("tmpbias")
                imgets (biasim//".fits[0]", "MPREPARE", >& "dev$null")
                if (imgets.value == "0") {
                    mprepare (biasim, rawpath="", outimage=tmpbias, \
                        verbose=no, fl_rescue=no,logfile=l_logfile)
                } else {
                    copy (biasim//".fits", tmpbias//".fits", verbose-)
                }

                imarith ("a"//source//"[1]", "-", tmpbias//"[SCI,1]", \
                    "a"//source//"[SCI,1,overwrite]")
                if (l_fl_std) {
                    imarith ("a"//reference//"[1]", "-", tmpbias//"[SCI,1]", \
                    "a"//reference//"[SCI,1,overwrite]") 
                }               

                imdelete (tmpbias, ver-, >& "dev$null")

            }     
            curpath = ""        # everything should now be in the current dir.
            
            imdelete ("wa"//source, verify-, >& "dev$null")
            nswavelength ("a"//source, outspectra="wa"//source, outprefix="",
                crval=INDEF, cdelt=INDEF, crpix=INDEF, dispaxis=1, database="",
                coordlist=l_linelist, fl_inter=l_fl_wavelength,
                nsappwavedb="gnirs$data/nsappwave.fits", fl_median=no, 
                sdist="", sdorder=4, 
                xorder=2, yorder=2, aptable="gnirs$data/apertures.fits",
                section="default", nsum=10, ftype="emission", fwidth=4.,
                cradius=5., threshold=0., minsep=2, match=-6,
                function="chebyshev", order=4, sample="*", niterate=10,
                low_reject=3., high_reject=3., grow=0., refit=yes, step=10,
                trace=no, nlost=3, fl_overwrite+, aiddebug="", fmatch=0.2,
                nfound=6, sigma=0.05, rms=0.1, logfile=l_logfile,
                verbose=l_verbose, debug=no)
            
            if (l_fl_std) {
                imdelete ("wa"//reference, verify-, >& "dev$null")
                nswavelength ("a"//reference, outspectra="wa"//reference,
                    outprefix="", crval=INDEF, cdelt=INDEF, crpix=INDEF,
                    dispaxis=1, database="", coordlist=l_linelist,
                    fl_inter=l_fl_wavelength, 
                    nsappwavedb="gnirs$data/nsappwave.fits", sdist="",
                    sdorder=4, xorder=2, yorder=2,
                    aptable="gnirs$data/apertures.fits", section="default",
                    nsum=10, ftype="emission", fwidth=4., cradius=5.,
                    threshold=0., minsep=2, match=-6, function="chebyshev",
                    order=4, sample="*", niterate=10, low_reject=3.,
                    high_reject=3., grow=0., refit=yes, step=10, trace=no,
                    nlost=3, fl_overwrite+, aiddebug="", fmatch=0.2, nfound=6,
                    sigma=0.05, rms=0.1, logfile=l_logfile, verbose=l_verbose,
                    debug=no)
            }
            
            # Find the wavelength transformation function across the
            # longslit sky spectrum.
            imdelete ("fr"//source, verify-, >& "dev$null")
            imdelete ("tfr"//source, verify-, >& "dev$null")
            nsfitcoords ("r"//source, outspectra="", outprefix="f",
              lamptransf="wa"//source, sdisttransf="", dispaxis=1,
              database="", fl_inter=l_fl_fitcoords, fl_align=no,
              function="chebyshev", lxorder=3, lyorder=3, sxorder=3,
              syorder=3, pixscale=1.,
              logfile=l_logfile, verbose=l_verbose, debug=no, force=no)
            if (nsfitcoords.status != 0) goto clean    
            nstransform (inimages='fr'//source, outspectra="tfr"//source,
              dispaxis=1, database="", fl_stripe=no,
              interptype="poly3", xlog=no, ylog=no, pixscale=1.,
              logfile=l_logfile, verbose=l_verbose, debug=no)              
            if (nstransform.status != 0) goto clean

            if (l_fl_std) {
                imdelete ("fr"//reference, verify-, >& "dev$null")
                imdelete ("tfr"//reference, verify-, >& "dev$null")
                nsfitcoords ("r"//reference, outspectra="",          
                  outpref="f", lamptransf="wa"//reference, sdisttransf="",
                  dispaxis=1, database="", fl_inter=l_fl_fitcoords,
                  fl_align=no, function="chebyshev", lxorder=3,
                  lyorder=3, sxorder=3, syorder=3, pixscale=1.,
                  logfile=l_logfile, verbose=l_verbose, debug=no, force=no)
                if (nsfitcoords.status != 0) goto clean
                nstransform (inimages='fr'//reference, 
                  outspectra="tfr"//reference,
                  outpref="", dispaxis=1, database="", fl_stripe=no,
                  interptype="poly3", xlog=no, ylog=no, pixscale=1.,
                  logfile=l_logfile, verbose=l_verbose, debug=no)
                if (nstransform.status != 0) goto clean

            }
        } # end if (no == l_fl_reextract)
        
        
        # Extract a wavelength calibrated spectrum from an aperture near line 
        # l_line in the target and standard observations, using the wavelength 
        # calibration from the reference images.
        
        imdelete ("xtfr"//source, verify-, >& "dev$null")
        nsextract ("tfr"//source, outspectra="xtfr"//source, outprefix="",
            dispaxis=1, database="", line=l_line, nsum=10, ylevel=0.1,
            upper=10, lower=-10, background="average", fl_vardq+, fl_addvar+, 
            fl_skylines+, fl_inter=l_fl_extract, fl_apall+,
            fl_trace=l_fl_retrace, aptable="gnirs$data/apertures.fits",
            fl_usetabap-, fl_flipped+, fl_project-, fl_findneg-, bgsample="*",
            trace="", tr_nsum=10, tr_step=10, tr_nlost=3,
            tr_function="legendre", tr_order=5, tr_sample="*", tr_naver=1,
            tr_niter=1, tr_lowrej=3., tr_highrej=3., tr_grow=0.,
            weights="variance", logfile=l_logfile, verbose=l_verbose)
        if (l_fl_std) {
            imdelete ("xtfr"//reference, verify-, >& "dev$null")
            nsextract ("tfr"//reference, outspectra="xtfr"//reference,
                outprefix="", dispaxis=1, database="", line=l_line, nsum=10,
                ylevel=0.1, upper=10, lower=-10, background="average",
                fl_vardq+, fl_addvar+, fl_skylines+, fl_inter=l_fl_extract,
                fl_apall+, fl_trace=l_fl_retrace,
                aptable="gnirs$data/apertures.fits", fl_usetabap-, fl_flipped+,
                fl_project-, fl_findneg-, bgsample="*", trace="", tr_nsum=10,
                tr_step=10, tr_nlost=3, tr_function="legendre", tr_order=5,
                tr_sample="*", tr_naver=1, tr_niter=1, tr_lowrej=3.,
                tr_highrej=3., tr_grow=0., weights="variance",
                logfile=l_logfile, verbose=l_verbose)
        }

        # The following is because "fl_findneg" does not appear to work for 
        # Michelle NOD mode spectra.  Thus I have to negate, extract, and sum 
        # by hand as it were.
        
        if (l_fl_negative) {
            imdelete ("ntfr"//source, verify-, >& "dev$null")
            imdelete ("ptfr"//source, verify-, >& "dev$null")
            if (l_fl_std) {
                imdelete ("ntfr"//reference, verify-, >& "dev$null")
                imdelete ("ptfr"//reference, verify-, >& "dev$null")
            }
            
            tmpneg = mktemp ("tmpneg")
            tmpwork = mktemp ("tmpwork")
            imarith ("0", "-", "tfr"//source//"[1]", tmpwork)
            wmef (tmpwork, tmpneg, extname="SCI", phu="tfr"//source//"[0]",
                verbose-, >& "dev$null")
            imdelete (tmpwork, verify-, >& "dev$null")
            
            # Effectively extract the negative spectrum
            imdelete ("ntfr"//source, verify-, >& "dev$null")
            nsextract (tmpneg, outspectra="ntfr"//source, outprefix="",
                dispaxis=1, database="", line=l_line, nsum=10, ylevel=0.1,
                upper=10, lower=-10, background="average", fl_vardq+,
                fl_addvar+, fl_skylines+, fl_inter=l_fl_extract, fl_apall+, 
                fl_trace=l_fl_retrace, aptable="gnirs$data/apertures.fits",
                fl_usetabap-, fl_flipped+, fl_project-, fl_findneg-,
                bgsample="*", trace="", tr_nsum=10, tr_step=10, tr_nlost=3, 
                tr_function="legendre", tr_order=5, tr_sample="*", tr_naver=1, 
                tr_niter=1, tr_lowrej=3., tr_highrej=3., tr_grow=0., 
                weights="variance", logfile=l_logfile, verbose=l_verbose)
            imdelete (tmpneg, verify-, >& "dev$null")
            
            if (l_fl_std) {
                tmpneg = mktemp ("tmpneg")
                tmpwork = mktemp ("tmpwork")
                imarith ("0", "-", "tfr"//reference//"[1]", tmpwork)
                wmef (tmpwork, tmpneg, extname="SCI",
                    phu="tfr"//reference//"[0]", verbose-, >& "dev$null")
                imdelete (tmpwork, ver-, >& "dev$null")
                
                # Effectively extract the negative spectrum
                imdelete ("ntfr"//reference, verify-, >& "dev$null")
                nsextract (tmpneg, outspectra="ntfr"//reference, outprefix="",
                    dispaxis=1, database="", line=l_line, nsum=10, ylevel=0.1,
                    upper=10, lower=-10, background="average", fl_vardq+, 
                    fl_addvar+, fl_skylines+, fl_inter=l_fl_extract, fl_apall+,
                    fl_trace=l_fl_retrace, aptable="gnirs$data/apertures.fits",
                    fl_usetabap-, fl_flipped+, fl_project-, fl_findneg-, 
                    bgsample="*", trace="", tr_nsum=10, tr_step=10, tr_nlost=3,
                    tr_function="legendre", tr_order=5, tr_sample="*",
                    tr_naver=1, tr_niter=1, tr_lowrej=3., tr_highrej=3.,
                    tr_grow=0., weights="variance", logfile=l_logfile, 
                    verbose=l_verbose)
                imdelete (tmpneg, verify-, >& "dev$null")
            }
            
            # Sum the extracted positive and the extracted negative
            tmpwork = mktemp ("tmpwork")
            imrename ("xtfr"//source, "ptfr"//source, ver-, >& "dev$null")
            imsum ("ptfr"//source//"[1]","ntfr"//source//"[1]", tmpwork, 
                title="", hparams="", pixtype="", calctype="", option="sum",
                low_reject=0., high_reject=0., verbose-, >& "dev$null")
            wmef (tmpwork, "xtfr"//source, extname="SCI",
                phu="ptfr"//source//"[0]", verbose-, >& "dev$null")
            imdelete (tmpwork, ver-, >& "dev$null")
            imdelete ("ptfr"//source, ver-, >& "dev$null")
            imdelete ("ntfr"//source, ver-, >& "dev$null")
            
            if (l_fl_std) {
                tmpwork = mktemp ("tmpwork")
                imrename ("xtfr"//reference, "ptfr"//reference, ver-,
                    >& "dev$null")
                imsum ("ptfr"//reference//"[1],ntfr"//reference//"[1]", 
                    tmpwork, title="", hparams="", pixtype="", calctype="",
                    option="sum", low_reject=0., high_reject=0., verbose-, 
                    >& "dev$null") 
                wmef (tmpwork, "xtfr"//reference, extname="SCI",
                    phu="ptfr"//source//"[0]", verbose-, >& "dev$null")
                imdelete (tmpwork, ver-, >& "dev$null")
                imdelete ("ptfr"//reference, ver-, >& "dev$null")
                imdelete ("ntfr"//reference, ver-, >& "dev$null")
            }
        }  # end if (l_fl_negative)
        
        # Defringe, if requested.
        #
        # The defringing is intrinsicly interactive, except possibly for the 
        # case of low resolution spectra.  Thus this step cannot be used if the
        # script is to be used for batch processing.
        
        prefix = "xtfr"
        sprefix = "xtfr"
        
        if (l_fl_defringe) {
            imdelete ("dxtfr"//source, verify-, >& "dev$null")
            msdefringe ("xtfr"//source, outspec="dxtfr"//source, outpref="",
                fl_lowres=l_fl_lowres, fl_zerocut=l_fl_zerocut,
                fl_interpolate=l_fl_dfinterp, fmin=l_fmin, fmax=l_fmax,
                fl_mef+, logfile=l_logfile, verbose=l_verbose) 
                           
            if (msdefringe.status == 0)
                prefix = "dxtfr"
            else {
                # prefix stay at 'xtfr'. Issue warning and move on.
                errmsg = "Source xtfr"//source//" not de-fringeds"
                glogprint (l_logfile, "msreduce", "status", type="warning",
                    errno=1, str=errmsg, verbose=l_verbose)
            }
            
            if (l_fl_std) {
                imdelete ("dxtfr"//reference, verify-, >& "dev$null")
                msdefringe ("xtfr"//reference, outspec="dxtfr"//reference,
                    outpref="", fl_lowres=l_fl_lowres, fl_zerocut=l_fl_zerocut,
                    fl_interpolate=l_fl_dfinterp, fmin=l_fmin, fmax=l_fmax,
                    fl_mef+, logfile=l_logfile, verbose=l_verbose)
                if (msdefringe.status == 0)
                    sprefix = "dxtfr"
                else {
                    # sprefix stay at 'xtfr'. Issue warning and move on.
                    errmsg = "Standard/telluric xtfr"//reference//" not \
                        de-fringed."
                    glogprint (l_logfile, "msreduce", "status", type="warning",
                        errno=1, str=errmsg, verbose=l_verbose)
                }
            }
        }
        l_fl_plots = yes
        if (l_fl_std) {
            if (l_fl_telluric == no) l_fl_plots = no
            if (l_fl_bbody) {
                # Finally, ratio the target spectrum and standard spectrum and 
                # then multiply by a normalized blackbody spectrum to produce a
                # "calibrated" target spectrum.  This spectrum should have the 
                # proper shape but does not have any overall normalization to 
                # the magnitude of the standard in N-band.
                #
                # "msabsflux" calls "telluric" from the noao.onespec package.
                #
                # The blackbody temperature used makes little difference at 
                # these wavelengths as long as it is sufficiently high (i.e. 
                # more than 3000 K).  A default value of 10000 K is 
                # recommended.
                
                imdelete ("a"//prefix//source, verify-, >& "dev$null")
                msabsflux (prefix//source, sprefix//reference, "",
                    outimage="a"//prefix//source, outpref="",
                    outtype=l_outtype, fl_bbody=yes, bbody=l_bbody, xcorr=yes,
                    lag=10, shift=0, scale=1., dshift=0.5, dscale=0.1,
                    threshold=0.01, fl_inter=l_fl_telluric, 
                    fl_plots=l_fl_plots, logfile=l_logfile, verbose=l_verbose)
            } else {
                # Finally, ratio the target spectrum and standard spectrum and 
                # then multiply by the spectrophotometric curve for the 
                # standard.

                imdelete ("a"//prefix//source, verify-, >& "dev$null")
                msabsflux (prefix//source, sprefix//reference, l_stdname,
                    outimage="a"//prefix//source, outpref="",
                    outtype=l_outtype, fl_bbody=no, bbody=1., xcorr=yes, 
                    lag=10, shift=0, scale=1., dshift=0.5, dscale=0.1,
                    threshold=0.01, fl_inter=l_fl_telluric, 
                    fl_plots=l_fl_plots, logfile=l_logfile, verbose=l_verbose)
            }
            
            # Inspect the output spectrum
            if (access ("a"//prefix//source//".fits") && l_fl_plots) {
                splot ("a"//prefix//source//"[1]", line=1, band=1, units="",
                    options="auto wreset", xmin=INDEF, xmax=INDEF, ymin=INDEF,
                    ymax=INDEF, save_file="splot.log", graphic="stdgraph",
                    cursor="", nerrsam=0, sigma0=INDEF, invgain=INDEF,
                    function="spline3", order=1, low_reject=2., high_reject=4.,
                    niterate=10, grow=1., markrej=yes, caldir=")_.caldir",
                    fnuzero=3.68e-20)
            }
            
            gemdate ()
            gemhedit ("a"//prefix//source//".fits[0]", "MSREDUCE",
                gemdate.outdate, "UT Time stamp for MSREDUCE", delete-)
            gemhedit ("a"//prefix//source//".fits[0]", "GEM-TLM", 
                gemdate.outdate, "UT Last modification with GEMINI", delete-)
        } else {
            gemdate ()
            gemhedit (prefix//source//".fits[0]", "MSREDUCE",
                gemdate.outdate, "UT Time stamp for MSREDUCE", delete-)
            gemhedit (prefix//source//".fits[0]", "GEM-TLM", gemdate.outdate,
                "UT Last modification with GEMINI", delete-)
        }

        ngood += 1
        
nextimage:
        glogprint (l_logfile, "msreduce", "visual", type="visual",
            vistype="shortdash", verbose=l_verbose)
    }


clean:
    delete (tmpfile, ver-, >& "dev$null")
    delete (tmpinimg, ver-, >& "dev$null")
    
    if (status == 0) {
        glogprint (l_logfile, "msreduce", "status", type="string",
            str="All "//str(nsources)//" images successfully processed.",
            verbose=l_verbose)
        glogclose (l_logfile, "msreduce", fl_success+, verbose=l_verbose)
    } else {
        glogprint (l_logfile, "msreduce", "status", type="string",
            str=str(ngood)//" out of "//str(ngood+nskipped)//" successfully \
            processed.",verbose=l_verbose)
        glogclose (l_logfile, "msreduce", fl_success-, verbose=l_verbose)
    }

exitnow:
    ;

end

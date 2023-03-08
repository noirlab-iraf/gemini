# Copyright(c) 2006-2009 Association of Universities for Research in Astronomy, Inc.

procedure msslice(inspec)

char    inspec          {prompt="Name of wavelength calibrated source file (output of nstransform)"}
char    outpref         {"x",prompt="Prefix for output slice file names"}
bool    fl_calibrate    {yes,prompt="Calibrate the sliced spectra"}
char    calpref         {"a",prompt="Prefix for output slice file names"}
real    scale           {1.0,prompt="Scale factor multiplying airmass ratio [telluric]"}
real    shift           {0.0,prompt="Shift of calibration spectrum (pixels) [telluric]"}
bool    fl_inter        {no,prompt="Interactive tweaking?  [telluric]"}
char    outtype         {"fnu",prompt="Type of output spectrum: fnu|flambda|lambda*flambda"}
char    rawpath         {"",prompt="Path for input images"}
char    std             {"",prompt="Name of extracted standard spectrum file (output of nsextract)"}
char    stdname         {"",prompt="Name of standard star"}
bool    fl_bbody        {no,prompt="Use blackbody option rather than spectrophotometry"}
real    bbody           {10000.0, prompt="Temperature of calibrator for black-body fit"}
int     nspectra        {2,prompt="Number of spectra to extract over the aperture"}
char    logfile         {"",prompt="Log file name"}
bool    verbose         {yes,prompt="verbose logging?"}
int     status          {0,prompt="Exit error status: (0=good, >0=bad)"}
struct *scanfile        {"",prompt="(Internal use only)"}

begin

    char    l_inspec = ""
    char    l_outpref = ""
    char    l_calpref = ""
    char    l_outtype = ""
    char    l_rawpath = ""
    char    l_std = ""
    char    l_stdname = ""
    char    l_logfile = ""
    int     l_nspectra
    real    l_scale, l_shift, l_bbody
    bool    l_fl_calibrate, l_fl_inter, l_fl_bbody, l_verbose

    char    outname, tmpname, tmplist, paramstr, lastchar
    char    fname
    int     junk
    int     i, j, k

    junk = fscan (inspec, l_inspec)
    junk = fscan (outpref, l_outpref)
    l_fl_calibrate = fl_calibrate
    junk = fscan (calpref, l_calpref)
    l_scale = scale
    l_shift = shift
    l_fl_inter = fl_inter
    junk = fscan (outtype, l_outtype)
    junk = fscan (rawpath, l_rawpath)
    junk = fscan (std, l_std)
    junk = fscan (stdname, l_stdname)
    l_fl_bbody = fl_bbody
    l_bbody = bbody
    l_nspectra = nspectra
    junk = fscan (logfile, l_logfile)
    l_verbose = verbose

    # Initialize
    status = 0

    # Generate temp file name
    tmpname = mktemp("tmpslice")
    tmplist = mktemp("tmpfilelist")

    # Add the trailing slash to rawpath, if missing
    lastchar = substr (l_rawpath, strlen(l_rawpath), strlen(l_rawpath))
    if ((l_rawpath != "") && (lastchar != "$") && (lastchar != "/"))
        l_rawpath = l_rawpath//"/"

    # Create the list of parameter/value pairs.  One pair per line.
    # All lines combined into one string.  Line delimiter is '\n'.
    paramstr =  "inspec         = "//inspec.p_value//"\n"
    paramstr += "outpref        = "//outpref.p_value//"\n"
    paramstr += "fl_calibrate   = "//fl_calibrate.p_value//"\n"
    paramstr += "calpref        = "//calpref.p_value//"\n"
    paramstr += "scale          = "//scale.p_value//"\n"
    paramstr += "shift          = "//shift.p_value//"\n"
    paramstr += "fl_inter       = "//fl_inter.p_value//"\n"
    paramstr += "outtype        = "//outtype.p_value//"\n"
    paramstr += "rawpath        = "//rawpath.p_value//"\n"
    paramstr += "std            = "//std.p_value//"\n"
    paramstr += "stdname        = "//stdname.p_value//"\n"
    paramstr += "fl_bbody       = "//fl_bbody.p_value//"\n"
    paramstr += "bbody          = "//bbody.p_value//"\n"
    paramstr += "nspectra       = "//nspectra.p_value//"\n"
    paramstr += "logfile        = "//logfile.p_value//"\n"
    paramstr += "verbose        = "//verbose.p_value

    # Assign a logfile name if not given.  Open logfile and start log.
    # Write parameter/value pairs ("paramstr") to log.
    gloginit (l_logfile, "msslice", "midir", paramstr, fl_append+,
        verbose=l_verbose)
    if (gloginit.status != 0) {
        status = gloginit.status
        goto exitnow
    }
    l_logfile = gloginit.logfile.p_value
    
    # KL - Should check if inspec is valid
    # KL - Also 'rawpath' should not be used.  The input spectrum is 
    #       never a raw image.
    
    # Check user inputs
    if (l_nspectra < 2) {
        status = 121
        goto clean
    }
    
    
    # Set output file name
    # KL - Should check if outname exists already.  Also should add
    #      parameter for output file name instead of forcing the use of
    #      a prefix.
    outname=outpref//l_inspec

    apsum (l_rawpath//l_inspec//"[1]", output=tmpname, apertures="",
        format="onedspec", references="", profiles="", interactive=yes,
        find=yes, recenter=yes, resize=yes, edit=yes, trace=yes, fittrace=yes,
        extract=yes, extras=no, review=yes, line=INDEF, nsum=10,
        background="none", weights="none", pfit="fit1d", clean=no, skybox=1,
        saturation=INDEF, readnoise="0.", lsigma=4., usigma=4.,
        nsubaps=l_nspectra)
    ls (tmpname//"*", > tmplist)

    if (l_fl_calibrate && (l_fl_bbody == no)) {
        if ((scale <= 0.) || (abs(shift) > 5.0)) {
            glogprint (l_logfile, "msslice", "status", type="error",
                errno=status, str="Bad scale/shift values were specified.",
                verbose=yes)
            status = 121
            goto clean
        }
        
        # KL - I really, really don't like the following 3 lines.  The
        #       telluric parameters should be set in the call itself.
        unlearn ("telluric")
        telluric.tweakrms = l_fl_inter
        telluric.interactive = l_fl_inter
    }

    i = 0
    j = 0
    scanfile = tmplist
    while (fscan(scanfile,fname) != EOF) {
        i = i + 1
        wmef (fname,outname//"_"//str(i), phu=l_rawpath//l_inspec//"[0]",
            verbose=no, extname="SCI")
        glogprint (l_logfile, "msslice", "status", type="string", 
            str="Sliced spectrum file "//outname//"_"//str(i)//\
            " has been written.", verbose=l_verbose)

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
            # more than 3000 K).  A default value of 10000 K is recommended.

            imdelete (l_calpref//outname//"_"//str(i), verify-, >& "dev$null")
            msabsflux (outname//"_"//str(i), l_rawpath//l_std, "",
                outimage=l_calpref//outname//"_"//str(i), outpref="",
                outtype=l_outtype, fl_bbody=yes, bbody=l_bbody, xcorr=yes,
                lag=10, shift=l_shift, scale=l_scale, dshift=0.5, dscale=0.1,
                threshold=0.01, fl_inter=l_fl_inter, fl_plots=no,
                logfile=l_logfile, verbose=l_verbose)
            if (msabsflux.status == 0) {
                j = j + 1
            } else {
                status = status + 1
            }
        } else {
            # Finally, ratio the target spectrum and standard spectrum and 
            # then multiply by the spectrophotometric curve for the 
            # standard.

            imdelete (l_calpref//outname//"_"//str(i), verify-, >& "dev$null")
            msabsflux (outname//"_"//str(i), l_rawpath//l_std, l_stdname,
                outimage=l_calpref//outname//"_"//str(i), outpref="",
                outtype=l_outtype, xcorr=yes, lag=10, shift=l_shift,
                scale=l_scale, dshift=0.5, dscale=0.1, threshold=0.01,
                fl_inter=l_fl_inter, fl_plots=no, logfile=l_logfile,
                verbose=l_verbose)
            if (msabsflux.status == 0) {
                j = j + 1
            } else {
                status = status + 1
            }
        }
        glogprint (l_logfile, "msslice", "status", type="string", 
            str="Calibrated spectrum file "//\
            l_calpref//outname//"_"//str(i)//" has been written.",
            verbose=l_verbose)

        # Inspect the output spectrum if interactive tweaking is requested
        if (l_fl_inter && access(l_calpref//outname//str(i)//".fits")) {
            splot (l_calpref//outname//"_"//str(i)//"[1]", line=1, band=1,
                units="", options="auto wreset", xmin=INDEF, xmax=INDEF,
                ymin=INDEF, ymax=INDEF, save_file="splot.log",
                graphic="stdgraph", cursor="", nerrsam=0, sigma0=INDEF,
                invgain=INDEF, function="spline3", order=1, low_reject=2.,
                high_reject=4., niterate=10, grow=1., markrej=yes,
                caldir=")_.caldir", fnuzero=3.68e-20)
        }

    }

clean:

    delete (tmplist, ver-, >& "dev$null")
    delete ("tmpslice*", ver-, >& "dev$null")

    if (status == 0) {
        glogprint (l_logfile, "msslice", "status", type="string",
            str="All "//str(nspectra)//" slices have been successfully \
            processed.", verbose=l_verbose)
        glogclose (l_logfile, "msslice", fl_success+, verbose=l_verbose)
    } else {
        glogprint (l_logfile, "msslice", "status", type="string",
            str=str(j)//" out of "//str(nspectra)//" slices were successfully \
            processed.",verbose=l_verbose)
        glogclose (l_logfile, "msslice", fl_success-, verbose=l_verbose)
    }

exitnow:
    ;

end

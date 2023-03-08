# Copyright(c) 2006-2009 Association of Universities for Research in Astronomy, Inc.
#
# Original Author: Kevin Volk  (March 2006)
#

procedure nfdispc (inimage, cmin, cmax)

char    inimage         {prompt="Input NIFS datacube file"} # OLDP-1-input-primary
real    cmin            {prompt="Starting spectral pixel or wavelength"}
real    cmax            {prompt="Ending spectral pixel or wavelength"}
bool    fl_wavelength   {no,prompt="Start/cmax values in wavelength"}
bool    fl_pause        {no,prompt="Pause and run imexam on each image"}
int     frame           {1,prompt="frame to be written into [display]"}
bool    zscale          {no,prompt="display range of greylevels near median [display]"}
bool    zrange          {yes,prompt="display full image intensity range [display]"}
real    z1              {0.,prompt="minimum greylevel to be displayed [display]"}
real    z2              {0.,prompt="maximum greylevel to be displayed [display]"}
char    ztrans          {"linear",prompt="greylevel transformation (linear|log|none|user) [display]"}
char    lutfile         {"",prompt="file containing user defined look up table [display]"}
char    logfile         {"",prompt="Log file name"}
bool    verbose         {yes,prompt="verbose logging?"}
int     status          {0,prompt="Exit error status: (0=good, >0=bad)"}

begin

    char    l_inimage, l_lutfile, l_ztrans, l_logfile
    real    l_cmin, l_cmax, l_z1, l_z2
    real    wavel0, delwl, pixel0, pixel1
    int     imin, imax, k, naxis[3], iaxis, l_frame
    bool    l_fl_pause, l_fl_wavelength, l_zscale, l_zrange, l_verbose
    char    paramstr, tmpstring

    l_inimage = inimage
    l_cmin = cmin
    l_cmax = cmax
    l_fl_wavelength = fl_wavelength
    l_fl_pause = fl_pause
    l_frame = frame
    l_zscale = zscale
    l_zrange = zrange
    l_z1 = z1
    l_z2 = z2
    l_ztrans = ztrans
    l_lutfile = lutfile
    l_logfile = logfile
    l_verbose = verbose

    # Create the list of parameter/value pairs.  One pair per line.
    # All lines combined into one string.  Line delimiter is '\n'.
    paramstr  =  "inimage        = "//inimage.p_value//"\n"
    paramstr += "cmin            = "//cmin.p_value//"\n"
    paramstr += "cmax            = "//cmax.p_value//"\n"
    paramstr += "fl_wavelength   = "//fl_wavelength.p_value//"\n"
    paramstr += "fl_pause        = "//fl_pause.p_value//"\n"
    paramstr += "zscale          = "//zscale.p_value//"\n"
    paramstr += "zrange          = "//zrange.p_value//"\n"
    paramstr += "z1              = "//z1.p_value//"\n"
    paramstr += "z2              = "//z2.p_value//"\n"
    paramstr += "ztrans          = "//ztrans.p_value//"\n"
    paramstr += "lutfile         = "//lutfile.p_value//"\n"
    paramstr += "verbose        = "//verbose.p_value

    status = 0

    # Assign a logfile name if not given.  Open logfile and start log.
    # Write parameter/value pairs ("paramstr") to log.
    gloginit (l_logfile, "nfdispc", "nifs", paramstr, fl_append+,
        verbose=l_verbose)
    if (gloginit.status != 0) {
        status = gloginit.status
        goto exitnow
    }
    l_logfile = gloginit.logfile.p_value

    if (no == imaccess(l_inimage)) {
        glogprint (l_logfile, "nfdispc", "status", type="error", errno=121,
            str="The input image was not found.", verbose+)
        status = 1
        goto clean
    }

    imgets (l_inimage//"[0]", "INSTRUME", >& "dev$null")
    if (imgets.value != "NIFS") {
        glogprint (l_logfile, "nfdispc", "status", type="error", errno=121,
            str="The instrument is not NIFS.", verbose+)
        status = 1
        goto clean
    }

    imgets (l_inimage//"[0]", "NSTRANSF", >& "dev$null")
    if (imgets.value == "") {
        glogprint (l_logfile, "nfdispc", "status", type="error", errno=121,
            str="The file is not an output file of NSTRANSFORM.", verbose+)
        status = 1
        goto clean
    }

    if (l_fl_wavelength) {

        imgets (l_inimage//"[1]", "CTYPE3", >& "dev$null")
        if (imgets.value != "LINEAR") {
            glogprint (l_logfile, "nfdispc", "status", type="error", errno=121,
                str="The spectrum is not linearized.", verbose+)
            status = 1 
            goto clean
        }

        imgets (l_inimage//"[1]", "CRVAL3", >& "dev$null")
        print (imgets.value) | scanf ("%f",wavel0)
        imgets (l_inimage//"[1]", "CD3_3", >& "dev$null")
        print (imgets.value) | scanf ("%f",delwl)
        imgets (l_inimage//"[1]", "CRPIX3", >& "dev$null")
        print (imgets.value) | scanf ("%f",pixel0)
        if ((wavel0 == 0.) || (delwl < 0.5) || (pixel0 < 1.)) {
            glogprint (l_logfile, "nfdispc", "status", type="error", errno=121,
                str="Could not get the wavelength parameters.", verbose+)
            status = 1
            goto clean
        }
        pixel1 = (l_cmin-wavel0) / delwl + pixel0
        imin = int(pixel1+0.5)
        pixel1 = (l_cmax-wavel0) / delwl + pixel0
        imax = int(pixel1+0.5)
        
    } else {
        imin = int(l_cmin)
        imax = int(l_cmax)
    }

    if (imin > imax) {
        k = imax
        imax = imin
        imin = k
    }
    imgets (l_inimage//"[1]", "i_naxis")
    iaxis = int(imgets.value)
    if (iaxis != 3) {
        tmpstring = "Image extension 1 has NAXIS = "//str(iaxis)//".\n\
            The expected number is 3."
        glogprint (l_logfile, "nfdispc", "status", type="error", errno=121,
            str=tmpstring,verbose+)
        status = 1
        goto clean
    }
    else {
        naxis[1] = 0
        naxis[2] = 0
        naxis[3] = 0

        for (k=1; k <= iaxis; k+=1) {
            tmpstring = "naxis"//str(k)
            imgets (l_inimage//"[1]", tmpstring)
            naxis[k] = int(imgets.value)
        }
    }

    if (imin > naxis[3] || imax < 1) {
        if (l_fl_wavelength) {
            tmpstring = "Specified wavelength range "//str(l_cmin)//" to "//\
                str(l_cmax)//" is out of range."
        } else {
            tmpstring = "Specified pixel range "//str(imin)//" to "//\
                str(imax)//" is out of range."
        }
        glogprint (l_logfile, "nfdispc", "status", type="error", errno=121,
            str=tmpstring, verbose+)
        status = 1
        goto clean
    }

    if (imax > naxis[3]) imax = naxis[3]
    if (imin < 1) imin = 1

    # loop to scan through planes of nifs cube and display them
    for (k=imin; k <= imax ; k+=1) {
        tmpstring = l_inimage//"[1][*,*,"//k//"]"
        if (l_verbose) {
            if (l_fl_wavelength) {
                pixel1 = wavel0 + (real(k)-pixel0) * delwl
                paramstr = "Displaying wavelength "//str(pixel1)//\
                    " (pixel "//str(k)//")"
            } else {
                paramstr = "Displaying pixel plane "//str(k)
            }
            glogprint (l_logfile, "nfdispc", "status", type="string", errno=121,
                str=paramstr, verbose+)
            display (tmpstring, l_frame, zscale=l_zscale, zrange=l_zrange,
                z1=l_z1, z2=l_z2, ztrans=l_ztrans, lutfile=l_lutfile, 
                erase-, border_erase-, nsample=5000, fill+)
        } else {
            display (tmpstring, l_frame, zscale=l_zscale, zrange=l_zrange,
                z1=l_z1, z2=l_z2, ztrans=l_ztrans, lutfile=l_lutfile, 
                erase-, border_erase-, nsample=5000, fill+, >& "dev$null")
        }
        if (l_fl_pause) {
            paramstr = "Image plane "//str(k)
            imexam (tmpstring, l_frame, paramstr, logfile="",
                display="display(image='$1',frame=$2)", use_display+)
        }  
    }

clean:

    if (status == 0) {
        glogclose (l_logfile, "nfdispc", fl_success+, verbose=l_verbose)
    } else {
        glogclose (l_logfile, "nfdispc", fl_success-, verbose=l_verbose)
    }

exitnow:
    ;

end

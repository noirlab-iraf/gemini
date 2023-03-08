# Copyright(c) 2005-2017 Association of Universities for Research in Astronomy, Inc.
#
# GMSKXY - Generate x,y MOS mask coordinates from a list of RA and Dec.
#
# Original author: Rachel Johnson
#
# 2006-05-29 added new transformation x(xi,eta) y(xi,eta)
# commented out transformdb
#
# 10/07 A. Wong: Removed hemisphere parameter and replaced with Instrument
#                Added placeholder variables for Flamingos II
#                Moved input ascii table check from gmskcreate to here
#                  - if not all optional entries present in table will not
#                    set all to the defaults
#                  - For objects with priority = 0 (acquisition objects,
#                    the script will set slit size = 2", and the tilt and
#                    position = 0
#
# 2017-03-27 (mischa):
#                    Transforming the 1x1 coordinate system to the new
#                    Hamamatsu and 2x2 binning, including an correction for a 6"
#					 offset	for GMOS-S that has been present for a long time.
# 2017-04-21 (mischa):
#                    Added pseudo-image support for F2
# 2017-05-16 (mischa):
#                    Updated pseudo-image transformations for
#                    F2, GMOS-N and GMOS-S




procedure gmskxy (indata, instrument, rafield, decfield, pa)
char    indata      {prompt="Input file containing info (id,ra,dec,mag + optional) for spectroscopy candidates"}
char    instrument  {enum="gmos-n|gmos-s|flamingos2", prompt="Instrument (gmos-n|gmos-s|flamingos2)"}
real    rafield     {prompt="RA value of field center"}
real    decfield    {prompt="Dec value of field center"}
real    pa          {prompt="PA of field if required by OT"}
char    iraunits    {"hours", enum="hours|degrees", prompt="Spectroscopy candidate RA units (hours|degrees)"}
char    fraunits    {"hours", enum="hours|degrees", prompt="Field centre RA units (hours|degrees)"}
char    outcoords   {"", prompt="Output file for x,y of spectroscopy candidates"}
char    outtab      {"", prompt="Output object table for gmmps"}
real    slitszx     {1.0,prompt="Default slit width in arcsecs"}
real    slitszy     {5.0,prompt="Default slit length in arcsecs"}
char    logfile     {"", prompt="Logfile name"}
char    glogpars    {"", prompt="Logging preferences"}
bool    verbose     {yes, prompt="Verbose"}
bool    fl_debug    {no, prompt="Print debugging information"}
int     status      {0, prompt="Exit status (0=good)"}
struct  *scanfile   {"", prompt = "Internal use only"}

begin

    char    l_indata = ""
    char    l_instrument = ""
    char    l_iraunits = ""
    char    l_fraunits = ""
    char    l_outcoords = ""
    char    l_outtab = ""
    char    l_logfile = ""
    real    l_rafield, l_decfield, l_pa
    bool    l_verbose
    bool    l_fl_debug

    char    paramstr, errmsg
    char    tmpcrds, tmpdata, tmpcd, tmpinput
    char    outcoofile, outtabfile
    int     i, j, k, n, junk
    real    pi
    real    l_slitszx, l_slitszy
    real    ra0, dec0, radiff, rac, decc
    real    xi, eta, pxi, peta
    real    xgemini, ygemini, theta
    struct  tmpstruct

    char    id, pri_char
    int     pri
    real    racand, deccand, mag, ssx, ssy, st, spy

    ###########################################################################
    ## these are the coefficients used for the transformation

#    real    xcfs_GMOS_S[15] = 3101.31, 0.000178737, -1.14711e-06, \
#                -3.23254e-10, 3.51214e-12, 2.82153, -3.38154e-07, \
#                9.21113e-09, -1.50887e-12, 1.83446e-06, -2.82184e-10, \
#                -5.98567e-13, 9.08167e-09, 1.59136e-12, -1.45280e-12
#    real    ycfs_GMOS_S[15] = 2367.49, 2.82072, -4.06813e-07, 8.77777e-09, \
#                1.06523e-12, 0.000198779, 6.17942e-08, -6.23632e-10, \
#                7.95186e-13, 1.25869e-06, 9.70030e-09, 6.03889e-13, \
#                -3.36675e-10, -2.67149e-13, -2.17115e-12

	real    xcfs_GMOS_S[15] = 1563.75, 0.00318698, 5.42125e-08, -9.7724e-12, \
	            0.0, -1.2872, 1.08287e-07, -4.65546e-09, 0.0, \
				-7.36554e-08, 3.00154e-10, 0.0, -3.93593e-09, 0.0, 0.0

	real    ycfs_GMOS_S[15] = 1043.44, -1.28762, 4.66584e-08, -4.02699e-09, \
	            0.0, -0.00338725, -3.20178e-08, -2.34703e-11, 0.0, \
				-1.57816e-07, -3.88237e-09, 0.0, 2.60055e-10, 0.0, 0.0
	
	
    # Coefficients which are equal to zero in the North list are necessary
    # because 3rd order polynomial gave best fit for Northern transform whereas
    # 4th order was best fit for South. Below we use the same transform
    # equation for both North and South.

#    real    xcfs_GMOS_N[15] = 3108.06, 0.00131027, -4.17610e-08, \
#                -2.32551e-09, 0.0, 2.82696, -3.79727e-07, 8.23253e-09, 0.0, \
#                3.99659e-07, -2.47552e-10, 0.0, 1.04171e-08, 0.0, 0.0
#    real    ycfs_GMOS_N[15] = 2304.06, 2.82767, -3.15178e-07, 6.56689e-09, \
#                0.0, -0.000864297, -6.18804e-07, 1.38929e-09, 0.0, \
#                -1.60831e-07, 1.00271e-08, 0.0, 1.92364e-09, 0.0, 0.0

    real    xcfs_GMOS_N[15] = 1569.4, 0.000151798, 1.30509e-08, 1.47095e-10, \
	            0.0, -1.2743, -3.31558e-07, -3.94706e-09, 0.0, \
				1.6597e-07, 6.71758e-10, 0.0, -3.59396e-09, 0.0, 0.0
    real    ycfs_GMOS_N[15] = 1043.15, -1.27425, -4.01144e-07, -3.58639e-09, \
	            0.0, -0.000207071, -2.80796e-08, -3.49208e-10, 0.0, \
	            1.39364e-07, -3.48446e-09, 0.0, -2.38642e-10, 0.0, 0.0
	
    # Flamingos II; a third order description is sufficient for residuals
	# of ~0.1 pixels, min / max +/-0.3 pixel
    real    xcfs_Flamingos[15] = 1023.98, 1.15006, -4.25093e-07, 1.90535e-09, \
	            0.0, -1.75298e-05, -4.3626e-07, 6.90963e-12, 0.0, \
			    -8.29324e-08, 2.07379e-09, 0.0, 1.47971e-11, 0.0, 0.0

    real    ycfs_Flamingos[15] = 1024.08, 0.000106509, 3.31497e-08, -1.53915e-10, \
	            0.0, 1.15023, -4.77258e-07, 1.67071e-09, 0.0, \
				-3.59315e-07, -3.965e-10, 0.0, 1.55948e-09, 0.0, 0.0
	
    real    xcfs[15]
    real    ycfs[15]

    ###########################################################################

    junk = fscan (indata, l_indata)
    junk = fscan (instrument, l_instrument)
    l_rafield = rafield
    l_decfield = decfield
    l_pa = pa
    junk = fscan (iraunits, l_iraunits)
    junk = fscan (fraunits, l_fraunits)
    junk = fscan (outcoords, l_outcoords)

	# GMOS and F2 have 90 deg different dispersion directions
	if (l_instrument == "gmos-n" || l_instrument == "gmos-s") {
	    l_slitszx = slitszx
        l_slitszy = slitszy
    } else if (l_instrument == "flamingos2") {
	    l_slitszx = slitszy
        l_slitszy = slitszx
    }
	
    junk = fscan (outtab, l_outtab)
    junk = fscan (logfile, l_logfile)
    l_verbose = verbose
    l_fl_debug = fl_debug

    # Initialize
    status = 0
    pi = 4. * atan2(1,1)

    # Create names for temp files
    tmpcrds = mktemp ("tmpcrds")
    tmpdata = mktemp ("tmpdata")
    tmpinput = mktemp("tmpinput")
    tmpcd = mktemp ("tmpcd")

    # Create the list of parameter/value pairs. One pair per line.
    # All lines combined into one string.  Line delimiter is '\n'.
    paramstr =  "indata         = "//indata.p_value//"\n"
    paramstr += "instrument     = "//instrument.p_value//"\n"
    paramstr += "rafield        = "//rafield.p_value//"\n"
    paramstr += "decfield       = "//decfield.p_value//"\n"
    paramstr += "pa             = "//pa.p_value//"\n"
    paramstr += "iraunits       = "//iraunits.p_value//"\n"
    paramstr += "fraunits       = "//fraunits.p_value//"\n"
    paramstr += "outcoords      = "//outcoords.p_value//"\n"
    paramstr += "outtab         = "//outtab.p_value//"\n"
    paramstr += "logfile        = "//logfile.p_value//"\n"
    paramstr += "verbose        = "//verbose.p_value//"\n"
    paramstr += "fl_debug       = "//fl_debug.p_value

    # Assign a logfile name if not given.  Open logfile and start log.
    # Write parameter/value pairs ("paramstr") to log.
    gloginit (l_logfile, "gmskxy", "mostools", paramstr, fl_append+,
        verbose=l_verbose)
    if (gloginit.status != 0) {
        status = gloginit.status
        goto exitnow
    }
    l_logfile = gloginit.logfile.p_value


    # Check input parameters

    # Input RA can be in hours or degrees. If in hours then fraunits is
    # 'hours', if in degrees is is 'degrees'
    # IRAF automatically convert hh:mm:ss.ss to a 'real'.
    # Therefore at this point RA is in format:
    #      hh.hhh (hours)
    #   or dd.ddd (degrees)
    # And Dec is in format
    #      dd.ddd (degrees)
    # RA and Dec must be transformed to radians before next stage

    glogprint (l_logfile, "gmskxy", "engineering", type="string",
        str="Requested field center: RA("//l_fraunits//") = "//l_rafield//\
        ", Dec = "//l_decfield, verbose=l_verbose)

    if (l_fraunits == "hours")
        ra0 = l_rafield * (360. / 24.) * (pi/180.0)
    else
        ra0 = l_rafield * (pi/180.0)

    dec0 = l_decfield * (pi/180.0)

    ################ check input data file ###################################
    # remove commented lines
    # remove lines with < 4 values
    # if 'priority', 'slitsize_x', 'slitsize_y', 'slittilt', 'slitpos_y'
    #    not set in input, then use default values

    scanfile = l_indata

    # for each line in the 'indata' file, check it has the minimum required
    # number of values, and write out the data to a tmpdata file. Each line
    # must contain id, ra, dec, mag.  If the entry does not have the minimum
    # number of values then it is not used.  The line can optionally also
    # contain priority, slitsize_x, slitsize_y, slittilt, slitpos_y.
    #
    # Check that the priority is in the range 0-3.
    #
    # If none of the optional values are set in the input file, then they are
    # set to default values here. If some but not all of the optional values
    # are set then none of them are used, and they are set to default values.
    # If all of the optional values are set then they are used.

    i = 0  # i counts the total number of lines in the input file,
           # including commented out ones
    j = 0  # j counts the total number of useful lines in the input file


    glogprint (l_logfile, "gmskxy", "engineering", type="string",
        str="Checking input data file for compliance", verbose=l_verbose)

    pri_char = INDEF
    while (fscan (scanfile,id,rac,decc,mag,pri_char,ssx,ssy,st,spy) != EOF) {
        i += 1
        if (substr(id,1,1) != "#") { # ignore commented out lines in input
            n = nscan()

            if (n < 4) {
                # if n < 4 then the required parameters (id,ra,dec,mag) are not
                # set and this line is ignored

                errmsg = "Problem with input file "//l_indata//\
                    ". Line number "//i//", object id = "//id//" does not \
                    contain at least 4 values (id, ra, dec, mag), so will not \
                    be included in the output files."
                status = 121
                glogprint (l_logfile, "gmskxy", "status", type="warning",
                    errno=status,  str=errmsg, verbose+)

            } else if (n > 9) {
                # if n >9  then too many parameters in the line have been
                # set and this line is ignored

                errmsg = "Problem with input file "//l_indata//\
                    ". Line number "//i//", object id = "//id//" contains \
                    more than the 9 columns of expected input, so will not \
                    be included in the output files."
                status = 121
                glogprint (l_logfile, "gmskxy", "status", type="warning",
                    errno=status,  str=errmsg, verbose+)

            } else {
                # Assumes that input data is order id, ra, dec, mag, priority,
                # slit size x, slit size y, slit tilt, and slit position y.
                # Based on the number of entries in each line, it will
                # enter default values for any missing parameters

                if (n == 4) {
                    pri = 1
                    ssx = l_slitszx
                    ssy = l_slitszy
                    st = 0
                    spy = 0
                    glogprint (l_logfile, "gmskxy", "status", type="string",
                        str="Setting priority, slit size x & y, slit tilt, \
                        and y-position in line "//i, verbose=l_verbose)

                } else if (n == 5) {
                    ssx = l_slitszx
                    ssy = l_slitszy
                    st = 0
                    spy = 0
                    glogprint (l_logfile, "gmskxy", "status", type="string",
                        str="Setting slit size x & y, slit tilt, and \
                        y-position in line "//i, verbose=l_verbose)
                } else if (n == 6) {
                    ssy = l_slitszy
                    st = 0
                    spy = 0
                    glogprint (l_logfile, "gmskxy", "status", type="string",
                        str="Setting slit size y, slit tilt, and "//\
                        y-position//"in line "//i, verbose=l_verbose)

                } else if (n == 7) {
                    st = 0
                    spy = 0
                    glogprint (l_logfile, "gmskxy", "status", type="string",
                        str="Setting slit tilt, and y-position in line "//i,
                        verbose=l_verbose)
                } else if (n == 8) {
                    spy = 0
                    glogprint (l_logfile, "gmskxy", "status", type="string",
                        str="Setting y-position in line "//i,
                        verbose=l_verbose)
                }

                # Convert pri_char to integer
                # Check for pri_char == x; else leave it to the pri check later
                # on to issue a warning. x or X is a valid character to exclude
                # a source.
                pri = INDEF
                if (isindef(pri_char)) {
                    glogprint (l_logfile, "gmskxy", "status", \
                        type="string", str="Priority is set to INDEF"//\
                        "for id = "//id//": setting to 1", verbose=l_verbose)
                    pri_char = "1"
                }
                print (pri_char) | scan (pri)
                if (strupr(pri_char) == "X") {
                    glogprint (l_logfile, "gmskxy", "status", \
                        type="string", str="Priority is set to "//pri_char//\
                        "for id = "//id, verbose=l_verbose)
                } else if ((pri != 0) && (pri != 1) && (pri != 2) && \
                        (pri != 3)) {
                    # check to see if priority is set to 0, 1, 2, or 3
                    # if not, automatically sets to one
                    errmsg = "Problem with input file "//l_indata//". The \
                        priority value (\""//pri_char//"\") in Line number "//\
                        i//", object \
                        id = "//id//" is outside the allowed range (0-3). The \
                        priority value for this line will be set to 1."
                    status = 121
                    glogprint (l_logfile, "gmskxy", "status", type="warning",
                        errno=status, str=errmsg, verbose+)
                    pri_char = str(1)
                    pri = 1
                }

                if (pri == 0) {
                    # objects with priority of 0 are acquisition objects and
                    # must be set to the following values to have the proper
                    # slit properties present when mask is cut and observed.

                    ssx = 2
                    ssy = 2
                    st = 0
                    spy = 0
                    glogprint (l_logfile, "gmskxy", "status", type="string",
                       str="Setting acquisition object to have default slit \
                       properties in line "//i, verbose=l_verbose)
                }

                print (id, " ", rac, decc, mag, pri_char, " ", ssx, ssy, st, \
                    spy, >> tmpinput)
                j += 1
            }
        }
        pri_char = INDEF
    }

    # reset status (the messages in the while loop above were warnings)
    status = 0

    # check that there are enough useful lines in the input data file
    # for gmskxy require >1 line

    glogprint (l_logfile, "gmskxy", "status", type="string",
        str="Number of useful lines in input data file = "//j,
        verbose=l_verbose)
    if (j == 0) {
        errmsg = "No useful lines in input file."
        status = 121
        glogprint (l_logfile, "gmskxy", "status", type="error",
            errno=status, str=errmsg, verbose+)
        goto exit
    }


    # Get angle for transformation.
    # Input angle is that used in OT
    # For GMOS, for PA=0 we have North Down and East Left (180-skyPA)
    # For F2, for PA=0 we have North Right and East Up 
    # Require skypa (theta) for transformation

	theta = l_pa
#	if (l_instrument == "gmos-n") {
#        theta = 180. + l_pa
#    } 

#	if (l_instrument == "gmos-s") {
#        theta = 180. + l_pa
#    } 

    if (instrument == "flamingos2") {
        # No idea why, but for F2 we must negate the angle to get the OT orientation
        theta = -1.0 * l_pa
    }

    if (theta >= 360.)
        theta = theta - 360.

    if (theta <= -360.)
        theta = theta + 360.

    theta = theta * (pi/180.)

    if (instrument == "gmos-n"){
        # The syntax 'xcfs = xcfs_north' does not work under ECL for some
        # reason.  The for-loop is less elegant but it works.
        for (k=1; k<16; k=k+1) {
            xcfs[k] = xcfs_GMOS_N[k]
            ycfs[k] = ycfs_GMOS_N[k]
        }
    } else if (instrument == "gmos-s") {
        for (k=1; k<16; k=k+1) {
            xcfs[k] = xcfs_GMOS_S[k]
            ycfs[k] = ycfs_GMOS_S[k]
        }
    } else if (instrument == "flamingos2") {
        for (k=1; k<16; k=k+1) {
            xcfs[k] = xcfs_Flamingos[k]
            ycfs[k] = ycfs_Flamingos[k]
        }
    }
    glogprint (l_logfile, "gmskxy", "status", type="string",
        str="Using transformations for "//l_instrument, verbose=l_verbose)

    ############## find gmos xy of spectroscopy candidates ##################
    # take RA Dec of each candidate and calculate standard coordinates

    glogprint (l_logfile, "gmskxy", "status", type="string",
        str="Finding "//l_instrument//" X,Y of spectroscopy candidates",
        verbose=l_verbose)
    glogprint (l_logfile, "gmskxy", "task", type="string",
        str="Output "//l_instrument//" X,Y to "//l_outcoords, \
        verbose=l_verbose)

    # l_indata contains id, racand, deccand, mag, pri, ssx, ssy, st, spy
    # the input info file is checked by gmskcreate and only good lines are
    # in input file here

    scanfile = tmpinput
    i = 1

    # First read the lines from the input info file into an array, info,
    # ignoring any lines which are commented out
    while (fscan (scanfile, id, racand, deccand, mag, pri_char, ssx, ssy, \
        st, spy) != EOF) {

        # IRAF automatically convert hh:mm:ss.ss to a 'real'.
        # Therefore at this point RA is in format:
        #      hh.hhh (hours)
        #   or dd.ddd (degrees)
        # And Dec is in format
        #      dd.ddd (degrees)
        # RA and Dec must be transformed to radians before next stage

        if (l_fl_debug) {
            printf ("Candidate #%d : RA(%s)=%.6f , Dec=%.6f\n", i, l_iraunits,
                racand, deccand) | scan (tmpstruct)
            glogprint (l_logfile, "gmskxy", "engineering", type="string",
                str=tmpstruct, verbose=yes)
        }

        if (l_iraunits == "hours")
            racand = racand * 360. / 24.

        racand = racand * (pi/180.0)
        deccand = deccand * (pi/180.0)

        radiff = racand - ra0

        # Transform to standard coordinates
        xi = (cos(deccand) * sin(radiff)) / \
            ((sin(dec0) * sin(deccand)) + \
            (cos(dec0) * cos(deccand) * cos(radiff)))
        eta = ((cos(dec0) * sin(deccand)) - \
            (sin(dec0) * cos(deccand) * cos(radiff))) / \
            ((sin(dec0) * sin(deccand)) + \
            (cos(dec0) * cos(deccand) * cos(radiff)))

        if (l_fl_debug) {
            printf ("Candidate #%d : xi=%.10f , eta=%.10f\n", i, xi, eta) |\
                scan (tmpstruct)
            glogprint (l_logfile, "gmskxy", "engineering", type="string",
                str=tmpstruct, verbose=yes)
        }
        xi = xi * 1e6
        eta = eta * 1e6

        # And rotate to axes at angle theta

        pxi = (xi * cos(theta)) - (eta * sin(theta))
        peta = (xi * sin(theta)) + (eta * cos(theta))

        if (l_fl_debug) {
            printf ("Candidate #%d : pxi=%.10f , peta=%.10f\n", i, pxi, \
                peta) | scan (tmpstruct)
            glogprint (l_logfile, "gmskxy", "engineering", type="string",
                str=tmpstruct, verbose=yes)
        }

        xgemini = xcfs[1]        + xcfs[2]*peta         + xcfs[3]*peta**2     \
              + xcfs[4]*peta**3  + xcfs[5]*peta**4      + xcfs[6]*pxi         \
              + xcfs[7]*pxi*peta + xcfs[8]*pxi*peta**2  + xcfs[9]*pxi*peta**3 \
              + xcfs[10]*pxi**2  + xcfs[11]*pxi**2*peta + xcfs[12]*pxi**2*peta**2 \
              + xcfs[13]*pxi**3  + xcfs[14]*pxi**3*peta + xcfs[15]*pxi**4

        ygemini = ycfs[1]        + ycfs[2]*peta         + ycfs[3]*peta**2     \
              + ycfs[4]*peta**3  + ycfs[5]*peta**4      + ycfs[6]*pxi         \
              + ycfs[7]*pxi*peta + ycfs[8]*pxi*peta**2  + ycfs[9]*pxi*peta**3 \
              + ycfs[10]*pxi**2  + ycfs[11]*pxi**2*peta + ycfs[12]*pxi**2*peta**2 \
              + ycfs[13]*pxi**3  + ycfs[14]*pxi**3*peta + ycfs[15]*pxi**4

        # EVIL HACK because we don't know how the transformation
        # above was determined.
        # Need to fold in the new pixel scale and detector geometry for
        # the Hamamatsus.
        # The constant offset ensures that the RA/DEC requested is
        # at the image center.

        # All of the above should be remeasured properly. It looks
        # like this could be a lot more accurate...
        # The numeric factor at the end is the old EEV pixscale divided by
        # the new Hamamatsu 2x2 pixscale
#        if (instrument == "gmos-s") {
#        	xgemini = 153.0 + xgemini * 0.0730 / 0.1600
#        	ygemini = -36.0 + ygemini * 0.0730 / 0.1600
#        }

#        if (instrument == "gmos-n") {
#        	xgemini = 170.0 + xgemini * 0.0727 / 0.1614
#        	ygemini =   3.0 + ygemini * 0.0727 / 0.1614
#        }

        if (l_fl_debug) {
            printf ("Candidate #%d : deltax=%.6f , deltay=%.6f\n", i,
                xgemini, ygemini) | scan (tmpstruct)
            glogprint (l_logfile, "gmskxy", "engineering", type="string",
                str=tmpstruct, verbose=yes)
        }

        printf ("%.3f %.3f\n", xgemini, ygemini, >> outcoords)

        # Write ID, RA, Dec, x, y, mag, slitsize_x, slitsize_y, slittilt
        # to output temp file for later writing to gmmps OT file
        # The rac has to be in decimal hours, decc in degrees

        rac = racand / (pi/180.0) / (360./24.)
        decc = deccand / (pi/180.0)

        print (id, " ", rac, decc, xgemini, ygemini, mag, pri_char, " ", ssx, \
            ssy, st, spy, >> tmpdata)
    }
    scanfile = ""

    ############### Turn into object table required for gmmps #################

    # Need a file containing ID, RA, Dec, x_ccd, y_ccd, MAG
    # and optionally slitsize_x, slitsize_y, slittilt, slitpos_y
    # This information is in tmpdata

    glogprint (l_logfile, "gmskxy", "status", type="string",
        str="Creating FITS object table "//l_outtab//" for use in gmmps",
        verbose=l_verbose)

    # Create table column file
    print ("ID ch*20 \"\" ## \nRA r \"\" H \nDEC r \"\" deg \n\
        x_ccd r \"\" pixels \ny_ccd r \"\" pixels \nMAG r \"\" mag \n\
        priority ch*1 \nslitsize_x r \"\" arcsec \n\
        slitsize_y r \"\" arcsec \nslittilt r \"\" degrees \n\
        slitpos_y r \"\" arcsec", > tmpcd)

    tcreate (l_outtab, tmpcd, tmpdata, uparfile="", nskip=0, nlines=0,
        nrows=0, hist=yes, extrapar=5, tbltype='default', extracol=0)

exit:
    scanfile = ""
    delete (tmpcd, verify-, >& "dev$null")
    delete (tmpdata, verify-, >& "dev$null")
    delete (tmpcrds, verify-, >& "dev$null")
    delete (tmpinput, verify-, >& "dev$null")

    if (status == 0)
        glogclose (l_logfile, "gmskxy", fl_success+, verbose=l_verbose)
    else
        glogclose (l_logfile, "gmskxy", fl_success-, verbose=l_verbose)

exitnow:
    ;

end

# Copyright(c) 2006-2017 Association of Universities for Research in Astronomy, Inc.
#
# GMSKCREATE - Produce the required inputs for gmmps from a non-Gemini image with 
#              good astrometry
#
# This cl script is a wrapper for the gmskyxy and gmskimg routines
#
# Original author: Rachel Johnson
#
# 10/07   A. Wong: Removed hemisphere and replaced with instrument
#                  Added support of Flamingos II w/ empty variable placeholders
#                  Moved WCS check to gmskimg
#                  Moved input ascii file check to gmskxy
#                  Moved header keyword adding to output image to gmskimg
# 2017-04-21 (mischa): Added pseudo-image support for F2
#

procedure gmskcreate (indata, gprgid, instrument, rafield, decfield, pa)

char    indata      {prompt="Input file containing information for spectroscopy candidates"}
char    gprgid      {prompt="Your Gemini program ID (e.g. GN-2007A-Q-4)"}
char    instrument  {enum="gmos-n|gmos-s|flamingos2", prompt="Instrument (gmos-n|gmos-s|flamingos2)"}
real    rafield     {prompt="RA value of field center"}
real    decfield    {prompt="Dec value of field center"}
real    pa          {prompt="PA of field if required by OT"}

bool    fl_getim    {yes, prompt="Do you require an output fake image ?"}
char    inimage     {"", prompt="Image to be transformed (if MEF, specify extension)"}
char    outimage    {"", prompt="Output transformed image"}

bool    fl_getxy    {yes, prompt="Do you require output x,y and GMMPS object table?"}
char    outcoords   {"", prompt="Output file for x,y of spectroscopy candidates"}
char    outtab      {"", prompt="Output object table for gmmps"}


char    iraunits    {"hours", enum="hours|degrees", prompt="Spectroscopy candidate RA units (hours|degrees)"} 
char    fraunits    {"hours", enum="hours|degrees", prompt="Field centre RA units (hours|degrees)"} 
real    slitszx     {1.0,prompt="Default slit width in arcsecs"}
real    slitszy     {5.0,prompt="Default slit length in arcsecs"}

char    logfile     {"", prompt="Logfile name"}
pset    glogpars    {"", prompt="Logging preferences"}
bool    fl_inter    {no, prompt="Interactive mode for geomap"}
bool    verbose     {yes, prompt="Verbose"}
bool    fl_debug    {no, prompt="Print debugging information"}
int     status      {0, prompt="Exit status (0=good)"}

struct  *scanfile   {"", prompt = "Internal use only"}

begin

    char    l_indata = ""
    char    l_gprgid = ""
    char    l_instrument = ""
    char    l_outcoords = ""
    char    l_outtab = ""
    char    l_inimage = ""
    char    l_outimage = ""
    char    l_iraunits = ""
    char    l_fraunits = ""
    char    l_logfile = ""
    real    l_rafield, l_decfield, l_pa
    bool    l_fl_getxy, l_fl_getim
    bool    l_fl_inter, l_verbose, l_fl_debug

    # other variable declarations here
    char    paramstr, errmsg
    char    outpref = "GMI"             # applied regardless of whether outimage is set or not
    real    l_slitszx, l_slitszy
    int     junk
    
    junk = fscan (indata, l_indata)
    junk = fscan (gprgid, l_gprgid)
    junk = fscan (instrument, l_instrument)
    l_rafield = rafield
    l_decfield = decfield
    l_pa = pa
    l_fl_getxy = fl_getxy
    junk = fscan (outcoords, l_outcoords)
    junk = fscan (outtab, l_outtab)
    l_fl_getim = fl_getim
    junk = fscan (inimage, l_inimage)
    junk = fscan (outimage, l_outimage)
    junk = fscan (iraunits, l_iraunits)
    junk = fscan (fraunits, l_fraunits)
    l_slitszx = slitszx
    l_slitszy = slitszy
    junk = fscan (logfile, l_logfile)
    l_fl_inter=fl_inter
    l_verbose = verbose
    l_fl_debug = fl_debug

    cache ("gloginit", "glogpars") 

    # Initialize
    status = 0

    # Create the list of parameter/value pairs.  One pair per line.
    # All lines combined into one string.  Line delimiter is '\n'.
    paramstr =  "indata         = "//indata.p_value//"\n"
    paramstr += "gprgid         = "//gprgid.p_value//"\n"
    paramstr += "rafield        = "//rafield.p_value//"\n"
    paramstr += "decfield       = "//decfield.p_value//"\n"
    paramstr += "pa             = "//pa.p_value//"\n"
    paramstr += "fl_getxy       = "//fl_getxy.p_value//"\n"
    paramstr += "outcoords      = "//outcoords.p_value//"\n"
    paramstr += "outtab         = "//outtab.p_value//"\n"
    paramstr += "fl_getim       = "//fl_getim.p_value//"\n"
    paramstr += "inimage        = "//inimage.p_value//"\n"
    paramstr += "outimage       = "//outimage.p_value//"\n"
    paramstr += "instrument     = "//instrument.p_value//"\n"
    paramstr += "iraunits       = "//iraunits.p_value//"\n"
    paramstr += "fraunits       = "//fraunits.p_value//"\n"
    paramstr += "slitszx        = "//slitszx.p_value//"\n"
    paramstr += "slitszy        = "//slitszy.p_value//"\n"
    paramstr += "logfile        = "//logfile.p_value//"\n"
    paramstr += "fl_inter       = "//fl_inter.p_value//"\n"
    paramstr += "verbose        = "//verbose.p_value//"\n"
    paramstr += "fl_debug       = "//fl_debug.p_value

    # Open the log file here then carry it though to the subroutines
    # Assign a logfile name if not given.  Open logfile and start log.
    # Write parameter/value pairs ("paramstr") to log.
    gloginit (l_logfile, "gmskcreate", "mostools", paramstr, fl_append+,
        verbose=l_verbose)
    if (gloginit.status != 0) {
        status = gloginit.status
        goto exitnow
    }
    l_logfile = gloginit.logfile.p_value

    ############# check input parameters are set correctly ################
    
    # if indata is not set then exit
    if (l_indata == "") {
        errmsg = "The  input file, indata, containing id, ra, dec, mag of \
            input objects, must be set"
        status = 121
    }	

    # if indata not found then exit
    if (no == access (l_indata)) {
        errmsg = "RA/Dec coordinate file "//l_indata//" not found."
        status = 101
    }

    if (status != 0) {
        glogprint (l_logfile, "gmskcreate", "status", type="error",
            errno=status, str=errmsg, verbose+)
        goto exit
    }

    # if neither l_fl_getxy or l_fl_getim is set then there is nothing to 
    # do and exit
    
    if ((l_fl_getxy == no) && (l_fl_getim == no)) {
        errmsg = "Neither fl_getxy nor fl_getim is set to yes, so there is \
            nothing for the script to do"
        status = 121
        glogprint (l_logfile, "gmskcreate", "status", type="error",
            errno=status,  str=errmsg, verbose+)
        goto exit
    }


    # if they want the output xy and gmmps OT (ie l_fl_getxy is set) then 
    # if l_outcoords is not defined then set it to prefix//incoords
    # if l_outtab is defined check whether it has .fits on the end, and 
    #       add it if not
    # if l_outtab is not defined then set it to prefix//incoords//fits
    # if either outcoords or outtab already exist then exit
    
    if (l_fl_getxy) {

        if (l_outcoords == "") {
            fparse (l_indata, ver-)
            l_outcoords = outpref//fparse.root//fparse.extension
        }

        if (l_outtab != "") {
            fparse (l_outtab, ver-)
            if (fparse.extension != ".fits")
                l_outtab = l_outtab//".fits"
        } else {
            fparse (l_indata, ver-)
            l_outtab = outpref//fparse.root//"_OT.fits"
        }

        if (yes == access(l_outcoords)) {
            errmsg = "Output file "//l_outcoords//" already exists."
            status = 102
            glogprint (l_logfile, "gmskcreate", "status", type="error",
                errno=status, str=errmsg, verbose+)
        }

        if (yes == access(l_outtab)) {
            errmsg = "Output file "//l_outtab//" already exists."
            status = 102
            glogprint (l_logfile, "gmskcreate", "status", type="error",
                errno=status, str=errmsg, verbose+)
        }
    }

    # if they want the output fake image then 
    # the input image must be defined
    # the input image must exist
    # if the output image name is not defined then set it to prefix//inimage
    # if the output image name is not defined and the prefix is not defined 
    #       then stop
    # if the output image exists then stop

    if (l_fl_getim) {
        if (l_inimage == "") {
            errmsg =  "The inimage parameter must be set in order to produce a \
                fake "//l_instrument//" image"
            status = 121
            glogprint (l_logfile, "gmkcreate", "status", type="error",
                errno=status, str=errmsg, verbose+)

        } else if (no == imaccess(l_inimage)) {
            errmsg = "Input image "//l_inimage//" not found"
            status = 101
            glogprint (l_logfile, "gmkcreate", "status", type="error",
                errno=status, str=errmsg, verbose+)
        }	    

        if (l_outimage == "")
            l_outimage = outpref//l_inimage
        else
            l_outimage = outpref//l_outimage

        if (yes == imaccess(l_outimage)) {
            errmsg = "Output image "//l_outimage//" already exists."
            status = 102
            glogprint (l_logfile, "gmkcreate", "status", type="error", 
                errno=status, str=errmsg, verbose+)
        } else {
            glogprint (l_logfile, "gmskcreate", "task", type="string",
                str="Output image = "//l_outimage, verbose=l_verbose)
       }
    }

    if (status != 0) {
        goto exit
    }

    ################### end check of input parameters ########################

    # GMSKXY
    if (l_fl_getxy == yes) {
        glogprint (l_logfile, "gmskcreate", "status", type="fork", 
            fork="forward", child="gmskxy", verbose=l_verbose)
        gmskxy (l_indata, l_instrument, l_rafield, l_decfield, l_pa,
            iraunits=l_iraunits, fraunits=l_fraunits, outcoords=l_outcoords,
            slitszx=l_slitszx, slitszy=l_slitszy, outtab=l_outtab,
            logfile=l_logfile, verbose=l_verbose, fl_debug=l_fl_debug,
            status=status)
        glogprint (l_logfile, "GMSKCREATE", "status", type="fork",
            fork="backward", child="gmskxy", verbose=l_verbose)
    }

    # GMSKIMG
    if (l_fl_getim == yes) {
        glogprint (l_logfile, "gmskcreate", "status", type="fork",
            fork="forward", child="gmskimg", verbose=l_verbose)
        gmskimg (l_inimage, l_gprgid, l_instrument, rafield=l_rafield,
            decfield=l_decfield, pa=l_pa, fraunits= l_fraunits,
            outimage=l_outimage, fl_inter=l_fl_inter, logfile=l_logfile, 
            verbose=l_verbose, fl_debug=l_fl_debug, status=status)	
        glogprint (l_logfile, "gmskcreate", "status", type="fork",
            fork="backward", child="gmskimg", verbose=l_verbose)
    }


exit:
    scanfile = ""

    if (status == 0)
        glogclose (l_logfile, "gmskcreate", fl_success+, verbose=l_verbose)
    else
        glogclose (l_logfile, "gmskcreate", fl_success-, verbose=l_verbose)

exitnow:
    ;

end

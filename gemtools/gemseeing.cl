# Copyright(c) 1997-2006 Inger Jorgensen
# Copyright(c) 2000-2015 Association of Universities for Research in Astronomy, Inc.

procedure gemseeing(image)

# Semi-automatic method for determination of image quality for 
# imaging data. 
# Output: direct fwhm, moffat fwhm, eed50, eed85, strehl ratio
#
# Suggested parameter settings:
#  Seeing  0.4-0.8 arcsec: fwhmin=0.6, radius=0, fwhmmin=0, fwhmmax=0
#  Seeing    2-4   arcsec: fwhmin=3, radius=0, fwhmmin=0, fwhmmax=0
# Default values in database gemtools$gemseeing.dat reflects expected
# seeing with various instruments.
# 
# Known bugs/limitations:
#  strehl ratio output is not reliable in crowded fields
#  a MEF file with valid extensions may have NEXTEND undefinded in PHU
#  this causes gemseeing to crash - and is not a valid FITS header
#  though it works for most tasks.
#
# Version  Sept 1, 2000  IJ
#          Nov 8,  2001  IJ  v1.3 release
#          Sept 20,2002  IJ  v1.4 release
#          Oct 31, 2002  IJ  work on raw GMOS images, output in current dir
#          Mar 20, 2003  BM  include GMOS-S, database table also updated
#          Aug 28, 2003  KL  IRAF2.12 - new parameters
#                              hedit: addonly
#                              imstat: nclip,lsigma,usigma,cache
#                              apphot.phot: wcsin,wcsout,cache
#                              apphot.daofind: wcsout,cache
#                            IRAF2.12 - replaced task
#                              txdump:  been replaced by pdump (same parameters)
#                            Bug fix: fix a typo in imexamine call
#                            Bug fix: initialize l_fl_ext to 'no' before using it.
#                            Workaround: added '.fits' extension to fxhead calls
#                                  to workaround a bug in fxhead.

char    image       {prompt="Image (w/ MEF extension)"}                 # OLDP-1-input-primary-update
char    coords      {"default",prompt="daofind output file"}            # OLDP-3
char    psffile     {"default",prompt="File for psfmeasure output"}     # OLDP-3
bool    fl_inter    {no,prompt="Interactive use"}                       # OLDP-3
bool    fl_useall   {no,prompt="Use all objects if N>maxnobj"}          # OLDP-3
int     maxnobj     {100,min=5,prompt="Approximate maximum number of objects"}  # OLDP-3
bool    fl_update   {yes,prompt="Update image keywords"}                # OLDP-3
bool    fl_keep     {yes,prompt="Keep output tables"}                   # OLDP-3
bool    fl_strehl   {yes,prompt="Derive Strehl ratio if defined for instrument"}   # OLDP-3
bool    fl_overwrite    {no,prompt="Overwrite results from previous use"}   # OLDP-2
bool    fl_cleed    {yes,prompt="Clean measurements of EED of outliers"}    # OLDP-3
char    key_inst    {"INSTRUME",prompt="Keyword for instrument"}            # OLDP-3
char    key_camera  {"CAMERA",prompt="Keyword for NIRI camera"}             # OLDP-3
char    key_ron     {"RDNOISE",prompt="Keyword for read out noise in electrons"}   # OLDP-3
char    key_gain    {"GAIN",prompt="Keyword for gain in electrons/ADU"} # OLDP-3
char    key_filter  {"FILTER",prompt="Keyword for filter"}              # OLDP-3
char    instrument  {"",prompt="Instrument"}                            # OLDP-3
char    camera      {"",prompt="NIRI camera"}                           # OLDP-3
real    ron         {1.,min=0.,prompt="Read out noise in electrons"}    # OLDP-3
real    gain        {1.,min=0.00001,prompt="Gain in e-/ADU"}            # OLDP-3
char    filter      {"",prompt="Filter name"}                           # OLDP-3
bool    fl_database {yes,prompt="Use database values for parameters"}   # OLDP-3
char    database    {"gemtools$gemseeing.dat",prompt="Database to use"} # OLDP-2
char    stdatabase  {"gemtools$gemseeing_strehl.dat",prompt="Database for Strehl ratio"} # OLDP-3
real    pixscale    {0.5,min=0.00001,prompt="Pixel scale arcsec/pixel [database]"} # OLDP-2
real    fwhmin      {1.,min=0.001,prompt="Input FWHM PSF in arcsec (daofind) [database]"} # OLDP-1
real    datamax     {40000.,prompt="Maximum valid data value [database]"}   # OLDP-2
real    radius      {0.,min=0.,prompt="Radius in pixels for fit [database]"}    # OLDP-3
real    buffer      {8.,min=0.,prompt="Buffer zone in pixels for sky [database]"}  # OLDP-3
real    width       {2.,min=1.,prompt="Width in pixels for sky [database]"} # OLDP-3
real    fwhmnonao   {0.5,min=0.05,prompt="App non-AO FWHM in arcsec if using AO"}  # OLDP-3
real    sigthres    {10.,min=0.,prompt="Threshold in sigma of sky (daofind)"}   # OLDP-2
real    threshold   {INDEF,min=0.,prompt="Threshold for daofind (daofind)"} # OLDP-2
int     niter       {3,min=1,prompt="Number of iterations (psfmeasure)"}    # OLDP-3
real    mdaofaint   {-2.0,max=0.,prompt="Faintest magnitude from daofind"}  # OLDP-3
real    fwhmmin     {0.,min=0.,prompt="Minimim fwhm in arcsec"}             # OLDP-3
real    fwhmmax     {0.,min=0.,prompt="Maximum fwhm in arcsec"}             # OLDP-3
real    sharplow    {0.3,prompt="Minimum sharp from daofind"}               # OLDP-3
real    sharphigh   {1.0,prompt="Maximum sharp from daofind"}               # OLDP-3
real    sharpsig    {3.,min=0.,prompt="Sigma selection on sharp from daofind"}     # OLDP-3
real    dmagsig     {2.,min=0.,prompt="Sigma selection on magnitude difference"}   # OLDP-3
real    fwhmsig     {2.5,min=0.,prompt="Sigma selection on FWHM"}       # OLDP-3
char    logfile     {"gemseeing.log",prompt="Logfile"}                  # OLDP-1
bool    verbose     {yes,prompt="Print intermediate results"}           # OLDP-4
int     status      {0,prompt="Exit status (0=good)"}                   # OLDP-4

begin

    char    l_image, l_coords, l_psffile, l_database, l_stdatabase
    char    l_key_ron, l_key_gain, l_key_inst, l_inst, l_camera, l_key_camera
    char    l_key_filter, l_filter
    char    l_logfile
    real    l_pixscale, l_fwhmin, l_ron, l_gain, l_threshold, l_sigthres
    real    l_fwhmmin, l_fwhmmax, l_radius, l_buffer, l_width, l_fwhmsig
    real    l_mdaofaint, l_datamax, l_sharpsig, l_dmagsig
    real    l_sharplow, l_sharphigh
    real    l_fwhmnonao
    bool    l_fl_inter, l_fl_overwrite, l_fl_cleed
    bool    l_fl_update, l_fl_useall, l_fl_keep, l_fl_database, l_fl_strehl

    int     l_niter, l_maxnobj, ii, next_inp
    real    out_npsf, out_fwhm, out_mfwhm, out_eed50, out_eed85, out_eps
    real    out_pos, out_strehl
    real    rms_npsf, rms_fwhm, rms_mfwhm, rms_eed50, rms_eed85, rms_eps
    real    rms_pos, rms_strehl, x_psf, y_psf

    char    tmpcoo, tmpdao, tmpime, tmpcd, tmpdat, tmpfit, tmpmom, tmpgra
    char    tmpsubset, junk
    char    l_expression, s_image, l_section, junk1, junk2, t_image
    char    l_sciext, l_save, l_dettype
    char    l_imagezero
    real    l_level, l_skysigma, l_npeak, l_wavelength, l_bin
    int     l_Xmax, l_Ymax, l_obj, n_i, l_ext, l_col, l_nsec
    bool    l_verbose, l_fl_ext, l_fl_disp, l_fl_one
    struct  l_struct

    l_verbose = verbose
    l_image=image ; l_coords=coords ; l_psffile=psffile
    
    cache ("tinfo")
    
    # Make image name without section if necessary, s_image is the name
    # without the section.  Cannot check access with one statement
    # if image section is present.
    # l_section contains MEF extension if present and section
    
    s_image=l_image ; t_image=l_image ; l_section=""
    if (stridx("[",l_image) != 0) { 
        s_image = substr(l_image, 1, (stridx("[",l_image)-1) )
        l_section = substr(l_image, stridx("[",l_image), strlen(l_image))
        t_image = s_image
    }

    #Initialize
    # KL : in some scenarios, l_fl_ext was used before it had been initialized
    #      resulting in a crash.  Same thing happens on IRAF2.11.3b.
    
    l_fl_ext = no

    status = 0
    # set logfile
    l_logfile = logfile
    if (l_logfile=="" || l_logfile==" ") {
        print ("WARNING - GEMSEEING: No logfile defined. Output to \
            screen only.")
        l_logfile = ""
    }

    date | scan(l_struct)
    printlog ("--------------------------------------------------------------\
        --------------", l_logfile, l_verbose)
    printlog ("GEMSEEING -- "//l_struct, l_logfile, l_verbose)
    printlog ("", l_logfile, l_verbose)

    # check access
    if (no==imaccess(s_image) && no==imaccess(s_image//"[0]") && no==imaccess(s_image//"[1]")) {
        printlog ("ERROR - GEMSEEING: Image "//l_image//" not found",
            l_logfile, yes)
        status = 1
        bye
    }
    # cut .imh off s_image if present 
    if (substr(s_image,strlen(s_image)-3,strlen(s_image)) == ".imh") 
        s_image = substr(s_image, 1, (strlen(s_image)-4))
        
    # cut .fits off s_image if present
    if (substr(s_image, strlen(s_image)-4, strlen(s_image)) == ".fits") 
        s_image = substr(s_image, 1, (strlen(s_image)-5))

    l_imagezero = s_image  # if simple fits or imh

    l_fl_update=fl_update ; l_fl_useall=fl_useall ; l_fl_keep=fl_keep
    l_fl_inter=fl_inter ; l_maxnobj=maxnobj

    cache ("imgets", "fparse", "gemdate")
    fparse (s_image)
    if (l_coords == "default")
        l_coords = fparse.root//"_coo"
    if (l_psffile == "default")
        l_psffile = fparse.root//"_see"

    # Check of whether this is a MEF file - if it is not a fits file
    # the check should not be made
    # Make image name with the MEF extension name if this is a MEF file
    # Use extension defined in input, otherwise use SCI extension if 
    # present, otherwise extension number 1
    # This can still crash in NEXTEND is not defined in header of MEF[0]
    
    if (access(s_image//".fits")) {
        imgets (s_image//"[0]", "EXTEND", >& "dev$null")
        if (imgets.value != "F") {
            imgets (s_image//"[0]", "NEXTEND", >& "dev$null")
            if (imgets.value != "0") {
                l_fl_ext = yes
                l_ext = 0 
            } else {
                # Most peculiar bug?, need the extra fxhead to make sure the 
                # lines get counted correctly???, IJ
                # KL: Also, the '.fits' extension had to be added as fxhead 
                #     sometimes won't work properly if the full name is not 
                #     given.  Strange, I know.
                
                fxhead (s_image//".fits", format_file="", long_head-,
                    count_lines-, >& "dev$null")
                fxhead (s_image//".fits", format_file="", long_head-,
                    count_lines-) | count ("STDIN") | scan (next_inp)
                fxhead (s_image//".fits", format_file="", long_head-,
                    count_lines-) | fields ("STDIN", "1", lines=str(next_inp),
                    quit_if_miss-, print_file-) | scan (next_inp)
                if (next_inp > 0) {
                    l_fl_ext=yes
                    l_ext=0
                }
            }
            if (l_fl_ext) {
                # If there is an extension specified then l_section!=""
                # l_section!="" may also be due to a section specified
                
                if (l_section != "") {
                    print (l_section) | match (",", stop-) | count | \
                        scan (l_nsec)
                    if (l_nsec != 0) {  # section defined, check for extension
                        print (l_section) | match ("\]\[", stop-) | count | \
                            scan (l_ext)
                        if (l_ext!=0) {  # section and extension defined
                            l_sciext = substr (l_section, stridx("[",l_section),
                                stridx("]",l_section))
                            l_section = substr (l_section, strlen(l_sciext)+1,
                                strlen(l_section))
                            print (l_section) | match (",", stop-) | count | \
                                scan (l_nsec)
                            if (l_nsec == 0) {
                                # just switch the variables if needed, rather 
                                # than exit
                                l_save = l_section
                                l_section = l_sciext
                                l_sciext = l_save
                            }
                        } else {  # only section defined, set default extension
                            if (imaccess(s_image//"[SCI]"))
                                l_sciext = "[SCI]" 
                            else
                                l_sciext = "[1]" 
                        }
                    } else {
                        l_ext = 1  # extension defined, no section
                        l_sciext = l_section
                        l_section = ""
                    }

                    l_imagezero = s_image//"[0]"
                    l_image = s_image//".fits"//l_sciext//l_section
                    t_image = s_image//".fits"//l_sciext
                    s_image = s_image//l_sciext//l_section
                } else {

                    if (imaccess(s_image//"[SCI]")) {
                        l_imagezero = s_image//"[0]"
                        l_image = s_image//".fits[SCI]"//l_section
                        t_image = s_image//".fits[SCI]"
                        s_image = s_image//"[SCI]"//l_section
                    } else {
                        l_imagezero = s_image//"[0]"
                        l_image = s_image//".fits[1]"//l_section
                        t_image = s_image//".fits[1]"
                        s_image = s_image//"[1]"//l_section
                    }
                }
            }
        }
    } else
        s_image = s_image//l_section

    # slightly redundant but catches all
    if (no == imaccess(l_image)) {
        printlog ("ERROR - GEMSEEING: Image "//l_image//" not found",
            l_logfile, yes)
        status = 1
        bye
    }

    l_fl_cleed = fl_cleed
    l_fl_overwrite = fl_overwrite 
    if (l_fl_overwrite) 
        delete (l_coords//","//l_psffile//","//l_psffile//".fits,\
            "//l_psffile//"_psf.fits", verify-, >& "dev$null")

    # check if final psffile exists
    if (access(l_psffile//"_psf.fits")) {
        printlog ("ERROR - GEMSEEING: PSF file "//l_psffile//"_psf.fits exists",
            l_logfile, yes)
        status = 1
        bye
    }

    # Get instrument (and camera), ron and gain from header if keywords are defined
    l_key_inst = key_inst
    if (l_key_inst!="" && l_key_inst!=" ") {
        imgets (l_imagezero, l_key_inst, >& "dev$null") 
        print (imgets.value) | scan (l_inst)
        if (l_inst == "0") 
            l_inst = instrument
    } else
        l_inst = instrument 

    if (l_inst=="NIRI") {
        l_key_camera = key_camera
        if (l_key_camera != "" && l_key_camera != " ") {
            imgets (l_imagezero, l_key_camera, >& "dev$null")
            l_camera = imgets.value
            if (l_camera == "0") 
                l_camera = camera
        } else
            l_camera = camera 

        l_inst = l_inst//"-"//l_camera
    }
    
    # disable database use if no valid instrument
    l_fl_database = fl_database
    if (l_inst == "" || l_inst == " ")
        l_fl_database = no

    l_key_ron = key_ron
    l_key_gain = key_gain
    if (l_key_ron != "" && l_key_ron != " ") {
        imgets (l_imagezero, l_key_ron, >& "dev$null")
        l_ron = real(imgets.value)
        if (l_ron == 0) 
            l_ron = ron
    } else
        l_ron = ron 

    if (l_key_gain != "" && l_key_gain != " ") {
        imgets (l_imagezero, l_key_gain, >& "dev$null")
        l_gain = real(imgets.value)
        if (l_gain == 0) 
            l_gain = gain
    } else
        l_gain = gain 

    # Use database if present
    l_database = database
    if (l_database=="" || l_database==" ") 
        l_fl_database = no
    if (no == access(l_database)) {
        l_fl_database = no
        printlog ("WARNING - GEMSEEING: Cannot access database "//l_database,
            l_logfile, yes)
        printlog ("                     Using parameter file values",
            l_logfile, yes)
    }

    l_pixscale = 0
    if (l_fl_database) {
        if ((l_inst == "GMOS-S") || (l_inst == "GMOS-N")) {
            keypar(l_imagezero, "DETTYPE", silent+)
            l_dettype = keypar.value
            match (l_dettype, files=database) | \
                fields ("STDIN", "1-7", lines="1-", quit_if_miss-, \
                    print_file_n-) | \
                    match (l_inst, stop-) | \
                        scan (l_inst, l_pixscale, l_datamax, l_fwhmin, \
                            l_radius, l_buffer, l_width)
        } else {
            fields (database, "1-7", lines="1-", quit_if_miss-, \
                print_file_n-) | \
                    match (l_inst, stop-) | \
                        scan (l_inst, l_pixscale, l_datamax, l_fwhmin, \
                            l_radius, l_buffer, l_width)
        }
        if (l_pixscale == 0) {
            printlog ("WARNING - GEMSEEING: Instrument "//l_inst//" not in \
                database", l_logfile, yes)
            printlog ("          Using parameter file values", l_logfile, yes)
        }
    } 
    
    if (l_pixscale == 0) {  #either fl_database was no, or instrument not in db
        l_pixscale = pixscale
        l_fwhmin = fwhmin
        l_radius = radius 
        l_buffer = buffer
        l_width = width
        l_datamax = datamax
    }

    # GMOS - adjust for binning
    if (l_inst=="GMOS-N" || l_inst=="GMOS-S") {
        imgets (l_image, "CCDSUM")
        print (imgets.value) | scan (l_bin)
        l_pixscale = l_pixscale * l_bin
    }

    # Set fwhmin to pixels rather than arcsec
    l_fwhmin = l_fwhmin / l_pixscale
    if (l_radius == 0)
        l_radius = l_fwhmin * 2.5  # changed from 2 to 3*1.1
        
    l_fwhmmin = fwhmmin / l_pixscale
    l_fwhmmax = fwhmmax / l_pixscale
    if (l_fwhmmin == 0 && l_fwhmmax == 0) {
        l_fwhmmin = l_fwhmin / 4.
        l_fwhmmax = l_fwhmin * 4.
    }
    if (l_fwhmmax == 0 || l_fwhmmax <= l_fwhmin) 
        l_fwhmmax = l_fwhmin * 4.

    l_fwhmnonao = fwhmnonao

    # Strehl ratio - check validity in database
    l_fl_strehl = fl_strehl
    if (l_fl_strehl) {
        l_stdatabase = stdatabase
        if (l_stdatabase=="" || l_stdatabase==" " || !access(l_stdatabase)) {
            print ("WARNING - GEMSEEING: Strehl ratio database not found, \
                fl_strehl=no", l_logfile, yes)
            l_fl_strehl = no
        }
    }

    if (l_fl_strehl) {
        #  l_stdatabase=stdatabase
        l_key_filter = key_filter
        if (l_key_filter != "" && l_key_filter != " ") {
            imgets (l_imagezero, l_key_filter, >& "dev$null")
            l_filter = imgets.value
            if (l_filter == "0") 
                l_filter = filter
        } else
            l_filter = filter
            
        if (l_filter == "" || l_filter == " " || l_filter == "0") {
            l_fl_strehl = no
            goto skipstrehl
        }
        # attach spaces to filter, so a simple match will find it
        l_filter = " "//l_filter//" "

        l_npeak = 0 ; l_wavelength = 0 ; junk = ""
        fields (l_stdatabase, "1-5", lines="1-", quit_if_miss-,
            print_file_n-) | match (l_inst, stop-) | match (l_filter, stop-) | \
            scan (junk, junk, l_wavelength, l_filter, l_npeak)
        if (l_npeak == 0)  {
            l_fl_strehl = no
            printlog ("WARNING - GEMSEEING: Filter/instrument combination "//\
                l_inst//" "//l_filter//" not in Strehl database",
                l_logfile, yes)
            printlog ("                     No calculation of Strehl ratio",
                l_logfile, yes)
            goto skipstrehl
        }
    }

skipstrehl:

    l_niter = niter 
    l_sharplow = sharplow ; l_sharphigh = sharphigh
    l_mdaofaint = mdaofaint ; l_sharpsig = sharpsig
    l_dmagsig = dmagsig ; l_fwhmsig = fwhmsig

    tmpcoo = mktemp("tmpcoo")
    tmpdao = mktemp("tmpdao")
    tmpsubset = mktemp("tmpsubset")
    tmpdat = mktemp("tmpdat")
    tmpime = mktemp("tmpime")
    tmpcd = mktemp("tmpcd")
    tmpfit = mktemp("tmpfit")
    tmpmom = mktemp("tmpmom")
    tmpgra = mktemp("tmpgra")

    # ---------------------------------------------------------------------
    # Find the objects using daofind or interact, if l_coords does not exist already
    # ---------------------------------------------------------------------

    l_fl_disp = no
    if (no == access(l_coords) ) {

        if (l_fl_inter) {
        
            # ----- Interactive finding -----
            
            print ("")
            print ("Use the image cursor to point to objects in the display")
            print ("   x  - select object")
            print ("   q  - quit, when done")
            print ("")
            # need to display the image explicitly to make the tvmark work
            display (l_image, 1, >& "dev$null")
            
            # KL - Fix a typo in imexamine call: 
            #    was ...(image'$1'...   ;  should be ...(image='$1'...
            
            imexamine (l_image, frame=1, logfile=l_coords, keeplog+,
                display="display(image='$1',frame=1)")
            l_fl_disp = yes
        } else {
        
            # ----- DAOFIND -----
            
            # set threshold  avoid extreme values
            imstat (l_image, fields="midpt", lower=-l_datamax, upper=l_datamax,
                nclip=0, lsigma=INDEF, usigma=INDEF, binwidth=0.1, format=no,
                cache=no) | scan(l_level)
            imstat (l_image, fields="midpt",
                lower=(l_level-5.*sqrt( (l_ron/l_gain)**2+l_level/l_gain )),
                upper=(l_level+5.*sqrt( (l_ron/l_gain)**2+l_level/l_gain )),
                nclip=0, lsigma=INDEF, usigma=INDEF, binwidth=0.1, format-,
                cache-) | scan (l_level)

            if (sigthres == INDEF) {
                l_threshold = threshold
                l_sigthres = l_threshold / \
                    sqrt( (l_ron/l_gain)**2 + l_level/l_gain )
            } else {
                l_sigthres = sigthres
            }
            n_i = strlen(s_image)
            printf ("Image  %-"//str(n_i)//"s   Sigma threshold %6.1f\n",
                s_image, l_sigthres) | scan (l_struct)
            printlog (l_struct, l_logfile, l_verbose)
            printf ("       %-"//str(n_i)//"s   Sky level [ADU] %6.1f\n"," ",
                l_level) | scan(l_struct)
            printlog (l_struct, l_logfile, l_verbose)

            # find objects try to get rid of cosmics here
            daofind (l_image, output=l_coords, starmap="", skymap="",
                datapars="", findpars="", boundary="nearest", constant=0.,
                interactive=no, icommands="", gcommands="", wcsout="logical",
                cache=no, verify=no, update=no,
                datapars.fwhmpsf=l_fwhmin, datapars.readnoise=l_ron,
                datapars.epadu=l_gain, datapars.ccdread="", datapars.gain="",
                datapars.noise="poisson", datapars.exposure="",
                datapars.airmass="", datapars.filter="", datapars.obstime="",
                datapars.itime=1,
                datapars.sigma=(sqrt( (l_ron/l_gain)**2+l_level/l_gain )),
                datapars.datamax=l_datamax, datapars.datamin=-l_datamax,
                datapars.scale=1., datapars.emission=yes,
                findpars.nsigma=1.5, findpars.threshold=l_sigthres,
                findpars.ratio=1., findpars.theta=0., findpars.sharplo=0.2,
                findpars.sharphi=1.0, findpars.roundlo=-1., findpars.roundhi=1.,
                findpars.mkdetections-)
        } # end of finding

    } else {
        printf ("Image  %s\n", s_image) | scan (l_struct)
        printlog (l_struct, l_logfile, l_verbose)
        print ("Coordinate file "//l_coords//" exists, daofind not run") | \
            scan (l_struct)
        printlog (l_struct, l_logfile, l_verbose)
    } # end of accessing coordinate file

    # get image size in any case
    imgets (l_image, "i_naxis1")
    l_Xmax = int(imgets.value)
    imgets (l_image, "i_naxis2")
    l_Ymax = int(imgets.value)

    # Do this check in any case
    # check if coordinate file is from interactive or daofind
    head (l_coords, nlines=1) | scan (junk1, junk2)
    if (junk2 != "IRAF" && !l_fl_inter) {
        l_fl_inter = yes
        printlog ("WARNING - GEMSEEING: Existing coordinate file is not from \
            daofind", l_logfile, yes)
    }
    if (junk2 == "IRAF" && l_fl_inter) {
        l_fl_inter = no
        printlog ("WARNING - GEMSEEING: Existing coordinate file is from \
            daofind", l_logfile, yes)
    }

    if (no == l_fl_inter)  # only if daofind has been run
        pdump (l_coords, "threshold", yes) | scan (l_sigthres)
    else
        l_sigthres = 10000.

    # ---------------------------------------------------------------------
    # Derive direct fwhm using psfmeasure, if l_psffile does not exist
    # ---------------------------------------------------------------------
    if (no == access(l_psffile)) {

        if (l_fl_inter) {
        
            # ---- Interactive -----
            
            fields (l_coords, "1,2", lines="1-", quit-, print-, > tmpcoo)
            count (tmpcoo) | scan (l_obj)
            if (l_obj == 0) {
                printlog ("ERROR - GEMSEEING: No objects in coordinate \
                    file "//l_coords, l_logfile, yes)
                delete (tmpcoo//","//l_coords, verify-, >& "dev$null")
                status = 1
                bye
            }
            printf ("c1\nc2\n", > tmpcd)
            tcreate (tmpcoo//".fits", tmpcd, tmpcoo, uparfile="", nskip=0, nlines=0,
                nrows=0, hist=yes, extrapar=5, tbltype="default", extracol=0)
            delete (tmpcd//","//tmpcoo, verify=no, >& "dev$null")
        } else {
        
            # ---- Non-interactive -----
            
            pdump (l_coords, "xcenter,ycenter,mag",
                "mag<-0.05 && sharpness<999", headers-, parameters+, > tmpdao)

            # Check number of objects
            count (tmpdao) | scan (l_obj)
            if (l_obj == 0) {
                printlog ("ERROR - GEMSEEING: No valid objects found in \
                    image "//l_image, l_logfile, yes)
                delete (tmpdao//","//l_coords, verify-)
                status = 1
                bye
            }

            # get rid of objects too close to the edge of the image
            tselect (tmpdao, tmpcoo//".fits",
                "c1>5 && c1<("//l_Xmax//"-5) && c2>5 && c2<("//l_Ymax//"-5)")
            delete (tmpdao, verify=no, >& "dev$null")

            # Count objects
            tinfo (tmpcoo//".fits", ttout=no)
            l_obj =  tinfo.nrows
            if (l_obj==0) {
                printlog ("ERROR - GEMSEEING: No valid objects found in \
                    image "//l_image, l_logfile, yes)
                delete (tmpdao//","//tmpcoo//".fits", verify-, >& "dev$null")
                status = 1
                bye
            }

            n_i = strlen(s_image)
            printf ("Image  %-"//str(n_i)//"s   Total Nobj      %6d\n",
                s_image, l_obj) | scan (l_struct)
            printlog (l_struct, l_logfile, l_verbose)

            # If fl_useall==no:  Select subset 
            if (no == l_fl_useall) {
                # If there are more than l_maxnobj objects cut 0.5 mag fainter 
                # than mdaofaint unless that is larger than zero
                if (l_obj > l_maxnobj && l_mdaofaint < -0.5) {
                    tselect (tmpcoo//".fits", tmpsubset//".fits", \
                        "c3 < "//str(l_mdaofaint+0.5))
                    delete (tmpcoo//".fits", verify=no)
                    rename (tmpsubset//".fits", tmpcoo//".fits", field="all")
                }
                # If there are still more than l_maxnobj objects select subset
                tinfo (tmpcoo//".fits", ttout=no)
                l_obj = tinfo.nrows
                if (l_obj>l_maxnobj ) {
                    l_obj = max(2, int(l_obj/real(l_maxnobj)+0.5))
                    tselect (tmpcoo//".fits", tmpsubset//".fits", \
                        "mod(row(),"//str(l_obj)//")==0")
                    delete (tmpcoo//".fits", verify=no)
                    rename (tmpsubset//".fits", tmpcoo//".fits", field="all")
                }

            } 

        } # end of non-interactive

        tinfo (tmpcoo//".fits", ttout=no)
        l_obj = tinfo.nrows
        n_i = strlen(s_image)
        printf ("       %-"//str(n_i)//"s   Fitting Nobj    %6d\n"," ",
            l_obj) | scan (l_struct)
        printlog (l_struct, l_logfile, l_verbose)
        printf ("       %-"//str(n_i)//"s   Radius [pixels] %6.2f\n"," ",
            l_radius) | scan (l_struct)
        printlog (l_struct, l_logfile, l_verbose)

        printf (":show "//l_psffile//"\n:size MFWHM\n:show "//l_psffile//\
            "\nq\n", > tmpgra)

        # if l_obj=1 need to make fake line in input file
        # this in fact works ...
        #  KL Oct 2012 - I'm commenting this out, this is silly.
        #                Hopefully this isn't necessary anymore now that we
        #                are working with a FITS tmpcoo.
        #  UPDATE: KL Nov 2012 - Actually we can't get away without it!
        #                       psfmeasure requires a second line with a number in it.
        #                       Even an empty line won't do.  It's got to be a number.
        #                       Ridiculous.   I'm taking care of this below, just before
        #                       the call to noao.osbutil.psfmeasure.
        #if (l_obj == 1) {
        #    type (tmpcoo, >> tmpcoo)
        #}

        # Use pixelscale 1, convert later ...
        # psfmeasure silently throws away any objects it cannot handle, without
        # putting an entry in the output file ....
        #
        # For debugging :
        # psfmeasure(l_image,coords="markall",wcs="logical",display-,frame=1,
        #   level=0.5,size="FWHM",scale=1.,radius=l_radius,sbuffer=l_buffer,
        #   swidth=l_width,saturation=l_datamax,ignore_sat+,iterations=l_niter,
        #   logfile=l_psffile,imagecur=tmpcoo,graphcur="")
        # bye
        # psfmeasure crashes if ignore_sat+ and only saturated input objects
        
        tdump (tmpcoo//".fits", cdfile="", pfile="", datafile=tmpcoo//".dat",
            columns="", rows="-", pwidth=-1)
        # KL Nov 2012:  Bug in noao.obsutil.psfmeasure; it barfs when there's only
        #               one set of coordinates.  Dumb.  So we add a line with a zero
        #               in it.  (an empty line won't work, double dumb)
        count (tmpcoo//".dat") | scan(l_obj)
        if (l_obj == 1) {
            printf ("0\n", >> tmpcoo//".dat")
        }
        noao.obsutil.psfmeasure (l_image, coords="markall", wcs="logical", 
            display-, frame=1, level=0.5, size="FWHM", beta=INDEF, scale=1.,
            radius=l_radius, sbuffer=l_buffer, swidth=l_width,
            saturation=l_datamax, ignore_sat-, iterations=l_niter, logfile="",
            xcenter=INDEF, ycenter=INDEF, imagecur=tmpcoo//".dat", graphcur=tmpgra, 
            >& "dev$null", >>G "dev$null")
        delete (tmpgra//","//tmpcoo//".fits,"//tmpcoo//".dat", verify-, >& "dev$null")
    } else {
        printlog ("PSF file "//l_psffile//" exists, fitting not run",
            l_logfile, l_verbose)
    }

    # ---------------------------------------------------------------------
    # Make psffile//".fits" if it does not already exist 
    # ---------------------------------------------------------------------

    if (no == access(l_psffile//".fits")) {

        if (l_fl_inter) {
            printf ("x_dao r f6.2\ny_dao r f6.2\n") | \
                tcreate (tmpdao//".fits", "STDIN", l_coords, upar="", nskip=0, 
                    nlines=0, nrows=0, hist-, extrapar=5, tbltype="default", 
                    extracol=0)
        } else {
            # make FITS table with output from psfmeasure and daofind
            pconvert (l_coords, tmpdao//".fits", "*", expr="yes", append-)
            tchcol (tmpdao//".fits", "xcenter", "x_dao", "", "", verbose-)
            tchcol (tmpdao//".fits", "ycenter", "y_dao", "", "", verbose-)
            tchcol (tmpdao//".fits", "mag", "m_dao", "", "", verbose-)
        }

        # make the psffile FITS table from psfmeasure output

        printf ("mfwhm r f7.3\nbeta r f6.2\n", > tmpcd)
        printf ("x_psf r f7.2\ny_psf r f7.2\nm_psf r f7.3\nfwhm r f7.3\n\
            eps r f7.3\npos r f7.3\n", >> tmpcd)
        print ("sat ch*1 %s", >> tmpcd)

        # Find the length of psffile; has nothing to do with number of 
        # input objects
        count (l_psffile) | scan (l_obj)
        l_obj = l_obj/2

        if (l_obj == 0) {
            printlog ("ERROR - GEMSEEING: No valid objects fitted in image "//\
                l_image, l_logfile, yes)
            delete (tmpdao//".fits,"//tmpcd, verify-, >& "dev$null")
            status = 1
            bye
        }

        # For each get rows:  l_obj/2 excl. first 3 and last 2
        # There are saturated objects in this because psfmeasure crashes
        # if only saturated input objects and run with ignore_sat+
        for (n_i=1; n_i<=2; n_i+=1) {
            if (n_i == 1) {
                fields (l_psffile, "2-8", lines=str((n_i-1)*l_obj+4), > tmpdat)
                fields (l_psffile, "1-7",
                    lines=str((n_i-1)*l_obj+5)//"-"//str((n_i*l_obj-2)),
                    >> tmpdat)
            } else {
                fields (l_psffile, "5,6", lines=str((n_i-1)*l_obj+4), > tmpgra)
                fields (l_psffile, "4,5",
                    lines=str((n_i-1)*l_obj+5)//"-"//str((n_i*l_obj-2)),
                    >> tmpgra)
                joinlines (tmpgra, tmpdat, output=tmpmom, delim=" ",
                    missing="Missing", maxchars=161, shortest+, verbose+)
                delete (tmpdat//","//tmpgra, verify-, >& "dev$null")
                rename (tmpmom, tmpdat, field="all")
            }
        }

        # need to change how saturated objects are marked since tselect cannot
        # handle this silly output from psfmeasure

        sed("-e 's/$/N/'", tmpdat," | sed -e 's/\*   N/Y/' ", >> tmpgra)

        tcreate (tmpime//".fits", tmpcd, tmpgra, upar="", nskip=0, nlines=1, nrows=0,
            hist-, extrapar=5, tbltype="default", extracol=0)
        delete (tmpgra//","//tmpdat//","//tmpcd, verify=no, >& "dev$null")

        # Get rid of objects that are not unique. There such objects because
        # daofind finds lots of spurious objects in AO images.
        # This does not take care of all of them, but a large part of them

        tproject (tmpime//".fits", tmpgra//".fits", "*", uniq+)

        # match the two tables
        if (l_fl_inter) {
            tmatch (tmpgra//".fits", tmpdao//".fits", l_psffile//".fits", 
                "x_psf,y_psf", "x_dao,y_dao", maxnorm=(l_fwhmin*3.0), incol1="", 
                incol2="", factor="", diagfile="", nmcol1="", nmcol2="", sphere-)
            tinfo (l_psffile//".fits", ttout=no)
            if (tinfo.nrows == 0) {
                printlog ("ERROR - GEMSEEING: No usable objects.",
                    l_logfile, yes)
                status = 1
                bye
            }
        } else {
            tmatch (tmpgra//".fits", tmpdao//".fits", l_psffile//".fits", 
                "x_psf,y_psf", "x_dao,y_dao", maxnorm=(l_fwhmin*0.5), incol1="", 
                incol2="", factor="", diagfile="", nmcol1="", nmcol2="", sphere-)
            tinfo (l_psffile//".fits", ttout=no)
            if (tinfo.nrows == 0) {
                printlog ("ERROR - GEMSEEING: No usable objects.",
                    l_logfile, yes)
                status = 1
                bye
            }
            # calculate magnitude difference
            tcalc (l_psffile//".fits", "dmag", "m_dao-m_psf", colfmt="f6.3",
                colunits="")
        }
        tdelete (tmpdao//".fits,"//tmpime//".fits,"//tmpgra//".fits", verify=no)

    } else {
        printlog ("PSF STSDAS table "//l_psffile//".fits exists.",
            l_logfile, l_verbose)
    }

    # -----------------------------------------------------------------
    # Select objecs from the table to be used for seeing determination
    # If there after the first selection is less than 5 objects left
    # the program goes directly to the final run of psfmeasure and the 
    # statistics. If there after any of the following selections is less
    # than 3 objects left, that selection is omitted and the final run
    # of psfmeasure and statistics is done on the sample before that last
    # selection.
    # -----------------------------------------------------------------
    
    cache ("tselect", "tstat", "tinfo")

    # check if l_psffile.fits has a dmag column 
    if (no == l_fl_inter) {
        l_col = 0
        tlcol (l_psffile//".fits", nlist=1) | match ("dmag", "STDIN") | \
            count | scan (l_col)
        if (l_col == 0) {
            print ("WARNING: Existing PSF file is from interactive use.")
            print ("         Selection turned off")
            l_fl_inter = yes
        }
    }

    # use all if interactive or if dmag column was missing
    if (l_fl_inter) {
        tcalc (l_psffile//".fits", "flag", "1", colfmt="i2", colunit="")
    } else {

        # select on saturated and magnitude
        tselect (l_psffile//".fits", tmpdao//".fits",
            "sat!='Y' && m_dao<"//l_mdaofaint//" && fwhm>="//l_fwhmmin//" && fwhm<"//l_fwhmmax )
        tinfo (tmpdao//".fits",ttout=no)
        if (tinfo.nrows < 3) {
            printlog ("WARNING - GEMSEEING: Too few objects, only saturated \
                objects and objects", l_logfile, yes)
            printf ("                     with fwhm < %7.4f pixels are \
                excluded\n", l_fwhmmin) | scan (l_struct)
            printlog (l_struct, l_logfile, yes)
            tdelete (tmpdao//".fits", verify=no)
            tselect (l_psffile//".fits", tmpdao//".fits", 
                "sat!='Y' && fwhm>="//l_fwhmmin)
            tinfo (tmpdao//".fits", ttout=no)
            if (tinfo.nrows < 1) {
                printlog ("ERROR - GEMSEEING: All objects saturated",
                    l_logfile, yes)
                tdelete (tmpdao//".fits", verify-)
                status = 1
                bye
            }
            goto finalselection
        }
        if (tinfo.nrows <= 5)
            goto finalselection

        tstat (tmpdao//".fits", "dmag", outtable="", lowlim=INDEF, highlim=INDEF,
            rows="-") 
        tselect (tmpdao//".fits", tmpime//".fits", 
            "m_psf<("//l_mdaofaint//"-"//tstat.median//")")
        tdelete (tmpdao//".fits", verify=no)

        tinfo (tmpime//".fits", ttout=no)
        if (tinfo.nrows < 3) {
            tdelete (tmpime//".fits", verify=no)
            goto finalselection
        }

        # select on sigma deviation in sharp
        tstat (tmpime//".fits", "sharpness", outtable="", lowlim=l_sharplow,
            highlim=l_sharphigh, rows="-") 
        if (tstat.stddev > 0) {
            l_expression = "sharpness>=("//tstat.median//"-"//l_sharpsig//"*"//tstat.stddev//") && sharpness<=("//tstat.median//"+"//l_sharpsig//"*"//tstat.stddev//")"
            tselect (tmpime//".fits", tmpdao//".fits", l_expression)
            tdelete (tmpime//".fits", verify=no)
        } else {
            printlog ("WARNING - GEMSEEING: No selection on sharpness",
                l_logfile, yes)
            rename (tmpime//".fits", tmpdao//".fits", field="all")
        }

        tinfo (tmpdao//".fits", ttout=no)
        if (tinfo.nrows < 3) {
            tdelete (tmpdao//".fits", verify=no)
            rename (tmpime//".fits", tmpdao//".fits", field="all")
            goto finalselection
        }

        # select on sigma deviation in dmag
        tstat (tmpdao//".fits", "dmag", outtable="", lowlim=INDEF, highlim=INDEF,
            rows="-")
        if (tstat.stddev > 0) {
            l_expression = "dmag>=("//tstat.median//"-"//l_dmagsig//"*"//tstat.stddev//") && dmag<=("//tstat.median//"+"//l_dmagsig//"*"//tstat.stddev//")"
            tselect (tmpdao//".fits", tmpime//".fits", l_expression)
            tdelete (tmpdao//".fits", verify=no)
        } else {
            printlog ("WARNING - GEMSEEING: No selection on rms dmag",
                l_logfile, yes)
            rename (tmpdao//".fits", tmpime//".fits", field="all")
        }

        tinfo (tmpime//".fits", ttout=no)
        if (tinfo.nrows < 3) {
            tdelete (tmpime//".fits", verify=no)
            goto finalselection
        }

        # select on sigma deviation in fwhm
        tstat (tmpime//".fits", "fwhm", outtable="", lowlim=INDEF, highlim=INDEF,
            rows="-")
        if (tstat.stddev > 0) {
            l_expression = "fwhm>=("//tstat.median//"-"//l_fwhmsig//"*"//tstat.stddev//") && fwhm<=("//tstat.median//"+"//l_fwhmsig//"*"//tstat.stddev//")"
            tselect (tmpime//".fits", tmpdao//".fits", l_expression)
            tdelete (tmpime//".fits", verify=no)
        } else {
            printlog ("WARNING - GEMSEEING: No selection on rms fwhm",
                l_logfile, yes)
            rename (tmpime//".fits", tmpdao//".fits", field="all")
        }

        tinfo (tmpdao//".fits", ttout=no)
        if (tinfo.nrows < 3) {
            tdelete (tmpdao//".fits", verify=no)
            rename (tmpime//".fits", tmpdao//".fits", field="all")
            goto finalselection
        }

finalselection:
        # Resulting direct FWHM 
        tstat (tmpdao//".fits", "fwhm", outtable="", lowlim=INDEF, highlim=INDEF,
            rows="-")
        n_i = strlen(s_image)
        printf ("Image  %-"//str(n_i)//"s   FWHM [arcsec]  %7.4f\n",
            s_image, tstat.median*l_pixscale) | scan(l_struct)
        printlog (l_struct, l_logfile, l_verbose)

        # mark objects in l_psffile//".fits" that have been used

        # objects not used
        tdiffer (l_psffile//".fits", tmpdao//".fits", tmpime//".fits", "id", "id")
        tinfo (tmpime//".fits",ttout=no)
        if (tinfo.nrows != 0)
            tcalc (tmpime//".fits", "flag", "0", colfmt="i2", colunit="")
        # objects used
        tcalc (tmpdao//".fits", "flag", "1", colfmt="i2", colunit="")
        tdelete (l_psffile//".fits", verify=no)

        tmerge (tmpdao//".fits,"//tmpime//".fits", l_psffile//".fits", "append", 
            allcols+, tbltype="default")
        tsort (l_psffile//".fits", "id", ascend+, casesens+)

        tdelete (tmpdao//".fits,"//tmpime//".fits", verify=no)
        delete (tmpcoo//".fits", verify=no, >& "dev$null")

    }  # end of selection

    # ----------------------------------------------------------------------
    # Rerun psfmeasure on much fewer objects to get 
    #  FWHM, EED50, EED85 and MFWHM
    # ----------------------------------------------------------------------

    tselect (l_psffile//".fits", tmpdao//".fits", "flag==1")
    tprint (tmpdao//".fits", prparam-, prdata+, pwidth=80, col="x_psf,y_psf", showr-,
        showh-, orig_row+, showu-, rows="-", option="plain", align+, sp_col="",
        lgroup=0, > tmpcoo)
    #print("q", >> tmpcoo)

    # count objects and make two lines if needed
    l_fl_one = no
    count (tmpcoo) | scan (l_obj)
    if (l_obj == 1) {
        type (tmpcoo, >> tmpcoo)
        l_fl_one = yes
    }

    # Change level before size or level will not get changed ...
    printf (":level 0.5\n:size Radius\n:show "//tmpime//".dat\n", >> tmpgra)
    printf (":level 0.85\n:size Radius\n:show "//tmpime//".dat\nq\n", >> tmpgra)

    # iterations=1 since radius should be fixed in this step if it is AO 
    # observations. The check on AO or not should be changed later to look 
    # for some Altair keyword
    if (l_inst=="NIRI-F/32" || l_inst=="Hokupaa/QUIRC" || \
        l_inst=="Hokupaa/QUIRC/Ded" || l_inst=="AO-SIM") {
        
        noao.obsutil.psfmeasure (l_image, coords="markall", wcs="logical", 
            display-, frame=1, level=0.5, size="FWHM", beta=INDEF, scale=1.,
            radius=(l_fwhmnonao*3./l_pixscale), 
            sbuffer=(l_fwhmnonao*2./l_pixscale), swidth=l_width,
            saturation=l_datamax, ignore_sat+, iterations=1, logfile="",
            xcenter=INDEF, ycenter=INDEF, imagecur=tmpcoo, graphcur=tmpgra,
            >& "dev$null", >G "dev$null")

    } else {

        noao.obsutil.psfmeasure (l_image, coords="markall", wcs="logical", 
            display-, frame=1, level=0.5, size="FWHM", beta=INDEF, scale=1.,
            radius=l_radius, sbuffer=l_buffer, swidth=l_width,
            saturation=l_datamax, ignore_sat+, iterations=l_niter, logfile="",
            xcenter=INDEF, ycenter=INDEF, imagecur=tmpcoo, graphcur=tmpgra, 
            >& "dev$null", >G "dev$null")
    }
    delete (tmpgra//","//tmpcoo, verify-, >& "dev$null")

    # make the psffile FITS table from psfmeasure output
    # need to get rid of the image name in the file ...
    printf ("x_psf r f7.2\ny_psf r f7.2\n", > tmpcd)
    printf ("eed50 r f7.3\need85 r f7.3\n", >> tmpcd)

    # Find the length of this psffile
    #tdump (tmpime//".fits", cdfile="", pfile="", datafile=tmpime//".dat",
    #    columns="", rows="-", pwidth=-1)
    count (tmpime//".dat") | scan (l_obj)
    l_obj = l_obj/2

    # For each get rows:  l_obj/2 excl. first 3 and last 2
    # There should be no saturated objects in this
    for (n_i=1; n_i<=2; n_i+=1) {
        if (n_i == 1) {
            fields (tmpime//".dat", "2-3,5", lines=str((n_i-1)*l_obj+4), > tmpdat)
            fields (tmpime//".dat", "1-2,4",
                lines=str((n_i-1)*l_obj+5)//"-"//str((n_i*l_obj-2)), >> tmpdat)
        } else {
            fields (tmpime//".dat", "5", lines=str((n_i-1)*l_obj+4), > tmpgra)
            fields (tmpime//".dat", "4",
                lines=str((n_i-1)*l_obj+5)//"-"//str((n_i*l_obj-2)), >> tmpgra)
            joinlines (tmpdat, tmpgra, output=tmpmom, delim=" ", \
                missing="Missing", maxchars=161, shortest+, verbose+)
            delete (tmpdat//","//tmpgra, verify-, >& "dev$null")
            rename (tmpmom, tmpdat, field="all")
        }
    }
    delete (tmpime//".dat", verify-, >& "dev$null")
    if (l_fl_one) {
        tcreate (tmpcoo//".fits", tmpcd, tmpdat, upar="", nskip=0, nlines=1, nrows=0,
            hist-, extrapar=5, tbltype="default", extracol=0)
        tselect (tmpcoo//".fits", tmpime//".fits", "row()==1")
        tdelete (tmpcoo//".fits", verify-)
    } else {
        tcreate (tmpime//".fits", tmpcd, tmpdat, upar="", nskip=0, nlines=1, nrows=0,
            hist-, extrapar=5, tbltype="default", extracol=0)
    }
    tmerge (tmpime//".fits,"//tmpdao//".fits",l_psffile//"_psf.fits", "merge", 
        allcols+, tbltype="default")
    delete (tmpcd//","//tmpdat//","//tmpdao//".fits,"//tmpime//".fits",
        verify-, >& "dev$null")

    # convert all measurements from pixels to arcsec, eed?? to diameters
    tcalc (l_psffile//"_psf.fits", "fwhm", "fwhm*"//str(l_pixscale),
        colfmt="f7.4", colunit="")
    tcalc (l_psffile//"_psf.fits", "mfwhm", "mfwhm*"//str(l_pixscale),
        colfmt="f7.4", colunit="")
    tcalc (l_psffile//"_psf.fits", "eed50", "eed50*2*"//str(l_pixscale),
        colfmt="f7.4", colunit="")
    tcalc (l_psffile//"_psf.fits", "eed85", "eed85*2*"//str(l_pixscale),
        colfmt="f7.4", colunit="")

    # Final determinations - median values
    if (l_fl_keep || l_fl_update)
        date | scan(l_struct)
        
    # put keywords in FITS table

    if (l_fl_keep) {
        parkey (l_struct, l_psffile//"_psf.fits", "GEMSEE", add=yes)
        if (l_fl_strehl)
            parkey ("gemseeing:  psfmeasure  direct/eed50/eed85/moffat/strehl",
                l_psffile//"_psf.fits", "SEEPRG", add=yes)
        else
            parkey ("gemseeing:  psfmeasure  direct/eed50/eed85/moffat",
                l_psffile//"_psf.fits", "SEEPRG", add=yes)

        parkey (l_image, l_psffile//"_psf.fits", "IMAGE", add=yes)
        parkey (l_image, l_psffile//"_psf.fits", "SECTION", add=yes)
        parkey (l_inst, l_psffile//"_psf.fits", "INSTRUME", add=yes)
        parkey (l_fwhmin, l_psffile//"_psf.fits", "FWHM_IN", add=yes)
        parkey (l_pixscale, l_psffile//"_psf.fits", "PIXSCALE", add=yes)
        parkey (l_ron, l_psffile//"_psf.fits", "RON", add=yes)
        parkey (l_gain, l_psffile//"_psf.fits", "EPADU", add=yes)
        parkey (l_sigthres, l_psffile//"_psf.fits", "SIGTHRES", add=yes)
        parkey ("psfmeasure", l_psffile//"_psf.fits", "FITPRG", add=yes)
        parkey ("direct", l_psffile//"_psf.fits", "FITTYP", add=yes)
        parkey (l_radius, l_psffile//"_psf.fits", "FITRAD", add=yes)
        parkey (l_buffer, l_psffile//"_psf.fits", "FITBUF", add=yes)
        parkey (l_width, l_psffile//"_psf.fits", "FITWID", add=yes)
        parkey (l_mdaofaint, l_psffile//"_psf.fits", "MDAOFAI", add=yes)
        parkey (l_datamax, l_psffile//"_psf.fits", "DATAMAX", add=yes)
        parkey (l_sharplow, l_psffile//"_psf.fits", "SHARPLO", add=yes)
        parkey (l_sharphigh, l_psffile//"_psf.fits", "SHARPHI", add=yes)
        parkey (l_sharpsig, l_psffile//"_psf.fits", "SHARPSIG", add=yes)
        parkey (l_dmagsig, l_psffile//"_psf.fits", "DMAGSIG", add=yes)
        parkey (l_fwhmsig, l_psffile//"_psf.fits", "FWHMSIG", add=yes)
    }

    tstat (l_psffile//"_psf.fits", "fwhm", outtable="", lowlim=INDEF,
        highlim=INDEF)
    out_fwhm = tstat.median
    out_npsf = tstat.nrows
    if (out_npsf > 2)
        rms_fwhm = tstat.stddev
    else
        rms_fwhm = 0.

    if (l_fl_keep) {
        parkey (tstat.median, l_psffile//"_psf.fits", "FWHMPSF", add=yes)
        parkey (tstat.mean, l_psffile//"_psf.fits", "FWHMMEAN", add=yes)
        parkey (tstat.median, l_psffile//"_psf.fits", "FWHMMED", add=yes)
        parkey (tstat.stddev, l_psffile//"_psf.fits", "FWHMRMS", add=yes)
        parkey (tstat.vmin, l_psffile//"_psf.fits", "FWHMMIN", add=yes)
        parkey (tstat.vmax, l_psffile//"_psf.fits", "FWHMMAX", add=yes)
        parkey (tstat.nrows, l_psffile//"_psf.fits", "NPSF", add=yes)
    }
    if (l_fl_update) {
        gemdate ()
        gemhedit (l_image, "GEMSEE", gemdate.outdate,
            "UT Time stamp for gemseeing", delete-)
        if (l_fl_strehl)
            hedit (l_image, "SEEPRG",
                "gemseeing:  psfmeasure  direct/eed50/eed85/moffat/strehl",
                add=yes, addonly=no, delete=no, verify=no, show-, update=yes)
        else
            hedit (l_image, "SEEPRG",
                "gemseeing:  psfmeasure  direct/eed50/eed85/moffat", add=yes,
                addonly=no, delete=no, verify=no, show-, update=yes)
        gemhedit (l_image, "FWHMPSF", tstat.median,
            "FWHM of PSF [arcsec], direct measurement")
        gemhedit (l_image, "NPSF", tstat.nrows, 
            "Number of stars used by gemseeing")
    }
    tstat (l_psffile//"_psf.fits", "mfwhm", outtable="", lowlim=INDEF,
        highlim=INDEF)
    out_mfwhm = tstat.median
    rms_mfwhm = tstat.stddev
    if (l_fl_keep) {
        parkey (tstat.median, l_psffile//"_psf.fits",  "MFWHM",add=yes)
        parkey (tstat.mean, l_psffile//"_psf.fits", "MFWHMMEA", add=yes)
        parkey (tstat.median, l_psffile//"_psf.fits", "MFWHMMED", add=yes)
        parkey (tstat.stddev, l_psffile//"_psf.fits", "MFWHMRMS", add=yes)
        parkey (tstat.vmin, l_psffile//"_psf.fits", "MFWHMMIN", add=yes)
        parkey (tstat.vmax, l_psffile//"_psf.fits", "MFWHMMAX", add=yes)
    }
    if (l_fl_update)
        gemhedit (l_image, "MFWHM", tstat.median,
            "FWHM of PSF [arcsec], Moffat fit")

    tstat (l_psffile//"_psf.fits", "eed50", outtable="", lowlim=INDEF,
        highlim=INDEF)
    # clean the measurement if enough objects
    if (tstat.nrows > 15 && l_fl_cleed) {
        tstat (l_psffile//"_psf.fits", "eed50", outtable="", lowlim=INDEF,
            highlim=(1.5*tstat.median))
        tstat (l_psffile//"_psf.fits", "eed50", outtable="", lowlim=INDEF,
            highlim=(1.5*tstat.median))
    }
    out_eed50 = tstat.median
    rms_eed50 = tstat.stddev
    if (l_fl_keep) {
        parkey (tstat.median, l_psffile//"_psf.fits", "EED50", add=yes)
        parkey (tstat.mean, l_psffile//"_psf.fits", "EED50MEA", add=yes)
        parkey (tstat.median, l_psffile//"_psf.fits", "EED50MED", add=yes)
        parkey (tstat.stddev, l_psffile//"_psf.fits", "EED50RMS", add=yes)
        parkey (tstat.vmin, l_psffile//"_psf.fits", "EED50MIN", add=yes)
        parkey (tstat.vmax, l_psffile//"_psf.fits", "EED50MAX", add=yes)
    }
    if (l_fl_update)
        gemhedit (l_image, "EED50", tstat.median,
            "50-percent encircled energy diameter [arcsec]")

    tstat (l_psffile//"_psf.fits", "eed85", outtable="", lowlim=INDEF,
        highlim=INDEF)
    # clean the measurement if enough objects
    if (tstat.nrows>15 && l_fl_cleed) {
        tstat (l_psffile//"_psf.fits", "eed85", outtable="", lowlim=INDEF,
            highlim=(1.5*tstat.median))
        tstat (l_psffile//"_psf.fits", "eed85", outtable="", lowlim=INDEF,
            highlim=(1.5*tstat.median))
    }
    out_eed85 = tstat.median
    rms_eed85 = tstat.stddev
    if (l_fl_keep) {
        parkey(tstat.median, l_psffile//"_psf.fits", "EED85", add=yes)
        parkey(tstat.mean, l_psffile//"_psf.fits", "EED85MEA", add=yes)
        parkey(tstat.median, l_psffile//"_psf.fits", "EED85MED", add=yes)
        parkey(tstat.stddev, l_psffile//"_psf.fits", "EED85RMS", add=yes)
        parkey(tstat.vmin, l_psffile//"_psf.fits", "EED85MIN", add=yes)
        parkey(tstat.vmax, l_psffile//"_psf.fits", "EED85MAX", add=yes)
    }
    if (l_fl_update)
        gemhedit (l_image, "EED85", tstat.median,
            "85-percent encircled energy diameter [arcsec]")

    tstat (l_psffile//"_psf.fits", "eps", outtable="", lowlim=INDEF,
        highlim=INDEF)
    out_eps = tstat.median
    rms_eps = tstat.stddev
    if (l_fl_keep) {
        parkey (tstat.median, l_psffile//"_psf.fits", "EPSPSF", add=yes)
        parkey (tstat.mean, l_psffile//"_psf.fits", "EPSMEAN", add=yes)
        parkey (tstat.median, l_psffile//"_psf.fits", "EPSMED", add=yes)
        parkey (tstat.stddev, l_psffile//"_psf.fits", "EPSRMS", add=yes)
        parkey (tstat.vmin, l_psffile//"_psf.fits", "EPSMIN", add=yes)
        parkey (tstat.vmax, l_psffile//"_psf.fits", "EPSMAX", add=yes)
    }
    if (l_fl_update)
        gemhedit (l_image, "EPSPSF", tstat.median, "Ellipticity of PSF")

    tstat (l_psffile//"_psf.fits", "pos", outtable="", lowlim=INDEF,
        highlim=INDEF)
    out_pos = tstat.median
    rms_pos = tstat.stddev
    if (l_fl_keep) {
        parkey (tstat.median, l_psffile//"_psf.fits", "POSPSF", add=yes)
        parkey (tstat.mean, l_psffile//"_psf.fits", "POSMEAN", add=yes)
        parkey (tstat.median, l_psffile//"_psf.fits", "POSMED", add=yes)
        parkey (tstat.stddev, l_psffile//"_psf.fits", "POSRMS", add=yes)
        parkey (tstat.vmin, l_psffile//"_psf.fits", "POSMIN", add=yes)
        parkey (tstat.vmax, l_psffile//"_psf.fits", "POSMAX", add=yes)
    }
    if (l_fl_update)
        gemhedit (l_image, "POSPSF", tstat.median,
            "Position angle of PSF, CCW from X-axis")

    # --------------------------------------------------------------------
    # Determination of approximative Strehl ratio if fl_strehl=yes
    # --------------------------------------------------------------------

    if (l_fl_strehl) {
        tprint (l_psffile//"_psf.fits", col="x_psf,y_psf", showr-, showh-,
            prdata+, prparam-, pwidth=80, plength=0, orig_row+, showunits-,
            rows="-", option="plain", align+, sp_col="", lgroup=0, > tmpcoo)
        phot (l_image, skyfile="", coords=tmpcoo, output=tmpdao, plotfile="",
            datapars="", centerpars="", fitskypars="", photpars="",
            interactive-, radplots-, icommands="", gcommands="",
            wcsin="logical", wcsout="logical", cache-, verify-, update-,
            verbose-,
            datapars.scale=1, datapars.emission+, datapars.datamin=-l_datamax,
            datapars.datamax=l_datamax, datapars.noise="poisson",
            datapars.ccdread="", datapars.gain="", datapars.readnoise=l_ron,
            datapars.epadu=l_gain, datapars.exposure="", datapars.airmass="",
            datapars.filter="", datapars.obstime="", datapars.itime=1,
            centerpars.calgorithm="none",
            fitskypars.salgorithm="gauss",
            fitskypars.annulus=(out_eed85/l_pixscale+l_buffer),
            fitskypars.dannulus=l_width, fitskypars.smaxiter=10,
            fitskypars.sloclip=0, fitskypars.shiclip=0, fitskypars.snreject=50,
            fitskypars.sloreject=3, fitskypars.shireject=3, fitskypars.khist=5,
            fitskypars.binsize=0.1, fitskypars.smooth-, fitskypars.rgrow=0,
            fitskypars.mksky-,
            photpars.weighting="constant",
            photpars.apertures=str(out_eed50/l_pixscale/2.),
            photpars.zmag=25., photpars.mkapert-)

        pconvert (tmpdao, tmpime//".fits", "flux,mag,merr,msky", expr="yes", append-)
        delete (tmpdao//","//tmpcoo, verify-, >& "dev$null")
        tmerge (l_psffile//"_psf.fits,"//tmpime//".fits", tmpdao//".fits", "merge")
        tdelete (l_psffile//"_psf.fits,"//tmpime//".fits", verify-)
        rename (tmpdao//".fits", l_psffile//"_psf.fits", field="all")

        # The approximate Strehl ratio
        # use the value in the central pixel plus a statistical correction
        tinfo (l_psffile//"_psf.fits", ttout=no)
        l_obj = tinfo.nrows
        for (ii=1; ii<=l_obj; ii+=1) {
            tprint (l_psffile//"_psf.fits", col="x_psf,y_psf", showr-, showh-,
                showu-, option="plain", row=str(ii), prpar-, prdat+) | \
                scan (x_psf, y_psf)
            listpix (t_image//"["//int(x_psf+0.5)//","//int(y_psf+0.5)//"]") | \
                fields ("STDIN", "2", lines="1", quit-, print-, >> tmpdao)
        }
        printf ("peak r f8.1\n") | tcreate (tmpdao//".fits", "STDIN", tmpdao,
            upar="", nskip=0, nlines=0, nrows=0, hist-, extrapar=5,
            tbltype="default", extracol=0)
        tmerge (l_psffile//"_psf.fits, "//tmpdao//".fits", tmpcoo//".fits",
            "merge", allcols+, tbltype="default", extracol=0)
        delete (l_psffile//"_psf.fits,"//tmpdao//".fits,"//tmpdao,
            verify-, >& "dev$null")
        rename (tmpcoo//".fits", l_psffile//"_psf.fits", field="all")
        
        # Statistical correction of 1.08 for objects not being centered on 
        # pixels - approx correct for Hokupaa/QUIRC, need database info later
        
        tcalc (l_psffile//"_psf.fits", "npeak", "1.08*(peak-msky)/(flux*2.)",
            colfmt="f7.4")
        tcalc (l_psffile//"_psf.fits", "strehl", "npeak/"//str(l_npeak),
            colfmt="f7.4")

        tstat (l_psffile//"_psf.fits", "strehl", outtable="", lowlim=INDEF,
            highlim=INDEF)
        out_strehl = tstat.median
        rms_strehl = tstat.stddev
        if (l_fl_keep) {
            parkey (tstat.median, l_psffile//"_psf.fits", "STREHLRA", add=yes)
            parkey (tstat.mean, l_psffile//"_psf.fits", "STRMEAN", add=yes)
            parkey (tstat.median, l_psffile//"_psf.fits", "STRMED", add=yes)
            parkey (tstat.stddev, l_psffile//"_psf.fits", "STRRMS", add=yes)
            parkey (tstat.vmin, l_psffile//"_psf.fits", "STRMIN", add=yes)
            parkey (tstat.vmax, l_psffile//"_psf.fits", "STRMAX", add=yes)
        }
        if (l_fl_update)
            gemhedit (l_image, "STREHLRA", tstat.median, "Strehl ratio")

    }

    # Output result
    printlog ("", l_logfile, yes)
    n_i = strlen(s_image)
    n_i = max(n_i,8)
    if (l_fl_strehl) {
        printf ("%-"//str(n_i)//"s  %4s %7s %7s %7s %7s %7s %6s %5s\n",
            "Image", "NPSF", "FWHM", "MFWHM", "EED50%", "EED85%", "STREHL",
            "EPS", "POS") | scan(l_struct)
        printlog (l_struct, l_logfile, yes)
        printf ("%-"//str(n_i)//"s  %4d %7.4f %7.4f %7.4f %7.4f %7.4f %6.3f %5.1f\n",
            s_image, out_npsf, out_fwhm, out_mfwhm, out_eed50, out_eed85,
            out_strehl, out_eps, out_pos) | scan(l_struct)
        printlog (l_struct, l_logfile, yes)
        printf ("%-"//str(n_i)//"s  %4s %7.4f %7.4f %7.4f %7.4f %7.4f %6.3f %5.1f\n",
            "[stddev]", "", rms_fwhm, rms_mfwhm, rms_eed50, rms_eed85,
            rms_strehl, rms_eps, rms_pos) | scan(l_struct)
        printlog (l_struct, l_logfile, yes)
    } else {
        printf ("%-"//str(n_i)//"s  %4s %7s %7s %7s %7s %6s %5s\n",
            "Image", "NPSF", "FWHM", "MFWHM", "EED50%", "EED85%", "EPS",
            "POS") | scan(l_struct)
        printlog (l_struct, l_logfile, yes)
        printf ("%-"//str(n_i)//"s  %4d %7.4f %7.4f %7.4f %7.4f %6.3f %5.1f\n",
            s_image, out_npsf, out_fwhm, out_mfwhm, out_eed50, out_eed85,
            out_eps, out_pos) | scan(l_struct)
        printlog (l_struct, l_logfile, yes)
        printf ("%-"//str(n_i)//"s  %4s %7.4f %7.4f %7.4f %7.4f %6.3f %5.1f\n",
            "[stddev]", "", rms_fwhm, rms_mfwhm, rms_eed50, rms_eed85, rms_eps,
            rms_pos) | scan(l_struct)
        printlog (l_struct, l_logfile, yes)
    }

    printlog ("", l_logfile, yes)
    printf ("%-"//str(n_i)//"s  %s %7.4f %s\n"," ", "Pixel scale = ",
        l_pixscale, "arcsec/pixel") | scan(l_struct)
    printlog (l_struct, l_logfile, yes)

    if (l_fl_inter) {
        print ("")
        tprint (l_psffile//"_psf.fits",
            col="x_psf,y_psf,fwhm,mfwhm,eed50,eed85,strehl,eps,pos", showr-,
            showh+, showunits-, rows="-", option="plain", align+, sp_col="",
            lgroup=0)
    }

    if (fl_inter) {
        if (no == l_fl_disp)    #it has not been displayed yet, so display now.
            display (l_image, 1, >& "dev$null")
        tprint (l_psffile//"_psf.fits", col="x_psf,y_psf,fwhm", showr-, showh-,
            showunits-, rows="-", option="plain", align+, sp_col="",
            lgroup=0) | tvmark (1, "STDIN", autolog-, mark="point", color=204,
            label=yes, nxoff=5, nyoff=5, txsize=1, interactive=no)
    }

    if (no == l_fl_keep)
        delete (l_coords//","//l_psffile//","//l_psffile//".fits,"//\
            l_psffile//"_psf.fits", verify-, >& "dev$null")


end

# Copyright(c) 2002-2016 Association of Universities for Research in Astronomy, Inc.

procedure gsstandard (input, sfile, sfunction)

# Determine flux calibration from reduced, extracted spectrum of a
# spectrophotometric standard stars. This task use STANDARD and
# SENSFUNCTION to establish the flux calibration.
#
# Version   Feb 28, 2002  RC  v1.3 release
#           May 23, 2002  BM  update for IFU, check calib directory
#           Sept 20, 2002    v1.4 release

string  input       {prompt="Input image(s)"}
string  sfile       {"std",prompt="Output flux file (used by SENSFUNC)"}
string  sfunction   {"sens",prompt="Output root sensitivity function image name"}
string  sci_ext     {"SCI",prompt="Name or number of science extension"}
string  var_ext     {"VAR",prompt="Name or number of variance extension"}
string  dq_ext      {"DQ",prompt="Name or number of data quality extension"}
string  key_airmass {"AIRMASS",prompt="Header keyword for airmass"}
string  key_exptime {"EXPTIME",prompt="Header keyword for exposure time"}
bool    fl_inter    {no,prompt="Run the task interactively"}

string  starname    {"",prompt="Standard star name(s) in calibration list"}
bool    samestar    {yes,prompt="Same star in all apertures"}
string  apertures   {"",prompt="Aperture selection list"}
bool    beamswitch  {no,prompt="Beam switch spectra"}
real    bandwidth   {INDEF,prompt="Bandpass width"}
real    bandsep     {INDEF,prompt="Bandpass separation"}
real    fnuzero     {3.68E-20,prompt="Absolute flux zero point"}
string  caldir      {"onedstds$spec50cal/",prompt="Directory containing calibration data"}
string  observatory {"Gemini-North",prompt="Observatory"}
string  mag         {"",prompt="Magnitude of stars"}
string  magband     {"",prompt="Magnitude types/bands (U|B|V|R|I|J|H|K|L|Lprime|M)"}
string  teff        {"",prompt="Effective temperature of spectral types"}

bool    ignoreaps   {yes,prompt="Ignore apertures and make one sensitivity function"}
string  extinction  {"",prompt="Extinction file"}
string  out_extinction  {"extinct.dat",prompt="Output revised extinction file"}
string  function    {"spline3",prompt="Fitting function",enum="spline3|chebyshev|legendre|spline1"}
int     order       {6,min=1,prompt="Order of fit"}
string  graphs      {"sr",prompt="Graphs per frame"}
string  marks       {"plus cross box",prompt="Data mark types (marks deleted added)"}
string  colors      {"2 1 3 4",prompt="Colors (lines marks deleted added)"}

bool    verbose     {yes,prompt="Verbose?"}
string  logfile     {"",prompt="Logfile name"}
int     status      {0,prompt="Exit status (0=good)"}
struct  *scanfile   {prompt="For internal use only"}

begin

    # Local variable definitions
    string  l_input, l_sname, l_sfunction, l_sci_ext, l_var_ext, l_dq_ext
    string  l_logfile, l_apertures, l_sfile, l_extinction, l_out_extinction
    string  l_caldir, l_observat, l_key_airmass, l_key_exptime, l_magband
    string  l_teff, l_mag, l_function, l_graphs, l_marks, l_colors, l_fl_answer
    bool    l_fl_inter, l_samestar, l_beamswitch, l_verbose,l_ignoreaps
    real    l_bandwidth, l_bandsep, l_fnuzero
    int     l_order

    # Other definitions

    string  l_key_qecorrim, l_key_qestate, l_comm_qestate
    string  tmplst, listtmp
    int     nblank, i, j, ninp, magnum, magbandnum, teffnum, nerror
    int     msktype, numext, mdfpos[200], nsciext[200], snum, sindex
    string  inimg[200],sname[200],stdtemp[200],stdmagband[200]
    string  suf, img, obsmode[200], inlst, observat[200]
    bool    mdf, inqecorr_state, first_inqecorr_state
    real    stdairmass[200], stdexptime[200], stdmag[200], tmpmag
    struct  sdate

    # Query parameters

    l_input=input; l_sname=starname; l_sfile=sfile
    l_sfunction=sfunction; l_sci_ext=sci_ext
    l_var_ext=var_ext; l_dq_ext=dq_ext
    l_samestar=samestar; l_beamswitch=beamswitch
    l_bandwidth=bandwidth; l_bandsep=bandsep; l_fnuzero=fnuzero
    l_extinction=extinction; l_out_extinction=out_extinction
    l_caldir=caldir; l_observat=observatory; l_fl_inter=fl_inter
    l_key_airmass=key_airmass; l_key_exptime=key_exptime
    l_mag=mag; l_magband=magband; l_teff=teff; l_apertures=apertures
    l_function=function; l_order=order; l_graphs=graphs
    l_marks=marks; l_colors=colors; l_verbose=verbose
    l_logfile=logfile; l_ignoreaps=ignoreaps

    status = 0
    l_fl_answer = "no"
    l_key_qecorrim = "QECORRIM"
    l_key_qestate = "QESTATE"
    l_comm_qestate = "Input Standard Data QE corrected?"
    inqecorr_state = no
    first_inqecorr_state = no

    # Make temp files

    tmplst = mktemp("tmplst")
    listtmp = mktemp("tmpchklst")

    # Cache some important tasks

    cache ("imgets", "gimverify", "gextverify", "gemhedit")
    cache ("specred.standard", "specred.sensfunc", "gemdate")

    # Set up date

    date | scan(sdate)

    #
    # Check the logfile
    #

    if (l_logfile == "STDOUT") {
        l_logfile = ""
    } else if (l_logfile == "" || stridx(" ",l_logfile)>0) {
        l_logfile = gmos.logfile
        if (l_logfile == "" || stridx(" ",l_logfile)>0) {
            l_logfile = "gmos.log"
            printlog ("WARNING - GSSTANDARD: both gsstandard.logfile and \
                gmos.logfile are empty.", l_logfile, l_verbose)
            printlog ("                      Using default file gmos.log.",
                l_logfile, l_verbose)
        }
    }

    # Logfile: what will be done:

    printlog ("---------------------------------------------------------------\
        -----------------", l_logfile, yes)
    printlog ("GSSTANDARD -- "//sdate, l_logfile, verbose+)
    printlog ("", l_logfile, yes)
    printlog ("input           = "//l_input, l_logfile, l_verbose)
    printlog ("starname        = "//l_sname, l_logfile, l_verbose)
    printlog ("sfile           = "//l_sfile, l_logfile, l_verbose)
    printlog ("sfunction       = "//l_sfunction, l_logfile, l_verbose)
    printlog ("sci_ext         = "//l_sci_ext, l_logfile, l_verbose)
    printlog ("var_ext         = "//l_var_ext, l_logfile, l_verbose)
    printlog ("dq_ext          = "//l_dq_ext, l_logfile, l_verbose)
    printlog ("samestar        = "//l_samestar, l_logfile, l_verbose)
    printlog ("apertures       = "//l_apertures, l_logfile, l_verbose)
    printlog ("beamswitch      = "//l_beamswitch, l_logfile, l_verbose)
    printlog ("fnuzero         = "//l_fnuzero, l_logfile, l_verbose)
    printlog ("extinction      = "//l_extinction, l_logfile, l_verbose)
    printlog ("out_extinction  = "//l_out_extinction, l_logfile, l_verbose)
    printlog ("caldir          = "//l_caldir, l_logfile, l_verbose)
    printlog ("observatory     = "//l_observat, l_logfile, l_verbose)
    printlog ("key_airmass     = "//l_key_airmass, l_logfile, l_verbose)
    printlog ("key_exptime     = "//l_key_exptime, l_logfile, l_verbose)
    printlog ("mag             = "//l_mag, l_logfile, l_verbose)
    printlog ("magband         = "//l_magband, l_logfile, l_verbose)
    printlog ("teff            = "//l_teff, l_logfile, l_verbose)
    printlog ("ignoreaps       = "//l_ignoreaps, l_logfile,l_verbose)
    printlog ("function        = "//l_function, l_logfile, l_verbose)
    printlog ("order           = "//l_order, l_logfile, l_verbose)
    printlog ("fl_inter        = "//l_fl_inter, l_logfile, l_verbose)
    printlog ("fl_answer       = "//l_fl_answer, l_logfile, l_verbose)
    printlog ("", l_logfile, verbose+)

    # Now, we start with a lot of verifications.

    nerror = 0

    # Check if the input file is not empty string

    if ((l_input=="") || (l_input==" ")) {
        printlog ("ERROR - GSSTANDARD: input image(s) or list are not \
            specified", l_logfile, verbose+)
        nerror = nerror+1
    }

    # Check existence of input list

    if (substr(l_input,1,1) == "@") {
        inlst = substr (l_input, 2, strlen(l_input))
        if (!access(inlst)) {
            printlog ("ERROR - GSSTANDARD: Input list "//inlst//" not found",
                l_logfile, verbose+)
            nerror = nerror+1
        }
    }

    # Check if input images exists

    if (substr(l_input,1,1)=="@") {
        inlst = substr(l_input,2,strlen(l_input))
        sections ("@"//inlst, > tmplst)
    } else {
        files (l_input, sort-, > tmplst)
    }

    scanfile = tmplst
    while(fscan(scanfile,img) != EOF) {
        gimverify (img)
        if (gimverify.status > 0) {
            printlog ("ERROR - GSSTANDARD: Input image "//img//" does not \
                exist.", l_logfile, verbose+)
            nerror = nerror+1
        }
    }

    scanfile = ""

    # Check that sfile is not an empty string

    if ((l_sfile=="") || (l_sfile==" ")) {
        printlog ("ERROR - GSSTANDARD: Output flux file is not specified",
            l_logfile, verbose+)
        nerror = nerror+1
    }

    # Check that sensitivity function is not an empty string

    if ((l_sfunction=="") || (l_sfunction==" ")) {
        printlog ("ERROR - GSSTANDARD: Output sensitivity function is not \
            specified", l_logfile, verbose+)
        nerror = nerror+1
    }

    # Determine if the sensitivity file and function already exist.

    if (access(l_sfile)) {
        printlog ("ERROR - GSSTANDARD: Output flux file "//l_sfile//\
            " already exist.", l_logfile,verbose+)
        nerror = nerror+1
    }

    gimverify (l_sfunction)
    # the sensitivity function is a simple FITS images,
    if (gimverify.status == 4) {
        printlog ("ERROR - GSSTANDARD: Output sensitivity function "//\
            l_sfunction//" already exist.", l_logfile, verbose+)
        nerror = nerror+1
    }

    # Check that the extension (SCI,VAR,DQ) are not empty.
    gextverify (l_sci_ext)
    l_sci_ext = gextverify.outext
    if (gextverify.status == 1) {
        printlog ("ERROR - GSSTANDARD: sci_ext is an empty string.",
            l_logfile, verbose+)
        nerror = nerror+1
    }

    gextverify (l_var_ext)
    l_var_ext = gextverify.outext
    if (gextverify.status == 1) {
        printlog ("ERROR - GSSTANDARD: var_ext is an empty string.",
            l_logfile, verbose+)
        nerror = nerror+1
    }

    gextverify (l_dq_ext)
    l_dq_ext = gextverify.outext
    if (gextverify.status == 1) {
        printlog ("ERROR - GSSTANDARD: dq_ext is an empty string.",
            l_logfile, verbose+)
        nerror = nerror+1
    }
    # Check the AIRMASS, EXPTIME and observatory parameters are not empty

    if ((l_key_airmass=="") || (l_key_airmass==" ")) {
        printlog ("ERROR - GSSTANDARD: airmass keyword parameter is not \
            specified", l_logfile, verbose+)
        nerror = nerror+1
    }

    if ((l_key_exptime=="") || (l_key_exptime==" ")) {
        printlog ("ERROR - GSSTANDARD: exptime keyword parameter is not \
            specified", l_logfile, verbose+)
        nerror = nerror+1
    }

    if (l_observat == "" || l_observat == " ") {
        printlog ("ERROR - GSSTANDARD: parameter observatory is not specified",
            l_logfile, verbose+)
        nerror = nerror+1
    }

    # Check the calibration directory. Should not be empty.

    if ((l_caldir == "") || (l_caldir==" ")) {
        printlog ("ERROR - GSSTANDARD: parameter caldir is not specified",
            l_logfile, verbose+)
        nerror = nerror+1
    }

    # If nerror!=0, go outside

    if (nerror > 0) {
        goto outerror
    }

    scanfile = ""
    nerror = 0

    # Check input images

    scanfile = tmplst
    #copy(tmplst,"quepasa.lst")
    i = 0
    mdf = no

    while (fscan(scanfile,img) != EOF) {

        i+=1
        suf = substr (img, strlen(img)-3, strlen(img))
        if (suf!="fits") {
            inimg[i] = img//".fits"
        } else {
            inimg[i] = img
        }
        # Check if the images are MEF
        imgets (inimg[i]//"[0]", "EXTEND", >& "dev$null")
        if (imgets.value=="F" || imgets.value=="0") {
            printlog ("ERROR - GSSTANDARD: image "//img//" is not a MEF file.",
                l_logfile, verbose+)
            nerror = nerror+1
        }

        # Check image type
        imgets (inimg[i]//"[0]", "MASKTYP", >& "dev$null")
        msktype = int(imgets.value)
        if (msktype != 1 && msktype != -1) {
            printlog ("ERROR - GSSTANDARD: "//img//" has MASKTYP other than \
                LONGSLIT, MOS, or IFU mode.", l_logfile, verbose+)
            nerror = nerror+1
        }
        # Check the MDF file.  Check where is it.
        imgets (inimg[i]//"[0]", "NEXTEND", >& "dev$null")
        numext = int(imgets.value)
        for (k=1; k<=numext; k+=1) {
            keypar (inimg[i]//"["//k//"]", keyword="EXTNAME", silent=yes)
            if (keypar.value == "MDF" && keypar.found){
                mdf = yes
                mdfpos[i] = k
            }
        }
        if (!mdf) {
            printlog ("ERROR - GSSTANDARD: Input file "//img//" does not have \
                an attached MDF table.", l_logfile, verbose+)
            nerror = nerror+1
        }
        # Get obsmode
        imgets (inimg[i]//"[0]", "OBSMODE", >& "dev$null")
        if (imgets.value == "" || imgets.value == " " || imgets.value == "0") {
            printlog ("ERROR - GSSTANDARD: Header keyword parameter OBSMODE \
                in "//img//" not found.",l_logfile,verbose+)
            nerror = nerror+1
        } else {
            obsmode[i] = imgets.value
        }
        # Check how many SCI extension we have
        imgets (inimg[i]//"[0]", "NSCIEXT", >& "dev$null")
        nsciext[i] = int(imgets.value)

        if (nsciext[i] == 0) {
            printlog ("ERROR - GSSTANDARD: Number of "//l_sci_ext//\
                " extensions unknown for image "//img, l_logfile, verbose+)
            nerror = nerror+1
        }
        for (j=1; j<=nsciext[i]; j+=1) {
            if (!imaccess(inimg[i]//"["//l_sci_ext//","//str(j)//"]")) {
                printlog ("ERROR - GSSTANDARD: can not access "//img//\
                    "["//l_sci_ext//","//str(j)//"]", l_logfile, verbose+)
                nerror = nerror+1
            }
        }

        # Check AIRMASS and EXPTIME keyword inside the header
        imgets (inimg[i]//"[0]", l_key_airmass, >& "dev$null")
        if (imgets.value=="" || imgets.value==" " || imgets.value=="0") {
            printlog ("ERROR - GSSTANDARD: Image header parameter "//\
                l_key_airmass//" not found in "//img, l_logfile, verbose+)
            nerror = nerror+1
        } else {
            stdairmass[i] = real(imgets.value)
            if (stdairmass[i] == 0.0) {
                printlog ("WARNING - GSSTANDARD: "//l_key_airmass//" for "//\
                    img//" is 0.0 Using airmass=1.", l_logfile, verbose+)
                stdairmass[i] = 1.0
            }
        }

        imgets (inimg[i]//"[0]", l_key_exptime, >& "dev$null")
        if (imgets.value=="" || imgets.value==" " || imgets.value=="0") {
            printlog ("ERROR - GSSTANDARD: Image header parameter "//\
                l_key_exptime//" not found in "//img, l_logfile, verbose+)
            nerror = nerror+1
        } else {
            stdexptime[i] = real(imgets.value)
        }

        # Now, check if the header parameter OBSERVAT exist, if not,
        # observatory[i]=l_observat

        imgets (inimg[i]//"[0]", "OBSERVAT", >& "dev$null")
        if (imgets.value=="" || imgets.value==" " || imgets.value=="0") {
            observat[i] = l_observat
        } else {
            observat[i] = imgets.value
        }

        # Now check if the standard stars have been GSEXTRACTed/GFAPSUMed.
        # If not, error.
        if (msktype==1) {
            imgets (inimg[i]//"[0]", "GSEXTRAC", >& "dev$null")
            if (imgets.value == "" || imgets.value == " " || \
                imgets.value=="0") {
                    printlog ("ERROR - GSSTANDARD: Spectra in "//img//\
                        " were not GSEXTRACTed. Run GSEXTRACT first.", \
                        l_logfile, verbose+)
                    nerror = nerror+1
            }
        } else if (msktype == -1) {
            imgets (inimg[i]//"[0]", "GFAPSUM", >& "dev$null")
            if (imgets.value == "" || imgets.value == " " || \
                imgets.value=="0") {
                    printlog ("ERROR - GSSTANDARD: Spectra in "//img//\
                        " were not GFAPSUMed. Run GFAPSUM first.", \
                        l_logfile, l_verbose)
                    nerror = nerror+1
            }
        }

        # Check the QE state of the input image
        keypar (inimg[i]//"[0]", l_key_qecorrim, silent+)
        if (keypar.found) {
            # Record the QE state of the inout images (default is no)
            inqecorr_state = yes
        } else {
            inqecorr_state = no
        }

        # Compare the QE state to all images
        if (i == 1) {
            first_inqecorr_state = inqecorr_state
        } else if (first_inqecorr_state != inqecorr_state) {
            # Cannot use data sets with differing QE correction states
            printlog ("ERROR - GSSTANDARD: "//img//" does not have the "//\
                "same QE corretion state as the previous images. ", \
                l_logfile, verbose+)
            nerror += 1
        }

    } #end of for loop over input images to check the parameters.

    ninp = i
    i = 0

    # Check if starname is not an empty string

    if ((l_sname=="") || (l_sname==" ")) {
        printlog ("ERROR - GSSTANDARD: Standard star name(s) starname or \
            list is not specified", l_logfile, verbose+)
        nerror = nerror+1
    }

    # Check existence of starname list

    if (substr(l_sname,1,1) == "@") {
        inlst = substr(l_sname,2,strlen(l_sname))
        if (!access(inlst)) {
            printlog ("ERROR - GSSTANDARD: starname list "//inlst//\
                " not found", l_logfile, verbose+)
            nerror = nerror+1
        }
    }

    if (nerror > 0) goto outerror

    delete (tmplst, verify-, >& "dev$null")

    # Generating a list with the starnames

    if (substr(l_sname,1,1) == "@") {
        inlst = substr(l_sname,2,strlen(l_sname))
        sections ("@"//inlst, > tmplst)
    } else {
        files (l_sname, sort-, > tmplst)
    }

    snum = 0
    scanfile = tmplst

    while (fscan(scanfile,img) != EOF) {
        if (img != "" || img != " ") {
            snum+=1
            sname[snum] = strlwr(img)
            img = ""
        }
    }
    scanfile = ""

    if (snum==0) {
        printlog ("ERROR - GSSTANDARD: Standard starnames in calibration \
            list are all an empty string.", l_logfile, verbose+)
        nerror = nerror+1
    }

    # imred.standard does not require a standards.men provided the user
        # gave correct starname. Task does not have to break if standards.men
        # does not exist.
    # Check that caldir$starname exists

    nerror = 0
    for (i=1; i<=snum; i=i+1) {
        if (access(l_caldir//sname[i]//".dat") == no) {
            printlog ("ERROR - GSSTANDARD: Starname "//sname[i]//" not found \
                in current user defined caldir "//l_caldir,
                l_logfile, verbose+)
            printlog ("ERROR - GSSTANDARD: specred.standard requires the \
                <starname>.dat file name to be lowercase. This could be the \
                problem if you are using a custom file.", l_logfile, verbose+)
            nerror = nerror + 1
        }
    }

    if (nerror > 0) goto outerror


    # Check if the number of input images is equal to the number of
    # standard star names

    if (snum != ninp && snum != 1) {
        printlog ("ERROR - GSSTANDARD: number of star names and input images \
            are different.", l_logfile, verbose+)
        nerror = nerror+1
    }

    # CALDIR PARAMETER = onedstds$blackbody/

    if (l_caldir == "onedstds$blackbody/") {
        if (l_mag=="" || l_mag==" ") {
            printlog ("ERROR - GSSTANDARD: parameter mag is empty string.",
                l_logfile, verbose+)
            nerror = nerror+1
        }
        scanfile = listtmp
        if (substr(l_mag,1,1)=="@") scanfile = substr(l_mag, 2, strlen(l_mag))
        else files(l_mag, sort-, > listtmp)
        magnum = 0
        while (fscan(scanfile,img) != EOF) {
            if (img!="" || img!=" ") {
                magnum+=1
                print (img) | scan (stdmag[magnum])
                img = ""
            }
        }
        delete (listtmp, verify-, >& "dev$null")
        if (magnum==0) {
            printlog ("ERROR - GSSTANDARD: parameter mag are all empty",
                l_logfile, verbose+)
            nerror = nerror+1
        } else if (magnum!=ninp) {
            printlog ("ERROR - GSSTANDARD: number of mag entries and standard \
                stars is different.", l_logfile, verbose+)
            nerror = nerror+1
        }
        # Check the MAGBAND parameter
        if (l_magband=="" || l_magband==" ") {
            printlog ("ERROR - GSSTANDARD: parameter magband is empty string.",
                l_logfile, verbose+)
            nerror = nerror+1
        }
        scanfile = listtmp
        if (substr(l_magband,1,1)=="@")
            scanfile = substr (l_magband, 2, strlen(l_magband))
        else
            files (l_magband, sort-, > listtmp)
        ii = 0
        while (fscan(scanfile,img) != EOF) {
            if (img!="" || img!=" ") {
                magbandnum+=1
                stdmagband[magbandnum] = img
                img = ""
            }
        }
        delete (listtmp, verify-, >& "dev$null")
        if (magbandnum==0) {
            printlog ("ERROR - GSSTANDARD: parameter magband are all empty \
                strings.", l_logfile, verbose+)
            nerror = nerror+1
        } else if (magbandnum!=ninp) {
            printlog ("ERROR - GSSTANDARD: number of magband entries and \
                input standard stars is different.", l_logfile, verbose+)
            nerror+=1
        }
        # Check the TEFF parameter
        if (l_teff=="" || l_teff==" ") {
            printlog ("ERROR - GSSTANDARD: parameter teff is empty string.",
                l_logfile, verbose+)
            nerror+=1
        }
        scanfile = listtmp
        if (substr(l_teff,1,1)=="@") scanfile = substr(l_teff,2,strlen(l_teff))
        else files(l_teff, sort-, > listtmp)
        teffnum = 0
        while (fscan(scanfile,img) != EOF) {
            if (img!="") {
                teffnum+=1
                stdtemp[teffnum] = img
                img = ""
            }
        }
        delete (listtmp, verify-, >& "dev$null")
        if (teffnum == 0) {
            printlog ("ERROR - GSSTANDARD: parameter teff are all empty \
                strings.", l_logfile, verbose+)
            nerror = nerror+1
        } else if (teffnum!=ninp) {
            printlog ("ERROR - GSSTANDARD: number of teff entries and "//\
                "standard stars is different", l_logfile, verbose+)
            nerror = nerror+1
        }
    } else if (l_caldir != "onedstds$blackbody/") {
        for (i=1; i<=ninp; i+=1) {
            stdmagband[i] = ""
            stdmag[i] = -99
            stdtemp[i] = ""
            printf("stdmag[i]="//stdmag[i])
        }
    }

    # If no error found, continue, else stop

    if (nerror > 0) {
        goto outerror
    }

    # Working in images
    delete (tmplst, verify-, >& "dev$null")

    # Define 'answer'. If fl_inter=yes, then answer="YES" else answer="NO"
    if (l_fl_inter) {
        l_fl_answer = "yes"
    } else {
        l_fl_answer = "no"
    }
    # Set specred task (standard and sensfunc)
    specred.logfile = l_logfile
    specred.standard.answer = l_fl_answer
    specred.sensfunc.answer = l_fl_answer

    printlog ("", l_logfile, l_verbose)
    printlog ("GSSTANDARD: STANDARD -- adding standard stars to sensitivity \
        file", l_logfile, l_verbose)
    printlog ("", l_logfile, l_verbose)
    for (i=1; i<=ninp; i+=1) {
        suf = substr (inimg[i], strlen(inimg[i])-3, strlen(inimg[i]))
        if (suf!="fits") {
            inimg[i] = inimg[i]//".fits"
        }
        printlog ("Working on image "//inimg[i], l_logfile, l_verbose)

        # Now we are ready to get the add the standard star to the sensitivity
        # file using the task standard in the specred package.
        # Working on the number of science extension (just in case)

        sindex = 1
        if (snum != 1) {
            sindex = i
        }
        for (j=1; j<=nsciext[i]; j+=1) {
            printlog ("Slit # "//j//"; star_name = "//sname[sindex],
                l_logfile, l_verbose)
            if (stdmag[i] < -98)
                tmpmag = INDEF
        
            if (stdmagband[i] != "") {
                specred.standard (inimg[i]//"["//l_sci_ext//","//j//"]",
                    l_sfile, sname[sindex], airmass=stdairmass[i],
                    exptime=stdexptime[i], mag=tmpmag, magband=stdmagband[i],
                    teff=stdtemp[i], interact=l_fl_inter, samestar=l_samestar,
                    beam_switch=l_beamswitch, apertures=l_apertures,
                    bandwidth=l_bandwidth, bandsep=l_bandsep, \
                    fnuzero=l_fnuzero,
                    extinction=l_extinction, caldir=l_caldir,
                    observatory=observat[i], graphics="stdgraph", cursor="")
            } else {
                specred.standard (inimg[i]//"["//l_sci_ext//","//j//"]",
                    l_sfile, sname[sindex], airmass=stdairmass[i],
                    exptime=stdexptime[i], mag=tmpmag, teff=stdtemp[i],
                    interact=l_fl_inter, samestar=l_samestar,
                    beam_switch=l_beamswitch, apertures=l_apertures,
                    bandwidth=l_bandwidth, bandsep=l_bandsep, \
                    fnuzero=l_fnuzero,
                    extinction=l_extinction, caldir=l_caldir,
                    observatory=observat[i], graphics="stdgraph", cursor="")
            }

            # Update the airmass keyword in the science and variance extensions
            # specred.standard writes them to the science extension; which gets
            # used by specred.calibrate. If missing from the header AIRMASS
            # gets set to INDEF which is not FITS standard - MS
            gemhedit (inimg[i]//"["//l_sci_ext//","//j//"]", l_key_airmass, \
                stdairmass[i], "Mean airmass for the observation", delete-, \
                upfile="")
        
            if (imaccess(inimg[i]//"["//l_var_ext//","//j//"]")) {
                gemhedit (inimg[i]//"["//l_var_ext//","//j//"]", \
                    l_key_airmass, stdairmass[i], \
                    "Mean airmass for the observation", \
                    delete-, upfile="")
            }
        
            if (imaccess(inimg[i]//"["//l_dq_ext//","//j//"]")) {
                gemhedit (inimg[i]//"["//l_dq_ext//","//j//"]", \
                    l_key_airmass, stdairmass[i], \
                    "Mean airmass for the observation", \
                    delete-, upfile="")
            }
        }
        printlog ("", l_logfile, l_verbose)
    }

    # Run sensfunc only when dosensfunc is yes
    # Determine the sensitivity and extinction functions with the task
    # sensfunc. The sensfunc appears to be the same in the longslit and
    # onedspec packages. Note all the objects are run with the same standard
    # star calibrations.

    printlog ("", l_logfile, l_verbose)
    printlog ("GSSTANDARD: SENSFUNC -- determining sensitivity and extinction \
        functions", l_logfile, l_verbose)
    printlog ("", l_logfile, l_verbose)

    specred.sensfunc (l_sfile, l_sfunction, apertures=l_apertures,
        ignoreaps=l_ignoreaps, logfile=l_logfile, extinction=l_extinction,
        newextinction=l_out_extinction, observatory=l_observat,
        function=l_function, order=l_order, interactive=l_fl_inter,
        graphs=l_graphs, marks=l_marks, colors=l_colors, cursor="",
        device="stdgraph")


    # Update the QE state of the sensitivity function
    # To get here all the QE states of the input spectra will be the same
    gemhedit (l_sfunction, l_key_qestate, first_inqecorr_state, \
        l_comm_qestate, delete-)

    gemdate ()
    printlog ("", l_logfile, l_verbose)
    gemhedit (l_sfunction, "GSSTAND", gemdate.outdate,
        "UT Time stamp for GSSTANDARD")
    gemhedit (l_sfunction, "GEM-TLM", gemdate.outdate,
        "UT Last modification with GEMINI")

    status = 0
    goto clean

    # Exit with error

outerror:
    delete (tmplst, verify-, >& "dev$null")
    status = 1
    goto clean

    # clean up

clean:
    scanfile = ""
    date | scan (sdate)
    if (status > 0) {
        printlog ("ERROR - GSSTANDARD: "//nerror//" error(s) found. Exiting.",
            l_logfile, verbose+)
        printlog ("", l_logfile, verbose+)
        printlog ("GSSTANDARD done. Exit status bad -- "//sdate,
            l_logfile, verbose+)
    } else {
        printlog ("GSSTANDARD done. Exit status good -- "//sdate,
            l_logfile, verbose+)
    }
    printlog ("-------------------------------------------------------------\
        ------------------", l_logfile, l_verbose)
    printlog ("", l_logfile, verbose+)

end

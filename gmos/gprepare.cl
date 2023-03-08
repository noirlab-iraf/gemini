# Copyright(c) 2002-2017 Association of Universities for Research in Astronomy, Inc.

procedure gprepare (inimages)

# Takes GMOS raw data (all modes), 3 SCI extensions
# Updated header keywords as necessary
# MDF may be attached
#
# Version   Feb 28, 2002  CW,BM,IJ  v1.3 release
#           Jul 7, 2002   IJ  N&S support, bugfix comma-sep lists
#           Aug 9, 2002   IJ add pixscale to headers to support imcoadd
#           Aug 12, 2002  IJ fixed logfile/l_logfile use
#           Aug 19, 2002  IJ more robust for incomplete headers/incorrect input
#                            warn if processing spec data w/ fl_addmdf-
#                            mirror-> MIRROR
#           Aug 26, 2002  IJ parameter encoding
#           Sept 20, 2002 IJ v1.4 release
#           Oct 14, 2002  IJ don't use ijk from within script
#           Mar 20, 2003  BM add different pixscale's for GMOS-N,S
#           May 9, 2003   IJ changed INSTRUME logic to support old GMOS-N data
#           Aug 26, 2003  KL IRAF2.12 - new parameter, addonly, in hedit
#           Oct  2, 2008  JT Automatic selection of IFU MDF for new data

char    inimages    {prompt="Input GMOS images or list"}  # OLDP-1-input-primary-single-prefix=g
char    rawpath     {"",prompt="Path for input raw images"}  # OLDP-4
char    outimages   {"",prompt="Output images or list"}      # OLDP-1-output
char    outpref     {"g",prompt="Prefix for output images"}  # OLDP-4
bool    fl_addmdf   {no,prompt="Add Mask Definition File (LONGSLIT/MOS/IFU modes)"} # OLDP-2
char    sci_ext     {"SCI",prompt="Name of science extension"}                  # OLDP-3
char    key_mdf     {"MASKNAME",prompt="Header keyword for the Mask Definition File"} # OLDP-3
char    mdffile     {"",prompt="MDF file to use if keyword not found"}          # OLDP-2-input
char    mdfdir      {"gmos$data/", prompt="MDF database directory"}             # OLDP-2
char    gaindb      {"default",prompt="Database with gain data"} # OLDP-3
char    key_ron     {"RDNOISE",prompt="Header keyword for readout noise"}       # OLDP-3
char    key_gain    {"GAIN",prompt="Header keyword for gain (e-/ADU)"}          # OLDP-3
real    ron         {3.5,prompt="Readout noise in electrons"}                   # OLDP-3
real    gain        {2.2,prompt="Gain in e-/ADU"}                               # OLDP-3
char    logfile     {"",prompt="Logfile"}                                       # OLDP-1
bool    verbose     {yes,prompt="Verbose?"}                                     # OLDP-4
int     status      {0,prompt="Exit status (0=good)"}                           # OLDP-4
struct  *scanfile   {"",prompt="Internal use only"}                             # OLDP-4

begin

    # Define local variables:
    char    l_inimages, l_outimages, l_outpref, l_logfile, l_sci_ext
    char    l_key_ron, l_key_gain, l_key_mdf
    char    l_mdffile,l_mdfdir, l_rawpath, l_temp
    char    l_gaindb = ""
    bool    l_verbose, l_fl_addmdf
    real    l_ron, l_gain, pixscale[2,3]
    int     maxfiles, maxscext, nbad
    struct  l_struct
    char    img, filelist, inimg[200], outimg[200], amatch
    char    tmpsci, dettype
    char    sci_tsec, suf, sci_dsec, tsec, mdfname[200]
    char    obsmode[200], gaindbname
    char    typ = ""
    string  mdfprefix, pathtest, obstype, mdf_suffix
    bool    addmdf[200]
    int     ninp, len, nout, i, next_inp, next_out, n_j, next_sci, n_i
    int     temp, l_test, ybin, dum, inst, iccd
    int     n_test, n_itest
    int     junk
    char    l_dir, rphend, ifu_mdf_stub

    # Make temporary files
    filelist = mktemp ("tmpfile")

    # Set local variable values:
    l_inimages = inimages
    l_outimages = outimages
    l_outpref = outpref
    l_logfile = logfile
    l_sci_ext = sci_ext
    l_key_ron = key_ron
    l_key_gain = key_gain
    l_key_mdf = key_mdf
    l_mdffile = mdffile
    l_mdfdir = mdfdir
    l_rawpath = rawpath
    l_verbose = verbose
    l_fl_addmdf = fl_addmdf
    l_ron = ron
    l_gain = gain
    junk = fscan (gaindb, l_gaindb)

    l_dir = ""

    status = 0
    maxfiles = 200
    maxscext = 60

    #pixscale[inst,iccd]*ybin will set the header keyword PIXSCALE #hcode
    #which is the 'Pixel scale in Y in arcsec/pixel'
    pixscale[1,1] = 0.0727  #GMOS-N EEV2
    pixscale[1,2] = 0.07288 #GMOS-N e2vDD CCDs
    pixscale[1,3] = 0.0807  #GMOS-N - New Hamamatsu Detector #hcode 
    pixscale[2,1] = 0.073   #GMOS-S
    pixscale[2,3] = 0.0800  #GMOS-S - New Hamamatsu Detector ##M PIXEL_SCALE

    # Keep task parameters from changing from the outside
    cache ("imgets", "gimverify", "fparse", "gemdate")

    # Test the logfile:
    if ((l_logfile == "") || (stridx(" ", l_logfile) > 0)) {
        l_logfile = gmos.logfile
        if ((l_logfile == "") || (stridx (" ", l_logfile) > 0)) {
            l_logfile = "gmos.log"
            printlog ("WARNING - GPREPARE: Both gprepare.logfile and \
                gmos.logfile fields are empty", l_logfile, verbose+)
            printlog ("                    Using default file gmos.log",
                l_logfile, verbose+)
        }
    }

    # Start logging
    date | scan(l_struct)
    printlog ("---------------------------------------------------------------\
        -------------", l_logfile, l_verbose)
    printlog ("GPREPARE -- "//l_struct, l_logfile, l_verbose)
    printlog ("", l_logfile, l_verbose)
    printlog ("Input list    = "//l_inimages, l_logfile, l_verbose)
    printlog ("Output list   = "//l_outimages, l_logfile, l_verbose)
    printlog ("Output prefix = "//l_outpref, l_logfile, l_verbose)
    printlog ("Raw path      = "//l_rawpath, l_logfile, l_verbose)
    printlog ("MDF dir       = "//l_mdfdir, l_logfile, l_verbose)
    printlog ("Add MDF       = "//l_fl_addmdf, l_logfile, l_verbose)
    if (l_fl_addmdf)
        printlog ("Input MDF in case header keyword not found = "//l_mdffile,
            l_logfile, l_verbose)
    printlog ("", l_logfile, l_verbose)

    # Check the rawpath name for a final "/"
    if (l_rawpath != "") {
        rphend = substr(l_rawpath,strlen(l_rawpath),strlen(l_rawpath))
        if (rphend == "$") {
            show (substr(l_rawpath,1,strlen(l_rawpath)-1)) | scan (pathtest)
            rphend = substr(pathtest,strlen(pathtest),strlen(pathtest))
        }
        if (rphend != "/") {
            l_rawpath = l_rawpath//"/"
        }
        if (!access(l_rawpath)) {
            printlog ("ERROR - GPREPARE: Cannot access rawpath: "//l_rawpath, \
                l_logfile, verbose+)
            goto clean
        }
    }

    # Check the mdfdir has a trainling "/"
    # Check that the rawpath has a trailing slash and is a valid entry
    if (l_mdfdir != "") {
        rphend = substr(l_mdfdir,strlen(l_mdfdir),strlen(l_mdfdir))
        if (rphend == "$") {
            show (substr(l_mdfdir,1,strlen(l_mdfdir)-1)) | scan (pathtest)
            rphend = substr(pathtest,strlen(pathtest),strlen(pathtest))
        }

        if (rphend != "/") {
            l_mdfdir = l_mdfdir//"/"
        }

        if (!access(l_mdfdir)) {
            printlog ("ERROR - GPREPARE: Cannot access mdfdir: "//l_mdfdir, \
                l_logfile, verbose+)
            goto crash
        }
    }

    # Test the SCI extension name:
    if ((l_sci_ext == "") || (stridx(" ",l_sci_ext) > 0)) {
        printlog ("ERROR - GPREPARE: Science extension name sci_ext is not \
            defined", l_logfile, verbose+)
        goto crash
    }

    # Test the MDF related parameters
    if (l_fl_addmdf) {
        if (((l_key_mdf == "") || (stridx(" ", l_key_mdf) > 0) ) && \
            ((l_mdffile == "") || (stridx(" ", l_mdffile) > 0))) {

            printlog ("ERROR - GPREPARE: Neither the MDF keyword key_mdf or \
                the", l_logfile, verbose+)
            printlog ("                  MDF filename mdffile are defined",
                l_logfile, verbose+)
            goto crash
        }
    }

    # Test gaindb
    gaindbname = ""    # required later on in the loop. gaindb does not have to
                       # be defined.
    if ((l_gaindb != "default") && (l_gaindb != "")) {
        gaindbname = osfn(l_gaindb)
        if (access(gaindbname) == no) {
            printlog ("ERROR - GPREPARE: Gain database not found",
                l_logfile, verbose+)
            goto crash
        }
    }

    #-------------------------------------------------------------------
    # Load up input name list: @list, * and ?, comma separated

    if ((l_inimages == "") || (stridx(" ", l_inimages) > 0)) {
        printlog ("ERROR - GPREPARE: Input file not specified", l_logfile,
            verbose+)
        goto crash
    }

    # Test for @filelist
    if ((substr(l_inimages,1,1) == "@") && \
        (access(substr(l_inimages,2,strlen(l_inimages))) == no) ) {

        printlog ("ERROR - GPREPARE: Input list "//\
            substr(l_inimages,2,strlen(l_inimages))//" does not exist",
            l_logfile, verbose+)
        goto crash
    }

    # parse wildcard and comma-separated lists
    if (substr(l_inimages,1,1) == "@") {
        scanfile = substr (l_inimages, 2, strlen(l_inimages))
        while (fscan(scanfile,l_temp) != EOF) {
            files (l_rawpath//l_temp, >> filelist)
        }
    } else {
        if (stridx(",",l_inimages) == 0)
            files (l_rawpath//l_inimages, > filelist)
        else {
            l_test = 9999
            while (l_test != 0) {
                l_test = stridx (",", l_inimages)
                if (l_test > 0)
                    files (l_rawpath//substr(l_inimages,1,l_test-1),
                        >> filelist)
                else
                    files (l_rawpath//l_inimages, >> filelist)
                l_inimages = substr (l_inimages, l_test+1, strlen(l_inimages))
            }
        }
    }
    scanfile = ""
    scanfile = filelist

    #--------------------------------------------------------------------
    # check the input images

    ninp = 0
    nbad = 0
    while (fscan(scanfile, img) != EOF) {
        # split off directory path
        fparse (img, verbose-)
        img = fparse.root
        l_dir = fparse.directory

        # Must be there, but cannot be GEIS or OIF
        gimverify (l_dir//img)
        if (gimverify.status == 0) {
            ninp = ninp+1
            if (ninp > maxfiles) {
                printlog ("ERROR - GPREPARE: Maximum number of input images \
                    exceeded", l_logfile, verbose+)
                goto crash
            }
            # now has path if relevant
            inimg[ninp] = gimverify.outname//".fits"
            imgets (l_dir//img//"[0]", "GPREPARE", >& "dev$null")
            if (imgets.value != "0") {
                printlog ("ERROR - GPREPARE: Image "//l_dir//img//" already \
                    GPREPAPREd", l_logfile, verbose+)
                nbad += 1
            }
        } else if (gimverify.status==1) {
            printlog ("ERROR - GPREPARE: Input image "//l_dir//img//" does \
                not exist", l_logfile, verbose+)
            nbad += 1
        } else {
            printlog ("ERROR - GPREPARE: Input image "//l_dir//img//" is not \
                MEF", l_logfile, verbose+)
            nbad += 1
        }
    } # end while

    if (nbad != 0) {
        printlog ("ERROR - GPREPARE: "//nbad//" input files do not exist, are \
            the wrong type, or", l_logfile, verbose+)
        printlog ("                  are already GPREPAREd", l_logfile,
            verbose+)
        goto crash
    }
    if (ninp == 0) {
        printlog ("ERROR - GPREPARE: No input images meet wildcard criteria",
            l_logfile, verbose+)
        goto crash
    }

    # inimg[ninp] now contains the input images incl. the directory path
    for (n_i=1; n_i<=ninp; n_i+=1) {
        fparse (inimg[n_i], verbose-)
        img = fparse.root
        l_dir = fparse.directory

        # Identify the OBSMODE and set the keyword value. Note that if it
        # cannot find the keyword will automatically set as IMAGE
        imgets (l_dir//img//"[0]", "MASKTYP", >& "dev$null")
        # Check whether it is BIAS or DARK - give themn an OBSMODE of IMAGE too
        # - MS
        obstype = ""
        keypar (l_dir//img//"[0]", "OBSTYPE", silent+)
        if (keypar.found) {
            obstype = str(keypar.value)
        } else {
            printlog ("ERROR - GPREPARE: OBSTYPE keyword not found", \
                l_logfile, verbose+)
            goto crash
        }

        if (imgets.value == "0" || obstype == "BIAS" || obstype == "DARK") {
            obsmode[n_i] = "IMAGE"
        } else if (imgets.value == "-1") {
            # Check that the value of the MASKNAME keyword starts with IFU to
            # ensure that the values of the MASKTYP and MASKNAME keywords are
            # consistent (see HD 3595)
            keypar (l_dir//img//"[0]", "MASKNAME", >& "dev$null")
            if (substr(keypar.value,1,3) == "IFU")
                obsmode[n_i] = "IFU"
            else {
                printlog ("ERROR - GPREPARE: The MASKTYP and MASKNAME \
                    keywords are inconsistent.", l_logfile, verbose+)
                printlog ("                  MASKTYP = -1, which implies that \
                    the data is IFU, but ", l_logfile, verbose+)
                printlog ("                  MASKNAME suggests that the data \
                    is not, since it does not ", l_logfile, verbose+)
                printlog ("                  start with 'IFU'", l_logfile,
                    verbose+)
                goto crash
            }
        } else if (imgets.value == "1") {
            # Check that the value of the MASKNAME keyword does not start with
            # IFU and is not equal to "None" to ensure that the values of the
            # MASKTYP and MASKNAME keywords are consistent (see HD 3595)
            keypar (l_dir//img//"[0]", "MASKNAME", >& "dev$null")
            if (substr(keypar.value,1,3) != "IFU" && keypar.value != "None")
                obsmode[n_i] = "MOS"
            else {
                if (substr(keypar.value,1,3) == "IFU")
                    typ = "'IFU'"
                if (keypar.value == "None")
                    typ = "'IMAGE'"
                printlog ("ERROR - GPREPARE: The MASKTYP and MASKNAME \
                    keywords are inconsistent.", l_logfile, verbose+)
                printlog ("                  MASKTYP = 1, which implies that \
                    the data is MOS, but ", l_logfile, verbose+)
                printlog ("                  MASKNAME suggests that the data \
                    is "//typ, l_logfile, verbose+)
                goto crash
            }
        } else {
            # If cannot determine what it is, assume it is IMAGE instead of
            # crashing
            printlog ("WARNING - GPREPARE: The MASKTYP keyword is not \
                recognizable - assuming OBSMODE=IMAGE", l_logfile, verbose+)
            obsmode[n_i] = "IMAGE"
        }

        # If OBSMODE = MOS, check the MASKNAME for longslit
        amatch = ""
        if (obsmode[n_i] == "MOS") {
            imgets (l_dir//img//"[0]", "MASKNAME", >>& "dev$null")
            print (imgets.value) | match ("arcsec","STDIN",stop-) | \
                match ("NS","STDIN",stop+) | scan (amatch)
            if (strlen(amatch) > 0)
                obsmode[n_i] = "LONGSLIT"
        }

        # If adding the MDF, test if the image is of a valid format to receive
        # an MDF (IMAGEs do not have MDFs)
        if (l_fl_addmdf) {
            if (obsmode[n_i] == "IMAGE")
                addmdf[n_i] = no
            else
                addmdf[n_i] = yes
        }
        # Warn if not adding MDF to a spectroscopic observation
        if (!l_fl_addmdf && (obsmode[n_i] != "IMAGE"))
            printlog ("WARNING - GPREPARE: Spectroscopic frame "//img//\
                " will have no MDF, use fl_addmdf=yes", l_logfile, l_verbose)

        # Double check the grating choice
        # Checking for it here will allow the MDF to be attached with
        # OBSMODE=IMAGE

        imgets (l_dir//img//"[0]", "GRATING", >& "dev$null")
        if ((imgets.value == "MIRROR") && (obsmode[n_i] != "IMAGE")) {
            printlog ("WARNING - GPREPARE: Mask or IFU used without grating, \
                setting OBSMODE=IMAGE", l_logfile, verbose+)
            obsmode[n_i] = "IMAGE"
        }

    } # end of for loop over input images -- while(img) loop

    scanfile = ""
    delete (filelist, verify-, >& "dev$null")

    #--------------------------------------------------------------------
    # Output name list: can be empty if prefix is defined, @list or
    # comma-separated list

    if (stridx(" ", l_outimages) > 0)
        l_outimages = ""

    if ((substr(l_outimages,1,1) == "@") && \
        (access(substr(l_outimages,2,strlen(l_outimages))) == no) ) {

        printlog ("ERROR - GPREPARE: File "//\
            substr(l_outimages,2,strlen(l_outimages))//" does not exist",
            l_logfile, verbose+)
        goto crash
    }

    files (l_outimages, sort-, > filelist)
    scanfile = filelist
    nout = 0

    # If empty string, prefix must be defined
    if (l_outimages == "") {
        if ((l_outpref == "") || (stridx(" ",l_outpref) > 0) ) {
            printlog ("ERROR - GPREPARE: Neither output name nor output \
                prefix is defined", l_logfile, verbose+)
            goto crash
        }
        nout = ninp
        n_i = 1

        while (n_i<=nout) {
            fparse (inimg[n_i], verbose=no)
            outimg[n_i] = l_outpref//fparse.root//".fits"
            gimverify (outimg[n_i])
            if (gimverify.status != 1) {
                printlog ("ERROR - GPREPARE: Output image "//outimg[n_i]//\
                    " already exists", l_logfile, verbose+)
                nbad += 1
            }
            n_i = n_i + 1
        } # end of while(n_i) loop

    } else {
        while (fscan(scanfile, img) != EOF) {
            gimverify (img)
            if (gimverify.status != 1) {
                printlog ("ERROR - GPREPARE: Output image "//img//" already \
                    exists", l_logfile, verbose+)
                nbad += 1
            }
            nout = nout+1
            if (nout > maxfiles) {
                printlog ("ERROR - GPREPARE: Maximum number of output images \
                    exceeded", l_logfile, verbose+)
                goto crash
            }
            outimg[nout] = gimverify.outname//".fits"
        } # end of while (img) loop

    } # end else

    scanfile = ""
    delete (filelist, verify-, >& "dev$null")

    if (nbad > 0)
        goto crash

    # Input and output number must be the same, if output names are defined
    if ((ninp != nout) && (l_outimages != "")) {
        printlog ("ERROR - GPREPARE: Number of input and output images are \
            not the same", l_logfile, verbose+)
        goto crash
    }
    #---------------------------------------------------------------
    # Start making the output images
    n_i = 1

    while (n_i<=ninp) {
        img = inimg[n_i]

        # Which instrument?
        inst = 1 # Default is GMOS-N, support for old data

        imgets (img//"[0]", "INSTRUME", >& "dev$null")
        if (imgets.value == "0") {
            printlog ("ERROR - GPREPARE: Instrument keyword not found.",
                l_logfile, verbose+)
            goto crash
        } else if (imgets.value == "GMOS-S") {
            inst=2 # GMOS-S
        }

        # Which detector type?
        keypar (img//"[0]", "DETTYPE", silent+)
        if (!keypar.found) {
            printlog ("ERROR - GPREPARE: DETTYPE not found in "//\
                img//"[0]", l_logfile, verbose+)
            goto crash
        }
        dettype = keypar.value

        # Default MDF extension to copy out to 1. Only gets reset for GMOS-S
        # Hammamatsu CCDs when using IFU MDFs in the package
        ifu_mdf_stub = ""
        if (dettype == "SDSU II CCD") {
            # CCDs GMOS-N and GMOS-S on 2011-07-24
            iccd = 1
        } else if (dettype == "SDSU II e2v DD CCD42-90") {
            # New e2v CCDs GMOS-N
            iccd = 2
        } else if (dettype == "S10892") {
            # GMOS-S Hamamatsu CCDs
            iccd = 3
            ifu_mdf_stub = "_HAM"
        } else if (dettype == "S10892-N") {
            iccd = 3
            ifu_mdf_stub = ""
        } else {
            printlog ("ERROR - GPREPARE: DETTYPE not recognised for "//\
                img, l_logfile, verbose+)
            goto crash
        }

        printlog ("Input  "//img//"   Output  "//outimg[n_i],
            l_logfile,l_verbose)

        # Raw images do not have NEXTEND in the header, so get the number of
        # extensions from fxheader
        # Most peculiar bug?, need the extra fxhead to make sure the lines get
        # counted correctly???, IJ
        fxhead (img, format_file="", long_head-, count_lines-, >& "dev$null")
        fxhead (img, format_file="", long_head-, count_lines-) | \
            count ("STDIN") | scan (next_inp)
        fxhead (img, format_file="", long_head-, count_lines-) | \
            fields ("STDIN", "1", lines=str(next_inp), quit_if_miss-,
            print_file-) | scan (next_inp)
        # Test if the MDF keyword exists in the header, check if it is the
        # same as the input file.
        if (l_fl_addmdf && addmdf[n_i]) {
            imgets (img//"[0]", l_key_mdf, >& "dev$null")
            if (imgets.value == "0" || imgets.value == "IFU") {
                printlog ("GPREPARE: Using MDF defined in the parameter \
                    list "//l_mdffile, l_logfile, verbose+)
                if (imgets.value=="IFU")
                  printlog ("  (header value 'IFU' in old data is not unique)",
                    l_logfile, verbose+)
                mdfname[n_i] = l_mdffile
            } else if (obsmode[n_i]=="IFU") {
              mdfname[n_i] = imgets.value
              printlog ("GPREPARE: Using MDF corresponding to "//mdfname[n_i]//
                        " in the header ", l_logfile, verbose+)
              if (inst==1)
                mdfprefix="gnifu_"
              else
                mdfprefix="gsifu_"
              if (mdfname[n_i] == "IFU-2")
                mdf_suffix = "slits_mdf"
              else if (mdfname[n_i] == "IFU-B")
                mdf_suffix = "slitb_mdf"
              else if (mdfname[n_i] == "IFU-R")
                mdf_suffix = "slitr_mdf"
              else if (mdfname[n_i] == "IFU-NS-2")
                mdf_suffix = "ns_slits_mdf"
              else if (mdfname[n_i] == "IFU-NS-B")
                mdf_suffix = "ns_slitb_mdf"
              else if (mdfname[n_i] == "IFU-NS-R")
                mdf_suffix = "ns_slitr_mdf"

              mdfname[n_i] = mdfprefix//mdf_suffix//ifu_mdf_stub//".fits"
              printlog ("  (-> "//mdfname[n_i]//")", l_logfile, verbose+)

            } else { # MDF in the header and obsmode is not IFU
                mdfname[n_i] = imgets.value
                printlog ("GPREPARE: Using MDF defined in the header "//\
                    mdfname[n_i], l_logfile, verbose+)
            }
            len = strlen (mdfname[n_i])
            suf = substr (mdfname[n_i], len-4, len)
            if (suf != ".fits" )
                mdfname[n_i] = mdfname[n_i]//".fits"
            if (access(mdfname[n_i]))
                printlog ("GPREPARE: Taking MDF from the current directory.",
                    l_logfile, l_verbose)
            else {
                printlog ("GPREPARE: Taking MDF from directory "//mdfdir,
                    l_logfile, l_verbose)
                mdfname[n_i] = l_mdfdir//mdfname[n_i]
            }
            if ((mdfname[n_i] == "") || (stridx(" ",mdfname[n_i]) > 0)) {
                printlog ("ERROR - GPREPARE: The MDF filename is not defined.",
                    l_logfile, verbose+)
                goto crash
            }
            if (access(mdfname[n_i]) == no) {
                printlog ("ERROR - GPREPARE: The MDF file does not exist.",
                    l_logfile, verbose+)
                goto crash
            }

        } # end of fl_addmdf

        #-------------------------
        #If no problems found, creates the new image and insert the MDF
        #if requested:
        copy (img, outimg[n_i], verbose-, >& "dev$null")

        # Update EXTNAME and EXTVER, assuming all extensions but the MDF
        # are SCI
        n_j = 1
        next_sci = 0
        while (n_j <= next_inp) {
            keypar (img//"["//n_j//"]", "EXTNAME", silent-, >& "dev$null")
            if (keypar.value == "MDF")
                addmdf[n_i] = no
            else {
                gemhedit (outimg[n_i]//"["//n_j//"]", "EXTNAME", l_sci_ext, \
                    "", delete-)
                gemhedit (outimg[n_i]//"["//n_j//"]", "EXTVER", n_j, "", \
                    delete-)

                # Set DISPAXIS if not in imaging mode
                if (obsmode[n_i] != "IMAGE") {
                    gemhedit (outimg[n_i]//"["//n_j//"]", "DISPAXIS", 1,
                        "Dispersion Axis", delete-)
                }

                next_sci = next_sci + 1
            }
            n_j = n_j+1
        }
        flpr # necessary to work around an IRAF bug

        if (next_sci > maxscext) {
            printlog ("ERROR - GPREPARE: Maximum number of "//l_sci_ext//\
                " extensions exceeded.", l_logfile, verbose+)
            imdelete (outimg[n_i], verify-, >& "dev$null")
            goto crash
        }

        if (l_fl_addmdf && addmdf[n_i]) {
            fxinsert (mdfname[n_i], outimg[n_i]//"["//next_inp//"]",
                groups="1", verbose-, >& "dev$null")
            next_inp = next_inp+1
            # Clean old MDFs with double EXTNAME
            n_test = 1
            tprint (outimg[n_i]//"["//next_inp//"]", prpar+, prdat-,
                option="plain") | match ("EXTNAME", "STDIN", stop-) | \
                count("STDIN") | scan(n_test)
            for (n_itest=1; n_itest<n_test; n_itest+=1)
                thedit (outimg[n_i]//"["//next_inp//"]", "EXTNAME", "",
                    delete+, show-)
            parkey ("MDF", outimg[n_i]//"["//next_inp//"]", "EXTNAME", add+)
            parkey ("1", outimg[n_i]//"["//next_inp//"]","EXTVER", add+)
        }
        next_out = next_inp

        # Update gain and readnoise
        gemhedit (outimg[n_i]//"[0]", "NSCIEXT", next_sci,
            "Number of science extensions", delete-)

        if (l_gaindb == "default") {
            ggdbhelper (outimg[n_i]//"[0]", logfile=l_logfile)
            if (ggdbhelper.status != 0)
                goto crash
            gaindbname = osfn(ggdbhelper.gaindb)
            if (access(gaindbname) == no) {
                printlog ("ERROR - GPREPARE: Gain database not found",
                    l_logfile, verbose+)
                goto crash
            }
        }

        if (gaindbname != "") {
            ggain (outimg[n_i], gaindb=gaindbname, logfile=l_logfile,
                key_gain=l_key_gain, key_ron=l_key_ron, gain=l_gain, ron=l_ron,
                fl_update+, fl_mult-, verbose=no, sci_ext=l_sci_ext,
                var_ext="VAR")  # because fl_mult- here, variance ext \
                                # name doesn't matter
            if (ggain.status != 0) {
                printlog ("ERROR - GPREPARE: error in GGAIN", l_logfile,
                    verbose+)
                goto crash
            }
        }

        # Now update the PHU
        # get the pixscale, use binning in Y
        ybin = 1
        keypar (outimg[n_i]//"["//l_sci_ext//",1]", "CCDSUM", silent+)
        if (keypar.found)
            print(keypar.value) | scan(dum,ybin)
        else
            printlog ("WARNING - GPREPARE: Cannot find CCDSUM, assume \
                PIXSCALE="//pixscale[inst,iccd], l_logfile, l_verbose)
        gemhedit (outimg[n_i]//"[0]", "PIXSCALE", (pixscale[inst,iccd]*ybin),
            "Pixel scale in Y in arcsec/pixel", delete-)

        gemdate ()
        gemhedit (outimg[n_i]//"[0]", "NEXTEND", next_out,
            "Number of extensions", delete-)
        gemhedit (outimg[n_i]//"[0]", "OBSMODE", obsmode[n_i],
            "Observing mode (IMAGE|IFU|MOS|LONGSLIT)", delete-)
        gemhedit (outimg[n_i]//"[0]", "GEM-TLM", gemdate.outdate,
            "UT Last modification with GEMINI", delete-)
        gemhedit (outimg[n_i]//"[0]", "GPREPARE", gemdate.outdate,
            "UT Time stamp for GPREPARE", delete-)

        # update only if keyword not defined
        if (l_fl_addmdf && addmdf[n_i] && (l_key_mdf == "")) {
            fparse (mdfname[n_i], verbose-)
            gemhedit (outimg[n_i]//"[0]", "MASKNAME", fparse.root, \
                "Input MDF", delete-)
        }
        n_i += 1

    } # end of while(i) loop

    goto clean

    #--------------------------------------------------------------------------
crash:
    # Exit with error subroutine
    status = 1
    goto clean

clean:
    # clean up
    delete (filelist, verify-, >& "dev$null")
    scanfile = ""

    #close log file
    if (status == 0) {
        printlog ("", l_logfile, l_verbose)
        printlog ("GPREPARE exit status: good.", l_logfile, l_verbose)
        printlog ("-----------------------------------------------------------\
            -----------------", l_logfile, l_verbose)
    } else {
        printlog ("", l_logfile, l_verbose)
        printlog ("GPREPARE exit status: error.", l_logfile, l_verbose)
        printlog ("-----------------------------------------------------------\
            -----------------", l_logfile, l_verbose)
    }

end

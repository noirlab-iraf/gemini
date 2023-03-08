# Copyright(c) 2002-2012 Association of Universities for Research in Astronomy, Inc.

procedure gfdisplay (input, frame)

# Display extracted GMOS IFU spectra using ldisplay
# The scale of the output image is 0.0125arcsec/pixel
#
# Version   Sept 20, 2002 BM  v1.4 release
#           Apr  03, 2004 BM  get ldisplay config files from MDF by default
#           Apr  07, 2004 BM  fix selection for sky-subtracted cubes
#           Aug  30, 2006 JT  fix for ldisplay loop in 1-slit mode

string  input   {prompt="Input GMOS IFU extracted spectra"}
int     frame   {1,prompt="Frame for display"}
real    z1      {0.0,prompt="Minimum greylevel to be displayed"}
real    z2      {0.0,prompt="Maximum greylevel to be displayed"}
string  ztrans  {"log",enum="linear|log",prompt="Greylevel transformation (linear|log)"}
string  output  {"",prompt="Name of output image"}
string  extname {"SCI",prompt="Which EXTNAME to display"}
string  version {"*",enum="1|2|*",prompt="Which EXTVER to display (1,2,*)"}
string  config  {"default",prompt="Ldisplay configuration file"}
string  deadfib {"default",prompt="Ldisplay dead fiber configuration file"}
int     status  {0,prompt="Exit status (0=good)"}

begin

    string  l_input,l_config, l_ver, l_extn, l_dead,l_output,l_ztrans
    int     l_frame
    real    l_z1,l_z2

    bool    rmout,rmconf
    int     nx1,nx2,len
    real    ref1,ref2,refpix
    string  tmpmsa,tmpmsb,msjoin,tmpout
    string  sec1,sec2,cmd
    string  mdfmatch,tmpgood,tmpbad

    l_input=input; l_frame=frame; l_extn=extname; l_ver=version
    l_config=config; l_dead=deadfib; l_output=output; l_z1=z1; l_z2=z2
    l_ztrans=ztrans

    cache ("imgets", "gemdate")

    tmpmsa = mktemp("tmpmsa")
    tmpmsb = mktemp("tmpmsb")
    msjoin = mktemp("tmpmsjoin")
    tmpout = mktemp("tmpldisplay")
    tmpgood = mktemp("tmpgood")//".fits"
    tmpbad = mktemp("tmpbad")//".fits"
    rmout = no
    rmconf = no
    status = 0

    # check that there are input files
    if (l_input == "" || l_input == " "){
        print ("ERROR - GFDISPLAY: input files not specified")
        goto error
    }

    # check existence of list file
    if (substr(l_input,1,1) == "@") {
        print ("ERROR - GFDISPLAY: lists are currently not supported")
        goto error
    }

    # check existence of input file
    gimverify (l_input)
    if (gimverify.status != 0) {
        print ("ERROR - GFDISPLAY: "//l_input//" does not exist or is \
            not a MEF")
        goto error
    }
    l_input = gimverify.outname//".fits"

    # check existence of config files
    if (l_config!="default" && !access(l_config)) {
        print ("ERROR - GFDISPLAY: "//l_config//" does not exist.")
        goto error
    }
    if (l_dead!="default" && !access(l_dead)) {
        print ("ERROR - GFDISPLAY: "//l_dead//" does not exist.")
        goto error
    }

    # if either l_config or l_dead==default, read MDF
    if (l_config=="default" || l_dead=="default") {
        #Check if file has a MDF
        mdfmatch = ""
        fxhead (l_input, format_file="", long-, count-) | \
            match("MDF","STDIN") | scan(mdfmatch)
        if (mdfmatch=="") {
            print ("ERROR - GFDISPLAY: "//l_input//" does not contain a MDF.  \
                Please specify CONFIG and DEADFIB.")
            goto error
        }
        # Check that MDF has the proper columns
        mdfmatch = ""
        tlcol (l_input//"[MDF]", nlist=1) | match("XLDIS","STDIN") | \
            scan(mdfmatch)
        if (mdfmatch=="") {
            print ("ERROR - GFDISPLAY: "//l_input//" contains an old or \
                improper MDF.  Please specify CONFIG and DEADFIB.")
            goto error
        }
        # If here then ok, make temporary config files
        rmconf = yes
        l_config = mktemp("tmpgfdconf")
        l_dead = mktemp("tmpgfddead")
        tselect (l_input//"[MDF]", tmpgood, "BEAM >= 0")
        tselect (l_input//"[MDF]", tmpbad, "BEAM == -1")
        tprint (tmpgood, showrow-, showhdr-, showunits-, col="XLDIS,YLDIS",\
            > l_config)
        tprint (tmpbad, showrow-, showhdr-, showunits-, col="XLDIS,YLDIS", \
            > l_dead)
        delete (tmpgood//","//tmpbad, verify-, >& "dev$null")
    }


    # Command for display
    cmd = "display frame="//str(l_frame)
    if (l_z1 != 0.0 || l_z2 != 0.0) {
        cmd = cmd//" zscale- zrange- z1="//str(l_z1)//" z2="//str(l_z2)
    }
    cmd = cmd//" ztrans="//l_ztrans//" image="

    # Temporary output output image if no output image given
    if (l_output=="" || l_output==" ") {
        l_output = tmpout
        rmout = yes
    }

    if (l_ver == "1" || l_ver=="2") {
        if (!imaccess(l_input//"["//l_extn//","//l_ver//"]")) {
            print ("ERROR - GFDISPLAY: "//l_input//"["//l_extn//","//l_ver//\
                "] does not exist.")
            goto error
        }

        # If true 1-slit case, use scopy as for 2-slit case, to renumber the 
        # APNUM header keywords, fixing a loop in ldisplay for 1 slit mode, 
        # Aug 2006, JT
        hselect (l_input//"["//l_extn//","//l_ver//"]", "NAXIS2", yes) | \
            scan (nx2)
        if (nx2 <= 750) { # 1-slit mode
            scopy (l_input//"["//l_extn//","//l_ver//"]",msjoin//".ms", renum+,
                merge+, clobber+, verbose-, rebin-)
        } else { # this is a pre-extracted 2-slit
            imcopy (l_input//"["//l_extn//","//l_ver//"]",msjoin//".ms", 
                verbose-)
            # Note: scopy won't work if nx2 > 1024
            #    but here the spectra have already been extracted so a simple
            #    imcopy will work just fine.
        }

    } else if (l_ver=="*") {
        if (!imaccess(l_input//"["//l_extn//",1]")) {
            print ("ERROR - GFDISPLAY: "//l_input//"["//l_extn//",1] does \
                not exist.")
            goto error
        }
        if (!imaccess(l_input//"["//l_extn//",2]")) {
            print ("ERROR - GFDISPLAY: "//l_input//"["//l_extn//",2] does \
                not exist.")
            goto error
        }
        imgets (l_input//"["//l_extn//",1]", "i_naxis1", >& "dev$null")
        nx1 = int(imgets.value)
        imgets (l_input//"["//l_extn//",1]", "refpix1", >& "dev$null")
        ref1 = real(imgets.value)
        imgets (l_input//"["//l_extn//",2]", "i_naxis1", >& "dev$null")
        nx2 = int(imgets.value)
        imgets (l_input//"["//l_extn//",2]", "refpix1", >& "dev$null")
        ref2 = real(imgets.value)

        len = min(nx1,nx2)
        refpix = min(ref1,ref2)
        sec1 = "["//(nint(ref1-refpix+1.))//":"//(nint(ref1-refpix+len))//",*]"
        sec2 = "["//(nint(ref2-refpix+1.))//":"//(nint(ref2-refpix+len))//",*]"
        #print(sec1, sec2)

        imcopy (l_input//"["//l_extn//",1]"//sec1,tmpmsa//".ms", verbose-)
        imcopy (l_input//"["//l_extn//",2]"//sec2,tmpmsb//".ms", verbose-)
        scopy (tmpmsa//".ms,"//tmpmsb//".ms",msjoin//".ms", renum+, merge+,
            clobber+, verbose-, rebin-)
    }
    count (l_config) | scan (nx2)
    imgets (msjoin//".ms.fits", "i_naxis2")
    if (nx2==int(imgets.value)) {
        ldisplay (msjoin//".ms", siz=1024, ldispdir="", lconf=l_config,
            lscale=0.2, deadfibs=l_dead, output=l_output, axisflip+,
            disqu=yes, display=cmd)
        if (!rmout) { # time stamps the output, if not tmp file
            gemdate ()
            gemhedit (l_output, "GFDISP", gemdate.outdate,
                "UT Time stamp for GFDISPLAY", delete-)
            gemhedit (l_output, "GEM-TLM", gemdate.outdate,
                "UT Last modification with GEMINI", delete-)
        }
    } else {
        print ("ERROR - GFDISPLAY: ldisplay configuration file does not \
            match image")
    }

    # clean up
    goto clean

error:
    status = 1
    goto clean

clean:
    imdelete (tmpmsa//".ms,"//tmpmsb//".ms,"//msjoin//".ms", verify-,
        >& "dev$null")
    if (rmout)
        imdelete (l_output, verify-, >& "dev$null")
    if (rmconf)
        delete (l_config//","//l_dead, verify-, >& "dev$null")

end

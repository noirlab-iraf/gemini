# Copyright(c) 2003-2011 Association of Universities for Research in Astronomy, Inc.

procedure tbackground(inimages)

#
# This routine automaticly finds bad frames in a T-ReCS raw image and writes
# the appropriate values to the extension headers.
#
# Version:  May  20, 2003:  KV routine created from BRs prototype
#           Sept  9, 2003;  KV changes to logging, made to work for all modes,
#                              works on a list of images
#           Sept 12, 2003:  KV added code to avoid duplicate badset numbers
#	    Oct   3, 2003:  KV small syntax changes to the parameters
#	    Oct  29, 2003:  KL IRAF 2.12 - new parameters
#			      imstat: nclip,lsigma,usigma,cache
#			      hedit: addonly
#           March 30, 2004:  KV bug fix (value of "m" can be undefined)
#                              I do not understand how this got past 
#                              the original testing of the script

char    inimages    {prompt="Input TReCS image(s)"}
char    outimages   {"",prompt="Output image(s)"}
char    outpref     {"b",prompt="Prefix for output image name(s)"}
char    rawpath     {"",prompt="Path for input raw images"}
real    sigma       {4.0,prompt="Sigma tolerance for bad frames"}
char    bsetfile    {"",prompt="Bad Frame list file"}
bool    writeps     {no,prompt="Write .ps file?"}
bool    sh_change   {no,prompt="Show changes to image headers?"}
char    logfile     {"",prompt="Logfile"}
bool    verbose     {yes,prompt="Verbose?"}
int     status      {0,prompt="Exit error status: (0=good, >0=bad)"}
struct  *scanfile   {"", prompt="Internal use only"}

begin

    char    l_image,l_logfile,l_bsetfile,l_inputimages
    char    l_filename,l_outputimages,l_prefix
    char    tmpfile,l_rawpath,tmplist
    char    in[100],out[100]
    char    in1,out1
    char    ref, fname, phu, shu, savesetstring
    char    tmpstata,tmpstatb, tmpmeana, tmpmeanb, tmpmeans
    char    tmpps,tmpbsets,tmpwork
    char    paramstr
    int     nnods, nsavesets, nnodsets, nextns, nframes
    int     i, j, k, l, m, n, maximages, nimages, modeflag, filetype
    int     noutputimages
    int     nbadsetorig,badsetorig[100],nbadset,badset[100]
    int     frmcoadd, chpcoadd, ncoadd, nref
    real    ADC_DARK,ADC_SAT
    real    pymin, pymax
    real    l_sigma
    bool    l_verbose, l_writeps, l_sh_change
    struct  l_struct


    tmpfile=mktemp("tmpfile")

    l_inputimages=inimages
    l_rawpath=rawpath
    l_outputimages=outimages
    l_sigma=sigma
    l_verbose=verbose
    l_writeps=writeps
    l_sh_change=sh_change
    l_logfile=logfile
    l_bsetfile=bsetfile
    l_prefix=outpref
    
    cache ("gemdate")

    status=0

    # Create the list of parameter/value pairs.  One pair per line.
    # All lines combined into one string.  Line delimiter is '\n'.
    paramstr =  "inimages       = "//inimages.p_value//"\n"
    paramstr += "outimages      = "//outimages.p_value//"\n"
    paramstr += "outpref        = "//outpref.p_value//"\n"
    paramstr += "rawpath        = "//rawpath.p_value//"\n"
    paramstr += "sigma          = "//sigma.p_value//"\n"
    paramstr += "bsetfile       = "//bsetfile.p_value//"\n"
    paramstr += "writeps        = "//writeps.p_value//"\n"
    paramstr += "sh_change      = "//sh_change.p_value//"\n"
    paramstr += "logfile        = "//logfile.p_value//"\n"
    paramstr += "verbose        = "//verbose.p_value

    # Assign a logfile name if not given.  Open logfile and start log.
    # Write parameter/value pairs ("paramstr") to log.
    gloginit( l_logfile, "tbackground", "midir", paramstr, fl_append+,
        verbose=l_verbose )
    if (gloginit.status != 0) {
        status = gloginit.status
        goto exitnow
    }
    l_logfile = gloginit.logfile.p_value

    # Check the rawpath name for a final /
    if (substr(l_rawpath,(strlen(l_rawpath)),(strlen(l_rawpath))) != "/")
        l_rawpath=l_rawpath//"/"
    if (l_rawpath=="/" || l_rawpath==" ")
        l_rawpath=""

    nimages=0
    maximages=100

    # Count the number of input images
    # First, generate the file list if needed

    if (stridx("*",l_inputimages) > 0) {
        files(l_rawpath//l_inputimages, > tmpfile)
        l_inputimages="@"//tmpfile
    }

    if (substr(l_inputimages,1,1)=="@") {
        scanfile=substr(l_inputimages,2,strlen(l_inputimages))
    } else {
        if (stridx(",",l_inputimages)==0) 
            files(l_inputimages, > tmpfile)
        else {
            j=9999
            while (j!=0) {
                j=stridx(",",l_inputimages)
                if (j>0)
                    files(substr(l_inputimages,1,j-1), >> tmpfile)
                else
                    files(l_inputimages, >> tmpfile)
                l_inputimages=substr(l_inputimages,j+1,strlen(l_inputimages))
            }
        }
        scanfile=tmpfile
    }

    i=0

    while (fscan(scanfile,l_filename) != EOF) {

        i=i+1

        if (substr(l_filename,strlen(l_filename)-4,strlen(l_filename))==".fits")
            l_filename=substr(l_filename,1,strlen(l_filename)-5)

        j=0
        if (stridx("/",l_rawpath) > 0 && stridx("/",l_filename) > 0) {
            j=stridx("/",l_filename)
            if (j > 0) {
                for (k=1; k < 100 && j > 0; k+=1) {
                    l_filename=substr(l_filename,j+1,strlen(l_filename))
                    j=stridx("/",l_filename)
                }
            }
        }

        if (!imaccess(l_rawpath//l_filename)) {
            glogprint( l_logfile, "tbackground", "status", type="error", 
                errno=101, str="Input image "//l_rawpath//l_filename//\
            " was not found.",verbose+)
        } else {
            nimages=nimages+1
            if (nimages > maximages) {
                glogprint( l_logfile, "tbackground", "status", type="error",
                    errno=121, str="Maximum number of input images \
                    ["//str(maximages)//"] has been exceeded.",verbose+)
                status=1
                goto exit
            }
            in[nimages]=l_rawpath//l_filename
            out[nimages]=l_filename
            j=stridx("/",out[nimages])
            if (j > 0) {
                for (k=1; k < 100 && j > 0; k+=1) {
                    out[nimages]=substr(out[nimages],j+1,strlen(out[nimages]))
                    j=stridx("/",out[nimages])
                }
            }
        }
    } #end of while-loop

    scanfile=""
    delete(tmpfile,ver-,>& "dev$null")

    if (nimages == 0) {
        glogprint( l_logfile, "tbackground", "status", type="error", errno=121,
            str="No input images were defined.",verbose+)
        status=1
        goto exit
    } else {
        glogprint( l_logfile, "tbackground", "status", type="string",
            str="Processing "//str(nimages)//" image(s).", verbose=l_verbose )
    }

    # Now, do the same counting for the output file
    tmpfile=mktemp("tmpfile")

    noutputimages=0
    if (l_outputimages != "" && l_outputimages != " ") {
        if (substr(l_outputimages,1,1) == "@") {
            scanfile=substr(l_outputimages,2,strlen(l_outputimages))
        } else {
            if (stridx("*",l_outputimages) > 0) {
                files(l_outputimages,sort-) | \
                    match(".hhd",stop+,print-,metach-, > tmpfile)
                scanfile=tmpfile
            } else {
                files(l_outputimages,sort-, > tmpfile)
                scanfile=tmpfile
            }
        }

        while (fscan(scanfile,l_filename) != EOF) {
            if (substr(l_filename,strlen(l_filename)-4,strlen(l_filename)) == ".fits")
                l_filename=substr(l_filename,1,strlen(l_filename)-5)

            noutputimages=noutputimages+1
            if (noutputimages > maximages) {
                glogprint( l_logfile, "tbackground", "status", type="error",
                    errno=121, str="Maximum number of output images "//\
                    str(maximages)//" exceeded.",verbose+)
                status=1
                goto exit
            }
            out[noutputimages]=l_filename
            if (imaccess(out[noutputimages])) {
                glogprint( l_logfile, "tbackground", "status", type="error",
                    errno=102, str="Output image "//l_filename//\
                    " already exists.",verbose+)
                status=1
                goto exit
            }
        }
        if (noutputimages != nimages) {
            glogprint( l_logfile, "tbackground", "status", type="error",
                errno=121, str="Different number of input ("//str(nimages)//") \
                and output ("//str(noutputimages)//" image names have been \
                specified.",verbose+)
            status=1
            goto exit
        }

        scanfile=""
        delete(tmpfile,ver-, >& "dev$null")
    } else {
        if (l_prefix == "" || l_prefix == " ")
            l_prefix="b"
        for (i=1; i <= nimages; i+=1) {
            out[i]=l_prefix//out[i]
            if (imaccess(out[i])) {
                glogprint( l_logfile, "tbackground", "status", type="error",
                    errno=102, str="Output image "//out[i]//" already exists.",
                    verbose+)
                status=1
                goto exit
            }
        }
    }

    if (l_sigma < 0.) {
        glogprint( l_logfile, "tbackground", "status", type="error", errno=121,
            str="The sigma parameter is negative.  Aborting.",verbose+) 
        status=1
        goto exit
    }

    if (l_sigma == 0.) {
        l_sigma=4.0
        glogprint( l_logfile, "tbackground", "science", type="warning",
            errno=0, str="The sigma parameter is zero.  Using the default \
            value of 4.0.",verbose+)
    }

    savesetstring=""

    #
    # The main loop: check each image in turn.
    #
    # Each image is copied from the "rawpath" directory to the local directory 
    # before the frames are examined.
    #

    l=1
    while (l <= nimages) {
        tmpstata=" "
        tmpmeana=" "
        tmpbsets=" "
        tmpwork=" "
        k=stridx(".",in[l])
        j=stridx(".",out[l])

        if (k == 0)
            in1=in[l]//".fits"
        else
            in1=in[l]

        if (j == 0)
            out1=out[l]//".fits"
        else
            out1=out[l]

        copy(in1,out1,verbose-)
        glogprint( l_logfile, "tbackground", "task", type="string",
            str="Copying image "//in[l]//" to "//out[l],verbose=l_verbose)
        l_image=out[l]

        #check if image exists
        if (!imaccess(l_image)) {
            glogprint( l_logfile, "tbackground", "status", type="error",
                errno=101, str="Copied image "//l_image//" is not found",
                verbose+)
            goto nextimage
        }

        glogprint( l_logfile, "tbackground", "status", type="string",
            str="Checking image "//l_image,verbose=l_verbose)
        glogprint( l_logfile, "tbackground", "status", type="fork",
            fork="forward", child="tcheckstructure", verbose=l_verbose )

        tcheckstructure(l_image,logfile=l_logfile,verbose=l_verbose)

        glogprint( l_logfile, "tbackground", "status", type="fork",
            fork="backward", child="tcheckstructure", verbose=l_verbose )

        j=tcheckstructure.status
        modeflag=tcheckstructure.modeflag
        if (modeflag > 4) {
            modeflag=modeflag-10
            filetype=2
        } else {
            filetype=1
        }

        if (j != 0) {
            glogprint( l_logfile, "tbackground", "status", type="error",
                errno=123, str="Image "//l_image//" does not have the \
                expected structure.",verbose+)
            status=1
            goto nextimage
        }
        if (filetype != 1) {
            glogprint( l_logfile, "tbackground", "status", type="error",
                errno=121, str="Image "//l_image//" appears to be a \
                'prepared' T-ReCS file.",verbose+)
            status=1
            goto nextimage
        }

        phu=l_image//"[0]"

        imgets(phu,"TBACKGRO",>& "dev$null")
        if (imgets.value != "0") {
            glogprint( l_logfile, "tbackground", "task", type="warning",
                errno=0, str="File "//in[i]//" has already been screened.",
                verbose=l_verbose)
            goto nextimage    
        }

        imgets(phu,"NNODS",>& "dev$null") ; nnods=int(imgets.value)
        imgets(phu,"NNODSETS",>& "dev$null") ; nnodsets=int(imgets.value)
        imgets(phu,"SAVESETS",>& "dev$null") ; nsavesets=int(imgets.value)
        imgets(phu,"NEXTEND",>& "dev$null") ; nextns=int(imgets.value)
        if (nextns == 0)
            nextns=nnods*nnodsets

        #check for 2 nod phases when nodding (can't handle anything else); 
        #assumes 2 chop phases
        if (nnods != 2) {
            glogprint( l_logfile, "tbackground", "status", type="error",
                errno=123, str="Image "//l_image//" has "//nnods//\
                " nod phases.",verbose+)
            goto nextimage
        }

        if (modeflag == 1) {
            nframes = nnodsets*nsavesets*4.
        } else {
            if (modeflag == 2 || modeflag == 3)
                nframes=nnodsets*nsavesets*2
            else
                nframes=nnodesets*nsavesets
        }

        imgets(phu,"frmcoadd",>& "dev$null") ; frmcoadd=int(imgets.value)
        imgets(phu,"chpcoadd",>& "dev$null") ; chpcoadd=int(imgets.value)
        ncoadd=frmcoadd*chpcoadd
        if (ncoadd <= 0) {
            glogprint( l_logfile, "tbackground", "status", type="error",
                errno=132, str="Header info missing or zero for FRMCOADD or \
                CHPCOADD",verbose+)
            status=1
            goto nextimage
        }

        #create tmpfiles
        tmpstata=mktemp("tmpstat"); tmpstatb=mktemp("tmpstat")
        tmpmeana=mktemp("tmpmean"); tmpmeanb=mktemp("tmpmean")
        tmpmeans=mktemp("tmpmeans")
        tmpbsets=mktemp("tmpbsets")
        tmpwork=mktemp("tmpwork")

        #CAN WE GET THESE FROM HEADER?? Or from TRECS package parameters
        # The following are just guesses at the moment...
        ADC_DARK=0.600e4
        ADC_SAT=65000.

        printf("%-25s %-3s %-3s %-3s %-3s %-3s\n","Image","n_savesets",
            "n_nodpos","n_nodsets","nextensions","total frames") | \
            scan(l_struct)
        glogprint( l_logfile, "tbackground", "science", type="string",
            str=l_struct, verbose=l_verbose )
        printf("%-25s   %-9d %-10d %-8d %-9d %-8d\n",l_image,nsavesets,nnods,
            nnodsets,nextns,nframes) | scan(l_struct)
        glogprint( l_logfile, "tbackground", "science", type="string",
            str=l_struct, verbose=l_verbose )

        # print to file all ref frames for each save & nod set: 
        # 	ref1=[nodA,chopB]; ref2=[nodB,chopA]
        # Save nodA and nodB to seperate lists
        #
        #FOLLOWING ASSUMES 1 NOD PHASE PER EXTENSION
        for (i=1;i<=nextns;i+=1) {
            shu=l_image//"["//str(i)//"]"
            for (j=1;j<101;j+=1) {
                badsetorig[j]=0
            }
            imgets(shu,"NBADSET",>& "dev$null")
            nbadsetorig=int(imgets.value)
            if (nbadsetorig > 0) {
                if (j < 10)
                    imgets(shu,"BADSET0"//str(j))
                else
                    imgets(shu,"BADSET"//str(j))
                badsetorig[j]=int(imgets.value)
            }
            if (modeflag == 1 || modeflag == 3) {
                imgets(shu,"NOD",>& "dev$null") 
                if (imgets.value == "A") {
                    ref = "2"; fname=tmpstata
                } else {
                    if (imgets.value == "B") {
                        ref = "1"; fname=tmpstatb
                    } else {
                        glogprint( l_logfile, "tbackground", "status", 
                            type="error", errno=132, str="NOD phase in \
                            extension header is not A or B",verbose=l_verbose)
                        goto nextimage
                    }
                }
            } else {
                if (modeflag == 4 || modeflag == 3)
                    ref=1
                else if (modeflag == 2)
                    ref=2
            }
            for (j=1;j<=nsavesets;j+=1) {
                print(l_image//"["//str(i)//"][*,*,"//ref//","//str(j)//"]",
                    >> fname)
            }
        }

        # stats do not exactly match f6bstat because f6bstat also demeans(?) 
        # result (subtracting 0th moment I think)
        # get means for each nod seperately for plotting
        imstat("@"//tmpstata,fields="midpt",lower=INDEF,upper=INDEF,nclip=0,
        lsigma=INDEF,usigma=INDEF,binwidth=0.1,format-,cache-,>tmpmeana)
        if (modeflag == 1 || modeflag == 3)
            imstat("@"//tmpstatb,fields="midpt",lower=INDEF,upper=INDEF,nclip=0,
                lsigma=INDEF,usigma=INDEF,binwidth=0.1,format-,cache-,>tmpmeanb)

        # change to %well: 
        tcalc(tmpmeana,"Well",
            "(c1/"//str(ncoadd)//"-"//str(ADC_DARK)//")/"//str(ADC_SAT-ADC_DARK),
            colfmt="f8.4")
        tcalc(tmpmeana,"Well","Well*100.",colfmt="f6.2")
        tcalc(tmpmeana,"Row","ROWNUM",colfmt="i5")
        if (modeflag == 1 || modeflag == 3) {
            tcalc(tmpmeanb,"Well",
                "(c1/"//str(ncoadd)//"-"//str(ADC_DARK)//")/"//str(ADC_SAT-ADC_DARK),
                colfmt="f8.4")
            tcalc(tmpmeanb,"Well","Well*100.",colfmt="f6.2")
            tcalc(tmpmeanb,"Row","ROWNUM",colfmt="i5")
            # join tables for statistics: 
            # final table = ref1(nodsets A) + ref2(nodsets B)
            tmerge(tmpmeana//","//tmpmeanb,tmpmeans,"append")
        } else
            tcopy(tmpmeana,tmpmeans)

        tstat(tmpmeans,"Well",outtable="STDOUT", >& "dev$null")
        glogprint( l_logfile, "tbackground", "science", type="string",
            str="Signal [percent full well] in reference frames:",
            verbose=l_verbose )
        glogprint( l_logfile, "tbackground", "visual", type="visual",
            vistype="empty", verbose=l_verbose )
        printf("                        Average  = %6.2f\n",tstat.mean) | \
            scan(l_struct)
        glogprint( l_logfile, "tbackground", "science", type="string",
            str=l_struct, verbose=l_verbose )
        printf("                        Stddev   = %8.4f\n",tstat.stddev) | \
            scan(l_struct)
        glogprint( l_logfile, "tbackground", "science", type="string",
            str=l_struct, verbose=l_verbose )
        printf("                        Minimum  = %6.2f\n",tstat.vmin) | \
            scan(l_struct)
        glogprint( l_logfile, "tbackground", "science", type="string",
            str=l_struct, verbose=l_verbose )
        printf("                        Maximum  = %6.2f\n",tstat.vmax) | \
            scan(l_struct)
        glogprint( l_logfile, "tbackground", "science", type="string",
            str=l_struct, verbose=l_verbose )

        #set x axis max value
        nref=nsavesets*nnodsets

        #set y axis values for plot to min/max or +/- 2 percentage points
        if (tstat.mean-2.0 < tstat.vmin)
            pymin=tstat.mean - 2.0
        else
            pymin=tstat.vmin - tstat.stddev

        if (tstat.mean+2.0 > tstat.vmax)
            pymax=tstat.mean + 2.0
        else
            pymax=tstat.vmax + tstat.stddev

        #plot results

        #get rid of "_" in image title
        printf(" %s Background\n",l_image) | scan(l_struct)

        #print("DEBUG: writing plot to file")
        if (l_writeps) {
            #plot to file first (takes awhile)
            mgograph (tmpmeana,3,2,rows="-",wx1=0,wx2=nref,wy1=pymin,wy2=pymax,
                excol="",eycol="",logx=no,logy=no,labelexp=1.5,boxexp=1.,
                xlabel="Reference Frame Number",ylabel="% of Full Well",
                title=l_struct, postitle="topcenter",append=no,pointmode=no,
                pattern="solid",crvstyle="straight",lweight=1,color=1,mkzero=no,
                device="psi_land",gkifile="mgo.gki")
            mgograph (tmpmeanb,3,2,append=yes,title="",pointmode=no,
                pattern="solid",crvstyle="straight",lweight=2,color=4,mkzero=no,
                device="psi_land",gkifile="mgo.gki")
            #add legend (position units are for normalized coords: 0-1 in x,y)
            igi (initcmd="DRELOCATE 0.35 0.85 LWEIGHT 1 COLOR 1 DDRAW 0.42 0.85 PUTLABEL 6 NodA ;DRELOCATE 0.6 0.85 LWEIGHT 2 COLOR 4 DDRAW 0.67 0.85 PUTLABEL 6 NodB; END",
                append=yes,device="psi_land", >&"dev$null")
            gflush
            #print("DEBUG: writing plot to screen")
            mgograph(tmpmeana,3,2,rows="-",wx1=0,wx2=nref,wy1=pymin,wy2=pymax,
                excol="",eycol="",logx=no,logy=no,labelexp=1.5,boxexp=1.,
                xlabel="Reference Frame Number",ylabel="% of Full Well",
                title=l_struct, postitle="topcenter",append=no,pointmode=no,
                pattern="solid",crvstyle="straight",lweight=1,color=1,
                mkzero=no,device="stdgraph",gkifile="mgo.gki")
            mgograph(tmpmeanb,3,2,append=yes,title="",pointmode=no,
                pattern="solid",crvstyle="straight",lweight=2,color=4,
                mkzero=no,device="stdgraph",gkifile="mgo.gki")
            #add legend
            igi(initcmd="DRELOCATE 0.35 0.85 LWEIGHT 1 COLOR 1 DDRAW 0.42 0.85 PUTLABEL 6 NodA ;DRELOCATE 0.6 0.85 LWEIGHT 2 COLOR 4 DDRAW 0.67 0.85 PUTLABEL 6 NodB; END",
                append=yes,device="stdgraph", >&"dev$null")
        }

        #print("DEBUG: Finding bad savesets")
        #Finding bad savesets: median values l_sigma(default 4) below/above 
        #average
        if ((tstat.vmin < (tstat.mean-tstat.stddev*l_sigma)) || (tstat.vmax > (tstat.mean+tstat.stddev*l_sigma))) {
            glogprint( l_logfile, "tbackground", "status", type="warning",
                errno=102, str="Bad frames (median outside "//l_sigma//\
                " sigma) exist.",verbose+)
            if ((l_bsetfile=="") || (l_bsetfile==" ")) {
                #if bad frame file not defined, default to image_name.bsets
                l_bsetfile=l_image//".bsets"
            }
            #check if file exists
            if (access(l_bsetfile))
                glogprint( l_logfile, "tbackground", "task", type="warning", 
                    errno=0, str="Appending to existing Bad frames \
                    file: "//l_bsetfile,verbose+)

            printf("\n# File: "//l_image//"  -  Save and Nod sets outside "//\
                l_sigma//" sigma of mean\n", >>l_bsetfile)
            tselect(tmpmeans,tmpbsets,"Well > ("//str(tstat.mean)//"+"//\
                str(tstat.stddev)//"*"//str(l_sigma)//") || Well < ("//\
                str(tstat.mean)//"-"//str(tstat.stddev)//"*"//str(l_sigma)//")")

            ######
            #NOTE: calculation of nodset saveset assumes all ref1's followed by ref2's:
            ######
            tcalc(tmpbsets,"Nodset","int((Row-1)/"//str(nsavesets)//")+1",
                datatype="int",colfmt="i5")
            tcalc(tmpbsets,"Saveset","Row-((Nodset-1)*"//str(nsavesets)//")",
                datatype="int",colfmt="i5")   
            tprint(tmpbsets,prdata+,showrow-,showhdr+,align+,>>l_bsetfile)
            glogprint( l_logfile, "tbackground", "task", type="string",
                str="Savesets outside "//l_sigma//" sigma of mean written \
                to "//l_bsetfile,verbose=l_verbose) 

            # Now, write the header parameters that mark the bad savesets in 
            # each nodset

            for (k=1;k<= nnodsets;k+=1) {
                if (access(tmpwork))
                    delete(tmpwork,verify-, >&"dev$null")    
                tselect(tmpbsets,tmpwork,"Nodset == "//str(k))
                tinfo(tmpwork,ttout=no)
                nbadset=int(tinfo.nrows)
                m=0
                if (nbadset > 0) {
                    tmplist=mktemp("tmplist")
                    tdump (tmpwork,cdfile="",pfile="",datafile=tmplist,
                        columns="Saveset")
                    scanfile=tmplist
                    for (i=1;i <= nbadset && fscan(scanfile,savesetstring) != EOF;i+=1) {
                        badset[i]=int(savesetstring)
                    }
                    m=nbadsetorig
                    n=0
                    for (i=1;i <= nbadset ; i+=1) {
                        for (j=1; j <= m && n == 0; j+=1) {
                            if (badsetorig[j] == badset[i])
                                n=1
                        }
                        if (n == 0) {
                            m=m+1
                            badsetorig[m]=badset[i]
                        }
                    }
                    delete(tmplist,ver-,>& "dev$null")
                }
                if (m > nbadsetorig) {
                    gemhedit (l_image//"["//str(2*k-1)//"]", "NBADSET", m,
                        "", delete-)
                    if (modeflag == 1) {
                        gemhedit (l_image//"["//str(2*k)//"]", "NBADSET", m,
                            "", delete-)
                    }
                    for (i=1;i <= m;i+=1) {
                        if (i < 10) {
                            gemhedit (l_image//"["//str(2*k-1)//"]",
                                "BADSET0"//str(i), str(badsetorig[i]), "",
                                delete-)
                            if (modeflag == 1)
                                gemhedit (l_image//"["//str(2*k)//"]",
                                    "BADSET0"//str(i), str(badsetorig[i]),
                                    "", delete-)
                        } else {
                            gemhedit (l_image//"["//str(2*k-1)//"]",
                                "BADSET"//str(i), str(badsetorig[i]),
                                "", delete-)
                            if (modeflag == 1)
                                gemhedit (l_image//"["//str(2*k)//"]",
                                    "BADSET"//str(i), str(badsetorig[i]),
                                    "", delete-)
                        }
                    }
                }
            } # end for-loop over nodsets
        } else {
            glogprint( l_logfile, "tbackground", "science", type="string",
                str="There are no frames with median outside "//l_sigma//\
                " sigma of mean.\n",verbose=l_verbose)
        }

        if (l_writeps) {
            # Move the PS image from temporary file
            tmpps=""
            files("ps*.eps",sort-) | scan(tmpps)
            if (tmpps!="") {
                if (access(l_image//"_ref.ps")) {
                    glogprint( l_logfile, "tbackground", "status", 
                        type="warning", errno=102, str="Overwriting \
                        previous .ps output image", verbose+)
                    delete(l_image//"_ref.ps",verify-, >&"dev$null")
                }
                rename(tmpps,l_image//"_ref.ps", field="all")
                glogprint( l_logfile, "tbackground", "task", type="string",
                    str="Postscript file of median sky level in ref \
                    frames: "//l_image//"_ref.ps",verbose=l_verbose)
            } else {
                glogprint( l_logfile, "tbackground", "status", type="error",
                    errno=101, str="Cannot find .ps output image",verbose+)
            }
        }

        # Time stamp the primary header
        #
        gemdate ()
        gemhedit (l_image//"[0]", "TBACKGRO", gemdate.outdate,
            "UT Time stamp for TBACKGROUND", delete-)
        gemhedit (l_image//"[0]", "GEM-TLM", gemdate.outdate,
            "UT Last modification with GEMINI", delete-)

        # jump to here if there is a problem
nextimage:

        #delete temporary files
        if (tmpstata != " ") 
            delete(tmpstata//","//tmpstatb,verify-, >& "dev$null")
        if (tmpmeana != " ") 
            delete(tmpmeana//","//tmpmeans//","//tmpmeanb,verify-,>& "dev$null")
        delete(tmpbsets, verify-, >& "dev$null")
        delete(tmpwork, verify-, >& "dev$null")

        l=l+1
    }

exit:
    delete(tmpfile,verify-, >&"dev$null")

    if (status == 0)
        glogclose( l_logfile, "tbackground", fl_success+, verbose=l_verbose )
    else
        glogclose( l_logfile, "tbackground", fl_success-, verbose=l_verbose )

exitnow:
    ;

end

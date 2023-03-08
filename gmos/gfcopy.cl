# Copyright(c) 2002-2004 Association of Universities for Research in Astronomy, Inc.

procedure gfcopy (inimages)

# Copy spectra into a new datacube
# 
# Version   2001/12/01  bmiller   created
#           Aug 25, 2003  KL  IRAF2.12 - new parameter, addonly, in hedit

string inimages     {"",prompt="Input images"}
string outimages    {"",prompt="Output images"}
string outprefix    {"p",prompt="Prefix for output images"}
string expr         {"XINST < 10.",prompt="Expression for selection of sky spectra"}
string logfile      {"gmos.log",prompt="Logfile name"}
bool   verbose      {yes,prompt="Verbose?"}
int    status       {0,prompt="Exit status (0=good)"}
struct *flist       {"",prompt="Internal use only"}
struct *slist       {"",prompt="Internal use only"}

begin

    # Local variables for input parameters
    string l_inimages,l_outimages,l_prefix,l_logfile
    string l_expr
    bool   l_verbose

    # Other local variables
    string inlist,outlist,temp1,temp2,temp3,mdf
    string img,outimg,suf,aplist,range,sdum,sky
    int    i,j,nbad,nin,nout,nsky,prev,rstart,nrange,dum,nextnd
    bool   useprefix
    struct sdate,sline

    # Initialize exit status
    status=0

    # cache some parameter files
    cache("imgets","gmos","hedit","gimverify","tinfo","tabpar")

    # Initialize local variables
    l_inimages=inimages ; l_outimages=outimages ; l_prefix=outprefix
    l_logfile=logfile ; l_verbose=verbose
    l_expr=expr

    # Start logfile
    if((l_logfile=="") || (l_logfile==" ")) {
	l_logfile="gmos.log"
	printlog("WARNING - GFSKYSUB:  logfile is empty, using gmos.log.",
	    l_logfile,l_verbose) 
    }

    printlog("-------------------------------------------------------------------------------",
        l_logfile,l_verbose)
    date | scan(sdate)
    printlog("GFSKYSUB -- "//sdate,l_logfile,l_verbose)
    printlog(" ",l_logfile,l_verbose)

    #The usual kind of checks to make sure everything we need is specified. 

    nbad=0

    if (l_inimages=="" || l_inimages==" ") {
        printlog("ERROR - GFSKYSUB: Input spectra is an empty string",
            l_logfile,yes)
        nbad=nbad+1 
    }

    if (substr(l_inimages,1,1)=="@") {
        inlist=substr(l_inimages,2,strlen(l_inimages))
        if (!access(inlist)) {
	    printlog("ERROR - GFSKYSUB: Input list "//inlist//" not found",
	        l_logfile,verbose+)
	    nbad=nbad+1
        }
    }

    if ((l_outimages=="" || l_outimages==" ") && (l_prefix=="" || l_prefix==" ")) {
	printlog("ERROR - GFSKYSUB: Neither the output spectra nor prefix is specified.",
	    l_logfile,yes)
	nbad=nbad+1
    } else if ((l_outimages!="" && l_outimages!=" ")) {
	useprefix=no
    } else {
	useprefix=yes
    }

    if (substr(l_outimages,1,1)=="@") {
        outlist=substr(l_outimages,2,strlen(l_outimages))
        if (!access(outlist)) {
	    printlog("ERROR - GFSKYSUB: Output list "//outlist//" not found",
	        l_logfile,yes)
	    nbad=nbad+1
        }
    }

    temp1=mktemp("tmpin")
    temp2=mktemp("tmpout") 
    temp3=mktemp("tmpfilelist") 

    files(l_inimages,sort-, > temp1)
    count(temp1) | scan(nin)
    if (!useprefix) { 
	files(l_outimages,sort-, > temp2)
	count(temp2) | scan(nout)
	if (nin != nout) {
            printlog("ERROR - GFSKYSUB: Different number of input and output spectra",
	        l_logfile,yes)
            nbad=nbad+1
	}
    } else {
       files(l_prefix//"//@"//temp1,sort-, > temp2)
    }

    flist=temp1
    while (fscan(flist,img) != EOF) { 
	gimverify(img)
	if (gimverify.status>0)
            nbad=nbad+1
     
    } #end of while loop over images.
    flist=""

    #If anything was wrong then exit. 
    if (nbad > 0) { 
	printlog("ERROR - GFSKYSUB: "//nbad//" errors found with input parameters. Exiting.",
	    l_logfile,yes) 
	goto error
    } 

    #If we are here then everything should be OK. 
    #Write all the relevant info to the logfile:
    #
    printlog("inimages = "//l_inimages,l_logfile,l_verbose) 
    printlog("outimages = "//l_outimages,l_logfile,l_verbose) 
    printlog("outprefix = "//l_prefix,l_logfile,l_verbose) 
    printlog(" ",l_logfile,l_verbose)

    # Make list file for main loop
    joinlines (temp1//","//temp2, output=temp3, delim=" ", missing="Missing", \
        maxchar=161, shortest=yes, verbose=no)
    flist=temp3 
    while (fscan(flist,img,outimg) != EOF) {
        # create tmp image/FITS file names used only within this loop
	mdf=mktemp("tmpmdf")
	sky=mktemp("tmpsky")
    
	# check .fits
	suf=substr(img,strlen(img)-4,strlen(img))
	if (suf!=".fits")
	  img=img//".fits"

	suf=substr(outimg,strlen(outimg)-4,strlen(outimg))
	if (suf!=".fits")
	  outimg=outimg//".fits"

	print(img," ",outimg)

	# copy MDF
	tselect(img//"[MDF]",mdf//".fits","BEAM >= 0")
	tcalc(mdf//".fits","LINE","rownum",datatype="int")

	# Select the spectra
	tcalc(mdf//".fits",
	    "BEAM","if BEAM != -1 && "//l_expr//" then 1 else -1")
	tselect(mdf//".fits",mdf//"sel.fits","BEAM > 0")
	wmef(mdf//"sel.fits",outimg,extnames="MDF",verb-,phu=img, >& "dev$null")
	delete(mdf//"*.fits",verify-)

	aplist=""
	tinfo(outimg//"[MDF]",ttout-)
	nsky=0
	prev=0
	range=""
	nrange=0
	for (i=1; i<=tinfo.nrows; i+=1) {
	    tabpar(outimg//"[MDF]","LINE",i)
	    j=int(tabpar.value)
	    if (prev==0) {
		prev=j
		rstart=j
		nrange=j
	    } else if (j-prev != 1 || i==tinfo.nrows) {
		range=str(rstart)
		if (nrange>1)
		    range=range//"-"//str(prev)
		if (aplist=="")
		    aplist=range
		else
		    aplist=aplist//","//range
		prev=j
		rstart=j
		range=""
		nrange=1
	    } else {
		nrange=nrange+1
		prev=j
	    }
	    nsky=nsky+1
	}
	printlog("Found "//nsky//" spectra",l_logfile,l_verbose)
	print (aplist)

	# Make the combined sky
	scopy(img//"[SCI,1]",sky,format="multispec",apertures=aplist,bands="",
	    beams="",renumber+,apmodulus=0,rebin+,verbose=l_verbose)


	# Put image in output
	imgets(outimg//"[0]","NEXTEND")
	nextnd=int(imgets.value)
	fxinsert(sky//".fits",l_outimage//"["//nextnd//"]","0",verbose-,
	    >& "dev$null")
	hedit(outimg//"[0]","NEXTEND",(nextnd+1),add-,addonly-,delete-,verify-,
	    show-,update+)

	delete(sky//".fits",verify-, >>& "dev$null")
    }
    flist=""

    goto clean

error:
    status=1

clean:
    delete(temp1//","//temp2//","//temp3,verify-, >& "dev$null") 
    printlog("",l_logfile,l_verbose)
    flist=""
    printlog("GFSKYSUB done",l_logfile,l_verbose)
    printlog("-------------------------------------------------------------------------------",
        l_logfile,l_verbose)

end

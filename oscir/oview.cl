# Copyright(c) 2001-2006 Association of Universities for Research in Astronomy, Inc.

procedure oview (image)

#View frames of raw OSCIR data according to specified type(default=sig);
# allows interactive (imexam) mode for sig frames only (too slow for
# other frame types which view generally hundreds of images
# Images are displayed in ximtool; nothing is written to logfile 
#
#NOTES: Need to have ximtool already running;
#	works on single image only
#
# Version: Sept 14, 2002 BR  Release v1.4

char image	{"",prompt="Image to display"}
char type	{"sig",prompt="ref1|ref2|src1|src2|dif1|dif2|sig|sigall"}
real delay	{0,prompt="update delay in seconds"}
bool manual	{no,prompt="Manual mode? (only for type=sig)"}
bool disp_fill	{no,prompt="Fill display?"}
bool test_mode	{no,prompt="Test mode?"}
bool verbose	{yes,prompt="Verbose?"}

begin

    char l_image,l_type
    real l_delay
    bool l_verbose,l_dispfill,l_testmode, l_manual

    char fname,tref1,tref2,tsrc1,tsrc2,tdif1,tdif2
    char tnod,tsig,tmpfile,tmpimg
    int  n_savesets,n_nodsets
    int  n_i,n_j,naxis3,naxis5

    l_image=image ; l_type=type ; l_delay=delay
    l_verbose=verbose; l_dispfill=disp_fill; 
    l_testmode=test_mode; l_manual=manual

    if (l_testmode) time

    #check if image exists
    if (!imaccess(l_image)) {
	print("ERROR - OVIEW: Image "//l_image//" not found")
	bye
    }

    #load needed packages, cache frequently used tasks
    cache("imgets","display")

    #get axes
    imgets(l_image,"i_naxis3") ; naxis3=int(imgets.value)
    imgets(l_image,"i_naxis4") ; n_savesets=int(imgets.value)
    imgets(l_image,"i_naxis5") ; naxis5=int(imgets.value)
    imgets(l_image,"i_naxis6") ; n_nodsets=int(imgets.value)

    #check if image contains chop-nod data
    if ((naxis3!=2)||(naxis5!=2)) {
	print("ERROR - OVIEW: Image "//l_image//" is not chop-nod data;")
	print("                     n_choppos= "//naxis3//", n_nodpos= "//naxis5)
	bye
    }

    if ((l_type != "sig") && (l_type != "sigall") && (l_type != "ref1") &&
    	(l_type != "ref2") && (l_type!="dif1") && (l_type!="dif2") &&
	(l_type!="src1")&&(l_type!="src2")) {
	print("ERROR - OVIEW: Image type invalid.  lpar oview for valid \
	    values.")
	bye
    }

    printf("oview: Image="//l_image//"   Type="//l_type)
    printf("         Savesets="//n_savesets//"   Nodsets="//n_nodsets//"\n")

    tref1=mktemp("tmpref1")
    tref2=mktemp("tmpref2")
    tsrc1=mktemp("tmpsrc1")
    tsrc2=mktemp("tmpsrc2")
    tdif1=mktemp("tmpdif1")
    tdif2=mktemp("tmpdif2") 
    tsig=mktemp("tmpsig")
    tnod=mktemp("tmpnod")
    tmpfile=mktemp("tmpfile")

    #create temporary images for each type (time-consuming step)
    printf("\n Extracting frames, this may take a moment...\n")

    imcopy(l_image//"[*,*,1,*,1,*]",tsrc1, >>& "dev$null")
    imcopy(l_image//"[*,*,2,*,1,*]",tref1, >>& "dev$null")
    imcopy(l_image//"[*,*,2,*,2,*]",tsrc2, >>& "dev$null")
    imcopy(l_image//"[*,*,1,*,2,*]",tref2, >>& "dev$null")
    if ((l_type == "sig")||(l_type=="sigall")||(l_type=="dif1")||(l_type=="dif2")) {
	imarith(tsrc1,"-",tref1,tdif1, >>& "dev$null")
	imarith(tsrc2,"-",tref2,tdif2, >>& "dev$null")
    }
    if ((l_type=="sigall")||(l_type=="sig"))
	imarith(tdif1,"+",tdif2,tsig, >>&"dev$null")

    if (l_testmode) time

    if (l_type == "sig") {
	imcopy(tsig//"[*,*,1,*]",tnod, >>& "dev$null")		#dummy 3-D file
	imarith(tnod,"*",0.0,tnod, >>& "dev$null")
	#sum dif files for each nodset
	for (n_j=1;n_j<=n_nodsets;n_j+=1) {
	    # Create tmp FITS file names used within this loop
	    tmpimg=mktemp("tmpimg")
	    
	    #if (l_testmode) print(l_type//": Nodset "//str(n_j))
	    for (n_i=1;n_i<=n_savesets;n_i+=1)
	  	print(tsig//"[*,*,"//str(n_i)//","//str(n_j)//"]", >> tmpfile)
	    imcomb("@"//tmpfile,tmpimg,combine="average",reject="none",
	        logfile="",outtype="double", >>& "dev$null")
	    imcopy(tmpimg,tnod//"[*,*,"//str(n_j)//"]", >>& "dev$null")
	    delete(tmpimg//".fits,"//tmpfile,verify-, >>&"dev$null")
	}

	if (l_testmode) time

	#display sig frames
	for (n_j=1;n_j<=n_nodsets;n_j+=1) {
	    if (l_verbose) print(l_type//": Nodset "//str(n_j))
		display(tnod//"[*,*,"//str(n_j)//"]",1,erase-,fill=l_dispfill,
	        >>& "dev$null")
		if (l_manual)
		    imexamine() 
		else
		    sleep(l_delay)
	}
    } else {
	#set filename for other types
	if (l_type == "src1")
	    fname = tsrc1
	else if (l_type == "src2")
	    fname = tsrc2
	else if (l_type == "ref1")
	    fname = tref1
	else if (l_type == "ref2")
	    fname = tref2
	else if (l_type == "dif1")
	    fname = tdif1
	else if (l_type == "dif2")
	    fname = tdif2
	else
	    fname = tsig

	#display 
	for (n_j=1;n_j<=n_nodsets;n_j+=1) {
	    for (n_i=1;n_i<=n_savesets;n_i+=1) {
		if (l_verbose)
		    print(l_type//": Saveset "//str(n_i)//"  Nodset "//str(n_j))
		display(fname//"[*,*,"//str(n_i)//","//str(n_j)//"]",1,
		    erase-,fill=l_dispfill, >>& "dev$null")
		sleep(l_delay)
	    }
	}
    }

    if (l_testmode) time

    print("Done.")

    #cleanup
    imdelete(tref1//","//tref2//","//tsrc1//","//tsrc2//","//tnod,verify-,
        >>& "dev$null")
    imdelete(tdif1//","//tdif2//","//tsig,verify-, >>& "dev$null")

end



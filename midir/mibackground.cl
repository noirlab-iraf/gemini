# Copyright(c) 2001-2004 Association of Universities for Research in Astronomy, Inc.

procedure miview(inimages)

#View frames of raw Michelle data: src, ref or diff (default=dif);
# interactive mode: step through frames, eventually will be able to
# select frames for removal from stack
# non-interactive: just display frames with specified delay
#
# Images are displayed in ximtool; nothing is written to logfile 
#
#NOTES: in current form, need to have ximtool already running;
#       works on single image only
#
# Version 06 Feb 2003 - BR
#
# History:
# 15 Nov 2001,BR- modeled after IDL UF f6movie task; 
#                 started work 22 Oct01 
# 25 Jul 2002, BR - copied from oview, modified to handle MEF
# 05 Feb 2003, BR - major rewrite for Michelle data (will work for T-ReCS too, after 'tprepare' is written)
# 21 Jun 2003, TB - edit to work with real michelle data - works.
# 05 Oct 2003  TB - edit to work with michelle data only.
# 29 Oct 2003  KL - IRAF 2.12 - new parameters
#                     imstat: nclip,lsigma,usigma,cache (no changes required)
#                     hedit: addonly
# 3  Dec 2003  KV - removed a 'delete("tmp*")' statement, went back 
#                     to specific deletes of tmp working files.
#

char  inimages     {"",prompt="Images to display"}
char  outimages    {"",prompt="Output images(s)"}               
char  outprefix    {"v",prompt="Prefix for out images(s)"} 
char  rawpath         {"",prompt="Path for input raw images"}
char type       {"dif",prompt="src|ref|dif|sig"}
real delay      {0,prompt="update delay in seconds"}
bool fl_inter   {no,prompt="Run interactively?"}
bool fl_disp_fill       {no,prompt="Fill display?"}
bool fl_test_mode       {yes,prompt="Test mode?"}
bool fl_verbose {yes,prompt="Verbose?"}
char  logfile      {"",prompt="Logfile"}                        
real  z1           {0.,prompt="Minimum level to be displayed"}
real  z2           {0., prompt="Maximum level to be displayed"}
bool   zscale   {yes,prompt="Auto set grayscale display range"}
bool   zrange   {yes,prompt="Auto set image intensity range"}
char   ztrans   {"linear",prompt="Greyscale transformation"}
struct* scanfile   {"",prompt="Internal use only"} 

begin

char l_inimages,l_type,cursinput,header, l_logfile, l_outimages,l_rawpath
char l_prefix, instrument, l_temp, tmphead,lztrans,check
char l_filename,in[100],out[100],tmpfile, tmpfile1,tmpfile2,tmpfile3,tmpin
real l_delay, lz1, lz2

bool l_verbose,l_dispfill,l_testmode, l_manual,lzscale, lzrange
int  n_nodsets,i,nextns,naxis3,n_a,n_b,status,nbad,nimages, noutimages
int  n_i,n_j,ext,ndim,wcs,nbadsets,npars,maximages,modeflag
real xpos,ypos
struct l_struct

tmpfile=mktemp("tmpfile")
tmpin=mktemp("tmpin")
tmphead=mktemp("tmphead")

l_inimages=inimages ; l_type=type ; l_delay=delay ;
l_verbose=fl_verbose ; l_dispfill=fl_disp_fill ; 
l_testmode=fl_test_mode ; l_manual=fl_inter
cursinput=""
l_logfile=logfile
l_outimages=outimages
l_rawpath=rawpath
l_prefix=outprefix
lz1=z1;lz2=z2
lztrans=ztrans; lzscale=zscale; lzrange=zrange

if (l_testmode) time

nimages=0
status=0
maximages=100

if((l_logfile=="") || (l_logfile==" ")) {
   l_logfile=midir.logfile
}

# Load up arrays of input name lists
# This version handles both *s, / and commas in l_inimages

# Check the rawpath name for a final /
if(substr(l_rawpath,(strlen(l_rawpath)),(strlen(l_rawpath))) != "/") {
   l_rawpath=l_rawpath//"/"
}
if(l_rawpath=="/" || l_rawpath==" ")
  l_rawpath=""

# check that list file exists
if(substr(l_inimages,1,1)=="@") {
  l_temp=substr(l_inimages,2,strlen(l_inimages))
  if(!access(l_temp) || !access(l_rawpath//l_temp)) {
    printlog("ERROR - MIVIEW:  Input file "//l_temp//" not found.",l_logfile,l_verbose)
    status=1
    goto clean
  }
}


# Count the number of in images
# First, generate the file list if needed


if (stridx("*",l_inimages) > 0) {
  files(l_inimages, > tmpfile)
  l_inimages="@"//tmpfile
}

if (substr(l_inimages,1,1)=="@") {
  scanfile=substr(l_inimages,2,strlen(l_inimages))
}
else {
  files(l_inimages,sort-, > tmpfile)
  scanfile=tmpfile
}

i=0


nimages=0
nbad=0
noutimages=0
npars=0

while (fscan(scanfile,l_filename) != EOF && i <= 100) {

  i=i+1

  if (substr(l_filename,strlen(l_filename)-4,strlen(l_filename)) == ".fits"){
    l_filename=substr(l_filename,1,strlen(l_filename)-5)
#    l_filename=l_filename//".fits"
  }

  l_temp=substr(l_inimages,2,strlen(l_inimages))
  if(!access(l_filename//".fits") && !access(l_rawpath//l_filename//".fits")) {
    printlog("ERROR - MIVIEW: Input image "//l_filename//" was not found.",l_logfile, l_verbose)
    status=1
    goto clean
  }  else {
    nimages=nimages+1
    if (nimages > maximages) {
      printlog("ERROR - MIVIEW: Maximum number of input images exceeded:"//maximages,l_logfile,l_verbose)
      status=1 
      goto clean
    }
}
if(l_rawpath=="" || l_rawpath==" "){
    in[nimages]=l_filename
} else {
    in[nimages]=l_rawpath//l_filename
  }
}

scanfile=""

delete(tmpfile//","//tmpin,ver-,>& "dev$null")

if (nimages == 0) {
  printlog("ERROR - MIVIEW: No input images defined.",l_logfile, verbose+)
  status=1
  goto clean
}

# If prefix is to be used instead of filename
if (l_manual) {
if(l_outimages=="" || l_outimages==" ") {
  print(l_prefix) | scan(l_prefix)
  if(l_prefix=="" || l_prefix==" ") {
     printlog("ERROR - MIVIEW: Neither output images name nor output prefix is defined.",l_logfile,verbose+)
     status=1
     goto clean
  }

if(l_outimages=="" || l_outimages==" ") {
noutimages = nimages
}
}

# Now, do the same counting for the out file

tmpfile=mktemp("tmpfile")

if (l_outimages != "") {
  if (substr(l_outimages,1,1) == "@") {
    scanfile=substr(l_outimages,2,strlen(l_outimages))
  }
  else {
    if (stridx("*",l_outimages) > 0) {
      files(l_outimages,sort-) | match(".hhd",stop+,print-,metach-, > tmpfile)
      scanfile=tmpfile
    }
    else {
      files(l_outimages,sort-, > tmpfile)
      scanfile=tmpfile
    }
  }


  while(fscan(scanfile,l_filename) != EOF) {
    noutimages=noutimages+1
    if (noutimages > maximages) {
      printlog("ERROR - MIVIEW: Maximum number of output images exceeded:"//maximages,l_logfile, verbose+)
      status=1
      goto clean
    }
    if (l_manual) {
  if (substr(l_filename,strlen(l_filename)-4,strlen(l_filename)) != ".fits"){
    out[noutimages]=l_filename//".fits"
  } else {
    out[noutimages]=l_filename
  }
    if (imaccess(out[noutimages])) {
      printlog("ERROR - MIVIEW: Output image "//l_filename//" already exists.",l_logfile,l_verbose)
      status = 1
      goto clean
    }
    }
  }

  if (noutimages != nimages) {
    printlog("ERROR - MIVIEW: Different number of in images ("//nimages//") and out images ("//noutimages//")",l_logfile, verbose+)
    status =1
    goto clean
  }
}
  scanfile=""
  delete(tmpfile,ver-, >& "dev$null")
}

  i=1

  while(i<=nimages) {

if(l_outimages=="" || l_outimages==" ") {
   if (l_manual) {
  if (substr(in[i],strlen(in[i])-4,strlen(in[i])) == ".fits"){
    out[i]=l_prefix//in[i]
  } else {
    out[i]=l_prefix//in[i]//".fits"
  }
    if(imaccess(out[i])) {
       printlog("ERROR - MIVIEW: Output image "//out[i]//" already exists.",
          l_logfile,verbose+)
       nbad+=1
    }
  }
  if (nbad > 0) {
    printlog("ERROR - MIVIEW: "//nbad//" image(s) already exist.",l_logfile,verbose+)
    status=1
    goto clean
  }
}
    i=i+1
}

  i=1

  while(i<=nimages) {



#check instrument
  header=in[i]//"[0]"

  imgets(header,"INSTRUMENT")
  instrument=imgets.value
  printlog("Instrument is:"//instrument,logfile=l_logfile, verbose+)

    imgets(header,"MISTACK", >& "dev$null")
  check=imgets.value

   if (check == "0") {

  if (instrument == "michelle") {
    imgets(header,"MPREPARE", >& "dev$null")
    if (imgets.value == "0") {
      printlog("ERROR - MIVIEW: Image "//in[i]//" not MPREPAREd.  Use MVIEW to display.",l_logfile,l_verbose)
      status=1
      goto clean
    } else {
    print("MIVIEW: Displaying michelle images.")
      imgets(header,"MODE", >& "dev$null")
      if (imgets.value == "0") {
        printlog("ERROR - MIVIEW: could not find the MODE from the primary header.",l_logfile,l_verbose)
        status=status+1
        goto clean
      }
      modeflag=0
      if (imgets.value == "chop-nod") modeflag=1
      if (imgets.value == "ndchop") modeflag=1
      if (imgets.value == "chop") modeflag=2
      if (imgets.value == "nod") modeflag=3
      if (imgets.value == "ndstare") modeflag=4
      if (imgets.value == "stare") modeflag=4
      if (modeflag == 0) {
        printlog("ERROR - MIVIEW: Unrecognised MODE ("//imgets.value//") in the primary header.",l_logfile,l_verbose)
        status=status+1
        goto clean
      }
}
}

  if (instrument == "TReCS") {
    imgets(header,"TPREPARE", >& "dev$null")
    if (imgets.value == "0") {
      printlog("ERROR - MIVIEW: Image "//in[i]//" not TPREPAREd.  Use MVIEW to display.",l_logfile,l_verbose)
      status=1
      goto clean
    } else {
    print("MIVIEW: Displaying TReCS images.")
      imgets(header,"OBSMODE", >& "dev$null")
      if (imgets.value == "0") {
        printlog("ERROR - MIVIEW: could not find the OBSMODE from the primary header.",l_logfile,l_verbose)
        status=status+1
        goto nextimage
      }
      modeflag=0
      if (imgets.value == "chop-nod") modeflag=1
      if (imgets.value == "chop") modeflag=2
      if (imgets.value == "nod") modeflag=3
      if (imgets.value == "stare") modeflag=4
      if (modeflag == 0) {
        printlog("ERROR - MIVIEW: Unrecognised MODE ("//imgets.value//") in the primary header.",l_logfile,l_verbose)
        status=status+1
        goto clean
      }
}
}


##############


if (modeflag == 1 || modeflag==2) {
#check type
if ((l_type!="sig")&&(l_type!="ref")&&(l_type!="dif")&&(l_type!="src")) {
        printlog("ERROR - MVIEW: Image type keyword invalid.  lpar mview for valid values.",l_logfile,l_verbose)    
      status=1
      goto nextimage
}
}
 
if (modeflag == 1 || modeflag==2) {

#set filename for display type
 if (l_type == "src") {
     ndim=1
 } else if (l_type == "ref") {
     ndim=2
 } else if (l_type == "dif") {
     ndim=3
 } else if (l_type == "sig") {
     printf("WARNING - MVIEW: sig type not implemented, using dif")
     ndim=3
 } 
}

if ((l_type!="sig")&&(l_type!="ref")&&(l_type!="dif")&&(l_type!="src")) {
        printlog("ERROR - MVIEW: Images type invalid.  lpar mview for valid values.",l_logfile,l_verbose)    
      status=1
      goto nextimage
}

if (modeflag ==3 || modeflag == 4) {
ndim=1
}

cache("imgets","display")

#get axes

imgets(in[i]//"[1]","i_naxis3") ; naxis3=int(imgets.value)

if (naxis3!=3) {
   printlog("ERROR - MVIEW: Images "//in[i]//" is not the correct format.",l_logfile,l_verbose)
   print("                     n_choppos= "//naxis3)
      status=1
      goto nextimage
}

#check extensions
#check number of extensions

  imgets(header,"NUMEXT", >& "dev$null")
  nextns=int(imgets.value)

              nbadsets=0
tmpfile2=mktemp("tmpfile2")
tmpfile3=mktemp("tmpfile3")


#display images in non-interactive mode:

n_i=1

    if (!l_manual) {
for (n_i=1;n_i<=nextns;n_i+=1) {
        if (l_verbose) print("Displaying "//l_type//": Nod "//str(n_i))
display(in[i]//"["//str(n_i)//"][*,*,"//str(ndim)//",1]",1,erase-,zscale=lzscale,z1=lz1,z2=lz2,ztrans=lztrans,zrange=lzrange,fill=l_dispfill, >& "dev$null")
}
}

#now, make the output file that can have the header edited:

print(out[i])

n_i=1
    if (l_manual) {
        if (l_verbose) print("MIVIEW: Copying image "//i//" into a new file.")
fxcopy(in[i]//".fits[0]",out[i])
        if (l_verbose) print(l_type//": Nod "//str(n_i))
for (n_i=1;n_i<=nextns;n_i+=1) {
tmpfile1=mktemp("tmpfile1")
    imcopy(in[i]//"["//n_i//"]",tmpfile1,verbose=no)
    fxinsert(tmpfile1//".fits",out[i]//"["//n_i-1//"]",groups="",ver-)
}
}

#Starting the header editing for interactive mode:

    if (l_manual) {

for (n_i=1;n_i<=nextns;n_i+=1) {

    if (npars != 2) {
    npars=0
}

        if (l_verbose) print("Displaying "//l_type//": Nod "//str(n_i))
display(out[i]//"["//str(n_i)//"][*,*,"//str(ndim)//"]",1,erase-,zscale=lzscale,z1=lz1,z2=lz2,ztrans=lztrans,zrange=lzrange,fill=l_dispfill, >& "dev$null")

    sleep(l_delay)

      while (npars == 0) {
        npars=1
            printf("MIVIEW:  Starting interactive cursor input mode for nod:"//n_i//".  Press h for help.\n")

        if (fscan(imcur,xpos,ypos,wcs,cursinput) != EOF) {

          if (cursinput == "q" || cursinput == "Q") {
            if (l_verbose) printf("MIVIEW: Exiting interactive mode for nod: "//n_i//".\n")
            npars=1
          }

          if (cursinput == "x" || cursinput == "X") {
            if (l_verbose) printf("MIVIEW: Exiting interactive mode.\n")
          npars=2
          }

#
# For the "h" command, loop to get another input value.  Any value aside 
# from the defined ones causes the next images to be displayed.
#
          if (cursinput == "h" || cursinput == "H") {
            printf("--------MIVIEW: INTERACTIVE HELP-------- \n ")
            printf("Available key commands:\n (h) print this help \n "//
            "(b) mark as a bad frame \n"//
            " (u) unmark as a bad frame \n "// 
           "(i) run imexamine on this image.  \n "//
            "(q) stop interactive mode and move on to next nod \n "//
            "(s) get images statistics  \n"// 
            " (x) exit interactive mode immediately.  \n")
  printf("Key commands can be entered in upper or lower case.  Any undefined keystroke \n will automatically advance the display to the next nod image.\n")
            npars=0
            printf("-------------------------------------- \n ")
          }
          if (cursinput == "s" || cursinput == "S") {
            imstat(out[i]//"["//str(n_i)//"][*,*,"//str(ndim)//",1,1]")
            npars=0
          }
#
# For the moment, disable marking bad frames....
#

          if (cursinput == "b" || cursinput == "B") {
              if (l_verbose) printlog("MIVIEW: Nod "//n_i//" marked as bad.\n",l_logfile,l_verbose)
              hedit(out[i]//"["//str(n_i)//"]","BADNOD","1",add=yes,addonly=no,
	          delete=no,verify=no,show=no,update=yes)
              nbadsets=nbadsets+1
            npars=0
          }
          if (cursinput == "u" || cursinput == "U") {
              if (l_verbose) printf("MIVIEW: Nod "//n_i//" unmarked as bad.\n")
              hedit(out[i]//"["//str(n_i)//"]","BADNOD","",add=no,addonly=no,
	          delete=yes,verify=no,show=no,update=yes)
             nbadsets=nbadsets-1
            npars=0
          }
          if (cursinput == "i" || cursinput == "I") {
            if (l_verbose) printf("MIVIEW: Entering imexam.\n")
            imexamine()
            if (l_verbose) printf("MIVIEW: Exiting imexam.\n")
            npars=0
          }
          if (npars == 1) {
            if (l_verbose) {
              printf("Going to next frame.\n")
              printf("")
           }
          }
        }
     } 
}
          sleep(l_delay) 
        }


}

if (check != "0") {

  imgets(header,"NUMEXT", >& "dev$null")
  nextns=int(imgets.value)


print(out[i])
print("MIVIEW:  Data has been already been stacked with MISTACK, displaying coadded frame.")

#display images in non-interactive mode:

n_i=1

    if (!l_manual) {
        if (l_verbose) print("Displaying "//in[i])
display(in[i]//"["//str(n_i)//"][*,*]",1,erase-,zscale=lzscale,z1=lz1,z2=lz2,ztrans=lztrans,zrange=lzrange,fill=l_dispfill, >& "dev$null")
}

    if (l_manual) {
              nbadsets=0
        if (l_verbose) print("MIVIEW: Copying image "//i//" into a new file.")
fxcopy(in[i]//".fits[0]",out[i])
tmpfile1=mktemp("tmpfile1")
    imcopy(in[i]//"["//n_i//"]",tmpfile1,verbose=no)
    fxinsert(tmpfile1//".fits",out[i]//"["//n_i-1//"]",groups="",ver-)
}


#Starting the header editing for interactive mode:

    if (l_manual) {

    if (npars != 2) {
    npars=0
}

        if (l_verbose) print("Displaying "//out[i])

display(out[i]//"["//str(n_i)//"]",1,erase-,zscale=lzscale,z1=lz1,z2=lz2,ztrans=lztrans,zrange=lzrange,fill=l_dispfill, >& "dev$null")

    sleep(l_delay)

     

      while (npars == 0) {
        npars=1
            printf("MIVIEW:  Starting interactive cursor input mode for nod:"//n_i//".  Press h for help.\n")

        if (fscan(imcur,xpos,ypos,wcs,cursinput) != EOF) {

          if (cursinput == "q" || cursinput == "Q") {
            if (l_verbose) printf("MIVIEW: Exiting interactive mode for nod: "//n_i//".\n")
            npars=1
          }

          if (cursinput == "x" || cursinput == "X") {
            if (l_verbose) printf("MIVIEW: Exiting interactive mode.\n")
          npars=2
          }

     

#
# For the "h" command, loop to get another input value.  Any value aside 
# from the defined ones causes the next images to be displayed.
#
          if (cursinput == "h" || cursinput == "H") {
            printf("--------MIVIEW: INTERACTIVE HELP-------- \n ")
            printf("Available key commands:\n (h) print this help \n "//
            "(b) mark as a bad frame \n"//
            " (u) unmark as a bad frame \n "// 
           "(i) run imexamine on this image.  \n "//
            "(q) stop interactive mode and move on to next nod \n "//
            "(s) get images statistics  \n"// 
            " (x) exit interactive mode immediately for all images.  \n")
  printf("Key commands can be entered in upper or lower case.  Any undefined keystroke \n will automatically advance the display to the next nod image.\n")
            npars=0
            printf("-------------------------------------- \n ")
          }
          if (cursinput == "s" || cursinput == "S") {
            imstat(out[i]//"["//str(n_i)//"][*,*,"//str(ndim)//",1,1]")
            npars=0
          }
#
# For the moment, disable marking bad frames....
#

          if (cursinput == "b" || cursinput == "B") {
              if (l_verbose) printlog("MIVIEW: Nod "//n_i//" marked as bad.\n",l_logfile,l_verbose)
              hedit(out[i]//"["//str(n_i)//"]","BADNOD","1",add=yes,addonly=no,
	          delete=no,verify=no,show=no,update=yes)
              nbadsets=nbadsets+1
            npars=0
          }
          if (cursinput == "u" || cursinput == "U") {
              if (l_verbose) printf("MIVIEW: Nod "//n_i//" unmarked as bad.\n")
              hedit(out[i]//"["//str(n_i)//"]","BADNOD","",add=no,addonly=no,
	          delete=yes,verify=no,show=no,update=yes)
             nbadsets=nbadsets-1
            npars=0
          }
          if (cursinput == "i" || cursinput == "I") {
            if (l_verbose) printf("MIVIEW: Entering imexam.\n")
            imexamine()
            if (l_verbose) printf("MIVIEW: Exiting imexam.\n")
            npars=0
          }

     

          if (npars == 1) {
            if (l_verbose) {
              printf("")
              printf("")
           }
          }
        }
    }
}

          sleep(l_delay) 
}
       

if(! l_manual) {
delete(out[i],ver-, >& "dev$null")
}

if( l_manual) {
if (nbadsets == 0) {
delete(out[i],ver-, >& "dev$null")
printf("MIVIEW: No bad nodsets identified.  Header not changed, so no output image has been written to disk. \n")
} else {
printf("MIVIEW: The number of nodsets marked as bad is:"//nbadsets//".\n")
  date | scan(l_struct)
  printf("%-8s= \'%-18s\' / %-s\n","GEM-TLM",l_struct,"Last modification with GEMINI IRAF", >> tmphead)
  printf("%-8s= \'%-18s\' / %-s\n","MIVIEW",l_struct,"Time stamp for MIVIEW", >> tmphead)
   mkheader(out[i]//"[0]",tmphead,append+,verbose-)
   delete(tmphead,verify-, >& "dev$null")
}
}

printf("MIVIEW: Now done displaying image:"//i//".\n")


   if (check == "0") {
npars=0
}

print("-----------------------------------------------------")
nextimage:

i=i+1

}

if (l_testmode) time

#cleanup

#   print("Done.")

print("-----------------------------------------------------")

clean:

  if(status==0) {
    printlog("MIVIEW exit status:  good.\n",l_logfile,l_verbose)
  }
  if(status!=0) {
    printlog("MIVIEW: Exited with errors. \n",l_logfile,l_verbose)
  }

scanfile="" 
delete("tmpfile*.fits",ver-, >& "dev$null")
delete("tmphead*.fits",ver-, >& "dev$null")
delete("tmpin*.fits",ver-, >& "dev$null")

end

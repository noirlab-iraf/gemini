# Copyright(c) 2003-2015 Association of Universities for Research in Astronomy , Inc.

procedure mdfplot (intable, instrume)

# Plot MDF
#
# Version Apr 23, 2003  IJ  First draft, minimal checking of input
#         Jul 15, 2003  IJ  updated for 2003B conventions, handle GS automatically
#         Aug 12, 2003  IJ  updated for the 0/1 semester convention
#      Oct 13, 2003  KL  moved to gmos$mostools
#         Dec 10, 2003  IJ  backwards compatable, make it work for DD/SV programs
#         Aug 14, 2006  KC  fix tilt angle bug, constrain (0 < slittilt < 180)
#         Jun 22, 2009  PG  added prompt for Instrument and added F2
#         Aug  7, 2009  JH  clean ups and pyraf fix

char      intable   {prompt="MDF FITS file be plotted"}
char      instrume  {prompt="Instrument to use: GMOS-N, GMOS-S, F2?"}
bool      fl_ps     {yes,prompt="Make PostScript on disk"}
char      barcode   {"default",prompt="Barcode for mask"}
bool      fl_inter  {yes,prompt="Examine mask plot interactively"}
bool      fl_over   {yes,prompt="Overwrite existing PostScript on disk"}
char      logfile   {"",prompt="Logfile"}
bool      verbose   {yes,prompt="Verbose?"}
int       status    {0,prompt="Exit status (0=good)"}
char      *scanfile {"",prompt="For internal use only"}

begin

    char l_intable, l_logfile, l_barcode, l_instrume
    bool l_verbose, l_fl_ps, l_fl_over, l_fl_inter

    char l_rootname, l_postscript, tmpgki, s_empty, tmpigi, l_nrootname
    char l_mdfid, estring, l_graphcap
    struct l_struct
    int l_ii, l_rows, l_istr, count
    real xpos, ypos, tilt, xsize, ysize, y1, y2, x11, x12, x21, x22
    bool l_fl_old

    # Local variables used for BARCODE definition
    # Q=0, C=1, LP=2, SV=8, DS/DD=9
    # type_programs and type_numbers must remain in order. Number of types must
    # match the number of type_programs
    int num_types=7
    char type_programs="Q","C","LP","FT","SV","DS","DD"
    char type_numbers="0","1","2","3","8","9","9"
    int ii
    char semster, mask_num, prog_num, program_type, tmpstr, tmp_type, semester
    int dash_loc, prg_type_offset, num_loc, max_len_types

    max_len_types = 2

    l_intable=intable ; l_logfile=logfile ; l_verbose=verbose ; l_fl_ps=fl_ps
    l_barcode=barcode ; l_fl_over=fl_over ; l_fl_inter=fl_inter
    l_instrume=instrume
    status=0

    tmpgki=mktemp("tmpgki")
    tmpigi=mktemp("tmpigi")

    cache("gimverify","imgets","parkey","tinfo","fparse")

    # Store the current graphcap setting
    show graphcap | scan (l_graphcap)

    # Define the name of the logfile
    s_empty=""; print(l_logfile) | scan(s_empty); l_logfile=s_empty
    if (l_logfile == "") {
        l_logfile = mostools.logfile
        if (l_logfile == "") {
            l_logfile = "gmos.log"
            printlog("WARNING - MDFPLOT: both mdfplot.logfile and \
                mostools.logfile are empty.", l_logfile, l_verbose)
            printlog("                   Using default file gmos.log.", \
                l_logfile, l_verbose)
        }
    }

    # Write to the logfile
    date | scan(l_struct)
    printlog("----------------------------------------------------------------\
    ------------",l_logfile,l_verbose)
    printlog("MDFPLOT -- "//l_struct,l_logfile,l_verbose)
    printlog("",l_logfile,l_verbose)

    if (!l_fl_ps && !l_fl_inter) {
        printlog("ERROR - MDFPLOT: Using fl_ps=no and fl_inter=no produces no \
            output", l_logfile, verbose+)
        goto crash
    }

    # Make sure input table has .fits suffix
    gimverify(l_intable) ; l_intable=gimverify.outname//".fits"
    if (!access(l_intable)) {
        printlog("ERROR - MDFPLOT: Input table "//l_intable//" not found",
            l_logfile, verbose+)
        goto crash
    }

    fparse(l_intable) ; l_rootname=fparse.root

    if (access(l_rootname//".ps") && l_fl_over==no) {
        printlog("ERROR - MDFPLOT: Output PostScript file exists", l_logfile,
            verbose+)
        goto crash
    } else
        delete(l_rootname//".ps",verify-, >>& "dev$null")

    # Decide on the barcode
    if (l_barcode == "default") {
        dash_loc = stridx("-", l_rootname)
        prg_type_offset = 1
        if (dash_loc > 8) {
            dash_loc = 8
            prg_type_offset = 0
	}
        tmpstr = substr(l_rootname, dash_loc - 1, strlen(l_rootname))
        semester = substr(tmpstr, 1, 1)
        if (semester == "A")
            l_barcode="10"
        else
            l_barcode="11"

        if (substr(tmpstr, 2, 2) == "-") {
             dash_loc = 2
	} else {
             dash_loc = 1
	}
        tmpstr = substr(tmpstr, dash_loc + 1, strlen(tmpstr))

        program_type = substr(tmpstr, 1, max_len_types)
        for (ii = 1; ii <= num_types; ii += 1) {

            tmp_type = type_programs[ii]
            if (tmp_type == substr(program_type, 1, strlen(tmp_type))) {
                l_barcode = l_barcode//type_numbers[ii]
                break
            } else if (ii == num_types) {
                printlog("WARNING - MDFPLOT: Cannot assign barcode", \
                         l_logfile, verbose+)
                l_barcode = "None"
                goto NO_BARCODE
            }
        }
        tmpstr = substr(tmpstr, strlen(program_type) + 1, strlen(tmpstr))
        if (substr(tmpstr, 1, 1) == "-") {
             dash_loc = 2
	} else {
             dash_loc = 1
	}

        prog_num = substr(tmpstr, dash_loc, strldx("-", tmpstr) - 1)
        printf("%03d\n", int(prog_num)) | scan (prog_num)
        l_barcode = l_barcode//prog_num

        tmpstr = substr(tmpstr, strldx("-", tmpstr) + 1, strlen(tmpstr))
        mask_num = tmpstr
        printf("%02d\n", int(mask_num)) | scan (mask_num)
        l_barcode = l_barcode//mask_num

    }
NO_BARCODE:

    # Make the plot

    # Assign plotting sequence

    if ((l_instrume == "GMOS-N") || (l_instrume == "GMOS-S")){
        mgograph(l_intable,"slitpos_mx","slitpos_my",pointe=0.1,marker="cross",
            wx1=-120,wx2=120,wy1=-120,wy2=120,xlab="X mask",ylab="Y mask",
            title=l_rootname//" "//l_barcode,postitle="topleft",labelexp=1.5,
            device="gkifile", gkifile=tmpgki)
        # Draw the frame
        type("mostools$mdfplot.igi") | igi(append+,debug-,device="stdgraph", \
            >>G tmpgki)
    }

    if (l_instrume == "F2"){
        # Make the plot
        mgograph(l_intable,"slitpos_mx","slitpos_my",pointe=0.1,marker="cross",
            wx1=-125,wx2=125,wy1=-125,wy2=125,xlab="X mask",ylab="Y mask",
            title=l_rootname//" "//l_barcode,postitle="topleft",labelexp=1.5,
            device="gkifile", gkifile=tmpgki)
        # Draw the frame
        type("mostools$mdfplotf2.igi") | igi(append+,debug-, \
            device="stdgraph", >>G tmpgki)
    }

    # Make igi command file that draws the slits with the correct size and tilt
    tinfo(l_intable,ttout-, >>& "dev$null")
    l_rows=tinfo.nrows
    for(l_ii=1;l_ii<=l_rows;l_ii+=1) {
        tprint(l_intable, showr-, row=str(l_ii), showh-, option="plain",
            prdata+, prparam-, col="slitpos_mx,slitpos_my,slittilt_m,\
            slitsize_mx,slitsize_my",pwidth=156) | \
            scan(xpos,ypos,tilt,xsize,ysize)

        # Constrain slit tilt angle between 0 < tilt < 180 (effects only plot)
        if (tilt < 0.0)
            tilt = tilt + 360.
        if (tilt > 180.0)
            tilt = tilt - 180.

        y1=ypos-ysize/2. ; y2=ypos+ysize/2.
        x11=xpos-xsize/2.+ysize/2*(-cos(tilt/180.*3.1416))
        x12=xpos+xsize/2.+ysize/2*(-cos(tilt/180.*3.1416))
        x21=xpos-xsize/2.-ysize/2*(-cos(tilt/180.*3.1416))
        x22=xpos+xsize/2.-ysize/2*(-cos(tilt/180.*3.1416))

        if (l_ii==1) {
            if ((l_instrume == "GMOS-N") || (l_instrume == "GMOS-S"))
                head("mostools$mdfplot.igi", nlines=3, > tmpigi)
            else if (l_instrume == "F2")
                head("mostools$mdfplotf2.igi", nlines=3, > tmpigi)
        }

        print("relo "//str(x11)//" "//str(y1), >> tmpigi)
        print("draw "//str(x12)//" "//str(y1)//" ; draw \
            "//str(x22)//" "//str(y2),>> tmpigi)
        print("draw "//str(x21)//" "//str(y2)//" ; draw \
            "//str(x11)//" "//str(y1),>> tmpigi)
    }

    type(tmpigi) | igi(append+,debug-,device="stdgraph", >>G tmpgki)

    if (l_fl_inter)
        gkimosaic(tmpgki, device="stdgraph", output="", nx=1, ny=1, \
            rotate=no, fill=no, interactive=yes, cursor="")

    if (l_fl_ps) {
        # Set the graphcap to mostools$graphcap, which contains a new device
        # that creates a postscript file called iraf.ps
        set graphcap="mostools$data/graphcap"

        # Create the postscript file
        gkimosaic(tmpgki, device="psfsq", output="", nx=1, ny=1, \
            rotate=no, fill=no, interactive=yes, cursor="")
        gflush

        # Set the graphcap back to the original setting
        set graphcap=str(l_graphcap)

        # Add a check here to make sure iraf.ps file exists before continuing
        # Sleep 3/5 is not enough ... odf2mdf RTF tests are failing randomly
        # with "Warning: Cannot open file (iraf.ps)" error
        count = 0
        while (!access("iraf.ps") && count < 3) {
            printlog("WARNING - MDFPLOT: Output from gkimosaic has not yet \
                been created, waiting ...", l_logfile, verbose+)
            sleep 3
            date
            count += 1
        }
        if (!access("iraf.ps")) {
            printlog("ERROR - MDFPLOT: Problem with gkimosaic", l_logfile, \
                verbose+)
            goto crash
        } else {
            printlog("MDFPLOT - Writing PostScript file "// l_rootname // \
                ".ps", l_logfile, l_verbose)
            rename("iraf.ps", l_rootname//".ps", field="all")
        }
    }

    printlog("MDFPLOT - MDF "//l_rootname//" Barcode "//l_barcode, l_logfile,
        l_verbose)

    goto clean

crash:
    # Set the graphcap back to the original setting
    set graphcap=str(l_graphcap)

    status=1
    goto theend

clean:
    status=0

theend:
    delete(tmpgki//","//tmpigi, verify-, >>& "dev$null")
    printlog("----------------------------------------------------------------\
    ------------",l_logfile, l_verbose)

end

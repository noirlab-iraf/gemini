# Copyright(c) 2002-2017 Association of Universities for Research in Astronomy, Inc.

procedure gfquick (image)

# Version  Sept 20, 2002  IH,IJ v1.4 releas
# Version Sept 27, 2005 PG now it deletes the *ifu* files with a warning
# Version Nov 15, 2005 PG added instrument identifier and definition of xpos1
#                         added MASKNAME identifier
# Version Dec 02, 2005 PG added "default' entries for xpos1, nslit, and ioffu
#                         as well as other cosmetic changes

string  image       {prompt = "Frame to be reconstructed" }
string  rawdir      {"adata$",prompt = "Directory containing raw data"}
string  xpos1       {"default",prompt = "X-position for slit 1 on chip 1, 1210 for GMOS-N"}
string  nslit       {"default",prompt = "1 or 2 slits?, default is 2 slits"}
bool    fl_off      {yes,prompt = "Calculate offsets on reconstructed image?"}
int     ioff        {1,prompt = "Use which image for offsets?"}
bool    fl_spos     {yes,prompt = "re-determine y-shift in fibre positions?"}
string  yoffu       {"default",prompt= "y-shift for IFU fibre positions. Offset is 12"}
char    issport     {"side",enum="up|side",prompt="ISS port (up|side)"}


begin

    char    l_image, l_rawdir, secstr, imlogfile, offim, tmpim, brun1x, brun3x
    char    l_issport, l_nslit, l_yoffu
    bool    l_fl_off, l_fl_spos
    int     xlo1, xhi1, xlo2, xhi2, l_ioff, xlo, xhi, ylo, yhi
    int     yref, yoff, l_xpos1, l_xpos2, inst, flag_0
    int     naxis1, naxis2, ampinteg, numslit, fnumslit, status
    real    xoffi, yoffi, xncen, yncen, xs, ys, vxpos1[2], yrefa[2]
    real    ysoffu
    char    sa_1, sa_2, sb_1, sb_2, sc_1, sc_2, sd_1, sd_2, se_1, se_2
    char    sf_1, sf_2, sg_1, sg_2
    char    sh_1, sh_2, si_1, si_2, sj_1, sj_2, s1_1, s1_2, s2_1, s2_2
    char    s3_1, s3_2, s4_1, s4_2
    char    s5_1, s5_2, s6_1, s6_2, s7_1, s7_2, s8_1, s8_2, s9_1, s9_2
    char    s10_1, s10_2
    char    s11_1, s11_2, s12_1, s12_2, s13_1, s13_2, s14_1, s14_2
    char    s15_1, s15_2, s16_1, s16_2
    char    s17_1, s17_2, s18_1, s18_2, s19_1, s19_2, s20_1, s20_2
    char    mslit, bslit, mslit1, mslit2
    char    t_xpos1
    char    tmpcoo, tmppref
    bool    ishamamatsu

    l_image = image
    l_rawdir = rawdir
    l_ioff = ioff
    l_issport = issport
    l_nslit = nslit
    t_xpos1 = xpos1
    l_fl_off = fl_off
    l_fl_spos = fl_spos
    l_yoffu = yoffu
    tmpcoo="uparm$"//mktemp("tmpcoo")
    flag_0 = 1
    status = 0

    tmppref = mktemp("tmpti")

    gtile (l_image, outimages="", outpref=tmppref, out_ccds="all", \
        ret_roi=yes, req_roi=0, fl_stats_only=no, fl_tile_det=no, \
        fl_app_rois=no, fl_pad=no, sci_pad=0., var_pad=0., \
        dq_pad=16, sci_fakeval=0., var_fakeval=0., dq_fakeval=16., \
        chipgap="default", sci_ext="", var_ext="VAR", dq_ext="DQ", \
        mdf_ext="MDF", key_detsec="DETSEC", key_ccdsec="CCDSEC", \
        key_datasec="DATASEC", key_biassec="BIASSEC", key_ccdsum="CCDSUM", \
        rawpath=l_rawdir, logfile="dev$null", fl_verbose=no)

    if (gtile.status != 0) {
        goto crash
    }

    l_image = tmppref//l_image

    # Due to gtile call
    l_rawdir = ""

    # Which instrument?
    imgets (l_rawdir//l_image//"[0]", "INSTRUME", >& "dev$null")
    if (imgets.value == "0") {
        print ("WARNING - GFQUICK: Instrument keyword not found.")
        inst = 1
    }
    if (imgets.value == "GMOS-N")
        inst = 1
    else
        inst = 2

    # read the proper xpos1
    # Relative to the Datasec
    vxpos1[1] = 1210. + 5. -32.0 # Adjusted becuase of gtile call
    vxpos1[2] = 1181.

    if ((t_xpos1 == "default") && (inst == 1)) {
        l_xpos1 = vxpos1[1]
    }
    if ((t_xpos1 == "default") && (inst == 2)) {
        l_xpos1 = vxpos1[2]
    }
    if (t_xpos1 != "default") {
        gemisnumber (t_xpos1, "decimal", verbose-)
        if (gemisnumber.fl_istype)
            l_xpos1 = real(t_xpos1)
        else {
            print ("ERROR - GFQUICK: xpos1 is not a number.")
            bye 
        }    
    }
    
    # Which IFU? 

    if (l_nslit != "default") {
        gemisnumber (l_nslit, "decimal", verbose-)
        if (gemisnumber.fl_istype)
            numslit = int(l_nslit)
        else {
            print ("ERROR - GFQUICK: nslit has a wrong value")
            bye 
        }    
    } else {
        imgets (l_rawdir//l_image//"[0]", "MASKNAME", >& "dev$null")
        if (imgets.value == "0") {
            print ("WARNING - GFQUICK: MASKNAME keyword not found.")
            numslit = 2
        }
        fnumslit = 20
        
        if (imgets.value == "IFU-1") {
            fnumslit = 1
        } else if (imgets.value == "IFU-R") {
            fnumslit = 3
        } else if (imgets.value == "IFU-B") {
            fnumslit = 4
        } else if (imgets.value == "IFU-2") {
            fnumslit = 2
        } else if (imgets.value == "IFU-NS-B") {
            fnumslit = 1
        } else if (imgets.value == "IFU-NS-R") {
            fnumslit = 1
        } else if (imgets.value == "IFU-NS-2") {
            fnumslit = 2
        } else if (imgets.value == "IFU") {
            print ("ERROR - GFQUICK: Unable to determine the IFU configuration")
            print ("                 Specify parameter 'nslit'")
            bye
        }
        switch (fnumslit) {
            case 1:
                numslit = 1
            case 3:
                numslit = 1
            case 4:
                numslit = 1
            case 2:
                numslit = 2
            default: {
                print ("ERROR - GFQUICK: invalid IFU type on MASKNAME. Maybe not \
                    an IFU image?")
                bye 
            }
        }
    }

    if ((numslit > 2) || (numslit < 1)) {
        print ("ERROR - GFQUICK: nslit has a value different than 1 or 2")
        bye 
    }   

    # this is really a y-position on the original frame!
    yrefa[1] = 2137
    yrefa[2] = 1909
    if (inst == 1) {
        keypar (l_image//"[0]", "DETTYPE", silent+)
        if (str(keypar.value) == "SDSU II e2v DD CCD42-90") {
            yrefa[1] = yrefa[1] - 22
            yrefa[2] = yrefa[2] - 22
	    }
        yref = yrefa[1]
    }

    ishamamatsu = no
    if (inst == 2) {
        yref = yrefa[2]
        keypar (l_image//"[0]", "DETTYPE", silent+)
        if (str(keypar.value) == "S10892" || str(keypar.value) == "S10892-N") {
            ishamamatsu = yes
        }
    }

    brun1x = "uparm$"//mktemp("brun1x")
    brun3x = "uparm$"//mktemp("brun3x")

    imlogfile = "uparm$"//mktemp("imlogfile")

    if (yes == imaccess(l_image//'_ifu1')) {
        print ("WARNING: Deleting "//l_image//"_ifu1")
        imdelete (image=l_image//'_ifu1' , verify=no, >& "dev$null")
    }

    if (yes == imaccess (l_image//'_ifu2')) {
        print ("WARNING: Deleting "//l_image//"_ifu2")
        imdelete (image=l_image//'_ifu2' , verify=no, >& "dev$null")
    }


    specred.dispaxis=2

    # set position of second output slit based on the position of the 
    # first slit
    l_xpos2 = l_xpos1-337 

    xlo1 = l_xpos1-8
    xhi1 = l_xpos1+9

    secstr = "[1]["//xlo1//":"//xhi1//",*]"

    imcopy (input=l_rawdir//l_image//secstr, out=brun1x, verbose=no)
    # next line is useful for debugging
    #display(image=brun1x, frame=1)
    #imexam(input=brun3x, frame=1)

    imtranspose (input=brun1x//"[-*,*]", output=brun1x)
    blkavg (input=brun1x, output='t1', b1=1, b2=18)
    imdelete (images=brun1x, verify=no)
    blkrep (input='t1', output=brun1x, b1=1, b2=4)
    imdelete (images='t1', verify=no, >& "dev$null")

    if (numslit == 2) {

        # these went with blkav=16
        #xlo2=xpos2-8
        #xhi2=xpos2+7
        # need a wider piece because fibres not so well aligned
        xlo2 = l_xpos2-8
        xhi2 = l_xpos2+9

        secstr = "[3]["//xlo2//":"//xhi2//",*]"

        imcopy (input=l_rawdir//l_image//secstr, out=brun3x, verbose=no)
        # next line is useful for debugging
        #display(image=brun3x, frame=1)
        #imexam(input=brun3x, frame=1)

        imtranspose (input=brun3x//"[-*,*]", output=brun3x)
        blkavg (input=brun3x,output='t3', b1=1, b2=18)
        imdelete (images=brun3x, verify=no, >& "dev$null")
        blkrep (input='t3',output=brun3x, b1=1, b2=4)
        imdelete (images='t3', verify=no, >& "dev$null")
    }

    # option of redetermining the offsets in y

    if (l_fl_spos == yes) {

        tmpim = "uparm$"//mktemp("tmpim")
        xlo = l_xpos1 - 8
        xhi = l_xpos1 + 9

        ylo = yref - 100
        yhi = yref + 100

        secstr = "[1]["//xlo//":"//xhi//","//ylo//":"//yhi//"]"
        imcopy (input=l_rawdir//l_image//secstr, out=tmpim, verbose=no)
        display (image=tmpim, frame=1, zs-, zr+)
        # This will display a circle around the top fiber of the lower set.
        print (str(10), " ", str(108), "Top", > tmpcoo)
        tvmark (1, tmpcoo, mark="circle", radii=10, color=204, label-)

        delete (tmpcoo, verify-, >& "dev$null")
        print ("**************** Determine fibre y-shift *********************")
        print ("Put cursor on top fibre of the lower set which should be")
        print ("roughly centered on red circle, press x then press q")
        print ("")
        print ("IF fibers are too faint, then then just press q and the")
        print ("script will run with default values")
        if (yes == access(imlogfile))
            delete(imlogfile, verify=no, >& "dev$null")

        imexam (input=tmpim, frame=1, keeplog=yes, logfile=imlogfile,
            >& "dev$null")
        xs = -100.      # dummy negative number
        ys = -100.      # dummy negative number
        fields (files=imlogfile, fields="1,2", lines="2-") | scan (xs, ys)
#        print(xs,ys)

        if ((xs >= -99.) && (ys >= -99.)) {

            yoff = int(ys) - 100 - 1

            l_xpos1 = int(xs) + xlo
            print ("*********************************************************\
                *****")
            print ("Y-offset measured=", yoff," pixels")
            if (yes == imaccess(tmpim))
                imdelete (tmpim, verify=no, >& "dev$null")

        } else {

            ysoffu = 12.
            yoff = 12.
            print ("**********************************************************\
                ****")
            print ("")
            print ("WARNING: xpos1 set to yes but using default values")
            print ("Y-offset assumed=", yoff," pixels")
            print ("X-position for slit 1 assumed=", l_xpos1," pixels")
            print ("**********************************************************\
                ****")
            flag_0 = 0
       
        }  

    } else {

        # xoff=6 for the R400
        # I estimate R831 offset should be R400+4
        # from flats, B600=R400+46pix
        #      imgets (l_rawdir//l_image//"[0]", "GRATING", >>& "dev$null")
        #      if (imgets.value =="B600+_G5303") {
        #         print("GFQUICK - Applying y-offset for B600 grating")
        #         yoff=52
        #      } else if (imgets.value == "R400+_G5305") {
        #         print("GFQUICK - Applying offset for R400 grating")
        #         yoff=6  
        #      } else if (imgets.value == "MIRROR") {
        #         print("GFQUICK - Applying offset for Mirror")
        ## Value used in commissioning
        ##         yoff=0
        ## Value needed Jan 2002
        #          yoff=15
        #      }

        if (l_yoffu == "default" ) {
            ysoffu=12.
        }

        if (l_yoffu != "default") {
            gemisnumber (l_yoffu, "decimal", verbose-)
            if (gemisnumber.fl_istype)
                ysoffu = real(l_yoffu)
            else {
                print ("ERROR - GFQUICK: yoffu is not a number.")
                bye 
            }    
        }

        yoff = ysoffu
    }

    

    yref = yref + yoff

    # new x value 
    imdelete (images=brun3x, verify=no, >& "dev$null")
    imdelete (images=brun1x, verify=no, >& "dev$null")
    if (flag_0 == 1) {
        print ("X-position for slit 1 measured=", l_xpos1," pixels")
    }
    l_xpos2 = l_xpos1 - 337 + 32

    xlo1 = l_xpos1-8
    xhi1 = l_xpos1+9
#   xlo1 = l_xpos1-4
 #  xhi1 = l_xpos1+5


    secstr = "[1]["//xlo1//":"//xhi1//",*]"

    imcopy (input=l_rawdir//l_image//secstr, out=brun1x, verbose=no)
    # next line is useful for debugging
    #display(image=brun1x, frame=1)
    #imexam(input=brun3x, frame=1)

    imtranspose (input=brun1x//"[-*,*]", output=brun1x)
    # 18 and 4
    blkavg (input=brun1x, output='t1', b1=1, b2=18)
    imdelete (images=brun1x, verify=no)
    blkrep (input='t1', output=brun1x, b1=1, b2=4)
    imdelete (images='t1', verify=no, >& "dev$null")

    if (numslit == 2) {

        # these went with blkav=16
        #xlo2=xpos2-8
        #xhi2=xpos2+7
        # need a wider piece because fibres not so well aligned
        xlo2 = l_xpos2 - 8
        xhi2 = l_xpos2 + 9

        secstr = "[3]["//xlo2//":"//xhi2//",*]"

        imcopy (input=l_rawdir//l_image//secstr, out=brun3x, verbose=no)
        # next line is useful for debugging
        #display(image=brun3x, frame=1)
        #imexam(input=brun3x, frame=1)

        imtranspose (input=brun3x//"[-*,*]", output=brun3x)
        blkavg (input=brun3x,output='t3', b1=1, b2=18)
        imdelete (images=brun3x, verify=no, >& "dev$null")
        blkrep (input='t3',output=brun3x, b1=1, b2=4)
        imdelete (images='t3', verify=no, >& "dev$null")
    }


    # subtract overscan level (roughly)
    hselect (l_rawdir//l_image//"[0]", "AMPINTEG", yes) | scan (ampinteg)
    if (ishamamatsu) {
        if (ampinteg == 11880) {
            # slow mode
            imarith (operand1=brun1x, op='-', operand2=756, result=brun1x)
            if (numslit == 2)
                imarith (operand1=brun3x, op='-', operand2=387, result=brun3x)

        } else if (ampinteg == 4000) {
            # fast mode
            imarith (operand1=brun1x, op='-', operand2=662, result=brun1x)
            if (numslit == 2)
                imarith (operand1=brun3x, op='-', operand2=337, result=brun3x)
        }
    
    } else {
        if (ampinteg == 5000) {
            # slow mode
            imarith (operand1=brun1x, op='-', operand2=756, result=brun1x)
            if (numslit == 2)
                imarith (operand1=brun3x, op='-', operand2=387, result=brun3x)

        } else if (ampinteg == 1000) {
            # fast mode
            imarith (operand1=brun1x, op='-', operand2=662, result=brun1x)
            if (numslit == 2)
                imarith (operand1=brun3x, op='-', operand2=337, result=brun3x)
        }
    }

    # background field
    imdelete (images='s*.fits', verify=no, >& "dev$null")
    print ("GFQUICK - Cutting out the sections")

    sa_1=mktemp("sa_1")
    sa_2=mktemp("sa_2")
    sb_1=mktemp("sb_1")
    sb_2=mktemp("sb_2")
    sc_1=mktemp("sc_1")
    sc_2=mktemp("sc_2")
    sd_1=mktemp("sd_1")
    sd_2=mktemp("sd_2")
    se_1=mktemp("se_1")
    se_2=mktemp("se_2")
    sf_1=mktemp("sf_1")
    sf_2=mktemp("sf_2")
    sg_1=mktemp("sg_1")
    sg_2=mktemp("sg_2")
    sh_1=mktemp("sh_1")
    sh_2=mktemp("sh_2")
    si_1=mktemp("si_1")
    si_2=mktemp("si_2")
    sj_1=mktemp("sj_1")
    sj_2=mktemp("sj_2")
    s1_1=mktemp("s1_1")
    s1_2=mktemp("s1_2")
    s2_1=mktemp("s2_1")
    s2_2=mktemp("s2_2")
    s3_1=mktemp("s3_1")
    s3_2=mktemp("s3_2")
    s4_1=mktemp("s4_1")
    s4_2=mktemp("s4_2")
    s5_1=mktemp("s5_1")
    s5_2=mktemp("s5_2")
    s6_1=mktemp("s6_1")
    s6_2=mktemp("s6_2")
    s7_1=mktemp("s7_1")
    s7_2=mktemp("s7_2")
    s8_1=mktemp("s8_1")
    s8_2=mktemp("s8_2")
    s9_1=mktemp("s9_1")
    s9_2=mktemp("s9_2")
    s10_1=mktemp("s10_1")
    s10_2=mktemp("s10_2")
    s11_1=mktemp("s11_1")
    s11_2=mktemp("s11_2")
    s12_1=mktemp("s12_1")
    s12_2=mktemp("s12_2")
    s13_1=mktemp("s13_1")
    s13_2=mktemp("s13_2")
    s14_1=mktemp("s14_1")
    s14_2=mktemp("s14_2")
    s15_1=mktemp("s15_1")
    s15_2=mktemp("s15_2")
    s16_1=mktemp("s16_1")
    s16_2=mktemp("s16_2")
    s17_1=mktemp("s17_1")
    s17_2=mktemp("s17_2")
    s18_1=mktemp("s18_1")
    s18_2=mktemp("s18_2")
    s19_1=mktemp("s19_1")
    s19_2=mktemp("s19_2")
    s20_1=mktemp("s20_1")
    s20_2=mktemp("s20_2")

    imcopy (input=brun1x//"["//  918+yref//":"// 1059+yref//",*]",output=sf_1,verbose=no)
    imcopy (input=brun1x//"["// 1201+yref//":"// 1060+yref//",*]",output=sf_2,verbose=no)
    imcopy (input=brun1x//"["//   19+yref//":"//  160+yref//",*]",output=sg_1,verbose=no)
    imcopy (input=brun1x//"["//  302+yref//":"//  161+yref//",*]",output=sg_2,verbose=no)
    imcopy (input=brun1x//"["// -878+yref//":"// -737+yref//",*]",output=sh_1,verbose=no)
    imcopy (input=brun1x//"["// -595+yref//":"// -736+yref//",*]",output=sh_2,verbose=no)
    imcopy (input=brun1x//"["//-2078+yref//":"//-1937+yref//",*]",output=si_1,verbose=no)
    imcopy (input=brun1x//"["//-1795+yref//":"//-1936+yref//",*]",output=si_2,verbose=no)
    imcopy (input=brun1x//"["// 2117+yref//":"// 2258+yref//",*]",output=sj_1,verbose=no)
    imcopy (input=brun1x//"["// 2400+yref//":"// 2259+yref//",*]",output=sj_2,verbose=no)

    if (numslit == 2) {
        imcopy (input=brun3x//"["//-2078+yref//":"//-1937+yref//",*]",output=sa_1,verbose=no)
        imcopy (input=brun3x//"["//-1795+yref//":"//-1936+yref//",*]",output=sa_2,verbose=no)
        imcopy (input=brun3x//"["// 2118+yref//":"// 2259+yref//",*]",output=sb_1,verbose=no)
        imcopy (input=brun3x//"["// 2401+yref//":"// 2260+yref//",*]",output=sb_2,verbose=no)
        imcopy (input=brun3x//"["//  918+yref//":"// 1059+yref//",*]",output=sc_1,verbose=no)
        imcopy (input=brun3x//"["// 1201+yref//":"// 1060+yref//",*]",output=sc_2,verbose=no)
        imcopy (input=brun3x//"["//   19+yref//":"//  160+yref//",*]",output=sd_1,verbose=no)
        imcopy (input=brun3x//"["//  302+yref//":"//  161+yref//",*]",output=sd_2,verbose=no)
        imcopy (input=brun3x//"["// -878+yref//":"// -737+yref//",*]",output=se_1,verbose=no)
        imcopy (input=brun3x//"["// -595+yref//":"// -736+yref//",*]",output=se_2,verbose=no)
    }

    imcopy (input=brun1x//"["//-1778+yref//":"//-1637+yref//",*]",output=s11_1,verbose=no)
    imcopy (input=brun1x//"["//-1495+yref//":"//-1636+yref//",*]",output=s11_2,verbose=no)
    imcopy (input=brun1x//"["//-1477+yref//":"//-1336+yref//",*]",output=s12_1,verbose=no)
    imcopy (input=brun1x//"["//-1194+yref//":"//-1335+yref//",*]",output=s12_2,verbose=no)
    imcopy (input=brun1x//"["//-1178+yref//":"//-1037+yref//",*]",output=s13_1,verbose=no)
    imcopy (input=brun1x//"["// -895+yref//":"//-1036+yref//",*]",output=s13_2,verbose=no)
    imcopy (input=brun1x//"["// -580+yref//":"// -439+yref//",*]",output=s14_1,verbose=no)
    imcopy (input=brun1x//"["// -297+yref//":"// -438+yref//",*]",output=s14_2,verbose=no)
    imcopy (input=brun1x//"["// -280+yref//":"// -139+yref//",*]",output=s15_1,verbose=no)
    imcopy (input=brun1x//"["//    2+yref//":"// -139+yref//",*]",output=s15_2,verbose=no)
    imcopy (input=brun1x//"["//  318+yref//":"//  459+yref//",*]",output=s16_1,verbose=no)
    imcopy (input=brun1x//"["//  600+yref//":"//  459+yref//",*]",output=s16_2,verbose=no)
    imcopy (input=brun1x//"["//  619+yref//":"//  760+yref//",*]",output=s17_1,verbose=no)
    imcopy (input=brun1x//"["//  901+yref//":"//  760+yref//",*]",output=s17_2,verbose=no)
    imcopy (input=brun1x//"["// 1214+yref//":"// 1355+yref//",*]",output=s18_1,verbose=no)
    imcopy (input=brun1x//"["// 1497+yref//":"// 1356+yref//",*]",output=s18_2,verbose=no)
    imcopy (input=brun1x//"["// 1516+yref//":"// 1657+yref//",*]",output=s19_1,verbose=no)
    imcopy (input=brun1x//"["// 1799+yref//":"// 1658+yref//",*]",output=s19_2,verbose=no)
    imcopy (input=brun1x//"["// 1818+yref//":"// 1959+yref//",*]",output=s20_1,verbose=no)
    imcopy (input=brun1x//"["// 2101+yref//":"// 1960+yref//",*]",output=s20_2,verbose=no)
    imdelete (images=brun1x, verify=no, >& "dev$null")

    if (numslit == 2) {
        imcopy (input=brun3x//"["//-1776+yref//":"//-1635+yref//",*]",output=s10_1,verbose=no)
        imcopy (input=brun3x//"["//-1493+yref//":"//-1634+yref//",*]",output=s10_2,verbose=no)
        imcopy (input=brun3x//"["//-1476+yref//":"//-1335+yref//",*]",output=s9_1,verbose=no)
        imcopy (input=brun3x//"["//-1193+yref//":"//-1334+yref//",*]",output=s9_2,verbose=no)
        imcopy (input=brun3x//"["//-1176+yref//":"//-1035+yref//",*]",output=s8_1,verbose=no)
        imcopy (input=brun3x//"["// -893+yref//":"//-1034+yref//",*]",output=s8_2,verbose=no)
        imcopy (input=brun3x//"["// -578+yref//":"// -437+yref//",*]",output=s7_1,verbose=no)
        imcopy (input=brun3x//"["// -295+yref//":"// -436+yref//",*]",output=s7_2,verbose=no)
        imcopy (input=brun3x//"["// -280+yref//":"// -139+yref//",*]",output=s6_1,verbose=no)
        imcopy (input=brun3x//"["//    3+yref//":"// -138+yref//",*]",output=s6_2,verbose=no)
        imcopy (input=brun3x//"["//  320+yref//":"//  461+yref//",*]",output=s5_1,verbose=no)
        imcopy (input=brun3x//"["//  603+yref//":"//  462+yref//",*]",output=s5_2,verbose=no)
        imcopy (input=brun3x//"["//  618+yref//":"//  759+yref//",*]",output=s4_1,verbose=no)
        imcopy (input=brun3x//"["//  901+yref//":"//  760+yref//",*]",output=s4_2,verbose=no)
        imcopy (input=brun3x//"["// 1218+yref//":"// 1359+yref//",*]",output=s3_1,verbose=no)
        imcopy (input=brun3x//"["// 1501+yref//":"// 1360+yref//",*]",output=s3_2,verbose=no)
        imcopy (input=brun3x//"["// 1518+yref//":"// 1659+yref//",*]",output=s2_1,verbose=no)
        imcopy (input=brun3x//"["// 1801+yref//":"// 1660+yref//",*]",output=s2_2,verbose=no)
        imcopy (input=brun3x//"["// 1817+yref//":"// 1958+yref//",*]",output=s1_1,verbose=no)
        imcopy (input=brun3x//"["// 2100+yref//":"// 1959+yref//",*]",output=s1_2,verbose=no)
        imdelete (images=brun3x, verify=no, >& "dev$null")
    }

    print("GFQUICK - Joining the pieces together")
    if (numslit == 2) {

        mslit = s20_2//","//s20_1//","//s19_2//","//s19_1//","//s18_2//","//\
            s18_1//","//s17_2//","//s17_1//","//s16_2//","//s16_1//","//\
            s15_2//","//s15_1//","//s14_2//","//s14_1//","//s13_2//","//\
            s13_1//","//s12_2//","//s12_1//","//s11_2//","//s11_1//","//\
            s10_1//","//s10_2//","//s9_1//","//s9_2//","//s8_1//","//\
            s8_2//","//s7_1//","//s7_2//","//s6_1//","//s6_2//","//\
            s5_1//","//s5_2//","//s4_1//","//s4_2//","//s3_1//","//\
            s3_2//","//s2_1//","//s2_2//","//s1_1//","//s1_2

        bslit = sj_2//","//sj_1//","//si_2//","//si_1//","//sh_2//","//\
            sh_1//","//sg_2//","//sg_1//","//sf_2//","//sf_1//","//\
            se_1//","//se_2//","//sd_1//","//sd_2//","//sc_1//","//\
            sc_2//","//sb_1//","//sb_2//","//sa_1//","//sa_2


    } else if (numslit == 1) {

        mslit = s20_2//","//s20_1//","//s19_2//","//s19_1//","//s18_2//","//\
            s18_1//","//s17_2//","//s17_1//","//s16_2//","//s16_1//","//\
            s15_2//","//s15_1//","//s14_2//","//s14_1//","//s13_2//","//\
            s13_1//","//s12_2//","//s12_1//","//s11_2//","//s11_1

        bslit = sj_2//","//sj_1//","//si_2//","//si_1//","//sh_2//","//\
            sh_1//","//sg_2//","//sg_1//","//sf_2//","//sf_1

    }

    imjoin (input=mslit, output=l_image//'_ifu1', join_dimension=2, verbose=no)
    imdelete (images=mslit,verify=no, >& "dev$null")
    imjoin (input=bslit, output=l_image//'_ifu2', join_dimension=2, verbose=no)
    imdelete (images=bslit, verify=no, >& "dev$null")

    # scale the output image to the right aspect ratio and flip 
    magnify (input=l_image//'_ifu1', output=l_image//'_ifu1', xmag=1.0,
        ymag=1.23, >& "dev$null")
    magnify (input=l_image//'_ifu2', output=l_image//'_ifu2', xmag=1.0,
        ymag=1.23, >& "dev$null")
    imtranspose (input=l_image//'_ifu1[*,-*]', output=l_image//'_ifu1',
        >& "dev$null")
    imtranspose (input=l_image//'_ifu2[*,-*]', output=l_image//'_ifu2',
        >& "dev$null")

    print ("++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
    print ("GFQUICK - Output target     image is ", l_image//'_ifu1')
    print ("GFQUICK - Output background image is ", l_image//'_ifu2')
    print ("++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
    print ("GFQUICK - Scale of output images is 0.035 arcsec/pix")
    print ("GFQUICK - Note that flux level is not preserved")

    if (l_fl_off == yes) {

        offim = l_image//'_ifu'//l_ioff

        display (image=offim, frame=1, zs-)


        print ("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
        print ("Put cursor on the object to be centred and press x, then press q")
        print ("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")

        if (yes == access(imlogfile)) 
            delete (imlogfile,verify=no, >& "dev$null")
        imexam (input=offim, frame=1, keeplog=yes, logfile=imlogfile,
            >& "dev$null")
        fields (files=imlogfile, fields="1,2", lines="2-") | scan (xs, ys)
        print ("Object coordinates (pix)=", xs, ys)

        hselect (offim, "i_naxis1,i_naxis2", yes) | scan (naxis1, naxis2)
        xncen = naxis1 / 2.0
        yncen = naxis2 / 2.0

        print ("Field centre (pix)=", xncen, yncen)

        xoffi = (xncen-xs) * 0.035    
        yoffi = (yncen-ys) * 0.035    

        print ("")
        print ("***************** RESULTS *****************")

        print ("Offsets derived for the "//l_issport//"-looking port")

        if (l_issport=="side")
            xoffi = -xoffi

        print ("Instrument offsets to apply (arcsec): P=",str(-xoffi),
            " Q=",str(yoffi))

        if (yes == access(imlogfile))
            delete (imlogfile, verify=no, >& "dev$null")
    }

    imdelete (l_image, verify=no, >& "dev$null")

goto clean

crash:
    status = 1

clean:

    if (status == 0) {
         print ("GFQUICK - Exit staus = GOOD")
    } else {
         print ("GFQUICK - Exit staus = ERROR")
    }

end

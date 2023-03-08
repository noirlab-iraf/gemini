# Copyright(c) 2004-2009 Association of Universities for Research in Astronomy, Inc.

include	<ctotok.h>
include	<imhdr.h>
include	<ctype.h>
include	<mach.h>
include	<imset.h>
include	<fset.h>
include	<lexnum.h>
include	<evvexpr.h>
include "../../../lib/mefio/mefio.h"
include	"gettok.h"
include "gemexpr.h"

define FULLOPNAM_LEN	32
define MAX_NDIM			100
# MIMEXPR -- Task procedure for the image expression evaluator.  This task
# generates an image by evaluating an arbitrary vector expression, which may
# reference other images as input operands.
#
# The input expression may be any legal EVVEXPR expression.  Input operands
# must be specified using the reserved names "a" through "z", hence there are
# a maximum of 26 input operands.  An input operand may be an image name or
# an image header parameter, a numeric constant, or the name
# of a builtin keyword.  Image header parameters are specified as, e.g.,
# "a.naxis1" where the operand "a" must be assigned to an input image.  The
# special image name "." refers to the output image generated in the last
# call to imexpr, making it easier to perform a sequence of operations.

procedure t_mimexpr()

double	dval
bool	verbose, rangecheck
pointer	out, st, sp, ie, dims, intype, outtype, ref_im
pointer	outim, fname, expr, xexpr, output, section, data, imname
pointer	oplist, opnam, opval, param, io, ip, op, o, im, ia, emsg
#GEMEXPR extravars
pointer phumef		# holds string name of MEF from which phu will be copied
pointer mep, extp, mepary, corary
pointer cluster					# temp to use to get cluster portions
pointer tmpimname,tmpimname2 	# used to copy PHU... and as needed
pointer tmpimp				 	# used to check for PHU, and as needed

int index, extveri
pointer extname
pointer parymalloc(), paryget(), gdrelate()
pointer mignls(),mignli(),mignll(),mignlr(),mignld()
pointer i2dget()
int  errget()

bool phucopied
int refimext
int extent[MAX_NDIM]
int ii
#GEMEXPR END

int	len_exprbuf, fd, nchars, noperands, dtype, status, i, j
int	ndim, npix, ch, percent, nlines, totlines, flags

double	imgetd()
int	imgftype(), ctod()
bool	clgetb(), imgetb(), streq(), strne()
pointer memap(), megetep()
int	imgnls(), imgnli(), imgnll(), imgnlr(), imgnld()
int	impnls(), impnli(), impnll(), impnlr(), impnld()
int	open(), getci(), ie_getops(), lexnum(), stridxs()
int	imgeti(), ctoi(), btoi(), strncmp()
pointer locpr()
pointer	ie_getexprdb(), ie_expandtext(), immap()
extern	ie_getop(), ie_fcn()
pointer	evvexpr()
long	fstatl()

char fullopnam[FULLOPNAM_LEN]
char tmpstr[SZ_LINE]
int tmpi

bool worktodo
bool bphu
int  htype
bool	ldebug

#for getting phu header vals
int megftype()
double megetd()
int megeti()
bool megetb()
int strlen()

string	s_nodata "bad image: no data"
string	s_badtype "unknown image type"

# goto labels
define	numeric_ 91
define	image_ 92
define  extiter_ 93


errchk megetep, megbnv, memaxver, meindbnv, countbn, mesetiter, 
errchk megnep, megntep, memap, evvexpr
begin
    ldebug = false 
 
    outim = 0
    
    phucopied = false  # phu is copied within the iterative loop if it does
    #	not exist already... 
    #   once done stop checking to see if it exists
	
    if (ldebug)
    {
	    call printf("t_mimexpr(): starting mimexpr\n");
        call flush(STDOUT)
        verbose = true # if debug output is off, this is verbose
    }

    # call memlog ("--------- START IMEXPR -----------")

    call smark (sp)
    call salloc (ie, LEN_IMEXPR, TY_STRUCT)
    call salloc (fname, SZ_PATHNAME, TY_CHAR)
    call salloc (output, SZ_PATHNAME, TY_CHAR)
    call salloc (imname, SZ_PATHNAME, TY_CHAR)
    call salloc (section, SZ_FNAME, TY_CHAR)
    call salloc (cluster, SZ_FNAME, TY_CHAR)
    call salloc (intype, SZ_FNAME, TY_CHAR)
	call salloc (outtype, SZ_FNAME, TY_CHAR)
	call salloc (oplist, SZ_LINE, TY_CHAR)
	call salloc (opval, SZ_LINE, TY_CHAR)
	call salloc (dims, SZ_LINE, TY_CHAR)
	call salloc (emsg, SZ_LINE, TY_CHAR)
	#GEMEXPR
	call salloc (extname, SZ_LINE, TY_CHAR)
	call salloc (tmpimname, SZ_FNAME, TY_CHAR)
	call salloc (tmpimname2, SZ_FNAME, TY_CHAR)
	call salloc (phumef, SZ_FNAME, TY_CHAR)
	Memc[phumef] = EOS
	
    
	call clgstr("mimexprpars.extname", Memc[extname], SZ_LINE)
	# Initialize the main imexpr descriptor.
	call aclri (Memi[ie], LEN_IMEXPR)

	verbose = clgetb ("verbose")
	rangecheck = clgetb ("rangecheck")

	# Load the expression database, if any.
	st = NULL
	call clgstr ("mimexprpars.exprdb", Memc[fname], SZ_PATHNAME)
	if (strne (Memc[fname], "none"))
	    st = ie_getexprdb (Memc[fname])
	IE_ST(ie) = st

	# Get the expression to be evaluated and expand any file inclusions
	# or macro references.
    
	len_exprbuf = SZ_COMMAND
	call malloc (expr, len_exprbuf, TY_CHAR)
	call clgstr ("mimexprpars.expr", Memc[expr], len_exprbuf)

	if (Memc[expr] == '@') {
	    fd = open (Memc[expr+1], READ_ONLY, TEXT_FILE)
	    nchars = fstatl (fd, F_FILESIZE)
	    if (nchars > len_exprbuf) {
		len_exprbuf = nchars
		call realloc (expr, len_exprbuf, TY_CHAR)
	    }
	    for (op=expr;  getci(fd,ch) != EOF;  op = op + 1) {
		if (ch == '\n')
		    Memc[op] = ' '
		else
		    Memc[op] = ch
	    }
	    Memc[op] = EOS
	    call close (fd)
	}

	if (st != NULL) {
	    xexpr = ie_expandtext (st, Memc[expr])
	    call mfree (expr, TY_CHAR)
	    expr = xexpr
	    if (ldebug) {
			call printf ("expr=%s\n")
		    call pargstr (Memc[expr])
		    call flush (STDOUT)
	    }
	}

	# Get output image name.
	call clgstr ("mimexprpars.output", Memc[output], SZ_PATHNAME)
	call imgimage (Memc[output], Memc[imname], SZ_PATHNAME)
	call strcat("[APPEND]", Memc[output],SZ_PATHNAME)

    # commented out as these parameters were removed on Inger's request
	IE_BWIDTH(ie) = 0 # clgeti ("bwidth")
	IE_BTYPE(ie)  = 2 # clgwrd ("btype", Memc[oplist], SZ_LINE, "|constant|nearest|reflect|wrap|project|")
    # note: oplist is just use as temporary in the commented out code above
    
	IE_BPIXVAL(ie) = 0 # clgetr ("bpixval")

	# Determine the minimum input operand type.
	call clgstr ("intype", Memc[intype], SZ_FNAME)

	if (strncmp (Memc[intype], "default", 4) == 0)
	    IE_INTYPE(ie) = 0
	else {
	    switch (Memc[intype]) {
	    case 'i', 'l':
		IE_INTYPE(ie) = TY_INT
	    case 'r':
		IE_INTYPE(ie) = TY_REAL
	    case 'd':
		IE_INTYPE(ie) = TY_DOUBLE
	    default:
		IE_INTYPE(ie) = 0
	    }
	}

	# NOTE: Below is code to fill the IE_ and IO_ imexpr and operand structs.
	# NOTE: This code has been duplicated in a function ie_fill that could
	# NOTE: be called here... however, this loop also sets the PHU
	# NOTE: source name... I need the function in t_gemexpr.x and will
	# NOTE: use it here later.  Note... it's not an exact copy and this
	# NOTE: loop below would require refactoring, not just replacement.

	# Parse the expression and generate a list of input operands.
	noperands = ie_getops (st, Memc[expr], Memc[oplist], SZ_LINE)
	IE_NOPERANDS(ie) = noperands

	
	# Process the list of input operands and initialize each operand.
	# This means fetch the value of the operand from the CL, determine
	# the operand type, and initialize the image operand descriptor.
	# The operand list is returned as a sequence of EOS delimited strings.

	opnam = oplist
	do i = 1, noperands {

	    io = IE_IMOP(ie,i)
	    if (Memc[opnam] == EOS) {
		    call error (GEERR_MALFORMEDOP, "malformed operand list")
        }

		call mkfullopname("mimexprpars", Memc[opnam], fullopnam, FULLOPNAM_LEN)
	    call clgstr (fullopnam, Memc[opval], SZ_LINE)
		if (ldebug)
		{
			call printf("opname,opval = %s,%s\n");
			call pargstr(Memc[opnam])
			call pargstr(Memc[opval])
			call flush(STDOUT)
		}
	    IO_OPNAME(io) = Memc[opnam]
	    ip = opval

        # Initialize the input operand; these values are overwritten below.
	    o = IO_OP(io)
	    call aclri (Memi[o], LEN_OPERAND)

	    if (Memc[ip] == '.' && (Memc[ip+1] == EOS || Memc[ip+1] == '[')) {
		# A "." is shorthand for the last output image.
		call strcpy (Memc[ip+1], Memc[section], SZ_FNAME)
		call clgstr ("lastout", Memc[opval], SZ_LINE)
		call strcat (Memc[section], Memc[opval], SZ_LINE)
		goto image_

	    } else if (IS_LOWER(Memc[ip]) && (Memc[ip+1] == '.' ||
                                ( Memc[ip+1] == '[' && Memc[ip+2] == '0' &&
                                   Memc[ip+3] == ']' && Memc[ip+4] == '.' ))) {
		# "a.foo" refers to parameter foo of image A.  Mark this as
		# a parameter operand for now, and patch it up later.
            
        # Extended to support a[0].foo which refers to the PHU value for "foo"

		IO_TYPE(io) = PARAMETER
        
 
		    # posterity --> imexpr bug: IO_DATA(io) = ip

            tmpi = strlen(Memc[ip])+1
            call salloc(IO_DATA(io),tmpi, TY_CHAR)
            call strcpy(Memc[ip], Memc[IO_DATA(io)], tmpi)


	    } else if (ctod (Memc, ip, dval) > 0) {
		if (Memc[ip] != EOS)
		    goto image_

		# A numeric constant.
numeric_	IO_TYPE(io) = NUMERIC

		ip = opval
		switch (lexnum (Memc, ip, nchars)) {
		case LEX_REAL:
		    dtype = TY_REAL
		    if (stridxs("dD",Memc[opval]) > 0 || nchars > NDIGITS_RP+3)
			dtype = TY_DOUBLE
		    O_TYPE(o) = dtype
		    if (dtype == TY_REAL)
			O_VALR(o) = dval
		    else
			O_VALD(o) = dval
		default:
		    O_TYPE(o) = TY_INT
		    O_LEN(o)  = 0
		    O_VALI(o) = int(dval)
		}

	    } else {
		# Anything else is assumed to be an image name.
image_
		ip = opval
		call imgimage (Memc[ip], Memc[fname], SZ_PATHNAME)
		if (streq (Memc[fname], Memc[imname]))
		    call error (GEERR_SAMEINOUT, "input and output images cannot be the same")
			
			#
            call clgstr("mimexprpars.refim", tmpstr, SZ_LINE)
			if (streq(Memc[opnam],tmpstr) && (Memc[phumef] == EOS)) {
				# if phumef is not set, this must be first image,
				# and therefore source of PHU for output MEF
				call imgcluster(Memc[ip], Memc[phumef], SZ_FNAME)
				if (ldebug) {
					call sprintf (tmpstr, SZ_LINE, "Will Use PHU from \"%s\" if needed")
					call pargstr(Memc[phumef])
                    call log_info(tmpstr)
				}
			}
# GEMEXPR
            call imgsection(Memc[ip], tmpstr, SZ_LINE)
            if ( strne(tmpstr,"") ) {
                # note: sections have some internal support, but 
                #  there are unresolved issue regarding how they
                #  should work, so currently, they are prohibited
                #  although you will see them "handled" in other
                #  parts of the code.  The prohibition is at this
                #  point as a gatekeeper.
                call sprintf(tmpstr,SZ_LINE, "FITS image sections are not supported,\n       MEF filenames only.\n       Operand was \"%s\"")
                call pargstr(Memc[ip])
                call error(MEERR_NOSECTIONS, tmpstr)
            }

			call imgimage(Memc[ip], Memc[cluster], SZ_FNAME)
            
            if  ((strncmp(Memc[ip+1],"[0].",4) == 0) 
                                    && (strne(Memc[ip], Memc[cluster], SZ_FNAME))) {
                call sprintf(tmpstr,SZ_LINE, "FITS kernel sections and other qualifiers\n       are not supported, MEF filenames only.\n       Operand was \"%s\"")
                call pargstr(Memc[ip])
                call error(MEERR_NOSECTIONS, tmpstr)
            }
            if (ldebug) {
                call printf("imgimage = %s\n")
                call pargstr(Memc[ip])
                call flush(STDOUT)
            }
            
			mep = memap (Memc[cluster])
            
            #extp = megbnv(mep, Memc[extname], 1, READ_ONLY)
			#extp = megetep(mep, 1, READ_ONLY)

			IO_MEF(io) = mep
            IO_TYPE(io) = IMAGE

            # call meepfunmap(extp, true)
	    }

	    # Get next operand name.
	    while (Memc[opnam] != EOS)
		opnam = opnam + 1
	    opnam = opnam + 1
	}

	# GEMEXPR WORK:
	# We go through the operands and make the empty frame relation table
    # by identifying which extensions are SCI
	mepary = parymalloc(noperands)
	do i = 1, noperands {
		index = i - 1
		io = IE_IMOP(ie,i)
		ip = IO_DATA(io)
		if (IO_TYPE(io) == IMAGE)
		{
			call paryset(mepary, index, IO_MEF(io))	
            # kludge alert... this conditional is here to tell fmrelate which
            #  operands are supposed to be SCI frames (i.e. when gemexpr is
            #  processing var_exp)  This is OK since it was decided not to 
            #  expose mimexpr as a task on it's own.  If it was exposed we'd
            #  have to go back to the generalization and not assume operands
            #  over "n" refer to SCI frames
            if (IO_OPNAME(io) >= 'n') {
                ME_OPTYPE(IO_MEF(io)) = OPT_SCI
            } else {
                ME_OPTYPE(IO_MEF(io)) = OPT_NORMAL
            }
		}
		else
		{
			call paryset(mepary, index, 0)
		}
	}
    
    iferr (	corary = gdrelate(mepary, noperands, Memc[extname])){
        tmpi = errget(tmpstr, SZ_LINE)
        call clputi("status", 1)
        call error(tmpi, tmpstr)
    }

    if(ldebug) {
        call i2dprint(corary)
    }
    if ((I2D_NROWS(corary) <= 0)) {
        call sprintf(tmpstr, SZ_LINE, "No work to do for EXTNAME=%s, not present")
		call pargstr(Memc[extname])
        call log_warn(tmpstr)
        call clputi("status", GEWARN_NOFRAMES)
		return
	}
    worktodo = false
	do i = 1, I2D_NROWS(corary) {
        if (i2dget(corary, 0 , i) == I2DROWVALID) {
            worktodo = true
            call flush(STDOUT)
            break
        }
    }
        
	if (worktodo == false) {
        call sprintf(tmpstr, SZ_LINE, "No work for EXTNAME=%s, not present in some image operands")
		call pargstr(Memc[extname])
        call log_warn(tmpstr)
		return
	}
	
# ----------------------------
# TOP OF EXTVER ITERATION LOOP 
# ----------------------------
	for (extveri = 1; extveri <= I2D_YW(corary); extveri = extveri + 1)
	{
		if (ldebug)
		{
			call printf("Start Processing extveri=#%d...\n")
			call pargi(extveri)
			call flush(STDOUT)
		}
        if (i2dget(corary,0,extveri) == I2DROWINVAL) {
            next
        }
# clear extent checking array for each extver pass
        call aclri(extent, MAX_NDIM)
	    do i = 1, noperands {
		    io = IE_IMOP(ie, i)
            o = IO_OP(io)
            ip = IO_DATA(io)
            
#            call printf("str%d\n")
#            call pargi(ip)
#            call pargstr(Memc[ip])
#            call flush(STDOUT)

            if (IO_TYPE(io) == IMAGE) {
# update current images operand structures
			    if (ldebug) {
					call printf("extp = megetep(%x, i2dget(%x, %d, %d)[==%d], READ_ONLY)\n")
					call pargi(IO_MEF(io))
					call pargi(corary)
					call pargi(i)
					call pargi(extveri)
					call pargi(i2dget(corary,i,extveri))
					call flush(STDOUT)
			    }
			
             	extp = megetep(IO_MEF(io),i2dget(corary, i, extveri), READ_ONLY)
			    IO_IM(io) = EXT_EXTP(extp)
			    im = IO_IM(io)
			    if(ldebug)	{
					call printf("extp=%x\tim=%x\n")
					call pargi(extp)
					call pargi(im)
				}
            #check extent matching requirements
                for (ii = 1; ii <= IM_NDIM(im); ii = ii+1) {
                    if (extent[ii]==0) {
                        extent[ii] = IM_LEN(im,ii)
                        if (ldebug) {
                            call printf("extent[%d]=%d, IM_LEN(im,%d) = %d\n")
                            call pargi(ii)
                            call pargi(extent[ii])
                            call pargi(ii)
                            call pargi(IM_LEN(im,ii))
                            call flush(STDOUT)
                        }
                    }
                    if ((extent[ii] != 0) && (extent[ii] != IM_LEN(im,ii))) {
                        call error(GEERR_EXTENTMISTMATCH, "The images are not the same size")
                    }
                }

			# Set any image options.
			if (IE_BWIDTH(ie) > 0) {
			call imseti (im, IM_NBNDRYPIX, IE_BWIDTH(ie))
		    	call imseti (im, IM_TYBNDRY, IE_BTYPE(ie))
		    	call imsetr (im, IM_BNDRYPIXVAL, IE_BPIXVAL(ie))
			}
			IO_TYPE(io) = IMAGE
			call amovkl (1, IO_V(io,1), IM_MAXDIM)
			IO_IM(io) = im		

			switch ( IM_PIXTYPE(im)) {
				case TY_SHORT, TY_INT, TY_LONG, TY_REAL, TY_DOUBLE:
				    O_TYPE(o) = IM_PIXTYPE(im)
				case TY_COMPLEX:
				    O_TYPE(o) = TY_REAL
				default:			# TY_USHORT
				    O_TYPE(o) = TY_INT
			}
                
                if (ldebug) {
                    call printf("intype %d - pixtype %d - otype %d\n")
                    call pargi(IE_INTYPE(ie))
                    call pargi(IM_PIXTYPE(im))
                    call pargi(O_TYPE(o))
                    call flush(STDOUT)
                }
                    

		O_TYPE(o) = max (IE_INTYPE(ie), O_TYPE(o))
		O_LEN(o) = IM_LEN(im,1)
		O_FLAGS(o) = 0

                if (ldebug) {
                    call printf("2- intype %d - pixtype %d - otype %d\n")
                    call pargi(IE_INTYPE(ie))
                    call pargi(IM_PIXTYPE(im))
                    call pargi(O_TYPE(o))
                    call flush(STDOUT)
                }
                
		# If one dimensional image read in data and be done with it.
		if (IM_NDIM(im) == 1) {
		    switch (O_TYPE(o)) {

		    case TY_SHORT:
			if (imgnls (im, IO_DATA(io), IO_V(io,1)) == EOF)
			    call error (GEERR_NODATA, s_nodata)

		    case TY_INT:
			if (imgnli (im, IO_DATA(io), IO_V(io,1)) == EOF)
			    call error (GEERR_NODATA, s_nodata)

		    case TY_LONG:
			if (imgnll (im, IO_DATA(io), IO_V(io,1)) == EOF)
			    call error (GEERR_NODATA, s_nodata)

		    case TY_REAL:
			if (imgnlr (im, IO_DATA(io), IO_V(io,1)) == EOF)
			    call error (GEERR_NODATA, s_nodata)

		    case TY_DOUBLE:
			if (imgnld (im, IO_DATA(io), IO_V(io,1)) == EOF)
			    call error (GEERR_NODATA, s_nodata)

		    default:
			    call error (GEERR_BADTYPE, s_badtype)
		    }
		}
		}
	}

        #parm updates, after IMAGE OPERAND PROCESSING BECAUSE we need to use
        # latest extension loaded for that operand on this extver iteration
        
        do i = 1, noperands {
	        io = IE_IMOP(ie,i)
	        ip = IO_DATA(io)
            
	        if (IO_TYPE(io) != PARAMETER)
		        next
            bphu = false
            
            # Locate referenced symbolic image operand (e.g. "a").
	        ia = NULL
	        do j = 1, noperands {
		        ia = IE_IMOP(ie,j)
		        if (IO_OPNAME(ia) == Memc[ip] && IO_TYPE(ia) == IMAGE)
		            break
		        ia = NULL
	        }

	        if (ia == NULL) {
		        call sprintf (Memc[emsg], SZ_LINE,
		        "bad images parameter reference %s")
		        call pargstr (Memc[ip])
			    call error (GEERR_BADIMPARMREF, Memc[emsg])
	        }

            mep = IO_MEF(ia)
            
            if (false) {
                call printf("\n%s --> %s\n")
                call pargstr(Memc[ip])
                call pargstr(IO_OPNAME(ia))
                call flush(STDOUT)
            }
            
            # check if this is a PHU parameter
            if (Memc[ip+1]  == '[' && Memc[ip+2] == '0' && 
                Memc[ip+3] == ']' && Memc[ip+4] == '.') {
                # the parameter is in the phu, true
                bphu = true
            }
            
            
            # Get the parameter value and set up operand struct.
	        im = IO_IM(ia)
            if (bphu) {
                param = ip + 5
            } else {
	            param = ip + 2
            }
            
            if (ldebug) {
                call printf("param name is %s (in %s) (mep=%x)\n")
                call pargstr(Memc[param])
                if (bphu) {
                    call pargstr("phu")
                } else {
                    call pargstr("extension header")
                }
                call pargi(mep)
                call flush(STDOUT)
            }
            #IO_TYPE(io) = NUMERIC
	        o = IO_OP(io)
	        O_LEN(o) = 0
            
            if (bphu) {
                htype = megftype(mep, Memc[param])
            } else {
                htype = imgftype (im, Memc[param])
            }
	        switch (htype) {
	            case TY_BOOL:
		        O_TYPE(o) = TY_BOOL
                if (bphu) {
                    O_VALI(o) = btoi (megetb(mep,Memc[param]))
                } else {
		            O_VALI(o) = btoi (imgetb (im, Memc[param]))
                }
                if (ldebug) {
                    call printf("b:%s == %b\n")
                    call pargstr(Memc[ip])
                    call pargi(O_VALI(o))
                    call flush(STDOUT)
                }
                
	            case TY_CHAR:
		        O_TYPE(o) = TY_CHAR
		        O_LEN(o)  = SZ_LINE
		        call malloc (O_VALP(o), SZ_LINE, TY_CHAR)
                if (bphu) {
                    call error(1,"string values from phu not supported")
                } else {
		            call imgstr (im, Memc[param], O_VALC(o), SZ_LINE)
                }
                if (ldebug) {
                    call printf("c:%s == %s\n")
                    call pargstr(Memc[ip])
                    call pargstr(O_VALC(o))
                    call flush(STDOUT)
                }
                
	            case TY_INT:
		        O_TYPE(o) = TY_INT
                if (bphu) {
                    O_VALI(o) = megeti (mep, Memc[param])
                } else {
                    O_VALI(o) = imgeti (im, Memc[param])
                }
                if (ldebug) {
                    call printf("i:%s == %d\n")
                    call pargstr(Memc[ip])
                    call pargi(O_VALI(o))
                    call flush(STDOUT)
                }
                
	            case TY_REAL:
		        O_TYPE(o) = TY_DOUBLE
                if (bphu) {
                    O_VALD(o) = megetd (mep, Memc[param])
                } else {
                    O_VALD(o) = imgetd (im, Memc[param])
                }
                if (ldebug) {
                    call printf("d:%s == %g\n")
                    call pargstr(Memc[ip])
                    call pargd(O_VALD(o))
                    call flush(STDOUT)
                }
                
	            default:
		        call sprintf (Memc[emsg], SZ_LINE, "param %s not found\n")
		        call pargstr (Memc[ip])
		        call error (GEERR_PARAMNOTFOUND, Memc[emsg])
	        }
            if (ldebug) {
                call printf("\n")
                call flush(STDOUT)
            }
	    }
        
        if (false) {
            call printf("I'm here\n")
            call flush(STDOUT)
        }
	# Determine the reference image from which we will inherit image
	# attributes such as the WCS.  If the user specifies this we use
	# the indicated image, otherwise we use the input image operand with
	# the highest dimension.

	call clgstr ("mimexprpars.refim", Memc[fname], SZ_PATHNAME)
	if (streq (Memc[fname], "default")) {
#	call printf("default\n")
#	call flush(STDOUT)
	    # Locate best reference image (highest dimension).
	        ndim = 0
	        ref_im = NULL
            
	        do i = 1, noperands {
		        io = IE_IMOP(ie,i)
		        if (IO_TYPE(io) != IMAGE || IO_IM(io) == NULL)
		            next
                
		        im = IO_IM(io)
		        if (IM_NDIM(im) > ndim) {
		            ref_im = im
                    refimext = i
		            ndim = IM_NDIM(im)
		        }
	        }
            # phumef should not be set yet in refim==default case
            if (Memc[phumef] == EOS) {
				# if phumef is not set, this must be first image,
				# and therefore source of PHU for output MEF
                if (refimext == 0) {
                    call error(GEERR_BADREFIM, "Reference Image Bad (PHU ERR)")
                }
				call imgcluster(ME_FILENAME(paryget(mepary, refimext-1)), Memc[phumef], SZ_FNAME)
            
                if (ldebug) {
				    call sprintf (tmpstr, SZ_LINE, "Will Use PHU from default \"%s\" if needed")
				    call pargstr(Memc[phumef])
                    call log_info(tmpstr)
			    }
            }
	    } else {
#	call printf("noauto\n")
#	call flush(STDOUT)
# Locate referenced symbolic image operand (e.g. "a").
	        io = NULL
	        do i = 1, noperands {
		        io = IE_IMOP(ie,i)
          
		        if (IO_OPNAME(io) == Memc[fname] && IO_TYPE(io) == IMAGE) {
                    refimext = i                    
		            break
                }
		        io = NULL
	        }
	        if (io == NULL) {
		        call sprintf (Memc[emsg], SZ_LINE, "bad wcsimage reference image %s")
		        call pargstr (Memc[fname])
		        call error (GEERR_BADWCS, Memc[emsg])
	        }
	        ref_im = IO_IM(io)
	    }

        if (ldebug) {
	        call printf("ref_im=(%s)%x-op#%d\n")
            call pargstr(Memc[fname])
	        call pargi(ref_im)
            call pargi(refimext)
	        call flush(STDOUT)
		}
# Copy PHU: IFF the image  does not exist (otherwise merely append to
# current file)
		#only do this once... probably should move this outside the loop (TODO::?)
		if (phucopied == false) {
			# NOTE:: this is really GDIO:: work... should move to MEFIO if that's
			#        where GDIO stuff goes.
			call imgcluster(Memc[output],Memc[tmpimname], SZ_FNAME)
			#open PHU, if it fails... we copy a PHU ot outname
			call strcat("[0,APPEND]", Memc[tmpimname], SZ_FNAME)

			if (ldebug)
			{
				call printf("out = %s\n")
				call pargstr (Memc[tmpimname])
				call flush(STDOUT)
			}
			
			iferr (	tmpimp = immap(Memc[tmpimname], READ_ONLY, NULL) )
			{
			# if here then the image does not exist (extension 0 doesn't exist so...)
			#loop over operands... find first image operand
			# Locate referenced symbolic image operand (e.g. "a").

				io = NULL
				do i = 1, noperands {
					io = IE_IMOP(ie,i)
					if (IO_TYPE(io) == IMAGE)
						break
					io = NULL
				}
				if (io == NULL)
				{
					call error(GEERR_NOIMOP, "ERROR: NO IMAGE OPERANDS\n")
				}
				call strcpy(Memc[phumef], Memc[tmpimname2], SZ_FNAME)
				call strcat("[0]", Memc[tmpimname2], SZ_FNAME)
# copy PHU from first image operand to output[0]
				if (verbose) {
					call sprintf(tmpstr, SZ_LINE,"PHU copy: %s --> %s)")
					call pargstr(Memc[tmpimname2])
					call pargstr(Memc[tmpimname])
					call flush(STDOUT)
                    call log_info(tmpstr)
				}
				call mimcopy(Memc[tmpimname2], Memc[tmpimname], verbose)
				if (ldebug)	{
					call printf("PHU copy...done\n")
					call flush(STDOUT)
				}
			} else {
			# file did exist, close it
				call imunmap(tmpimp)
				if (ldebug)	{
					call printf("PHU already present in %s\n")
					call pargstr(Memc[tmpimname])
					call flush(STDOUT)
				}
			}
			phucopied = true
		}

	# Determine the dimension and size of the output image.  If the "dims"
	# parameter is set this determines the image dimension, otherwise we
	# determine the best output image dimension and size from the input
	# images.  The exception is the line length, which is determined by
	# the image line operand returned when the first line of the image
	# is evaluated.

	call clgstr ("dims", Memc[dims], SZ_LINE)
	if (streq (Memc[dims], "default")) {
	    # Determine the output image dimensions from the input images.
	    call amovki (1, IE_AXLEN(ie,2), IM_MAXDIM-1)
	    IE_AXLEN(ie,1) = 0
	    ndim = 1

	    do i = 1, noperands {
		io = IE_IMOP(ie,i)
		im = IO_IM(io)
		if (IO_TYPE(io) != IMAGE || im == NULL)
		    next

		ndim = max (ndim, IM_NDIM(im))
		do j = 2, IM_NDIM(im) {
		    npix = IM_LEN(im,j)
		    if (npix > 1) {
			if (IE_AXLEN(ie,j) <= 1)
			    IE_AXLEN(ie,j) = npix
			else
			    IE_AXLEN(ie,j) = min (IE_AXLEN(ie,j), npix)
		    }
		}
	    }
	    IE_NDIM(ie) = ndim
	} else {
	    # Use user specified output image dimensions.
	    ndim = 0
	    for (ip=dims;  ctoi(Memc,ip,npix) > 0;  ) {
		ndim = ndim + 1
		IE_AXLEN(ie,ndim) = npix
		for (ch=Memc[ip];  IS_WHITE(ch) || ch == ',';  ch=Memc[ip])
		    ip = ip + 1
	    }
	    IE_NDIM(ie) = ndim
	}

	# Determine the pixel type of the output image.
	call clgstr ("outtype", Memc[outtype], SZ_FNAME)

#		call printf("outtype= %s\n")
#		call pargstr(Memc[outtype])
#		call flush(STDOUT)
		
	if (strncmp (Memc[outtype], "default", 7) == 0) {
	    IE_OUTTYPE(ie) = 0
	} else if (strncmp (Memc[outtype], "ref", 3) == 0) {
	    if (ref_im != NULL)
		IE_OUTTYPE(ie) = IM_PIXTYPE(ref_im)
	    else
		IE_OUTTYPE(ie) = 0
	} else {
	    switch (Memc[outtype]) {
	    case 'u':
		IE_OUTTYPE(ie) = TY_USHORT
	    case 's':
		IE_OUTTYPE(ie) = TY_SHORT
	    case 'i':
		IE_OUTTYPE(ie) = TY_INT
	    case 'l':
		IE_OUTTYPE(ie) = TY_LONG
	    case 'r':
		IE_OUTTYPE(ie) = TY_REAL
	    case 'd':
		IE_OUTTYPE(ie) = TY_DOUBLE
	    default:
		call error (GEERR_BADOUTTYPE, "bad outtype")
	    }
	}

	# Open the output image.  If the output image name has a section we
	# are writing to a section of an existing image.

	    call imgsection (Memc[output], Memc[section], SZ_FNAME)
        call imgcluster(Memc[output],tmpstr, SZ_LINE)
        call sprintf(Memc[tmpimname], SZ_FNAME, "%s[%s,%d,APPEND]")
        call pargstr(tmpstr)
        call pargstr(Memc[extname])
        call pargi(extveri)
        
	if (Memc[section] != EOS) {
            # this isn't allowed and should be prohibited by the 
            # calling routines.  however, we may bring sections back in... I'll emit error
            # for now but not remove old code
            call error(GEERR_GENERR, "image coordinate sections not allowed in output")
            if (ldebug) {
                call printf("outim=%x\n")
                call pargi(outim)
                call flush(STDOUT)
            }
	        outim = immap (Memc[output], READ_WRITE, 0)
	        IE_AXLEN(ie,1) = IM_LEN(outim,1)
	} else {
	    if (ref_im != NULL) {
				if (ldebug)	{
					call printf("NEW_COPY extension to %s (%s)\n")
					call pargstr(Memc[output])
                    call pargstr(Memc[tmpimname])
                    call printf("pointer outim=%x\n")
                    call pargi(outim)
                    call flush(STDOUT)
            	}	
                
				iferr (outim = immap (Memc[tmpimname], NEW_COPY, ref_im)) {
                    tmpi = errget(tmpstr, SZ_LINE)
                    if (ldebug) {
                        call printf("immap failed: %d - %s\n")
                        call pargi(tmpi)
                        call pargstr(tmpstr)
                        call flush(STDOUT)
                    }
                    call log_err(GEERR_BADCOPY, "extension creation failed\n(ensure refim image header is well formed/valid)")
                    call error(GEERR_BADCOPY, "immap NEW_COPY failed")
                }
			} else {
				if (ldebug)
				{
					call printf("NEW_IMAGE extension to %s\n")
					call pargstr(Memc[output])
					call flush(STDOUT)
				}	
				iferr (outim = immap (Memc[output], NEW_IMAGE, 0)) {
                       call error(GEERR_BADCOPY, "immap NEW_IMAGE failed")
                }
			}

	        IM_LEN(outim,1) = 0
	        call amovl (IE_AXLEN(ie,2), IM_LEN(outim,2), IM_MAXDIM-1)
	        IM_NDIM(outim) = IE_NDIM(ie)
	        IM_PIXTYPE(outim) = 0

	}

	# Initialize output image line pointer.
	call amovkl (1, IE_V(ie,1), IM_MAXDIM)

	percent = 0
	nlines = 0
	totlines = 1
	do i = 2, IM_NDIM(outim)
	    totlines = totlines * IM_LEN(outim,i)

	# Generate the pixel data for the output image line by line, 
	# evaluating the user supplied expression to produce each image
	# line.  Images may be any dimension, datatype, or size.

	# call memlog ("--------- PROCESS IMAGE -----------")

	out = NULL
	repeat {
	    # call memlog1 ("--------- line %d ----------", nlines + 1)
	    # Output image line generated by last iteration.
	    if (out != NULL) {
		op = data
		if (O_LEN(out) == 0) {
		    # Output image line is a scalar.

		    switch (O_TYPE(out)) {
		    case TY_BOOL:
			Memi[op] = O_VALI(out)

		    case TY_SHORT:
			Mems[op] = O_VALS(out)

		    case TY_INT:
			Memi[op] = O_VALI(out)

		    case TY_LONG:
			Meml[op] = O_VALL(out)

		    case TY_REAL:
			Memr[op] = O_VALR(out)

		    case TY_DOUBLE:
			Memd[op] = O_VALD(out)

		    }

		} else {
		    # Output image line is a vector.
		    npix = min (O_LEN(out), IM_LEN(outim,1))
		    ip = O_VALP(out)
		    switch (O_TYPE(out)) {
		    case TY_BOOL:
			call amovi (Memi[ip], Memi[op], npix)

		    case TY_SHORT:
			call amovs (Mems[ip], Mems[op], npix)

		    case TY_INT:
			call amovi (Memi[ip], Memi[op], npix)

		    case TY_LONG:
			call amovl (Meml[ip], Meml[op], npix)

		    case TY_REAL:
			call amovr (Memr[ip], Memr[op], npix)

		    case TY_DOUBLE:
			call amovd (Memd[ip], Memd[op], npix)

		    }
		}

		call evvfree (out)
		out = NULL
	    }

	    # Get the next line in all input images.  If EOF is seen on the
	    # image we merely rewind and keep going.  This allows a vector,
	    # plane, etc. to be applied to each line, band, etc. of a higher
	    # dimensioned image.

	    do i = 1, noperands {
		io = IE_IMOP(ie,i)
		if (IO_TYPE(io) != IMAGE || IO_IM(io) == NULL)
		    next

		mep = IO_MEF(io)
		extp = megetep(mep, i2dget(corary, i, extveri), READ_ONLY)
		IO_EXT(io) = extp
		im = EXT_EXTP(extp)
		IO_IM(io) = im
		#im = IO_IM(io)
		o  = IO_OP(io)

		# Data for a 1D image was read in above.
		if (IM_NDIM(im) == 1)
		    next
		switch (O_TYPE(o)) {

		case TY_SHORT:
		    if (mignls (extp, IO_DATA(io), IO_V(io,1)) == EOF) {
			call amovkl (1, IO_V(io,1), IM_MAXDIM)
			if (mignls (extp, IO_DATA(io), IO_V(io,1)) == EOF)
			    call error (GEERR_NODATA, s_nodata)
		    }

		case TY_INT:
		    if (mignli (extp, IO_DATA(io), IO_V(io,1)) == EOF) {
			call amovkl (1, IO_V(io,1), IM_MAXDIM)
			if (mignli (extp, IO_DATA(io), IO_V(io,1)) == EOF)
			    call error (GEERR_NODATA, s_nodata)
		    }

		case TY_LONG:
		    if (mignll (extp, IO_DATA(io), IO_V(io,1)) == EOF) {
			call amovkl (1, IO_V(io,1), IM_MAXDIM)
			if (mignll (extp, IO_DATA(io), IO_V(io,1)) == EOF)
			    call error (GEERR_NODATA, s_nodata)
		    }

		case TY_REAL:
		    if (mignlr (extp, IO_DATA(io), IO_V(io,1)) == EOF) {
			call amovkl (1, IO_V(io,1), IM_MAXDIM)
			if (mignlr (extp, IO_DATA(io), IO_V(io,1)) == EOF)
			    call error (GEERR_NODATA, s_nodata)
		    }

		case TY_DOUBLE:
		    if (mignld (extp, IO_DATA(io), IO_V(io,1)) == EOF) {
			call amovkl (1, IO_V(io,1), IM_MAXDIM)
			if (mignld (extp, IO_DATA(io), IO_V(io,1)) == EOF)
			    call error (GEERR_NODATA, s_nodata)
		    }

		default:
		    call error (GEERR_BADTYPE, s_badtype)
		}

	    }

	    # call memlog (".......... enter evvexpr ..........")

	    # This is it!  Evaluate the vector expression.
	    flags = 0
	    if (rangecheck)
		flags = or (flags, EV_RNGCHK)
#REMOVE
#call printf("it's here \n")
#call flush(STDOUT)

    	    out = evvexpr (Memc[expr],
		    locpr(ie_getop), ie, locpr(ie_fcn), ie, flags)

#REMOVE
#call printf("it's there \n")
#call flush(STDOUT)
	    # call memlog (".......... exit evvexpr ..........")

	    # If the pixel type and line length of the output image are
	    # still undetermined set them to match the output operand.

	    if (IM_PIXTYPE(outim) == 0) {
		if (IE_OUTTYPE(ie) == 0) {
		    if (O_TYPE(out) == TY_BOOL)
			IE_OUTTYPE(ie) = TY_INT
		    else
			IE_OUTTYPE(ie) = O_TYPE(out)
		    IM_PIXTYPE(outim) = IE_OUTTYPE(ie)
		} else
		    IM_PIXTYPE(outim) = IE_OUTTYPE(ie)
	    }
	    if (IM_LEN(outim,1) == 0) {
		if (IE_AXLEN(ie,1) == 0) {
		    if (O_LEN(out) == 0) {
			IE_AXLEN(ie,1) = 1
			IM_LEN(outim,1) = 1
		    } else {
			IE_AXLEN(ie,1) = O_LEN(out)
			IM_LEN(outim,1) = O_LEN(out)
		    }
		} else
		    IM_LEN(outim,1) = IE_AXLEN(ie,1)
	    }

	    # Print percent done.
	    if (ldebug) {
		        nlines = nlines + 1
		        if (nlines * 100 / totlines >= percent + 10) {
		            percent = percent + 10
		            call printf ("%2d%% ")
			        call pargi (percent)
		            call flush (STDOUT)
		        }
	    }

	    switch (O_TYPE(out)) {
	    case TY_BOOL:
		status = impnli (outim, data, IE_V(ie,1))

	    case TY_SHORT:
		status = impnls (outim, data, IE_V(ie,1))

	    case TY_INT:
		status = impnli (outim, data, IE_V(ie,1))

	    case TY_LONG:
		status = impnll (outim, data, IE_V(ie,1))

	    case TY_REAL:
		status = impnlr (outim, data, IE_V(ie,1))

	    case TY_DOUBLE:
		status = impnld (outim, data, IE_V(ie,1))

	    default:
		call error (GEERR_INCOMPATTYPE, "expression type incompatible with image")
	    }
	} until (status == EOF)
		# newline after percentage
	if (ldebug) { call printf("\n") }
	    if (ldebug) {
		    call printf("... Finished Processing extveri=#%d\n")
		    call pargi(extveri)
		    call flush(STDOUT)
	    }

	# clean up for next pass
	#  Unmap images...
		call imunmap (outim)
        outim = 0
		do i = 1, noperands {
			io = IE_IMOP(ie,i)
			if (IO_TYPE(io) == IMAGE && IO_EXT(io) != NULL)
				call meepfunmap (IO_EXT(io), true)
		}
		
	}  # end of extver iteration loop
	
	# call memlog ("--------- DONE PROCESSING IMAGE -----------")

	if (ldebug) {
	    call printf ("- done\n")
	    call flush (STDOUT)
	}
 

	# Clean up.
	do i = 1, noperands {
	    io = IE_IMOP(ie,i)
	    o = IO_OP(io)
	    if (O_TYPE(o) == TY_CHAR)
		call mfree (O_VALP(o), TY_CHAR)
	}
    
        
    call evvfree (out)
	call mfree (expr, TY_CHAR)
	if (st != NULL)
	    call stclose (st)

	call clpstr ("lastout", Memc[output])
	call sfree (sp)
end


# IE_GETOP -- Called by evvexpr to fetch an input image operand.

procedure ie_getop (ie, opname, o)

pointer	ie			#I imexpr descriptor
char	opname[ARB]		#I operand name
pointer	o			#I output operand to be filled in

int	axis, i
pointer	param, data
pointer	sp, im, io, v, mep
char tmpstr[SZ_LINE]
bool	imgetb()
int	imgeti()
double	imgetd()
int	imgftype(), btoi(), megftype()
errchk	malloc
bool ldebug
double megetd()
int megeti()
bool megetb()
bool bphu

define	err_ 91
#include "glogcommon.h"
begin
    bphu = false
	call smark (sp)
    ldebug = false
    
    if (ldebug) {
        call sprintf(tmpstr, SZ_LINE, "opname %s\n")
        call pargstr(opname)
        call printf(tmpstr)
    }
    
	if (IS_LOWER(opname[1]) && opname[2] == EOS) {
	    # Image operand.
#REMOVE
#call printf("it's in a %s\n")
#call pargstr(opname)
#call flush(STDOUT)

	    io = NULL
	    do i = 1, IE_NOPERANDS(ie) {
		io = IE_IMOP(ie,i)
		if (IO_OPNAME(io) == opname[1])
		    break
		io = NULL
	    }

	    if (io == NULL)
		goto err_
	    else
		v = IO_OP(io)

	    call amovi (Memi[v], Memi[o], LEN_OPERAND)
	    if (IO_TYPE(io) == IMAGE) {
		O_VALP(o) = IO_DATA(io)
		O_FLAGS(o) = 0
	    }

	    call sfree (sp)
	    return

	} else if (IS_LOWER(opname[1]) && (opname[2] == '.' || 
                (opname[2]  == '[' && opname[3] == '0' && 
                        opname[4] == ']' && opname[5] == '.'))) {
	    # Image parameter reference, e.g., "a.foo" or phu parm (a[0].foo)
        if (opname[2]  == '[' && opname[3] == '0' && 
                opname[4] == ']' && opname[5] == '.') {
            # the parameter is in the phu, true
            bphu = true
        }
	    call salloc (param, SZ_FNAME, TY_CHAR)
#REMOVE
#call printf("it's in here B\n")
#call flush(STDOUT)

	    # Locate referenced symbolic image operand (e.g. "a").
	    io = NULL
	    do i = 1, IE_NOPERANDS(ie) {
		io = IE_IMOP(ie,i)
		if (IO_OPNAME(io) == opname[1] && IO_TYPE(io) == IMAGE)
		    break
		io = NULL
	    }
	    if (io == NULL)
		goto err_

	    # Get the parameter value and set up operand struct.
        if (bphu) {
            call strcpy (opname[6], Memc[param], SZ_FNAME)
        } else {
            call strcpy (opname[3], Memc[param], SZ_FNAME)
        }
	    im = IO_IM(io)
        mep = IO_MEF(io)
            

        if (bphu) {
            iferr (O_TYPE(o) = megftype (mep, Memc[param])) {
                goto err_
            }
        } else {
            iferr (O_TYPE(o) = imgftype (im, Memc[param])) {
                goto err_
            }
        }

	    switch (O_TYPE(o)) {
	    case TY_BOOL:
            if (bphu) {
                iferr (O_VALI(o) = btoi (megetb(mep,Memc[param]))) {
                    goto err_
                }
            } else {
                iferr (O_VALI(o) = btoi (imgetb (im, Memc[param]))) {
		            goto err_
                }
            }

	    case TY_CHAR:
		O_LEN(o) = SZ_LINE
		O_FLAGS(o) = O_FREEVAL
		iferr {
                if (bphu) {
                    #not supported from phu syntax
                    call error(1,"string parms from phu not supported")
                } else {
                    call malloc (O_VALP(o), SZ_LINE, TY_CHAR)
		            call imgstr (im, Memc[param], O_VALC(o), SZ_LINE)
                }
		    } then {
		        goto err_
            }
	    case TY_INT:
            if (bphu) {
                iferr (O_VALI(o) = megeti (mep, Memc[param])) {
                    goto err_
                }
            } else {
		        iferr (O_VALI(o) = imgeti (im, Memc[param])) {
		            goto err_
                }
            } 
	    case TY_REAL:
		    O_TYPE(o) = TY_DOUBLE
            
            if (bphu) {
                iferr(O_VALD(o) = megetd (mep, Memc[param])) {
                    goto err_
                }
            } else {
                iferr (O_VALD(o) = imgetd (im, Memc[param])) {
                    goto err_
                }
            }

	    default:
		goto err_
	    }

	    call sfree (sp)
	    return

	} else if (IS_UPPER(opname[1]) && opname[2] == EOS) {

#REMOVE
#call printf("it's in here C\n")
#call flush(STDOUT)

		# The current pixel coordinate [I,J,K,...].  The line coordinate

	    # is a special case since the image is computed a line at a time.
	    # If "I" is requested return a vector where v[i] = i.  For J, K,
	    # etc. just return the scalar index value.
	    axis = opname[1] - 'I' + 1
	    if (axis == 1) {
		O_TYPE(o) = TY_INT
		if (IE_AXLEN(ie,1) > 0)
		    O_LEN(o) = IE_AXLEN(ie,1)
		else {
		    # Line length not known yet.
		    O_LEN(o) = DEF_LINELEN
		}
		call malloc (data, O_LEN(o), TY_INT)
		do i = 1, O_LEN(o)
		    Memi[data+i-1] = i
		O_VALP(o) = data
		O_FLAGS(o) = O_FREEVAL
	    } else {
		O_TYPE(o) = TY_INT
		#O_LEN(o) = 0
		#if (axis < 1 || axis > IM_MAXDIM)
		    #O_VALI(o) = 1
		#else
		    #O_VALI(o) = IE_V(ie,axis)
		#O_FLAGS(o) = 0
		if (IE_AXLEN(ie,1) > 0)
		    O_LEN(o) = IE_AXLEN(ie,1)
		else 
		    # Line length not known yet.
		    O_LEN(o) = DEF_LINELEN
		call malloc (data, O_LEN(o), TY_INT)
		if (axis < 1 || axis > IM_MAXDIM)
		    call amovki (1, Memi[data], O_LEN(o))
		else
		    call amovki (IE_V(ie,axis), Memi[data], O_LEN(o))
		O_VALP(o) = data
		O_FLAGS(o) = O_FREEVAL
	    }

	    call sfree (sp)
	    return
	}

err_
	O_TYPE(o) = ERR
	call sfree (sp)
end


# IE_FCN -- Called by evvexpr to execute an imexpr special function.

procedure ie_fcn (ie, fcn, args, nargs, o)

pointer	ie			#I imexpr descriptor
char	fcn[ARB]		#I function name
pointer	args[ARB]		#I input arguments
int	nargs			#I number of input arguments
pointer	o			#I output operand to be filled in

begin
	# No functions yet.
	O_TYPE(o) = ERR
end


# IE_GETEXPRDB -- Read the expression database into a symbol table.  The
# input file has the following structure:
#
#	<symbol>['(' arg-list ')'][':'|'=']	replacement-text
#		
# Symbols must be at the beginning of a line.  The expression text is
# terminated by a nonempty, noncomment line with no leading whitespace.

pointer procedure ie_getexprdb (fname)

char	fname[ARB]		#I file to be read

pointer	sym, sp, lbuf, st, a_st, ip, symname, tokbuf, text
int	tok, fd, line, nargs, op, token, buflen, offset, stpos, n
errchk	open, getlline, stopen, stenter, ie_puttok
int	open(), getlline(), ctotok(), stpstr()
pointer	stopen(), stenter()
define	skip_ 91

begin
	call smark (sp)
	call salloc (lbuf, SZ_COMMAND, TY_CHAR)
	call salloc (text, SZ_COMMAND, TY_CHAR)
	call salloc (tokbuf, SZ_COMMAND, TY_CHAR)
	call salloc (symname, SZ_FNAME, TY_CHAR)

	fd = open (fname, READ_ONLY, TEXT_FILE)
	st = stopen ("imexpr", DEF_LENINDEX, DEF_LENSTAB, DEF_LENSBUF)
	a_st = stopen ("args", DEF_LENINDEX, DEF_LENSTAB, DEF_LENSBUF)
	line = 0

	while (getlline (fd, Memc[lbuf], SZ_COMMAND) != EOF) {
	    line = line + 1
	    ip = lbuf

	    # Skip comments and blank lines.
	    while (IS_WHITE(Memc[ip]))
		ip = ip + 1
	    if (Memc[ip] == '\n' || Memc[ip] == '#')
		next
	
	    # Get symbol name.
	    if (ctotok (Memc,ip,Memc[symname],SZ_FNAME) != TOK_IDENTIFIER) {
		call eprintf ("exprdb: expected identifier at line %d\n")
		    call pargi (line)
skip_		while (getlline (fd, Memc[lbuf], SZ_COMMAND) != EOF) {
		    line = line + 1
		    if (Memc[lbuf] == '\n')
			break
		}
	    }

	    call stmark (a_st, stpos)

	    # Check for the optional argument-symbol list.  Allow only a
	    # single space between the symbol name and its argument list,
	    # otherwise we can't tell the difference between an argument
	    # list and the parenthesized expression which follows.

	    if (Memc[ip] == ' ')
		ip = ip + 1

	    if (Memc[ip] == '(') {
		ip = ip + 1
		n = 0
		repeat {
		    tok = ctotok (Memc, ip, Memc[tokbuf], SZ_FNAME)
		    if (tok == TOK_IDENTIFIER) {
			sym = stenter (a_st, Memc[tokbuf], LEN_ARGSYM)
			n = n + 1
			ARGNO(sym) = n
		    } else if (Memc[tokbuf] == ',') {
			;
		    } else if (Memc[tokbuf] != ')') {
			call eprintf ("exprdb: bad arglist at line %d\n")
			    call pargi (line)
			call stfree (a_st, stpos)
			goto skip_
		    }
		} until (Memc[tokbuf] == ')')
	    }

	    # Check for the optional ":" or "=".
	    while (IS_WHITE(Memc[ip]))
		ip = ip + 1
	    if (Memc[ip] == ':' || Memc[ip] == '=')
		ip = ip + 1

	    # Accumulate the expression text.
	    buflen = SZ_COMMAND
	    op = 1
	    
	    repeat {
		repeat {
		    token = ctotok (Memc, ip, Memc[tokbuf], SZ_COMMAND)
		    if (Memc[tokbuf] == '#')
			break
		    else if (token != TOK_EOS && token != TOK_NEWLINE)
			call ie_puttok (a_st, text, op, buflen, Memc[tokbuf])
		} until (token == TOK_EOS)

		if (getlline (fd, Memc[lbuf], SZ_COMMAND) == EOF)
		    break
		else
		    line = line + 1

		for (ip=lbuf;  IS_WHITE(Memc[ip]);  ip=ip+1)
		    ;
		if (ip == lbuf) {
		    call ungetline (fd, Memc[lbuf])
		    line = line - 1
		    break
		}
	    }

	    # Free any argument list symbols.
	    call stfree (a_st, stpos)

	    # Scan the expression text and count the number of $N arguments.
	    nargs = 0
	    for (ip=text;  Memc[ip] != EOS;  ip=ip+1)
		if (Memc[ip] == '$' && IS_DIGIT(Memc[ip+1])) {
		    nargs = max (nargs, TO_INTEG(Memc[ip+1]))
		    ip = ip + 1
		}

	    # Enter symbol in table.
	    sym = stenter (st, Memc[symname], LEN_SYM)
	    offset = stpstr (st, Memc[text], 0)
	    SYM_TEXT(sym) = offset
	    SYM_NARGS(sym) = nargs
	}

	call stclose (a_st)
	call sfree (sp)

	return (st)
end


# IE_PUTTOK -- Append a token string to a text buffer.

procedure ie_puttok (a_st, text, op, buflen, token)

pointer	a_st			#I argument-symbol table
pointer	text			#U text buffer
int	op			#U output pointer
int	buflen			#U buffer length, chars
char	token[ARB]		#I token string

pointer	sym
int	ip, ch1, ch2
pointer	stfind()
errchk	realloc

begin
	# Replace any symbolic arguments by "$N".
	if (a_st != NULL && IS_ALPHA(token[1])) {
	    sym = stfind (a_st, token)
	    if (sym != NULL) {
		token[1] = '$'
		token[2] = TO_DIGIT(ARGNO(sym))
		token[3] = EOS
	    }
	}

	# Append the token string to the text buffer.
	for (ip=1;  token[ip] != EOS;  ip=ip+1) {
	    if (op + 1 > buflen) {
		buflen = buflen + SZ_COMMAND
		call realloc (text, buflen, TY_CHAR)
	    }

	    # The following is necessary because ctotok parses tokens such as
	    # "$N", "==", "!=", etc.  as two tokens.  We need to rejoin these
	    # characters to make one token.

	    if (op > 1 && token[ip+1] == EOS) {
		ch1 = Memc[text+op-3]
		ch2 = token[ip]

		if (ch1 == '$' && IS_DIGIT(ch2))
		    op = op - 1
		else if (ch1 == '*' && ch2 == '*')
		    op = op - 1
		else if (ch1 == '/' && ch2 == '/')
		    op = op - 1
		else if (ch1 == '<' && ch2 == '=')
		    op = op - 1
		else if (ch1 == '>' && ch2 == '=')
		    op = op - 1
		else if (ch1 == '=' && ch2 == '=')
		    op = op - 1
		else if (ch1 == '!' && ch2 == '=')
		    op = op - 1
		else if (ch1 == '?' && ch2 == '=')
		    op = op - 1
		else if (ch1 == '&' && ch2 == '&')
		    op = op - 1
		else if (ch1 == '|' && ch2 == '|')
		    op = op - 1
	    }

	    Memc[text+op-1] = token[ip]
	    op = op + 1
	}

	# Append a space to ensure that tokens are delimited.
	Memc[text+op-1] = ' '
	op = op + 1

	Memc[text+op-1] = EOS
end



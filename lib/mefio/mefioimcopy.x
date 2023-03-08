# Copyright(c) 1986-2005 Association of Universities for Research in Astronomy Inc.

include	<imhdr.h>
include <time.h>
include "mefio.h"


# MIMCOPY -- Copy an image.  Use sequential routines to permit copying
# images of any dimension.  Perform pixel i/o in the datatype of the image,
# to avoid unnecessary type conversion.

procedure mimcopy (image1, image2, verbose)

char	image1[ARB]			# Input image
char	image2[ARB]			# Output image
bool	verbose				# Print the operation

int	npix, junk
pointer	buf1, buf2, im1, im2
pointer	sp, root1, root2, imtemp, section
long	v1[IM_MAXDIM], v2[IM_MAXDIM]
char tmpstr[SZ_LINE]
long tmpl,tmpu
int tmpi
int tm[LEN_TMSTRUCT]

bool	strne()
int	imgnls(), imgnll(), imgnlr(), imgnld(), imgnlx()
int	impnls(), impnll(), impnlr(), impnld(), impnlx()
pointer	immap()
long clktime(), lsttogmt()

long tmp

bool ldebug

int errget()

errchk immap()

begin
	ldebug = false
	call smark (sp)
	call salloc (root1, SZ_PATHNAME, TY_CHAR)
	call salloc (root2, SZ_PATHNAME, TY_CHAR)
	call salloc (imtemp, SZ_PATHNAME, TY_CHAR)
	call salloc (section, SZ_FNAME, TY_CHAR)

	# If verbose print the operation.
	if (ldebug) {
	    call printf ("%s -> %s\n")
		call pargstr (image1)
		call pargstr (image2)
	    call flush (STDOUT)
	}

	# Get the input and output root names and the output section.
	call imgimage (image1, Memc[root1], SZ_PATHNAME)
	call imgimage (image2, Memc[root2], SZ_PATHNAME)
	call imgsection (image2, Memc[section], SZ_FNAME)

	if (ldebug) {
	    call printf ("r: %s -> %s\n")
		call pargstr (Memc[root1])
		call pargstr (Memc[root2])
	    call flush (STDOUT)
	}
	
	
    if (ldebug) {
        call printf("Parsed Names\n")
        call flush(STDOUT)
    }
    
	# Map the input image.
	im1 = immap (image1, READ_ONLY, 0)

    if (ldebug) {
        call printf("mapped im1\n")
        call flush(STDOUT)
    }
    
	# If the output has a section appended we are writing to a
	# section of an existing image.  Otherwise get a temporary
	# output image name and map it as a copy of the input image.
	# Copy the input image to the temporary output image and unmap
	# the images.  Release the temporary image name.

	if (strne (Memc[root1], Memc[root2]) && Memc[section] != EOS) {
	    call strcpy (image2, Memc[imtemp], SZ_PATHNAME)
        if (ldebug) {
            call printf("w/section image2 = %s\n")
            call pargstr(image2)
            call flush(STDOUT)
        }
	    im2 = immap (image2, READ_WRITE, 0)
	} else {
	    call xt_mkimtemp (image1, image2, Memc[imtemp], SZ_PATHNAME)
        if (ldebug) {
            call printf("wo/section image2 = %s\n")
            call pargstr(image2)
            call flush(STDOUT)
        }
		
	    iferr (im2 = immap (image2, NEW_COPY, im1)){
			tmpi = errget(tmpstr, SZ_LINE)
			if (ldebug) {
				call printf("err #%d, %s\n")
				call pargi(tmpi)
				call pargstr(tmpstr)
				call flush(STDOUT)
			}
			call error(MEERR_COPYFAILED, tmpstr)
		}		
	}
    
    if (ldebug) {
        call printf("mapped im2\n")
        call flush(STDOUT)
    }

    # put modification time in copy (note... the way this is used 
    # in mefio/gemarith/gemexpr this becomes the Time of Last Modification in the PHU
    #tmpl = 0
    tmp = 0
    while (tmpl < 400000) {
        tmpl = clktime(tmpl)
        tmpu = lsttogmt(tmpl)
#call cnvtime(tmpl, tmpstr, SZ_LINE)
        call brktime(tmpu, tm)
        call sprintf(tmpstr,SZ_LINE, "%4d-%02d-%02dT%02d:%02d:%02d")
        call pargi(TM_YEAR(tm))
        call pargi(TM_MONTH(tm))
        call pargi(TM_MDAY(tm))
        call pargi(TM_HOUR(tm))
        call pargi(TM_MIN(tm))
        call pargi(TM_SEC(tm))
        
        if (ldebug) {
            call printf("time string: %s\n(%d)(%d)\n")
            call pargstr(tmpstr)
            call pargi(tmpl)
            call pargi(tmpu)
            call flush(STDOUT)
        }

    }
    #call imastr(im2, "GEM-TLM", tmpstr) 
    call imastrc(im2, "GEM-TLM", tmpstr, "UT Last modification with GEMINI")
    
	# Setup start vector for sequential reads and writes.

	call amovkl (long(1), v1, IM_MAXDIM)
	call amovkl (long(1), v2, IM_MAXDIM)

	# Copy the image.

	npix = IM_LEN(im1, 1)
	switch (IM_PIXTYPE(im1)) {
	case TY_SHORT:
	    while (imgnls (im1, buf1, v1) != EOF) {
		junk = impnls (im2, buf2, v2)
		call amovs (Mems[buf1], Mems[buf2], npix)
	    }
	case TY_USHORT, TY_INT, TY_LONG:
	    while (imgnll (im1, buf1, v1) != EOF) {
		junk = impnll (im2, buf2, v2)
		call amovl (Meml[buf1], Meml[buf2], npix)
	    }
	case TY_REAL:
	    while (imgnlr (im1, buf1, v1) != EOF) {
		junk = impnlr (im2, buf2, v2)
		call amovr (Memr[buf1], Memr[buf2], npix)
	    }
	case TY_DOUBLE:
	    while (imgnld (im1, buf1, v1) != EOF) {
		junk = impnld (im2, buf2, v2)
		call amovd (Memd[buf1], Memd[buf2], npix)
	    }
	case TY_COMPLEX:
	    while (imgnlx (im1, buf1, v1) != EOF) {
	        junk = impnlx (im2, buf2, v2)
		call amovx (Memx[buf1], Memx[buf2], npix)
	    }
	default:
	    call error (1, "unknown pixel datatype")
	}

	# Unmap the images.

	call imunmap (im2)
	call imunmap (im1)
	call xt_delimtemp (image2, Memc[imtemp])
	call sfree (sp)
end

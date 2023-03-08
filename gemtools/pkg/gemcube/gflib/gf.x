# Copyright(c) 2006 Association of Universities for Research in Astronomy, Inc.
# GF -- Interface routines to the geometry function drivers.


include	"gf.h"


# GF_OPEN -- Geometry function open procedure.

procedure gf_open (geofunc, gf, images, index)

char	geofunc[ARB]			#I Function
pointer	gf				#O Function pointer
pointer	images				#I List of images
int	index				#I Index of image to open

int	i
pointer	ptr1, ptr2

bool	streq()

include	"gf.com"

begin
	gf = NULL

	# Initialize geometry function table.
	call gf_load ()

	# Find function pointer.  Note that we assume that once the
	# function table is initialized the function table pointer
	# is never changed; # i.e. reallocated.

	ptr1 = NULL
	do i = 0, gf_n-1 {
	    ptr2 = gfs + i * GF_LEN
	    if (streq (geofunc, GF_NAME(ptr2))) {
		ptr1 = ptr2
		break
	    }
	}

	if (ptr1 == NULL)
	    call error (1, "Geometry function not found")

	# Open function.
	call zcall4 (GF_OPEN(ptr1), ptr2, images, index, geofunc)

	# Allocate return structure.
	call malloc (gf, 2, TY_STRUCT)
	Memi[gf] = ptr1
	Memi[gf+1] = ptr2
end


# GF_OUT -- Geometry function output procedure.

procedure gf_out (geofunc, images, refmw, mw, imlen)

char	geofunc[ARB]			#I Function
pointer	images[ARB]			#I List of images
pointer	refmw				#I Reference WCS
pointer	mw				#O Suggested WCS
int	imlen[3]			#0 Suggested size

int	i
pointer	ptr1, ptr2

bool	streq()

include	"gf.com"

begin
	# Initialize geometry function table.
	call gf_load ()

	# Find function pointer.  Note that we assume that once the
	# function table is initialized the function table pointer
	# is never changed; # i.e. reallocated.

	ptr1 = NULL
	do i = 0, gf_n-1 {
	    ptr2 = gfs + i * GF_LEN
	    if (streq (geofunc, GF_NAME(ptr2))) {
		ptr1 = ptr2
		break
	    }
	}

	if (ptr1 == NULL)
	    call error (1, "Geometry function not found")

	# Open function.
	call zcall5 (GF_OUT(ptr1), images, refmw, geofunc, mw, imlen)
end


# GF_CLOSE -- Geometry function close procedure.

procedure gf_close (gf)

pointer	gf				#U Function pointer

begin
	if (gf == NULL)
	    return
	call zcall1 (GF_CLOSE(Memi[gf]), Memi[gf+1])
	call mfree (gf, TY_STRUCT)
end


# GF_PIXEL -- Geometry function pixel procedure.

procedure gf_pixel (gf, pixel, world)

pointer	gf				#I Function pointer
double	pixel[3]			#I Pixel coordinate
double	world[3]			#0 World coordinate

begin
	call zcall3 (GF_PIXEL(Memi[gf]), Memi[gf+1], pixel, world)
end


# GF_GEOM -- Geometry function geometry procedure.

procedure gf_geom (gf, pixel, shape, axmap)

pointer	gf				#I Function pointer
double	pixel[3]			#I Pixel coordinate
char	shape[1024]			#0 Shape string
int	axmap[3]			#O Axis map

begin
	call zcall4 (GF_GEOM(Memi[gf]), Memi[gf+1], pixel, shape, axmap)
end


# GF_DEFINE -- Define a function in the function table structure.

procedure gf_define (geofunc, open, out, close, pixel, geom)

char	geofunc[ARB]			#I Name
int	open				#I Open procedure
int	out				#I Out procedure
int	close				#I Close procedure
int	pixel				#I Pixel procedure
int	geom				#I Geometry procedure

int	n
data	n /0/
pointer	gf
include	"gf.com"

begin
	# Get a new slot.
	if (n == 0) 
	    call malloc (gfs, GF_NALLOC*GF_LEN, TY_STRUCT)
	else if (mod (n, GF_NALLOC) == 0)
	    call realloc (gfs, (n+GF_NALLOC)*GF_LEN, TY_STRUCT)
	gf = gfs + n * GF_LEN
	n = n + 1

	# Load the function.
	call strcpy (geofunc, GF_NAME(gf), GF_LENNAME)
	GF_OPEN(gf) = open
	GF_OUT(gf) = out
	GF_CLOSE(gf) = close
	GF_PIXEL(gf) = pixel
	GF_GEOM(gf) = geom
	gf_n = n
end


# GF_LIST -- List geometry functions.

procedure gf_list (fd)

int	fd				#I File descriptor

int	i
pointer	gf

include	"gf.com"

begin
	call gf_load ()

	do i = 0, gf_n-1 {
	    gf = gfs + i * GF_LEN
	    call fprintf (fd, "%s\n")
	        call pargstr (GF_NAME(gf))
	}
	call flush (fd)
end

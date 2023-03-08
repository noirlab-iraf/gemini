# Copyright(c) 2006 Association of Universities for Research in Astronomy, Inc.

include	<fset.h>


# T_TRANSCUBE -- Cube transform procedure.

procedure t_transcube ()

pointer	input				# Input image list
pointer	output				# Output image list
pointer	masks				# Output mask list
pointer	weights				# Output weights list
pointer	logfiles			# Output logfile list
pointer	bpm				# Input bad pixel list
pointer	scale				# Input scale list
pointer	wt				# Input weight list
pointer	wcsreference			# WCS reference
pointer	wttype				# Weighting type
real	drizscale[3]			# Drizzle scale factor
real	blank				# Blank value
pointer	geofunc				# Geometry function
real	memalloc			# Memory to alloc (Mb)

int	i, nlogs
pointer	sp, str, logs 

bool	streq()
int	imtlen(), imtgetim(), open(), nscan()
real	clgetr()
pointer	imtopenp()

begin
	call smark (sp)
	call salloc (wcsreference, SZ_FNAME, TY_CHAR)
	call salloc (wttype, SZ_FNAME, TY_CHAR)
	call salloc (geofunc, SZ_FNAME, TY_CHAR)
	call salloc (str, SZ_LINE, TY_CHAR)

	input = imtopenp ("input")
	output = imtopenp ("output")
	masks = imtopenp ("masks")
	weights = imtopenp ("weights")
	logfiles = imtopenp ("logfiles")
	bpm = imtopenp ("bpm")
	scale = imtopenp ("scale")
	wt = imtopenp ("wt")
	call clgstr ("wcsreference", Memc[wcsreference], SZ_FNAME)
	call clgstr ("wttype", Memc[wttype], SZ_FNAME)
	call clgstr ("drizscale", Memc[str], SZ_FNAME)
	blank = clgetr ("blank")
	call clgstr ("geofunc", Memc[geofunc], SZ_FNAME)
	memalloc = clgetr ("memalloc")

	# Parse the drizscale string.
	call sscan (Memc[str])
	call gargr (drizscale[1])
	call gargr (drizscale[2])
	call gargr (drizscale[3])
	switch (nscan()) {
	case 0:
	    call error (1, "Syntax error in drizscale parameter")
	case 1:
	    drizscale[2] = drizscale[1]
	    drizscale[3] = drizscale[2]
	case 2:
	    drizscale[3] = drizscale[2]
	}

	# Open logfiles.
	nlogs = imtlen (logfiles)
	if (nlogs > 0) {
	    call salloc (logs, nlogs, TY_INT)
	    nlogs = 0
	    while (imtgetim (logfiles, Memc[str], SZ_LINE) != EOF) {
		Memi[logs+nlogs] = open (Memc[str], APPEND, TEXT_FILE)
		nlogs = nlogs + 1

		call fseti (Memi[logs+nlogs-1], F_FLUSHNL, YES)
		call sysid (Memc[str], SZ_LINE)
		call fprintf (Memi[logs+nlogs-1], "\n%s: %s\n")
		    call pargstr ("GEMCUBE")
		    call pargstr (Memc[str])
	    }
	}

	# Check whether to just list the geometry functions.
	if (streq ("list", Memc[geofunc]))
	    call gf_list (STDOUT)
	else
	    call transcube_list (input, output, masks, weights, bpm, scale, wt,
	        Memc[wcsreference], Memc[wttype], drizscale, blank,
		Memc[geofunc], memalloc, Memi[logs], nlogs)

	do i = 0, nlogs-1
	    call close (Memi[logs+i])
	call imtclose (wt)
	call imtclose (scale)
	call imtclose (bpm)
	call imtclose (logfiles)
	call imtclose (weights)
	call imtclose (masks)
	call imtclose (output)
	call imtclose (input)

	call sfree (sp)
end

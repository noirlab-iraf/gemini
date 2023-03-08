# Copyright(c) 2006 Association of Universities for Research in Astronomy, Inc.
# GF_LOAD -- Load geometry functions.

procedure gf_load ()

# Add driver functions here.
extern	gfwcs_open(),
	gfwcs_out(),
	gfwcs_close(),
	gfwcs_pixel(),
	gfwcs_geom()
extern	gfstack_open(),
	gfstack_out()

bool	first_time
data	first_time /true/
errchk	gf_define
int	locpr()

begin
	# Only do this the first time.
	if (!first_time)
	    return
	first_time = false

	# Define driver functions here.
	call gf_define ("gfwcs",
	    locpr(gfwcs_open),
	    locpr(gfwcs_out),
	    locpr(gfwcs_close),
	    locpr (gfwcs_pixel),
	    locpr (gfwcs_geom))
	call gf_define ("gfstack",
	    locpr(gfstack_open),
	    locpr(gfstack_out),
	    locpr(gfwcs_close),
	    locpr (gfwcs_pixel),
	    locpr (gfwcs_geom))
end

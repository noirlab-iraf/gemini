# Copyright(c) 2006 Association of Universities for Research in Astronomy, Inc.
# GF common.

int	gf_n				# Number of functions
pointer	gfs				# Function table structure

common	/gfcom/ gf_n, gfs

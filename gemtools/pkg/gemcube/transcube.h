# Copyright(c) 2006 Association of Universities for Research in Astronomy, Inc.


define	WTTYPES	"|nearest|drizzle|linear|"
define	WT_NEAREST	1		# Nearest pixel weighting
define	WT_DRIZ		2		# Drizzle weighting
define	WT_LIN		3		# Linear interpolation weighting
define	WT_ERR		99		# Error code for weight type error

define	SHAPES	"|rectangle|cylinder|hexagon|"

#define	NSAMPLEX	21
#define	NSAMPLEY	21
#define	NSAMPLEZ	21
define	NSAMPLEX	5
define	NSAMPLEY	5
define	NSAMPLEZ	5

define	NSUBPIX		100		# Number of subpixels

define	TOL1		1E-3		# Tolerance on equality of pixel size
define	TOL2		1		# Tolerance on equality ofposition angle

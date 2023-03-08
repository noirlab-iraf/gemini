# Copyright(c) 2006 Association of Universities for Research in Astronomy, Inc.
#
# Function table structure definitions.
define	GF_LENNAME	31		#I Maximum length of function name
define	GF_NALLOC	10		#I Number of functions per realloc
define	GF_LEN		21		#I Structure length per function
define	GF_NAME		Memc[P2C($1)]	#I Function name
define	GF_OPEN		Memi[$1+16]	#I Open function address
define	GF_OUT		Memi[$1+17]	#I Open function address
define	GF_CLOSE	Memi[$1+18]	#I Close function address
define	GF_PIXEL	Memi[$1+19]	#I Pixel function address
define	GF_GEOM		Memi[$1+20]	#I Pixel function address

# Copyright(c) 2004-2006 Association of Universities for Research in Astronomy, Inc.
#
# See imexplode.x for main documentation

# In the case of foo/bar.baz[2/3][SCI,1,goo=gar,gaz][1,2]
#
# Memc[IMX_PATH(imx)] = "foo/"
# Memc[IMX_FILE(imx)] = "bar"
# Memc[IMX_EXTENSION(imx)] = "baz"
# IMX_CLINDEX(imx) = 2
# IMX_CLSIZE(imx) = 3
# Memc[IMX_EXTNAME(imx)] = "SCI"
# IMX_EXTVERSION(imx) = 1
# Memc[IMX_IKPARAMS(imx)] = "goo=gar,gaz"
# Memc[IMX_SECTION(imx)] = "[1,2]"

define	LEN_IMEXPLODE 9

define	IMX_PATH	Memi[$1]
define	IMX_FILE	Memi[$1+1]
define	IMX_EXTENSION	Memi[$1+2]
define	IMX_CLINDEX	Memi[$1+3]
define	IMX_CLSIZE	Memi[$1+4]
define	IMX_EXTNAME	Memi[$1+5]
define	IMX_EXTVERSION	Memi[$1+6]
define	IMX_IKPARAMS	Memi[$1+7]
define	IMX_SECTION	Memi[$1+8]

define	EXTNAME		"extname"
define	EXTVERSION	"extver"

define	NO_INDEX	-1

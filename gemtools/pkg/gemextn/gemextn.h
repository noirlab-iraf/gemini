# Copyright(c) 2004-2006 Association of Universities for Research in Astronomy, Inc.
#

define  NO_MATCH	0
define	DEBUG		true	# Select debug info
define	MAX_RANGES	100

define	OMIT		"|path|extension|index|name|version|params|kernel|section|"
define	OMIT_PATH	1	# File path
define	OMIT_EXTENSION	2	# File extension (".fits" etc)
define	OMIT_INDEX	3	# Extension index
define	OMIT_NAME	4	# Extension index
define	OMIT_VERSION	5	# Extension version
define	OMIT_PARAMS	6	# Kernel parameters
define	OMIT_KERNEL	7	# Entire kernel section
define	OMIT_SECTION	8	# Image section

define	IMG_CHECKS	"|absent|exists|mef|write|"
define	EXT_CHECKS	"|absent|empty|exists|image|table|"
define	AUTO_CHECKS	"|absent|empty|exists|image|mef|table|write|force|"
define	IMG_ABSENT	1	# Image does not exist
define	IMG_EXISTS	2	# Image exists
define	IMG_MEF		3	# EXTEND=T in primary header
define	IMG_WRITE	4	# Write access to image possible
define	EXT_ABSENT	1	# Extension does not exist
define	EXT_EMPTY	2	# Extension exists, but contains no data
define	EXT_EXISTS	3	# Extension exists
define	EXT_IMAGE	4	# Extension has XTENSION=IMAGE
define	EXT_TABLE	5	# Extension has XTENSION=BINTABLE
define	AUTO_ABSENT	1	# Extension does not exist
define	AUTO_EMPTY	2	# Extension exists, but contains no data
define	AUTO_EXISTS	3	# Extension exists
define	AUTO_IMAGE	4	# Extension has XTENSION=IMAGE
define	AUTO_MEF	5	# EXTEND=T in primary header
define	AUTO_TABLE	6	# Extension has XTENSION=BINTABLE
define	AUTO_WRITE	7	# Write access to image possible
define	AUTO_FORCE	8	# Omit error checking for inconsistent checks

define	LEN_GXN		15

define	GXN_P_IMG_CHECK		Memi[$1]
define	GXN_P_EXT_CHECK		Memi[$1+1]
define	GXN_P_AUTO_CHECK	Memi[$1+2]
define	GXN_P_INDEX		Memi[$1+3]
define	GXN_P_EXTNAME		Memi[$1+4]
define	GXN_P_EXTVER		Memi[$1+5]
define	GXN_P_IKPARAMS		Memi[$1+6]
define	GXN_P_OMIT		Memi[$1+7]
define	GXN_P_REPLACE		Memi[$1+8]
define	GXN_FDOUT		Memi[$1+9]
define	GXN_COUNT		Memi[$1+10]	# index conversion?
define	GXN_FAILCOUNT		Memi[$1+11]
define	GXN_NARGS		Memi[$1+12]
define	GXN_P_GL		Memi[$1+13]	# GEMLOG GL structure
define	GXN_P_OP		Memi[$1+14]	# GEMLOG OP structure

# direct access to contents
define	GXN_IMG_CHECK	Memi[GXN_P_IMG_CHECK[$1]]
define	GXN_EXT_CHECK	Memi[GXN_P_EXT_CHECK[$1]]
define	GXN_AUTO_CHECK	Memi[GXN_P_AUTO_CHECK[$1]]
define	GXN_INDEX	Memc[GXN_P_INDEX[$1]]
define	GXN_EXTNAME	Memc[GXN_P_EXTNAME[$1]]
define	GXN_EXTVER	Memc[GXN_P_EXTVER[$1]]
define	GXN_IKPARAMS	Memc[GXN_P_IKPARAMS[$1]]
define	GXN_OMIT	Memi[GXN_P_OMIT[$1]]
define	GXN_REPLACE	Memc[GXN_P_REPLACE[$1]]


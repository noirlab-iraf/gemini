# Copyright(c) 2004-2005 Association of Universities for Research in Astronomy, Inc.

# ToDo List (aka "caveats" until they get done!) NOTE: "TODO::" is for easy greping 
# TODO::() mef filename validation (e.g. no extension numbers.. allowed section parsing)
# TODO::() determine final lock policy, mefs are not kept open after memap(..)
# TODO::() validate PHU and MEF in general on open (which conditions are MEFIO memap errors?
# TODO::() [DONE] store EXTVER in MEXT structure?
# TODO::() NAMING: should EXTCOUNTS be REFCOUNTS?

#extension types
define ET_ANY	-1
define ET_OTHER	0
define ET_IMAGE 1
define ET_TABLE 2

#extension optypes
define OPT_NULL   0
define OPT_NORMAL 1
define OPT_SCI    2
define OPT_VAR    3
define OPT_DQ     4    

# MEFIO structure types
define MST_MEXT		1
define MST_EXT		2

# structures
#NOTE: due to the way these structures have to be defined in SPP I
# 	have just been adding elements to the bottom, I will be
#   rearranging them when the library is stable

# MEXT STRUCTURE
# REVISIT:: TODO:: I got an access violation freeing this structure when I
# incremented 	by 1, which can noy be true, 
define LEN_MEXT		30 			# length of MEXTP "structure"
define ME_MST		Memi[$1]	# MEFIO structure type
define ME_NAMES		Memi[$1+2]	# pointer to array of ext names (char arrays)
define ME_TYPES		Memi[$1+4]	# pointer to array of types (ints, see ET_xxx enum)
define ME_NUMEXTS	Memi[$1+6]	# number of extension physically counted including PHU 
define ME_EPS 		Memi[$1+8]	# pointer to array of extpointers
define ME_NEXTEND	Memi[$1+10]	# NEXTEND value from PHU 
define ME_CUREXTI 	Memi[$1+12]     # Current extension for iterator funcs
define ME_PFILENAME	Memi[$1+14]	# pointer to filename
define ME_PSECSTR	Memi[$1+16]	# pointer to section string
define ME_FILENAME	Memc[ME_PFILENAME[$1]]	# character ref for filename pointer
define ME_SECSTR	Memc[ME_PSECSTR[$1]]	# character ref for filename pointer
define ME_EXTVERS	Memi[$1+18]     # pointer to array of EXTVERS for each extension
define ME_EXTCOUNTS	Memi[$1+20]     # pointer to array of extension access counts
define ME_REFCOUNTS	Memi[$1+20]     # pointer to array of extension access counts
define ME_NUMNAMEDEXTS	Memi[$1+21]     # number of named extensions
define ME_OPTYPE        Memi[$1+22]     # additional type related to operator type when used in expressions
define ME_PHU           Memi[$1+23]     # pointer for PHU


# EXT STRUCTURE (TABLE and IMAGE types supported)
define LEN_EXT		16
define EXT_MST		Memi[$1]	# MEFIO structure type
define EXT_MEXTP	Memi[$1+2]	# pointer to parent MEXTP 
define EXT_INDEX	Memi[$1+4]	# absolute extension position (also position in MEXT arrays)
define EXT_EXTTYPE	Memi[$1+6]	# ET_ type 
define EXT_EXTP		Memi[$1+8]	# pointer to mapped in image or table pointer
define EXT_OPENMODE	Memi[$1+9]	# the mode used to open the extension (READ_WRITE, READ_ONLY)

# I2DARY
define I2D_ALLVALID -1			# OLDNROWS value if all rows are valid
define LEN_I2D	4
define I2D_XW	Memi[$1]
define I2D_YW	Memi[$1+1]
define I2D_NCOLS Memi[$1]
define I2D_NROWS Memi[$1+1]
define I2D_REALNROWS Memi[$1+2]
define I2D_BUF	Memi[$1+3]

define I2DROWVALID    0
define I2DROWINVAL    1

#FM related: frame correlation
define FM_SCALEROPERAND -1

#errors
define MEERR_NOFILE	201
define MEERR_BADSTRUCT  202
define MEERR_OUTOFRANGE 203
define MEERR_INVAL	204
define MEERR_GDRELATE   205
define MEERR_BADTYPE    206
define MEERR_NONEXTEND  207
define MEERR_BADOPERAND 208
define MEERR_NOSECTIONS 209
define MEERR_MEMAPFAIL  210
define MEERR_NOTEXIST   211
define MEERR_BUFFERERR  212
define MEERR_EXTNOTAVAILABLE 213
define MEERR_VALNOTFOUND 214
define MEERR_COPYFAILED	215

#warnings
define MEWRN_UNKNOWN_EXTTYPE 250

#return values for memap
define ME_FAILURE       0
define ME_SUCCESS 	1

#special OPENMODEs

define ANY_MODE 100

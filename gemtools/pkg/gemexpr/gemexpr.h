# Copyright(c) 2004-2006 Association of Universities for Research in Astronomy, Inc.

# IMEXPR.X -- Image expression evaluator.

define	MAX_OPERANDS	26
define	MAX_ALIASES	10
define	DEF_LENINDEX	97
define	DEF_LENSTAB	1024
define	DEF_LENSBUF	8192
define	DEF_LINELEN	8192

# Input image operands.
define	LEN_IMOPERAND	22
define	IO_OPNAME	Memi[$1]		# symbolic operand name
define	IO_TYPE		Memi[$1+1]		# operand type
define	IO_MEF		Memi[$1+2]		# mef pointer from mefio
define  IO_EXT		Memi[$1+3]		# extension pointer
define	IO_IM		Memi[$1+4]		# image pointer if image
define	IO_V		Memi[$1+6+($2)-1]	# image i/o pointer
define	IO_DATA		Memi[$1+14]		# current image line
			# align
define	IO_OP		($1+16)			# pointer to evvexpr operand

# Image operand types (IO_TYPE). 
define	IMAGE		1			# image (vector) operand
define	NUMERIC		2			# numeric constant
define	PARAMETER	3			# image parameter reference

# Main imexpr descriptor.
define	LEN_IMEXPR	(24+LEN_IMOPERAND*MAX_OPERANDS)
define	IE_ST		Memi[$1]		# symbol table
define	IE_IM		Memi[$1+1]		# output image
define	IE_NDIM		Memi[$1+2]		# dimension of output image
define	IE_AXLEN	Memi[$1+3+($2)-1]	# dimensions of output image
define	IE_INTYPE	Memi[$1+10]		# minimum input operand type
define	IE_OUTTYPE	Memi[$1+11]		# datatype of output image
define	IE_BWIDTH	Memi[$1+12]		# npixels boundary extension
define	IE_BTYPE	Memi[$1+13]		# type of boundary extension
define	IE_BPIXVAL	Memr[$1+14]		# boundary pixel value
define	IE_V		Memi[$1+15+($2)-1]	# position in output image
define	IE_NOPERANDS	Memi[$1+22]		# number of input operands
			# align
define	IE_IMOP		($1+24+(($2)-1)*LEN_IMOPERAND)	# image operand array

# Expression database symbol.
define	LEN_SYM		2
define	SYM_TEXT	Memi[$1]
define	SYM_NARGS	Memi[$1+1]

# Argument list symbol
define	LEN_ARGSYM	1
define	ARGNO		Memi[$1]

# gemexpr errors

define GEERR_GENERR           300
define GEERR_INVALOPERANDS    301
define GEERR_OUTPUTEXISTS     302
define GEERR_W_GEMEXPRFAILED  303
define GEERR_EXPRPROC         304
define GEERR_NEEDIMFORVAR     305
define GEERR_OPERANDRANGE     306
define GEERR_NOIMOP           307
define GEERR_MALFORMEDOP      308
define GEERR_SAMEINOUT        309
define GEERR_NODATA           310
define GEERR_BADTYPE          311
define GEERR_BADIMPARMREF     312
define GEERR_PARAMNOTFOUND    313
define GEERR_BADWCS           314
define GEERR_BADOUTTYPE       315
define GEERR_INCOMPATTYPE     316
define GEERR_MALFORMEDOP      317
define GEERR_BADREFIM         318
define GEERR_BADCOPY          319
define GEERR_NOSCIENCE        320
define GEERR_PARMCOLLIDE      321
define GEERR_BADOPERATOR      322
define GEERR_EXTENTMISTMATCH  323
define GEERR_INVALIDEXPR      324
define GEERR_BADPARM          325

define GEWARN_NOFRAMES        400




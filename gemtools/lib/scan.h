# Copyright(c) 2004-2006 Association of Universities for Research in Astronomy, Inc.
#
# A structure for storing string information as it is scanned

define	LEN_SSCAN 4

define	GT_SC_IP	Memi[$1]
define	GT_SC_NTOKENS	Memi[$1+1]
define	GT_SC_STOPSCAN	Memb[$1+2]
define	GT_P_SC_SCANBUF	Memi[$1+3]

define	GT_SC_SCANBUF	Memc[GT_P_SC_SCANBUF($1)+$2-1]	# unit indexing


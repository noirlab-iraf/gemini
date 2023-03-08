# Copyright(c) 2004-2006 Association of Universities for Research in Astronomy, Inc.

# Author: Kathleen Labrie  3-May-2004

# GLOG.H -- GEMLOG internal definitions.

###################\
# Sizes and Lengths \
#####################\
#
define G_SZ_LABEL	   4
define G_SZ_SSTR	   15	#Size of a short string
define G_MAX_PARAM	   100
define G_SZ_PARAMSTR	   G_MAX_PARAM*SZ_LINE
define G_MAX_LINES	   10	#Max number of lines to be written at once
define G_MAX_BLK	   100	#Maximum number of block to retrieve
define G_MAX_ENTRY_PER_BLK 1000	#Max number of lines to retrieve in a block	

#######################\
# Structure definitions \
#########################\
#
#### GL Structure ####
#
define LEN_GL		11		#Size of the GL structure

define GL_FD		Meml[$1]	#File descriptor for the logfile
define GL_REQSTAT	Memi[($1)+2]	#Are status level entries requested?
define GL_REQSCI	Memi[($1)+3]	#Are science level entries requested?
define GL_REQENG	Memi[($1)+4]	#Are engineering level entries requested?
define GL_REQVIS	Memi[($1)+5]	#Are visual level entries requested?
define GL_REQTSK	Memi[($1)+6]	#Are task level entries requested?
define GL_VERBOSE	Memi[($1)+7]	#Verbose?
define GL_LOG_P		Memi[($1)+8]	#Logfile pointer
define GL_CPKG_P	Memi[($1)+9]	#Current package name pointer
define GL_CTASK_P	Memi[($1)+10]	#Current task name pointer
define GL_LOGFILE	Memc[GL_LOG_P($1)]	#Name of the logfile
define GL_CURPACK	Memc[GL_CPKG_P($1)]	#Name of current package
define GL_CURTASK	Memc[GL_CTASK_P($1)]	#Name of current task
#Next macro: define GL_???  Mem?[($1)+11]

#
#### SL Structure ####
# (Block Selection Structure)
#
define LEN_SL		8+G_MAX_BLK	#Size of the SEL structure

define SL_TSK_P		Memi[$1]	#Task name pointer
define SL_LTIME		Memi[($1)+1]    #Lower limit on time range
define SL_UTIME		Memi[($1)+2]	#Upper limit on time range
define SL_CHILD		Memi[($1)+3]	#Retrieve the logs for child processes?
define SL_NCHILD	Memi[($1)+4]    #Maximum number of subprocesses
define SL_NBLK		Memi[($1)+5]	#Number of blocks to retrieve
define SL_BLKS		Memi[($1)+5+($2)]  #Array of the block position to retrieve
define SL_BPOS		Memi[($1)+6+G_MAX_BLK] #Current block position
define SL_TSKNAME	Memc[SL_TSK_P($1)]	#Name of the task to retrieve
#Next macro: define SL_???  Mem?[($1)+7+G_MAX_BLK]

#
#### OP Structure ####
# (API options)
#
define LEN_OP		9		#Size of the OP structure

define OP_FL_APPEND	Memi[$1]	#Append to logfile?
define OP_FORCE_APPEND	Memi[($1)+1]	#File must already exists?
define OP_DEFLOG	Memi[($1)+2]	#Using default logfile name
define OP_VERBOSE	Memi[($1)+3]	#Verbose?
define OP_STATUS	Memi[($1)+4]	#Exit status (0=good)
define OP_VISTYPE	Memi[($1)+5]	#Type of visual enhancement
define OP_ERRNO		Memi[($1)+6]	#Gemini error code
define OP_FORK		Memi[($1)+7]	#Fork to or back from child process?
define OP_CHILD_P	Memi[($1)+8]	#Child name pointer
define OP_CHILD		Memc[OP_CHILD_P($1)]	#Name of child process
#Next macro: define OP_??? Mem?[($1)+9]

#############################\
# Log entry level definitions \
###############################\
#
define LEN_LEVEL_STR	4		#Length of level tag string

define IGNORE_LEVEL	0		#Ignore level (eg. str already has tag)
define NO_LEVEL		1		#Leave level section empty
define STAT_LEVEL	2		#Status level
define SCI_LEVEL	4		#Science level
define ENG_LEVEL	8		#Engineering level
define VIS_LEVEL	16		#Visual level (to improve readability)
define TSK_LEVEL	32		#Task level

##################################\
# Log type definitions (GLOGPRINT) \
####################################\
#
define G_STR_LOG	1		#Normal string
define G_ERR_LOG	2		#Error message
define G_WARN_LOG	4		#Warning message
define G_FORK_LOG	8		#Fork to or back from child process
define G_FILE_LOG	16		#Entries taken from a file
define G_VIS_LOG	32		#Visual aid

##############################\
# Error codes (from gemerrmsg) \
################################\
#
define	G_INTERNAL_ERROR	99
define  G_OPEN_FILE_ERR		100
define	G_FILE_NOT_ACCESSIBLE	101
define  G_FILE_EXIST		102

define	G_USING_DEFAULT		120
define	G_INPUT_ERROR		121
define	G_OP_UNRECOGNIZED	122
define	G_WRONG_IMG_FORMAT	123

define	G_KEY_NOT_FOUND		131
define	G_HDR_ERR		132

###################\
# Exit status codes \
#####################\
#
define	G_SUCCESS	YES
define	G_FAILURE	NO

#####################\
# Log tag definitions \
#######################\
#
define BEGIN_TAG	1		#Tag beginning of logs for current task
define END_TAG		2		#Tag end of logs for current task

################################################\
# Type of visual (readability improvement) lines \
##################################################\
#
define G_EMPTY		1	#An empty line
define G_LONG_DASH	2	#A long (80 ch) line of dashes
define G_SHRT_DASH	4	#A short (20 ch) line of dashes

###############\
# Miscellaneous \
#################\
#
define G_INDEF		-1

define G_FORWARD	0	#Fork to child
define G_BACKWARD	1	#Come back to parent

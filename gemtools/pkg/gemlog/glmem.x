# Copyright(c) 2004-2006 Association of Universities for Research in Astronomy, Inc.

# Author: Kathleen Labrie  1-Jul-2004

include "glog.h"

.help
.nf
This file contains procedures to allocate and free memory for the structures
defined in 'glog.h'.

   glalloc - Allocate memory for a GL structure
    glfree - Free memory attached to a GL structure

   opalloc - Allocate memory for an OP structure
    opfree - Free memory attached to an OP structure

   slalloc - Allocate memory for an SL structure
    slfree - Free memory attached to an SL structure

.fi
.endhelp

#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
#--------------------------------------------------------------------------

# GLALLOC -- Allocate memory for a GL structure (defined in glog.h)
#	call glalloc ( gl )
#
#	gl		: GL structure  [output, (GL)]

procedure glalloc (gl)

pointer	gl			#O  GL structure pointer

begin
	call malloc (gl, LEN_GL, TY_STRUCT)

	call malloc (GL_LOG_P(gl), SZ_FNAME, TY_CHAR)
	call malloc (GL_CPKG_P(gl), SZ_FNAME, TY_CHAR)
	call malloc (GL_CTASK_P(gl), SZ_FNAME, TY_CHAR)

	return
end

#--------------------------------------------------------------------------

# GLFREE -- Free memory attached to a GL structure (defined in glog.h)
#	call glfree ( gl )
#
#	gl		: GL structure  [input/output, (GL)]

procedure glfree (gl)

pointer gl			#IO  GL structure pointer

begin
	call mfree (GL_CTASK_P(gl), TY_CHAR)
	call mfree (GL_CPKG_P(gl), TY_CHAR)
	call mfree (GL_LOG_P(gl), TY_CHAR)

	call mfree (gl, TY_STRUCT)
	gl = NULL

	return
end

#--------------------------------------------------------------------------

# OPALLOC -- Allocate memory for an OP structure (defined in glog.h)
#	call opalloc ( op )
#
#	op		: OP structure  [output, (OP)]

procedure opalloc (op)

pointer	op		#O  OP structure pointer

begin
	call malloc (op, LEN_OP, TY_STRUCT)
	call malloc (OP_CHILD_P(op), SZ_FNAME, TY_CHAR)

	return
end

#--------------------------------------------------------------------------

# OPFREE -- Free memory attached to an OP structure (defined in glog.h)
#	call opfree ( op )
#
#	op		: OP structure  [input/output, (OP)]

procedure opfree (op)

pointer op			#IO  OP structure pointer

begin
	call mfree (OP_CHILD_P(op), TY_CHAR)
	call mfree (op, TY_STRUCT)
	op = NULL

	return
end

#--------------------------------------------------------------------------

# SLALLOC -- Allocate memory for an SL structure (defined in glog.h)
#	call slalloc ( sl )
#
#	sl		: SL structure  [output, (SL)]

procedure slalloc (sl)

pointer	sl			#O  SL structure pointer

begin
	call malloc (sl, LEN_SL, TY_STRUCT)
	call malloc (SL_TSK_P(sl), SZ_FNAME, TY_CHAR)

	return
end

#--------------------------------------------------------------------------

# SLFREE -- Free memory attached to an SL structure (defined in glog.h)
#	call slfree ( sl )
#
#	sl		: SL structure  [input/output, (SL)]

procedure slfree (sl)

pointer sl			#IO  SL structure pointer

begin
	call mfree (SL_TSK_P(sl), TY_CHAR)
	call mfree (sl, TY_STRUCT)
	sl = NULL

	return
end

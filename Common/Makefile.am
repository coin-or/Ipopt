# Copyright (C) 2004, International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Common Public License.

## $Id$

# Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

AUTOMAKE_OPTIONS = foreign

noinst_LIBRARIES = libcommon.a

libcommon_a_SOURCES = \
	IpCachedResults.hpp \
	IpCompoundMatrix.cpp IpCompoundMatrix.hpp \
	IpCompoundSymMatrix.cpp IpCompoundSymMatrix.hpp \
	IpCompoundVector.cpp IpCompoundVector.hpp \
	IpDebug.cpp IpDebug.hpp \
	IpDiagMatrix.cpp IpDiagMatrix.hpp \
	IpException.hpp \
	IpIdentityMatrix.cpp IpIdentityMatrix.hpp \
	IpIpoptType.hpp IpIpoptType.cpp \
	IpJournalist.cpp IpJournalist.hpp \
	IpMatrix.cpp \
	IpObserver.cpp IpObserver.hpp \
	IpOptionsList.cpp IpOptionsList.hpp \
	IpRegOptions.cpp IpRegOptions.hpp \
	IpReferenced.cpp IpReferenced.hpp \
	IpSmartPtr.hpp \
	IpSumMatrix.cpp IpSumMatrix.hpp \
	IpSumSymMatrix.cpp IpSumSymMatrix.hpp \
	IpSymMatrix.cpp \
	IpTaggedObject.cpp IpTaggedObject.hpp \
	IpTypes.hpp \
	IpUtils.hpp \
	IpVector.cpp \
	IpZeroMatrix.cpp IpZeroMatrix.hpp

CXXFLAGS += -I$(srcdir)/../Interfaces

# Astyle stuff

ASTYLE_SOURCES = \
	IpCachedResults.hppbak \
	IpCompoundMatrix.cppbak IpCompoundMatrix.hppbak \
	IpCompoundSymMatrix.cppbak IpCompoundSymMatrix.hppbak \
	IpCompoundVector.cppbak IpCompoundVector.hppbak \
	IpDebug.cppbak IpDebug.hppbak \
	IpDiagMatrix.cppbak IpDiagMatrix.hppbak \
	IpException.hppbak \
	IpIdentityMatrix.cppbak IpIdentityMatrix.hppbak \
	IpIpoptType.hppbak IpIpoptType.cppbak \
	IpJournalist.cppbak IpJournalist.hppbak \
	IpMatrix.cppbak \
	IpObserver.cppbak IpObserver.hppbak \
	IpOptionsList.hppbak \
	IpRegOptions.cppbak IpRegOptions.hppbak \
	IpReferenced.cppbak IpReferenced.hppbak \
	IpSmartPtr.hppbak \
	IpSumMatrix.cppbak IpSumMatrix.hppbak \
	IpSumSymMatrix.cppbak IpSumSymMatrix.hppbak \
	IpSymMatrix.cppbak \
	IpTaggedObject.cppbak IpTaggedObject.hppbak \
	IpTypes.hppbak \
	IpUtils.hppbak \
	IpVector.cppbak \
	IpZeroMatrix.cppbak IpZeroMatrix.hppbak

ASTYLE = @ASTYLE@
ASTYLEFLAGS = @ASTYLEFLAGS@

CLEANFILES = $(ASTYLE_SOURCES)

SUFFIXES = .cppbak .hppbak

astyle: $(ASTYLE_SOURCES)

.hpp.hppbak:
	mv $< $@
	$(ASTYLE) $(ASTYLEFLAGS) < $@ > $<
	touch $@

.cpp.cppbak:
	mv $< $@
	$(ASTYLE) $(ASTYLEFLAGS) < $@ > $<
	touch $@

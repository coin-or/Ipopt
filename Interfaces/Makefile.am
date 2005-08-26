# Copyright (C) 2004, International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Common Public License.

## $Id$

# Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

AUTOMAKE_OPTIONS = foreign

includeipoptdir = $(includedir)/ipopt
includeipopt_HEADERS = \
	IpAlgTypes.hpp \
	IpIpoptApplication.hpp \
	IpNLP.hpp \
	IpReturnCodes.h \
	IpReturnCodes_inc.h \
	IpReturnCodes.inc \
	IpSolveStatistics.hpp \
	IpStdCInterface.h \
	IpTNLP.hpp

noinst_LIBRARIES = libinterfaces.a

libinterfaces_a_SOURCES = \
	IpAlgTypes.hpp \
	IpInterfacesRegOp.cpp IpInterfacesRegOp.hpp \
	IpIpoptApplication.cpp IpIpoptApplication.hpp \
	IpNLP.hpp \
	IpReturnCodes.h \
	IpSolveStatistics.cpp IpSolveStatistics.hpp \
	IpStdCInterface.cpp IpStdCInterface.h \
	IpStdFInterface.c \
	IpStdInterfaceTNLP.cpp IpStdInterfaceTNLP.hpp \
	IpTNLP.hpp \
	IpTNLPAdapter.cpp IpTNLPAdapter.hpp

CXXFLAGS += -I$(srcdir)/../Common \
	-I$(srcdir)/../LinAlg \
	-I$(srcdir)/../Algorithm/LinearSolvers \
	-I$(srcdir)/../LinAlg/TMatrices \
	-I$(srcdir)/../Algorithm

CFLAGS += -I$(srcdir)/../Common

# Astyle stuff

ASTYLE_FILES = \
	IpAlgTypes.hppbak \
	IpInterfacesRegOp.cppbak IpInterfacesRegOp.hppbak \
	IpIpoptApplication.cppbak IpIpoptApplication.hppbak \
	IpNLP.hppbak \
	IpSolveStatistics.cppbak IpSolveStatistics.hppbak \
	IpStdCInterface.cppbak \
	IpStdInterfaceTNLP.cppbak IpStdInterfaceTNLP.hppbak \
	IpTNLP.hppbak \
	IpTNLPAdapter.cppbak IpTNLPAdapter.hppbak

ASTYLE = @ASTYLE@
ASTYLEFLAGS = @ASTYLEFLAGS@

CLEANFILES = $(ASTYLE_FILES)

SUFFIXES = .cppbak .hppbak

astyle: $(ASTYLE_FILES)

.hpp.hppbak:
	mv $< $@
	$(ASTYLE) $(ASTYLEFLAGS) < $@ > $<
	touch $@

.cpp.cppbak:
	mv $< $@
	$(ASTYLE) $(ASTYLEFLAGS) < $@ > $<
	touch $@

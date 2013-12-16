# Copyright (C) 2004, 2008 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.

## $Id$

# Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

AUTOMAKE_OPTIONS = foreign

includeipoptdir = $(includedir)/coin
includeipopt_HEADERS = \
	IpAlgTypes.hpp \
	IpIpoptApplication.hpp \
	IpNLP.hpp \
	IpReturnCodes.h \
	IpReturnCodes.hpp \
	IpReturnCodes_inc.h \
	IpReturnCodes.inc \
	IpSolveStatistics.hpp \
	IpStdCInterface.h \
	IpTNLP.hpp \
	IpTNLPAdapter.hpp \
	IpTNLPReducer.hpp

lib_LTLIBRARIES = libipopt.la

libipopt_la_SOURCES = \
	IpAlgTypes.hpp \
	IpInterfacesRegOp.cpp IpInterfacesRegOp.hpp \
	IpIpoptApplication.cpp IpIpoptApplication.hpp \
	IpNLP.hpp \
	IpReturnCodes.h \
	IpReturnCodes.hpp \
	IpSolveStatistics.cpp IpSolveStatistics.hpp \
	IpStdCInterface.cpp IpStdCInterface.h \
	IpStdFInterface.c \
	IpStdInterfaceTNLP.cpp IpStdInterfaceTNLP.hpp \
	IpTNLP.hpp \
	IpTNLPAdapter.cpp IpTNLPAdapter.hpp \
	IpTNLPReducer.cpp IpTNLPReducer.hpp

libipopt_la_LIBADD = $(IPALLLIBS)
if DEPENDENCY_LINKING
libipopt_la_LIBADD += $(IPOPTLIB_LIBS)
endif

libipopt_la_DEPENDENCIES = $(IPALLLIBS)

libipopt_la_LDFLAGS = $(LT_LDFLAGS)

CXXFLAGS += -I`$(CYGPATH_W) $(srcdir)/../Common` \
	-I`$(CYGPATH_W) $(srcdir)/../LinAlg` \
	-I`$(CYGPATH_W) $(srcdir)/../Algorithm/LinearSolvers` \
	-I`$(CYGPATH_W) $(srcdir)/../LinAlg/TMatrices` \
	-I`$(CYGPATH_W) $(srcdir)/../Algorithm` \
	-I`$(CYGPATH_W) $(srcdir)/../Algorithm/Inexact` \
	-I`$(CYGPATH_W) $(srcdir)/../contrib/CGPenalty` \
	-I`$(CYGPATH_W) $(srcdir)/../contrib/LinearSolverLoader` \
	$(HSL_CFLAGS)

CFLAGS += -I`$(CYGPATH_W) $(srcdir)/../Common` \
	$(HSL_CFLAGS)

# This line is necessary to allow VPATH compilation
DEFAULT_INCLUDES = -I. -I`$(CYGPATH_W) $(srcdir)` -I$(top_builddir)/src/Common

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
	IpTNLPAdapter.cppbak IpTNLPAdapter.hppbak \
	IpTNLPReducer.cppbak IpTNLPReducer.hppbak

ASTYLE = @ASTYLE@
ASTYLEFLAGS = @ASTYLEFLAGS@

CLEANFILES =

DISTCLEANFILES = $(ASTYLE_FILES)

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

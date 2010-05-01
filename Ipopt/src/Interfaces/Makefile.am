# Copyright (C) 2004, 2008 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Common Public License.

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
	IpTNLPReducer.hpp

coinlibdir = $(libdir)/coin
coinlib_LTLIBRARIES = libipopt.la

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

libipopt_la_DEPENDENCIES = $(IPALLLIBS)

libipopt_la_LDFLAGS = $(LT_LDFLAGS)

addlibsdir = $(datadir)/coin/doc/Ipopt
addlibs_DATA = ipopt_addlibs_f.txt

if COIN_HAS_PKGCONFIG
addlibs_DATA += ipopt_addlibs_cpp.txt

ipopt_addlibs_cpp.txt : ../../ipopt.pc
	PKG_CONFIG_PATH=../..:@COIN_PKG_CONFIG_PATH@ $(PKG_CONFIG) --libs ipopt > ipopt_addlibs_cpp.txt
endif

CXXFLAGS += -I`$(CYGPATH_W) $(srcdir)/../Common` \
	-I`$(CYGPATH_W) $(srcdir)/../LinAlg` \
	-I`$(CYGPATH_W) $(srcdir)/../Algorithm/LinearSolvers` \
	-I`$(CYGPATH_W) $(srcdir)/../LinAlg/TMatrices` \
	-I`$(CYGPATH_W) $(srcdir)/../Algorithm` \
	-I`$(CYGPATH_W) $(srcdir)/../Algorithm/Inexact` \
	-I`$(CYGPATH_W) $(srcdir)/../contrib/CGPenalty`

CFLAGS += -I`$(CYGPATH_W) $(srcdir)/../Common`

if COIN_HAS_HSL
CXXFLAGS += $(HSL_CFLAGS)
CFLAGS   += $(HSL_CFLAGS)
endif

# This line is necessary to allow VPATH compilation with MS compilers
# on Cygwin
DEFAULT_INCLUDES = -I. -I`$(CYGPATH_W) $(srcdir)` -I$(top_builddir)/inc

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

CLEANFILES = ipopt_addlibs_cpp.txt

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

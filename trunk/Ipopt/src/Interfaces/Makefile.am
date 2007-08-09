# Copyright (C) 2004, 2006 International Business Machines and others.
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
	IpReturnCodes.hpp \
	IpReturnCodes_inc.h \
	IpReturnCodes.inc \
	IpSolveStatistics.hpp \
	IpStdCInterface.h \
	IpTNLP.hpp

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
	IpTNLPAdapter.cpp IpTNLPAdapter.hpp

libipopt_la_LIBADD = $(IPALLLIBS)

libipopt_la_DEPENDENCIES = $(IPALLLIBS)

libipopt_la_LDFLAGS = $(LT_LDFLAGS)

ADDLIBS_FILES = \
	ipopt_addlibs_cpp.txt \
	ipopt_addlibs_f.txt

CXXFLAGS += -I`$(CYGPATH_W) $(srcdir)/../Common` \
	-I`$(CYGPATH_W) $(srcdir)/../LinAlg` \
	-I`$(CYGPATH_W) $(srcdir)/../Algorithm/LinearSolvers` \
	-I`$(CYGPATH_W) $(srcdir)/../LinAlg/TMatrices` \
	-I`$(CYGPATH_W) $(srcdir)/../Algorithm` \
	-I`$(CYGPATH_W) $(srcdir)/../contrib/CGPenalty`

CFLAGS += -I`$(CYGPATH_W) $(srcdir)/../Common`

# This line is necessary to allow VPATH compilation with MS compilers
# on Cygwin
DEFAULT_INCLUDES = -I. -I`$(CYGPATH_W) $(srcdir)` -I$(top_builddir)/inc

# Make sure the addlibs files are installed in lib
install-exec-local:
	for i in $(ADDLIBS_FILES); do \
	  $(install_sh_DATA) $$i $(DESTDIR)$(libdir)/$$i; \
	done

uninstall-local:
	for i in $(ADDLIBS_FILES); do \
	  rm -f $(DESTDIR)$(libdir)/$$i; \
	done

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

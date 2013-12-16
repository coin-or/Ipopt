# Copyright (C) 2004, 2006 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.

## $Id$

# Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

AUTOMAKE_OPTIONS = foreign

lib_LTLIBRARIES = libipoptamplinterface.la
bin_PROGRAMS = ipopt

libipoptamplinterface_la_SOURCES = \
	AmplTNLP.cpp AmplTNLP.hpp

if DEPENDENCY_LINKING
libipoptamplinterface_la_LIBADD = ../../Interfaces/libipopt.la \
	$(IPOPTAMPLINTERFACELIB_LIBS)
endif

libipoptamplinterface_la_LDFLAGS = $(LT_LDFLAGS)

ipopt_SOURCES = \
	ampl_ipopt.cpp

ipopt_LDADD = libipoptamplinterface.la ../../Interfaces/libipopt.la \
	$(IPOPTAMPLINTERFACELIB_LIBS) $(IPOPTLIB_LIBS)

ipopt_DEPENDENCIES = libipoptamplinterface.la ../../Interfaces/libipopt.la \
	$(IPOPTAMPLINTERFACELIB_DEPENDENCIES) $(IPOPTLIB_DEPENDENCIES)

AM_CPPFLAGS = -I`$(CYGPATH_W) $(srcdir)/../../Common` \
	-I`$(CYGPATH_W) $(srcdir)/../../LinAlg` \
	-I`$(CYGPATH_W) $(srcdir)/../../LinAlg/TMatrices` \
	-I`$(CYGPATH_W) $(srcdir)/../../Algorithm` \
	-I`$(CYGPATH_W) $(srcdir)/../../Interfaces`

AM_CPPFLAGS += $(ASL_CFLAGS)

includeipoptdir = $(includedir)/coin
includeipopt_HEADERS = AmplTNLP.hpp

pkgconfiglibdir = $(libdir)/pkgconfig
pkgconfiglib_DATA = ../../../ipoptamplinterface.pc

# This line is necessary to allow VPATH compilation
DEFAULT_INCLUDES = -I. -I`$(CYGPATH_W) $(srcdir)` -I$(top_builddir)/src/Common

if COIN_CXX_IS_CL
install-exec-hook:
	test -s ipopt.exe.manifest && \
	mt -manifest ipopt.exe.manifest -outputresource:ipopt.exe;\
	cp ipopt.exe $(bindir)/ipopt.exe
endif

# Astyle stuff

ASTYLE_FILES = \
	AmplTNLP.cppbak AmplTNLP.hppbak \
	ampl_ipopt.cppbak

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

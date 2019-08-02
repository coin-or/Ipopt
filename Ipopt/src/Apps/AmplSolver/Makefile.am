# Copyright (C) 2004, 2006 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.

# Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

AUTOMAKE_OPTIONS = foreign

lib_LTLIBRARIES = libipoptamplinterface.la
bin_PROGRAMS = ipopt

libipoptamplinterface_la_SOURCES = \
	AmplTNLP.cpp AmplTNLP.hpp

libipoptamplinterface_la_LIBADD = ../../Interfaces/libipopt.la \
	$(IPOPTAMPLINTERFACELIB_LIBS)

libipoptamplinterface_la_LDFLAGS = $(LT_LDFLAGS)

ipopt_SOURCES = \
	ampl_ipopt.cpp

ipopt_LDADD = libipoptamplinterface.la ../../Interfaces/libipopt.la \
	$(IPOPTAMPLINTERFACELIB_LIBS) $(IPOPTLIB_LIBS)

AM_CPPFLAGS = -DIPOPT_BUILD \
	-I$(srcdir)/../../Common \
	-I$(srcdir)/../../LinAlg \
	-I$(srcdir)/../../LinAlg/TMatrices \
	-I$(srcdir)/../../Algorithm \
	-I$(srcdir)/../../Interfaces \
	$(ASL_CFLAGS)

includeipoptdir = $(includedir)/coin
includeipopt_HEADERS = AmplTNLP.hpp

pkgconfiglibdir = $(libdir)/pkgconfig
pkgconfiglib_DATA = ../../../ipoptamplinterface.pc

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
	$(ASTYLE) $(ASTYLEFLAGS) < $@ | sed -e 's/  *$$//' > $<
	touch $@

.cpp.cppbak:
	mv $< $@
	$(ASTYLE) $(ASTYLEFLAGS) < $@ | sed -e 's/  *$$//' > $<
	touch $@

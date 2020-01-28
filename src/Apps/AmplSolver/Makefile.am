# Copyright (C) 2004, 2006 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.

# Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

AUTOMAKE_OPTIONS = foreign

lib_LTLIBRARIES = libipoptamplinterface.la
bin_PROGRAMS = ipopt

libipoptamplinterface_la_SOURCES = \
	AmplTNLP.cpp

libipoptamplinterface_la_LIBADD = ../../Interfaces/libipopt.la \
	$(IPOPTAMPLINTERFACELIB_LFLAGS)

ipopt_SOURCES = \
	ampl_ipopt.cpp

ipopt_LDADD = libipoptamplinterface.la ../../Interfaces/libipopt.la

LDFLAGS += $(LT_LDFLAGS)

AM_CPPFLAGS = -DIPOPTAMPL_BUILD \
	-I$(srcdir)/../../Common \
	-I$(srcdir)/../../LinAlg \
	-I$(srcdir)/../../LinAlg/TMatrices \
	-I$(srcdir)/../../Algorithm \
	-I$(srcdir)/../../Interfaces \
	$(IPOPTAMPLINTERFACELIB_CFLAGS)

includeipoptdir = $(includedir)/coin-or
includeipopt_HEADERS = AmplTNLP.hpp

pkgconfiglibdir = $(libdir)/pkgconfig
pkgconfiglib_DATA = ../../../ipoptamplinterface.pc

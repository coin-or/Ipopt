# Copyright (C) 2004, 2006 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.

# Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

AUTOMAKE_OPTIONS = foreign

includeipoptdir = $(includedir)/coin
includeipopt_HEADERS = \
	IpGenTMatrix.hpp \
	IpSymTMatrix.hpp \
	IpTripletHelper.hpp
noinst_LTLIBRARIES = libtmatrices.la

libtmatrices_la_SOURCES = \
	IpGenTMatrix.cpp \
	IpSymTMatrix.cpp \
	IpTripletHelper.cpp

libtmatrices_la_LDFLAGS = $(LT_LDFLAGS)

AM_CPPFLAGS = -DIPOPT_BUILD \
	-I$(srcdir)/../../Common \
	-I$(srcdir)/..

# Astyle stuff

ASTYLE_FILES = \
	IpGenTMatrix.cppbak IpGenTMatrix.hppbak \
	IpSymTMatrix.cppbak IpSymTMatrix.hppbak \
	IpTripletHelper.cppbak IpTripletHelper.hppbak

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

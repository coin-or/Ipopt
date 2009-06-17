# Copyright (C) 2009 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Common Public License.

## $Id$

# Authors:  Andreas Waechter               IBM    2009-05-28

AUTOMAKE_OPTIONS = foreign

#includeipoptdir = $(includedir)/coin
#includeipopt_HEADERS = 

noinst_LTLIBRARIES = liblinalgpar.la

liblinalgpar_la_SOURCES = \
	IpParExpansionMatrix.cpp IpParExpansionMatrix.hpp \
	IpParGenMatrix.cpp IpParGenMatrix.hpp \
	IpParSymMatrix.cpp IpSymGenMatrix.hpp \
	IpParVector.cpp IpParVector.hpp

liblinalgpar_la_LDFLAGS = $(LT_LDFLAGS)

AM_CPPFLAGS = -I`$(CYGPATH_W) $(srcdir)/..` \
	-I`$(CYGPATH_W) $(srcdir)/../TMatrices` \
	-I`$(CYGPATH_W) $(srcdir)/../../Common`

# This line is necessary to allow VPATH compilation with MS compilers
# on Cygwin
DEFAULT_INCLUDES = -I. -I`$(CYGPATH_W) $(srcdir)` -I$(top_builddir)/inc

# Astyle stuff

ASTYLE_FILES = \
	IpParExpansionMatrix.cppbak IpParExpansionMatrix.hppbak \
	IpParGenMatrix.cppbak IpParGenMatrix.hppbak \
	IpParSymMatrix.cppbak IpSymGenMatrix.hppbak \
	IpParVector.cppbak IpParVector.hppbak

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

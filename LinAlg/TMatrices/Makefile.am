# Copyright (C) 2004, 2006 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Common Public License.

## $Id$

# Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

AUTOMAKE_OPTIONS = foreign

noinst_LIBRARIES = libtmatrices.a

libtmatrices_a_SOURCES = \
	IpGenTMatrix.cpp IpGenTMatrix.hpp \
	IpSymTMatrix.cpp IpSymTMatrix.hpp \
	IpTripletHelper.cpp IpTripletHelper.hpp

CXXFLAGS += -I`$(CYGPATH_W) $(srcdir)/../../Common` \
	-I`$(CYGPATH_W) $(srcdir)/..`

# Astyle stuff

ASTYLE_FILES = \
	IpGenTMatrix.cppbak IpGenTMatrix.hppbak \
	IpSymTMatrix.cppbak IpSymTMatrix.hppbak \
	IpTripletHelper.cppbak IpTripletHelper.hppbak

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

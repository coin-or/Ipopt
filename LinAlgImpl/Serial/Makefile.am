# Copyright (C) 2004, International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Common Public License.

## $Id$

# Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

AUTOMAKE_OPTIONS = foreign

noinst_LIBRARIES = liblaserial.a

liblaserial_a_SOURCES = \
	IpBlas.cpp IpBlas.hpp \
	IpDenseVector.cpp IpDenseVector.hpp \
	IpExpansionMatrix.hpp IpExpansionMatrix.cpp \
	IpGenTMatrix.hpp IpGenTMatrix.cpp \
	IpSymTMatrix.hpp IpSymTMatrix.cpp \
	IpTNLPAdapter.hpp IpTNLPAdapter.cpp

CXXFLAGS += -I$(srcdir)/../../Common \
	-I$(srcdir)/../../Interfaces

# Astyle stuff

ASTYLE_SOURCES = IpBlas.cppbak IpBlas.hppbak \
	IpDenseVector.cppbak IpDenseVector.hppbak \
	IpExpansionMatrix.hppbak IpExpansionMatrix.cppbak \
	IpGenTMatrix.hppbak IpGenTMatrix.cppbak \
	IpSymTMatrix.hppbak IpSymTMatrix.cppbak \
	IpTNLPAdapter.hppbak IpTNLPAdapter.cppbak

ASTYLE = @ASTYLE@
ASTYLEFLAGS = @ASTYLEFLAGS@

CLEANFILES = $(ASTYLE_SOURCES)

SUFFIXES = .cppbak .hppbak

astyle: $(ASTYLE_SOURCES)

.hpp.hppbak:
	mv $< $@
	$(ASTYLE) $(ASTYLEFLAGS) < $@ > $<
	touch $@

.cpp.cppbak:
	mv $< $@
	$(ASTYLE) $(ASTYLEFLAGS) < $@ > $<
	touch $@

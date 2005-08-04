# Copyright (C) 2004, International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Common Public License.

## $Id$

# Authors:  Carl Laird, Andreas Waechter     IBM    2004-09-02

AUTOMAKE_OPTIONS = foreign

noinst_LIBRARIES = libstdinterface.a

libstdinterface_a_SOURCES = \
	IpInterfaceTypes.hpp \
	IpStdInterfaceTNLP.hpp IpStdInterfaceTNLP.cpp \
	IpStdCInterface.cpp IpStdFInterface.c \
	IpIpoptApplication.cpp IpIpoptApplication.hpp

CXXFLAGS += -I$(srcdir)/../../Common
CXXFLAGS += -I$(srcdir)/../../Interfaces
CXXFLAGS += -I$(srcdir)/../../LinAlgImpl/Serial
CXXFLAGS += -I$(srcdir)/../../Algorithm
CXXFLAGS += -I$(srcdir)/../../PDSystemImpl/PDFullSpace

CFLAGS += -I$(srcdir)/../../Interfaces

# Astyle stuff

ASTYLE_FILES = \
	IpInterfaceTypes.hppbak \
	IpStdInterfaceTNLP.hppbak IpStdInterfaceTNLP.cppbak \
	IpStdCInterface.cppbak \
	IpIpoptApplication.cppbak IpIpoptApplication.hppbak

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

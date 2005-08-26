# Copyright (C) 2004, International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Common Public License.

## $Id$

# Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

AUTOMAKE_OPTIONS = foreign

noinst_LIBRARIES = libamplinterface.a
bin_PROGRAMS = ipopt

libamplinterface_a_SOURCES = \
	AmplTNLP.cpp AmplTNLP.hpp

ipopt_SOURCES = \
	ampl_ipopt.cpp

IPLIBS = \
	../../Interfaces/libinterfaces.a \
	../../Algorithm/libipoptalg.a \
	../../Algorithm/LinearSolvers/liblinsolvers.a \
	../../LinAlg/TMatrices/libtmatrices.a \
	../../LinAlg/liblinalg.a \
	../../Common/libcommon.a \
	libamplinterface.a \
	$(ASLLIB)

ipopt_LDADD = $(IPLIBS)

CXXFLAGS += -I$(srcdir)/../../Common \
	-I$(srcdir)/../../LinAlg \
	-I$(srcdir)/../../LinAlg/TMatrices \
	-I$(srcdir)/../../Algorithm \
	-I$(srcdir)/../../Interfaces

CXXFLAGS += $(AMPL_CPPFLAGS)

# Astyle stuff

ASTYLE_FILES = \
	AmplTNLP.cppbak AmplTNLP.hppbak \
	ampl_ipopt.cppbak

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

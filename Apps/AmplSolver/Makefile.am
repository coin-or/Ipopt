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
	../../Apps/StdInterface/libstdinterface.a \
	../../Algorithm/libipoptalg.a \
	../../LinAlgImpl/Serial/liblaserial.a \
	../../PDSystemImpl/PDFullSpace/libpdfullspace.a \
	../../Algorithm/libipoptalg.a \
	../../Common/libcommon.a \
	libamplinterface.a \
	../../PDSystemImpl/PDFullSpace/libpdfullspace.a \
	../../Algorithm/libipoptalg.a \
	../../Common/libcommon.a \
	../../LinAlgImpl/Serial/liblaserial.a \
	$(ASLLIB)

ipopt_LDADD = $(IPLIBS)

CXXFLAGS += -I$(srcdir)/../../Common
CXXFLAGS += -I$(srcdir)/../../Interfaces
CXXFLAGS += -I$(srcdir)/../../LinAlgImpl/Serial
CXXFLAGS += -I$(srcdir)/../../Algorithm
CXXFLAGS += -I$(srcdir)/../../PDSystemImpl/PDFullSpace
CXXFLAGS += -I$(srcdir)/../../Apps/StdInterface
CXXFLAGS += -I$(srcdir)/../../Apps/AmplSolver \
	    -I$(srcdir)/../../OTHERS/ASL/ \
	    -I$(srcdir)/../../OTHERS/ASL/solvers
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

# Copyright (C) 2004, International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Common Public License.

## $Id$

# Authors:  Carl Laird, Andreas Waechter     IBM    2004-09-02

AUTOMAKE_OPTIONS = foreign

noinst_LIBRARIES = libcompositenlp.a
bin_PROGRAMS = composite_ampl_ipopt

libcompositenlp_a_SOURCES = \
	CompositeNLP.hpp CompositeNLP.cpp

composite_ampl_ipopt_SOURCES = \
	comp_ampl_ipopt.cpp

IPLIBS = ../../Algorithm/libipoptalg.a \
	../../LinAlgImpl/Serial/liblaserial.a \
	../../PDSystemImpl/PDFullSpace/libpdfullspace.a \
	../../Common/libcommon.a \
	../../Apps/AmplSolver/libamplinterface.a \
	$(ASLLIB)

composite_ampl_ipopt_LDADD = $(IPLIBS) \
	libcompositenlp.a


CXXFLAGS += -I$(srcdir)/../../Common
CXXFLAGS += -I$(srcdir)/../../Interfaces
CXXFLAGS += -I$(srcdir)/../../LinAlgImpl/Serial
CXXFLAGS += -I$(srcdir)/../../Algorithm
CXXFLAGS += -I$(srcdir)/../../PDSystemImpl/PDFullSpace
CXXFLAGS += -I$(srcdir)/../../Apps/AmplSolver \
	    -I$(srcdir)/../../OTHERS/ASL/ \
	    -I$(srcdir)/../../OTHERS/ASL/solvers
CXXFLAGS += $(AMPL_CPPFLAGS)

# Astyle stuff

ASTYLE_FILES = \
	CompositeNLP.hppbak CompositeNLP.cppbak \
	comp_ampl_ipopt.cppbak

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

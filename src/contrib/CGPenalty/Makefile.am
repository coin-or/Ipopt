# Copyright (C) 2007 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.

# Authors:  Andreas Waechter        IBM       2007-06-01

AUTOMAKE_OPTIONS = foreign

noinst_LTLIBRARIES = libcgpenalty.la

libcgpenalty_la_SOURCES = \
	IpCGPenaltyCq.cpp \
	IpCGPenaltyData.cpp \
	IpCGPenaltyLSAcceptor.cpp \
	IpCGPenaltyRegOp.cpp \
	IpCGPerturbationHandler.cpp \
	IpCGSearchDirCalc.cpp \
	IpPiecewisePenalty.cpp

libcgpenalty_la_LDFLAGS = $(LT_LDFLAGS)

AM_CPPFLAGS = -DIPOPT_BUILD \
	-I$(srcdir)/../../Common \
	-I$(srcdir)/../../LinAlg \
	-I$(srcdir)/../../LinAlg/TMatrices \
	-I$(srcdir)/../../Interfaces \
	-I$(srcdir)/../../Algorithm

# Astyle stuff

ASTYLE_FILES = \
	IpCGPenaltyCq.cppbak IpCGPenaltyCq.hppbak \
	IpCGPenaltyData.cppbak IpCGPenaltyData.hppbak \
	IpCGPenaltyLSAcceptor.cppbak IpCGPenaltyLSAcceptor.hppbak \
	IpCGPenaltyRegOp.cppbak IpCGPenaltyRegOp.hppbak \
	IpCGPerturbationHandler.cppbak IpCGPerturbationHandler.hppbak \
	IpCGSearchDirCalc.cppbak IpCGSearchDirCalc.hppbak \
	IpPiecewisePenalty.cppbak IpPiecewisePenalty.hppbak

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

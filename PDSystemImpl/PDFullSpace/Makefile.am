# Copyright (C) 2004, International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Common Public License.

## $Id$

# Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

AUTOMAKE_OPTIONS = foreign

noinst_LIBRARIES = libpdfullspace.a

libpdfullspace_a_SOURCES = \
	IpAugRestoSystemSolver.cpp IpAugRestoSystemSolver.hpp \
	IpAugSystemSolver.hpp \
	IpStdAugSystemSolver.hpp IpStdAugSystemSolver.cpp \
	IpAugTSystemSolver.cpp IpAugTSystemSolver.hpp \
	IpLeastSquareMults.hpp IpLeastSquareMults.cpp \
	IpMa27SymLinearSolver.hpp IpMa27SymLinearSolver.cpp \
	IpMc19SymTScalingMethod.hpp IpMc19SymTScalingMethod.cpp \
	IpPDFullSpaceSolver.cpp	IpPDFullSpaceSolver.hpp \
	IpTripletHelper.hpp IpTripletHelper.cpp \
	IpSymLinearSolver.hpp 


CXXFLAGS += -I$(srcdir)/../../Common \
	-I$(srcdir)/../../Interfaces \
	-I$(srcdir)/../../LinAlgImpl/Serial \
	-I$(srcdir)/../../Algorithm

# Astyle stuff

ASTYLE_SOURCES = \
	IpAugRestoSystemSolver.cppbak IpAugRestoSystemSolver.hppbak \
	IpAugSystemSolver.hppbak \
	IpStdAugSystemSolver.hppbak IpStdAugSystemSolver.cppbak \
	IpAugTSystemSolver.cppbak IpAugTSystemSolver.hppbak \
	IpLeastSquareMults.hppbak IpLeastSquareMults.cppbak \
	IpMa27SymLinearSolver.hppbak IpMa27SymLinearSolver.cppbak \
	IpMc19SymTScalingMethod.hppbak IpMc19SymTScalingMethod.cppbak \
	IpPDFullSpaceSolver.cppbak IpPDFullSpaceSolver.hppbak \
	IpTripletHelper.hppbak IpTripletHelper.cppbak \
	IpSymLinearSolver.hppbak

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

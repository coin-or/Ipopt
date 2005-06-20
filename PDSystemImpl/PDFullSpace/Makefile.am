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
	IpLeastSquareMults.hpp IpLeastSquareMults.cpp \
	IpMc19TSymScalingMethod.hpp IpMc19TSymScalingMethod.cpp \
	IpPDFullSpaceSolver.cpp	IpPDFullSpaceSolver.hpp \
	IpTripletHelper.hpp IpTripletHelper.cpp \
	IpSymLinearSolver.hpp \
	IpTSymScalingMethod.hpp \
	IpTSymLinearSolver.hpp \
	IpTSymLinearSolver.cpp \
	IpSparseSymLinearSolverInterface.hpp \
	IpMa27TSolverInterface.cpp IpMa27TSolverInterface.hpp \
	IpTripletToCSRConverter.cpp IpTripletToCSRConverter.hpp

if HAVE_PARDISO
  libpdfullspace_a_SOURCES += \
	IpPardisoSolverInterface.cpp IpPardisoSolverInterface.hpp
endif

if HAVE_TAUCS
  libpdfullspace_a_SOURCES += \
	IpTAUCSSolverInterface.hpp IpTAUCSSolverInterface.cpp
endif


CXXFLAGS += -I$(srcdir)/../../Common \
	-I$(srcdir)/../../Interfaces \
	-I$(srcdir)/../../LinAlgImpl/Serial \
	-I$(srcdir)/../../Algorithm

# Astyle stuff

ASTYLE_FILES = \
	IpAugRestoSystemSolver.cppbak IpAugRestoSystemSolver.hppbak \
	IpAugSystemSolver.hppbak \
	IpStdAugSystemSolver.hppbak IpStdAugSystemSolver.cppbak \
	IpLeastSquareMults.hppbak IpLeastSquareMults.cppbak \
	IpMc19TSymScalingMethod.hppbak IpMc19TSymScalingMethod.cppbak \
	IpPDFullSpaceSolver.cppbak IpPDFullSpaceSolver.hppbak \
	IpTripletHelper.hppbak IpTripletHelper.cppbak \
	IpSymLinearSolver.hppbak \
	IpTSymScalingMethod.hppbak \
	IpTSymLinearSolver.hppbak \
	IpTSymLinearSolver.cppbak \
	IpSparseSymLinearSolverInterface.hppbak \
	IpMa27TSolverInterface.cppbak IpMa27TSolverInterface.hppbak \
	IpTripletToCSRConverter.cppbak IpTripletToCSRConverter.hppbak \
	IpPardisoSolverInterface.cppbak IpPardisoSolverInterface.hppbak \
	IpTAUCSSolverInterface.hppbak IpTAUCSSolverInterface.cppbak

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

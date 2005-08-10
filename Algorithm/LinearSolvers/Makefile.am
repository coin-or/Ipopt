# Copyright (C) 2004, International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Common Public License.

## $Id$

# Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

AUTOMAKE_OPTIONS = foreign

noinst_LIBRARIES = liblinsolvers.a

liblinsolvers_a_SOURCES = \
	IpMa27TSolverInterface.cpp IpMa27TSolverInterface.hpp \
	IpMc19TSymScalingMethod.cpp IpMc19TSymScalingMethod.hpp \
	IpSparseSymLinearSolverInterface.hpp \
	IpSymLinearSolver.hpp \
	IpTripletToCSRConverter.cpp IpTripletToCSRConverter.hpp \
	IpTSymLinearSolver.cpp IpTSymLinearSolver.hpp \
	IpTSymScalingMethod.hpp

if HAVE_PARDISO
  liblinsolvers_a_SOURCES += \
	IpPardisoSolverInterface.cpp IpPardisoSolverInterface.hpp
endif

if HAVE_TAUCS
  liblinsolvers_a_SOURCES += \
	IpTAUCSSolverInterface.cpp IpTAUCSSolverInterface.hpp
endif


CXXFLAGS += -I$(srcdir)/../../Common \
	-I$(srcdir)/../../LinAlg \
	-I$(srcdir)/../../LinAlg/TMatrices \
	-I$(srcdir)/.. \
	-I$(srcdir)/../../Interfaces

# Astyle stuff

ASTYLE_FILES = \
	IpMa27TSolverInterface.cppbak IpMa27TSolverInterface.hppbak \
	IpMc19TSymScalingMethod.cppbak IpMc19TSymScalingMethod.hppbak \
	IpSparseSymLinearSolverInterface.hppbak \
	IpSymLinearSolver.hppbak \
	IpTripletToCSRConverter.cppbak IpTripletToCSRConverter.hppbak \
	IpTSymLinearSolver.cppbak IpTSymLinearSolver.hppbak \
	IpTSymScalingMethod.hppbak
	IpPardisoSolverInterface.cppbak IpPardisoSolverInterface.hppbak
	IpTAUCSSolverInterface.cppbak IpTAUCSSolverInterface.hppbak


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

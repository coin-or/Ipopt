# Copyright (C) 2004, 2011 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.

# Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

AUTOMAKE_OPTIONS = foreign

includeipoptdir = $(includedir)/coin
includeipopt_HEADERS = IpSymLinearSolver.hpp

noinst_LTLIBRARIES = liblinsolvers.la

liblinsolvers_la_SOURCES = \
	IpLinearSolversRegOp.cpp \
	IpSlackBasedTSymScalingMethod.cpp \
	IpTripletToCSRConverter.cpp \
	IpTSymDependencyDetector.cpp \
	IpTSymLinearSolver.cpp

liblinsolvers_la_SOURCES += \
	IpMa27TSolverInterface.cpp \
	IpMa57TSolverInterface.cpp \
	IpMa86SolverInterface.cpp \
	IpMa97SolverInterface.cpp \
	IpMc19TSymScalingMethod.cpp \
	IpMa77SolverInterface.cpp

if BUILD_LINEARSOLVERLOADER
  liblinsolvers_la_SOURCES += IpPardisoSolverInterface.cpp
else
if HAVE_PARDISO
  liblinsolvers_la_SOURCES += IpPardisoSolverInterface.cpp
endif
endif # BUILD_LINEARSOLVERLOADER

if HAVE_MA28
  liblinsolvers_la_SOURCES += IpMa28TDependencyDetector.cpp
  liblinsolvers_la_SOURCES += IpMa28Partition.F
endif

if HAVE_WSMP
  liblinsolvers_la_SOURCES += IpWsmpSolverInterface.cpp IpIterativeWsmpSolverInterface.cpp
endif

if COIN_HAS_MUMPS
  liblinsolvers_la_SOURCES += IpMumpsSolverInterface.cpp
endif

liblinsolvers_la_LDFLAGS = $(LT_LDFLAGS)

AM_CPPFLAGS = -DIPOPT_BUILD \
	-I$(srcdir)/../../Common \
	-I$(srcdir)/../../LinAlg \
	-I$(srcdir)/../../LinAlg/TMatrices \
	-I$(srcdir)/.. \
	-I$(srcdir)/../../Interfaces \
	-I$(srcdir)/../../contrib/CGPenalty \
	$(IPOPTLIB_CFLAGS)

# Astyle stuff

ASTYLE_FILES = \
	IpGenKKTSolverInterface.hppbak \
	IpIterativeWsmpSolverInterface.cppbak \
	IpIterativeWsmpSolverInterface.hppbak \
	IpLinearSolversRegOp.cppbak IpLinearSolversRegOp.hppbak \
	IpMa27TSolverInterface.cppbak IpMa27TSolverInterface.hppbak \
	IpMa28TDependencyDetector.cppbak IpMa28TDependencyDetector.hppbak \
	IpMa57TSolverInterface.cppbak IpMa57TSolverInterface.hppbak \
	IpMa77SolverInterface.cppbak IpMa77SolverInterface.hppbak \
	IpMa86SolverInterface.cppbak IpMa86SolverInterface.hppbak \
	IpMa97SolverInterface.cppbak IpMa97SolverInterface.hppbak \
	IpMc19TSymScalingMethod.cppbak IpMc19TSymScalingMethod.hppbak \
	IpMumpsSolverInterface.cppbak IpMumpsSolverInterface.hppbak \
	IpSlackBasedTSymScalingMethod.cppbak IpSlackBasedTSymScalingMethod.hppbak \
	IpSparseSymLinearSolverInterface.hppbak \
	IpSymLinearSolver.hppbak \
	IpTDependencyDetector.hppbak \
	IpTripletToCSRConverter.cppbak IpTripletToCSRConverter.hppbak \
	IpTSymDependencyDetector.cppbak IpTSymDependencyDetector.hppbak \
	IpTSymLinearSolver.cppbak IpTSymLinearSolver.hppbak \
	IpTSymScalingMethod.hppbak \
	IpPardisoSolverInterface.cppbak IpPardisoSolverInterface.hppbak \
	IpWsmpSolverInterface.cppbak IpWsmpSolverInterface.hppbak


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

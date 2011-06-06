# Copyright (C) 2008 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.

## $Id$

# Authors:  Andreas Waechter        IBM       2008-08-31

AUTOMAKE_OPTIONS = foreign

if BUILD_INEXACT
  noinst_LTLIBRARIES = libinexact.la
endif

libinexact_la_SOURCES = \
	IpInexactAlgBuilder.cpp IpInexactAlgBuilder.hpp \
	IpInexactCq.cpp IpInexactCq.hpp \
	IpInexactData.cpp IpInexactData.hpp \
	IpInexactDoglegNormal.cpp IpInexactDoglegNormal.hpp \
	IpInexactLSAcceptor.cpp IpInexactLSAcceptor.hpp \
	IpInexactNewtonNormal.cpp IpInexactNewtonNormal.hpp \
	IpInexactNormalStepCalc.hpp \
	IpInexactNormalTerminationTester.cpp IpInexactNormalTerminationTester.hpp \
	IpInexactPDSolver.cpp IpInexactPDSolver.hpp \
	IpInexactPDTerminationTester.cpp IpInexactPDTerminationTester.hpp \
	IpInexactRegOp.cpp IpInexactRegOp.hpp \
	IpInexactSearchDirCalc.cpp IpInexactSearchDirCalc.hpp \
	IpInexactTSymScalingMethod.cpp IpInexactTSymScalingMethod.hpp \
	IpIterativePardisoSolverInterface.cpp IpIterativePardisoSolverInterface.hpp \
	IpIterativeSolverTerminationTester.cpp IpIterativeSolverTerminationTester.hpp

libinexact_la_LDFLAGS = $(LT_LDFLAGS)

AM_CPPFLAGS = -I`$(CYGPATH_W) $(srcdir)/../../Common` \
	-I`$(CYGPATH_W) $(srcdir)/../../LinAlg` \
	-I`$(CYGPATH_W) $(srcdir)/../../LinAlg/TMatrices` \
	-I`$(CYGPATH_W) $(srcdir)/../../Interfaces` \
	-I`$(CYGPATH_W) $(srcdir)/../` \
	-I`$(CYGPATH_W) $(srcdir)/../LinearSolvers` \
	-I`$(CYGPATH_W) $(srcdir)/../../contrib/LinearSolverLoader`

if COIN_HAS_HSL
AM_CPPFLAGS += $(HSL_CFLAGS)
endif

# This line is necessary to allow VPATH compilation
DEFAULT_INCLUDES = -I. -I`$(CYGPATH_W) $(srcdir)` -I$(top_builddir)/src/Common

# Astyle stuff

ASTYLE_FILES = \
	IpInexactAlgBuilder.cppbak IpInexactAlgBuilder.hppbak \
	IpInexactCq.cppbak IpInexactCq.hppbak \
	IpInexactData.cppbak IpInexactData.hppbak \
	IpInexactDoglegNormal.cppbak IpInexactDoglegNormal.hppbak \
	IpInexactLSAcceptor.cppbak IpInexactLSAcceptor.hppbak \
	IpInexactNewtonNormal.cppbak IpInexactNewtonNormal.hppbak \
	IpInexactNormalStepCalc.hppbak \
	IpInexactNormalTerminationTester.cppbak IpInexactNormalTerminationTester.hppbak \
	IpInexactPDSolver.cppbak IpInexactPDSolver.hppbak \
	IpInexactPDTerminationTester.cppbak IpInexactPDTerminationTester.hppbak \
	IpInexactRegOp.cppbak IpInexactRegOp.hppbak \
	IpInexactSearchDirCalc.cppbak IpInexactSearchDirCalc.hppbak \
	IpInexactTSymScalingMethod.cppbak IpInexactTSymScalingMethod.hppbak \
	IpIterativePardisoSolverInterface.cppbak IpIterativePardisoSolverInterface.hppbak \
	IpIterativeSolverTerminationTester.cppbak IpIterativeSolverTerminationTester.hppbak

ASTYLE = @ASTYLE@
ASTYLEFLAGS = @ASTYLEFLAGS@

DISTCLEANFILES = $(ASTYLE_FILES)

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

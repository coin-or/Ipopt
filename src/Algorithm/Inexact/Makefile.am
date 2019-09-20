# Copyright (C) 2008 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.

# Authors:  Andreas Waechter        IBM       2008-08-31

AUTOMAKE_OPTIONS = foreign

if BUILD_INEXACT
  noinst_LTLIBRARIES = libinexact.la
endif

libinexact_la_SOURCES = \
	IpInexactAlgBuilder.cpp \
	IpInexactCq.cpp \
	IpInexactData.cpp \
	IpInexactDoglegNormal.cpp \
	IpInexactLSAcceptor.cpp \
	IpInexactNewtonNormal.cpp \
	IpInexactNormalTerminationTester.cpp \
	IpInexactPDSolver.cpp \
	IpInexactPDTerminationTester.cpp \
	IpInexactRegOp.cpp \
	IpInexactSearchDirCalc.cpp \
	IpInexactTSymScalingMethod.cpp \
	IpIterativePardisoSolverInterface.cpp \
	IpIterativeSolverTerminationTester.cpp

#libinexact_la_LDFLAGS = $(LT_LDFLAGS)

AM_CPPFLAGS = -DIPOPT_BUILD \
	-I$(srcdir)/../../Common \
	-I$(srcdir)/../../LinAlg \
	-I$(srcdir)/../../LinAlg/TMatrices \
	-I$(srcdir)/../../Interfaces \
	-I$(srcdir)/.. \
	-I$(srcdir)/../LinearSolvers \
	-I$(srcdir)/../../contrib/LinearSolverLoader \
	$(IPOPTLIB_CFLAGS)


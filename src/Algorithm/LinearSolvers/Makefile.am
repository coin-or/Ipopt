# Copyright (C) 2004, 2011 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.

# Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

AUTOMAKE_OPTIONS = foreign

includeipoptdir = $(includedir)/coin-or
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

#liblinsolvers_la_LDFLAGS = $(LT_LDFLAGS)

AM_CPPFLAGS = -DIPOPT_BUILD \
	-I$(srcdir)/../../Common \
	-I$(srcdir)/../../LinAlg \
	-I$(srcdir)/../../LinAlg/TMatrices \
	-I$(srcdir)/.. \
	-I$(srcdir)/../../Interfaces \
	-I$(srcdir)/../../contrib/CGPenalty \
	$(IPOPTLIB_CFLAGS)

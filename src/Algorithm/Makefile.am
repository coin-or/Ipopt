# Copyright (C) 2004, 2010 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.

# Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

AUTOMAKE_OPTIONS = foreign

SUBDIRS = LinearSolvers Inexact

includeipoptdir = $(includedir)/coin-or
includeipopt_HEADERS = \
	IpAlgBuilder.hpp \
	IpAlgStrategy.hpp \
	IpAugSystemSolver.hpp \
	IpConvCheck.hpp \
	IpEqMultCalculator.hpp \
	IpHessianUpdater.hpp \
	IpIpoptAlg.hpp \
	IpIpoptCalculatedQuantities.hpp \
	IpIpoptData.hpp \
	IpIpoptNLP.hpp \
	IpIterateInitializer.hpp \
	IpIteratesVector.hpp \
	IpIterationOutput.hpp \
	IpOrigIpoptNLP.hpp \
	IpLineSearch.hpp \
	IpMuUpdate.hpp \
	IpNLPScaling.hpp \
	IpPDSystemSolver.hpp \
	IpSearchDirCalculator.hpp \
	IpTimingStatistics.hpp

noinst_LTLIBRARIES = libipoptalg.la

libipoptalg_la_SOURCES = \
	IpAdaptiveMuUpdate.cpp \
	IpAlgBuilder.cpp \
	IpAlgorithmRegOp.cpp \
	IpAugRestoSystemSolver.cpp \
	IpBacktrackingLineSearch.cpp \
	IpDefaultIterateInitializer.cpp \
	IpEquilibrationScaling.cpp \
	IpExactHessianUpdater.cpp \
	IpFilter.cpp \
	IpFilterLSAcceptor.cpp \
	IpGenAugSystemSolver.cpp \
	IpGradientScaling.cpp \
	IpIpoptAlg.cpp \
	IpIpoptCalculatedQuantities.cpp \
	IpIpoptData.cpp \
	IpIteratesVector.cpp \
	IpLeastSquareMults.cpp \
	IpLimMemQuasiNewtonUpdater.cpp \
	IpLoqoMuOracle.cpp \
	IpLowRankAugSystemSolver.cpp \
	IpLowRankSSAugSystemSolver.cpp \
	IpMonotoneMuUpdate.cpp \
	IpNLPBoundsRemover.cpp \
	IpNLPScaling.cpp \
	IpOptErrorConvCheck.cpp \
	IpOrigIpoptNLP.cpp \
	IpOrigIterationOutput.cpp \
	IpPDFullSpaceSolver.cpp \
	IpPDPerturbationHandler.cpp \
	IpPDSearchDirCalc.cpp \
	IpPenaltyLSAcceptor.cpp \
	IpProbingMuOracle.cpp \
	IpQualityFunctionMuOracle.cpp \
	IpRestoConvCheck.cpp \
	IpRestoFilterConvCheck.cpp \
	IpRestoIpoptNLP.cpp \
	IpRestoIterateInitializer.cpp \
	IpRestoIterationOutput.cpp \
	IpRestoMinC_1Nrm.cpp \
	IpRestoPenaltyConvCheck.cpp \
	IpRestoRestoPhase.cpp \
	IpStdAugSystemSolver.cpp \
	IpTimingStatistics.cpp \
	IpUserScaling.cpp \
	IpWarmStartIterateInitializer.cpp

#libipoptalg_la_LDFLAGS = $(LT_LDFLAGS)

AM_CPPFLAGS = -DIPOPT_BUILD \
	-I$(srcdir)/../Common \
	-I$(srcdir)/../LinAlg \
	-I$(srcdir)/../LinAlg/TMatrices \
	-I$(srcdir)/../Interfaces \
	-I$(srcdir)/LinearSolvers \
	-I$(srcdir)/../contrib/CGPenalty \
	-I$(srcdir)/../contrib/LinearSolverLoader \
	$(IPOPTLIB_CFLAGS)

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

libipoptalg_la_LDFLAGS = $(LT_LDFLAGS)

AM_CPPFLAGS = -DIPOPT_BUILD \
	-I$(srcdir)/../Common \
	-I$(srcdir)/../LinAlg \
	-I$(srcdir)/../LinAlg/TMatrices \
	-I$(srcdir)/../Interfaces \
	-I$(srcdir)/LinearSolvers \
	-I$(srcdir)/../contrib/CGPenalty \
	-I$(srcdir)/../contrib/LinearSolverLoader \
	$(IPOPTLIB_CFLAGS)

# Astyle stuff

ASTYLE_FILES = \
	IpAdaptiveMuUpdate.cppbak IpAdaptiveMuUpdate.hppbak \
	IpAlgBuilder.cppbak IpAlgBuilder.hppbak \
	IpAlgorithmRegOp.cppbak IpAlgorithmRegOp.hppbak \
	IpAlgStrategy.hppbak \
	IpAugRestoSystemSolver.cppbak IpAugRestoSystemSolver.hppbak \
	IpAugSystemSolver.hppbak \
	IpBacktrackingLSAcceptor.hppbak \
	IpBacktrackingLineSearch.cppbak IpBacktrackingLineSearch.hppbak \
	IpConvCheck.hppbak \
	IpDefaultIterateInitializer.cppbak IpDefaultIterateInitializer.hppbak \
	IpEqMultCalculator.hppbak \
	IpEquilibrationScaling.cppbak IpEquilibrationScaling.hppbak \
	IpExactHessianUpdater.cppbak IpExactHessianUpdater.hppbak \
	IpFilter.cppbak IpFilter.hppbak \
	IpFilterLSAcceptor.cppbak IpFilterLSAcceptor.hppbak \
	IpGenAugSystemSolver.cppbak IpGenAugSystemSolver.hppbak \
	IpGradientScaling.cppbak IpGradientScaling.hppbak \
	IpHessianUpdater.hppbak \
	IpIpoptAlg.cppbak IpIpoptAlg.hppbak \
	IpIpoptCalculatedQuantities.cppbak IpIpoptCalculatedQuantities.hppbak \
	IpIpoptData.cppbak IpIpoptData.hppbak \
	IpIpoptNLP.hppbak \
	IpIterateInitializer.hppbak \
	IpIteratesVector.cppbak IpIteratesVector.hppbak \
	IpIterationOutput.hppbak \
	IpLeastSquareMults.cppbak IpLeastSquareMults.hppbak \
	IpLimMemQuasiNewtonUpdater.cppbak IpLimMemQuasiNewtonUpdater.hppbak \
	IpLineSearch.hppbak \
	IpLoqoMuOracle.cppbak IpLoqoMuOracle.hppbak \
	IpLowRankAugSystemSolver.cppbak IpLowRankAugSystemSolver.hppbak \
	IpLowRankSSAugSystemSolver.cppbak IpLowRankSSAugSystemSolver.hppbak \
	IpMonotoneMuUpdate.cppbak IpMonotoneMuUpdate.hppbak \
	IpMuOracle.hppbak \
	IpMuUpdate.hppbak \
	IpNLPBoundsRemover.cppbak IpNLPBoundsRemover.hppbak \
	IpNLPScaling.cppbak IpNLPScaling.hppbak \
	IpOptErrorConvCheck.cppbak IpOptErrorConvCheck.hppbak \
	IpOrigIpoptNLP.cppbak IpOrigIpoptNLP.hppbak \
	IpOrigIterationOutput.cppbak IpOrigIterationOutput.hppbak \
	IpPDFullSpaceSolver.cppbak IpPDFullSpaceSolver.hppbak \
	IpPDPerturbationHandler.cppbak IpPDPerturbationHandler.hppbak \
	IpPDSearchDirCalc.cppbak IpPDSearchDirCalc.hppbak \
	IpPDSystemSolver.hppbak \
	IpPenaltyLSAcceptor.cppbak IpPenaltyLSAcceptor.hppbak \
	IpProbingMuOracle.cppbak IpProbingMuOracle.hppbak \
	IpQualityFunctionMuOracle.cppbak IpQualityFunctionMuOracle.hppbak \
	IpRestoConvCheck.cppbak IpRestoConvCheck.hppbak \
	IpRestoFilterConvCheck.cppbak IpRestoFilterConvCheck.hppbak \
	IpRestoIpoptNLP.cppbak IpRestoIpoptNLP.hppbak \
	IpRestoIterateInitializer.cppbak IpRestoIterateInitializer.hppbak \
	IpRestoIterationOutput.cppbak IpRestoIterationOutput.hppbak \
	IpRestoMinC_1Nrm.cppbak IpRestoMinC_1Nrm.hppbak \
	IpRestoPenaltyConvCheck.cppbak IpRestoPenaltyConvCheck.hppbak \
	IpRestoPhase.hppbak \
	IpRestoRestoPhase.cppbak IpRestoRestoPhase.hppbak \
	IpSearchDirCalculator.hppbak \
	IpStdAugSystemSolver.cppbak IpStdAugSystemSolver.hppbak \
	IpTimingStatistics.cppbak IpTimingStatistics.hppbak \
	IpUserScaling.cppbak IpUserScaling.hppbak \
	IpWarmStartIterateInitializer.cppbak IpWarmStartIterateInitializer.hppbak

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

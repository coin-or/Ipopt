# Copyright (C) 2004, 2010 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.

## $Id$

# Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

AUTOMAKE_OPTIONS = foreign

SUBDIRS = LinearSolvers Inexact

includeipoptdir = $(includedir)/coin
includeipopt_HEADERS = IpIpoptCalculatedQuantities.hpp \
	IpIpoptData.hpp IpIteratesVector.hpp IpTimingStatistics.hpp \
	IpIpoptNLP.hpp IpOrigIpoptNLP.hpp IpNLPScaling.hpp \
	IpAlgBuilder.hpp IpIpoptAlg.hpp IpAlgStrategy.hpp \
	IpSearchDirCalculator.hpp IpLineSearch.hpp IpMuUpdate.hpp \
	IpConvCheck.hpp IpIterateInitializer.hpp IpIterationOutput.hpp \
	IpHessianUpdater.hpp IpEqMultCalculator.hpp IpAugSystemSolver.hpp \
	IpPDSystemSolver.hpp

noinst_LTLIBRARIES = libipoptalg.la

libipoptalg_la_SOURCES = \
	IpAdaptiveMuUpdate.cpp IpAdaptiveMuUpdate.hpp \
	IpAlgBuilder.cpp IpAlgBuilder.hpp \
	IpAlgorithmRegOp.cpp IpAlgorithmRegOp.hpp \
	IpAlgStrategy.hpp \
	IpAugRestoSystemSolver.cpp IpAugRestoSystemSolver.hpp \
	IpAugSystemSolver.hpp \
	IpBacktrackingLSAcceptor.hpp \
	IpBacktrackingLineSearch.cpp IpBacktrackingLineSearch.hpp \
	IpConvCheck.hpp \
	IpDefaultIterateInitializer.cpp IpDefaultIterateInitializer.hpp \
	IpEqMultCalculator.hpp \
	IpEquilibrationScaling.cpp IpEquilibrationScaling.hpp \
	IpExactHessianUpdater.cpp IpExactHessianUpdater.hpp \
	IpFilter.cpp IpFilter.hpp \
	IpFilterLSAcceptor.cpp IpFilterLSAcceptor.hpp \
	IpGenAugSystemSolver.cpp IpGenAugSystemSolver.hpp \
	IpGradientScaling.cpp IpGradientScaling.hpp \
	IpHessianUpdater.hpp \
	IpIpoptAlg.cpp IpIpoptAlg.hpp \
	IpIpoptCalculatedQuantities.cpp IpIpoptCalculatedQuantities.hpp \
	IpIpoptData.cpp IpIpoptData.hpp \
	IpIpoptNLP.hpp \
	IpIterateInitializer.hpp \
	IpIteratesVector.cpp IpIteratesVector.hpp \
	IpIterationOutput.hpp \
	IpLeastSquareMults.cpp IpLeastSquareMults.hpp \
	IpLimMemQuasiNewtonUpdater.cpp IpLimMemQuasiNewtonUpdater.hpp \
	IpLineSearch.hpp \
	IpLoqoMuOracle.cpp IpLoqoMuOracle.hpp \
	IpLowRankAugSystemSolver.cpp IpLowRankAugSystemSolver.hpp \
	IpLowRankSSAugSystemSolver.cpp IpLowRankSSAugSystemSolver.hpp \
	IpMonotoneMuUpdate.cpp IpMonotoneMuUpdate.hpp \
	IpMuOracle.hpp \
	IpMuUpdate.hpp \
	IpNLPBoundsRemover.cpp IpNLPBoundsRemover.hpp \
	IpNLPScaling.cpp IpNLPScaling.hpp \
	IpOptErrorConvCheck.cpp IpOptErrorConvCheck.hpp \
	IpOrigIpoptNLP.cpp IpOrigIpoptNLP.hpp \
	IpOrigIterationOutput.cpp IpOrigIterationOutput.hpp \
	IpPDFullSpaceSolver.cpp IpPDFullSpaceSolver.hpp \
	IpPDPerturbationHandler.cpp IpPDPerturbationHandler.hpp \
	IpPDSearchDirCalc.cpp IpPDSearchDirCalc.hpp \
	IpPDSystemSolver.hpp \
	IpPenaltyLSAcceptor.cpp IpPenaltyLSAcceptor.hpp \
	IpProbingMuOracle.cpp IpProbingMuOracle.hpp \
	IpQualityFunctionMuOracle.cpp IpQualityFunctionMuOracle.hpp \
	IpRestoConvCheck.cpp IpRestoConvCheck.hpp \
	IpRestoFilterConvCheck.cpp IpRestoFilterConvCheck.hpp \
	IpRestoIpoptNLP.cpp IpRestoIpoptNLP.hpp \
	IpRestoIterateInitializer.cpp IpRestoIterateInitializer.hpp \
	IpRestoIterationOutput.cpp IpRestoIterationOutput.hpp \
	IpRestoMinC_1Nrm.cpp IpRestoMinC_1Nrm.hpp \
	IpRestoPenaltyConvCheck.cpp IpRestoPenaltyConvCheck.hpp \
	IpRestoPhase.hpp \
	IpRestoRestoPhase.cpp IpRestoRestoPhase.hpp \
	IpSearchDirCalculator.hpp \
	IpStdAugSystemSolver.cpp IpStdAugSystemSolver.hpp \
	IpTimingStatistics.cpp IpTimingStatistics.hpp \
	IpUserScaling.cpp IpUserScaling.hpp \
	IpWarmStartIterateInitializer.cpp IpWarmStartIterateInitializer.hpp

libipoptalg_la_LDFLAGS = $(LT_LDFLAGS)

AM_CPPFLAGS = -I`$(CYGPATH_W) $(srcdir)/../Common` \
	-I`$(CYGPATH_W) $(srcdir)/../LinAlg` \
	-I`$(CYGPATH_W) $(srcdir)/../LinAlg/TMatrices` \
	-I`$(CYGPATH_W) $(srcdir)/../Interfaces` \
	-I`$(CYGPATH_W) $(srcdir)/LinearSolvers` \
	-I`$(CYGPATH_W) $(srcdir)/../contrib/CGPenalty` \
	-I`$(CYGPATH_W) $(srcdir)/../contrib/LinearSolverLoader`

#if COIN_HAS_MUMPS
#  AM_CPPFLAGS += $(MUMPS_CFLAGS)
#endif

if COIN_HAS_HSL
AM_CPPFLAGS += $(HSL_CFLAGS)
endif

# This line is necessary to allow VPATH compilation
DEFAULT_INCLUDES = -I. -I`$(CYGPATH_W) $(srcdir)` -I$(top_builddir)/src/Common

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
	IpLowRankAugSystemSolver.cppbak IpLowRankAugSystemSolver.cppbak \
	IpLowRankSSAugSystemSolver.cppbak IpLowRankSSAugSystemSolver.cppbak \
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
	$(ASTYLE) $(ASTYLEFLAGS) < $@ > $<
	touch $@

.cpp.cppbak:
	mv $< $@
	$(ASTYLE) $(ASTYLEFLAGS) < $@ > $<
	touch $@

# Copyright (C) 2004, 2006 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Common Public License.

## $Id$

# Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

AUTOMAKE_OPTIONS = foreign

SUBDIRS = LinearSolvers

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
	IpLineSearch.hpp \
	IpLoqoMuOracle.cpp IpLoqoMuOracle.hpp \
	IpLowRankAugSystemSolver.cpp IpLowRankAugSystemSolver.hpp \
	IpMonotoneMuUpdate.cpp IpMonotoneMuUpdate.hpp \
	IpMuOracle.hpp \
	IpMuUpdate.hpp \
	IpNLPScaling.cpp IpNLPScaling.hpp \
	IpOptErrorConvCheck.cpp IpOptErrorConvCheck.hpp \
	IpOrigIpoptNLP.cpp IpOrigIpoptNLP.hpp \
	IpOrigIterationOutput.cpp IpOrigIterationOutput.hpp \
	IpLimMemQuasiNewtonUpdater.cpp IpLimMemQuasiNewtonUpdater.hpp \
	IpPDFullSpaceSolver.cpp IpPDFullSpaceSolver.hpp \
	IpPDPerturbationHandler.cpp IpPDPerturbationHandler.hpp \
	IpPDSystemSolver.hpp \
	IpProbingMuOracle.cpp IpProbingMuOracle.hpp \
	IpQualityFunctionMuOracle.cpp IpQualityFunctionMuOracle.hpp \
	IpRestoFilterConvCheck.cpp IpRestoFilterConvCheck.hpp \
	IpRestoIpoptNLP.cpp IpRestoIpoptNLP.hpp \
	IpRestoIterateInitializer.cpp IpRestoIterateInitializer.hpp \
	IpRestoIterationOutput.cpp IpRestoIterationOutput.hpp \
	IpRestoMinC_1Nrm.cpp IpRestoMinC_1Nrm.hpp \
	IpRestoPhase.hpp \
	IpRestoRestoPhase.cpp IpRestoRestoPhase.hpp \
	IpStdAugSystemSolver.cpp IpStdAugSystemSolver.hpp \
	IpTimingStatistics.cpp IpTimingStatistics.hpp \
	IpUserScaling.cpp IpUserScaling.hpp \
	IpWarmStartIterateInitializer.cpp IpWarmStartIterateInitializer.hpp

libipoptalg_la_LDFLAGS = $(LT_LDFLAGS)

AM_CPPFLAGS = -I`$(CYGPATH_W) $(srcdir)/../Common` \
	-I`$(CYGPATH_W) $(srcdir)/../LinAlg` \
	-I`$(CYGPATH_W) $(srcdir)/../LinAlg/TMatrices` \
	-I`$(CYGPATH_W) $(srcdir)/../Interfaces` \
	-I`$(CYGPATH_W) $(srcdir)/LinearSolvers`

if COIN_HAS_MUMPS
  AM_CPPFLAGS += $(MUMPS_INCFLAGS)
endif

# This line is necessary to allow VPATH compilation with MS compilers
# on Cygwin
DEFAULT_INCLUDES = -I. -I`$(CYGPATH_W) $(srcdir)` -I$(top_builddir)/inc

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
	IpLineSearch.hppbak \
	IpLoqoMuOracle.cppbak IpLoqoMuOracle.hppbak \
	IpLowRankAugSystemSolver.cppbak IpLowRankAugSystemSolver.cppbak \
	IpMonotoneMuUpdate.cppbak IpMonotoneMuUpdate.hppbak \
	IpMuOracle.hppbak \
	IpMuUpdate.hppbak \
	IpNLPScaling.cppbak IpNLPScaling.hppbak \
	IpOptErrorConvCheck.cppbak IpOptErrorConvCheck.hppbak \
	IpOrigIpoptNLP.cppbak IpOrigIpoptNLP.hppbak \
	IpOrigIterationOutput.cppbak IpOrigIterationOutput.hppbak \
	IpLimMemQuasiNewtonUpdater.cppbak IpLimMemQuasiNewtonUpdater.hppbak \
	IpPDFullSpaceSolver.cppbak IpPDFullSpaceSolver.hppbak \
	IpPDPerturbationHandler.cppbak IpPDPerturbationHandler.hppbak \
	IpPDSystemSolver.hppbak \
	IpProbingMuOracle.cppbak IpProbingMuOracle.hppbak \
	IpQualityFunctionMuOracle.cppbak IpQualityFunctionMuOracle.hppbak \
	IpRestoFilterConvCheck.cppbak IpRestoFilterConvCheck.hppbak \
	IpRestoIpoptNLP.cppbak IpRestoIpoptNLP.hppbak \
	IpRestoIterateInitializer.cppbak IpRestoIterateInitializer.hppbak \
	IpRestoIterationOutput.cppbak IpRestoIterationOutput.hppbak \
	IpRestoMinC_1Nrm.cppbak IpRestoMinC_1Nrm.hppbak \
	IpRestoPhase.hppbak \
	IpRestoRestoPhase.cppbak IpRestoRestoPhase.hppbak \
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

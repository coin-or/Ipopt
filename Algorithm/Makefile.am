# Copyright (C) 2004, 2005 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Common Public License.

## $Id$

# Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

AUTOMAKE_OPTIONS = foreign

SUBDIRS = LinearSolvers

noinst_LIBRARIES = libipoptalg.a

libipoptalg_a_SOURCES = \
	IpAdaptiveMuUpdate.cpp IpAdaptiveMuUpdate.hpp \
	IpAlgBuilder.cpp IpAlgBuilder.hpp \
	IpAlgorithmRegOp.cpp IpAlgorithmRegOp.hpp \
	IpAlgStrategy.hpp \
	IpAugRestoSystemSolver.cpp IpAugRestoSystemSolver.hpp \
	IpAugSystemSolver.hpp \
	IpConvCheck.hpp \
	IpDefaultIterateInitializer.cpp IpDefaultIterateInitializer.hpp \
	IpEqMultCalculator.hpp \
	IpFilter.cpp IpFilter.hpp \
	IpFilterLineSearch.cpp IpFilterLineSearch.hpp \
	IpGradientScaling.cpp IpGradientScaling.hpp \
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
	IpMonotoneMuUpdate.cpp IpMonotoneMuUpdate.hpp \
	IpMuOracle.hpp \
	IpMuUpdate.hpp \
	IpNLPScaling.cpp IpNLPScaling.hpp \
	IpOptErrorConvCheck.cpp IpOptErrorConvCheck.hpp \
	IpOrigIpoptNLP.cpp IpOrigIpoptNLP.hpp \
	IpOrigIterationOutput.cpp IpOrigIterationOutput.hpp \
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
	IpUserScaling.cpp IpUserScaling.hpp \
	IpWarmStartIterateInitializer.cpp IpWarmStartIterateInitializer.hpp

CXXFLAGS += -I$(srcdir)/../Common \
	-I$(srcdir)/../LinAlg \
	-I$(srcdir)/../LinAlg/TMatrices \
	-I$(srcdir)/../Interfaces \
	-I$(srcdir)/LinearSolvers

# Astyle stuff

ASTYLE_FILES = \
	IpAdaptiveMuUpdate.cppbak IpAdaptiveMuUpdate.hppbak \
	IpAlgBuilder.cppbak IpAlgBuilder.hppbak \
	IpAlgorithmRegOp.cppbak IpAlgorithmRegOp.hppbak \
	IpAlgStrategy.hppbak \
	IpAugRestoSystemSolver.cppbak IpAugRestoSystemSolver.hppbak \
	IpAugSystemSolver.hppbak \
	IpConvCheck.hppbak \
	IpDefaultIterateInitializer.cppbak IpDefaultIterateInitializer.hppbak \
	IpEqMultCalculator.hppbak \
	IpFilter.cppbak IpFilter.hppbak \
	IpFilterLineSearch.cppbak IpFilterLineSearch.hppbak \
	IpGradientScaling.cppbak IpGradientScaling.hppbak \
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
	IpMonotoneMuUpdate.cppbak IpMonotoneMuUpdate.hppbak \
	IpMuOracle.hppbak \
	IpMuUpdate.hppbak \
	IpNLPScaling.cppbak IpNLPScaling.hppbak \
	IpOptErrorConvCheck.cppbak IpOptErrorConvCheck.hppbak \
	IpOrigIpoptNLP.cppbak IpOrigIpoptNLP.hppbak \
	IpOrigIterationOutput.cppbak IpOrigIterationOutput.hppbak \
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
	IpUserScaling.cppbak IpUserScaling.hppbak \
	IpWarmStartIterateInitializer.cppbak IpWarmStartIterateInitializer.hppbak

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

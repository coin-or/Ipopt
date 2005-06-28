# Copyright (C) 2004, International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Common Public License.

## $Id$

# Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

AUTOMAKE_OPTIONS = foreign

noinst_LIBRARIES = libipoptalg.a

libipoptalg_a_SOURCES = \
	IpAlgBuilder.hpp IpAlgBuilder.cpp \
	IpAlgStrategy.hpp \
	IpConvCheck.hpp \
	IpDefaultIterateInitializer.hpp IpDefaultIterateInitializer.cpp \
	IpWarmStartIterateInitializer.hpp IpWarmStartIterateInitializer.cpp \
	IpFilter.cpp IpFilter.hpp \
	IpFilterLineSearch.cpp IpFilterLineSearch.hpp \
	IpIpoptAlg.cpp IpIpoptAlg.hpp \
	IpIpoptCalculatedQuantities.cpp IpIpoptCalculatedQuantities.hpp \
	IpIpoptData.cpp IpIpoptData.hpp \
	IpIpoptNLP.hpp  \
	IpIterateInitializer.hpp \
	IpIteratesVector.hpp IpIteratesVector.cpp \
	IpIterationOutput.hpp \
	IpLineSearch.hpp \
	IpLoqoMuOracle.hpp IpLoqoMuOracle.cpp \
	IpMonotoneMuUpdate.hpp IpMonotoneMuUpdate.cpp \
	IpMuOracle.hpp \
	IpMuUpdate.hpp IpMuUpdate.cpp \
	IpNonmonotoneMuUpdate.hpp IpNonmonotoneMuUpdate.cpp \
	IpNLPScaling.hpp IpNLPScaling.cpp \
	IpOptErrorConvCheck.cpp IpOptErrorConvCheck.hpp \
	IpOrigIpoptNLP.cpp IpOrigIpoptNLP.hpp \
	IpOrigIterationOutput.hpp IpOrigIterationOutput.cpp \
	IpProbingMuOracle.hpp IpProbingMuOracle.cpp \
	IpQualityFunctionMuOracle.hpp IpQualityFunctionMuOracle.cpp \
	IpRestoFilterConvCheck.cpp IpRestoFilterConvCheck.hpp \
	IpRestoIpoptNLP.hpp IpRestoIpoptNLP.cpp \
	IpRestoIterateInitializer.hpp IpRestoIterateInitializer.cpp \
	IpRestoIterationOutput.hpp IpRestoIterationOutput.cpp \
	IpRestoMinC_1Nrm.cpp IpRestoMinC_1Nrm.hpp \
	IpRestoPhase.hpp \
	IpRestoRestoPhase.hpp IpRestoRestoPhase.cpp \
	IpUserScaling.cpp IpUserScaling.hpp

CXXFLAGS += -I$(srcdir)/../Common
CXXFLAGS += -I$(srcdir)/../Interfaces
CXXFLAGS += -I$(srcdir)/../PDSystemImpl/PDFullSpace
CXXFLAGS += -I$(srcdir)/../LinAlgImpl/Serial

# Astyle stuff

ASTYLE_FILES = \
	IpAlgBuilder.hppbak IpAlgBuilder.cppbak \
	IpAlgStrategy.hppbak \
	IpConvCheck.hppbak \
	IpDefaultIterateInitializer.hppbak IpDefaultIterateInitializer.cppbak \
	IpWarmStartIterateInitializer.hppbak IpWarmStartIterateInitializer.cppbak \
	IpFilter.cppbak IpFilter.hppbak \
	IpFilterLineSearch.cppbak IpFilterLineSearch.hppbak \
	IpIpoptAlg.cppbak IpIpoptAlg.hppbak \
	IpIpoptCalculatedQuantities.cppbak IpIpoptCalculatedQuantities.hppbak \
	IpIpoptData.cppbak IpIpoptData.hppbak \
	IpIpoptNLP.hppbak \
	IpIterateInitializer.hppbak \
	IpIteratesVector.hppbak IpIteratesVector.cppbak \
	IpIterationOutput.hppbak \
	IpLineSearch.hppbak \
	IpLoqoMuOracle.hppbak IpLoqoMuOracle.cppbak \
	IpMonotoneMuUpdate.hppbak IpMonotoneMuUpdate.cppbak \
	IpMuOracle.hppbak \
	IpMuUpdate.hppbak \
	IpNonmonotoneMuUpdate.hppbak IpNonmonotoneMuUpdate.cppbak \
	IpNLPScaling.hppbak IpNLPScaling.cppbak \
	IpOptErrorConvCheck.cppbak IpOptErrorConvCheck.hppbak \
	IpOrigIpoptNLP.cppbak IpOrigIpoptNLP.hppbak \
	IpOrigIterationOutput.hppbak IpOrigIterationOutput.cppbak \
	IpProbingMuOracle.hppbak IpProbingMuOracle.cppbak \
	IpQualityFunctionMuOracle.hppbak IpQualityFunctionMuOracle.cppbak \
	IpRestoFilterConvCheck.cppbak IpRestoFilterConvCheck.hppbak \
	IpRestoIpoptNLP.hppbak IpRestoIpoptNLP.cppbak \
	IpRestoIterateInitializer.hppbak IpRestoIterateInitializer.cppbak \
	IpRestoIterationOutput.hppbak IpRestoIterationOutput.cppbak \
	IpRestoMinC_1Nrm.cppbak IpRestoMinC_1Nrm.hppbak \
	IpRestoPhase.hppbak \
	IpRestoRestoPhase.hppbak IpRestoRestoPhase.cppbak \
	IpUserScaling.cppbak IpUserScaling.hppbak

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

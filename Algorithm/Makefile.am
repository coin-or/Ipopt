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
	IpFilter.cpp IpFilter.hpp \
	IpFilterLineSearch.cpp IpFilterLineSearch.hpp \
	IpIpoptAlg.cpp IpIpoptAlg.hpp \
	IpIpoptCalculatedQuantities.cpp IpIpoptCalculatedQuantities.hpp \
	IpIpoptData.cpp IpIpoptData.hpp \
	IpIpoptNLP.hpp \
	IpIterateInitializer.hpp \
	IpLineSearch.hpp \
	IpLoqoMuOracle.hpp IpLoqoMuOracle.cpp \
	IpMonotoneMuUpdate.hpp IpMonotoneMuUpdate.cpp \
	IpMuOracle.hpp \
	IpMuUpdate.hpp \
	IpNonmonotoneMuUpdate.hpp IpNonmonotoneMuUpdate.cpp \
	IpOptErrorConvCheck.cpp IpOptErrorConvCheck.hpp \
	IpOrigIpoptNLP.cpp IpOrigIpoptNLP.hpp \
	IpIterationOutput.hpp \
	IpOrigIterationOutput.hpp IpOrigIterationOutput.cpp \
	IpProbingMuOracle.hpp IpProbingMuOracle.cpp \
	IpRestoFilterConvCheck.cpp IpRestoFilterConvCheck.hpp \
	IpRestoIpoptNLP.hpp IpRestoIpoptNLP.cpp \
	IpRestoIterateInitializer.hpp IpRestoIterateInitializer.cpp \
	IpRestoIterationOutput.hpp IpRestoIterationOutput.cpp \
	IpRestoMinC_1Nrm.cpp IpRestoMinC_1Nrm.hpp \
	IpRestoPhase.hpp

CXXFLAGS += -I$(srcdir)/../Common
CXXFLAGS += -I$(srcdir)/../Interfaces
CXXFLAGS += -I$(srcdir)/../PDSystemImpl/PDFullSpace
CXXFLAGS += -I$(srcdir)/../LinAlgImpl/Serial

# Astyle stuff

ASTYLE_SOURCES = \
	IpAlgBuilder.hppbak IpAlgBuilder.cppbak \
	IpAlgStrategy.hppbak \
	IpConvCheck.hppbak \
	IpDefaultIterateInitializer.hppbak IpDefaultIterateInitializer.cppbak \
	IpFilter.cppbak IpFilter.hppbak \
	IpFilterLineSearch.cppbak IpFilterLineSearch.hppbak \
	IpIpoptAlg.cppbak IpIpoptAlg.hppbak \
	IpIpoptCalculatedQuantities.cppbak IpIpoptCalculatedQuantities.hppbak \
	IpIpoptData.cppbak IpIpoptData.hppbak \
	IpIpoptNLP.hppbak \
	IpIterateInitializer.hppbak \
	IpLineSearch.hppbak \
	IpLoqoMuOracle.hppbak IpLoqoMuOracle.cppbak \
	IpMonotoneMuUpdate.hppbak IpMonotoneMuUpdate.cppbak \
	IpMuOracle.hppbak \
	IpMuUpdate.hppbak \
	IpNonmonotoneMuUpdate.hppbak IpNonmonotoneMuUpdate.cppbak \
	IpOptErrorConvCheck.cppbak IpOptErrorConvCheck.hppbak \
	IpOrigIpoptNLP.cppbak IpOrigIpoptNLP.hppbak \
	IpIterationOutput.hppbak \
	IpOrigIterationOutput.hppbak IpOrigIterationOutput.cppbak \
	IpProbingMuOracle.hppbak IpProbingMuOracle.cppbak \
	IpRestoFilterConvCheck.cppbak IpRestoFilterConvCheck.hppbak \
	IpRestoIpoptNLP.hppbak IpRestoIpoptNLP.cppbak \
	IpRestoIterateInitializer.hppbak IpRestoIterateInitializer.cppbak \
	IpRestoIterationOutput.hppbak IpRestoIterationOutput.cppbak \
	IpRestoMinC_1Nrm.cppbak IpRestoMinC_1Nrm.hppbak \
	IpRestoPhase.hppbak

ASTYLE = @ASTYLE@
ASTYLEFLAGS = @ASTYLEFLAGS@

CLEANFILES = $(ASTYLE_SOURCES)

SUFFIXES = .cppbak .hppbak

astyle: $(ASTYLE_SOURCES)

.hpp.hppbak:
	mv $< $@
	$(ASTYLE) $(ASTYLEFLAGS) < $@ > $<
	touch $@

.cpp.cppbak:
	mv $< $@
	$(ASTYLE) $(ASTYLEFLAGS) < $@ > $<
	touch $@

# Copyright (C) 2004, 2008 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.


# Authors:  Carl Laird, Andreas Waechter, Hans Pirnay     IBM    2009-12-19

AUTOMAKE_OPTIONS = foreign

includesipoptdir = $(includedir)/coin

lib_LTLIBRARIES = libsipopt.la

libsipopt_la_SOURCES = \
	SensBacksolver.hpp \
	SensAlgorithm.cpp \
	SensAlgorithm.hpp \
	SensRegOp.cpp \
	SensRegOp.hpp \
	SensDenseGenSchurDriver.cpp \
	SensDenseGenSchurDriver.hpp \
	SensIndexPCalculator.cpp \
	SensIndexPCalculator.hpp \
	SensIndexSchurData.cpp \
	SensIndexSchurData.hpp \
	SensMeasurement.hpp \
	SensMetadataMeasurement.cpp \
	SensMetadataMeasurement.hpp \
	SensApplication.cpp \
	SensApplication.hpp \
	SensUtils.cpp \
	SensUtils.hpp \
	SensPCalculator.hpp \
	SensReducedHessianCalculator.cpp \
	SensReducedHessianCalculator.hpp \
	SensBuilder.cpp \
	SensBuilder.hpp \
	SensSchurData.hpp \
	SensSchurDriver.hpp \
	SensStepCalc.hpp \
	SensSimpleBacksolver.cpp \
	SensSimpleBacksolver.hpp \
	SensStdStepCalc.cpp \
	SensStdStepCalc.hpp \
	SensSuffixHandler.hpp

libsipopt_la_LDFLAGS = $(LT_LDFLAGS)

addlibsdir = $(datadir)/coin/doc/Ipopt

CXXFLAGS += -I`$(CYGPATH_W) $(srcdir)/../../../src/Common` \
-I`$(CYGPATH_W) $(srcdir)/../../../src/Interfaces` \
-I`$(CYGPATH_W) $(srcdir)/../../../src/LinAlg` \
-I`$(CYGPATH_W) $(srcdir)/../../../src/Algorithm/LinearSolvers` \
-I`$(CYGPATH_W) $(srcdir)/../../../src/LinAlg/TMatrices` \
-I`$(CYGPATH_W) $(srcdir)/../../../src/Algorithm` \
-I`$(CYGPATH_W) $(srcdir)/../../../src/Algorithm/Inexact` \
-I`$(CYGPATH_W) $(srcdir)/../../../src/contrib/CGPenalty`

CFLAGS += -I`$(CYGPATH_W) $(srcdir)/../Common`

# This line is necessary to allow VPATH compilation
DEFAULT_INCLUDES = -I. -I`$(CYGPATH_W) $(srcdir)` -I$(top_builddir)/src/Common

# Astyle stuff

ASTYLE_FILES = \
	SensBacksolver.hppbak \
	SensAlgorithm.cppbak \
	SensAlgorithm.hppbak \
	SensRegOp.cppbak \
	SensRegOp.hppbak \
	SensDenseGenSchurDriver.cppbak \
	SensDenseGenSchurDriver.hppbak \
	SensIndexPCalculator.cppbak \
	SensIndexPCalculator.hppbak \
	SensIndexSchurData.cppbak \
	SensIndexSchurData.hppbak \
	SensMeasurement.hppbak \
	SensMetadataMeasurement.cppbak \
	SensMetadataMeasurement.hppbak \
	SensApplication.cppbak \
	SensApplication.hppbak \
	SensUtils.cppbak \
	SensUtils.hppbak \
	SensPCalculator.hppbak \
	SensReducedHessianCalculator.cppbak \
	SensReducedHessianCalculator.hppbak \
	SensBuilder.cppbak \
	SensBuilder.hppbak \
	SensSchurData.hppbak \
	SensSchurDriver.hppbak \
	SensStepCalc.hppbak \
	SensSimpleBacksolver.cppbak \
	SensSimpleBacksolver.hppbak \
	SensStdStepCalc.cppbak \
	SensStdStepCalc.hppbak \
	SensSuffixHandler.hppbak

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


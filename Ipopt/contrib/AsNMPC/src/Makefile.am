# Copyright (C) 2004, 2008 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.


# Authors:  Carl Laird, Andreas Waechter, Hans Pirnay     IBM    2009-12-19

AUTOMAKE_OPTIONS = foreign

includeasnmpcdir = $(includedir)/coin

coinlibdir = $(libdir)/coin
coinlib_LTLIBRARIES = libasnmpc.la

libasnmpc_la_SOURCES = \
	AsAsBacksolver.hpp \
	AsAsNmpController.cpp \
	AsAsNmpController.hpp \
	AsAsNMPCRegOp.cpp \
	AsAsNMPCRegOp.hpp \
	AsDenseGenSchurDriver.cpp \
	AsDenseGenSchurDriver.hpp \
	AsIFTSchurDriver.cpp \
	AsIFTSchurDriver.hpp \
	AsIndexPCalculator.cpp \
	AsIndexPCalculator.hpp \
	AsIndexSchurData.cpp \
	AsIndexSchurData.hpp \
	AsMeasurement.hpp \
	AsMetadataMeasurement.cpp \
	AsMetadataMeasurement.hpp \
	AsNMPCApplication.cpp \
	AsNMPCApplication.hpp \
	AsNmpcUtils.cpp \
	AsNmpcUtils.hpp \
	AsPCalculator.hpp \
	AsReducedHessianCalculator.cpp \
	AsReducedHessianCalculator.hpp \
	AsSchurBuilder.cpp \
	AsSchurBuilder.hpp \
	AsSchurData.hpp \
	AsSchurDriver.hpp \
	AsSensStepCalc.hpp \
	AsSimpleBacksolver.cpp \
	AsSimpleBacksolver.hpp \
	AsStdPCalculator.cpp \
	AsStdPCalculator.hpp \
	AsStdSchurData.cpp \
	AsStdSchurData.hpp \
	AsStdStepCalc.cpp \
	AsStdStepCalc.hpp \
	AsSuffixHandler.hpp

#libasnmpc_la_LIBADD = 

#libasnmpc_la_DEPENDENCIES = $(IPALLLIBS)

libasnmpc_la_LDFLAGS = $(LT_LDFLAGS)

addlibsdir = $(datadir)/coin/doc/Ipopt
#addlibs_DATA = asnmpc_addlibs_cpp.txt asnmpc_addlibs_f.txt

CXXFLAGS += -I`$(CYGPATH_W) $(srcdir)/../../../src/Common` \
-I`$(CYGPATH_W) $(srcdir)/../../../src/Interfaces` \
-I`$(CYGPATH_W) $(srcdir)/../../../src/LinAlg` \
-I`$(CYGPATH_W) $(srcdir)/../../../src/Algorithm/LinearSolvers` \
-I`$(CYGPATH_W) $(srcdir)/../../../src/LinAlg/TMatrices` \
-I`$(CYGPATH_W) $(srcdir)/../../../src/Algorithm` \
-I`$(CYGPATH_W) $(srcdir)/../../../src/Algorithm/Inexact` \
-I`$(CYGPATH_W) $(srcdir)/../../../src/contrib/CGPenalty`

CFLAGS += -I`$(CYGPATH_W) $(srcdir)/../Common`

# This line is necessary to allow VPATH compilation with MS compilers
# on Cygwin
DEFAULT_INCLUDES = -I. -I`$(CYGPATH_W) $(srcdir)` -I$(top_builddir)/inc

# Astyle stuff

ASTYLE_FILES = \
	AsAsBacksolver.hppbak \
	AsAsNmpController.cppbak \
	AsAsNmpController.hppbak \
	AsAsNMPCRegOp.cppbak \
	AsAsNMPCRegOp.hppbak \
	AsDenseGenSchurDriver.cppbak \
	AsDenseGenSchurDriver.hppbak \
	AsIFTSchurDriver.cppbak \
	AsIFTSchurDriver.hppbak \
	AsIndexPCalculator.cppbak \
	AsIndexPCalculator.hppbak \
	AsIndexSchurData.cppbak \
	AsIndexSchurData.hppbak \
	AsMeasurement.hppbak \
	AsMetadataMeasurement.cppbak \
	AsMetadataMeasurement.hppbak \
	AsNMPCAppbaklication.cppbak \
	AsNMPCAppbaklication.hppbak \
	AsNmpcUtils.cppbak \
	AsNmpcUtils.hppbak \
	AsPCalculator.hppbak \
	AsReducedHessianCalculator.cppbak \
	AsReducedHessianCalculator.hppbak \
	AsSchurBuilder.cppbak \
	AsSchurBuilder.hppbak \
	AsSchurData.hppbak \
	AsSchurDriver.hppbak \
	AsSensStepCalc.hppbak \
	AsSimpleBacksolver.cppbak \
	AsSimpleBacksolver.hppbak \
	AsStdPCalculator.cppbak \
	AsStdPCalculator.hppbak \
	AsStdSchurData.cppbak \
	AsStdSchurData.hppbak \
	AsStdStepCalc.cppbak \
	AsStdStepCalc.hppbak \
	AsSuffixHandler.hppbak

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


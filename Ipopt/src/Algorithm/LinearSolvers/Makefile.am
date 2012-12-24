# Copyright (C) 2004, 2011 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.

## $Id$

# Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

AUTOMAKE_OPTIONS = foreign

noinst_LTLIBRARIES = liblinsolvers.la

liblinsolvers_la_SOURCES = \
	IpGenKKTSolverInterface.hpp \
	IpLinearSolversRegOp.cpp IpLinearSolversRegOp.hpp \
	IpSlackBasedTSymScalingMethod.cpp IpSlackBasedTSymScalingMethod.hpp \
	IpSparseSymLinearSolverInterface.hpp \
	IpSymLinearSolver.hpp \
	IpTDependencyDetector.hpp \
	IpTripletToCSRConverter.cpp IpTripletToCSRConverter.hpp \
	IpTSymDependencyDetector.cpp IpTSymDependencyDetector.hpp \
	IpTSymLinearSolver.cpp IpTSymLinearSolver.hpp \
	IpTSymScalingMethod.hpp

liblinsolvers_la_SOURCES += \
	IpMa27TSolverInterface.cpp IpMa27TSolverInterface.hpp \
	IpMa57TSolverInterface.cpp IpMa57TSolverInterface.hpp \
	IpMa86SolverInterface.cpp IpMa86SolverInterface.hpp \
	IpMa97SolverInterface.cpp IpMa97SolverInterface.hpp \
	IpMc19TSymScalingMethod.cpp IpMc19TSymScalingMethod.hpp \
	IpMa28TDependencyDetector.cpp IpMa28TDependencyDetector.hpp \
	IpMa77SolverInterface.cpp IpMa77SolverInterface.hpp

if BUILD_LINEARSOLVERLOADER
  liblinsolvers_la_SOURCES += \
	IpPardisoSolverInterface.cpp IpPardisoSolverInterface.hpp \
	IpMa28Partition.F
else
if HAVE_PARDISO
  liblinsolvers_la_SOURCES += \
	IpPardisoSolverInterface.cpp IpPardisoSolverInterface.hpp
endif
if HAVE_MA28
  liblinsolvers_la_SOURCES += IpMa28Partition.F
endif
endif # BUILD_LINEARSOLVERLOADER

if HAVE_WSMP
  liblinsolvers_la_SOURCES += \
	IpWsmpSolverInterface.cpp IpWsmpSolverInterface.hpp \
	IpIterativeWsmpSolverInterface.cpp IpIterativeWsmpSolverInterface.hpp
endif

if COIN_HAS_MUMPS
liblinsolvers_la_SOURCES += IpMumpsSolverInterface.cpp IpMumpsSolverInterface.hpp
endif

liblinsolvers_la_LDFLAGS = $(LT_LDFLAGS)

AM_CPPFLAGS = -I`$(CYGPATH_W) $(srcdir)/../../Common` \
	-I`$(CYGPATH_W) $(srcdir)/../../LinAlg` \
	-I`$(CYGPATH_W) $(srcdir)/../../LinAlg/TMatrices` \
	-I`$(CYGPATH_W) $(srcdir)/..` \
	-I`$(CYGPATH_W) $(srcdir)/../../Interfaces` \
	-I`$(CYGPATH_W) $(srcdir)/../../contrib/CGPenalty` \
	$(MUMPS_CFLAGS) $(HSL_CFLAGS)

# This line is necessary to allow VPATH compilation
DEFAULT_INCLUDES = -I. -I`$(CYGPATH_W) $(srcdir)` -I$(top_builddir)/src/Common

# Workaround a bug in the AIX C++ compiler
if AIX_FIX_COMPILER_BUG
IpTSymLinearSolver.lo: IpTSymLinearSolver.cpp
@am__fastdepCXX_TRUE@	if $(LTCXXCOMPILE) -qnooptimize -MT $@ -MD -MP -MF "$(DEPDIR)/$*.Tpo" -c -o $@ $<; \
@am__fastdepCXX_TRUE@	then mv -f "$(DEPDIR)/$*.Tpo" "$(DEPDIR)/$*.Plo"; else rm -f "$(DEPDIR)/$*.Tpo"; exit 1; fi
@AMDEP_TRUE@@am__fastdepCXX_FALSE@	source='$<' object='$@' libtool=yes @AMDEPBACKSLASH@
@AMDEP_TRUE@@am__fastdepCXX_FALSE@	DEPDIR=$(DEPDIR) $(CXXDEPMODE) $(depcomp) @AMDEPBACKSLASH@
@am__fastdepCXX_FALSE@	$(LTCXXCOMPILE) -qnooptimize -c -o $@ $<
endif

# Need to overwrite automake definition since otherwise libtool does not
# understand the fortran file (need to add --tag=F77)
AM_FCFLAGS = $(HSL_F90_INCFLAGS)
LTFCCOMPILE = $(LIBTOOL) --tag=F77 --mode=compile $(FC) $(AM_FCFLAGS) $(FCFLAGS)

# Astyle stuff

ASTYLE_FILES = \
	IpGenKKTSolverInterface.hppbak \
	IpIterativeWsmpSolverInterface.cppbak \
	IpIterativeWsmpSolverInterface.hppbak \
	IpLinearSolversRegOp.cppbak IpLinearSolversRegOp.hppbak \
	IpMa27TSolverInterface.cppbak IpMa27TSolverInterface.hppbak \
	IpMa28TDependencyDetector.cppbak IpMa28TDependencyDetector.hppbak \
	IpMa57TSolverInterface.cppbak IpMa57TSolverInterface.hppbak \
	IpMa86SolverInterface.cppbak IpMa86SolverInterface.hppbak \
	IpMa97SolverInterface.cppbak IpMa97SolverInterface.hppbak \
	IpMc19TSymScalingMethod.cppbak IpMc19TSymScalingMethod.hppbak \
	IpMumpsSolverInterface.cppbak IpMumpsSolverInterface.hppbak \
	IpSlackBasedTSymScalingMethod.cppbak IpSlackBasedTSymScalingMethod.hppbak \
	IpSparseSymLinearSolverInterface.hppbak \
	IpSymLinearSolver.hppbak \
	IpTDependencyDetector.hppbak \
	IpTripletToCSRConverter.cppbak IpTripletToCSRConverter.hppbak \
	IpTSymDependencyDetector.cppbak IpTSymDependencyDetector.hppbak \
	IpTSymLinearSolver.cppbak IpTSymLinearSolver.hppbak \
	IpTSymScalingMethod.hppbak \
	IpPardisoSolverInterface.cppbak IpPardisoSolverInterface.hppbak \
	IpWsmpSolverInterface.cppbak IpWsmpSolverInterface.hppbak


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

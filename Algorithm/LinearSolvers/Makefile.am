# Copyright (C) 2004, 2006 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Common Public License.

## $Id$

# Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

AUTOMAKE_OPTIONS = foreign

noinst_LIBRARIES = liblinsolvers.a

liblinsolvers_a_SOURCES = \
	IpLinearSolversRegOp.cpp IpLinearSolversRegOp.hpp \
	IpSparseSymLinearSolverInterface.hpp \
	IpSymLinearSolver.hpp \
	IpTripletToCSRConverter.cpp IpTripletToCSRConverter.hpp \
	IpTSymLinearSolver.cpp IpTSymLinearSolver.hpp \
	IpTSymScalingMethod.hpp

if HAVE_MA27
  liblinsolvers_a_SOURCES += \
	IpMa27TSolverInterface.cpp IpMa27TSolverInterface.hpp
endif

if HAVE_MA57
  liblinsolvers_a_SOURCES += \
	IpMa57TSolverInterface.cpp IpMa57TSolverInterface.hpp
endif

if HAVE_MC19
  liblinsolvers_a_SOURCES += \
	IpMc19TSymScalingMethod.cpp IpMc19TSymScalingMethod.hpp
endif

if HAVE_PARDISO
  liblinsolvers_a_SOURCES += \
	IpPardisoSolverInterface.cpp IpPardisoSolverInterface.hpp
endif

if HAVE_TAUCS
  liblinsolvers_a_SOURCES += \
	IpTAUCSSolverInterface.cpp IpTAUCSSolverInterface.hpp
endif

if HAVE_WSMP
  liblinsolvers_a_SOURCES += \
	IpWsmpSolverInterface.cpp IpWsmpSolverInterface.hpp
endif


CXXFLAGS += -I`$(CYGPATH_W) $(srcdir)/../../Common` \
	-I`$(CYGPATH_W) $(srcdir)/../../LinAlg` \
	-I`$(CYGPATH_W) $(srcdir)/../../LinAlg/TMatrices` \
	-I`$(CYGPATH_W) $(srcdir)/..` \
	-I`$(CYGPATH_W) $(srcdir)/../../Interfaces`

if HAVE_MUMPS
liblinsolvers_a_SOURCES += \
	IpMumpsSolverInterface.cpp IpMumpsSolverInterface.hpp
CXXFLAGS += $(MUMPS_INCFLAGS)
endif

# Astyle stuff

ASTYLE_FILES = \
	IpLinearSolversRegOp.cppbak IpLinearSolversRegOp.hppbak \
	IpMa27TSolverInterface.cppbak IpMa27TSolverInterface.hppbak \
	IpMa57TSolverInterface.cppbak IpMa57TSolverInterface.hppbak \
	IpMc19TSymScalingMethod.cppbak IpMc19TSymScalingMethod.hppbak \
	IpMumpsSolverInterface.cppbak IpMumpsSolverInterface.hppbak \
	IpSparseSymLinearSolverInterface.hppbak \
	IpSymLinearSolver.hppbak \
	IpTripletToCSRConverter.cppbak IpTripletToCSRConverter.hppbak \
	IpTSymLinearSolver.cppbak IpTSymLinearSolver.hppbak \
	IpTSymScalingMethod.hppbak \
	IpPardisoSolverInterface.cppbak IpPardisoSolverInterface.hppbak \
	IpTAUCSSolverInterface.cppbak IpTAUCSSolverInterface.hppbak \
	IpWsmpSolverInterface.cppbak IpWsmpSolverInterface.hppbak


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

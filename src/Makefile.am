# Copyright (C) 2021 COIN-OR Foundation
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.
#
# Authors: Stefan Vigerske

includeipoptdir = $(includedir)/coin-or
includeipopt_HEADERS = \
  Common/IpCachedResults.hpp \
  Common/IpDebug.hpp \
  Common/IpException.hpp \
  Common/IpJournalist.hpp \
  Common/IpLibraryLoader.hpp \
  Common/IpObserver.hpp \
  Common/IpOptionsList.hpp \
  Common/IpReferenced.hpp \
  Common/IpRegOptions.hpp \
  Common/IpSmartPtr.hpp \
  Common/IpTaggedObject.hpp \
  Common/IpTimedTask.hpp \
  Common/IpTypes.hpp \
  Common/IpTypes.h \
  Common/IpUtils.hpp \
  LinAlg/IpBlas.hpp \
  LinAlg/IpCompoundMatrix.hpp \
  LinAlg/IpCompoundSymMatrix.hpp \
  LinAlg/IpCompoundVector.hpp \
  LinAlg/IpDenseVector.hpp \
  LinAlg/IpDiagMatrix.hpp \
  LinAlg/IpExpansionMatrix.hpp \
  LinAlg/IpIdentityMatrix.hpp \
  LinAlg/IpLapack.hpp \
  LinAlg/IpMatrix.hpp \
  LinAlg/IpScaledMatrix.hpp \
  LinAlg/IpSumSymMatrix.hpp \
  LinAlg/IpSymMatrix.hpp \
  LinAlg/IpSymScaledMatrix.hpp \
  LinAlg/IpVector.hpp \
  LinAlg/IpZeroSymMatrix.hpp \
  LinAlg/TMatrices/IpGenTMatrix.hpp \
  LinAlg/TMatrices/IpSymTMatrix.hpp \
  LinAlg/TMatrices/IpTripletHelper.hpp \
  Algorithm/IpAlgBuilder.hpp \
  Algorithm/IpAlgStrategy.hpp \
  Algorithm/IpAugSystemSolver.hpp \
  Algorithm/IpConvCheck.hpp \
  Algorithm/IpEqMultCalculator.hpp \
  Algorithm/IpHessianUpdater.hpp \
  Algorithm/IpIpoptAlg.hpp \
  Algorithm/IpIpoptCalculatedQuantities.hpp \
  Algorithm/IpIpoptData.hpp \
  Algorithm/IpIpoptNLP.hpp \
  Algorithm/IpIterateInitializer.hpp \
  Algorithm/IpIteratesVector.hpp \
  Algorithm/IpIterationOutput.hpp \
  Algorithm/IpOrigIpoptNLP.hpp \
  Algorithm/IpLineSearch.hpp \
  Algorithm/IpMuUpdate.hpp \
  Algorithm/IpNLPScaling.hpp \
  Algorithm/IpPDSystemSolver.hpp \
  Algorithm/IpSearchDirCalculator.hpp \
  Algorithm/IpStdAugSystemSolver.hpp \
  Algorithm/IpTimingStatistics.hpp \
  Algorithm/LinearSolvers/IpLinearSolvers.h \
  Algorithm/LinearSolvers/IpSlackBasedTSymScalingMethod.hpp \
  Algorithm/LinearSolvers/IpSparseSymLinearSolverInterface.hpp \
  Algorithm/LinearSolvers/IpSymLinearSolver.hpp \
  Algorithm/LinearSolvers/IpTripletToCSRConverter.hpp \
  Algorithm/LinearSolvers/IpTSymLinearSolver.hpp \
  Algorithm/LinearSolvers/IpTSymScalingMethod.hpp \
  Interfaces/IpAlgTypes.hpp \
  Interfaces/IpIpoptApplication.hpp \
  Interfaces/IpNLP.hpp \
  Interfaces/IpReturnCodes.h \
  Interfaces/IpReturnCodes.hpp \
  Interfaces/IpReturnCodes_inc.h \
  Interfaces/IpReturnCodes.inc \
  Interfaces/IpSolveStatistics.hpp \
  Interfaces/IpStdCInterface.h \
  Interfaces/IpTNLP.hpp \
  Interfaces/IpTNLPAdapter.hpp \
  Interfaces/IpTNLPReducer.hpp

lib_LTLIBRARIES = libipopt.la
libipopt_la_SOURCES = \
  Common/IpDebug.cpp \
  Common/IpJournalist.cpp \
  Common/IpObserver.cpp \
  Common/IpOptionsList.cpp \
  Common/IpRegOptions.cpp \
  Common/IpTaggedObject.cpp \
  Common/IpUtils.cpp \
  Common/IpLibraryLoader.cpp \
  LinAlg/IpBlas.cpp \
  LinAlg/IpCompoundMatrix.cpp \
  LinAlg/IpCompoundSymMatrix.cpp \
  LinAlg/IpCompoundVector.cpp \
  LinAlg/IpDenseGenMatrix.cpp \
  LinAlg/IpDenseSymMatrix.cpp \
  LinAlg/IpDenseVector.cpp \
  LinAlg/IpDiagMatrix.cpp \
  LinAlg/IpExpandedMultiVectorMatrix.cpp \
  LinAlg/IpExpansionMatrix.cpp \
  LinAlg/IpIdentityMatrix.cpp \
  LinAlg/IpLapack.cpp \
  LinAlg/IpLowRankUpdateSymMatrix.cpp \
  LinAlg/IpMatrix.cpp \
  LinAlg/IpMultiVectorMatrix.cpp \
  LinAlg/IpScaledMatrix.cpp \
  LinAlg/IpSumMatrix.cpp \
  LinAlg/IpSumSymMatrix.cpp \
  LinAlg/IpSymScaledMatrix.cpp \
  LinAlg/IpTransposeMatrix.cpp \
  LinAlg/IpVector.cpp \
  LinAlg/IpZeroMatrix.cpp \
  LinAlg/IpZeroSymMatrix.cpp \
  LinAlg/TMatrices/IpGenTMatrix.cpp \
  LinAlg/TMatrices/IpSymTMatrix.cpp \
  LinAlg/TMatrices/IpTripletHelper.cpp \
  Algorithm/IpAdaptiveMuUpdate.cpp \
  Algorithm/IpAlgBuilder.cpp \
  Algorithm/IpAlgorithmRegOp.cpp \
  Algorithm/IpAugRestoSystemSolver.cpp \
  Algorithm/IpBacktrackingLineSearch.cpp \
  Algorithm/IpDefaultIterateInitializer.cpp \
  Algorithm/IpEquilibrationScaling.cpp \
  Algorithm/IpExactHessianUpdater.cpp \
  Algorithm/IpFilter.cpp \
  Algorithm/IpFilterLSAcceptor.cpp \
  Algorithm/IpGenAugSystemSolver.cpp \
  Algorithm/IpGradientScaling.cpp \
  Algorithm/IpIpoptAlg.cpp \
  Algorithm/IpIpoptCalculatedQuantities.cpp \
  Algorithm/IpIpoptData.cpp \
  Algorithm/IpIteratesVector.cpp \
  Algorithm/IpLeastSquareMults.cpp \
  Algorithm/IpLimMemQuasiNewtonUpdater.cpp \
  Algorithm/IpLoqoMuOracle.cpp \
  Algorithm/IpLowRankAugSystemSolver.cpp \
  Algorithm/IpLowRankSSAugSystemSolver.cpp \
  Algorithm/IpMonotoneMuUpdate.cpp \
  Algorithm/IpNLPBoundsRemover.cpp \
  Algorithm/IpNLPScaling.cpp \
  Algorithm/IpOptErrorConvCheck.cpp \
  Algorithm/IpOrigIpoptNLP.cpp \
  Algorithm/IpOrigIterationOutput.cpp \
  Algorithm/IpPDFullSpaceSolver.cpp \
  Algorithm/IpPDPerturbationHandler.cpp \
  Algorithm/IpPDSearchDirCalc.cpp \
  Algorithm/IpPenaltyLSAcceptor.cpp \
  Algorithm/IpProbingMuOracle.cpp \
  Algorithm/IpQualityFunctionMuOracle.cpp \
  Algorithm/IpRestoConvCheck.cpp \
  Algorithm/IpRestoFilterConvCheck.cpp \
  Algorithm/IpRestoIpoptNLP.cpp \
  Algorithm/IpRestoIterateInitializer.cpp \
  Algorithm/IpRestoIterationOutput.cpp \
  Algorithm/IpRestoMinC_1Nrm.cpp \
  Algorithm/IpRestoPenaltyConvCheck.cpp \
  Algorithm/IpRestoRestoPhase.cpp \
  Algorithm/IpStdAugSystemSolver.cpp \
  Algorithm/IpTimingStatistics.cpp \
  Algorithm/IpUserScaling.cpp \
  Algorithm/IpWarmStartIterateInitializer.cpp \
  Algorithm/LinearSolvers/IpLinearSolversRegOp.cpp \
  Algorithm/LinearSolvers/IpLinearSolvers.c \
  Algorithm/LinearSolvers/IpSlackBasedTSymScalingMethod.cpp \
  Algorithm/LinearSolvers/IpTripletToCSRConverter.cpp \
  Algorithm/LinearSolvers/IpTSymDependencyDetector.cpp \
  Algorithm/LinearSolvers/IpTSymLinearSolver.cpp \
  contrib/CGPenalty/IpCGPenaltyCq.cpp \
  contrib/CGPenalty/IpCGPenaltyData.cpp \
  contrib/CGPenalty/IpCGPenaltyLSAcceptor.cpp \
  contrib/CGPenalty/IpCGPenaltyRegOp.cpp \
  contrib/CGPenalty/IpCGPerturbationHandler.cpp \
  contrib/CGPenalty/IpCGSearchDirCalc.cpp \
  contrib/CGPenalty/IpPiecewisePenalty.cpp \
  Interfaces/IpInterfacesRegOp.cpp \
  Interfaces/IpIpoptApplication.cpp \
  Interfaces/IpSolveStatistics.cpp \
  Interfaces/IpStdCInterface.cpp \
  Interfaces/IpStdInterfaceTNLP.cpp \
  Interfaces/IpStdFInterface.c \
  Interfaces/IpTNLP.cpp \
  Interfaces/IpTNLPAdapter.cpp \
  Interfaces/IpTNLPReducer.cpp

if HAVE_PARDISO_MKL
  libipopt_la_SOURCES += Algorithm/LinearSolvers/IpPardisoMKLSolverInterface.cpp
endif

if !IPOPT_INT64
libipopt_la_SOURCES += \
  Algorithm/LinearSolvers/IpMc19TSymScalingMethod.cpp \
  Algorithm/LinearSolvers/IpMa27TSolverInterface.cpp \
  Algorithm/LinearSolvers/IpMa57TSolverInterface.cpp \
  Algorithm/LinearSolvers/IpMa77SolverInterface.cpp \
  Algorithm/LinearSolvers/IpMa86SolverInterface.cpp \
  Algorithm/LinearSolvers/IpMa97SolverInterface.cpp \
  Algorithm/LinearSolvers/IpPardisoSolverInterface.cpp

if HAVE_MA28
  libipopt_la_SOURCES += Algorithm/LinearSolvers/IpMa28TDependencyDetector.cpp
if IPOPT_SINGLE
  libipopt_la_SOURCES += Algorithm/LinearSolvers/IpMa28sPartition.F
else
  libipopt_la_SOURCES += Algorithm/LinearSolvers/IpMa28Partition.F
endif  # IPOPT_SINGLE
endif  # HAVE_MA28

endif  # IPOPT_INT64

if HAVE_WSMP
  libipopt_la_SOURCES += Algorithm/LinearSolvers/IpWsmpSolverInterface.cpp Algorithm/LinearSolvers/IpIterativeWsmpSolverInterface.cpp
endif

if COIN_HAS_MUMPS
  libipopt_la_SOURCES += Algorithm/LinearSolvers/IpMumpsSolverInterface.cpp
endif

if COIN_HAS_SPRAL
  libipopt_la_SOURCES += Algorithm/LinearSolvers/IpSpralSolverInterface.cpp
endif

if BUILD_INEXACT
libipopt_la_SOURCES += \
  Algorithm/Inexact/IpInexactAlgBuilder.cpp \
  Algorithm/Inexact/IpInexactCq.cpp \
  Algorithm/Inexact/IpInexactData.cpp \
  Algorithm/Inexact/IpInexactDoglegNormal.cpp \
  Algorithm/Inexact/IpInexactLSAcceptor.cpp \
  Algorithm/Inexact/IpInexactNewtonNormal.cpp \
  Algorithm/Inexact/IpInexactNormalTerminationTester.cpp \
  Algorithm/Inexact/IpInexactPDSolver.cpp \
  Algorithm/Inexact/IpInexactPDTerminationTester.cpp \
  Algorithm/Inexact/IpInexactRegOp.cpp \
  Algorithm/Inexact/IpInexactSearchDirCalc.cpp \
  Algorithm/Inexact/IpInexactTSymScalingMethod.cpp \
  Algorithm/Inexact/IpIterativePardisoSolverInterface.cpp \
  Algorithm/Inexact/IpIterativeSolverTerminationTester.cpp
endif

if BUILD_JAVA
  libipopt_la_SOURCES += Interfaces/IpStdJInterface.cpp org_coinor_Ipopt.h
endif

AM_CPPFLAGS = \
  -I$(srcdir)/Common \
  -I$(srcdir)/LinAlg \
  -I$(srcdir)/LinAlg/TMatrices \
  -I$(srcdir)/Algorithm \
  -I$(srcdir)/Algorithm/LinearSolvers \
  -I$(srcdir)/Algorithm/Inexact \
  -I$(srcdir)/Interfaces \
  -I$(srcdir)/contrib/CGPenalty \
  $(IPOPTLIB_CFLAGS)

AM_LDFLAGS = $(LT_LDFLAGS)

libipopt_la_LIBADD = $(IPOPTLIB_LFLAGS)

if BUILD_JAVA
BUILT_SOURCES = org_coinor_Ipopt.h

org.coinor.ipopt.jar : org/coinor/Ipopt.class
	$(JAR) cf $@ $<

if IPOPT_SINGLE
# replace double by float in Ipopt.java; filename needs to be Ipopt.java, so put into subdir
Ipopt.java : $(srcdir)/Interfaces/Ipopt.java
	sed -e s/double/float/g $< > $@

org/coinor/Ipopt.class org_coinor_Ipopt.h : Ipopt.java
	$(JAVAC) $< -d $(PWD) -h $(PWD)
else
# had  -target 1.8 -source 1.8, but that produces a warning if Java 8 RT is not setup
org/coinor/Ipopt.class org_coinor_Ipopt.h : Interfaces/Ipopt.java
	$(JAVAC) $< -d $(PWD) -h $(PWD)
endif

all-local : org.coinor.ipopt.jar

CLEANFILES = org.coinor.ipopt.jar org/coinor/Ipopt.class org_coinor_Ipopt.h
endif

install-exec-local:
	$(install_sh_DATA) Common/config_ipopt.h $(DESTDIR)$(includeipoptdir)/IpoptConfig.h
if BUILD_JAVA
	$(MKDIR_P) "$(DESTDIR)$(datadir)/java"
	$(LIBTOOL) $(AM_LIBTOOLFLAGS) $(LIBTOOLFLAGS) --mode=install $(INSTALL) $(INSTALL_STRIP_FLAG) org.coinor.ipopt.jar "$(DESTDIR)$(datadir)/java"
endif

uninstall-local:
	rm -f $(DESTDIR)$(includeipoptdir)/IpoptConfig.h
if BUILD_JAVA
	$(LIBTOOL) $(AM_LIBTOOLFLAGS) $(LIBTOOLFLAGS) --mode=uninstall rm -f "$(DESTDIR)$(datadir)/java/org.coinor.ipopt.jar"
endif

pkgconfiglibdir = $(libdir)/pkgconfig
pkgconfiglib_DATA = ipopt.pc

if COIN_HAS_ASL
  SUBDIRS = . Apps/AmplSolver
endif

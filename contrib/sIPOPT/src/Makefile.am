# Copyright (C) 2004, 2008 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.
#
# Authors:  Carl Laird, Andreas Waechter, Hans Pirnay     IBM    2009-12-19

includesipoptdir = $(includedir)/coin-or
includesipopt_HEADERS = \
  SensAlgorithm.hpp \
  SensApplication.hpp \
  SensBacksolver.hpp \
  SensMeasurement.hpp \
  SensPCalculator.hpp \
  SensRegOp.hpp \
  SensSchurData.hpp \
  SensSchurDriver.hpp \
  SensSimpleBacksolver.hpp \
  SensStepCalc.hpp \
  SensUtils.hpp

lib_LTLIBRARIES = libsipopt.la

libsipopt_la_SOURCES = \
  SensAlgorithm.cpp \
  SensRegOp.cpp \
  SensDenseGenSchurDriver.cpp \
  SensIndexPCalculator.cpp \
  SensIndexSchurData.cpp \
  SensMetadataMeasurement.cpp \
  SensApplication.cpp \
  SensUtils.cpp \
  SensReducedHessianCalculator.cpp \
  SensBuilder.cpp \
  SensSimpleBacksolver.cpp \
  SensStdStepCalc.cpp

libsipopt_la_LIBADD = ../../../src/libipopt.la

AM_LDFLAGS = $(LT_LDFLAGS)

AM_CPPFLAGS = -DSIPOPTLIB_BUILD \
  -I$(srcdir)/../../../src/Common \
  -I$(srcdir)/../../../src/Interfaces \
  -I$(srcdir)/../../../src/LinAlg \
  -I$(srcdir)/../../../src/Algorithm/LinearSolvers \
  -I$(srcdir)/../../../src/LinAlg/TMatrices \
  -I$(srcdir)/../../../src/Algorithm \
  -I$(srcdir)/../../../src/Algorithm/Inexact \
  -I$(srcdir)/../../../src/contrib/CGPenalty

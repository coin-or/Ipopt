# Copyright (C) 2019 COIN-OR
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.

bin_PROGRAMS = ipopt_sens

ipopt_sens_SOURCES = ampl_sipopt.cpp SensAmplTNLP.cpp

ipopt_sens_LDADD = ../src/libsipopt.la \
  ../../../src/Apps/AmplSolver/libipoptamplinterface.la \
  $(SIPOPTAMPLINTERFACELIB_LFLAGS)

AM_LDFLAGS = $(LT_LDFLAGS)

AM_CPPFLAGS = \
  -I$(srcdir)/../src \
  -I$(srcdir)/../../../src/Common \
  -I$(srcdir)/../../../src/Interfaces \
  -I$(srcdir)/../../../src/LinAlg \
  -I$(srcdir)/../../../src/Algorithm \
  -I$(srcdir)/../../../src/Apps/AmplSolver \
  $(SIPOPTAMPLINTERFACELIB_CFLAGS)

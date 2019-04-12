# Copyright (C) 2003, International Business Machines and others
# All Rights Reserved.
# This file is distributed under the Common Public License.

## $Id: Makefile.am 632 2004-07-20 23:52:15Z andreasw $

iplibdir = $(top_builddir)/lib

AUTOMAKE_OPTIONS = foreign

noinst_LIBRARIES = libiphsl.a

nodist_libiphsl_a_SOURCES =

if BUILDFD05
  nodist_libiphsl_a_SOURCES += fd05ad.f
endif

if BUILDMA27
  nodist_libiphsl_a_SOURCES += ma27ad.f
endif

if BUILDMA28
  nodist_libiphsl_a_SOURCES += ma28ad.f ma30ad.f mc22ad.f mc23ad.f mc24ad.f
endif

if BUILDMA47
  nodist_libiphsl_a_SOURCES += ma47ad.f
endif

if BUILDMA48
  nodist_libiphsl_a_SOURCES += ma48ad.f ma50ad.f ma51ad.f mc41ad.f
endif

if BUILDMA57
  nodist_libiphsl_a_SOURCES += ma57ad.f mc34ad.f mc47ad.f mc59ad.f mc71ad.f
endif

if BUILDMC13
  nodist_libiphsl_a_SOURCES += mc13dd.f
endif

if BUILDMC19
  nodist_libiphsl_a_SOURCES += mc19ad.f
endif

if BUILDMC20
  nodist_libiphsl_a_SOURCES += mc20ad.f
endif

if BUILDMC21
  nodist_libiphsl_a_SOURCES += mc21ad.f
endif

if BUILDMC29
  nodist_libiphsl_a_SOURCES += mc29ad.f
endif

if BUILDMC30
  nodist_libiphsl_a_SOURCES += mc30ad.f
endif

if BUILDMC35
  nodist_libiphsl_a_SOURCES += mc35ad.f
endif

if BUILDMC39
  nodist_libiphsl_a_SOURCES += mc39ad.f
endif

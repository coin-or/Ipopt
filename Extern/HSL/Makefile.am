# Copyright (C) 2003, 2006 International Business Machines and others
# All Rights Reserved.
# This file is distributed under the Common Public License.

## $Id$

AUTOMAKE_OPTIONS = foreign

EXTRA_DIST = INSTALL.HSL

iplibdir = $(top_builddir)/lib

noinst_LTLIBRARIES = libiphsl.la

nodist_libiphsl_la_SOURCES =

if BUILDMA27
  nodist_libiphsl_la_SOURCES += ma27ad.f
endif

if BUILDMA57
  nodist_libiphsl_la_SOURCES += ma57ad.f
endif

if BUILDMC19
  nodist_libiphsl_la_SOURCES += mc19ad.f
endif

libiphsl_la_LDFLAGS = -no-undefined


# Copyright (C) 2003, 2006 International Business Machines and others
# All Rights Reserved.
# This file is distributed under the Common Public License.

## $Id$

AUTOMAKE_OPTIONS = foreign

EXTRA_DIST = INSTALL.HSL

iplibdir = $(top_builddir)/lib

noinst_LIBRARIES = libiphsl.a

nodist_libiphsl_a_SOURCES =

if BUILDMA27
  nodist_libiphsl_a_SOURCES += ma27ad.f
endif

if BUILDMA57
  nodist_libiphsl_a_SOURCES += ma57ad.f
endif

if BUILDMC19
  nodist_libiphsl_a_SOURCES += mc19ad.f
endif

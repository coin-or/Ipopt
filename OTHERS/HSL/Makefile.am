# Copyright (C) 2003, International Business Machines and others
# All Rights Reserved.
# This file is distributed under the Common Public License.

## $Id$

iplibdir = $(top_builddir)/lib

AUTOMAKE_OPTIONS = foreign

noinst_LIBRARIES = libiphsl.a

nodist_libiphsl_a_SOURCES =

if BUILDMA27
  nodist_libiphsl_a_SOURCES += ma27ad.f
endif

if BUILDMC19
  nodist_libiphsl_a_SOURCES += mc19ad.f
endif

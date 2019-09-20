# Copyright (C) 2008, 2011    GAMS Development
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.
#
# Author: Stefan Vigerske

AUTOMAKE_OPTIONS = foreign

EXTRA_DIST = 

noinst_LTLIBRARIES = libLinearSolverLoader.la

libLinearSolverLoader_la_SOURCES = \
  LibraryHandler.c \
  HSLLoader.c

includeipoptdir = $(includedir)/coin-or
includeipopt_HEADERS = HSLLoader.h

AM_CPPFLAGS = -DIPOPT_BUILD \
  -I$(srcdir)/../../Common \
  -I$(srcdir)/../../Algorithm/LinearSolvers \
  $(IPOPTLIB_CFLAGS)

if !HAVE_PARDISO
  libLinearSolverLoader_la_SOURCES += PardisoLoader.c
  includeipopt_HEADERS += PardisoLoader.h
endif

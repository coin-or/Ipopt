# Copyright (C) 2008     GAMS Development
# All Rights Reserved.
# This file is distributed under the Common Public License.
#
## $Id: Makefile.am 340 2008-02-14 16:02:10Z stefan $
#
# Author: Stefan Vigerske

AUTOMAKE_OPTIONS = foreign

EXTRA_DIST = 

noinst_LTLIBRARIES = libLinearSolverLoader.la

libLinearSolverLoader_la_SOURCES = LibraryHandler.c LibraryHandler.h

if BUILD_HSLLOADER
  libLinearSolverLoader_la_SOURCES += HSLLoader.c HSLLoader.h
endif

if BUILD_PARDISOLOADER
  libLinearSolverLoader_la_SOURCES += PardisoLoader.c PardisoLoader.h
endif

AM_CPPFLAGS = 

libLinearSolverLoader_la_LIBADD = 
libLinearSolverLoader_la_DEPENDENCIES =  

DEFAULT_INCLUDES = -I. -I`$(CYGPATH_W) $(srcdir)` -I$(top_builddir)/inc

DISTCLEANFILES = 

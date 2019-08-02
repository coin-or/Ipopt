# Copyright (C) 2008, 2011    GAMS Development
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.
#
# Author: Stefan Vigerske

AUTOMAKE_OPTIONS = foreign

EXTRA_DIST = 

noinst_LTLIBRARIES = libLinearSolverLoader.la

libLinearSolverLoader_la_SOURCES = \
  LibraryHandler.c LibraryHandler.h \
  HSLLoader.c HSLLoader.h

AM_CPPFLAGS = -DIPOPT_BUILD \
  -I$(srcdir)/../../Common \
  -I$(srcdir)/../../Algorithm/LinearSolvers \
  $(HSL_CFLAGS)

if BUILD_PARDISOLOADER
  libLinearSolverLoader_la_SOURCES += PardisoLoader.c PardisoLoader.h
endif

libLinearSolverLoader_la_LIBADD = 

includeipoptdir = $(includedir)/coin
includeipopt_HEADERS = HSLLoader.h PardisoLoader.h

# Astyle stuff

ASTYLE_FILES = \
    LibraryHandler.cbak LibraryHandler.hbak \
    HSLLoader.cbak HSLLoader.hbak

ASTYLE = @ASTYLE@
ASTYLEFLAGS = @ASTYLEFLAGS@

DISTCLEANFILES = $(ASTYLE_FILES)

SUFFIXES = .cppbak .hppbak

astyle: $(ASTYLE_FILES)

.h.hbak:
	mv $< $@
	$(ASTYLE) $(ASTYLEFLAGS) < $@ | sed -e 's/  *$$//' > $<
	touch $@

.c.cbak:
	mv $< $@
	$(ASTYLE) $(ASTYLEFLAGS) < $@ | sed -e 's/  *$$//' > $<
	touch $@

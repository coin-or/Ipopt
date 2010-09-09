# Copyright (C) 2004, 2006 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Common Public License.

## $Id$

# Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

AUTOMAKE_OPTIONS = foreign

noinst_LTLIBRARIES = libamplinterface.la
bin_PROGRAMS = ipopt

libamplinterface_la_SOURCES = \
	AmplTNLP.cpp AmplTNLP.hpp

ipopt_SOURCES = \
	ampl_ipopt.cpp

ipopt_LDADD = libamplinterface.la ../../Interfaces/libipopt.la
	
LIBS = $(ASL_LIBS) $(MUMPS_LIBS) $(HSL_LIBS) $(ADDLIBS) $(LAPACK_LIBS) $(BLAS_LIBS)

AM_CPPFLAGS = -I`$(CYGPATH_W) $(srcdir)/../../Common` \
	-I`$(CYGPATH_W) $(srcdir)/../../LinAlg` \
	-I`$(CYGPATH_W) $(srcdir)/../../LinAlg/TMatrices` \
	-I`$(CYGPATH_W) $(srcdir)/../../Algorithm` \
	-I`$(CYGPATH_W) $(srcdir)/../../Interfaces`

AM_CPPFLAGS += $(ASL_CFLAGS)

# This line is necessary to allow VPATH compilation with MS compilers
# on Cygwin
DEFAULT_INCLUDES = -I. -I`$(CYGPATH_W) $(srcdir)` -I$(top_builddir)/inc

# Astyle stuff

ASTYLE_FILES = \
	AmplTNLP.cppbak AmplTNLP.hppbak \
	ampl_ipopt.cppbak

ASTYLE = @ASTYLE@
ASTYLEFLAGS = @ASTYLEFLAGS@

DISTCLEANFILES = $(ASTYLE_FILES)

SUFFIXES = .cppbak .hppbak

astyle: $(ASTYLE_FILES)

.hpp.hppbak:
	mv $< $@
	$(ASTYLE) $(ASTYLEFLAGS) < $@ > $<
	touch $@

.cpp.cppbak:
	mv $< $@
	$(ASTYLE) $(ASTYLEFLAGS) < $@ > $<
	touch $@

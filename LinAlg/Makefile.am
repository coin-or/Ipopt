# Copyright (C) 2004, International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Common Public License.

## $Id$

# Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

AUTOMAKE_OPTIONS = foreign

SUBDIRS = TMatrices

includeipoptdir = $(includedir)/ipopt
includeipopt_HEADERS = \
	IpMatrix.hpp \
	IpSymMatrix.hpp \
	IpVector.hpp

noinst_LIBRARIES = liblinalg.a

liblinalg_a_SOURCES = \
	IpBlas.cpp IpBlas.hpp \
	IpCompoundMatrix.cpp IpCompoundMatrix.hpp \
	IpCompoundSymMatrix.cpp IpCompoundSymMatrix.hpp \
	IpCompoundVector.cpp IpCompoundVector.hpp \
	IpDenseVector.cpp IpDenseVector.hpp \
	IpDiagMatrix.cpp IpDiagMatrix.hpp \
	IpExpansionMatrix.cpp IpExpansionMatrix.hpp \
	IpIdentityMatrix.cpp IpIdentityMatrix.hpp \
	IpMatrix.cpp IpMatrix.hpp \
	IpScaledMatrix.cpp IpScaledMatrix.hpp \
	IpSumMatrix.cpp IpSumMatrix.hpp \
	IpSumSymMatrix.cpp IpSumSymMatrix.hpp \
	IpSymMatrix.cpp IpSymMatrix.hpp \
	IpSymScaledMatrix.cpp IpSymScaledMatrix.hpp \
	IpVector.cpp IpVector.hpp \
	IpZeroMatrix.cpp IpZeroMatrix.hpp

CXXFLAGS += -I$(srcdir)/../Common

# Astyle stuff

ASTYLE_FILES = IpBlas.cppbak IpBlas.hppbak \
	IpCompoundMatrix.cppbak IpCompoundMatrix.hppbak \
	IpCompoundSymMatrix.cppbak IpCompoundSymMatrix.hppbak \
	IpCompoundVector.cppbak IpCompoundVector.hppbak \
	IpDenseVector.cppbak IpDenseVector.hppbak \
	IpDiagMatrix.cppbak IpDiagMatrix.hppbak \
	IpExpansionMatrix.cppbak IpExpansionMatrix.hppbak \
	IpIdentityMatrix.cppbak IpIdentityMatrix.hppbak \
	IpMatrix.cppbak IpMatrix.hppbak \
	IpScaledMatrix.cppbak IpScaledMatrix.hppbak \
	IpSumMatrix.cppbak IpSumMatrix.hppbak \
	IpSumSymMatrix.cppbak IpSumSymMatrix.hppbak \
	IpSymMatrix.cppbak IpSymMatrix.hppbak \
	IpSymScaledMatrix.cppbak IpSymScaledMatrix.hppbak \
	IpVector.cppbak IpVector.hppbak \
	IpZeroMatrix.cppbak IpZeroMatrix.hppbak

ASTYLE = @ASTYLE@
ASTYLEFLAGS = @ASTYLEFLAGS@

CLEANFILES = $(ASTYLE_FILES)

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

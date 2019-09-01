# Copyright (C) 2004, 2010 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.

# Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

AUTOMAKE_OPTIONS = foreign

SUBDIRS = TMatrices

includeipoptdir = $(includedir)/coin-or
includeipopt_HEADERS = \
	IpBlas.hpp \
	IpCompoundMatrix.hpp \
	IpCompoundSymMatrix.hpp \
	IpCompoundVector.hpp \
	IpDenseVector.hpp \
	IpDiagMatrix.hpp \
	IpExpansionMatrix.hpp \
	IpIdentityMatrix.hpp \
	IpLapack.hpp \
	IpMatrix.hpp \
	IpScaledMatrix.hpp \
	IpSumSymMatrix.hpp \
	IpSymMatrix.hpp \
	IpSymScaledMatrix.hpp \
	IpVector.hpp \
	IpZeroSymMatrix.hpp

noinst_LTLIBRARIES = liblinalg.la

liblinalg_la_SOURCES = \
	IpBlas.cpp \
	IpCompoundMatrix.cpp \
	IpCompoundSymMatrix.cpp \
	IpCompoundVector.cpp \
	IpDenseGenMatrix.cpp \
	IpDenseSymMatrix.cpp \
	IpDenseVector.cpp \
	IpDiagMatrix.cpp \
	IpExpandedMultiVectorMatrix.cpp \
	IpExpansionMatrix.cpp \
	IpIdentityMatrix.cpp \
	IpLapack.cpp \
	IpLowRankUpdateSymMatrix.cpp \
	IpMatrix.cpp \
	IpMultiVectorMatrix.cpp \
	IpScaledMatrix.cpp \
	IpSumMatrix.cpp \
	IpSumSymMatrix.cpp \
	IpSymScaledMatrix.cpp \
	IpTransposeMatrix.cpp \
	IpVector.cpp \
	IpZeroMatrix.cpp \
	IpZeroSymMatrix.cpp

liblinalg_la_LDFLAGS = $(LT_LDFLAGS)

AM_CPPFLAGS = -DIPOPT_BUILD \
	-I$(srcdir)/../Common

# Astyle stuff

ASTYLE_FILES = IpBlas.cppbak IpBlas.hppbak \
	IpCompoundMatrix.cppbak IpCompoundMatrix.hppbak \
	IpCompoundSymMatrix.cppbak IpCompoundSymMatrix.hppbak \
	IpCompoundVector.cppbak IpCompoundVector.hppbak \
	IpDenseGenMatrix.cppbak IpDenseGenMatrix.hppbak \
	IpDenseSymMatrix.cppbak IpDenseSymMatrix.hppbak \
	IpDenseVector.cppbak IpDenseVector.hppbak \
	IpDiagMatrix.cppbak IpDiagMatrix.hppbak \
	IpExpandedMultiVectorMatrix.cppbak IpExpandedMultiVectorMatrix.hppbak \
	IpExpansionMatrix.cppbak IpExpansionMatrix.hppbak \
	IpIdentityMatrix.cppbak IpIdentityMatrix.hppbak \
	IpLapack.cppbak IpLapack.hppbak \
	IpLowRankUpdateSymMatrix.cppbak IpLowRankUpdateSymMatrix.hppbak \
	IpMatrix.cppbak IpMatrix.hppbak \
	IpMultiVectorMatrix.cppbak IpMultiVectorMatrix.hppbak \
	IpScaledMatrix.cppbak IpScaledMatrix.hppbak \
	IpSumMatrix.cppbak IpSumMatrix.hppbak \
	IpSumSymMatrix.cppbak IpSumSymMatrix.hppbak \
	IpSymMatrix.hppbak \
	IpSymScaledMatrix.cppbak IpSymScaledMatrix.hppbak \
	IpTransposeMatrix.cppbak IpTransposeMatrix.hppbak \
	IpVector.cppbak IpVector.hppbak \
	IpZeroMatrix.cppbak IpZeroMatrix.hppbak \
	IpZeroSymMatrix.cppbak IpZeroSymMatrix.hppbak

ASTYLE = @ASTYLE@
ASTYLEFLAGS = @ASTYLEFLAGS@

DISTCLEANFILES = $(ASTYLE_FILES)

SUFFIXES = .cppbak .hppbak

astyle: $(ASTYLE_FILES)

.hpp.hppbak:
	mv $< $@
	$(ASTYLE) $(ASTYLEFLAGS) < $@ | sed -e 's/  *$$//' > $<
	touch $@

.cpp.cppbak:
	mv $< $@
	$(ASTYLE) $(ASTYLEFLAGS) < $@ | sed -e 's/  *$$//' > $<
	touch $@

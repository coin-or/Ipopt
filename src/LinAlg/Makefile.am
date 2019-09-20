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

#liblinalg_la_LDFLAGS = $(LT_LDFLAGS)

AM_CPPFLAGS = -DIPOPT_BUILD \
	-I$(srcdir)/../Common

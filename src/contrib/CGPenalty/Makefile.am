# Copyright (C) 2007 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.
#
# Authors:  Andreas Waechter        IBM       2007-06-01

noinst_LTLIBRARIES = libcgpenalty.la

libcgpenalty_la_SOURCES = \
	IpCGPenaltyCq.cpp \
	IpCGPenaltyData.cpp \
	IpCGPenaltyLSAcceptor.cpp \
	IpCGPenaltyRegOp.cpp \
	IpCGPerturbationHandler.cpp \
	IpCGSearchDirCalc.cpp \
	IpPiecewisePenalty.cpp

AM_CPPFLAGS = \
	-I$(srcdir)/../../Common \
	-I$(srcdir)/../../LinAlg \
	-I$(srcdir)/../../LinAlg/TMatrices \
	-I$(srcdir)/../../Interfaces \
	-I$(srcdir)/../../Algorithm \
	$(IPOPTLIB_CFLAGS)

# Copyright (C) 2003, 2006 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Common Public License.

## $Id$

AUTOMAKE_OPTIONS = foreign

EXTRA_DIST = get.blas INSTALL.BLAS

iplibdir = $(top_builddir)/lib

noinst_LIBRARIES = libipblas.a

nodist_libipblas_a_SOURCES =

if BUILDBLAS
nodist_libipblas_a_SOURCES += dasum.f \
	daxpy.f \
	dcopy.f \
	ddot.f \
	dgbmv.f \
	dgemm.f \
	dgemv.f \
	dger.f \
	dnrm2.f \
	drot.f \
	drotg.f \
	drotm.f \
	drotmg.f \
	dsbmv.f \
	dscal.f \
	dsdot.f \
	dspmv.f \
	dspr.f \
	dspr2.f \
	dswap.f \
	dsymm.f \
	dsymv.f \
	dsyr.f \
	dsyr2.f \
	dsyr2k.f \
	dsyrk.f \
	dtbmv.f \
	dtbsv.f \
	dtpmv.f \
	dtpsv.f \
	dtrmm.f \
	dtrmv.f \
	dtrsm.f \
	dtrsv.f \
	idamax.f \
	lsame.f \
	xerbla.f
endif

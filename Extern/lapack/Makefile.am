# Copyright (C) 2006 International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Common Public License.

## $Id$

AUTOMAKE_OPTIONS = foreign

EXTRA_DIST = get.lapack INSTALL.LAPACK

iplibdir = $(top_builddir)/lib

noinst_LIBRARIES = libiplapack.a

nodist_libiplapack_a_SOURCES =

LAPACK_FILES = \
	dlae2.f \
	dlaev2.f \
	dlamch.f \
	dlanst.f \
	dlansy.f \
	dlapy2.f \
	dlarf.f \
	dlarfb.f \
	dlarfg.f \
	dlarft.f \
	dlartg.f \
	dlascl.f \
	dlaset.f \
	dlasr.f \
	dlasrt.f \
	dlassq.f \
	dlatrd.f \
	dorg2l.f \
	dorg2r.f \
	dorgql.f \
	dorgqr.f \
	dorgtr.f \
	dpotf2.f \
	dpotrf.f \
	dpotrs.f \
	dsteqr.f \
	dsterf.f \
	dsyev.f \
	dsytd2.f \
	dsytrd.f \
	ieeeck.f \
	ilaenv.f

if BUILDLAPACK
nodist_libiplapack_a_SOURCES += $(LAPACK_FILES)
endif

DISTCLEANFILES = $(LAPACK_FILES)


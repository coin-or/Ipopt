# Copyright (C) 2003, International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Common Public License.

## $Id$

AUTOMAKE_OPTIONS = foreign

noinst_LIBRARIES = libiplapack.a

lapacksrc = dgesv.f \
	dgetf2.f \
	dgetrf.f \
	dgetrs.f \
	dlabad.f \
	dlacon.f \
	dlacpy.f \
	dlae2.f \
	dlaebz.f \
	dlaev2.f \
	dlagtf.f \
	dlagts.f \
	dlamch.f \
	dlansp.f \
	dlanst.f \
	dlapy2.f \
	dlaqsp.f \
	dlarf.f \
	dlarfg.f \
	dlarnv.f \
	dlartg.f \
	dlaruv.f \
	dlascl.f \
	dlaset.f \
	dlasr.f \
	dlasrt.f \
	dlassq.f \
	dlaswp.f \
	dlatps.f \
	dopgtr.f \
	dopmtr.f \
	dorg2l.f \
	dorg2r.f \
	dppcon.f \
	dppequ.f \
	dpprfs.f \
	dppsv.f \
	dppsvx.f \
	dpptrf.f \
	dpptrs.f \
	drscl.f \
	dspev.f \
	dspevx.f \
	dsptrd.f \
	dstebz.f \
	dstein.f \
	dsteqr.f \
	dsterf.f \
	ieeeck.f \
	ilaenv.f

nodist_libiplapack_a_SOURCES = $(lapacksrc)

distclean-local:
	rm -f $(lapacksrc)

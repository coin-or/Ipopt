# Copyright (C) 2003, International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Common Public License.

## $Id: Makefile.am 632 2004-07-20 23:52:15Z andreasw $

AUTOMAKE_OPTIONS = foreign

noinst_LIBRARIES = libiptron.a

patchsrc = dtron_p.f dspcg_p.f dcauchy_p.f dtrpcg_p.f dgpstep_p.f dprsrch_p.f

origsrc = dtrqsol.f dmid.f dicfs.f dssyax.f dstrsol.f dicf.f ihsort.f \
	insort.f dsel2.f dgpnrm2.f dbreakpt.f

EXTRA_DIST = dcauchy_p.patch dgpstep_p.patch dprsrch_p.patch dspcg_p.patch \
	dtron_p.patch dtrpcg_p.patch

PATCH=@PATCH@

tronsrcdir = $(srcdir)/tron/src/tron/

nodist_libiptron_a_SOURCES = $(origsrc) $(patchsrc)

$(patchsrc):
	cp $(tronsrcdir)/`echo $@ | sed s/_p.f/.f/` $@
	$(PATCH) $@ < $(srcdir)/`echo $@ | sed s/_p.f/_p.patch/`

FFLAGS += -I$(srcdir)/../../IPOPT/ipopt

CLEANFILES = $(patchsrc)

distclean-local:
	rm -f $(origsrc)

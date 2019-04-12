# Copyright (C) 2003, International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Common Public License.

## $Id: Makefile.am 632 2004-07-20 23:52:15Z andreasw $

AUTOMAKE_OPTIONS = foreign

IPINCLUDEDIR = @IPINCLUDEDIR@
FPP = @FPP@
IPADDLIBS = @IPADDLIBS@

EXTRA_DIST = $(preproc_sources) CUTE.INC

ipoptma_files = CUTEr_interface_main.f \
	eval_a.f \
	eval_c.f \
	eval_f.f \
	eval_g.f \
	eval_h.f \
	eval_hesslag_v.f \
	cute_init.f \
	eval_hesscon_v.f \
	eval_hessobj_v.f

preproc_sources = CUTEr_interface_main.F cute_init.F eval_a.F eval_c.F \
	eval_f.F eval_g.F eval_h.F eval_hesscon_v.F eval_hesslag_v.F \
	eval_hessobj_v.F

CUTEr_interface_main.f: CUTEr_interface_main.F

#all-local: ipoptma.f

ipoptma.f: $(ipoptma_files)
	cat $(ipoptma_files) >$@

CLEANFILES = $(ipoptma_files) ipoptma.f

$(preproc_sources): $(IPINCLUDEDIR)/config.h

.F.f:
	$(FPP) $(DEFS) -I$(IPINCLUDEDIR) $< >$@

cuter:	ipoptma.f
	if test -z "$$MYCUTER" -o -z "$$CUTER"; then \
	  echo CUTEr environment variables not set - cannot install ipoptma.f; \
	else \
	  echo "-L$(libdir) -lipopt $(IPADDLIBS)" >$$MYCUTER/double/lib/ipopt.liblist ;\
	  cp ipoptma.f $$CUTER/common/src/tools/ ;\
	  cd $$MYCUTER; $(MAKE) ;\
	fi

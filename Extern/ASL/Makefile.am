# Copyright (C) 2003, International Business Machines and others.
# All Rights Reserved.
# This file is distributed under the Common Public License.

## $Id$

AUTOMAKE_OPTIONS = foreign

compdir = compdir

libampl   = amplsolver.a

LN_S=@LN_S@

ASLMAKEFLAGS=@ASLMAKEFLAGS@
ASLMAKEFILECMDS=@ASLMAKEFILECMDS@

all-local: $(libampl) stdio1.h arith.h funcadd0.o

$(libampl):
	set -x; \
	rm -rf $(compdir) ; \
	$(mkinstalldirs) $(compdir); \
	files=`cd ASLsolvers; ls -d *`; \
	cd $(compdir); \
	for file in $$files; \
	do \
	  $(LN_S) ../ASLsolvers/$$file $$file ;\
	done; \
	rm -f makefile; \
	$(ASLMAKEFILECMDS); \
	$(MAKE) $(ASLMAKEFLAGS) CC="$(CC)" clean ;\
	$(MAKE) $(ASLMAKEFLAGS) CC="$(CC)" ;\
	mv amplsolver.a ../$(libampl)
	-cd $(compdir); mv stdio1.h arith.h funcadd0.o ..
	rm -rf $(compdir)

stdio1.h arith.h funcadd0.o: $(libampl)

CLEANFILES = $(libampl) stdio1.h arith.h funcadd0.o

# ASLsolvers

distclean-local:
	rm -f ASLsolvers

# Copyright (C) 2002, Carnegie Mellon University and others.
# All Rights Reserved.
# This code is published under the Common Public License.

# $Id: Makefile 531 2004-03-11 01:31:07Z andreasw $
THIS_DIST_DIR   := DYNOPT/dae2nlp
#
#  The next lines automatically determine the root dir for the objects
#
slash := /
empty :=
space := $(empty) $(empty)
dirs  := $(subst /, ,$(THIS_DIST_DIR))
dots  := $(foreach dir,$(dirs), .. )
OBJROOTDIR := $(subst $(space),$(slash),$(strip $(dots)))

MAKEFILEDEFS = $(OBJROOTDIR)/Makefile.defs
MAKEFILETRGS = $(OBJROOTDIR)/Makefile.trgs

include $(MAKEFILEDEFS)

LIBOBJS = \
	addcon.o \
	appsln.o \
	basics.o \
	colrhs.o \
	dae2nlp.o \
	estimnnz.o \
	fact_c.o \
	fact_fw.o \
	fillblok.o \
	filljac.o \
	lagvec.o \
	mult_c.o \
	mult_n.o \
	multcw.o \
	multmat.o \
	partele.o \
	solve_c.o \
	solve_fw.o \
	structure.o

LIB = libdae2nlp.a

ifdef SOLIBDIR
  SOLIB := $(patsubst %.a,%.so,$(LIB))
endif

CLEANOBJS = $(LIBOBJS)
CLEANLIBS = $(LIB) $(patsubst %.a,%.so,$(LIB))

all: $(LIB) $(SOLIB)


#
#  Dependencies
#
addcon.o \
appsln.o \
basics.o \
colrhs.o \
fillblok.o \
lagvec.o \
multcw.o : DAE2NLP.INC

include $(MAKEFILETRGS)

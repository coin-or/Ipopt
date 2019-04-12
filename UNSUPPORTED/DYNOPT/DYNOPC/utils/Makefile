# Copyright (C) 2002, Carnegie Mellon University and others.
# All Rights Reserved.
# This code is published under the Common Public License.

# $Id: Makefile 531 2004-03-11 01:31:07Z andreasw $
THIS_DIST_DIR   := DYNOPT/DYNOPC/utils
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

OBJS = \
	calc_nx.o \
	cprintf.o \
	daspk_call.o \
	eval_f.o \
	eval_g.o \
	eval_hessobj_v.o \
	get_filename.o \
	get_time.o \
	main.o \
	read_adb.o \
	read_cmp.o \
	read_ojf.o \
	read_stp.o \
	write_stp.o


CLEANOBJS = $(OBJS)

all: ;    # don't do anything - create object files only when needed

#
#  Dependencies
#
incdirs    = $(SRCROOTDIR)/$(DIST_DYNOPT_AUXIL)  \
	$(SRCROOTDIR)/$(DIST_DYNOPT_DAE2NLP)
FINCFLAGS  = $(patsubst %, -I%, $(incdirs))

#FFLAGS += -DOLDSTYLE

eval_f.o eval_g.o read_adb.o read_cmp.o read_ojf.o read_stp.o write_stp.o: \
	DYNOPC.INC

eval_f.o eval_g.o daspk_call.o read_adb.o read_cmp.o read_ojf.o \
read_stp.o write_stp.o: \
	$(SRCROOTDIR)/$(DIST_DYNOPT_AUXIL)/DYNAUX.INC

daspk_call.o: \
	$(SRCROOTDIR)/$(DIST_DYNOPT_DAE2NLP)/DAE2NLP.INC


include $(MAKEFILETRGS)

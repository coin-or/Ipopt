# Copyright (C) 2002, Carnegie Mellon University and others.
# All Rights Reserved.
# This code is published under the Common Public License.

# $Id: Makefile 531 2004-03-11 01:31:07Z andreasw $
THIS_DIST_DIR   := DYNOPT/auxil
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
	constr.o \
	daemodel_df.o \
	daemodel_df_struc.o \
	daemodel_f.o \
	daemodel_hv.o \
	get_iele4t.o \
	get_start.o \
	ip2dae_check_derivs.o \
	ip2dae_main.o \
	model_hv_old_dummy.o \
	model_hv_dummy.o \
	model_init_dummy.o \
        myalloc.o \
	readsol.o \
	refinesol.o

CLEANOBJS = $(OBJS)

all: ;    # don't do anything - create object files only when needed

#
#  Dependencies
#
incdirs    = $(SRCROOTDIR)/$(DIST_DYNOPT_DAE2NLP) \
	$(SRCROOTDIR)/$(DIST_IPOPT_IPOPT)
FINCFLAGS  = $(patsubst %, -I%, $(incdirs))

constr.o: $(SRCROOTDIR)/$(DIST_IPOPT_IPOPT)/IPOPT.INC

constr.o \
daemodel_df.o \
daemodel_df_struc.o \
daemodel_f.o \
daemodel_hv.o \
get_iele4t.o \
get_start.o \
ip2dae_check_derivs.o \
ip2dae_main.o \
model_hv_dummy.o \
readsol.o \
refinesol.o: DYNAUX.INC

readsol.o refinesol.o: $(SRCROOTDIR)/$(DIST_DYNOPT_DAE2NLP)/DAE2NLP.INC

include $(MAKEFILETRGS)

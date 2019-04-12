# Copyright (C) 2002, Carnegie Mellon University and others.
# All Rights Reserved.
# This code is published under the Common Public License.

# $Id: Makefile 531 2004-03-11 01:31:07Z andreasw $
THIS_DIST_DIR   := DYNOPT/Fortran_example
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

# Make sure that the files in $SRCROOTDIR/DYNOPT/auxil have been compiled
# without the -DOLDADOLC option
ifneq (,$(findstring -DOLDADOLC,$(DEFS)))
  dummy := $(error You need to compile this example without the -DOLDADOLC definition)
endif
#
#  Need to redefine search path for source code
#
model    = tank
exe      = $(model)_Fortran
execheck = $(exe)_check
modelinc = $(datadir)/TANK_CONSTS.INC

mainobj  = $(OBJROOTDIR)/$(DIST_DYNOPT_AUXIL)/ip2dae_main.o 
maincheckobj = $(OBJROOTDIR)/$(DIST_DYNOPT_AUXIL)/ip2dae_check_derivs.o

hereobjs =  $(model)_f.o $(model)_df.o $(model)_df_struc.o $(model)_hv.o \
	eval.o

ipoptauxilobjfiles = \
	c_out.o \
	iter_out.o \
	eval_h_dummy.o \
	eval_hesslag_v_dummy.o
auxilobjs1 = $(addprefix  $(OBJROOTDIR)/$(DIST_IPOPT_AUXIL)/, \
	$(ipoptauxilobjfiles))

dynoptauxilobjfiles = \
	constr.o \
	daemodel_df.o \
	daemodel_df_struc.o \
	daemodel_f.o \
	daemodel_hv.o \
	get_iele4t.o \
	get_start.o \
	model_init_dummy.o \
	model_hv_old_dummy.o
auxilobjs2 = $(addprefix  $(OBJROOTDIR)/$(DIST_DYNOPT_AUXIL)/, \
	$(dynoptauxilobjfiles))

dynoptauxilobjfilescheck = \
	model_init_dummy.o 
auxilobjscheck = $(addprefix  $(OBJROOTDIR)/$(DIST_DYNOPT_AUXIL)/, \
	$(dynoptauxilobjfilescheck))

AUXILOBJS = $(mainobj) $(maincheckobj) $(auxilobjs1) $(auxilobjs2)

data    = $(notdir $(wildcard $(datadir)/*.DAT))
datadir = $(SRCROOTDIR)/$(DIST_DYNOPT_EXAMPLEDATA)

CLEANOBJS = $(hereobjs) $(mainobj) $(maincheckobj)
CLEANEXES = $(exe) $(execheck)
CLEANDATA = $(data) IPOPT.OUT X.sol

objs    = $(hereobjs) $(auxilobjs1) $(auxilobjs2) $(auxilobjs3)
objscheck = $(hereobjs) $(auxilobjscheck)
LIBS    = $(IPOPTLIB) $(TRONLIB) $(IPOPTLIB) $(DAE2NLPLIB) $(LINALGLIBS)

ifdef SOLIBDIR
  lnklibcheck = -L$(SOLIBDIR) -lblas
else
  lnklibcheck = $(BLASLIB)
endif

all: $(exe) $(execheck) $(data)

$(exe): $(mainobj) $(objs) $(DEPLIBS)
	$(F77) $(LNFLAGS) -o $@ $(mainobj) $(objs) \
	$(LNKLIBS) $(LNLIBS)

$(execheck): $(maincheckobj) $(objscheck)
	$(F77) $(LNFLAGS) -o $@ $(maincheckobj) $(objscheck) \
	$(lnklibcheck) $(LNLIBS)

$(data):
	cp $(datadir)/$@ .
#
#  Mainly for development:
#  Force check of libraries and object files in other directories
#
$(patsubst %.a,%.so, $(DAE2NLPLIB) $(IPOPTLIB)) \
$(AUXILOBJS) : FORCE

#
#  Dependencies
#
incdirs    = $(SRCROOTDIR)/$(DIST_DYNOPT_AUXIL) $(datadir)
FINCFLAGS  = $(patsubst %, -I%, $(incdirs))

eval.o: $(SRCROOTDIR)/$(DIST_DYNOPT_AUXIL)/DYNAUX.INC

$(model)_f.o \
$(model)_df.o \
$(model)_df_struc.o \
$(model)_hv.o: $(modelinc)

#
#  Need to establish some order for parallelism in case FLAGS changed
#
break := $(firstword $(auxilobj1))
$(break): $(filter-out $(break), $(auxilobj1))

$(mainobj): $(auxilobj2) # Need to point these in order to avoid 

include $(MAKEFILETRGS)

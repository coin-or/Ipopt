# Copyright (C) 2002, Carnegie Mellon University and others.
# All Rights Reserved.
# This code is published under the Common Public License.

# $Id: Makefile 531 2004-03-11 01:31:07Z andreasw $
THIS_DIST_DIR   := DYNOPT/DYNOPC/optimization
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

## Uncomment following line to compile ESO version
##ESOVERSION := 1

#
#  Need to redefine search path for source code
#
exe      = DYNOPC_optimization

ifndef PROCESS_DATA_DIR
  PROCESS_DATA_DIR = .
endif

# modellist := $(shell cd $(PROCESS_DATA_DIR); ls -d

model    = $(MODEL)
modelf   = $(PROCESS_DATA_DIR)/$(MODEL)/$(model).f
modelc   = $(model).c
modelcxx = $(model).cc
modelobj = $(model).o
ifdef SOLIBDIR
  modelso  = $(PROCESS_DATA_DIR)/$(MODEL)/$(model).so
else
  modelso = $(error You need to specify SOLIBDIR for DYNOPC.)
endif

get_prof = get_prof

mainobj  = $(OBJROOTDIR)/$(DIST_DYNOPC_UTILS)/main.o

hereobjs =  opt_fmain.o	update_stp.o push_bnd.o set_bnds.o #eval.o

#dynopcutilsobjfiles = \
#	calc_nx.o \
#	cprintf.o \
#	get_filename.o \
#	get_time.o \
#	read_adb.o \
#	read_cmp.o \
#	read_ojf.o \
#	read_stp.o \
#	write_stp.o

dynopcutilsobjfiles = \
	calc_nx.o \
	cprintf.o \
	eval_f.o \
	eval_g.o \
	eval_hessobj_v.o \
	get_filename.o \
	get_time.o \
	read_adb.o \
	read_cmp.o \
	read_ojf.o \
	read_stp.o \
	write_stp.o

auxilobjs1 = $(addprefix  $(OBJROOTDIR)/$(DIST_DYNOPC_UTILS)/, \
	$(dynopcutilsobjfiles))

dynoptauxilobjfiles = \
	constr.o \
        daemodel_df.o \
        daemodel_df_struc.o \
        daemodel_f.o \
        daemodel_hv.o \
	get_iele4t.o \
	myalloc.o

auxilobjs2 = $(addprefix  $(OBJROOTDIR)/$(DIST_DYNOPT_AUXIL)/, \
        $(dynoptauxilobjfiles))

dynoptadolc2daefiles = \
	cprogs.o \
	fprogs.o
auxilobjs3 = $(addprefix  $(OBJROOTDIR)/$(DIST_DYNOPT_ADOLC2DAE)/, \
	$(dynoptadolc2daefiles))

modelinitobj = $(OBJROOTDIR)/$(DIST_DYNOPT_ADOLC2DAE)/model_init.o
modelinitaddobj = $(OBJROOTDIR)/$(DIST_DYNOPT_ADOLC2DAE)/model_init_add.o

auxilobjs4 = $(modelinitobj) $(modelinitaddobj)

ipoptauxilobjfiles = \
	c_out.o \
	iter_out.o \
	eval_h_dummy.o \
	eval_hesslag_v_dummy.o
auxilobjs5 = $(addprefix  $(OBJROOTDIR)/$(DIST_IPOPT_AUXIL)/, \
	$(ipoptauxilobjfiles))

AUXILOBJS = $(auxilobjs1) $(auxilobjs2) $(auxilobjs3) $(auxilobjs4) \
	$(auxilobjs5) $(mainobj)

CLEANOBJS = $(hereobjs) $(modelc) $(modelcxx) $(modelobj) \
	$(wildcard *.o) $(wildcard *.c) $(wildcard *.cc)
CLEANEXES = $(exe)
CLEANDATA = $(data)

ifndef ESOVERSION
objs    = $(hereobjs) $(auxilobjs1) $(auxilobjs2) $(auxilobjs3) \
	$(auxilobjs5) $(NLEQ1SOBJ)
LIBS    = $(DASPKLIB) $(IPOPTLIB) $(TRONLIB) $(IPOPTLIB) $(DAE2NLPLIB) \
	$(ADOLCLIB) $(LINALGLIBS) $(DLLIB)
checkobjs = $(checkhereobjs) $(auxilobjs1) $(auxilobjs2) $(auxilobjs3) \
	$(auxilobjs5) $(NLEQ1SOBJ)
ifdef SOLIBDIR
  objs += $(modelinitaddobj)
  checkobjs += $(modelinitaddobj)
else
  objs += $(modelinitobj) $(modelobj)
  checkobjs += $(modelinitobj) $(modelobj)
endif
else
dynoptauxilobjfiles += model_hv_dummy.o
objs    = $(hereobjs) $(auxilobjs1) $(auxilobjs2) \
	$(auxilobjs5) $(NLEQ1SOBJ)
LIBS    = $(DASPKLIB) $(IPOPTLIB) $(DAE2NLPLIB) $(ESO2DAELIB) \
	$(IDLGPROMSLIB) $(LINALGLIBS) $(CORBALIBS)
checkobjs = $(checkhereobjs) $(auxilobjs1) $(auxilobjs2) \
	$(auxilobjs5) $(NLEQ1SOBJ)
endif

#
# Stuff for the get_prof program:
#

getobjs = get_prof.o update_stp.o $(auxilobjs1) $(auxilobjs2) \
	$(auxilobjs5) 

all: $(exe)

$(exe): $(mainobj) $(objs) $(DEPLIBS)
	$(CXX) -pg  $(CXXLNFLAGS) -o $@ $(mainobj) $(objs) $(LNKLIBS) $(CXXLNLIBS)

$(data):
	cp $(datadir)/$@ .

$(get_prof): $(mainobj) $(getobjs) $(DEPLIBS)
	$(CXX)  $(CXXLNFLAGS) -o $@ $(mainobj) $(getobjs) $(LNKLIBS) $(CXXLNLIBS)

model: $(modelso)
#
#  Special Building rules using f2c to obtain C++ model from Fortran
#
.INTERMEDIATE: $(modelc) $(modelcxx)
$(modelc): $(modelf) $(modelinc)
	$(F2C) $(F2CFLAGS) $(FINCFLAGS) $<

$(modelcxx): $(modelc) $(F2CPOST)
	$(F2CPOST) $< $@

$(modelobj): $(modelcxx) $(f2cadolc_h)
	$(CXX) -c $(CXXFLAGS) $(CXXINCFLAGS) $< -o $@

$(modelso): $(modelinitobj) $(modelobj)
	$(SOAR) $(SOARFLAGS) $@ $(modelinitobj) $(modelobj)

ifdef SOLIBDIR
.INTERMEDIATE: $(modelobj) #$(modelcxx)
endif
#
#  Mainly for development:
#  Force check of libraries and object files in other directories
#
ifdef SOLIBDIR
$(patsubst %.a,%.so, $(DAE2NLPLIB) $(IPOPTLIB) $(ESO2DAELIB) $(IDLGPROMSLIB)) \
$(AUXILOBJS) $(F2CPOST): FORCE
else
$(DAE2NLPLIB) $(IPOPTLIB) $(ESO2DAELIB) $(IDLGPROMSLIB) \
$(AUXILOBJS) $(F2CPOST): FORCE
endif

#
#  Dependencies
#
incdirs    = $(SRCROOTDIR)/$(DIST_DYNOPC_UTILS) \
	$(SRCROOTDIR)/$(DIST_DYNOPT_AUXIL)  \
	$(SRCROOTDIR)/$(DIST_DYNOPT_DAE2NLP)
FINCFLAGS  = $(patsubst %, -I%, $(incdirs))

opt_fmain.o push_bnd.o set_bnds.o update_stp.o: \
	$(SRCROOTDIR)/$(DIST_DYNOPC_UTILS)/DYNOPC.INC

opt_fmain.o push_bnd.o set_bnds.o update_stp.o: \
	$(SRCROOTDIR)/$(DIST_DYNOPT_AUXIL)/DYNAUX.INC

update_stp.o: \
	$(SRCROOTDIR)/$(DIST_DYNOPT_DAE2NLP)/DAE2NLP.INC

#
f2cadol_hdir = $(SRCROOTDIR)/$(DIST_DYNOPT_ADOLC2DAE)
cxxincdirs   = $(f2cadol_hdir) $(ADOLCSRCDIR)
CXXINCFLAGS  = $(patsubst %, -I%, $(cxxincdirs))

f2cadolc_h = $(f2cadol_hdir)/f2c_adolc.h
$(modelcc): $(f2cadolc_h)

#
#  Need to establish some order for parallelism in case FLAGS changed
#
break := $(firstword $(auxilobjs))
$(break): $(filter-out $(break), $(auxilobjs))

$(mainobj): $(hereobjs) # Need to point these in order to avoid 

include $(MAKEFILETRGS)

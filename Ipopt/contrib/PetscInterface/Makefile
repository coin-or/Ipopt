# The location of the mesh library
LIBMESH_DIR ?= $(HOME)/Develop/libmesh/libmesh-trunk

# include the library options determined by configure.  This will
# set the variables INCLUDE and LIBS that we will need to build and
# link with the library.
include $(LIBMESH_DIR)/Make.common

# Add flags for compilation with Ipopt
IPOPTDIR ?= $(HOME)/COIN-svn/OBJgcc-debug/Ipopt/branches/parallel-inexact
libmesh_INCLUDE += -I$(IPOPTDIR)/include/coin

libmesh_LIBS += -L$(IPOPTDIR)/lib/coin/ -lipopt 
libmesh_LIBS += -L$(MUMPS_DIR)/lib -ldmumps -lmumps_common -L$(SCALPACK_DIR) -lscalapack -lblacsF77all
libmesh_LIBS += -L$(METIS_DIR) -lmetis -L$(MUMPS_DIR)/PORD/lib -lpord
libmesh_LIBS += -L$(HOME)/Develop/ATLAS/build/lib -llapack -lcblas -lf77blas -latlas 
libmesh_LIBS += -L$(IPOPTLIBDIR) -L$(IPOPTLIBDIR)/ThirdParty -L$(IPOPTLIBDIR) -lipopt -L/usr/lib/atlas -llapack -lblas -L/usr/lib/atlas -llapack -lblas -L$(MUMPS_DIR)/lib -ldmumps -lmumps_common -L$(METIS_DIR) -lmetis -L$(MUMPS_DIR)/PORD/lib -lpord -L$(IPOPTLIBDIR)/ThirdParty/HSL/libcoinhsl.la -lm   -L/usr/lib -L/usr/lib/gcc/x86_64-linux-gnu/4.4.1 -L/usr/lib/gcc/x86_64-linux-gnu/4.4.1/../../../../lib -L/lib/../lib -L/usr/lib/../lib -L/usr/lib/gcc/x86_64-linux-gnu/4.4.1/../../.. -lmpich -lopa -lpthread -lrt -lgfortranbegin -lgfortran -lm -lgcc_s  

libmesh_LIBS += -L$(HOME)/Code/lib -lscalapack -lblacsF77all -L$(HOME)/lib -lpwsmpGNU64 $(HOME)/lib/libpardiso.so -lpthread -lgomp -lm  -ldl

#libmesh_CXXFLAGS += -DCHECK_SMARTPTR
# libmesh_CXXFLAGS = -DDEBUG  -g -ansi -pedantic -DPETSC_USE_DYNAMIC_LIBRARIES -Wno-long-long  -D_GLIBCXX_DEBUG -D_GLIBCXX_DEBUG_PEDANTIC
# -felide-constructors -llapack_netlib -lblas_netlib  

###############################################################################
# File management.  This is where the source, header, and object files are
# defined

#
# source files
srcfiles 	:= Airflow.C IpLibMeshPDE.C IpLibMeshPDENLP.C ProblemGeometry.C

#
# object files
objects		:= $(patsubst %.C, %.$(obj-suffix), $(srcfiles))
###############################################################################

.PHONY: clean clobber distclean

###############################################################################
# Target:
#
target 	   := ./Airflow


all:: $(target)

# Production rules:  how to make the target - depends on library configuration
$(target): $(objects)
	@echo "Linking "$@"..."
	$(libmesh_CXX) $(libmesh_CXXFLAGS) $(objects) -o $@ $(libmesh_LIBS) $(libmesh_LDFLAGS)


# Useful rules.
clean:
	rm -f $(objects) *~
	rm -f $(target)
clobber:
	@$(MAKE) clean
	@rm -f $(target) out_*

distclean:
	@$(MAKE) clobber
	@rm -f *.o *.g.o *.pg.o

# include the dependency list
include .depend


#
# Dependencies
#
.depend:
	@$(perl) $(LIBMESH_DIR)/contrib/bin/make_dependencies.pl -I. $(foreach i, $(wildcard $(LIBMESH_DIR)/include/*), -I$(i)) "-S\$$(obj-suffix)" $(srcfiles) > .depend

###############################################################################

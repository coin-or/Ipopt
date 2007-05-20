# Linux settings.
MEX           = mex
MEXSUFFIX     = mexglx
MATLAB_HOME   = /cs/local/generic/lib/pkg/matlab-7.3
CXX           = g++-3.4.5
F77           = g77-3.4.5
CFLAGS        = -g -ggdb -fPIC -pthread 
IPOPT_INCLUDE = -I$(HOME)/include
LIBS          = -L$(HOME)/lib/gcc-3.4.5 -lipopt -llapack -lblas -lhsl \
                -lg2c -lm

# Mac OS X settings.
MEX           = mex
MEXSUFFIX     = mexmac
MATLAB_HOME   = /Applications/MATLAB72
CXX           = g++
F77           = g77 
CFLAGS        = -g -ggdb -fPIC -fno-common -fexceptions -no-cpp-precomp \
                -DMWINDEXISINT
IPOPT_INCLUDE = -I$(HOME)/include
LIBS          = -L$(HOME)/lib -lipopt -llapack -lhsl -lg2c -lm

TARGET = ipopt.$(MEXSUFFIX)
OBJS   = matlabexception.o matlabscalar.o matlabstring.o matlabjournal.o \
         matlabmatrix.o arrayofmatrices.o sparsematrix.o matlabprogram.o \
         matlaboption.o ipopt.o

MATLAB_INCLUDE = -I$(MATLAB_HOME)/extern/include
INCLUDE        = $(IPOPT_INCLUDE) $(MATLAB_INCLUDE) 
CFLAGS        += -Wall -ansi -DMATLAB_MEXFILE
MEXFLAGS       = -v -g -cxx CXX=$(CXX) CC=$(CXX) FC=$(FCC) LD=$(CXX)

all: $(TARGET)

%.o: %.cpp
	$(CXX) $(CFLAGS) $(INCLUDE) -o $@ -c $^

$(TARGET): $(OBJS)
	$(MEX) $(MEXFLAGS) $(LIBS) -output $@ $^

# Make the tar archive for public release.
ipoptmatlab.tar.gz: ipoptmatlab.tar
	gzip $^

ipoptmatlab.tar:
	tar -c -v -f $@ Makefile *.h *.hpp *.cpp *.m

clean:
	rm -f *.o $(TARGET) ipoptmatlab.tar ipoptmatlab.tar.gz


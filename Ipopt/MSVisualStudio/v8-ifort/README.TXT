This file describes how to use the Microsoft Visual Studio 2005 and 2010 
project files included in the subdirectories in order to compile a DLL 
version of Ipopt both for MSVC 2005 and 2010.

These DLLs include the BLAS, LAPACK and MUMPS solvers which can be freely
distributed. Alternatively, the solution contains a configuration where the 
resulting DLL uses the (optimized, multi-threaded) Intel MKL versions of 
the BLAS and LAPACK libraries rather than the netlib versions. The current 
configuration requires MSVC 2005 and the Intel Fortran compiler to compile
the Fortran 90 sources of the MUMPS solver. However, the resulting code DLL 
can be used both by MSVC 2005/2008 and MSVC 2010.

The non-free HSL and PARDISO solvers can be dynamically added at runtime by 
compiling these solvers into separate DLLs that are dynamically loaded by the 
Ipopt DLL when these solvers are selected. A separate project to compile the 
HSL dll is included in the Visual Studio solution.

The resulting DLLs (which are available for download from the COIN-OR site), can 
be used in any Visual Studio project without the need for a Fortran compiler 
to compile the Ipopt solver components. This project itself relies on the Intel 
Fortran compiler to compile the third party Fortran components (especially 
since the MUMPS solver uses F90). The solution contains project files for both
the MSVC C++ compiler or the Intel C++ compiler. The binary distribution only 
contains the optimized MKL version of the Fortran solver DLL IpoptFSS.dll, which 
can be re-used by all configurations of the Ipopt solver project itself. The latter
only require a working MSVC compiler.

Below, BASDIR refers to the basic directory of your downloaded Ipopt
package.  This is the directory where you find subdirectories like
"BuildTools"

Here are the steps:

1. Download the Blas, Lapack, and MUMPS and Metis source code as usual; read the
   INSTALL files in the appropriate subdirectories of

   BASDIR\ThirdParty

2. Edit the file

   BASDIR\Ipopt\MSVisualStudio\v8-ifort\Ipopt\IntelPaths.vsprops
   
   to reflect the base path to the Intel compiler on your computer to allow the
   link step succeed in finding the appropriate libraries. You will only need to set these
   if you want to compile IpoptFSS.dll yourself.
   
3. Now you can open the solution file

   BASDIR\Ipopt\MSVisualStudio\v8-ifort\IpOpt-ifort.sln

   After selecting a build configuration (win32/x64 release/MKL release/debug), you can build the 
   Ipopt project to obtain the Ipopt DLL for that configuration. The project dependencies 
   will make sure the all third party components are built as well. Expect the compile step 
   of the MUMPS Fortran project to fail for a couple times until all interdependencies 
   between the F90 files are resolved properly.

4. If you want to use the HSL solver, then you have to obtain the HSL source code as
   usual, and build either the libhsl project if you have both the sources for the MA27 and MA57
   solvers, or libhsl-no-MA57 if you just have access to the MA27 solver. Both projects result in 
   libhsl.dll, which you have to place into your path to make it available to Ipopt. 

5. When linking against Ipopt.dll in your own project, use the IpoptApplicationFactory method
   exported by the DLL to obtain an IpoptApplication pointer rather than using the new operator.
   From there on, you can call the methods of the IpoptApplication object as usual. Note, however,
   that the virtual interfaces exposed (recursively) by the DLL do not provide access to the more 
   exotic interfaces or member functions in the Ipopt library.
   
   Alternatively, you can use the C interface of Ipopt which is also exported by the DLL. 

Author:  Marcel Roelofs (Marcel.Roelofs at aimms.com)
	 Paragon Decision Technology 
	 2013-05-24

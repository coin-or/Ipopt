This file describes how to use the Microsoft Visual Studio 2005 project files
included in the subdirectories in order to compile a DLL version of Ipopt. 
This DLL includes the BLAS, LAPACK and MUMPS solvers which can be freely
distributed. The non-free HSL and PARDISO solvers can be dynamically added 
at runtime by compiling these solvers into separate DLLs that are dynamically 
loaded by the Ipopt DLL when these solvers are selected. A separate project to 
compile the HSL dll is included in the Visual Studio solution.

The resulting DLL (which is available for download from the COIN-OR site), can 
be used in any Visual Studio project without the need for a Fortran compiler 
to compile the Ipopt solver components. This project itself relies on the Intel 
Fortran compiler to compile the third party Fortran components (especially 
since the MUMPS solver uses F90). The solution contains project files for both
the MSVC C++ compiler or the Intel C++ compiler.

Below, BASDIR refers to the basic directory of your downloaded Ipopt
package.  This is the directory where you find subdirectories like
"BuildTools"

Here are the steps:

1. Download the Blas, Lapack, and MUMPS and Metis source code as usual; read the
   INSTALL files in the appropriate subdirectories of

   BASDIR\ThirdParty

2. Now you can open the solution file

   BASDIR\Ipopt\MSVisualStudio\v8-ifort\IpOpt-ifort.sln

   After selecting a build configuration (win32/x64 release/debug), you can build the 
   Ipopt project to obtain the Ipopt DLL for that configuration. The project dependencies 
   will make sure the all third party components are built as well. Expect the compile step 
   of the MUMPS Fortran project to fail for a couple times until all interdependencies 
   between the F90 files are resolved properly.

3. If you want to use the HSL solver, then you have to obtain the HSL source code as
   usual, and build the libhsl project. This results in libhsl.dll, which you have to place
   into your path to make it available to Ipopt. 

4. When linking against Ipopt.dll in your own project, use the IpoptApplicationFactory method
   exported by the DLL to obtain an IpoptApplication pointer rather than using the new operator.
   From there on, you can call the methods of the IpoptApplication object as usual. 
   
   Alternatively, you can use the C interface of Ipopt which is also exported by the DLL. 

   
Author:  Marcel Roelofs (Marcel.Roelofs at aimms.com)
	 Paragon Decision Technology 
	 2009-10-25

The DLLs in this directory export the C interface of the Ipopt solver, as well as its C++ interface 
(which can be accessed via the IpoptApplicationFactory method). By linking against the import library 
you can use these DLLs from within any code compiled using the MSVC 2005/2008/2010 and/or Intel compilers 
for the Windows platform. 

The DLLs include the MUMPS solver by default, and can dynamically load a separate DLL containing the 
(non-free) HSL solvers if available. If you have access to the HSL solver sources, you have to compile 
the DLL for the HSL solvers yourself. The v8-ifort directory in the Ipopt source distribution contains 
a separate libhsl project to create the DLL using the Intel fortran compiler. Alternatively, you can compile
libhsl using f2c and C++ (cf. the v8 MSVC project).
 
To use these DLLs, you need to add the include\coin directory to the include path of your MSVC project.
In addition, you have to make sure that you have HAVE_CONFIG_H defined in your project, and that the Ipopt 
import library for your platform is added to list of  linker inputs of your project. To run your program, 
you have to copy the correct Ipopt DLL and all necessary runtime DLLs for your platform and configuration 
to the directory where your executable is located. 

Please note that the debug and release versions of the Ipopt DLL are *not* binary compatible, because 
the MSVC compiler unfortunately treats std::string objects slightly differently in both configurations. 
Using a release DLL in debug code and vice versa will therefore propably lead to stack corruption, runtime 
errors, or other hard-to-explain crashes.

This binary distribution only contains the optimized MKL version of the Fortran solver DLL IpoptFSS.dll. You 
can re-use this for any configuration of the Ipopt solver itself (Debug/Release, MSVC 2005/2008/2010). 

It also contains just the Release versions of the Ipopt solver itself plus all necessary MSVC runtime DLLs. Because 
the Debug runtime DLLs of MSVC are not redistributable, you have to build the Debug configurations of the Ipopt solver 
yourself using the v8-ifort MSVC solution, if you are developing Debug versions of your own software using 
the Ipopt solver. When doing so, you can safely link against the MKL release version of IpoptFSS.dll.

Author:  Marcel Roelofs (Marcel.Roelofs at aimms.com)
	 Paragon Decision Technology 
	 2013-05-24

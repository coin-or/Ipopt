To compile this file using MSVC, retrieve the appropriate Ipopt-3.x.y-win32-dll.zip file from the
binary distribution of Ipopt, and unpack its contents into the BinaryDLL-Link-Example directory
of the Ipopt distribution.

This will allow you to compile and link the example program using MSVC. To run it you should copy the 
correct version of the Ipopt dll in the directory where the compiled executable is located. 

Please note that the debug and release versions of the Ipopt dll are *not* binary compatible, because 
the MSVC compiler treats std::string arguments of methods in the Ipopt C++ interface slightly differently.
Using a release dll in debug code and vice versa will propably lead to stack corruption, runtime errors, 
or other hard-to-explain crashes.





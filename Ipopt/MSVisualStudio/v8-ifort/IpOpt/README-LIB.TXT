The DLLs in this directory export the C interface of the Ipopt solver, as well as its C++ interface 
(via the IpoptApplicationFactory method). By linking against the import library you can use these DLLs
from within any code compiled using the MSVC and/or Intel compilers for the Windows platform.

The DLLs are statically linked, so do not require any other runtime DLLs. 

Please note that the debug and release versions of the Ipopt dll are *not* binary compatible, because 
the MSVC compiler treats std::string arguments of methods in the Ipopt C++ interface slightly differently.
Using a release dll in debug code and vice versa will propably lead to stack corruption, runtime errors, 
or even hard-to-explain crashes.


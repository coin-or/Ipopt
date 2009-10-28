The DLLs in this directory export the C interface of the Ipopt solver, as well as its C++ interface 
(via the IpoptApplicationFactory method). By linking against the import library you can use these DLLs
from within any code compiled using the MSVC and/or Intel compilers for the Windows platform.

The DLLs are statically linked, so do not require any other runtime DLLs. 

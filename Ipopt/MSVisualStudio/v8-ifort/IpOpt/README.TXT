The DLLs in this directory export the C interface of the Ipopt solver, as well as its C++ interface 
(via the IpoptApplicationFactory method). The C interface can be linked against from within any compiler, 
the C++ interface can only be used from within MSVC generated code. 

The DLLs are statically linked, so do not require any other runtime DLLs. 

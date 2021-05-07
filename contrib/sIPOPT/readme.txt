                            sIPOPT Toolbox for IPOPT
                            ========================

This is the Sensitivity with IPOPT toolbox. Its purpose is to compute
fast approximate solutions when parameters in the NLP change. For
more information on the project please see the implementation paper,
or the project website.

Documentation, tutorials and test examples can be found in the IPOPT
documentation as well as on the project documentation and website.

Just like IPOPT, the sIPOPT code is separated into a library
"libsipopt" that holds the main algorithm, and an executable that acts
as a solver for AMPL. By default, the library is installed in the same
directory as libipopt, and the executable is installed in the same
directory as IPOPT's ampl executable. They are build in the same build
process as for main Ipopt, unless --disable-sipopt has been specified
for configure.

Contact: The sIPOPT code was developed by Hans Pirnay (RWTH-Aachen),
Rodrigo Lopez-Negrete (CMU), and Prof. Lorenz Biegler (CMU). Any
questions / problems / bugs may be sent to the Hans Pirnay or the IPOPT
mailing list.

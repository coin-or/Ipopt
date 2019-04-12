var index =
[
    [ "Overview", "index.html#Overview", [
      [ "Mathematical Background", "index.html#MATHBACKGROUND", null ],
      [ "Availability", "index.html#AVAILABILITY", null ],
      [ "Prerequisites", "index.html#PREREQUISITES", null ],
      [ "How to use Ipopt", "index.html#HOWTOUSE", null ],
      [ "More Information and Contributions", "index.html#MOREINFO", null ],
      [ "History of Ipopt", "index.html#HISTORY_IPOPT", null ],
      [ "History of this document", "index.html#HISTORY_DOC", null ]
    ] ],
    [ "Installing Ipopt", "INSTALL.html", [
      [ "Getting System Packages (Compilers, ...)", "INSTALL.html#SYSTEMPACKAGES", null ],
      [ "Getting the Ipopt Code", "INSTALL.html#GETIPOPT", null ],
      [ "Getting the Ipopt code via git", "INSTALL.html#GETIPOPT_GIT", null ],
      [ "Getting the Ipopt code as a tarball", "INSTALL.html#GETIPOPT_TAR", null ],
      [ "Download External Code", "INSTALL.html#EXTERNALCODE", [
        [ "Download BLAS, LAPACK and ASL", "INSTALL.html#DOWNLOAD_LINALGASL", null ],
        [ "Download HSL Subroutines", "INSTALL.html#DOWNLOAD_HSL", null ],
        [ "Obtaining the MUMPS Linear Solver", "INSTALL.html#DOWNLOAD_MUMPS", null ],
        [ "Obtaining the Linear Solver Pardiso", "INSTALL.html#DOWNLOAD_PARDISO", null ],
        [ "Obtaining the Linear Solver WSMP", "INSTALL.html#DOWNLOAD_WSMP", null ],
        [ "Using the Linear Solver Loader", "INSTALL.html#LINEARSOLVERLOADER", null ],
        [ "Obtaining METIS", "INSTALL.html#DOWNLOAD_METIS", null ]
      ] ],
      [ "Compiling and Installing Ipopt", "INSTALL.html#COMPILEINSTALL", null ],
      [ "Installation on Windows", "INSTALL.html#INSTALL_WINDOWS", [
        [ "Installation with Cygwin using GNU compilers", "INSTALL.html#INSTALL_CYGWIN", null ],
        [ "Installation with Cygwin using the MSVC++ compiler", "INSTALL.html#INSTALL_CYGWINNATIVE", null ],
        [ "Installation with MSYS/MinGW", "INSTALL.html#INSTALL_MINGW", null ]
      ] ],
      [ "Compiling and Installing the Java Interface JIpopt", "INSTALL.html#INSTALL_JAVA", null ],
      [ "Compiling and Installing the R Interface ipoptr", "INSTALL.html#INSTALL_R", null ],
      [ "Compiling and Installing the MATLAB interface", "INSTALL.html#INSTALL_MATLAB", [
        [ "Setting up mex", "INSTALL.html#INSTALL_MATLAB_MEX", null ],
        [ "Adjusting configuration and build of Ipopt", "INSTALL.html#INSTALL_MATLAB_IPOPT", null ],
        [ "Building the MATLAB interface", "INSTALL.html#INSTALL_MATLAB_BUILD", null ],
        [ "Making MATLAB aware of the mex file", "INSTALL.html#INSTALL_MATLAB_ADDPATH", null ],
        [ "Additional notes", "INSTALL.html#INSTALL_MATLAB_NOTES", null ],
        [ "Troubleshooting", "INSTALL.html#INSTALL_MATLAB_TROUBLESHOOT", null ]
      ] ],
      [ "Expert Installation Options for Ipopt", "INSTALL.html#EXPERT_INSTALL", null ]
    ] ],
    [ "Interfacing your NLP to Ipopt", "INTERFACES.html", [
      [ "Using Ipopt through AMPL", "INTERFACES.html#INTERFACE_AMPL", [
        [ "Using Ipopt from the command line", "INTERFACES.html#INTERFACE_AMPL_CL", null ]
      ] ],
      [ "Interfacing with Ipopt through code", "INTERFACES.html#INTERFACE_CODE", [
        [ "The C++ Interface", "INTERFACES.html#INTERFACE_CPP", [
          [ "Coding the Problem Representation", "INTERFACES.html#INTERFACE_CPP_NLP", null ],
          [ "Coding the Executable", "INTERFACES.html#INTERFACE_CPP_MAIN", null ],
          [ "Compiling and Testing the Example", "INTERFACES.html#INTERFACE_CPP_COMPILE", null ],
          [ "Additional methods in TNLP", "INTERFACES.html#INTERFACE_CPP_ADDITIONAL", null ]
        ] ],
        [ "The C Interface", "INTERFACES.html#INTERFACE_C", null ],
        [ "The Fortran Interface", "INTERFACES.html#INTERFACE_FORTRAN", null ],
        [ "The Java Interface JIpopt", "INTERFACES.html#INTERFACE_JAVA", null ],
        [ "The R Interface ipoptr", "INTERFACES.html#INTERFACE_R", null ],
        [ "The MATLAB Interface", "INTERFACES.html#INTERFACE_MATLAB", null ]
      ] ]
    ] ],
    [ "Special Features", "SPECIALS.html", [
      [ "Derivative Checker", "SPECIALS.html#DERIVCHECK", null ],
      [ "Quasi-Newton Approximation of Second Derivatives", "SPECIALS.html#QUASI_NEWTON", null ],
      [ "Warm-Starting Capabilities via AMPL", "SPECIALS.html#AMPL_WARMSTART", null ],
      [ "sIpopt: Optimal Sensitivity Based on Ipopt", "SPECIALS.html#SIPOPT", null ],
      [ "Inertia-Free Curvature Test", "SPECIALS.html#INERTIAFREE_CURVTEST", null ]
    ] ],
    [ "Ipopt Options", "OPTIONS.html", [
      [ "Options Reference", "OPTIONS.html#OPTIONS_REF", [
        [ "Output", "OPTIONS.html#OPT_Output", null ],
        [ "Termination", "OPTIONS.html#OPT_Termination", null ],
        [ "NLP Scaling", "OPTIONS.html#OPT_NLP_Scaling", null ],
        [ "NLP", "OPTIONS.html#OPT_NLP", null ],
        [ "Initialization", "OPTIONS.html#OPT_Initialization", null ],
        [ "Barrier Parameter", "OPTIONS.html#OPT_Barrier_Parameter", null ],
        [ "Multiplier Updates", "OPTIONS.html#OPT_Multiplier_Updates", null ],
        [ "Line Search", "OPTIONS.html#OPT_Line_Search", null ],
        [ "Warm Start", "OPTIONS.html#OPT_Warm_Start", null ],
        [ "Restoration Phase", "OPTIONS.html#OPT_Restoration_Phase", null ],
        [ "Linear Solver", "OPTIONS.html#OPT_Linear_Solver", null ],
        [ "Hessian Perturbation", "OPTIONS.html#OPT_Hessian_Perturbation", null ],
        [ "Quasi-Newton", "OPTIONS.html#OPT_Quasi-Newton", null ],
        [ "Derivative Test", "OPTIONS.html#OPT_Derivative_Test", null ],
        [ "MA27 Linear Solver", "OPTIONS.html#OPT_MA27_Linear_Solver", null ],
        [ "MA57 Linear Solver", "OPTIONS.html#OPT_MA57_Linear_Solver", null ],
        [ "MA77 Linear Solver", "OPTIONS.html#OPT_MA77_Linear_Solver", null ],
        [ "MA86 Linear Solver", "OPTIONS.html#OPT_MA86_Linear_Solver", null ],
        [ "MA97 Linear Solver", "OPTIONS.html#OPT_MA97_Linear_Solver", null ],
        [ "MUMPS Linear Solver", "OPTIONS.html#OPT_MUMPS_Linear_Solver", null ],
        [ "Pardiso Linear Solver", "OPTIONS.html#OPT_Pardiso_Linear_Solver", null ]
      ] ],
      [ "Options available via the AMPL Interface", "OPTIONS.html#OPTIONS_AMPL", null ]
    ] ],
    [ "Ipopt Output", "OUTPUT.html", null ],
    [ "Implementation Details", "IMPL.html", [
      [ "Triplet Format for Sparse Matrices", "IMPL.html#TRIPLET", null ],
      [ "The Smart Pointer Implementation: SmartPtr<T>", "IMPL.html#SMARTPTR", null ]
    ] ],
    [ "Frequenty Asked Questions", "FAQ.html", null ]
];
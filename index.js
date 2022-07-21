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
      [ "Getting System Packages (Compilers, ...)", "INSTALL.html#SYSTEMPACKAGES", [
        [ "Linux distributions", "INSTALL.html#SYSTEMPACKAGES_LINUX", null ],
        [ "macOS", "INSTALL.html#SYSTEMPACKAGES_MAC", null ],
        [ "Windows with MSYS2/MinGW", "INSTALL.html#SYSTEMPACKAGES_MSYS", null ]
      ] ],
      [ "Download, build, and install dependencies", "INSTALL.html#EXTERNALCODE", [
        [ "ASL (Ampl Solver Library)", "INSTALL.html#EXTERNALCODE_ASL", null ],
        [ "BLAS and LAPACK", "INSTALL.html#EXTERNALCODE_LINALG", null ],
        [ "HSL (Harwell Subroutines Library)", "INSTALL.html#DOWNLOAD_HSL", [
          [ "Providing a HSL library at runtime", "INSTALL.html#LINEARSOLVERLOADER", null ]
        ] ],
        [ "MUMPS Linear Solver", "INSTALL.html#DOWNLOAD_MUMPS", null ],
        [ "Pardiso (Parallel Sparse Direct Linear Solver) from Pardiso Project", "INSTALL.html#DOWNLOAD_PARDISO", null ],
        [ "Pardiso from Intel MKL", "INSTALL.html#DOWNLOAD_PARDISOMKL", null ],
        [ "SPRAL (Sparse Parallel Robust Algorithms Library)", "INSTALL.html#DOWNLOAD_SPRAL", null ],
        [ "WSMP (Watson Sparse Matrix Package)", "INSTALL.html#DOWNLOAD_WSMP", null ]
      ] ],
      [ "Getting the Ipopt Code", "INSTALL.html#GETIPOPT", [
        [ "Getting the Ipopt code via git", "INSTALL.html#GETIPOPT_GIT", null ],
        [ "Getting the Ipopt code as a tarball", "INSTALL.html#GETIPOPT_TAR", null ]
      ] ],
      [ "Compiling and Installing Ipopt", "INSTALL.html#COMPILEINSTALL", [
        [ "Flags to configure", "INSTALL.html#CONFIGURE_FLAGS", null ],
        [ "Additional Flags for Compiler Preprocessor", "INSTALL.html#CPP_FLAGS", null ]
      ] ],
      [ "Using CoinBrew", "INSTALL.html#COINBREW", null ],
      [ "Compiling and Installing the R Interface ipoptr", "INSTALL.html#INSTALL_R", null ],
      [ "Building for single-precision floating-point arithmetic", "INSTALL.html#SINGLEPRECISION_BUILD", null ],
      [ "Building for 64-bit integers", "INSTALL.html#INT64_BUILD", null ]
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
        [ "The R Interface ipoptr", "INTERFACES.html#INTERFACE_R", null ]
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
      [ "Options in the AMPL Interface", "OPTIONS.html#OPTIONS_AMPL", null ],
      [ "Options Reference", "OPTIONS.html#OPTIONS_REF", [
        [ "Termination", "OPTIONS.html#OPT_Termination", null ],
        [ "Output", "OPTIONS.html#OPT_Output", null ],
        [ "NLP", "OPTIONS.html#OPT_NLP", null ],
        [ "NLP Scaling", "OPTIONS.html#OPT_NLP_Scaling", null ],
        [ "Initialization", "OPTIONS.html#OPT_Initialization", null ],
        [ "Warm Start", "OPTIONS.html#OPT_Warm_Start", null ],
        [ "Miscellaneous", "OPTIONS.html#OPT_Miscellaneous", null ],
        [ "Barrier Parameter Update", "OPTIONS.html#OPT_Barrier_Parameter_Update", null ],
        [ "Line Search", "OPTIONS.html#OPT_Line_Search", null ],
        [ "Linear Solver", "OPTIONS.html#OPT_Linear_Solver", null ],
        [ "Step Calculation", "OPTIONS.html#OPT_Step_Calculation", null ],
        [ "Restoration Phase", "OPTIONS.html#OPT_Restoration_Phase", null ],
        [ "Hessian Approximation", "OPTIONS.html#OPT_Hessian_Approximation", null ],
        [ "Derivative Checker", "OPTIONS.html#OPT_Derivative_Checker", null ],
        [ "MA27 Linear Solver", "OPTIONS.html#OPT_MA27_Linear_Solver", null ],
        [ "MA57 Linear Solver", "OPTIONS.html#OPT_MA57_Linear_Solver", null ],
        [ "MA77 Linear Solver", "OPTIONS.html#OPT_MA77_Linear_Solver", null ],
        [ "MA86 Linear Solver", "OPTIONS.html#OPT_MA86_Linear_Solver", null ],
        [ "MA97 Linear Solver", "OPTIONS.html#OPT_MA97_Linear_Solver", null ],
        [ "Pardiso (pardiso-project.org) Linear Solver", "OPTIONS.html#OPT_Pardiso__pardiso_project_org__Linear_Solver", null ],
        [ "Pardiso (MKL) Linear Solver", "OPTIONS.html#OPT_Pardiso__MKL__Linear_Solver", null ],
        [ "SPRAL Linear Solver", "OPTIONS.html#OPT_SPRAL_Linear_Solver", null ],
        [ "WSMP Linear Solver", "OPTIONS.html#OPT_WSMP_Linear_Solver", null ],
        [ "Mumps Linear Solver", "OPTIONS.html#OPT_Mumps_Linear_Solver", null ],
        [ "MA28 Linear Solver", "OPTIONS.html#OPT_MA28_Linear_Solver", null ]
      ] ]
    ] ],
    [ "Ipopt Output", "OUTPUT.html", null ],
    [ "Implementation Details", "IMPL.html", [
      [ "Triplet Format for Sparse Matrices", "IMPL.html#TRIPLET", null ],
      [ "The Smart Pointer Implementation: SmartPtr<T>", "IMPL.html#SMARTPTR", null ]
    ] ],
    [ "Frequenty Asked Questions", "FAQ.html", null ],
    [ "ChangeLog", "md_ChangeLog.html", [
      [ "3.14", "md_ChangeLog.html#autotoc_md10", [
        [ "What is Ipopt?", "FAQ.html#autotoc_md0", null ],
        [ "How do I use Ipopt?", "FAQ.html#autotoc_md1", null ],
        [ "What license is Ipopt released under?", "FAQ.html#autotoc_md2", null ],
        [ "What do I need to build Ipopt?", "FAQ.html#autotoc_md3", null ],
        [ "On what operating systems can Ipopt be used?", "FAQ.html#autotoc_md4", null ],
        [ "Is Ipopt thread-safe?", "FAQ.html#autotoc_md5", null ],
        [ "What is the method behind Ipopt?", "FAQ.html#autotoc_md6", null ],
        [ "Who do I contact with questions about Ipopt?", "FAQ.html#autotoc_md7", null ],
        [ "What is the difference between the Fortran version and the C++ version?", "FAQ.html#autotoc_md8", null ],
        [ "3.14.9 (2022-07-21)", "md_ChangeLog.html#autotoc_md11", null ],
        [ "3.14.8 (2022-07-13)", "md_ChangeLog.html#autotoc_md12", null ],
        [ "3.14.7 (2022-06-24)", "md_ChangeLog.html#autotoc_md13", null ],
        [ "3.14.6 (2022-05-02)", "md_ChangeLog.html#autotoc_md14", null ],
        [ "3.14.5 (2022-02-09)", "md_ChangeLog.html#autotoc_md15", null ],
        [ "3.14.4 (2021-09-20)", "md_ChangeLog.html#autotoc_md16", null ],
        [ "3.14.3 (2021-09-03)", "md_ChangeLog.html#autotoc_md17", null ],
        [ "3.14.2 (2021-07-21)", "md_ChangeLog.html#autotoc_md18", null ],
        [ "3.14.1 (2021-06-25)", "md_ChangeLog.html#autotoc_md19", null ],
        [ "3.14.0 (2021-06-15)", "md_ChangeLog.html#autotoc_md20", null ]
      ] ],
      [ "3.13", "md_ChangeLog.html#autotoc_md21", [
        [ "3.13.4 (2021-02-24)", "md_ChangeLog.html#autotoc_md22", null ],
        [ "3.13.3 (2020-10-16)", "md_ChangeLog.html#autotoc_md23", null ],
        [ "3.13.2 (2020-04-30)", "md_ChangeLog.html#autotoc_md24", null ],
        [ "3.13.1 (2020-03-11)", "md_ChangeLog.html#autotoc_md25", null ],
        [ "3.13.0 (2019-10-19)", "md_ChangeLog.html#autotoc_md26", null ]
      ] ],
      [ "3.12", "md_ChangeLog.html#autotoc_md27", [
        [ "3.12.13 (2019-04-08)", "md_ChangeLog.html#autotoc_md28", null ],
        [ "3.12.12 (2018-11-17)", "md_ChangeLog.html#autotoc_md29", null ],
        [ "3.12.11 (2018-09-16)", "md_ChangeLog.html#autotoc_md30", null ],
        [ "3.12.10 (2018-06-02)", "md_ChangeLog.html#autotoc_md31", null ],
        [ "3.12.9 (2018-01-15)", "md_ChangeLog.html#autotoc_md32", null ],
        [ "3.12.8 (2017-06-12)", "md_ChangeLog.html#autotoc_md33", null ],
        [ "3.12.7 (2017-02-25)", "md_ChangeLog.html#autotoc_md34", null ],
        [ "3.12.6 (2016-07-20)", "md_ChangeLog.html#autotoc_md35", null ],
        [ "3.12.5 (2016-04-30)", "md_ChangeLog.html#autotoc_md36", null ],
        [ "3.12.4 (2015-08-09)", "md_ChangeLog.html#autotoc_md37", null ],
        [ "3.12.3 and 3.11.11 (2015-04-15)", "md_ChangeLog.html#autotoc_md38", null ],
        [ "3.12.2 (2015-04-04)", "md_ChangeLog.html#autotoc_md39", null ],
        [ "3.12.1 (2015-02-13)", "md_ChangeLog.html#autotoc_md40", null ],
        [ "3.12.0 (2015-01-23)", "md_ChangeLog.html#autotoc_md41", null ]
      ] ],
      [ "3.11", "md_ChangeLog.html#autotoc_md42", [
        [ "3.11.10 (2015-01-18)", "md_ChangeLog.html#autotoc_md43", null ],
        [ "3.11.9 (2014-08-16)", "md_ChangeLog.html#autotoc_md44", null ],
        [ "3.11.8 (2014-04-08)", "md_ChangeLog.html#autotoc_md45", null ],
        [ "3.11.7 (2013-12-18)", "md_ChangeLog.html#autotoc_md46", null ],
        [ "3.11.6 (2013-11-16)", "md_ChangeLog.html#autotoc_md47", null ],
        [ "3.11.5 (2013-10-26)", "md_ChangeLog.html#autotoc_md48", null ],
        [ "3.11.4 (2013-09-12)", "md_ChangeLog.html#autotoc_md49", null ],
        [ "3.11.3 (2013-08-08)", "md_ChangeLog.html#autotoc_md50", null ],
        [ "3.11.2 (2013-07-01)", "md_ChangeLog.html#autotoc_md51", null ],
        [ "3.11.1 (2013-06-14)", "md_ChangeLog.html#autotoc_md52", null ],
        [ "3.11.0 (2013-05-07)", "md_ChangeLog.html#autotoc_md53", null ]
      ] ],
      [ "3.10", "md_ChangeLog.html#autotoc_md54", [
        [ "3.10.4 (2013-05-05)", "md_ChangeLog.html#autotoc_md55", null ],
        [ "3.10.3 (2012-11-19)", "md_ChangeLog.html#autotoc_md56", null ],
        [ "3.10.2 (2012-02-12)", "md_ChangeLog.html#autotoc_md57", null ],
        [ "3.10.1 (2011-09-20)", "md_ChangeLog.html#autotoc_md58", null ],
        [ "3.10.0 (2011-06-20)", "md_ChangeLog.html#autotoc_md59", null ]
      ] ],
      [ "3.9", "md_ChangeLog.html#autotoc_md60", [
        [ "3.9.3 (2011-04-07)", "md_ChangeLog.html#autotoc_md61", null ],
        [ "3.9.2 (2010-12-22)", "md_ChangeLog.html#autotoc_md62", null ],
        [ "3.9.1 (2010-11-26)", "md_ChangeLog.html#autotoc_md63", null ],
        [ "3.9.0 (2010-11-05)", "md_ChangeLog.html#autotoc_md64", null ]
      ] ],
      [ "3.8", "md_ChangeLog.html#autotoc_md65", [
        [ "3.8.3 (2010-06-29)", "md_ChangeLog.html#autotoc_md66", null ],
        [ "3.8.2 (2010-06-16)", "md_ChangeLog.html#autotoc_md67", null ],
        [ "3.8.1 (2009-10-30)", "md_ChangeLog.html#autotoc_md68", null ],
        [ "3.8.0 (2009-10-30)", "md_ChangeLog.html#autotoc_md69", null ]
      ] ],
      [ "3.7", "md_ChangeLog.html#autotoc_md70", [
        [ "3.7.1 (2009-10-06)", "md_ChangeLog.html#autotoc_md71", null ],
        [ "3.7.0 (2009-07-16)", "md_ChangeLog.html#autotoc_md72", null ]
      ] ],
      [ "3.6", "md_ChangeLog.html#autotoc_md73", [
        [ "3.6.1 (2009-05-01)", "md_ChangeLog.html#autotoc_md74", null ],
        [ "3.6.0 (2009-04-29)", "md_ChangeLog.html#autotoc_md75", null ]
      ] ],
      [ "3.5", "md_ChangeLog.html#autotoc_md76", [
        [ "3.5.5 (2009-01-13)", "md_ChangeLog.html#autotoc_md77", null ],
        [ "3.5.4 (2008-09-29)", "md_ChangeLog.html#autotoc_md78", null ],
        [ "3.5.3 (2008-09-19)", "md_ChangeLog.html#autotoc_md79", null ],
        [ "3.5.2 (2008-09-18)", "md_ChangeLog.html#autotoc_md80", null ],
        [ "3.5.1 (2008-08-26)", "md_ChangeLog.html#autotoc_md81", null ],
        [ "3.5.0 (2008-08-25)", "md_ChangeLog.html#autotoc_md82", null ]
      ] ],
      [ "3.4", "md_ChangeLog.html#autotoc_md83", [
        [ "3.4.2 (2008-07-18)", "md_ChangeLog.html#autotoc_md84", null ],
        [ "3.4.1 (2008-05-30)", "md_ChangeLog.html#autotoc_md85", null ],
        [ "3.4.0 (2008-04-25)", "md_ChangeLog.html#autotoc_md86", null ]
      ] ],
      [ "3.3", "md_ChangeLog.html#autotoc_md87", [
        [ "3.3.5 (2008-02-28)", "md_ChangeLog.html#autotoc_md88", null ],
        [ "3.3.4 (2007-12-27)", "md_ChangeLog.html#autotoc_md89", null ],
        [ "3.3.3 (2007-09-25)", "md_ChangeLog.html#autotoc_md90", null ],
        [ "3.3.1 (2007-06-20)", "md_ChangeLog.html#autotoc_md91", null ]
      ] ],
      [ "3.2", "md_ChangeLog.html#autotoc_md92", [
        [ "3.2.4 (2007-04-24)", "md_ChangeLog.html#autotoc_md93", null ],
        [ "3.2.3 (2006-11-29)", "md_ChangeLog.html#autotoc_md94", null ],
        [ "3.2 r795 (2006-10-11)", "md_ChangeLog.html#autotoc_md95", null ],
        [ "3.2.1 (2006-07-14) - dev release number 764", "md_ChangeLog.html#autotoc_md96", null ],
        [ "3.2.0 (2006-07-07) - dev release number 757", "md_ChangeLog.html#autotoc_md97", null ]
      ] ],
      [ "3.1", "md_ChangeLog.html#autotoc_md98", [
        [ "3.1.0 (2006-04-08) - dev release number 714", "md_ChangeLog.html#autotoc_md99", null ]
      ] ],
      [ "3.0", "md_ChangeLog.html#autotoc_md100", [
        [ "3.0.1 (2005-12-04)", "md_ChangeLog.html#autotoc_md101", null ],
        [ "3.0.0 (2005-08-26) - dev release number 510", "md_ChangeLog.html#autotoc_md102", null ]
      ] ],
      [ "2.2", "md_ChangeLog.html#autotoc_md103", [
        [ "no new release (2005-08-19)", "md_ChangeLog.html#autotoc_md104", null ],
        [ "2.2.1e (2005-05-30)", "md_ChangeLog.html#autotoc_md105", null ],
        [ "no new release (2005-01-07)", "md_ChangeLog.html#autotoc_md106", null ],
        [ "2.2.1d (2004-10-05)", "md_ChangeLog.html#autotoc_md107", null ],
        [ "2.2.1c (2004-07-20)", "md_ChangeLog.html#autotoc_md108", null ],
        [ "2.2.1b (2004-05-21)", "md_ChangeLog.html#autotoc_md109", null ],
        [ "2.2.1a (2004-05-13)", "md_ChangeLog.html#autotoc_md110", null ],
        [ "no new release (2004-04-28)", "md_ChangeLog.html#autotoc_md111", null ],
        [ "2.2.1 (2004-04-25)", "md_ChangeLog.html#autotoc_md112", null ],
        [ "2.2.0 (2004-03-10)", "md_ChangeLog.html#autotoc_md113", null ]
      ] ]
    ] ],
    [ "Authors and Contributors", "AUTHORS.html", null ],
    [ "License", "LICENSE.html", null ]
];
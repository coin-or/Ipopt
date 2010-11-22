# Function that takes as input the directory where Ipopt was built
# and the directory where the source code of ipoptr is located. 
# This function reads a Makefile from the Ipopt examples that has
# been configured for your system and uses this to create a Makevars
# file to compile ipoptr.
#
# example usage:
# On windows:
# configure_ipoptr( ipopt.build.dir = "C:/tools/CoinIpopt/build", 
#                   ipoptr.dir = "C:/work/projects/Rpackages/ipoptr" )
#
# On Linux:
# configure_ipoptr( ipopt.build.dir = "/cvos/shared/apps/Ipopt/3.5.5-source",
#                   ipoptr.dir = "~/Ipopt/contrib/RInterface" )

configure_ipoptr <- function( ipopt.build.dir, ipoptr.dir ) {

    inputfile <- readLines( paste( ipopt.build.dir, '/Ipopt/examples/hs071_cpp/Makefile', sep='' ) )

    # replace these variables to agree with Makevars.in
    inputfile <- gsub( "IPOPTINCDIR", "IPOPT_INCDIR", inputfile )
    inputfile <- gsub( "IPOPTLIBDIR", "IPOPT_LIBDIR", inputfile )

    # grep returns line number n
    # then we want the n-th element of inputfile
    # in this line we want to replace our search string with empty space to be left with the flags
    cxx.str <- "^CXX[ ]*="
    cxx = gsub( cxx.str, '', inputfile[ grep( cxx.str , inputfile ) ] )

    cxxflags.str <- "^CXXFLAGS[ ]*="
    cxxflags = gsub( cxxflags.str, '', inputfile[ grep( cxxflags.str , inputfile ) ] )

    cxxlinkflags.str <- "^CXXLINKFLAGS[ ]*="
    cxxlinkflags = gsub( cxxlinkflags.str, '', inputfile[ grep( cxxlinkflags.str , inputfile ) ] )

    # extract prefix
    prefix.str <- "^prefix[ ]*="
    prefix = gsub( prefix.str, '', inputfile[ grep( prefix.str , inputfile ) ] )
    if ( .Platform$OS.type == "windows" ) {
        # if we are on Windows, change if it begins with /t -> t:, /c -> c:, etc.
        prefix <- gsub( "^ /([a-z])/", " \\1:/", prefix )
    }

    # extract exec_prefix
    exec_prefix.str <- "^exec_prefix[ ]*="
    exec_prefix = gsub( exec_prefix.str, '', inputfile[ grep( exec_prefix.str , inputfile ) ] )

    ipopt_incdir.str <- "^IPOPT_INCDIR[ ]*="
    ipopt_incdir = gsub( ipopt_incdir.str, '', inputfile[ grep( ipopt_incdir.str , inputfile ) ] )

    ipopt_libdir.str <- "^IPOPT_LIBDIR[ ]*="
    ipopt_libdir = gsub( ipopt_libdir.str, '', inputfile[ grep( ipopt_libdir.str , inputfile ) ] )

    libs.str <- "^LIBS[ ]*="
    libs = gsub( libs.str, '', inputfile[ grep( libs.str , inputfile ) ] )

    incl.str <- "^INCL[ ]*="
    incl = gsub( incl.str, '', inputfile[ grep( incl.str , inputfile ) ] )

    cygpath.str <- "^CYGPATH_W[ ]*="
    cygpath = gsub( cygpath.str, '', inputfile[ grep( cygpath.str , inputfile ) ] )


    # load Makevars.in and replace the variables there
    outputfile <- readLines( paste( ipoptr.dir, "/src/Makevars.in", sep='' ) )

    outputfile <- gsub( "@CXX@", cxx, outputfile )
    outputfile <- gsub( "@CXXFLAGS@", cxxflags, outputfile )
    outputfile <- gsub( "@CXXLINKFLAGS@", cxxlinkflags, outputfile )
    outputfile <- gsub( "@includedir@", ipopt_incdir, outputfile )
    outputfile <- gsub( "@libdir@", ipopt_libdir, outputfile )
    outputfile <- gsub( "@exec_prefix@", exec_prefix, outputfile )
    outputfile <- gsub( "@prefix@", prefix, outputfile )
    outputfile <- gsub( "@ipoptlib@", libs, outputfile )
    outputfile <- gsub( "@ipoptinc@", incl, outputfile )
    outputfile <- gsub( "@CYGPATH_W@", cygpath, outputfile )

    # save as Makevars or Makevars.win
    if ( .Platform$OS.type == "windows" ) {
        writeLines( outputfile, paste( ipoptr.dir, "/src/Makevars.win", sep='' ) )
    } else {
        writeLines( outputfile, paste( ipoptr.dir, "/src/Makevars", sep='' ) )
    }

}

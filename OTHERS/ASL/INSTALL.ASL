Installation of the Ampl Solver Library (ASL) for IPOPT

**********************************************************************
Note: It is YOUR RESPONSIBILITY to ensure that you are entitled to
      download and use this third party package.
**********************************************************************

If you are working under Unix/Linux-like system and you have `wget'
installed, you can downloaded the required source code for the ASL
by executing the included script `get.ASL' in this directory.

Otherwise, follow those instructions:

1. Obtain the ASL source code, for example from
   http://www.netlib.org/ampl/solvers/ .
   Check http://www.netlib.org/ampl/solvers/README for instructions on
   how to obtain the entire directory 'solvers' in a convenient manner.

   For example:

   In this directory enter:

   ftp www.netlib.org

   and log in as 'anonymous' and enter your email address as password

   In FTP, enter

   ftp> cd ampl
   ftp> bin
   ftp> get solvers.tar
   ftp> quit

2. Unpack this package, so that the directory solvers is now located
   in $SRCROOTDIR/OTHERS/ASL , i.e. you have a directory
   $SRCROOTDIR/OTHERS/ASL/solvers containing the source code for the
   ASL.

   For example:

   If you obtained the file 'solvers.tar' as described above and it is
   in this directory, enter

   tar xf solvers.tar

   Since the files are compressed when obtained as described above,
   you need to uncompress them, e.g. by typing

   gunzip -fr solvers

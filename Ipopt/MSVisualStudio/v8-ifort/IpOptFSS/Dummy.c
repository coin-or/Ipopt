// Copyright (C) 2011 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Authors:  Marcel Roelofs     Paragon Decision Technology B.V.    2011-03-14

// Some trickery to get MSVC to export symbols from the system libraries and 
// libraries that are generated in dependent projects. Note that dummy_call is 
// exported nor called, but is only present as a trick for the linker to make the 
// the symbols to be exported available for exporting. Sigh.

// Correct prototypes not needed, so use the most simple ones. 
// Only necessary to get the linker to pull the symbols in.

extern int MA28PART();
extern int MPI_Finalize();
extern int MPI_Init();
extern int MPI_Comm_rank();
extern int dmumps_c();
extern int DAXPY();
extern int DDOT();
extern int DSCAL();
extern int DTRSM();
extern int DSYMV();
extern int DASUM();
extern int DGEMM();
extern int DSYRK();
extern int DNRM2();
extern int IDAMAX();
extern int DCOPY();
extern int DGEMV();
extern int DPOTRS();
extern int DGETRS();
extern int DGETRF();
extern int DPOTRF();
extern int DSYEV();

void dummy_call()
{
	MA28PART();
	MPI_Finalize();
	MPI_Init();
	MPI_Comm_rank();
	dmumps_c();
	DAXPY();
	DDOT();
	DSCAL();
	DTRSM();
	DSYMV();
	DASUM();
	DGEMM();
	DSYRK();
	DNRM2();
	IDAMAX();
	DCOPY();
	DGEMV();
	DPOTRS();
	DGETRS();
	DGETRF();
	DPOTRF();
	DSYEV();
}

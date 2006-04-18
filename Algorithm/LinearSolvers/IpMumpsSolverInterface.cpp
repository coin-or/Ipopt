// Copyright (C) 2006 Damien Hocking, KBC Advanced Technologies
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Authors: Damien Hocking                 KBC    2006-03-20
//        (included his original contribution into Ipopt package on 2006-03-25)
//          Andreas Waechter               IBM    2006-03-25
//           (minor changes and corrections)

#include "IpMumpsSolverInterface.hpp"

extern "C"
{
#include "mpi.h"
}

namespace Ipopt
{
#ifdef IP_DEBUG
  static const Index dbg_verbosity = 0;
#endif

#define JOB_INIT -1
#define JOB_END -2
#define USE_COMM_WORLD -987654

  MumpsSolverInterface::MumpsSolverInterface()
      :
      n(0),
      nz(0),
      a(NULL),
      irn_(NULL),
      jcn_(NULL),
      negevals(-1)
  {
    DBG_START_METH("MumpsSolverInterface::MumpsSolverInterface()", dbg_verbosity);
    mumps_data.a = 0;
    mumps_data.n = 0;
    mumps_data.nz = 0;
  }


  MumpsSolverInterface::~MumpsSolverInterface()
  {
    DBG_START_METH("MumpsSolverInterface::~MumpsSolverInterface()", dbg_verbosity);
    mumps_data.job = JOB_END;
    dmumps_c(&mumps_data); /* Terminate instance */
    MPI_Finalize();

    delete [] a;
    delete [] irn_;
    delete [] jcn_;
  }

  void MumpsSolverInterface::RegisterOptions(SmartPtr<RegisteredOptions> roptions)
  {}

  bool MumpsSolverInterface::InitializeImpl(const OptionsList& options,
      const std::string& prefix)
  {
    return true;
  }

  ESymSolverStatus MumpsSolverInterface::MultiSolve(bool new_matrix, const Index* ia, const Index* ja,
      Index nrhs, double* rhs_vals, bool check_NegEVals,
      Index numberOfNegEVals)
  {
    DBG_START_METH("MumpsSolverInterface::MultiSolve", dbg_verbosity);
    DBG_ASSERT(!check_NegEVals || ProvidesInertia());
    //DBG_ASSERT(initialized_);

    // check if a factorization has to be done
    // perform the factorization
    if (new_matrix) {
      ESymSolverStatus retval = Factorization(ia, ja, check_NegEVals, numberOfNegEVals);
      if (retval != SYMSOLVER_SUCCESS)  {
        DBG_PRINT((1, "FACTORIZATION FAILED!\n"));
        return retval;  // Matrix singular or error occurred
      }
    }
    // do the solve
    return Solve(ia, ja, nrhs, rhs_vals);
  }


  double* MumpsSolverInterface::GetValuesArrayPtr()
  {
    return a;
  }

  /** Initialize the local copy of the positions of the nonzero
      elements */
  ESymSolverStatus MumpsSolverInterface::InitializeStructure(Index dim, Index nonzeros,
      const Index* ia, const Index* ja)
  {
    ESymSolverStatus retval = SYMSOLVER_SUCCESS;
    DBG_START_METH("MumpsSolverInterface::InitializeStructure", dbg_verbosity);
    if (a) {
      delete [] a;
    }
    n = dim;
    nz = nonzeros;
    delete [] a;
    delete [] irn_;
    delete [] jcn_;
    a = new double[nz];
    irn_ = new int[nz];
    jcn_ = new int[nz];
    for (Index i=0; i<nz; i++) {
      irn_[i] = ia[i];
      jcn_[i] = ja[i];
    }

    int argc=1;
    //const char* name = "ipopt";
    char ** argv = 0;
    int myid, ierr;
    ierr = MPI_Init(&argc, &argv);
    ierr = MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    mumps_data.job = JOB_INIT;
    mumps_data.par = 1;//working host
    mumps_data.sym = 2;//general symmetric
    mumps_data.comm_fortran = USE_COMM_WORLD;
    dmumps_c(&mumps_data);

    mumps_data.a = a;
    mumps_data.irn = irn_;
    mumps_data.jcn = jcn_;
    mumps_data.n = n;
    mumps_data.nz = nz;

    return retval;
  }


  ESymSolverStatus MumpsSolverInterface::Factorization(const Index* ia, const Index* ja,
      bool check_NegEVals, Index numberOfNegEVals)
  {
    DBG_START_METH("MumpsSolverInterface::Factorization", dbg_verbosity);
    ESymSolverStatus retval = SYMSOLVER_SUCCESS;
    //Do symbolic pass
    mumps_data.job = 1;//symbolic ordering pass
    mumps_data.icntl[1] = 0;
    mumps_data.icntl[2] = 0;//QUIETLY!
    mumps_data.icntl[3] = 0;
    //mumps_data.icntl[5] = 0;//no column permutation
    //mumps_data.icntl[6] = 0;//AMD ordering
    mumps_data.icntl[7] = 2;//MC29 scaling
    mumps_data.icntl[9] = 3;//Iterative refinement iterations
    //mumps_data.icntl[13] = 1000.0;//Allowed % increase in workspace
    //mumps_data.cntl[0] = 0.0001;//Pivot tolerance
    mumps_data.cntl[0] = 0.01;//Pivot tolerance TODO: make option and flexible
    dmumps_c(&mumps_data);
    int error = mumps_data.info[0];
    if (error < 0) {
      Jnlst().Printf(J_ERROR, J_LINEAR_ALGEBRA,
                     "Error=%d returned from MUMPS in Factorization.\n",
                     error);
      retval = SYMSOLVER_FATAL_ERROR;
    }
    else {
      retval = SYMSOLVER_SUCCESS;
    }
    //int sizeMAXS = mumps_data.info[6];
    //int sizeMAXIS = mumps_data.info[7];
    //mumps_data.info[6] = 10*sizeMAXS;
    //mumps_data.info[7] = 10*sizeMAXIS;
    mumps_data.job = 2;//numerical factorisation
    int trycount = 0;
    while (trycount < 5) {
      dmumps_c(&mumps_data);
      error = mumps_data.info[0];
      if (error < 0) {
        if (error == -8 || error == -9)//not enough memory
        {
          Jnlst().Printf(J_WARNING, J_LINEAR_ALGEBRA,
                         "MUMPS requires more memory, reallocating.  Attempt %d\n",
                         trycount+1);
          double mem_percent = mumps_data.icntl[13];
          mumps_data.icntl[13] = (Index)(2.0 * mem_percent);
          trycount++;
          continue;
        }
      }
      else {
        break;
      }
    }

    if (trycount == 5) {
      Jnlst().Printf(J_ERROR, J_LINEAR_ALGEBRA,
                     "MUMPS was not able to obtain enough memory.\n");
      retval = SYMSOLVER_FATAL_ERROR;
    }
    else {
      negevals = mumps_data.infog[11];
      Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                     "MUMPS determined %d negative eigenvalues.\n",
                     negevals);
      if (check_NegEVals && (numberOfNegEVals!=negevals)) {
        retval = SYMSOLVER_WRONG_INERTIA;
      }
    }
    return retval;
  }

  ESymSolverStatus MumpsSolverInterface::Solve(const Index* ia, const Index* ja, Index nrhs, double *rhs_vals)
  {
    DBG_START_METH("MumpsSolverInterface::Solve", dbg_verbosity);
    ESymSolverStatus retval = SYMSOLVER_SUCCESS;
    IpData().TimingStats().LinearSystemBackSolve().Start();
    for (Index i = 0; i < nrhs; i++) {
      Index offset = i * mumps_data.n;
      mumps_data.rhs = &(rhs_vals[offset]);
      mumps_data.job = 3;//solve
      dmumps_c(&mumps_data);
      int error = mumps_data.info[0];
      if (error < 0) {
        Jnlst().Printf(J_ERROR, J_LINEAR_ALGEBRA,
                       "Error=%d returned from MUMPS in Solve.\n",
                       error);
        retval = SYMSOLVER_FATAL_ERROR;
      }
    }
    IpData().TimingStats().LinearSystemBackSolve().End();
    return retval;
  }

  Index MumpsSolverInterface::NumberOfNegEVals() const
  {
    DBG_START_METH("MumpsSolverInterface::NumberOfNegEVals", dbg_verbosity);
    DBG_ASSERT(negevals >= 0);
    return negevals;
  }

  bool MumpsSolverInterface::IncreaseQuality()
  {
    return false;
  }

}//end Ipopt namespace




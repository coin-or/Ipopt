// Copyright (C) 2011, Science and Technology Facilities Council.
// Copyright (C) 2009, Jonathan Hogg <jdh41.at.cantab.net>.
// Copyright (C) 2004, 2007 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// $Id$
//
// Authors: Jonathan Hogg                    STFC   2011-03-15
//          Jonathan Hogg                           2009-07-29
//          Carl Laird, Andreas Waechter     IBM    2004-03-17

#include "IpMa86SolverInterface.hpp"
#include <iostream>
#include <stdio.h>
#include <cmath>
using namespace std;

extern "C"
{
  /*
   * Easier to just have our own definition than include the full metis.h
   */
  extern void METIS_NodeND(int *n, int *xadj, int *adjncy,
                             int *numflag, int *options, int *perm, int *iperm);

  /*#include "hsl_mc64d.h"*/
}

namespace Ipopt
{

  Ma86SolverInterface::~Ma86SolverInterface()
  {
    delete [] val_;

    ma86_finalise(&keep_, &control_);
  }

  void Ma86SolverInterface::RegisterOptions(SmartPtr<RegisteredOptions> roptions)
  {
    roptions->AddIntegerOption(
      "ma86_print_level",
      "Debug printing level for the linear solver MA86",
      0,
      "Meep");
    /*
    "<0 no printing.\n"
    "0  Error and warning messages only.\n"
    "=1 Limited diagnostic printing.\n"
    ">1 Additional diagnostic printing.");
    */
    roptions->AddLowerBoundedIntegerOption(
      "ma86_nemin",
      "Node Amalgamation parameter",
      1, 32,
      "Two nodes in elimination tree are merged if result has fewer than "
      "ma86_nemin variables.");
    roptions->AddLowerBoundedNumberOption(
      "ma86_small",
      "Zero Pivot Threshold",
      0.0, false, 1e-20,
      "Any pivot less than ma86_small is treated as zero.");
    roptions->AddLowerBoundedNumberOption(
      "ma86_static",
      "Static Pivoting Threshold",
      0.0, false, 0.0,
      "See MA86 documentation. Either ma86_static=0.0 or "
      "ma86_static>ma86_small. ma86_static=0.0 disables static pivoting.");
    roptions->AddBoundedNumberOption(
      "ma86_u",
      "Pivoting Threshold",
      0.0, false, 0.5, false, 1e-8,
      "See MA86 documentation.");
    roptions->AddBoundedNumberOption(
      "ma86_umax",
      "Maximum Pivoting Threshold",
      0.0, false, 0.5, false, 1e-4,
      "See MA86 documentation.");
    /*roptions->AddStringOption2(
      "ma86_mc64_scaling",
      "Controls scaling of matrix using HSL_MC64",
      "yes",
      "no", "Do not scale the linear system matrix",
      "yes", "Scale the linear system matrix",
      "This option controls scaling for the solve HSL_MA86.");*/
  }

  bool Ma86SolverInterface::InitializeImpl(const OptionsList& options,
      const std::string& prefix)
  {
    ma86_default_control(&control_);
    control_.f_arrays = 1; // Use Fortran numbering (faster)

    options.GetIntegerValue("ma86_print_level", control_.diagnostics_level,
                            prefix);
    options.GetIntegerValue("ma86_nemin", control_.nemin, prefix);
    options.GetNumericValue("ma86_small", control_.small, prefix);
    options.GetNumericValue("ma86_static", control_.static_, prefix);
    options.GetNumericValue("ma86_u", control_.u, prefix);
    options.GetNumericValue("ma86_u", umax_, prefix);
    /*options.GetBoolValue("ma86_mc64_scaling", scale_, prefix);*/

    return true; // All is well
  }

  /*  Method for initializing internal stuctures.  Here, ndim gives
   *  the number of rows and columns of the matrix, nonzeros give
   *  the number of nonzero elements, and ia and ja give the
   *  positions of the nonzero elements, given in the matrix format
   *  determined by MatrixFormat.
   */
  ESymSolverStatus Ma86SolverInterface::InitializeStructure(Index dim,
      Index nonzeros, const Index* ia, const Index* ja)
  {
    struct ma86_info info;

    // Store size for later use
    ndim_ = dim;

    // Determine an ordering
    order_ = new Index[dim];
    MetisOrder(dim, ia, ja, order_);
    //for(int i=0; i<dim; i++) perm[i] = i+1;

    // Setup memory for values
    if (val_!=NULL) delete[] val_;
    val_ = new double[nonzeros];

    // Setup memory for scaling
    //if(scale_) scaling_ = new double[2*dim]; //size m+n for mc64 call

    if (HaveIpData()) {
      IpData().TimingStats().LinearSystemSymbolicFactorization().Start();
    }

    // perform analyse
    ma86_analyse(dim, ia, ja, order_, &keep_, &control_, &info);

    if (HaveIpData()) {
      IpData().TimingStats().LinearSystemSymbolicFactorization().End();
    }

    if (info.flag>=0) {
      return SYMSOLVER_SUCCESS;
    }
    else {
      return SYMSOLVER_FATAL_ERROR;
    }
  }

  /*  Solve operation for multiple right hand sides.  Solves the
   *  linear system A * x = b with multiple right hand sides, where
   *  A is the symmtric indefinite matrix.  Here, ia and ja give the
   *  positions of the values (in the required matrix data format).
   *  The actual values of the matrix will have been given to this
   *  object by copying them into the array provided by
   *  GetValuesArrayPtr. ia and ja are identical to the ones given
   *  to InitializeStructure.  The flag new_matrix is set to true,
   *  if the values of the matrix has changed, and a refactorzation
   *  is required.
   *
   *  The return code is SYMSOLV_SUCCESS if the factorization and
   *  solves were successful, SYMSOLV_SINGULAR if the linear system
   *  is singular, and SYMSOLV_WRONG_INERTIA if check_NegEVals is
   *  true and the number of negative eigenvalues in the matrix does
   *  not match numberOfNegEVals.  If SYMSOLV_CALL_AGAIN is
   *  returned, then the calling function will request the pointer
   *  for the array for storing a again (with GetValuesPtr), write
   *  the values of the nonzero elements into it, and call this
   *  MultiSolve method again with the same right-hand sides.  (This
   *  can be done, for example, if the linear solver realized it
   *  does not have sufficient memory and needs to redo the
   *  factorization; e.g., for MA27.)
   *
   *  The number of right-hand sides is given by nrhs, the values of
   *  the right-hand sides are given in rhs_vals (one full right-hand
   *  side stored immediately after the other), and solutions are
   *  to be returned in the same array.
   *
   *  check_NegEVals will not be chosen true, if ProvidesInertia()
   *  returns false.
   */
  ESymSolverStatus Ma86SolverInterface::MultiSolve(bool new_matrix,
      const Index* ia, const Index* ja, Index nrhs, double* rhs_vals,
      bool check_NegEVals, Index numberOfNegEVals)
  {
    struct ma86_info info;
    /*struct mc64_control control64;
    struct mc64_info info64;
    int *perm64;*/

    if (new_matrix || pivtol_changed_) {
      /*if(scale_) {
         if (HaveIpData()) {
            IpData().TimingStats().LinearSystemScaling().Start();
         }
         control64.checking=1; //disabled
         perm64 = new int[2*ndim_];
         mc64_matching(5, 4, ndim_, ndim_, ia, ja, val_, &control64, &info64,
            perm64, scaling_);
         delete[] perm64;
         for(int i=0; i<ndim_; i++) {
            scaling_[i] = exp(scaling_[i]);
            if(scaling_[i]>1e10) scaling_[i] = 1.0;
         }
         if (HaveIpData()) {
            IpData().TimingStats().LinearSystemScaling().End();
         }
      }*/

      if (HaveIpData()) {
        IpData().TimingStats().LinearSystemFactorization().Start();
      }
      //ma86_factor(ndim_, ia, ja, val_, order_, &keep_, &control_, &info);
      ma86_factor_solve(ndim_, ia, ja, val_, order_, &keep_, &control_, &info,
                        1, ndim_, rhs_vals, scaling_);
      if (HaveIpData()) {
        IpData().TimingStats().LinearSystemFactorization().End();
      }
      if (info.flag<0) return SYMSOLVER_FATAL_ERROR;
      if (check_NegEVals && info.num_neg!=numberOfNegEVals)
        return SYMSOLVER_WRONG_INERTIA;

      /*if (HaveIpData()) {
         IpData().TimingStats().LinearSystemBackSolve().Start();
      }
      ma86_solve(0, 1, ndim_, rhs_vals, order_, &keep_, &control_, &info);
      if (HaveIpData()) {
         IpData().TimingStats().LinearSystemBackSolve().End();
      }*/

      numneg_ = info.num_neg;
      pivtol_changed_ = false;
    }
    else {
      if (HaveIpData()) {
        IpData().TimingStats().LinearSystemBackSolve().Start();
      }
      ma86_solve(0, 1, ndim_, rhs_vals, order_, &keep_, &control_, &info,
                 scaling_);
      if (HaveIpData()) {
        IpData().TimingStats().LinearSystemBackSolve().End();
      }
    }

    return SYMSOLVER_SUCCESS;
  }

  /*
   * Call metis_NodeND to perform ordering on the graph, return it in perm
   */
  void Ma86SolverInterface::MetisOrder(const int ndim, const Index *ptr,
                                       const Index *row, Index *perm)
  {
    int options[8];
    options[0] = 0; // Defaults
    int numflag = 1;
    int ndim_nc = ndim;

    // Metis requires the full matrix, sans diagonal entries
    // As we only have the lower triangle we must expand
    Index *ptr_tmp = new Index[ndim+3]; // Two bigger than final ptr
    Index *row_tmp = new Index[2*(ptr[ndim]-1)];

    // First pass counts number of entries in each column
    // We do this at an offset of two so we can do everything in place
    for (int i=0; i<ndim+2; i++) ptr_tmp[i] = 0;
    for (int i=0; i<ndim; i++) {
      for (int j=ptr[i]-1; j<ptr[i+1]-1; j++) {
        if (i==row[j]-1) continue; // Skip diagonals
        ptr_tmp[i+2]++;
        ptr_tmp[row[j]-1+2]++;
      }
    }
    // Set ptr_tmp up so that ptr_tmp[i+1] is start of row i
    // This allows us to use ptr_tmp[i+1] as insert location for row i in
    // second pass of matrix data
    ptr_tmp[0] = 1;
    ptr_tmp[1]=1;
    for (int i=1;i<ndim;i++)
      ptr_tmp[i+1] = ptr_tmp[i] + ptr_tmp[i+1];
    // Second pass drops entries into correct locations
    ptr_tmp[0] = 1;
    for (int i=0; i<ndim; i++) {
      for (int j=ptr[i]-1; j<ptr[i+1]-1; j++) {
        if (i==row[j]-1) continue; // Skip diagonals
        int k = row[j]-1;
        row_tmp[ptr_tmp[i+1]-1] = row[j];
        row_tmp[ptr_tmp[k+1]-1] = i+1;
        ptr_tmp[i+1]++;
        ptr_tmp[k+1]++;
      }
    }

    // Note that MeTiS's iperm is our perm and vice-versa
    Index *iperm = new Index[ndim];
    //METIS_NodeND(&ndim_nc, ptr_tmp, row_tmp, &numflag, options, iperm, perm);
    THROW_EXCEPTION(INTERNAL_ABORT,
                    "Code in the MA86 currently needs to be changed.  The above line require Metis, but not all Ipopt builts have it.");
    delete[] iperm;
    delete[] row_tmp;
    delete[] ptr_tmp;
  }

  bool Ma86SolverInterface::IncreaseQuality()
  {
    if (control_.u >= umax_) {
      return false;
    }
    pivtol_changed_ = true;

    Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                   "Indreasing pivot tolerance for HSL_MA86 from %7.2e ",
                   control_.u);
    control_.u = Min(umax_, pow(control_.u,0.75));
    Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                   "to %7.2e.\n",
                   control_.u);
    return true;
  }

} // namespace Ipopt

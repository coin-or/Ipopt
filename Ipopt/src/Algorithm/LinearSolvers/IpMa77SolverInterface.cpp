// Copyright (C) 2013, Science and Technology Facilities Council.
// Copyright (C) 2009, Jonathan Hogg <jdh41.at.cantab.net>
// Copyright (C) 2004, 2007 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// $Id$
//
// Authors: Jonathan Hogg                    STFC   2012-02-14
//          Jonathan Hogg                           2009-07-29
//          Carl Laird, Andreas Waechter     IBM    2004-03-17

#include "IpoptConfig.h"

#ifdef COIN_HAS_HSL
#include "CoinHslConfig.h"
#endif

// if we do not have HSL_MA77 in HSL or the linear solver loader, then we want to build the MA77 interface
#if defined(COINHSL_HAS_MA77) || defined(HAVE_LINEARSOLVERLOADER)

#include "IpMa77SolverInterface.hpp"
#include <iostream>
#include <cmath>
using namespace std;

extern "C"
{
#include "hsl_mc68i.h"
}

namespace Ipopt
{

  Ma77SolverInterface::~Ma77SolverInterface()
  {
    delete [] val_;

    struct ma77_info info;
    if (keep_) ma77_finalise(&keep_, &control_, &info);
  }

  void Ma77SolverInterface::RegisterOptions(SmartPtr<RegisteredOptions> roptions)
  {
    roptions->AddIntegerOption(
      "ma77_print_level",
      "Debug printing level for the linear solver MA77",
      -1,
      "");
    /*
    "<0 no printing.\n"
    "0  Error and warning messages only.\n"
    "=1 Limited diagnostic printing.\n"
    ">1 Additional diagnostic printing.");
    */
    roptions->AddLowerBoundedIntegerOption(
      "ma77_buffer_lpage",
      "Number of scalars per MA77 buffer page",
      1, 4096,
      "Number of scalars per an in-core buffer in the out-of-core solver "
      "MA77. Must be at most ma77_file_size.");
    roptions->AddLowerBoundedIntegerOption(
      "ma77_buffer_npage",
      "Number of pages that make up MA77 buffer",
      1, 1600,
      "Number of pages of size buffer_lpage that exist in-core for the "
      "out-of-core solver MA77.");
    roptions->AddLowerBoundedIntegerOption(
      "ma77_file_size",
      "Target size of each temporary file for MA77, scalars per type",
      1, 2097152,
      "MA77 uses many temporary files, this option controls the size of "
      "each one. It is measured in the number of entries (int or double), "
      "NOT bytes.");
    roptions->AddLowerBoundedIntegerOption(
      "ma77_maxstore",
      "Maximum storage size for MA77 in-core mode",
      0, 0,
      "If greater than zero, the maximum size of factors stored in core "
      "before out-of-core mode is invoked.");
    roptions->AddLowerBoundedIntegerOption(
      "ma77_nemin",
      "Node Amalgamation parameter",
      1, 8,
      "Two nodes in elimination tree are merged if result has fewer than "
      "ma77_nemin variables.");
    roptions->AddLowerBoundedNumberOption(
      "ma77_small",
      "Zero Pivot Threshold",
      0.0, false, 1e-20,
      "Any pivot less than ma77_small is treated as zero.");
    roptions->AddLowerBoundedNumberOption(
      "ma77_static",
      "Static Pivoting Threshold",
      0.0, false, 0.0,
      "See MA77 documentation. Either ma77_static=0.0 or "
      "ma77_static>ma77_small. ma77_static=0.0 disables static pivoting.");
    roptions->AddBoundedNumberOption(
      "ma77_u",
      "Pivoting Threshold",
      0.0, false, 0.5, false, 1e-8,
      "See MA77 documentation.");
    roptions->AddBoundedNumberOption(
      "ma77_umax",
      "Maximum Pivoting Threshold",
      0.0, false, 0.5, false, 1e-4,
      "Maximum value to which u will be increased to improve quality.");
    roptions->AddStringOption2(
      "ma77_order",
      "Controls type of ordering used by HSL_MA77",
#ifdef COINHSL_HAS_METIS
      "metis",
#else
      "amd",
#endif
      "amd", "Use the HSL_MC68 approximate minimum degree algorithm",
      "metis", "Use the MeTiS nested dissection algorithm (if available)",
      "This option controls ordering for the solver HSL_MA77.");
  }

  bool Ma77SolverInterface::InitializeImpl(const OptionsList& options,
      const std::string& prefix)
  {
    ma77_default_control(&control_);
    control_.f_arrays = 1; // Use Fortran numbering (faster)
    control_.bits=32;
    // FIXME: HSL_MA77 should be updated to allow a matrix with new
    // values to be refactorized after a -11 (singular) error.
    //control_.action = 0; // false, should exit with error on singularity

    options.GetIntegerValue("ma77_print_level", control_.print_level, prefix);
    options.GetIntegerValue("ma77_buffer_lpage", control_.buffer_lpage[0], prefix);
    options.GetIntegerValue("ma77_buffer_lpage", control_.buffer_lpage[1], prefix);
    options.GetIntegerValue("ma77_buffer_npage", control_.buffer_npage[0], prefix);
    options.GetIntegerValue("ma77_buffer_npage", control_.buffer_npage[1], prefix);
    int temp;
    options.GetIntegerValue("ma77_file_size", temp, prefix);
    control_.file_size = temp;
    options.GetIntegerValue("ma77_maxstore", temp, prefix);
    control_.maxstore = temp;
    options.GetIntegerValue("ma77_nemin", control_.nemin, prefix);
    options.GetNumericValue("ma77_small", control_.small, prefix);
    options.GetNumericValue("ma77_static", control_.static_, prefix);
    options.GetNumericValue("ma77_u", control_.u, prefix);
    options.GetNumericValue("ma77_umax", umax_, prefix);
    std::string order_method;
    options.GetStringValue("ma77_order", order_method, prefix);
    if (order_method == "metis") {
      ordering_ = ORDER_METIS;
    } else {
      ordering_ = ORDER_AMD;
    }
    return true; // All is well
  }

  /*  Method for initializing internal structures.  Here, ndim gives
   *  the number of rows and columns of the matrix, nonzeros give
   *  the number of nonzero elements, and ia and ja give the
   *  positions of the nonzero elements, given in the matrix format
   *  determined by MatrixFormat.
   */
  ESymSolverStatus Ma77SolverInterface::InitializeStructure(Index dim, 
      Index nonzeros, const Index* ia, const Index* ja)
  {
    struct ma77_info info;
    struct mc68_control control68;
    struct mc68_info info68;

    // Store size for later use
    ndim_ = dim;

    if (HaveIpData()) {
      IpData().TimingStats().LinearSystemSymbolicFactorization().Start();
    }

    // mc68 requires a half matrix. A future version will support full
    // matrix entry, and this code should be removed when it is available.
    Index *ia_half = new Index[dim+1];
    Index *ja_half = new Index[ia[dim]-1];
    {
       int k = 0;
       for(int i=0; i<dim; i++) {
          ia_half[i] = k+1;
          for(int j=ia[i]-1; j<ia[i+1]-1; j++)
             if(ja[j]-1 >= i)
                ja_half[k++] = ja[j];
       }
       ia_half[dim] = k+1;
    }

    // Determine an ordering
    mc68_default_control(&control68);
    control68.f_array_in = 1; // Use Fortran numbering (faster)
    control68.f_array_out = 1; // Use Fortran numbering (faster)
    Index *perm = new Index[dim];
    if(ordering_ == ORDER_METIS) {
      mc68_order(3, dim, ia_half, ja_half, perm, &control68, &info68); /* MeTiS */
      if(info68.flag == -5) {
         // MeTiS not available
         ordering_ = ORDER_AMD;
      } else if(info68.flag<0) {
         delete[] ia_half;
         delete[] ja_half;
         return SYMSOLVER_FATAL_ERROR;
      }
    }
    if(ordering_ == ORDER_AMD) {
      mc68_order(1, dim, ia_half, ja_half, perm, &control68, &info68); /* AMD */
      if(info68.flag<0) {
         delete[] ia_half;
         delete[] ja_half;
         return SYMSOLVER_FATAL_ERROR;
      }
    }
    delete[] ia_half;
    delete[] ja_half;

    // Open files
    ma77_open(ndim_, "ma77_int", "ma77_real", "ma77_work", "ma77_delay", &keep_,
      &control_, &info);
    if (info.flag < 0) return SYMSOLVER_FATAL_ERROR;

    // Store data into files
    for(int i=0; i<dim; i++) {
      ma77_input_vars(i+1, ia[i+1]-ia[i], &(ja[ia[i]-1]), &keep_,
        &control_, &info);
      if (info.flag < 0) return SYMSOLVER_FATAL_ERROR;
    }

    // Perform analyse
    ma77_analyse(perm, &keep_, &control_, &info);
    delete[] perm; // Done with order

    if (HaveIpData())
      IpData().TimingStats().LinearSystemSymbolicFactorization().End();

    // Setup memory for values
    if (val_!=NULL) delete[] val_;
    val_ = new double[nonzeros];

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
  ESymSolverStatus Ma77SolverInterface::MultiSolve(bool new_matrix,
      const Index* ia, const Index* ja, Index nrhs, double* rhs_vals,
      bool check_NegEVals, Index numberOfNegEVals)
  {
    struct ma77_info info;

    if (new_matrix || pivtol_changed_)
    {
      for(int i=0; i<ndim_; i++) {
         ma77_input_reals(i+1, ia[i+1]-ia[i], &(val_[ia[i]-1]), &keep_,
            &control_, &info);
         if (info.flag < 0) return SYMSOLVER_FATAL_ERROR;
      }

      if (HaveIpData()) {
        IpData().TimingStats().LinearSystemFactorization().Start();
      }
      //ma77_factor(0, &keep_, &control_, &info, NULL);
      ma77_factor_solve(0, &keep_, &control_, &info, NULL,
                        nrhs, ndim_, rhs_vals);
      if (HaveIpData()) {
        IpData().TimingStats().LinearSystemFactorization().End();
      }
      if (info.flag==4 || info.flag==-11) return SYMSOLVER_SINGULAR;
      if (info.flag<0) return SYMSOLVER_FATAL_ERROR;
      if (check_NegEVals && info.num_neg!=numberOfNegEVals)
        return SYMSOLVER_WRONG_INERTIA;

      numneg_ = info.num_neg;
      pivtol_changed_ = false;
    }
    else {
      if (HaveIpData()) {
        IpData().TimingStats().LinearSystemBackSolve().Start();
      }
      ma77_solve(0, nrhs, ndim_, rhs_vals, &keep_, &control_, &info, NULL);
      if (HaveIpData()) {
        IpData().TimingStats().LinearSystemBackSolve().End();
      }
    }

    return SYMSOLVER_SUCCESS;
  }

  bool Ma77SolverInterface::IncreaseQuality()
  {
    if (control_.u >= umax_) {
      return false;
    }
    pivtol_changed_ = true;

    Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                   "Increasing pivot tolerance for HSL_MA77 from %7.2e ",
                   control_.u);
    control_.u = Min(umax_, pow(control_.u,0.75));
    Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                   "to %7.2e.\n",
                   control_.u);
    return true;
  }

} // namespace Ipopt

#endif /* COINHSL_HAS_MA77 or HAVE_LINEARSOLVERLOADER */

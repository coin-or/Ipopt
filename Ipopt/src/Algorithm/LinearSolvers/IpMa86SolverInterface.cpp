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

#include "IpoptConfig.h"

#ifdef COIN_HAS_HSL
#include "CoinHslConfig.h"
#endif

// if we do not have MA86 in HSL or the linear solver loader, then we want to build the MA86 interface
#if defined(COINHSL_HAS_MA86) || defined(HAVE_LINEARSOLVERLOADER)

#include "IpMa86SolverInterface.hpp"
#include <iostream>
#include <cmath>
using namespace std;

extern "C"
{
#include "hsl_mc68i.h"
}

namespace Ipopt
{

  Ma86SolverInterface::~Ma86SolverInterface()
  {
    delete [] val_;

    if(keep_) ma86_finalise(&keep_, &control_);
  }

  void Ma86SolverInterface::RegisterOptions(SmartPtr<RegisteredOptions> roptions)
  {
    roptions->AddIntegerOption(
      "ma86_print_level",
      "Debug printing level for the linear solver MA86",
      -1,
      "");
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
      "Maximum value to which u will be increased to improve quality.");
    roptions->AddStringOption3(
      "ma86_scaling",
      "Controls scaling of matrix",
      "mc64",
      "none", "Do not scale the linear system matrix",
      "mc64", "Scale linear system matrix using MC64",
      "mc77", "Scale linear system matrix using MC77 [1,3,0]",
      "This option controls scaling for the solver HSL_MA86.");
    roptions->AddStringOption3(
      "ma86_order",
      "Controls type of ordering used by HSL_MA86",
#ifdef COINHSL_HAS_METIS
      "auto",
#else
      "amd",
#endif
      "auto", "Try both AMD and MeTiS, pick best",
      "amd", "Use the HSL_MC68 approximate minimum degree algorithm",
      "metis", "Use the MeTiS nested dissection algorithm (if available)",
      "This option controls ordering for the solver HSL_MA86.");
  }

  bool Ma86SolverInterface::InitializeImpl(const OptionsList& options,
      const std::string& prefix)
  {
    ma86_default_control(&control_);
    control_.f_arrays = 1; // Use Fortran numbering (faster)
    /* Note: we can't set control_.action = false as we need to know the
     * intertia. (Otherwise we just enter the restoration phase and fail) */

    options.GetIntegerValue("ma86_print_level", control_.diagnostics_level,
                            prefix);
    options.GetIntegerValue("ma86_nemin", control_.nemin, prefix);
    options.GetNumericValue("ma86_small", control_.small_, prefix);
    options.GetNumericValue("ma86_static", control_.static_, prefix);
    options.GetNumericValue("ma86_u", control_.u, prefix);
    options.GetNumericValue("ma86_umax", umax_, prefix);
    std::string order_method, scaling_method;
    options.GetStringValue("ma86_order", order_method, prefix);
    if(order_method == "metis") {
      ordering_ = ORDER_METIS;
    } else if(order_method == "amd") {
      ordering_ = ORDER_AMD;
    } else {
      ordering_ = ORDER_AUTO;
    }
    options.GetStringValue("ma86_scaling", scaling_method, prefix);
    if(scaling_method == "mc64") {
      control_.scaling = 1;
    } else if(scaling_method == "mc77") {
      control_.scaling = 2;
    } else {
      control_.scaling = 0;
    }

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
    struct ma86_info info, info2;
    struct mc68_control control68;
    struct mc68_info info68;
    Index *order_amd, *order_metis;
    void *keep_amd, *keep_metis;

    // Store size for later use
    ndim_ = dim;

    // Determine an ordering
    mc68_default_control(&control68);
    control68.f_array_in = 1; // Use Fortran numbering (faster)
    control68.f_array_out = 1; // Use Fortran numbering (faster)
    order_amd = NULL; order_metis = NULL;
    if(ordering_ == ORDER_METIS || ordering_ == ORDER_AUTO) {
      order_metis = new Index[dim];
      mc68_order(3, dim, ia, ja, order_metis, &control68, &info68); /* MeTiS */
      if(info68.flag == -5) {
         // MeTiS not available
         ordering_ = ORDER_AMD;
         delete[] order_metis;
      } else if(info68.flag<0) {
         return SYMSOLVER_FATAL_ERROR;
      }
    }
    if(ordering_ == ORDER_AMD || ordering_ == ORDER_AUTO) {
      order_amd = new Index[dim];
      mc68_order(1, dim, ia, ja, order_amd, &control68, &info68); /* AMD */
    }
    if(info68.flag<0) return SYMSOLVER_FATAL_ERROR;

    if (HaveIpData()) {
      IpData().TimingStats().LinearSystemSymbolicFactorization().Start();
    }

    // perform analyse
    if(ordering_ == ORDER_AUTO) {
      ma86_analyse(dim, ia, ja, order_amd, &keep_amd, &control_, &info2);
      if (info2.flag<0) return SYMSOLVER_FATAL_ERROR;
      ma86_analyse(dim, ia, ja, order_metis, &keep_metis, &control_, &info);
      if (info.flag<0) return SYMSOLVER_FATAL_ERROR;
      if(info.num_flops > info2.num_flops) {
         // Use AMD
         //cout << "Choose AMD\n";
         order_ = order_amd;
         keep_ = keep_amd;
         delete[] order_metis;
         ma86_finalise(&keep_metis, &control_);
      } else {
         // Use MeTiS
         //cout << "Choose MeTiS\n";
         order_ = order_metis;
         keep_ = keep_metis;
         delete[] order_amd;
         ma86_finalise(&keep_amd, &control_);
      }
    } else {
      if(ordering_ == ORDER_AMD) order_ = order_amd;
      if(ordering_ == ORDER_METIS) order_ = order_metis;
      ma86_analyse(dim, ia, ja, order_, &keep_, &control_, &info);
    }

    if (HaveIpData()) {
      IpData().TimingStats().LinearSystemSymbolicFactorization().End();
    }

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
  ESymSolverStatus Ma86SolverInterface::MultiSolve(bool new_matrix,
      const Index* ia, const Index* ja, Index nrhs, double* rhs_vals,
      bool check_NegEVals, Index numberOfNegEVals)
  {
    struct ma86_info info;

    if (new_matrix || pivtol_changed_) {


      if (HaveIpData()) {
        IpData().TimingStats().LinearSystemFactorization().Start();
      }
      //ma86_factor(ndim_, ia, ja, val_, order_, &keep_, &control_, &info);
      ma86_factor_solve(ndim_, ia, ja, val_, order_, &keep_, &control_, &info,
                        nrhs, ndim_, rhs_vals, NULL);
      if (HaveIpData()) {
        IpData().TimingStats().LinearSystemFactorization().End();
      }
      if (info.flag<0) return SYMSOLVER_FATAL_ERROR;
      if (info.flag==2 || info.flag==-3) return SYMSOLVER_SINGULAR;
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
      ma86_solve(0, nrhs, ndim_, rhs_vals, order_, &keep_, &control_, &info,
                 NULL);
      if (HaveIpData()) {
        IpData().TimingStats().LinearSystemBackSolve().End();
      }
    }

    return SYMSOLVER_SUCCESS;
  }


  bool Ma86SolverInterface::IncreaseQuality()
  {
    if (control_.u >= umax_) {
      return false;
    }
    pivtol_changed_ = true;

    Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                   "Increasing pivot tolerance for HSL_MA86 from %7.2e ",
                   control_.u);
    control_.u = Min(umax_, pow(control_.u,0.75));
    Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                   "to %7.2e.\n",
                   control_.u);
    return true;
  }

} // namespace Ipopt

#endif /* COINHSL_HAS_MA86 or HAVE_LINEARSOLVERLOADER */

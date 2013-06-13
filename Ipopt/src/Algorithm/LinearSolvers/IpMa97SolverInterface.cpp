// Copyright (C) 2012, The Science and Technology Facilities Council.
// Copyright (C) 2009, Jonathan Hogg <jhogg41.at.gmail.com>.
// Copyright (C) 2004, 2007 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// $Id: IpMa97SolverInterface.cpp 2003 2011-06-03 03:53:44Z andreasw $
//
// Authors: Jonathan Hogg                    STFC   2012-12-21
//          Jonathan Hogg                           2009-07-29
//          Carl Laird, Andreas Waechter     IBM    2004-03-17

#include "IpoptConfig.h"

#ifdef COIN_HAS_HSL
#include "CoinHslConfig.h"
#endif

// if we have MA97 in HSL or the linear solver loader, then we want to build the MA97 interface
#if defined(COINHSL_HAS_MA97) || defined(HAVE_LINEARSOLVERLOADER)

#include "IpMa97SolverInterface.hpp"
#include <iostream>
#include <stdio.h>
#include <cmath>
using namespace std;

/* Uncomment the following line to enable the ma97_dump_matrix option.
 * This option requires a version of the coinhsl that supports this function.
 */
//#define MA97_DUMP_MATRIX

#ifdef MA97_DUMP_MATRIX
extern "C" {
  extern void F77_FUNC (dump_mat_csc, DUMP_MAT_CSC) (
      const ipfint    *factidx,
      const ipfint    *n,
      const ipfint    *ptr,
      const ipfint    *row,
      const double    *a);
}
#endif

namespace Ipopt
{

  Ma97SolverInterface::~Ma97SolverInterface()
  {
    delete [] val_;
    if(scaling_) delete [] scaling_;

    ma97_finalise(&akeep_, &fkeep_);
  }

  void Ma97SolverInterface::RegisterOptions(SmartPtr<RegisteredOptions> roptions)
  {
    roptions->AddIntegerOption(
      "ma97_print_level",
      "Debug printing level for the linear solver MA97",
      0, ""
      /*
      "<0 no printing.\n"
      "0  Error and warning messages only.\n"
      "=1 Limited diagnostic printing.\n"
      ">1 Additional diagnostic printing."*/);
    roptions->AddLowerBoundedIntegerOption(
      "ma97_nemin",
      "Node Amalgamation parameter",
      1, 8,
      "Two nodes in elimination tree are merged if result has fewer than "
      "ma97_nemin variables.");
    roptions->AddLowerBoundedNumberOption(
      "ma97_small",
      "Zero Pivot Threshold",
      0.0, false, 1e-20,
      "Any pivot less than ma97_small is treated as zero.");
    roptions->AddBoundedNumberOption(
      "ma97_u",
      "Pivoting Threshold",
      0.0, false, 0.5, false, 1e-8,
      "See MA97 documentation.");
    roptions->AddBoundedNumberOption(
      "ma97_umax",
      "Maximum Pivoting Threshold",
      0.0, false, 0.5, false, 1e-4,
      "See MA97 documentation.");
    roptions->AddStringOption5(
      "ma97_scaling",
      "Specifies strategy for scaling in HSL_MA97 linear solver",
      "dynamic",
      "none", "Do not scale the linear system matrix",
      "mc30", "Scale all linear system matrices using MC30",
      "mc64", "Scale all linear system matrices using MC64",
      "mc77", "Scale all linear system matrices using MC77 [1,3,0]",
      "dynamic", "Dynamically select scaling according to rules specified by ma97_scalingX and ma97_switchX options.",
      "");
    roptions->AddStringOption4(
      "ma97_scaling1",
      "First scaling.",
      "mc64",
      "none", "No scaling",
      "mc30", "Scale linear system matrix using MC30",
      "mc64", "Scale linear system matrix using MC64",
      "mc77", "Scale linear system matrix using MC77 [1,3,0]",
      "If ma97_scaling=dynamic, this scaling is used according to the trigger "
      "ma97_switch1. If ma97_switch2 is triggered it is disabled.");
    roptions->AddStringOption9(
      "ma97_switch1",
      "First switch, determine when ma97_scaling1 is enabled.",
      "od_hd_reuse",
      "never", "Scaling is never enabled.",
      "at_start", "Scaling to be used from the very start.",
      "at_start_reuse", "Scaling to be used on first iteration, then reused thereafter.",
      "on_demand", "Scaling to be used after Ipopt request improved solution (i.e. iterative refinement has failed).",
      "on_demand_reuse", "As on_demand, but reuse scaling from previous itr",
      "high_delay", "Scaling to be used after more than 0.05*n delays are present",
      "high_delay_reuse", "Scaling to be used only when previous itr created more that 0.05*n additional delays, otherwise reuse scaling from previous itr",
      "od_hd", "Combination of on_demand and high_delay",
      "od_hd_reuse", "Combination of on_demand_reuse and high_delay_reuse",
      "If ma97_scaling=dynamic, ma97_scaling1 is enabled according to this condition. If ma97_switch2 occurs this option is henceforth ignored.");
    roptions->AddStringOption4(
      "ma97_scaling2",
      "Second scaling.",
      "mc64",
      "none", "No scaling",
      "mc30", "Scale linear system matrix using MC30",
      "mc64", "Scale linear system matrix using MC64",
      "mc77", "Scale linear system matrix using MC77 [1,3,0]",
      "If ma97_scaling=dynamic, this scaling is used according to the trigger "
      "ma97_switch2. If ma97_switch3 is triggered it is disabled.");
    roptions->AddStringOption9(
      "ma97_switch2",
      "Second switch, determine when ma97_scaling2 is enabled.",
      "never",
      "never", "Scaling is never enabled.",
      "at_start", "Scaling to be used from the very start.",
      "at_start_reuse", "Scaling to be used on first iteration, then reused thereafter.",
      "on_demand", "Scaling to be used after Ipopt request improved solution (i.e. iterative refinement has failed).",
      "on_demand_reuse", "As on_demand, but reuse scaling from previous itr",
      "high_delay", "Scaling to be used after more than 0.05*n delays are present",
      "high_delay_reuse", "Scaling to be used only when previous itr created more that 0.05*n additional delays, otherwise reuse scaling from previous itr",
      "od_hd", "Combination of on_demand and high_delay",
      "od_hd_reuse", "Combination of on_demand_reuse and high_delay_reuse",
      "If ma97_scaling=dynamic, ma97_scaling2 is enabled according to this condition. If ma97_switch3 occurs this option is henceforth ignored.");
    roptions->AddStringOption4(
      "ma97_scaling3",
      "Third scaling.",
      "mc64",
      "none", "No scaling",
      "mc30", "Scale linear system matrix using MC30",
      "mc64", "Scale linear system matrix using MC64",
      "mc77", "Scale linear system matrix using MC77 [1,3,0]",
      "If ma97_scaling=dynamic, this scaling is used according to the trigger "
      "ma97_switch3.");
    roptions->AddStringOption9(
      "ma97_switch3",
      "Third switch, determine when ma97_scaling3 is enabled.",
      "never",
      "never", "Scaling is never enabled.",
      "at_start", "Scaling to be used from the very start.",
      "at_start_reuse", "Scaling to be used on first iteration, then reused thereafter.",
      "on_demand", "Scaling to be used after Ipopt request improved solution (i.e. iterative refinement has failed).",
      "on_demand_reuse", "As on_demand, but reuse scaling from previous itr",
      "high_delay", "Scaling to be used after more than 0.05*n delays are present",
      "high_delay_reuse", "Scaling to be used only when previous itr created more that 0.05*n additional delays, otherwise reuse scaling from previous itr",
      "od_hd", "Combination of on_demand and high_delay",
      "od_hd_reuse", "Combination of on_demand_reuse and high_delay_reuse",
      "If ma97_scaling=dynamic, ma97_scaling3 is enabled according to this condition.");
    roptions->AddStringOption7(
      "ma97_order",
      "Controls type of ordering used by HSL_MA97",
      "auto",
      "auto", "Use HSL_MA97 heuristic to guess best of AMD and METIS",
      "best", "Try both AMD and MeTiS, pick best",
      "amd", "Use the HSL_MC68 approximate minimum degree algorithm",
      "metis", "Use the MeTiS nested dissection algorithm",
      "matched-auto", "Use the HSL_MC80 matching with heuristic choice of AMD or METIS",
      "matched-metis", "Use the HSL_MC80 matching based ordering with METIS",
      "matched-amd", "Use the HSL_MC80 matching based ordering with AMD",
      "");
#ifdef MA97_DUMP_MATRIX
    roptions->AddStringOption2(
      "ma97_dump_matrix",
      "Controls whether HSL_MA97 dumps each matrix to a file",
      "no",
      "no", "Do not dump matrix",
      "yes", "Do dump matrix",
      "");
#endif
    roptions->AddStringOption2(
      "ma97_solve_blas3",
      "Controls if blas2 or blas3 routines are used for solve",
      "no",
      "no", "Use BLAS2 (faster, some implementations bit incompatible)",
      "yes", "Use BLAS3 (slower)",
      "");
  }

  int Ma97SolverInterface::ScaleNameToNum(const std::string& name) {
     if(name=="none") return 0;
     if(name=="mc64") return 1;
     if(name=="mc77") return 2;
     if(name=="mc30") return 4;
     assert(0);
     return -1;
  }

  bool Ma97SolverInterface::InitializeImpl(const OptionsList& options,
      const std::string& prefix)
  {
    ma97_default_control(&control_);
    control_.f_arrays = 1; // Use Fortran numbering (faster)
    control_.action = 0; // false, shuold exit with error on singularity

    options.GetIntegerValue("ma97_print_level", control_.print_level,
                            prefix);
    options.GetIntegerValue("ma97_nemin", control_.nemin, prefix);
    options.GetNumericValue("ma97_small", control_.small, prefix);
    options.GetNumericValue("ma97_u", control_.u, prefix);
    options.GetNumericValue("ma97_umax", umax_, prefix);
    std::string order_method, scaling_method, rescale_strategy;
    options.GetStringValue("ma97_order", order_method, prefix);
    if(order_method == "metis") {
      ordering_ = ORDER_METIS;
    } else if(order_method == "amd") {
      ordering_ = ORDER_AMD;
    } else if(order_method == "best") {
       ordering_ = ORDER_BEST;
    } else if(order_method == "matched-metis") {
       ordering_ = ORDER_MATCHED_METIS;
    } else if(order_method == "matched-amd") {
       ordering_ = ORDER_MATCHED_AMD;
    } else if(order_method == "matched-auto") {
       ordering_ = ORDER_MATCHED_AUTO;
    } else {
      ordering_ = ORDER_AUTO;
    }
    options.GetStringValue("ma97_scaling", scaling_method, prefix);
    current_level_ = 0;
    if(scaling_method == "dynamic") {
       scaling_type_ = 0;
       string switch_val[3], scale_val[3];
       options.GetStringValue("ma97_switch1", switch_val[0], prefix);
       options.GetStringValue("ma97_scaling1", scale_val[0], prefix);
       options.GetStringValue("ma97_switch2", switch_val[1], prefix);
       options.GetStringValue("ma97_scaling2", scale_val[1], prefix);
       options.GetStringValue("ma97_switch3", switch_val[2], prefix);
       options.GetStringValue("ma97_scaling3", scale_val[2], prefix);
       for(int i=0; i<3; i++) {
          scaling_val_[i] = ScaleNameToNum(scale_val[i]);
          if(switch_val[i]=="never")
             switch_[i] = SWITCH_NEVER;
          else if(switch_val[i]=="at_start") {
             switch_[i] = SWITCH_AT_START;
             scaling_type_ = scaling_val_[i];
             current_level_ = i;
             Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
               "HSL_MA97: Enabled scaling level %d on initialization\n", current_level_); 
          }
          else if(switch_val[i]=="at_start_reuse") {
             switch_[i] = SWITCH_AT_START_REUSE;
             scaling_type_ = scaling_val_[i];
             current_level_ = i;
             Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
               "HSL_MA97: Enabled scaling level %d on initialization\n", current_level_); 
          }
          else if(switch_val[i]=="on_demand")
             switch_[i] = SWITCH_ON_DEMAND;
          else if(switch_val[i]=="on_demand_reuse")
             switch_[i] = SWITCH_ON_DEMAND_REUSE;
          else if(switch_val[i]=="high_delay")
             switch_[i] = SWITCH_NDELAY;
          else if(switch_val[i]=="high_delay_reuse")
             switch_[i] = SWITCH_NDELAY_REUSE;
          else if(switch_val[i]=="od_hd")
             switch_[i] = SWITCH_OD_ND;
          else if(switch_val[i]=="od_hd_reuse")
             switch_[i] = SWITCH_OD_ND_REUSE;
       }
    } else {
       switch_[0] = SWITCH_AT_START;
       switch_[1] = SWITCH_NEVER;
       switch_[2] = SWITCH_NEVER;
       scaling_type_ = ScaleNameToNum(scaling_method);
    }
#ifdef MA97_DUMP_MATRIX
    options.GetBoolValue("ma97_dump_matrix", dump_, prefix);
#endif
    bool solve_blas3;
    options.GetBoolValue("ma97_solve_blas3", solve_blas3, prefix);
    control_.solve_blas3 = solve_blas3 ? 1 : 0;

    // Set whether we scale on first iteration or not
    switch(switch_[current_level_]) {
       case SWITCH_NEVER:
       case SWITCH_ON_DEMAND:
       case SWITCH_ON_DEMAND_REUSE:
       case SWITCH_NDELAY:
       case SWITCH_NDELAY_REUSE:
       case SWITCH_OD_ND:
       case SWITCH_OD_ND_REUSE:
          rescale_ = false;
          break;
       case SWITCH_AT_START:
       case SWITCH_AT_START_REUSE:
          rescale_ = true;
          break;
    }
    // Set scaling
    control_.scaling = scaling_type_;

    return true; // All is well
  }

  /*  Method for initializing internal stuctures.  Here, ndim gives
   *  the number of rows and columns of the matrix, nonzeros give
   *  the number of nonzero elements, and ia and ja give the
   *  positions of the nonzero elements, given in the matrix format
   *  determined by MatrixFormat.
   */
  ESymSolverStatus Ma97SolverInterface::InitializeStructure(Index dim,
      Index nonzeros, const Index* ia, const Index* ja)
  {
    struct ma97_info info, info2;
    void *akeep_amd, *akeep_metis;

    // Store size for later use
    ndim_ = dim;

    // Setup memory for values
    if (val_!=NULL) delete[] val_;
    val_ = new double[nonzeros];

    // Check if analyse needs to be postponed
    if(ordering_ == ORDER_MATCHED_AMD || ordering_ == ORDER_MATCHED_METIS) {
      // Ordering requires values. Just signal success and return
       Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                     "HSL_MA97: Delaying analyse until values are available\n"); 
      switch(ordering_)
      {
        case ORDER_MATCHED_AMD:
          control_.ordering = 7; // HSL_MC80 with AMD
          break;
        case ORDER_MATCHED_METIS:
          control_.ordering = 8; // HSL_MC80 with METIS
          break;
      }
       return SYMSOLVER_SUCCESS;
    }

    if (HaveIpData()) {
      IpData().TimingStats().LinearSystemSymbolicFactorization().Start();
    }

    // perform analyse
    if(ordering_ == ORDER_BEST) {
      Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                    "HSL_MA97: Use best of AMD or MeTiS:\n");
      control_.ordering = 1; // AMD
      ma97_analyse(0, dim, ia, ja, NULL, &akeep_amd, &control_, &info2, NULL);
      if (info2.flag<0) return SYMSOLVER_FATAL_ERROR;
      Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                    "AMD   nfactor = %d, nflops = %d:\n",
                    info2.num_factor, info2.num_flops);
      control_.ordering = 3; // METIS
      ma97_analyse(0, dim, ia, ja, NULL, &akeep_metis, &control_, &info, NULL);
      if (info.flag<0) return SYMSOLVER_FATAL_ERROR;
      Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                    "MeTiS nfactor = %d, nflops = %d:\n",
                    info.num_factor, info.num_flops);
      if(info.num_flops > info2.num_flops) {
        // Use AMD
        Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                      "HSL_MA97: Choose AMD\n");
        akeep_ = akeep_amd;
        ma97_free_akeep(&akeep_metis);
        info = info2;
      } else {
         // Use MeTiS
        Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                      "HSL_MA97: Choose MeTiS\n");
         akeep_ = akeep_metis;
         ma97_free_akeep(&akeep_amd);
      }
    } else {
      switch(ordering_)
      {
        case ORDER_AMD:
        case ORDER_MATCHED_AMD:
          control_.ordering = 1; // AMD
          break;
        case ORDER_METIS:
        case ORDER_MATCHED_METIS:
          control_.ordering = 3; // METIS
          break;
        case ORDER_AUTO:
        case ORDER_MATCHED_AUTO:
          Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
            "HSL_MA97: Make heuristic choice of AMD or MeTiS\n");
          control_.ordering = 5; // Use heuristic to pick which to use
      }
      ma97_analyse(0, dim, ia, ja, NULL, &akeep_, &control_, &info, NULL);
      switch(info.ordering) {
        case 1:
          Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                        "HSL_MA97: Used AMD\n");
          if(ordering_ == ORDER_MATCHED_AUTO)
             ordering_ = ORDER_MATCHED_AMD;
          break;
        case 3:
          Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                        "HSL_MA97: Used MeTiS\n");
          if(ordering_ == ORDER_MATCHED_AUTO)
             ordering_ = ORDER_MATCHED_METIS;
          break;
        default:
          Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                        "HSL_MA97: Used ordering %d\n", info.ordering);
          break;
      }
    }

    Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                  "HSL_MA97: PREDICTED nfactor %d, maxfront %d\n",
                  info.num_factor, info.maxfront); 

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
  ESymSolverStatus Ma97SolverInterface::MultiSolve(bool new_matrix,
      const Index* ia, const Index* ja, Index nrhs, double* rhs_vals,
      bool check_NegEVals, Index numberOfNegEVals)
  {
    struct ma97_info info;
    Number t1=0,t2;

    if (new_matrix || pivtol_changed_) {

#ifdef MA97_DUMP_MATRIX
      if(dump_) {
        Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                      "Dumping matrix %d\n", fctidx_);
        F77_FUNC (dump_mat_csc, DUMP_MAT_CSC)
        (&fctidx_, &ndim_, ia, ja, val_);
        fctidx_++;
      }
#endif

      // Set scaling option
      if(rescale_) {
         control_.scaling = scaling_type_;
         if(scaling_type_!=0 && scaling_==NULL)
            scaling_ = new double[ndim_]; // alloc if not already
      } else {
         control_.scaling = 0; // None or user (depends if scaling_ is alloc'd)
      }

      if((ordering_ == ORDER_MATCHED_AMD || ordering_ == ORDER_MATCHED_METIS) &&  rescale_) {
         /*
          * Perform delayed analyse
          */
         if (HaveIpData()) {
           IpData().TimingStats().LinearSystemSymbolicFactorization().Start();
         }

         switch(ordering_)
         {
           case ORDER_MATCHED_AMD:
             control_.ordering = 7; // HSL_MC80 with AMD
             break;
           case ORDER_MATCHED_METIS:
             control_.ordering = 8; // HSL_MC80 with METIS
             break;
         }

         ma97_analyse(0, ndim_, ia, ja, val_, &akeep_, &control_, &info, NULL);
         if(scaling_type_ == 1)
            control_.scaling = 3; // use mc64 from ordering

         Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                        "HSL_MA97: PREDICTED nfactor %d, maxfront %d\n",
                        info.num_factor, info.maxfront); 

         if (HaveIpData()) {
           IpData().TimingStats().LinearSystemSymbolicFactorization().End();
         }

         if (info.flag==6 || info.flag==-7) {
           Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                         "In Ma97SolverInterface::Factorization: "
                         "Singular system, estimated rank %d of %d\n",
                         info.matrix_rank, ndim_); 
            return SYMSOLVER_SINGULAR;
         }
         if (info.flag<0) return SYMSOLVER_FATAL_ERROR;

      }

      if (HaveIpData()) {
        t1 = IpData().TimingStats().LinearSystemFactorization().TotalWallclockTime();
        IpData().TimingStats().LinearSystemFactorization().Start();
      }
      ma97_factor(4, ia, ja, val_, &akeep_, &fkeep_, &control_, &info, scaling_);
      //ma97_factor_solve(4, ia, ja, val_, 1, rhs_vals, ndim_, &akeep_, &fkeep_,
      //                  &control_, &info, scaling_);
      Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                    "HSL_MA97: delays %d, nfactor %d, nflops %ld, maxfront %d\n",
                    info.num_delay, info.num_factor, info.num_flops, info.maxfront); 
      if (HaveIpData()) {
        IpData().TimingStats().LinearSystemFactorization().End();
        t2 = IpData().TimingStats().LinearSystemFactorization().TotalWallclockTime();
        Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                      "Ma97SolverInterface::Factorization: "
                      "ma97_factor_solve took %10.3f\n",
                      t2-t1); 
      }
      if (info.flag==7 || info.flag==-7) {
        Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                      "In Ma97SolverInterface::Factorization: "
                      "Singular system, estimated rank %d of %d\n",
                      info.matrix_rank, ndim_); 
         return SYMSOLVER_SINGULAR;
      }
      for(int i=current_level_; i<3; i++) {
         switch(switch_[i]) {
            case SWITCH_NEVER:
            case SWITCH_AT_START:
            case SWITCH_ON_DEMAND:
               // Nothing to do here
               break;
            case SWITCH_AT_START_REUSE:
               rescale_ = false; // Scaled exactly once, never changed again
               break;
            case SWITCH_ON_DEMAND_REUSE:
               if(i==current_level_ && rescale_) rescale_ = false;
               break;
            case SWITCH_NDELAY_REUSE:
            case SWITCH_OD_ND_REUSE:
               if(rescale_) numdelay_ = info.num_delay; // Need to do this before we reset rescale_
               if(i==current_level_ && rescale_) rescale_ = false;
               // Falls through to:
            case SWITCH_NDELAY:
            case SWITCH_OD_ND:
               if(rescale_) numdelay_ = info.num_delay;
               if(info.num_delay - numdelay_ > 0.05*ndim_) {
                  // number of delays has signficantly increased, so trigger
                  current_level_ = i;
                  scaling_type_ = scaling_val_[i];
                  Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                     "HSL_MA97: Enabling scaling %d due to excess delays\n", i); 
                  rescale_ = true;
               }
               break;
         }
      }
      if (info.flag<0) {
        Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                      "In Ma97SolverInterface::Factorization: "
                      "Unhandled error. info.flag = %d\n",
                      info.flag); 
         return SYMSOLVER_FATAL_ERROR;
      }
      if (check_NegEVals && info.num_neg!=numberOfNegEVals) {
        Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                      "In Ma97SolverInterface::Factorization: "
                      "info.num_neg = %d, but numberOfNegEVals = %d\n",
                      info.num_neg, numberOfNegEVals); 
        return SYMSOLVER_WRONG_INERTIA;
      }

      if (HaveIpData()) {
         IpData().TimingStats().LinearSystemBackSolve().Start();
      }
      ma97_solve(0, 1, rhs_vals, ndim_, &akeep_, &fkeep_, &control_, &info);
      if (HaveIpData()) {
         IpData().TimingStats().LinearSystemBackSolve().End();
      }

      numneg_ = info.num_neg;
      pivtol_changed_ = false;
    }
    else {
      if (HaveIpData()) {
        IpData().TimingStats().LinearSystemBackSolve().Start();
      }
      ma97_solve(0, 1, rhs_vals, ndim_, &akeep_, &fkeep_, &control_, &info);
      if (HaveIpData()) {
        IpData().TimingStats().LinearSystemBackSolve().End();
      }
    }

    return SYMSOLVER_SUCCESS;
  }


  bool Ma97SolverInterface::IncreaseQuality()
  {
    for(int i=current_level_; i<3; i++) {
      switch(switch_[i]) {
        case SWITCH_ON_DEMAND:
        case SWITCH_ON_DEMAND_REUSE:
        case SWITCH_OD_ND:
        case SWITCH_OD_ND_REUSE:
          rescale_ = true;
          current_level_ = i;
          scaling_type_ = scaling_val_[i];
          Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
            "HSL_MA97: Enabling scaling %d due to failure of iterative refinement\n", current_level_); 
          break;
        default: ;
      }
    }

    if (control_.u >= umax_) {
      return false;
    }
    pivtol_changed_ = true;
    Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                   "Indreasing pivot tolerance for HSL_MA97 from %7.2e ",
                   control_.u);
    control_.u = Min(umax_, pow(control_.u,0.75));
    Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                   "to %7.2e.\n",
                   control_.u);
    return true;
  }

} // namespace Ipopt

#endif /* COINHSL_HAS_MA97 or HAVE_LINEARSOLVERLOADER */

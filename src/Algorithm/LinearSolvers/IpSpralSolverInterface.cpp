// Copyright (C) 2012, The Science and Technology Facilities Council.
// Copyright (C) 2009, Jonathan Hogg <jhogg41.at.gmail.com>.
// Copyright (C) 2004, 2007 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Authors: Byron Tasseff                    LANL   2020-03-21
//          Jonathan Hogg                    STFC   2012-12-21
//          Jonathan Hogg                           2009-07-29
//          Carl Laird, Andreas Waechter     IBM    2004-03-17

#include "IpoptConfig.h"
#include "IpSpralSolverInterface.hpp"

#include <cassert>
#include <cmath>
#include <cinttypes>

using namespace std;

namespace Ipopt
{

SpralSolverInterface::~SpralSolverInterface()
{
   delete[] val_;
   delete[] scaling_;

   spral_ssids_free(&akeep_, &fkeep_);
}

void SpralSolverInterface::RegisterOptions(
   SmartPtr<RegisteredOptions> roptions
)
{
   roptions->AddLowerBoundedIntegerOption(
      "spral_cpu_block_size",
      "CPU Parallelization Block Size",
      1,
      256,
      "Block size to use for parallelization of large nodes on CPU resources.");

   roptions->AddLowerBoundedNumberOption(
      "spral_gpu_perf_coeff",
      "GPU Performance Coefficient",
      0.0, true,
      1.0,
      "How many times faster a GPU is than a CPU at factoring a subtree.");

   roptions->AddStringOption2(
      "spral_ignore_numa", "Non-uniform memory access (NUMA) region setting.",
      "yes", "no", "Do not treat CPUs and GPUs as belonging to a single NUMA region.",
      "yes", "Treat CPUs and GPUs as belonging to a single NUMA region.");

   roptions->AddLowerBoundedNumberOption(
      "spral_max_load_inbalance",
      "Maximum Permissible Load",
      1.0, true,
      1.2,
      "Maximum permissible load inbalance for leaf subtree allocations.");

   roptions->AddLowerBoundedNumberOption(
      "spral_min_gpu_work",
      "Minimum GPU Work",
      0.0, false,
      5.0e9,
      "Minimum number of FLOPS in subtree before scheduling on GPU.");

   roptions->AddLowerBoundedIntegerOption(
      "spral_nemin",
      "Node Amalgamation Parameter",
      1,
      32,
      "Two nodes in the elimination tree are merged if the result has fewer than spral_nemin variables.");

   roptions->AddStringOption2(
      "spral_order",
      "Controls type of ordering used by SPRAL",
      "matching",
      "metis", "Use METIS with default settings.",
      "matching", "Use matching-based elimination ordering.");

   roptions->AddStringOption3(
      "spral_pivot_method",
      "Specifies strategy for scaling in SPRAL linear solver.",
      "block",
      "aggressive", "Aggressive a posteori pivoting.",
      "block", "Block a posteori pivoting.",
      "threshold", "Threshold partial pivoting (not parallel).");

   roptions->AddIntegerOption(
      "spral_print_level",
      "Print level for the linear solver SPRAL",
      -1,
      "<0: no printing, 0: errors and warning messages, 1: limited diagnostics, >1: additional diagnostics");

   roptions->AddStringOption6(
      "spral_scaling",
      "Specifies strategy for scaling in SPRAL linear solver.",
      "matching",
      "none", "Do not scale the linear system matrix.",
      "mc64", "Scale using weighted bipartite matching (MC64).",
      "auction", "Scale using the auction algorithm.",
      "matching", "Scale using the matching-based ordering.",
      "ruiz", "Scale using the norm-equilibration algorithm of Ruiz (MC77).",
      "dynamic", "Dynamically select scaling according to switch options.");

   roptions->AddStringOption5(
      "spral_scaling_1",
      "First scaling strategy.",
      "matching",
      "none", "Do not scale the linear system matrix.",
      "mc64", "Scale using weighted bipartite matching (MC64).",
      "auction", "Scale using the auction algorithm.",
      "matching", "Scale using the matching-based ordering.",
      "ruiz", "Scale using the norm-equilibration algorithm of Ruiz (MC77).",
      "If spral_scaling = dynamic, this scaling is used according to the trigger "
      "spral_switch_1. If spral_switch_2 is triggered, it is disabled.",
      true);

   roptions->AddStringOption5(
      "spral_scaling_2",
      "Second scaling strategy.",
      "mc64",
      "none", "Do not scale the linear system matrix.",
      "mc64", "Scale using weighted bipartite matching (MC64).",
      "auction", "Scale using the auction algorithm.",
      "matching", "Scale using the matching-based ordering.",
      "ruiz", "Scale using the norm-equilibration algorithm of Ruiz (MC77).",
      "If spral_scaling = dynamic, this scaling is used according to the trigger "
      "spral_switch_2. If spral_switch_3 is triggered, it is disabled.",
      true);

   roptions->AddStringOption5(
      "spral_scaling_3",
      "Third scaling strategy.",
      "none",
      "none", "Do not scale the linear system matrix.",
      "mc64", "Scale using weighted bipartite matching (MC64).",
      "auction", "Scale using the auction algorithm.",
      "matching", "Scale using the matching-based ordering.",
      "ruiz", "Scale using the norm-equilibration algorithm of Ruiz (MC77).",
      "If spral_scaling = dynamic, this scaling is used according to the trigger spral_switch_3.",
      true);

   roptions->AddStringOption9(
      "spral_switch_1",
      "First switch, determining when spral_scaling_1 is enabled.",
      "at_start",
      "never", "Scaling is never enabled.",
      "at_start", "Scaling is used from the very start.",
      "at_start_reuse", "Scaling is used on the first iteration, then reused thereafter.",
      "on_demand", "Scaling is used when iterative refinement has failed.",
      "on_demand_reuse", "As on_demand, but scaling from previous iteration is reused.",
      "high_delay", "Scaling is used after more than 0.05*n delays are present.",
      "high_delay_reuse", "Scaling is used only when previous iteration created "
      "more that 0.05*n additional delays; otherwise, reuse scaling from the previous iteration.",
      "od_hd", "Combination of on_demand and high_delay.",
      "od_hd_reuse", "Combination of on_demand_reuse and high_delay_reuse",
      "If spral_scaling = dynamic, spral_scaling_1 is enabled according to this "
      "condition. If spral_switch_2 occurs, this option is henceforth ignored.",
      true);

   roptions->AddStringOption9(
      "spral_switch_2",
      "Second switch, determining when spral_scaling_2 is enabled.",
      "on_demand",
      "never", "Scaling is never enabled.",
      "at_start", "Scaling is used from the very start.",
      "at_start_reuse", "Scaling is used on the first iteration, then reused thereafter.",
      "on_demand", "Scaling is used when iterative refinement has failed.",
      "on_demand_reuse", "As on_demand, but scaling from previous iteration is reused.",
      "high_delay", "Scaling is used after more than 0.05*n delays are present.",
      "high_delay_reuse", "Scaling is used only when previous iteration created "
      "more that 0.05*n additional delays; otherwise, reuse scaling from the previous iteration.",
      "od_hd", "Combination of on_demand and high_delay.",
      "od_hd_reuse", "Combination of on_demand_reuse and high_delay_reuse",
      "If spral_scaling = dynamic, spral_scaling_2 is enabled according to this "
      "condition. If spral_switch_3 occurs, this option is henceforth ignored.",
      true);

   roptions->AddStringOption9(
      "spral_switch_3",
      "Third switch, determining when spral_scaling_3 is enabled.",
      "never",
      "never", "Scaling is never enabled.",
      "at_start", "Scaling is used from the very start.",
      "at_start_reuse", "Scaling is used on the first iteration, then reused thereafter.",
      "on_demand", "Scaling is used when iterative refinement has failed.",
      "on_demand_reuse", "As on_demand, but scaling from previous iteration is reused.",
      "high_delay", "Scaling is used after more than 0.05*n delays are present.",
      "high_delay_reuse", "Scaling is used only when previous iteration created "
      "more that 0.05*n additional delays; otherwise, reuse scaling from the previous iteration.",
      "od_hd", "Combination of on_demand and high_delay.",
      "od_hd_reuse", "Combination of on_demand_reuse and high_delay_reuse",
      "If spral_scaling = dynamic, spral_scaling_3 is enabled according to this condition.",
      true);

   roptions->AddLowerBoundedNumberOption(
      "spral_small",
      "Zero Pivot Threshold",
      0.0, true,
      1.0e-20,
      "Any pivot less than spral_small is treated as zero.");

   roptions->AddLowerBoundedNumberOption(
      "spral_small_subtree_threshold",
      "Small Subtree Threshold",
      0.0, true,
      4.0e6,
      "Maximum number of FLOPS in a subtree treated as a single task.");

   roptions->AddBoundedNumberOption(
      "spral_u",
      "Pivoting Threshold",
      0.0, true,
      0.5, false,
      1.0e-8,
      "Relative pivot threshold used in symmetric indefinite case.");

   roptions->AddBoundedNumberOption(
      "spral_umax",
      "Maximum Pivoting Threshold",
      0.0, true,
      0.5, false,
      1.0e-4,
      "See SPRAL documentation.");

   roptions->AddBoolOption(
      "spral_use_gpu",
      "Specifies whether or not graphics processing units (GPUs) are used by the SPRAL linear solver if present.",
      true);
}

int SpralSolverInterface::PivotMethodNameToNum(
   const std::string& name
)
{
   if( name == "aggressive" )
   {
      return 0;
   }
   else if( name == "block" )
   {
      return 1;
   }
   else if( name == "threshold" )
   {
      return 2;
   }
   else
   {
      assert(0);
      return -1;
   }
}

int SpralSolverInterface::ScaleNameToNum(
   const std::string& name
)
{
   if( name == "none" )
   {
      return 0;
   }
   else if( name == "mc64" )
   {
      return 1;
   }
   else if( name == "auction" )
   {
      return 2;
   }
   else if( name == "matching" )
   {
      return 3;
   }
   else if( name == "ruiz" )
   {
      return 4;
   }
   else
   {
      assert(0);
      return -1;
   }
}

bool SpralSolverInterface::InitializeImpl(
   const OptionsList& options,
   const std::string& prefix
)
{
   spral_ssids_default_options(&control_);
   control_.array_base = 1; // Use Fortran numbering
   control_.action = true; // Continue factorization on discovery of a zero pivot.
   /* Note: we can't set control_.action = false as we need to know the
    * inertia. (Otherwise we just enter the restoration phase and fail.) */

   options.GetBoolValue("spral_ignore_numa", control_.ignore_numa, prefix);
   options.GetBoolValue("spral_use_gpu", control_.use_gpu, prefix);
   options.GetIntegerValue("spral_cpu_block_size", control_.cpu_block_size, prefix);
   options.GetIntegerValue("spral_nemin", control_.nemin, prefix);
   options.GetIntegerValue("spral_print_level", control_.print_level, prefix);
   options.GetNumericValue("spral_small", control_.small, prefix);
   options.GetNumericValue("spral_u", control_.u, prefix);
   options.GetNumericValue("spral_umax", umax_, prefix);

   // Set gpu_perf_coeff.
   double gpu_perf_coeff_tmp = 1.0;
   options.GetNumericValue("spral_gpu_perf_coeff", gpu_perf_coeff_tmp, prefix);
   control_.gpu_perf_coeff = (float)gpu_perf_coeff_tmp;

   // Set max_load_inbalance.
   double max_load_inbalance_tmp = 1.2;
   options.GetNumericValue("spral_max_load_inbalance", max_load_inbalance_tmp, prefix);
   control_.max_load_inbalance = (float)max_load_inbalance_tmp;

   // Set min_gpu_work.
   double min_gpu_work_tmp = 5.0e9;
   options.GetNumericValue("spral_min_gpu_work", min_gpu_work_tmp, prefix);
   control_.min_gpu_work = (int64_t)min_gpu_work_tmp;

   // Set the pivot method.
   std::string pivot_method;
   options.GetStringValue("spral_pivot_method", pivot_method, prefix);  // TODO use GetEnumValue?
   control_.pivot_method = PivotMethodNameToNum(pivot_method);

   // Set small_subtree_threshold.
   double small_subtree_threshold_tmp = 4.0e6;
   options.GetNumericValue("spral_small_subtree_threshold", small_subtree_threshold_tmp, prefix);
   control_.small_subtree_threshold = (int64_t)small_subtree_threshold_tmp;

   // Reset all private data.
   pivtol_changed_ = false;

   std::string order_method;
   options.GetStringValue("spral_order", order_method, prefix);

   if( order_method == "metis" )
   {
      control_.ordering = 1;
   }
   else if( order_method == "matching" )
   {
      control_.ordering = 2;
   }

   std::string scaling_method;
   options.GetStringValue("spral_scaling", scaling_method, prefix);
   current_level_ = 0;

   if( scaling_method == "dynamic" )
   {
      scaling_type_ = 0;
      std::string switch_val[3], scaling_val[3];

      options.GetStringValue("spral_switch_1",   switch_val[0], prefix);
      options.GetStringValue("spral_scaling_1", scaling_val[0], prefix);
      options.GetStringValue("spral_switch_2",   switch_val[1], prefix);
      options.GetStringValue("spral_scaling_2", scaling_val[1], prefix);
      options.GetStringValue("spral_switch_3",   switch_val[2], prefix);
      options.GetStringValue("spral_scaling_3", scaling_val[2], prefix);

      for( int i = 0; i < 3; i++ )
      {
         scaling_val_[i] = ScaleNameToNum(scaling_val[i]);

         if( switch_val[i] == "never" )
         {
            switch_[i] = SWITCH_NEVER;
         }
         else if( switch_val[i] == "at_start" )
         {
            switch_[i] = SWITCH_AT_START;
            scaling_type_ = scaling_val_[i];
            current_level_ = i;
            Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA, "SPRAL: Enabled "
                           "scaling level %d on initialization\n", current_level_);
         }
         else if( switch_val[i] == "at_start_reuse" )
         {
            switch_[i] = SWITCH_AT_START_REUSE;
            scaling_type_ = scaling_val_[i];
            current_level_ = i;
            Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA, "SPRAL: Enabled "
                           "scaling level %d on initialization\n", current_level_);
         }
         else if( switch_val[i] == "on_demand" )
         {
            switch_[i] = SWITCH_ON_DEMAND;
         }
         else if( switch_val[i] == "on_demand_reuse" )
         {
            switch_[i] = SWITCH_ON_DEMAND_REUSE;
         }
         else if( switch_val[i] == "high_delay" )
         {
            switch_[i] = SWITCH_NDELAY;
         }
         else if( switch_val[i] == "high_delay_reuse" )
         {
            switch_[i] = SWITCH_NDELAY_REUSE;
         }
         else if( switch_val[i] == "od_hd" )
         {
            switch_[i] = SWITCH_OD_ND;
         }
         else if( switch_val[i] == "od_hd_reuse" )
         {
            switch_[i] = SWITCH_OD_ND_REUSE;
         }
      }
   }
   else
   {
      switch_[0] = SWITCH_AT_START;
      switch_[1] = SWITCH_NEVER;
      switch_[2] = SWITCH_NEVER;
      scaling_type_ = ScaleNameToNum(scaling_method);
   }

   // Set whether we scale on first iteration or not
   switch ( switch_[current_level_] )
   {
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

   // Set scaling and ordering.
   control_.scaling = scaling_type_;
   control_.ordering = scaling_type_ != 3 ? 1 : 2;

   return true; // All is well.
}

ESymSolverStatus SpralSolverInterface::InitializeStructure(
   Index        dim,
   Index        nonzeros,
   const Index* ia,
   const Index* ja
)
{
   struct spral_ssids_inform info;

   // Store size for later use
   ndim_ = dim;

   // Setup memory for values
   delete[] val_;
   val_ = new double[nonzeros];

   // Correct scaling and ordering if necessary.
   if( control_.ordering == 2 && control_.scaling != 3 )
   {
      Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA, "In SpralSolverInterface, "
                     "matching-based ordering was used, but matching-based scaling was "
                     "not. Setting scaling using the matching-based ordering.\n");
      control_.scaling = scaling_type_ = 3;
   }

   if( control_.ordering != 2 && control_.scaling == 3 )
   {
      Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA, "In SpralSolverInterface, "
                     "matching-based scaling was used, but matching-based ordering was "
                     "not. Setting ordering using the matching-based algorithm.\n");
      control_.ordering = 2;
   }

   // Perform analyse.
   if( !( control_.ordering == 2 && control_.scaling == 3 ) )
   {
      if( HaveIpData() )
      {
         IpData().TimingStats().LinearSystemSymbolicFactorization().Start();
      }

      spral_ssids_analyse_ptr32(false, ndim_, NULL, ia, ja, NULL, &akeep_, &control_, &info);

      Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA, "nfactor = %" PRId64 ", nflops = %" PRId64 ":\n",
                     info.num_factor, info.num_flops);

      if( HaveIpData() )
      {
         IpData().TimingStats().LinearSystemSymbolicFactorization().End();
      }

      if( info.flag < 0 )
      {
         Jnlst().Printf(J_STRONGWARNING, J_LINEAR_ALGEBRA, "In SpralSolverInterface::InitializeStructure: "
                        "Unhandled error. info.flag = %d.\n", info.flag);
         return SYMSOLVER_FATAL_ERROR;
      }
   }

   return SYMSOLVER_SUCCESS;
}

ESymSolverStatus SpralSolverInterface::MultiSolve(
   bool         new_matrix,
   const Index* ia,
   const Index* ja,
   Index        nrhs,
   double*      rhs_vals,
   bool         check_NegEVals,
   Index        numberOfNegEVals
)
{
   struct spral_ssids_inform info;

   if( new_matrix || pivtol_changed_ )
   {
      // Set scaling option
      if( rescale_ )
      {
         control_.scaling = scaling_type_;
         control_.ordering = scaling_type_ != 3 ? 1 : 2;

         if( scaling_type_ != 0 && scaling_ == NULL )
         {
            scaling_ = new double[ndim_];
         }
      }
      else
      {
         control_.scaling = 0; // None or user (depends if scaling_ is allocated).
      }

      if( control_.ordering == 2 && control_.scaling == 3 && rescale_ )
      {
         if( HaveIpData() )
         {
            IpData().TimingStats().LinearSystemSymbolicFactorization().Start();
         }

         spral_ssids_analyse_ptr32(false, ndim_, NULL, ia, ja, val_, &akeep_, &control_, &info);

         Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA, "nfactor = %" PRId64 ", nflops = %" PRId64 ":\n",
                        info.num_factor, info.num_flops);

         if( HaveIpData() )
         {
            IpData().TimingStats().LinearSystemSymbolicFactorization().End();
         }

         if ( info.flag == 6 || info.flag == -5 )
         {
            Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA, "In SpralSolverInterface::Factorization: "
                           "Singular system, estimated rank %d of %d\n", info.matrix_rank, ndim_);
            return SYMSOLVER_SINGULAR;
         }

         if ( info.flag < 0 )
         {
            Jnlst().Printf(J_STRONGWARNING, J_LINEAR_ALGEBRA, "In SpralSolverInterface::Factorization: "
                           "Unhandled error. info.flag = %d.\n", info.flag);
            return SYMSOLVER_FATAL_ERROR;
         }
      }

      Number t1 = 0;
      if( HaveIpData() )
      {
         t1 = IpData().TimingStats().LinearSystemFactorization().TotalWallclockTime();
         IpData().TimingStats().LinearSystemFactorization().Start();
      }

      spral_ssids_factor_ptr32(false, ia, ja, val_, scaling_, akeep_, &fkeep_, &control_, &info);

      Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA, "SPRAL: delays %d, nfactor %" PRId64 ", "
                     "nflops %" PRId64 ", maxfront %d\n", info.num_delay, info.num_factor, info.num_flops,
                     info.maxfront);

      if( HaveIpData() )
      {
         IpData().TimingStats().LinearSystemFactorization().End();
         Number t2 = IpData().TimingStats().LinearSystemFactorization().TotalWallclockTime();
         Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA, "SpralSolverInterface::Factorization: "
                        "spral_factor_solve took %10.3f\n", t2 - t1);
      }

      if( info.flag == 7 || info.flag == 6 || info.flag == -5 )
      {
         Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA, "In SpralSolverInterface::Factorization: "
                        "Singular system, estimated rank %d of %d\n", info.matrix_rank, ndim_);
         return SYMSOLVER_SINGULAR;
      }

      if( info.flag < 0 )
      {
         Jnlst().Printf(J_STRONGWARNING, J_LINEAR_ALGEBRA, "In SpralSolverInterface::Factorization: "
                        "Unhandled error. info.flag = %d.\n", info.flag);
         if( info.flag == -53 )
         {
            Jnlst().Printf(J_STRONGWARNING, J_LINEAR_ALGEBRA, "Maybe one forgot to set environment variable OMP_CANCELLATION to TRUE.\n");
         }
         if( control_.print_level < 0 )
         {
            Jnlst().Printf(J_STRONGWARNING, J_LINEAR_ALGEBRA, "Set spral_print_level=0 to see more details.\n");
         }
         return SYMSOLVER_FATAL_ERROR;
      }

      for( int i = current_level_; i < 3; i++ )
      {
         switch( switch_[i] )
         {
            case SWITCH_NEVER:
            case SWITCH_AT_START:
            case SWITCH_ON_DEMAND:
               // Nothing to do here.
               break;
            case SWITCH_AT_START_REUSE:
               // Scaled exactly once, never changed again.
               rescale_ = false;
               break;
            case SWITCH_ON_DEMAND_REUSE:
               if( i == current_level_ && rescale_ )
               {
                  rescale_ = false;
               }

               break;
            case SWITCH_NDELAY_REUSE:
            case SWITCH_OD_ND_REUSE:
               if( rescale_ )
               {
                  // Need to do this before we reset rescale.
                  numdelay_ = info.num_delay;
               }

               if( i == current_level_ && rescale_ )
               {
                  rescale_ = false;
               }
            // fall through
            case SWITCH_NDELAY:
            case SWITCH_OD_ND:
               if( rescale_ )
               {
                  numdelay_ = info.num_delay;
               }

               if( info.num_delay - numdelay_ > 0.05 * ndim_ )
               {
                  // Number of delays has signficantly increased, so trigger.
                  current_level_ = i;
                  scaling_type_ = scaling_val_[i];
                  Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA, "SPRAL: "
                                 "Enabling scaling %d due to excess delays\n", i);
                  rescale_ = true;
               }

               break;
         }
      }

      if( check_NegEVals && info.num_neg != numberOfNegEVals )
      {
         Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA, "In SpralSolverInterface::Factorization: "
                        "info.num_neg = %d, but numberOfNegEVals = %" IPOPT_INDEX_FORMAT "\n", info.num_neg, numberOfNegEVals);
         return SYMSOLVER_WRONG_INERTIA;
      }

      if( HaveIpData() )
      {
         IpData().TimingStats().LinearSystemBackSolve().Start();
      }

      spral_ssids_solve(0, nrhs, rhs_vals, ndim_, akeep_, fkeep_, &control_, &info);

      if( HaveIpData() )
      {
         IpData().TimingStats().LinearSystemBackSolve().End();
      }

      numneg_ = info.num_neg;

      pivtol_changed_ = false;
   }
   else
   {
      if( HaveIpData() )
      {
         IpData().TimingStats().LinearSystemBackSolve().Start();
      }

      spral_ssids_solve(0, nrhs, rhs_vals, ndim_, akeep_, fkeep_, &control_, &info);

      if( HaveIpData() )
      {
         IpData().TimingStats().LinearSystemBackSolve().End();
      }
   }

   return SYMSOLVER_SUCCESS;
}

bool SpralSolverInterface::IncreaseQuality()
{
   for( int i = current_level_; i < 3; i++ )
   {
      switch( switch_[i] )
      {
         case SWITCH_ON_DEMAND:
         case SWITCH_ON_DEMAND_REUSE:
         case SWITCH_OD_ND:
         case SWITCH_OD_ND_REUSE:
            rescale_ = true;
            current_level_ = i;
            scaling_type_ = scaling_val_[i];
            Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA, "SPRAL: Enabling scaling "
                           "%d due to failure of iterative refinement\n", current_level_);
            break;
         default:
            ;
      }
   }

   if( control_.u >= umax_ )
   {
      return false;
   }

   pivtol_changed_ = true;
   Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA, "Increasing pivot tolerance "
                  "for SPRAL from %7.2e ", control_.u);
   control_.u = Min(umax_, std::pow(control_.u, 0.75));
   Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA, "to %7.2e.\n", control_.u);
   return true;
}

} // namespace Ipopt

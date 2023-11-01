// Copyright (C) 2013, Science and Technology Facilities Council.
// Copyright (C) 2009, Jonathan Hogg <jdh41.at.cantab.net>
// Copyright (C) 2004, 2007 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Authors: Jonathan Hogg                    STFC   2012-02-14
//          Jonathan Hogg                           2009-07-29
//          Carl Laird, Andreas Waechter     IBM    2004-03-17

#include "IpoptConfig.h"
#include "IpMa77SolverInterface.hpp"

#include <iostream>
#include <cmath>

#ifdef IPOPT_HAS_HSL
#include "CoinHslConfig.h"
#endif

using namespace std;

namespace Ipopt
{
static IPOPT_DECL_MA77_DEFAULT_CONTROL(*user_ma77_default_control) = NULL;
static IPOPT_DECL_MA77_OPEN_NELT(*user_ma77_open_nelt) = NULL;
static IPOPT_DECL_MA77_OPEN(*user_ma77_open) = NULL;
static IPOPT_DECL_MA77_INPUT_VARS(*user_ma77_input_vars) = NULL;
static IPOPT_DECL_MA77_INPUT_REALS(*user_ma77_input_reals) = NULL;
static IPOPT_DECL_MA77_ANALYSE(*user_ma77_analyse) = NULL;
static IPOPT_DECL_MA77_FACTOR(*user_ma77_factor) = NULL;
static IPOPT_DECL_MA77_FACTOR_SOLVE(*user_ma77_factor_solve) = NULL;
static IPOPT_DECL_MA77_SOLVE(*user_ma77_solve) = NULL;
static IPOPT_DECL_MA77_RESID(*user_ma77_resid) = NULL;
static IPOPT_DECL_MA77_SCALE(*user_ma77_scale) = NULL;
static IPOPT_DECL_MA77_ENQUIRE_POSDEF(*user_ma77_enquire_posdef) = NULL;
static IPOPT_DECL_MA77_ENQUIRE_INDEF(*user_ma77_enquire_indef) = NULL;
static IPOPT_DECL_MA77_ALTER(*user_ma77_alter) = NULL;
static IPOPT_DECL_MA77_RESTART(*user_ma77_restart) = NULL;
static IPOPT_DECL_MA77_FINALISE(*user_ma77_finalise) = NULL;
static IPOPT_DECL_MC68_DEFAULT_CONTROL(*user_mc68_default_control) = NULL;
static IPOPT_DECL_MC68_ORDER(*user_mc68_order) = NULL;

Ma77SolverInterface::~Ma77SolverInterface()
{
   delete[] val_;

   struct ma77_info info;
   if( keep_ )
   {
      ma77_finalise(&keep_, &control_, &info);
   }
}

void Ma77SolverInterface::RegisterOptions(
   SmartPtr<RegisteredOptions> roptions
)
{
   roptions->AddIntegerOption(
      "ma77_print_level",
      "Debug printing level for the linear solver MA77",
      -1,
      "<0: no printing; 0: Error and warning messages only; 1: Limited diagnostic printing; >1 Additional diagnostic printing.");
   roptions->AddLowerBoundedIntegerOption(
      "ma77_buffer_lpage",
      "Number of scalars per MA77 in-core buffer page in the out-of-core solver MA77",
      1,
      4096,
      "Must be at most ma77_file_size.");
   roptions->AddLowerBoundedIntegerOption(
      "ma77_buffer_npage",
      "Number of pages that make up MA77 buffer",
      1,
      1600,
      "Number of pages of size buffer_lpage that exist in-core for the out-of-core solver MA77.");
   roptions->AddLowerBoundedIntegerOption(
      "ma77_file_size",
      "Target size of each temporary file for MA77, scalars per type",
      1,
      2097152,
      "MA77 uses many temporary files, this option controls the size of each one. "
      "It is measured in the number of entries (int or double), NOT bytes.");
   roptions->AddLowerBoundedIntegerOption(
      "ma77_maxstore",
      "Maximum storage size for MA77 in-core mode",
      0,
      0,
      "If greater than zero, the maximum size of factors stored in core before out-of-core mode is invoked.");
   roptions->AddLowerBoundedIntegerOption(
      "ma77_nemin",
      "Node Amalgamation parameter",
      1,
      8,
      "Two nodes in elimination tree are merged if result has fewer than ma77_nemin variables.");
   roptions->AddLowerBoundedNumberOption(
      "ma77_small",
      "Zero Pivot Threshold",
      0.0, false,
      1e-20,
      "Any pivot less than ma77_small is treated as zero.");
   roptions->AddLowerBoundedNumberOption(
      "ma77_static",
      "Static Pivoting Threshold", 0.0, false,
      0.0,
      "See MA77 documentation. "
      "Either ma77_static=0.0 or ma77_static>ma77_small. "
      "ma77_static=0.0 disables static pivoting.");
   roptions->AddBoundedNumberOption(
      "ma77_u",
      "Pivoting Threshold",
      0.0, false,
      0.5, false,
      1e-8,
      "See MA77 documentation.");
   roptions->AddBoundedNumberOption(
      "ma77_umax",
      "Maximum Pivoting Threshold",
      0.0, false,
      0.5, false,
      1e-4,
      "Maximum value to which u will be increased to improve quality.");
   roptions->AddStringOption2(
      "ma77_order",
      "Controls type of ordering used by MA77",
#ifdef COINHSL_HAS_METIS
      "metis",
#else
      "amd",
#endif
      "amd", "Use the HSL_MC68 approximate minimum degree algorithm",
      "metis", "Use the MeTiS nested dissection algorithm (if available)");
}

/// set MA77 functions to use for every instantiation of this class
void Ma77SolverInterface::SetFunctions(
   IPOPT_DECL_MA77_DEFAULT_CONTROL(*ma77_default_control),
   IPOPT_DECL_MA77_OPEN_NELT(*ma77_open_nelt),
   IPOPT_DECL_MA77_OPEN(*ma77_open),
   IPOPT_DECL_MA77_INPUT_VARS(*ma77_input_vars),
   IPOPT_DECL_MA77_INPUT_REALS(*ma77_input_reals),
   IPOPT_DECL_MA77_ANALYSE(*ma77_analyse),
   IPOPT_DECL_MA77_FACTOR(*ma77_factor),
   IPOPT_DECL_MA77_FACTOR_SOLVE(*ma77_factor_solve),
   IPOPT_DECL_MA77_SOLVE(*ma77_solve),
   IPOPT_DECL_MA77_RESID(*ma77_resid),
   IPOPT_DECL_MA77_SCALE(*ma77_scale),
   IPOPT_DECL_MA77_ENQUIRE_POSDEF(*ma77_enquire_posdef),
   IPOPT_DECL_MA77_ENQUIRE_INDEF(*ma77_enquire_indef),
   IPOPT_DECL_MA77_ALTER(*ma77_alter),
   IPOPT_DECL_MA77_RESTART(*ma77_restart),
   IPOPT_DECL_MA77_FINALISE(*ma77_finalise),
   IPOPT_DECL_MC68_DEFAULT_CONTROL(*mc68_default_control),
   IPOPT_DECL_MC68_ORDER(*mc68_order)
)
{
   DBG_ASSERT(ma77_default_control != NULL);
   DBG_ASSERT(ma77_open_nelt != NULL);
   DBG_ASSERT(ma77_open != NULL);
   DBG_ASSERT(ma77_input_vars != NULL);
   DBG_ASSERT(ma77_input_reals != NULL);
   DBG_ASSERT(ma77_analyse != NULL);
   DBG_ASSERT(ma77_factor != NULL);
   DBG_ASSERT(ma77_factor_solve != NULL);
   DBG_ASSERT(ma77_solve != NULL);
   DBG_ASSERT(ma77_resid != NULL);
   DBG_ASSERT(ma77_scale != NULL);
   DBG_ASSERT(ma77_enquire_posdef != NULL);
   DBG_ASSERT(ma77_enquire_indef != NULL);
   DBG_ASSERT(ma77_alter != NULL);
   DBG_ASSERT(ma77_restart != NULL);
   DBG_ASSERT(ma77_finalise != NULL);
   DBG_ASSERT(mc68_default_control != NULL);
   DBG_ASSERT(mc68_order != NULL);

   user_ma77_default_control = ma77_default_control;
   user_ma77_open_nelt = ma77_open_nelt;
   user_ma77_open = ma77_open;
   user_ma77_input_vars = ma77_input_vars;
   user_ma77_input_reals = ma77_input_reals;
   user_ma77_analyse = ma77_analyse;
   user_ma77_factor = ma77_factor;
   user_ma77_factor_solve = ma77_factor_solve;
   user_ma77_solve = ma77_solve;
   user_ma77_resid = ma77_resid;
   user_ma77_scale = ma77_scale;
   user_ma77_enquire_posdef = ma77_enquire_posdef;
   user_ma77_enquire_indef = ma77_enquire_indef;
   user_ma77_alter = ma77_alter;
   user_ma77_restart = ma77_restart;
   user_ma77_finalise = ma77_finalise;
   user_mc68_default_control = mc68_default_control;
   user_mc68_order = mc68_order;
}

bool Ma77SolverInterface::InitializeImpl(
   const OptionsList& options,
   const std::string& prefix
)
{
   if( user_ma77_default_control != NULL )
   {
      ma77_default_control = user_ma77_default_control;
      ma77_open_nelt = user_ma77_open_nelt;
      ma77_open = user_ma77_open;
      ma77_input_vars = user_ma77_input_vars;
      ma77_input_reals = user_ma77_input_reals;
      ma77_analyse = user_ma77_analyse;
      ma77_factor = user_ma77_factor;
      ma77_factor_solve = user_ma77_factor_solve;
      ma77_solve = user_ma77_solve;
      ma77_resid = user_ma77_resid;
      ma77_scale = user_ma77_scale;
      ma77_enquire_posdef = user_ma77_enquire_posdef;
      ma77_enquire_indef = user_ma77_enquire_indef;
      ma77_alter = user_ma77_alter;
      ma77_restart = user_ma77_restart;
      ma77_finalise = user_ma77_finalise;
      mc68_default_control = user_mc68_default_control;
      mc68_order = user_mc68_order;
   }
   else
   {
#if (defined(COINHSL_HAS_MA77) && !defined(IPOPT_SINGLE)) || (defined(COINHSL_HAS_MA77S) && defined(IPOPT_SINGLE))
      // use HSL functions that should be available in linked HSL library
      ma77_default_control = &::ma77_default_control;
      ma77_open_nelt = &::ma77_open_nelt;
      ma77_open = &::ma77_open;
      ma77_input_vars = &::ma77_input_vars;
      ma77_input_reals = &::ma77_input_reals;
      ma77_analyse = &::ma77_analyse;
      ma77_factor = &::ma77_factor;
      ma77_factor_solve = &::ma77_factor_solve;
      ma77_solve = &::ma77_solve;
      ma77_resid = &::ma77_resid;
      ma77_scale = &::ma77_scale;
      ma77_enquire_posdef = &::ma77_enquire_posdef;
      ma77_enquire_indef = &::ma77_enquire_indef;
      ma77_alter = &::ma77_alter;
      ma77_restart = &::ma77_restart;
      ma77_finalise = &::ma77_finalise;
      mc68_default_control = &::mc68_default_control;
      mc68_order = &::mc68_order;
#else
      // try to load HSL functions from a shared library at runtime
      DBG_ASSERT(IsValid(hslloader));

#define STR2(x) #x
#define STR(x) STR2(x)
      ma77_default_control = (IPOPT_DECL_MA77_DEFAULT_CONTROL(*))hslloader->loadSymbol(STR(ma77_default_control));
      ma77_open_nelt = (IPOPT_DECL_MA77_OPEN_NELT(*))hslloader->loadSymbol(STR(ma77_open_nelt));
      ma77_open = (IPOPT_DECL_MA77_OPEN(*))hslloader->loadSymbol(STR(ma77_open));
      ma77_input_vars = (IPOPT_DECL_MA77_INPUT_VARS(*))hslloader->loadSymbol(STR(ma77_input_vars));
      ma77_input_reals = (IPOPT_DECL_MA77_INPUT_REALS(*))hslloader->loadSymbol(STR(ma77_input_reals));
      ma77_analyse = (IPOPT_DECL_MA77_ANALYSE(*))hslloader->loadSymbol(STR(ma77_analyse));
      ma77_factor = (IPOPT_DECL_MA77_FACTOR(*))hslloader->loadSymbol(STR(ma77_factor));
      ma77_factor_solve = (IPOPT_DECL_MA77_FACTOR_SOLVE(*))hslloader->loadSymbol(STR(ma77_factor_solve));
      ma77_solve = (IPOPT_DECL_MA77_SOLVE(*))hslloader->loadSymbol(STR(ma77_solve));
      ma77_resid = (IPOPT_DECL_MA77_RESID(*))hslloader->loadSymbol(STR(ma77_resid));
      ma77_scale = (IPOPT_DECL_MA77_SCALE(*))hslloader->loadSymbol(STR(ma77_scale));
      ma77_enquire_posdef = (IPOPT_DECL_MA77_ENQUIRE_POSDEF(*))hslloader->loadSymbol(STR(ma77_enquire_posdef));
      ma77_enquire_indef = (IPOPT_DECL_MA77_ENQUIRE_INDEF(*))hslloader->loadSymbol(STR(ma77_enquire_indef));
      ma77_alter = (IPOPT_DECL_MA77_ALTER(*))hslloader->loadSymbol(STR(ma77_alter));
      ma77_restart = (IPOPT_DECL_MA77_RESTART(*))hslloader->loadSymbol(STR(ma77_restart));
      ma77_finalise = (IPOPT_DECL_MA77_FINALISE(*))hslloader->loadSymbol(STR(ma77_finalise));
      mc68_default_control = (IPOPT_DECL_MC68_DEFAULT_CONTROL(*))hslloader->loadSymbol(STR(mc68_default_control));
      mc68_order = (IPOPT_DECL_MC68_ORDER(*))hslloader->loadSymbol(STR(mc68_order));
#endif
   }

   DBG_ASSERT(ma77_default_control != NULL);
   DBG_ASSERT(ma77_open_nelt != NULL);
   DBG_ASSERT(ma77_open != NULL);
   DBG_ASSERT(ma77_input_vars != NULL);
   DBG_ASSERT(ma77_input_reals != NULL);
   DBG_ASSERT(ma77_analyse != NULL);
   DBG_ASSERT(ma77_factor != NULL);
   DBG_ASSERT(ma77_factor_solve != NULL);
   DBG_ASSERT(ma77_solve != NULL);
   DBG_ASSERT(ma77_resid != NULL);
   DBG_ASSERT(ma77_scale != NULL);
   DBG_ASSERT(ma77_enquire_posdef != NULL);
   DBG_ASSERT(ma77_enquire_indef != NULL);
   DBG_ASSERT(ma77_alter != NULL);
   DBG_ASSERT(ma77_restart != NULL);
   DBG_ASSERT(ma77_finalise != NULL);
   DBG_ASSERT(mc68_default_control != NULL);
   DBG_ASSERT(mc68_order != NULL);

   ma77_default_control(&control_);
   control_.f_arrays = 1; // Use Fortran numbering (faster)
   control_.bits = 32;
   // FIXME: HSL_MA77 should be updated to allow a matrix with new
   // values to be refactorized after a -11 (singular) error.
   //control_.action = 0; // false, should exit with error on singularity

   Index temp;
   options.GetIntegerValue("ma77_print_level", temp, prefix);
   control_.print_level = temp;
   options.GetIntegerValue("ma77_buffer_lpage", temp, prefix);
   control_.buffer_lpage[0] = temp;
   options.GetIntegerValue("ma77_buffer_lpage", temp, prefix);
   control_.buffer_lpage[1] = temp;
   options.GetIntegerValue("ma77_buffer_npage", temp, prefix);
   control_.buffer_npage[0] = temp;
   options.GetIntegerValue("ma77_buffer_npage", temp, prefix);
   control_.buffer_npage[1] = temp;
   options.GetIntegerValue("ma77_file_size", temp, prefix);
   control_.file_size = temp;
   options.GetIntegerValue("ma77_maxstore", temp, prefix);
   control_.maxstore = temp;
   options.GetIntegerValue("ma77_nemin", temp, prefix);
   control_.nemin = temp;
   options.GetNumericValue("ma77_small", control_.small, prefix);
   options.GetNumericValue("ma77_static", control_.static_, prefix);
   options.GetNumericValue("ma77_u", control_.u, prefix);
   options.GetNumericValue("ma77_umax", umax_, prefix);

   std::string order_method;
   options.GetStringValue("ma77_order", order_method, prefix);
   if( order_method == "metis" )
   {
      ordering_ = ORDER_METIS;
   }
   else
   {
      ordering_ = ORDER_AMD;
   }

   return true; // All is well
}

ESymSolverStatus Ma77SolverInterface::InitializeStructure(
   Index        dim,
   Index        nonzeros,
   const Index* ia,
   const Index* ja
)
{
   struct ma77_info info;
   struct mc68_control control68;
   struct mc68_info info68;

   // Store size for later use
   ndim_ = dim;

   if( HaveIpData() )
   {
      IpData().TimingStats().LinearSystemSymbolicFactorization().Start();
   }

   // mc68 requires a half matrix. A future version will support full
   // matrix entry, and this code should be removed when it is available.
   int* ia_half = new int[dim + 1];
   int* ja_half = new int[ia[dim] - 1];
   {
      int k = 0;
      for( int i = 0; i < dim; i++ )
      {
         ia_half[i] = k + 1;
         for( int j = ia[i] - 1; j < ia[i + 1] - 1; j++ )
            if( ja[j] - 1 >= i )
            {
               ja_half[k++] = ja[j];
            }
      }
      ia_half[dim] = k + 1;
   }

   // Determine an ordering
   mc68_default_control(&control68);
   control68.f_array_in = 1; // Use Fortran numbering (faster)
   control68.f_array_out = 1; // Use Fortran numbering (faster)
   int* perm = new int[dim];
   if( ordering_ == ORDER_METIS )
   {
      mc68_order(3, dim, ia_half, ja_half, perm, &control68, &info68); /* MeTiS */
      if( info68.flag == -5 )
      {
         // MeTiS not available
         ordering_ = ORDER_AMD;
      }
      else if( info68.flag < 0 )
      {
         delete[] ia_half;
         delete[] ja_half;
         delete[] perm;
         return SYMSOLVER_FATAL_ERROR;
      }
   }
   if( ordering_ == ORDER_AMD )
   {
      mc68_order(1, dim, ia_half, ja_half, perm, &control68, &info68); /* AMD */
      if( info68.flag < 0 )
      {
         delete[] ia_half;
         delete[] ja_half;
         delete[] perm;
         return SYMSOLVER_FATAL_ERROR;
      }
   }
   delete[] ia_half;
   delete[] ja_half;

   // Open files
   ma77_open(ndim_, "ma77_int", "ma77_real", "ma77_work", "ma77_delay", &keep_, &control_, &info);
   if( info.flag < 0 )
   {
      delete[] perm;
      return SYMSOLVER_FATAL_ERROR;
   }

   // Store data into files
   for( int i = 0; i < dim; i++ )
   {
      ma77_input_vars(i + 1, ia[i + 1] - ia[i], &(ja[ia[i] - 1]), &keep_, &control_, &info);
      if( info.flag < 0 )
      {
         delete[] perm;
         return SYMSOLVER_FATAL_ERROR;
      }
   }

   // Perform analyse
   ma77_analyse(perm, &keep_, &control_, &info);
   delete[] perm; // Done with order

   if( HaveIpData() )
   {
      IpData().TimingStats().LinearSystemSymbolicFactorization().End();
   }

   // Setup memory for values
   if( val_ != NULL )
   {
      delete[] val_;
   }
   val_ = new Number[nonzeros];

   if( info.flag >= 0 )
   {
      return SYMSOLVER_SUCCESS;
   }
   else
   {
      return SYMSOLVER_FATAL_ERROR;
   }
}

ESymSolverStatus Ma77SolverInterface::MultiSolve(
   bool         new_matrix,
   const Index* ia,
   const Index* /*ja*/,
   Index        nrhs,
   Number*      rhs_vals,
   bool         check_NegEVals,
   Index        numberOfNegEVals
)
{
   struct ma77_info info;

   if( new_matrix || pivtol_changed_ )
   {
      for( int i = 0; i < ndim_; i++ )
      {
         ma77_input_reals(i + 1, ia[i + 1] - ia[i], &(val_[ia[i] - 1]), &keep_, &control_, &info);
         if( info.flag < 0 )
         {
            return SYMSOLVER_FATAL_ERROR;
         }
      }

      if( HaveIpData() )
      {
         IpData().TimingStats().LinearSystemFactorization().Start();
      }
      //ma77_factor(0, &keep_, &control_, &info, NULL);
      ma77_factor_solve(0, &keep_, &control_, &info, NULL, nrhs, ndim_, rhs_vals);
      if( HaveIpData() )
      {
         IpData().TimingStats().LinearSystemFactorization().End();
      }
      if( info.flag == 4 || info.flag == -11 )
      {
         return SYMSOLVER_SINGULAR;
      }
      if( info.flag < 0 )
      {
         return SYMSOLVER_FATAL_ERROR;
      }
      if( check_NegEVals && info.num_neg != numberOfNegEVals )
      {
         return SYMSOLVER_WRONG_INERTIA;
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
      ma77_solve(0, nrhs, ndim_, rhs_vals, &keep_, &control_, &info, NULL);
      if( HaveIpData() )
      {
         IpData().TimingStats().LinearSystemBackSolve().End();
      }
   }

   return SYMSOLVER_SUCCESS;
}

bool Ma77SolverInterface::IncreaseQuality()
{
   if( control_.u >= umax_ )
   {
      return false;
   }
   pivtol_changed_ = true;

   Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                  "Increasing pivot tolerance for HSL_MA77 from %7.2e ", control_.u);
   control_.u = Min(umax_, std::pow(control_.u, Number(0.75)));
   Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                  "to %7.2e.\n", control_.u);
   return true;
}

} // namespace Ipopt

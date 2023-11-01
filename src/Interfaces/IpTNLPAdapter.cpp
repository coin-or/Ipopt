// Copyright (C) 2004, 2010 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#include "IpoptConfig.h"
#include "IpTNLPAdapter.hpp"
#include "IpBlas.hpp"
#include "IpIpoptData.hpp"
#include "IpIpoptCalculatedQuantities.hpp"
#include "IpDenseVector.hpp"
#include "IpExpansionMatrix.hpp"
#include "IpGenTMatrix.hpp"
#include "IpSymTMatrix.hpp"
#include "IpTDependencyDetector.hpp"
#include "IpTSymDependencyDetector.hpp"
#include "IpTripletToCSRConverter.hpp"

#ifdef IPOPT_HAS_HSL
#include "CoinHslConfig.h"
#include "IpMa28TDependencyDetector.hpp"
#endif

#ifdef IPOPT_HAS_MUMPS
# include "IpMumpsSolverInterface.hpp"
#endif
#ifdef IPOPT_HAS_WSMP
# include "IpWsmpSolverInterface.hpp"
#endif

#include <cmath>
#include <cstdio>
#include <cstring>

namespace Ipopt
{
#if IPOPT_VERBOSITY > 0
static const Index dbg_verbosity = 0;
#endif

TNLPAdapter::TNLPAdapter(
   const SmartPtr<TNLP>             tnlp,
   const SmartPtr<const Journalist> jnlst /* = NULL */
)
   : tnlp_(tnlp),
     jnlst_(jnlst),
     full_x_(NULL),
     full_lambda_(NULL),
     full_g_(NULL),
     jac_g_(NULL),
     c_rhs_(NULL),
     x_tag_for_iterates_(0),
     y_c_tag_for_iterates_(0),
     y_d_tag_for_iterates_(0),
     x_tag_for_g_(0),
     x_tag_for_jac_g_(0),
     jac_idx_map_(NULL),
     h_idx_map_(NULL),
     x_fixed_map_(NULL),
     findiff_jac_ia_(NULL),
     findiff_jac_ja_(NULL),
     findiff_jac_postriplet_(NULL),
     findiff_x_l_(NULL),
     findiff_x_u_(NULL)
{
   ASSERT_EXCEPTION(IsValid(tnlp_), INVALID_TNLP, "The TNLP passed to TNLPAdapter is NULL. This MUST be a valid TNLP!");
}

TNLPAdapter::~TNLPAdapter()
{
   delete[] full_x_;
   delete[] full_lambda_;
   delete[] full_g_;
   delete[] jac_g_;
   delete[] c_rhs_;
   delete[] jac_idx_map_;
   delete[] h_idx_map_;
   delete[] x_fixed_map_;
   delete[] findiff_jac_ia_;
   delete[] findiff_jac_ja_;
   delete[] findiff_jac_postriplet_;
   delete[] findiff_x_l_;
   delete[] findiff_x_u_;
}

void TNLPAdapter::RegisterOptions(
   SmartPtr<RegisteredOptions> roptions
)
{
   roptions->SetRegisteringCategory("NLP");
   roptions->AddNumberOption(
      "nlp_lower_bound_inf",
      "any bound less or equal this value will be considered -inf (i.e. not lower bounded).",
      -1e19);
   roptions->AddNumberOption(
      "nlp_upper_bound_inf",
      "any bound greater or this value will be considered +inf (i.e. not upper bounded).",
      1e19);
   roptions->AddStringOption4(
      "fixed_variable_treatment",
      "Determines how fixed variables should be handled.",
      "make_parameter",
      "make_parameter", "Remove fixed variable from optimization variables",
      "make_parameter_nodual", "Remove fixed variable from optimization variables and do not compute bound multipliers for fixed variables",
      "make_constraint", "Add equality constraints fixing variables",
      "relax_bounds", "Relax fixing bound constraints",
      "The main difference between those options is that the starting point in the \"make_constraint\" case still "
      "has the fixed variables at their given values, whereas in the case \"make_parameter(_nodual)\" the functions are always "
      "evaluated with the fixed values for those variables.  "
      "Also, for \"relax_bounds\", the fixing bound constraints are relaxed (according to\" bound_relax_factor\"). "
      "For all but \"make_parameter_nodual\", bound multipliers are computed for the fixed variables.");

   std::vector<std::string> options;
   std::vector<std::string> descrs;
   options.push_back("none");
   descrs.push_back("don't check; no extra work at beginning");
#ifdef IPOPT_HAS_MUMPS
   options.push_back("mumps");
   descrs.push_back("use MUMPS");
#endif
#ifdef IPOPT_HAS_WSMP
   options.push_back("wsmp");
   descrs.push_back("use WSMP");
#endif
#if ((defined(COINHSL_HAS_MA28) && !defined(IPOPT_SINGLE)) || (defined(COINHSL_HAS_MA28S) && defined(IPOPT_SINGLE))) && defined(F77_FUNC) && !defined(IPOPT_INT64)
   options.push_back("ma28");
   descrs.push_back("use MA28");
#endif
   roptions->AddStringOption(
      "dependency_detector",
      "Indicates which linear solver should be used to detect linearly dependent equality constraints.",
      "none",
      options,
      descrs,
      "This is experimental and does not work well.",
      true);
   roptions->AddBoolOption(
      "dependency_detection_with_rhs",
      "Indicates if the right hand sides of the constraints should be considered in addition to gradients during dependency detection",
      false,
      "",
      true);
   roptions->AddLowerBoundedIntegerOption(
      "num_linear_variables",
      "Number of linear variables",
      0,
      0,
      "When the Hessian is approximated, it is assumed that the first num_linear_variables variables are linear. "
      "The Hessian is then not approximated in this space. "
      "If the get_number_of_nonlinear_variables method in the TNLP is implemented, this option is ignored.",
      true);
   roptions->AddStringOption2(
      "jacobian_approximation",
      "Specifies technique to compute constraint Jacobian",
      "exact",
      "exact", "user-provided derivatives",
      "finite-difference-values", "user-provided structure, values by finite differences",
      "",
      true);
   roptions->AddStringOption2(
      "gradient_approximation",
      "Specifies technique to compute objective Gradient",
      "exact",
      "exact", "user-provided gradient",
      "finite-difference-values", "values by finite differences",
      "",
      true);
   roptions->AddLowerBoundedNumberOption(
      "findiff_perturbation",
      "Size of the finite difference perturbation for derivative approximation.",
      0., true,
      1e-7,
      "This determines the relative perturbation of the variable entries.",
      true);

   roptions->SetRegisteringCategory("Derivative Checker");
   roptions->AddStringOption4(
      "derivative_test",
      "Enable derivative checker",
      "none",
      "none", "do not perform derivative test",
      "first-order", "perform test of first derivatives at starting point",
      "second-order", "perform test of first and second derivatives at starting point",
      "only-second-order", "perform test of second derivatives at starting point",
      "If this option is enabled, a (slow!) derivative test will be performed before the optimization. "
      "The test is performed at the user provided starting point and marks derivative values that seem suspicious");
   roptions->AddLowerBoundedIntegerOption(
      "derivative_test_first_index",
      "Index of first quantity to be checked by derivative checker",
      -2,
      -2,
      "If this is set to -2, then all derivatives are checked. "
      "Otherwise, for the first derivative test it specifies the first variable for which the test is done (counting starts at 0). "
      "For second derivatives, it specifies the first constraint for which the test is done; "
      "counting of constraint indices starts at 0, and -1 refers to the objective function Hessian.");
   roptions->AddLowerBoundedNumberOption(
      "derivative_test_perturbation",
      "Size of the finite difference perturbation in derivative test.",
      0., true,
      1e-8,
      "This determines the relative perturbation of the variable entries.");
   roptions->AddLowerBoundedNumberOption(
      "derivative_test_tol",
      "Threshold for indicating wrong derivative.",
      0., true,
      1e-4,
      "If the relative deviation of the estimated derivative from the given one is larger than this value, "
      "the corresponding derivative is marked as wrong.");
   roptions->AddBoolOption(
      "derivative_test_print_all",
      "Indicates whether information for all estimated derivatives should be printed.",
      false,
      "Determines verbosity of derivative checker.");
   roptions->AddLowerBoundedNumberOption(
      "point_perturbation_radius",
      "Maximal perturbation of an evaluation point.",
      0., false,
      10.,
      "If a random perturbation of a points is required, this number indicates the maximal perturbation. "
      "This is for example used when determining the center point at which the finite difference derivative test is executed.");
}

bool TNLPAdapter::ProcessOptions(
   const OptionsList& options,
   const std::string& prefix
)
{
   DBG_START_METH("TNLPAdapter::ProcessOptions", dbg_verbosity);
   options.GetNumericValue("nlp_lower_bound_inf", nlp_lower_bound_inf_, prefix);
   options.GetNumericValue("nlp_upper_bound_inf", nlp_upper_bound_inf_, prefix);

   ASSERT_EXCEPTION(nlp_lower_bound_inf_ < nlp_upper_bound_inf_, OPTION_INVALID,
                    "Option \"nlp_lower_bound_inf\" must be smaller than \"nlp_upper_bound_inf\".");

   // Registered in IpOrigIpoptNLP
   options.GetNumericValue("bound_relax_factor", bound_relax_factor_, prefix);

   Index enum_int;
   options.GetEnumValue("fixed_variable_treatment", enum_int, prefix);
   fixed_variable_treatment_ = FixedVariableTreatmentEnum(enum_int);
   options.GetEnumValue("derivative_test", enum_int, prefix);
   derivative_test_ = DerivativeTestEnum(enum_int);
   options.GetNumericValue("derivative_test_perturbation", derivative_test_perturbation_, prefix);
   options.GetNumericValue("derivative_test_tol", derivative_test_tol_, prefix);
   options.GetBoolValue("derivative_test_print_all", derivative_test_print_all_, prefix);
   options.GetIntegerValue("derivative_test_first_index", derivative_test_first_index_, prefix);

   // The option warm_start_same_structure is registered by OrigIpoptNLP
   options.GetBoolValue("warm_start_same_structure", warm_start_same_structure_, prefix);
   // The following is registered in OrigIpoptNLP
   options.GetEnumValue("hessian_approximation", enum_int, prefix);
   hessian_approximation_ = HessianApproximationType(enum_int);
   options.GetIntegerValue("num_linear_variables", num_linear_variables_, prefix);

   options.GetEnumValue("jacobian_approximation", enum_int, prefix);
   jacobian_approximation_ = JacobianApproxEnum(enum_int);
   options.GetEnumValue("gradient_approximation", enum_int, prefix);
   gradient_approximation_ = GradientApproxEnum(enum_int);
   options.GetNumericValue("findiff_perturbation", findiff_perturbation_, prefix);

   options.GetNumericValue("point_perturbation_radius", point_perturbation_radius_, prefix);

   options.GetNumericValue("tol", tol_, prefix);

   options.GetBoolValue("dependency_detection_with_rhs", dependency_detection_with_rhs_, prefix);
   std::string dependency_detector;
   options.GetStringValue("dependency_detector", dependency_detector, prefix);
#ifdef IPOPT_HAS_MUMPS
   if( dependency_detector == "mumps" )
   {
      SmartPtr<SparseSymLinearSolverInterface> SolverInterface;
      SolverInterface = new MumpsSolverInterface();
      SmartPtr<TSymLinearSolver> ScaledSolver =
         new TSymLinearSolver(SolverInterface, NULL);
      dependency_detector_ = new TSymDependencyDetector(*ScaledSolver);
   }
#endif

#ifdef IPOPT_HAS_WSMP
   if( dependency_detector == "wsmp" )
   {
      SmartPtr<SparseSymLinearSolverInterface> SolverInterface;
#ifdef PARDISO_MATCHING_PREPROCESS
      std::string libname;
      options.GetStringValue("pardisolib", libname, prefix);
      SolverInterface = new WsmpSolverInterface(new LibraryLoader(libname));
#else
      SolverInterface = new WsmpSolverInterface();
#endif
      SmartPtr<TSymLinearSolver> ScaledSolver =
         new TSymLinearSolver(SolverInterface, NULL);
      dependency_detector_ = new TSymDependencyDetector(*ScaledSolver);
   }
#endif

#if ((defined(COINHSL_HAS_MA28) && !defined(IPOPT_SINGLE)) || (defined(COINHSL_HAS_MA28S) && defined(IPOPT_SINGLE))) && defined(F77_FUNC) && !defined(IPOPT_INT64)
   if( dependency_detector == "ma28" )
   {
      dependency_detector_ = new Ma28TDependencyDetector();
   }
#endif

   if( IsValid(dependency_detector_) )
      if( !dependency_detector_->ReducedInitialize(*jnlst_, options, prefix) )
      {
         return false;
      }

   return true;
}

bool TNLPAdapter::GetSpaces(
   SmartPtr<const VectorSpace>&    x_space,
   SmartPtr<const VectorSpace>&    c_space,
   SmartPtr<const VectorSpace>&    d_space,
   SmartPtr<const VectorSpace>&    x_l_space,
   SmartPtr<const MatrixSpace>&    px_l_space,
   SmartPtr<const VectorSpace>&    x_u_space,
   SmartPtr<const MatrixSpace>&    px_u_space,
   SmartPtr<const VectorSpace>&    d_l_space,
   SmartPtr<const MatrixSpace>&    pd_l_space,
   SmartPtr<const VectorSpace>&    d_u_space,
   SmartPtr<const MatrixSpace>&    pd_u_space,
   SmartPtr<const MatrixSpace>&    Jac_c_space,
   SmartPtr<const MatrixSpace>&    Jac_d_space,
   SmartPtr<const SymMatrixSpace>& Hess_lagrangian_space
)
{
   DBG_START_METH("TNLPAdapter::GetSpaces", dbg_verbosity);

   // First, if required, perform derivative test
   if( derivative_test_ != NO_TEST )
   {
      bool retval = CheckDerivatives(derivative_test_, derivative_test_first_index_);
      if( !retval )
      {
         return retval;
      }
   }

   if( warm_start_same_structure_ )
   {
      ASSERT_EXCEPTION(full_x_, INVALID_WARMSTART,
                       "warm_start_same_structure chosen, but TNLPAdapter is called for the first time.");
      if( IsValid(jnlst_) )
      {
         jnlst_->Printf(J_DETAILED, J_INITIALIZATION, "Reusing previous information for warm start in TNLPAdapter.\n");
      }
   }
   else
   {
      // In case the Adapter has been used before, but this is not a
      // warm start, make sure we delete all previously allocated
      // memory
      delete[] full_x_;
      full_x_ = NULL;
      delete[] full_lambda_;
      full_lambda_ = NULL;
      delete[] full_g_;
      full_g_ = NULL;
      delete[] jac_g_;
      jac_g_ = NULL;
      delete[] c_rhs_;
      c_rhs_ = NULL;
      delete[] jac_idx_map_;
      jac_idx_map_ = NULL;
      delete[] h_idx_map_;
      h_idx_map_ = NULL;
      delete[] x_fixed_map_;
      x_fixed_map_ = NULL;
      jac_fixed_idx_map_.clear();
      jac_fixed_iRow_.clear();
      jac_fixed_jCol_.clear();
#if __cplusplus >= 201103L
      jac_fixed_idx_map_.shrink_to_fit();
      jac_fixed_iRow_.shrink_to_fit();
      jac_fixed_jCol_.shrink_to_fit();
#endif
   }

   // Get the full dimensions of the problem
   Index n_full_x, n_full_g, nz_full_jac_g, nz_full_h;
   bool retval = tnlp_->get_nlp_info(n_full_x, n_full_g, nz_full_jac_g, nz_full_h, index_style_);
   ASSERT_EXCEPTION(retval, INVALID_TNLP, "get_nlp_info returned false");
   ASSERT_EXCEPTION(n_full_x >= 0, INVALID_TNLP, "number of variables negative");
   ASSERT_EXCEPTION(n_full_g >= 0, INVALID_TNLP, "number of constraints negative");
   ASSERT_EXCEPTION(nz_full_jac_g >= 0, INVALID_TNLP, "number of nonzeros in Jacobian negative");
   ASSERT_EXCEPTION(nz_full_h >= 0, INVALID_TNLP, "number of nonzeros in Hessian negative");
   ASSERT_EXCEPTION(
      !warm_start_same_structure_
      || (n_full_x == n_full_x_ && n_full_g == n_full_g_ && nz_full_jac_g == nz_full_jac_g_
          && nz_full_h == nz_full_h_), INVALID_WARMSTART,
      "warm_start_same_structure chosen, but problem dimensions are different.");
   n_full_x_ = n_full_x;
   n_full_g_ = n_full_g;
   nz_full_jac_g_ = nz_full_jac_g;
   nz_full_h_ = nz_full_h;

   if( !warm_start_same_structure_ )
   {
      // create space to store vectors that are the full length of x
      full_x_ = new Number[n_full_x_];

      // create space to store vectors that area the full length of lambda
      full_lambda_ = new Number[n_full_g_];

      // create space to store vectors that are the full length of g
      full_g_ = new Number[n_full_g_];
      // check if there is any meta data for the variables and constraints
      StringMetaDataMapType var_string_md;
      IntegerMetaDataMapType var_integer_md;
      NumericMetaDataMapType var_numeric_md;
      StringMetaDataMapType con_string_md;
      IntegerMetaDataMapType con_integer_md;
      NumericMetaDataMapType con_numeric_md;
      if( !tnlp_->get_var_con_metadata(n_full_x_, var_string_md, var_integer_md, var_numeric_md, n_full_g_,
                                       con_string_md, con_integer_md, con_numeric_md) )
      {
         var_string_md.clear();
         var_integer_md.clear();
         var_numeric_md.clear();
         con_string_md.clear();
         con_integer_md.clear();
         con_numeric_md.clear();
      }

      // allocate internal space to store the full jacobian
      jac_g_ = new Number[nz_full_jac_g_];

      /* Spaces for bounds. We need to remove the fixed variables
       * and find out which bounds do not exist. */
      Number* x_l = new Number[n_full_x_];
      Number* x_u = new Number[n_full_x_];
      Number* g_l = new Number[n_full_g_];
      Number* g_u = new Number[n_full_g_];
      retval = tnlp_->get_bounds_info(n_full_x_, x_l, x_u, n_full_g_, g_l, g_u);
      ASSERT_EXCEPTION(retval, INVALID_TNLP, "get_bounds_info returned false in GetSpaces");

      //*********************************************************
      // Create the spaces and permutation spaces
      //*********************************************************

      Index n_x_var;
      Index n_x_l;
      Index n_x_u;
      Index* x_not_fixed_map = new Index[n_full_x_];
      Index* x_l_map = new Index[n_full_x_];
      Index* x_u_map = new Index[n_full_x_];

      Index n_c;
      Index n_d;
      Index n_d_l;
      Index n_d_u;
      Index* c_map = new Index[n_full_g_]; // we do not know n_c yet!
      Index* d_map = new Index[n_full_g_]; // we do not know n_d yet!
      Index* d_l_map = new Index[n_full_g_]; // "
      Index* d_u_map = new Index[n_full_g_]; // "

      bool done = false;
      // We might have to do the following twice: If we detect that we
      // don't have enough degrees of freedom, we simply redo
      // everything with fixed_variable_treatment to set RELAX_BOUNDS
      while( !done )
      {
         n_x_var = 0;
         n_x_l = 0;
         n_x_u = 0;
         n_x_fixed_ = 0;
         Index* x_fixed_map_tmp = new Index[n_full_x_];

         for( Index i = 0; i < n_full_x_; i++ )
         {
            Number lower_bound = x_l[i];
            Number upper_bound = x_u[i];
            if( lower_bound == upper_bound )
            {
               switch( fixed_variable_treatment_ )
               {
                  case MAKE_PARAMETER:
                  case MAKE_PARAMETER_NODUAL:
                     // Variable is fixed, remove it from the problem
                     full_x_[i] = lower_bound;
                     x_fixed_map_tmp[n_x_fixed_] = i;
                     n_x_fixed_++;
                     break;
                  case MAKE_CONSTRAINT:
                     x_fixed_map_tmp[n_x_fixed_] = i; // don't really need this
                     // array then
                     n_x_fixed_++;
                     x_not_fixed_map[n_x_var] = i;
                     n_x_var++;
                     break;
                  case RELAX_BOUNDS:
                     x_l_map[n_x_l] = n_x_var;
                     n_x_l++;
                     x_u_map[n_x_u] = n_x_var;
                     n_x_u++;
                     n_x_var++;
                     break;
                  default:
                     DBG_ASSERT(false && "invalid fixed_variable_treatment_");
               }
            }
            else if( lower_bound > upper_bound )
            {
               char string[128];
               Snprintf(string, 127,
                        "There are inconsistent bounds on variable %" IPOPT_INDEX_FORMAT ": lower = %25.16e and upper = %25.16e.", i, lower_bound,
                        upper_bound);
               delete[] x_l;
               delete[] x_u;
               delete[] g_l;
               delete[] g_u;
               delete[] x_not_fixed_map;
               delete[] x_fixed_map_tmp;
               delete[] x_l_map;
               delete[] x_u_map;
               delete[] c_map;
               delete[] d_map;
               delete[] d_l_map;
               delete[] d_u_map;
               THROW_EXCEPTION(INCONSISTENT_BOUNDS, string);
            }
            else
            {
               x_not_fixed_map[n_x_var] = i;
               if( lower_bound > nlp_lower_bound_inf_ )
               {
                  x_l_map[n_x_l] = n_x_var;
                  n_x_l++;
               }

               if( upper_bound < nlp_upper_bound_inf_ )
               {
                  x_u_map[n_x_u] = n_x_var;
                  n_x_u++;
               }
               n_x_var++;
            }
         }

         // If there are fixed variables, we keep their position around
         // for a possible warm start later or if fixed variables are
         // treated by added equality constraints
         if( n_x_fixed_ > 0 )
         {
            delete[] x_fixed_map_;
            x_fixed_map_ = NULL;
            x_fixed_map_ = new Index[n_x_fixed_];
            for( Index i = 0; i < n_x_fixed_; i++ )
            {
               x_fixed_map_[i] = x_fixed_map_tmp[i];
            }
         }
         else
         {
            delete[] x_fixed_map_;
            x_fixed_map_ = NULL;
         }
         delete[] x_fixed_map_tmp;

         // Create the spaces for c and d
         // - includes the internal permutation matrices for
         //  full_g to c and d
         // - includes the permutation matrices for d_l and d_u
         // c(x) = (P_c)T * g(x)
         // d(x) = (P_d)T * g(x)
         // d_L = (P_d_L)T * (P_d)T * g_l
         // d_U = (P_d_U)T * (P_d)T * g_u
         n_c = 0;
         n_d = 0;
         n_d_l = 0;
         n_d_u = 0;

         for( Index i = 0; i < n_full_g_; i++ )
         {
            Number lower_bound = g_l[i];
            Number upper_bound = g_u[i];
            if( lower_bound == upper_bound )
            {
               // equality constraint
               c_map[n_c] = i;
               n_c++;
            }
            else if( lower_bound > upper_bound )
            {
               delete[] x_l;
               delete[] x_u;
               delete[] g_l;
               delete[] g_u;
               delete[] x_not_fixed_map;
               delete[] x_l_map;
               delete[] x_u_map;
               delete[] c_map;
               delete[] d_map;
               delete[] d_l_map;
               delete[] d_u_map;
               char string[128];
               Snprintf(string, 127,
                        "There are inconsistent bounds on constraint function %" IPOPT_INDEX_FORMAT ": lower = %25.16e and upper = %25.16e.", i,
                        lower_bound, upper_bound);
               THROW_EXCEPTION(INCONSISTENT_BOUNDS, string);
            }
            else
            {
               // inequality constraint
               d_map[n_d] = i;
               if( lower_bound > nlp_lower_bound_inf_ )
               {
                  d_l_map[n_d_l] = n_d;
                  n_d_l++;
               }
               if( upper_bound < nlp_upper_bound_inf_ )
               {
                  d_u_map[n_d_u] = n_d;
                  n_d_u++;
               }
               n_d++;
            }
         }

         if( fixed_variable_treatment_ == RELAX_BOUNDS || n_x_fixed_ == 0 || n_x_var >= n_c || n_x_var == 0 )
         {
            done = true;
         }
         else
         {
            fixed_variable_treatment_ = RELAX_BOUNDS;
            jnlst_->Printf(J_WARNING, J_INITIALIZATION,
                           "Too few degrees of freedom (n_x = %" IPOPT_INDEX_FORMAT ", n_c = %" IPOPT_INDEX_FORMAT ").\n  Trying fixed_variable_treatment = RELAX_BOUNDS\n\n",
                           n_x_var, n_c);
         }
      } // while (!done)

      // If requested, check if there are linearly dependent equality
      // constraints
      if( n_c > 0 && IsValid(dependency_detector_) )
      {
         std::list<Index> c_deps;
         if( !DetermineDependentConstraints(n_x_var, x_not_fixed_map, x_l, x_u, g_l, g_u, n_c, c_map, c_deps) )
         {
            jnlst_->Printf(J_WARNING, J_INITIALIZATION,
                           "Dependent constraint detector had a problem, assume full rank.\n");
         }
         c_deps.sort();
         if( !c_deps.empty() )
         {
            jnlst_->Printf(J_WARNING, J_INITIALIZATION,
                           "\nDetected %zd linearly dependent equality constraints; taking those out.\n\n", c_deps.size());
         }
         else
         {
            jnlst_->Printf(J_DETAILED, J_INITIALIZATION, "\nNo dependent constraints detected.\n\n");
         }
         if( jnlst_->ProduceOutput(J_DETAILED, J_INITIALIZATION) )
         {
            jnlst_->Printf(J_DETAILED, J_INITIALIZATION, "\nList of indices of dependent constraints:\n");
            Index count = 0;
            for( std::list<Index>::iterator i = c_deps.begin(); i != c_deps.end(); ++i )
            {
               jnlst_->Printf(J_DETAILED, J_INITIALIZATION, "c_dep[%" IPOPT_INDEX_FORMAT "] = %" IPOPT_INDEX_FORMAT "\n", count++, *i);
            }
            jnlst_->Printf(J_DETAILED, J_INITIALIZATION, "\n");
         }
         if( !c_deps.empty() )
         {
            // Take the dependent constraints out.
            // We assume that the list in c_dep is sorted
            std::list<Index>::iterator idep = c_deps.begin();
            Index new_n_c = *idep;
            for( Index i = *idep; i < n_c; i++ )
            {
               if( i == *idep )
               {
                  ++idep;
               }
               else
               {
                  c_map[new_n_c] = c_map[i];
                  new_n_c++;
               }
               if( idep == c_deps.end() )
               {
                  // just copy the rest and done
                  for( Index j = i + 1; j < n_c; j++ )
                  {
                     c_map[new_n_c++] = c_map[j];
                  }
                  break;
               }
            }
            n_c = new_n_c;
            // We also need to set the multipliers for the constraints
            // to zero (could do only for dependent ones... was too lazy
            // right now)
            const Number zero = 0.;
            IpBlasCopy(n_full_g_, &zero, 0, full_lambda_, 1);
         }
      }
      delete[] x_l;
      x_l = NULL;
      delete[] x_u;
      x_u = NULL;

      // create x spaces
      SmartPtr<DenseVectorSpace> dv_x_space = new DenseVectorSpace(n_x_var);
      x_space_ = GetRawPtr(dv_x_space);
      SmartPtr<DenseVectorSpace> dv_x_l_space = new DenseVectorSpace(n_x_l);
      x_l_space_ = GetRawPtr(dv_x_l_space);
      SmartPtr<DenseVectorSpace> dv_x_u_space = new DenseVectorSpace(n_x_u);
      x_u_space_ = GetRawPtr(dv_x_u_space);

      if( n_x_fixed_ > 0 && (fixed_variable_treatment_ == MAKE_PARAMETER || fixed_variable_treatment_ == MAKE_PARAMETER_NODUAL) )
      {
         P_x_full_x_space_ = new ExpansionMatrixSpace(n_full_x_, n_x_var, x_not_fixed_map);
         P_x_full_x_ = P_x_full_x_space_->MakeNewExpansionMatrix();
      }
      else
      {
         P_x_full_x_space_ = NULL;
         P_x_full_x_ = NULL;
      }

      P_x_x_L_space_ = new ExpansionMatrixSpace(n_x_var, n_x_l, x_l_map);
      px_l_space_ = GetRawPtr(P_x_x_L_space_);
      P_x_x_L_ = P_x_x_L_space_->MakeNewExpansionMatrix();
      P_x_x_U_space_ = new ExpansionMatrixSpace(n_x_var, n_x_u, x_u_map);
      px_u_space_ = GetRawPtr(P_x_x_U_space_);
      P_x_x_U_ = P_x_x_U_space_->MakeNewExpansionMatrix();

      // setup the variable meta data if present
      if( var_string_md.size() > 0 )
      {
         for( StringMetaDataMapType::iterator iter = var_string_md.begin(); iter != var_string_md.end(); ++iter )
         {
            std::vector<std::string> string_md(n_x_var);
            const Index* pos_idx = NULL;
            if( IsValid(P_x_full_x_space_) )
            {
               pos_idx = P_x_full_x_->ExpandedPosIndices();
               for( Index i = 0; i < n_x_var; i++ )
               {
                  string_md[i] = iter->second[pos_idx[i]];
               }
            }
            else
            {
               for( Index i = 0; i < n_x_var; i++ )
               {
                  string_md[i] = iter->second[i];
               }
            }
            dv_x_space->SetStringMetaData(iter->first, string_md);

            string_md.clear();
            string_md.resize(n_x_l);
            const Index* pos_idxL = P_x_x_L_space_->ExpandedPosIndices();
            for( Index i = 0; i < n_x_l; i++ )
            {
               if( pos_idx != NULL )
               {
                  string_md[i] = iter->second[pos_idx[pos_idxL[i]]];
               }
               else
               {
                  string_md[i] = iter->second[pos_idxL[i]];
               }
            }
            dv_x_l_space->SetStringMetaData(iter->first, string_md);

            string_md.clear();
            string_md.resize(n_x_u);
            const Index* pos_idxU = P_x_x_U_space_->ExpandedPosIndices();
            for( Index i = 0; i < n_x_u; i++ )
            {
               if( pos_idx != NULL )
               {
                  string_md[i] = iter->second[pos_idx[pos_idxU[i]]];
               }
               else
               {
                  string_md[i] = iter->second[pos_idxU[i]];
               }
            }
            dv_x_u_space->SetStringMetaData(iter->first, string_md);
         }
      }

      if( var_integer_md.size() > 0 )
      {
         for( IntegerMetaDataMapType::iterator iter = var_integer_md.begin(); iter != var_integer_md.end(); ++iter )
         {
            std::vector<Index> integer_md(n_x_var);
            const Index* pos_idx = NULL;
            if( IsValid(P_x_full_x_space_) )
            {
               pos_idx = P_x_full_x_->ExpandedPosIndices();
               for( Index i = 0; i < n_x_var; i++ )
               {
                  integer_md[i] = iter->second[pos_idx[i]];
               }
            }
            else
            {
               for( Index i = 0; i < n_x_var; i++ )
               {
                  integer_md[i] = iter->second[i];
               }
            }
            dv_x_space->SetIntegerMetaData(iter->first, integer_md);

            integer_md.clear();
            integer_md.resize(n_x_l);
            const Index* pos_idxL = P_x_x_L_space_->ExpandedPosIndices();
            for( Index i = 0; i < n_x_l; i++ )
            {
               if( pos_idx != NULL )
               {
                  integer_md[i] = iter->second[pos_idx[pos_idxL[i]]];
               }
               else
               {
                  integer_md[i] = iter->second[pos_idxL[i]];
               }
            }
            dv_x_l_space->SetIntegerMetaData(iter->first, integer_md);

            integer_md.clear();
            integer_md.resize(n_x_u);
            const Index* pos_idxU = P_x_x_U_space_->ExpandedPosIndices();
            for( Index i = 0; i < n_x_u; i++ )
            {
               if( pos_idx != NULL )
               {
                  integer_md[i] = iter->second[pos_idx[pos_idxU[i]]];
               }
               else
               {
                  integer_md[i] = iter->second[pos_idxU[i]];
               }
            }
            dv_x_u_space->SetIntegerMetaData(iter->first, integer_md);
         }
      }

      if( var_numeric_md.size() > 0 )
      {
         for( NumericMetaDataMapType::iterator iter = var_numeric_md.begin(); iter != var_numeric_md.end(); ++iter )
         {
            std::vector<Number> numeric_md(n_x_var);
            const Index* pos_idx = NULL;
            if( IsValid(P_x_full_x_space_) )
            {
               pos_idx = P_x_full_x_->ExpandedPosIndices();
               for( Index i = 0; i < n_x_var; i++ )
               {
                  numeric_md[i] = iter->second[pos_idx[i]];
               }
            }
            else
            {
               for( Index i = 0; i < n_x_var; i++ )
               {
                  numeric_md[i] = iter->second[i];
               }
            }
            dv_x_space->SetNumericMetaData(iter->first, numeric_md);

            numeric_md.clear();
            numeric_md.resize(n_x_l);
            const Index* pos_idxL = P_x_x_L_space_->ExpandedPosIndices();
            for( Index i = 0; i < n_x_l; i++ )
            {
               if( pos_idx != NULL )
               {
                  numeric_md[i] = iter->second[pos_idx[pos_idxL[i]]];
               }
               else
               {
                  numeric_md[i] = iter->second[pos_idxL[i]];
               }
            }
            dv_x_l_space->SetNumericMetaData(iter->first, numeric_md);

            numeric_md.clear();
            numeric_md.resize(n_x_u);
            const Index* pos_idxU = P_x_x_U_space_->ExpandedPosIndices();
            for( Index i = 0; i < n_x_u; i++ )
            {
               if( pos_idx != NULL )
               {
                  numeric_md[i] = iter->second[pos_idx[pos_idxU[i]]];
               }
               else
               {
                  numeric_md[i] = iter->second[pos_idxU[i]];
               }
            }
            dv_x_u_space->SetNumericMetaData(iter->first, numeric_md);
         }
      }

      delete[] x_not_fixed_map;
      x_not_fixed_map = NULL;
      delete[] x_l_map;
      x_l_map = NULL;
      delete[] x_u_map;
      x_u_map = NULL;

      // create the required c_space

      SmartPtr<DenseVectorSpace> dc_space;
      if( n_x_fixed_ == 0 || fixed_variable_treatment_ == MAKE_PARAMETER || fixed_variable_treatment_ == MAKE_PARAMETER_NODUAL )
      {
         dc_space = new DenseVectorSpace(n_c);
      }
      else
      {
         dc_space = new DenseVectorSpace(n_c + n_x_fixed_);
      }
      c_space_ = GetRawPtr(dc_space);
      c_rhs_ = new Number[dc_space->Dim()];

      // create the internal expansion matrix for c to g
      P_c_g_space_ = new ExpansionMatrixSpace(n_full_g_, n_c, c_map);
      P_c_g_ = P_c_g_space_->MakeNewExpansionMatrix();
      delete[] c_map;
      c_map = NULL;

      // create the required d_space
      SmartPtr<DenseVectorSpace> dv_d_space = new DenseVectorSpace(n_d);
      d_space_ = GetRawPtr(dv_d_space);
      // create the internal expansion matrix for d to g
      P_d_g_space_ = new ExpansionMatrixSpace(n_full_g_, n_d, d_map);
      P_d_g_ = P_d_g_space_->MakeNewExpansionMatrix();
      delete[] d_map;
      d_map = NULL;

      // create the required d_l space
      SmartPtr<DenseVectorSpace> dv_d_l_space = new DenseVectorSpace(n_d_l);
      d_l_space_ = GetRawPtr(dv_d_l_space);
      // create the required expansion matrix for d_L to d_L_exp
      SmartPtr<ExpansionMatrixSpace> P_d_l_space = new ExpansionMatrixSpace(n_d, n_d_l, d_l_map);
      pd_l_space_ = GetRawPtr(P_d_l_space);
      delete[] d_l_map;
      d_l_map = NULL;

      // create the required d_u space
      SmartPtr<DenseVectorSpace> dv_d_u_space = new DenseVectorSpace(n_d_u);
      d_u_space_ = GetRawPtr(dv_d_u_space);
      // create the required expansion matrix for d_U to d_U_exp
      SmartPtr<ExpansionMatrixSpace> P_d_u_space = new ExpansionMatrixSpace(n_d, n_d_u, d_u_map);
      pd_u_space_ = GetRawPtr(P_d_u_space);
      delete[] d_u_map;
      d_u_map = NULL;

      delete[] g_l;
      g_l = NULL;
      delete[] g_u;
      g_u = NULL;

      // set the constraint meta data if present
      if( con_string_md.size() > 0 )
      {
         for( StringMetaDataMapType::iterator iter = con_string_md.begin(); iter != con_string_md.end(); ++iter )
         {
            std::vector<std::string> string_md(n_c);
            const Index* pos_idx = P_c_g_space_->ExpandedPosIndices();
            for( Index i = 0; i < n_c; i++ )
            {
               string_md[i] = iter->second[pos_idx[i]];
            }
            dc_space->SetStringMetaData(iter->first, string_md);

            string_md.clear();
            string_md.resize(n_d);
            pos_idx = P_d_g_space_->ExpandedPosIndices();
            for( Index i = 0; i < n_d; i++ )
            {
               string_md[i] = iter->second[pos_idx[i]];
            }
            dv_d_space->SetStringMetaData(iter->first, string_md);

            string_md.clear();
            string_md.resize(n_d_l);
            const Index* d_pos_idx = P_d_l_space->ExpandedPosIndices();
            for( Index i = 0; i < n_d_l; i++ )
            {
               string_md[i] = iter->second[pos_idx[d_pos_idx[i]]];
            }
            dv_d_l_space->SetStringMetaData(iter->first, string_md);

            string_md.clear();
            string_md.resize(n_d_u);
            d_pos_idx = P_d_u_space->ExpandedPosIndices();
            for( Index i = 0; i < n_d_u; i++ )
            {
               string_md[i] = iter->second[pos_idx[d_pos_idx[i]]];
            }
            dv_d_u_space->SetStringMetaData(iter->first, string_md);
         }
      }

      if( con_integer_md.size() > 0 )
      {
         for( IntegerMetaDataMapType::iterator iter = con_integer_md.begin(); iter != con_integer_md.end(); ++iter )
         {
            std::vector<Index> integer_md(n_c);
            const Index* pos_idx = P_c_g_space_->ExpandedPosIndices();
            for( Index i = 0; i < n_c; i++ )
            {
               integer_md[i] = iter->second[pos_idx[i]];
            }
            dc_space->SetIntegerMetaData(iter->first, integer_md);

            integer_md.clear();
            integer_md.resize(n_d);
            pos_idx = P_d_g_space_->ExpandedPosIndices();
            for( Index i = 0; i < n_d; i++ )
            {
               integer_md[i] = iter->second[pos_idx[i]];
            }
            dv_d_space->SetIntegerMetaData(iter->first, integer_md);

            integer_md.clear();
            integer_md.resize(n_d_l);
            const Index* d_pos_idx = P_d_l_space->ExpandedPosIndices();
            for( Index i = 0; i < n_d_l; i++ )
            {
               integer_md[i] = iter->second[pos_idx[d_pos_idx[i]]];
            }
            dv_d_l_space->SetIntegerMetaData(iter->first, integer_md);

            integer_md.clear();
            integer_md.resize(n_d_u);
            d_pos_idx = P_d_u_space->ExpandedPosIndices();
            for( Index i = 0; i < n_d_u; i++ )
            {
               integer_md[i] = iter->second[pos_idx[d_pos_idx[i]]];
            }
            dv_d_u_space->SetIntegerMetaData(iter->first, integer_md);
         }
      }

      if( con_numeric_md.size() > 0 )
      {
         for( NumericMetaDataMapType::iterator iter = con_numeric_md.begin(); iter != con_numeric_md.end(); ++iter )
         {
            std::vector<Number> numeric_md(n_c);
            const Index* pos_idx = P_c_g_space_->ExpandedPosIndices();
            for( Index i = 0; i < n_c; i++ )
            {
               numeric_md[i] = iter->second[pos_idx[i]];
            }
            dc_space->SetNumericMetaData(iter->first, numeric_md);

            numeric_md.clear();
            numeric_md.resize(n_d);
            pos_idx = P_d_g_space_->ExpandedPosIndices();
            for( Index i = 0; i < n_d; i++ )
            {
               numeric_md[i] = iter->second[pos_idx[i]];
            }
            dv_d_space->SetNumericMetaData(iter->first, numeric_md);

            numeric_md.clear();
            numeric_md.resize(n_d_l);
            const Index* d_pos_idx = P_d_l_space->ExpandedPosIndices();
            for( Index i = 0; i < n_d_l; i++ )
            {
               numeric_md[i] = iter->second[pos_idx[d_pos_idx[i]]];
            }
            dv_d_l_space->SetNumericMetaData(iter->first, numeric_md);

            numeric_md.clear();
            numeric_md.resize(n_d_u);
            d_pos_idx = P_d_u_space->ExpandedPosIndices();
            for( Index i = 0; i < n_d_u; i++ )
            {
               numeric_md[i] = iter->second[pos_idx[d_pos_idx[i]]];
            }
            dv_d_u_space->SetNumericMetaData(iter->first, numeric_md);
         }
      }

      /** Create the matrix space for the jacobians
       */
      // Get the non zero structure
      Index* g_iRow = new Index[nz_full_jac_g_];
      Index* g_jCol = new Index[nz_full_jac_g_];

      bool retval = tnlp_->eval_jac_g(n_full_x_, NULL, false, n_full_g_, nz_full_jac_g_, g_iRow, g_jCol, NULL);

      if( !retval )
      {
         delete[] g_iRow;
         delete[] g_jCol;

         THROW_EXCEPTION(INVALID_TNLP, "eval_jac_g returned false in GetSpaces");
      }

      if( index_style_ != TNLP::FORTRAN_STYLE )
      {
         for( Index i = 0; i < nz_full_jac_g_; i++ )
         {
            DBG_ASSERT(g_iRow[i] >= 0);
            DBG_ASSERT(g_iRow[i] < n_full_g_);
            DBG_ASSERT(g_jCol[i] >= 0);
            DBG_ASSERT(g_jCol[i] < n_full_x_);
            g_iRow[i] += 1;
            g_jCol[i] += 1;
         }
      }
#if IPOPT_CHECKLEVEL > 0
      else
      {
         for( Index i = 0; i < nz_full_jac_g_; i++ )
         {
            DBG_ASSERT(g_iRow[i] > 0);
            DBG_ASSERT(g_iRow[i] <= n_full_g_);
            DBG_ASSERT(g_jCol[i] > 0);
            DBG_ASSERT(g_jCol[i] <= n_full_x_);
         }
      }
#endif

      if( jacobian_approximation_ == JAC_FINDIFF_VALUES )
      {
         initialize_findiff_jac(g_iRow, g_jCol);
      }

      // ... build the non-zero structure for jac_c
      // ... (the permutation from rows in jac_g to jac_c is
      // ...  the same as P_c_g_)
      Index nz_jac_all;
      if( fixed_variable_treatment_ == MAKE_PARAMETER || fixed_variable_treatment_ == MAKE_PARAMETER_NODUAL )
      {
         nz_jac_all = nz_full_jac_g_;
      }
      else
      {
         nz_jac_all = nz_full_jac_g_ + n_x_fixed_;
      }
      jac_idx_map_ = new Index[nz_jac_all];
      Index* jac_c_iRow = new Index[nz_jac_all];
      Index* jac_c_jCol = new Index[nz_jac_all];
      Index current_nz = 0;

      const Index* c_row_pos = P_c_g_->CompressedPosIndices();
      if( IsValid(P_x_full_x_) )
      {
         // there are missing variables x
         const Index* c_col_pos = P_x_full_x_->CompressedPosIndices();
         for( Index i = 0; i < nz_full_jac_g_; i++ )
         {
            const Index& c_row = c_row_pos[g_iRow[i] - 1];
            const Index& c_col = c_col_pos[g_jCol[i] - 1];
            if( c_col != -1 && c_row != -1 )
            {
               jac_idx_map_[current_nz] = i;
               jac_c_iRow[current_nz] = c_row + 1;
               jac_c_jCol[current_nz] = c_col + 1;
               current_nz++;
            }
            else if( c_col == -1 && fixed_variable_treatment_ == MAKE_PARAMETER )
            {
               // c_col == -1 should mean a fixed variables
               // c_row == -1 should mean a row in d(x)  (but the distinction into c(x) and d(x) isn't relevant for us here)
               jac_fixed_idx_map_.push_back(i);
               jac_fixed_iRow_.push_back(g_iRow[i]);
               jac_fixed_jCol_.push_back(g_jCol[i]);
            }
         }
      }
      else
      {
         for( Index i = 0; i < nz_full_jac_g_; i++ )
         {
            const Index& c_row = c_row_pos[g_iRow[i] - 1];
            const Index& c_col = g_jCol[i] - 1;
            if( c_row != -1 )
            {
               jac_idx_map_[current_nz] = i;
               jac_c_iRow[current_nz] = c_row + 1;
               jac_c_jCol[current_nz] = c_col + 1;
               current_nz++;
            }
         }
      }
      nz_jac_c_no_extra_ = current_nz;
      Index n_added_constr;
      if( fixed_variable_treatment_ == MAKE_PARAMETER || fixed_variable_treatment_ == MAKE_PARAMETER_NODUAL )
      {
         nz_jac_c_ = nz_jac_c_no_extra_;
         n_added_constr = 0;
      }
      else
      {
         nz_jac_c_ = nz_jac_c_no_extra_ + n_x_fixed_;
         for( Index i = 0; i < n_x_fixed_; i++ )
         {
            jac_c_iRow[current_nz] = n_c + i + 1;
            jac_c_jCol[current_nz] = x_fixed_map_[i] + 1;
            current_nz++;
         }
         n_added_constr = n_x_fixed_;
      }

      Jac_c_space_ = new GenTMatrixSpace(n_c + n_added_constr, n_x_var, nz_jac_c_, jac_c_iRow, jac_c_jCol);
      delete[] jac_c_iRow;
      jac_c_iRow = NULL;
      delete[] jac_c_jCol;
      jac_c_jCol = NULL;

      // ... build the nonzero structure for jac_d
      // ... (the permutation from rows in jac_g to jac_c is the
      // ...  the same as P_d_g_)
      Index* jac_d_iRow = new Index[nz_full_jac_g_];
      Index* jac_d_jCol = new Index[nz_full_jac_g_];
      current_nz = 0;
      const Index* d_row_pos = P_d_g_->CompressedPosIndices();
      if( IsValid(P_x_full_x_) )
      {
         const Index* d_col_pos = P_x_full_x_->CompressedPosIndices();
         for( Index i = 0; i < nz_full_jac_g_; i++ )
         {
            const Index& d_row = d_row_pos[g_iRow[i] - 1];
            const Index& d_col = d_col_pos[g_jCol[i] - 1];
            if( d_col != -1 && d_row != -1 )
            {
               jac_idx_map_[current_nz + nz_jac_c_no_extra_] = i;
               jac_d_iRow[current_nz] = d_row + 1;
               jac_d_jCol[current_nz] = d_col + 1;
               current_nz++;
            }
         }
      }
      else
      {
         for( Index i = 0; i < nz_full_jac_g_; i++ )
         {
            const Index& d_row = d_row_pos[g_iRow[i] - 1];
            const Index& d_col = g_jCol[i] - 1;
            if( d_row != -1 )
            {
               jac_idx_map_[current_nz + nz_jac_c_no_extra_] = i;
               jac_d_iRow[current_nz] = d_row + 1;
               jac_d_jCol[current_nz] = d_col + 1;
               current_nz++;
            }
         }
      }
      nz_jac_d_ = current_nz;
      Jac_d_space_ = new GenTMatrixSpace(n_d, n_x_var, nz_jac_d_, jac_d_iRow, jac_d_jCol);
      delete[] jac_d_iRow;
      jac_d_iRow = NULL;
      delete[] jac_d_jCol;
      jac_d_jCol = NULL;

      delete[] g_iRow;
      g_iRow = NULL;
      delete[] g_jCol;
      g_jCol = NULL;

      if( hessian_approximation_ == EXACT )
      {
         /** Create the matrix space for the hessian of the lagrangian */
         Index* full_h_iRow = new Index[nz_full_h_];
         Index* full_h_jCol = new Index[nz_full_h_];
         Index* h_iRow = new Index[nz_full_h_];
         Index* h_jCol = new Index[nz_full_h_];
         retval = tnlp_->eval_h(n_full_x_, NULL, false, 0, n_full_g_,
                                NULL, false, nz_full_h_, full_h_iRow, full_h_jCol, NULL);
         if( !retval )
         {
            delete[] full_h_iRow;
            delete[] full_h_jCol;
            delete[] h_iRow;
            delete[] h_jCol;
            jnlst_->Printf(J_ERROR, J_INITIALIZATION,
                           "Option \"hessian_approximation\" is not chosen as \"limited-memory\", but eval_h returns false.\n");
            THROW_EXCEPTION(OPTION_INVALID, "eval_h is called but has not been implemented");
         }

         if( index_style_ != TNLP::FORTRAN_STYLE )
         {
            for( Index i = 0; i < nz_full_h_; i++ )
            {
               DBG_ASSERT(full_h_iRow[i] >= 0);
               DBG_ASSERT(full_h_iRow[i] < n_full_x);
               DBG_ASSERT(full_h_jCol[i] >= 0);
               DBG_ASSERT(full_h_jCol[i] < n_full_x);
               full_h_iRow[i] += 1;
               full_h_jCol[i] += 1;
            }
         }
#if IPOPT_CHECKLEVEL > 0
         else
         {
            for( Index i = 0; i < nz_full_h_; i++ )
            {
               DBG_ASSERT(full_h_iRow[i] > 0);
               DBG_ASSERT(full_h_iRow[i] <= n_full_x);
               DBG_ASSERT(full_h_jCol[i] > 0);
               DBG_ASSERT(full_h_jCol[i] <= n_full_x);
            }
         }
#endif

         current_nz = 0;
         if( IsValid(P_x_full_x_) )
         {
            h_idx_map_ = new Index[nz_full_h_];
            const Index* h_pos = P_x_full_x_->CompressedPosIndices();
            for( Index i = 0; i < nz_full_h_; i++ )
            {
               const Index& h_row = h_pos[full_h_iRow[i] - 1];
               const Index& h_col = h_pos[full_h_jCol[i] - 1];
               if( h_row != -1 && h_col != -1 )
               {
                  h_idx_map_[current_nz] = i;
                  h_iRow[current_nz] = h_row + 1;
                  h_jCol[current_nz] = h_col + 1;
                  current_nz++;
               }
            }
         }
         else
         {
            h_idx_map_ = NULL;
            for( Index i = 0; i < nz_full_h_; i++ )
            {
               const Index& h_row = full_h_iRow[i] - 1;
               const Index& h_col = full_h_jCol[i] - 1;
               h_iRow[i] = h_row + 1;
               h_jCol[i] = h_col + 1;
               current_nz++;
            }
            current_nz = nz_full_h_;
         }
         nz_h_ = current_nz;
         Hess_lagrangian_space_ = new SymTMatrixSpace(n_x_var, nz_h_, h_iRow, h_jCol);
         delete[] full_h_iRow;
         full_h_iRow = NULL;
         delete[] full_h_jCol;
         full_h_jCol = NULL;
         delete[] h_iRow;
         h_iRow = NULL;
         delete[] h_jCol;
         h_jCol = NULL;
      }
      else
      {
         nz_h_ = 0;
         Hess_lagrangian_space_ = NULL;
      }
   } /* if (warm_start_same_structure_) { */

   // Assign the spaces to the returned pointers
   x_space = x_space_;
   c_space = c_space_;
   d_space = d_space_;
   x_l_space = x_l_space_;
   px_l_space = px_l_space_;
   x_u_space = x_u_space_;
   px_u_space = px_u_space_;
   d_l_space = d_l_space_;
   pd_l_space = pd_l_space_;
   d_u_space = d_u_space_;
   pd_u_space = pd_u_space_;
   Jac_c_space = Jac_c_space_;
   Jac_d_space = Jac_d_space_;
   Hess_lagrangian_space = Hess_lagrangian_space_;

   if( IsValid(jnlst_) )
   {
      jnlst_->Printf(J_ITERSUMMARY, J_STATISTICS, "Number of nonzeros in equality constraint Jacobian...:%9" IPOPT_INDEX_FORMAT "\n",
                     nz_jac_c_);
      jnlst_->Printf(J_ITERSUMMARY, J_STATISTICS, "Number of nonzeros in inequality constraint Jacobian.:%9" IPOPT_INDEX_FORMAT "\n",
                     nz_jac_d_);
      jnlst_->Printf(J_ITERSUMMARY, J_STATISTICS, "Number of nonzeros in Lagrangian Hessian.............:%9" IPOPT_INDEX_FORMAT "\n\n",
                     nz_h_);
   }

   return true;
}

bool TNLPAdapter::GetBoundsInformation(
   const Matrix& Px_L,
   Vector&       x_L,
   const Matrix& Px_U,
   Vector&       x_U,
   const Matrix& Pd_L,
   Vector&       d_L,
   const Matrix& Pd_U,
   Vector&       d_U
)
{
   // This could be done more efficiently, I have already called this method
   // once to setup the structure for the problem, I could store the values
   // and use them here ?
   // Actually, this is better for a warm start
   Number* x_l = new Number[n_full_x_];
   Number* x_u = new Number[n_full_x_];
   Number* g_l = new Number[n_full_g_];
   Number* g_u = new Number[n_full_g_];
   bool retval = tnlp_->get_bounds_info(n_full_x_, x_l, x_u, n_full_g_, g_l, g_u);
   ASSERT_EXCEPTION(retval, INVALID_TNLP, "get_bounds_info returned false in GetBoundsInformation");

   if( fixed_variable_treatment_ == MAKE_PARAMETER || fixed_variable_treatment_ == MAKE_PARAMETER_NODUAL )
   {
      // Set the values of fixed variables
      for( Index i = 0; i < n_x_fixed_; i++ )
      {
         DBG_ASSERT(x_l[x_fixed_map_[i]] == x_u[x_fixed_map_[i]]);
         full_x_[x_fixed_map_[i]] = x_l[x_fixed_map_[i]];
      }
   }
   else if( fixed_variable_treatment_ == RELAX_BOUNDS )
   {
      // Relax the bounds for fixed variables
      const Number bound_relax = Max(Number(1e-8), bound_relax_factor_);
      for( Index i = 0; i < n_x_fixed_; i++ )
      {
         if( x_l[i] == x_u[i] )
         {
            x_l[i] -= bound_relax * Max(Number(1.), std::abs(x_l[i]));
            x_u[i] += bound_relax * Max(Number(1.), std::abs(x_u[i]));
         }
      }
   }

   // Set the bounds values for x
   DenseVector* dx_L = static_cast<DenseVector*>(&x_L);
   DBG_ASSERT(dynamic_cast<DenseVector*>(&x_L));
   Number* values = dx_L->Values();
   const ExpansionMatrix* em_Px_L = static_cast<const ExpansionMatrix*>(&Px_L);
   DBG_ASSERT(dynamic_cast<const ExpansionMatrix*>(&Px_L));
   if( IsValid(P_x_full_x_) )
   {
      for( Index i = 0; i < Px_L.NCols(); i++ )
      {
         const Index& ipopt_idx = em_Px_L->ExpandedPosIndices()[i];
         const Index& full_idx = P_x_full_x_->ExpandedPosIndices()[ipopt_idx];
         const Number& lower_bound = x_l[full_idx];
         values[i] = lower_bound;
      }
   }
   else
   {
      for( Index i = 0; i < Px_L.NCols(); i++ )
      {
         const Index& ipopt_idx = em_Px_L->ExpandedPosIndices()[i];
         const Number& lower_bound = x_l[ipopt_idx];
         values[i] = lower_bound;
      }
   }

   DenseVector* dx_U = static_cast<DenseVector*>(&x_U);
   DBG_ASSERT(dynamic_cast<DenseVector*>(&x_U));
   values = dx_U->Values();
   const ExpansionMatrix* em_Px_U = static_cast<const ExpansionMatrix*>(&Px_U);
   DBG_ASSERT(dynamic_cast<const ExpansionMatrix*>(&Px_U));
   if( IsValid(P_x_full_x_) )
   {
      for( Index i = 0; i < Px_U.NCols(); i++ )
      {
         const Index& ipopt_idx = em_Px_U->ExpandedPosIndices()[i];
         const Index& full_idx = P_x_full_x_->ExpandedPosIndices()[ipopt_idx];
         const Number& upper_bound = x_u[full_idx];
         values[i] = upper_bound;
      }
   }
   else
   {
      for( Index i = 0; i < Px_U.NCols(); i++ )
      {
         const Index& ipopt_idx = em_Px_U->ExpandedPosIndices()[i];
         const Number& upper_bound = x_u[ipopt_idx];
         values[i] = upper_bound;
      }
   }

   // get the bounds values (rhs values to subtract) for c
   // i.e. if gL == gU, then we actually have g(x) = gL = gU,
   // since we solve c(x) = 0, we actually need c(x) - gL = 0
   for( Index i = 0; i < P_c_g_->NCols(); i++ )
   {
      Index full_idx = P_c_g_->ExpandedPosIndices()[i];
      Number rhs = g_l[full_idx];
      c_rhs_[i] = rhs;
   }
   // similarly, if we have fixed variables, consider them here
   if( fixed_variable_treatment_ == MAKE_CONSTRAINT )
   {
      Index n_c_no_fixed = P_c_g_->NCols();
      for( Index i = 0; i < n_x_fixed_; i++ )
      {
         DBG_ASSERT(x_l[x_fixed_map_[i]] == x_u[x_fixed_map_[i]]);
         c_rhs_[n_c_no_fixed + i] = x_l[x_fixed_map_[i]];
      }
   }

   // get the bounds values for d
   DenseVector* dd_L = static_cast<DenseVector*>(&d_L);
   DBG_ASSERT(dynamic_cast<DenseVector*>(&d_L));
   values = dd_L->Values();
   const ExpansionMatrix* em_Pd_L = static_cast<const ExpansionMatrix*>(&Pd_L);
   DBG_ASSERT(dynamic_cast<const ExpansionMatrix*>(&Pd_L));
   for( Index i = 0; i < Pd_L.NCols(); i++ )
   {
      Index d_exp_idx = em_Pd_L->ExpandedPosIndices()[i];
      Index full_idx = P_d_g_->ExpandedPosIndices()[d_exp_idx];
      Number lower_bound = g_l[full_idx];
      values[i] = lower_bound;
   }

   DenseVector* dd_U = static_cast<DenseVector*>(&d_U);
   DBG_ASSERT(dynamic_cast<DenseVector*>(&d_U));
   values = dd_U->Values();
   const ExpansionMatrix* em_Pd_U = static_cast<const ExpansionMatrix*>(&Pd_U);
   DBG_ASSERT(dynamic_cast<const ExpansionMatrix*>(&Pd_U));
   for( Index i = 0; i < Pd_U.NCols(); i++ )
   {
      Index d_exp_idx = em_Pd_U->ExpandedPosIndices()[i];
      Index full_idx = P_d_g_->ExpandedPosIndices()[d_exp_idx];
      Number upper_bound = g_u[full_idx];
      values[i] = upper_bound;
   }

   // In case we are doing finite differences, keep a copy of the bounds
   if( jacobian_approximation_ != JAC_EXACT || gradient_approximation_ != OBJGRAD_EXACT )
   {
      delete[] findiff_x_l_;
      delete[] findiff_x_u_;
      findiff_x_l_ = x_l;
      findiff_x_u_ = x_u;
      x_l = NULL;
      x_u = NULL;
   }

   delete[] x_l;
   x_l = NULL;
   delete[] x_u;
   x_u = NULL;
   delete[] g_l;
   g_l = NULL;
   delete[] g_u;
   g_u = NULL;

   return true;
}

bool TNLPAdapter::GetStartingPoint(
   SmartPtr<Vector> x,
   bool             need_x,
   SmartPtr<Vector> y_c,
   bool             need_y_c,
   SmartPtr<Vector> y_d,
   bool             need_y_d,
   SmartPtr<Vector> z_L,
   bool             need_z_L,
   SmartPtr<Vector> z_U,
   bool             need_z_U
)
{
   Number* full_x = new Number[n_full_x_];
   Number* full_z_l = new Number[n_full_x_];
   Number* full_z_u = new Number[n_full_x_];
   Number* full_lambda = new Number[n_full_g_];
   bool init_x = need_x;
   bool init_z = need_z_L || need_z_U;
   bool init_lambda = need_y_c || need_y_d;

   bool retvalue = tnlp_->get_starting_point(n_full_x_, init_x, full_x, init_z, full_z_l, full_z_u, n_full_g_,
                   init_lambda, full_lambda);

   if( !retvalue )
   {
      delete[] full_x;
      delete[] full_z_l;
      delete[] full_z_u;
      delete[] full_lambda;
      return false;
   }

   if( need_x )
   {
      DenseVector* dx = static_cast<DenseVector*>(GetRawPtr(x));
      DBG_ASSERT(dynamic_cast<DenseVector*>(GetRawPtr(x)));
      Number* values = dx->Values();
      const Index& n_x_var = x->Dim();
      if( IsValid(P_x_full_x_) )
      {
         const Index* x_pos = P_x_full_x_->ExpandedPosIndices();
         for( Index i = 0; i < n_x_var; i++ )
         {
            values[i] = full_x[x_pos[i]];
         }
      }
      else
      {
         IpBlasCopy(n_x_var, full_x, 1, values, 1);
      }
   }

   if( need_y_c )
   {
      DenseVector* dy_c = static_cast<DenseVector*>(GetRawPtr(y_c));
      DBG_ASSERT(dynamic_cast<DenseVector*>(GetRawPtr(y_c)));
      Number* values = dy_c->Values();
      const Index* y_c_pos = P_c_g_->ExpandedPosIndices();
      for( Index i = 0; i < P_c_g_->NCols(); i++ )
      {
         values[i] = full_lambda[y_c_pos[i]];
      }
      if( fixed_variable_treatment_ == MAKE_CONSTRAINT )
      {
         // ToDo maybe use info from z_L and Z_U here?
         const Number zero = 0.;
         IpBlasCopy(n_x_fixed_, &zero, 0, &values[P_c_g_->NCols()], 1);
      }
   }

   if( need_y_d )
   {
      DenseVector* dy_d = static_cast<DenseVector*>(GetRawPtr(y_d));
      DBG_ASSERT(dynamic_cast<DenseVector*>(GetRawPtr(y_d)));
      Number* values = dy_d->Values();
      const Index* y_d_pos = P_d_g_->ExpandedPosIndices();
      for( Index i = 0; i < y_d->Dim(); i++ )
      {
         values[i] = full_lambda[y_d_pos[i]];
      }
   }

   if( need_z_L )
   {
      DenseVector* dz_l = static_cast<DenseVector*>(GetRawPtr(z_L));
      DBG_ASSERT(dynamic_cast<DenseVector*>(GetRawPtr(z_L)));
      Number* values = dz_l->Values();
      const Index& n_z_l = z_L->Dim();
      const Index* z_l_pos = P_x_x_L_->ExpandedPosIndices();
      if( IsValid(P_x_full_x_) )
      {
         const Index* x_pos = P_x_full_x_->ExpandedPosIndices();
         for( Index i = 0; i < n_z_l; i++ )
         {
            Index idx = z_l_pos[i]; // convert from x_L to x (ipopt)
            idx = x_pos[idx]; // convert from x (ipopt) to x_full
            values[i] = full_z_l[idx];
         }
      }
      else
      {
         for( Index i = 0; i < n_z_l; i++ )
         {
            Index idx = z_l_pos[i]; // convert from x_L to x (ipopt)
            values[i] = full_z_l[idx];
         }
      }
   }

   if( need_z_U )
   {
      DenseVector* dz_u = static_cast<DenseVector*>(GetRawPtr(z_U));
      DBG_ASSERT(dynamic_cast<DenseVector*>(GetRawPtr(z_U)));
      Number* values = dz_u->Values();
      const Index* z_u_pos = P_x_x_U_->ExpandedPosIndices();
      if( IsValid(P_x_full_x_) )
      {
         const Index* x_pos = P_x_full_x_->ExpandedPosIndices();
         for( Index i = 0; i < z_U->Dim(); i++ )
         {
            Index idx = z_u_pos[i]; // convert from x_u to x (ipopt)
            idx = x_pos[idx]; // convert from x (ipopt) to x_full
            values[i] = full_z_u[idx];
         }
      }
      else
      {
         for( Index i = 0; i < z_U->Dim(); i++ )
         {
            Index idx = z_u_pos[i]; // convert from x_u to x (ipopt)
            values[i] = full_z_u[idx];
         }
      }
   }

   delete[] full_x;
   full_x = NULL;
   delete[] full_z_l;
   full_z_l = NULL;
   delete[] full_z_u;
   full_z_u = NULL;
   delete[] full_lambda;
   full_lambda = NULL;

   return true;
}

bool TNLPAdapter::GetWarmStartIterate(
   IteratesVector& warm_start_iterate
)
{
   return tnlp_->get_warm_start_iterate(warm_start_iterate);
}

bool TNLPAdapter::Eval_f(
   const Vector& x,
   Number&       f
)
{
   bool new_x = false;
   if( update_local_x(x) )
   {
      new_x = true;
   }
   return tnlp_->eval_f(n_full_x_, full_x_, new_x, f);
}

bool TNLPAdapter::Eval_grad_f(
   const Vector& x,
   Vector&       g_f
)
{
   bool retvalue = false;
   bool new_x = false;
   if( update_local_x(x) )
   {
      new_x = true;
   }

   DenseVector* dg_f = static_cast<DenseVector*>(&g_f);
   DBG_ASSERT(dynamic_cast<DenseVector*>(&g_f));
   Number* values = dg_f->Values();

   if( gradient_approximation_ == OBJGRAD_EXACT )
   {
      if( IsValid(P_x_full_x_) )
      {
         Number* full_grad_f = new Number[n_full_x_];
         if( tnlp_->eval_grad_f(n_full_x_, full_x_, new_x, full_grad_f) )
         {
            const Index* x_pos = P_x_full_x_->ExpandedPosIndices();
            for( Index i = 0; i < g_f.Dim(); i++ )
            {
               values[i] = full_grad_f[x_pos[i]];
            }
            retvalue = true;
         }
         delete[] full_grad_f;
      }
      else
      {
         retvalue = tnlp_->eval_grad_f(n_full_x_, full_x_, new_x, values);
      }
   }
   else
   {
      // make sure we have the value of the objective at the point
      Number f;
      retvalue = tnlp_->eval_f(n_full_x_, full_x_, new_x, f);
      if( retvalue )
      {
         Number* full_x_pert = new Number[n_full_x_];
         IpBlasCopy(n_full_x_, full_x_, 1, full_x_pert, 1);
         const Index* x_pos = NULL;
         if( IsValid(P_x_full_x_) )
         {
            x_pos = P_x_full_x_->ExpandedPosIndices();
         }

         // Compute the finite difference objective
         for( Index i = 0; i < g_f.Dim(); i++ )
         {
            Index ivar = x_pos != NULL ? x_pos[i] : i;
            if( findiff_x_l_[ivar] < findiff_x_u_[ivar] )
            {
               const Number xorig = full_x_pert[ivar];
               Number this_perturbation = findiff_perturbation_ * Max(Number(1.), std::abs(full_x_[ivar]));
               full_x_pert[ivar] += this_perturbation;
               if( full_x_pert[ivar] > findiff_x_u_[ivar] )
               {
                  // if at upper bound, then change direction towards lower bound
                  this_perturbation = -this_perturbation;
                  full_x_pert[ivar] = xorig + this_perturbation;
               }
               Number f_pert;
               retvalue = tnlp_->eval_f(n_full_x_, full_x_pert, true, f_pert);
               if( !retvalue )
               {
                  break;
               }

               values[i] = (f_pert - f) / this_perturbation;

               full_x_pert[ivar] = xorig;
            }
         }

         delete[] full_x_pert;
      }
   }

   return retvalue;
}

bool TNLPAdapter::Eval_c(
   const Vector& x,
   Vector&       c
)
{
   bool new_x = false;
   if( update_local_x(x) )
   {
      new_x = true;
   }

   if( internal_eval_g(new_x) )
   {
      DenseVector* dc = static_cast<DenseVector*>(&c);
      DBG_ASSERT(dynamic_cast<DenseVector*>(&c));
      Number* values = dc->Values();
      const Index* c_pos = P_c_g_->ExpandedPosIndices();
      Index n_c_no_fixed = P_c_g_->NCols();
      for( Index i = 0; i < n_c_no_fixed; i++ )
      {
         values[i] = full_g_[c_pos[i]];
         values[i] -= c_rhs_[i];
      }
      if( fixed_variable_treatment_ == MAKE_CONSTRAINT )
      {
         for( Index i = 0; i < n_x_fixed_; i++ )
         {
            values[n_c_no_fixed + i] = full_x_[x_fixed_map_[i]] - c_rhs_[n_c_no_fixed + i];
         }
      }
      return true;
   }

   return false;
}

bool TNLPAdapter::Eval_jac_c(
   const Vector& x,
   Matrix&       jac_c
)
{
   bool new_x = false;
   if( update_local_x(x) )
   {
      new_x = true;
   }

   if( internal_eval_jac_g(new_x) )
   {
      GenTMatrix* gt_jac_c = static_cast<GenTMatrix*>(&jac_c);
      DBG_ASSERT(dynamic_cast<GenTMatrix*>(&jac_c));
      Number* values = gt_jac_c->Values();

      for( Index i = 0; i < nz_jac_c_no_extra_; i++ )
      {
         // Assume the same structure as initially given
         values[i] = jac_g_[jac_idx_map_[i]];
      }
      if( fixed_variable_treatment_ == MAKE_CONSTRAINT )
      {
         const Number one = 1.;
         IpBlasCopy(n_x_fixed_, &one, 0, &values[nz_jac_c_no_extra_], 1);
      }
      return true;
   }
   return false;
}

bool TNLPAdapter::Eval_d(
   const Vector& x,
   Vector&       d
)
{
   bool new_x = false;
   if( update_local_x(x) )
   {
      new_x = true;
   }

   DenseVector* dd = static_cast<DenseVector*>(&d);
   DBG_ASSERT(dynamic_cast<DenseVector*>(&d));
   Number* values = dd->Values();
   if( internal_eval_g(new_x) )
   {
      const Index* d_pos = P_d_g_->ExpandedPosIndices();
      for( Index i = 0; i < d.Dim(); i++ )
      {
         values[i] = full_g_[d_pos[i]];
      }
      return true;
   }

   return false;
}

bool TNLPAdapter::Eval_jac_d(
   const Vector& x,
   Matrix&       jac_d
)
{
   bool new_x = false;
   if( update_local_x(x) )
   {
      new_x = true;
   }

   if( internal_eval_jac_g(new_x) )
   {
      GenTMatrix* gt_jac_d = static_cast<GenTMatrix*>(&jac_d);
      DBG_ASSERT(dynamic_cast<GenTMatrix*>(&jac_d));
      Number* values = gt_jac_d->Values();

      for( Index i = 0; i < nz_jac_d_; i++ )
      {
         // Assume the same structure as initially given
         values[i] = jac_g_[jac_idx_map_[nz_jac_c_no_extra_ + i]];
      }
      return true;
   }
   return false;
}

bool TNLPAdapter::Eval_h(
   const Vector& x,
   Number        obj_factor,
   const Vector& yc,
   const Vector& yd,
   SymMatrix&    h
)
{
   // First see if all weights are set to zero (for example, when
   // computing the least square multiplier estimates, this is what
   // we do).  In that case, there is no need to compute values, just
   // set them to zero.
   if( obj_factor == 0. && yc.Asum() == 0. && yd.Asum() == 0. )
   {
      SymTMatrix* st_h = static_cast<SymTMatrix*>(&h);
      DBG_ASSERT(dynamic_cast<SymTMatrix*>(&h));
      Number* values = st_h->Values();
      for( Index i = 0; i < nz_h_; i++ )
      {
         values[i] = 0.;
      }
      return true;
   }

   bool retval = false;
   bool new_x = false;
   if( update_local_x(x) )
   {
      new_x = true;
   }
   bool new_y = false;
   if( update_local_lambda(yc, yd) )
   {
      new_y = true;
   }

   SymTMatrix* st_h = static_cast<SymTMatrix*>(&h);
   DBG_ASSERT(dynamic_cast<SymTMatrix*>(&h));
   Number* values = st_h->Values();

   if( h_idx_map_ )
   {
      Number* full_h = new Number[nz_full_h_];

      if( tnlp_->eval_h(n_full_x_, full_x_, new_x, obj_factor, n_full_g_, full_lambda_, new_y, nz_full_h_, NULL, NULL,
                        full_h) )
      {
         for( Index i = 0; i < nz_h_; i++ )
         {
            values[i] = full_h[h_idx_map_[i]];
         }
         retval = true;
      }
      delete[] full_h;
   }
   else
   {
      retval = tnlp_->eval_h(n_full_x_, full_x_, new_x, obj_factor, n_full_g_, full_lambda_, new_y, nz_full_h_, NULL,
                             NULL, values);
   }

   return retval;
}

void TNLPAdapter::GetScalingParameters(
   const SmartPtr<const VectorSpace> x_space,
   const SmartPtr<const VectorSpace> c_space,
   const SmartPtr<const VectorSpace> d_space,
   Number&                           obj_scaling,
   SmartPtr<Vector>&                 x_scaling,
   SmartPtr<Vector>&                 c_scaling,
   SmartPtr<Vector>&                 d_scaling
) const
{
   x_scaling = x_space->MakeNew();
   c_scaling = c_space->MakeNew();
   d_scaling = d_space->MakeNew();
   DBG_ASSERT(c_scaling->Dim() + d_scaling->Dim() == n_full_g_ + (fixed_variable_treatment_ == MAKE_CONSTRAINT ? n_x_fixed_ : 0));

   DenseVector* dx = static_cast<DenseVector*>(GetRawPtr(x_scaling));
   DBG_ASSERT(dynamic_cast<DenseVector*>(GetRawPtr(x_scaling)));
   DenseVector* dc = static_cast<DenseVector*>(GetRawPtr(c_scaling));
   DBG_ASSERT(dynamic_cast<DenseVector*>(GetRawPtr(c_scaling)));
   DenseVector* dd = static_cast<DenseVector*>(GetRawPtr(d_scaling));
   DBG_ASSERT(dynamic_cast<DenseVector*>(GetRawPtr(d_scaling)));
   Number* dx_values = dx->Values();
   Number* dc_values = dc->Values();
   Number* dd_values = dd->Values();

   Number* full_g_scaling = new Number[n_full_g_];
   bool use_x_scaling = true;
   bool use_g_scaling = true;

   if( IsValid(P_x_full_x_) )
   {
      Number* full_x_scaling = new Number[n_full_x_];
      bool retval = tnlp_->get_scaling_parameters(obj_scaling, use_x_scaling, n_full_x_, full_x_scaling, use_g_scaling,
                    n_full_g_, full_g_scaling);
      if( !retval )
      {
         delete[] full_x_scaling;
         jnlst_->Printf(J_ERROR, J_INITIALIZATION,
                        "Option nlp_scaling_method selected as user-scaling, but no user-scaling available, or it cannot be computed.\n");
         THROW_EXCEPTION(OPTION_INVALID, "User scaling chosen, but get_scaling_parameters returned false.");
      }

      if( use_x_scaling )
      {
         const Index* x_pos = P_x_full_x_->ExpandedPosIndices();
         for( Index i = 0; i < dx->Dim(); i++ )
         {
            dx_values[i] = full_x_scaling[x_pos[i]];
         }
      }
      delete[] full_x_scaling;
   }
   else
   {
      bool retval = tnlp_->get_scaling_parameters(obj_scaling, use_x_scaling, n_full_x_, dx_values, use_g_scaling,
                    n_full_g_, full_g_scaling);
      if( !retval )
      {
         jnlst_->Printf(J_ERROR, J_INITIALIZATION,
                        "Option nlp_scaling_method selected as user-scaling, but no user-scaling available, or it cannot be computed.\n");
         THROW_EXCEPTION(OPTION_INVALID, "User scaling chosen, but get_scaling_parameters returned false.");
      }
   }

   if( !use_x_scaling )
   {
      x_scaling = NULL;
   }

   if( use_g_scaling )
   {
      const Index* c_pos = P_c_g_->ExpandedPosIndices();
      for( Index i = 0; i < P_c_g_->NCols(); i++ )
      {
         dc_values[i] = full_g_scaling[c_pos[i]];
      }
      if( fixed_variable_treatment_ == MAKE_CONSTRAINT )
      {
         const Number one = 1.;
         IpBlasCopy(n_x_fixed_, &one, 0, &dc_values[P_c_g_->NCols()], 1);
      }

      const Index* d_pos = P_d_g_->ExpandedPosIndices();
      for( Index i = 0; i < dd->Dim(); i++ )
      {
         dd_values[i] = full_g_scaling[d_pos[i]];
      }
   }
   else
   {
      c_scaling = NULL;
      d_scaling = NULL;
   }

   delete[] full_g_scaling;
}

void TNLPAdapter::FinalizeSolution(
   SolverReturn               status,
   const Vector&              x,
   const Vector&              z_L,
   const Vector&              z_U,
   const Vector&              c,
   const Vector&              d,
   const Vector&              y_c,
   const Vector&              y_d,
   Number                     obj_value,
   const IpoptData*           ip_data,
   IpoptCalculatedQuantities* ip_cq
)
{
   DBG_START_METH("TNLPAdapter::FinalizeSolution", dbg_verbosity);

   update_local_x(x);
   update_local_lambda(y_c, y_d);

   ResortX(x, full_x_);

   StringMetaDataMapType var_string_md;
   IntegerMetaDataMapType var_integer_md;
   NumericMetaDataMapType var_numeric_md;
   StringMetaDataMapType con_string_md;
   IntegerMetaDataMapType con_integer_md;
   NumericMetaDataMapType con_numeric_md;

   SmartPtr<const DenseVectorSpace> x_space = dynamic_cast<const DenseVectorSpace*>(GetRawPtr(x.OwnerSpace()));
   const NumericMetaDataMapType x_meta = x_space->GetNumericMetaData();
   NumericMetaDataMapType::const_iterator x_meta_iter;
   for( x_meta_iter = x_meta.begin(); x_meta_iter != x_meta.end(); ++x_meta_iter )
   {
      if( (Index) x_meta_iter->second.size() == x.Dim() )
      {
         std::vector<Number> new_meta_data;
         new_meta_data.resize(n_full_x_);
         SmartPtr<DenseVector> x_meta_vector = x_space->MakeNewDenseVector();
         x_meta_vector->SetValues(&(x_meta_iter->second)[0]);
         ResortX(*x_meta_vector, &new_meta_data[0]);
         var_numeric_md[x_meta_iter->first] = new_meta_data;
      }
   }

   ResortG(y_c, y_d, full_lambda_);

   SmartPtr<const DenseVectorSpace> y_c_space = dynamic_cast<const DenseVectorSpace*>(GetRawPtr(y_c.OwnerSpace()));
   SmartPtr<const DenseVectorSpace> y_d_space = dynamic_cast<const DenseVectorSpace*>(GetRawPtr(y_d.OwnerSpace()));
   const NumericMetaDataMapType y_c_meta = y_c_space->GetNumericMetaData();
   const NumericMetaDataMapType y_d_meta = y_d_space->GetNumericMetaData();
   NumericMetaDataMapType::const_iterator y_c_meta_iter;
   for( y_c_meta_iter = y_c_meta.begin(); y_c_meta_iter != y_c_meta.end(); ++y_c_meta_iter )
   {
      if( (Index) y_c_meta_iter->second.size() == y_c.Dim() )
      {
         if( y_d_space->HasNumericMetaData(y_c_meta_iter->first.c_str())    // There exists a corresponding y_d metadata
             && ((Index) (y_d_meta.find(y_c_meta_iter->first)->second).size()) == y_d.Dim() ) // and y_d metadata vector has size y_d.Dim()
         {
            std::vector<Number> y_d_second = y_d_space->GetNumericMetaData(y_c_meta_iter->first);
            std::vector<Number> new_g_meta_data;
            new_g_meta_data.resize(n_full_g_);
            SmartPtr<DenseVector> y_c_meta_vector = y_c_space->MakeNewDenseVector();
            SmartPtr<DenseVector> y_d_meta_vector = y_d_space->MakeNewDenseVector();
            y_c_meta_vector->SetValues(&(y_c_meta_iter->second)[0]);
            if( y_d_second.size() > 0 )
            {
               y_d_meta_vector->SetValues(&y_d_second[0]);
            }
            ResortG(*y_c_meta_vector, *y_d_meta_vector, &new_g_meta_data[0]);
            con_numeric_md[y_c_meta_iter->first] = new_g_meta_data;
         }
      }
   }

   Number* full_g = new Number[n_full_g_];
   // TODO:
   if( c.Dim() + d.Dim() < n_full_g_ )
   {
      const Number zero = 0.;
      IpBlasCopy(n_full_g_, &zero, 0, full_g, 1);
   }
   ResortG(c, d, full_g, true);

   Number* full_z_L = new Number[n_full_x_];
   Number* full_z_U = new Number[n_full_x_];
   switch( status )
   {
      case SUCCESS:
      case MAXITER_EXCEEDED:
      case STOP_AT_TINY_STEP:
      case STOP_AT_ACCEPTABLE_POINT:
      case LOCAL_INFEASIBILITY:
      case USER_REQUESTED_STOP:
      case FEASIBLE_POINT_FOUND:
      case DIVERGING_ITERATES:
      case RESTORATION_FAILURE:
      case ERROR_IN_STEP_COMPUTATION:
         // only for these status codes, IpoptApplication calls without all vectors 0
         if( !ResortBoundMultipliers(x, y_c, y_d, z_L, full_z_L, z_U, full_z_U) )
         {
            jnlst_->Printf(J_WARNING, J_INITIALIZATION, "Failed to evaluate gradient of objective or constraints when computing bound multipliers for fixed variables.\n");
         }
         break;
      default:
         // if IpoptApplication doesn't provide an actual solution, do not bother to setup good multipliers for fixed variables
         ResortBounds(z_L, full_z_L, z_U, full_z_U);
         break;
   }

   SmartPtr<const DenseVectorSpace> z_L_space = dynamic_cast<const DenseVectorSpace*>(GetRawPtr(z_L.OwnerSpace()));
   SmartPtr<const DenseVectorSpace> z_U_space = dynamic_cast<const DenseVectorSpace*>(GetRawPtr(z_U.OwnerSpace()));
   const NumericMetaDataMapType z_L_meta = z_L_space->GetNumericMetaData();
   const NumericMetaDataMapType z_U_meta = z_U_space->GetNumericMetaData();
   NumericMetaDataMapType::const_iterator z_L_meta_iter;
   for( z_L_meta_iter = z_L_meta.begin(); z_L_meta_iter != z_L_meta.end(); ++z_L_meta_iter )
   {
      if( (Index) z_L_meta_iter->second.size() == z_L.Dim() )
      {
         if( z_U_space->HasNumericMetaData(z_L_meta_iter->first.c_str())
             && (Index) z_U_meta.find(z_L_meta_iter->first.c_str())->second.size() == z_U.Dim() )
         {
            std::vector<Number> z_U_second = z_U_space->GetNumericMetaData(z_L_meta_iter->first);
            SmartPtr<DenseVector> z_L_meta_vector = z_L_space->MakeNewDenseVector();
            SmartPtr<DenseVector> z_U_meta_vector = z_U_space->MakeNewDenseVector();
            z_L_meta_vector->SetValues(&(z_L_meta_iter->second)[0]);
            if( z_U_second.size() > 0 )
            {
               z_U_meta_vector->SetValues(&z_U_second[0]);
            }
            std::vector<Number> new_z_L_meta_data(n_full_x_, 0.0);
            std::vector<Number> new_z_U_meta_data(n_full_x_, 0.0);
            ResortBounds(*z_L_meta_vector, &new_z_L_meta_data[0], *z_U_meta_vector, &new_z_U_meta_data[0]);
            std::string z_L_meta_data_tag = z_L_meta_iter->first;
            std::string z_U_meta_data_tag = z_L_meta_iter->first;
            z_L_meta_data_tag += "_z_L";
            z_U_meta_data_tag += "_z_U";
            var_numeric_md[z_L_meta_data_tag] = new_z_L_meta_data;
            var_numeric_md[z_U_meta_data_tag] = new_z_U_meta_data;
         }
      }
   }

   tnlp_->finalize_metadata(n_full_x_, var_string_md, var_integer_md, var_numeric_md, n_full_g_, con_string_md,
                            con_integer_md, con_numeric_md);

   tnlp_->finalize_solution(status, n_full_x_, full_x_, full_z_L, full_z_U, n_full_g_, full_g, full_lambda_, obj_value,
                            ip_data, ip_cq);

   delete[] full_z_L;
   full_z_L = NULL;
   delete[] full_z_U;
   full_z_U = NULL;
   delete[] full_g;
   full_g = NULL;

   if( c.Dim() + d.Dim() < n_full_g_ )
   {
      // Temporary: Check if we really have a feasible point:
      Number max_viol = 0.;
      Number* x_l = new Number[n_full_x_];
      Number* x_u = new Number[n_full_x_];
      Number* g_l = new Number[n_full_g_];
      Number* g_u = new Number[n_full_g_];
      bool retval = tnlp_->get_bounds_info(n_full_x_, x_l, x_u, n_full_g_, g_l, g_u);
      ASSERT_EXCEPTION(retval, INVALID_TNLP, "get_bounds_info returned false in FinalizeSolution");
      for( Index i = 0; i < n_full_g_; i++ )
      {
         max_viol = Max(max_viol, full_g_[i] - g_u[i], g_l[i] - full_g_[i]);
      }
      jnlst_->Printf(J_ITERSUMMARY, J_INITIALIZATION, "Constraint violation for ALL constraints is %e.\n", max_viol);
      delete[] x_l;
      delete[] x_u;
      delete[] g_l;
      delete[] g_u;
   }
}

bool TNLPAdapter::IntermediateCallBack(
   AlgorithmMode              mode,
   Index                      iter,
   Number                     obj_value,
   Number                     inf_pr,
   Number                     inf_du,
   Number                     mu,
   Number                     d_norm,
   Number                     regularization_size,
   Number                     alpha_du,
   Number                     alpha_pr,
   Index                      ls_trials,
   const IpoptData*           ip_data,
   IpoptCalculatedQuantities* ip_cq
)
{
   return tnlp_->intermediate_callback(mode, iter, obj_value, inf_pr, inf_du, mu, d_norm, regularization_size, alpha_du,
                                       alpha_pr, ls_trials, ip_data, ip_cq);
}

void TNLPAdapter::GetQuasiNewtonApproximationSpaces(
   SmartPtr<VectorSpace>& approx_space,
   SmartPtr<Matrix>&      P_approx
)
{
   Index num_nonlin_vars = tnlp_->get_number_of_nonlinear_variables();

   if( num_nonlin_vars < 0 && num_linear_variables_ == 0 )
   {
      approx_space = NULL;
      P_approx = NULL;
      return;
   }

   Index* pos_nonlin_vars = NULL;
   if( num_nonlin_vars < 0 )
   {
      num_nonlin_vars = n_full_x_ - num_linear_variables_;
      pos_nonlin_vars = new Index[num_nonlin_vars];
      Index ii = 0;
      for( Index i = num_linear_variables_; i < n_full_x_; i++ )
      {
         pos_nonlin_vars[ii++] = i;
      }
   }
   else if( num_nonlin_vars > 0 )
   {
      pos_nonlin_vars = new Index[num_nonlin_vars];
      bool retval = tnlp_->get_list_of_nonlinear_variables(num_nonlin_vars, pos_nonlin_vars);
      if( !retval )
      {
         delete[] pos_nonlin_vars;
         jnlst_->Printf(J_ERROR, J_INITIALIZATION,
                        "TNLP's get_number_of_nonlinear_variables returns non-negative number, but get_list_of_nonlinear_variables returns false.\n");
         THROW_EXCEPTION(INVALID_TNLP, "get_list_of_nonlinear_variables has not been overwritten");
      }
      // Correct indices in case user starts counting variables at 1
      // and not 0
      if( index_style_ == TNLP::FORTRAN_STYLE )
      {
         for( Index i = 0; i < num_nonlin_vars; i++ )
         {
            DBG_ASSERT(pos_nonlin_vars[i] > 0);
            DBG_ASSERT(pos_nonlin_vars[i] <= n_full_x_);
            pos_nonlin_vars[i]--;
         }
      }
#if IPOPT_CHECKLEVEL > 0
      else
      {
         for( Index i = 0; i < num_nonlin_vars; i++ )
         {
            DBG_ASSERT(pos_nonlin_vars[i] >= 0);
            DBG_ASSERT(pos_nonlin_vars[i] < n_full_x_);
         }
      }
#endif
   }

   if( IsNull(P_x_full_x_) )
   {
      if( num_nonlin_vars == n_full_x_ )
      {
         approx_space = NULL;
         P_approx = NULL;
      }
      else
      {
         SmartPtr<ExpansionMatrixSpace> ex_sp = new ExpansionMatrixSpace(n_full_x_, num_nonlin_vars, pos_nonlin_vars);
         P_approx = ex_sp->MakeNew();
         approx_space = new DenseVectorSpace(num_nonlin_vars);
      }
   }
   else
   {
      const Index* compr_pos = P_x_full_x_->CompressedPosIndices();
      Index* nonfixed_pos_nonlin_vars = new Index[num_nonlin_vars];

      Index nonfixed_nonlin_vars = 0;
      for( Index i = 0; i < num_nonlin_vars; i++ )
      {
         Index full_pos = pos_nonlin_vars[i];
         Index nonfixed_pos = compr_pos[full_pos];
         if( nonfixed_pos >= 0 )
         {
            nonfixed_pos_nonlin_vars[nonfixed_nonlin_vars] = nonfixed_pos;
            nonfixed_nonlin_vars++;
         }
      }

      const Index n_x_free = n_full_x_ - n_x_fixed_;
      if( nonfixed_nonlin_vars == n_x_free )
      {
         approx_space = NULL;
         P_approx = NULL;
      }
      else
      {
         SmartPtr<ExpansionMatrixSpace> ex_sp = new ExpansionMatrixSpace(n_x_free, nonfixed_nonlin_vars,
               nonfixed_pos_nonlin_vars);
         P_approx = ex_sp->MakeNew();
         approx_space = new DenseVectorSpace(nonfixed_nonlin_vars);
      }

      delete[] nonfixed_pos_nonlin_vars;
   }
   delete[] pos_nonlin_vars;
}

void TNLPAdapter::ResortX(
   const Vector& x,
   Number*       x_orig,
   bool          usefixedvals
)
{
   const DenseVector* dx = static_cast<const DenseVector*>(&x);
   DBG_ASSERT(dynamic_cast<const DenseVector*>(&x));

   if( IsValid(P_x_full_x_) )
   {
      const Index* x_pos = P_x_full_x_->CompressedPosIndices();

      if( dx->IsHomogeneous() )
      {
         const Number& scalar = dx->Scalar();
         for( Index i = 0; i < n_full_x_; i++ )
         {
            Index idx = x_pos[i];
            if( idx != -1 )
            {
               x_orig[i] = scalar;
            }
            else if( usefixedvals )
            {
               x_orig[i] = full_x_[i];
            }
            else
            {
               x_orig[i] = 0.0;
            }
         }
      }
      else
      {
         const Number* x_values = dx->Values();
         for( Index i = 0; i < n_full_x_; i++ )
         {
            Index idx = x_pos[i];
            if( idx != -1 )
            {
               x_orig[i] = x_values[idx];
            }
            else if( usefixedvals )
            {
               x_orig[i] = full_x_[i];
            }
            else
            {
               x_orig[i] = 0.0;
            }
         }
      }
   }
   else
   {
      if( dx->IsHomogeneous() )
      {
         const Number& scalar = dx->Scalar();
         IpBlasCopy(n_full_x_, &scalar, 0, x_orig, 1);
      }
      else
      {
         IpBlasCopy(n_full_x_, dx->Values(), 1, x_orig, 1);
      }
   }
}

void TNLPAdapter::ResortG(
   const Vector& c,
   const Vector& d,
   Number*       g_orig,
   bool          correctrhs
)
{
   const DenseVector* dc = static_cast<const DenseVector*>(&c);
   DBG_ASSERT(dynamic_cast<const DenseVector*>(&c));

   const Index* c_pos = P_c_g_->ExpandedPosIndices();
   if( dc->IsHomogeneous() )
   {
      Number scalar = dc->Scalar();
      for( Index i = 0; i < P_c_g_->NCols(); i++ )
      {
         g_orig[c_pos[i]] = scalar;
         if( correctrhs )
         {
            g_orig[c_pos[i]] += c_rhs_[i];
         }
      }
   }
   else
   {
      const Number* c_values = dc->Values();
      for( Index i = 0; i < P_c_g_->NCols(); i++ )
      {
         g_orig[c_pos[i]] = c_values[i];
         if( correctrhs )
         {
            g_orig[c_pos[i]] += c_rhs_[i];
         }
      }
   }

   const DenseVector* dd = static_cast<const DenseVector*>(&d);
   DBG_ASSERT(dynamic_cast<const DenseVector*>(&d));

   const Index* d_pos = P_d_g_->ExpandedPosIndices();
   if( dd->IsHomogeneous() )
   {
      Number scalar = dd->Scalar();
      for( Index i = 0; i < d.Dim(); i++ )
      {
         g_orig[d_pos[i]] = scalar;
      }
   }
   else
   {
      const Number* d_values = dd->Values();
      for( Index i = 0; i < d.Dim(); i++ )
      {
         g_orig[d_pos[i]] = d_values[i];
      }
   }
}

void TNLPAdapter::ResortBounds(
   const Vector& x_L,
   Number*       x_L_orig,
   const Vector& x_U,
   Number*       x_U_orig
)
{
   if( x_L_orig )
   {
      const DenseVector* dx_L = static_cast<const DenseVector*>(&x_L);
      DBG_ASSERT(dynamic_cast<const DenseVector*>(&x_L));

      const Index* bnds_pos_not_fixed = P_x_x_L_->ExpandedPosIndices();
      const Index& n_xL = x_L.Dim();

      if( n_xL < n_full_x_ )
      {
         memset(x_L_orig, 0, n_full_x_ * sizeof(Number));
      }

      if( IsValid(P_x_full_x_) )
      {
         const Index* bnds_pos_full = P_x_full_x_->ExpandedPosIndices();
         if( dx_L->IsHomogeneous() )
         {
            Number scalar = dx_L->Scalar();
            for( Index i = 0; i < n_xL; i++ )
            {
               Index idx = bnds_pos_not_fixed[i];
               idx = bnds_pos_full[idx];
               x_L_orig[idx] = scalar;
            }
         }
         else
         {
            const Number* x_L_values = dx_L->Values();
            for( Index i = 0; i < n_xL; i++ )
            {
               Index idx = bnds_pos_not_fixed[i];
               idx = bnds_pos_full[idx];
               x_L_orig[idx] = x_L_values[i];
            }
         }
      }
      else
      {
         if( dx_L->IsHomogeneous() )
         {
            Number scalar = dx_L->Scalar();
            for( Index i = 0; i < n_xL; i++ )
            {
               Index idx = bnds_pos_not_fixed[i];
               x_L_orig[idx] = scalar;
            }
         }
         else
         {
            const Number* x_L_values = dx_L->Values();
            for( Index i = 0; i < n_xL; i++ )
            {
               Index idx = bnds_pos_not_fixed[i];
               x_L_orig[idx] = x_L_values[i];
            }
         }
      }
   }

   if( x_U_orig )
   {
      const DenseVector* dx_U = static_cast<const DenseVector*>(&x_U);
      DBG_ASSERT(dynamic_cast<const DenseVector*>(&x_U));

      if( x_U.Dim() < n_full_x_ )
      {
         memset(x_U_orig, 0, n_full_x_ * sizeof(Number));
      }

      const Index* bnds_pos_not_fixed = P_x_x_U_->ExpandedPosIndices();

      if( IsValid(P_x_full_x_) )
      {
         const Index* bnds_pos_full = P_x_full_x_->ExpandedPosIndices();
         if( dx_U->IsHomogeneous() )
         {
            Number scalar = dx_U->Scalar();
            for( Index i = 0; i < x_U.Dim(); i++ )
            {
               Index idx = bnds_pos_not_fixed[i];
               idx = bnds_pos_full[idx];
               x_U_orig[idx] = scalar;
            }
         }
         else
         {
            const Number* x_U_values = dx_U->Values();
            for( Index i = 0; i < x_U.Dim(); i++ )
            {
               Index idx = bnds_pos_not_fixed[i];
               idx = bnds_pos_full[idx];
               x_U_orig[idx] = x_U_values[i];
            }
         }
      }
      else
      {
         if( dx_U->IsHomogeneous() )
         {
            Number scalar = dx_U->Scalar();
            for( Index i = 0; i < x_U.Dim(); i++ )
            {
               Index idx = bnds_pos_not_fixed[i];
               x_U_orig[idx] = scalar;
            }
         }
         else
         {
            const Number* x_U_values = dx_U->Values();
            for( Index i = 0; i < x_U.Dim(); i++ )
            {
               Index idx = bnds_pos_not_fixed[i];
               x_U_orig[idx] = x_U_values[i];
            }
         }
      }
   }
}

bool TNLPAdapter::ResortBoundMultipliers(
   const Vector& x,
   const Vector& y_c,
   const Vector& y_d,
   const Vector& z_L,
   Number*       z_L_orig,
   const Vector& z_U,
   Number*       z_U_orig
)
{
   ResortBounds(z_L, z_L_orig, z_U, z_U_orig);

   if( n_x_fixed_ == 0 )
   {
      return true;
   }

   // recover the bound multipliers for fixed variables
   if( fixed_variable_treatment_ == MAKE_CONSTRAINT )
   {
      const DenseVector* dy_c = static_cast<const DenseVector*>(&y_c);
      DBG_ASSERT(dynamic_cast<const DenseVector*>(&y_c));
      Index n_c_no_fixed = y_c.Dim() - n_x_fixed_;
      if( !dy_c->IsHomogeneous() )
      {
         const Number* values = dy_c->Values();
         for( Index i = 0; i < n_x_fixed_; i++ )
         {
            if( z_L_orig != NULL )
            {
               z_L_orig[x_fixed_map_[i]] = Max(Number(0.), -values[n_c_no_fixed + i]);
            }
            if( z_U_orig != NULL )
            {
               z_U_orig[x_fixed_map_[i]] = Max(Number(0.), values[n_c_no_fixed + i]);
            }
         }
      }
      else
      {
         Number value = dy_c->Scalar();
         for( Index i = 0; i < n_x_fixed_; i++ )
         {
            if( z_L_orig != NULL )
            {
               z_L_orig[x_fixed_map_[i]] = Max(Number(0.), -value);
            }
            if( z_U_orig != NULL )
            {
               z_U_orig[x_fixed_map_[i]] = Max(Number(0.), value);
            }
         }
      }
   }

   if( fixed_variable_treatment_ == MAKE_PARAMETER )
   {
      // Lagrangian should be grad_f + lambda jac - z_L + z_U  == 0
      // so fixed variables get z_L - z_U = grad_f + lambda_jac
      Number* mult = new Number[n_full_x_];
      memset(mult, 0, sizeof(Number) * n_full_x_);

      bool new_x = update_local_x(x);
      if( !tnlp_->eval_grad_f(n_full_x_, full_x_, new_x, mult) )
      {
         delete[] mult;
         return false;
      }

      if( !jac_fixed_idx_map_.empty() )
      {
         if( !internal_eval_jac_g(false) )
         {
            delete[] mult;
            return false;
         }
         DBG_ASSERT(dynamic_cast<const DenseVector*>(&y_c));
         DBG_ASSERT(dynamic_cast<const DenseVector*>(&y_d));
         const DenseVector* dy_c = static_cast<const DenseVector*>(&y_c);
         const DenseVector* dy_d = static_cast<const DenseVector*>(&y_d);
         // mappings from full g() indices to index in c() and d()
         const Index* c_row_pos = P_c_g_->CompressedPosIndices();
         const Index* d_row_pos = P_d_g_->CompressedPosIndices();
         for( size_t i = 0; i < jac_fixed_idx_map_.size(); i++ )
         {
            // Assume the same structure as initially given
            // correct for 1-based indexing in jac_fixed_iRow_ and jac_fixed_jCol_
            Index row = jac_fixed_iRow_[i] - 1;
            Index col = jac_fixed_jCol_[i] - 1;
            Number val = jac_g_[jac_fixed_idx_map_[i]];
            DBG_ASSERT(row >= 0);
            DBG_ASSERT(row < n_full_g_);
            DBG_ASSERT(col >= 0);
            DBG_ASSERT(col < n_full_x_);

            Number lambda = 0.0;
            if( c_row_pos[row] != -1 )
            {
               if( dy_c->IsHomogeneous() )
               {
                  lambda = dy_c->Scalar();
               }
               else
               {
                  lambda = dy_c->Values()[c_row_pos[row]];
               }
            }
            else if( d_row_pos[row] != -1 )
            {
               if( dy_d->IsHomogeneous() )
               {
                  lambda = dy_d->Scalar();
               }
               else
               {
                  lambda = dy_d->Values()[d_row_pos[row]];
               }
            }
            // else: a constraint that is neither in c() nor d(), so assuming lambda=0 seems ok

            // add lambda*jac part
            mult[col] += lambda * val;
         }
      }

      // set z_L = max(0,mult), z_U = max(0,-mult)
      for( Index i = 0; i < n_x_fixed_; ++i )
      {
         Index xidx = x_fixed_map_[i];
         if( z_L_orig != NULL )
         {
            z_L_orig[xidx] = Max(Number(0.0),  mult[xidx]);
         }
         if( z_U_orig != NULL )
         {
            z_U_orig[xidx] = Max(Number(0.0), -mult[xidx]);
         }
      }

      delete[] mult;
   }

   return true;
}

bool TNLPAdapter::update_local_x(
   const Vector& x
)
{
   if( x.GetTag() == x_tag_for_iterates_ )
   {
      return false;
   }

   ResortX(x, full_x_);

   x_tag_for_iterates_ = x.GetTag();

   return true;
}

bool TNLPAdapter::update_local_lambda(
   const Vector& y_c,
   const Vector& y_d
)
{
   if( y_c.GetTag() == y_c_tag_for_iterates_ && y_d.GetTag() == y_d_tag_for_iterates_ )
   {
      return false;
   }

   ResortG(y_c, y_d, full_lambda_);

   y_c_tag_for_iterates_ = y_c.GetTag();
   y_d_tag_for_iterates_ = y_d.GetTag();

   return true;
}

bool TNLPAdapter::internal_eval_g(
   bool new_x
)
{
   if( x_tag_for_g_ == x_tag_for_iterates_ )
   {
      // already calculated!
      return true;
   }

   x_tag_for_g_ = x_tag_for_iterates_;

   bool retval = tnlp_->eval_g(n_full_x_, full_x_, new_x, n_full_g_, full_g_);

   if( !retval )
   {
      x_tag_for_jac_g_ = 0;
   }

   return retval;
}

bool TNLPAdapter::internal_eval_jac_g(
   bool new_x
)
{
   if( x_tag_for_jac_g_ == x_tag_for_iterates_ )
   {
      // already calculated!
      return true;
   }

   x_tag_for_jac_g_ = x_tag_for_iterates_;

   bool retval;
   if( jacobian_approximation_ == JAC_EXACT )
   {
      retval = tnlp_->eval_jac_g(n_full_x_, full_x_, new_x, n_full_g_, nz_full_jac_g_, NULL, NULL, jac_g_);
   }
   else
   {
      // make sure we have the value of the constraints at the point
      retval = internal_eval_g(new_x);
      if( retval )
      {
         Number* full_g_pert = new Number[n_full_g_];
         Number* full_x_pert = new Number[n_full_x_];
         IpBlasCopy(n_full_x_, full_x_, 1, full_x_pert, 1);
         // Compute the finite difference Jacobian
         for( Index ivar = 0; ivar < n_full_x_; ivar++ )
         {
            if( findiff_x_l_[ivar] < findiff_x_u_[ivar] )
            {
               const Number xorig = full_x_pert[ivar];
               Number this_perturbation = findiff_perturbation_ * Max(Number(1.), std::abs(full_x_[ivar]));
               full_x_pert[ivar] += this_perturbation;
               if( full_x_pert[ivar] > findiff_x_u_[ivar] )
               {
                  // if at upper bound, then change direction towards lower bound
                  this_perturbation = -this_perturbation;
                  full_x_pert[ivar] = xorig + this_perturbation;
               }
               retval = tnlp_->eval_g(n_full_x_, full_x_pert, true, n_full_g_, full_g_pert);
               if( !retval )
               {
                  break;
               }
               for( Index i = findiff_jac_ia_[ivar]; i < findiff_jac_ia_[ivar + 1]; i++ )
               {
                  const Index& icon = findiff_jac_ja_[i];
                  const Index& ipos = findiff_jac_postriplet_[i];
                  jac_g_[ipos] = (full_g_pert[icon] - full_g_[icon]) / this_perturbation;
               }
               full_x_pert[ivar] = xorig;
            }
         }
         delete[] full_g_pert;
         delete[] full_x_pert;
      }
   }

   if( !retval )
   {
      x_tag_for_jac_g_ = 0;
   }

   return retval;
}

void TNLPAdapter::initialize_findiff_jac(
   const Index* iRow,
   const Index* jCol
)
{

   SmartPtr<TripletToCSRConverter> findiff_jac_converter = new TripletToCSRConverter(0);
   // construct structure of a sparse matrix with only Jacobian in it
   // TODO: This could be done more efficiently without detour via
   // symmetric matrix
   Index* airn = new Index[nz_full_jac_g_];
   Index* ajcn = new Index[nz_full_jac_g_];
   for( Index i = 0; i < nz_full_jac_g_; i++ )
   {
      airn[i] = jCol[i];
      ajcn[i] = iRow[i] + n_full_x_;
   }
   // Get the column ordered sparse representation
   findiff_jac_nnz_ = findiff_jac_converter->InitializeConverter(n_full_g_ + n_full_x_, nz_full_jac_g_, airn, ajcn);
   delete[] airn;
   delete[] ajcn;
   if( findiff_jac_nnz_ != nz_full_jac_g_ )
   {
      THROW_EXCEPTION(INVALID_TNLP,
                      "Sparsity structure of Jacobian has multiple occurrences of the same position.  This is not allowed for finite differences.");
   }

   // Finally, get the right numbers out of the converter object
   delete[] findiff_jac_ia_;
   delete[] findiff_jac_ja_;
   delete[] findiff_jac_postriplet_;
   findiff_jac_ia_ = NULL;
   findiff_jac_ja_ = NULL;
   findiff_jac_postriplet_ = NULL;
   findiff_jac_ia_ = new Index[n_full_x_ + 1];
   findiff_jac_ja_ = new Index[findiff_jac_nnz_];
   findiff_jac_postriplet_ = new Index[findiff_jac_nnz_];
   const Index* ia = findiff_jac_converter->IA();
   if( ia != NULL )
   {
      for( Index i = 0; i < n_full_x_ + 1; i++ )
      {
         findiff_jac_ia_[i] = ia[i];
      }
   }
   else
   {
      memset(findiff_jac_ia_, 0, n_full_x_ * sizeof(Index));
   }
   const Index* ja = findiff_jac_converter->JA();
   for( Index i = 0; i < findiff_jac_nnz_; i++ )
   {
      findiff_jac_ja_[i] = ja[i] - n_full_x_;
   }
   const Index* postrip = findiff_jac_converter->iPosFirst();
   for( Index i = 0; i < findiff_jac_nnz_; i++ )
   {
      findiff_jac_postriplet_[i] = postrip[i];
   }

}

bool TNLPAdapter::CheckDerivatives(
   TNLPAdapter::DerivativeTestEnum deriv_test,
   Index                           deriv_test_start_index
)
{
   if( deriv_test == NO_TEST )
   {
      return true;
   }

   Index nerrors = 0;

   ASSERT_EXCEPTION(IsValid(jnlst_), ERROR_IN_TNLP_DERIVATIVE_TEST,
                    "No Journalist given to TNLPAdapter.  Need Journalist, otherwise can't produce any output in DerivativeChecker!");

   bool retval = true;
   // Since this method should be independent of all other internal
   // data (so that it can be called indpenendent of GetSpace etc),
   // we are not using any internal fields

   // Obtain the problem size
   Index nx; // number of variables
   Index ng; // number of constriants
   Index nz_jac_g; // number of nonzeros in constraint Jacobian
   Index nz_hess_lag; // number of nonzeros in Lagrangian Hessian
   TNLP::IndexStyleEnum index_style;
   retval = tnlp_->get_nlp_info(nx, ng, nz_jac_g, nz_hess_lag, index_style);
   ASSERT_EXCEPTION(retval, INVALID_TNLP, "get_nlp_info returned false for derivative checker");
   ASSERT_EXCEPTION(nx >= 0, INVALID_TNLP, "number of variables negative");
   ASSERT_EXCEPTION(ng >= 0, INVALID_TNLP, "number of constraints negative");
   ASSERT_EXCEPTION(nz_jac_g >= 0, INVALID_TNLP, "number of nonzeros in Jacobian negative");
   ASSERT_EXCEPTION(nz_hess_lag >= 0, INVALID_TNLP, "number of nonzeros in Hessian negative");

   // Obtain starting point as reference point at which derivative
   // test should be performed
   Number* xref = new Number[nx];
   tnlp_->get_starting_point(nx, true, xref, false, NULL, NULL, ng, false, NULL);

   // Perform a random perturbation.  We need the bounds to make sure
   // they are not violated
   Number* x_l = new Number[nx];
   Number* x_u = new Number[nx];
   Number* g_l = new Number[ng];
   Number* g_u = new Number[ng];
   retval = tnlp_->get_bounds_info(nx, x_l, x_u, ng, g_l, g_u);
   ASSERT_EXCEPTION(retval, INVALID_TNLP, "get_bounds_info returned false in derivative checker");
   IpResetRandom01();
   for( Index i = 0; i < nx; i++ )
   {
      const Number lower = Max(x_l[i], xref[i] - point_perturbation_radius_);
      const Number upper = Min(x_u[i], xref[i] + point_perturbation_radius_);
      const Number interval = upper - lower;
      const Number random_number = IpRandom01();
      xref[i] = lower + random_number * interval;
   }
   delete[] x_l;
   delete[] x_u;
   delete[] g_l;
   delete[] g_u;

   // Obtain value of objective and constraints at reference point
   bool new_x = true;
   Number fref;
   Number* gref = NULL;
   if( ng > 0 )
   {
      gref = new Number[ng];
   }
   retval = tnlp_->eval_f(nx, xref, new_x, fref);
   ASSERT_EXCEPTION(retval, ERROR_IN_TNLP_DERIVATIVE_TEST,
                    "In TNLP derivative test: f could not be evaluated at reference point.");
   new_x = false;
   if( ng > 0 )
   {
      retval = tnlp_->eval_g(nx, xref, new_x, ng, gref);
      ASSERT_EXCEPTION(retval, ERROR_IN_TNLP_DERIVATIVE_TEST,
                       "In TNLP derivative test: g could not be evaluated at reference point.");
   }

   // Obtain gradient of objective function at reference point
   Number* grad_f = new Number[nx];
   retval = tnlp_->eval_grad_f(nx, xref, true, grad_f);
   ASSERT_EXCEPTION(retval, ERROR_IN_TNLP_DERIVATIVE_TEST,
                    "In TNLP derivative test: grad_f could not be evaluated at reference point.");

   Index* g_iRow = NULL;
   Index* g_jCol = NULL;
   Number* jac_g = NULL;
   if( ng > 0 )
   {
      // Obtain constraint Jacobian at reference point (including structure)
      g_iRow = new Index[nz_jac_g];
      g_jCol = new Index[nz_jac_g];
      retval = tnlp_->eval_jac_g(nx, NULL, false, ng, nz_jac_g, g_iRow, g_jCol, NULL);
      ASSERT_EXCEPTION(retval, ERROR_IN_TNLP_DERIVATIVE_TEST,
                       "In TNLP derivative test: Jacobian structure could not be evaluated.");
      // Correct counting if required to C-style
      if( index_style == TNLP::FORTRAN_STYLE )
      {
         for( Index i = 0; i < nz_jac_g; i++ )
         {
            DBG_ASSERT(g_iRow[i] > 0);
            DBG_ASSERT(g_iRow[i] <= ng);
            DBG_ASSERT(g_jCol[i] > 0);
            DBG_ASSERT(g_jCol[i] <= nx);
            g_iRow[i] -= 1;
            g_jCol[i] -= 1;
         }
      }
#if IPOPT_CHECKLEVEL > 0
      else
      {
         for( Index i = 0; i < nz_jac_g; i++ )
         {
            DBG_ASSERT(g_iRow[i] >= 0);
            DBG_ASSERT(g_iRow[i] < ng);
            DBG_ASSERT(g_jCol[i] >= 0);
            DBG_ASSERT(g_jCol[i] < nx);
         }
      }
#endif
      // Obtain values at reference point
      jac_g = new Number[nz_jac_g];
      retval = tnlp_->eval_jac_g(nx, xref, new_x, ng, nz_jac_g, NULL, NULL, jac_g);
      ASSERT_EXCEPTION(retval, ERROR_IN_TNLP_DERIVATIVE_TEST,
                       "In TNLP derivative test: Jacobian values could not be evaluated at reference point.");
   }

   // Space for the perturbed point
   Number* xpert = new Number[nx];
   IpBlasCopy(nx, xref, 1, xpert, 1);

   // Space for constraints at perturbed point
   Number* gpert = NULL;
   if( ng > 0 )
   {
      gpert = new Number[ng];
   }

   Index index_correction = 0;
   if( index_style == TNLP::FORTRAN_STYLE )
   {
      index_correction = 1;
   }

   if( deriv_test == FIRST_ORDER_TEST || deriv_test == SECOND_ORDER_TEST )
   {
      jnlst_->Printf(J_SUMMARY, J_NLP, "Starting derivative checker for first derivatives.\n\n");

      // Now go through all variables and check the partial derivatives
      const Index ivar_first = Max(Index(0), deriv_test_start_index);
      for( Index ivar = ivar_first; ivar < nx; ivar++ )
      {
         Number this_perturbation = derivative_test_perturbation_ * Max(Number(1.), std::abs(xref[ivar]));
         xpert[ivar] = xref[ivar] + this_perturbation;

         Number fpert;
         new_x = true;
         retval = tnlp_->eval_f(nx, xpert, new_x, fpert);
         ASSERT_EXCEPTION(retval, ERROR_IN_TNLP_DERIVATIVE_TEST,
                          "In TNLP derivative test: f could not be evaluated at perturbed point.");
         new_x = false;

         Number deriv_approx = (fpert - fref) / this_perturbation;
         Number deriv_exact = grad_f[ivar];
         Number rel_error = std::abs(deriv_approx - deriv_exact) / Max(std::abs(deriv_approx), derivative_test_tol_);
         char cflag = ' ';
         if( rel_error >= derivative_test_tol_ )
         {
            cflag = '*';
            nerrors++;
         }
         if( cflag != ' ' || derivative_test_print_all_ )
         {
            jnlst_->Printf(J_WARNING, J_NLP, "%c grad_f[      %5" IPOPT_INDEX_FORMAT "] = %23.16e    ~ %23.16e  [%10.3e]\n", cflag,
                           ivar + index_correction, deriv_exact, deriv_approx, rel_error);
         }

         if( ng > 0 )
         {
            retval = tnlp_->eval_g(nx, xpert, new_x, ng, gpert);
            ASSERT_EXCEPTION(retval, ERROR_IN_TNLP_DERIVATIVE_TEST,
                             "In TNLP derivative test: g could not be evaluated at reference point.");
            for( Index icon = 0; icon < ng; icon++ )
            {
               deriv_approx = (gpert[icon] - gref[icon]) / this_perturbation;
               deriv_exact = 0.;
               bool found = false;
               for( Index i = 0; i < nz_jac_g; i++ )
               {
                  if( g_iRow[i] == icon && g_jCol[i] == ivar )
                  {
                     found = true;
                     deriv_exact += jac_g[i];
                  }
               }
               rel_error = std::abs(deriv_approx - deriv_exact) / Max(std::abs(deriv_approx), derivative_test_tol_);
               cflag = ' ';
               if( rel_error >= derivative_test_tol_ )
               {
                  cflag = '*';
                  nerrors++;
               }
               char sflag = ' ';
               if( found )
               {
                  sflag = 'v';
               }
               if( cflag != ' ' || derivative_test_print_all_ )
               {
                  jnlst_->Printf(J_WARNING, J_NLP, "%c jac_g [%5" IPOPT_INDEX_FORMAT ",%5" IPOPT_INDEX_FORMAT "] = %23.16e %c  ~ %23.16e  [%10.3e]\n", cflag,
                                 icon + index_correction, ivar + index_correction, deriv_exact, sflag, deriv_approx, rel_error);
               }
            }
         }

         xpert[ivar] = xref[ivar];
      }

   }
   const Number zero = 0.;
   if( deriv_test == SECOND_ORDER_TEST || deriv_test == ONLY_SECOND_ORDER_TEST )
   {
      jnlst_->Printf(J_SUMMARY, J_NLP, "Starting derivative checker for second derivatives with obj_factor or lambda[i] set to 1.5.\n\n");

      // Get sparsity structure of Hessian
      Index* h_iRow = new Index[nz_hess_lag];
      Index* h_jCol = new Index[nz_hess_lag];
      retval = tnlp_->eval_h(nx, NULL, false, 0., ng, NULL, false, nz_hess_lag, h_iRow, h_jCol, NULL);
      ASSERT_EXCEPTION(retval, ERROR_IN_TNLP_DERIVATIVE_TEST,
                       "In TNLP derivative test: Hessian structure could not be evaluated.");

      if( index_style == TNLP::FORTRAN_STYLE )
      {
         for( Index i = 0; i < nz_hess_lag; i++ )
         {
            DBG_ASSERT(h_iRow[i] > 0);
            DBG_ASSERT(h_iRow[i] <= nx);
            DBG_ASSERT(h_jCol[i] > 0);
            DBG_ASSERT(h_jCol[i] <= nx);
            h_iRow[i] -= 1;
            h_jCol[i] -= 1;
         }
      }
#if IPOPT_CHECKLEVEL > 0
      else
      {
         for( Index i = 0; i < nz_hess_lag; i++ )
         {
            DBG_ASSERT(h_iRow[i] >= 0);
            DBG_ASSERT(h_iRow[i] < nx);
            DBG_ASSERT(h_jCol[i] >= 0);
            DBG_ASSERT(h_jCol[i] < nx);
         }
      }
#endif
      Number* h_values = new Number[nz_hess_lag];

      Number* lambda = NULL;
      if( ng > 0 )
      {
         lambda = new Number[ng];
         IpBlasCopy(ng, &zero, 0, lambda, 1);
      }
      Number* gradref = new Number[nx]; // gradient of objective or constraint at reference point
      Number* gradpert = new Number[nx]; // gradient of objective or constraint at perturbed point
      Number* jacpert = new Number[nz_jac_g];

      // Check all Hessians
      const Index icon_first = Max(Index(-1), deriv_test_start_index);
      for( Index icon = icon_first; icon < ng; icon++ )
      {
         Number objfact = 0.;
         if( icon == -1 )
         {
            objfact = 1.5;
            IpBlasCopy(nx, grad_f, 1, gradref, 1);
         }
         else
         {
            DBG_ASSERT(lambda != NULL);
            lambda[icon] = 1.5;
            IpBlasCopy(nx, &zero, 0, gradref, 1);
            for( Index i = 0; i < nz_jac_g; i++ )
            {
               if( g_iRow[i] == icon )
               {
                  gradref[g_jCol[i]] += jac_g[i];
               }
            }
         }
         // Hessian at reference point
         new_x = true;
         bool new_y = true;
         retval = tnlp_->eval_h(nx, xref, new_x, objfact, ng, lambda, new_y, nz_hess_lag, NULL, NULL, h_values);
         ASSERT_EXCEPTION(retval, ERROR_IN_TNLP_DERIVATIVE_TEST,
                          "In TNLP derivative test: Hessian could not be evaluated at reference point.");

         for( Index ivar = 0; ivar < nx; ivar++ )
         {
            Number this_perturbation = derivative_test_perturbation_ * Max(Number(1.), std::abs(xref[ivar]));
            xpert[ivar] = xref[ivar] + this_perturbation;

            new_x = true;
            if( icon == -1 )
            {
               // we are looking at the objective function
               retval = tnlp_->eval_grad_f(nx, xpert, new_x, gradpert);
               ASSERT_EXCEPTION(retval, ERROR_IN_TNLP_DERIVATIVE_TEST,
                                "In TNLP derivative test: grad_f could not be evaluated at perturbed point.");
            }
            else
            {
               // this is the icon-th constraint
               retval = tnlp_->eval_jac_g(nx, xpert, new_x, ng, nz_jac_g, NULL, NULL, jacpert);
               ASSERT_EXCEPTION(retval, ERROR_IN_TNLP_DERIVATIVE_TEST,
                                "In TNLP derivative test: Jacobian values could not be evaluated at reference point.");
               // ok, now we need to filter the gradient of the icon-th constraint
               IpBlasCopy(nx, &zero, 0, gradpert, 1);
               IpBlasCopy(nx, &zero, 0, gradref, 1);
               for( Index i = 0; i < nz_jac_g; i++ )
               {
                  if( g_iRow[i] == icon )
                  {
                     gradpert[g_jCol[i]] += jacpert[i];
                     gradref[g_jCol[i]] += jac_g[i];
                  }
               }
            }

            for( Index ivar2 = 0; ivar2 < nx; ivar2++ )
            {
               Number deriv_approx = 1.5 * (gradpert[ivar2] - gradref[ivar2]) / this_perturbation;
               Number deriv_exact = 0.;
               bool found = false;
               for( Index i = 0; i < nz_hess_lag; i++ )
               {
                  if( (h_iRow[i] == ivar && h_jCol[i] == ivar2) || (h_jCol[i] == ivar && h_iRow[i] == ivar2) )
                  {
                     deriv_exact += h_values[i];
                     found = true;
                  }
               }
               Number rel_error = std::abs(deriv_approx - deriv_exact) / Max(std::abs(deriv_approx), derivative_test_tol_);
               char cflag = ' ';
               if( rel_error >= derivative_test_tol_ )
               {
                  cflag = '*';
                  nerrors++;
               }
               char sflag = ' ';
               if( found )
               {
                  sflag = 'v';
               }
               if( cflag != ' ' || derivative_test_print_all_ )
               {
                  if( icon == -1 )
                  {
                     jnlst_->Printf(J_WARNING, J_NLP,
                                    "%c             obj_hess[%5" IPOPT_INDEX_FORMAT ",%5" IPOPT_INDEX_FORMAT "] = %23.16e %c  ~ %23.16e  [%10.3e]\n", cflag,
                                    ivar + index_correction, ivar2 + index_correction, deriv_exact, sflag, deriv_approx, rel_error);
                  }
                  else
                  {
                     jnlst_->Printf(J_WARNING, J_NLP,
                                    "%c %5" IPOPT_INDEX_FORMAT "-th constr_hess[%5" IPOPT_INDEX_FORMAT ",%5" IPOPT_INDEX_FORMAT "] = %23.16e %c  ~ %23.16e  [%10.3e]\n", cflag,
                                    icon + index_correction, ivar + index_correction, ivar2 + index_correction, deriv_exact, sflag,
                                    deriv_approx, rel_error);
                  }
               }

            }

            xpert[ivar] = xref[ivar];
         }

         if( icon >= 0 )
         {
            lambda[icon] = 0.;
         }
      }

      delete[] h_iRow;
      delete[] h_jCol;
      delete[] h_values;
      delete[] lambda;
      delete[] gradref;
      delete[] gradpert;
      delete[] jacpert;
   }

   delete[] xref;
   delete[] gref;
   delete[] grad_f;
   delete[] xpert;
   delete[] g_iRow;
   delete[] g_jCol;
   delete[] jac_g;
   delete[] gpert;

   if( nerrors == 0 )
   {
      jnlst_->Printf(J_SUMMARY, J_NLP, "\nNo errors detected by derivative checker.\n\n");
   }
   else
   {
      jnlst_->Printf(J_WARNING, J_NLP, "\nDerivative checker detected %" IPOPT_INDEX_FORMAT " error(s).\n\n", nerrors);
   }

   return retval;
}

bool TNLPAdapter::DetermineDependentConstraints(
   Index             n_x_var,
   const Index*      x_not_fixed_map,
   const Number*     x_l,
   const Number*     x_u,
   const Number*     g_l,
   const Number*     /*g_u*/,
   Index             n_c,
   const Index*      c_map,
   std::list<Index>& c_deps
)
{
   // First get a temporary expansion matrix for getting the equality
   // constraints
   SmartPtr<ExpansionMatrixSpace> P_c_g_space = new ExpansionMatrixSpace(n_full_g_, n_c, c_map);
   SmartPtr<ExpansionMatrix> P_c_g = P_c_g_space->MakeNewExpansionMatrix();

   // Get the structure of the big Jacobian of g and get the map for
   // the equality constraints entries
   Index* g_iRow = new Index[nz_full_jac_g_];
   Index* g_jCol = new Index[nz_full_jac_g_];
   if( !tnlp_->eval_jac_g(n_full_x_, NULL, false, n_full_g_, nz_full_jac_g_, g_iRow, g_jCol, NULL) )
   {
      delete[] g_iRow;
      delete[] g_jCol;
      return false;
   }
   if( index_style_ == TNLP::FORTRAN_STYLE )
   {
      for( Index i = 0; i < nz_full_jac_g_; i++ )
      {
         DBG_ASSERT(g_iRow[i] > 0);
         DBG_ASSERT(g_iRow[i] <= n_full_g_);
         DBG_ASSERT(g_jCol[i] > 0);
         DBG_ASSERT(g_jCol[i] <= n_full_x_);
         g_iRow[i] -= 1;
         g_jCol[i] -= 1;
      }
   }
#if IPOPT_CHECKLEVEL > 0
   else
   {
      for( Index i = 0; i < nz_full_jac_g_; i++ )
      {
         DBG_ASSERT(g_iRow[i] >= 0);
         DBG_ASSERT(g_iRow[i] < n_full_g_);
         DBG_ASSERT(g_jCol[i] >= 0);
         DBG_ASSERT(g_jCol[i] < n_full_x_);
      }
   }
#endif
   // TODO: Here we don't handle
   // fixed_variable_treatment_==MAKE_PARAMETER correctly (yet?)
   // Include space for the RHS
   Index* jac_c_map = new Index[nz_full_jac_g_];
   Index* jac_c_iRow = new Index[nz_full_jac_g_ + n_c];
   Index* jac_c_jCol = new Index[nz_full_jac_g_ + n_c];
   Index nz_jac_c = 0;
   const Index* c_row_pos = P_c_g->CompressedPosIndices();
   Index n_fixed = n_full_x_ - n_x_var;
   if( n_fixed > 0 )
   {
      // Get the reverse map for the fixed variables
      Index* c_col_pos = new Index[n_full_x_];
      for( Index i = 0; i < n_full_x_; i++ )
      {
         c_col_pos[i] = -1;
      }
      for( Index i = 0; i < n_x_var; i++ )
      {
         c_col_pos[x_not_fixed_map[i]] = i;
      }
      for( Index i = 0; i < nz_full_jac_g_; i++ )
      {
         const Index& c_row = c_row_pos[g_iRow[i]];
         const Index& c_col = c_col_pos[g_jCol[i]];
         if( c_col != -1 && c_row != -1 )
         {
            jac_c_map[nz_jac_c] = i;
            jac_c_iRow[nz_jac_c] = c_row + 1;
            jac_c_jCol[nz_jac_c] = c_col + 1;
            nz_jac_c++;
         }
      }
      delete[] c_col_pos;
   }
   else
   {
      for( Index i = 0; i < nz_full_jac_g_; i++ )
      {
         const Index& c_row = c_row_pos[g_iRow[i]];
         const Index& c_col = g_jCol[i];
         if( c_row != -1 )
         {
            jac_c_map[nz_jac_c] = i;
            jac_c_iRow[nz_jac_c] = c_row + 1;
            jac_c_jCol[nz_jac_c] = c_col + 1;
            nz_jac_c++;
         }
      }
   }
   delete[] g_iRow;
   delete[] g_jCol;

   // First we evaluate the equality constraint Jacobian at the
   // starting point with some random perturbation (projected into bounds)
   if( !tnlp_->get_starting_point(n_full_x_, true, full_x_, false, NULL,
                                  NULL, n_full_g_, false, NULL) )
   {
      delete[] jac_c_iRow;
      delete[] jac_c_jCol;
      delete[] jac_c_map;
      return false;
   }
   // Here we reset the random number generator
   IpResetRandom01();
   for( Index i = 0; i < n_full_x_; i++ )
   {
      const Number lower = Max(x_l[i], full_x_[i] - point_perturbation_radius_);
      const Number upper = Min(x_u[i], full_x_[i] + point_perturbation_radius_);
      const Number interval = upper - lower;
      const Number random_number = IpRandom01();
      full_x_[i] = lower + random_number * interval;
   }
   Number* g_vals = NULL;
   if( dependency_detection_with_rhs_ )
   {
      g_vals = new Number[n_full_g_];
      if( !tnlp_->eval_g(n_full_x_, full_x_, true, n_full_g_, g_vals) )
      {
         delete[] jac_c_iRow;
         delete[] jac_c_jCol;
         delete[] jac_c_map;
         delete[] g_vals;
         return false;
      }
   }
   if( !tnlp_->eval_jac_g(n_full_x_, full_x_, !dependency_detection_with_rhs_, n_full_g_, nz_full_jac_g_, NULL, NULL,
                          jac_g_) )
   {
      delete[] jac_c_iRow;
      delete[] jac_c_jCol;
      delete[] jac_c_map;
      delete[] g_vals;
      return false;
   }

   // Get the equality constraint Jacobian out
   Number* jac_c_vals = new Number[nz_jac_c + n_c];
   for( Index i = 0; i < nz_jac_c; i++ )
   {
      jac_c_vals[i] = jac_g_[jac_c_map[i]];
   }
   if( dependency_detection_with_rhs_ )
   {
      // Add the right hand side column
      const Index* c_row_pos2 = P_c_g->ExpandedPosIndices();
      for( Index i = 0; i < n_c; i++ )
      {
         jac_c_iRow[nz_jac_c + i] = i + 1;
         jac_c_jCol[nz_jac_c + i] = n_x_var + 1;
         jac_c_vals[nz_jac_c + i] = g_vals[c_row_pos2[i]] - g_l[c_row_pos2[i]];
      }
      n_x_var += 1;
      nz_jac_c += n_c;
   }

   ASSERT_EXCEPTION(IsValid(dependency_detector_), OPTION_INVALID,
                    "No dependency_detector_ object available in TNLPAdapter::DetermineDependentConstraints");

   bool retval = dependency_detector_->DetermineDependentRows(n_c, n_x_var, nz_jac_c, jac_c_vals, jac_c_iRow,
                 jac_c_jCol, c_deps);

   // For now, we just get rid of the dependency_detector_ object, in
   // order to save memory.  Maybe we need to add a clean method at
   // some point if we think that this is actually used more than
   // once...
   dependency_detector_ = NULL;

   delete[] jac_c_iRow;
   delete[] jac_c_jCol;
   delete[] jac_c_map;
   delete[] jac_c_vals;
   delete[] g_vals;

   return retval;
}

} // namespace Ipopt

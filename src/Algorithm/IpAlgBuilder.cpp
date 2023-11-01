// Copyright (C) 2004, 2011 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-09-29

#include "IpoptConfig.h"
#include "IpAlgBuilder.hpp"

#include "IpOrigIpoptNLP.hpp"
#include "IpIpoptData.hpp"
#include "IpIpoptCalculatedQuantities.hpp"
#include "IpCGPenaltyData.hpp"
#include "IpCGPenaltyCq.hpp"

#include "IpStdAugSystemSolver.hpp"
#include "IpAugRestoSystemSolver.hpp"
#include "IpPDFullSpaceSolver.hpp"
#include "IpPDPerturbationHandler.hpp"
#include "IpCGPerturbationHandler.hpp"
#include "IpOptErrorConvCheck.hpp"
#include "IpBacktrackingLineSearch.hpp"
#include "IpFilterLSAcceptor.hpp"
#include "IpCGPenaltyLSAcceptor.hpp"
#include "IpPenaltyLSAcceptor.hpp"
#include "IpPDSearchDirCalc.hpp"
#include "IpCGSearchDirCalc.hpp"
#include "IpMonotoneMuUpdate.hpp"
#include "IpAdaptiveMuUpdate.hpp"
#include "IpLoqoMuOracle.hpp"
#include "IpProbingMuOracle.hpp"
#include "IpQualityFunctionMuOracle.hpp"
#include "IpRestoMinC_1Nrm.hpp"
#include "IpLeastSquareMults.hpp"
#include "IpDefaultIterateInitializer.hpp"
#include "IpWarmStartIterateInitializer.hpp"
#include "IpOrigIterationOutput.hpp"
#include "IpLimMemQuasiNewtonUpdater.hpp"
#include "IpOrigIpoptNLP.hpp"
#include "IpLowRankAugSystemSolver.hpp"
#include "IpLowRankSSAugSystemSolver.hpp"
#include "IpRestoIterationOutput.hpp"
#include "IpRestoFilterConvCheck.hpp"
#include "IpRestoIterateInitializer.hpp"
#include "IpRestoPenaltyConvCheck.hpp"
#include "IpRestoRestoPhase.hpp"
#include "IpTSymLinearSolver.hpp"
#include "IpUserScaling.hpp"
#include "IpGradientScaling.hpp"
#include "IpEquilibrationScaling.hpp"
#include "IpExactHessianUpdater.hpp"
#include "IpSlackBasedTSymScalingMethod.hpp"

#include "IpLinearSolvers.h"
#include "IpMa27TSolverInterface.hpp"
#include "IpMa57TSolverInterface.hpp"
#include "IpMa77SolverInterface.hpp"
#include "IpMa86SolverInterface.hpp"
#include "IpMa97SolverInterface.hpp"
#include "IpMc19TSymScalingMethod.hpp"
#include "IpPardisoSolverInterface.hpp"
#ifdef IPOPT_HAS_PARDISO_MKL
# include "IpPardisoMKLSolverInterface.hpp"
#endif
#ifdef IPOPT_HAS_SPRAL
# include "IpSpralSolverInterface.hpp"
#endif
#ifdef IPOPT_HAS_WSMP
# include "IpWsmpSolverInterface.hpp"
# include "IpIterativeWsmpSolverInterface.hpp"
#endif
#ifdef IPOPT_HAS_MUMPS
# include "IpMumpsSolverInterface.hpp"
#endif

namespace Ipopt
{
#if IPOPT_VERBOSITY > 0
static const Index dbg_verbosity = 0;
#endif

AlgorithmBuilder::AlgorithmBuilder(
   SmartPtr<AugSystemSolver> custom_solver, /*=NULL*/
   const std::string& custom_solver_name    /*=std::string()*/
)
   : custom_solver_(custom_solver),
     custom_solver_name_(custom_solver_name)
{ }

void AlgorithmBuilder::RegisterOptions(
   SmartPtr<RegisteredOptions> roptions
)
{
   IpoptLinearSolver availablesolvers = IpoptGetAvailableLinearSolvers(false);
   IpoptLinearSolver availablesolverslinked = IpoptGetAvailableLinearSolvers(true);

   std::vector<std::string> options;
   std::vector<std::string> descrs;
   options.reserve(10);
   descrs.reserve(10);

   if( availablesolvers & IPOPTLINEARSOLVER_MA27 )
   {
      options.push_back("ma27");
      if( availablesolverslinked & IPOPTLINEARSOLVER_MA27 )
      {
         descrs.push_back("use the Harwell routine MA27");
      }
      else
      {
         descrs.push_back("load the Harwell routine MA27 from library at runtime");
      }
   }

   if( availablesolvers & IPOPTLINEARSOLVER_MA57 )
   {
      options.push_back("ma57");
      if( availablesolverslinked & IPOPTLINEARSOLVER_MA57 )
      {
         descrs.push_back("use the Harwell routine MA57");
      }
      else
      {
         descrs.push_back("load the Harwell routine MA57 from library at runtime");
      }
   }

   if( availablesolvers & IPOPTLINEARSOLVER_MA77 )
   {
      options.push_back("ma77");
      if( availablesolverslinked & IPOPTLINEARSOLVER_MA77 )
      {
         descrs.push_back("use the Harwell routine HSL_MA77");
      }
      else
      {
         descrs.push_back("load the Harwell routine HSL_MA77 from library at runtime");
      }
   }

   if( availablesolvers & IPOPTLINEARSOLVER_MA86 )
   {
      options.push_back("ma86");
      if( availablesolverslinked & IPOPTLINEARSOLVER_MA86 )
      {
         descrs.push_back("use the Harwell routine HSL_MA86");
      }
      else
      {
         descrs.push_back("load the Harwell routine MA86 from library at runtime");
      }
   }

   if( availablesolvers & IPOPTLINEARSOLVER_MA97 )
   {
      options.push_back("ma97");
      if( availablesolverslinked & IPOPTLINEARSOLVER_MA97 )
      {
         descrs.push_back("use the Harwell routine HSL_MA97");
      }
      else
      {
         descrs.push_back("load the Harwell routine MA97 from library at runtime");
      }
   }

   if( availablesolvers & IPOPTLINEARSOLVER_PARDISO )
   {
      options.push_back("pardiso");
      if( availablesolverslinked & IPOPTLINEARSOLVER_PARDISO )
      {
         descrs.push_back("use the Pardiso package from pardiso-project.org");
      }
      else
      {
         descrs.push_back("load the Pardiso package from pardiso-project.org from user-provided library at runtime");
      }
   }

   if( availablesolvers & IPOPTLINEARSOLVER_PARDISOMKL )
   {
      options.push_back("pardisomkl");
      descrs.push_back("use the Pardiso package from Intel MKL");
   }

   if( availablesolvers & IPOPTLINEARSOLVER_SPRAL )
   {
      options.push_back("spral");
      descrs.push_back("use the Spral package");
   }

   if( availablesolvers & IPOPTLINEARSOLVER_WSMP )
   {
      options.push_back("wsmp");
      descrs.push_back("use the Wsmp package");
   }

   if( availablesolvers & IPOPTLINEARSOLVER_MUMPS )
   {
      options.push_back("mumps");
      descrs.push_back("use the Mumps package");
   }

   options.push_back("custom");
   descrs.push_back("use custom linear solver (expert use)");

   std::string defaultsolver;
   if( availablesolverslinked & IPOPTLINEARSOLVER_MA27 )
   {
      defaultsolver = "ma27";
   }
   else if( availablesolverslinked & IPOPTLINEARSOLVER_MA57 )
   {
      defaultsolver = "ma57";
   }
   else if( availablesolverslinked & IPOPTLINEARSOLVER_MA97 )
   {
      defaultsolver = "ma97";
   }
   else if( availablesolverslinked & IPOPTLINEARSOLVER_MA86 )
   {
      defaultsolver = "ma86";
   }
   else if( availablesolverslinked & IPOPTLINEARSOLVER_PARDISO )
   {
      defaultsolver = "pardiso";
   }
   else if( availablesolverslinked & IPOPTLINEARSOLVER_WSMP )
   {
      defaultsolver = "wsmp";
   }
   else if( availablesolverslinked & IPOPTLINEARSOLVER_MUMPS )
   {
      defaultsolver = "mumps";
   }
   else if( availablesolverslinked & IPOPTLINEARSOLVER_PARDISOMKL )
   {
      defaultsolver = "pardisomkl";
   }
   else if( availablesolverslinked & IPOPTLINEARSOLVER_SPRAL )
   {
      defaultsolver = "spral";
   }
   else if( availablesolverslinked & IPOPTLINEARSOLVER_MA77 )
   {
      defaultsolver = "ma77";
   }
   else if( availablesolvers & IPOPTLINEARSOLVER_MA27 )
   {
      defaultsolver = "ma27";
   }
   else
   {
      defaultsolver = "custom";
   }

   roptions->SetRegisteringCategory("Linear Solver");
   roptions->AddStringOption(
      "linear_solver",
      "Linear solver used for step computations.",
      defaultsolver,
      options,
      descrs,
      "Determines which linear algebra package is to be used for the solution of the augmented linear system (for obtaining the search directions).");

   options.clear();
   descrs.clear();

   std::string longdescr =
      "Determines the method used to compute symmetric scaling factors for the augmented system "
      "(see also the \"linear_scaling_on_demand\" option). "
      "This scaling is independent of the NLP problem scaling.";

   options.push_back("none");
   descrs.push_back("no scaling will be performed");
   defaultsolver = "none";

   if( availablesolvers & IPOPTLINEARSOLVER_MC19 )
   {
      options.push_back("mc19");

      if( availablesolverslinked & IPOPTLINEARSOLVER_MC19 )
      {
         descrs.push_back("use the Harwell routine MC19");
         defaultsolver = "mc19";
         longdescr += " The default is MC19 only if MA27, MA57, MA77, or MA86 are selected as linear solvers. Otherwise it is 'none'.";
      }
      else
      {
         descrs.push_back("load the Harwell routine MC19 from library at runtime");
      }
   }

   options.push_back("slack-based");
   descrs.push_back("use the slack values");

   roptions->AddStringOption(
      "linear_system_scaling", "Method for scaling the linear system.",
      defaultsolver,
      options, descrs,
      longdescr);

   // have hsllib option if some HSL solvers are not linked but can be loaded
   if( (availablesolverslinked ^ availablesolvers) & IPOPTLINEARSOLVER_ALLHSL )
      roptions->AddStringOption1(
         "hsllib", "Name of library containing HSL routines for load at runtime",
         "libhsl." IPOPT_SHAREDLIBEXT,
         "*", "Any acceptable filename (may contain path, too)");

   roptions->AddStringOption1(
      "pardisolib", "Name of library containing Pardiso routines (from pardiso-project.org) for load at runtime",
#ifdef PARDISO_LIB
      PARDISO_LIB,
#else
      "libpardiso." IPOPT_SHAREDLIBEXT,
#endif
      "*", "Any acceptable filename (may contain path, too)");

   roptions->SetRegisteringCategory("NLP Scaling");

   options.clear();
   descrs.clear();
   options.push_back("none");
   descrs.push_back("no problem scaling will be performed");
   options.push_back("user-scaling");
   descrs.push_back("scaling parameters will come from the user");
   options.push_back("gradient-based");
   descrs.push_back("scale the problem so the maximum gradient at the starting point is nlp_scaling_max_gradient");

   if( availablesolvers & IPOPTLINEARSOLVER_MC19 )
   {
      options.push_back("equilibration-based");
      descrs.push_back("scale the problem so that first derivatives are of order 1 at random points");
      if( availablesolverslinked & IPOPTLINEARSOLVER_MC19 )
      {
         descrs.back() += " (uses Harwell routine MC19)";
      }
      else
      {
         descrs.back() += " (load the Harwell routine MC19 from library at runtime)";
      }
   }
   roptions->AddStringOption(
      "nlp_scaling_method", "Select the technique used for scaling the NLP.",
      "gradient-based",
      options,
      descrs,
      "Selects the technique used for scaling the problem internally before it is solved. "
      "For user-scaling, the parameters come from the NLP."
#ifdef IPOPT_HAS_ASL
      " If you are using AMPL, they can be specified through suffixes (\"scaling_factor\")"
#endif
   );

   roptions->SetRegisteringCategory("Barrier Parameter Update");
   roptions->AddStringOption2(
      "mu_strategy",
      "Update strategy for barrier parameter.",
      "monotone",
      "monotone", "use the monotone (Fiacco-McCormick) strategy",
      "adaptive", "use the adaptive update strategy",
      "Determines which barrier parameter update strategy is to be used.");
   roptions->AddStringOption3(
      "mu_oracle",
      "Oracle for a new barrier parameter in the adaptive strategy.",
      "quality-function",
      "probing", "Mehrotra's probing heuristic",
      "loqo", "LOQO's centrality rule",
      "quality-function", "minimize a quality function",
      "Determines how a new barrier parameter is computed in each \"free-mode\" iteration of the adaptive barrier parameter strategy. "
      "(Only considered if \"adaptive\" is selected for option \"mu_strategy\").");
   roptions->AddStringOption4(
      "fixed_mu_oracle",
      "Oracle for the barrier parameter when switching to fixed mode.",
      "average_compl",
      "probing", "Mehrotra's probing heuristic",
      "loqo", "LOQO's centrality rule",
      "quality-function", "minimize a quality function",
      "average_compl", "base on current average complementarity",
      "Determines how the first value of the barrier parameter should be computed when switching to the \"monotone mode\" in the adaptive strategy. "
      "(Only considered if \"adaptive\" is selected for option \"mu_strategy\".)");

   roptions->SetRegisteringCategory("Hessian Approximation");
   roptions->AddStringOption2(
      "limited_memory_aug_solver",
      "Strategy for solving the augmented system for low-rank Hessian.",
      "sherman-morrison",
      "sherman-morrison", "use Sherman-Morrison formula",
      "extended", "use an extended augmented system",
      "",
      true);

   roptions->SetRegisteringCategory("Line Search");
   roptions->AddStringOption3(
      "line_search_method",
      "Globalization method used in backtracking line search",
      "filter",
      "filter", "Filter method",
      "cg-penalty", "Chen-Goldfarb penalty function",
      "penalty", "Standard penalty function",
      "Only the \"filter\" choice is officially supported. "
      "But sometimes, good results might be obtained with the other choices.",
      true);
   roptions->SetRegisteringCategory("Undocumented");
   roptions->AddBoolOption(
      "wsmp_iterative",
      "Switches to use iterative instead of direct solver in WSMP.",
      false,
      "EXPERIMENTAL!",
      true);
}

SmartPtr<SymLinearSolver> AlgorithmBuilder::GetSymLinearSolver(
   const Journalist&     jnlst,
   const OptionsList&    options,
   const std::string&    prefix
)
{
   if( IsNull(SymSolver_) )
   {
      SymSolver_ = SymLinearSolverFactory(jnlst, options, prefix);
   }
   DBG_ASSERT(IsValid(SymSolver_));
   return SymSolver_;
}

SmartPtr<SymLinearSolver> AlgorithmBuilder::SymLinearSolverFactory(
   const Journalist&     /*jnlst*/,
   const OptionsList&    options,
   const std::string&    prefix
)
{
   SmartPtr<SparseSymLinearSolverInterface> SolverInterface;
   options.GetStringValue("linear_solver", linear_solver, prefix);

   if( false ) ;

#ifndef IPOPT_INT64
   else if( linear_solver == "ma27" )
   {
      SolverInterface = new Ma27TSolverInterface(GetHSLLoader(options, prefix));
   }

   else if( linear_solver == "ma57" )
   {
      SolverInterface = new Ma57TSolverInterface(GetHSLLoader(options, prefix));
   }

   else if( linear_solver == "ma77" )
   {
      SolverInterface = new Ma77SolverInterface(GetHSLLoader(options, prefix));
   }

   else if( linear_solver == "ma86" )
   {
      SolverInterface = new Ma86SolverInterface(GetHSLLoader(options, prefix));
   }

   else if( linear_solver == "ma97" )
   {
      SolverInterface = new Ma97SolverInterface(GetHSLLoader(options, prefix));
   }

   else if( linear_solver == "pardiso" )
   {
      SolverInterface = new PardisoSolverInterface(GetPardisoLoader(options, prefix));
   }
#endif

#ifdef IPOPT_HAS_PARDISO_MKL
   else if( linear_solver == "pardisomkl" )
   {
      SolverInterface = new PardisoMKLSolverInterface();
   }
#endif

#ifdef IPOPT_HAS_SPRAL
   else if( linear_solver == "spral" )
   {
      SolverInterface = new SpralSolverInterface();
   }
#endif

#ifdef IPOPT_HAS_WSMP
   else if( linear_solver == "wsmp" )
   {
      bool wsmp_iterative;
      options.GetBoolValue("wsmp_iterative", wsmp_iterative, prefix);
      if( wsmp_iterative )
      {
         SolverInterface = new IterativeWsmpSolverInterface();
      }
      else
      {
#ifdef PARDISO_MATCHING_PREPROCESS
         SolverInterface = new WsmpSolverInterface(GetPardisoLoader(options, prefix));
#else
         SolverInterface = new WsmpSolverInterface();
#endif
      }
      int V, R, M;
      WsmpSolverInterface::GetVersion(V, R, M);
      char buffer[100];
      Snprintf(buffer, 100, "WSMP %d.%d.%d\n", V, R, M);
      linear_solver = buffer;
   }
#endif

#ifdef IPOPT_HAS_MUMPS
   else if( linear_solver == "mumps" )
   {
      SolverInterface = new MumpsSolverInterface();
      linear_solver = MumpsSolverInterface::GetName();
   }
#endif

   else if( linear_solver == "custom" )
   {
      SolverInterface = NULL;
   }

   else
   {
      // this should have been checked earlier
      THROW_EXCEPTION(OPTION_INVALID, "Invalid value selected for option linear_solver");
   }

   SmartPtr<TSymScalingMethod> ScalingMethod;
   std::string linear_system_scaling;
   if( !options.GetStringValue("linear_system_scaling", linear_system_scaling, prefix) )
   {
      // By default, don't use mc19 for non-HSL solvers, or HSL_MA97
      if( linear_solver != "ma27" && linear_solver != "ma57" && linear_solver != "ma77" && linear_solver != "ma86" )
      {
         linear_system_scaling = "none";
      }
   }

   if( linear_system_scaling == "slack-based" )
   {
      ScalingMethod = new SlackBasedTSymScalingMethod();
   }
#ifndef IPOPT_INT64
   else if( linear_system_scaling == "mc19" )
   {
      ScalingMethod = new Mc19TSymScalingMethod(GetHSLLoader(options, prefix));
   }
#endif

   SmartPtr<SymLinearSolver> ScaledSolver = new TSymLinearSolver(SolverInterface, ScalingMethod);
   return ScaledSolver;
}

SmartPtr<AugSystemSolver> AlgorithmBuilder::GetAugSystemSolver(
   const Journalist&     jnlst,
   const OptionsList&    options,
   const std::string&    prefix
)
{
   if( IsNull(AugSolver_) )
   {
      AugSolver_ = AugSystemSolverFactory(jnlst, options, prefix);
   }
   DBG_ASSERT(IsValid(AugSolver_));
   return AugSolver_;
}

SmartPtr<AugSystemSolver> AlgorithmBuilder::AugSystemSolverFactory(
   const Journalist&     jnlst,
   const OptionsList&    options,
   const std::string&    prefix
)
{
   SmartPtr<AugSystemSolver> AugSolver;
   options.GetStringValue("linear_solver", linear_solver, prefix);
   if( linear_solver == "custom" )
   {
      ASSERT_EXCEPTION(IsValid(custom_solver_), OPTION_INVALID, "Selected linear solver CUSTOM not available.");
      AugSolver = custom_solver_;
      if( !custom_solver_name_.empty() )
      {
         linear_solver = custom_solver_name_;
      }
   }
   else
   {
      AugSolver = new StdAugSystemSolver(*GetSymLinearSolver(jnlst, options, prefix));
   }

   Index enum_int;
   options.GetEnumValue("hessian_approximation", enum_int, prefix);
   HessianApproximationType hessian_approximation = HessianApproximationType(enum_int);
   if( hessian_approximation == LIMITED_MEMORY )
   {
      std::string lm_aug_solver;
      options.GetStringValue("limited_memory_aug_solver", lm_aug_solver, prefix);
      if( lm_aug_solver == "sherman-morrison" )
      {
         AugSolver = new LowRankAugSystemSolver(*AugSolver);
      }
      else if( lm_aug_solver == "extended" )
      {
         Index lm_history;
         options.GetIntegerValue("limited_memory_max_history", lm_history, prefix);
         Index max_rank;
         std::string lm_type;
         options.GetStringValue("limited_memory_update_type", lm_type, prefix);
         if( lm_type == "bfgs" )
         {
            max_rank = 2 * lm_history;
         }
         else if( lm_type == "sr1" )
         {
            max_rank = lm_history;
         }
         else
         {
            THROW_EXCEPTION(OPTION_INVALID, "Unknown value for option \"limited_memory_update_type\".");
         }
         AugSolver = new LowRankSSAugSystemSolver(*AugSolver, max_rank);
      }
      else
      {
         THROW_EXCEPTION(OPTION_INVALID, "Unknown value for option \"limited_memory_aug_solver\".");
      }
   }
   return AugSolver;
}

SmartPtr<PDSystemSolver> AlgorithmBuilder::GetPDSystemSolver(
   const Journalist&     jnlst,
   const OptionsList&    options,
   const std::string&    prefix
)
{
   if( IsNull(PDSolver_) )
   {
      PDSolver_ = PDSystemSolverFactory(jnlst, options, prefix);
   }
   DBG_ASSERT(IsValid(PDSolver_));
   return PDSolver_;
}

SmartPtr<PDSystemSolver> AlgorithmBuilder::PDSystemSolverFactory(
   const Journalist&     jnlst,
   const OptionsList&    options,
   const std::string&    prefix
)
{
   SmartPtr<PDPerturbationHandler> pertHandler;
   std::string lsmethod;
   options.GetStringValue("line_search_method", lsmethod, prefix);
   if( lsmethod == "cg-penalty" )
   {
      pertHandler = new CGPerturbationHandler();
   }
   else
   {
      pertHandler = new PDPerturbationHandler();
   }

   SmartPtr<PDSystemSolver> PDSolver = new PDFullSpaceSolver(*GetAugSystemSolver(jnlst, options, prefix), *pertHandler);
   return PDSolver;
}

void AlgorithmBuilder::BuildIpoptObjects(
   const Journalist&     jnlst,
   const OptionsList&    options,
   const std::string&    prefix,
   const SmartPtr<NLP>&  nlp,
   SmartPtr<IpoptNLP>&   ip_nlp,
   SmartPtr<IpoptData>&  ip_data,
   SmartPtr<IpoptCalculatedQuantities>& ip_cq
)
{
   DBG_ASSERT(prefix == "");

   SmartPtr<NLPScalingObject> nlp_scaling;
   std::string nlp_scaling_method;
   options.GetStringValue("nlp_scaling_method", nlp_scaling_method, "");
   if( nlp_scaling_method == "user-scaling" )
   {
      nlp_scaling = new UserScaling(ConstPtr(nlp));
   }
   else if( nlp_scaling_method == "gradient-based" )
   {
      nlp_scaling = new GradientScaling(nlp);
   }
   else if( nlp_scaling_method == "equilibration-based" )
   {
      nlp_scaling = new EquilibrationScaling(nlp, GetHSLLoader(options, prefix));
   }
   else
   {
      nlp_scaling = new NoNLPScalingObject();
   }

   // Create the IpoptData.  Check if there is additional data that
   // is needed
   std::string lsmethod;
   SmartPtr<IpoptAdditionalData> add_data;
   options.GetStringValue("line_search_method", lsmethod, prefix);
   if( lsmethod == "cg-penalty" )
   {
      add_data = new CGPenaltyData();
   }
   ip_data = new IpoptData(add_data);

   ip_nlp = new OrigIpoptNLP(&jnlst, GetRawPtr(nlp), nlp_scaling, ip_data->TimingStats());

   // Create the IpoptCalculators.  Check if there are additional
   // calculated quantities that are needed
   ip_cq = new IpoptCalculatedQuantities(ip_nlp, ip_data);
   if( lsmethod == "cg-penalty" )
   {
      SmartPtr<IpoptAdditionalCq> add_cq = new CGPenaltyCq(GetRawPtr(ip_nlp), GetRawPtr(ip_data), GetRawPtr(ip_cq));
      ip_cq->SetAddCq(add_cq);
   }
}

SmartPtr<IpoptAlgorithm> AlgorithmBuilder::BuildBasicAlgorithm(
   const Journalist&     jnlst,
   const OptionsList&    options,
   const std::string&    prefix
)
{

   /* These three objects don't have any dependencies on other
    * other major components of the algorithm */
   IterOutput_ = BuildIterationOutput(jnlst, options, prefix);
   HessUpdater_ = BuildHessianUpdater(jnlst, options, prefix);
   ConvCheck_ = BuildConvergenceCheck(jnlst, options, prefix);

   /* All solver type factory methods are likely
    * first invoked here */
   SearchDirCalc_ = BuildSearchDirectionCalculator(jnlst, options, prefix);
   EqMultCalculator_ = BuildEqMultiplierCalculator(jnlst, options, prefix);

   /* Requires:
    *   -> EqMultCalculator_
    */
   IterInitializer_ = BuildIterateInitializer(jnlst, options, prefix);

   /* Requires:
    *   -> ConvCheck_
    *   -> EqMultCalculator_
    */
   LineSearch_ = BuildLineSearch(jnlst, options, prefix);

   /* Requires:
    *   -> LineSearch_
    *      -> ConvCheck_
    *      -> EqMultCalculator_
    */
   MuUpdate_ = BuildMuUpdate(jnlst, options, prefix);

   SmartPtr<IpoptAlgorithm> alg = new IpoptAlgorithm(SearchDirCalc_, LineSearch_, MuUpdate_, ConvCheck_,
         IterInitializer_, IterOutput_, HessUpdater_, EqMultCalculator_, linear_solver);

   return alg;
}

SmartPtr<IterationOutput> AlgorithmBuilder::BuildIterationOutput(
   const Journalist&     /*jnlst*/,
   const OptionsList&    /*options*/,
   const std::string&    /*prefix*/
)
{
   // Create the object for the iteration output
   SmartPtr<IterationOutput> IterOutput = new OrigIterationOutput();
   return IterOutput;
}

SmartPtr<HessianUpdater> AlgorithmBuilder::BuildHessianUpdater(
   const Journalist&     /*jnlst*/,
   const OptionsList&    options,
   const std::string&    prefix
)
{
   // Get the Hessian updater for the main algorithm
   SmartPtr<HessianUpdater> HessUpdater;
   Index enum_int;
   options.GetEnumValue("hessian_approximation", enum_int, prefix);
   HessianApproximationType hessian_approximation = HessianApproximationType(enum_int);
   switch( hessian_approximation )
   {
      case EXACT:
         HessUpdater = new ExactHessianUpdater();
         break;
      case LIMITED_MEMORY:
         // ToDo This needs to be replaced!
         HessUpdater = new LimMemQuasiNewtonUpdater(false);
         break;
   }
   return HessUpdater;
}

SmartPtr<ConvergenceCheck> AlgorithmBuilder::BuildConvergenceCheck(
   const Journalist&     /*jnlst*/,
   const OptionsList&    /*options*/,
   const std::string&    /*prefix*/
)
{
   SmartPtr<ConvergenceCheck> ConvCheck = new OptimalityErrorConvergenceCheck();
   return ConvCheck;
}

SmartPtr<SearchDirectionCalculator> AlgorithmBuilder::BuildSearchDirectionCalculator(
   const Journalist&     jnlst,
   const OptionsList&    options,
   const std::string&    prefix
)
{
   std::string lsmethod;
   options.GetStringValue("line_search_method", lsmethod, prefix);
   SmartPtr<SearchDirectionCalculator> SearchDirCalc;
   if( lsmethod == "cg-penalty" )
   {
      SearchDirCalc = new CGSearchDirCalculator(GetRawPtr(GetPDSystemSolver(jnlst, options, prefix)));
   }
   else
   {
      SearchDirCalc = new PDSearchDirCalculator(GetRawPtr(GetPDSystemSolver(jnlst, options, prefix)));
   }
   return SearchDirCalc;
}

SmartPtr<EqMultiplierCalculator> AlgorithmBuilder::BuildEqMultiplierCalculator(
   const Journalist&     jnlst,
   const OptionsList&    options,
   const std::string&    prefix
)
{
   SmartPtr<EqMultiplierCalculator> EqMultCalculator = new LeastSquareMultipliers(
      *GetAugSystemSolver(jnlst, options, prefix));
   return EqMultCalculator;
}

SmartPtr<IterateInitializer> AlgorithmBuilder::BuildIterateInitializer(
   const Journalist&     jnlst,
   const OptionsList&    options,
   const std::string&    prefix
)
{
   DBG_ASSERT(IsValid(EqMultCalculator_));

   // Create the object for initializing the iterates Initialization
   // object.  We include both the warm start and the default
   // initializer, so that the warm start options can be activated
   // without having to rebuild the algorithm
   SmartPtr<IterateInitializer> WarmStartInitializer = new WarmStartIterateInitializer();

   SmartPtr<IterateInitializer> IterInitializer = new DefaultIterateInitializer(EqMultCalculator_, WarmStartInitializer,
         GetAugSystemSolver(jnlst, options, prefix));
   return IterInitializer;
}

SmartPtr<LineSearch> AlgorithmBuilder::BuildLineSearch(
   const Journalist&     jnlst,
   const OptionsList&    options,
   const std::string&    prefix
)
{
   DBG_ASSERT(IsValid(ConvCheck_));
   DBG_ASSERT(IsValid(EqMultCalculator_));

   Index enum_int;
   options.GetEnumValue("hessian_approximation", enum_int, prefix);
   HessianApproximationType hessian_approximation = HessianApproximationType(enum_int);

   SmartPtr<RestorationPhase> resto_phase;
   SmartPtr<RestoConvergenceCheck> resto_convCheck;

   // We only need a restoration phase object if we use the filter
   // line search
   std::string lsmethod;
   options.GetStringValue("line_search_method", lsmethod, prefix);
   if( lsmethod == "filter" || lsmethod == "penalty" )
   {
      // Solver for the restoration phase
      SmartPtr<AugSystemSolver> resto_AugSolver = new AugRestoSystemSolver(*GetAugSystemSolver(jnlst, options, prefix));
      SmartPtr<PDPerturbationHandler> resto_pertHandler = new PDPerturbationHandler();
      SmartPtr<PDSystemSolver> resto_PDSolver = new PDFullSpaceSolver(*resto_AugSolver, *resto_pertHandler);

      // Convergence check in the restoration phase
      if( lsmethod == "filter" )
      {
         resto_convCheck = new RestoFilterConvergenceCheck();
      }
      else if( lsmethod == "penalty" )
      {
         resto_convCheck = new RestoPenaltyConvergenceCheck();
      }

      // Line search method for the restoration phase
      SmartPtr<RestoRestorationPhase> resto_resto = new RestoRestorationPhase();

      SmartPtr<BacktrackingLSAcceptor> resto_LSacceptor;
      std::string resto_lsacceptor;
      options.GetStringValue("line_search_method", resto_lsacceptor, "resto." + prefix);
      if( resto_lsacceptor == "filter" )
      {
         resto_LSacceptor = new FilterLSAcceptor(GetRawPtr(resto_PDSolver));
      }
      else if( resto_lsacceptor == "cg-penalty" )
      {
         resto_LSacceptor = new CGPenaltyLSAcceptor(GetRawPtr(resto_PDSolver));
      }
      else if( resto_lsacceptor == "penalty" )
      {
         resto_LSacceptor = new PenaltyLSAcceptor(GetRawPtr(resto_PDSolver));
      }
      SmartPtr<LineSearch> resto_LineSearch = new BacktrackingLineSearch(resto_LSacceptor, GetRawPtr(resto_resto),
            GetRawPtr(resto_convCheck));

      // Create the mu update that will be used by the restoration phase
      // algorithm
      SmartPtr<MuUpdate> resto_MuUpdate;
      std::string resto_smuupdate;
      if( !options.GetStringValue("mu_strategy", resto_smuupdate, "resto." + prefix) )
      {
         // Change default for quasi-Newton option (then we use adaptive)
         if( hessian_approximation == LIMITED_MEMORY )
         {
            resto_smuupdate = "adaptive";
         }
      }

      std::string resto_smuoracle;
      std::string resto_sfixmuoracle;
      if( resto_smuupdate == "adaptive" )
      {
         options.GetStringValue("mu_oracle", resto_smuoracle, "resto." + prefix);
         options.GetStringValue("fixed_mu_oracle", resto_sfixmuoracle, "resto." + prefix);
      }

      if( resto_smuupdate == "monotone" )
      {
         resto_MuUpdate = new MonotoneMuUpdate(GetRawPtr(resto_LineSearch));
      }
      else if( resto_smuupdate == "adaptive" )
      {
         SmartPtr<MuOracle> resto_MuOracle;
         if( resto_smuoracle == "loqo" )
         {
            resto_MuOracle = new LoqoMuOracle();
         }
         else if( resto_smuoracle == "probing" )
         {
            resto_MuOracle = new ProbingMuOracle(resto_PDSolver);
         }
         else if( resto_smuoracle == "quality-function" )
         {
            resto_MuOracle = new QualityFunctionMuOracle(resto_PDSolver);
         }
         SmartPtr<MuOracle> resto_FixMuOracle;
         if( resto_sfixmuoracle == "loqo" )
         {
            resto_FixMuOracle = new LoqoMuOracle();
         }
         else if( resto_sfixmuoracle == "probing" )
         {
            resto_FixMuOracle = new ProbingMuOracle(resto_PDSolver);
         }
         else if( resto_sfixmuoracle == "quality-function" )
         {
            resto_FixMuOracle = new QualityFunctionMuOracle(resto_PDSolver);
         }
         else
         {
            resto_FixMuOracle = NULL;
         }
         resto_MuUpdate = new AdaptiveMuUpdate(GetRawPtr(resto_LineSearch), resto_MuOracle, resto_FixMuOracle);
      }

      // Initialization of the iterates for the restoration phase
      SmartPtr<EqMultiplierCalculator> resto_EqMultCalculator = new LeastSquareMultipliers(*resto_AugSolver);
      SmartPtr<IterateInitializer> resto_IterInitializer = new RestoIterateInitializer(resto_EqMultCalculator);

      // Create the object for the iteration output during restoration
      SmartPtr<OrigIterationOutput> resto_OrigIterOutput = NULL;
      //   new OrigIterationOutput();
      SmartPtr<IterationOutput> resto_IterOutput = new RestoIterationOutput(resto_OrigIterOutput);

      // Get the Hessian updater for the restoration phase
      SmartPtr<HessianUpdater> resto_HessUpdater;
      switch( hessian_approximation )
      {
         case EXACT:
            resto_HessUpdater = new ExactHessianUpdater();
            break;
         case LIMITED_MEMORY:
            // ToDo This needs to be replaced!
            resto_HessUpdater = new LimMemQuasiNewtonUpdater(true);
            break;
      }

      // Put together the overall restoration phase IP algorithm
      SmartPtr<SearchDirectionCalculator> resto_SearchDirCalc;
      if( resto_lsacceptor == "cg-penalty" )
      {
         resto_SearchDirCalc = new CGSearchDirCalculator(GetRawPtr(resto_PDSolver));
      }
      else
      {
         resto_SearchDirCalc = new PDSearchDirCalculator(GetRawPtr(resto_PDSolver));
      }

      SmartPtr<IpoptAlgorithm> resto_alg = new IpoptAlgorithm(resto_SearchDirCalc, GetRawPtr(resto_LineSearch),
            GetRawPtr(resto_MuUpdate), GetRawPtr(resto_convCheck), resto_IterInitializer, resto_IterOutput,
            resto_HessUpdater, resto_EqMultCalculator, linear_solver);

      // Set the restoration phase
      resto_phase = new MinC_1NrmRestorationPhase(*resto_alg, EqMultCalculator_);
   }

   // Create the line search to be used by the main algorithm
   SmartPtr<BacktrackingLSAcceptor> LSacceptor;
   if( lsmethod == "filter" )
   {
      LSacceptor = new FilterLSAcceptor(GetRawPtr(GetPDSystemSolver(jnlst, options, prefix)));
   }
   else if( lsmethod == "cg-penalty" )
   {
      LSacceptor = new CGPenaltyLSAcceptor(GetRawPtr(GetPDSystemSolver(jnlst, options, prefix)));
   }
   else if( lsmethod == "penalty" )
   {
      LSacceptor = new PenaltyLSAcceptor(GetRawPtr(GetPDSystemSolver(jnlst, options, prefix)));
   }
   SmartPtr<LineSearch> LineSearch = new BacktrackingLineSearch(LSacceptor, GetRawPtr(resto_phase), ConvCheck_);

   // The following cross reference is not good: We have to store a
   // pointer to the LineSearch_ object in resto_convCheck as a
   // non-SmartPtr to make sure that things are properly deleted when
   // the IpoptAlgorithm returned by the Builder is destructed.
   if( IsValid(resto_convCheck) )
   {
      resto_convCheck->SetOrigLSAcceptor(*LSacceptor);
   }

   return LineSearch;
}

SmartPtr<MuUpdate> AlgorithmBuilder::BuildMuUpdate(
   const Journalist&     jnlst,
   const OptionsList&    options,
   const std::string&    prefix
)
{
   DBG_ASSERT(IsValid(LineSearch_));

   bool mehrotra_algorithm;
   options.GetBoolValue("mehrotra_algorithm", mehrotra_algorithm, prefix);

   // Create the mu update that will be used by the main algorithm
   SmartPtr<MuUpdate> MuUpdate;
   std::string smuupdate;
   if( !options.GetStringValue("mu_strategy", smuupdate, prefix) )
   {
      // Change default for quasi-Newton option (then we use adaptive)
      Index enum_int;
      if( options.GetEnumValue("hessian_approximation", enum_int, prefix) )
      {
         HessianApproximationType hessian_approximation = HessianApproximationType(enum_int);
         if( hessian_approximation == LIMITED_MEMORY )
         {
            smuupdate = "adaptive";
         }
      }
      if( mehrotra_algorithm )
      {
         smuupdate = "adaptive";
      }
   }
   ASSERT_EXCEPTION(!mehrotra_algorithm || smuupdate == "adaptive", OPTION_INVALID,
                    "If mehrotra_algorithm=yes, mu_strategy must be \"adaptive\".");
   std::string smuoracle;
   std::string sfixmuoracle;
   if( smuupdate == "adaptive" )
   {
      if( !options.GetStringValue("mu_oracle", smuoracle, prefix) )
      {
         if( mehrotra_algorithm )
         {
            smuoracle = "probing";
         }
      }
      options.GetStringValue("fixed_mu_oracle", sfixmuoracle, prefix);
      ASSERT_EXCEPTION(!mehrotra_algorithm || smuoracle == "probing", OPTION_INVALID,
                       "If mehrotra_algorithm=yes, mu_oracle must be \"probing\".");
   }

   if( smuupdate == "monotone" )
   {
      MuUpdate = new MonotoneMuUpdate(GetRawPtr(LineSearch_));
   }
   else if( smuupdate == "adaptive" )
   {
      SmartPtr<MuOracle> muOracle;
      if( smuoracle == "loqo" )
      {
         muOracle = new LoqoMuOracle();
      }
      else if( smuoracle == "probing" )
      {
         muOracle = new ProbingMuOracle(GetPDSystemSolver(jnlst, options, prefix));
      }
      else if( smuoracle == "quality-function" )
      {
         muOracle = new QualityFunctionMuOracle(GetPDSystemSolver(jnlst, options, prefix));
      }
      SmartPtr<MuOracle> FixMuOracle;
      if( sfixmuoracle == "loqo" )
      {
         FixMuOracle = new LoqoMuOracle();
      }
      else if( sfixmuoracle == "probing" )
      {
         FixMuOracle = new ProbingMuOracle(GetPDSystemSolver(jnlst, options, prefix));
      }
      else if( sfixmuoracle == "quality-function" )
      {
         FixMuOracle = new QualityFunctionMuOracle(GetPDSystemSolver(jnlst, options, prefix));
      }
      else
      {
         FixMuOracle = NULL;
      }
      MuUpdate = new AdaptiveMuUpdate(GetRawPtr(LineSearch_), muOracle, FixMuOracle);
   }
   return MuUpdate;
}

SmartPtr<LibraryLoader> AlgorithmBuilder::GetHSLLoader(
   const OptionsList& options,
   const std::string& prefix
)
{
   if( !IsValid(hslloader) )
   {
      IpoptLinearSolver availablesolvers = IpoptGetAvailableLinearSolvers(false);
      IpoptLinearSolver availablesolverslinked = IpoptGetAvailableLinearSolvers(true);

      // we don't have the hsllib option if linked against all hsl routines
      // but then we also don't use the hslloader, so can return NULL
      if( (availablesolverslinked ^ availablesolvers) & IPOPTLINEARSOLVER_ALLHSL )
      {
         std::string libname;
         options.GetStringValue("hsllib", libname, prefix);
         hslloader = new LibraryLoader(libname);
      }
   }

   return hslloader;
}

SmartPtr<LibraryLoader> AlgorithmBuilder::GetPardisoLoader(
   const OptionsList& options,
   const std::string& prefix
)
{
   if( !IsValid(pardisoloader) )
   {
      std::string libname;
      options.GetStringValue("pardisolib", libname, prefix);
      pardisoloader = new LibraryLoader(libname);
   }

   return pardisoloader;
}

} // namespace

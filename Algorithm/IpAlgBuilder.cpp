// Copyright (C) 2004, International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id$
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-09-29

#include "IpAlgBuilder.hpp"

#include "IpStdAugSystemSolver.hpp"
#include "IpAugRestoSystemSolver.hpp"
#include "IpPDFullSpaceSolver.hpp"
#include "IpOptErrorConvCheck.hpp"
#include "IpFilterLineSearch.hpp"
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
#include "IpRestoIterationOutput.hpp"
#include "IpRestoFilterConvCheck.hpp"
#include "IpRestoIterateInitializer.hpp"
#include "IpRestoRestoPhase.hpp"
#include "IpTSymLinearSolver.hpp"
#include "IpMa27TSolverInterface.hpp"
#include "IpMc19TSymScalingMethod.hpp"
#ifdef HAVE_PARDISO
# include "IpPardisoSolverInterface.hpp"
#endif
#ifdef HAVE_TAUCS
# include "IpTAUCSSolverInterface.hpp"
#endif

namespace Ipopt
{
  DBG_SET_VERBOSITY(0);

  DefineIpoptType(AlgorithmBuilder);

  void AlgorithmBuilder::RegisterOptions(SmartPtr<RegisteredOptions> roptions)
  {
    roptions->AddStringOption3(
      "linear_solver",
      "Linear solver used for step computations.",
      "ma27",
      "ma27", "use the Harwell routine MA27",
      "pardiso", "use the Pardiso package",
      "taucs", "use TAUCS package",
      "Determines which linear algebra package is to be used for the "
      "solution of the linear system from which the search directions is "
      "obtained.  Note that depending on your Ipopt installation, not all "
      "options might be available.");
    roptions->AddStringOption2(
      "linear_system_scaling",
      "Method for scaling the linear system.",
      "none",
      "none", "no scaling will be performed",
      "mc19", "use the Harwell routine mc19",
      "Determines which method should be use to compute symmetric scaling "
      "factors for the augmented system.");
    roptions->AddStringOption2(
      "mu_strategy",
      "Update strategy for barrier parameter.",
      "monotone",
      "monotone", "use the monotone (Fiacco-McCormick) strategy",
      "adaptive", "use the adaptive update strategy",
      "Determines which barrier parameter strategy is to be used.");
    roptions->AddStringOption3(
      "mu_oracle",
      "Oracle for a new barrier parameters in the adaptive strategy",
      "probing",
      "probing", "Mehrotra's probing heuristic",
      "loqo", "LOQO's centrality rule",
      "quality_function", "minimize a quality function",
      "Determines how a new barrier parameter is computed in each "
      "\"free-mode\" iteration of the adaptive barrier parameter "
      "strategy. (Only considered if \"adaptive\" is selected for "
      "option \"mu_strategy\".");
    roptions->AddStringOption4(
      "fixed_mu_oracle",
      "Oracle for the barrier parameter when switching to fixed mode.",
      "average_compl",
      "probing", "Mehrotra's probing heuristic",
      "loqo", "LOQO's centrality rule",
      "quality_function", "minimize a quality function",
      "average_compl", "base on current average complementarity",
      "Determines how the first value of the barrier parameter should be "
      "computed when switching to the \"monotone mode\" in the adaptive "
      "strategy. (Only considered if \"adaptive\" is selected for option "
      "\"mu_strategy\".");
    roptions->AddStringOption2(
      "warm_start_init_point",
      "Warm-start for initial point", "no",
      "no", "do not use the warm start initialization",
      "yes", "use the warm start initialization",
      "Indicates whether this optimization should use a warm start "
      "initialization, where values of primal and dual variables are "
      "given (e.g., from a previous optimization of a related problem.)");
  }

  SmartPtr<IpoptAlgorithm>
  AlgorithmBuilder::BuildBasicAlgorithm(const Journalist& jnlst,
                                        const OptionsList& options,
                                        const std::string& prefix)
  {
    DBG_START_FUN("AlgorithmBuilder::BuildBasicAlgorithm",
                  dbg_verbosity);
    // Create the convergence check
    SmartPtr<ConvergenceCheck> convCheck =
      new OptimalityErrorConvergenceCheck();

    // Create the solvers that will be used by the main algorithm
    SmartPtr<TSymScalingMethod> ScalingMethod;
    std::string linear_system_scaling;
    options.GetValue("linear_system_scaling",
                     linear_system_scaling, prefix);
    if (linear_system_scaling=="mc19") {
      ScalingMethod = new Mc19TSymScalingMethod();
    }

    SmartPtr<SparseSymLinearSolverInterface> SolverInterface;
    std::string linear_solver;
    options.GetValue("linear_solver", linear_solver, prefix);
    if (linear_solver=="ma27") {
      SolverInterface = new Ma27TSolverInterface();
    }
    else if (linear_solver=="pardiso") {
#ifdef HAVE_PARDISO
      SolverInterface = new PardisoSolverInterface();
#else

      THROW_EXCEPTION(OptionsList::OPTION_OUT_OF_RANGE,
                      "Selected solver Pardiso not available.");
#endif

    }
    else if (linear_solver=="taucs") {
#ifdef HAVE_TAUCS
      SolverInterface = new TAUCSSolverInterface();
#else

      ASSERT_EXCEPTION(false,
                       OptionsList::OPTION_OUT_OF_RANGE,
                       "Selected solver TAUCS not available.");
#endif

    }

    SmartPtr<SymLinearSolver> ScaledSolver =
      new TSymLinearSolver(SolverInterface, ScalingMethod);

    SmartPtr<AugSystemSolver> AugSolver =
      //        = new AugTSystemSolver(*Ma27Solver);
      new StdAugSystemSolver(*ScaledSolver);
    SmartPtr<PDSystemSolver> PDSolver =
      new PDFullSpaceSolver(*AugSolver);

    // Create the object for initializing the iterates
    // Initialization object
    SmartPtr<EqMultiplierCalculator> EqMultCalculator =
      new LeastSquareMultipliers(*AugSolver);
    SmartPtr<IterateInitializer> IterInitializer;
    bool warm_start_init_point;
    std::string warm_start_option;
    options.GetValue("warm_start_init_point", warm_start_option, prefix);
    warm_start_init_point = (warm_start_option == "yes");

    if (warm_start_init_point) {
      IterInitializer = new WarmStartIterateInitializer();
    }
    else {
      IterInitializer = new DefaultIterateInitializer(EqMultCalculator);
    }

    // Solver for the restoration phase
    SmartPtr<AugSystemSolver> resto_AugSolver =
      new AugRestoSystemSolver(*AugSolver);
    SmartPtr<PDSystemSolver> resto_PDSolver =
      new PDFullSpaceSolver(*resto_AugSolver);

    // Convergence check in the restoration phase
    SmartPtr<RestoFilterConvergenceCheck> resto_convCheck =
      new RestoFilterConvergenceCheck();

    // Line search method for the restoration phase
    SmartPtr<RestoRestorationPhase> resto_resto =
      new RestoRestorationPhase();
    SmartPtr<FilterLineSearch> resto_LineSearch =
      new FilterLineSearch(GetRawPtr(resto_resto), GetRawPtr(resto_PDSolver),
                           GetRawPtr(resto_convCheck));

    // Create the mu update that will be used by the restoration phase
    // algorithm
    SmartPtr<MuUpdate> resto_MuUpdate;
    std::string resto_smuupdate;
    options.GetValue("mu_strategy", resto_smuupdate, "resto."+prefix);

    std::string resto_smuoracle;
    std::string resto_sfixmuoracle;
    if (resto_smuupdate=="adaptive" ) {
      options.GetValue("mu_oracle", resto_smuoracle, "resto."+prefix);
      options.GetValue("fixed_mu_oracle", resto_sfixmuoracle, "resto."+prefix);
    }

    if (resto_smuupdate=="monotone" ) {
      resto_MuUpdate = new MonotoneMuUpdate(GetRawPtr(resto_LineSearch));
    }
    else if (resto_smuupdate=="adaptive") {
      SmartPtr<MuOracle> resto_MuOracle;
      if (resto_smuoracle=="loqo") {
        resto_MuOracle = new LoqoMuOracle();
      }
      else if (resto_smuoracle=="probing") {
        resto_MuOracle = new ProbingMuOracle(resto_PDSolver);
      }
      else if (resto_smuoracle=="quality_function") {
        resto_MuOracle = new QualityFunctionMuOracle(resto_PDSolver);
      }
      SmartPtr<MuOracle> resto_FixMuOracle;
      if (resto_sfixmuoracle=="loqo") {
        resto_FixMuOracle = new LoqoMuOracle();
      }
      else if (resto_sfixmuoracle=="probing") {
        resto_FixMuOracle = new ProbingMuOracle(resto_PDSolver);
      }
      else if (resto_sfixmuoracle=="quality_function") {
        resto_FixMuOracle = new QualityFunctionMuOracle(resto_PDSolver);
      }
      else {
        resto_FixMuOracle = NULL;
      }
      resto_MuUpdate =
        new AdaptiveMuUpdate(GetRawPtr(resto_LineSearch),
                             resto_MuOracle, resto_FixMuOracle);
    }

    // Initialization of the iterates for the restoration phase
    SmartPtr<EqMultiplierCalculator> resto_EqMultCalculator =
      new LeastSquareMultipliers(*resto_AugSolver);
    SmartPtr<IterateInitializer> resto_IterInitializer =
      new RestoIterateInitializer(resto_EqMultCalculator);

    // Create the object for the iteration output during restoration
    SmartPtr<OrigIterationOutput> resto_OrigIterOutput = NULL;
    //   new OrigIterationOutput();
    SmartPtr<IterationOutput> resto_IterOutput =
      new RestoIterationOutput(resto_OrigIterOutput);

    // Put together the overall restoration phase IP algorithm
    SmartPtr<IpoptAlgorithm> resto_alg =
      new IpoptAlgorithm(resto_PDSolver,
                         GetRawPtr(resto_LineSearch),
                         GetRawPtr(resto_MuUpdate),
                         GetRawPtr(resto_convCheck),
                         resto_IterInitializer,
                         resto_IterOutput);

    // Set the restoration phase
    SmartPtr<RestorationPhase> resto_phase =
      new MinC_1NrmRestorationPhase(*resto_alg, EqMultCalculator);

    // Create the line search to be used by the main algorithm
    SmartPtr<FilterLineSearch> lineSearch =
      new FilterLineSearch(GetRawPtr(resto_phase), GetRawPtr(PDSolver),
                           convCheck);

    // The following cross reference is not good: We have to store a
    // pointer to the lineSearch object in resto_convCheck as a
    // non-SmartPtr to make sure that things are properly deleted when
    // the IpoptAlgorithm return by the Builder is destructed.
    resto_convCheck->SetOrigFilterLineSearch(*lineSearch);

    // Create the mu update that will be used by the main algorithm
    SmartPtr<MuUpdate> MuUpdate;
    std::string smuupdate;
    options.GetValue("mu_strategy", smuupdate, prefix);
    std::string smuoracle;
    std::string sfixmuoracle;
    if (smuupdate=="adaptive" ) {
      options.GetValue("mu_oracle", smuoracle, prefix);
      options.GetValue("fixed_mu_oracle", sfixmuoracle, prefix);
    }

    if (smuupdate=="monotone" ) {
      MuUpdate = new MonotoneMuUpdate(GetRawPtr(lineSearch));
    }
    else if (smuupdate=="adaptive") {
      SmartPtr<MuOracle> muOracle;
      if (smuoracle=="loqo") {
        muOracle = new LoqoMuOracle();
      }
      else if (smuoracle=="probing") {
        muOracle = new ProbingMuOracle(PDSolver);
      }
      else if (smuoracle=="quality_function") {
        muOracle = new QualityFunctionMuOracle(PDSolver);
      }
      SmartPtr<MuOracle> FixMuOracle;
      if (sfixmuoracle=="loqo") {
        FixMuOracle = new LoqoMuOracle();
      }
      else if (sfixmuoracle=="probing") {
        FixMuOracle = new ProbingMuOracle(PDSolver);
      }
      else if (sfixmuoracle=="quality_function") {
        FixMuOracle = new QualityFunctionMuOracle(PDSolver);
      }
      else {
        FixMuOracle = NULL;
      }
      MuUpdate = new AdaptiveMuUpdate(GetRawPtr(lineSearch),
                                      muOracle, FixMuOracle);
    }

    // Create the object for the iteration output
    SmartPtr<IterationOutput> IterOutput =
      new OrigIterationOutput();

    // Create the main algorithm
    SmartPtr<IpoptAlgorithm> alg =
      new IpoptAlgorithm(PDSolver,
                         GetRawPtr(lineSearch), MuUpdate,
                         convCheck, IterInitializer, IterOutput);

    return alg;
  }

} // namespace

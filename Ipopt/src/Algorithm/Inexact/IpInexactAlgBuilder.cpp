// Copyright (C) 2008 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// $Id$
//
// Authors:  Andreas Waechter     IBM                  2008-09-05
//            based on IpAlgBuilder.cpp (rev 1311)

#include "IpoptConfig.h"
#include "IpInexactAlgBuilder.hpp"
#include "IpInexactData.hpp"
#include "IpInexactCq.hpp"
#include "IpNLPBoundsRemover.hpp"

#include "IpOptErrorConvCheck.hpp"
#include "IpStdAugSystemSolver.hpp"
#include "IpTSymLinearSolver.hpp"
#include "IpPDPerturbationHandler.hpp"

#include "IpBacktrackingLineSearch.hpp"
#include "IpInexactLSAcceptor.hpp"

#include "IpMonotoneMuUpdate.hpp"
#include "IpAdaptiveMuUpdate.hpp"
#include "IpLoqoMuOracle.hpp"
#include "IpDefaultIterateInitializer.hpp"
#include "IpWarmStartIterateInitializer.hpp"
#include "IpOrigIterationOutput.hpp"
#include "IpUserScaling.hpp"
#include "IpGradientScaling.hpp"
#include "IpEquilibrationScaling.hpp"
#include "IpExactHessianUpdater.hpp"

#include "IpInexactDoglegNormal.hpp"
#include "IpInexactSearchDirCalc.hpp"
#include "IpInexactNewtonNormal.hpp"
#include "IpInexactPDSolver.hpp"

#include "IpMa27TSolverInterface.hpp"
#include "IpMa57TSolverInterface.hpp"
#include "IpMc19TSymScalingMethod.hpp"
#include "IpInexactTSymScalingMethod.hpp"
#include "IpIterativePardisoSolverInterface.hpp"
#include "IpInexactNormalTerminationTester.hpp"
#include "IpInexactPDTerminationTester.hpp"

#ifdef COIN_HAS_HSL
# include "CoinHslConfig.h"
#endif

#ifdef HAVE_WSMP
# include "IpWsmpSolverInterface.hpp"
#endif
#ifdef COIN_HAS_MUMPS
# include "IpMumpsSolverInterface.hpp"
#endif

#ifdef HAVE_LINEARSOLVERLOADER
# include "HSLLoader.h"
# include "PardisoLoader.h"
#endif

namespace Ipopt
{
#if COIN_IPOPT_VERBOSITY > 0
  static const Index dbg_verbosity = 0;
#endif

  InexactAlgorithmBuilder::InexactAlgorithmBuilder()
      :
      AlgorithmBuilder()
  {}

  void InexactAlgorithmBuilder::BuildIpoptObjects(const Journalist& jnlst,
      const OptionsList& options,
      const std::string& prefix,
      const SmartPtr<NLP>& nlp,
      SmartPtr<IpoptNLP>& ip_nlp,
      SmartPtr<IpoptData>& ip_data,
      SmartPtr<IpoptCalculatedQuantities>& ip_cq)
  {
    DBG_ASSERT(prefix == "");

    DBG_ASSERT(dynamic_cast<NLPBoundsRemover*>(GetRawPtr(nlp)));

    // use the original method to get the basic quantites
    AlgorithmBuilder::BuildIpoptObjects(jnlst, options, prefix, nlp,
                                        ip_nlp, ip_data, ip_cq);

    // Now add the objects specific for the inexact step version
    if (ip_data->HaveAddData()) {
      THROW_EXCEPTION(OPTION_INVALID, "The Inexact step computation of Ipopt has been chosen, but some option has been set that requires additional Ipopt data beside the one for the chosen inexact step computation");
    }
    ip_data->SetAddData(new InexactData());

    if (ip_cq->HaveAddCq()) {
      THROW_EXCEPTION(OPTION_INVALID, "The Inexact step computation of Ipopt has been chosen, but some option has been set that requires additional Ipopt calculated quantities beside the one for the chosen inexact step computation");
    }
    ip_cq->SetAddCq(new InexactCq(GetRawPtr(ip_nlp), GetRawPtr(ip_data), GetRawPtr(ip_cq)));
  }

  void InexactAlgorithmBuilder::RegisterOptions(SmartPtr<RegisteredOptions> roptions)
  {
    roptions->SetRegisteringCategory("Linear Solver");
    roptions->AddStringOption2(
      "inexact_linear_system_scaling",
      "Method for scaling the linear system for the inexact approach",
      "slack-based",
      "none", "no scaling will be performed",
      "slack-based", "scale the linear system as in paper",
      "");
  }

  SmartPtr<IpoptAlgorithm>
  InexactAlgorithmBuilder::BuildBasicAlgorithm(const Journalist& jnlst,
      const OptionsList& options,
      const std::string& prefix)
  {
    DBG_START_FUN("InexactAlgorithmBuilder::BuildBasicAlgorithm",
                  dbg_verbosity);

    // Create the convergence check
    SmartPtr<ConvergenceCheck> convCheck =
      new OptimalityErrorConvergenceCheck();

    SmartPtr<InexactNormalTerminationTester> NormalTester;
    SmartPtr<SparseSymLinearSolverInterface> SolverInterface;
    std::string linear_solver;
    options.GetStringValue("linear_solver", linear_solver, prefix);
    if (linear_solver=="ma27") {
#ifndef COINHSL_HAS_MA27
# ifdef HAVE_LINEARSOLVERLOADER
      SolverInterface = new Ma27TSolverInterface();
      char buf[256];
      int rc = LSL_loadHSL(NULL, buf, 255);
      if (rc) {
        std::string errmsg;
        errmsg = "Selected linear solver MA27 not available.\nTried to obtain MA27 from shared library \"";
        errmsg += LSL_HSLLibraryName();
        errmsg += "\", but the following error occured:\n";
        errmsg += buf;
        THROW_EXCEPTION(OPTION_INVALID, errmsg.c_str());
      }
# else
      THROW_EXCEPTION(OPTION_INVALID, "Support for MA27 has not been compiled into Ipopt.");
# endif
#else
      SolverInterface = new Ma27TSolverInterface();
#endif

    }
    else if (linear_solver=="ma57") {
#ifndef COINHSL_HAS_MA57
# ifdef HAVE_LINEARSOLVERLOADER
      SolverInterface = new Ma57TSolverInterface();
      char buf[256];
      int rc = LSL_loadHSL(NULL, buf, 255);
      if (rc) {
        std::string errmsg;
        errmsg = "Selected linear solver MA57 not available.\nTried to obtain MA57 from shared library \"";
        errmsg += LSL_HSLLibraryName();
        errmsg += "\", but the following error occured:\n";
        errmsg += buf;
        THROW_EXCEPTION(OPTION_INVALID, errmsg.c_str());
      }
# else
      THROW_EXCEPTION(OPTION_INVALID, "Support for MA57 has not been compiled into Ipopt.");
# endif
#else
      SolverInterface = new Ma57TSolverInterface();
#endif

    }
    else if (linear_solver=="pardiso") {
      NormalTester = new InexactNormalTerminationTester();
      SmartPtr<IterativeSolverTerminationTester> pd_tester =
        new InexactPDTerminationTester();
#ifndef HAVE_PARDISO
# ifdef HAVE_LINEARSOLVERLOADER
      SolverInterface = new IterativePardisoSolverInterface(*NormalTester, *pd_tester);
      char buf[256];
      int rc = LSL_loadPardisoLib(NULL, buf, 255);
      if (rc) {
        std::string errmsg;
        errmsg = "Selected linear solver Pardiso not available.\nTried to obtain Pardiso from shared library \"";
        errmsg += LSL_PardisoLibraryName();
        errmsg += "\", but the following error occured:\n";
        errmsg += buf;
        THROW_EXCEPTION(OPTION_INVALID, errmsg.c_str());
      }
# else
      THROW_EXCEPTION(OPTION_INVALID, "Support for Pardiso has not been compiled into Ipopt.");
# endif
#else
      SolverInterface = new IterativePardisoSolverInterface(*NormalTester, *pd_tester);
#endif

    }
    else if (linear_solver=="wsmp") {
#ifdef HAVE_WSMP
      SolverInterface = new WsmpSolverInterface();
#else

      THROW_EXCEPTION(OPTION_INVALID,
                      "Selected linear solver WSMP not available.");
#endif

    }
    else if (linear_solver=="mumps") {
#ifdef COIN_HAS_MUMPS
      SolverInterface = new MumpsSolverInterface();
#else

      THROW_EXCEPTION(OPTION_INVALID,
                      "Selected linear solver MUMPS not available.");
#endif

    }
    else {
      THROW_EXCEPTION(OPTION_INVALID,
                      "Inexact version not available for this selection of linear solver.");
    }

    SmartPtr<TSymScalingMethod> ScalingMethod;

    std::string inexact_linear_system_scaling;
    options.GetStringValue("inexact_linear_system_scaling",
                           inexact_linear_system_scaling, prefix);
    if (inexact_linear_system_scaling=="slack-based") {
      ScalingMethod = new InexactTSymScalingMethod();
    }

    SmartPtr<SymLinearSolver> ScaledSolver =
      new TSymLinearSolver(SolverInterface, ScalingMethod);

    SmartPtr<AugSystemSolver> AugSolver =
      new StdAugSystemSolver(*ScaledSolver);

    // Create the object for initializing the iterates Initialization
    // object.  We include both the warm start and the defaut
    // initializer, so that the warm start options can be activated
    // without having to rebuild the algorithm
    SmartPtr<IterateInitializer> WarmStartInitializer =
      new WarmStartIterateInitializer();
    SmartPtr<IterateInitializer> IterInitializer =
      new DefaultIterateInitializer(NULL, WarmStartInitializer, NULL);

    // Create the line search to be used by the main algorithm
    SmartPtr<BacktrackingLSAcceptor> LSacceptor =
      new InexactLSAcceptor();
    SmartPtr<LineSearch> lineSearch =
      new BacktrackingLineSearch(LSacceptor, NULL, convCheck);


    // Create the mu update that will be used by the main algorithm
    SmartPtr<MuUpdate> MuUpdate;
    std::string smuupdate;
    options.GetStringValue("mu_strategy", smuupdate, prefix);
    if (smuupdate=="monotone" ) {
      MuUpdate = new MonotoneMuUpdate(GetRawPtr(lineSearch));
    }
    else if (smuupdate=="adaptive") {
      // for now, we only allow Loqo oracle since it does not require
      // linear system solve
      SmartPtr<MuOracle> muOracle = new LoqoMuOracle();
      SmartPtr<MuOracle> FixMuOracle= new LoqoMuOracle();
      MuUpdate = new AdaptiveMuUpdate(GetRawPtr(lineSearch),
                                      muOracle, FixMuOracle);
    }

    // Create the object for the iteration output
    SmartPtr<IterationOutput> IterOutput =
      new OrigIterationOutput();

    // Get the Hessian updater for the main algorithm
    SmartPtr<HessianUpdater> HessUpdater = new ExactHessianUpdater();

    SmartPtr<InexactNewtonNormalStep> NewtonNormalStep =
      new InexactNewtonNormalStep(AugSolver);

    SmartPtr<InexactNormalStepCalculator> normal_step_calculator =
      new InexactDoglegNormalStep(NewtonNormalStep, NormalTester);

    SmartPtr<PDPerturbationHandler> perturbHandler =
      new PDPerturbationHandler();

    SmartPtr<InexactPDSolver> inexact_pd_solver =
      new InexactPDSolver(*AugSolver, *perturbHandler);

    SmartPtr<SearchDirectionCalculator> SearchDirCalc =
      new InexactSearchDirCalculator(normal_step_calculator, inexact_pd_solver);

    // Create the main algorithm
    SmartPtr<IpoptAlgorithm> alg =
      new IpoptAlgorithm(SearchDirCalc,
                         GetRawPtr(lineSearch), MuUpdate,
                         convCheck, IterInitializer, IterOutput,
                         HessUpdater);

    return alg;
  }

  void
  AddInexactDefaultOptions(OptionsList& options_list)
  {
    options_list.SetIntegerValueIfUnset("max_soc", 0);
    options_list.SetStringValueIfUnset("constraint_violation_norm_type",
                                       "2-norm");
    options_list.SetNumericValueIfUnset("constr_mult_init_max", 0.);

    // TODO: Find out about the following:
    //options_list.SetNumericValueIfUnset("bound_relax_factor", 0.);
    options_list.SetNumericValueIfUnset("kappa_d", 0.);
    options_list.SetStringValueIfUnset("linear_solver", "pardiso");
    options_list.SetStringValue("linear_scaling_on_demand", "no");
    options_list.SetStringValue("replace_bounds", "yes");
  }
} // namespace

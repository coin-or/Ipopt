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
#include "IpNonmonotoneMuUpdate.hpp"
#include "IpLoqoMuOracle.hpp"
#include "IpProbingMuOracle.hpp"
#include "IpOptProbingMuOracle.hpp"
#include "IpRestoMinC_1Nrm.hpp"
#include "IpLeastSquareMults.hpp"
#include "IpDefaultIterateInitializer.hpp"
#include "IpOrigIterationOutput.hpp"
#include "IpRestoIterationOutput.hpp"
#include "IpRestoFilterConvCheck.hpp"
#include "IpRestoIterateInitializer.hpp"
#include "IpRestoRestoPhase.hpp"
#include "IpTSymLinearSolver.hpp"
#include "IpMa27TSolverInterface.hpp"
#include "IpMc19TSymScalingMethod.hpp"
#include "IpPardisoSolverInterface.hpp"

namespace Ipopt
{
  DBG_SET_VERBOSITY(0);

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
    std::string scaling_method;
    if (options.GetValue("scaling_method", scaling_method, prefix)) {
      ASSERT_EXCEPTION(scaling_method=="mc19" || scaling_method=="none",
                       OptionsList::OPTION_OUT_OF_RANGE,
                       "Option \"scaling_method\" has invalid value.");
    }
    else {
      scaling_method = "mc19";
    }
    if (scaling_method=="mc19") {
      ScalingMethod = new Mc19TSymScalingMethod();
    }

    SmartPtr<SparseSymLinearSolverInterface> SolverInterface;
    std::string linear_solver;
    if (options.GetValue("linear_solver", linear_solver, prefix)) {
      ASSERT_EXCEPTION(linear_solver=="ma27" || linear_solver=="pardiso",
                       OptionsList::OPTION_OUT_OF_RANGE,
                       "Option \"linear_solver\" has invalid value.");
    }
    else {
      linear_solver = "ma27";
    }
    if (linear_solver=="ma27") {
      SolverInterface = new Ma27TSolverInterface();
    }
    else if (linear_solver=="pardiso") {
#ifdef HAVE_PARDISO
      SolverInterface = new PardisoSolverInterface();
#else

      ASSERT_EXCEPTION(false,
                       OptionsList::OPTION_OUT_OF_RANGE,
                       "Selected solver Pardiso not available.");
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
    SmartPtr<IterateInitializer> IterInitializer =
      new DefaultIterateInitializer(EqMultCalculator);

    // Solver for the restoration phase
    SmartPtr<AugSystemSolver> resto_AugSolver =
      new AugRestoSystemSolver(*AugSolver);
    SmartPtr<PDSystemSolver> resto_PDSolver =
      new PDFullSpaceSolver(*resto_AugSolver);

    // Line search method for the restoration phase
    SmartPtr<RestoRestorationPhase> resto_resto =
      new RestoRestorationPhase();
    SmartPtr<FilterLineSearch> resto_LineSearch =
      new FilterLineSearch(GetRawPtr(resto_resto), GetRawPtr(resto_PDSolver));

    // Create the mu update that will be used by the restoration phase
    // algorithm
    SmartPtr<MuUpdate> resto_MuUpdate;
    std::string resto_smuupdate;
    if (options.GetValue("muupdate", resto_smuupdate, "resto."+prefix)) {
      ASSERT_EXCEPTION(resto_smuupdate=="monotone" || resto_smuupdate=="nonmonotone",
                       OptionsList::OPTION_OUT_OF_RANGE,
                       "Option \"resto_muupdate\" has invalid value.");
    }
    else {
      resto_smuupdate = "monotone";
    }
    std::string resto_smuoracle;
    std::string resto_sfixmuoracle;
    if (resto_smuupdate=="nonmonotone" ) {
      if (options.GetValue("muoracle", resto_smuoracle, "resto."+prefix)) {
        ASSERT_EXCEPTION(resto_smuoracle=="loqo" || resto_smuoracle=="probing" || resto_smuoracle=="optprobing",
                         OptionsList::OPTION_OUT_OF_RANGE,
                         "Option \"resto_muoracle\" has invalid value.");
      }
      else {
        resto_smuoracle = "probing";
      }

      if (options.GetValue("fixmuoracle", resto_sfixmuoracle, "resto."+prefix)) {
        ASSERT_EXCEPTION(resto_sfixmuoracle=="loqo" || resto_sfixmuoracle=="probing" || resto_sfixmuoracle=="optprobing" || resto_sfixmuoracle=="average_compl",
                         OptionsList::OPTION_OUT_OF_RANGE,
                         "Option \"resto_fixmuoracle\" has invalid value.");
      }
      else {
        resto_sfixmuoracle = "average_compl";
      }
    }

    if (resto_smuupdate=="monotone" ) {
      resto_MuUpdate = new MonotoneMuUpdate(GetRawPtr(resto_LineSearch));
    }
    else if (resto_smuupdate=="nonmonotone") {
      SmartPtr<MuOracle> resto_MuOracle;
      if (resto_smuoracle=="loqo") {
        resto_MuOracle = new LoqoMuOracle();
      }
      else if (resto_smuoracle=="probing") {
        resto_MuOracle = new ProbingMuOracle(resto_PDSolver);
      }
      else if (resto_smuoracle=="optprobing") {
        resto_MuOracle = new OptProbingMuOracle(resto_PDSolver);
      }
      SmartPtr<MuOracle> resto_FixMuOracle;
      if (resto_sfixmuoracle=="loqo") {
        resto_FixMuOracle = new LoqoMuOracle();
      }
      else if (resto_sfixmuoracle=="probing") {
        resto_FixMuOracle = new ProbingMuOracle(resto_PDSolver);
      }
      else if (resto_sfixmuoracle=="optprobing") {
        resto_FixMuOracle = new OptProbingMuOracle(resto_PDSolver);
      }
      else {
        resto_FixMuOracle = NULL;
      }
      resto_MuUpdate =
        new NonmonotoneMuUpdate(GetRawPtr(resto_LineSearch),
                                resto_MuOracle, resto_FixMuOracle);
    }

    // Convergence check in the restoration phase
    SmartPtr<RestoFilterConvergenceCheck> resto_convCheck =
      new RestoFilterConvergenceCheck();

    // Initialization of the iterates for the restoration phase
    SmartPtr<EqMultiplierCalculator> resto_EqMultCalculator =
      new LeastSquareMultipliers(*resto_AugSolver);
    SmartPtr<IterateInitializer> resto_IterInitializer =
      //      new DefaultIterateInitializer(resto_EqMultCalculator); //TODO
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
      new FilterLineSearch(GetRawPtr(resto_phase), GetRawPtr(PDSolver));

    // The following cross reference is not good: We have to store a
    // pointer to the lineSearch object in resto_convCheck as a
    // non-SmartPtr to make sure that things are properly deleted when
    // the IpoptAlgorithm return by the Builder is destructed.
    resto_convCheck->SetOrigFilterLineSearch(*lineSearch);

    // Create the mu update that will be used by the main algorithm
    SmartPtr<MuUpdate> MuUpdate;
    std::string smuupdate;
    if (options.GetValue("muupdate", smuupdate, prefix)) {
      ASSERT_EXCEPTION(smuupdate=="monotone" || smuupdate=="nonmonotone",
                       OptionsList::OPTION_OUT_OF_RANGE,
                       "Option \"muupdate\" has invalid value.");
    }
    else {
      smuupdate = "monotone";
    }
    std::string smuoracle;
    std::string sfixmuoracle;
    if (smuupdate=="nonmonotone" ) {
      if (options.GetValue("muoracle", smuoracle, prefix)) {
        ASSERT_EXCEPTION(smuoracle=="loqo" || smuoracle=="probing" || smuoracle=="optprobing",
                         OptionsList::OPTION_OUT_OF_RANGE,
                         "Option \"muoracle\" has invalid value.");
      }
      else {
        smuoracle = "probing";
      }

      if (options.GetValue("fixmuoracle", sfixmuoracle, prefix)) {
        ASSERT_EXCEPTION(sfixmuoracle=="loqo" || sfixmuoracle=="probing" || sfixmuoracle=="optprobing" || sfixmuoracle=="average_compl",
                         OptionsList::OPTION_OUT_OF_RANGE,
                         "Option \"fixmuoracle\" has invalid value.");
      }
      else {
        sfixmuoracle = "average_compl";
      }
    }

    if (smuupdate=="monotone" ) {
      MuUpdate = new MonotoneMuUpdate(GetRawPtr(lineSearch));
    }
    else if (smuupdate=="nonmonotone") {
      SmartPtr<MuOracle> muOracle;
      if (smuoracle=="loqo") {
        muOracle = new LoqoMuOracle();
      }
      else if (smuoracle=="probing") {
        muOracle = new ProbingMuOracle(PDSolver);
      }
      else if (smuoracle=="optprobing") {
        muOracle = new OptProbingMuOracle(PDSolver);
      }
      SmartPtr<MuOracle> FixMuOracle;
      if (sfixmuoracle=="loqo") {
        FixMuOracle = new LoqoMuOracle();
      }
      else if (sfixmuoracle=="probing") {
        FixMuOracle = new ProbingMuOracle(PDSolver);
      }
      else if (sfixmuoracle=="optprobing") {
        FixMuOracle = new OptProbingMuOracle(PDSolver);
      }
      else {
        FixMuOracle = NULL;
      }
      MuUpdate = new NonmonotoneMuUpdate(GetRawPtr(lineSearch),
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

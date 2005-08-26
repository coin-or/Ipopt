// Copyright (C) 2005, International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id$
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2005-08-16

#include "IpAlgorithmRegOp.hpp"
#include "IpRegOptions.hpp"

#include "IpAdaptiveMuUpdate.hpp"
#include "IpAlgBuilder.hpp"
#include "IpDefaultIterateInitializer.hpp"
#include "IpFilterLineSearch.hpp"
#include "IpGradientScaling.hpp"
#include "IpIpoptAlg.hpp"
#include "IpIpoptCalculatedQuantities.hpp"
#include "IpIpoptData.hpp"
#include "IpMonotoneMuUpdate.hpp"
#include "IpNLPScaling.hpp"
#include "IpOptErrorConvCheck.hpp"
#include "IpOrigIpoptNLP.hpp"
#include "IpOrigIterationOutput.hpp"
#include "IpPDFullSpaceSolver.hpp"
#include "IpPDPerturbationHandler.hpp"
#include "IpProbingMuOracle.hpp"
#include "IpQualityFunctionMuOracle.hpp"
#include "IpRestoFilterConvCheck.hpp"
#include "IpRestoIpoptNLP.hpp"
#include "IpRestoMinC_1Nrm.hpp"
#include "IpWarmStartIterateInitializer.hpp"


namespace Ipopt
{

  void RegisterOptions_Algorithm(const SmartPtr<RegisteredOptions>& roptions)
  {
    roptions->SetRegisteringCategory("Mu Update");
    AdaptiveMuUpdate::RegisterOptions(roptions);
    roptions->SetRegisteringCategory("Initialization");
    DefaultIterateInitializer::RegisterOptions(roptions);
    roptions->SetRegisteringCategory("Main Algorithm");
    AlgorithmBuilder::RegisterOptions(roptions);
    roptions->SetRegisteringCategory("Line Search");
    FilterLineSearch::RegisterOptions(roptions);
    roptions->SetRegisteringCategory("NLP Scaling");
    StandardScalingBase::RegisterOptions(roptions);
    roptions->SetRegisteringCategory("NLP Scaling");
    GradientScaling::RegisterOptions(roptions);
    roptions->SetRegisteringCategory("Uncategorized");
    IpoptAlgorithm::RegisterOptions(roptions);
    roptions->SetRegisteringCategory("Uncategorized");
    IpoptData::RegisterOptions(roptions);
    roptions->SetRegisteringCategory("Uncategorized");
    IpoptCalculatedQuantities::RegisterOptions(roptions);
    roptions->SetRegisteringCategory("Mu Update");
    MonotoneMuUpdate::RegisterOptions(roptions);
    roptions->SetRegisteringCategory("Convergence");
    OptimalityErrorConvergenceCheck::RegisterOptions(roptions);
    roptions->SetRegisteringCategory("NLP");
    OrigIpoptNLP::RegisterOptions(roptions);
    roptions->SetRegisteringCategory("Output");
    OrigIterationOutput::RegisterOptions(roptions);
    roptions->SetRegisteringCategory("Step Calculation");
    PDFullSpaceSolver::RegisterOptions(roptions);
    roptions->SetRegisteringCategory("Step Calculation");
    PDPerturbationHandler::RegisterOptions(roptions);
    roptions->SetRegisteringCategory("Mu Update");
    ProbingMuOracle::RegisterOptions(roptions);
    roptions->SetRegisteringCategory("Mu Update");
    QualityFunctionMuOracle::RegisterOptions(roptions);
    roptions->SetRegisteringCategory("Restoration");
    RestoFilterConvergenceCheck::RegisterOptions(roptions);
    roptions->SetRegisteringCategory("Restoration");
    RestoIpoptNLP::RegisterOptions(roptions);
    roptions->SetRegisteringCategory("Uncategorized");
    roptions->SetRegisteringCategory("Restoration");
    MinC_1NrmRestorationPhase::RegisterOptions(roptions);
    roptions->SetRegisteringCategory("Warm Start");
    WarmStartIterateInitializer::RegisterOptions(roptions);
  }

} // namespace Ipopt

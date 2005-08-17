// Copyright (C) 2004, International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id: IpMa27TSolverInterface.hpp 430 2005-08-10 00:19:54Z andreasw $
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
    AdaptiveMuUpdate::RegisterOptions(roptions);
    AlgorithmBuilder::RegisterOptions(roptions);
    DefaultIterateInitializer::RegisterOptions(roptions);
    FilterLineSearch::RegisterOptions(roptions);
    GradientScaling::RegisterOptions(roptions);
    IpoptAlgorithm::RegisterOptions(roptions);
    IpoptCalculatedQuantities::RegisterOptions(roptions);
    IpoptData::RegisterOptions(roptions);
    MonotoneMuUpdate::RegisterOptions(roptions);
    StandardScalingBase::RegisterOptions(roptions);
    OptimalityErrorConvergenceCheck::RegisterOptions(roptions);
    OrigIpoptNLP::RegisterOptions(roptions);
    PDFullSpaceSolver::RegisterOptions(roptions);
    PDPerturbationHandler::RegisterOptions(roptions);
    ProbingMuOracle::RegisterOptions(roptions);
    QualityFunctionMuOracle::RegisterOptions(roptions);
    RestoFilterConvergenceCheck::RegisterOptions(roptions);
    RestoIpoptNLP::RegisterOptions(roptions);
    MinC_1NrmRestorationPhase::RegisterOptions(roptions);
    WarmStartIterateInitializer::RegisterOptions(roptions);
  }

} // namespace Ipopt

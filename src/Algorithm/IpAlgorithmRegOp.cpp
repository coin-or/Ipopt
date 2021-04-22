// Copyright (C) 2005, 2008 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2005-08-16

#include "IpAlgorithmRegOp.hpp"
#include "IpRegOptions.hpp"

#include "IpAdaptiveMuUpdate.hpp"
#include "IpAlgBuilder.hpp"
#include "IpDefaultIterateInitializer.hpp"
#include "IpBacktrackingLineSearch.hpp"
#include "IpFilterLSAcceptor.hpp"
#include "IpGradientScaling.hpp"
#include "IpEquilibrationScaling.hpp"
#include "IpIpoptAlg.hpp"
#include "IpIpoptCalculatedQuantities.hpp"
#include "IpIpoptData.hpp"
#include "IpMonotoneMuUpdate.hpp"
#include "IpNLPScaling.hpp"
#include "IpOptErrorConvCheck.hpp"
#include "IpOrigIpoptNLP.hpp"
#include "IpOrigIterationOutput.hpp"
#include "IpLimMemQuasiNewtonUpdater.hpp"
#include "IpPDFullSpaceSolver.hpp"
#include "IpPDPerturbationHandler.hpp"
#include "IpPDSearchDirCalc.hpp"
#include "IpPenaltyLSAcceptor.hpp"
#include "IpProbingMuOracle.hpp"
#include "IpQualityFunctionMuOracle.hpp"
#include "IpRestoConvCheck.hpp"
#include "IpRestoFilterConvCheck.hpp"
#include "IpRestoIpoptNLP.hpp"
#include "IpRestoMinC_1Nrm.hpp"
#include "IpRestoPenaltyConvCheck.hpp"
#include "IpWarmStartIterateInitializer.hpp"

namespace Ipopt
{

void RegisterOptions_Algorithm(
   const SmartPtr<RegisteredOptions>& roptions
)
{
   roptions->SetRegisteringCategory("Barrier Parameter Update", 390000);
   AdaptiveMuUpdate::RegisterOptions(roptions);
   roptions->SetRegisteringCategory("Initialization", 460000);
   DefaultIterateInitializer::RegisterOptions(roptions);
   roptions->SetRegisteringCategory("Main Algorithm", -50000);
   AlgorithmBuilder::RegisterOptions(roptions);
   roptions->SetRegisteringCategory("Line Search", 380000);
   BacktrackingLineSearch::RegisterOptions(roptions);
   roptions->SetRegisteringCategory("Line Search");
   FilterLSAcceptor::RegisterOptions(roptions);
   roptions->SetRegisteringCategory("Line Search");
   PenaltyLSAcceptor::RegisterOptions(roptions);
   roptions->SetRegisteringCategory("NLP Scaling", 480000);
   StandardScalingBase::RegisterOptions(roptions);
   roptions->SetRegisteringCategory("NLP Scaling");
   GradientScaling::RegisterOptions(roptions);
   roptions->SetRegisteringCategory("NLP Scaling");
   EquilibrationScaling::RegisterOptions(roptions);
   roptions->SetRegisteringCategory("");
   IpoptAlgorithm::RegisterOptions(roptions);
   roptions->SetRegisteringCategory("");
   IpoptData::RegisterOptions(roptions);
   roptions->SetRegisteringCategory("");
   IpoptCalculatedQuantities::RegisterOptions(roptions);
   roptions->SetRegisteringCategory("Hessian Approximation", 280000);
   LimMemQuasiNewtonUpdater::RegisterOptions(roptions);
   roptions->SetRegisteringCategory("Barrier Parameter Update", 390000);
   MonotoneMuUpdate::RegisterOptions(roptions);
   roptions->SetRegisteringCategory("Termination", 490000);
   OptimalityErrorConvergenceCheck::RegisterOptions(roptions);
   roptions->SetRegisteringCategory("NLP", 470000);
   OrigIpoptNLP::RegisterOptions(roptions);
   roptions->SetRegisteringCategory("Output", 900000);
   OrigIterationOutput::RegisterOptions(roptions);
   roptions->SetRegisteringCategory("Step Calculation", 350000);
   PDSearchDirCalculator::RegisterOptions(roptions);
   roptions->SetRegisteringCategory("Linear Solver", 360000);
   PDFullSpaceSolver::RegisterOptions(roptions);
   roptions->SetRegisteringCategory("Step Calculation");
   PDPerturbationHandler::RegisterOptions(roptions);
   roptions->SetRegisteringCategory("Barrier Parameter Update", 390000);
   ProbingMuOracle::RegisterOptions(roptions);
   roptions->SetRegisteringCategory("Barrier Parameter Update");
   QualityFunctionMuOracle::RegisterOptions(roptions);
   roptions->SetRegisteringCategory("Restoration Phase", 340000);
   RestoConvergenceCheck::RegisterOptions(roptions);
   roptions->SetRegisteringCategory("Restoration Phase");
   RestoFilterConvergenceCheck::RegisterOptions(roptions);
   roptions->SetRegisteringCategory("Restoration Phase");
   RestoIpoptNLP::RegisterOptions(roptions);
   roptions->SetRegisteringCategory("Restoration Phase");
   RestoPenaltyConvergenceCheck::RegisterOptions(roptions);
   roptions->SetRegisteringCategory("Restoration Phase");
   MinC_1NrmRestorationPhase::RegisterOptions(roptions);
   roptions->SetRegisteringCategory("Warm Start", 370000);
   WarmStartIterateInitializer::RegisterOptions(roptions);
}

} // namespace Ipopt

// Copyright (C) 2005 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id: IpSolveStatistics.cpp 507 2005-08-26 18:33:43Z andreasw $
//
// Authors:  Andreas Waechter          IBM    2005-09-19

#include "IpTimingStatistics.hpp"

namespace Ipopt
{
  void
  TimingStatistics::PrintAllTimingStatistics(
    Journalist& jnlst,
    EJournalLevel level,
    EJournalCategory category) const
  {
    if (!jnlst.ProduceOutput(level, category))
      return;

    jnlst.Printf(level, category,
                 "OverallAlgorithm....................: %10.3f\n",
                 OverallAlgorithm_.TotalTime());
    jnlst.Printf(level, category,
                 " PrintProblemStatistics.............: %10.3f\n",
                 PrintProblemStatistics_.TotalTime());
    jnlst.Printf(level, category,
                 " InitializeIterates.................: %10.3f\n",
                 InitializeIterates_.TotalTime());
    jnlst.Printf(level, category,
                 " ActualizeHessian...................: %10.3f\n",
                 ActualizeHessian_.TotalTime());
    jnlst.Printf(level, category,
                 " OutputIteration....................: %10.3f\n",
                 OutputIteration_.TotalTime());
    jnlst.Printf(level, category,
                 " UpdateBarrierParameter.............: %10.3f\n",
                 UpdateBarrierParameter_.TotalTime());
    jnlst.Printf(level, category,
                 " ComputeSearchDirection.............: %10.3f\n",
                 ComputeSearchDirection_.TotalTime());
    jnlst.Printf(level, category,
                 " ComputeAcceptableTrialPoint........: %10.3f\n",
                 ComputeAcceptableTrialPoint_.TotalTime());
    jnlst.Printf(level, category,
                 " AcceptTrialPoint...................: %10.3f\n",
                 AcceptTrialPoint_.TotalTime());
    jnlst.Printf(level, category,
                 " CheckConvergence...................: %10.3f\n",
                 CheckConvergence_.TotalTime());

    jnlst.Printf(level, category,
                 "PDSystemSolverTotal.................: %10.3f\n",
                 PDSystemSolverTotal_.TotalTime());
    jnlst.Printf(level, category,
                 " PDSystemSolverSolveOnce............: %10.3f\n",
                 PDSystemSolverSolveOnce_.TotalTime());
    jnlst.Printf(level, category,
                 " ComputeResiduals...................: %10.3f\n",
                 ComputeResiduals_.TotalTime());
    jnlst.Printf(level, category,
                 " LinearSystemScaling................: %10.3f\n",
                 LinearSystemScaling_.TotalTime());
    jnlst.Printf(level, category,
                 " LinearSystemSymbolicFactorization..: %10.3f\n",
                 LinearSystemSymbolicFactorization_.TotalTime());
    jnlst.Printf(level, category,
                 " LinearSystemFactorization..........: %10.3f\n",
                 LinearSystemFactorization_.TotalTime());
    jnlst.Printf(level, category,
                 " LinearSystemBackSolve..............: %10.3f\n",
                 LinearSystemBackSolve_.TotalTime());
    jnlst.Printf(level, category,
                 "QualityFunctionSearch...............: %10.3f\n",
                 QualityFunctionSearch_.TotalTime());
    jnlst.Printf(level, category,
                 "TryCorrector........................: %10.3f\n",
                 TryCorrector_.TotalTime());
    jnlst.Printf(level, category,
                 "Task1...............................: %10.3f\n",
                 Task1_.TotalTime());
    jnlst.Printf(level, category,
                 "Task2...............................: %10.3f\n",
                 Task2_.TotalTime());
    jnlst.Printf(level, category,
                 "Task3...............................: %10.3f\n",
                 Task3_.TotalTime());
    jnlst.Printf(level, category,
                 "Task4...............................: %10.3f\n",
                 Task4_.TotalTime());
    jnlst.Printf(level, category,
                 "Task5...............................: %10.3f\n",
                 Task5_.TotalTime());
  }
} // namespace Ipopt

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
                 "OverallAlgorithm....................: %10.2f\n",
                 OverallAlgorithm.TotalTime());
    jnlst.Printf(level, category,
                 " PrintProblemStatistics.............: %10.2f\n",
                 PrintProblemStatistics.TotalTime());
    jnlst.Printf(level, category,
                 " InitializeIterates.................: %10.2f\n",
                 InitializeIterates.TotalTime());
    jnlst.Printf(level, category,
                 " ActualizeHessian...................: %10.2f\n",
                 ActualizeHessian.TotalTime());
    jnlst.Printf(level, category,
                 " OutputIteration....................: %10.2f\n",
                 OutputIteration.TotalTime());
    jnlst.Printf(level, category,
                 " UpdateBarrierParameter.............: %10.2f\n",
                 UpdateBarrierParameter.TotalTime());
    jnlst.Printf(level, category,
                 " ComputeSearchDirection.............: %10.2f\n",
                 ComputeSearchDirection.TotalTime());
    jnlst.Printf(level, category,
                 " ComputeAcceptableTrialPoint........: %10.2f\n",
                 ComputeAcceptableTrialPoint.TotalTime());
    jnlst.Printf(level, category,
                 " AcceptTrialPoint...................: %10.2f\n",
                 AcceptTrialPoint.TotalTime());
    jnlst.Printf(level, category,
                 " CheckConvergence...................: %10.2f\n",
                 CheckConvergence.TotalTime());

    jnlst.Printf(level, category,
                 "PDSystemSolverTotal.................: %10.2f\n",
                 PDSystemSolverTotal.TotalTime());
    jnlst.Printf(level, category,
                 " PDSystemSolverSolveOnce............: %10.2f\n",
                 PDSystemSolverSolveOnce.TotalTime());
    jnlst.Printf(level, category,
                 " LinearSystemScaling................: %10.2f\n",
                 LinearSystemScaling.TotalTime());
    jnlst.Printf(level, category,
                 " LinearSystemSymbolicFactorization..: %10.2f\n",
                 LinearSystemSymbolicFactorization.TotalTime());
    jnlst.Printf(level, category,
                 " LinearSystemFactorization..........: %10.2f\n",
                 LinearSystemFactorization.TotalTime());
    jnlst.Printf(level, category,
                 " LinearSystemBackSolve..............: %10.2f\n",
                 LinearSystemBackSolve.TotalTime());
    jnlst.Printf(level, category,
                 "QualityFunctionSearch...............: %10.2f\n",
                 QualityFunctionSearch.TotalTime());
    jnlst.Printf(level, category,
                 "TryCorrector........................: %10.2f\n",
                 TryCorrector.TotalTime());
    jnlst.Printf(level, category,
                 "Task1...............................: %10.2f\n",
                 Task1.TotalTime());
    jnlst.Printf(level, category,
                 "Task2...............................: %10.2f\n",
                 Task2.TotalTime());
    jnlst.Printf(level, category,
                 "Task3...............................: %10.2f\n",
                 Task3.TotalTime());
    jnlst.Printf(level, category,
                 "Task4...............................: %10.2f\n",
                 Task4.TotalTime());
    jnlst.Printf(level, category,
                 "Task5...............................: %10.2f\n",
                 Task5.TotalTime());
  }
} // namespace Ipopt

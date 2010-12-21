// Copyright (C) 2005, 2009 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// $Id$
//
// Authors:  Andreas Waechter          IBM    2005-09-19

#include "IpTimingStatistics.hpp"

namespace Ipopt
{
  void
  TimingStatistics::ResetTimes()
  {
    OverallAlgorithm_.Reset();
    PrintProblemStatistics_.Reset();
    InitializeIterates_.Reset();
    UpdateHessian_.Reset();
    OutputIteration_.Reset();
    UpdateBarrierParameter_.Reset();
    ComputeSearchDirection_.Reset();
    ComputeAcceptableTrialPoint_.Reset();
    AcceptTrialPoint_.Reset();
    CheckConvergence_.Reset();
    PDSystemSolverTotal_.Reset();
    PDSystemSolverSolveOnce_.Reset();
    ComputeResiduals_.Reset();
    StdAugSystemSolverMultiSolve_.Reset();
    LinearSystemScaling_.Reset();
    LinearSystemSymbolicFactorization_.Reset();
    LinearSystemFactorization_.Reset();
    LinearSystemBackSolve_.Reset();
    LinearSystemStructureConverter_.Reset();
    LinearSystemStructureConverterInit_.Reset();
    QualityFunctionSearch_.Reset();
    TryCorrector_.Reset();
    Task1_.Reset();
    Task2_.Reset();
    Task3_.Reset();
    Task4_.Reset();
    Task5_.Reset();
    Task6_.Reset();
  }

  void
  TimingStatistics::PrintAllTimingStatistics(
    Journalist& jnlst,
    EJournalLevel level,
    EJournalCategory category) const
  {
    if (!jnlst.ProduceOutput(level, category))
      return;

    jnlst.Printf(level, category,
                 "OverallAlgorithm....................: %10.3f (sys: %10.3f wall: %10.3f)\n",
                 OverallAlgorithm_.TotalCpuTime(),
                 OverallAlgorithm_.TotalSysTime(),
                 OverallAlgorithm_.TotalWallclockTime());
    jnlst.Printf(level, category,
                 " PrintProblemStatistics.............: %10.3f (sys: %10.3f wall: %10.3f)\n",
                 PrintProblemStatistics_.TotalCpuTime(),
                 PrintProblemStatistics_.TotalSysTime(),
                 PrintProblemStatistics_.TotalWallclockTime());
    jnlst.Printf(level, category,
                 " InitializeIterates.................: %10.3f (sys: %10.3f wall: %10.3f)\n",
                 InitializeIterates_.TotalCpuTime(),
                 InitializeIterates_.TotalSysTime(),
                 InitializeIterates_.TotalWallclockTime());
    jnlst.Printf(level, category,
                 " UpdateHessian......................: %10.3f (sys: %10.3f wall: %10.3f)\n",
                 UpdateHessian_.TotalCpuTime(),
                 UpdateHessian_.TotalSysTime(),
                 UpdateHessian_.TotalWallclockTime());
    jnlst.Printf(level, category,
                 " OutputIteration....................: %10.3f (sys: %10.3f wall: %10.3f)\n",
                 OutputIteration_.TotalCpuTime(),
                 OutputIteration_.TotalSysTime(),
                 OutputIteration_.TotalWallclockTime());
    jnlst.Printf(level, category,
                 " UpdateBarrierParameter.............: %10.3f (sys: %10.3f wall: %10.3f)\n",
                 UpdateBarrierParameter_.TotalCpuTime(),
                 UpdateBarrierParameter_.TotalSysTime(),
                 UpdateBarrierParameter_.TotalWallclockTime());
    jnlst.Printf(level, category,
                 " ComputeSearchDirection.............: %10.3f (sys: %10.3f wall: %10.3f)\n",
                 ComputeSearchDirection_.TotalCpuTime(),
                 ComputeSearchDirection_.TotalSysTime(),
                 ComputeSearchDirection_.TotalWallclockTime());
    jnlst.Printf(level, category,
                 " ComputeAcceptableTrialPoint........: %10.3f (sys: %10.3f wall: %10.3f)\n",
                 ComputeAcceptableTrialPoint_.TotalCpuTime(),
                 ComputeAcceptableTrialPoint_.TotalSysTime(),
                 ComputeAcceptableTrialPoint_.TotalWallclockTime());
    jnlst.Printf(level, category,
                 " AcceptTrialPoint...................: %10.3f (sys: %10.3f wall: %10.3f)\n",
                 AcceptTrialPoint_.TotalCpuTime(),
                 AcceptTrialPoint_.TotalSysTime(),
                 AcceptTrialPoint_.TotalWallclockTime());
    jnlst.Printf(level, category,
                 " CheckConvergence...................: %10.3f (sys: %10.3f wall: %10.3f)\n",
                 CheckConvergence_.TotalCpuTime(),
                 CheckConvergence_.TotalSysTime(),
                 CheckConvergence_.TotalWallclockTime());

    jnlst.Printf(level, category,
                 "PDSystemSolverTotal.................: %10.3f (sys: %10.3f wall: %10.3f)\n",
                 PDSystemSolverTotal_.TotalCpuTime(),
                 PDSystemSolverTotal_.TotalSysTime(),
                 PDSystemSolverTotal_.TotalWallclockTime());
    jnlst.Printf(level, category,
                 " PDSystemSolverSolveOnce............: %10.3f (sys: %10.3f wall: %10.3f)\n",
                 PDSystemSolverSolveOnce_.TotalCpuTime(),
                 PDSystemSolverSolveOnce_.TotalSysTime(),
                 PDSystemSolverSolveOnce_.TotalWallclockTime());
    jnlst.Printf(level, category,
                 " ComputeResiduals...................: %10.3f (sys: %10.3f wall: %10.3f)\n",
                 ComputeResiduals_.TotalCpuTime(),
                 ComputeResiduals_.TotalSysTime(),
                 ComputeResiduals_.TotalWallclockTime());
    jnlst.Printf(level, category,
                 " StdAugSystemSolverMultiSolve.......: %10.3f (sys: %10.3f wall: %10.3f)\n",
                 StdAugSystemSolverMultiSolve_.TotalCpuTime(),
                 StdAugSystemSolverMultiSolve_.TotalSysTime(),
                 StdAugSystemSolverMultiSolve_.TotalWallclockTime());
    jnlst.Printf(level, category,
                 " LinearSystemScaling................: %10.3f (sys: %10.3f wall: %10.3f)\n",
                 LinearSystemScaling_.TotalCpuTime(),
                 LinearSystemScaling_.TotalSysTime(),
                 LinearSystemScaling_.TotalWallclockTime());
    jnlst.Printf(level, category,
                 " LinearSystemSymbolicFactorization..: %10.3f (sys: %10.3f wall: %10.3f)\n",
                 LinearSystemSymbolicFactorization_.TotalCpuTime(),
                 LinearSystemSymbolicFactorization_.TotalSysTime(),
                 LinearSystemSymbolicFactorization_.TotalWallclockTime());
    jnlst.Printf(level, category,
                 " LinearSystemFactorization..........: %10.3f (sys: %10.3f wall: %10.3f)\n",
                 LinearSystemFactorization_.TotalCpuTime(),
                 LinearSystemFactorization_.TotalSysTime(),
                 LinearSystemFactorization_.TotalWallclockTime());
    jnlst.Printf(level, category,
                 " LinearSystemBackSolve..............: %10.3f (sys: %10.3f wall: %10.3f)\n",
                 LinearSystemBackSolve_.TotalCpuTime(),
                 LinearSystemBackSolve_.TotalSysTime(),
                 LinearSystemBackSolve_.TotalWallclockTime());
    jnlst.Printf(level, category,
                 " LinearSystemStructureConverter.....: %10.3f (sys: %10.3f wall: %10.3f)\n",
                 LinearSystemStructureConverter_.TotalCpuTime(),
                 LinearSystemStructureConverter_.TotalSysTime(),
                 LinearSystemStructureConverter_.TotalWallclockTime());
    jnlst.Printf(level, category,
                 "  LinearSystemStructureConverterInit: %10.3f (sys: %10.3f wall: %10.3f)\n",
                 LinearSystemStructureConverterInit_.TotalCpuTime(),
                 LinearSystemStructureConverterInit_.TotalSysTime(),
                 LinearSystemStructureConverterInit_.TotalWallclockTime());
    jnlst.Printf(level, category,
                 "QualityFunctionSearch...............: %10.3f (sys: %10.3f wall: %10.3f)\n",
                 QualityFunctionSearch_.TotalCpuTime(),
                 QualityFunctionSearch_.TotalSysTime(),
                 QualityFunctionSearch_.TotalWallclockTime());
    jnlst.Printf(level, category,
                 "TryCorrector........................: %10.3f (sys: %10.3f wall: %10.3f)\n",
                 TryCorrector_.TotalCpuTime(),
                 TryCorrector_.TotalSysTime(),
                 TryCorrector_.TotalWallclockTime());
    jnlst.Printf(level, category,
                 "Task1...............................: %10.3f (sys: %10.3f wall: %10.3f)\n",
                 Task1_.TotalCpuTime(),
                 Task1_.TotalSysTime(),
                 Task1_.TotalWallclockTime());
    jnlst.Printf(level, category,
                 "Task2...............................: %10.3f (sys: %10.3f wall: %10.3f)\n",
                 Task2_.TotalCpuTime(),
                 Task2_.TotalSysTime(),
                 Task2_.TotalWallclockTime());
    jnlst.Printf(level, category,
                 "Task3...............................: %10.3f (sys: %10.3f wall: %10.3f)\n",
                 Task3_.TotalCpuTime(),
                 Task3_.TotalSysTime(),
                 Task3_.TotalWallclockTime());
    jnlst.Printf(level, category,
                 "Task4...............................: %10.3f (sys: %10.3f wall: %10.3f)\n",
                 Task4_.TotalCpuTime(),
                 Task4_.TotalSysTime(),
                 Task4_.TotalWallclockTime());
    jnlst.Printf(level, category,
                 "Task5...............................: %10.3f (sys: %10.3f wall: %10.3f)\n",
                 Task5_.TotalCpuTime(),
                 Task5_.TotalSysTime(),
                 Task5_.TotalWallclockTime());
  }
} // namespace Ipopt

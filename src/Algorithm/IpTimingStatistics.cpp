// Copyright (C) 2005, 2009 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Authors:  Andreas Waechter          IBM    2005-09-19

#include "IpTimingStatistics.hpp"

namespace Ipopt
{

bool TimingStatistics::IsFunctionEvaluationTimeEnabled() const
{
   return f_eval_time_.IsEnabled() || grad_f_eval_time_.IsEnabled() || c_eval_time_.IsEnabled()
          || d_eval_time_.IsEnabled() || jac_c_eval_time_.IsEnabled() || jac_d_eval_time_.IsEnabled()
          || h_eval_time_.IsEnabled();
}

Number TimingStatistics::TotalFunctionEvaluationCpuTime() const
{
   return f_eval_time_.TotalCpuTime() + grad_f_eval_time_.TotalCpuTime() + c_eval_time_.TotalCpuTime()
          + d_eval_time_.TotalCpuTime() + jac_c_eval_time_.TotalCpuTime() + jac_d_eval_time_.TotalCpuTime()
          + h_eval_time_.TotalCpuTime();
}

Number TimingStatistics::TotalFunctionEvaluationSysTime() const
{
   return f_eval_time_.TotalSysTime() + grad_f_eval_time_.TotalSysTime() + c_eval_time_.TotalSysTime()
          + d_eval_time_.TotalSysTime() + jac_c_eval_time_.TotalSysTime() + jac_d_eval_time_.TotalSysTime()
          + h_eval_time_.TotalSysTime();
}

Number TimingStatistics::TotalFunctionEvaluationWallclockTime() const
{
   return f_eval_time_.TotalWallclockTime() + grad_f_eval_time_.TotalWallclockTime() + c_eval_time_.TotalWallclockTime()
          + d_eval_time_.TotalWallclockTime() + jac_c_eval_time_.TotalWallclockTime()
          + jac_d_eval_time_.TotalWallclockTime() + h_eval_time_.TotalWallclockTime();
}

/** Method for enabling all timed tasked. */
void TimingStatistics::EnableTimes()
{
   OverallAlgorithm_.IsEnabled();
   PrintProblemStatistics_.Enable();
   InitializeIterates_.Enable();
   UpdateHessian_.Enable();
   OutputIteration_.Enable();
   UpdateBarrierParameter_.Enable();
   ComputeSearchDirection_.Enable();
   ComputeAcceptableTrialPoint_.Enable();
   AcceptTrialPoint_.Enable();
   CheckConvergence_.Enable();
   PDSystemSolverTotal_.Enable();
   PDSystemSolverSolveOnce_.Enable();
   ComputeResiduals_.Enable();
   StdAugSystemSolverMultiSolve_.Enable();
   LinearSystemScaling_.Enable();
   LinearSystemSymbolicFactorization_.Enable();
   LinearSystemFactorization_.Enable();
   LinearSystemBackSolve_.Enable();
   LinearSystemStructureConverter_.Enable();
   LinearSystemStructureConverterInit_.Enable();
   QualityFunctionSearch_.Enable();
   TryCorrector_.Enable();
   Task1_.Enable();
   Task2_.Enable();
   Task3_.Enable();
   Task4_.Enable();
   Task5_.Enable();
   Task6_.Enable();
   f_eval_time_.Enable();
   grad_f_eval_time_.Enable();
   c_eval_time_.Enable();
   d_eval_time_.Enable();
   jac_c_eval_time_.Enable();
   jac_d_eval_time_.Enable();
   h_eval_time_.Enable();
}

/** Method for disabling all timed tasks except for OverallAlgorithm */
void TimingStatistics::DisableTimes()
{
   PrintProblemStatistics_.Disable();
   InitializeIterates_.Disable();
   UpdateHessian_.Disable();
   OutputIteration_.Disable();
   UpdateBarrierParameter_.Disable();
   ComputeSearchDirection_.Disable();
   ComputeAcceptableTrialPoint_.Disable();
   AcceptTrialPoint_.Disable();
   CheckConvergence_.Disable();
   PDSystemSolverTotal_.Disable();
   PDSystemSolverSolveOnce_.Disable();
   ComputeResiduals_.Disable();
   StdAugSystemSolverMultiSolve_.Disable();
   LinearSystemScaling_.Disable();
   LinearSystemSymbolicFactorization_.Disable();
   LinearSystemFactorization_.Disable();
   LinearSystemBackSolve_.Disable();
   LinearSystemStructureConverter_.Disable();
   LinearSystemStructureConverterInit_.Disable();
   QualityFunctionSearch_.Disable();
   TryCorrector_.Disable();
   Task1_.Disable();
   Task2_.Disable();
   Task3_.Disable();
   Task4_.Disable();
   Task5_.Disable();
   Task6_.Disable();
   f_eval_time_.Disable();
   grad_f_eval_time_.Disable();
   c_eval_time_.Disable();
   d_eval_time_.Disable();
   jac_c_eval_time_.Disable();
   jac_d_eval_time_.Disable();
   h_eval_time_.Disable();
}

void TimingStatistics::ResetTimes()
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
   f_eval_time_.Reset();
   grad_f_eval_time_.Reset();
   c_eval_time_.Reset();
   d_eval_time_.Reset();
   jac_c_eval_time_.Reset();
   jac_d_eval_time_.Reset();
   h_eval_time_.Reset();
}

void TimingStatistics::PrintAllTimingStatistics(
   const Journalist& jnlst,
   EJournalLevel     level,
   EJournalCategory  category
) const
{
   if( !jnlst.ProduceOutput(level, category) )
   {
      return;
   }

   if( OverallAlgorithm_.IsEnabled() )
      jnlst.Printf(level, category,
                   "OverallAlgorithm....................: %10.3f (sys: %10.3f wall: %10.3f)\n", OverallAlgorithm_.TotalCpuTime(), OverallAlgorithm_.TotalSysTime(), OverallAlgorithm_.TotalWallclockTime());
   else
   {
      jnlst.Printf(level, category, "OverallAlgorithm\n");
   }
   if( PrintProblemStatistics_.IsEnabled() )
      jnlst.Printf(level, category,
                   " PrintProblemStatistics.............: %10.3f (sys: %10.3f wall: %10.3f)\n", PrintProblemStatistics_.TotalCpuTime(), PrintProblemStatistics_.TotalSysTime(), PrintProblemStatistics_.TotalWallclockTime());
   if( InitializeIterates_.IsEnabled() )
      jnlst.Printf(level, category,
                   " InitializeIterates.................: %10.3f (sys: %10.3f wall: %10.3f)\n", InitializeIterates_.TotalCpuTime(), InitializeIterates_.TotalSysTime(), InitializeIterates_.TotalWallclockTime());
   if( UpdateHessian_.IsEnabled() )
      jnlst.Printf(level, category,
                   " UpdateHessian......................: %10.3f (sys: %10.3f wall: %10.3f)\n", UpdateHessian_.TotalCpuTime(), UpdateHessian_.TotalSysTime(), UpdateHessian_.TotalWallclockTime());
   if( OutputIteration_.IsEnabled() )
      jnlst.Printf(level, category,
                   " OutputIteration....................: %10.3f (sys: %10.3f wall: %10.3f)\n", OutputIteration_.TotalCpuTime(), OutputIteration_.TotalSysTime(), OutputIteration_.TotalWallclockTime());
   if( UpdateBarrierParameter_.IsEnabled() )
      jnlst.Printf(level, category,
                   " UpdateBarrierParameter.............: %10.3f (sys: %10.3f wall: %10.3f)\n", UpdateBarrierParameter_.TotalCpuTime(), UpdateBarrierParameter_.TotalSysTime(), UpdateBarrierParameter_.TotalWallclockTime());
   if( ComputeSearchDirection_.IsEnabled() )
      jnlst.Printf(level, category,
                   " ComputeSearchDirection.............: %10.3f (sys: %10.3f wall: %10.3f)\n", ComputeSearchDirection_.TotalCpuTime(), ComputeSearchDirection_.TotalSysTime(), ComputeSearchDirection_.TotalWallclockTime());
   if( ComputeAcceptableTrialPoint_.IsEnabled() )
      jnlst.Printf(level, category,
                   " ComputeAcceptableTrialPoint........: %10.3f (sys: %10.3f wall: %10.3f)\n", ComputeAcceptableTrialPoint_.TotalCpuTime(), ComputeAcceptableTrialPoint_.TotalSysTime(), ComputeAcceptableTrialPoint_.TotalWallclockTime());
   if( AcceptTrialPoint_.IsEnabled() )
      jnlst.Printf(level, category,
                   " AcceptTrialPoint...................: %10.3f (sys: %10.3f wall: %10.3f)\n", AcceptTrialPoint_.TotalCpuTime(), AcceptTrialPoint_.TotalSysTime(), AcceptTrialPoint_.TotalWallclockTime());
   if( CheckConvergence_.IsEnabled() )
      jnlst.Printf(level, category,
                   " CheckConvergence...................: %10.3f (sys: %10.3f wall: %10.3f)\n", CheckConvergence_.TotalCpuTime(), CheckConvergence_.TotalSysTime(), CheckConvergence_.TotalWallclockTime());

   if( PDSystemSolverTotal_.IsEnabled() )
      jnlst.Printf(level, category,
                   "PDSystemSolverTotal.................: %10.3f (sys: %10.3f wall: %10.3f)\n", PDSystemSolverTotal_.TotalCpuTime(), PDSystemSolverTotal_.TotalSysTime(), PDSystemSolverTotal_.TotalWallclockTime());
   else if( PDSystemSolverSolveOnce_.IsEnabled() || ComputeResiduals_.IsEnabled() || StdAugSystemSolverMultiSolve_.IsEnabled() || LinearSystemScaling_.IsEnabled() || LinearSystemSymbolicFactorization_.IsEnabled() || LinearSystemFactorization_.IsEnabled() || LinearSystemBackSolve_.IsEnabled() || LinearSystemStructureConverter_.IsEnabled() || LinearSystemStructureConverterInit_.IsEnabled() )
   {
      jnlst.Printf(level, category, "PDSystemSolverTotal\n");
   }
   if( PDSystemSolverSolveOnce_.IsEnabled() )
      jnlst.Printf(level, category,
                   " PDSystemSolverSolveOnce............: %10.3f (sys: %10.3f wall: %10.3f)\n", PDSystemSolverSolveOnce_.TotalCpuTime(), PDSystemSolverSolveOnce_.TotalSysTime(), PDSystemSolverSolveOnce_.TotalWallclockTime());
   if( ComputeResiduals_.IsEnabled() )
      jnlst.Printf(level, category,
                   " ComputeResiduals...................: %10.3f (sys: %10.3f wall: %10.3f)\n", ComputeResiduals_.TotalCpuTime(), ComputeResiduals_.TotalSysTime(), ComputeResiduals_.TotalWallclockTime());
   if( StdAugSystemSolverMultiSolve_.IsEnabled() )
      jnlst.Printf(level, category,
                   " StdAugSystemSolverMultiSolve.......: %10.3f (sys: %10.3f wall: %10.3f)\n", StdAugSystemSolverMultiSolve_.TotalCpuTime(), StdAugSystemSolverMultiSolve_.TotalSysTime(), StdAugSystemSolverMultiSolve_.TotalWallclockTime());
   if( LinearSystemScaling_.IsEnabled() )
      jnlst.Printf(level, category,
                   " LinearSystemScaling................: %10.3f (sys: %10.3f wall: %10.3f)\n", LinearSystemScaling_.TotalCpuTime(), LinearSystemScaling_.TotalSysTime(), LinearSystemScaling_.TotalWallclockTime());
   if( LinearSystemSymbolicFactorization_.IsEnabled() )
      jnlst.Printf(level, category,
                   " LinearSystemSymbolicFactorization..: %10.3f (sys: %10.3f wall: %10.3f)\n", LinearSystemSymbolicFactorization_.TotalCpuTime(), LinearSystemSymbolicFactorization_.TotalSysTime(), LinearSystemSymbolicFactorization_.TotalWallclockTime());
   if( LinearSystemFactorization_.IsEnabled() )
      jnlst.Printf(level, category,
                   " LinearSystemFactorization..........: %10.3f (sys: %10.3f wall: %10.3f)\n", LinearSystemFactorization_.TotalCpuTime(), LinearSystemFactorization_.TotalSysTime(), LinearSystemFactorization_.TotalWallclockTime());
   if( LinearSystemBackSolve_.IsEnabled() )
      jnlst.Printf(level, category,
                   " LinearSystemBackSolve..............: %10.3f (sys: %10.3f wall: %10.3f)\n", LinearSystemBackSolve_.TotalCpuTime(), LinearSystemBackSolve_.TotalSysTime(), LinearSystemBackSolve_.TotalWallclockTime());
   if( LinearSystemStructureConverter_.IsEnabled() )
      jnlst.Printf(level, category,
                   " LinearSystemStructureConverter.....: %10.3f (sys: %10.3f wall: %10.3f)\n", LinearSystemStructureConverter_.TotalCpuTime(), LinearSystemStructureConverter_.TotalSysTime(), LinearSystemStructureConverter_.TotalWallclockTime());
   if( LinearSystemStructureConverterInit_.IsEnabled() )
      jnlst.Printf(level, category,
                   "  LinearSystemStructureConverterInit: %10.3f (sys: %10.3f wall: %10.3f)\n", LinearSystemStructureConverterInit_.TotalCpuTime(), LinearSystemStructureConverterInit_.TotalSysTime(), LinearSystemStructureConverterInit_.TotalWallclockTime());

   if( QualityFunctionSearch_.IsEnabled() )
      jnlst.Printf(level, category,
                   "QualityFunctionSearch...............: %10.3f (sys: %10.3f wall: %10.3f)\n", QualityFunctionSearch_.TotalCpuTime(), QualityFunctionSearch_.TotalSysTime(), QualityFunctionSearch_.TotalWallclockTime());
   if( TryCorrector_.IsEnabled() )
      jnlst.Printf(level, category,
                   "TryCorrector........................: %10.3f (sys: %10.3f wall: %10.3f)\n", TryCorrector_.TotalCpuTime(), TryCorrector_.TotalSysTime(), TryCorrector_.TotalWallclockTime());
   if( Task1_.IsEnabled() )
      jnlst.Printf(level, category,
                   "Task1...............................: %10.3f (sys: %10.3f wall: %10.3f)\n", Task1_.TotalCpuTime(), Task1_.TotalSysTime(), Task1_.TotalWallclockTime());
   if( Task2_.IsEnabled() )
      jnlst.Printf(level, category,
                   "Task2...............................: %10.3f (sys: %10.3f wall: %10.3f)\n", Task2_.TotalCpuTime(), Task2_.TotalSysTime(), Task2_.TotalWallclockTime());
   if( Task3_.IsEnabled() )
      jnlst.Printf(level, category,
                   "Task3...............................: %10.3f (sys: %10.3f wall: %10.3f)\n", Task3_.TotalCpuTime(), Task3_.TotalSysTime(), Task3_.TotalWallclockTime());
   if( Task4_.IsEnabled() )
      jnlst.Printf(level, category,
                   "Task4...............................: %10.3f (sys: %10.3f wall: %10.3f)\n", Task4_.TotalCpuTime(), Task4_.TotalSysTime(), Task4_.TotalWallclockTime());
   if( Task5_.IsEnabled() )
      jnlst.Printf(level, category,
                   "Task5...............................: %10.3f (sys: %10.3f wall: %10.3f)\n", Task5_.TotalCpuTime(), Task5_.TotalSysTime(), Task5_.TotalWallclockTime());
   if( Task6_.IsEnabled() )
      jnlst.Printf(level, category,
                   "Task6...............................: %10.3f (sys: %10.3f wall: %10.3f)\n", Task6_.TotalCpuTime(), Task6_.TotalSysTime(), Task6_.TotalWallclockTime());

   if( IsFunctionEvaluationTimeEnabled() )
      jnlst.Printf(level, category,
                   "Function Evaluations................: %10.3f (sys: %10.3f wall: %10.3f)\n", TotalFunctionEvaluationCpuTime(), TotalFunctionEvaluationSysTime(), TotalFunctionEvaluationWallclockTime());
   if( f_eval_time_.IsEnabled() )
      jnlst.Printf(level, category,
                   " Objective function.................: %10.3f (sys: %10.3f wall: %10.3f)\n", f_eval_time_.TotalCpuTime(), f_eval_time_.TotalSysTime(), f_eval_time_.TotalWallclockTime());
   if( grad_f_eval_time_.IsEnabled() )
      jnlst.Printf(level, category,
                   " Objective function gradient........: %10.3f (sys: %10.3f wall: %10.3f)\n", grad_f_eval_time_.TotalCpuTime(), grad_f_eval_time_.TotalSysTime(), grad_f_eval_time_.TotalWallclockTime());
   if( c_eval_time_.IsEnabled() )
      jnlst.Printf(level, category,
                   " Equality constraints...............: %10.3f (sys: %10.3f wall: %10.3f)\n", c_eval_time_.TotalCpuTime(), c_eval_time_.TotalSysTime(), c_eval_time_.TotalWallclockTime());
   if( d_eval_time_.IsEnabled() )
      jnlst.Printf(level, category,
                   " Inequality constraints.............: %10.3f (sys: %10.3f wall: %10.3f)\n", d_eval_time_.TotalCpuTime(), d_eval_time_.TotalSysTime(), d_eval_time_.TotalWallclockTime());
   if( jac_c_eval_time_.IsEnabled() )
      jnlst.Printf(level, category,
                   " Equality constraint Jacobian.......: %10.3f (sys: %10.3f wall: %10.3f)\n", jac_c_eval_time_.TotalCpuTime(), jac_c_eval_time_.TotalSysTime(), jac_c_eval_time_.TotalWallclockTime());
   if( jac_d_eval_time_.IsEnabled() )
      jnlst.Printf(level, category,
                   " Inequality constraint Jacobian.....: %10.3f (sys: %10.3f wall: %10.3f)\n", jac_d_eval_time_.TotalCpuTime(), jac_d_eval_time_.TotalSysTime(), jac_d_eval_time_.TotalWallclockTime());
   if( h_eval_time_.IsEnabled() )
      jnlst.Printf(level, category,
                   " Lagrangian Hessian.................: %10.3f (sys: %10.3f wall: %10.3f)\n", h_eval_time_.TotalCpuTime(), h_eval_time_.TotalSysTime(), h_eval_time_.TotalWallclockTime());
}

} // namespace Ipopt

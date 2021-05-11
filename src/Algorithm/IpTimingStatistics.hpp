// Copyright (C) 2005, 2008 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Authors:  Andreas Waechter               IBM    2005-09-19

#ifndef __IPTIMINGSTATISTICS_HPP__
#define __IPTIMINGSTATISTICS_HPP__

#include "IpReferenced.hpp"
#include "IpJournalist.hpp"
#include "IpTimedTask.hpp"

namespace Ipopt
{
/** This class collects all timing statistics for Ipopt.
 */
class IPOPTLIB_EXPORT TimingStatistics: public ReferencedObject
{
public:
   /**@name Constructors/Destructors */
   ///@{
   /** Default constructor. */
   TimingStatistics()
   { }

   /** Destructor */
   virtual ~TimingStatistics()
   { }
   ///@}

   /// Whether timing of function evaluation has been enabled
   /// @since 3.14.0
   bool IsFunctionEvaluationTimeEnabled() const;
   /// total CPU time spend in function evaluation
   /// @since 3.14.0
   Number TotalFunctionEvaluationCpuTime() const;
   /// total system time spend in function evaluation
   /// @since 3.14.0
   Number TotalFunctionEvaluationSysTime() const;
   /// total wall-clock time spend in function evaluation
   /// @since 3.14.0
   Number TotalFunctionEvaluationWallclockTime() const;

   /** Method for resetting all times. */
   void ResetTimes();

   /** Method for enabling all timed tasked.
    * @since 3.14.0
    */
   void EnableTimes();

   /** Method for disabling all timed tasks except for OverallAlgorithm
    * @since 3.14.0
    */
   void DisableTimes();

   /** Method for printing all timing information */
   void PrintAllTimingStatistics(
      const Journalist& jnlst,
      EJournalLevel     level,
      EJournalCategory  category
   ) const;

   /**@name Accessor methods to all timed tasks. */
   ///@{
   TimedTask& OverallAlgorithm()
   {
      return OverallAlgorithm_;
   }

   /// @since 3.14.0
   const TimedTask& OverallAlgorithm() const
   {
      return OverallAlgorithm_;
   }

   TimedTask& PrintProblemStatistics()
   {
      return PrintProblemStatistics_;
   }

   TimedTask& InitializeIterates()
   {
      return InitializeIterates_;
   }

   TimedTask& UpdateHessian()
   {
      return UpdateHessian_;
   }

   TimedTask& OutputIteration()
   {
      return OutputIteration_;
   }

   TimedTask& UpdateBarrierParameter()
   {
      return UpdateBarrierParameter_;
   }

   TimedTask& ComputeSearchDirection()
   {
      return ComputeSearchDirection_;
   }

   TimedTask& ComputeAcceptableTrialPoint()
   {
      return ComputeAcceptableTrialPoint_;
   }

   TimedTask& AcceptTrialPoint()
   {
      return AcceptTrialPoint_;
   }

   TimedTask& CheckConvergence()
   {
      return CheckConvergence_;
   }

   TimedTask& PDSystemSolverTotal()
   {
      return PDSystemSolverTotal_;
   }

   TimedTask& PDSystemSolverSolveOnce()
   {
      return PDSystemSolverSolveOnce_;
   }

   TimedTask& ComputeResiduals()
   {
      return ComputeResiduals_;
   }

   TimedTask& StdAugSystemSolverMultiSolve()
   {
      return StdAugSystemSolverMultiSolve_;
   }

   TimedTask& LinearSystemScaling()
   {
      return LinearSystemScaling_;
   }

   TimedTask& LinearSystemSymbolicFactorization()
   {
      return LinearSystemSymbolicFactorization_;
   }

   TimedTask& LinearSystemFactorization()
   {
      return LinearSystemFactorization_;
   }

   TimedTask& LinearSystemBackSolve()
   {
      return LinearSystemBackSolve_;
   }

   TimedTask& LinearSystemStructureConverter()
   {
      return LinearSystemStructureConverter_;
   }

   TimedTask& LinearSystemStructureConverterInit()
   {
      return LinearSystemStructureConverterInit_;
   }

   TimedTask& QualityFunctionSearch()
   {
      return QualityFunctionSearch_;
   }

   TimedTask& TryCorrector()
   {
      return TryCorrector_;
   }

   TimedTask& Task1()
   {
      return Task1_;
   }

   TimedTask& Task2()
   {
      return Task2_;
   }

   TimedTask& Task3()
   {
      return Task3_;
   }

   TimedTask& Task4()
   {
      return Task4_;
   }

   TimedTask& Task5()
   {
      return Task5_;
   }

   TimedTask& Task6()
   {
      return Task6_;
   }

   /// @since 3.14.0
   TimedTask& f_eval_time()
   {
      return f_eval_time_;
   }

   /// @since 3.14.0
   TimedTask& grad_f_eval_time()
   {
      return grad_f_eval_time_;
   }

   /// @since 3.14.0
   TimedTask& c_eval_time()
   {
      return c_eval_time_;
   }

   /// @since 3.14.0
   TimedTask& jac_c_eval_time()
   {
      return jac_c_eval_time_;
   }

   /// @since 3.14.0
   TimedTask& d_eval_time()
   {
      return d_eval_time_;
   }

   /// @since 3.14.0
   TimedTask& jac_d_eval_time()
   {
      return jac_d_eval_time_;
   }

   /// @since 3.14.0
   TimedTask& h_eval_time()
   {
      return h_eval_time_;
   }
   ///@}

private:
   /**@name Default Compiler Generated Methods
    * (Hidden to avoid implicit creation/calling).
    *
    * These methods are not implemented and
    * we do not want the compiler to implement
    * them for us, so we declare them private
    * and do not define them. This ensures that
    * they will not be implicitly created/called.
    */
   ///@{
   /** Copy Constructor */
   TimingStatistics(
      const TimingStatistics&
   );

   /** Default Assignment Operator */
   void operator=(
      const TimingStatistics&
   );
   ///@}

   /**@name All timed tasks. */
   ///@{
   TimedTask OverallAlgorithm_;
   TimedTask PrintProblemStatistics_;
   TimedTask InitializeIterates_;
   TimedTask UpdateHessian_;
   TimedTask OutputIteration_;
   TimedTask UpdateBarrierParameter_;
   TimedTask ComputeSearchDirection_;
   TimedTask ComputeAcceptableTrialPoint_;
   TimedTask AcceptTrialPoint_;
   TimedTask CheckConvergence_;

   TimedTask PDSystemSolverTotal_;
   TimedTask PDSystemSolverSolveOnce_;
   TimedTask ComputeResiduals_;
   TimedTask StdAugSystemSolverMultiSolve_;
   TimedTask LinearSystemScaling_;
   TimedTask LinearSystemSymbolicFactorization_;
   TimedTask LinearSystemFactorization_;
   TimedTask LinearSystemBackSolve_;
   TimedTask LinearSystemStructureConverter_;
   TimedTask LinearSystemStructureConverterInit_;
   TimedTask QualityFunctionSearch_;
   TimedTask TryCorrector_;

   TimedTask Task1_;
   TimedTask Task2_;
   TimedTask Task3_;
   TimedTask Task4_;
   TimedTask Task5_;
   TimedTask Task6_;

   TimedTask f_eval_time_;
   TimedTask grad_f_eval_time_;
   TimedTask c_eval_time_;
   TimedTask jac_c_eval_time_;
   TimedTask d_eval_time_;
   TimedTask jac_d_eval_time_;
   TimedTask h_eval_time_;
   ///@}
};

} // namespace Ipopt

#endif

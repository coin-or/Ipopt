// Copyright (C) 2005, International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id: IpSolveStatistics.hpp 507 2005-08-26 18:33:43Z andreasw $
//
// Authors:  Andreas Waechter               IBM    2005-09-19

#ifndef __IPTIMINGSTATISTICS_HPP__
#define __IPTIMINGSTATISTICS_HPP__

#include "IpReferenced.hpp"
#include "IpJournalist.hpp"

#ifdef HAVE_CTIME
# include <ctime>
#else
# ifdef HAVE_TIME_H
#  include <time.h>
# else
#  error "don't have header file for time"
# endif
#endif

// The following lines are copied from CoinTime.hpp
// We should probably make some more tests here
#if defined(_MSC_VER)
// Turn off compiler warning about long names
#  pragma warning(disable:4786)
#else
// MacOS-X and FreeBSD needs sys/time.h
#if defined(__MACH__) || defined (__FreeBSD__)
#include <sys/time.h>
#endif
#include <sys/resource.h>
#endif

namespace Ipopt
{
  /** This class is used to collect timing information for a
   *  particular task. */
  class TimedTask
  {
  public:
    /**@name Constructors/Destructors */
    //@{
    /** Default constructor. */
    TimedTask()
        :
        total_time_(0.),
        start_called_(false),
        end_called_(true)
    {}

    /** Default destructor */
    ~TimedTask()
    {}
    //@}

    /** Method that is called before execution of the task. */
    void Start()
    {
      DBG_ASSERT(end_called_);
      DBG_ASSERT(!start_called_);
      end_called_ = false;
      start_called_ = true;
      start_time_ = CpuTime();
    }

    /** Method that is called after execution of the task. */
    void End()
    {
      DBG_ASSERT(!end_called_);
      DBG_ASSERT(start_called_);
      end_called_ = true;
      start_called_ = false;
      total_time_ += CpuTime() - start_time_;
    }

    /** Method that is called after execution of the task for which
     *  timing might have been started.  This only updates the timing
     *  if the timing has indeed been conducted. This is useful to
     *  stop timing after catching exceptions. */
    void EndIfStarted()
    {
      if (start_called_) {
        end_called_ = true;
        start_called_ = false;
        total_time_ += CpuTime() - start_time_;
      }
      DBG_ASSERT(end_called_);
    }

    /** Method returning total time spend for task so far. */
    Number TotalTime() const
    {
      DBG_ASSERT(end_called_);
      return total_time_;
    }

  private:
    /**@name Default Compiler Generated Methods (Hidden to avoid
     * implicit creation/calling).  These methods are not
     * implemented and we do not want the compiler to implement them
     * for us, so we declare them private and do not define
     * them. This ensures that they will not be implicitly
     * created/called. */
    //@{
    /** Copy Constructor */
    TimedTask(const TimedTask&);

    /** Overloaded Equals Operator */
    void operator=(const TimedTask&);
    //@}

    /** Time at beginning of task. */
    Number start_time_;
    /** Total time for task measured so far. */
    Number total_time_;

    /** @name fields for debugging */
    //@{
    bool start_called_;
    bool end_called_;
    //@}

    // The following lines were taken from CoinTime.hpp in COIN/Coin
    /** method determining CPU executed since start of program */
    static inline Number CpuTime()
    {
      double cpu_temp;
#if defined(_MSC_VER)

      unsigned int ticksnow;        /* clock_t is same as int */

      ticksnow = (unsigned int)clock();

      cpu_temp = (double)((double)ticksnow/CLOCKS_PER_SEC);
#else

      struct rusage usage;
      getrusage(RUSAGE_SELF,&usage);
      cpu_temp = usage.ru_utime.tv_sec;
      cpu_temp += 1.0e-6*((double) usage.ru_utime.tv_usec);
#endif

      return cpu_temp;
    }
  };

  /** This class collects all timing statistics for Ipopt.
   */
  class TimingStatistics : public ReferencedObject
  {
  public:
    /**@name Constructors/Destructors */
    //@{
    /** Default constructor. */
    TimingStatistics()
    {}

    /** Default destructor */
    virtual ~TimingStatistics()
    {}
    //@}

    /** Method for printing all timing information */
    void PrintAllTimingStatistics(Journalist& jnlst,
                                  EJournalLevel level,
                                  EJournalCategory category) const;

    /**@name All timed tasks. */
    //@{
    TimedTask OverallAlgorithm;
    TimedTask PrintProblemStatistics;
    TimedTask InitializeIterates;
    TimedTask ActualizeHessian;
    TimedTask OutputIteration;
    TimedTask UpdateBarrierParameter;
    TimedTask ComputeSearchDirection;
    TimedTask ComputeAcceptableTrialPoint;
    TimedTask AcceptTrialPoint;
    TimedTask CheckConvergence;

    TimedTask PDSystemSolverTotal;
    TimedTask PDSystemSolverSolveOnce;
    TimedTask LinearSystemScaling;
    TimedTask LinearSystemSymbolicFactorization;
    TimedTask LinearSystemFactorization;
    TimedTask LinearSystemBackSolve;
    TimedTask QualityFunctionSearch;
    TimedTask TryCorrector;

    TimedTask Task1;
    TimedTask Task2;
    TimedTask Task3;
    TimedTask Task4;
    TimedTask Task5;
    TimedTask Task6;
    //@}

  private:
    /**@name Default Compiler Generated Methods
     * (Hidden to avoid implicit creation/calling).
     * These methods are not implemented and 
     * we do not want the compiler to implement
     * them for us, so we declare them private
     * and do not define them. This ensures that
     * they will not be implicitly created/called. */
    //@{
    /** Copy Constructor */
    TimingStatistics(const TimingStatistics&);

    /** Overloaded Equals Operator */
    void operator=(const TimingStatistics&);
    //@}

  };

} // namespace Ipopt

#endif

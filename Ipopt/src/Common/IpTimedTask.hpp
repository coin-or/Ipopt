// Copyright (C) 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id$
//
// Authors:  Andreas Waechter               IBM    2005-09-19

#ifndef __IPTIMEDTASK_HPP__
#define __IPTIMEDTASK_HPP__

#include "IpUtils.hpp"

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

    /** Method for resetting time to zero. */
    void Reset()
    {
      total_time_ = 0.;
      start_called_ = false;
      end_called_ = true;
    }

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

  };
} // namespace Ipopt

#endif

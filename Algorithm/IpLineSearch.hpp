// Copyright (C) 2004, International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id$
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#ifndef __IPLINESEARCH_HPP__
#define __IPLINESEARCH_HPP__

#include "IpUtils.hpp"
#include "IpAlgStrategy.hpp"
#include "IpIpoptCalculatedQuantities.hpp"

namespace Ipopt
{

  /** Base class for line search objects.
   */
  class LineSearch : public AlgorithmStrategyObject
  {
  public:
    /**@name Constructors/Destructors */
    //@{
    /** Default Constructor */
    LineSearch()
    {}
    ;

    /** Default destructor */
    virtual ~LineSearch()
    {}
    ;
    //@}

    /** Perform the line search.  As search direction the delta
     *  in the data object is used
     */
    virtual void FindAcceptableTrialPoint() = 0;

    /** Resest the line search.
     *  This function should be called if all previous information
     *  should be discarded when the line search is performed the
     *  next time.  For example, this method should be calledm if
     *  the barrier parameter is changed.
     */
    virtual void Reset() = 0;


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
    LineSearch(const LineSearch&);

    /** Overloaded Equals Operator */
    void operator=(const LineSearch&);
    //@}

  };

} // namespace Ipopt

#endif

// Copyright (C) 2004,2005 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id$
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-09-29

#ifndef __IPALGBUILDER_HPP__
#define __IPALGBUILDER_HPP__

#include "IpIpoptAlg.hpp"

namespace Ipopt
{

  /** Builder to create a complete IpoptAlg object.  This object
   *  contains all subelements (such as line search objects etc).  How
   *  the resulting IpoptAlg object is build can be influenced by the
   *  options. */
  class AlgorithmBuilder
  {
  public:
    /**@name Static methods for building Algorithm objects.
    *		Most users will want to use BuildAlgorithmFromTNLP. 
    *     The other methods are for more specialized algorithms
    */
    //@{
    //	static SmartPtr<IpoptAlgorithm> BuildAlgorithmFromTNLP(...)

    static SmartPtr<IpoptAlgorithm> BuildBasicAlgorithm(const Journalist& jnlst,
        const OptionsList& options,
        const std::string& prefix);
    //@}

    /** Methods for IpoptTypeInfo */
    //@{
    /** register the options used by the algorithm builder */
    static void RegisterOptions(SmartPtr<RegisteredOptions> roptions);
    //@}
  };
} // namespace Ipopt

#endif

// Copyright (C) 2004, 2007 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// $Id$
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-09-29

#ifndef __IPALGBUILDER_HPP__
#define __IPALGBUILDER_HPP__

#include "IpIpoptAlg.hpp"
#include "IpReferenced.hpp"
#include "IpAugSystemSolver.hpp"

namespace Ipopt
{

  /** Builder to create a complete IpoptAlg object.  This object
   *  contains all subelements (such as line search objects etc).  How
   *  the resulting IpoptAlg object is built can be influenced by the
   *  options.
   *
   *  The optional argument custom_solver allows the expert user to
   *  provide a specialized linear solver (e.g., of the type
   *  GenAugSystemSolver), possibly for selfmade matrix objects.
   *
   * TODO: Currently, this is a basic implementation with everything
   *  in one method that can be overloaded. This will need to be expanded
   *  to allow customization of different parts without recoding everything. 
   */
  class AlgorithmBuilder : public ReferencedObject
  {
  public:
    /**@name Constructors/Destructors */
    //@{
    /** Constructor */
    AlgorithmBuilder(SmartPtr<AugSystemSolver> custom_solver=NULL);

    /** Destructor */
    virtual ~AlgorithmBuilder()
    {}

    //@}

    /** @name Methods to build parts of the algorithm */
    //@{
    virtual void BuildIpoptObjects(const Journalist& jnlst,
                                   const OptionsList& options,
                                   const std::string& prefix,
                                   const SmartPtr<NLP>& nlp,
                                   SmartPtr<IpoptNLP>& ip_nlp,
                                   SmartPtr<IpoptData>& ip_data,
                                   SmartPtr<IpoptCalculatedQuantities>& ip_cq);

    virtual SmartPtr<IpoptAlgorithm> BuildBasicAlgorithm(const Journalist& jnlst,
        const OptionsList& options,
        const std::string& prefix);
    //@}

    /** Methods for IpoptTypeInfo */
    //@{
    /** register the options used by the algorithm builder */
    static void RegisterOptions(SmartPtr<RegisteredOptions> roptions);
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
    /** Default Constructor */
    //AlgorithmBuilder();

    /** Copy Constructor */
    AlgorithmBuilder(const AlgorithmBuilder&);

    /** Overloaded Equals Operator */
    void operator=(const AlgorithmBuilder&);
    //@}

    /** Optional pointer to AugSystemSolver.  If this is set in the
     *  contructor, we will use this to solver the linear systems if
     *  the option linear_solver=custerm is chosen. */
    SmartPtr<AugSystemSolver> custom_solver_;

  };
} // namespace Ipopt

#endif

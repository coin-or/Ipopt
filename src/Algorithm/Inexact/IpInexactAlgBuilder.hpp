// Copyright (C) 2008 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Authors:  Andreas Waechter     IBM                  2008-09-05
//            based on IpAlgBuilder.hpp (rev 913)

#ifndef __IPINEXACTALGBUILDER_HPP__
#define __IPINEXACTALGBUILDER_HPP__

#include "IpAlgBuilder.hpp"

namespace Ipopt
{

/** Builder to create a complete IpoptAlg object for the inexact
 *  step computation version.
 *
 * TODO: The AlorithmBuilder base class has been reorganized to
 *       allow for easier customization. This class could be
 *       reimplemented to take advantage of that. In particular, a
 *       substantial amount code for generating the SymLinearSolver
 *       and AugSystemSolver is available for reuse.
 */
class InexactAlgorithmBuilder: public AlgorithmBuilder
{
public:
   /**@name Constructors/Destructors */
   ///@{
   /** Constructor */
   InexactAlgorithmBuilder();

   /** Destructor */
   virtual ~InexactAlgorithmBuilder()
   { }
   ///@}

   /** @name Methods to build parts of the algorithm */
   ///@{
   virtual void BuildIpoptObjects(
      const Journalist&                    jnlst,
      const OptionsList&                   options,
      const std::string&                   prefix,
      const SmartPtr<NLP>&                 nlp,
      SmartPtr<IpoptNLP>&                  ip_nlp,
      SmartPtr<IpoptData>&                 ip_data,
      SmartPtr<IpoptCalculatedQuantities>& ip_cq
   );

   virtual SmartPtr<IpoptAlgorithm> BuildBasicAlgorithm(
      const Journalist&  jnlst,
      const OptionsList& options,
      const std::string& prefix
   );
   ///@}

   /** register the options used by the algorithm builder */
   static void RegisterOptions(
      SmartPtr<RegisteredOptions> roptions
   );

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
   InexactAlgorithmBuilder(
      const InexactAlgorithmBuilder&
   );

   /** Overloaded Assignment Operator */
   void operator=(
      const InexactAlgorithmBuilder&
   );
   ///@}

   /** Optional pointer to AugSystemSolver.
    *
    *  If this is set in the constructor, we will use this to solve
    *  the linear systems if the option linear_solver=custom is chosen.
    */
   SmartPtr<AugSystemSolver> custom_solver_;

};

/** Function for setting options who's default is different for the
 *  inexact algorithm compared to the defaults for the regular Ipopt
 *  algorithm.
 *
 *  The options_list is augmented by the different default values,
 *  but only if the corresponding option has not yet been set.
 */
void AddInexactDefaultOptions(
   OptionsList& options_list
);

} // namespace Ipopt

#endif

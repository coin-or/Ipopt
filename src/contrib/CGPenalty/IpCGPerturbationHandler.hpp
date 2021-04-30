// Copyright (C) 2005, 2008 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Authors:  Carl Laird, Andreas Waechter              IBM    2005-08-04

#ifndef __IPCGPERTURBATIONHANDLER_HPP__
#define __IPCGPERTURBATIONHANDLER_HPP__

#include "IpPDPerturbationHandler.hpp"
#include "IpCGPenaltyCq.hpp"

namespace Ipopt
{

/** Class for handling the perturbation factors delta_x, delta_s,
 *  delta_c, and delta_d in the primal dual system.
 *
 *  This class is
 *  used by the PDFullSpaceSolver to handle the cases where the
 *  primal-dual system is singular or has the wrong inertia.  The
 *  perturbation factors are obtained based on simple heuristics,
 *  taking into account the size of previous perturbations.
 */
class CGPerturbationHandler: public PDPerturbationHandler
{
public:
   /**@name Constructors/Destructors */
   ///@{
   /** Default Constructor */
   CGPerturbationHandler();

   /** Destructor */
   virtual ~CGPerturbationHandler()
   { }
   ///@}

   /* overloaded from AlgorithmStrategyObject */
   virtual bool InitializeImpl(
      const OptionsList& options,
      const std::string& prefix
   );

   /** This method must be called for each new matrix, and before any
    *  other method for generating perturbation factors.
    *
    *  Usually,
    *  the returned perturbation factors are zero, but if the system
    *  is thought to be structurally singular, they might be
    *  positive.  If the return value is false, no suitable
    *  perturbation could be found.
    */
   bool ConsiderNewSystem(
      Number& delta_x,
      Number& delta_s,
      Number& delta_c,
      Number& delta_d
   );

   /** This method returns perturbation factors for the case when the
    *  most recent factorization resulted in a singular matrix.
    *
    *  @return false, if no suitable perturbation could be found
    */
   bool PerturbForSingularity(
      Number& delta_x,
      Number& delta_s,
      Number& delta_c,
      Number& delta_d
   );

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
   CGPerturbationHandler(
      const CGPerturbationHandler&
   );

   /** Default Assignment Operator */
   void operator=(
      const CGPerturbationHandler&
   );
   ///@}

   /** Method to easily access CGPenalty data */
   CGPenaltyData& CGPenData()
   {
      CGPenaltyData& cg_pen_data = static_cast<CGPenaltyData&>(IpData().AdditionalData());
      DBG_ASSERT(dynamic_cast<CGPenaltyData*>(&IpData().AdditionalData()));
      return cg_pen_data;
   }

   /** Method to easily access CGPenalty calculated quantities */
   CGPenaltyCq& CGPenCq()
   {
      CGPenaltyCq& cg_pen_cq = static_cast<CGPenaltyCq&>(IpCq().AdditionalCq());
      DBG_ASSERT(dynamic_cast<CGPenaltyCq*>(&IpCq().AdditionalCq()));
      return cg_pen_cq;
   }

   /** The max reference value for scaling the penalty parameter */
   Number penalty_max_;
   /** Feasibility for perturbation in pure Newton method*/
   Number mult_diverg_feasibility_tol_;

};

} // namespace Ipopt

#endif

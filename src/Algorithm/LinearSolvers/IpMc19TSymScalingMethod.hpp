// Copyright (C) 2004, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-03-17

#ifndef __IPMC19TSYMSCALINGMETHOD_HPP__
#define __IPMC19TSYMSCALINGMETHOD_HPP__

#include "IpUtils.hpp"
#include "IpTSymScalingMethod.hpp"
#include "IpLibraryLoader.hpp"
#include "IpTypes.h"

// note that R,C,W are single-precision also in the double-precision version of MC19 (MC19AD)
// here we assume that float corresponds to Fortran's single precision
/// @since 3.14.0
#define IPOPT_DECL_MC19A(x) void (x)( \
   const ipindex* N,   \
   const ipindex* NZ,  \
   ipnumber*      A,   \
   ipindex*       IRN, \
   ipindex*       ICN, \
   float*         R,   \
   float*         C,   \
   float*         W    \
)

namespace Ipopt
{

/** Class for the method for computing scaling factors for symmetric
 *  matrices in triplet format, using MC19.
 */
class Mc19TSymScalingMethod: public TSymScalingMethod
{
public:
   /** @name Constructor/Destructor */
   ///@{
   Mc19TSymScalingMethod(
      SmartPtr<LibraryLoader> hslloader_  ///< @since 3.14.0
   ) : hslloader(hslloader_),
      mc19a(NULL)
   { }

   virtual ~Mc19TSymScalingMethod()
   { }
   ///@}

   virtual bool InitializeImpl(
      const OptionsList& options,
      const std::string& prefix
   );

   /** Method for computing the symmetric scaling factors, given the
    *  symmetric matrix in triplet (MA27) format.
    */
   virtual bool ComputeSymTScalingFactors(
      Index         n,
      Index         nnz,
      const Index*  airn,
      const Index*  ajcn,
      const Number* a,
      Number*       scaling_factors
   );

   /// set MC19 function to use for every instantiation of this class
   /// @since 3.14.0
   static void SetFunctions(
      IPOPT_DECL_MC19A(*mc19a)
   );

   /// get MC19A function that has been set via SetFunctions
   ///
   /// this does not return a MC19A that has been linked in or loaded from a library at runtime
   /// @since 3.14.0
   static IPOPT_DECL_MC19A(*GetMC19A());

private:
   /**@name Default Compiler Generated Methods (Hidden to avoid
    * implicit creation/calling).  These methods are not implemented
    * and we do not want the compiler to implement them for us, so we
    * declare them private and do not define them. This ensures that
    * they will not be implicitly created/called. */
   ///@{
   /** Copy Constructor */
   Mc19TSymScalingMethod(
      const Mc19TSymScalingMethod&
   );

   /** Default Assignment Operator */
   void operator=(
      const Mc19TSymScalingMethod&
   );

   /**@name MC19 function pointer
    * @{
    */
   SmartPtr<LibraryLoader> hslloader;

   IPOPT_DECL_MC19A(*mc19a);
   /**@} */
};

} // namespace Ipopt

#endif

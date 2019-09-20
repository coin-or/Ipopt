// Copyright 2009, 2011 Hans Pirnay
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Date   : 2009-05-06

#ifndef __ASPCALCULATOR_HPP__
#define __ASPCALCULATOR_HPP__

#include "IpAlgStrategy.hpp"
#include "SensSimpleBacksolver.hpp"
#include "SensSchurData.hpp"

namespace Ipopt
{

/** This class is the interface for implementations of any class that calculates the matrix \f$P =K^{-1}A\f$
 *  of the following matrix:
 *  \f[
 *  \left(\begin{array}{cc}
 *  K & A\\
 *  B & 0
 *  \end{array}\right)
 *  \f]
 */
class SIPOPTLIB_EXPORT PCalculator: public AlgorithmStrategyObject
{
public:
   PCalculator(
      SmartPtr<SensBacksolver> backsolver,
      SmartPtr<SchurData>      A_data
   )
      : backsolver_(backsolver),
        data_A_init(ConstPtr(A_data->MakeNewSchurDataCopy())),
        data_A_(A_data)
   { }

   virtual ~PCalculator()
   { }

   /* Overloaded from AlgorithmStrategyObject */
   virtual bool InitializeImpl(
      const OptionsList& /*options*/,
      const std::string& /*prefix*/
   )
   {
      return true;
   }

   /** Function to start the computation of  P from E_0 and KKT*/
   virtual bool ComputeP() = 0;

   /** Function to extract a SchurMatrix corresponding to $B K^{-1} A$.
    *
    * If B==NULL, use A=B.
    */
   virtual bool GetSchurMatrix(
      const SmartPtr<const SchurData>& B,
      SmartPtr<Matrix>&                S
   ) = 0;

   virtual void PrintImpl(
      const Journalist&  jnlst,
      EJournalLevel      level,
      EJournalCategory   category,
      const std::string& name,
      Index              indent,
      const std::string& prefix
   ) const = 0;

   void Print(
      const Journalist&  jnlst,
      EJournalLevel      level,
      EJournalCategory   category,
      const std::string& name,
      Index              indent = 0,
      const std::string& prefix = ""
   ) const
   {
      if( jnlst.ProduceOutput(level, category) )
      {
         PrintImpl(jnlst, level, category, name, indent, prefix);
      }
   }

   void Print(
      SmartPtr<const Journalist> jnlst,
      EJournalLevel              level,
      EJournalCategory           category,
      const std::string&         name,
      Index                      indent,
      const std::string&         prefix
   ) const
   {
      if( IsValid(jnlst) && jnlst->ProduceOutput(level, category) )
      {
         PrintImpl(*jnlst, level, category, name, indent, prefix);
      }
   }

   /** Accessor methods for data and backsolver.
    *
    *  This unconstness seems
    *  kind of dangerous but I don't think there is a way around it. Anyway,
    *  there is no difference between this and the IpData() method of AlgStrategy.
    */
   SmartPtr<SensBacksolver> Solver() const
   {
      return backsolver_;
   }

   SmartPtr<const SchurData> data_A() const
   {
      return ConstPtr(data_A_);
   }

   SmartPtr<SchurData> data_A_nonconst() const
   {
      return data_A_;
   }

   void reset_data_A()
   {
      data_A_ = data_A_init->MakeNewSchurDataCopy();
   }

private:

   SmartPtr<SensBacksolver> backsolver_;

   SmartPtr<const SchurData> data_A_init;
   SmartPtr<SchurData> data_A_;
};

}

#endif

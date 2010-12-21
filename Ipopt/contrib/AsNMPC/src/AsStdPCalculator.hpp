// Copyright 2009, 2010 Hans Pirnay
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Date   : 2009-05-06


#ifndef __ASSTDPCALCULATOR_HPP__
#define __ASSTDPCALCULATOR_HPP__

#include "AsPCalculator.hpp"
#include "AsSchurData.hpp"
#include "AsSimpleBacksolver.hpp"


namespace Ipopt
{

  class StdPCalculator : public PCalculator
  {
    /** This is the standard implementation: single processor, SimpleBacksolver als linear solver */

  public:
    StdPCalculator(SmartPtr<AsBacksolver> backsolver,
		   SmartPtr<SchurData> A_data);

    virtual ~StdPCalculator();

    /** Overloaded from AlgorithmStrategyObject */
    virtual bool InitializeImpl(const OptionsList& options,
                                const std::string& prefix);

    /** Virtual functions overloaded from PCalculator */
    virtual bool ComputeP();

    virtual bool GetSchurMatrix(const SmartPtr<const SchurData>& B, SmartPtr<Matrix>& S);

    virtual void PrintImpl(const Journalist& jnlst,
			   EJournalLevel level,
			   EJournalCategory category,
			   const std::string& name,
			   Index indent,
			   const std::string& prefix) const;

  private:
    /** standard constructor defined here so it can't be called */
    StdPCalculator();
    
    /** Rows of P = Rows of KKT */
    Index nrows_; 

    /** Cols of P = data_A_->GetNRowsAdded() */
    Index ncols_;
    
    /** Stores entries of P columnwise */
    Number* P_;

  };

}


#endif 

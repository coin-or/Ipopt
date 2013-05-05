// Copyright 2009, 2011 Hans Pirnay
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Date   : 2009-05-06

#ifndef __ASINDEXPCALCULATOR_HPP__
#define __ASINDEXPCALCULATOR_HPP__

#include "SensPCalculator.hpp"

namespace Ipopt
{
  /* Forward declarations */
  class PColumn;

  class IndexPCalculator : public PCalculator
  {
    /** This class is the implementation of the PCalculator that corresponds
     *  to IndexSchurData. It expects to be used with a kind of IndexSchurData. */

  public:

    IndexPCalculator(SmartPtr<SensBacksolver> backsolver,
		     SmartPtr<SchurData> A_data);

    virtual ~IndexPCalculator();

    /** Overloaded from PCalculator */
    virtual bool InitializeImpl(const OptionsList& options,
                                const std::string& prefix);

    virtual bool ComputeP();

    virtual bool GetSchurMatrix(const SmartPtr<const SchurData>& B, SmartPtr<Matrix>& S);

    virtual void PrintImpl(const Journalist& jnlst,
			   EJournalLevel level,
			   EJournalCategory category,
			   const std::string& name,
			   Index indent,
			   const std::string& prefix) const;

  private:

    /** Rows of P = Rows of KKT */
    Index nrows_;

    /** Cols of P */
    Index ncols_;

    std::map< Index, SmartPtr<PColumn> > cols_;

  };

  class PColumn : public ReferencedObject
  {
    /** This class provides an easy interface for PCalculators with data where columns are
     *  not necessarily in adjacent parts of memory. */

  public:
    PColumn(Index nrows, Number* values);

    virtual ~PColumn();

    virtual void GetSchurMatrixRows(const std::vector<Index>* row_idx_B, Number* S) const;

    virtual const Number* Values() const;

  private:

    Index nrows_;
    Number* val_;
  };

}

#endif

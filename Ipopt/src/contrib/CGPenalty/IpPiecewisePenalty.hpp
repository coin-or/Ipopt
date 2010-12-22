// Copyright (C) 2007 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// $Id$
//
// Authors:  Lifeng Chen/Zaiwen Wen      Columbia Univ

#ifndef __IPPIECEWISEPENALTY_HPP__
#define __IPPIECEWISEPENALTY_HPP__

#include "IpJournalist.hpp"
#include "IpDebug.hpp"
#include "IpOptionsList.hpp"
#include "IpIpoptCalculatedQuantities.hpp"
#include "IpBacktrackingLSAcceptor.hpp"
#include "IpPDSystemSolver.hpp"
#include <list>
#include <vector>

namespace Ipopt
{

  /** struct for one Piecewise Penalty entry. */
  typedef struct
        PiecewisePenEntry
  {
    Number pen_r;
    Number barrier_obj;
    Number infeasi;
  }
  PiecewisePenEntry;


  /** Class for the Piecewise Penalty.  This class contains all
   *  Piecewise Penalty entries.
   *  The entries are stored as the corner point, including the
   *  margin. */
  class PiecewisePenalty
  {
  public:
    /**@name Constructors/Destructors */
    //@{
    /** Default Constructor */
    PiecewisePenalty(Index dim);
    /** Default Destructor */
    ~PiecewisePenalty()
    {
      //ToDo figure out if that here is necessary
      //    Clear();
    }
    //@}


    //@{
    // Initialize Piecewise Penalty list
    bool IsPiecewisePenaltyListEmpty()
    {
      return PiecewisePenalty_list_.empty();
    }

    void InitPiecewisePenaltyList(Number pen_r, Number barrier_obj, Number infeasi)
    {
      AddEntry( pen_r, barrier_obj,  infeasi);
    }

    /** Check acceptability of given coordinates with respect
     *  to the Piecewise Penalty.  Returns true, if pair is acceptable
     */
    bool Acceptable(Number Fzconst, Number Fzlin);

    /** Get the value of the biggest barrier function so far*/
    Number BiggestBarr();

    /** Update Piecewise Penalty entry for given coordinates.
     */
    void UpdateEntry(Number barrier_obj, Number infeasi);

    /** Add a entry to the list
     */
    void AddEntry(Number pen_r, Number barrier_obj, Number infeasi)
    {
      PiecewisePenEntry TmpEntry;
      if (IsPiecewisePenaltyListEmpty()) {
        TmpEntry.pen_r = 0.0;
      }
      else {
        TmpEntry.pen_r = pen_r;
      }
      TmpEntry.barrier_obj = barrier_obj;
      TmpEntry.infeasi = infeasi;
      PiecewisePenalty_list_.push_back(TmpEntry);
    }

    /** Clear and reset the piecewise penalty list
     */
    void ResetList(Number pen_r, Number barrier_obj, Number infeasi)
    {
      PiecewisePenalty_list_.clear();
      AddEntry( pen_r, barrier_obj,  infeasi);
    }

    //@}

    /** Delete all Piecewise Penalty entries */
    void Clear()
    {
      PiecewisePenalty_list_.clear();
    }

    /** Print current Piecewise Penalty entries */
    void Print(const Journalist& jnlst);

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
    PiecewisePenalty();
    /** Copy Constructor */
    PiecewisePenalty(const PiecewisePenalty&);

    /** Overloaded Equals Operator */
    void operator=(const PiecewisePenalty&);
    //@}

    /** Dimension of the Piecewise Penalty (number of coordinates per entry) */
    Index dim_;

    /** The min penalty value for the piecewise penalty list */
    Number min_piece_penalty_;

    /** The max number of the break points in the piecewise penalty list */
    Index max_piece_number_;

    /** vector storing the Piecewise Penalty entries */
    std::vector<PiecewisePenEntry> PiecewisePenalty_list_;

  };







} // namespace Ipopt

#endif

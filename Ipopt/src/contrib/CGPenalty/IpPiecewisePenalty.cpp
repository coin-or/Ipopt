// Copyright (C) 2007, 2009 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// $Id$
//
// Authors:  Lifeng Chen/Zaiwen Wen      Columbia Univ

#include "IpPiecewisePenalty.hpp"
#include "IpJournalist.hpp"
#include "IpRestoPhase.hpp"
#include "IpAlgTypes.hpp"


namespace Ipopt
{

#if COIN_IPOPT_VERBOSITY > 0
  static const Index dbg_verbosity = 0;
#endif

  PiecewisePenalty::PiecewisePenalty(Index dim)
      :
      dim_(dim),
      min_piece_penalty_(0),  // make it regular option here or elsewhere?
      max_piece_number_(100)  // make it regular option here or elsewhere?
  {}

  bool PiecewisePenalty::Acceptable(Number Fzconst, Number Fzlin)
  {
    DBG_START_METH("PiebcewisePenalty::Acceptable", dbg_verbosity);
    DBG_ASSERT(!IsPiecewisePenaltyListEmpty());
    bool acceptable = false;
    std::vector<PiecewisePenEntry>::iterator iter;
    // Avoid the entries of the piecewise penalty list becoming too many.
    // Here, we require the entry number is less than or equal to max_piece_number,
    // unless the entries are accepted by the regular Armijo conditions.
    Index size = (Index)PiecewisePenalty_list_.size();
    if (size >= max_piece_number_) {
      Number trial_inf = Fzlin;
      Number trial_barrier = Fzconst;
      // First check the starting entry of the list.
      iter = PiecewisePenalty_list_.begin();
      Number value = iter->barrier_obj + iter->pen_r * iter->infeasi
                     - trial_barrier - iter->pen_r * trial_inf;
      if (value >= 0.) {
        iter++;
        value = iter->barrier_obj + iter->pen_r * iter->infeasi
                - trial_barrier - iter->pen_r * trial_inf;
        if (value <= 0.) {
          return false;
        }
      }
      // Then check the ending entry of the list.
      iter = PiecewisePenalty_list_.end();
      iter--;
      value = iter->barrier_obj + iter->pen_r * iter->infeasi
              - trial_barrier - iter->pen_r * trial_inf;
      if (value <= 0. && trial_inf <= iter->infeasi) {
        return false;
      }
      // Check the next to the ending entry.
      if (value >= 0. && trial_inf >= iter->infeasi) {
        iter-=1;
        value = iter->barrier_obj + iter->pen_r * iter->infeasi
                - trial_barrier - iter->pen_r * trial_inf;
        if (value <= 0.) {
          return false;
        }
      }
      // Finally, check the middle entries of the list.
      Number value_left, value_mid, value_right;
      for (iter = PiecewisePenalty_list_.begin() +1; iter != PiecewisePenalty_list_.end();
           iter++) {
        value_mid = iter->barrier_obj + iter->pen_r * iter->infeasi
                    - trial_barrier - iter->pen_r * trial_inf;
        iter++;
        value_right = iter->barrier_obj + iter->pen_r * iter->infeasi
                      - trial_barrier - iter->pen_r * trial_inf;
        iter-=2;
        value_left = iter->barrier_obj + iter->pen_r * iter->infeasi
                     - trial_barrier - iter->pen_r * trial_inf;
        iter++;
        if (value_left <= 0. && value_mid >= 0. && value_right <= 0.) {
          return false;
        }
      }
    }
    // Check if the trial point is acceptable to the piecewise list
    Number Fz;
    for (iter = PiecewisePenalty_list_.begin(); iter != PiecewisePenalty_list_.end();
         iter++) {
      Fz = Fzconst + iter->pen_r * (Fzlin - iter->infeasi) - iter->barrier_obj ;
      if (Fz < 0.) {
        acceptable = true;
        break;
      }
    }
    iter = PiecewisePenalty_list_.end() -1;
    if (acceptable == false && Fzlin < iter->infeasi) {
      acceptable = true;
    }
    return acceptable;
  }

  Number PiecewisePenalty::BiggestBarr()
  {
    DBG_START_METH("PiecewisePenalty::BiggestBarr", dbg_verbosity);
    DBG_ASSERT(!IsPiecewisePenaltyListEmpty());
    Number value = -1e20;
    if (PiecewisePenalty_list_.size() > 0) {
      std::vector<PiecewisePenEntry>::iterator iter;
      iter = PiecewisePenalty_list_.end();
      iter--;
      value = iter->barrier_obj;
    }
    return value;
  }

  void PiecewisePenalty::UpdateEntry(Number barrier_obj, Number infeasi)
  {

    DBG_START_METH("PiecewisePenalty::UpdateEntry", dbg_verbosity);
    DBG_ASSERT(!IsPiecewisePenaltyListEmpty());

    Number Gzi1, Gzi2;
    Number epsM = 0.0; //1e-20;
    Number TmpPen = 0.0;
    // construt a temp list, copy current list to the temp list
    std::vector<PiecewisePenEntry> TmpList(PiecewisePenalty_list_);
    // Erases the elements of current list
    PiecewisePenalty_list_.clear();
    std::vector<PiecewisePenEntry>::iterator iter = TmpList.begin(), iter2;
    Gzi1 = barrier_obj + iter->pen_r * ( infeasi - iter->infeasi) - iter->barrier_obj;
    for (; iter <= TmpList.end()-1; iter++) {
      // Be careful about this
      if ( TmpList.size() > 1 &&  iter <= TmpList.end()-2 ) {
        iter2 = iter+1;
        Gzi2 = barrier_obj + iter2->pen_r * ( infeasi - iter2->infeasi) - iter2->barrier_obj;
      }
      else {
        Gzi2 = infeasi - iter->infeasi;
      }
      if ( Gzi1 < -epsM && Gzi2 >= epsM ) {
        if ( IsPiecewisePenaltyListEmpty() ) {
          AddEntry(TmpPen, barrier_obj, infeasi);
        }
        if (Gzi2 > epsM) {
          TmpPen = (iter->barrier_obj - barrier_obj)/( infeasi - iter->infeasi );
          AddEntry(TmpPen, iter->barrier_obj, iter->infeasi);
        }
      }
      if (Gzi1 >= epsM && Gzi2 < -epsM) {
        if (Gzi1 > epsM) {
          AddEntry(iter->pen_r, iter->barrier_obj, iter->infeasi);
        }
        TmpPen = (iter->barrier_obj - barrier_obj)/(infeasi - iter->infeasi);
        AddEntry(TmpPen, barrier_obj, infeasi);
      }
      if (Gzi1 >= epsM && Gzi2 >= epsM) {
        AddEntry(iter->pen_r, iter->barrier_obj, iter->infeasi);
      }
      // handle the last point
      if ( iter == TmpList.end()-1 ) {
        if ( Gzi1 < - epsM && Gzi2 < - epsM ) {
          if ( IsPiecewisePenaltyListEmpty() ) {
            AddEntry(0.0, barrier_obj, infeasi);
          }
        }
      }
      Gzi1 = Gzi2;
    }
    dim_ = (Index)PiecewisePenalty_list_.size();
  }

  /*
  void PiecewisePenalty::Clear()
  {
    
    DBG_START_METH("PenaltyLineSearch::Filter::Clear", dbg_verbosity);
    while (!PiecewisePenalty_list_.empty()) {
      PiecewisePenEntry entry = PiecewisePenalty_list_.back();
      PiecewisePenalty_list_.pop_back();
      delete entry;
    }
    
  }
  */

  void PiecewisePenalty::Print(const Journalist& jnlst)
  {
    // DBG_START_METH("FilterLineSearch::Filter::Print", dbg_verbosity);
    jnlst.Printf(J_DETAILED, J_LINE_SEARCH,
                 "The current piecewise penalty has %d entries.\n",PiecewisePenalty_list_.size());
    jnlst.Printf(J_DETAILED, J_LINE_SEARCH,
                 "We only allow %d entries.\n", max_piece_number_);
    jnlst.Printf(J_DETAILED, J_LINE_SEARCH,
                 "The min piecewise penalty is %d .\n", min_piece_penalty_);
    if (!jnlst.ProduceOutput(J_DETAILED, J_LINE_SEARCH)) {
      return;
    }
    std::vector<PiecewisePenEntry>::iterator iter;
    Index count = 0;
    for (iter = PiecewisePenalty_list_.begin(); iter != PiecewisePenalty_list_.end();
         iter++) {
      if (count % 10 == 0) {
        jnlst.Printf(J_DETAILED, J_LINE_SEARCH,
                     "                pen_r                    barrier_obj            infeasi\n");
      }
      count++;
      jnlst.Printf(J_DETAILED, J_LINE_SEARCH, "%5d ", count);
      jnlst.Printf(J_DETAILED, J_LINE_SEARCH, "%23.16e %23.16e  %23.16e \n", iter->pen_r,
                   iter->barrier_obj, iter->infeasi);
    }
  }

} // namespace Ipopt

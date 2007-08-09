// Copyright (C) 2007 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id: IpIpoptData.hpp 988 2007-06-01 21:57:27Z andreasw $
//
// Authors:  Andreas Waechter           IBM     2007-06-04
//                   based on IpIpoptData.hpp

#ifndef __IPCGPENALTYDATA_HPP__
#define __IPCGPENALTYDATA_HPP__

#include "IpIteratesVector.hpp"
#include "IpOptionsList.hpp"

namespace Ipopt
{

  /** Class to organize all the additional data required by the
   *  Chen-Goldfarb penalty function algorithm. */
  class CGPenaltyData : public ReferencedObject
  {
  public:
    /**@name Constructors/Destructors */
    //@{
    /** Constructor */
    CGPenaltyData();

    /** Default destructor */
    ~CGPenaltyData();
    //@}

    /** This method must be called to initialize the global
     *  algorithmic parameters.  The parameters are taken from the
     *  OptionsList object. */    
    bool Initialize(const Journalist& jnlst,
		    const OptionsList& options,
		    const std::string& prefix);

    /** Initialize Data Structures */
    bool InitializeDataStructures();

    /** Delta for the Chen-Goldfarb search direction */
    SmartPtr<const IteratesVector> delta_cgpen() const;

    /** Set the delta_cgpen - like the trial point, this method copies
     *  the pointer for efficiency (no copy and to keep cache tags the
     *  same) so after you call set, you cannot modify the data.
     */
    void set_delta_cgpen(SmartPtr<IteratesVector>& delta_pen);

    /** Set the delta_cgpen - like the trial point, this method copies
     *  the pointer for efficiency (no copy and to keep cache tags the
     *  same) so after you call set, you cannot modify the data.  This
     *  is the version that is happy with a pointer to const
     *  IteratesVector.
     */
    void set_delta_cgpen(SmartPtr<const IteratesVector>& delta_pen);

    /** Delta for the fast Chen-Goldfarb search direction */
    SmartPtr<const IteratesVector> delta_cgfast() const;

    /** Set the delta_cgpen - like the trial point, this method copies
     *  the pointer for efficiency (no copy and to keep cache tags the
     *  same) so after you call set, you cannot modify the data.
     */
    void set_delta_cgfast(SmartPtr<IteratesVector>& delta_fast);

    /** @name Chen-Goldfarb step2.  Those fields can be used to store
     *  directions related to the Chen-Goldfarb algorithm */
    //@{
    bool HaveCgPenDeltas() const
    {
      return have_cgpen_deltas_;
    }
    void SetHaveCgPenDeltas(bool have_cgpen_deltas)
    {
      have_cgpen_deltas_ = have_cgpen_deltas;
    }

    bool HaveCgFastDeltas() const
    {
      return have_cgfast_deltas_;
    }
    void SetHaveCgFastDeltas(bool have_cgfast_deltas)
    {
      have_cgfast_deltas_ = have_cgfast_deltas;
    }
    //@}

    /** @name Public Methods for updating iterates */
    //@{
    /** Set the current iterate values from the
     *  trial values. */
    void AcceptTrialPoint();
    //@}

    Number curr_penalty() const
    {
      DBG_ASSERT(penalty_initialized_);
      return curr_penalty_;
    }
    void Set_penalty(Number penalty)
    {
      curr_penalty_ = penalty;
      penalty_initialized_ = true;
    }
    void SetPenaltyUninitialized()
    {
      penalty_initialized_ = false;
    }
    bool PenaltyInitialized() const
    {
      return penalty_initialized_;
    }

  private:

    /** @name Pure Chen-Goldfarb step for the penatly function.  This
     *  used to transfer the information about the step from the
     *  computation of the overall search direction to the line
     *  search. */
    //@{
    SmartPtr<const IteratesVector> delta_cgpen_;
    /** The following flag is set to true, if some other part of the
     *  algorithm has already computed the Chen-Goldfarb step.  This
     *  flag is reset when the AcceptTrialPoint method is called.
     * ToDo: we could cue off of a null delta_cgpen_;
     */
    bool have_cgpen_deltas_;
    //@}

    /** @name Fast Chen-Goldfarb step for the penatly function.  This
     *  used to transfer the information about the step from the
     *  computation of the overall search direction to the line
     *  search. */
    //@{
    SmartPtr<const IteratesVector> delta_cgfast_;
    /** The following flag is set to true, if some other part of the
     *  algorithm has already computed the fast Chen-Goldfarb step.
     *  This flag is reset when the AcceptTrialPoint method is called.
     *  * ToDo: we could cue off of a null delta_cgfast_;
     */
    bool have_cgfast_deltas_;
    //@}

    /**@name current penalty parameter */
    //@{
    Number curr_penalty_;
    bool penalty_initialized_;
    //@}

    /** flag indicating if Initialize method has been called (for
     *  debugging) */
    bool initialize_called_;

    /**@name Default Compiler Generated Methods
     * (Hidden to avoid implicit creation/calling).
     * These methods are not implemented and 
     * we do not want the compiler to implement
     * them for us, so we declare them private
     * and do not define them. This ensures that
     * they will not be implicitly created/called. */
    //@{
    /** Copy Constructor */
    CGPenaltyData(const CGPenaltyData&);

    /** Overloaded Equals Operator */
    void operator=(const CGPenaltyData&);
    //@}

#if COIN_IPOPT_CHECKLEVEL > 0
    /** Some debug flags to make sure vectors are not changed
     *  behind the CGPenaltyData's back
     */
    //@{
    TaggedObject::Tag debug_delta_cgpen_tag_;
    TaggedObject::Tag debug_delta_cgfast_tag_;
    TaggedObject::Tag debug_delta_cgpen_tag_sum_;
    TaggedObject::Tag debug_delta_cgfast_tag_sum_;
    //@}
#endif

  };

  inline
  SmartPtr<const IteratesVector> CGPenaltyData::delta_cgpen() const
  {
    DBG_ASSERT(IsNull(delta_cgpen_) || (delta_cgpen_->GetTag() == debug_delta_cgpen_tag_ && delta_cgpen_->GetTagSum() == debug_delta_cgpen_tag_sum_) );

    return delta_cgpen_;
  }

  inline
  SmartPtr<const IteratesVector> CGPenaltyData::delta_cgfast() const
  {
    DBG_ASSERT(IsNull(delta_cgfast_) || (delta_cgfast_->GetTag() == debug_delta_cgfast_tag_ && delta_cgfast_->GetTagSum() == debug_delta_cgfast_tag_sum_) );

    return delta_cgfast_;
  }

  inline
  void CGPenaltyData::set_delta_cgpen(SmartPtr<IteratesVector>& delta_cgpen)
  {
    delta_cgpen_ = ConstPtr(delta_cgpen);
#if COIN_IPOPT_CHECKLEVEL > 0

    if (IsValid(delta_cgpen)) {
      debug_delta_cgpen_tag_ = delta_cgpen->GetTag();
      debug_delta_cgpen_tag_sum_ = delta_cgpen->GetTagSum();
    }
    else {
      debug_delta_cgpen_tag_ = 0;
      debug_delta_cgpen_tag_sum_ = delta_cgpen->GetTagSum();
    }
#endif

    delta_cgpen = NULL;
  }

  inline
  void CGPenaltyData::set_delta_cgpen(SmartPtr<const IteratesVector>& delta_cgpen)
  {
    delta_cgpen_ = delta_cgpen;
#if COIN_IPOPT_CHECKLEVEL > 0

    if (IsValid(delta_cgpen)) {
      debug_delta_cgpen_tag_ = delta_cgpen->GetTag();
      debug_delta_cgpen_tag_sum_ = delta_cgpen->GetTagSum();
    }
    else {
      debug_delta_cgpen_tag_ = 0;
      debug_delta_cgpen_tag_sum_ = delta_cgpen->GetTagSum();
    }
#endif

    delta_cgpen = NULL;
  }

  inline
  void CGPenaltyData::set_delta_cgfast(SmartPtr<IteratesVector>& delta_cgfast)
  {
    delta_cgfast_ = ConstPtr(delta_cgfast);
#if COIN_IPOPT_CHECKLEVEL > 0

    if (IsValid(delta_cgfast)) {
      debug_delta_cgfast_tag_ = delta_cgfast->GetTag();
      debug_delta_cgfast_tag_sum_ = delta_cgfast->GetTagSum();
    }
    else {
      debug_delta_cgfast_tag_ = 0;
      debug_delta_cgfast_tag_sum_ = delta_cgfast->GetTagSum();
    }
#endif

    delta_cgfast = NULL;
  }

} // namespace Ipopt

#endif

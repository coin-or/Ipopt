// Copyright (C) 2004, 2005 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id$
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#include "IpVector.hpp"

namespace Ipopt
{

  void Vector::Print(SmartPtr<const Journalist> jnlst,
                     EJournalLevel level,
                     EJournalCategory category,
                     const std::string& name,
                     Index indent,
                     const std::string& prefix) const
  {
    if (IsValid(jnlst) && jnlst->ProduceOutput(level, category)) {
      PrintImpl(*jnlst, level, category, name, indent, prefix);
    }
  }

  void Vector::Print(const Journalist& jnlst,
                     EJournalLevel level,
                     EJournalCategory category,
                     const std::string& name,
                     Index indent,
                     const std::string& prefix) const
  {
    if (jnlst.ProduceOutput(level, category)) {
      PrintImpl(jnlst, level, category, name, indent, prefix);
    }
  }

  /* Prototype implementation for specialized functions */
  void Vector::AddTwoVectorsImpl(Number a, const Vector& v1,
                                 Number b, const Vector& v2, Number c)
  {
    if (c==0.) {
      if (a==1.) {
        Copy(v1);
        if (b!=0.) {
          Axpy(b, v2);
        }
      }
      else if (a==0.) {
        if (b==0.) {
          Set(0.);
        }
        else {
          Copy(v2);
          if (b!=1.) {
            Scal(b);
          }
        }
      }
      else {
        if (b==1.) {
          Copy(v2);
          Axpy(a, v1);
        }
        else if (b==0.) {
          Copy(v1);
          Scal(a);
        }
        else {
          Copy(v1);
          Scal(a);
          Axpy(b, v2);
        }
      }
    }
    else { /* c==0. */
      if (c!=1.) {
        Scal(c);
      }
      if (a!=0.) {
        Axpy(a, v1);
      }
      if (b!=0.) {
        Axpy(b, v2);
      }
    }
  }

  Number Vector::FracToBoundImpl(const Vector& delta, Number tau) const
  {
    DBG_ASSERT(tau>=0.);
    DBG_ASSERT(Dim() == delta.Dim());
    if (Dim() == 0 && delta.Dim() == 0) {
      return 1.0;
    }

    SmartPtr<Vector> inv_alpha_bar = MakeNew();
    inv_alpha_bar->AddOneVector(-1.0/tau, delta, 0.);
    inv_alpha_bar->ElementWiseDivide(*this);

    Number alpha = inv_alpha_bar->Max();
    if (alpha > 0) {
      alpha = Ipopt::Min(1.0/alpha, 1.0);
    }
    else {
      alpha = 1.0;
    }

    return alpha;
  }

  void Vector::AddVectorQuotientImpl(Number a, const Vector& z,
                                     const Vector& s, Number c)
  {
    DBG_ASSERT(Dim() == z.Dim());
    DBG_ASSERT(Dim() == s.Dim());

    if (c==0.) {
      AddOneVector(a, z, 0.);
      ElementWiseDivide(s);
    }
    else {
      SmartPtr<Vector> tmp = MakeNew();
      tmp->Copy(z);
      tmp->ElementWiseDivide(s);
      AddOneVector(a, *tmp, c);
    }
  }

} // namespace Ipopt

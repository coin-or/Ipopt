// Copyright (C) 2008 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// $Id$
//
// Authors:  Andreas Waechter            IBM    2008-09-19

#include "IpIterativeSolverTerminationTester.hpp"
#include "IpTripletHelper.hpp"

namespace Ipopt
{

#if COIN_IPOPT_VERBOSITY > 0
  static const Index dbg_verbosity = 0;
#endif

  void
  IterativeSolverTerminationTester::GetVectors(Index ndim, const Number* array,
      SmartPtr<const Vector>& comp_x,
      SmartPtr<const Vector>& comp_s,
      SmartPtr<const Vector>& comp_c,
      SmartPtr<const Vector>& comp_d)
  {
    DBG_ASSERT(ndim==IpData().curr()->x()->Dim()+IpData().curr()->s()->Dim()+
               IpData().curr()->y_c()->Dim()+IpData().curr()->y_d()->Dim());

    // x
    SmartPtr<Vector> sol_x = IpData().curr()->x()->MakeNew();
    Index dim = sol_x->Dim();
    TripletHelper::PutValuesInVector(dim, array, *sol_x);
    array += dim;
    comp_x = ConstPtr(sol_x);

    // s
    SmartPtr<Vector> sol_s = IpData().curr()->s()->MakeNew();
    dim = sol_s->Dim();
    TripletHelper::PutValuesInVector(dim, array, *sol_s);
    array += dim;
    comp_s = ConstPtr(sol_s);

    // c
    SmartPtr<Vector> sol_c = IpData().curr()->y_c()->MakeNew();
    dim = sol_c->Dim();
    TripletHelper::PutValuesInVector(dim, array, *sol_c);
    array += dim;
    comp_c = ConstPtr(sol_c);

    // d
    SmartPtr<Vector> sol_d = IpData().curr()->y_d()->MakeNew();
    dim = sol_d->Dim();
    TripletHelper::PutValuesInVector(dim, array, *sol_d);
    array += dim;
    comp_d = ConstPtr(sol_d);
  }

} // namespace Ipopt

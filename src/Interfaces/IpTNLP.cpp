// Copyright (C) 2021 COIN-OR Foundation
// All Rights Reserved.
// This code is published under the Eclipse Public License.

#include "IpoptConfig.h"
#include "IpTNLP.hpp"
#include "IpIpoptData.hpp"
#include "IpIpoptCalculatedQuantities.hpp"
#include "IpOrigIpoptNLP.hpp"
#include "IpRestoIpoptNLP.hpp"
#include "IpTNLPAdapter.hpp"

namespace Ipopt
{

bool TNLP::get_curr_iterate(
   const IpoptData*           ip_data,
   IpoptCalculatedQuantities* ip_cq,
   Index                      n,
   Number*                    x,
   Number*                    z_L,
   Number*                    z_U,
   Index                      m,
   Number*                    lambda
   ) const
{
   if( ip_data == NULL || !IsValid(ip_data->curr()) )
      return false;
   if( ip_cq == NULL )
      return false;

   Ipopt::OrigIpoptNLP* orignlp;
   TNLPAdapter* tnlp_adapter;
   Index n_full;
   Index m_full;

   // check whether we use a OrigIpoptNLP
   orignlp = dynamic_cast<OrigIpoptNLP*>(GetRawPtr(ip_cq->GetIpoptNLP()));
   if( orignlp != NULL )
   {
      tnlp_adapter = dynamic_cast<TNLPAdapter*>(GetRawPtr(orignlp->nlp()));
      if( tnlp_adapter == NULL )
         return false;

      tnlp_adapter->GetFullDimensions(n_full, m_full);
      if( n != n_full && (x != NULL || z_L != NULL || z_U != NULL) )
         THROW_EXCEPTION(IpoptException, "Incorrect dimension of x given to TNLP::get_curr_iterate().\n");
      if( m != m_full && lambda != NULL )
         THROW_EXCEPTION(IpoptException, "Incorrect dimension of g(x) given to TNLP::get_curr_iterate().\n");

      // resort Ipopt-internal x to TNLP-version of x, i.e., reinsert fixed variables
      if( x != NULL )
         tnlp_adapter->ResortX(*ip_data->curr()->x(), x);

      // resort Ipopt-internal variable duals to TNLP-version
      if( z_L != NULL || z_U != NULL )
         tnlp_adapter->ResortBnds(*ip_data->curr()->z_L(), z_L, *ip_data->curr()->z_U(), z_U, true);

      // resort Ipopt-internal constraint duals to TNLP-version
      if( lambda != NULL )
         tnlp_adapter->ResortG(*ip_data->curr()->y_c(), *ip_data->curr()->y_d(), lambda);

      return true;
   }

   // check whether we are in restoration phase, so get a RestoIpoptNLP
   Ipopt::RestoIpoptNLP* restonlp = dynamic_cast<RestoIpoptNLP*>(GetRawPtr(ip_cq->GetIpoptNLP()));
   if( restonlp != NULL )
   {
      if( (orignlp = dynamic_cast<OrigIpoptNLP*>(&restonlp->OrigIpNLP())) == NULL )
         return false;

      tnlp_adapter = dynamic_cast<TNLPAdapter*>(GetRawPtr(orignlp->nlp()));
      if( tnlp_adapter == NULL )
         return false;

      tnlp_adapter->GetFullDimensions(n_full, m_full);
      if( n != n_full && (x != NULL || z_L != NULL || z_U != NULL) )
         THROW_EXCEPTION(IpoptException, "Incorrect dimension of x given to TNLP::get_curr_iterate().\n");
      if( m != m_full && lambda != NULL )
         THROW_EXCEPTION(IpoptException, "Incorrect dimension of g(x) given to TNLP::get_curr_iterate().\n");

      const CompoundVector* c_vec;

      // restoration phase NLP optimizes over x and slack p and n
      if( x != NULL )
      {
         // get x from the compound vector (x,p,n)
         DBG_ASSERT(dynamic_cast<const CompoundVector*>(GetRawPtr(ip_data->curr()->x())) != NULL);
         c_vec = static_cast<const CompoundVector*>(GetRawPtr(ip_data->curr()->x()));
         SmartPtr<const Vector> x_only = c_vec->GetComp(0);
         DBG_ASSERT(IsValid(x_only));

         // resort Ipopt-internal x to TNLP-version of x, i.e., reinsert fixed variables
         tnlp_adapter->ResortX(*x_only, x);
      }

      if( z_L != NULL || z_U != NULL )
      {
         // get duals for x from the compound vector for the duals of (x,p,n)
         DBG_ASSERT(dynamic_cast<const CompoundVector*>(GetRawPtr(ip_data->curr()->z_L())) != NULL);
         c_vec = static_cast<const CompoundVector*>(GetRawPtr(ip_data->curr()->z_L()));
         SmartPtr<const Vector> z_L_only = c_vec->GetComp(0);
         DBG_ASSERT(IsValid(z_L_only));

         DBG_ASSERT(dynamic_cast<const CompoundVector*>(GetRawPtr(ip_data->curr()->z_U())) != NULL);
         c_vec = static_cast<const CompoundVector*>(GetRawPtr(ip_data->curr()->z_U()));
         SmartPtr<const Vector> z_U_only = c_vec->GetComp(0);
         DBG_ASSERT(IsValid(z_U_only));

         tnlp_adapter->ResortBnds(*z_L_only, z_L, *z_U_only, z_U, true);
      }

      if( lambda != NULL )
      {
         // yc and yd should be trivial compound vectors, see RestoIpoptNLP::h()
         DBG_ASSERT(dynamic_cast<const CompoundVector*>(GetRawPtr(ip_data->curr()->y_c())) != NULL);
         c_vec = static_cast<const CompoundVector*>(GetRawPtr(ip_data->curr()->y_c()));
         DBG_ASSERT(c_vec->NComps() == 1);
         SmartPtr<const Vector> yc_only = c_vec->GetComp(0);
         DBG_ASSERT(IsValid(yc_only));

         DBG_ASSERT(dynamic_cast<const CompoundVector*>(GetRawPtr(ip_data->curr()->y_d())) != NULL);
         c_vec = static_cast<const CompoundVector*>(GetRawPtr(ip_data->curr()->y_d()));
         DBG_ASSERT(c_vec->NComps() == 1);
         SmartPtr<const Vector> yd_only = c_vec->GetComp(0);
         DBG_ASSERT(IsValid(yd_only));

         // resort Ipopt-internal to TNLP-version
         tnlp_adapter->ResortG(*yc_only, *yd_only, lambda);
      }

      return true;
   }

   // if it is neither OrigIpoptNLP nor RestoIpoptNLP, then we don't know how to retrieve x
   return false;
}

} // namespace Ipopt

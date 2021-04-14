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
#include "IpDenseVector.hpp"
#include "IpBlas.hpp"

namespace Ipopt
{

// TODO return g(x) as well
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

         DBG_ASSERT(dynamic_cast<const CompoundVector*>(GetRawPtr(ip_data->curr()->y_c())) != NULL);
         c_vec = static_cast<const CompoundVector*>(GetRawPtr(ip_data->curr()->y_c()));
         DBG_ASSERT(c_vec->NComps() == 1);
         SmartPtr<const Vector> yc_only = c_vec->GetComp(0);
         DBG_ASSERT(IsValid(yc_only));

         tnlp_adapter->ResortBoundMultipliers(*yc_only, *z_L_only, z_L, *z_U_only, z_U);
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

bool TNLP::get_curr_violations(
   const IpoptData*           ip_data,
   IpoptCalculatedQuantities* ip_cq,
   bool                       scaled,
   Index                      n,
   Number*                    compl_x_L,
   Number*                    compl_x_U,
   Number*                    grad_lag_x,
   Index                      m,
   Number*                    nlp_constraint_violation,
   Number*                    compl_g
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
      if( n != n_full && (compl_x_L != NULL || compl_x_U != NULL || grad_lag_x != NULL) )
         THROW_EXCEPTION(IpoptException, "Incorrect dimension of x given to TNLP::get_curr_violations().\n");
      if( m != m_full && (nlp_constraint_violation != NULL || compl_g != NULL) )
         THROW_EXCEPTION(IpoptException, "Incorrect dimension of g(x) given to TNLP::get_curr_violations().\n");

      // TODO handle fixed variable treatment being make_constraint below

      if( compl_x_L != NULL || compl_x_U != NULL )
      {
         // this should give XZe from (5)
         tnlp_adapter->ResortBnds(*ip_cq->curr_compl_x_L(), compl_x_L, *ip_cq->curr_compl_x_U(), compl_x_U, true);

         if( !scaled )
         {
            // IpoptCalculatedQuantities::unscaled_curr_complementarity() calls unapply_obj_scaling() on norm of complementarity vector
            // so we multiply here each entry of the complementarity vector with the same factor
            Number obj_unscal = orignlp->NLP_scaling()->unapply_obj_scaling(1.0);
            if( compl_x_L != NULL )
               IpBlasScal(n, obj_unscal, compl_x_L, 1);
            if( compl_x_U != NULL )
               IpBlasScal(n, obj_unscal, compl_x_U, 1);
         }
      }

      if( grad_lag_x != NULL )
      {
         if( scaled )
         {
            tnlp_adapter->ResortX(*ip_cq->curr_grad_lag_x(), grad_lag_x, false);
         }
         else
         {
            // adapted from IpoptCalculatedQuantities::unscaled_curr_dual_infeasibility()
            SmartPtr<const Vector> intern_grad_lag_x = orignlp->NLP_scaling()->unapply_grad_obj_scaling(ip_cq->curr_grad_lag_x());
            tnlp_adapter->ResortX(*intern_grad_lag_x, grad_lag_x, false);
         }
      }

      if( nlp_constraint_violation != NULL || compl_g != NULL )
      {
         // adapted from IpoptCalculatedQuantities::curr_(unscaled_)nlp_constraint_violation

         // violation of d_L <= d(x)
         SmartPtr<Vector> d_viol_L;
         if( orignlp->d_L()->Dim() > 0 )
         {
            d_viol_L = ip_cq->curr_d()->MakeNewCopy();
            orignlp->Pd_L()->MultVector(1., *orignlp->d_L(), -1., *d_viol_L);  // d_L - d
            if( !scaled )
               d_viol_L = orignlp->NLP_scaling()->unapply_vector_scaling_d_NonConst(ConstPtr(d_viol_L));
         }
         else
         {
            d_viol_L = ip_cq->curr_d()->MakeNew();
            d_viol_L->Set(0.);
         }

         // violation of d(x) <= d_U
         SmartPtr<Vector> d_viol_U;
         if( orignlp->d_U()->Dim() > 0 )
         {
            d_viol_U = ip_cq->curr_d()->MakeNewCopy();
            orignlp->Pd_U()->MultVector(-1., *orignlp->d_U(), 1., *d_viol_U);  // d - d_U
            if( !scaled )
               d_viol_U = orignlp->NLP_scaling()->unapply_vector_scaling_d_NonConst(ConstPtr(d_viol_U));
         }
         else
         {
            d_viol_U = ip_cq->curr_d()->MakeNew();
            d_viol_U->Set(0.);
         }

         // c(x) = 0, d_L <= d(x) <= d_U should result in complementarity constraints
         // y_c*c(x), (d(x)-d_L)*y_d^+, (d_U-d(x))*y_d^-,  where y_d^+ = max(0,y_d) ~= v_L, y_d^- = -min(0,y_d) ~= v_U   (I hope these are the correct signs)
         // we will merge the latter two into one vector, taking the nonzero entries, i.e., (d(x)-d_L)*y_d^+ + (d_U-d(x))*y_d^-
         // this is then also consistent with 0 <= c(x) <= 0, since c(x)*y_c^+ + (-c(x))*y_c^- = c(x)*(y_c^+ - y_c^-) = c(x)*y_c
         if( compl_g != NULL )
         {
            SmartPtr<Vector> yd_pos = ip_data->curr()->y_d()->MakeNewCopy();
            SmartPtr<Vector> yd_neg = ip_data->curr()->y_d()->MakeNewCopy();
            SmartPtr<Vector> zero = yd_pos->MakeNew();
            zero->Set(0.);
            yd_pos->ElementWiseMax(*zero);
            yd_neg->ElementWiseMin(*zero);  // -y_d^-

            yd_pos->ElementWiseMultiply(*d_viol_L);  // (d_L-d(x)) * y_d^+
            yd_neg->ElementWiseMultiply(*d_viol_U);  // (d(x)-d_U) * (-y_d^-)

            yd_pos->Scal(-1.0);          // (d(x)-d_L) * y_d^+
            yd_pos->Axpy(1.0, *yd_neg);  // (d(x)-d_L) * y_d^+ + (d_U-d(x)) * y_d^-

            SmartPtr<Vector> c_compl;
            if( scaled )
               c_compl = ip_cq->curr_c()->MakeNewCopy();
            else
               c_compl = ip_cq->unscaled_curr_c()->MakeNewCopy();
            c_compl->ElementWiseMultiply(*ip_data->curr()->y_c());      // c(x)*y_c

            tnlp_adapter->ResortG(*c_compl, *yd_pos, compl_g);
         }

         if( nlp_constraint_violation != NULL )
         {
            // violation of c(x) = 0
            SmartPtr<Vector> c_viol;
            if( scaled )
               c_viol = ip_cq->curr_c()->MakeNewCopy();
            else
               c_viol = ip_cq->unscaled_curr_c()->MakeNewCopy();
            c_viol->ElementWiseAbs();  // |c(x)|

            // violation of d_L <= d(x) <= d_U:   d_viol_L := max(d_viol_L, d_viol_U, 0)
            d_viol_L->ElementWiseMax(*d_viol_U);
            SmartPtr<Vector> tmp = d_viol_L->MakeNew();
            tmp->Set(0.);
            d_viol_L->ElementWiseMax(*tmp);

            tnlp_adapter->ResortG(*c_viol, *d_viol_L, nlp_constraint_violation);
         }
      }

      return true;
   }

   return false;
}


} // namespace Ipopt

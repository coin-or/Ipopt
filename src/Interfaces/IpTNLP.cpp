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

#include <cstring>

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
   Number*                    g,
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
      if( m != m_full && (lambda != NULL || g != NULL) )
         THROW_EXCEPTION(IpoptException, "Incorrect dimension of g(x) given to TNLP::get_curr_iterate().\n");

      // resort Ipopt-internal x to TNLP-version of x, i.e., reinsert fixed variables
      if( x != NULL )
         tnlp_adapter->ResortX(*ip_data->curr()->x(), x);

      // resort Ipopt-internal variable duals to TNLP-version
      if( z_L != NULL || z_U != NULL )
         tnlp_adapter->ResortBnds(*ip_data->curr()->z_L(), z_L, *ip_data->curr()->z_U(), z_U, true);

      // resort Ipopt-interval constraint activity to TNLP-version
      if( g != NULL )
         tnlp_adapter->ResortG(*ip_cq->unscaled_curr_c(), *ip_cq->unscaled_curr_d(), g, true);

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
      if( m != m_full && (g != NULL || lambda != NULL) )
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

      if( g != NULL )
      {
         // get nc, pc, nd, pd from the compound vector (x,nc,pc,nd,pd,...)
         DBG_ASSERT(dynamic_cast<const CompoundVector*>(GetRawPtr(ip_data->curr()->x())) != NULL);
         c_vec = static_cast<const CompoundVector*>(GetRawPtr(ip_data->curr()->x()));
         SmartPtr<const Vector> nc_only = c_vec->GetComp(1);
         SmartPtr<const Vector> pc_only = c_vec->GetComp(2);
         SmartPtr<const Vector> nd_only = c_vec->GetComp(3);
         SmartPtr<const Vector> pd_only = c_vec->GetComp(4);
         DBG_ASSERT(IsValid(nc_only));
         DBG_ASSERT(IsValid(pc_only));
         DBG_ASSERT(IsValid(nd_only));
         DBG_ASSERT(IsValid(pd_only));

         // get scaled c from restonlp
         DBG_ASSERT(dynamic_cast<const CompoundVector*>(GetRawPtr(ip_cq->curr_c())) != NULL);
         c_vec = static_cast<const CompoundVector*>(GetRawPtr(ip_cq->curr_c()));
         SmartPtr<Vector> c_resto = c_vec->GetComp(0)->MakeNewCopy();

         // get scaled d from restonlp
         DBG_ASSERT(dynamic_cast<const CompoundVector*>(GetRawPtr(ip_cq->curr_d())) != NULL);
         c_vec = static_cast<const CompoundVector*>(GetRawPtr(ip_cq->curr_d()));
         SmartPtr<Vector> d_resto = c_vec->GetComp(0)->MakeNewCopy();

         // undo addition of slacks nc-pc
         c_resto->Axpy(-1.0, *nc_only);
         c_resto->Axpy(1.0, *pc_only);

         // undo addition of slacks nd-pd
         d_resto->Axpy(-1.0, *nd_only);
         d_resto->Axpy(1.0, *pd_only);

         // unscale using scaling in original NLP
         c_resto = orignlp->NLP_scaling()->unapply_vector_scaling_c_NonConst(c_resto);
         d_resto = orignlp->NLP_scaling()->unapply_vector_scaling_d_NonConst(d_resto);

         tnlp_adapter->ResortG(*c_resto, *d_resto, g, true);
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

      SmartPtr<Vector> c_compl;

      int n_x_fixed;
      Index* x_fixed_map;
      TNLPAdapter::FixedVariableTreatmentEnum fixed_variable_treatment;
      tnlp_adapter->GetFixedVariables(n_x_fixed, x_fixed_map, fixed_variable_treatment);

      if( compl_x_L != NULL || compl_x_U != NULL )
      {
         // this should give XZe from (5)

         if( n_x_fixed > 0 && fixed_variable_treatment == TNLPAdapter::MAKE_CONSTRAINT )
         {
            // if fixed variables are treated as constraints, we have equality constraints x - x_L = 0 at the end of c(x)=0
            // we can then use c(x)*y_c as complementarity for these variables

            if( scaled )
               c_compl = ip_cq->curr_c()->MakeNewCopy();
            else
               c_compl = ip_cq->unscaled_curr_c()->MakeNewCopy();
            c_compl->ElementWiseMultiply(*ip_data->curr()->y_c());      // c(x)*y_c

            tnlp_adapter->ResortBoundMultipliers(*c_compl, *ip_cq->curr_compl_x_L(), compl_x_L, *ip_cq->curr_compl_x_U(), compl_x_U);
         }
         else
         {
            tnlp_adapter->ResortBnds(*ip_cq->curr_compl_x_L(), compl_x_L, *ip_cq->curr_compl_x_U(), compl_x_U, true);
         }

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
         // this will set the derivative of the Lagrangian w.r.t. fixed variables to 0 (for any fixed_variables_treatment)
         // the actual values are nowhere stored within Ipopt Data or CQ, since TNLPAdapter does not seem to pass them on
         // we may have to reevaluate and capture them in TNLPAdapater
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

         if( n_x_fixed > 0 && fixed_variable_treatment == TNLPAdapter::MAKE_CONSTRAINT )
         {
            // if fixed_variable_treatment is make_constraint, then fixed variable contribute y_c*x to the Lagrangian
            // so we should subtract y_c for these entries; let's hope we don't need to deal with scaling
            const DenseVector* y_c = static_cast<const DenseVector*>(GetRawPtr(ip_data->curr()->y_c()));
            DBG_ASSERT(dynamic_cast<const DenseVector*>(GetRawPtr(ip_data->curr()->y_c())) != NULL);
            DBG_ASSERT(y_c->Dim() >= n_x_fixed);
            if( y_c->IsHomogeneous() )
               for( Index i = 0; i < n_x_fixed; ++i )
                  grad_lag_x[x_fixed_map[i]] -= y_c->Scalar();
            else
               for( Index i = 0; i < n_x_fixed; ++i )
                  grad_lag_x[x_fixed_map[i]] -= y_c->Values()[y_c->Dim()-n_x_fixed+i];
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

         // c(x) = 0, d_L <= d(x) <= d_U should result in complementarities
         // y_c*c(x), (d(x)-d_L)*y_d^-, (d_U-d(x))*y_d^+,  where y_d^+ = max(0,y_d), y_d^- = max(0,-y_d)   (I took the signs from TNLPAdapter::ResortBoundMultipliers)
         // we will merge the latter two into one vector, taking the nonzero entries, i.e., (d(x)-d_L)*y_d^- + (d_U-d(x))*y_d^+
         // to be consistent, it looks like we need to negate for 0 <= c(x) <= 0, since c(x)*y_c^- + (-c(x))*y_c^+ = c(x)*(y_c^- - y_c^+) = -c(x)*y_c
         if( compl_g != NULL )
         {
            SmartPtr<Vector> yd_pos = ip_data->curr()->y_d()->MakeNewCopy();
            SmartPtr<Vector> yd_neg = ip_data->curr()->y_d()->MakeNewCopy();
            SmartPtr<Vector> zero = yd_pos->MakeNew();
            zero->Set(0.);
            yd_pos->ElementWiseMax(*zero);
            yd_neg->ElementWiseMin(*zero);  // -y_d^-

            yd_pos->ElementWiseMultiply(*d_viol_U);  // (d(x)-d_U) * y_d^+
            yd_neg->ElementWiseMultiply(*d_viol_L);  // (d_L-d(x)) * (-y_d^-)

            yd_neg->Axpy(-1.0, *yd_pos); // (d(x)-d_L) * y_d^- + (d_U-d(x)) * y_d^+

            if( !IsValid(c_compl) )
            {
               if( scaled )
                  c_compl = ip_cq->curr_c()->MakeNewCopy();
               else
                  c_compl = ip_cq->unscaled_curr_c()->MakeNewCopy();
               c_compl->ElementWiseMultiply(*ip_data->curr()->y_c());      // c(x)*y_c
            }
            c_compl->Scal(-1.0);  // -c(x)*y_c

            tnlp_adapter->ResortG(*c_compl, *yd_neg, compl_g);
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
      if( n != n_full && (compl_x_L != NULL || compl_x_U != NULL || grad_lag_x != NULL) )
         THROW_EXCEPTION(IpoptException, "Incorrect dimension of x given to TNLP::get_curr_violations().\n");
      if( m != m_full && (nlp_constraint_violation != NULL || compl_g != NULL) )
         THROW_EXCEPTION(IpoptException, "Incorrect dimension of g(x) given to TNLP::get_curr_violations().\n");

      const CompoundVector* c_vec;
      SmartPtr<Vector> c;        // will hold c(x) of original NLP
      SmartPtr<Vector> c_compl;  // will hold c(x)*y_c of original NLP

      // if fixed variables are treated as constraints, we have equality constraints x - x_L = 0 at the end of c(x)=0
      // we can then use c(x)*y_c as complementarity for these variables
      // but for restoration nlp, we actually have c(x)-pc+nc=0, so need to undo this, too

      // get activity for c(x)
      DBG_ASSERT(dynamic_cast<const CompoundVector*>(GetRawPtr(ip_cq->curr_c())) != NULL);
      c_vec = static_cast<const CompoundVector*>(GetRawPtr(ip_cq->curr_c()));
      c = c_vec->GetComp(0)->MakeNewCopy();

      // get nc, pc
      DBG_ASSERT(dynamic_cast<const CompoundVector*>(GetRawPtr(ip_data->curr()->x())) != NULL);
      c_vec = static_cast<const CompoundVector*>(GetRawPtr(ip_data->curr()->x()));
      SmartPtr<const Vector> nc_only = c_vec->GetComp(1);
      SmartPtr<const Vector> pc_only = c_vec->GetComp(2);
      DBG_ASSERT(IsValid(nc_only));
      DBG_ASSERT(IsValid(pc_only));

      // undo addition of slacks nc-pc
      c->Axpy(-1.0, *nc_only);
      c->Axpy(1.0, *pc_only);

      // unscale using scaling in original NLP
      if( !scaled )
         c = orignlp->NLP_scaling()->unapply_vector_scaling_c_NonConst(c);

      // get duals for c(x)=0
      DBG_ASSERT(dynamic_cast<const CompoundVector*>(GetRawPtr(ip_data->curr()->y_c())) != NULL);
      c_vec = static_cast<const CompoundVector*>(GetRawPtr(ip_data->curr()->y_c()));
      DBG_ASSERT(c_vec->NComps() == 1);
      SmartPtr<const Vector> yc_only = c_vec->GetComp(0);
      DBG_ASSERT(IsValid(yc_only));

      c_compl = c->MakeNewCopy();
      c_compl->ElementWiseMultiply(*yc_only);      // c(x)*y_c

      if( compl_x_L != NULL || compl_x_U != NULL )
      {
         // get duals z_L and z_U for x from the compound vector for the duals of (x,p,n)
         DBG_ASSERT(dynamic_cast<const CompoundVector*>(GetRawPtr(ip_data->curr()->z_L())) != NULL);
         c_vec = static_cast<const CompoundVector*>(GetRawPtr(ip_data->curr()->z_L()));
         SmartPtr<const Vector> z_L_only = c_vec->GetComp(0);
         DBG_ASSERT(IsValid(z_L_only));

         DBG_ASSERT(dynamic_cast<const CompoundVector*>(GetRawPtr(ip_data->curr()->z_U())) != NULL);
         c_vec = static_cast<const CompoundVector*>(GetRawPtr(ip_data->curr()->z_U()));
         SmartPtr<const Vector> z_U_only = c_vec->GetComp(0);
         DBG_ASSERT(IsValid(z_U_only));

         // get slacks for x
         DBG_ASSERT(dynamic_cast<const CompoundVector*>(GetRawPtr(ip_cq->curr_slack_x_L())) != NULL);
         c_vec = static_cast<const CompoundVector*>(GetRawPtr(ip_cq->curr_slack_x_L()));
         SmartPtr<const Vector> slack_x_L = c_vec->GetComp(0);
         DBG_ASSERT(IsValid(slack_x_L));

         DBG_ASSERT(dynamic_cast<const CompoundVector*>(GetRawPtr(ip_cq->curr_slack_x_U())) != NULL);
         c_vec = static_cast<const CompoundVector*>(GetRawPtr(ip_cq->curr_slack_x_U()));
         SmartPtr<const Vector> slack_x_U = c_vec->GetComp(0);
         DBG_ASSERT(IsValid(slack_x_U));

         // calculate complementarity for x_L and x_U
         SmartPtr<Vector> compl_x_L_v = slack_x_L->MakeNewCopy();
         compl_x_L_v->ElementWiseMultiply(*z_L_only);
         SmartPtr<Vector> compl_x_U_v = slack_x_U->MakeNewCopy();
         compl_x_U_v->ElementWiseMultiply(*z_U_only);

         tnlp_adapter->ResortBoundMultipliers(*c_compl, *compl_x_L_v, compl_x_L, *compl_x_U_v, compl_x_U);

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
         // this looks like it would require reevaluating to get the gradient for the Lagrangian in the original NLP
         memset(grad_lag_x, 0, m*sizeof(Number));
      }

      if( nlp_constraint_violation != NULL || compl_g != NULL )
      {
         // adapted from IpoptCalculatedQuantities::curr_(unscaled_)nlp_constraint_violation

         // get nd, pd from the compound vector (x,nc,pc,nd,pd,...)
         DBG_ASSERT(dynamic_cast<const CompoundVector*>(GetRawPtr(ip_data->curr()->x())) != NULL);
         c_vec = static_cast<const CompoundVector*>(GetRawPtr(ip_data->curr()->x()));
         SmartPtr<const Vector> nd_only = c_vec->GetComp(3);
         SmartPtr<const Vector> pd_only = c_vec->GetComp(4);
         DBG_ASSERT(IsValid(nd_only));
         DBG_ASSERT(IsValid(pd_only));

         // get scaled d from restonlp
         DBG_ASSERT(dynamic_cast<const CompoundVector*>(GetRawPtr(ip_cq->curr_d())) != NULL);
         c_vec = static_cast<const CompoundVector*>(GetRawPtr(ip_cq->curr_d()));
         SmartPtr<Vector> d_resto = c_vec->GetComp(0)->MakeNewCopy();

         // undo addition of slacks nd-pd
         d_resto->Axpy(-1.0, *nd_only);
         d_resto->Axpy(1.0, *pd_only);

         // violation of d_L <= d(x)
         SmartPtr<Vector> d_viol_L;
         if( orignlp->d_L()->Dim() > 0 )
         {
            d_viol_L = d_resto->MakeNewCopy();
            orignlp->Pd_L()->MultVector(1., *orignlp->d_L(), -1., *d_viol_L);  // d_L - d
            if( !scaled )
               d_viol_L = orignlp->NLP_scaling()->unapply_vector_scaling_d_NonConst(ConstPtr(d_viol_L));
         }
         else
         {
            d_viol_L = d_resto->MakeNew();
            d_viol_L->Set(0.);
         }

         // violation of d(x) <= d_U
         SmartPtr<Vector> d_viol_U;
         if( orignlp->d_U()->Dim() > 0 )
         {
            d_viol_U = d_resto->MakeNewCopy();
            orignlp->Pd_U()->MultVector(-1., *orignlp->d_U(), 1., *d_viol_U);  // d - d_U
            if( !scaled )
               d_viol_U = orignlp->NLP_scaling()->unapply_vector_scaling_d_NonConst(ConstPtr(d_viol_U));
         }
         else
         {
            d_viol_U = d_resto->MakeNew();
            d_viol_U->Set(0.);
         }

         // c(x) = 0, d_L <= d(x) <= d_U should result in complementarities
         // y_c*c(x), (d(x)-d_L)*y_d^-, (d_U-d(x))*y_d^+,  where y_d^+ = max(0,y_d), y_d^- = max(0,-y_d)   (I took the signs from TNLPAdapter::ResortBoundMultipliers)
         // we will merge the latter two into one vector, taking the nonzero entries, i.e., (d(x)-d_L)*y_d^- + (d_U-d(x))*y_d^+
         // to be consistent, it looks like we need to negate for 0 <= c(x) <= 0, since c(x)*y_c^- + (-c(x))*y_c^+ = c(x)*(y_c^- - y_c^+) = -c(x)*y_c
         if( compl_g != NULL )
         {
            DBG_ASSERT(dynamic_cast<const CompoundVector*>(GetRawPtr(ip_data->curr()->y_d())) != NULL);
            c_vec = static_cast<const CompoundVector*>(GetRawPtr(ip_data->curr()->y_d()));
            DBG_ASSERT(c_vec->NComps() == 1);
            SmartPtr<const Vector> yd_only = c_vec->GetComp(0);
            DBG_ASSERT(IsValid(yd_only));

            SmartPtr<Vector> yd_pos = yd_only->MakeNewCopy();
            SmartPtr<Vector> yd_neg = yd_only->MakeNewCopy();
            SmartPtr<Vector> zero = yd_pos->MakeNew();
            zero->Set(0.);
            yd_pos->ElementWiseMax(*zero);
            yd_neg->ElementWiseMin(*zero);  // -y_d^-

            yd_pos->ElementWiseMultiply(*d_viol_U);  // (d(x)-d_U) * y_d^+
            yd_neg->ElementWiseMultiply(*d_viol_L);  // (d_L-d(x)) * (-y_d^-)

            yd_neg->Axpy(-1.0, *yd_pos); // (d(x)-d_L) * y_d^- + (d_U-d(x)) * y_d^+

            c_compl->Scal(-1.0);  // -c(x)*y_c

            tnlp_adapter->ResortG(*c_compl, *yd_neg, compl_g);
         }

         if( nlp_constraint_violation != NULL )
         {
            c->ElementWiseAbs();  // |c(x)|

            // violation of d_L <= d(x) <= d_U:   d_viol_L := max(d_viol_L, d_viol_U, 0)
            d_viol_L->ElementWiseMax(*d_viol_U);
            SmartPtr<Vector> tmp = d_viol_L->MakeNew();
            tmp->Set(0.);
            d_viol_L->ElementWiseMax(*tmp);

            tnlp_adapter->ResortG(*c, *d_viol_L, nlp_constraint_violation);
         }
      }

      return true;
   }

   return false;
}


} // namespace Ipopt

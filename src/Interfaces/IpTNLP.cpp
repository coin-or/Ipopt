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

static
SmartPtr<const DenseVector> curr_x(
   const IpoptData*           ip_data,
   IpoptCalculatedQuantities* /* ip_cq */,
   OrigIpoptNLP*              orignlp,
   RestoIpoptNLP*             restonlp,
   bool                       scaled
   )
{
   SmartPtr<const Vector> x;

   if( restonlp == NULL )
   {
      x = ip_data->curr()->x();
   }
   else
   {
      // get x from the compound vector (x,p,n)
      DBG_ASSERT(dynamic_cast<const CompoundVector*>(GetRawPtr(ip_data->curr()->x())) != NULL);
      const CompoundVector* c_vec = static_cast<const CompoundVector*>(GetRawPtr(ip_data->curr()->x()));
      x = c_vec->GetComp(0);
   }
   DBG_ASSERT(IsValid(x));

   if( !scaled && orignlp->NLP_scaling()->have_x_scaling() )
   {
      x = orignlp->NLP_scaling()->unapply_vector_scaling_x(x);
   }

   DBG_ASSERT(dynamic_cast<const DenseVector*>(GetRawPtr(x)) != NULL);
   return static_cast<const DenseVector*>(GetRawPtr(x));
}

static
SmartPtr<const DenseVector> curr_z_L(
   const IpoptData*           ip_data,
   IpoptCalculatedQuantities* ip_cq,
   OrigIpoptNLP*              orignlp,
   RestoIpoptNLP*             restonlp,
   bool                       scaled
   )
{
   SmartPtr<const Vector> z_L;

   if( restonlp == NULL )
   {
      z_L = ip_data->curr()->z_L();
   }
   else
   {
      DBG_ASSERT(dynamic_cast<const CompoundVector*>(GetRawPtr(ip_data->curr()->z_L())) != NULL);
      const CompoundVector* c_vec = static_cast<const CompoundVector*>(GetRawPtr(ip_data->curr()->z_L()));
      z_L = c_vec->GetComp(0);
   }
   DBG_ASSERT(IsValid(z_L));

   if( !scaled )
   {
      Number obj_unscale_factor = orignlp->NLP_scaling()->unapply_obj_scaling(1.);
      if( orignlp->NLP_scaling()->have_x_scaling() )
      {
         // get copy with x scaling unapplied
         int x_dim = curr_x(ip_data, ip_cq, orignlp, restonlp, true)->Dim();
         SmartPtr<Vector> tmp = orignlp->NLP_scaling()->apply_vector_scaling_x_LU_NonConst(*orignlp->Px_L(), z_L, *new DenseVectorSpace(x_dim));
         // unapply obj scaling
         tmp->Scal(obj_unscale_factor);
         z_L = ConstPtr(tmp);
      }
      else if( obj_unscale_factor != 1. )
      {
         // make copy and unapply obj scaling
         SmartPtr<Vector> tmp = z_L->MakeNewCopy();
         tmp->Scal(obj_unscale_factor);
         z_L = ConstPtr(tmp);
      }
   }

   DBG_ASSERT(dynamic_cast<const DenseVector*>(GetRawPtr(z_L)) != NULL);
   return static_cast<const DenseVector*>(GetRawPtr(z_L));
}

static
SmartPtr<const DenseVector> curr_z_U(
   const IpoptData*           ip_data,
   IpoptCalculatedQuantities* ip_cq,
   OrigIpoptNLP*              orignlp,
   RestoIpoptNLP*             restonlp,
   bool                       scaled
   )
{
   SmartPtr<const Vector> z_U;

   if( restonlp == NULL )
   {
      z_U = ip_data->curr()->z_U();
   }
   else
   {
      DBG_ASSERT(dynamic_cast<const CompoundVector*>(GetRawPtr(ip_data->curr()->z_U())) != NULL);
      const CompoundVector* c_vec = static_cast<const CompoundVector*>(GetRawPtr(ip_data->curr()->z_U()));
      z_U = c_vec->GetComp(0);
   }
   DBG_ASSERT(IsValid(z_U));

   if( !scaled )
   {
      Number obj_unscale_factor = orignlp->NLP_scaling()->unapply_obj_scaling(1.);
      if( orignlp->NLP_scaling()->have_x_scaling() )
      {
         // get copy with x scaling unapplied
         int x_dim = curr_x(ip_data, ip_cq, orignlp, restonlp, true)->Dim();
         SmartPtr<Vector> tmp = orignlp->NLP_scaling()->apply_vector_scaling_x_LU_NonConst(*orignlp->Px_U(), z_U, *new DenseVectorSpace(x_dim));
         // unapply obj scaling
         tmp->Scal(obj_unscale_factor);
         z_U = ConstPtr(tmp);
      }
      else if( obj_unscale_factor != 1. )
      {
         // make copy and unapply obj scaling
         SmartPtr<Vector> tmp = z_U->MakeNewCopy();
         tmp->Scal(obj_unscale_factor);
         z_U = ConstPtr(tmp);
      }
   }

   DBG_ASSERT(dynamic_cast<const DenseVector*>(GetRawPtr(z_U)) != NULL);
   return static_cast<const DenseVector*>(GetRawPtr(z_U));
}

static
SmartPtr<const DenseVector> curr_grad_lag_x(
   const IpoptData*           ip_data,
   IpoptCalculatedQuantities* ip_cq,
   OrigIpoptNLP*              orignlp,
   RestoIpoptNLP*             restonlp,
   bool                       scaled
   )
{
   SmartPtr<const Vector> grad;

   if( restonlp == NULL )
   {
      if( scaled )
      {
         grad = ip_cq->curr_grad_lag_x();
      }
      else
      {
         // adapted from IpoptCalculatedQuantities::unscaled_curr_dual_infeasibility()
         // TODO this doesn't look as if it does something correct if the NLP is scaled
         grad = orignlp->NLP_scaling()->unapply_grad_obj_scaling(ip_cq->curr_grad_lag_x());
      }
   }
   else
   {
      // this looks like it would require reevaluating to get the gradient for the Lagrangian in the original NLP
      // we don't do this, but just return a 0-vector
      // get x from the compound vector (x,p,n)
      DBG_ASSERT(dynamic_cast<const CompoundVector*>(GetRawPtr(ip_data->curr()->x())) != NULL);
      const CompoundVector* c_vec = static_cast<const CompoundVector*>(GetRawPtr(ip_data->curr()->x()));
      SmartPtr<Vector> zero = c_vec->GetComp(0)->MakeNew();
      zero->Set(0.);
      grad = zero;
   }
   DBG_ASSERT(IsValid(grad));

   DBG_ASSERT(dynamic_cast<const DenseVector*>(GetRawPtr(grad)) != NULL);
   return static_cast<const DenseVector*>(GetRawPtr(grad));
}

static
SmartPtr<const DenseVector> curr_c(
   const IpoptData*           ip_data,
   IpoptCalculatedQuantities* ip_cq,
   OrigIpoptNLP*              orignlp,
   RestoIpoptNLP*             restonlp,
   bool                       scaled
   )
{
   SmartPtr<const Vector> c;

   if( restonlp == NULL )
   {
      c = ip_cq->curr_c();
   }
   else
   {
      // get nc, pc from the compound vector (x,nc,pc,nd,pd,...)
      DBG_ASSERT(dynamic_cast<const CompoundVector*>(GetRawPtr(ip_data->curr()->x())) != NULL);
      const CompoundVector* c_vec = static_cast<const CompoundVector*>(GetRawPtr(ip_data->curr()->x()));
      SmartPtr<const Vector> nc_only = c_vec->GetComp(1);
      SmartPtr<const Vector> pc_only = c_vec->GetComp(2);
      DBG_ASSERT(IsValid(nc_only));
      DBG_ASSERT(IsValid(pc_only));

      // get scaled c from restonlp
      DBG_ASSERT(dynamic_cast<const CompoundVector*>(GetRawPtr(ip_cq->curr_c())) != NULL);
      c_vec = static_cast<const CompoundVector*>(GetRawPtr(ip_cq->curr_c()));
      SmartPtr<Vector> c_resto = c_vec->GetComp(0)->MakeNewCopy();

      // undo addition of slacks nc-pc
      c_resto->Axpy(-1.0, *nc_only);
      c_resto->Axpy(1.0, *pc_only);

      c = c_resto;
   }
   DBG_ASSERT(IsValid(c));

   if( !scaled )
      c = orignlp->NLP_scaling()->unapply_vector_scaling_c(c);

   DBG_ASSERT(dynamic_cast<const DenseVector*>(GetRawPtr(c)) != NULL);
   return static_cast<const DenseVector*>(GetRawPtr(c));
}

static
SmartPtr<const DenseVector> curr_y_c(
   const IpoptData*           ip_data,
   IpoptCalculatedQuantities* /* ip_cq */,
   OrigIpoptNLP*              orignlp,
   RestoIpoptNLP*             restonlp,
   bool                       scaled
   )
{
   SmartPtr<const Vector> y_c;

   if( restonlp == NULL )
   {
      y_c = ip_data->curr()->y_c();
   }
   else
   {
      DBG_ASSERT(dynamic_cast<const CompoundVector*>(GetRawPtr(ip_data->curr()->y_c())) != NULL);
      const CompoundVector* c_vec = static_cast<const CompoundVector*>(GetRawPtr(ip_data->curr()->y_c()));
      DBG_ASSERT(c_vec->NComps() == 1);
      y_c = c_vec->GetComp(0);
   }
   DBG_ASSERT(IsValid(y_c));

   if( !scaled )
   {
      Number obj_unscale_factor = orignlp->NLP_scaling()->unapply_obj_scaling(1.);

      if( orignlp->NLP_scaling()->have_c_scaling() )
      {
         SmartPtr<Vector> tmp = orignlp->NLP_scaling()->apply_vector_scaling_c_NonConst(y_c);
         tmp->Scal(obj_unscale_factor);
         y_c = ConstPtr(tmp);
      }
      else if( obj_unscale_factor != 1. )
      {
         // make copy and unapply obj scaling
         SmartPtr<Vector> tmp = y_c->MakeNewCopy();
         tmp->Scal(obj_unscale_factor);
         y_c = ConstPtr(tmp);
      }
   }

   DBG_ASSERT(dynamic_cast<const DenseVector*>(GetRawPtr(y_c)) != NULL);
   return static_cast<const DenseVector*>(GetRawPtr(y_c));
}

static
SmartPtr<const DenseVector> curr_d(
   const IpoptData*           ip_data,
   IpoptCalculatedQuantities* ip_cq,
   OrigIpoptNLP*              orignlp,
   RestoIpoptNLP*             restonlp,
   bool                       scaled
   )
{
   SmartPtr<const Vector> d;

   if( restonlp == NULL )
   {
      d = ip_cq->curr_d();
   }
   else
   {
      // get nd, pd from the compound vector (x,nc,pc,nd,pd,...)
      DBG_ASSERT(dynamic_cast<const CompoundVector*>(GetRawPtr(ip_data->curr()->x())) != NULL);
      const CompoundVector* c_vec = static_cast<const CompoundVector*>(GetRawPtr(ip_data->curr()->x()));
      SmartPtr<const Vector> nd_only = c_vec->GetComp(3);
      SmartPtr<const Vector> pd_only = c_vec->GetComp(4);
      DBG_ASSERT(IsValid(nd_only));
      DBG_ASSERT(IsValid(pd_only));

      // get scaled d from restonlp
      DBG_ASSERT(dynamic_cast<const CompoundVector*>(GetRawPtr(ip_cq->curr_d())) != NULL);
      c_vec = static_cast<const CompoundVector*>(GetRawPtr(ip_cq->curr_d()));
      SmartPtr<Vector> d_resto = c_vec->GetComp(0)->MakeNewCopy();

      // undo addition of slacks nc-pc
      d_resto->Axpy(-1.0, *nd_only);
      d_resto->Axpy(1.0, *pd_only);

      d = d_resto;
   }
   DBG_ASSERT(IsValid(d));

   if( !scaled )
      d = orignlp->NLP_scaling()->unapply_vector_scaling_d(d);

   DBG_ASSERT(dynamic_cast<const DenseVector*>(GetRawPtr(d)) != NULL);
   return static_cast<const DenseVector*>(GetRawPtr(d));
}

static
SmartPtr<const DenseVector> curr_y_d(
   const IpoptData*           ip_data,
   IpoptCalculatedQuantities* /* ip_cq */,
   OrigIpoptNLP*              orignlp,
   RestoIpoptNLP*             restonlp,
   bool                       scaled
   )
{
   SmartPtr<const Vector> y_d;

   if( restonlp == NULL )
   {
      y_d = ip_data->curr()->y_d();
   }
   else
   {
      DBG_ASSERT(dynamic_cast<const CompoundVector*>(GetRawPtr(ip_data->curr()->y_d())) != NULL);
      const CompoundVector* c_vec = static_cast<const CompoundVector*>(GetRawPtr(ip_data->curr()->y_d()));
      DBG_ASSERT(c_vec->NComps() == 1);
      y_d = c_vec->GetComp(0);
   }
   DBG_ASSERT(IsValid(y_d));

   if( !scaled )
   {
      Number obj_unscale_factor = orignlp->NLP_scaling()->unapply_obj_scaling(1.);

      if( orignlp->NLP_scaling()->have_d_scaling() )
      {
         SmartPtr<Vector> tmp = orignlp->NLP_scaling()->apply_vector_scaling_d_NonConst(y_d);
         tmp->Scal(obj_unscale_factor);
         y_d = ConstPtr(tmp);
      }
      else if( obj_unscale_factor != 1. )
      {
         // make copy and unapply obj scaling
         SmartPtr<Vector> tmp = y_d->MakeNewCopy();
         tmp->Scal(obj_unscale_factor);
         y_d = ConstPtr(tmp);
      }
   }

   DBG_ASSERT(dynamic_cast<const DenseVector*>(GetRawPtr(y_d)) != NULL);
   return static_cast<const DenseVector*>(GetRawPtr(y_d));
}

static
SmartPtr<const DenseVector> curr_compl_x_L(
   const IpoptData*           ip_data,
   IpoptCalculatedQuantities* ip_cq,
   OrigIpoptNLP*              orignlp,
   RestoIpoptNLP*             restonlp,
   bool                       scaled
   )
{
   SmartPtr<const Vector> compl_x_L;

   Number obj_unscal = scaled ? 1.0 : orignlp->NLP_scaling()->unapply_obj_scaling(1.0);

   if( restonlp == NULL )
   {
      compl_x_L = ip_cq->curr_compl_x_L();
      if( obj_unscal != 1.0 )
      {
         SmartPtr<Vector> tmp = compl_x_L->MakeNewCopy();
         tmp->Scal(obj_unscal);
         compl_x_L = tmp;
      }
   }
   else
   {
      // get duals z_L for x from the compound vector for the duals of (x,nc,pc,nd,pd,...)
      DBG_ASSERT(dynamic_cast<const CompoundVector*>(GetRawPtr(ip_data->curr()->z_L())) != NULL);
      const CompoundVector* c_vec = static_cast<const CompoundVector*>(GetRawPtr(ip_data->curr()->z_L()));
      SmartPtr<const Vector> z_L_only = c_vec->GetComp(0);
      DBG_ASSERT(IsValid(z_L_only));

      // get slacks for x w.r.t. x_L
      DBG_ASSERT(dynamic_cast<const CompoundVector*>(GetRawPtr(ip_cq->curr_slack_x_L())) != NULL);
      c_vec = static_cast<const CompoundVector*>(GetRawPtr(ip_cq->curr_slack_x_L()));
      SmartPtr<const Vector> slack_x_L = c_vec->GetComp(0);
      DBG_ASSERT(IsValid(slack_x_L));

      // calculate complementarity for x_L
      SmartPtr<Vector> compl_x_L_v = slack_x_L->MakeNewCopy();
      compl_x_L_v->ElementWiseMultiply(*z_L_only);

      // unscale, if desired
      compl_x_L_v->Scal(obj_unscal);

      compl_x_L = compl_x_L_v;
   }
   DBG_ASSERT(IsValid(compl_x_L));

   DBG_ASSERT(dynamic_cast<const DenseVector*>(GetRawPtr(compl_x_L)) != NULL);
   return static_cast<const DenseVector*>(GetRawPtr(compl_x_L));
}

static
SmartPtr<const DenseVector> curr_compl_x_U(
   const IpoptData*           ip_data,
   IpoptCalculatedQuantities* ip_cq,
   OrigIpoptNLP*              orignlp,
   RestoIpoptNLP*             restonlp,
   bool                       scaled
   )
{
   SmartPtr<const Vector> compl_x_U;

   Number obj_unscal = scaled ? 1.0 : orignlp->NLP_scaling()->unapply_obj_scaling(1.0);

   if( restonlp == NULL )
   {
      compl_x_U = ip_cq->curr_compl_x_U();
      if( obj_unscal != 1.0 )
      {
         SmartPtr<Vector> tmp = compl_x_U->MakeNewCopy();
         tmp->Scal(obj_unscal);
         compl_x_U = tmp;
      }
   }
   else
   {
      // get duals z_U for x from the compound vector for the duals of (x,nc,pc,nd,pd,...)
      DBG_ASSERT(dynamic_cast<const CompoundVector*>(GetRawPtr(ip_data->curr()->z_U())) != NULL);
      const CompoundVector* c_vec = static_cast<const CompoundVector*>(GetRawPtr(ip_data->curr()->z_U()));
      SmartPtr<const Vector> z_U_only = c_vec->GetComp(0);
      DBG_ASSERT(IsValid(z_U_only));

      // get slacks for x w.r.t. x_U
      DBG_ASSERT(dynamic_cast<const CompoundVector*>(GetRawPtr(ip_cq->curr_slack_x_U())) != NULL);
      c_vec = static_cast<const CompoundVector*>(GetRawPtr(ip_cq->curr_slack_x_U()));
      SmartPtr<const Vector> slack_x_U = c_vec->GetComp(0);
      DBG_ASSERT(IsValid(slack_x_U));

      // calculate complementarity for x_U
      SmartPtr<Vector> compl_x_U_v = slack_x_U->MakeNewCopy();
      compl_x_U_v->ElementWiseMultiply(*z_U_only);

      // unscale, if desired
      compl_x_U_v->Scal(obj_unscal);

      compl_x_U = compl_x_U_v;
   }
   DBG_ASSERT(IsValid(compl_x_U));

   DBG_ASSERT(dynamic_cast<const DenseVector*>(GetRawPtr(compl_x_U)) != NULL);
   return static_cast<const DenseVector*>(GetRawPtr(compl_x_U));
}

bool TNLP::get_curr_iterate(
   const IpoptData*           ip_data,
   IpoptCalculatedQuantities* ip_cq,
   bool                       scaled,
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
   Ipopt::RestoIpoptNLP* restonlp = NULL;
   TNLPAdapter* tnlp_adapter;
   Index n_full;
   Index m_full;

   // check whether we use a OrigIpoptNLP
   orignlp = dynamic_cast<OrigIpoptNLP*>(GetRawPtr(ip_cq->GetIpoptNLP()));
   if( orignlp == NULL )
   {
      // check whether we are in restoration phase, so get a RestoIpoptNLP
      restonlp = dynamic_cast<RestoIpoptNLP*>(GetRawPtr(ip_cq->GetIpoptNLP()));

      if( restonlp == NULL )
         // if it is neither OrigIpoptNLP nor RestoIpoptNLP, then we don't know how to retrieve x
         return false;

      if( (orignlp = dynamic_cast<OrigIpoptNLP*>(&restonlp->OrigIpNLP())) == NULL )
         return false;
   }

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
      tnlp_adapter->ResortX(*curr_x(ip_data, ip_cq, orignlp, restonlp, scaled), x);

   // resort Ipopt-internal variable duals to TNLP-version
   if( z_L != NULL || z_U != NULL )
      tnlp_adapter->ResortBoundMultipliers(
         *curr_y_c(ip_data, ip_cq, orignlp, restonlp, scaled),
         *curr_z_L(ip_data, ip_cq, orignlp, restonlp, scaled), z_L,
         *curr_z_U(ip_data, ip_cq, orignlp, restonlp, scaled), z_U);

   // resort Ipopt-interval constraint activity to TNLP-version
   if( g != NULL )
   {
      if( !scaled || !orignlp->NLP_scaling()->have_c_scaling() )
      {
         tnlp_adapter->ResortG(*curr_c(ip_data, ip_cq, orignlp, restonlp, scaled), *curr_d(ip_data, ip_cq, orignlp, restonlp, scaled), g, true);
      }
      else
      {
         // scaled: add c(x) + c_rhs here, so we can scale c_rhs first
         SmartPtr<const DenseVector> c_scaled = curr_c(ip_data, ip_cq, orignlp, restonlp, true);

         SmartPtr<DenseVector> c_rhs = new DenseVector(new DenseVectorSpace(c_scaled->Dim()));
         c_rhs->SetValues(tnlp_adapter->GetC_Rhs());
         SmartPtr<Vector> c_rhs_scaled = orignlp->NLP_scaling()->apply_vector_scaling_c_NonConst(c_rhs);

         c_rhs_scaled->Axpy(1.0, *c_scaled);  // c(x) + c_rhs = g(x)  (scaled)

         tnlp_adapter->ResortG(*c_rhs_scaled, *curr_d(ip_data, ip_cq, orignlp, restonlp, scaled), g);
      }
   }

   // resort Ipopt-internal constraint duals to TNLP-version
   if( lambda != NULL )
      tnlp_adapter->ResortG(*curr_y_c(ip_data, ip_cq, orignlp, restonlp, scaled), *curr_y_d(ip_data, ip_cq, orignlp, restonlp, scaled), lambda);

   return true;
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
   Ipopt::RestoIpoptNLP* restonlp = NULL;
   TNLPAdapter* tnlp_adapter;
   Index n_full;
   Index m_full;

   // check whether we use a OrigIpoptNLP
   orignlp = dynamic_cast<OrigIpoptNLP*>(GetRawPtr(ip_cq->GetIpoptNLP()));
   if( orignlp == NULL )
   {
      // check whether we are in restoration phase, so get a RestoIpoptNLP
      restonlp = dynamic_cast<RestoIpoptNLP*>(GetRawPtr(ip_cq->GetIpoptNLP()));

      if( restonlp == NULL )
         // if it is neither OrigIpoptNLP nor RestoIpoptNLP, then we don't know how to retrieve x
         return false;

      if( (orignlp = dynamic_cast<OrigIpoptNLP*>(&restonlp->OrigIpNLP())) == NULL )
         return false;
   }

   tnlp_adapter = dynamic_cast<TNLPAdapter*>(GetRawPtr(orignlp->nlp()));
   if( tnlp_adapter == NULL )
      return false;

   tnlp_adapter->GetFullDimensions(n_full, m_full);
   if( n != n_full && (compl_x_L != NULL || compl_x_U != NULL || grad_lag_x != NULL) )
      THROW_EXCEPTION(IpoptException, "Incorrect dimension of x given to TNLP::get_curr_violations().\n");
   if( m != m_full && (nlp_constraint_violation != NULL || compl_g != NULL) )
      THROW_EXCEPTION(IpoptException, "Incorrect dimension of g(x) given to TNLP::get_curr_violations().\n");

   int n_x_fixed;
   Index* x_fixed_map;
   TNLPAdapter::FixedVariableTreatmentEnum fixed_variable_treatment;
   tnlp_adapter->GetFixedVariables(n_x_fixed, x_fixed_map, fixed_variable_treatment);

   if( compl_x_L != NULL || compl_x_U != NULL )
   {
      // this should give XZe from (5)
      tnlp_adapter->ResortBnds(*curr_compl_x_L(ip_data, ip_cq, orignlp, restonlp, scaled), compl_x_L,
         *curr_compl_x_U(ip_data, ip_cq, orignlp, restonlp, scaled), compl_x_U);

      if( n_x_fixed > 0 && fixed_variable_treatment == TNLPAdapter::MAKE_CONSTRAINT )
      {
         // compl_x_L = z_L * c(x) = y_c^- * c(x)
         // compl_x_U = z_U * (-c(x)) = -y_c^+ * c(x)

         SmartPtr<const Vector> y_c = curr_y_c(ip_data, ip_cq, orignlp, restonlp, scaled);
         SmartPtr<Vector> yc_pos = y_c->MakeNewCopy();
         SmartPtr<Vector> yc_neg = y_c->MakeNewCopy();
         SmartPtr<Vector> zero = yc_pos->MakeNew();
         zero->Set(0.);
         yc_pos->ElementWiseMax(*zero);  //  y_c^+
         yc_neg->ElementWiseMin(*zero);  // -y_c^-

         SmartPtr<const Vector> c = curr_c(ip_data, ip_cq, orignlp, restonlp, scaled);
         yc_pos->ElementWiseMultiply(*c);  // y_c^+ * c(x)
         yc_neg->ElementWiseMultiply(*c);  // -y_c^- * c(x)

         DBG_ASSERT(dynamic_cast<DenseVector*>(GetRawPtr(yc_pos)) != NULL);
         DBG_ASSERT(dynamic_cast<DenseVector*>(GetRawPtr(yc_neg)) != NULL);
         Number* yc_pos_val = static_cast<DenseVector*>(GetRawPtr(yc_pos))->ExpandedValues();
         Number* yc_neg_val = static_cast<DenseVector*>(GetRawPtr(yc_neg))->ExpandedValues();
         for( int i = 0; i < n_x_fixed; ++i )
         {
            compl_x_L[x_fixed_map[i]] = -yc_neg_val[yc_neg->Dim()-n_x_fixed+i];
            compl_x_U[x_fixed_map[i]] = -yc_pos_val[yc_pos->Dim()-n_x_fixed+i];
         }
      }
   }

   if( grad_lag_x != NULL )
   {
      // this will set the derivative of the Lagrangian w.r.t. fixed variables to 0 (for any fixed_variables_treatment)
      // the actual values are nowhere stored within Ipopt Data or CQ, since TNLPAdapter does not seem to pass them on
      // we may have to reevaluate and capture them in TNLPAdapater
      tnlp_adapter->ResortX(*curr_grad_lag_x(ip_data, ip_cq, orignlp, restonlp, scaled), grad_lag_x, false);

      // if fixed_variable_treatment is make_constraint, then fixed variable contribute y_c*x to the Lagrangian
      // so we should subtract y_c for these entries; FIXME deal with scaling
      // but if we are in restoration phase, then curr_grad_lag_x() just returns a zero vector, so don't change this
      if( restonlp == NULL && n_x_fixed > 0 && fixed_variable_treatment == TNLPAdapter::MAKE_CONSTRAINT )
      {
         SmartPtr<const DenseVector> y_c = curr_y_c(ip_data, ip_cq, orignlp, restonlp, scaled);
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

      SmartPtr<const DenseVector> c = curr_c(ip_data, ip_cq, orignlp, restonlp, scaled);
      SmartPtr<const DenseVector> d = curr_d(ip_data, ip_cq, orignlp, restonlp, true);

      // violation of d_L <= d(x)
      SmartPtr<Vector> d_viol_L;
      if( orignlp->d_L()->Dim() > 0 )
      {
         d_viol_L = d->MakeNewCopy();
         orignlp->Pd_L()->MultVector(1., *orignlp->d_L(), -1., *d_viol_L);  // d_L - d, scaled
         if( !scaled && orignlp->NLP_scaling()->have_d_scaling() )
            d_viol_L = orignlp->NLP_scaling()->unapply_vector_scaling_d_NonConst(ConstPtr(d_viol_L));
      }
      else
      {
         d_viol_L = d->MakeNew();
         d_viol_L->Set(0.);
      }

      // violation of d(x) <= d_U
      SmartPtr<Vector> d_viol_U;
      if( orignlp->d_U()->Dim() > 0 )
      {
         d_viol_U = d->MakeNewCopy();
         orignlp->Pd_U()->MultVector(-1., *orignlp->d_U(), 1., *d_viol_U);  // d - d_U, scaled
         if( !scaled && orignlp->NLP_scaling()->have_d_scaling() )
            d_viol_U = orignlp->NLP_scaling()->unapply_vector_scaling_d_NonConst(ConstPtr(d_viol_U));
      }
      else
      {
         d_viol_U = d->MakeNew();
         d_viol_U->Set(0.);
      }

      // c(x) = 0, d_L <= d(x) <= d_U should result in complementarities
      // y_c*c(x), (d(x)-d_L)*y_d^-, (d_U-d(x))*y_d^+,  where y_d^+ = max(0,y_d), y_d^- = max(0,-y_d)   (I took the signs from TNLPAdapter::ResortBoundMultipliers)
      // we will merge the latter two into one vector, taking the nonzero entries, i.e., (d(x)-d_L)*y_d^- + (d_U-d(x))*y_d^+
      // to be consistent, it looks like we need to negate for 0 <= c(x) <= 0, since c(x)*y_c^- + (-c(x))*y_c^+ = c(x)*(y_c^- - y_c^+) = -c(x)*y_c
      if( compl_g != NULL )
      {
         SmartPtr<const DenseVector> y_d = curr_y_d(ip_data, ip_cq, orignlp, restonlp, scaled);

         SmartPtr<Vector> yd_pos = y_d->MakeNewCopy();
         SmartPtr<Vector> yd_neg = y_d->MakeNewCopy();
         SmartPtr<Vector> zero = yd_pos->MakeNew();
         zero->Set(0.);
         yd_pos->ElementWiseMax(*zero);
         yd_neg->ElementWiseMin(*zero);  // -y_d^-

         yd_pos->ElementWiseMultiply(*d_viol_U);  // (d(x)-d_U) * y_d^+
         yd_neg->ElementWiseMultiply(*d_viol_L);  // (d_L-d(x)) * (-y_d^-)

         yd_neg->Axpy(-1.0, *yd_pos); // (d(x)-d_L) * y_d^- + (d_U-d(x)) * y_d^+

         SmartPtr<const DenseVector> y_c = curr_y_c(ip_data, ip_cq, orignlp, restonlp, scaled);

         SmartPtr<Vector> c_compl = c->MakeNewCopy();
         c_compl->ElementWiseMultiply(*y_c);      // c(x)*y_c
         c_compl->Scal(-1.0);  // -c(x)*y_c

         tnlp_adapter->ResortG(*c_compl, *yd_neg, compl_g);
      }

      if( nlp_constraint_violation != NULL )
      {
         // violation of c(x) = 0
         SmartPtr<Vector> c_viol = c->MakeNewCopy();
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

} // namespace Ipopt

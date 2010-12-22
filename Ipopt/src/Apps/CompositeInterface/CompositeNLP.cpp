// Copyright (C) 2004, International Business Machines and others.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// $Id$
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#include "CompositeNLP.hpp"
#include "IpCompoundVector.hpp"
#include "IpCompoundMatrix.hpp"
#include "IpCompoundSymMatrix.hpp"
#include "IpZeroMatrix.hpp"
#include "IpMa27TSolverInterface.hpp"
#include "IpBDSymLinearSolver.hpp"
#include "IpTSymLinearSolver.hpp"
#include "IpMPSymLinearSolver.hpp"

namespace Ipopt
{

  CompositeNLP::CompositeNLP(std::vector<SmartPtr<NLP> > nlps, SmartPtr<VectorSpace> q_space,
                             std::vector<SmartPtr<VectorSpace> > linking_eqn_c_spaces,
                             std::vector<SmartPtr<Matrix> > Jx_linking_eqns,
                             std::vector<SmartPtr<Matrix> > Jq_linking_eqns)
      :
      nlps_(nlps),
      q_space_(q_space),
      linking_eqn_c_spaces_(linking_eqn_c_spaces),
      Jx_linking_eqns_(Jx_linking_eqns),
      Jq_linking_eqns_(Jq_linking_eqns)
  {
    DBG_ASSERT(nlps.size() == linking_eqn_c_spaces.size());
    DBG_ASSERT(nlps.size() == Jx_linking_eqns.size());
    DBG_ASSERT(nlps.size() == Jq_linking_eqns.size());
  }

  CompositeNLP::~CompositeNLP()
  {}

  bool CompositeNLP::ProcessOptions(const OptionsList& options,
                                    const std::string& prefix)
  {
    bool ret = true;
    // need to pass the processing down to the other nlps
    for (Index i=0; i<(Index)nlps_.size(); i++) {
      if (!nlps_[i]->ProcessOptions(options, prefix)) {
	ret = false;
      }
    }
    return ret;
  }

  bool CompositeNLP::GetSpaces(SmartPtr<const VectorSpace>& x_space,
                               SmartPtr<const VectorSpace>& c_space,
                               SmartPtr<const VectorSpace>& d_space,
                               SmartPtr<const VectorSpace>& x_l_space,
                               SmartPtr<const MatrixSpace>& px_l_space,
                               SmartPtr<const VectorSpace>& x_u_space,
                               SmartPtr<const MatrixSpace>& px_u_space,
                               SmartPtr<const VectorSpace>& d_l_space,
                               SmartPtr<const MatrixSpace>& pd_l_space,
                               SmartPtr<const VectorSpace>& d_u_space,
                               SmartPtr<const MatrixSpace>& pd_u_space,
                               SmartPtr<const MatrixSpace>& Jac_c_space,
                               SmartPtr<const MatrixSpace>& Jac_d_space,
                               SmartPtr<const SymMatrixSpace>& Hess_lagrangian_space)
  {
    Index n_nlps = nlps_.size();
    DBG_ASSERT(n_nlps > 0);
    Index q_dim = Jq_linking_eqns_[0]->NCols();
    Index linking_eqn_dim = 0;

    Index x_dim=0;
    Index c_dim=0;
    Index d_dim=0;
    Index x_l_dim = 0;
    Index px_l_cols = 0;
    Index px_l_rows = 0;
    Index x_u_dim = 0;
    Index px_u_cols = 0;
    Index px_u_rows = 0;
    Index d_l_dim = 0;
    Index pd_l_cols = 0;
    Index pd_l_rows = 0;
    Index d_u_dim = 0;
    Index pd_u_cols = 0;
    Index pd_u_rows = 0;
    Index jac_c_cols = 0;
    Index jac_c_rows = 0;
    Index jac_d_cols = 0;
    Index jac_d_rows = 0;
    Index h_cols = 0;
    Index h_rows = 0;

    // retrieve the necessary spaces from the individual NLPS
    std::vector<SmartPtr<const VectorSpace> > x_spaces;
    std::vector<SmartPtr<const VectorSpace> > c_spaces;
    std::vector<SmartPtr<const VectorSpace> > d_spaces;
    std::vector<SmartPtr<const VectorSpace> > x_l_spaces;
    std::vector<SmartPtr<const MatrixSpace> > px_l_spaces;
    std::vector<SmartPtr<const VectorSpace> > x_u_spaces;
    std::vector<SmartPtr<const MatrixSpace> > px_u_spaces;
    std::vector<SmartPtr<const VectorSpace> > d_l_spaces;
    std::vector<SmartPtr<const MatrixSpace> > pd_l_spaces;
    std::vector<SmartPtr<const VectorSpace> > d_u_spaces;
    std::vector<SmartPtr<const MatrixSpace> > pd_u_spaces;
    std::vector<SmartPtr<const MatrixSpace> > jac_c_spaces;
    std::vector<SmartPtr<const MatrixSpace> > jac_d_spaces;
    std::vector<SmartPtr<const SymMatrixSpace> > h_spaces;
    for (Index i=0; i<n_nlps; i++) {
      DBG_ASSERT(IsValid(nlps_[i]));
      ASSERT_EXCEPTION(Jq_linking_eqns_[i]->NCols() == q_dim, INVALID_JACOBIAN_DIMENSION_FOR_LINKING_EQUATIONS,
                       "The # of columns in the Jq_linking_eqn must be the same for all nlp periods, i.e. Jc_linking_eqns[i] == q_dim for all i");

      linking_eqn_dim += Jx_linking_eqns_[i]->NRows();
      ASSERT_EXCEPTION(Jx_linking_eqns_[i]->NRows() == Jq_linking_eqns_[i]->NRows() , INVALID_JACOBIAN_DIMENSION_FOR_LINKING_EQUATIONS,
                       "The # of rows in the Jx_linking_eqns[i] must be the same as the number of rows in the Jq_linking_eqns[i]");

      SmartPtr<const VectorSpace> x_space_i;
      SmartPtr<const VectorSpace> c_space_i;
      SmartPtr<const VectorSpace> d_space_i;
      SmartPtr<const VectorSpace> x_l_space_i;
      SmartPtr<const MatrixSpace> px_l_space_i;
      SmartPtr<const VectorSpace> x_u_space_i;
      SmartPtr<const MatrixSpace> px_u_space_i;
      SmartPtr<const VectorSpace> d_l_space_i;
      SmartPtr<const MatrixSpace> pd_l_space_i;
      SmartPtr<const VectorSpace> d_u_space_i;
      SmartPtr<const MatrixSpace> pd_u_space_i;
      SmartPtr<const MatrixSpace> jac_c_space_i;
      SmartPtr<const MatrixSpace> jac_d_space_i;
      SmartPtr<const SymMatrixSpace> h_space_i;

      bool retvalue = nlps_[i]->GetSpaces(x_space_i,
                                          c_space_i,
                                          d_space_i,
                                          x_l_space_i,
                                          px_l_space_i,
                                          x_u_space_i,
                                          px_u_space_i,
                                          d_l_space_i,
                                          pd_l_space_i,
                                          d_u_space_i,
                                          pd_u_space_i,
                                          jac_c_space_i,
                                          jac_d_space_i,
                                          h_space_i);

      if (!retvalue) {
        return false;
      }

      x_spaces.push_back(x_space_i);
      x_dim += x_space_i->Dim();
      ASSERT_EXCEPTION(Jx_linking_eqns_[i]->NCols() == x_space_i->Dim(), INVALID_JACOBIAN_DIMENSION_FOR_LINKING_EQUATIONS,
                       "The # of columns in the Jx_linking_eqn must be the same as the dimension of x for that NLP period.");

      c_spaces.push_back(c_space_i);
      c_dim += c_space_i->Dim();

      d_spaces.push_back(d_space_i);
      d_dim += d_space_i->Dim();

      x_l_spaces.push_back(x_l_space_i);
      x_l_dim += x_l_space_i->Dim();

      px_l_spaces.push_back(px_l_space_i);
      px_l_cols += px_l_space_i->NCols();
      px_l_rows += px_l_space_i->NRows();

      x_u_spaces.push_back(x_u_space_i);
      x_u_dim += x_u_space_i->Dim();

      px_u_spaces.push_back(px_u_space_i);
      px_u_cols += px_u_space_i->NCols();
      px_u_rows += px_u_space_i->NRows();

      d_l_spaces.push_back(d_l_space_i);
      d_l_dim += d_l_space_i->Dim();

      pd_l_spaces.push_back(pd_l_space_i);
      pd_l_cols += pd_l_space_i->NCols();
      pd_l_rows += pd_l_space_i->NRows();

      d_u_spaces.push_back(d_u_space_i);
      d_u_dim += d_u_space_i->Dim();

      pd_u_spaces.push_back(pd_u_space_i);
      pd_u_cols += pd_u_space_i->NCols();
      pd_u_rows += pd_u_space_i->NRows();

      jac_c_spaces.push_back(jac_c_space_i);
      jac_c_cols += jac_c_space_i->NCols();
      jac_c_rows += jac_c_space_i->NRows();

      jac_d_spaces.push_back(jac_d_space_i);
      jac_d_cols += jac_d_space_i->NCols();
      jac_d_rows += jac_d_space_i->NRows();

      h_spaces.push_back(h_space_i);
      h_cols += h_space_i->NCols();
      h_rows += h_space_i->NRows();
    }

    // Create the compound vector spaces for the composite nlps
    // Need space for each nlp + one for the common variables
    SmartPtr<CompoundVectorSpace> Cx_space = new CompoundVectorSpace(n_nlps + 1, x_dim + q_dim);
    SmartPtr<CompoundVectorSpace> Cc_space = new CompoundVectorSpace(2*n_nlps, c_dim + linking_eqn_dim);
    SmartPtr<CompoundVectorSpace> Cd_space = new CompoundVectorSpace(n_nlps, d_dim);
    SmartPtr<CompoundVectorSpace> Cx_l_space = new CompoundVectorSpace(n_nlps, x_l_dim);
    SmartPtr<CompoundMatrixSpace> Cpx_l_space = new CompoundMatrixSpace(n_nlps + 1, n_nlps, px_l_rows + q_dim, px_l_cols);
    SmartPtr<CompoundVectorSpace> Cx_u_space = new CompoundVectorSpace(n_nlps, x_u_dim);
    SmartPtr<CompoundMatrixSpace> Cpx_u_space = new CompoundMatrixSpace(n_nlps + 1, n_nlps, px_u_rows + q_dim, px_u_cols);
    SmartPtr<CompoundVectorSpace> Cd_l_space = new CompoundVectorSpace(n_nlps, d_l_dim);
    SmartPtr<CompoundMatrixSpace> Cpd_l_space = new CompoundMatrixSpace(n_nlps, n_nlps, pd_l_rows, pd_l_cols);
    SmartPtr<CompoundVectorSpace> Cd_u_space = new CompoundVectorSpace(n_nlps, d_u_dim);
    SmartPtr<CompoundMatrixSpace> Cpd_u_space = new CompoundMatrixSpace(n_nlps, n_nlps, pd_u_rows, pd_u_cols);
    DBG_ASSERT(jac_c_cols == x_dim && jac_c_rows == c_dim);
    SmartPtr<CompoundMatrixSpace> CJac_c_space = new CompoundMatrixSpace(2*n_nlps, n_nlps+1, c_dim + linking_eqn_dim, x_dim + q_dim);
    DBG_ASSERT(jac_d_cols == x_dim && jac_d_rows == d_dim);
    SmartPtr<CompoundMatrixSpace> CJac_d_space = new CompoundMatrixSpace(n_nlps, n_nlps+1, d_dim, x_dim + q_dim);
    SmartPtr<CompoundSymMatrixSpace> CHess_lagrangian_space = new CompoundSymMatrixSpace(n_nlps + 1, x_dim + q_dim);

    // Set the compound space dimensions
    for (Index i=0; i<n_nlps; i++) {
      Cpx_l_space->SetBlockRows(i, x_spaces[i]->Dim());
      Cpx_l_space->SetBlockCols(i, x_l_spaces[i]->Dim());
      Cpx_u_space->SetBlockRows(i, x_spaces[i]->Dim());
      Cpx_u_space->SetBlockCols(i, x_u_spaces[i]->Dim());

      Cpd_l_space->SetBlockRows(i, d_spaces[i]->Dim());
      Cpd_l_space->SetBlockCols(i, d_l_spaces[i]->Dim());
      Cpd_u_space->SetBlockRows(i, d_spaces[i]->Dim());
      Cpd_u_space->SetBlockCols(i, d_u_spaces[i]->Dim());

      CJac_c_space->SetBlockRows(i, jac_c_spaces[i]->NRows());
      CJac_c_space->SetBlockCols(i, x_spaces[i]->Dim());
      CJac_c_space->SetBlockRows(n_nlps + i, Jx_linking_eqns_[i]->NRows());

      CJac_d_space->SetBlockRows(i, jac_d_spaces[i]->NRows());
      CJac_d_space->SetBlockCols(i, x_spaces[i]->Dim());

      CHess_lagrangian_space->SetBlockDim(i, x_spaces[i]->Dim());
    }

    Cpx_l_space->SetBlockRows(n_nlps, q_dim);
    Cpx_u_space->SetBlockRows(n_nlps, q_dim);
    CJac_c_space->SetBlockCols(n_nlps, q_dim);
    CJac_d_space->SetBlockCols(n_nlps, q_dim);
    CHess_lagrangian_space->SetBlockDim(n_nlps, q_dim);

    // Set the compound spaces
    x_space = GetRawPtr(Cx_space);
    c_space = GetRawPtr(Cc_space);
    d_space = GetRawPtr(Cd_space);
    x_l_space = GetRawPtr(Cx_l_space);
    px_l_space = GetRawPtr(Cpx_l_space);
    x_u_space = GetRawPtr(Cx_u_space);
    px_u_space = GetRawPtr(Cpx_u_space);
    d_l_space = GetRawPtr(Cd_l_space);
    pd_l_space = GetRawPtr(Cpd_l_space);
    d_u_space = GetRawPtr(Cd_u_space);
    pd_u_space = GetRawPtr(Cpd_u_space);
    Jac_c_space = GetRawPtr(CJac_c_space);
    Jac_d_space = GetRawPtr(CJac_d_space);
    Hess_lagrangian_space = GetRawPtr(CHess_lagrangian_space);


    for (Index i=0; i<n_nlps; i++) {
      // Assign the individual vector spaces to the compound vectors
      Cx_space->SetCompSpace(i, *x_spaces[i]);
      Cc_space->SetCompSpace(i, *c_spaces[i]);
      Cc_space->SetCompSpace(i + n_nlps, *linking_eqn_c_spaces_[i]);
      Cd_space->SetCompSpace(i, *d_spaces[i]);
      Cx_l_space->SetCompSpace(i, *x_l_spaces[i]);
      bool auto_allocate = true;
      Cpx_l_space->SetCompSpace(i, i, *px_l_spaces[i], auto_allocate);
      Cx_u_space->SetCompSpace(i, *x_u_spaces[i]);
      Cpx_u_space->SetCompSpace(i, i, *px_u_spaces[i], auto_allocate);
      Cd_l_space->SetCompSpace(i, *d_l_spaces[i]);
      Cpd_l_space->SetCompSpace(i, i, *pd_l_spaces[i], auto_allocate);
      Cd_u_space->SetCompSpace(i, *d_u_spaces[i]);
      Cpd_u_space->SetCompSpace(i, i, *pd_u_spaces[i], auto_allocate);
      CJac_c_space->SetCompSpace(i, i, *jac_c_spaces[i], true);
      CJac_c_space->SetCompSpace(n_nlps + i, i, *Jx_linking_eqns_[i]->OwnerSpace());
      CJac_c_space->SetCompSpace(n_nlps + i, n_nlps, *Jq_linking_eqns_[i]->OwnerSpace());
      CJac_d_space->SetCompSpace(i, i, *jac_d_spaces[i], true);
      //      CJac_d_space->SetCompSpace(i, n_nlps, *(new ZeroMatrixSpace(jac_d_spaces[i]->NRows(), q_dim)), true);
      CHess_lagrangian_space->SetCompSpace(i, i, *h_spaces[i], true);
      //      h_space->SetCompSpace(n_nlps, i, new ZeroMatrixSpace(q_dim, jac_c_spaces[i]->NCols()));
    }

    Cx_space->SetCompSpace(n_nlps, *q_space_);
    //    Cpx_l_space->SetCompSpace(n_nlps, n_nlps - 1, *(new ZeroMatrixSpace(q_dim, px_l_spaces[n_nlps-1]->NCols())), true);
    //    Cpx_u_space->SetCompSpace(n_nlps, n_nlps - 1, *(new ZeroMatrixSpace(q_dim, px_u_spaces[n_nlps-1]->NCols())), true);
    //    CHess_lagrangian_space->SetCompSpace(n_nlps, n_nlps, *(new ZeroMatrixSpace(q_dim, q_dim)));

    x_space = GetRawPtr(Cx_space);
    c_space = GetRawPtr(Cc_space);
    d_space = GetRawPtr(Cd_space);
    x_l_space = GetRawPtr(Cx_l_space);
    px_l_space = GetRawPtr(Cpx_l_space);
    x_u_space = GetRawPtr(Cx_u_space);
    px_u_space = GetRawPtr(Cpx_u_space);
    d_l_space = GetRawPtr(Cd_l_space);
    pd_l_space = GetRawPtr(Cpd_l_space);
    d_u_space = GetRawPtr(Cd_u_space);
    pd_u_space = GetRawPtr(Cpd_u_space);
    Jac_c_space = GetRawPtr(CJac_c_space);
    Jac_d_space = GetRawPtr(CJac_d_space);
    Hess_lagrangian_space = GetRawPtr(CHess_lagrangian_space);

    return true;
  }

  bool CompositeNLP::GetBoundsInformation(const Matrix& Px_L,
                                          Vector& x_L,
                                          const Matrix& Px_U,
                                          Vector& x_U,
                                          const Matrix& Pd_L,
                                          Vector& d_L,
                                          const Matrix& Pd_U,
                                          Vector& d_U)
  {
    CompoundVector* cx_l = dynamic_cast<CompoundVector*>(&x_L);
    const CompoundMatrix* cpx_l = dynamic_cast<const CompoundMatrix*>(&Px_L);
    CompoundVector* cx_u = dynamic_cast<CompoundVector*>(&x_U);
    const CompoundMatrix* cpx_u= dynamic_cast<const CompoundMatrix*>(&Px_U);

    CompoundVector* cd_l = dynamic_cast<CompoundVector*>(&d_L);
    const CompoundMatrix* cpd_l = dynamic_cast<const CompoundMatrix*>(&Pd_L);
    CompoundVector* cd_u = dynamic_cast<CompoundVector*>(&d_U);
    const CompoundMatrix* cpd_u = dynamic_cast<const CompoundMatrix*>(&Pd_U);

    DBG_ASSERT(cx_l && cpx_l && cx_u && cpx_u && cd_l && cpd_l && cd_u && cpd_u);

    Index n_nlps = nlps_.size();
    for (Index i=0; i<n_nlps; i++) {
      SmartPtr<Vector> x_l = cx_l->GetCompNonConst(i);
      SmartPtr<Vector> x_u = cx_u->GetCompNonConst(i);
      SmartPtr<Vector> d_l = cd_l->GetCompNonConst(i);
      SmartPtr<Vector> d_u = cd_u->GetCompNonConst(i);

      SmartPtr<const Matrix> px_l = cpx_l->GetComp(i, i);
      SmartPtr<const Matrix> px_u = cpx_u->GetComp(i, i);
      SmartPtr<const Matrix> pd_l = cpd_l->GetComp(i, i);
      SmartPtr<const Matrix> pd_u = cpd_u->GetComp(i, i);

      DBG_ASSERT(IsValid(x_l) && IsValid(x_u) && IsValid(d_l) && IsValid(d_u)
                 && IsValid(px_l) && IsValid(px_u) && IsValid(pd_l) && IsValid(pd_u));

      if (!nlps_[i]->GetBoundsInformation(*px_l, *x_l, *px_u, *x_u,
                                          *pd_l, *d_l, *pd_u, *d_u) ) {
        return false;
      }

    }
    return true;
  }

  bool CompositeNLP::GetStartingPoint(
    SmartPtr<Vector> x,
    bool need_x,
    SmartPtr<Vector> y_c,
    bool need_y_c,
    SmartPtr<Vector> y_d,
    bool need_y_d,
    SmartPtr<Vector> z_L,
    bool need_z_L,
    SmartPtr<Vector> z_U,
    bool need_z_U)
  {

    CompoundVector* cx = dynamic_cast<CompoundVector*>(GetRawPtr(x));
    CompoundVector* cy_c = dynamic_cast<CompoundVector*>(GetRawPtr(y_c));
    CompoundVector* cy_d = dynamic_cast<CompoundVector*>(GetRawPtr(y_d));
    CompoundVector* cz_l = dynamic_cast<CompoundVector*>(GetRawPtr(z_L));
    CompoundVector* cz_u = dynamic_cast<CompoundVector*>(GetRawPtr(z_U));

    Index n_nlps = nlps_.size();
    for (Index i=0; i<n_nlps; i++) {
      SmartPtr<Vector> xi = (need_x) ? cx->GetCompNonConst(i) : NULL;
      SmartPtr<Vector> y_ci = (need_y_c) ? cy_c->GetCompNonConst(i) : NULL;
      SmartPtr<Vector> y_di = (need_y_d) ? cy_d->GetCompNonConst(i) : NULL;
      SmartPtr<Vector> z_li = (need_z_L) ? cz_l->GetCompNonConst(i) : NULL;
      SmartPtr<Vector> z_ui = (need_z_U) ? cz_u->GetCompNonConst(i) : NULL;

      if (!nlps_[i]->GetStartingPoint(xi, need_x,
                                      y_ci, need_y_c, y_di, need_y_d,
                                      z_li, need_z_L, z_ui, need_z_U) ) {
        return false;
      }
    }

    // Don't forget to initialize q
    if (need_x) {
      SmartPtr<Vector> q = cx->GetCompNonConst(n_nlps);
      q->Set(0.0);
    }

    return true;
  }

  bool CompositeNLP::Eval_f(const Vector& x, Number& f)
  {
    // evaluate the f of all the nlps and sum them
    f = 0;
    const CompoundVector* cx = dynamic_cast<const CompoundVector*>(&x);
    DBG_ASSERT(cx);
    Index n_nlps = nlps_.size();
    for (Index i=0; i<n_nlps; i++) {
      SmartPtr<const Vector> x_i = cx->GetComp(i);
      DBG_ASSERT(IsValid(x_i));
      Number f_i = 0;
      if (!nlps_[i]->Eval_f(*x_i, f_i)) {
        return false;
      }
      f += f_i;
    }
    return true;
  }

  bool CompositeNLP::Eval_grad_f(const Vector& x, Vector& g_f)
  {
    const CompoundVector* cx = dynamic_cast<const CompoundVector*>(&x);
    DBG_ASSERT(cx);
    CompoundVector* cg_f = dynamic_cast<CompoundVector*>(&g_f);
    DBG_ASSERT(cg_f);

    Index n_nlps = nlps_.size();
    for (Index i=0; i<n_nlps; i++) {
      SmartPtr<const Vector> x_i = cx->GetComp(i);
      SmartPtr<Vector> g_f_i = cg_f->GetCompNonConst(i);
      DBG_ASSERT(IsValid(x_i) && IsValid(g_f_i));
      if (!nlps_[i]->Eval_grad_f(*x_i, *g_f_i)) {
        return false;
      }
    }

    // gradient of obj w.r.t q is zero
    SmartPtr<Vector> g_f_q = cg_f->GetCompNonConst(n_nlps);
    DBG_ASSERT(IsValid(g_f_q));
    g_f_q->Set(0.0);

    return true;
  }

  bool CompositeNLP::Eval_c(const Vector& x, Vector& c)
  {
    const CompoundVector* cx = dynamic_cast<const CompoundVector*>(&x);
    DBG_ASSERT(cx);
    CompoundVector* cc = dynamic_cast<CompoundVector*>(&c);
    DBG_ASSERT(cc);

    Index n_nlps = nlps_.size();
    SmartPtr<const Vector> q = cx->GetComp(n_nlps);
    for (Index i=0; i<n_nlps; i++) {
      SmartPtr<const Vector> x_i = cx->GetComp(i);
      SmartPtr<Vector> c_i = cc->GetCompNonConst(i);
      DBG_ASSERT(IsValid(x_i) && IsValid(c_i));
      if (!nlps_[i]->Eval_c(*x_i, *c_i)) {
        return false;
      }

      // Now calculate the part for the linking equations
      c_i = cc->GetCompNonConst(n_nlps + i);
      DBG_ASSERT(IsValid(c_i));
      Jx_linking_eqns_[i]->MultVector(1.0, *x_i, 0.0, *c_i);
      Jq_linking_eqns_[i]->MultVector(1.0, *q, 1.0, *c_i);
    }

    return true;
  }

  bool CompositeNLP::Eval_d(const Vector& x, Vector& d)
  {
    const CompoundVector* cx = dynamic_cast<const CompoundVector*>(&x);
    DBG_ASSERT(cx);
    CompoundVector* cd = dynamic_cast<CompoundVector*>(&d);
    DBG_ASSERT(cd);

    Index n_nlps = nlps_.size();
    SmartPtr<const Vector> q = cx->GetComp(n_nlps);
    for (Index i=0; i<n_nlps; i++) {
      SmartPtr<const Vector> x_i = cx->GetComp(i);
      SmartPtr<Vector> d_i = cd->GetCompNonConst(i);
      DBG_ASSERT(IsValid(x_i) && IsValid(d_i));
      if (!nlps_[i]->Eval_d(*x_i, *d_i)) {
        return false;
      }
    }

    return true;
  }

  bool CompositeNLP::Eval_jac_c(const Vector& x, Matrix& jac_c)
  {
    const CompoundVector* cx = dynamic_cast<const CompoundVector*>(&x);
    DBG_ASSERT(cx);
    CompoundMatrix* cj = dynamic_cast<CompoundMatrix*>(&jac_c);
    DBG_ASSERT(cj);

    Index n_nlps = nlps_.size();
    SmartPtr<const Vector> q = cx->GetComp(n_nlps);
    for (Index i=0; i<n_nlps; i++) {
      SmartPtr<const Vector> x_i = cx->GetComp(i);
      SmartPtr<Matrix> cj_i = cj->GetCompNonConst(i,i);
      DBG_ASSERT(IsValid(x_i) && IsValid(cj_i));
      if (!nlps_[i]->Eval_jac_c(*x_i, *cj_i)) {
        return false;
      }

      DBG_ASSERT(IsValid(Jx_linking_eqns_[i]) && IsValid(Jq_linking_eqns_[i]));
      cj->SetComp(n_nlps+i,i,*Jx_linking_eqns_[i]);
      cj->SetComp(n_nlps+i,n_nlps,*Jq_linking_eqns_[i]);
    }

    return true;
  }

  bool CompositeNLP::Eval_jac_d(const Vector& x, Matrix& jac_d)
  {
    const CompoundVector* cx = dynamic_cast<const CompoundVector*>(&x);
    DBG_ASSERT(cx);
    CompoundMatrix* cj = dynamic_cast<CompoundMatrix*>(&jac_d);
    DBG_ASSERT(cj);

    Index n_nlps = nlps_.size();
    SmartPtr<const Vector> q = cx->GetComp(n_nlps);
    for (Index i=0; i<n_nlps; i++) {
      SmartPtr<const Vector> x_i = cx->GetComp(i);
      SmartPtr<Matrix> cj_i = cj->GetCompNonConst(i,i);
      DBG_ASSERT(IsValid(x_i) && IsValid(cj_i));
      if (!nlps_[i]->Eval_jac_d(*x_i, *cj_i)) {
        return false;
      }
    }

    return true;
  }

  bool CompositeNLP::Eval_h(const Vector& x,
                            Number obj_factor,
                            const Vector& yc,
                            const Vector& yd,
                            SymMatrix& h)
  {
    const CompoundVector* cx = dynamic_cast<const CompoundVector*>(&x);
    DBG_ASSERT(cx);
    const CompoundVector* cyc = dynamic_cast<const CompoundVector*>(&yc);
    DBG_ASSERT(cyc);
    const CompoundVector* cyd = dynamic_cast<const CompoundVector*>(&yd);
    DBG_ASSERT(cyd);
    CompoundSymMatrix* ch = dynamic_cast<CompoundSymMatrix*>(&h);
    DBG_ASSERT(ch);

    Index n_nlps = nlps_.size();
    SmartPtr<const Vector> q = cx->GetComp(n_nlps);
    for (Index i=0; i<n_nlps; i++) {
      SmartPtr<const Vector> x_i = cx->GetComp(i);
      SmartPtr<const Vector> yc_i = cyc->GetComp(i);
      SmartPtr<const Vector> yd_i = cyd->GetComp(i);
      SmartPtr<Matrix> ch_i = ch->GetCompNonConst(i,i);
      SymMatrix* sh_i = dynamic_cast<SymMatrix*>(GetRawPtr(ch_i));
      DBG_ASSERT(IsValid(x_i) && IsValid(yc_i) && IsValid(yd_i) && IsValid(ch_i) && sh_i);
      if (!nlps_[i]->Eval_h(*x_i, obj_factor, *yc_i, *yd_i, *sh_i)) {
        return false;
      }
    }

    return true;
  }

  SmartPtr<SymLinearSolver> CompositeNLP::CreateLinearSolver()
  {
    // create the block solvers
    std::vector< SmartPtr<SymLinearSolver> > block_solvers;
    for (Index i=0; i<(Index)nlps_.size(); i++) {
      SmartPtr<SparseSymLinearSolverInterface> ma27_i 
	= new Ma27TSolverInterface();
      SmartPtr<SymLinearSolver> solver_i
	= new TSymLinearSolver(ma27_i, NULL);
      block_solvers.push_back(solver_i);
    }

    /*
    // create the permutations
    std::vector< SmartPtr<CompoundMatrixPermuter> > permuters;
    for (Index i=0; i<(Index)nlps_.size(); i++) {
      SmartPtr<CompoundMatrixPermuter> permuter
	= new CompoundMatrixPermuter();

      // Create the various permutations
      // H_i
      SmartPtr<SingleCompoundPermutation> permutation
	= new SingleCompoundPermutation(0,0);
      permutation->AddSourceDesignation(0, 0,0);
      permutation->AddSourceDesignation(1, i,i);
      permuter->AddPermutation(permutation);
      
      // D_s_i
      permutation = new SingleCompoundPermutation(1,1);
      permutation->AddSourceDesignation(0, 1,1);
      permutation->AddSourceDesignation(1, i,i);
      permuter->AddPermutation(permutation);

      // J_c_i
      permutation = new SingleCompoundPermutation(2,0);
      permutation->AddSourceDesignation(0, 2,0);
      permutation->AddSourceDesignation(1, i*2,i);
      permuter->AddPermutation(permutation);

      // D_c_i
      permutation = new SingleCompoundPermutation(2,2);
      permutation->AddSourceDesignation(0, 2,2);
      permutation->AddSourceDesignation(1, i,i);
      permuter->AddPermutation(permutation);

      // Lx_i
      permutation = new SingleCompoundPermutation(3,0);
      permutation->AddSourceDesignation(0, 2,0);
      permutation->AddSourceDesignation(1, i*2+1,i);
      permuter->AddPermutation(permutation);

      // D_L_i
      permutation = new SingleCompoundPermutation(3,3);
      permutation->AddSourceDesignation(0, 2,2);
      permutation->AddSourceDesignation(1, i*2+1,i*2+1);
      permuter->AddPermutation(permutation);

      // J_d_i
      permutation = new SingleCompoundPermutation(4,0);
      permutation->AddSourceDesignation(0, 3,0);
      permutation->AddSourceDesignation(1, i,i);
      permuter->AddPermutation(permutation);

      // J_d_ident_i
      permutation = new SingleCompoundPermutation(4,1);
      permutation->AddSourceDesignation(0, 3,1);
      permutation->AddSourceDesignation(1, i,i);
      permuter->AddPermutation(permutation);

      // D_d_i
      permutation = new SingleCompoundPermutation(4,4);
      permutation->AddSourceDesignation(0, 3,3);
      permutation->AddSourceDesignation(1, i,i);
      permuter->AddPermutation(permutation);

      permuters.push_back(permuter);
    }

    */
    SmartPtr<SymLinearSolver> mp_solver
      = new MPSymLinearSolver(block_solvers);

    return mp_solver;
  }
} // namespace Ipopt

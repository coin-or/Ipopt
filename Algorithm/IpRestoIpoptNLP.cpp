// Copyright (C) 2004, International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id$
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#include "IpRestoIpoptNLP.hpp"
#include "IpIdentityMatrix.hpp"
#include "IpSumSymMatrix.hpp"
#include "IpSumMatrix.hpp"

#ifdef OLD_C_HEADERS
#include <math.h>
#else
#include <cmath>
#endif

namespace Ipopt
{
  DBG_SET_VERBOSITY(0);

  RestoIpoptNLP::RestoIpoptNLP(IpoptNLP& orig_ip_nlp,
                               IpoptData& orig_ip_data,
                               IpoptCalculatedQuantities& orig_ip_cq,
                               IpoptData& curr_ip_data)
      :
      IpoptNLP(),
      orig_ip_nlp_(&orig_ip_nlp),
      orig_ip_data_(&orig_ip_data),
      orig_ip_cq_(&orig_ip_cq),
      ip_data_(&curr_ip_data),
      eta_factor_(1.0),
      eta_mu_exponent_(0.5),
      rho_(1000.)
  {}

  RestoIpoptNLP::~RestoIpoptNLP()
  {}

  bool RestoIpoptNLP::Initialize(const Journalist& jnlst,
                                 const OptionsList& options,
                                 const std::string& prefix)
  {
    initialized_ = true;
    return true;
  }

  bool RestoIpoptNLP::InitializeStructures(SmartPtr<Vector>& x,
      bool init_x,
      SmartPtr<Vector>& y_c,
      bool init_y_c,
      SmartPtr<Vector>& y_d,
      bool init_y_d,
      SmartPtr<Vector>& z_L,
      bool init_z_L,
      SmartPtr<Vector>& z_U,
      bool init_z_U,
      SmartPtr<Vector>& v_L,
      bool init_v_L,
      SmartPtr<Vector>& v_U,
      bool init_v_U
                                          )
  {
    DBG_START_METH("RestoIpoptNLP::InitializeStructures", 0);
    DBG_ASSERT(initialized_);
    ///////////////////////////////////////////////////////////
    // Get the vector/matrix spaces for the original problem //
    ///////////////////////////////////////////////////////////

    SmartPtr<const VectorSpace> orig_x_space;
    SmartPtr<const VectorSpace> orig_c_space;
    SmartPtr<const VectorSpace> orig_d_space;
    SmartPtr<const VectorSpace> orig_x_l_space;
    SmartPtr<const MatrixSpace> orig_px_l_space;
    SmartPtr<const VectorSpace> orig_x_u_space;
    SmartPtr<const MatrixSpace> orig_px_u_space;
    SmartPtr<const VectorSpace> orig_d_l_space;
    SmartPtr<const MatrixSpace> orig_pd_l_space;
    SmartPtr<const VectorSpace> orig_d_u_space;
    SmartPtr<const MatrixSpace> orig_pd_u_space;
    SmartPtr<const MatrixSpace> orig_jac_c_space;
    SmartPtr<const MatrixSpace> orig_jac_d_space;
    SmartPtr<const SymMatrixSpace> orig_h_space;

    orig_ip_nlp_->GetSpaces(orig_x_space, orig_c_space, orig_d_space,
                            orig_x_l_space, orig_px_l_space,
                            orig_x_u_space, orig_px_u_space,
                            orig_d_l_space, orig_pd_l_space,
                            orig_d_u_space, orig_pd_u_space,
                            orig_jac_c_space, orig_jac_d_space,
                            orig_h_space);

    // Create the restoration phase problem vector/matrix spaces, based
    // on the original spaces (pretty inconvenient with all the
    // matrix spaces, isn't it?!?)

    DBG_PRINT((1, "Creating the x_space_\n"));
    // vector x
    Index total_dim = orig_x_space->Dim() + 2*orig_c_space->Dim()
                      //orig  + orig_d_l_space->Dim() + orig_d_u_space->Dim();
                      + 2*orig_d_space->Dim();
    x_space_ = new CompoundVectorSpace(5, total_dim);
    x_space_->SetCompSpace(0, *orig_x_space);
    x_space_->SetCompSpace(1, *orig_c_space); // n_c
    x_space_->SetCompSpace(2, *orig_c_space); // p_c
    x_space_->SetCompSpace(3, *orig_d_space); // n_d
    x_space_->SetCompSpace(4, *orig_d_space); // p_d
    //orig    x_space_->SetCompSpace(3, *orig_d_l_space); // n_d
    //orig    x_space_->SetCompSpace(4, *orig_d_u_space); // p_d

    DBG_PRINT((1, "Setting the c_space_\n"));
    // vector c
    c_space_ = orig_c_space;

    DBG_PRINT((1, "Setting the d_space_\n"));
    // vector d
    d_space_ = orig_d_space;

    DBG_PRINT((1, "Creating the x_l_space_\n"));
    // vector x_L
    total_dim = orig_x_l_space->Dim() + 2*orig_c_space->Dim()
                //orig      + orig_d_l_space->Dim() + orig_d_u_space->Dim();
                + 2*orig_d_space->Dim();
    x_l_space_ = new CompoundVectorSpace(5, total_dim);
    x_l_space_->SetCompSpace(0, *orig_x_l_space);
    x_l_space_->SetCompSpace(1, *orig_c_space); // n_c >=0
    x_l_space_->SetCompSpace(2, *orig_c_space); // p_c >=0
    x_l_space_->SetCompSpace(3, *orig_d_space); // n_d >=0
    x_l_space_->SetCompSpace(4, *orig_d_space); // p_d >=0
    //orig    x_l_space_->SetCompSpace(3, *orig_d_l_space); // n_d >=0
    //orig    x_l_space_->SetCompSpace(4, *orig_d_u_space); // p_d >=0

    DBG_PRINT((1, "Setting the x_u_space_\n"));
    // vector x_U
    x_u_space_ = orig_x_u_space;

    DBG_PRINT((1, "Creating the px_l_space_\n"));
    // matrix px_l
    Index total_rows = orig_x_space->Dim() + 2*orig_c_space->Dim()
                       //orig      + orig_d_l_space->Dim() + orig_d_u_space->Dim();
                       + 2*orig_d_space->Dim();
    Index total_cols = orig_x_l_space->Dim() + 2*orig_c_space->Dim()
                       //orig      + orig_d_l_space->Dim() + orig_d_u_space->Dim();
                       + 2*orig_d_space->Dim();
    px_l_space_ = new CompoundMatrixSpace(5, 5, total_rows, total_cols);
    px_l_space_->SetBlockRows(0, orig_x_space->Dim());
    px_l_space_->SetBlockRows(1, orig_c_space->Dim());
    px_l_space_->SetBlockRows(2, orig_c_space->Dim());
    px_l_space_->SetBlockRows(3, orig_d_space->Dim());
    px_l_space_->SetBlockRows(4, orig_d_space->Dim());
    //orig    px_l_space_->SetBlockRows(3, orig_d_l_space->Dim());
    //orig    px_l_space_->SetBlockRows(4, orig_d_u_space->Dim());
    px_l_space_->SetBlockCols(0, orig_x_l_space->Dim());
    px_l_space_->SetBlockCols(1, orig_c_space->Dim());
    px_l_space_->SetBlockCols(2, orig_c_space->Dim());
    px_l_space_->SetBlockCols(3, orig_d_space->Dim());
    px_l_space_->SetBlockCols(4, orig_d_space->Dim());
    //orig    px_l_space_->SetBlockCols(3, orig_d_l_space->Dim());
    //orig    px_l_space_->SetBlockCols(4, orig_d_u_space->Dim());

    px_l_space_->SetCompSpace(0, 0, *orig_px_l_space);
    // now setup the identity matrix
    // This could be changed to be something like...
    // px_l_space_->SetBlockToIdentity(1,1,1.0);
    // px_l_space_->SetBlockToIdentity(2,2,other_factor);
    // ... etc with some simple changes to the CompoundMatrixSpace
    // to allow this (space should auto create the matrices)
    //
    // for now, we use the new feature and set the true flag for this block
    // to say that the matrices should be auto_allocated
    SmartPtr<const MatrixSpace> identity_mat_space_nc
    = new IdentityMatrixSpace(orig_c_space->Dim());
    px_l_space_->SetCompSpace(1, 1, *identity_mat_space_nc, true);
    px_l_space_->SetCompSpace(2, 2, *identity_mat_space_nc, true);
    //orig    SmartPtr<const MatrixSpace> identity_mat_space_nd_l
    //orig    = new IdentityMatrixSpace(orig_d_l_space->Dim());
    SmartPtr<const MatrixSpace> identity_mat_space_nd
    = new IdentityMatrixSpace(orig_d_space->Dim());
    //orig    px_l_space_->SetCompSpace(3, 3, *identity_mat_space_nd_l, true);
    px_l_space_->SetCompSpace(3, 3, *identity_mat_space_nd, true);
    //orig    SmartPtr<const MatrixSpace> identity_mat_space_nd_u
    //orig      = new IdentityMatrixSpace(orig_d_u_space->Dim());
    //orig    px_l_space_->SetCompSpace(4, 4, *identity_mat_space_nd_u, true);
    px_l_space_->SetCompSpace(4, 4, *identity_mat_space_nd, true);

    DBG_PRINT((1, "Creating the px_u_space_\n"));
    // matrix px_u    px_u_space_->SetBlockRows(0, orig_x_space->Dim());

    total_rows = orig_x_space->Dim() + 2*orig_c_space->Dim()
                 //orig      + orig_d_l_space->Dim() + orig_d_u_space->Dim();
                 + 2*orig_d_space->Dim();
    total_cols = orig_x_u_space->Dim();
    DBG_PRINT((1, "total_rows = %d, total_cols = %d\n",total_rows, total_cols));
    px_u_space_ = new CompoundMatrixSpace(5, 1, total_rows, total_cols);
    px_u_space_->SetBlockRows(0, orig_x_space->Dim());
    px_u_space_->SetBlockRows(1, orig_c_space->Dim());
    px_u_space_->SetBlockRows(2, orig_c_space->Dim());
    px_u_space_->SetBlockRows(3, orig_d_space->Dim());
    px_u_space_->SetBlockRows(4, orig_d_space->Dim());
    //orig    px_u_space_->SetBlockRows(3, orig_d_l_space->Dim());
    //orig    px_u_space_->SetBlockRows(4, orig_d_u_space->Dim());
    px_u_space_->SetBlockCols(0, orig_x_u_space->Dim());

    px_u_space_->SetCompSpace(0, 0, *orig_px_u_space);
    // other matrices are zero'ed out

    // vector d_L
    d_l_space_ = orig_d_l_space;

    // vector d_U
    d_u_space_ = orig_d_u_space;

    // matrix pd_L
    pd_l_space_ = orig_pd_l_space;

    // matrix pd_U
    pd_u_space_ = orig_pd_u_space;

    DBG_PRINT((1, "Creating the jac_c_space_\n"));
    // matrix jac_c
    total_rows = orig_c_space->Dim();
    total_cols = orig_x_space->Dim() + 2*orig_c_space->Dim()
                 + 2*orig_d_space->Dim();
    //orig      + orig_d_l_space->Dim() + orig_d_u_space->Dim();
    jac_c_space_ = new CompoundMatrixSpace(1, 5, total_rows, total_cols);
    jac_c_space_->SetBlockRows(0, orig_c_space->Dim());
    jac_c_space_->SetBlockCols(0, orig_x_space->Dim());
    jac_c_space_->SetBlockCols(1, orig_c_space->Dim());
    jac_c_space_->SetBlockCols(2, orig_c_space->Dim());
    jac_c_space_->SetBlockCols(3, orig_d_space->Dim());
    jac_c_space_->SetBlockCols(4, orig_d_space->Dim());
    //orig    jac_c_space_->SetBlockCols(3, orig_d_l_space->Dim());
    //orig    jac_c_space_->SetBlockCols(4, orig_d_u_space->Dim());

    jac_c_space_->SetCompSpace(0, 0, *orig_jac_c_space);
    jac_c_space_->SetCompSpace(0, 1, *identity_mat_space_nc, true);
    jac_c_space_->SetCompSpace(0, 2, *identity_mat_space_nc, true);
    // remaining blocks are zero'ed

    DBG_PRINT((1, "Creating the jac_d_space_\n"));
    // matrix jac_d
    total_rows = orig_d_space->Dim();
    total_cols = orig_x_space->Dim() + 2*orig_c_space->Dim()
                 + 2*orig_d_space->Dim();
    //orig      + orig_d_l_space->Dim() + orig_d_u_space->Dim();
    jac_d_space_ = new CompoundMatrixSpace(1, 5, total_rows, total_cols);
    jac_d_space_->SetBlockRows(0, orig_d_space->Dim());
    jac_d_space_->SetBlockCols(0, orig_x_space->Dim());
    jac_d_space_->SetBlockCols(1, orig_c_space->Dim());
    jac_d_space_->SetBlockCols(2, orig_c_space->Dim());
    jac_d_space_->SetBlockCols(3, orig_d_space->Dim());
    jac_d_space_->SetBlockCols(4, orig_d_space->Dim());
    //orig    jac_d_space_->SetBlockCols(3, orig_d_l_space->Dim());
    //orig    jac_d_space_->SetBlockCols(4, orig_d_u_space->Dim());

    jac_d_space_->SetCompSpace(0, 0, *orig_jac_d_space);
    // Blocks (0,1) and (0,2) are zero'ed out
    jac_d_space_->SetCompSpace(0, 3, *identity_mat_space_nd, true);
    jac_d_space_->SetCompSpace(0, 4, *identity_mat_space_nd, true);
    //orig    jac_d_space_->SetCompSpace(0, 3, *orig_pd_l_space, true);
    //orig    SmartPtr<SumMatrixSpace> sum_pd_u
    //orig    = new SumMatrixSpace(orig_d_space->Dim(), orig_d_u_space->Dim(), 1);
    //orig    jac_d_space_->SetCompSpace(0, 4, *sum_pd_u, true);

    DBG_PRINT((1, "Creating the h_space_\n"));
    // matrix h
    total_dim = orig_x_space->Dim() + 2*orig_c_space->Dim()
                + 2*orig_d_space->Dim();
    //orig      + orig_d_l_space->Dim() + orig_d_u_space->Dim();
    h_space_ = new CompoundSymMatrixSpace(5, total_dim);
    h_space_->SetBlockDim(0, orig_x_space->Dim());
    h_space_->SetBlockDim(1, orig_c_space->Dim());
    h_space_->SetBlockDim(2, orig_c_space->Dim());
    h_space_->SetBlockDim(3, orig_d_space->Dim());
    h_space_->SetBlockDim(4, orig_d_space->Dim());
    //orig    h_space_->SetBlockDim(3, orig_d_l_space->Dim());
    //orig    h_space_->SetBlockDim(4, orig_d_u_space->Dim());

    SmartPtr<const MatrixSpace> sumsym_mat_space =
      new SumSymMatrixSpace(orig_x_space->Dim(), 2);
    h_space_->SetCompSpace(0, 0, *sumsym_mat_space, true);
    // All remaining blocks are zero'ed out

    ///////////////////////////
    // Create the bound data //
    ///////////////////////////

    // x_L
    x_L_ = x_l_space_->MakeNewCompoundVector();
    x_L_->SetComp(0, *orig_ip_nlp_->x_L()); // x >= x_L
    x_L_->GetCompNonConst(1)->Set(0.0); // n_c >= 0
    x_L_->GetCompNonConst(2)->Set(0.0); // p_c >= 0
    x_L_->GetCompNonConst(3)->Set(0.0); // n_d >= 0
    x_L_->GetCompNonConst(4)->Set(0.0); // p_d >= 0
    DBG_PRINT_VECTOR(2,"resto_x_L", *x_L_);

    // x_U
    x_U_ = orig_ip_nlp_->x_U();

    // d_L
    d_L_ = orig_ip_nlp_->d_L();

    // d_U
    d_U_ = orig_ip_nlp_->d_U();

    // Px_L
    Px_L_ = px_l_space_->MakeNewCompoundMatrix();
    Px_L_->SetComp(0, 0, *orig_ip_nlp_->Px_L());
    // Identities are auto-created (true flag passed into SetCompSpace)

    // Px_U
    Px_U_ = px_u_space_->MakeNewCompoundMatrix();
    Px_U_->SetComp(0, 0, *orig_ip_nlp_->Px_U());
    // Remaining matrices will be zero'ed out

    // Pd_L
    Pd_L_ = orig_ip_nlp_->Pd_L();

    // Pd_U
    Pd_U_ = orig_ip_nlp_->Pd_U();

    /////////////////////////////////////////////////////////////////////////
    // Create and initialize the vectors for the restoration phase problem //
    /////////////////////////////////////////////////////////////////////////

    // Vector x
    SmartPtr<CompoundVector> comp_x = x_space_->MakeNewCompoundVector();
    if (init_x) {
      comp_x->GetCompNonConst(0)->Copy(*orig_ip_data_->curr_x());
      comp_x->GetCompNonConst(1)->Set(1.0);
      comp_x->GetCompNonConst(2)->Set(1.0);
      comp_x->GetCompNonConst(3)->Set(1.0);
      comp_x->GetCompNonConst(4)->Set(1.0);
    }
    x = GetRawPtr(comp_x);

    // Vector y_c
    y_c = c_space_->MakeNew();
    if (init_y_c) {
      y_c->Set(0.0);  // ToDo
    }

    // Vector y_d
    y_d = d_space_->MakeNew();
    if (init_y_d) {
      y_d->Set(0.0);
    }

    // Vector z_L
    z_L = x_l_space_->MakeNew();
    if (init_z_L) {
      z_L->Set(1.0);
    }

    // Vector z_U
    z_U = x_u_space_->MakeNew();
    if (init_z_U) {
      z_U->Set(1.0);
    }

    // Vector v_L
    v_L = d_l_space_->MakeNew();

    // Vector v_U
    v_U = d_u_space_->MakeNew();

    // Initialize other data needed by the restoration nlp.  x_ref is
    // the point to reference to which we based the regularization
    // term
    x_ref_ = orig_x_space->MakeNew();
    x_ref_->Copy(*orig_ip_data_->curr_x());

    SmartPtr<DiagMatrixSpace> DR_x_space
    = new DiagMatrixSpace(orig_x_space->Dim());
    dr_x_ = orig_x_space->MakeNew();
    dr_x_->Set(1.0);
    SmartPtr<Vector> tmp = dr_x_->MakeNew();
    tmp->Copy(*x_ref_);
    dr_x_->ElementWiseMax(*tmp);
    tmp->Scal(-1.);
    dr_x_->ElementWiseMax(*tmp);
    dr_x_->ElementWiseReciprocal();
    DBG_PRINT_VECTOR(2, "dr_x_", *dr_x_);
    DR_x_ = DR_x_space->MakeNewDiagMatrix();
    DR_x_->SetDiag(*dr_x_);

    return true;
  }

  Number RestoIpoptNLP::f(const Vector& x)
  {
    DBG_START_METH("RestoIpoptNLP::f",
                   dbg_verbosity);
    Number ret = 0.0;
    // rho*(pcTe + ncTe + pdT*e + ndT*e) + eta/2*||Dr*(x-xr)||_2^2
    const CompoundVector* c_vec = dynamic_cast<const CompoundVector*>(&x);
    DBG_ASSERT(c_vec);
    SmartPtr<const Vector> x_only = c_vec->GetComp(0);
    ret = x.Sum() - x_only->Sum();
    DBG_PRINT((1,"xdiff sum = %e\n",ret));
    ret = rho_ * ret;
    DBG_PRINT((1,"rho_ = %e\n",rho_));

    SmartPtr<Vector> x_diff = x_only->MakeNew();
    x_diff->Copy(*x_only);
    x_diff->Axpy(-1.0, *x_ref_);
    DBG_PRINT_VECTOR(2,"x_ref",*x_ref_);
    x_diff->ElementWiseMultiply(*dr_x_);
    Number ret2 = x_diff->Nrm2();
    DBG_PRINT((1,"Eta = %e\n",Eta()));
    ret2 = Eta()/2.0*ret2*ret2;

    ret += ret2;
    return ret;
  }

  SmartPtr<const Vector> RestoIpoptNLP::grad_f(const Vector& x)
  {
    SmartPtr<Vector> retPtr = x.MakeNew();
    // Scale the p's and n's by rho (Scale all, take out the x part later)
    retPtr->Set(rho_);

    const CompoundVector* c_vec_in = dynamic_cast<const CompoundVector*>(&x);
    SmartPtr<const Vector> x_only_in = c_vec_in->GetComp(0);

    CompoundVector* c_vec = dynamic_cast<CompoundVector*>(GetRawPtr(retPtr));
    DBG_ASSERT(c_vec);
    SmartPtr<Vector> x_only = c_vec->GetCompNonConst(0);
    x_only->Copy(*x_only_in);
    x_only->Axpy(-1.0, *x_ref_);
    x_only->ElementWiseMultiply(*dr_x_);
    x_only->Scal(Eta());

    return ConstPtr(retPtr);
  }

  SmartPtr<const Vector> RestoIpoptNLP::c(const Vector& x)
  {
    const CompoundVector* c_vec = dynamic_cast<const CompoundVector*>(&x);
    SmartPtr<const Vector> x_only = c_vec->GetComp(0);
    SmartPtr<const Vector> nc_only = c_vec->GetComp(1);
    SmartPtr<const Vector> pc_only = c_vec->GetComp(2);

    SmartPtr<const Vector> orig_c = orig_ip_nlp_->c(*x_only);
    SmartPtr<Vector> retPtr = c_space_->MakeNew();
    retPtr->Copy(*orig_c);
    retPtr->Axpy(1.0, *nc_only);
    retPtr->Axpy(-1.0, *pc_only);

    return GetRawPtr(retPtr);
  }


  SmartPtr<const Vector> RestoIpoptNLP::d(const Vector& x)
  {
    const CompoundVector* c_vec = dynamic_cast<const CompoundVector*>(&x);
    SmartPtr<const Vector> x_only = c_vec->GetComp(0);
    SmartPtr<const Vector> nd_only = c_vec->GetComp(3);
    SmartPtr<const Vector> pd_only = c_vec->GetComp(4);

    SmartPtr<const Vector> orig_d = orig_ip_nlp_->d(*x_only);
    SmartPtr<Vector> retPtr = d_space_->MakeNew();
    retPtr->Copy(*orig_d);
    retPtr->Axpy(1., *nd_only);
    retPtr->Axpy(-1., *pd_only);
#ifdef orig

    SmartPtr<Vector> tmp = orig_d->MakeNew();
    orig_ip_nlp_->Pd_L()->MultVector(1.0, *nd_only, 0.0, *tmp);
    retPtr->Axpy(1.0, *tmp);
    orig_ip_nlp_->Pd_U()->MultVector(1.0, *pd_only, 0.0, *tmp);
    retPtr->Axpy(-1.0, *tmp);
#endif

    return GetRawPtr(retPtr);
  }

  SmartPtr<const Matrix> RestoIpoptNLP::jac_c(const Vector& x)
  {

    // Here, we set the (0,0) block with the values from the
    // original jac_c and set the factor for the -I (jac w.r.t. p_c)

    // get out the x_only part
    const CompoundVector* c_vec = dynamic_cast<const CompoundVector*>(&x);
    DBG_ASSERT(c_vec);
    SmartPtr<const Vector> x_only = c_vec->GetComp(0);

    // calculate the jacobian for the original problem
    SmartPtr<const Matrix> jac_c_only = orig_ip_nlp_->jac_c(*x_only);

    // Create the new compound matrix
    // The zero parts remain NULL, the identities are created from the matrix
    // space (since auto_allocate was set to true in SetCompSpace)
    SmartPtr<CompoundMatrix> retPtr = jac_c_space_->MakeNewCompoundMatrix();

    // set the (0,0) block to the original jacobian
    retPtr->SetComp(0,0,*jac_c_only);

    // we currently do not have a default factor in the matrix spaces
    // so we need to set the factor on the identity (jacobian of the
    // restoration c w.r.t. p_c is -I)
    // This could easily be changed to include special processing
    // for identities in the CompoundMatrixSpace (and a factor)
    SmartPtr<Matrix> jac_c_pc_mat = retPtr->GetCompNonConst(0,2);
    IdentityMatrix* jac_c_pc = dynamic_cast<IdentityMatrix*>(GetRawPtr(jac_c_pc_mat));
    DBG_ASSERT(jac_c_pc);
    jac_c_pc->SetFactor(-1.0);

    return GetRawPtr(retPtr);
  }

  SmartPtr<const Matrix> RestoIpoptNLP::jac_d(const Vector& x)
  {
    // Here, we set the (0,0) block with the values from the
    // original jac_d and set the factor for the -I (jac w.r.t. p_d)

    // get out the x_only part
    const CompoundVector* c_vec = dynamic_cast<const CompoundVector*>(&x);
    DBG_ASSERT(c_vec);
    SmartPtr<const Vector> x_only = c_vec->GetComp(0);

    // calculate the jacobian for the original problem
    SmartPtr<const Matrix> jac_d_only = orig_ip_nlp_->jac_d(*x_only);

    // Create the new compound matrix
    // The zero parts remain NULL, the identities are created from the matrix
    // space (since auto_allocate was set to true in SetCompSpace)
    SmartPtr<CompoundMatrix> retPtr = jac_d_space_->MakeNewCompoundMatrix();

    // Set the block for the original jacobian
    retPtr->SetComp(0,0,*jac_d_only);

    // (0,1) and (0,2) blocks are zero (NULL)

    // set the factor for the identity matrix for the pd variables
    // (likr in jac_c)
    SmartPtr<Matrix> jac_d_pd_mat = retPtr->GetCompNonConst(0,4);
    IdentityMatrix* jac_d_pd = dynamic_cast<IdentityMatrix*>(GetRawPtr(jac_d_pd_mat));
    DBG_ASSERT(jac_d_pd);
    jac_d_pd->SetFactor(-1.0);

#ifdef orig
    // Jacobian of resto d w.r.t. n_d is Pd_L
    retPtr->SetComp(0, 3, *orig_ip_nlp_->Pd_L());

    // Change the matrix factor to -1 for Pd_U
    // Jacobian of the resto d w.r.t. p_d is -Pd_U
    SmartPtr<Matrix> jac_d_pd_mat = retPtr->GetCompNonConst(0,4);
    SumMatrix* jac_d_pd_sum = dynamic_cast<SumMatrix*>(GetRawPtr(jac_d_pd_mat));
    DBG_ASSERT(jac_d_pd_sum);
    jac_d_pd_sum->SetTerm(0, -1.0, *orig_ip_nlp_->Pd_U());
#endif

    return GetRawPtr(retPtr);
  }

  SmartPtr<const SymMatrix> RestoIpoptNLP::h(const Vector& x,
      Number obj_factor,
      const Vector& yc,
      const Vector& yd)
  {
    // Here, we use a SumSymMatrix for the (0,0) block of the
    // Hessian. We need to set this to the hessian of the restoration
    // problem, which is the hessian of the objective from the restoration
    // problem + the constraint only part of the hessian from the original
    // problem
    // All other blocks are zero'ed (NULL)

    // get the x_only part
    const CompoundVector* c_vec = dynamic_cast<const CompoundVector*>(&x);
    DBG_ASSERT(c_vec);
    SmartPtr<const Vector> x_only = c_vec->GetComp(0);

    // yc and yd should not be compound vectors

    // calculate the original hessian
    SmartPtr<const SymMatrix> h_con_orig = orig_ip_nlp_->h(*x_only, 0.0, yc, yd);

    // Create the new compound matrix
    // The SumSymMatrix is auto_allocated
    SmartPtr<CompoundSymMatrix> retPtr = h_space_->MakeNewCompoundSymMatrix();

    // Set the entries in the SumSymMatrix
    SmartPtr<Matrix> h_sum_mat = retPtr->GetCompNonConst(0,0);
    SmartPtr<SumSymMatrix> h_sum = dynamic_cast<SumSymMatrix*>(GetRawPtr(h_sum_mat));
    h_sum->SetTerm(0, 1.0, *h_con_orig);
    h_sum->SetTerm(1, obj_factor*Eta(), *DR_x_);

    return GetRawPtr(retPtr);
  }

  void RestoIpoptNLP::GetSpaces(SmartPtr<const VectorSpace>& x_space,
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
    x_space = GetRawPtr(x_space_);
    c_space = GetRawPtr(c_space_);
    d_space = GetRawPtr(d_space_);
    x_l_space = GetRawPtr(x_l_space_);
    px_l_space = GetRawPtr(px_l_space_);
    x_u_space = GetRawPtr(x_u_space_);
    px_u_space = GetRawPtr(px_u_space_);
    d_l_space = GetRawPtr(d_l_space_);
    pd_l_space = GetRawPtr(pd_l_space_);
    d_u_space = GetRawPtr(d_u_space_);
    pd_u_space = GetRawPtr(pd_u_space_);
    Jac_c_space = GetRawPtr(jac_c_space_);
    Jac_d_space = GetRawPtr(jac_d_space_);
    Hess_lagrangian_space = GetRawPtr(h_space_);
  }

  Number RestoIpoptNLP::Eta() const
  {
    return eta_factor_ * pow(ip_data_->curr_mu(), eta_mu_exponent_);
  }

  void RestoIpoptNLP::AdjustVariableBounds(const Vector& new_x_L, const Vector& new_x_U,
      const Vector& new_d_L, const Vector& new_d_U)
  {

    const CompoundVector* comp_new_x_L =
      dynamic_cast<const CompoundVector*>(&new_x_L);
    DBG_ASSERT(comp_new_x_L);

    SmartPtr<const Vector> new_orig_x_L = comp_new_x_L->GetComp(0);

    // adapt bounds for the original NLP
    orig_ip_nlp_->AdjustVariableBounds(*new_orig_x_L, new_x_U, new_d_L, new_d_U);

    // adapt bounds for the p and n variables
    SmartPtr<const Vector> new_nc_L = comp_new_x_L->GetComp(1);
    SmartPtr<const Vector> new_pc_L = comp_new_x_L->GetComp(2);
    SmartPtr<const Vector> new_nd_L = comp_new_x_L->GetComp(3);
    SmartPtr<const Vector> new_pd_L = comp_new_x_L->GetComp(4);

    x_L_->GetCompNonConst(1)->Copy(*new_nc_L);
    x_L_->GetCompNonConst(2)->Copy(*new_pc_L);
    x_L_->GetCompNonConst(3)->Copy(*new_nd_L);
    x_L_->GetCompNonConst(4)->Copy(*new_pd_L);

  }

} // namespace Ipopt

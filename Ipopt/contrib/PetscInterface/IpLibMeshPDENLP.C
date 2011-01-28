// Copyright (C) 2010 International Business Machines.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// $Id$
//
// Authors:  Johannes Huber, Andreas Waechter     IBM        2010-09-03
#include "mpi.h"
#include "IpLibMeshPDENLP.hpp"
#include "IpJournalist.hpp"

#include "petscmat.h"
#include "sparse_matrix.h"
#include "petsc_matrix.h"
#include "petsc_vector.h"
#include "dense_vector.h"

int GetProcID();

//#define DBG_PRINT(s) {std::cout << GetProcID() << ": " << s << std::endl;}
#define DBG_PRINT(s) {}

using namespace Ipopt;
using namespace libMesh;

typedef libMesh::Number lm_Number;

LibMeshPDENLP::LibMeshPDENLP(LibMeshPDEBase& libmeshPDE,
                             Journalist& jnlst)
    :
    libmeshPDE_(&libmeshPDE),
    jnlst_(&jnlst)
{}

bool
LibMeshPDENLP::
get_nlp_info(Index num_proc, Index proc_id,
             Index& n, Index& n_first, Index& n_last,
             Index& m, Index& m_first, Index& m_last,
             Index& nnz_jac_g_part, Index& nnz_h_lag_part,
             IndexStyleEnum& index_style)
{
  DBG_PRINT("LibMeshDPENLP::get_nlp_info called")

//  NumericVector<lm_Number>& lm_control_vec = libmeshPDE_->getControlVector();
//  ImplicitSystem& lm_sys = libmeshPDE_->getImplicitSystem();
  int m_pde;
  int m_aux;
  int n_state;
  int n_control;
  // Distribution of constraints
  {
    PetscInt tmp_first;
    PetscInt tmp_end;
    PetscVector<lm_Number> *p_vec = dynamic_cast<PetscVector<lm_Number>*>(&(libmeshPDE_->getAuxConstrVector()));
    Vec petsc_vec = p_vec->vec();
    VecGetOwnershipRange(petsc_vec,&tmp_first,&tmp_end);
    aux_first_ = (int)tmp_first;
    aux_last_ = (int)(tmp_end-1);
    VecGetSize(petsc_vec,&m_aux);
    p_vec = dynamic_cast<PetscVector<lm_Number>*>(&(libmeshPDE_->getPDEResVector()));
    petsc_vec = p_vec->vec();
    VecGetOwnershipRange(petsc_vec,&tmp_first,&tmp_end);
    pde_first_ = (int)tmp_first;
    pde_last_ = (int)(tmp_end-1);
    VecGetSize(petsc_vec,&m_pde);

    int m_local = aux_last_-aux_first_+1 + pde_last_-pde_first_+1;
    int m_glob_end;
    MPI_Scan(&m_local, &m_glob_end, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    m_first_ = m_first = m_glob_end - m_local;
    m_last_ = m_last = m_glob_end-1;
    m = m_pde+m_aux;
  }

  // Distribution of variables (state & controls)
  {
    // states distribution
    PetscVector<lm_Number>* p_vec = dynamic_cast<PetscVector<lm_Number>*>(&(libmeshPDE_->getStateVector()));

    PetscInt p_first;
    PetscInt p_last;
    Vec v = p_vec->vec();
    VecGetOwnershipRange(v, &p_first, &p_last);
    state_first_ = p_first;
    state_last_ = p_last-1;
    VecGetSize(v,&n_state);
    // control distribution
    p_vec = dynamic_cast<PetscVector<lm_Number>*>(&(libmeshPDE_->getControlVector()));
    v = p_vec->vec();
    VecGetOwnershipRange(v, &p_first, &p_last);
    control_first_ = p_first;
    control_last_ = p_last-1;
    VecGetSize(v,&n_control);
    // Get the local offset in large x vector (ParNLP distribution)
    int n_local = state_last_-state_first_+1 + control_last_-control_first_+1;
    int n_glob_end;
    MPI_Scan(&n_local, &n_glob_end, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    n_first_ = n_first = n_glob_end - n_local;
    n_last_ = n_last = n_glob_end-1;
    n = n_state+n_control;
  }

  libmeshPDE_->calc_type_ = LibMeshPDEBase::StructureOnly;
  SparseMatrix<lm_Number>* jac_state;
  SparseMatrix<lm_Number>* jac_control;
  libmeshPDE_->calcPDE_jacobians(jac_state,jac_control);

  PetscMatrix<lm_Number>* p_jac_state = dynamic_cast<PetscMatrix<lm_Number>*>(jac_state);
  Mat mat = p_jac_state->mat();
  MatInfo info;
  MatGetInfo(mat,MAT_LOCAL,&info);
  nnz_jac_g_part = info.nz_used;
  PetscMatrix<lm_Number>* p_jac_control = dynamic_cast<PetscMatrix<lm_Number>*>(jac_control);
  mat = p_jac_control->mat();
  MatGetInfo(mat,MAT_LOCAL,&info);
  nnz_jac_g_part += info.nz_used;

  SparseMatrix<lm_Number>* jac_aux_state;
  SparseMatrix<lm_Number>* jac_aux_control;
  libmeshPDE_->calcAux_jacobians(jac_aux_state,jac_aux_control);
  PetscMatrix<lm_Number>* p_aux_jac_state = dynamic_cast<PetscMatrix<lm_Number>*>(jac_aux_state);
  mat = p_aux_jac_state->mat();
  MatGetInfo(mat,MAT_LOCAL,&info);
  nnz_jac_g_part += info.nz_used;
  PetscMatrix<lm_Number>* p_aux_jac_control = dynamic_cast<PetscMatrix<lm_Number>*>(jac_aux_control);
  mat = p_aux_jac_control->mat();
  MatGetInfo(mat,MAT_LOCAL,&info);
  nnz_jac_g_part += info.nz_used;

  nnz_h_lag_part = 0;

  SparseMatrix<lm_Number>* hss;
  SparseMatrix<lm_Number>* hcc;
  SparseMatrix<lm_Number>* hcs;
  libMesh::DenseVector<libMesh::Number> lambda_pde;
  libMesh::DenseVector<lm_Number> lambda_aux;

  lambda_pde.resize(m_pde);
  for (int i=0.0;i<lambda_pde.size();++i)
    lambda_pde(i) = 1.0;
  lambda_aux.resize(m_aux);
  for (int i=0.0;i<lambda_aux.size();++i)
    lambda_aux(i) = 1.0;
  libmeshPDE_->calc_hessians(1.0,lambda_pde,lambda_aux,hcc,hcs,hss);
  PetscMatrix<lm_Number>* lm_mat_h = dynamic_cast<PetscMatrix<lm_Number>*>(hcc);
  mat = lm_mat_h->mat();
  MatGetInfo(mat,MAT_LOCAL,&info);
  nnz_h_lag_part += info.nz_used;
  lm_mat_h = dynamic_cast<PetscMatrix<lm_Number>*>(hcs);
  mat = lm_mat_h->mat();
  MatGetInfo(mat,MAT_LOCAL,&info);
  nnz_h_lag_part += info.nz_used;
  lm_mat_h = dynamic_cast<PetscMatrix<lm_Number>*>(hss);
  mat = lm_mat_h->mat();
  MatGetInfo(mat,MAT_LOCAL,&info);
  nnz_h_lag_part += info.nz_used;

  index_style = C_STYLE;
  DBG_PRINT( "LibMeshPDENLP::get_nlp_info finished" )
  libmeshPDE_->calc_type_ = LibMeshPDEBase::Values;
  return true;
}

void LibMeshPDENLP::optim_var_libMesh2local(NumericVector<lm_Number>& state, NumericVector<lm_Number>& control, Number* plocalDest)
{ // creates [local_state local_control]
  DBG_PRINT( "LibMeshPDENLP::optim_var_libMesh2local called" )
  int iVal = 0;
  int iOffset = 0;
  PetscVector<lm_Number>* p_vec = dynamic_cast<PetscVector<lm_Number>*>(&state);
  Vec v = p_vec->vec();
  PetscScalar *Vals;
  VecGetArray(v,&Vals);
  for ( iVal=0; iVal<=(state_last_-state_first_); ++iVal, ++iOffset)
    *(plocalDest+iOffset) = Vals[iVal];
  VecRestoreArray(v,&Vals);
  p_vec = dynamic_cast<PetscVector<lm_Number>*>(&control);
  v = p_vec->vec();
  VecGetArray(v,&Vals);
  for ( iVal=0; iVal<=(control_last_-control_first_); ++iVal, ++iOffset)
    *(plocalDest+iOffset) = Vals[iVal];
  VecRestoreArray(v,&Vals);
  DBG_PRINT( "LibMeshPDENLP::optim_var_libMesh2local finished" );
}

void LibMeshPDENLP::optim_var_global2libMesh(const Number* pglobal, NumericVector<lm_Number>& state, NumericVector<lm_Number>& control)
{ // [state_proc0 contrl_proc0 state_proc1 contrl_proc1 .... state_procN contrl_procNn] -> local_state, local_cntrl
  DBG_PRINT( "LibMeshPDENLP::optim_var_global2libMesh called" )
  int iVal = 0;
  int iOffset = n_first_;
  PetscVector<lm_Number>* p_vec = dynamic_cast<PetscVector<lm_Number>*>(&state);
  Vec v = p_vec->vec();
  PetscScalar *Vals;
  VecGetArray(v,&Vals);
  for ( iVal=0; iVal<=(state_last_-state_first_); ++iVal, ++iOffset)
    Vals[iVal] = *(pglobal+iOffset);
  VecRestoreArray(v,&Vals);
  p_vec = dynamic_cast<PetscVector<lm_Number>*>(&control);
  v = p_vec->vec();
  VecGetArray(v,&Vals);
  for ( iVal=0; iVal<=(control_last_-control_first_); ++iVal, ++iOffset)
    Vals[iVal] = *(pglobal+iOffset);
  VecRestoreArray(v,&Vals);
  assert(iOffset == n_last_+1);
  DBG_PRINT( "LibMeshPDENLP::optim_var_global2libMesh finished" );
}

bool LibMeshPDENLP::get_bounds_info(Index num_proc, Index proc_id,
                                    Index n, Index n_first, Index n_last,
                                    Number* x_l_part, Number* x_u_part,
                                    Index m, Index m_first, Index m_last,
                                    Number* g_l_part, Number* g_u_part)
{
  DBG_PRINT("LibMeshPDENLP::get_bounds_info called")
  AutoPtr<  NumericVector< lm_Number > > state_l = libmeshPDE_->getStateVector().clone();
  AutoPtr<  NumericVector< lm_Number > > state_u = libmeshPDE_->getStateVector().clone();
  AutoPtr<  NumericVector< lm_Number > > control_l = libmeshPDE_->getControlVector().clone();
  AutoPtr<  NumericVector< lm_Number > > control_u = libmeshPDE_->getControlVector().clone();
  AutoPtr<  NumericVector< lm_Number > > aux_constr_l = libmeshPDE_->getAuxConstrVector().clone();
  AutoPtr<  NumericVector< lm_Number > > aux_constr_u = libmeshPDE_->getAuxConstrVector().clone();

  libmeshPDE_->get_bounds(*state_l,*state_u,*control_l,*control_u,*aux_constr_l,*aux_constr_u);
  optim_var_libMesh2local(*state_l, *control_l, x_l_part);
  optim_var_libMesh2local(*state_u, *control_u, x_u_part);
  libMesh::NumericVector<Number>* pResidual;
  libmeshPDE_->calcPDE_residual(pResidual);
  int iVal=0;
  for (iVal=0;iVal<=(pde_last_-pde_first_);++iVal) {
    g_l_part[iVal] = 0.0;
    g_u_part[iVal] = 0.0;
  }
  int iAuxVal=0;
  PetscVector<lm_Number>* lm_vec_aux_l = dynamic_cast<PetscVector<lm_Number>*>(aux_constr_l.get());
  PetscVector<lm_Number>* lm_vec_aux_u = dynamic_cast<PetscVector<lm_Number>*>(aux_constr_u.get());
  Vec vec_aux_l = lm_vec_aux_l->vec();
  Vec vec_aux_u = lm_vec_aux_u->vec();
  PetscScalar *Vals_l(NULL),*Vals_u(NULL);
  VecGetArray(vec_aux_l,&Vals_l);
  VecGetArray(vec_aux_u,&Vals_u);
  for (iAuxVal=0;iAuxVal<=(aux_last_-aux_first_);++iVal,++iAuxVal) {
    g_l_part[iVal] = Vals_l[iAuxVal];
    g_u_part[iVal] = Vals_u[iAuxVal];
  }
  VecRestoreArray(vec_aux_l,&Vals_l);
  VecRestoreArray(vec_aux_u,&Vals_u);
  DBG_PRINT( "LibMeshPDENLP::get_bounds_info finished" );
  return true;
}


bool
LibMeshPDENLP::get_starting_point(Index num_proc, Index proc_id,
                                  Index n, Index n_first, Index n_last,
                                  bool init_x, Number* x_part,
                                  bool init_z, Number* z_L_part, Number* z_U_part,
                                  Index m, Index m_first, Index m_last,
                                  bool init_lambda, Number* lambda_part)
{
  DBG_PRINT( "LibMeshPDENLP::get_starting_point called" );
  assert(!init_z && !init_lambda);

  AutoPtr<  NumericVector< lm_Number > > state = libmeshPDE_->getStateVector().clone();
  AutoPtr<  NumericVector< lm_Number > > control = libmeshPDE_->getControlVector().clone();

  libMesh::NumericVector<libMesh::Number>* state_lb_mults(NULL);
  libMesh::NumericVector<libMesh::Number>* state_ub_mults(NULL);
  libMesh::NumericVector<libMesh::Number>* control_lb_mults(NULL);
  libMesh::NumericVector<libMesh::Number>* control_ub_mults(NULL);
  libMesh::NumericVector<libMesh::Number>* pde_residual_mults(NULL);
  libMesh::NumericVector<libMesh::Number>* aux_constr_mults(NULL);

  libmeshPDE_->get_starting_point(*state, *control, init_x, state_lb_mults,
				  state_ub_mults, control_lb_mults,
				  control_ub_mults, init_lambda,
				  pde_residual_mults, aux_constr_mults);
  optim_var_libMesh2local(*state, *control, x_part);
  if (init_z) {
    optim_var_libMesh2local(*state_lb_mults, *control_lb_mults, z_L_part);
    optim_var_libMesh2local(*state_ub_mults, *control_ub_mults, z_U_part);
  }
  if (init_lambda) {
    PetscVector<lm_Number>* p_vec =
      dynamic_cast<PetscVector<lm_Number>*>(pde_residual_mults);
    Vec v = p_vec->vec();
    PetscScalar *Vals = NULL;
    VecGetArray(v,&Vals);
    int iVal;
    for (iVal=0; iVal<=(pde_last_-pde_first_); ++iVal) {
      lambda_part[iVal] = Vals[iVal];
    }
    VecRestoreArray(v,&Vals);

    p_vec = dynamic_cast<PetscVector<lm_Number>*>(aux_constr_mults);
    v = p_vec->vec();
    VecGetArray(v,&Vals);
    int iAuxVal;
    for (int iAuxVal=0;iAuxVal<=(aux_last_-aux_first_);++iAuxVal,++iVal) {
      lambda_part[iVal] = Vals[iAuxVal];
    }
    VecRestoreArray(v,&Vals);
  }

  DBG_PRINT( "LibMeshPDENLP::get_starting_point finished" );
  return true;
}

void LibMeshPDENLP::update_x(const Number* x)
{
  DBG_PRINT( "LibMeshPDENLP::update_x called" );
  optim_var_global2libMesh(x, libmeshPDE_->getStateVector(), libmeshPDE_->getControlVector());

  // Need to make sure that local values are up-to-date
  libmeshPDE_->update();
  DBG_PRINT( "LibMeshPDENLP::update_x finished" );
}

bool
LibMeshPDENLP::eval_f(Index num_proc, Index proc_id,
                      Index n, Index n_first, Index n_last,
                      const Number* x, bool new_x, Number& obj_value)
{
  DBG_PRINT( "LibMeshPDENLP::eval_f called" );
  if (new_x) {
    update_x(x);
  }

  libmeshPDE_->calc_objective_part(obj_value);

  DBG_PRINT( "LibMeshPDENLP::eval_f finished" );
  return true;
}


bool
LibMeshPDENLP::eval_grad_f(Index num_proc, Index proc_id,
                           Index n,  Index n_first, Index n_last,
                           const Number* x, bool new_x,
                           Number* grad_f_part)
{
  DBG_PRINT( "LibMeshPDENLP::eval_grad_f called" );
  if (new_x) {
    update_x(x);
  }

  AutoPtr<  NumericVector< lm_Number > > state = libmeshPDE_->getStateVector().clone();
  AutoPtr<  NumericVector< lm_Number > > control = libmeshPDE_->getControlVector().clone();
  libmeshPDE_->calc_objective_gradient(*state, *control);
  optim_var_libMesh2local(*state, *control, grad_f_part);

  DBG_PRINT( "LibMeshPDENLP::eval_grad_f called" );
  return true;
}

bool
LibMeshPDENLP::eval_g(Index num_proc, Index proc_id,
                      Index n, const Number* x, bool new_x,
                      Index m, Index m_first, Index m_last,
                      Number* g_part)
{
  // g = [Aux_Proc1 PDE_Proc_1 Aux_Proc2 PDE_Proc_2 .... Aux_ProcN PDE_ProcN]
  DBG_PRINT( "LibMeshPDENLP::eval_g called" );
  if (new_x) {
    update_x(x);
  }

  // PDE resuduals
  NumericVector<lm_Number>* residual;
  libmeshPDE_->calcPDE_residual(residual);
  PetscVector<lm_Number>* p_vec = dynamic_cast<PetscVector<lm_Number>*>(residual);
  Vec v = p_vec->vec();
  PetscScalar *Vals = NULL;
  VecGetArray(v,&Vals);
  int iVal;
  for (iVal=0; iVal<=(pde_last_-pde_first_); ++iVal) {
    g_part[iVal] = Vals[iVal];
  }
  VecRestoreArray(v,&Vals);

  // auxiliary constraints
  libMesh::NumericVector<lm_Number>* constr;
  libmeshPDE_->calcAux_constr(constr);
  p_vec = dynamic_cast<PetscVector<lm_Number>*>(constr);
  v = p_vec->vec();
  VecGetArray(v,&Vals);
  int iAuxVal(0);
  for (int iAuxVal=0;iAuxVal<=(aux_last_-aux_first_);++iAuxVal,++iVal) {
    g_part[iVal] = Vals[iAuxVal];
  }
  VecRestoreArray(v,&Vals);

  DBG_PRINT( "LibMeshPDENLP::eval_g finished" );

  return true;
}

bool
LibMeshPDENLP::eval_jac_g(Index num_proc, Index proc_id,
                          Index n, const Number* x, bool new_x,
                          Index m, Index m_first, Index m_last,
                          Index nele_jac_part, Index* iRow_part,
                          Index *jCol_part, Number* values_part)
{
  DBG_PRINT( "LibMeshPDENLP::eval_jac_g called" );
  if (new_x) {
    update_x(x);
  }
  if (values_part==NULL) {
    libmeshPDE_->calc_type_ = LibMeshPDEBase::StructureOnly;
    assert(iRow_part && jCol_part);
  }
  else {
    libmeshPDE_->calc_type_ = LibMeshPDEBase::Values;
    assert(!iRow_part && !jCol_part);
  }

  SparseMatrix<lm_Number>* jac_state;
  SparseMatrix<lm_Number>* jac_control;
  libmeshPDE_->calcPDE_jacobians(jac_state, jac_control);
  PetscMatrix<lm_Number>* p_jac_state = dynamic_cast<PetscMatrix<lm_Number>*>(jac_state);
  Mat mat_state = p_jac_state->mat();
  PetscMatrix<lm_Number>* p_jac_control = dynamic_cast<PetscMatrix<lm_Number>*>(jac_control);
  Mat mat_control = p_jac_control->mat();
  SparseMatrix<lm_Number>* jac_aux_state;
  SparseMatrix<lm_Number>* jac_aux_control;
  libmeshPDE_->calcAux_jacobians(jac_aux_state, jac_aux_control);
  PetscMatrix<lm_Number>* p_aux_jac_state = dynamic_cast<PetscMatrix<lm_Number>*>(jac_aux_state);
  Mat mat_aux_state = p_aux_jac_state->mat();
  PetscMatrix<lm_Number>* p_aux_jac_control = dynamic_cast<PetscMatrix<lm_Number>*>(jac_aux_control);
  Mat mat_aux_control = p_aux_jac_control->mat();

  int iOffset = 0;
  if (values_part==NULL) {
    int *p_state_first = new int[num_proc];
    int *p_state_last = new int[num_proc];
    int *p_control_first = new int[num_proc];
    int *p_control_last = new int[num_proc];
    int *p_optvar_first = new int[num_proc];
    int *p_loc_control_offset = new int[num_proc];
    MPI_Allgather(&state_first_,1,MPI_INT,p_state_first,1,MPI_INT,MPI_COMM_WORLD);
    MPI_Allgather(&state_last_ ,1,MPI_INT,p_state_last ,1,MPI_INT,MPI_COMM_WORLD);
    MPI_Allgather(&control_first_,1,MPI_INT,p_control_first,1,MPI_INT,MPI_COMM_WORLD);
    MPI_Allgather(&control_last_ ,1,MPI_INT,p_control_last ,1,MPI_INT,MPI_COMM_WORLD);

    DBG_PRINT( "p_state_first[0]:"<<p_state_first[0] );
    DBG_PRINT( "p_state_last[0]:"<<p_state_last[0] );
    DBG_PRINT( "p_control_first[0]:"<<p_control_first[0] );
    DBG_PRINT( "p_control_last[0]:"<<p_control_last[0] );

    int i_var(0);
    for (int i_proc=0;i_proc<num_proc;++i_proc) {
      p_optvar_first[i_proc] = i_var;
      i_var += (p_state_last[i_proc]-p_state_first[i_proc]+1) +
               (p_control_last[i_proc]-p_control_first[i_proc]+1);
      p_loc_control_offset[i_proc] = p_state_last[i_proc] - p_state_first[i_proc] + 1;
    }

    int prev_proc(0);
    // PDE
    const int n_first_control = n_first_ + (state_last_-state_first_+1);
    const int n_first_state = n_first_;
    PetscInt irow;
    for (irow = pde_first_; irow<=pde_last_; ++irow) {
      PetscInt ncols;
      const PetscInt* cols;
      // state
      MatGetRow(mat_state, irow, &ncols, &cols, PETSC_NULL);
      for (int i=0; i<ncols; ++i) {
        iRow_part[iOffset] = (irow-pde_first_);
        unsigned int i_var_glob = GlobStateControlIdx2OptVarIdx(cols[i],num_proc,p_state_first,p_state_last,
                                  p_optvar_first,NULL,prev_proc);
        jCol_part[iOffset] = i_var_glob;//n_first_state + (cols[i]-state_first_);
        iOffset++;
      }
      MatRestoreRow(mat_state,irow,&ncols,&cols,PETSC_NULL);
      // control
      MatGetRow(mat_control, irow, &ncols, &cols, PETSC_NULL);
      for (int i=0; i<ncols; ++i) {
        iRow_part[iOffset] = (irow-pde_first_);
        unsigned int i_var_glob = GlobStateControlIdx2OptVarIdx(cols[i],num_proc,
                                  p_control_first,p_control_last,p_optvar_first,p_loc_control_offset,prev_proc);
        jCol_part[iOffset] = i_var_glob;//n_first_control + (cols[i]-control_first_);
        iOffset++;
      }
      MatRestoreRow(mat_control,irow,&ncols,&cols,PETSC_NULL);
    }
    // auxiliary constraints
    int iAuxConstr(0);
    PetscInt msz,nsz;
    MatGetSize(mat_aux_state,&msz,&nsz);
    const int m_loc_aux_offset = pde_last_-pde_first_+1;
    for (irow=aux_first_; irow<=aux_last_; ++irow) {
      PetscInt ncols;
      const PetscInt* cols;
      // state
      MatGetRow(mat_aux_state, irow, &ncols, &cols, PETSC_NULL);
      for (int i=0; i<ncols; ++i) {
        iRow_part[iOffset] = m_loc_aux_offset + (irow-aux_first_);
        unsigned int i_var_glob = GlobStateControlIdx2OptVarIdx(cols[i],num_proc,p_state_first,p_state_last,
                                  p_optvar_first,NULL,prev_proc);
        jCol_part[iOffset] = i_var_glob;//n_first_state + (cols[i]-state_first_);
        iOffset++;
      }
      MatRestoreRow(mat_aux_state,irow,&ncols,&cols,PETSC_NULL);
      // control
      MatGetRow(mat_aux_control, irow, &ncols, &cols, PETSC_NULL);
      for (int i=0; i<ncols; ++i) {
        iRow_part[iOffset] = m_loc_aux_offset + (irow-aux_first_);
        unsigned int i_var_glob = GlobStateControlIdx2OptVarIdx(cols[i],num_proc,
                                  p_control_first,p_control_last,p_optvar_first,p_loc_control_offset,prev_proc);
        jCol_part[iOffset] = i_var_glob;//n_first_control + (cols[i]-control_first_);
        iOffset++;
      }
      MatRestoreRow(mat_aux_control,irow,&ncols,&cols,PETSC_NULL);
    }
    delete [] p_state_first;
    delete [] p_state_last;
    delete [] p_control_first;
    delete [] p_control_last;
    delete [] p_optvar_first;
    delete [] p_loc_control_offset;
  }
  else {
    // PDE constraints
    PetscInt irow;
    for (irow = pde_first_; irow<=pde_last_; ++irow) {
      PetscInt ncols;
      const PetscScalar* vals;
      // state
      MatGetRow(mat_state, irow, &ncols, PETSC_NULL, &vals);
      for (int i=0; i<ncols; ++i) {
        values_part[iOffset] = vals[i];
        iOffset++;
      }
      MatRestoreRow(mat_state,irow,&ncols,PETSC_NULL,&vals);
      // control
      MatGetRow(mat_control, irow, &ncols, PETSC_NULL, &vals);
      for (int i=0; i<ncols; ++i) {
        values_part[iOffset] = vals[i];
        iOffset++;
      }
      MatRestoreRow(mat_control,irow,&ncols,PETSC_NULL,&vals);
    }

    // auxiliary constraints
    int iAuxConstr(0);
    for (irow=aux_first_; irow<=aux_last_; ++irow) {
      PetscInt ncols;
      const PetscScalar* vals;
      // state
      MatGetRow(mat_aux_state, irow, &ncols, PETSC_NULL, &vals);
      for (int i=0; i<ncols; ++i) {
        values_part[iOffset] = vals[i];
        iOffset++;
      }
      MatRestoreRow(mat_aux_state,irow,&ncols,PETSC_NULL,&vals);
      // control
      MatGetRow(mat_aux_control, irow, &ncols, PETSC_NULL, &vals);
      for (int i=0; i<ncols; ++i) {
        values_part[iOffset] = vals[i];
        iOffset++;
      }
      MatRestoreRow(mat_aux_control,irow,&ncols,PETSC_NULL,&vals);
    }
  }
  if (nele_jac_part != iOffset) {
    DBG_PRINT( "n_ele_jac_part=" << nele_jac_part << "!=" << iOffset << "=iOffset" );
    assert(nele_jac_part == iOffset);
  }
  DBG_PRINT( "LibMeshPDENLP::eval_jac_g finished" );
  return true;
}

int LibMeshPDENLP::GlobStateControlIdx2OptVarIdx(int i_var, int num_of_procs,
    const int* proc_first_arr, const int* proc_last_arr,
    const int* proc_var_first, const int* proc_loc_offset, int & prev_proc)
{
  // start with previous proc, in hope to accelerate search, if distribution is contigious
  int rv;
  for (int i_proc=prev_proc;i_proc<num_of_procs; ++i_proc) {
    if ( (proc_first_arr[i_proc]<=i_var)
         && (proc_last_arr[i_proc]>=i_var) ) {
      prev_proc = i_proc;
      rv = proc_var_first[i_proc] + (i_var-proc_first_arr[i_proc]);
      if (proc_loc_offset)
        rv += proc_loc_offset[i_proc];
      return rv;
    }
  }

  for (int i_proc=0;i_proc<prev_proc; ++i_proc) {
    if ( (proc_first_arr[i_proc]<=i_var)
         && (proc_last_arr[i_proc]>=i_var) ) {
      prev_proc = i_proc;
      rv = proc_var_first[i_proc] + (i_var-proc_first_arr[i_proc]);
      if (proc_loc_offset)
        rv += proc_loc_offset[i_proc];
      return rv;
    }
  }

  int sz;
  MPI_Comm_size(MPI_COMM_WORLD,&sz);
  DBG_PRINT( "LibMeshPDENLP::GlobStateControlIdx2OptVarIdx(" << i_var << ")=? (no proc found)" );
  std::cout << "proc_first_arr:";
  for (int i=0;i<sz;i++) std::cout << ", " << proc_first_arr[i];
  std::cout << std::endl;
  std::cout << "proc_last_arr:";
  for (int i=0;i<sz;i++) std::cout << ", " << proc_last_arr[i];
  std::cout << std::endl;

  assert(false);
}

bool
LibMeshPDENLP::eval_h(Index num_proc, Index proc_id,
                      Index n, Index n_first, Index n_last,
                      const Number* x, bool new_x, Number obj_factor,
                      Index m, Index m_first, Index m_last,
                      const Number* lambda,
                      bool new_lambda, Index nele_hess_part,
                      Index* iRow_part, Index* jCol_part,
                      Number* values_part)
{
  DBG_PRINT( "LibMeshPDENLP::eval_h called" );
  if (new_x) {
    update_x(x);
  }
  if (values_part==NULL) {
    libmeshPDE_->calc_type_ = LibMeshPDEBase::StructureOnly;
    assert(iRow_part && jCol_part);
  }
  else {
    libmeshPDE_->calc_type_ = LibMeshPDEBase::Values;
    assert(!iRow_part && !jCol_part);
  }

  libMesh::DenseVector<libMesh::Number> lambda_pde;
  libMesh::DenseVector<lm_Number> lambda_aux;
  int m_pde = libmeshPDE_->getPDEResVector().size();
  int m_aux = libmeshPDE_->getAuxConstrVector().size();
  lambda_pde.resize(m_pde);
  lambda_aux.resize(m_aux);

  int m_aux_first =m_first+(pde_last_-pde_first_+1);
  if (NULL==values_part) {
    for (int i_constr=0;i_constr<m_pde;++i_constr) {
      lambda_pde(i_constr) = 1.0;
    }
    for ( int i_constr=0;i_constr<m_aux;++i_constr) {
      lambda_aux(i_constr) = 1.0;
    }
  }
  else {
    int *p_pde_first = new int[num_proc];
    int *p_pde_last = new int[num_proc];
    int *p_aux_first = new int[num_proc];
    int *p_aux_last = new int[num_proc];
    MPI_Allgather(&pde_first_,1,MPI_INT,p_pde_first,1,MPI_INT,MPI_COMM_WORLD);
    MPI_Allgather(&pde_last_ ,1,MPI_INT,p_pde_last ,1,MPI_INT,MPI_COMM_WORLD);
    MPI_Allgather(&aux_first_,1,MPI_INT,p_aux_first,1,MPI_INT,MPI_COMM_WORLD);
    MPI_Allgather(&aux_last_ ,1,MPI_INT,p_aux_last ,1,MPI_INT,MPI_COMM_WORLD);

    int iVal(0);
    for (int i_proc=0;i_proc<num_proc;i_proc++) {
      for (int i_constr=0;i_constr<=(p_pde_last[i_proc]-p_pde_first[i_proc]);++i_constr) {
        lambda_pde(i_constr+p_pde_first[i_proc]) = lambda[iVal++];
      }
      for (int i_constr=0;i_constr<=(p_aux_last[i_proc]-p_aux_first[i_proc]);++i_constr) {
        lambda_aux(i_constr+p_aux_first[i_proc]) = lambda[iVal++];
      }
    }
    delete[] p_pde_first;
    delete[] p_pde_last;
    delete[] p_aux_first;
    delete[] p_aux_last;
  }

  libMesh::SparseMatrix<Number> *Hcc, *Hcs, *Hss;
  libmeshPDE_->calc_hessians(obj_factor, lambda_pde, lambda_aux, Hcc, Hcs, Hss);
  unsigned int n_control = Hcc->m();
  unsigned int n_state = Hss->m();
  Mat mat_hcc = (dynamic_cast<PetscMatrix<lm_Number>*>(Hcc))->mat();
  Mat mat_hcs = (dynamic_cast<PetscMatrix<lm_Number>*>(Hcs))->mat();
  Mat mat_hss = (dynamic_cast<PetscMatrix<lm_Number>*>(Hss))->mat();
  //MatView(mat_hss,PETSC_VIEWER_STDOUT_WORLD);
  PetscInt ncols;
  const PetscInt *cols;
  const PetscScalar* vals;
  unsigned long iOffset(0);
  if (NULL==values_part) {
    int *p_state_first = new int[num_proc];
    int *p_state_last = new int[num_proc];
    int *p_control_first = new int[num_proc];
    int *p_control_last = new int[num_proc];
    int *p_optvar_first = new int[num_proc];
    int *p_loc_control_offset = new int[num_proc];
    MPI_Allgather(&state_first_,1,MPI_INT,p_state_first,1,MPI_INT,MPI_COMM_WORLD);
    MPI_Allgather(&state_last_ ,1,MPI_INT,p_state_last ,1,MPI_INT,MPI_COMM_WORLD);
    MPI_Allgather(&control_first_,1,MPI_INT,p_control_first,1,MPI_INT,MPI_COMM_WORLD);
    MPI_Allgather(&control_last_ ,1,MPI_INT,p_control_last ,1,MPI_INT,MPI_COMM_WORLD);
    int i_var(0);
    for (int i_proc=0;i_proc<num_proc;++i_proc) {
      p_optvar_first[i_proc] = i_var;
      i_var += (p_state_last[i_proc]-p_state_first[i_proc]+1) +
               (p_control_last[i_proc]-p_control_first[i_proc]+1);
      p_loc_control_offset[i_proc] = p_state_last[i_proc] - p_state_first[i_proc] + 1;
    }

    int prev_proc_row(0);
    int prev_proc_col(0);
    int i_loc_state_offset = state_last_-state_first_+1;
    // state_state
    for (unsigned int i_row=0;i_row<n_state;++i_row) {
      MatGetRow(mat_hss, i_row, &ncols, &cols, PETSC_NULL);
      unsigned int i_var_row = GlobStateControlIdx2OptVarIdx(i_row,num_proc,
                               p_state_first,p_state_last,p_optvar_first,NULL,prev_proc_row);
      for (unsigned int i_col=0; i_col<ncols; ++i_col) {
        unsigned int i_var_col = GlobStateControlIdx2OptVarIdx(cols[i_col],num_proc,
                                 p_state_first,p_state_last,p_optvar_first,NULL,prev_proc_col);
        iRow_part[iOffset] = i_var_row;
        jCol_part[iOffset] = i_var_col;
        ++iOffset;
      }
      MatRestoreRow(mat_hss, i_row, &ncols, &cols, PETSC_NULL);
    }

    // control control
    prev_proc_row=0;
    prev_proc_col=0;
    for (unsigned int i_row=0;i_row<n_control;++i_row) {
      MatGetRow(mat_hcc, i_row, &ncols, &cols, PETSC_NULL);
      unsigned int i_var_row = GlobStateControlIdx2OptVarIdx(i_row,num_proc,
                               p_control_first,p_control_last,p_optvar_first,p_loc_control_offset,prev_proc_row);
      for (unsigned int i_col=0; i_col<ncols; ++i_col) {
        unsigned int i_var_col = GlobStateControlIdx2OptVarIdx(cols[i_col],num_proc,p_control_first,p_control_last,
                                 p_optvar_first,p_loc_control_offset,prev_proc_col);
        iRow_part[iOffset] = i_var_row;
        jCol_part[iOffset] = i_var_col;
        ++iOffset;
      }
      MatRestoreRow(mat_hcc, i_row, &ncols, &cols, PETSC_NULL);
    }

    // control state
    prev_proc_row=0;
    prev_proc_col=0;
    for (unsigned int i_row=0;i_row<n_control;++i_row) {
      MatGetRow(mat_hcs, i_row, &ncols, &cols, PETSC_NULL);
      unsigned int i_var_row = GlobStateControlIdx2OptVarIdx(i_row,num_proc,
                               p_control_first,p_control_last,p_optvar_first,p_loc_control_offset,prev_proc_row);
      for (unsigned int i_col=0; i_col<ncols; ++i_col) {
        unsigned int i_var_col = GlobStateControlIdx2OptVarIdx(cols[i_col],num_proc,p_state_first,p_state_last,
                                 p_optvar_first,NULL,prev_proc_col);
        iRow_part[iOffset] = i_var_row;
        jCol_part[iOffset] = i_var_col;
        ++iOffset;
      }
      MatRestoreRow(mat_hcs, i_row, &ncols, &cols, PETSC_NULL);
    }

    assert(iOffset==nele_hess_part);

    delete[] p_state_first;
    delete[] p_state_last;
    delete[] p_control_first;
    delete[] p_control_last;
    delete[] p_optvar_first;
    delete[] p_loc_control_offset;
  }
  else {
    // state_state
    for (unsigned int i_row=0;i_row<n_state;++i_row) {
      MatGetRow(mat_hss, i_row, &ncols, PETSC_NULL, &vals);
      for (unsigned int i_col=0; i_col<ncols; ++i_col) {
        values_part[iOffset] = vals[i_col];
        ++iOffset;
      }
      MatRestoreRow(mat_hss, i_row, &ncols, &cols, PETSC_NULL);
    }

    // control control
    for (unsigned int i_row=0;i_row<n_control;++i_row) {
      MatGetRow(mat_hcc, i_row, &ncols, &cols, PETSC_NULL);
      for (unsigned int i_col=0; i_col<ncols; ++i_col) {
        values_part[iOffset] = vals[i_col];
        ++iOffset;
      }
      MatRestoreRow(mat_hcc, i_row, &ncols, &cols, PETSC_NULL);
    }

    // control state
    for (unsigned int i_row=0;i_row<n_control;++i_row) {
      MatGetRow(mat_hcs, i_row, &ncols, &cols, PETSC_NULL);
      for (unsigned int i_col=0; i_col<ncols; ++i_col) {
        values_part[iOffset] = vals[i_col];
        ++iOffset;
      }
      MatRestoreRow(mat_hcs, i_row, &ncols, &cols, PETSC_NULL);
    }
    assert(iOffset==nele_hess_part);
  }

  DBG_PRINT( "LibMeshPDENLP::eval_h finished" );
  return true;
}

void LibMeshPDENLP::finalize_solution(SolverReturn status,
                                      Index n, const Number* x, const Number* z_L, const Number* z_U,
                                      Index m, const Number* g, const Number* lambda,
                                      Number obj_value,
                                      const IpoptData* ip_data,
                                      IpoptCalculatedQuantities* ip_cq)
{
  update_x(x);
  libmeshPDE_->Write2File("Sol");

  libMesh::NumericVector<libMesh::Number>* lm_state_lb_mults;
  libMesh::NumericVector<libMesh::Number>* lm_state_ub_mults;
  libMesh::NumericVector<libMesh::Number>* lm_control_lb_mults;
  libMesh::NumericVector<libMesh::Number>* lm_control_ub_mults;
  libMesh::NumericVector<libMesh::Number>* lm_pde_residual_mults;
  libMesh::NumericVector<libMesh::Number>* lm_aux_constr_mults;

  libmeshPDE_->get_finalize_vectors(lm_state_lb_mults, lm_state_ub_mults,
				    lm_control_lb_mults, lm_control_ub_mults,
				    lm_pde_residual_mults, lm_aux_constr_mults);

  if (lm_state_lb_mults) {
    assert(lm_control_lb_mults);
    optim_var_global2libMesh(z_L, *lm_state_lb_mults, *lm_control_lb_mults);
  }
  if (lm_state_ub_mults) {
    assert(lm_control_ub_mults);
    optim_var_global2libMesh(z_U, *lm_state_ub_mults, *lm_control_ub_mults);
  }
  if (lm_pde_residual_mults) {
    PetscVector<lm_Number>* p_vec =
      dynamic_cast<PetscVector<lm_Number>*>(lm_pde_residual_mults);
    Vec v = p_vec->vec();
    PetscScalar *Vals = NULL;
    VecGetArray(v,&Vals);
    int iVal;
    int iLam = m_first_;
      for (iVal=0; iVal<=(pde_last_-pde_first_); ++iVal, ++iLam) {
      Vals[iVal] = lambda[iLam];
    }
    VecRestoreArray(v,&Vals);
  }
  if (lm_aux_constr_mults) {
    PetscVector<lm_Number>* p_vec =
      dynamic_cast<PetscVector<lm_Number>*>(lm_aux_constr_mults);
    Vec v = p_vec->vec();
    PetscScalar *Vals = NULL;
    VecGetArray(v,&Vals);
    int iVal;
    int iLam = m_first_ + (pde_last_-pde_first_)+1;
    for (iVal=0; iVal<=(aux_last_-aux_first_); ++iVal, ++iLam) {
      Vals[iVal] = lambda[iLam];
    }
    VecRestoreArray(v,&Vals);
  }

}

// Copyright (C) 2010 International Business Machines.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// $Id$
//
// Authors:  Johannes Huber, Andreas Waechter     IBM        2010-09-03
#include "mpi.h"
#include "IpLibMeshPDE.hpp"
#include "petsc.h"

#include "string_to_enum.h"
#include "fe.h"
#include "quadrature_gauss.h"
#include "dof_map.h"
#include "sparse_matrix.h"
#include "petsc_matrix.h"
#include "petscvec.h"
#include "numeric_vector.h"
#include "dense_matrix.h"
#include "dense_vector.h"
#include "elem.h"
#include "linear_implicit_system.h"
#include "nonlinear_implicit_system.h"
#include "vtk_io.h"
#include "exodusII_io.h"
#include "error_vector.h"
#include "kelly_error_estimator.h"
#include "fourth_error_estimators.h"
#include "patch_recovery_error_estimator.h"
#include "mesh_refinement.h"


#include <fstream>
#include <set>

#define SCALE_AUX_BOUNDS
#define USE_NEW_DIRICHLET
//#define PHI_IN_OBJECTIVE

int GetProcID();

//#define DBG_PRINT(s) {std::cout << GetProcID() << ":" << s << std::endl;}
#define DBG_PRINT(s) {}

using namespace libMesh;

LibMeshPDEBase::LibMeshPDEBase() :
    calc_type_(Values),
    simulation_mode_(false),
    jac_control_(NULL),
    jac_aux_state_(NULL),
    jac_aux_control_(NULL),
    hess_control_control_(NULL),
    hess_control_state_(NULL),
    hess_state_state_(NULL),
    lm_eqn_sys_(NULL),
    lm_control_vec_(NULL),
    lm_pde_residual_vec_(NULL),
    lm_aux_constr_vec_(NULL),
    lm_aux_constr_vec_low_bd_(NULL),
    lm_control_lb_mults_(NULL),
    lm_control_ub_mults_(NULL),
    lm_aux_constr_mults_(NULL),
    min_airflow(1.0),
    first_aux_constr_(0),
    node_id_for_phi0_(-1)  
{}

// clear all matrices and vectors, but not problem geomatry data
// used e.g. after refinement (since problem size changes)
void LibMeshPDEBase::clear_math_obj()
{
  AuxConstrBoundMarkerList_.clear();
  DetroySelfOwnedLibMeshPetscVector(lm_control_vec_);
  DetroySelfOwnedLibMeshPetscVector(lm_pde_residual_vec_);
  DetroySelfOwnedLibMeshPetscVector(lm_aux_constr_vec_);
  DetroySelfOwnedLibMeshPetscVector(lm_aux_constr_vec_low_bd_);

  /*
  DetroySelfOwnedLibMeshPetscVector(lm_state_lb_mults_);
  DetroySelfOwnedLibMeshPetscVector(lm_state_ub_mults_);
  DetroySelfOwnedLibMeshPetscVector(lm_control_lb_mults_);
  DetroySelfOwnedLibMeshPetscVector(lm_control_ub_mults_);
  DetroySelfOwnedLibMeshPetscVector(lm_pde_residual_mults_);
  DetroySelfOwnedLibMeshPetscVector(lm_aux_constr_mults_);
  */

  DetroySelfOwnedLibMeshPetscMatrix(jac_control_);
  DetroySelfOwnedLibMeshPetscMatrix(jac_aux_state_);
  DetroySelfOwnedLibMeshPetscMatrix(jac_aux_control_);
  DetroySelfOwnedLibMeshPetscMatrix(hess_control_control_);
  DetroySelfOwnedLibMeshPetscMatrix(hess_control_state_);
  DetroySelfOwnedLibMeshPetscMatrix(hess_state_state_);

  if (lm_eqn_sys_)
    delete lm_eqn_sys_;
  lm_eqn_sys_ = NULL;
}

LibMeshPDEBase::~LibMeshPDEBase()
{
  clear_math_obj();
}

void LibMeshPDEBase::InitProblemData(std::istream& is)
{
  DBG_PRINT("ibMeshPDEBase::InitProblemData called");
  PG_.ReadFromStream(is);
  PG_.CreateMesh(&mesh_);

  WriteNodeFile(mesh_, "MeshGen.node");
  WriteEleFile(mesh_, "MeshGen.ele");

//  RefineMesh();

//  WriteNodeFile(mesh_, "MeshGen1.node");
//  WriteEleFile(mesh_, "MeshGen1.ele");
  DBG_PRINT("ibMeshPDEBase::InitProblemData finished");
}

void LibMeshPDEBase::reinit()
{
  libMesh::LinearImplicitSystem* lm_sys;
  if (!lm_eqn_sys_) {
    // We are running this for the first time

    // initalize all members that store vectors and matrices
    clear_math_obj();
    lm_eqn_sys_ = new EquationSystems(mesh_);

    std::string order("FIRST");
    std::string family("LAGRANGE");

    lm_sys = &lm_eqn_sys_->add_system<LinearImplicitSystem>("PDE");
    lm_sys->add_variable("Phi", Utility::string_to_enum<Order>(order), Utility::string_to_enum<FEFamily>(family));
    lm_sys->attach_assemble_function(LibMeshPDEBase::assemble_Phi_PDE);
    lm_eqn_sys_->parameters.set<LibMeshPDEBase*>("LibMeshPDEBase") = this;
    lm_eqn_sys_->parameters.set<bool>("b_struct_only") =false;  // true -> 0/1 values

    // Add vectors that we want to prolong later
    lm_sys->add_vector("state_lb_mults");
    lm_sys->add_vector("state_ub_mults");
    lm_sys->add_vector("pde_residual_mults");
    lm_sys->add_vector("full_aux_constr_mults");

    lm_eqn_sys_->init(); // here the matrix and the vector for solution of the PDE is generated
  }
  else {
    DetroySelfOwnedLibMeshPetscVector(lm_aux_constr_vec_);
    DetroySelfOwnedLibMeshPetscVector(lm_aux_constr_vec_low_bd_);
    DetroySelfOwnedLibMeshPetscMatrix(jac_control_);
    DetroySelfOwnedLibMeshPetscMatrix(jac_aux_state_);
    DetroySelfOwnedLibMeshPetscMatrix(jac_aux_control_);
    DetroySelfOwnedLibMeshPetscMatrix(hess_control_control_);
    DetroySelfOwnedLibMeshPetscMatrix(hess_control_state_);
    DetroySelfOwnedLibMeshPetscMatrix(hess_state_state_);
    lm_sys = &lm_eqn_sys_->get_system<LinearImplicitSystem>("PDE");
  }

  // determine local and global ranges of state
  int n_state_global, n_state_local;
  int m_pde_constr_global, m_pde_constr_local;
  {
    //nstate = lm_sys->matrix->n();
    n_state_global = lm_sys->matrix->n();
    PetscMatrix<Number> *pPetscMat = dynamic_cast<PetscMatrix<Number>*>(lm_sys->matrix);
    Mat mat = pPetscMat->mat();
    int start, end;
    MatGetOwnershipRangeColumn(mat, &start, &end );
    n_state_local = end-start;
    //int mlocal = lm_sys->matrix->row_stop() - lm_sys->matrix->row_start();
    m_pde_constr_global = lm_sys->matrix->m();
    m_pde_constr_local = lm_sys->matrix->row_stop()-lm_sys->matrix->row_start();
  }

  // determine local and global ranges of control
  int n_control_global, n_control_local;
  {
    n_control_global = PG_._ParamIdx2BCParam.size();
    n_control_local = (0==GetProcID()) ? n_control_global : 0;
  }

  // create vectro and matrices for auxilliary constraints
  int m_aux_constr_global, m_aux_constr_local;
  std::list<Number> LocIneqFactList;
  InitAuxConstr(&m_aux_constr_local, &m_aux_constr_global, &LocIneqFactList);

  std::cout <<  GetProcID() << ": n_state_global: " << n_state_global <<
  ", n_state_local: " << n_state_local <<
  ", n_control_global: " << n_control_global <<
  ", n_control_local: " << n_control_local <<
  ", m_pde_constr_global: " << m_pde_constr_global <<
  ", m_pde_constr_local: " << m_pde_constr_local <<
  ", m_aux_constr_global: " << m_aux_constr_global <<
  ", m_aux_constr_local: " << m_aux_constr_local << std::endl;

  PetscErrorCode ierr;
  /// TODO: The following block nca probably be deleted (see tw blocks down)
  if (!lm_control_vec_) {
    // create Petsc and libmesh vector for controls
    //lm_control_vec_ = new PetscVector<Number>::PetscVector(n_control_global,n_control_local); does NOT work, later Petsc-calls never terminate (?)
    Vec petsc_vec;
    ierr = VecCreateMPI(PETSC_COMM_WORLD,n_control_local,PETSC_DETERMINE,&petsc_vec);
    CHKERRV(ierr);
    lm_control_vec_ = new PetscVector<Number>::PetscVector(petsc_vec);
    *lm_control_vec_ = 1.0;
    lm_control_vec_->close();
  }
  { // TODO: analyze nonzero structure more closely to pass tight numbers to init
    // create Jacobi matrix of PDE constraints for control variables part
    Mat petsc_mat;
    ierr = MatCreateMPIAIJ(PETSC_COMM_WORLD,m_pde_constr_local,n_control_local,PETSC_DETERMINE,PETSC_DETERMINE,64,PETSC_NULL,16,PETSC_NULL,&petsc_mat); // alloc 64 entries per row on diagonal block and 16 on the off-diagonal block
    CHKERRV(ierr);
    jac_control_ = new PetscMatrix<Number>::PetscMatrix(petsc_mat);
  }
  {
    Vec petsc_vec;
    ierr = VecCreateMPI(PETSC_COMM_WORLD,m_pde_constr_local,PETSC_DETERMINE,&petsc_vec);
    CHKERRV(ierr);
    lm_pde_residual_vec_ = new PetscVector<Number>::PetscVector(petsc_vec);
    lm_pde_residual_vec_->close();
  }
  //  lm_pde_residual_vec_ = new PetscVector<Number>::PetscVector(m_pde_constr_global,m_pde_constr_local);
  // compute actual Jacobian at Some(?) point to get distribution of rows and columns in parallel
  SparseMatrix<Number> *pde_jac_state, *pde_jac_control;
  calcPDE_jacobians(pde_jac_state, pde_jac_control); // required to have something initialized (check)
  {
    // constraint values
    Vec petsc_vec;
    ierr = VecCreateMPI(PETSC_COMM_WORLD,m_aux_constr_local,PETSC_DETERMINE,&petsc_vec);
    CHKERRV(ierr);
    lm_aux_constr_vec_ = new PetscVector<Number>::PetscVector(petsc_vec);
    lm_aux_constr_vec_->close();

    // constrains boundaries
    Vec petsc_vec1;
    ierr = VecDuplicate(petsc_vec, &petsc_vec1);
    CHKERRV(ierr);
    PetscScalar *Vals;
    VecGetArray(petsc_vec1,&Vals);
#ifdef SCALE_AUX_BOUNDS
    std::list<Number>::iterator it = LocIneqFactList.begin();
    std::cout << "LocIneqFactList.size=" << LocIneqFactList.size() << std::endl;
    for ( unsigned int iVal=0; iVal<LocIneqFactList.size(); ++iVal, ++it) {
      Vals[iVal] = min_airflow * (*it);
    }
#else
    for ( unsigned int iVal=0; iVal<LocIneqFactList.size(); ++iVal) {
      Vals[iVal] = min_airflow;
    }
#endif
    std::cout << "Checkpoint 1" << std::endl;
    VecRestoreArray(petsc_vec1,&Vals);
    lm_aux_constr_vec_low_bd_ = new PetscVector<Number>::PetscVector(petsc_vec1);
    lm_aux_constr_vec_low_bd_->close();
  }
  //lm_aux_constr_vec_ = new PetscVector<Number>::PetscVector(m_aux_constr_global,m_aux_constr_local);
  //*lm_aux_constr_vec_ = 0.0;
  {
    Mat petsc_mat;
    MatCreateMPIAIJ(PETSC_COMM_WORLD,m_aux_constr_local,n_state_local,m_aux_constr_global,n_state_global,64,PETSC_NULL,16,PETSC_NULL,&petsc_mat); // alloc 64 entries per row on diagonal block and 16 on the off-diagonal block
    //MatCreateSeqAIJ(PETSC_COMM_SELF,m_aux_constr_local,n_state_global,4,PETSC_NULL,&petsc_mat);
    jac_aux_state_ = new PetscMatrix<Number>::PetscMatrix(petsc_mat);
  }
  {
    Mat petsc_mat;
    MatCreateMPIAIJ(PETSC_COMM_WORLD,m_aux_constr_local,n_control_local,m_aux_constr_global,n_control_global,64,PETSC_NULL,16,PETSC_NULL,&petsc_mat); // alloc 64 entries per row on diagonal block and 16 on the off-diagonal block
    //MatCreateSeqAIJ(PETSC_COMM_SELF,m_aux_constr_local,n_control_global,4,PETSC_NULL,&petsc_mat);
    jac_aux_control_ = new PetscMatrix<Number>(petsc_mat);
  }
  {
    // each processor has hole matrix, but sets only part of local constraints
    Mat petsc_mat;
    //MatCreateMPIAIJ(PETSC_COMM_WORLD,n_state_global,n_state_global,n_state_global,n_state_global,8,PETSC_NULL,0,PETSC_NULL,&petsc_mat);
    MatCreateSeqAIJ(PETSC_COMM_SELF,n_state_global,n_state_global,8,PETSC_NULL,&petsc_mat);
    hess_state_state_ = new PetscMatrix<Number>(petsc_mat);
  }
  {
    // each processor has hole matrix, but sets only part of local constraints
    Mat petsc_mat;
    //MatCreateMPIAIJ(PETSC_COMM_WORLD,n_control_global,n_state_global,n_control_global,n_state_global,8,PETSC_NULL,0,PETSC_NULL,&petsc_mat);
    MatCreateSeqAIJ(PETSC_COMM_SELF,n_control_global,n_state_global,8,PETSC_NULL,&petsc_mat);
    hess_control_state_ = new PetscMatrix<Number>(petsc_mat);
  }
  {
    // each processor has hole matrix, but sets only part of local constraints
    Mat petsc_mat;
    //MatCreateMPIAIJ(PETSC_COMM_WORLD,n_control_global,n_control_global,n_control_global,n_control_global,8,PETSC_NULL,0,PETSC_NULL,&petsc_mat);
    MatCreateSeqAIJ(PETSC_COMM_SELF,n_control_global,n_control_global,8,PETSC_NULL,&petsc_mat);
    hess_control_control_ = new PetscMatrix<Number>(petsc_mat);
  }
#ifndef PHI_IN_OBJECTIVE
#ifdef EXHAUST_AS_CONTROL
  // determine a node ID of an interior node that is used to pinn down Phi
  if (GetProcID()==0) {
    const MeshBase& mesh = lm_eqn_sys_->get_mesh();
#if 0
    // first mark all local nodes in an element that is at a boundary
    std::set<unsigned int> outer_nodes;
    MeshBase::const_element_iterator itCurEl = mesh.active_local_elements_begin();
    const MeshBase::const_element_iterator itEndEl = mesh.active_local_elements_end();
    for ( ; itCurEl != itEndEl; ++itCurEl) {
      const Elem* CurElem = *itCurEl;
      for (unsigned int side=0; side<CurElem->n_sides(); side++) {
	if (CurElem->neighbor(side) == NULL) {
	  for (unsigned int i=0;i<CurElem->n_nodes();i++) {
	    printf("CurElem->node(i) = %d\n",CurElem->node(i));
	    outer_nodes.insert(CurElem->node(i));
	  }
	  break;
	}
      }
    }
    node_id_for_phi0_ = -1;
    // now go through the local nodes to find the first that is not in a boundary element
    MeshBase::const_node_iterator itCurNode = mesh.local_nodes_begin();
    const MeshBase::const_node_iterator itEndNode = mesh.local_nodes_end();
    for ( ; itCurNode != itEndNode; ++itCurNode) {
      const Node* CurNode = *itCurNode;
      const unsigned int nodeID = CurNode->id();
	    printf("nodeID = %d\n",nodeID);
      std::set<unsigned int>::const_iterator it = outer_nodes.find(nodeID);
      if (it == outer_nodes.end()) {
	node_id_for_phi0_ = nodeID;
	break;
      }
    }
    assert(node_id_for_phi0_ != -1);
#endif
    node_id_for_phi0_ = mesh.n_nodes()-1;
    printf("Node ID for pinning down Phi = %d\n", node_id_for_phi0_);
    // TODO: broadcast
  }
#endif
#endif
}

void LibMeshPDEBase::DetroySelfOwnedLibMeshPetscMatrix(SparseMatrix<Number>*& matrix)
{
  if (NULL==matrix)
    return;
  PetscMatrix<Number>* lm_petsc_matrix = dynamic_cast<PetscMatrix<Number>*>(matrix);
  assert(NULL!=lm_petsc_matrix);  // Problem, if not a libMesh::PetscMatrix<Number>
  Mat petsc_matrix = lm_petsc_matrix->mat();
  MatDestroy(petsc_matrix);
  delete matrix;
  matrix = NULL;
}

void LibMeshPDEBase::DetroySelfOwnedLibMeshPetscVector(NumericVector<Number>*& vector)
{
  if (NULL==vector)
    return;
  PetscVector<Number>* lm_petsc_vector = dynamic_cast<PetscVector<Number>*>(vector);
  assert(NULL!=lm_petsc_vector);  // Problem, if not a libMesh::PetscMatrix<Number>
  Vec petsc_vector = lm_petsc_vector->vec();
  //TODO: VecDestroy(petsc_vector); does not work
  VecDestroy(petsc_vector);
  delete vector;
  vector = NULL;
}

// call this to transfor control vector data (e.g., from Ipopt) to the PG_ data structure that is used in assembly function
void LibMeshPDEBase::ConvertControl2PGData()
{
//  DBG_PRINT("LibMeshPDEBase::ConvertControl2PGData called");
  // first create vector, where all control variables are locally accessable
  PetscVector<Number>* lm_vec = dynamic_cast<PetscVector<Number>*>(lm_control_vec_);
  assert(NULL!=lm_vec);  // Problem, if not a libMesh::PetscMatrix<Number>
  Vec vec_distr = lm_vec->vec();

  VecScatter vscat;
  Vec vec_gathrd;
  PetscErrorCode ierr;
  ierr = VecScatterCreateToAll(vec_distr,&vscat,&vec_gathrd); CHKERRV(ierr);
  ierr = VecScatterBegin(vscat,vec_distr,vec_gathrd,INSERT_VALUES,SCATTER_FORWARD); CHKERRV(ierr);
  ierr = VecScatterEnd(vscat,vec_distr,vec_gathrd,INSERT_VALUES,SCATTER_FORWARD); CHKERRV(ierr);
  ierr = VecScatterDestroy(vscat); CHKERRV(ierr);
  PetscInt sz;
  VecGetSize(vec_gathrd,&sz);
  PetscScalar *vals;
  ierr=VecGetArray(vec_gathrd,&vals);CHKERRV(ierr);
  int iControl = 0;
#ifndef EXHAUST_AS_CONTROL
  assert(PG_._AC.size()==(unsigned int)sz);
  for (int iAC=0;iAC<sz;++iAC) {
    int BoundaryMarker = PG_._AC[iAC].BoundaryMarker[0];
    BoundaryConditionSquarePhiRhs* pBC = dynamic_cast<BoundaryConditionSquarePhiRhs*>(PG_._BoundCond[BoundaryMarker]);
    assert(pBC);
    pBC->_PhiRhsScale = vals[iControl++];
  }
#else
  assert(PG_._AC.size()+PG_._Exh.size()==(unsigned int)sz);
  for (int iAC=0;iAC<PG_._AC.size();++iAC) {
    int BoundaryMarker = PG_._AC[iAC].BoundaryMarker[0];
    BoundaryConditionSquarePhiRhs* pBC = dynamic_cast<BoundaryConditionSquarePhiRhs*>(PG_._BoundCond[BoundaryMarker]);
    assert(pBC);
    pBC->_PhiRhsScale = vals[iControl++];
    //printf("DEBUG: iAC = %d BoundaryMarker=%d PhiRhs = %e\n", iAC, BoundaryMarker, PG_._BoundCond[BoundaryMarker].PhiRhs);
  }
  for (int iExh=0;iExh<PG_._Exh.size();++iExh) {
    int BoundaryMarker = PG_._Exh[iExh].BoundaryMarker[0];
    BoundaryConditionConstValues* pBC = dynamic_cast<BoundaryConditionConstValues*>(PG_._BoundCond[BoundaryMarker]);
    assert(pBC);
    pBC->_PhiRhs = vals[iControl++];
    //printf("DEBUG: iExh = %d BoundaryMarker=%d PhiRhs = %e\n", iExh, BoundaryMarker, PG_._BoundCond[BoundaryMarker].PhiRhs);
  }
#endif
  ierr=VecRestoreArray(vec_gathrd,&vals);CHKERRV(ierr);
  VecDestroy(vec_gathrd);
//  DBG_PRINT("LibMeshPDEBase::ConvertControl2PGData finished");
}

void LibMeshPDEBase::calc_objective_part(Number& Val)
{
  Val = 0.0;
  ConvertControl2PGData();
  if (GetProcID()==0) {
    for (unsigned int iAC=0;iAC<PG_._AC.size();iAC++) {
      int BoundaryMarker = PG_._AC[iAC].BoundaryMarker[0];
      BoundaryConditionSquarePhiRhs* pBC = dynamic_cast<BoundaryConditionSquarePhiRhs*>(PG_._BoundCond[BoundaryMarker]);
      assert(pBC);
#ifdef QUAD_OBJ_FUNC
      Val += pBC->_PhiRhsScale*pBC->_PhiRhsScale;
#else
      Val += pBC->_PhiRhsScale;
#endif //QUAD_OBJ_FUNC
    }
#ifdef EXHAUST_AS_CONTROL
# ifdef PHI_IN_OBJECTIVE
    // quadratic objective to make the first node in the first element zero
    // this is to pin down the potential
    const MeshBase& mesh = lm_eqn_sys_->get_mesh();
    LinearImplicitSystem& system = lm_eqn_sys_->get_system<LinearImplicitSystem>("PDE");
    const DofMap& dof_map = system.get_dof_map();
    MeshBase::const_element_iterator itCurEl = mesh.active_local_elements_begin();
    const Elem* CurElem = *itCurEl;
    std::vector<unsigned int> dof_indices;
    dof_map.dof_indices(CurElem, dof_indices);
    Number Phi_first_node = system.current_local_solution->el(dof_indices[0]);
    Val += Phi_first_node*Phi_first_node;
# endif
#endif
  }
}

void LibMeshPDEBase::calc_objective_gradient(libMesh::NumericVector<libMesh::Number>& grad_state,
    libMesh::NumericVector<libMesh::Number>& grad_control)
{
  DBG_PRINT("LibMeshPDEBase::calc_objective_gradient called");
  grad_state.zero();
  grad_control.zero();
  ConvertControl2PGData();
  if (GetProcID()==0) {
    for (unsigned int iAC=0;iAC<PG_._AC.size();iAC++) {
#ifdef QUAD_OBJ_FUNC
      int BoundaryMarker = PG_._AC[iAC].BoundaryMarker[0];
      grad_control.set(iAC,2.0*PG_._BoundCond[BoundaryMarker].PhiRhs);
#else
      grad_control.set(iAC,1.0);
#endif //QUAD_OBJ_FUNC
    }
#ifdef EXHAUST_AS_CONTROL
# ifdef PHI_IN_OBJECTIVE
    // quadratic objective to make the first node in the first element zero
    // this is to pin down the potential
    const MeshBase& mesh = lm_eqn_sys_->get_mesh();
    LinearImplicitSystem& system = lm_eqn_sys_->get_system<LinearImplicitSystem>("PDE");
    const DofMap& dof_map = system.get_dof_map();
    MeshBase::const_element_iterator itCurEl = mesh.active_local_elements_begin();
    const Elem* CurElem = *itCurEl;
    std::vector<unsigned int> dof_indices;
    dof_map.dof_indices(CurElem, dof_indices);
    Number Phi_first_node = system.current_local_solution->el(dof_indices[0]);
    grad_state.set(dof_indices[0], 2.*Phi_first_node);
# endif
#endif
  }
  grad_state.close();
  grad_control.close();
  DBG_PRINT("LibMeshPDEBase::calc_objective_gradient finished");
}

void LibMeshPDEBase::calcPDE_residual(libMesh::NumericVector<libMesh::Number>*& residual)
{
  DBG_PRINT("LibMeshPDEBase::calcPDE_residual called");
  ConvertControl2PGData();

  const double eps = 1e-8;
  lm_pde_residual_vec_->zero();
  std::vector<BoundaryConditionBase*>& BCs=PG_._BoundCond;
  const MeshBase& mesh = lm_eqn_sys_->get_mesh();
  const unsigned int dim = mesh.mesh_dimension();
  LinearImplicitSystem& system = lm_eqn_sys_->get_system<LinearImplicitSystem>("PDE");
  const DofMap& dof_map = system.get_dof_map();
  FEType fe_type = dof_map.variable_type(0);

  AutoPtr<FEBase> fe (FEBase::build(dim, fe_type)); // this object will by the current volume object, actualy holding values for current phi and dphi
  QGauss qrule (dim, FIFTH);
  fe->attach_quadrature_rule (&qrule);
  AutoPtr<FEBase> fe_face (FEBase::build(dim, fe_type)); // surface object
  QGauss qface(dim-1, FIFTH);
  fe_face->attach_quadrature_rule (&qface);
  const std::vector<Real>& JxW = fe->get_JxW();   // References to values hold by fe (i.e. the chang, ones fe changes
  const std::vector<std::vector<Real> >& phi = fe->get_phi();
  const std::vector<std::vector<RealGradient> >& dphi = fe->get_dphi();

  DenseVector<Number> ElemRes;
  std::vector<unsigned int> dof_indices; // mapping local index -> global index
  std::vector<Number> LocalSolution;   // solution values of current element
  MeshBase::const_element_iterator itCurEl = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator itEndEl = mesh.active_local_elements_end();
  for ( ; itCurEl != itEndEl; ++itCurEl) {
    const Elem* CurElem = *itCurEl;
    dof_map.dof_indices(CurElem, dof_indices);  // setup local->global mapping in form of array
    fe->reinit(CurElem); // reinit fe to fit the current element, i.e. recalulate JxW and dPhi
    ElemRes.resize(dof_indices.size());
    LocalSolution.resize(dof_indices.size());
    for (unsigned int iCurSol=0;iCurSol<LocalSolution.size();iCurSol++) {
      LocalSolution[iCurSol] = system.current_local_solution->el(dof_indices[iCurSol]);
    }
    // Point Center = (CurElem->point(0)+CurElem->point(1)+CurElem->point(2))/3.0;
    // std::cout << "Element at " << Center << std::endl;
    unsigned int qp,i,j;
    for (qp=0; qp<qrule.n_points(); qp++)
      for (i=0; i<phi.size(); i++)
        for (j=0; j<phi.size(); j++)
          ElemRes(i) += JxW[qp]*(dphi[i][qp]*dphi[j][qp])*LocalSolution[j];

    {
      for (unsigned int side=0; side<CurElem->n_sides(); side++) {
        if (CurElem->neighbor(side) == NULL) {
          // std::cout << side << std::endl;
          short int bc_id = mesh.boundary_info->boundary_id (CurElem,side);
          assert(bc_id!=BoundaryInfo::invalid_id);
          double NeumCoef, DiriCoef, Rhs;
          // std::cout << "bc_id=" << bc_id << std::endl;
          const std::vector<std::vector<Real> >&  phi_face = fe_face->get_phi();
          const std::vector<Real>& JxW_face = fe_face->get_JxW();
          fe_face->reinit(CurElem, side);
          const std::vector<Point>& q_point = fe_face->get_xyz(); // location of quadrature points in physical space
          /*Point Center;
          for(unsigned int iNode=0;iNode<CurElem->n_nodes();iNode++)
           if( CurElem->is_node_on_side(iNode,side) )
            Center += CurElem->point(iNode);
          Center = Center/dim;
          std::cout << "BC at " << Center << " and : ";
          std::cout << NeumCoef << " dPhi/dn = " << DiriCoef << " Phi + " << Rhs << std::endl;*/
          NeumCoef = BCs[bc_id]->PhiNeumannCoef(q_point[0]);  // PhiNeumCoef is constant anyway, so it doesn't matter at which pq we evaluate
          if ( fabs(NeumCoef)>eps) { // handle non-Dirichlet boundary conditions
            for (unsigned int qp=0; qp<qface.n_points(); qp++) {
              NeumCoef = BCs[bc_id]->PhiNeumannCoef(q_point[qp]);
              DiriCoef = BCs[bc_id]->PhiDirichletCoef(q_point[qp]);
              Rhs = BCs[bc_id]->PhiRhs(q_point[qp]);
              for (unsigned int i=0; i<phi_face.size(); i++)
                for (unsigned int j=0; j<phi_face.size(); j++)
                  ElemRes(i) -= (DiriCoef/NeumCoef)*JxW_face[qp]*phi_face[i][qp]*LocalSolution[j];
              for (unsigned int i=0; i<phi_face.size(); i++)
                ElemRes(i) -= (Rhs/NeumCoef)*JxW_face[qp]*phi_face[i][qp];
            }
          }
#ifdef USE_NEW_DIRICHLET
	        else {
	          // loop over all nodes and find the ones corresponding to this side
	          for (unsigned int i=0; i<CurElem->n_nodes(); i++) {
	            if (CurElem->is_node_on_side(i, side)) {
		            ElemRes(i) += BCs[bc_id]->PhiDirichletCoef(CurElem->node(i))*LocalSolution[i];
		            ElemRes(i) += BCs[bc_id]->PhiRhs(CurElem->node(i));
	            }
	          }
	        }
#endif
        } // end if CurElem->neigbor(side)==NULL
      } // endfor side
#ifndef USE_NEW_DIRICHLET
#if 0
      for (unsigned int side=0; side<CurElem->n_sides(); side++) {
	if (CurElem->neighbor(side) == NULL) {
          short int bc_id = mesh.boundary_info->boundary_id (CurElem,side);
          assert(bc_id!=BoundaryInfo::invalid_id);
          double NeumCoef, DiriCoef, Rhs;
          // std::cout << "bc_id=" << bc_id << std::endl;
          NeumCoef = BCs[bc_id].PhiNeumannCoef;
          DiriCoef = BCs[bc_id].PhiDirichletCoef;
          Rhs = BCs[bc_id].PhiRhs;
          if ( fabs(NeumCoef)<=eps) { // handle Dirichlet boundary conditions
	    assert(DiriCoef != 0.);
	    // loop over all nodes and find the ones corresponding to this side
	    for (unsigned int i=0; i<CurElem->n_nodes(); i++) {
	      if (CurElem->is_node_on_side(i, side)) {
		ElemRes(i) = DiriCoef*LocalSolution[i];
		ElemRes(i) += Rhs;
	      }
	    }
	  }	  
	}
      }
#endif
      for (unsigned int i=0;i<CurElem->n_nodes();i++) { // handle Dirichlet BC for nodes: set whole line to 0, set DiriCoef on main diagonal and Rhs in vector
        std::vector<short int> bc_id = mesh.boundary_info->boundary_ids(CurElem->get_node(i));
        if (bc_id.empty())
          continue;
        if (bc_id[0]==BoundaryInfo::invalid_id)
          continue;
	      double diricoeff = 1.;
	      double rhs = 0.;
	      if (bc_id[0]!=666666) {
	        diricoeff = BCs[bc_id[0]]->PhiDirichletCoef(CurElem->get_node(i));
	        rhs = BCs[bc_id[0]].PhiRhs(CurElem->get_node(i));
	        if ( fabs( BCs[bc_id[0]].PhiNeumannCoef(CurElem->get_node(i)) ) > eps )
	          continue;
	      }
        // std::cout << "Dirichlet: bc_id=" << bc_id[0] << std::endl;
        ElemRes(i) = diricoeff*LocalSolution[i];
        ElemRes(i) += rhs;
      }
#endif
    }
#ifdef EXHAUST_AS_CONTROL
#if 0
# ifndef PHI_IN_OBJECTIVE
    // add first node as diriclet to pin down velocity potential
    if (GetProcID()==0 && itCurEl == mesh.active_local_elements_begin()) {
      ElemRes(0) += 1e20*LocalSolution[0];
    }
# endif
#endif
#endif
    /*std::cout << "Assembled ElemMat: " << std::endl;
    std::cout << ElemMat;
    std::cout << "Element points: " << std::endl;
    std::cout << CurElem->point(0);
    std::cout << CurElem->point(1);
    std::cout << CurElem->point(2);
    std::cout.flush();*/
    dof_map.constrain_element_vector(ElemRes, dof_indices); // Add constrains for hanging nodes (refinement)
    //std::cout << "Constrained ElemMat: " << std::endl;
    //std::cout << ElemMat;
#ifdef EXHAUST_AS_CONTROL
#ifndef PHI_IN_OBJECTIVE
    for (unsigned int i=0;i<CurElem->n_nodes();i++) { // handle Dirichlet BC for nodes: set whole line to 0, set DiriCoef on main diagonal and Rhs in vector
      if(CurElem->node(i)==node_id_for_phi0_) {
        ElemRes(i) = LocalSolution[i];
      }
    }
#endif
#endif
    lm_pde_residual_vec_->add_vector(ElemRes, dof_indices);
  }
  lm_pde_residual_vec_->close();
  residual = lm_pde_residual_vec_;
  DBG_PRINT("LibMeshPDEBase::calcPDE_residual finished");
}

void LibMeshPDEBase::calcPDE_jacobians(libMesh::SparseMatrix<libMesh::Number>*& jac_state, libMesh::SparseMatrix<libMesh::Number>*& jac_control)
{
  DBG_PRINT( "LibMeshPDEBase::calcPDE_jacobians called" );
  calcPDE_jacobian_state(jac_state);
  calcPDE_jacobian_control(jac_control);
  DBG_PRINT( "LibMeshPDEBase::calcPDE_jacobians finished" );
}

void LibMeshPDEBase::calcPDE_jacobian_state(libMesh::SparseMatrix<libMesh::Number>*& jac_state)
{
  ConvertControl2PGData();
  DBG_PRINT( "LibMeshPDEBase::calcPDE_jacobian_state called" );
  ImplicitSystem& sys=lm_eqn_sys_->get_system<ImplicitSystem>("PDE");
  sys.assemble();
  sys.matrix->close();
  sys.rhs->close();
  jac_state = sys.matrix;
  DBG_PRINT( "LibMeshPDEBase::calcPDE_jacobian_state finished" );
}

void LibMeshPDEBase::calcPDE_jacobian_control(libMesh::SparseMatrix<libMesh::Number>*& jac_control)
{
  DBG_PRINT( "LibMeshPDEBase::calcPDE_jacobian_control called" );
  ConvertControl2PGData();
  jac_control_->zero();

  const double eps = 1e-8;

  const unsigned int dim = mesh_.mesh_dimension();
  LinearImplicitSystem& system = lm_eqn_sys_->get_system<LinearImplicitSystem>("PDE");
  const DofMap& dof_map = system.get_dof_map();
  FEType fe_type = dof_map.variable_type(0);
  std::vector<BoundaryConditionBase*>& BCs = PG_._BoundCond; // Mapping boudary info (index) -> boundary condition, set up in Problem geometry class

  AutoPtr<FEBase> fe (FEBase::build(dim, fe_type)); // this object will by the current volume object, actualy holding values for current phi and dphi
  QGauss qrule (dim, FIFTH);
  fe->attach_quadrature_rule (&qrule);
  AutoPtr<FEBase> fe_face (FEBase::build(dim, fe_type)); // surface object
  QGauss qface(dim-1, FIFTH);
  fe_face->attach_quadrature_rule (&qface);

  DenseMatrix<Number> ElemMat;  // Matrices for current element
  std::vector<unsigned int> dof_indices; // mapping local index -> global index
  std::vector<Number> LocalSolution;   // solution values of current element

  MeshBase::const_element_iterator itCurEl = mesh_.active_local_elements_begin();
  const MeshBase::const_element_iterator itEndEl = mesh_.active_local_elements_end();
  for ( ; itCurEl != itEndEl; ++itCurEl) {
    const Elem* CurElem = *itCurEl;
    dof_map.dof_indices(CurElem, dof_indices);  // setup local->global mapping in form of array
    fe->reinit(CurElem); // reinit fe to fit the current element, i.e. recalulate JxW and dPhi
    ElemMat.resize (dof_indices.size(),dof_indices.size());
    // Point Center = (CurElem->point(0)+CurElem->point(1)+CurElem->point(2))/3.0;
    // std::cout << "Element at " << Center << std::endl;
    for (unsigned int side=0; side<CurElem->n_sides(); side++) {
      if (CurElem->neighbor(side) == NULL) {
        fe_face->reinit(CurElem,side);        // Missing before, check
        const std::vector<Point>& q_point = fe_face->get_xyz();
        // std::cout << side << std::endl;
        short int bc_id = mesh_.boundary_info->boundary_id (CurElem,side);
        assert(bc_id!=BoundaryInfo::invalid_id);
        for (unsigned int icontrol=0;icontrol<lm_control_vec_->size();icontrol++) {
          int BoundCntrl = PG_._ParamIdx2BCParam[icontrol].BoundaryMarker;
          if (BoundCntrl!=bc_id)
            continue;
          LocalSolution.resize(dof_indices.size());
          for (unsigned int iCurSol=0;iCurSol<LocalSolution.size();iCurSol++) {
            LocalSolution[iCurSol] = system.current_local_solution->el(dof_indices[iCurSol]);
          }
          double NeumCoef, DiriCoef, Rhs;
          // std::cout << "bc_id=" << bc_id << std::endl;
          NeumCoef = BCs[bc_id]->PhiNeumannCoef(q_point[0]);
          DiriCoef = BCs[bc_id]->PhiDirichletCoef(q_point[0]);
          if(fabs(NeumCoef)>eps) { // handle non-Dirichlet boundary conditions
            const std::vector<std::vector<Real> >&  phi_face = fe_face->get_phi();
            const std::vector<Real>& JxW_face = fe_face->get_JxW();

            switch (PG_._ParamIdx2BCParam[icontrol].BCParameter) { // 0:PhiDiriCoef, 1: PhiNeumCoef, 2: PhiRhs, 3:TDiriCoef, 4: TNeumCoef, 5: TiRhs
            case 0:
              if (calc_type_ == StructureOnly) {
                for (unsigned int qp=0; qp<qface.n_points(); qp++)
                  for (unsigned int i=0; i<phi_face.size(); i++)
                    for (unsigned int j=0; j<phi_face.size(); j++) {
                      jac_control_->add(dof_indices[i], icontrol, 1.0);
                    }
              }
              else {
                for (unsigned int qp=0; qp<qface.n_points(); qp++) {
                  for (unsigned int i=0; i<phi_face.size(); i++)
                    for (unsigned int j=0; j<phi_face.size(); j++) {
                      jac_control_->add(dof_indices[i], icontrol, -(1.0/NeumCoef)*JxW_face[qp]*phi_face[i][qp]*LocalSolution[j]);
                    }
                }
              }
              break;
            case 1:
              if (calc_type_ == StructureOnly) {
                for (unsigned int qp=0; qp<qface.n_points(); qp++) {
                  for (unsigned int i=0; i<phi_face.size(); i++)
                    for (unsigned int j=0; j<phi_face.size(); j++)
                      jac_control_->add(dof_indices[i], icontrol,1.0);
                  for (unsigned int i=0; i<phi_face.size(); i++)
                    jac_control_->add(dof_indices[i], icontrol,1.0);
                }
              }
              else {
                for (unsigned int qp=0; qp<qface.n_points(); qp++) {
                  for (unsigned int i=0; i<phi_face.size(); i++)
                    for (unsigned int j=0; j<phi_face.size(); j++)
                      jac_control_->add(dof_indices[i], icontrol,(DiriCoef/(NeumCoef*NeumCoef))*JxW_face[qp]*phi_face[i][qp]*LocalSolution[j]);
                  for (unsigned int i=0; i<phi_face.size(); i++)
                    jac_control_->add(dof_indices[i], icontrol,(Rhs/(NeumCoef*NeumCoef))*JxW_face[qp]*phi_face[i][qp]);
                }
              }
              break;
            case 2:
              if (calc_type_ == StructureOnly) {
                for (unsigned int qp=0; qp<qface.n_points(); qp++) {
                  for (unsigned int i=0; i<phi_face.size(); i++)
                    jac_control_->add(dof_indices[i], icontrol,1.0);
                }
              }
              else {
                for (unsigned int qp=0; qp<qface.n_points(); qp++) {
                  for (unsigned int i=0; i<phi_face.size(); i++) {
                    BoundaryConditionSquarePhiRhs* pBC=dynamic_cast<BoundaryConditionSquarePhiRhs*>(BCs[bc_id]);
                    assert(pBC);
                    double dPhiRhs_dCntrl = BCs[bc_id]->PhiRhs(q_point[qp])/pBC->_PhiRhsScale;
                    jac_control_->add(dof_indices[i], icontrol,-(dPhiRhs_dCntrl/NeumCoef)*JxW_face[qp]*phi_face[i][qp]);
                  }
                }
              }
              break;
            }
	        } // end if(fabs(NeumCoef)>eps)
#ifdef USE_NEW_DIRICHLET
	        else {
	          switch (PG_._ParamIdx2BCParam[icontrol].BCParameter) { // 0:PhiDiriCoef, 1: PhiNeumCoef, 2: PhiRhs, 3:TDiriCoef, 4: TNeumCoef, 5: TiRhs
	          case 0:
	            if (calc_type_ == StructureOnly)
		            for (unsigned int i=0; i<CurElem->n_nodes(); i++) {
		              if (CurElem->is_node_on_side(i, side)) {
		                jac_control_->add(dof_indices[i], icontrol,1.0);
		              }
		            }
	            else
		            for (unsigned int i=0; i<CurElem->n_nodes(); i++) {
		              if (CurElem->is_node_on_side(i, side)) {
		                jac_control_->add(dof_indices[i], icontrol,system.current_local_solution->el(dof_indices[i]));
		              }
		            }
	            break;
	          case 1:
	            assert(false); // derivative would be ... / NeumCoef^2, but NeumCoef < 1e-8
	          case 2:
	            if (calc_type_ == StructureOnly)
		            for (unsigned int i=0; i<CurElem->n_nodes(); i++) {
		              if (CurElem->is_node_on_side(i, side)) {
		                jac_control_->add(dof_indices[i], icontrol,1.0);
		              }
		            }
	            else
		            for (unsigned int i=0; i<CurElem->n_nodes(); i++) {
		              if (CurElem->is_node_on_side(i, side)) {
		                jac_control_->add(dof_indices[i], icontrol,1.0);
		              }
		            }
	          }
	        } // end else(fabs(NeumCoef)>eps)
#endif
        } // end Loop over Controls
      } // end if boundary side
    } // end side loop
#ifndef USE_NEW_DIRICHLET
    for (unsigned int i=0;i<CurElem->n_nodes();i++) { // handle Dirichlet BC for nodes: set whole line to 0, set DiriCoef on main diagonal and Rhs in vector
      std::vector<short int> bc_id = mesh_.boundary_info->boundary_ids(CurElem->get_node(i));
      if (bc_id.empty())
        continue;
      if (bc_id[0]==BoundaryInfo::invalid_id)
        continue;
      if ( fabs(BCs[bc_id[0]].PhiNeumannCoef) > eps )
        continue;
      for (unsigned int iControl=0;iControl<lm_control_vec_->size();iControl++) {
        int BoundCntrl = PG_._ParamIdx2BCParam[iControl].BoundaryMarker;
        if (BoundCntrl!=bc_id[0])
          continue;
        switch (PG_._ParamIdx2BCParam[iControl].BCParameter) { // 0:PhiDiriCoef, 1: PhiNeumCoef, 2: PhiRhs, 3:TDiriCoef, 4: TNeumCoef, 5: TiRhs
        case 0:
          if (calc_type_ == StructureOnly)
            jac_control_->set(dof_indices[i], iControl,1.0);
          else
            jac_control_->set(dof_indices[i], iControl,system.current_local_solution->el(dof_indices[i]));
          break;
        case 1:
          assert(false); // derivative would be ... / NeumCoef^2, but NeumCoef < 1e-8
        case 2:
          if (calc_type_ == StructureOnly)
            jac_control_->set(dof_indices[i], iControl,1.0);
          else
            jac_control_->set(dof_indices[i], iControl,1.0);
        }
      }
    }
#endif
    // dof_map.constrain_element_matrix(ElemMat, dof_indices); //There are no hangin nodes on the boundary, thus, no constraints needed
  }
  jac_control_->close();
  jac_control = jac_control_;
  DBG_PRINT( "LibMeshPDEBase::calcPDE_jacobian_control finished" );
}

// Function to assemble the PDE for Phi, i.e.
// \delta Phi = 0, a dPhi/dn = b Phi + c
// Weak formulation: \int_O \delta Phi v dx = - \int_O \nabla Phi \nabla v + \int_bd dPhi/dn v do
// => \int_O \nabla Phi \nabla v - \int_bd b/a Phi v do = +\int_bd c/a v do
// last eqation is implemented, Dirichlet-Cond are set by line with diagonal entry and rhs
void LibMeshPDEBase::assemble_Phi_PDE(EquationSystems& es, const std::string& system_name)
{
  DBG_PRINT( "LibMeshPDEBase::assemble_Phi_PDE called" );
  const double eps = 1e-8;
  libmesh_assert (system_name == "PDE");

  const MeshBase& mesh = es.get_mesh();
  const unsigned int dim = mesh.mesh_dimension();
  LinearImplicitSystem& system = es.get_system<LinearImplicitSystem>("PDE");
  const DofMap& dof_map = system.get_dof_map();
  FEType fe_type = dof_map.variable_type(0);
  LibMeshPDEBase* pData = es.parameters.get< LibMeshPDEBase* >("LibMeshPDEBase");
  std::vector<BoundaryConditionBase*>& BCs = pData->PG_._BoundCond; // Mapping boudary info (index) -> boundary condition, set up in Problem geometry class
  CalculationModeType calc_type=pData->calc_type_;

  AutoPtr<FEBase> fe (FEBase::build(dim, fe_type)); // this object will by the current volume object, actualy holding values for current phi and dphi
  QGauss qrule (dim, FIFTH);
  fe->attach_quadrature_rule (&qrule);
  AutoPtr<FEBase> fe_face (FEBase::build(dim, fe_type)); // surface object
  QGauss qface(dim-1, FIFTH);
  fe_face->attach_quadrature_rule (&qface);
  const std::vector<Real>& JxW = fe->get_JxW();   // References to values hold by fe (i.e. the chang, ones fe changes
  const std::vector<std::vector<Real> >& phi = fe->get_phi();
  const std::vector<std::vector<RealGradient> >& dphi = fe->get_dphi();

  DenseMatrix<Number> ElemMat;  // Matrices for current element
  std::vector<unsigned int> dof_indices; // mapping local index -> global index

  MeshBase::const_element_iterator itCurEl = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator itEndEl = mesh.active_local_elements_end();
  for ( ; itCurEl != itEndEl; ++itCurEl) {
    const Elem* CurElem = *itCurEl;
    dof_map.dof_indices(CurElem, dof_indices);  // setup local->global mapping in form of array
    fe->reinit (CurElem); // reinit fe to fit the current element, i.e. recalulate JxW and dPhi
    ElemMat.resize (dof_indices.size(),dof_indices.size());
    // Point Center = (CurElem->point(0)+CurElem->point(1)+CurElem->point(2))/3.0;
    // std::cout << "Element at " << Center << std::endl;
    unsigned int qp,i,j;
    for (qp=0; qp<qrule.n_points(); qp++)
      for (i=0; i<phi.size(); i++)
        for (j=0; j<phi.size(); j++)
          if (calc_type==StructureOnly)
            ElemMat(i,j) += 1.0;
          else
            ElemMat(i,j) += JxW[qp]*(dphi[i][qp]*dphi[j][qp]);

    for (unsigned int side=0; side<CurElem->n_sides(); side++) {
      if (CurElem->neighbor(side) == NULL) {
        fe_face->reinit(CurElem, side);
        const std::vector<Point> q_point = fe_face->get_xyz();
        // std::cout << side << std::endl;
        short int bc_id = mesh.boundary_info->boundary_id (CurElem,side);
        assert(bc_id!=BoundaryInfo::invalid_id);
        double NeumCoef, DiriCoef, Rhs;
        // std::cout << "bc_id=" << bc_id << std::endl;
        NeumCoef = BCs[bc_id]->PhiNeumannCoef(q_point[0]);
        DiriCoef = BCs[bc_id]->PhiDirichletCoef(q_point[0]);
        const std::vector<std::vector<Real> >&  phi_face = fe_face->get_phi();
        const std::vector<Real>& JxW_face = fe_face->get_JxW();
        /*Point Center;
        for(unsigned int iNode=0;iNode<CurElem->n_nodes();iNode++)
         if( CurElem->is_node_on_side(iNode,side) )
          Center += CurElem->point(iNode);
        Center = Center/dim;
        std::cout << "BC at " << Center << " and : ";
        std::cout << NeumCoef << " dPhi/dn = " << DiriCoef << " Phi + " << Rhs << std::endl;*/
        if ( fabs(NeumCoef)>eps) { // handle non-Dirichlet boundary conditions
          for (unsigned int qp=0; qp<qface.n_points(); qp++) {
            for (unsigned int i=0; i<phi_face.size(); i++)
              for (unsigned int j=0; j<phi_face.size(); j++)
                if (calc_type==StructureOnly)
                  ElemMat(i,j) += 1.0;
                else
                  ElemMat(i,j) -= (DiriCoef/NeumCoef)*JxW_face[qp]*phi_face[i][qp];
          }
        }
#ifdef USE_NEW_DIRICHLET
	else {
	  assert(DiriCoef != 0.);
	  for (unsigned int i=0; i<CurElem->n_nodes(); i++) {
	    if (CurElem->is_node_on_side(i, side)) {
	      //for (j=0;j<CurElem->n_nodes();j++)
	      //ElemMat(i,j) = 0;
	      if (calc_type==StructureOnly)
      		ElemMat(i,i) += 1.0;
	      else
      		ElemMat(i,i) += DiriCoef;
	    }
	  }
	}
#endif
      } // end if CurElem->neigbor(side)==NULL
    } // endfor side
#ifndef USE_NEW_DIRICHLET
    for (unsigned int i=0;i<CurElem->n_nodes();i++) { // handle Dirichlet BC for nodes: set whole line to 0, set DiriCoef on main diagonal and Rhs in vector
      std::vector<short int> bc_id = mesh.boundary_info->boundary_ids(CurElem->get_node(i));
      if (bc_id.empty())
        continue;
      if (bc_id[0]==BoundaryInfo::invalid_id)
        continue;
      //if ( fabs(BCs[bc_id[0]].PhiNeumannCoef) > eps )
      //continue;
      // std::cout << "Dirichlet: bc_id=" << bc_id[0] << std::endl;
      double diricoeff = 1.;
      if (bc_id[0]!=666666) {
	diricoeff = +BCs[bc_id[0]].PhiDirichletCoef;
	if ( fabs(BCs[bc_id[0]].PhiNeumannCoef) > eps )
	  continue;
      }
      //printf("bc_id[0] = %d\n", bc_id[0]);
      for (j=0;j<CurElem->n_nodes();j++)
        ElemMat(i,j)=0;
      if (calc_type==StructureOnly)
        ElemMat(i,i) = 1.0;
      else
        ElemMat(i,i) = diricoeff;
    }
#endif
#ifdef EXHAUST_AS_CONTROL
#if 0
# ifndef PHI_IN_OBJECTIVE
    // add first node as diriclet to pin down velocity potential
    if (GetProcID()==0 && itCurEl == mesh.active_local_elements_begin()) {
      ElemMat(0,0) += 1.0e20;
    }
# endif
#endif
#endif
    dof_map.constrain_element_matrix(ElemMat, dof_indices); // Add constrains for hanging nodes (refinement)
#ifdef EXHAUST_AS_CONTROL
#ifndef PHI_IN_OBJECTIVE
    for (unsigned int i=0;i<CurElem->n_nodes();i++) { // handle Dirichlet BC for nodes: set whole line to 0, set DiriCoef on main diagonal and Rhs in vector
      if(CurElem->node(i)==pData->node_id_for_phi0_) {
	for (j=0;j<CurElem->n_nodes();j++)
	  ElemMat(i,j)=0;
        ElemMat(i,i) = 1.0;
      }
    }
#endif
#endif
    system.matrix->add_matrix(ElemMat, dof_indices);
  }
  DBG_PRINT( "LibMeshPDEBase::assemble_Phi_PDE finished" );
}

// sigma and lambda_pde not used, avoid warning
void LibMeshPDEBase::calc_hessians(Number sigma, libMesh::DenseVector<Number>& lambda_pde,
                                   libMesh::DenseVector<Number>& lambda_aux,
                                   libMesh::SparseMatrix<Number>*& Hcc,
                                   libMesh::SparseMatrix<Number>*& Hcs,
                                   libMesh::SparseMatrix<Number>*& Hss)
{
  DBG_PRINT( "LibMeshPDEBase::calc_hessians called" );
  hess_control_control_->zero();
  hess_state_state_->zero();
  hess_control_state_->zero();
  ConvertControl2PGData();

  // The only functions, that have an non-zero hessian are the inequality constraiants (aux_constr)
  int i_aux_constr(first_aux_constr_);
  const MeshBase& mesh = lm_eqn_sys_->get_mesh();
  const unsigned int dim = mesh.mesh_dimension();
  LinearImplicitSystem& system = lm_eqn_sys_->get_system<LinearImplicitSystem>("PDE");
  const DofMap& dof_map = system.get_dof_map();
  FEType fe_type = dof_map.variable_type(0);
  AutoPtr<FEBase> fe_face (FEBase::build(dim, fe_type));  // surface object
  QGauss qface(dim-1, FIFTH);
  fe_face->attach_quadrature_rule (&qface);
  const std::vector<std::vector<RealGradient> >&  dphi_face = fe_face->get_dphi();
  const std::vector<Real>& JxW_face = fe_face->get_JxW();
  std::vector<unsigned int> dof_indices;  // mapping local index -> global index
  MeshBase::const_element_iterator itCurEl = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator itEndEl = mesh.active_local_elements_end();

  for ( ; itCurEl != itEndEl; ++itCurEl) {
    const Elem* CurElem = *itCurEl;
    for (unsigned int side=0; side<CurElem->n_sides(); side++) {
      if (CurElem->neighbor(side) == NULL) {
        short int bc_id = mesh.boundary_info->boundary_id (CurElem,side);
        assert(bc_id!=BoundaryInfo::invalid_id);
        if ( AuxConstrBoundMarkerList_.find(bc_id)!=AuxConstrBoundMarkerList_.end() ) { // Heating boundary, apply min. velocity-constraint
          fe_face->reinit(CurElem, side);
          dof_map.dof_indices(CurElem, dof_indices);   // setup local->global mapping in form of array
          DenseMatrix<Number> loc_l2dphi;
          loc_l2dphi.resize(dof_indices.size(),dof_indices.size());
          if (calc_type_==StructureOnly) {
            for (unsigned short inode=0;inode<CurElem->n_nodes();++inode) {
              for (unsigned short jnode=0;jnode<inode;++jnode) {
                loc_l2dphi(inode,jnode) = 1.0;
              }
            }
            hess_state_state_->add_matrix(loc_l2dphi,dof_indices);
          }
          else {
            DenseMatrix<Number> loc_l2dphi;
            loc_l2dphi.resize(dof_indices.size(),dof_indices.size());
            DenseVector<Number> tmp;
#ifndef SCALE_AUX_BOUNDS
	    Number SideFact = 0.0;
	    for (unsigned int qp=0; qp<qface.n_points(); qp++) {
	      SideFact += JxW_face[qp];
	    }
#endif
            for (unsigned int qp=0; qp<qface.n_points(); qp++) {
              for (unsigned short inode=0;inode<CurElem->n_nodes();++inode) {
                for ( unsigned short jnode=0;jnode<=inode;++jnode) { // only lower triangle
                  loc_l2dphi(inode,jnode) = dphi_face[inode][qp]*dphi_face[jnode][qp];
                }
              }
#ifdef SCALE_AUX_BOUNDS
              loc_l2dphi *= 2.0*lambda_aux(i_aux_constr)*JxW_face[qp];
#else
              loc_l2dphi *= 2.0*lambda_aux(i_aux_constr)*JxW_face[qp]/SideFact;
#endif
              hess_state_state_->add_matrix(loc_l2dphi,dof_indices);
            }
          }
          i_aux_constr++;
        }
      }
    }
  }
#ifdef EXHAUST_AS_CONTROL
# ifdef PHI_IN_OBJECTIVE
  if (GetProcID()==0 && sigma!=0.0) {
    // quadratic objective to make the first node in the first element zero
    // this is to pin down the potential
    MeshBase::const_element_iterator itCurEl = mesh.active_local_elements_begin();
    const Elem* CurElem = *itCurEl;
    dof_map.dof_indices(CurElem, dof_indices);
    hess_state_state_->add(dof_indices[0], dof_indices[0], 2.*sigma);
  }
# endif
#endif

  hess_control_control_->close();
  hess_state_state_->close();
  hess_control_state_->close();
  Hcc = hess_control_control_;
  Hss = hess_state_state_;
  Hcs = hess_control_state_;
  DBG_PRINT( "LibMeshPDEBase::calc_hessians finished" );
}

void LibMeshPDEBase::get_bounds(libMesh::NumericVector<libMesh::Number>& state_l,
                                libMesh::NumericVector<libMesh::Number>& state_u,
                                libMesh::NumericVector<libMesh::Number>& control_l,
                                libMesh::NumericVector<libMesh::Number>& control_u,
                                libMesh::NumericVector<libMesh::Number>& aux_constr_l,
                                libMesh::NumericVector<libMesh::Number>& aux_constr_u)
{
  DBG_PRINT( "LibMeshPDE::get_bounds finished" );
  double Inf = 1e20;
  state_l=-Inf;
  state_u=Inf;
  if (simulation_mode_) {
    control_l=*lm_control_vec_;
    control_u=*lm_control_vec_;
    aux_constr_l = -100;
#ifdef EXHAUST_AS_CONTROL
    if (GetProcID()==0) {
      int iControl=PG_._AC.size();
      for (int iExh=0;iExh<PG_._Exh.size();++iExh) {
	control_l.set(iControl,-1e50);
	control_u.set(iControl,1e50);
	iControl++;
      }
    }
#endif    
  }
  else {
    control_l=0.0;
    control_u=Inf;
    aux_constr_l = *lm_aux_constr_vec_low_bd_; // lower bound for tangetial flow on equipment boundary
  }
  aux_constr_u = Inf;
  DBG_PRINT( "LibMeshPDE::get_bounds finished" );
}

void LibMeshPDEBase::
get_starting_point(libMesh::NumericVector<libMesh::Number>& state,
		   libMesh::NumericVector<libMesh::Number>& control,
		   bool init_z,
		   libMesh::NumericVector<libMesh::Number>* state_lb_mults,
		   libMesh::NumericVector<libMesh::Number>* state_ub_mults,
		   libMesh::NumericVector<libMesh::Number>* control_lb_mults,
		   libMesh::NumericVector<libMesh::Number>* control_ub_mults,
		   bool init_lambda,
		   libMesh::NumericVector<libMesh::Number>* pde_residual_mults,
		   libMesh::NumericVector<libMesh::Number>* aux_constr_mults)
{
  DBG_PRINT( "LibMeshPDE::get_starting_point called" );
  state = getStateVector();
  control = getControlVector();

  libMesh::LinearImplicitSystem& lm_sys = lm_eqn_sys_->get_system<LinearImplicitSystem>("PDE");

  if (init_z) {
    state_lb_mults = &lm_sys.get_vector("state_lb_mults");
    state_ub_mults = &lm_sys.get_vector("state_ub_mults");
    control_lb_mults = &*lm_control_lb_mults_;
    control_ub_mults = &*lm_control_ub_mults_;
  }

  if (init_lambda) {
    pde_residual_mults = &lm_sys.get_vector("pde_residual_mults");
    aux_constr_mults = &*lm_aux_constr_mults_;
  }

  DBG_PRINT( "LibMeshPDE::get_starting_point finished" );
}

void LibMeshPDEBase::WriteNodes(const MeshBase& mesh, std::ostream& os)
{
  os << "# Nodes:" << std::endl;
  unsigned short iDim, Dim = mesh.mesh_dimension();
  unsigned int NumNodes = mesh.n_nodes();
  os << NumNodes << " " << Dim << " 0 0" << std::endl;

  MeshBase::const_node_iterator CurNode = mesh.nodes_begin();
  MeshBase::const_node_iterator EndNode = mesh.nodes_end();
  int iNode = 1;
  for (;CurNode!=EndNode;CurNode++) {
    os << "  " << iNode++ << " ";
    for (iDim=0;iDim<Dim;iDim++)
      os << (*(*CurNode))(iDim) << " ";
    os << std::endl;
  }
}

void LibMeshPDEBase::WriteNodeFile(const MeshBase& mesh, const std::string& Filename)
{
  std::ofstream file(Filename.c_str(),std::ios::out);
  WriteNodes(mesh, file);
}

void LibMeshPDEBase::WriteElems(const MeshBase& mesh, std::ostream& os)
{
  unsigned short Dim = mesh.mesh_dimension(), iNode;
  os << "Elements:" << std::endl;
  os << mesh.n_active_elem() << " " << Dim+1 << " 0" << std::endl;
  MeshBase::const_element_iterator CurEl = mesh.active_elements_begin();
  MeshBase::const_element_iterator EndEl = mesh.active_elements_end();

  int iEl = 1;
  for (;CurEl!=EndEl;CurEl++) {
    os << "  " << iEl++ << " ";
    for (iNode=0;iNode<((*CurEl)->n_nodes());iNode++)
      os << (*CurEl)->node(iNode)+1 << " ";
    os << std::endl;
  }
}

void LibMeshPDEBase::WriteEleFile(const MeshBase& mesh, const std::string& Filename)
{
  std::ofstream file(Filename.c_str(),std::ios::out);
  WriteElems(mesh, file);
}

void LibMeshPDEBase::Write2File( const std::string& pre_filename)
{
  DBG_PRINT( "LibMeshPDEBase::Write2File called" );
  std::string my_pre_filename = pre_filename;
  if (simulation_mode_)
    my_pre_filename += "Sim";

  std::string filename = my_pre_filename + "Cntrl.dat";
  std::ofstream f;
  f.open(filename.c_str(),std::ios::out);
  lm_control_vec_->print(f);
  f.close();

  std::vector<libMesh::Number> State;
  lm_eqn_sys_->build_solution_vector(State);

  filename = my_pre_filename + "State.dat";
  f.open(filename.c_str(),std::ios::out);
  for (unsigned int iVal=0;iVal<State.size();iVal++)
    f << State[iVal] << std::endl;
  f.close();

  filename = my_pre_filename + "State";
  //VTKIO vtkio(mesh_);
  //vtkio.write_equation_systems(filename,*lm_eqn_sys_);

  /* Creates Memory leak
  ExodusII_IO exoio(mesh_);
  exoio.write_equation_systems(my_pre_filename + "State.ex2", *lm_eqn_sys_);
  */

#if 0
  filename = my_pre_filename + "StatePot.csv";
  WritePotentialCSV(filename);

  WriteAirflowCSVs(my_pre_filename + "StateVolFlow.csv", my_pre_filename + "StateSurfFlow.csv");
  DBG_PRINT( "LibMeshPDEBase::Write2File finished" );

  WriteAirflowTKVs("VolFlow", "SurfFlow");
#endif
}

void LibMeshPDEBase::WritePotentialCSV(const std::string& Filename)
{
  DBG_PRINT( "LibMeshPDEBase::WritePotentialCSV called" );
  std::ofstream f;
  f.open(Filename.c_str(),std::ios::out);

  std::vector<libMesh::Number> State;
  lm_eqn_sys_->build_solution_vector(State);

  const unsigned int dim = mesh_.mesh_dimension();
  MeshBase::node_iterator it_cur_node = mesh_.active_nodes_begin();
  const MeshBase::node_iterator end_node = mesh_.active_nodes_end();
  for (;it_cur_node!=end_node;it_cur_node++) {
    const Node* cur_node = *it_cur_node;
    for (int i=0;i<dim;i++)
      f << (*cur_node)(i) << ", ";
    f << State[cur_node->id()] << std::endl;
  }
  f.close();
  DBG_PRINT( "LibMeshPDEBase::WritePotentialCSV finished" );
}

void LibMeshPDEBase::WriteAirflowTKVs(const std::string& VolumeFilename, const std::string& SurfFilename)
{
  DBG_PRINT( "LibMeshPDEBase::WriteAirflowCSVs called" );

  const unsigned int dim = mesh_.mesh_dimension();
  LinearImplicitSystem& system = lm_eqn_sys_->get_system<LinearImplicitSystem>("PDE");
  AutoPtr< NumericVector<Number> > pvx = system.solution->zero_clone(); // x component of velocity
  AutoPtr< NumericVector<Number> > pvy = system.solution->zero_clone(); // y component of velocity
  AutoPtr< NumericVector<Number> > pvz = system.solution->zero_clone(); // z component of velocity
  AutoPtr< NumericVector<Number> > pvnorm = system.solution->zero_clone(); // norm of velocity
  AutoPtr< NumericVector<Number> > pvn = system.solution->zero_clone(); // number of neighboring elements (to know by how much to divide when computing average)
  const DofMap& dof_map = system.get_dof_map();
  FEType fe_type = dof_map.variable_type(0);
  AutoPtr<FEBase> fe_vol (FEBase::build(dim, fe_type));  // volume object
  QGauss qvol(dim, FIFTH);
  fe_vol->attach_quadrature_rule (&qvol);
  const std::vector<std::vector<RealGradient> >&  dphi_vol = fe_vol->get_dphi();
  AutoPtr<FEBase> fe_face (FEBase::build(dim, fe_type));  // surface object
  QGauss qface(dim-1, FIFTH);
  fe_face->attach_quadrature_rule (&qface);
  const std::vector<std::vector<RealGradient> >&  dphi_face = fe_face->get_dphi(); // gradeint of basis function phi at quadarture points
  std::vector<unsigned int> dof_indices;  // mapping local index -> global index
  std::vector<Number> tmp;
  std::vector<Number> One;
  MeshBase::const_element_iterator itCurEl = mesh_.active_elements_begin();
  const MeshBase::const_element_iterator itEndEl = mesh_.active_elements_end();
  RealGradient CurGrad;
  for ( ; itCurEl != itEndEl; ++itCurEl) {
    // volume air flow
    const Elem* CurElem = *itCurEl;
    fe_vol->reinit(CurElem); // here: Phi and DPhi (basis function) are computed for this element
    dof_map.dof_indices(CurElem, dof_indices);   // setup local->global mapping in form of array
    CurGrad = 0.0;
    for (unsigned short inode=0;inode<CurElem->n_nodes();++inode) {
      Number node_val = system.current_local_solution->el(dof_indices[inode]);
      printf("nodeval = %e\n", node_val);
      CurGrad += node_val*dphi_vol[inode][0];  // assume linear FE -> grad = const on one element
    }
    tmp.resize(dof_indices.size(),CurGrad(0));
    pvx->add_vector(tmp,dof_indices);
    tmp.resize(dof_indices.size(),CurGrad(1));
    pvy->add_vector(tmp,dof_indices);
    tmp.resize(dof_indices.size(),CurGrad(2));
    pvz->add_vector(tmp,dof_indices);
    One.resize(dof_indices.size(),1.0);
    pvn->add_vector(One,dof_indices);
  }

  pvx->close();
  pvy->close();
  pvz->close();
  pvn->close();
  pvnorm->close();
  MeshBase::node_iterator it_cur_node = mesh_.active_nodes_begin();
  const MeshBase::node_iterator end_node = mesh_.active_nodes_end();
  for (;it_cur_node!=end_node;it_cur_node++) {
    const Node* cur_node = *it_cur_node;
    unsigned int id = cur_node->id();
    double ncount = pvn->el(id);
    double vx = pvx->el(id)/ncount;
    double vy = pvy->el(id)/ncount;
    double vz = pvz->el(id)/ncount;
    printf("%d pvx = %e pvy = %e pvz = %e\n", id, vx, vy ,vz);
    pvx->set( vx, id );
    pvy->set( vy, id );
    pvz->set( vz, id );
    printf("pvnorm[%d] = %e\n", id, std::sqrt(vx*vx + vy*vy + vz*vz));
    pvnorm->set( sqrt(vx*vx + vy*vy + vz*vz), id );
  }
  pvx->close();
  pvy->close();
  pvz->close();
  pvn->close();
  AutoPtr< NumericVector<Number> > tmpSol = system.solution->clone();
  system.solution = pvx;
  update();

  std::string VolumeFilenameX = VolumeFilename + "Velox";
  VTKIO vtkio(mesh_);
  vtkio.write_equation_systems(VolumeFilenameX,*lm_eqn_sys_);

  system.solution = pvnorm;
  update();
  VTKIO vtkio2(mesh_);
  vtkio2.write_equation_systems(VolumeFilename + "VeloNorm",*lm_eqn_sys_);

  system.solution = tmpSol;
  update();
  VTKIO vtkio3(mesh_);
  vtkio3.write_equation_systems(VolumeFilename + "Phi",*lm_eqn_sys_);
}

void LibMeshPDEBase::WriteAirflowCSVs(const std::string& VolumeFilename, const std::string& SurfFilename)
{
  DBG_PRINT( "LibMeshPDEBase::WriteAirflowCSVs called" );
  std::ofstream VolFile;
  VolFile.open(VolumeFilename.c_str(),std::ios::out);
  std::ofstream SurFile;
  SurFile.open(SurfFilename.c_str(),std::ios::out);

  const unsigned int dim = mesh_.mesh_dimension();
  LinearImplicitSystem& system = lm_eqn_sys_->get_system<LinearImplicitSystem>("PDE");
  const DofMap& dof_map = system.get_dof_map();
  FEType fe_type = dof_map.variable_type(0);
  AutoPtr<FEBase> fe_vol (FEBase::build(dim, fe_type));  // volume object
  QGauss qvol(dim, FIFTH);
  fe_vol->attach_quadrature_rule (&qvol);
  const std::vector<std::vector<RealGradient> >&  dphi_vol = fe_vol->get_dphi();
  AutoPtr<FEBase> fe_face (FEBase::build(dim, fe_type));  // surface object
  QGauss qface(dim-1, FIFTH);
  fe_face->attach_quadrature_rule (&qface);
  const std::vector<std::vector<RealGradient> >&  dphi_face = fe_face->get_dphi();
  std::vector<unsigned int> dof_indices;  // mapping local index -> global index
  MeshBase::const_element_iterator itCurEl = mesh_.active_elements_begin();
  const MeshBase::const_element_iterator itEndEl = mesh_.active_elements_end();
  RealGradient CurGrad;
  for ( ; itCurEl != itEndEl; ++itCurEl) {
    // volume air flow
    const Elem* CurElem = *itCurEl;
    fe_vol->reinit(CurElem);
    dof_map.dof_indices(CurElem, dof_indices);   // local->global mapping in form of array
    CurGrad = 0.0;
    for (unsigned short inode=0;inode<CurElem->n_nodes();++inode) {
      Number node_val = system.current_local_solution->el(dof_indices[inode]);
      CurGrad += node_val*dphi_vol[inode][0];  // assume linear FE -> grad = const on one element
    }
    Point Center = CurElem->centroid();
    for (int i=0;i<dim;i++)
      VolFile << Center(i) << ", ";
    for (int i=0;i<dim-1;i++)
      VolFile << CurGrad(i) << ", ";
    VolFile << CurGrad(dim-1) << std::endl;

    // air flow at boundary
    for (unsigned int side=0; side<CurElem->n_sides(); side++) {
      if (CurElem->neighbor(side) == NULL) {
        fe_face->reinit(CurElem, side);
        DenseVector<Number> loc_sol;
        DenseMatrix<Number> loc_l2dphi;
        loc_sol.resize(dof_indices.size());
        loc_l2dphi.resize(dof_indices.size(),dof_indices.size());
        double GradL2=0.0;
        for (unsigned short inode=0;inode<CurElem->n_nodes();++inode) {
          loc_sol(inode) = system.current_local_solution->el(dof_indices[inode]);
          for ( unsigned short jnode=0;jnode<CurElem->n_nodes();++jnode) {
            loc_l2dphi(inode,jnode) = dphi_face[inode][0]*dphi_face[jnode][0];
          }
        }
        DenseVector<Number> tmp;
        loc_l2dphi.vector_mult(tmp,loc_sol);
        GradL2 = loc_sol.dot(tmp);
        AutoPtr<Elem> Side = CurElem->build_side(side);
        Point Center = Side->centroid();
        for (int i=0;i<dim;i++)
          SurFile << Center(i) << ", ";
        SurFile << GradL2 << std::endl;
      }
    }
  }
  DBG_PRINT( "LibMeshPDEBase::WriteAirflowCSVs finished" );
}

void LibMeshPDEBase::InitAuxConstr(int *plocal, int *pglobal, std::list<Number>* pFactList)
{
  DBG_PRINT( "LibMeshPDEBase::InitAuxConstr called" );
  const double eps = 1e-8;
  std::vector<ProblemGeometry::Item> _Exh;

  int nEquipWalls = 4;
  if (mesh_.mesh_dimension()==3)
    nEquipWalls = 5;
  int BoundaryMarker;
  AuxConstrBoundMarkerList_.clear();
  pFactList->clear();
  std::set<int> EquipHeatExchangeBdList;
  for (unsigned int iEquip=0; iEquip<PG_._Equip.size();++iEquip)  {
    PG_.GetHeatExchangeBoundaryMarkers(iEquip, &EquipHeatExchangeBdList);
    AuxConstrBoundMarkerList_.insert(EquipHeatExchangeBdList.begin(),EquipHeatExchangeBdList.end());
  }

  LinearImplicitSystem& system = lm_eqn_sys_->get_system<LinearImplicitSystem>("PDE");
  const DofMap& dof_map = system.get_dof_map();
  FEType fe_type = dof_map.variable_type(0);
  const unsigned int dim = mesh_.mesh_dimension();
  AutoPtr<FEBase> fe_face (FEBase::build(dim, fe_type));  // surface object
  QGauss qface(dim-1, FIFTH);
  fe_face->attach_quadrature_rule (&qface);
  const std::vector<Real>& JxW_face = fe_face->get_JxW();
  std::vector<unsigned int> dof_indices;  // mapping local index -> global index

  first_aux_constr_=0;
  assert(plocal);
  (*plocal) = 0;
  MeshBase::const_element_iterator itCurEl = mesh_.active_local_elements_begin();
  const MeshBase::const_element_iterator itEndEl = mesh_.active_local_elements_end();
  for ( ; itCurEl != itEndEl; ++itCurEl) {
    const Elem* CurElem = *itCurEl;
    for (unsigned int side=0; side<CurElem->n_sides(); side++) {
      if (CurElem->neighbor(side) == NULL) {
        int bc_id = mesh_.boundary_info->boundary_id (CurElem,side);
        if ( AuxConstrBoundMarkerList_.find(bc_id)!=AuxConstrBoundMarkerList_.end() ) {
          ++(*plocal);

          fe_face->reinit(CurElem, side);
          dof_map.dof_indices(CurElem, dof_indices);   // setup local->global mapping in form of array
          DenseMatrix<Number> loc_l2dphi;
          loc_l2dphi.resize(dof_indices.size(),dof_indices.size());
          Number SideFact = 0.0;
          for (unsigned int qp=0; qp<qface.n_points(); qp++) {
            // JxW is Jacobian determinant times quadrature weights -> sum is quadrature area
            SideFact += JxW_face[qp];
          }
          pFactList->push_back(SideFact);
        }
      }
    }
  }

  assert(pglobal);
  MPI_Allreduce ( plocal, pglobal, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD );
  MPI_Scan(plocal,&first_aux_constr_,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
  first_aux_constr_ -= (*plocal);
  DBG_PRINT( "LibMeshPDEBase::InitAuxConstr finished" );
}

void LibMeshPDEBase::calcAux_constr(libMesh::NumericVector<libMesh::Number>*& constr)
{
  DBG_PRINT( "LibMeshPDEBase::calcAux_constr called" );
  lm_aux_constr_vec_->zero();
  int i_aux_constr(first_aux_constr_);
  const unsigned int dim = mesh_.mesh_dimension();
  LinearImplicitSystem& system = lm_eqn_sys_->get_system<LinearImplicitSystem>("PDE");
  const DofMap& dof_map = system.get_dof_map();
  FEType fe_type = dof_map.variable_type(0);
  AutoPtr<FEBase> fe_face (FEBase::build(dim, fe_type));  // surface object
  QGauss qface(dim-1, FIFTH);
  fe_face->attach_quadrature_rule (&qface);
  const std::vector<std::vector<RealGradient> >&  dphi_face = fe_face->get_dphi();
  const std::vector<Real>& JxW_face = fe_face->get_JxW();
  std::vector<unsigned int> dof_indices;  // mapping local index -> global index
  MeshBase::const_element_iterator itCurEl = mesh_.active_local_elements_begin();
  const MeshBase::const_element_iterator itEndEl = mesh_.active_local_elements_end();
  RealGradient CurGrad;
  for ( ; itCurEl != itEndEl; ++itCurEl) {
    const Elem* CurElem = *itCurEl;
    for (unsigned int side=0; side<CurElem->n_sides(); side++) {
      if (CurElem->neighbor(side) == NULL) {
        short int bc_id = mesh_.boundary_info->boundary_id (CurElem,side);
        assert(bc_id!=BoundaryInfo::invalid_id);
        if ( AuxConstrBoundMarkerList_.find(bc_id)!=AuxConstrBoundMarkerList_.end() ) { // Heating boundary, apply min. velocity-constraint
          fe_face->reinit(CurElem, side);
          dof_map.dof_indices(CurElem, dof_indices);   // setup local->global mapping in form of array
          DenseVector<Number> loc_sol;
          DenseMatrix<Number> loc_l2dphi;
          loc_sol.resize(dof_indices.size());
          loc_l2dphi.resize(dof_indices.size(),dof_indices.size());
          double GradL2=0.0;
	  Number SideFact = 0.0;
	  for (unsigned int qp=0; qp<qface.n_points(); qp++) {
            for (unsigned short inode=0;inode<CurElem->n_nodes();++inode) {
              loc_sol(inode) = system.current_local_solution->el(dof_indices[inode]);
              for ( unsigned short jnode=0;jnode<CurElem->n_nodes();++jnode) {
                loc_l2dphi(inode,jnode) = dphi_face[inode][qp]*dphi_face[jnode][qp];
              }
            }
            DenseVector<Number> tmp;
            loc_l2dphi.vector_mult(tmp,loc_sol);
            GradL2 += JxW_face[qp]*loc_sol.dot(tmp);
	    SideFact += JxW_face[qp];
          }
#ifdef SCALE_AUX_BOUNDS
          lm_aux_constr_vec_->set(i_aux_constr,GradL2);
#else
          lm_aux_constr_vec_->set(i_aux_constr,GradL2/SideFact);
#endif
          ++i_aux_constr;
        }
      }
    }
  }
  lm_aux_constr_vec_->close();
  constr = lm_aux_constr_vec_;
  DBG_PRINT( "LibMeshPDEBase::calcAux_constr finished" );
}

void LibMeshPDEBase::calcAux_jacobians(libMesh::SparseMatrix<libMesh::Number>*& jac_state,
                                       libMesh::SparseMatrix<libMesh::Number>*& jac_control)
{
  DBG_PRINT( "LibMeshPDEBase::calcAux_jacobians called" );
  calcAux_jacobian_state(jac_state);
  calcAux_jacobian_control(jac_control);
  DBG_PRINT( "LibMeshPDEBase::calcAux_jacobians finished" );
}

void LibMeshPDEBase::calcAux_jacobian_state(libMesh::SparseMatrix<libMesh::Number>*& jac_state)
{
  DBG_PRINT( "LibMeshPDEBase::calcAux_jacobian_state called" );
  jac_aux_state_->zero();

  int i_aux_constr(first_aux_constr_);
  const MeshBase& mesh = lm_eqn_sys_->get_mesh();
  const unsigned int dim = mesh.mesh_dimension();
  LinearImplicitSystem& system = lm_eqn_sys_->get_system<LinearImplicitSystem>("PDE");
  const DofMap& dof_map = system.get_dof_map();
  FEType fe_type = dof_map.variable_type(0);
  AutoPtr<FEBase> fe_face (FEBase::build(dim, fe_type));  // surface object
  QGauss qface(dim-1, FIFTH);
  fe_face->attach_quadrature_rule (&qface);
  const std::vector<std::vector<RealGradient> >&  dphi_face = fe_face->get_dphi();
  const std::vector<Real>& JxW_face = fe_face->get_JxW();
  std::vector<unsigned int> dof_indices;  // mapping local index -> global index
  MeshBase::const_element_iterator itCurEl = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator itEndEl = mesh.active_local_elements_end();
  for ( ; itCurEl != itEndEl; ++itCurEl) {
    const Elem* CurElem = *itCurEl;
    for (unsigned int side=0; side<CurElem->n_sides(); side++) {
      if (CurElem->neighbor(side) == NULL) {
        short int bc_id = mesh.boundary_info->boundary_id (CurElem,side);
        assert(bc_id!=BoundaryInfo::invalid_id);
        if ( AuxConstrBoundMarkerList_.find(bc_id)!=AuxConstrBoundMarkerList_.end() ) { // Heating boundary, apply min. velocity-constraint
          fe_face->reinit(CurElem, side);
          dof_map.dof_indices(CurElem, dof_indices);   // setup local->global mapping in form of array
          if (calc_type_==StructureOnly)
            for (unsigned short inode=0;inode<CurElem->n_nodes();++inode) {
              jac_aux_state_->add(i_aux_constr,dof_indices[inode],1.0);
            }
          else {
            DenseVector<Number> loc_sol;
            DenseMatrix<Number> loc_l2dphi;
            loc_sol.resize(dof_indices.size());
            loc_l2dphi.resize(dof_indices.size(),dof_indices.size());
            DenseVector<Number> tmp;
#ifndef SCALE_AUX_BOUNDS
	    Number SideFact = 0.0;
	    for (unsigned int qp=0; qp<qface.n_points(); qp++) {
	      SideFact += JxW_face[qp];
	    }
#endif
            for (unsigned int qp=0; qp<qface.n_points(); qp++) {
              for (unsigned short inode=0;inode<CurElem->n_nodes();++inode) {
                loc_sol(inode) = system.current_local_solution->el(dof_indices[inode]);
                for ( unsigned short jnode=0;jnode<CurElem->n_nodes();++jnode) {
                  loc_l2dphi(inode,jnode) = dphi_face[inode][qp]*dphi_face[jnode][qp];
                }
              }
              loc_l2dphi.vector_mult(tmp,loc_sol);
              for (unsigned short inode=0;inode<CurElem->n_nodes();++inode) {
#ifdef SCALE_AUX_BOUNDS
                jac_aux_state_->add(i_aux_constr,dof_indices[inode],2.0*JxW_face[qp]*tmp(inode));
#else
                jac_aux_state_->add(i_aux_constr,dof_indices[inode],2.0*JxW_face[qp]*tmp(inode)/SideFact);
#endif
              }
            }
          }
          i_aux_constr++;
        }
      }
    }
  }
  jac_aux_state_->close();
  jac_state = jac_aux_state_;
  DBG_PRINT( "LibMeshPDEBase::calcAux_jacobian_state finished" );
}

void LibMeshPDEBase::calcAux_jacobian_control(libMesh::SparseMatrix<libMesh::Number>*& jac_control)
{
  DBG_PRINT( "LibMeshPDEBase::calcAux_jacobian_control called" );
  jac_aux_control_->zero();
  jac_aux_control_->close();
  jac_control = jac_aux_control_;
  DBG_PRINT( "LibMeshPDEBase::calcAux_jacobian_control finished" );
}

void LibMeshPDEBase::RefineMesh(int iter)
{
  // here we print the mesh boundary before refinement
  if (0) {
    MeshBase::const_element_iterator itCurEl = mesh_.active_local_elements_begin();
    const MeshBase::const_element_iterator itEndEl = mesh_.active_local_elements_end();
    for ( ; itCurEl != itEndEl; ++itCurEl) {
      const Elem* CurElem = *itCurEl;
      const unsigned int n_nodes = CurElem->n_nodes();
      for (unsigned int side=0; side<CurElem->n_sides(); side++) {
	if (CurElem->neighbor(side) == NULL) {
	  short int bc_id = mesh_.boundary_info->boundary_id (CurElem,side);
	  printf("side: %d bc_id: %2d ", side, bc_id);
	  for (unsigned int node=0; node<n_nodes; node++) {
	    if (CurElem->is_node_on_side(node, side)) {
	      const Point point = CurElem->point(node);
	      printf("x=(%d)%e y(%d)=%e ", node, point(0), node, point(1));
	    }
	  }
	  printf("\n");
	}
      }
    }
    itCurEl = mesh_.active_local_elements_begin();
    for ( ; itCurEl != itEndEl; ++itCurEl) {
      const Elem* CurElem = *itCurEl;
      for (unsigned int i=0;i<CurElem->n_nodes();i++) { // handle Dirichlet BC for nodes: set whole line to 0, set DiriCoef on main diagonal and Rhs in vector
        std::vector<short int> bc_id = mesh_.boundary_info->boundary_ids(CurElem->get_node(i));
	if (bc_id.size()>0) {
	  printf("node: %d ", i);
	  for (int j=0; j<bc_id.size(); j++) {
	    DBG_PRINT("bc_id[" << j << "] = " << bc_id[j] );
	  }
	}
      }
    }    
    fflush(stdout);
  }

  DBG_PRINT( "LibMeshPDEBase::RefineMesh called" );
  //MeshRefinement rf(mesh_);
  //rf.uniformly_refine(1);
  MeshRefinement mesh_refinement(mesh_);
  mesh_refinement.refine_fraction() = 1.0;
  mesh_refinement.coarsen_fraction() = 0.0;
  mesh_refinement.max_h_level() = 100;
  mesh_refinement.absolute_global_tolerance() = 1e-1;
  
  ErrorVector error;
  KellyErrorEstimator error_estimator;
  //LaplacianErrorEstimator error_estimator;
  //PatchRecoveryErrorEstimator error_estimator;
  ImplicitSystem& sys = lm_eqn_sys_->get_system<ImplicitSystem>("PDE");
  error_estimator.estimate_error(sys, error);

  if (0) {
    // plot the error
    char buffer[255];
    sprintf(buffer, "error_mesh-%d.ex2", iter);
    std::string fname(buffer);
    error.plot_error(fname, mesh_);
  }

#if 0
  //expand aux_constr_mults onto the full mesh
  libMesh::NumericVector<libMesh::Number>& full_aux_constr_mults =
    sys.get_vector("full_aux_constr_mults");
  full_aux_constr_mults.zero();
  
  int i_aux_constr = first_aux_constr_;
  MeshBase::const_element_iterator itCurEl = mesh_.active_local_elements_begin();
  const MeshBase::const_element_iterator itEndEl = mesh_.active_local_elements_end();
  DenseVector<Number> ElemRes;
  for ( ; itCurEl != itEndEl; ++itCurEl) {
    const Elem* CurElem = *itCurEl;
    dof_map.dof_indices(CurElem, dof_indices);   // setup local->global mapping in form of array
    ElemRes.resize(dof_indices.size());
    for (unsigned int side=0; side<CurElem->n_sides(); side++) {
      if (CurElem->neighbor(side) == NULL) {
        short int bc_id = mesh_.boundary_info->boundary_id (CurElem,side);
        assert(bc_id!=BoundaryInfo::invalid_id);
        if ( AuxConstrBoundMarkerList_.find(bc_id)!=AuxConstrBoundMarkerList_.end() ) {
	  

	  DenseVector<Number> loc_sol;
          DenseMatrix<Number> loc_l2dphi;
          loc_sol.resize(dof_indices.size());
          loc_l2dphi.resize(dof_indices.size(),dof_indices.size());
          double GradL2=0.0;
	  Number SideFact = 0.0;
	  for (unsigned int qp=0; qp<qface.n_points(); qp++) {
            for (unsigned short inode=0;inode<CurElem->n_nodes();++inode) {
              loc_sol(inode) = system.current_local_solution->el(dof_indices[inode]);
              for ( unsigned short jnode=0;jnode<CurElem->n_nodes();++jnode) {
                loc_l2dphi(inode,jnode) = dphi_face[inode][qp]*dphi_face[jnode][qp];
              }
            }
            DenseVector<Number> tmp;
            loc_l2dphi.vector_mult(tmp,loc_sol);
            GradL2 += JxW_face[qp]*loc_sol.dot(tmp);
	    SideFact += JxW_face[qp];
          }
#ifdef SCALE_AUX_BOUNDS
          lm_aux_constr_vec_->set(i_aux_constr,GradL2);
#else
          lm_aux_constr_vec_->set(i_aux_constr,GradL2/SideFact);
#endif
          ++i_aux_constr;
        }
      }
    }
  }
#endif


  //mesh_refinement.flag_elements_by_error_fraction(error);
  mesh_refinement.flag_elements_by_error_tolerance(error);
  mesh_refinement.refine_and_coarsen_elements();
  lm_eqn_sys_->reinit();

  std::cout << "RefineMesh" << std::endl;
  reinit();
  DBG_PRINT( "LibMeshPDEBase::RefineMesh finished" );
}

void LibMeshPDEBase::
get_finalize_vectors(libMesh::NumericVector<libMesh::Number>*& lm_state_lb_mults,
		     libMesh::NumericVector<libMesh::Number>*& lm_state_ub_mults,
		     libMesh::NumericVector<libMesh::Number>*& lm_control_lb_mults,
		     libMesh::NumericVector<libMesh::Number>*& lm_control_ub_mults,
		     libMesh::NumericVector<libMesh::Number>*& lm_pde_residual_mults,
		     libMesh::NumericVector<libMesh::Number>*& lm_aux_constr_mults)
{
  lm_control_lb_mults_ = getControlVector().clone();
  lm_control_ub_mults_ = getControlVector().clone();

  lm_aux_constr_mults_ = lm_aux_constr_vec_->clone();

  libMesh::LinearImplicitSystem& lm_sys = lm_eqn_sys_->get_system<LinearImplicitSystem>("PDE");
  lm_state_lb_mults = &lm_sys.get_vector("state_lb_mults");
  lm_state_ub_mults = &lm_sys.get_vector("state_ub_mults");
  lm_pde_residual_mults = &lm_sys.get_vector("pde_residual_mults");

  lm_control_lb_mults = &*lm_control_lb_mults_;
  lm_control_ub_mults = &*lm_control_ub_mults_;
  lm_aux_constr_mults = &*lm_aux_constr_mults_;
}

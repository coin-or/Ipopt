// Copyright (C) 2010 International Business Machines.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id$
//
// Authors:  Johannes Huber, Andreas Waechter     IBM        2010-09-03
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

#include <fstream>

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
  min_airflow(1.0),
  first_aux_constr_(0)
{
}

// clear all matrices and vectors, but not problem geomatry data
// used e.g. after refinement (since problem size changes)
void LibMeshPDEBase::clear_math_obj()
{
  AuxConstrBoundMarkerList_.clear();
  DetroySelfOwnedLibMeshPetscVector(lm_control_vec_);
  DetroySelfOwnedLibMeshPetscVector(lm_pde_residual_vec_);
  DetroySelfOwnedLibMeshPetscVector(lm_aux_constr_vec_);
  DetroySelfOwnedLibMeshPetscVector(lm_aux_constr_vec_low_bd_);

  DetroySelfOwnedLibMeshPetscMatrix(jac_control_);
  DetroySelfOwnedLibMeshPetscMatrix(jac_aux_state_);
  DetroySelfOwnedLibMeshPetscMatrix(jac_aux_control_);
  DetroySelfOwnedLibMeshPetscMatrix(hess_control_control_);
  DetroySelfOwnedLibMeshPetscMatrix(hess_control_state_);
  DetroySelfOwnedLibMeshPetscMatrix(hess_state_state_);

  if(lm_eqn_sys_)
    delete lm_eqn_sys_;
  lm_eqn_sys_ = NULL;
}

LibMeshPDEBase::~LibMeshPDEBase()
{
  clear_math_obj();
}

// is not used yet
void LibMeshPDEBase::InitProblemData(std::istream& is)
{
  PG_.ReadFromStream(is);
#if 0  // Todo: read Data from file
  double RoomXSize, RoomYSize, RoomZSize;
  int example=5;
  PG_._h=0.1;
  std::vector<double> p1,p2;
  std::cout << "Running example " << example << " with h=" << PG_._h << std::endl; 
  p1.resize(2);
  p2.resize(2);
  switch(example)
  {
    case 0:
      p1[0]=RoomXSize=1.0; p1[1]=RoomYSize=1.0;
      PG_._RoomSize = p1;
      p1[0]=0.0*RoomXSize; p1[1]=0.3*RoomYSize;
      p2[0]=0.0*RoomXSize; p2[1]=0.7*RoomYSize;
      PG_.AddAC(p1,p2,3,10);
      p1[0]=0.0*RoomXSize; p1[1]=0.7*RoomYSize;
      p2[0]=0.0*RoomXSize; p2[1]=0.9*RoomYSize;
      //PG_.AddAC(p1,p2,1,10);
      p1[0]=1.0*RoomXSize; p1[1]=0.4*RoomYSize;
      p2[0]=1.0*RoomXSize; p2[1]=0.6*RoomYSize;
      PG_.AddExhaust(p1,p2);
      p1[0]=1.0*RoomXSize; p1[1]=0.7*RoomYSize;
      p2[0]=1.0*RoomXSize; p2[1]=0.9*RoomYSize;
      //PG_.AddExhaust(p1,p2);
      p1[0]=0.4*RoomXSize; p1[1]=0.4*RoomYSize;
      p2[0]=0.6*RoomXSize; p2[1]=0.5*RoomYSize;
      PG_.AddEquipment(p1,p2,35,1);
      p1[0]=0.4*RoomXSize; p1[1]=0.2*RoomYSize;
      p2[0]=0.8*RoomXSize; p2[1]=0.4*RoomYSize;
      //PG_.AddEquipment(p1,p2,35,1);
      break;
    case 1:   // Place a lot of AC on right wall -> the ones in the cornes are most important -> ignore the middle ones
      p1[0]=RoomXSize=2.0; p1[1]=RoomYSize=1.0;
      PG_._RoomSize = p1;
      p1[0]=0.0*RoomXSize; p2[0]=0.0*RoomXSize;
      for(int i=0;i<10;i++) {
        p1[1]=(0.025+0.1*i)*RoomYSize; p2[1]=p1[1] + 0.05*RoomYSize;
        PG_.AddAC(p1,p2,3,10);
      }
      p1[0]=1.0*RoomXSize; p1[1]=0.4*RoomYSize; p2[0]=1.0*RoomXSize; p2[1]=0.6*RoomYSize;
      PG_.AddExhaust(p1,p2);
      for(int i=0;i<2;i++) {
        p1[0]=(0.2+0.4*i)*RoomXSize; p2[0]=p1[0] + 0.2*RoomXSize;
        for(int j=0;j<2;j++) {
          p1[1]=(0.2+0.4*j)*RoomYSize; p2[1]=p1[1] + 0.2*RoomYSize;
          PG_.AddEquipment(p1,p2,35,1);
        }
      }
      break;
    case 2: // Test only corner clostest two AC's -> Result doesn't change to much
      p1[0]=RoomXSize=1.0; p1[1]=RoomYSize=0.2;
      PG_._RoomSize = p1;
      p1[0]=0.0*RoomXSize; p2[0]=0.0*RoomXSize; p1[1]=0.025*RoomYSize; p2[1]=p1[1] + 0.05*RoomYSize;
      PG_.AddAC(p1,p2,3,10);
      p1[0]=0.0*RoomXSize; p2[0]=0.0*RoomXSize; p1[1]=0.925*RoomYSize; p2[1]=p1[1] + 0.05*RoomYSize;
      PG_.AddAC(p1,p2,3,10);
      p1[0]=1.0*RoomXSize; p1[1]=0.4*RoomYSize; p2[0]=1.0*RoomXSize; p2[1]=0.6*RoomYSize;
      PG_.AddExhaust(p1,p2);
      for(int i=0;i<2;i++) {
        p1[0]=(0.2+0.4*i)*RoomXSize; p2[0]=p1[0] + 0.2*RoomXSize;
        for(int j=0;j<2;j++) {
          p1[1]=(0.2+0.4*j)*RoomYSize; p2[1]=p1[1] + 0.2*RoomYSize;
          PG_.AddEquipment(p1,p2,35,1);
        }
      }
      break;
    case 3: // What, if server room isn't fully occupied? -> Turn only the AC on, that is closest to the free place 
      p1[0]=RoomXSize=1.0; p1[1]=RoomYSize=0.2;
      PG_._RoomSize = p1;
      p1[0]=0.0*RoomXSize; p2[0]=0.0*RoomXSize; p1[1]=0.025*RoomYSize; p2[1]=p1[1] + 0.05*RoomYSize;
      PG_.AddAC(p1,p2,3,10);
      p1[0]=0.0*RoomXSize; p2[0]=0.0*RoomXSize; p1[1]=0.925*RoomYSize; p2[1]=p1[1] + 0.05*RoomYSize;
      PG_.AddAC(p1,p2,3,10);
      p1[0]=1.0*RoomXSize; p1[1]=0.4*RoomYSize; p2[0]=1.0*RoomXSize; p2[1]=0.6*RoomYSize;
      PG_.AddExhaust(p1,p2);
      for(int i=0;i<2;i++) {
        p1[0]=(0.2+0.4*i)*RoomXSize; p2[0]=p1[0] + 0.2*RoomXSize;
        for(int j=0;j<=i;j++) {
          p1[1]=(0.2+0.4*j)*RoomYSize; p2[1]=p1[1] + 0.2*RoomYSize;
          PG_.AddEquipment(p1,p2,35,1);
        }
      }
      break;
    case 4: // What, if I do it the other way arround? in fact higher velocity needed
      p1[0]=RoomXSize=1.0; p1[1]=RoomYSize=0.2;
      PG_._RoomSize = p1;
      p1[0]=0.0*RoomXSize; p2[0]=0.0*RoomXSize; p1[1]=0.025*RoomYSize; p2[1]=p1[1] + 0.05*RoomYSize;
      PG_.AddAC(p1,p2,3,10);
      p1[0]=1.0*RoomXSize; p1[1]=0.4*RoomYSize; p2[0]=1.0*RoomXSize; p2[1]=0.6*RoomYSize;
      PG_.AddExhaust(p1,p2);
      for(int i=0;i<2;i++) {
        p1[0]=(0.2+0.4*i)*RoomXSize; p2[0]=p1[0] + 0.2*RoomXSize;
        for(int j=0;j<=i;j++) {
          p1[1]=(0.2+0.4*j)*RoomYSize; p2[1]=p1[1] + 0.2*RoomYSize;
          PG_.AddEquipment(p1,p2,35,1);
        }
      }
      break;
    case 5:
      p1.resize(3); p2.resize(3);
      p1[0]=RoomXSize=1.0; p1[1]=RoomYSize=1.0; p1[2]=RoomZSize=1.0;
      PG_._RoomSize = p1;
      p1[0]=0.0*RoomXSize; p1[1]=0.1*RoomYSize; p1[2]=0.2*RoomZSize;
      p2[0]=0.0*RoomXSize; p2[1]=0.3*RoomYSize; p2[2]=0.4*RoomZSize;
      PG_.AddAC(p1,p2,3,10);
      p1[0]=0.0*RoomXSize; p1[1]=0.7*RoomYSize; p1[2]=0.2*RoomZSize;
      p2[0]=0.0*RoomXSize; p2[1]=0.9*RoomYSize; p2[2]=0.3*RoomZSize;
      //PG_.AddAC(p1,p2,1,10);
      p1[0]=1.0*RoomXSize; p1[1]=0.1*RoomYSize; p1[2]=0.6*RoomZSize;
      p2[0]=1.0*RoomXSize; p2[1]=0.3*RoomYSize; p2[2]=0.8*RoomZSize;
      PG_.AddExhaust(p1,p2);
      p1[0]=1.0*RoomXSize; p1[1]=0.7*RoomYSize; p1[2]=0.6*RoomZSize;
      p2[0]=1.0*RoomXSize; p2[1]=0.9*RoomYSize; p2[2]=0.8*RoomZSize;
      //PG_.AddExhaust(p1,p2);
      p1[0]=0.6*RoomXSize; p1[1]=0.2*RoomYSize; p1[2]=0.0*RoomZSize;
      p2[0]=0.8*RoomXSize; p2[1]=0.4*RoomYSize; p2[2]=0.6*RoomZSize;
      //PG_.AddEquipment(p1,p2,35,1);
      p1[0]=0.2*RoomXSize; p1[1]=0.5*RoomYSize; p1[2]=0.0*RoomZSize;
      p2[0]=0.7*RoomXSize; p2[1]=0.85*RoomYSize; p2[2]=0.8*RoomZSize;
      PG_.AddEquipment(p1,p2,35,1);
      break;
    default: std::cout << "Error: unkown constellation" << std::endl; exit(0);
  }
#endif
  PG_.CreateMesh(&mesh_);

  WriteNodeFile(mesh_, "MeshGen.node");
  WriteEleFile(mesh_, "MeshGen.ele");
}

void LibMeshPDEBase::reinit()
{
  clear_math_obj();
  lm_eqn_sys_ = new EquationSystems(mesh_);
  
  std::string order("FIRST");
  std::string family("LAGRANGE"); 

  libMesh::LinearImplicitSystem& lm_sys = lm_eqn_sys_->add_system<LinearImplicitSystem>("PDE");
  lm_sys.add_variable("Phi", Utility::string_to_enum<Order>(order), Utility::string_to_enum<FEFamily>(family));
  lm_sys.attach_assemble_function(LibMeshPDEBase::assemble_Phi_PDE);
  lm_eqn_sys_->parameters.set<LibMeshPDEBase*>("LibMeshPDEBase") = this;
  lm_eqn_sys_->parameters.set<bool>("b_struct_only") =false;
  lm_eqn_sys_->init();

  int n_state_global, n_state_local;
  int m_pde_constr_global, m_pde_constr_local;
  {
    //nstate = lm_sys.matrix->n();
    n_state_global = lm_sys.matrix->n();
    PetscMatrix<Number> *pPetscMat = dynamic_cast<PetscMatrix<Number>*>(lm_sys.matrix);
    Mat mat = pPetscMat->mat();
    int start, end;
    MatGetOwnershipRangeColumn(mat, &start, &end );
    n_state_local = end-start;
    //int mlocal = lm_sys.matrix->row_stop() - lm_sys.matrix->row_start();
    m_pde_constr_global = lm_sys.matrix->m();
    m_pde_constr_local = lm_sys.matrix->row_stop()-lm_sys.matrix->row_start();
  }

  int n_control_global, n_control_local;
  {
    n_control_global = PG_._ParamIdx2BCParam.size();
    n_control_local = (0==GetProcID()) ? n_control_global : 0;
  }

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
  {
    //lm_control_vec_ = new PetscVector<Number>::PetscVector(n_control_global,n_control_local); does NOT work, later Petsc-calls never terminate (?)
    Vec petsc_vec;
    ierr = VecCreateMPI(PETSC_COMM_WORLD,n_control_local,PETSC_DETERMINE,&petsc_vec);
    CHKERRV(ierr);
    lm_control_vec_ = new PetscVector<Number>::PetscVector(petsc_vec);
    *lm_control_vec_ = 1.0;
    lm_control_vec_->close();
  }
  { // TODO: analyze nonzero structure more closely to pass tight numbers to init
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
    lm_control_vec_->close();
  }
  //  lm_pde_residual_vec_ = new PetscVector<Number>::PetscVector(m_pde_constr_global,m_pde_constr_local);
  SparseMatrix<Number> *pde_jac_state, *pde_jac_control;
  calcPDE_jacobians(pde_jac_state, pde_jac_control);
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
    std::list<Number>::iterator it = LocIneqFactList.begin();
    std::cout << "LocIneqFactList.size=" << LocIneqFactList.size() << std::endl;
    for( unsigned int iVal=0; iVal<LocIneqFactList.size(); ++iVal, ++it) {
      Vals[iVal] = min_airflow * (*it);
    }
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
}

void LibMeshPDEBase::DetroySelfOwnedLibMeshPetscMatrix(SparseMatrix<Number>*& matrix)
{
  if(NULL==matrix)
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
  if(NULL==vector)
    return;
  PetscVector<Number>* lm_petsc_vector = dynamic_cast<PetscVector<Number>*>(vector);
  assert(NULL!=lm_petsc_vector);  // Problem, if not a libMesh::PetscMatrix<Number>
  //Vec petsc_vector = lm_petsc_vector->vec();
  //TODO: VecDestroy(petsc_vector); does not work
  delete vector;
  vector = NULL;
}

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
  
  PetscInt sz;
  VecGetSize(vec_gathrd,&sz);
  PetscScalar *vals;
  ierr=VecGetArray(vec_gathrd,&vals);CHKERRV(ierr);
	int iControl = 0;
  assert(PG_._AC.size()==(unsigned int)sz);
	for(int iAC=0;iAC<sz;++iAC) {
		int BoundaryMarker = PG_._AC[iAC].BoundaryMarker[0];
		PG_._BoundCond[BoundaryMarker].PhiRhs = vals[iControl++]; 
	}
  ierr=VecRestoreArray(vec_gathrd,&vals);CHKERRV(ierr);
//  DBG_PRINT("LibMeshPDEBase::ConvertControl2PGData finished");
}

void LibMeshPDEBase::calc_objective_part(Number& Val)
{
	Val = 0.0;
	ConvertControl2PGData();
	if(GetProcID()==0)
	  for(unsigned int iAC=0;iAC<PG_._AC.size();iAC++) {
		  int BoundaryMarker = PG_._AC[iAC].BoundaryMarker[0];
      #ifdef QUAD_OBJ_FUNC
		    Val += (PG_._BoundCond[BoundaryMarker].PhiRhs)*(PG_._BoundCond[BoundaryMarker].PhiRhs);
      #else
        Val += PG_._BoundCond[BoundaryMarker].PhiRhs;
      #endif //QUAD_OBJ_FUNC
	  }
}

void LibMeshPDEBase::calc_objective_gradient(libMesh::NumericVector<libMesh::Number>& grad_state,
					     libMesh::NumericVector<libMesh::Number>& grad_control)
{
  DBG_PRINT("LibMeshPDEBase::calc_objective_gradient called");
	grad_state.zero();
	ConvertControl2PGData();
	if(GetProcID()==0)
	{
		for(unsigned int iAC=0;iAC<PG_._AC.size();iAC++) {
			#ifdef QUAD_OBJ_FUNC
        int BoundaryMarker = PG_._AC[iAC].BoundaryMarker[0];
			  grad_control.set(iAC,2.0*PG_._BoundCond[BoundaryMarker].PhiRhs);
      #else
        grad_control.set(iAC,1.0);
      #endif //QUAD_OBJ_FUNC
		}
	}
  DBG_PRINT("LibMeshPDEBase::calc_objective_gradient finished");
}

void LibMeshPDEBase::calcPDE_residual(libMesh::NumericVector<libMesh::Number>*& residual)
{
  DBG_PRINT("LibMeshPDEBase::calcPDE_residual called");
  ConvertControl2PGData();

	const double eps = 1e-8;
  lm_pde_residual_vec_->zero();
  std::vector<BoundaryCondition>& BCs=PG_._BoundCond;
  const MeshBase& mesh = lm_eqn_sys_->get_mesh();
  const unsigned int dim = mesh.mesh_dimension();
  LinearImplicitSystem& system = lm_eqn_sys_->get_system<LinearImplicitSystem>("PDE");
  const DofMap& dof_map = system.get_dof_map();
  FEType fe_type = dof_map.variable_type(0);

  AutoPtr<FEBase> fe (FEBase::build(dim, fe_type));	// this object will by the current volume object, actualy holding values for current phi and dphi
  QGauss qrule (dim, FIFTH);
  fe->attach_quadrature_rule (&qrule);
  AutoPtr<FEBase> fe_face (FEBase::build(dim, fe_type));	// surface object
  QGauss qface(dim-1, FIFTH);
  fe_face->attach_quadrature_rule (&qface);
  const std::vector<Real>& JxW = fe->get_JxW();			// References to values hold by fe (i.e. the chang, ones fe changes
  const std::vector<std::vector<Real> >& phi = fe->get_phi();
  const std::vector<std::vector<RealGradient> >& dphi = fe->get_dphi();

  DenseVector<Number> ElemRes;
  std::vector<unsigned int> dof_indices;	// mapping local index -> global index
	std::vector<Number> LocalSolution;			// solution values of current element
  MeshBase::const_element_iterator itCurEl = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator itEndEl = mesh.active_local_elements_end();
  for ( ; itCurEl != itEndEl; ++itCurEl)
    {
      const Elem* CurElem = *itCurEl;
      dof_map.dof_indices(CurElem, dof_indices);		// setup local->global mapping in form of array
      fe->reinit(CurElem);	// reinit fe to fit the current element, i.e. recalulate JxW and dPhi 
      ElemRes.resize(dof_indices.size());
			LocalSolution.resize(dof_indices.size());
			for(unsigned int iCurSol=0;iCurSol<LocalSolution.size();iCurSol++) {
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
        for (unsigned int side=0; side<CurElem->n_sides(); side++)
				{
					if (CurElem->neighbor(side) == NULL)
					{
						// std::cout << side << std::endl;
						short int bc_id = mesh.boundary_info->boundary_id (CurElem,side);
						assert(bc_id!=BoundaryInfo::invalid_id);
						double NeumCoef, DiriCoef, Rhs;
						// std::cout << "bc_id=" << bc_id << std::endl; 
						NeumCoef = BCs[bc_id].PhiNeumannCoef;
						DiriCoef = BCs[bc_id].PhiDiricheltCoef;
						Rhs = BCs[bc_id].PhiRhs;
						const std::vector<std::vector<Real> >&  phi_face = fe_face->get_phi();
						const std::vector<Real>& JxW_face = fe_face->get_JxW();
						fe_face->reinit(CurElem, side);
						/*Point Center;
						for(unsigned int iNode=0;iNode<CurElem->n_nodes();iNode++)
							if( CurElem->is_node_on_side(iNode,side) )
								Center += CurElem->point(iNode);
						Center = Center/dim;
						std::cout << "BC at " << Center << " and : ";
						std::cout << NeumCoef << " dPhi/dn = " << DiriCoef << " Phi + " << Rhs << std::endl;*/
						if( fabs(NeumCoef)>eps) {	// handle non-Dirichlet boundary conditions
							for (unsigned int qp=0; qp<qface.n_points(); qp++) {
									for (unsigned int i=0; i<phi_face.size(); i++)
										for (unsigned int j=0; j<phi_face.size(); j++)
											ElemRes(i) -= (DiriCoef/NeumCoef)*JxW_face[qp]*phi_face[i][qp]*LocalSolution[j];
									for (unsigned int i=0; i<phi_face.size(); i++)
										ElemRes(i) -= (Rhs/NeumCoef)*JxW_face[qp]*phi_face[i][qp];
							}
						}
					} // end if CurElem->neigbor(side)==NULL
				} // endfor side
				for(unsigned int i=0;i<CurElem->n_nodes();i++) { // handle Dirichlet BC for nodes: set whole line to 0, set DiriCoef on main diagonal and Rhs in vector
					std::vector<short int> bc_id = mesh.boundary_info->boundary_ids(CurElem->get_node(i));
					if(bc_id.empty())
						continue;
					if(bc_id[0]==BoundaryInfo::invalid_id)
						continue;
					if( fabs(BCs[bc_id[0]].PhiNeumannCoef) > eps )
						continue;
					// std::cout << "Dirichlet: bc_id=" << bc_id[0] << std::endl; 
					ElemRes(i) = BCs[bc_id[0]].PhiDiricheltCoef*LocalSolution[i];
					ElemRes(i) += BCs[bc_id[0]].PhiRhs;
				}
			}
			/*std::cout << "Assembled ElemMat: " << std::endl;
			std::cout << ElemMat;
			std::cout << "Element points: " << std::endl;
			std::cout << CurElem->point(0);
			std::cout << CurElem->point(1);
			std::cout << CurElem->point(2);
			std::cout.flush();*/
			dof_map.constrain_element_vector(ElemRes, dof_indices);	// Add constrains for hanging nodes (refinement)
			//std::cout << "Constrained ElemMat: " << std::endl;
			//std::cout << ElemMat;

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
	std::vector<BoundaryCondition>& BCs = PG_._BoundCond;	// Mapping boudary info (index) -> boundary condition, set up in Problem geometry class

  AutoPtr<FEBase> fe (FEBase::build(dim, fe_type));	// this object will by the current volume object, actualy holding values for current phi and dphi
  QGauss qrule (dim, FIFTH);
  fe->attach_quadrature_rule (&qrule);
  AutoPtr<FEBase> fe_face (FEBase::build(dim, fe_type));	// surface object
  QGauss qface(dim-1, FIFTH);
  fe_face->attach_quadrature_rule (&qface);

  DenseMatrix<Number> ElemMat;		// Matrices for current element
  std::vector<unsigned int> dof_indices;	// mapping local index -> global index
	std::vector<Number> LocalSolution;			// solution values of current element

  MeshBase::const_element_iterator itCurEl = mesh_.active_local_elements_begin();
  const MeshBase::const_element_iterator itEndEl = mesh_.active_local_elements_end();
  for ( ; itCurEl != itEndEl; ++itCurEl)
	{
		const Elem* CurElem = *itCurEl;
		dof_map.dof_indices(CurElem, dof_indices);		// setup local->global mapping in form of array
		fe->reinit (CurElem);	// reinit fe to fit the current element, i.e. recalulate JxW and dPhi 
		ElemMat.resize (dof_indices.size(),dof_indices.size());
		// Point Center = (CurElem->point(0)+CurElem->point(1)+CurElem->point(2))/3.0;
		// std::cout << "Element at " << Center << std::endl;
		for (unsigned int side=0; side<CurElem->n_sides(); side++)
		{
			if (CurElem->neighbor(side) == NULL)
			{
				// std::cout << side << std::endl;
				short int bc_id = mesh_.boundary_info->boundary_id (CurElem,side);
				assert(bc_id!=BoundaryInfo::invalid_id);
				for(unsigned int icontrol=0;icontrol<lm_control_vec_->size();icontrol++)
				{
					int BoundCntrl = PG_._ParamIdx2BCParam[icontrol].BoundaryMarker;
					if(BoundCntrl!=bc_id)
						continue;
					LocalSolution.resize(dof_indices.size());
					for(unsigned int iCurSol=0;iCurSol<LocalSolution.size();iCurSol++) {
            LocalSolution[iCurSol] = system.current_local_solution->el(dof_indices[iCurSol]);
					}
					double NeumCoef, DiriCoef, Rhs;
					// std::cout << "bc_id=" << bc_id << std::endl; 
					NeumCoef = BCs[bc_id].PhiNeumannCoef;
					DiriCoef = BCs[bc_id].PhiDiricheltCoef;
					Rhs = BCs[bc_id].PhiRhs;
					if(fabs(NeumCoef)>eps) {	// handle non-Dirichlet boundary conditions
						const std::vector<std::vector<Real> >&  phi_face = fe_face->get_phi();
						const std::vector<Real>& JxW_face = fe_face->get_JxW();
						fe_face->reinit(CurElem, side);
						switch(PG_._ParamIdx2BCParam[icontrol].BCParameter)
						{	// 0:PhiDiriCoef, 1: PhiNeumCoef, 2: PhiRhs, 3:TDiriCoef, 4: TNeumCoef, 5: TiRhs
							case 0:	
                if(calc_type_ == StructureOnly) {
								  for (unsigned int qp=0; qp<qface.n_points(); qp++)
									  for (unsigned int i=0; i<phi_face.size(); i++)
										  for (unsigned int j=0; j<phi_face.size(); j++){
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
                if(calc_type_ == StructureOnly) {
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
                if(calc_type_ == StructureOnly) {
                  for (unsigned int qp=0; qp<qface.n_points(); qp++) {
                    for (unsigned int i=0; i<phi_face.size(); i++)
                      jac_control_->add(dof_indices[i], icontrol,1.0);
                  }
                }
                else {
                  for (unsigned int qp=0; qp<qface.n_points(); qp++) {
                    for (unsigned int i=0; i<phi_face.size(); i++)
                      jac_control_->add(dof_indices[i], icontrol,-(1.0/NeumCoef)*JxW_face[qp]*phi_face[i][qp]);
                  }
                }
								break;
						}
					} // end if(fabs(NeumCoef)>eps)
				} // end Loop over Controls
			} // end if boundary side
		} // end side loop
		for(unsigned int i=0;i<CurElem->n_nodes();i++) { // handle Dirichlet BC for nodes: set whole line to 0, set DiriCoef on main diagonal and Rhs in vector
			std::vector<short int> bc_id = mesh_.boundary_info->boundary_ids(CurElem->get_node(i));
			if(bc_id.empty())
				continue;
			if(bc_id[0]==BoundaryInfo::invalid_id)
				continue;
			if( fabs(BCs[bc_id[0]].PhiNeumannCoef) > eps )
				continue;
			for(unsigned int iControl=0;iControl<lm_control_vec_->size();iControl++)
			{
				int BoundCntrl = PG_._ParamIdx2BCParam[iControl].BoundaryMarker;
				if(BoundCntrl!=bc_id[0])
					continue;
				switch(PG_._ParamIdx2BCParam[iControl].BCParameter)
				{	// 0:PhiDiriCoef, 1: PhiNeumCoef, 2: PhiRhs, 3:TDiriCoef, 4: TNeumCoef, 5: TiRhs
					case 0:	
            if(calc_type_ == StructureOnly)
              jac_control_->set(dof_indices[i], iControl,1.0);
            else
						  jac_control_->set(dof_indices[i], iControl,system.current_local_solution->el(dof_indices[i]));
						break;
					case 1:
						assert(false); // derivative would be ... / NeumCoef^2, but NeumCoef < 1e-8
					case 2:
            if(calc_type_ == StructureOnly)
              jac_control_->set(dof_indices[i], iControl,1.0);
            else
						  jac_control_->set(dof_indices[i], iControl,1.0);
				}
			}
		}
	// dof_map.constrain_element_matrix(ElemMat, dof_indices);	//There are no hangin nodes on the boundary, thus, no constraints needed 
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
	std::vector<BoundaryCondition>& BCs = pData->PG_._BoundCond;	// Mapping boudary info (index) -> boundary condition, set up in Problem geometry class
  CalculationModeType calc_type=pData->calc_type_;

  AutoPtr<FEBase> fe (FEBase::build(dim, fe_type));	// this object will by the current volume object, actualy holding values for current phi and dphi
  QGauss qrule (dim, FIFTH);
  fe->attach_quadrature_rule (&qrule);
  AutoPtr<FEBase> fe_face (FEBase::build(dim, fe_type));	// surface object
  QGauss qface(dim-1, FIFTH);
  fe_face->attach_quadrature_rule (&qface);
  const std::vector<Real>& JxW = fe->get_JxW();			// References to values hold by fe (i.e. the chang, ones fe changes
  const std::vector<std::vector<Real> >& phi = fe->get_phi();
  const std::vector<std::vector<RealGradient> >& dphi = fe->get_dphi();

  DenseMatrix<Number> ElemMat;		// Matrices for current element
  std::vector<unsigned int> dof_indices;	// mapping local index -> global index

  MeshBase::const_element_iterator itCurEl = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator itEndEl = mesh.active_local_elements_end();
  for ( ; itCurEl != itEndEl; ++itCurEl)
	{
		const Elem* CurElem = *itCurEl;
		dof_map.dof_indices(CurElem, dof_indices);		// setup local->global mapping in form of array
		fe->reinit (CurElem);	// reinit fe to fit the current element, i.e. recalulate JxW and dPhi 
		ElemMat.resize (dof_indices.size(),dof_indices.size());
		// Point Center = (CurElem->point(0)+CurElem->point(1)+CurElem->point(2))/3.0;
		// std::cout << "Element at " << Center << std::endl;
		unsigned int qp,i,j;
		for (qp=0; qp<qrule.n_points(); qp++)
			for (i=0; i<phi.size(); i++)
				for (j=0; j<phi.size(); j++)
          if(calc_type==StructureOnly)
            ElemMat(i,j) += 1.0;
          else
					  ElemMat(i,j) += JxW[qp]*(dphi[i][qp]*dphi[j][qp]);
		
		for (unsigned int side=0; side<CurElem->n_sides(); side++)
		{
			if (CurElem->neighbor(side) == NULL)
			{
				// std::cout << side << std::endl;
				short int bc_id = mesh.boundary_info->boundary_id (CurElem,side);
				assert(bc_id!=BoundaryInfo::invalid_id);
				double NeumCoef, DiriCoef, Rhs;
				// std::cout << "bc_id=" << bc_id << std::endl; 
				NeumCoef = BCs[bc_id].PhiNeumannCoef;
				DiriCoef = BCs[bc_id].PhiDiricheltCoef;
				Rhs = BCs[bc_id].PhiRhs;
				const std::vector<std::vector<Real> >&  phi_face = fe_face->get_phi();
				const std::vector<Real>& JxW_face = fe_face->get_JxW();
				fe_face->reinit(CurElem, side);
				/*Point Center;
				for(unsigned int iNode=0;iNode<CurElem->n_nodes();iNode++)
					if( CurElem->is_node_on_side(iNode,side) )
						Center += CurElem->point(iNode);
				Center = Center/dim;
				std::cout << "BC at " << Center << " and : ";
				std::cout << NeumCoef << " dPhi/dn = " << DiriCoef << " Phi + " << Rhs << std::endl;*/
				if( fabs(NeumCoef)>eps) {	// handle non-Dirichlet boundary conditions
					for (unsigned int qp=0; qp<qface.n_points(); qp++) {
							for (unsigned int i=0; i<phi_face.size(); i++)
								for (unsigned int j=0; j<phi_face.size(); j++)
                  if(calc_type==StructureOnly)
                    ElemMat(i,j) += 1.0;
                  else
									  ElemMat(i,j) -= (DiriCoef/NeumCoef)*JxW_face[qp]*phi_face[i][qp];
					}
				}
			} // end if CurElem->neigbor(side)==NULL
		} // endfor side
		for(unsigned int i=0;i<CurElem->n_nodes();i++) { // handle Dirichlet BC for nodes: set whole line to 0, set DiriCoef on main diagonal and Rhs in vector
			std::vector<short int> bc_id = mesh.boundary_info->boundary_ids(CurElem->get_node(i));
			if(bc_id.empty())
				continue;
			if(bc_id[0]==BoundaryInfo::invalid_id)
				continue;
			if( fabs(BCs[bc_id[0]].PhiNeumannCoef) > eps )
				continue;
			// std::cout << "Dirichlet: bc_id=" << bc_id[0] << std::endl; 
			for(j=0;j<CurElem->n_nodes();j++)
				ElemMat(i,j)=0;
      if(calc_type==StructureOnly)
        ElemMat(i,i) = 1.0;
      else
  			ElemMat(i,i) = +BCs[bc_id[0]].PhiDiricheltCoef;
		}
		dof_map.constrain_element_matrix(ElemMat, dof_indices);	// Add constrains for hanging nodes (refinement)
		system.matrix->add_matrix(ElemMat, dof_indices);
	}
  DBG_PRINT( "LibMeshPDEBase::assemble_Phi_PDE finished" );
}

// sigma and lambda_pde not used, avoid warning
void LibMeshPDEBase::calc_hessians(Number /*sigma*/, libMesh::DenseVector<Number>& /*lambda_pde*/,
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
        if( AuxConstrBoundMarkerList_.find(bc_id)!=AuxConstrBoundMarkerList_.end() ) {  // Heating boundary, apply min. velocity-constraint
          fe_face->reinit(CurElem, side);
          dof_map.dof_indices(CurElem, dof_indices);   // setup local->global mapping in form of array
          DenseMatrix<Number> loc_l2dphi;
          loc_l2dphi.resize(dof_indices.size(),dof_indices.size());
          if(calc_type_==StructureOnly) {
            for(unsigned short inode=0;inode<CurElem->n_nodes();++inode) {
              for(unsigned short jnode=0;jnode<inode;++jnode) {
                loc_l2dphi(inode,jnode) = 1.0;
              }
            }
            hess_state_state_->add_matrix(loc_l2dphi,dof_indices);
          }
          else {
            DenseMatrix<Number> loc_l2dphi;
            loc_l2dphi.resize(dof_indices.size(),dof_indices.size());
            DenseVector<Number> tmp;
            for (unsigned int qp=0; qp<qface.n_points(); qp++) {
              for(unsigned short inode=0;inode<CurElem->n_nodes();++inode) {
                for( unsigned short jnode=0;jnode<=inode;++jnode) { // only lower triangle 
                  loc_l2dphi(inode,jnode) = dphi_face[inode][qp]*dphi_face[jnode][qp];
                }
              }
              loc_l2dphi *= 2.0*lambda_aux(i_aux_constr)*JxW_face[qp];
              hess_state_state_->add_matrix(loc_l2dphi,dof_indices);
            }
          }
          i_aux_constr++;
        }
      }
    }
  }

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
  if(simulation_mode_) {
    control_l=*lm_control_vec_;
    control_u=*lm_control_vec_;
    aux_constr_l = -100;
  }
  else {
    control_l=0.0;
    control_u=Inf;
    aux_constr_l = *lm_aux_constr_vec_low_bd_; // lower bound for tangetial flow on equipment boundary
  }
  aux_constr_u = Inf;
  DBG_PRINT( "LibMeshPDE::get_bounds finished" );
}

void LibMeshPDEBase::get_starting_point(libMesh::NumericVector<libMesh::Number>& state,
        libMesh::NumericVector<libMesh::Number>& control)
{
  DBG_PRINT( "LibMeshPDE::get_starting_point called" );
  state = getStateVector();
  control = getControlVector();
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
  for(;CurNode!=EndNode;CurNode++)
  {
    os << "  " << iNode++ << " "; 
    for(iDim=0;iDim<Dim;iDim++)
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
  for(;CurEl!=EndEl;CurEl++)
  {
    os << "  " << iEl++ << " "; 
    for(iNode=0;iNode<((*CurEl)->n_nodes());iNode++)
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
  std::string my_pre_filename = pre_filename;
  if(simulation_mode_)
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
  for(unsigned int iVal=0;iVal<State.size();iVal++)
    f << State[iVal] << std::endl;
  f.close();
}

void LibMeshPDEBase::InitAuxConstr(int *plocal, int *pglobal, std::list<Number>* pFactList)
{
  DBG_PRINT( "LibMeshPDEBase::InitAuxConstr called" );
  const double eps = 1e-8;
  std::vector<ProblemGeometry::Item> _Exh;

  int nEquipWalls = 4;
  if(mesh_.mesh_dimension()==3)
    nEquipWalls = 5;
  int BoundaryMarker;
  AuxConstrBoundMarkerList_.clear();
  pFactList->clear();
  for(unsigned int iEquip=0; iEquip<PG_._Equip.size();++iEquip)  {
    for(int iWall=0; iWall<nEquipWalls; iWall++)  {
      BoundaryMarker = PG_._Equip[iEquip].BoundaryMarker[iWall];
      // We notice, that ineq. const old here, by checking, if we have a real robin boundary condition for T, i.e. all coefficients are nonzero
      if(BoundaryMarker<0)  // in 3D, floor of equipment
        continue;
      if(fabs(PG_._BoundCond[BoundaryMarker].TDiricheltCoef
              *PG_._BoundCond[BoundaryMarker].TNeumannCoef
              *PG_._BoundCond[BoundaryMarker].TRhs) > eps) {
        AuxConstrBoundMarkerList_.insert(BoundaryMarker);
      }
    }
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
      if (CurElem->neighbor(side) == NULL)
      {
        int bc_id = mesh_.boundary_info->boundary_id (CurElem,side);
        if( AuxConstrBoundMarkerList_.find(bc_id)!=AuxConstrBoundMarkerList_.end() ) {
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
        if( AuxConstrBoundMarkerList_.find(bc_id)!=AuxConstrBoundMarkerList_.end() ) {  // Heating boundary, apply min. velocity-constraint
          fe_face->reinit(CurElem, side);
          dof_map.dof_indices(CurElem, dof_indices);   // setup local->global mapping in form of array
          DenseVector<Number> loc_sol;
          DenseMatrix<Number> loc_l2dphi;
          loc_sol.resize(dof_indices.size());
          loc_l2dphi.resize(dof_indices.size(),dof_indices.size());
          double GradL2=0.0;
          for (unsigned int qp=0; qp<qface.n_points(); qp++) {
            for(unsigned short inode=0;inode<CurElem->n_nodes();++inode) {
              loc_sol(inode) = system.current_local_solution->el(dof_indices[inode]);
              for( unsigned short jnode=0;jnode<CurElem->n_nodes();++jnode) {
                loc_l2dphi(inode,jnode) = dphi_face[inode][qp]*dphi_face[jnode][qp];
              }
            }
            DenseVector<Number> tmp;
            loc_l2dphi.vector_mult(tmp,loc_sol);
            GradL2 += JxW_face[qp]*loc_sol.dot(tmp);
          }
          lm_aux_constr_vec_->set(i_aux_constr,GradL2);
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
        if( AuxConstrBoundMarkerList_.find(bc_id)!=AuxConstrBoundMarkerList_.end() ) {  // Heating boundary, apply min. velocity-constraint
          fe_face->reinit(CurElem, side);
          dof_map.dof_indices(CurElem, dof_indices);   // setup local->global mapping in form of array
          if(calc_type_==StructureOnly)
            for(unsigned short inode=0;inode<CurElem->n_nodes();++inode) {
              jac_aux_state_->add(i_aux_constr,dof_indices[inode],1.0);
          }
          else {
            DenseVector<Number> loc_sol;
            DenseMatrix<Number> loc_l2dphi;
            loc_sol.resize(dof_indices.size());
            loc_l2dphi.resize(dof_indices.size(),dof_indices.size());
            DenseVector<Number> tmp;
            for (unsigned int qp=0; qp<qface.n_points(); qp++) {
              for(unsigned short inode=0;inode<CurElem->n_nodes();++inode) {
                loc_sol(inode) = system.current_local_solution->el(dof_indices[inode]);
                for( unsigned short jnode=0;jnode<CurElem->n_nodes();++jnode) {
                  loc_l2dphi(inode,jnode) = dphi_face[inode][qp]*dphi_face[jnode][qp];
                }
              }
              loc_l2dphi.vector_mult(tmp,loc_sol);
              for(unsigned short inode=0;inode<CurElem->n_nodes();++inode) {
                jac_aux_state_->add(i_aux_constr,dof_indices[inode],2.0*JxW_face[qp]*tmp(inode));
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
  DBG_PRINT( "LibMeshPDEBase::calcAux_jacobian_finished called" );
}

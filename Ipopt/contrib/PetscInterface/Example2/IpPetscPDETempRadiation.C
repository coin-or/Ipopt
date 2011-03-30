// Copyright (C) 2010 University Basel
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Authors:  Johannes Huber

#include "mpi.h"
#include "IpPetscPDETempRadiation.hpp"
#include "petsc.h"

#include "string_to_enum.h"
#include "fe.h"
#include "quadrature_gauss.h"
#include "dof_map.h"
#include "petsc_matrix.h"
#include "petscvec.h"
#include "numeric_vector.h"
#include "dense_matrix.h"
#include "dense_vector.h"
#include "elem.h"
#include "boundary_info.h"
#include "point_locator_tree.h"
#include "fe_interface.h"
#include "mesh_generation.h"
#include "exodusII_io.h"

#include <fstream>

int GetProcID();

// PDE boundary condition: dy/dn = alpha * (y^4-u^4) -> CONTROL_EXPONENT=4
// control = u^2 -> dy/dn = alpha * (y^4-u^4) = alpha (y^4-Control^2) -> CONTROL_EXPONENT=2
// control = u^4 -> dy/dn = alpha * (y^4-u^4) = alpha (y^4-Control  ) -> CONTROL_EXPONENT=1

#define CONTROL_EXPONENT 1

#define PRINT_LEVEL 10
//#define MY_MY_DBG_PRINT(s) {std::cout << GetProcID() << ": " << s << std::endl;}
#define MY_DBG_PRINT(s) {}
#define START_FUNCTION { if(PRINT_LEVEL>10) {std::cout << "[" << GetProcID() << "]" << __FILE__ << ":" << __LINE__ <<":" << "called " << __func__ << std::endl;} }
#define END_FUNCTION { if(PRINT_LEVEL>10) {std::cout << "[" << GetProcID() << "]" << __FILE__ << ":" << __LINE__ <<":" << "finished " << __func__ << std::endl;} }
#define CHECKPOINT { std::cout << GetProcID() << ", " << __func__ << ", Line " << __LINE__ << std::endl; }

using namespace libMesh;

PetscPDETempRadiation::PetscPDETempRadiation() :
    m_SolutionMode(Optimize),
    m_Control(NULL),
    m_State(NULL),
    m_PDEConstr(NULL),
    m_AuxConstr(NULL),
    m_GradObjControl(NULL),
    m_GradObjState(NULL),
    m_JacPDEState(NULL),
    m_JacPDEControl(NULL),
    m_JacAuxState(NULL),
    m_JacAuxControl(NULL),
    m_HessControlControl(NULL),
    m_HessControlState(NULL),
    m_HessStateState(NULL),
    m_Dim(3),
    m_StateMesh(m_Dim),
    m_FEOrder(FIRST),
    m_StateDofMap(NULL),
    m_PDEConstrScale(1.0),
    m_Alpha(1.0),
    m_HotMinTemp1(2.0),
    m_HotMinTemp2(1.0),
    m_ReguStateParam(0.0)
{
  m_Omega1[0] = 0.2;
  m_Omega1[1] = 0.3;
  m_Omega1[2] = m_Dim==3 ? 0.4 : 0.0;
  m_Omega1[3] = 0.5;
  m_Omega1[4] = 0.6;
  m_Omega1[5] = m_Dim==3 ? 0.7 : 0.0;
  m_Omega2[0] = 0.1;
  m_Omega2[1] = 0.2;
  m_Omega2[2] = 0.3;
  m_Omega2[3] = 0.6;
  m_Omega2[4] = 0.7;
  m_Omega2[5] = 0.8;  
}


PetscPDETempRadiation::~PetscPDETempRadiation()
{
  START_FUNCTION
  clear_math_obj();
  
  END_FUNCTION
}

// clear all matrices and vectors, but not problem geomatry data
// used e.g. after refinement (since problem size changes)
void PetscPDETempRadiation::clear_math_obj()
{
  START_FUNCTION
  DetroySelfOwnedLibMeshPetscVector(m_Control);
  DetroySelfOwnedLibMeshPetscVector(m_State);
  DetroySelfOwnedLibMeshPetscVector(m_PDEConstr);
  DetroySelfOwnedLibMeshPetscVector(m_AuxConstr);
  DetroySelfOwnedLibMeshPetscVector(m_GradObjControl);
  DetroySelfOwnedLibMeshPetscVector(m_GradObjState);
  DetroySelfOwnedLibMeshPetscMatrix(m_JacPDEState);
  DetroySelfOwnedLibMeshPetscMatrix(m_JacPDEControl);
  DetroySelfOwnedLibMeshPetscMatrix(m_JacAuxState);
  DetroySelfOwnedLibMeshPetscMatrix(m_JacAuxControl);
  DetroySelfOwnedLibMeshPetscMatrix(m_HessControlControl);
  DetroySelfOwnedLibMeshPetscMatrix(m_HessControlState);
  DetroySelfOwnedLibMeshPetscMatrix(m_HessStateState);
  END_FUNCTION
}

void PetscPDETempRadiation::DetroySelfOwnedLibMeshPetscMatrix(AutoPtr<PetscMatrix<Number> >& matrix)
{
  START_FUNCTION
  if (NULL==matrix.get())
    return;
  // Mat petsc_matrix = matrix->mat();
  // MatDestroy(petsc_matrix);
  matrix.reset(NULL);
  END_FUNCTION
}

void PetscPDETempRadiation::DetroySelfOwnedLibMeshPetscVector(AutoPtr<PetscVector<Number> >& vector)
{
  START_FUNCTION
  if (NULL==vector.get())
    return;
  vector->close();
  // Vec petsc_vector = vector->vec();
  // VecDestroy(petsc_vector);
  vector.reset(NULL);
  END_FUNCTION
}

double PetscPDETempRadiation::GetHotMinTemp1() {return m_HotMinTemp1;}
double PetscPDETempRadiation::GetHotMinTemp2() {return m_HotMinTemp2;}

void PetscPDETempRadiation::Init(const std::string filename)
{
  START_FUNCTION
  
  std::ifstream f("Problem.dat",std::ios::in);
  if (f.fail()) {
    std::cout << "Can't open file Problem.dat" << std::endl;
    exit(1);
  }

  double h(0.1);
  m_Dim = 3;
  {
    static char Buf[1024];
    int n;
    double Vals[8];
    while (!f.eof()) {
      f.getline(Buf,1024);
      if (strlen(Buf)<1)
        continue;
      if (Buf[0]=='#')
        continue;
      n = sscanf(Buf,"h=%lf",Vals);
      if (n==1) {
        h = Vals[0];
      }
      else {
        n = sscanf(Buf,"Omega1=%lf,%lf,%lf,%lf,%lf,%lf",Vals+0,Vals+1,Vals+2,Vals+3,Vals+4,Vals+5);
        if(n>0) {
          if(n==4) {
            if((m_Dim==0) || (m_Dim==2))
              m_Dim = 2;
            else {
              std::cout << "Dimension of Omega1 and Omega2 not constsitent" << std::endl;
              std::cout << "Problem interpreting line:" << std::endl;
              std::cout << Buf << std::endl;
              exit(1);
            }
            m_Omega1[0] = Vals[0];
            m_Omega1[1] = Vals[1];
            m_Omega1[2] = 0.0;
            m_Omega1[3] = Vals[2];
            m_Omega1[4] = Vals[3];
            m_Omega1[5] = 0.0;
          }
          else if(n==6) {
            if((m_Dim==0) || (m_Dim==3))
              m_Dim = 3;
            else {
              std::cout << "Dimension of Omega1 and Omega2 not constsitent" << std::endl;
              std::cout << "Problem interpreting line:" << std::endl;
              std::cout << Buf << std::endl;
              exit(1);
            }
            m_Omega1[0] = Vals[0];
            m_Omega1[1] = Vals[1];
            m_Omega1[2] = Vals[2];
            m_Omega1[3] = Vals[3];
            m_Omega1[4] = Vals[4];
            m_Omega1[5] = Vals[5];
          }
          else {
            std::cout << "Dimension not implemented, must be 2 or 3" << std::endl;
            std::cout << "Problem interpreting line:" << std::endl;
            std::cout << Buf << std::endl;
            exit(1);
          }
        }
        else {
          n = sscanf(Buf,"Omega2=%lf,%lf,%lf,%lf,%lf,%lf",Vals+0,Vals+1,Vals+2,Vals+3,Vals+4,Vals+5);
          if(n>0) {
            if(n==4) {
              if((m_Dim==0) || (m_Dim==2))
                m_Dim = 2;
              else {
                std::cout << "Dimension of Omega1 and Omega2 not constsitent" << std::endl;
                std::cout << "Problem interpreting line:" << std::endl;
                std::cout << Buf << std::endl;
                exit(1);
              }
              m_Omega2[0] = Vals[0];
              m_Omega2[1] = Vals[1];
              m_Omega2[2] = 0.0;
              m_Omega2[3] = Vals[2];
              m_Omega2[4] = Vals[3];
              m_Omega2[5] = 0.0;
            }
            else if(n==6) {
              if((m_Dim==0) || (m_Dim==3))
                m_Dim = 3;
              else {
                std::cout << "Dimension of Omega1 and Omega2 not constsitent" << std::endl;
                std::cout << "Problem interpreting line:" << std::endl;
                std::cout << Buf << std::endl;
                exit(1);
              }
              m_Omega2[0] = Vals[0];
              m_Omega2[1] = Vals[1];
              m_Omega2[2] = Vals[2];
              m_Omega2[3] = Vals[3];
              m_Omega2[4] = Vals[4];
              m_Omega2[5] = Vals[5];
            }
            else {
              std::cout << "Dimension not implemented, must be 2 or 3" << std::endl;
              std::cout << "Problem interpreting line:" << std::endl;
              std::cout << Buf << std::endl;
              exit(1);
            }
          }
          else
          {
            n = sscanf(Buf,"MinTemp1=%lf",Vals);
            if(1==n) {
              m_HotMinTemp1 = Vals[0];
            }
            else {
              n = sscanf(Buf,"MinTemp2=%lf",Vals);
              if(1==n) {
                m_HotMinTemp2 = Vals[0];
              }
              else {
                n = sscanf(Buf,"StateReg=%lf",Vals);
                if(1==n) {
                  m_ReguStateParam = Vals[0];
                }
                else {
                  std::cout << "Can't interprete line:" << std::endl;
                  std::cout << Buf << std::endl;
                  exit(1);
                }
              }
            }
          }
        }
      }
    }
  }

  {
    // Create Mesh for state variable
    int n = (int)(1/h);
    if(m_Dim==2) {
      MeshTools::Generation::build_cube (m_StateMesh,n,n, 0,
                                         0, 1, 0, 1, 0, 0, 
                                         m_FEOrder==FIRST ? TRI3 : TRI6);
    } else {
      MeshTools::Generation::build_cube (m_StateMesh,n,n,n,
                                         0, 1, 0, 1, 0, 1, 
                                         m_FEOrder==FIRST ? TET4 : TET10);
    }
    // Setup State degrees of freedom mapping
    MeshBase::node_iterator       node_it  = m_StateMesh.nodes_begin();
    const MeshBase::node_iterator node_end = m_StateMesh.nodes_end();
    for ( ; node_it != node_end; ++node_it)
      (*node_it)->add_system();
    MeshBase::element_iterator       elem_it  = m_StateMesh.elements_begin();
    const MeshBase::element_iterator elem_end = m_StateMesh.elements_end();
    for ( ; elem_it != elem_end; ++elem_it)
      (*elem_it)->add_system();
    m_StateDofMap = AutoPtr<DofMap>(new DofMap(0));
    std::string Varname("State");
    m_StateDofMap->add_variable(Variable(Varname, 0, FEType(m_FEOrder)));
    m_StateDofMap->distribute_dofs(m_StateMesh);
    m_StateDofMap->create_dof_constraints(m_StateMesh);
    m_StateDofMap->process_constraints();
    m_StateDofMap->prepare_send_list();
  }  
  unsigned int iControlDof(0);
  {
    m_ControlNodeIDToControlDOF.clear();
    std::vector<bool> bNodeDofAlreadySet;
    bNodeDofAlreadySet.resize(m_StateMesh.n_nodes(),false);
    MeshBase::element_iterator       elem_it  = m_StateMesh.elements_begin();
    const MeshBase::element_iterator elem_end = m_StateMesh.elements_end();
    for ( ; elem_it != elem_end; ++elem_it) {
      for (unsigned int iside=0; iside<(*elem_it)->n_sides(); iside++) {
        if ((*elem_it)->neighbor(iside) == NULL) {
          AutoPtr<libMesh::Elem> Side = (*elem_it)->build_side(iside);
          for(unsigned short iNode=0; iNode<Side->n_nodes(); iNode++) {
            if(!bNodeDofAlreadySet[Side->node(iNode)]) {
              bNodeDofAlreadySet[Side->node(iNode)] = true;
              m_ControlNodeIDToControlDOF.insert(std::pair<unsigned int,unsigned int>(Side->node(iNode),iControlDof++));
            }
          }
        }
      }
    }
  }
  
  int NumGlobStates(m_StateDofMap->n_dofs()), NumLocStates(m_StateDofMap->n_local_dofs());
  int NumGlobPDEConstr(NumGlobStates), NumLocPDEConstr(NumLocStates);
  int NumGlobControls(iControlDof), NumLocControls(iControlDof);
  int NumGlobAuxConstr(0), NumLocAuxConstr(0);
  
  std::cout <<  GetProcID() << ": order: " << m_FEOrder << std::endl << 
  ", # States:      " << NumGlobStates    << ", local: " << NumLocStates << std::endl <<
  ", # Controls:    " << NumGlobControls  << ", local: " << NumLocControls << std::endl <<
  ", # PDE constr:  " << NumGlobPDEConstr << ", local: " << NumLocPDEConstr << std::endl <<
  ", # Aux. constr: " << NumGlobAuxConstr << ", local: " << NumLocAuxConstr << std::endl;

  m_State =   AutoPtr<PetscVector<Number> >(new PetscVector<Number>(NumGlobStates, NumLocStates, m_StateDofMap->get_send_list(), GHOSTED));
  m_Control = AutoPtr<PetscVector<Number> >(new PetscVector<Number>(NumGlobControls, NumLocControls));
  m_PDEConstrScale = 1/pow(h,m_StateMesh.mesh_dimension()-2);
  m_PDEConstr = AutoPtr<PetscVector<Number> >(new PetscVector<Number>(NumGlobPDEConstr, NumLocPDEConstr));
  m_AuxConstr = AutoPtr<PetscVector<Number> >(new PetscVector<Number>(NumGlobAuxConstr, NumLocAuxConstr));
  Vec tmp = ((PetscVector<double>*)m_Control.get())->vec();
  int sz(9);
  VecGetSize(tmp,&sz);
  *m_Control = 0;
  m_Control->close();
  {
    Vec tmp;
    VecDuplicate(m_Control->vec(),&tmp);
    m_GradObjControl = AutoPtr<PetscVector<Number> >(new PetscVector<Number>(tmp));
  }
  {
    Vec tmp;
    VecDuplicate(m_State->vec(),&tmp);
    m_GradObjState = AutoPtr<PetscVector<Number> >(new PetscVector<Number>(tmp));
  }
  // Setup PDE jacobian including constraint distribution
  m_JacPDEState = AutoPtr<PetscMatrix<Number> >(new PetscMatrix<Number>);
  m_StateDofMap->attach_matrix(*m_JacPDEState);
  m_StateDofMap->compute_sparsity(m_StateMesh);
  m_JacPDEState->init();
  m_JacPDEControl = AutoPtr<PetscMatrix<Number> >(new PetscMatrix<Number>);
  m_JacPDEControl->init(NumGlobPDEConstr,NumGlobControls,NumLocPDEConstr,NumLocControls);
  m_JacAuxState = AutoPtr<PetscMatrix<Number> >(new PetscMatrix<Number>);
  m_JacAuxState->init(NumGlobAuxConstr,NumGlobStates,NumLocAuxConstr,NumLocStates);
  m_JacAuxControl = AutoPtr<PetscMatrix<Number> >(new PetscMatrix<Number>);
  m_JacAuxControl->init(NumGlobAuxConstr,NumGlobControls,NumLocAuxConstr,NumLocControls);
  {
    // each processor has hole matrix, but sets only part of local constraints
    Mat petsc_mat;
    MatCreateSeqAIJ(PETSC_COMM_SELF,NumGlobStates,NumGlobStates,8,PETSC_NULL,&petsc_mat);
    m_HessStateState = AutoPtr<PetscMatrix<Number> >(new PetscMatrix<Number>(petsc_mat));
  }
  {
    // each processor has hole matrix, but sets only part of local constraints
    Mat petsc_mat;  MY_DBG_PRINT("LibMeshPDEBase::calc_objective_gradient called");

    MatCreateSeqAIJ(PETSC_COMM_SELF,NumGlobControls,NumGlobStates,8,PETSC_NULL,&petsc_mat);
    m_HessControlState = AutoPtr<PetscMatrix<Number> >(new PetscMatrix<Number>(petsc_mat));
  }
  {
    // each processor has hole matrix, but sets only part of local constraints
    Mat petsc_mat;
    MatCreateSeqAIJ(PETSC_COMM_SELF,NumGlobControls,NumGlobControls,8,PETSC_NULL,&petsc_mat);
    m_HessControlControl = AutoPtr<PetscMatrix<Number> >(new PetscMatrix<Number>(petsc_mat));
  }
  END_FUNCTION
}


void PetscPDETempRadiation::get_bounds(Vec& state_l,
                          Vec& state_u,
                          Vec& control_l,
                          Vec& control_u,
                          Vec& aux_constr_l,
                          Vec& aux_constr_u)
{
  START_FUNCTION
  double Inf = 1e30;

  if( m_SolutionMode==Simulate ) {
    VecSet(state_l,0.0);
    VecSet(state_u,Inf);
    VecSet(control_l,1.0);
    VecSet(control_u,1.0);
  }
  else {
    VecSet(state_l,0.0);
    PetscVector<Number> lm_state_u(state_u);
    PetscVector<Number> lm_state_l(state_l);
    lm_state_u = Inf;
    lm_state_l = -Inf;
    MeshBase::node_iterator       node_it  = m_StateMesh.local_nodes_begin();
    const MeshBase::node_iterator node_end = m_StateMesh.local_nodes_end();
    for ( ; node_it != node_end; ++node_it) {
      if(IsInOmega1(**node_it)) {
        lm_state_l.set((*node_it)->dof_number(0,0,0),m_HotMinTemp1);
      }
      if(IsInOmega2(**node_it)) {
        lm_state_l.set((*node_it)->dof_number(0,0,0),m_HotMinTemp2);
      }
    }
    lm_state_u.close();
    
    VecSet(control_l,0.0);
    VecSet(control_u,Inf);
  }
  VecSet(aux_constr_l,0.0);
  VecSet(aux_constr_u,0.0);
  END_FUNCTION
}

void PetscPDETempRadiation::get_starting_point(Vec& state,
                                  Vec& control,
				                          bool init_z,
				                          Vec& state_lb_mults,
				                          Vec& state_ub_mults,
				                          Vec& control_lb_mults,
				                          Vec& control_ub_mults,
				                          bool init_lambda,
				                          Vec& pde_residual_mults,
				                          Vec& aux_constr_mults)
{
  START_FUNCTION
  state = m_State->vec();
  control = m_Control->vec();
  if(init_z) {
    VecSet(state_lb_mults,0.0);
    VecSet(state_ub_mults,0.0);
    VecSet(control_lb_mults,0.0);
    VecSet(control_ub_mults,0.0);
  }
  if(init_lambda) {
    VecSet(pde_residual_mults,0.0);
    VecSet(aux_constr_mults,0.0);
  }
  END_FUNCTION
}

void PetscPDETempRadiation::calc_objective_part(Number& Val)
{
  START_FUNCTION
  Val = 0.0;
  const unsigned int Dim = m_StateMesh.mesh_dimension();
  FEType fe_type = m_StateDofMap->variable_type(0);

  AutoPtr<FEBase> FEVol(FEBase::build(Dim, fe_type)); // this object will by the current volume object, actualy holding values for current phi and dphi
  AutoPtr<QBase> QuadVol=fe_type.default_quadrature_rule(Dim,0);
  FEVol->attach_quadrature_rule(QuadVol.get());
  const std::vector<Real>& JxWVol = FEVol->get_JxW();   // References to values hold by fe (i.e. they change, ones fe changes)
  const std::vector<std::vector<Number> >& VolPhi = FEVol->get_phi();
  std::vector<unsigned int> DofIdx; // mapping local index -> global index

  AutoPtr<FEBase> FEFace(FEBase::build(Dim, fe_type)); // surface object
  AutoPtr<QBase> QuadFace=fe_type.default_quadrature_rule(Dim-1,0);
  FEFace->attach_quadrature_rule(QuadFace.get());
  const std::vector<Real>& JxWFace = FEFace->get_JxW();
  const std::vector<std::vector<Real> >& FacePhi = FEFace->get_phi();
  std::vector<Number> LocalControl;
  std::vector<Number> LocalState;   // solution values of current element
  std::vector<int> FaceToVol;       // index face node -> index elem node
  double ControlVal(0.0), StateVal(0.0);
  
  MeshBase::const_element_iterator itCurEl = m_StateMesh.active_local_elements_begin();
  const MeshBase::const_element_iterator itEndEl = m_StateMesh.active_local_elements_end();
  for ( ; itCurEl != itEndEl; ++itCurEl) {
    const Elem* CurElem = *itCurEl;
    FEVol->reinit(CurElem);
    m_StateDofMap->dof_indices(CurElem, DofIdx);  // setup local->global mapping in form of array
    LocalState.resize(DofIdx.size());
    for (unsigned int i=0;i<LocalState.size();i++) {
      LocalState[i] = m_State->el(DofIdx[i]);
    }

    for (unsigned short qp=0; qp<QuadVol->n_points(); qp++) {
      StateVal=0.0;
      for(unsigned int iNode=0; iNode<CurElem->n_nodes(); iNode++){
        StateVal += LocalState[iNode] * VolPhi[iNode][qp];
      }
      Val += JxWVol[qp]*m_ReguStateParam*StateVal*StateVal;
    }

    for (unsigned int side=0; side<CurElem->n_sides(); side++) {
      if (CurElem->neighbor(side) == NULL) {
        FEFace->reinit(CurElem, side);
        AutoPtr<libMesh::Elem> Side = CurElem->build_side(side);
        FaceToVol.resize(Side->n_nodes());
        for(unsigned short iNode=0; iNode<Side->n_nodes(); iNode++) {
          for(unsigned short jNode=0; jNode<CurElem->n_nodes(); jNode++) {
            if(Side->node(iNode)==CurElem->node(jNode))
              FaceToVol[iNode] = jNode;
          }
        }
        LocalControl.resize(Side->n_nodes());
        for(unsigned short iNode=0; iNode<Side->n_nodes(); iNode++) {
          LocalControl[iNode] = m_Control->el(m_ControlNodeIDToControlDOF[Side->node(iNode)]);
        }
        for (unsigned short qp=0; qp<QuadFace->n_points(); qp++) {
          ControlVal=0.0;
          for(unsigned int iNode=0; iNode<Side->n_nodes(); iNode++){
            ControlVal += LocalControl[iNode] * FacePhi[FaceToVol[iNode]][qp];
          }
          Val += JxWFace[qp]*ControlVal*ControlVal;
        }
      }
    } // End side loop
  }
  END_FUNCTION
}

void PetscPDETempRadiation::calc_objective_gradient(Vec& grad_state, Vec& grad_control)
{
  START_FUNCTION
  m_GradObjState->zero();
  m_GradObjControl->zero();
  const unsigned int Dim = m_StateMesh.mesh_dimension();
  FEType fe_type = m_StateDofMap->variable_type(0);

  AutoPtr<FEBase> FEVol(FEBase::build(Dim, fe_type)); // this object will by the current volume object, actualy holding values for current phi and dphi
  AutoPtr<QBase> QuadVol=fe_type.default_quadrature_rule(Dim,0);
  FEVol->attach_quadrature_rule(QuadVol.get());
  const std::vector<Real>& JxWVol = FEVol->get_JxW();   // References to values hold by fe (i.e. they change, ones fe changes)
  const std::vector<std::vector<Number> >& VolPhi = FEVol->get_phi();
  std::vector<unsigned int> DofIdx; // mapping local index -> global index

  AutoPtr<FEBase> FEFace(FEBase::build(Dim, fe_type)); // surface object
  AutoPtr<QBase> QuadFace=fe_type.default_quadrature_rule(Dim-1,0);
  FEFace->attach_quadrature_rule(QuadFace.get());
  const std::vector<Real>& JxWFace = FEFace->get_JxW();
  const std::vector<std::vector<Real> >& FacePhi = FEFace->get_phi();
  std::vector<Number> LocalControl, LocalState;
  std::vector<int> FaceToVol;       // index face node -> index elem node
  double ControlVal(0.0), StateVal(0.0);

  MeshBase::const_element_iterator itCurEl = m_StateMesh.active_local_elements_begin();
  const MeshBase::const_element_iterator itEndEl = m_StateMesh.active_local_elements_end();
  for ( ; itCurEl != itEndEl; ++itCurEl) {
    const Elem* CurElem = *itCurEl;
    FEVol->reinit(CurElem);
    m_StateDofMap->dof_indices(CurElem, DofIdx);  // setup local->global mapping in form of array
    LocalState.resize(DofIdx.size());
    for (unsigned int i=0;i<LocalState.size();i++) {
      LocalState[i] = m_State->el(DofIdx[i]);
    }

    for (unsigned short qp=0; qp<QuadVol->n_points(); qp++) {
      StateVal=0.0;
      for(unsigned int iNode=0; iNode<CurElem->n_nodes(); iNode++){
        StateVal += LocalState[iNode] * VolPhi[iNode][qp];
      }
      for(unsigned int iNode=0; iNode<CurElem->n_nodes(); iNode++){
        m_GradObjState->add(DofIdx[iNode],2.0*m_ReguStateParam*JxWVol[qp]*StateVal*VolPhi[iNode][qp]);
      }
    }

    for (unsigned int side=0; side<CurElem->n_sides(); side++) {
      if (CurElem->neighbor(side) == NULL) {
        FEFace->reinit(CurElem, side);
        AutoPtr<libMesh::Elem> Side = CurElem->build_side(side);
        FaceToVol.resize(Side->n_nodes());
        for(unsigned short iNode=0; iNode<Side->n_nodes(); iNode++) {
          for(unsigned short jNode=0; jNode<CurElem->n_nodes(); jNode++) {
            if(Side->node(iNode)==CurElem->node(jNode))
              FaceToVol[iNode] = jNode;
          }
        }
        LocalControl.resize(Side->n_nodes());
        for(unsigned short iNode=0; iNode<Side->n_nodes(); iNode++) {
          LocalControl[iNode] = m_Control->el(m_ControlNodeIDToControlDOF[Side->node(iNode)]);
        }
        for (unsigned short qp=0; qp<QuadFace->n_points(); qp++) {
          ControlVal=0.0;
          for(unsigned short iNode=0; iNode<Side->n_nodes(); iNode++){
            ControlVal += LocalControl[iNode] * FacePhi[FaceToVol[iNode]][qp];
          }
          for (unsigned short iNode=0; iNode<Side->n_nodes(); iNode++) {
            m_GradObjControl->add(m_ControlNodeIDToControlDOF[Side->node(iNode)],
                                  2.0*JxWFace[qp]*ControlVal*FacePhi[FaceToVol[iNode]][qp]);
          }
        }
      }
    }
  }

  m_GradObjState->close();
  m_GradObjControl->close();
  grad_state = m_GradObjState->vec();
  grad_control = m_GradObjControl->vec();
  END_FUNCTION
}

void PetscPDETempRadiation::calcPDE_residual(Vec& residual)
{
  START_FUNCTION
  m_PDEConstr->zero();

  double ControlVal(0.0), StateVal(0.0);
  const unsigned int Dim = m_StateMesh.mesh_dimension();
  FEType fe_type = m_StateDofMap->variable_type(0);
  AutoPtr<FEBase> FEVol(FEBase::build(Dim, fe_type)); // this object will by the current volume object, actualy holding values for current phi and dphi
  AutoPtr<QBase> QuadVol=fe_type.default_quadrature_rule(Dim,0);
  FEVol->attach_quadrature_rule(QuadVol.get());
  const std::vector<Real>& JxWVol = FEVol->get_JxW();   // References to values hold by fe (i.e. they change, ones fe changes)
  const std::vector<std::vector<RealGradient> >& VolDPhi = FEVol->get_dphi();

  AutoPtr<FEBase> FEFace(FEBase::build(Dim, fe_type)); // surface object
  AutoPtr<QBase> QuadFace=fe_type.default_quadrature_rule(Dim-1,0);
  FEFace->attach_quadrature_rule(QuadFace.get());
  const std::vector<Real>& JxWFace = FEFace->get_JxW();
  const std::vector<std::vector<Real> >& FacePhi = FEFace->get_phi();
  
  DenseVector<Number> ElemRes;
  std::vector<unsigned int> DofIdx; // mapping local index -> global index
  std::vector<Number> LocalState;   // solution values of current element
  std::vector<Number> LocalBdState;   // solution values of current element
  std::vector<Number> LocalControl;
  std::vector<int> FaceToVol;       // index face node -> index elem node
  MeshBase::const_element_iterator itCurEl = m_StateMesh.active_local_elements_begin();
  const MeshBase::const_element_iterator itEndEl = m_StateMesh.active_local_elements_end();
  for ( ; itCurEl != itEndEl; ++itCurEl) {
    const Elem* CurElem = *itCurEl;
    m_StateDofMap->dof_indices(CurElem, DofIdx);  // setup local->global mapping in form of array
    FEVol->reinit(CurElem); // reinit fe to fit the current element, i.e. recalulate JxW and dPhi
    ElemRes.resize(DofIdx.size());
    LocalState.resize(DofIdx.size());
    for (unsigned int iCurSol=0;iCurSol<LocalState.size();iCurSol++) {
      LocalState[iCurSol] = m_State->el(DofIdx[iCurSol]);
    }
    unsigned int qp,i,j;
    for (qp=0; qp<QuadVol->n_points(); qp++) {
      for (i=0; i<VolDPhi.size(); i++) {
        for (j=0; j<VolDPhi.size(); j++) {
          ElemRes(i) += JxWVol[qp]*(VolDPhi[i][qp]*VolDPhi[j][qp])*LocalState[j];
        }
      }
    }

    for (unsigned int side=0; side<CurElem->n_sides(); side++) {
      if (CurElem->neighbor(side) == NULL) {
        FEFace->reinit(CurElem, side);
        AutoPtr<libMesh::Elem> Side = CurElem->build_side(side);
        FaceToVol.resize(Side->n_nodes());
        for(unsigned short iNode=0; iNode<Side->n_nodes(); iNode++) {
          for(unsigned short jNode=0; jNode<CurElem->n_nodes(); jNode++) {
            if(Side->node(iNode)==CurElem->node(jNode))
              FaceToVol[iNode] = jNode;
          }
        }
        LocalControl.resize(Side->n_nodes());
        LocalState.resize(Side->n_nodes());
        for(unsigned short iNode=0; iNode<Side->n_nodes(); iNode++) {
          LocalControl[iNode] = m_Control->el(m_ControlNodeIDToControlDOF[Side->node(iNode)]);
          LocalState[iNode] = m_State->el(DofIdx[FaceToVol[iNode]]);
        }
        unsigned int qp,i;
        for (qp=0; qp<QuadFace->n_points(); qp++) {
          ControlVal=0.0;
          StateVal=0.0;
          for(unsigned int iNode=0; iNode<Side->n_nodes(); iNode++){
            ControlVal += LocalControl[iNode] * FacePhi[FaceToVol[iNode]][qp];
            StateVal += LocalState[iNode]*FacePhi[FaceToVol[iNode]][qp];
          }
          for (i=0; i<Side->n_nodes(); i++) {
            ElemRes(FaceToVol[i]) -= JxWFace[qp]*m_Alpha* (pow(ControlVal,CONTROL_EXPONENT) - pow(StateVal,4.0))*FacePhi[FaceToVol[i]][qp];
          }
        }
      }
    }
    m_StateDofMap->constrain_element_vector(ElemRes, DofIdx); // Add constraints for hanging nodes (refinement)
    ElemRes *= m_PDEConstrScale;
    m_PDEConstr->add_vector(ElemRes, DofIdx);
  }
  m_PDEConstr->close();
  residual = m_PDEConstr->vec();
  END_FUNCTION
}

void PetscPDETempRadiation::calcPDE_jacobians(Mat& jac_state, Mat& jac_control)
{
  START_FUNCTION
  m_JacPDEState->zero();
  m_JacPDEControl->zero();

  double ControlVal, StateVal;
  const unsigned int Dim = m_StateMesh.mesh_dimension();
  FEType fe_type = m_StateDofMap->variable_type(0);

  AutoPtr<FEBase> FEVol(FEBase::build(Dim, fe_type)); // this object will by the current volume object, actualy holding values for current phi and dphi
  AutoPtr<QBase> QuadVol=fe_type.default_quadrature_rule(Dim,0);
  FEVol->attach_quadrature_rule(QuadVol.get());
  const std::vector<Real>& JxWVol = FEVol->get_JxW();   // References to values hold by fe (i.e. they change, ones fe changes)
  const std::vector<std::vector<RealGradient> >& VolDPhi = FEVol->get_dphi();

  AutoPtr<FEBase> FEFace(FEBase::build(Dim, fe_type)); // surface object
  AutoPtr<QBase> QuadFace=fe_type.default_quadrature_rule(Dim-1,0);
  FEFace->attach_quadrature_rule(QuadFace.get());
  const std::vector<Real>& JxWFace = FEFace->get_JxW();
  const std::vector<std::vector<Real> >& FacePhi = FEFace->get_phi();

  DenseMatrix<Number> ElemJacPDEState;
  std::vector<unsigned int> DofIdx; // mapping local index -> global index
  std::vector<Number> LocalState;   // solution values of current element
  std::vector<Number> LocalControl; // control values of current element
  std::vector<int> FaceToVol;       // index face node -> index elem node
  MeshBase::const_element_iterator itCurEl = m_StateMesh.active_local_elements_begin();
  const MeshBase::const_element_iterator itEndEl = m_StateMesh.active_local_elements_end();
  for ( ; itCurEl != itEndEl; ++itCurEl) {
    const Elem* CurElem = *itCurEl;
    m_StateDofMap->dof_indices(CurElem, DofIdx);  // setup local->global mapping in form of array
    FEVol->reinit(CurElem); // reinit fe to fit the current element, i.e. recalulate JxW and dPhi
    ElemJacPDEState.resize(DofIdx.size(),DofIdx.size());
    LocalState.resize(DofIdx.size());
    for (unsigned int iCurSol=0;iCurSol<LocalState.size();iCurSol++) {
      LocalState[iCurSol] = m_State->el(DofIdx[iCurSol]);
    }
    unsigned int qp,i,j;
    for (qp=0; qp<QuadVol->n_points(); qp++) {
      for (i=0; i<VolDPhi.size(); i++) {
        for (j=0; j<VolDPhi.size(); j++) {
          if(m_CalcType==Values) {
            ElemJacPDEState(i,j) += JxWVol[qp]*(VolDPhi[i][qp]*VolDPhi[j][qp]);
          }
          else {
            ElemJacPDEState(i,j) += 1.0;
          }
        }
      }
    }

    for (unsigned int side=0; side<CurElem->n_sides(); side++) {
      if (CurElem->neighbor(side) == NULL) {
        FEFace->reinit(CurElem, side);
        AutoPtr<libMesh::Elem> Side = CurElem->build_side(side);
        FaceToVol.resize(Side->n_nodes());
        for(unsigned short iNode=0; iNode<Side->n_nodes(); iNode++) {
          for(unsigned short jNode=0; jNode<CurElem->n_nodes(); jNode++) {
            if(Side->node(iNode)==CurElem->node(jNode))
              FaceToVol[iNode] = jNode;
          }
        }
        LocalControl.resize(Side->n_nodes());
        LocalState.resize(Side->n_nodes());
        for(unsigned short iNode=0; iNode<Side->n_nodes(); iNode++) {
          LocalControl[iNode] = m_Control->el(m_ControlNodeIDToControlDOF[Side->node(iNode)]);
          LocalState[iNode] = m_State->el(DofIdx[FaceToVol[iNode]]);
        }
        unsigned int qp,i;
        for (qp=0; qp<QuadFace->n_points(); qp++) {
          ControlVal=0.0;
          StateVal=0.0;
          for(unsigned int iNode=0; iNode<Side->n_nodes(); iNode++){
            ControlVal += LocalControl[iNode] * FacePhi[FaceToVol[iNode]][qp];
            StateVal += LocalState[iNode]*FacePhi[FaceToVol[iNode]][qp];
          }
          for (i=0; i<Side->n_nodes(); i++) {
            for (j=0; j<Side->n_nodes(); j++) {
              if(m_CalcType==Values) {
                ElemJacPDEState(FaceToVol[i],FaceToVol[j]) += JxWFace[qp]*m_Alpha*4.0*pow(StateVal,3.0)*FacePhi[FaceToVol[j]][qp]*FacePhi[FaceToVol[i]][qp];
                m_JacPDEControl->add(DofIdx[FaceToVol[i]],m_ControlNodeIDToControlDOF[Side->node(j)],-m_PDEConstrScale*JxWFace[qp]*m_Alpha*CONTROL_EXPONENT*pow(ControlVal,CONTROL_EXPONENT-1.0)*FacePhi[FaceToVol[j]][qp]*FacePhi[FaceToVol[i]][qp]);
              }
              else {
                ElemJacPDEState(FaceToVol[i],FaceToVol[j]) += 1.0;
                m_JacPDEControl->add(DofIdx[FaceToVol[i]],m_ControlNodeIDToControlDOF[Side->node(j)],1.0);
              }
            }
          }
        }
      }
    }
    m_StateDofMap->constrain_element_matrix(ElemJacPDEState, DofIdx); // Add constraints for hanging nodes (refinement)
    ElemJacPDEState *= m_PDEConstrScale;
    m_JacPDEState->add_matrix(ElemJacPDEState, DofIdx);
  }

  m_JacPDEState->close();
  m_JacPDEControl->close();
  jac_state = m_JacPDEState->mat();
  jac_control = m_JacPDEControl->mat();
  END_FUNCTION
}

void PetscPDETempRadiation::calc_hessians(Number sigma, Vec& lambda_loc_pde, Vec& lambda_loc_aux,
                                                 Mat& Hcc, Mat& Hcs, Mat& Hss)
{
  START_FUNCTION
  m_HessControlControl->zero();
  m_HessControlState->zero();
  m_HessStateState->zero();
  
  double ControlVal(0.0);
  const unsigned int Dim = m_StateMesh.mesh_dimension();
  FEType fe_type = m_StateDofMap->variable_type(0);

  AutoPtr<FEBase> FEVol(FEBase::build(Dim, fe_type)); // this object will by the current volume object, actualy holding values for current phi and dphi
  AutoPtr<QBase> QuadVol=fe_type.default_quadrature_rule(Dim,0);
  FEVol->attach_quadrature_rule(QuadVol.get());
  const std::vector<Number>& JxWVol = FEVol->get_JxW();   // References to values hold by fe (i.e. they change, ones fe changes)
  const std::vector<std::vector<Number> >& VolPhi = FEVol->get_phi();
  const std::vector<Point>& xyzVol = FEVol->get_xyz();   // References to values hold by fe (i.e. they change, ones fe changes)

  AutoPtr<FEBase> FEFace(FEBase::build(Dim, fe_type)); // surface object
  AutoPtr<QBase> QuadFace=fe_type.default_quadrature_rule(Dim-1,0);
  FEFace->attach_quadrature_rule(QuadFace.get());
  const std::vector<Real>& JxWFace = FEFace->get_JxW();
  const std::vector<std::vector<Real> >& FacePhi = FEFace->get_phi();

  DenseMatrix<Number> ElemHessStateState;
  std::vector<unsigned int> StateDofIdx; // mapping local index -> global index
  std::vector<Number> LocalState;   // solution values of current element
  std::vector<Number> LocalControl; // control values of current element
  std::vector<int> FaceToVol;       // index face node -> index elem node
  MeshBase::const_element_iterator itCurEl = m_StateMesh.active_local_elements_begin();
  const MeshBase::const_element_iterator itEndEl = m_StateMesh.active_local_elements_end();

  PetscVector<Number> lambda_pde(lambda_loc_pde);

  for ( ; itCurEl != itEndEl; ++itCurEl) {
    const Elem* CurElem = *itCurEl;
    FEVol->reinit(CurElem); // reinit fe to fit the current element, i.e. recalulate JxW and dPhi
    m_StateDofMap->dof_indices(CurElem, StateDofIdx);  // setup local->global mapping in form of array
    ElemHessStateState.resize(StateDofIdx.size(),StateDofIdx.size());

    for (unsigned short qp=0; qp<QuadVol->n_points(); qp++) {
      for (unsigned short i=0; i<CurElem->n_nodes(); i++) {
        for (unsigned short j=0; j<=i; j++) {
          ElemHessStateState(i,j) += JxWVol[qp]*sigma*m_ReguStateParam*VolPhi[i][qp]*VolPhi[j][qp];
        }
      }
    }

    unsigned int qp,i,j;
    for (unsigned int side=0; side<CurElem->n_sides(); side++) {
      if (CurElem->neighbor(side) == NULL) {
        FEFace->reinit(CurElem, side);
        AutoPtr<libMesh::Elem> Side = CurElem->build_side(side);
        FaceToVol.resize(Side->n_nodes());
        for(unsigned short iNode=0; iNode<Side->n_nodes(); iNode++) {
          for(unsigned short jNode=0; jNode<CurElem->n_nodes(); jNode++) {
            if(Side->node(iNode)==CurElem->node(jNode))
              FaceToVol[iNode] = jNode;
          }
        }
        LocalControl.resize(Side->n_nodes());
        LocalState.resize(Side->n_nodes());
        for(unsigned short iNode=0; iNode<Side->n_nodes(); iNode++) {
          LocalControl[iNode] = m_Control->el(m_ControlNodeIDToControlDOF[Side->node(iNode)]);
          LocalState[iNode] = m_State->el(StateDofIdx[FaceToVol[iNode]]);
        }
        unsigned int qp,i,j,k;
        double StateVal;
        for (qp=0; qp<QuadFace->n_points(); qp++) {
          ControlVal=0.0;
          StateVal=0.0;
          for(unsigned int iNode=0; iNode<Side->n_nodes(); iNode++){
            ControlVal += LocalControl[iNode] * FacePhi[FaceToVol[iNode]][qp];
            StateVal += LocalState[iNode]*FacePhi[FaceToVol[iNode]][qp];
          }
          for (i=0; i<Side->n_nodes(); i++) {
            for (j=0; j<=i; j++) {
              m_HessControlControl->add(m_ControlNodeIDToControlDOF[Side->node(i)],m_ControlNodeIDToControlDOF[Side->node(j)],JxWFace[qp]*2.0*sigma*FacePhi[FaceToVol[j]][qp]*FacePhi[FaceToVol[i]][qp]);
              for (k=0; k<Side->n_nodes(); k++) {
                if(m_CalcType==Values) {
                  ElemHessStateState(FaceToVol[i],FaceToVol[j]) += lambda_pde.el(StateDofIdx[FaceToVol[k]])*JxWFace[qp]*m_Alpha*12.0*pow(StateVal  ,2.0)*FacePhi[FaceToVol[j]][qp]*FacePhi[FaceToVol[i]][qp]*FacePhi[FaceToVol[k]][qp];
                  if(CONTROL_EXPONENT>1) {
                    m_HessControlControl->add(m_ControlNodeIDToControlDOF[Side->node(i)],m_ControlNodeIDToControlDOF[Side->node(j)],-lambda_pde.el(StateDofIdx[FaceToVol[k]])*JxWFace[qp]*m_Alpha*CONTROL_EXPONENT*(CONTROL_EXPONENT-1)*pow(ControlVal,CONTROL_EXPONENT-2.0)*FacePhi[FaceToVol[j]][qp]*FacePhi[FaceToVol[i]][qp]*FacePhi[FaceToVol[k]][qp]*m_PDEConstrScale);
                  }
                }
                else {
                  ElemHessStateState(FaceToVol[i],FaceToVol[j]) += 1.0;
                  m_HessControlControl->add(m_ControlNodeIDToControlDOF[Side->node(i)],m_ControlNodeIDToControlDOF[Side->node(j)], 1.0);
                }
              }
            }
          }
        }
      }
    }
    m_StateDofMap->constrain_element_matrix(ElemHessStateState, StateDofIdx); // Add constraints for hanging nodes (refinement)
    
    ElemHessStateState *= m_PDEConstrScale;
    m_HessStateState->add_matrix(ElemHessStateState, StateDofIdx);
  }

  m_HessControlControl->close();
  m_HessControlState->close();
  m_HessStateState->close();
  Hcc = m_HessControlControl->mat();
  Hcs = m_HessControlState->mat();
  Hss = m_HessStateState->mat();
  END_FUNCTION
}

void PetscPDETempRadiation::calcAux_constr(Vec& constr)
{
  START_FUNCTION
  m_AuxConstr->zero();
  m_AuxConstr->close();
  constr = m_AuxConstr->vec();
  END_FUNCTION
}

void PetscPDETempRadiation::calcAux_jacobians(Mat& jac_state,
                                 Mat& jac_control)
{
  START_FUNCTION
  m_JacAuxState->zero();
  m_JacAuxState->close();
  jac_state = m_JacAuxState->mat();
  m_JacAuxControl->zero();
  m_JacAuxControl->close();
  jac_control = m_JacAuxControl->mat();
  END_FUNCTION
}

void PetscPDETempRadiation::Write2File( const std::string& pre_filename)
{
  START_FUNCTION
  std::vector<Number> State;
  m_State->close();
  m_State->localize_to_one(State);
  std::string filename(pre_filename);
  filename.append("State.dat");
  std::ofstream f;
  f.open(filename.c_str(),std::ios::out);
  f.precision(10);
  const unsigned int dim = m_StateMesh.mesh_dimension();
  MeshBase::node_iterator it_cur_node = m_StateMesh.active_nodes_begin();
  const MeshBase::node_iterator end_node = m_StateMesh.active_nodes_end();

  for (;it_cur_node!=end_node;it_cur_node++) {
    const Node* cur_node = *it_cur_node;
    for (unsigned int i=0;i<dim;i++)
      f << (*cur_node)(i) << ", ";
    f << State[cur_node->dof_number(0,0,0)] << std::endl;
  }
  f.close();
  END_FUNCTION
}

void PetscPDETempRadiation::set_finalize_vectors(Vec& lm_state_lb_mults,
				                            Vec& lm_state_ub_mults,
				                            Vec& lm_control_lb_mults,
				                            Vec& lm_control_ub_mults,
				                            Vec& lm_pde_residual_mults,
				                            Vec& lm_aux_constr_mults)
{
  START_FUNCTION

  ExodusII_IO exoio(m_StateMesh);
  std::vector<double> Var;
  m_State->localize_to_one(Var);
  std::vector< std::string > VarNames;
  VarNames.push_back("State");
  exoio.write_nodal_data (std::string("State.ex2"), Var, VarNames);

  ExodusII_IO exoio2(m_StateMesh);
  MeshBase::node_iterator it_cur_node = m_StateMesh.active_nodes_begin();
  const MeshBase::node_iterator end_node = m_StateMesh.active_nodes_end();
  for (;it_cur_node!=end_node;it_cur_node++) {
    const Node* cur_node = *it_cur_node;
    if(IsOnBoundary(*cur_node)) {
      Var[cur_node->dof_number(0,0,0)] = m_Control->el(m_ControlNodeIDToControlDOF[cur_node->id()]);
    }
    else {
      Var[cur_node->dof_number(0,0,0)] = 0;
    }
  }

  VarNames[0] = "Control";
  exoio2.write_nodal_data (std::string("Control.ex2"), Var, VarNames);

  END_FUNCTION
}

bool PetscPDETempRadiation::IsOnBoundary(const libMesh::Point& pt)
{
  START_FUNCTION
  double eps = 1e-8;
  if( (fabs(pt(0))>eps) && (fabs(pt(0)-1)>eps) &&
      (fabs(pt(1))>eps) && (fabs(pt(1)-1)>eps) &&
      (fabs(pt(2))>eps) && (fabs(pt(2)-1)>eps) ) 
    return false;
  else
    return true;
  END_FUNCTION
}

bool PetscPDETempRadiation::IsInOmega1(const Point& pt)
{
  START_FUNCTION
  // is in the rectangle defined by OmegaInner?
  if(pt(0)<m_Omega1[0])
    return false;
  if(pt(1)<m_Omega1[1])
    return false;
  if( (m_Dim==3) && (pt(2)<m_Omega1[2]) )
    return false;
  if(pt(0)>m_Omega1[3])
    return false;
  if(pt(1)>m_Omega1[4])
    return false;
  if( (m_Dim==3) && (pt(2)>m_Omega1[5]) )
    return false;
  return true;
  END_FUNCTION
}

bool PetscPDETempRadiation::IsInOmega2(const Point& pt)
{
  START_FUNCTION
  if(pt(0)<m_Omega2[0])
    return false;
  if(pt(1)<m_Omega2[1])
    return false;
  if( (m_Dim==3) && (pt(2)<m_Omega2[2]) )
    return false;
  if(pt(0)>m_Omega2[3])
    return false;
  if(pt(1)>m_Omega2[4])
    return false;
  if( (m_Dim==3) && (pt(2)>m_Omega2[5]) )
    return false;
  return true;
  END_FUNCTION
}

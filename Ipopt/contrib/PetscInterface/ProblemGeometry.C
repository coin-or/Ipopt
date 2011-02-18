// Copyright (C) 2010 International Business Machines.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// $Id$
//
// Authors:  Johannes Huber     IBM        2010-09-03

#include "ProblemGeometry.h"
#include <stdexcept>
#include "mesh_generation.h"


extern int GetProcID();

//#define MY_DBG_PRINT(s) {std::cout << GetProcID() << __FILE__ << ":" << __LINE__ <<":" << s << std::endl;}
#define MY_DBG_PRINT(s) {}


///////////////////////////////////////////////////////////////////////////////////
//  ProblemGeometry::Item
ProblemGeometry::Item::Item(const std::vector<double>& min, const std::vector<double>& max)
{
  assert(min.size()==max.size());
  min_=min;
  max_=max;
  ValidateMinMax();
  BoundaryMarker.resize(2*min_.size(),-1); // those are for floor and ceiling
}

void ProblemGeometry::Item::ValidateMinMax()
{
  for (int i=0;i<min_.size();i++) {
    if (min_[i]>max_[i]) {
      double tmp=min_[i];
      min_[i]=max_[i];
      max_[i]=tmp;
    }
  }
}

bool ProblemGeometry::Item::IsInArea(const std::vector<double>& pt, double epsilon/*=1e-8*/)
{
  assert(pt.size()==min_.size());
  int dim = min_.size();
  for ( int i=0; i<dim; i++ ) {
    if (min_[i]-epsilon>pt[i])
      return false;
    if (max_[i]+epsilon<pt[i])
      return false;
  }
  return true;
}

int ProblemGeometry::Item::GetClosestBoundary(const std::vector<double>& pt)
{
  assert(pt.size()==min_.size());
  int dim = min_.size();
  double min_dist=fabs(pt[0]-min_[0]);
  int cur_min_dir=0;
  for (int i=0; i<dim; i++) {
    if (fabs(pt[i]-min_[i])<min_dist) {
      min_dist = fabs(pt[i]-min_[i]);
      cur_min_dir = 2*i;
    }
    if (fabs(max_[i]-pt[i])<min_dist) {
      min_dist = fabs(max_[i]-pt[i]);
      cur_min_dir = 2*i+1;
    }
  }
  return cur_min_dir;
}

int ProblemGeometry::Item::GetBoundaryMarker(const std::vector<double>& pt)
{
  assert(pt.size()==min_.size());
  if (!IsInArea(pt)) return -1;
  return BoundaryMarker[GetClosestBoundary(pt)];
}

//////////////////////////////////////////////////////////////////////////////////////
// ProblemGeometry
ProblemGeometry::ProblemGeometry()
{
  _BoundCond.push_back(new BoundaryConditionConstValues(0.0,1.0,0.0,0.0,1.0,0.0)); //   Add one Dummy to start from 1 (0:Interior Point)
  _BoundCond.push_back(new BoundaryConditionConstValues(0.0,1.0,0.0,0.0,1.0,0.0)); //   homogeneous boundary conditions on wall
  NextFreeBoundaryMarker = 2;
  _h = 0.0;
}

ProblemGeometry::~ProblemGeometry()
{
  for(int iBC=0;iBC<_BoundCond.size();iBC++) {
    delete _BoundCond[iBC];
    _BoundCond[iBC] = NULL;
  }
}

void ProblemGeometry::AddEquipment(const std::vector<double>& min, const std::vector<double>& max, double TempEquip, double kappa)
{
  assert(min.size()==GetDim());
  assert(max.size()==GetDim());
  Item item(min,max);
  assert(_BoundCond.size()==NextFreeBoundaryMarker);
  int BoundMark = NextFreeBoundaryMarker;
  // 1.st & 3rd Boundary: Head exchange
  item.BoundaryMarker[2] = BoundMark;
  item.BoundaryMarker[3] = BoundMark;
  // Phi: hom. Neum. T: dT/dn = Kappa(T-Teq) => 1*dT/dn = Kappa*T -Kappa*Teq
  _BoundCond.push_back(new BoundaryConditionConstValues(0.0,1.0,0.0,kappa,-1.0,-kappa*TempEquip));
  // No control parameters
  // 2.nd Boundary: Isolation
  item.BoundaryMarker[0] = BoundMark+1;
  item.BoundaryMarker[1] = BoundMark+1;
  // Phi: hom. Neum. T: hom. Neum
  _BoundCond.push_back(new BoundaryConditionConstValues(0.0,1.0,0.0,0.0,1.0,0.0));
  // No control parameters
  NextFreeBoundaryMarker+=2;
  if (GetDim()>2) { // Equipmentceiling: heat exchange (same as on side)
    item.BoundaryMarker[4] = -1;   // should never be accessed, since Equipment is standing on floor, so no boundary
    int BoundMark = NextFreeBoundaryMarker++;
    // 1.st & 3rd Boundary: Head exchange
    item.BoundaryMarker[5] = BoundMark;
    // Phi: hom. Neum. T: dT/dn = Kappa(T-Teq)
    _BoundCond.push_back(new BoundaryConditionConstValues(0.0,1.0,0.0,kappa,1.0,kappa*TempEquip));
  }
  _Equip.push_back(item);
}

void ProblemGeometry::GetHeatExchangeBoundaryMarkers(int iEquip, std::set<int>* pVals)
{
  assert(pVals);
  pVals->clear();
  // Only one boundary marker for both walls for heatr exchange (see AddEquipment)
  pVals->insert(_Equip[iEquip].BoundaryMarker[2]);
  pVals->insert(_Equip[iEquip].BoundaryMarker[3]);
}

void ProblemGeometry::AddAC(const std::vector<double>& min, const std::vector<double>& max, double vAc, double TempAc)
{
  assert(min.size()==GetDim());
  assert(max.size()==GetDim());
  _ItemAtWall.clear();  // invalidate item at wall and rebuild it, when the mesh is created
  Item item(min,max);
  int BoundMark = NextFreeBoundaryMarker++;
  for (unsigned int iBd=0;iBd<2*GetDim();iBd++)
    item.BoundaryMarker[iBd] = BoundMark;
  assert(_BoundCond.size()==BoundMark);
  _BoundCond.push_back(new BoundaryConditionSquarePhiRhs(0.0,-1.0,vAc,1.0,0.0,-TempAc,min,max)); // dPhi/dn = -vAc,T = TAc, i.e. 0 dT/dn = 1*T-TAc
  _ParamIdx2BCParam.push_back(ControlParameter(BoundMark,2));
  //_ParamIdx2BCParam.push_back(ControlParameter(BoundMark,5));
  _AC.push_back(item);
}

#ifndef EXHAUST_AS_CONTROL
void ProblemGeometry::AddExhaust(const std::vector<double>& min, const std::vector<double>& max)
{
  assert(min.size()==GetDim());
  assert(max.size()==GetDim());
  _ItemAtWall.clear();
  Item item(min,max);
  unsigned int BoundMark = NextFreeBoundaryMarker++;
  for ( int iBd=0;iBd<2*GetDim();iBd++)
    item.BoundaryMarker[iBd] = BoundMark;
  assert(_BoundCond.size()==BoundMark);
  _BoundCond.push_back(new BoundaryConditionConstValues(1.0,0.0,0.0,0.0,1.0,0.0));         // Phi = 0, dT/dn = 0
  //std::cout << "Adding BoundaryConditionConstValues at " << _BoundCond[_BoundCond.size()-1] << std::endl;
  //_BoundCond[_BoundCond.size()-1]->print();
  // No Control Parameter
  _Exh.push_back(item);
}
#else
void ProblemGeometry::AddExhaust(const std::vector<double>& min, const std::vector<double>& max)
{
  assert(min.size()==GetDim());
  assert(max.size()==GetDim());
  _ItemAtWall.clear();
  Item item(min,max);
  unsigned int BoundMark = NextFreeBoundaryMarker++;
  for ( int iBd=0;iBd<2*GetDim();iBd++)
    item.BoundaryMarker[iBd] = BoundMark;
  assert(_BoundCond.size()==BoundMark);
  // TODO: Retried vExh from input
  double vExh = 1.0;
  double TempExh = 1e30;
  _BoundCond.push_back(new BoundaryConditionConstValues(0.0,+1.0,vExh,1.0,0.0,-TempExh)); // dPhi/dn = -vAc,T = TAc, i.e. 0 dT/dn = 1*T-TAc
  _ParamIdx2BCParam.push_back(ControlParameter(BoundMark,2));
  //_ParamIdx2BCParam.push_back(ControlParameter(BoundMark,5));
  _Exh.push_back(item);
}
#endif

void ProblemGeometry::ReadFromStream(std::istream& is)
{
  static char Buf[1024];
  int n;
  double Vals[8];
  _h=0;
  while (!is.eof()) {
    is >> Buf;
    if (strlen(Buf)<1)
      continue;
    if (Buf[0]=='#')
      continue;
    n = sscanf(Buf,"h=%lf",Vals);
    if ( (n==1) ) {
      _h = Vals[0];
      break;
    }
    else {
      std::cout << "WARNING while reading input file and looking for \"h=*\": Line in wrong position?" << std::endl;
      std::cout << Buf << std::endl;
      continue;
    }
  }

  _RoomSize.clear();
  while (!is.eof()) {
    is >> Buf;
    if (strlen(Buf)<1)
      continue;
    if (Buf[0]=='#')
      continue;
    n = sscanf(Buf,"Roomsize=%lf,%lf,%lf",Vals,Vals+1,Vals+2);
    if ( (n==2) || (n==3) ) {
      for (int i=0;i<n;i++)
        _RoomSize.push_back(Vals[i]);
      break;
    }
    else {
      std::cout << "WARNING while reading input file and looking for \"Roomsize=*,*[,*]\": Line in wrong position?" << std::endl;
      std::cout << Buf << std::endl;
      continue;
    }
  }
  if (0==_RoomSize.size()) {
    throw std::runtime_error("Can't find line \"Roomsize=*,*[,*]\"");
  }
  else {
    if(GetProcID()==0)
      std::cout << "Roomsize=" << _RoomSize[0] <<"," << _RoomSize[1];
    if (_RoomSize.size()>2)
      std::cout << "," << _RoomSize[2];
    std::cout << std::endl;
  }

  _Equip.clear();
  _AC.clear();
  _Exh.clear();
  std::vector<double> tmp_min;
  std::vector<double> tmp_max;
  while (!is.eof()) {
    tmp_max.clear();
    tmp_min.clear();
    is >> Buf;
    if(GetProcID()==0)
      std::cout << Buf << "line read" << std::endl;
    if (strlen(Buf)<1)
      continue;
    if (Buf[0]=='#')
      continue;

    n = sscanf(Buf,"AC:Min=%lf,%lf,%lf;",Vals,Vals+1,Vals+2);
    if (n>1) {
      if (n==2) {
        n = sscanf(Buf,"AC:Min=%lf,%lf;Max=%lf,%lf;vAC=%lf;TAC=%lf",Vals,Vals+1,Vals+2,Vals+3,Vals+4,Vals+5);
        if (n!=6) {
          std::string str("Can't interpret line as AC in 2D:");
          str += Buf;
          throw std::runtime_error(str);
        }
        tmp_min.push_back(Vals[0]);
        tmp_min.push_back(Vals[1]);
        tmp_max.push_back(Vals[2]);
        tmp_max.push_back(Vals[3]);
        AddAC(tmp_min, tmp_max, Vals[4], Vals[5]);
      }
      else {
        n = sscanf(Buf,"AC:Min=%lf,%lf,%lf;Max=%lf,%lf,%lf;vAC=%lf;TAC=%lf",Vals,Vals+1,Vals+2,Vals+3,Vals+4,Vals+5,Vals+6,Vals+7);
        if (n!=8) {
          std::string str("Can't interpret line as AC in 3D:");
          str += Buf;
          throw std::runtime_error(str);
        }
        tmp_min.push_back(Vals[0]);
        tmp_min.push_back(Vals[1]);
        tmp_min.push_back(Vals[2]);
        tmp_max.push_back(Vals[3]);
        tmp_max.push_back(Vals[4]);
        tmp_max.push_back(Vals[5]);
        AddAC(tmp_min, tmp_max, Vals[6], Vals[7]);
      }
    }
    else {
      n = sscanf(Buf,"Exh:Min=%lf,%lf,%lf;",Vals,Vals+1,Vals+2);
      if (n>1) {
        if (n==2) {
          n = sscanf(Buf,"Exh:Min=%lf,%lf;Max=%lf,%lf",Vals,Vals+1,Vals+2,Vals+3);
          if (n!=4) {
            std::string str("Can't interpret line as Exh in 2D:");
            str += Buf;
            throw std::runtime_error(str);
          }
          tmp_min.push_back(Vals[0]);
          tmp_min.push_back(Vals[1]);
          tmp_max.push_back(Vals[2]);
          tmp_max.push_back(Vals[3]);
          AddExhaust(tmp_min, tmp_max);
        }
        else {
          n = sscanf(Buf,"Exh:Min=%lf,%lf,%lf;Max=%lf,%lf,%lf",Vals,Vals+1,Vals+2,Vals+3,Vals+4,Vals+5);
          if (n!=6) {
            std::string str("Can't interpret line as Exh in 3D:");
            str += Buf;
            throw std::runtime_error(str);
          }
          tmp_min.push_back(Vals[0]);
          tmp_min.push_back(Vals[1]);
          tmp_min.push_back(Vals[2]);
          tmp_max.push_back(Vals[3]);
          tmp_max.push_back(Vals[4]);
          tmp_max.push_back(Vals[5]);
          AddExhaust(tmp_min, tmp_max);
        }
      }
      else {
        n = sscanf(Buf,"Equip:Min=%lf,%lf,%lf;",Vals,Vals+1,Vals+2);
        if (n>1) {
          if (n==2) {
            n = sscanf(Buf,"Equip:Min=%lf,%lf;Max=%lf,%lf;TEquip=%lf;kappa=%lf",Vals,Vals+1,Vals+2,Vals+3,Vals+4,Vals+5);
            if (n!=6) {
              std::string str("Can't interpret line as Equip in 2D:");
              str += Buf;
              throw std::runtime_error(str);
            }
            tmp_min.push_back(Vals[0]);
            tmp_min.push_back(Vals[1]);
            tmp_max.push_back(Vals[2]);
            tmp_max.push_back(Vals[3]);
            AddEquipment(tmp_min, tmp_max, Vals[4], Vals[5]);
          }
          else {
            n = sscanf(Buf,"Equip:Min=%lf,%lf,%lf;Max=%lf,%lf,%lf;TEquip=%lf;kappa=%lf",Vals,Vals+1,Vals+2,Vals+3,Vals+4,Vals+5,Vals+6,Vals+7);
            if (n!=8) {
              std::string str("Can't interpret line as Equip in 3D:");
              str += Buf;
              throw std::runtime_error(str);
            }
            tmp_min.push_back(Vals[0]);
            tmp_min.push_back(Vals[1]);
            tmp_min.push_back(Vals[2]);
            tmp_max.push_back(Vals[3]);
            tmp_max.push_back(Vals[4]);
            tmp_max.push_back(Vals[5]);
            AddEquipment(tmp_min, tmp_max, Vals[6], Vals[7]);
          }
        }
        else {
          std::string str("Can't interpret line:");
          str += Buf;
          throw std::runtime_error(str);
        }
      }
    }
  }
  if(GetProcID()==0)
    std::cout << "Input file read" << std::endl;
}

int ProblemGeometry::GetBoundaryMarker(const std::vector<double>& pt)
{
  int iWall = GetWall(pt);
  int Tmp;
  if (iWall<0)  {
    for (int iEq=0;iEq<_Equip.size();iEq++) {
      Tmp=_Equip[iEq].GetBoundaryMarker(pt);
      if (Tmp>=0)
        return Tmp;
    }
    return -1;
  }
  else {
    if (iWall>3) // floor and ceiling
      return 1;
    ValidateItemAtWall();
    std::list< Item* >::iterator ppItem;
    ppItem = _ItemAtWall[iWall].begin();
    for (;ppItem!=_ItemAtWall[iWall].end();ppItem++) {
      Tmp=(*ppItem)->GetBoundaryMarker(pt);
      if (Tmp>=0)
        return Tmp;
    }
    return 1;
  }
}

int ProblemGeometry::GetWall(const std::vector<double>& pt)
{
  assert(pt.size()==GetDim());
  std::vector<double> Null;
  Null.resize(GetDim(),0.0);

  Item room(Null,_RoomSize);
  if (room.IsInArea(pt,-1e-8)) // Check if it's inside
    return -1;
  if (room.IsInArea(pt))
    return room.GetClosestBoundary(pt);
  else
    return -1;
}

// build up _TimeAtWall member for each 4 walls
void ProblemGeometry::ValidateItemAtWall()
{
  // _ItemAtWall: Index: 0:x0=0, 1:x0=Max, 2:x1=0, 3:x1=Max
  double epsilon = 1e-8;
  if (_ItemAtWall.size()>0)
    return;
  std::list< Item* > EmptyList;
  _ItemAtWall.resize(4,EmptyList);
  Item* pItem;
  unsigned int iItem;
  for (iItem=0;iItem<_AC.size();iItem++) {
    pItem = &(_AC[iItem]);
    for (int iDim=0;iDim<2;iDim++) // No check for iDim=2, thats floor or ceiling, special case
    {
      if ( fabs(pItem->min_[iDim])<epsilon && fabs(pItem->max_[iDim])<epsilon )
      {
        _ItemAtWall[2*iDim].push_back(pItem);
        continue;
      }
      if ( fabs(pItem->min_[iDim]-_RoomSize[iDim])<epsilon && fabs(pItem->max_[iDim]-_RoomSize[iDim])<epsilon )
      {
        _ItemAtWall[2*iDim+1].push_back(pItem);
        continue;
      }
    }
  }

  for (iItem=0;iItem<_Exh.size();iItem++) {
    pItem = &(_Exh[iItem]);
    for (int iDim=0;iDim<2;iDim++) // No check for iDim=2, thats floor or ceiling, special case
    {
      if ( fabs(pItem->min_[iDim])<epsilon && fabs(pItem->max_[iDim])<epsilon )
      {
        _ItemAtWall[2*iDim].push_back(pItem);
        continue;
      }
      if ( fabs(pItem->min_[iDim]-_RoomSize[iDim])<epsilon && fabs(pItem->max_[iDim]-_RoomSize[iDim])<epsilon )
      {
        _ItemAtWall[2*iDim+1].push_back(pItem);
        continue;
      }
    }
  }
}


// return index of wall at which item is located
int ProblemGeometry::GetWallItemWall(const Item& item)
{
  // Returnvalue: Min x1 -> 0, Max x1->1, Min x2 -> 2 Max x2 -> 3, Min x3 -> 4 Max x3 -> 5
  double epsilon = 1e-8;
  for (int iDim=0;iDim<GetDim();iDim++) {
    if ( (fabs(item.min_[iDim])<epsilon) && (fabs(item.max_[iDim])<epsilon) )
      return 2*iDim;
    if ( (fabs(item.min_[iDim]-_RoomSize[iDim])<epsilon) && (fabs(item.max_[iDim]-_RoomSize[iDim])<epsilon) )
      return 2*iDim+1;
  }
  // We should never get here
  assert(false);
  return -1;
}

// Update p_mesh to include BoundaryMarkers for all nodes that have Dirichlet
void ProblemGeometry::SetBoundaryInfo(libMesh::Mesh* p_mesh)
{
  MY_DBG_PRINT("");
  const double eps = 1e-8;
//  p_mesh->boundary_info->clear();

  //mark element sides
  std::vector<double> Pos;
  Pos.resize(GetDim());
  MeshBase::const_element_iterator el = p_mesh->active_local_elements_begin();
  const MeshBase::const_element_iterator end_el = p_mesh->active_local_elements_end();
  // loop over all active local elements
  for ( ; el != end_el; ++el) {
    for (unsigned int side=0; side<(*el)->n_sides(); side++) {
      if ((*el)->neighbor(side) == NULL) {
        // this is element side on boundary
        libMesh::Point Center;
        for (unsigned int iPt=0; iPt<(*el)->n_nodes();iPt++)
          if ((*el)->is_node_on_side(iPt,side))
            Center += (*el)->point(iPt);
        Center = (1.0/GetDim()) * Center;
        for (unsigned int iDim=0;iDim<GetDim();iDim++)
          Pos[iDim] = Center(iDim);
        int BoundaryMarker = GetBoundaryMarker(Pos);
        assert(BoundaryMarker!=-1);
//        p_mesh->boundary_info->add_side(*el,side,BoundaryMarker);

        // mark nodes if Dirichlet cond.: Those are set at points (not sides) by putting diagonals into the matrix
        if ( _BoundCond[BoundaryMarker]->IsPhiDirichlet() || _BoundCond[BoundaryMarker]->IsTDirichlet() ) {
          // mark side points
          for (unsigned int iPt=0; iPt<(*el)->n_nodes();iPt++)
            if ((*el)->is_node_on_side(iPt,side)) {
              MY_DBG_PRINT("added Dirichlet Boundary marker" << (*el)->get_node(iPt));
              p_mesh->boundary_info->add_node((*el)->get_node(iPt),BoundaryMarker);
            }
        }
      }
    }
  }
}

/////////////////////////////////////////////////////////////////////
//  3D
void ProblemGeometry::SetFacette(tetgenio::facet *f, int Base, int Offset0, int Offset1, int Offset2, int Offset3, int nHoles)
{
  tetgenio::init(f);
  f->numberofholes=nHoles;
  if (nHoles>0)
    f->holelist = new REAL[f->numberofholes * 3];
  f->numberofpolygons=1+nHoles;
  f->polygonlist = new tetgenio::polygon[f->numberofpolygons];
  tetgenio::polygon *p = f->polygonlist;
  tetgenio::init(p);
  p->numberofvertices = 4;
  p->vertexlist = new int[p->numberofvertices];
  p->vertexlist[0]=Base+Offset0;
  p->vertexlist[1]=Base+Offset1;
  p->vertexlist[2]=Base+Offset2;
  p->vertexlist[3]=Base+Offset3;
}

void ProblemGeometry::SetPoint(REAL* pointlist,int idx,REAL x,REAL y,REAL z)
{
  pointlist[3*idx+0]=x;
  pointlist[3*idx+1]=y;
  pointlist[3*idx+2]=z;
}


void ProblemGeometry::Tetgen2Mesh(const tetgenio& tet, libMesh::UnstructuredMesh* p_mesh)
{
  const double eps(1e-8);
  assert(GetDim()==3);
  int dim = GetDim();
  p_mesh->set_mesh_dimension(3);
  // Add Nodes:
  unsigned long int iPt;
  int BoundaryMarker;
  for (iPt=0;iPt<tet.numberofpoints;iPt++)
    p_mesh->add_point(libMesh::Point(tet.pointlist[3*iPt+0],tet.pointlist[3*iPt+1],tet.pointlist[3*iPt+2]) );
  // elements
  libMesh::Elem *pElem;
  for (unsigned int iEl=0; iEl<tet.numberoftetrahedra; iEl++) {
    pElem = new libMesh::Tet4;
    for (iPt=0;iPt<4;iPt++) {
      pElem->set_node(iPt) = p_mesh->node_ptr(tet.tetrahedronlist[4*iEl+iPt]);
    }
    // Finally, add this element to the mesh
    p_mesh->add_elem(pElem);
  }

  if (tet.neighborlist) {
    for (unsigned int iEl=0; iEl<tet.numberoftetrahedra; iEl++) {
      for (int iNeig=0;iNeig<dim+1;++iNeig) {
        if (tet.neighborlist[4*iEl+iNeig]==-1) {
          libMesh::Tet4 *pElem = static_cast<libMesh::Tet4*>(p_mesh->elem(iEl));
          int iBdNode1 = tet.tetrahedronlist[(dim+1)*iEl+((iNeig+1)%(dim+1))];
          int iBdNode2 = tet.tetrahedronlist[(dim+1)*iEl+((iNeig+2)%(dim+1))];
          int iBdNode3 = tet.tetrahedronlist[(dim+1)*iEl+((iNeig+3)%(dim+1))];
          int sum = iBdNode1+iBdNode2+iBdNode3;
          for (int iSide=0;iSide<dim+1;iSide++) {
            // sum is correct -> side is correct (holds only, if side has one node less than elem)
            if ( (pElem->node(pElem->side_nodes_map[iSide][0])
                  +pElem->node(pElem->side_nodes_map[iSide][1])
                  +pElem->node(pElem->side_nodes_map[iSide][2])) == sum ) {
              std::vector<double> Center;
              Center.resize(dim);
              for (int iDim=0;iDim<dim;iDim++)
                Center[iDim] = (tet.pointlist[dim*iBdNode1+iDim] + tet.pointlist[dim*iBdNode2+iDim]+ tet.pointlist[dim*iBdNode3+iDim])/dim;
              int BoundaryMarker = GetBoundaryMarker(Center);
              assert(BoundaryMarker!=-1);
              p_mesh->boundary_info->add_side(iEl,iSide,BoundaryMarker);
              if ( (_BoundCond[BoundaryMarker]->IsPhiDirichlet()) ||
                   (_BoundCond[BoundaryMarker]->IsTDirichlet() ) ) {
                p_mesh->boundary_info->add_node(iBdNode1,BoundaryMarker);
                p_mesh->boundary_info->add_node(iBdNode2,BoundaryMarker);
                p_mesh->boundary_info->add_node(iBdNode3,BoundaryMarker);
                // std::cout << "Dirichlet at " << iBdNode1 << ", " << iBdNode2 << ", " << iBdNode3 << ", " << "Center: " << Center << std::endl;
              }
            }
          }
        }
      }
    }
  }
}

void ProblemGeometry::CreateMesh(libMesh::UnstructuredMesh* p_mesh)
{
  switch (GetDim()) {
  case 2:
    CreateMesh2D(p_mesh);
    break;
  case 3:
    CreateMesh3D(p_mesh);
    break;
  default:
    std::cout << "Dimension not implemented" << std::endl;
    exit(0);
  }
}

void ProblemGeometry::CreateMesh3D(libMesh::UnstructuredMesh* p_mesh)
{
  tetgenio tetgen_in;
  //tetgenmesh      tetgen_mesh;

  tetgen_in.mesh_dim = 3;  // Three-dimemsional accoordinates.
  tetgen_in.numberofpointattributes = 0;  // no point attribute.
  tetgen_in.firstnumber = 0;
  tetgen_in.numberofcorners = 1;
  tetgenio::facet *FacWall;
  tetgenio::polygon *p;

  tetgen_in.numberofpoints = 8+8*_Equip.size()+4*(_AC.size()+_Exh.size());
  tetgen_in.pointlist = new REAL[3*tetgen_in.numberofpoints];
  int iPtInList = 0;
  int iFacInList = 0;

  SetPoint(tetgen_in.pointlist,iPtInList++,0             ,0             ,0             );
  SetPoint(tetgen_in.pointlist,iPtInList++,0             ,0             ,_RoomSize[2]);
  SetPoint(tetgen_in.pointlist,iPtInList++,0             ,_RoomSize[1],0             );
  SetPoint(tetgen_in.pointlist,iPtInList++,0             ,_RoomSize[1],_RoomSize[2]);
  SetPoint(tetgen_in.pointlist,iPtInList++,_RoomSize[0],0             ,0             );
  SetPoint(tetgen_in.pointlist,iPtInList++,_RoomSize[0],0             ,_RoomSize[2]);
  SetPoint(tetgen_in.pointlist,iPtInList++,_RoomSize[0],_RoomSize[1],0             );
  SetPoint(tetgen_in.pointlist,iPtInList++,_RoomSize[0],_RoomSize[1],_RoomSize[2]);

  // Where are the AC's and the Exhausts?
  ValidateItemAtWall();
  std::list<Item*>::iterator ppItem;

  int NumOfWallItems[4] = {0,0,0,0};
  unsigned int iItem;
  for (iItem=0;iItem<_AC.size();iItem++)
    NumOfWallItems[GetWallItemWall(_AC[iItem])]++;

  for (iItem=0;iItem<_Exh.size();iItem++)
    NumOfWallItems[GetWallItemWall(_Exh[iItem])]++;

  // facet marker:0: normal boundary condition, rest increasing
  int iFacMakerUsed = 1;
  tetgen_in.numberoffacets = 6+5*_Equip.size()+_AC.size()+_Exh.size();
  tetgen_in.facetlist = new tetgenio::facet[tetgen_in.numberoffacets];
  tetgen_in.facetmarkerlist = new int[tetgen_in.numberoffacets];
  memset(tetgen_in.facetmarkerlist,(int)0,tetgen_in.numberoffacets*sizeof(int));
  int iHoleInList;      // #Hole for this Facette
  int iFirstItemPt;
  // floor
  iHoleInList = 0;
  tetgen_in.facetmarkerlist[iFacInList] = 0;
  FacWall = tetgen_in.facetlist+(iFacInList++);
  SetFacette(FacWall, 0, 0, 2, 6, 4,_Equip.size()); // Set Facet of wall

  for (iItem=0;iItem<_Equip.size();iItem++) {
    iFirstItemPt = iPtInList;
    // Add Points of Item
    SetPoint(tetgen_in.pointlist,iPtInList++,_Equip[iItem].min_[0],_Equip[iItem].min_[1],_Equip[iItem].min_[2]);
    SetPoint(tetgen_in.pointlist,iPtInList++,_Equip[iItem].min_[0],_Equip[iItem].min_[1],_Equip[iItem].max_[2]);
    SetPoint(tetgen_in.pointlist,iPtInList++,_Equip[iItem].min_[0],_Equip[iItem].max_[1],_Equip[iItem].min_[2]);
    SetPoint(tetgen_in.pointlist,iPtInList++,_Equip[iItem].min_[0],_Equip[iItem].max_[1],_Equip[iItem].max_[2]);
    SetPoint(tetgen_in.pointlist,iPtInList++,_Equip[iItem].max_[0],_Equip[iItem].min_[1],_Equip[iItem].min_[2]);
    SetPoint(tetgen_in.pointlist,iPtInList++,_Equip[iItem].max_[0],_Equip[iItem].min_[1],_Equip[iItem].max_[2]);
    SetPoint(tetgen_in.pointlist,iPtInList++,_Equip[iItem].max_[0],_Equip[iItem].max_[1],_Equip[iItem].min_[2]);
    SetPoint(tetgen_in.pointlist,iPtInList++,_Equip[iItem].max_[0],_Equip[iItem].max_[1],_Equip[iItem].max_[2]);

    // Set Hole in Wall Facet
    p = FacWall->polygonlist+1+iHoleInList;
    p->numberofvertices = 4;
    p->vertexlist = new int[p->numberofvertices];
    p->vertexlist[0]=iFirstItemPt+0;
    p->vertexlist[1]=iFirstItemPt+2;
    p->vertexlist[2]=iFirstItemPt+6;
    p->vertexlist[3]=iFirstItemPt+4;
    FacWall->holelist[3*iHoleInList+0] = 0.5*(_Equip[iItem].min_[0]+_Equip[iItem].max_[0]);
    FacWall->holelist[3*iHoleInList+1] = 0.5*(_Equip[iItem].min_[1]+_Equip[iItem].max_[1]);
    FacWall->holelist[3*iHoleInList+2] = 0.0;
    iHoleInList++;

    // Add inner Facets (Equipment Boundary)
    SetFacette(tetgen_in.facetlist+(iFacInList++), iFirstItemPt, 1, 3, 7, 5);
    SetFacette(tetgen_in.facetlist+(iFacInList++), iFirstItemPt, 0, 1, 3, 2);
    SetFacette(tetgen_in.facetlist+(iFacInList++), iFirstItemPt, 4, 5, 7, 6);
    SetFacette(tetgen_in.facetlist+(iFacInList++), iFirstItemPt, 0, 1, 5, 4);
    SetFacette(tetgen_in.facetlist+(iFacInList++), iFirstItemPt, 2, 3, 7, 6);
  }

  // ceiling
  iHoleInList = 0;
  tetgen_in.facetmarkerlist[iFacInList] = 0;
  FacWall = tetgen_in.facetlist+(iFacInList++);
  tetgenio::init(FacWall);
  FacWall->numberofholes=0;
  FacWall->numberofpolygons=1+FacWall->numberofholes;
  FacWall->polygonlist = new tetgenio::polygon[FacWall->numberofpolygons];
  p = FacWall->polygonlist;
  tetgenio::init(p);
  p->numberofvertices = 4;
  p->vertexlist = new int[p->numberofvertices];
  p->vertexlist[0]=1;
  p->vertexlist[1]=3;
  p->vertexlist[2]=7;
  p->vertexlist[3]=5;

  Item* pItem;

  // x0==0
  // Add wall facet
  int iWall(0);
  iHoleInList = 0;
  tetgen_in.facetmarkerlist[iFacInList] = 0;
  FacWall = tetgen_in.facetlist+(iFacInList++);
  SetFacette(FacWall, 0, 0, 1, 3, 2,NumOfWallItems[iWall]);
  // Handle items, i.e. add hole in wall, add points, add facet
  ppItem=_ItemAtWall[iWall].begin();
  for (;ppItem!=_ItemAtWall[iWall].end();ppItem++) {
    pItem = *ppItem;
    iFirstItemPt = iPtInList;
    SetPoint(tetgen_in.pointlist,iPtInList++,pItem->min_[0],pItem->min_[1],pItem->min_[2]);
    SetPoint(tetgen_in.pointlist,iPtInList++,pItem->min_[0],pItem->min_[1],pItem->max_[2]);
    SetPoint(tetgen_in.pointlist,iPtInList++,pItem->min_[0],pItem->max_[1],pItem->min_[2]);
    SetPoint(tetgen_in.pointlist,iPtInList++,pItem->min_[0],pItem->max_[1],pItem->max_[2]);
    p = FacWall->polygonlist+1+iHoleInList;
    p->numberofvertices = 4;
    p->vertexlist = new int[p->numberofvertices];
    p->vertexlist[0]=iFirstItemPt+0;
    p->vertexlist[1]=iFirstItemPt+1;
    p->vertexlist[2]=iFirstItemPt+3;
    p->vertexlist[3]=iFirstItemPt+2;
    FacWall->holelist[3*iHoleInList+0] = 0.0;
    FacWall->holelist[3*iHoleInList+1] = 0.5*(pItem->min_[1]+pItem->max_[1]);
    FacWall->holelist[3*iHoleInList+2] = 0.5*(pItem->min_[2]+pItem->max_[2]);
    iHoleInList++;
    tetgen_in.facetmarkerlist[iFacInList] = iFacMakerUsed++;
    SetFacette(tetgen_in.facetlist+(iFacInList++), iFirstItemPt, 0, 1, 3, 2);
  }

  // x0==_RoomSize[0]
  iWall++;
  iHoleInList = 0;
  tetgen_in.facetmarkerlist[iFacInList] = 0;
  FacWall = tetgen_in.facetlist+(iFacInList++);
  SetFacette(FacWall, 0, 4, 5, 7, 6,NumOfWallItems[iWall]);
  ppItem=_ItemAtWall[iWall].begin();
  for (;ppItem!=_ItemAtWall[iWall].end();ppItem++) {
    pItem = *ppItem;
    iFirstItemPt = iPtInList;
    SetPoint(tetgen_in.pointlist,iPtInList++,pItem->max_[0],pItem->min_[1],pItem->min_[2]);
    SetPoint(tetgen_in.pointlist,iPtInList++,pItem->max_[0],pItem->min_[1],pItem->max_[2]);
    SetPoint(tetgen_in.pointlist,iPtInList++,pItem->max_[0],pItem->max_[1],pItem->min_[2]);
    SetPoint(tetgen_in.pointlist,iPtInList++,pItem->max_[0],pItem->max_[1],pItem->max_[2]);
    p = FacWall->polygonlist+1+iHoleInList;
    p->numberofvertices = 4;
    p->vertexlist = new int[p->numberofvertices];
    p->vertexlist[0]=iFirstItemPt+0;
    p->vertexlist[1]=iFirstItemPt+1;
    p->vertexlist[2]=iFirstItemPt+3;
    p->vertexlist[3]=iFirstItemPt+2;
    FacWall->holelist[3*iHoleInList+0] = _RoomSize[0];
    FacWall->holelist[3*iHoleInList+1] = 0.5*(pItem->min_[1]+pItem->max_[1]);
    FacWall->holelist[3*iHoleInList+2] = 0.5*(pItem->min_[2]+pItem->max_[2]);
    iHoleInList++;
    tetgen_in.facetmarkerlist[iFacInList] = iFacMakerUsed++;
    SetFacette(tetgen_in.facetlist+(iFacInList++), iFirstItemPt, 0, 1, 3, 2);
  }

  // x1==0
  iWall++;
  iHoleInList = 0;
  tetgen_in.facetmarkerlist[iFacInList] = 0;
  FacWall = tetgen_in.facetlist+(iFacInList++);
  SetFacette(FacWall, 0, 0, 1, 5, 4, NumOfWallItems[iWall]);
  ppItem=_ItemAtWall[iWall].begin();
  for (;ppItem!=_ItemAtWall[iWall].end();ppItem++) {
    pItem = *ppItem;
    iFirstItemPt = iPtInList;
    SetPoint(tetgen_in.pointlist,iPtInList++,pItem->min_[0],pItem->min_[1],pItem->min_[2]);
    SetPoint(tetgen_in.pointlist,iPtInList++,pItem->min_[0],pItem->min_[1],pItem->max_[2]);
    SetPoint(tetgen_in.pointlist,iPtInList++,pItem->max_[0],pItem->min_[1],pItem->min_[2]);
    SetPoint(tetgen_in.pointlist,iPtInList++,pItem->max_[0],pItem->min_[1],pItem->max_[2]);
    p = FacWall->polygonlist+1+iHoleInList;
    p->numberofvertices = 4;
    p->vertexlist = new int[p->numberofvertices];
    p->vertexlist[0]=iFirstItemPt+0;
    p->vertexlist[1]=iFirstItemPt+1;
    p->vertexlist[2]=iFirstItemPt+3;
    p->vertexlist[3]=iFirstItemPt+2;
    FacWall->holelist[3*iHoleInList+0] = 0.5*(pItem->min_[0]+pItem->max_[0]);
    FacWall->holelist[3*iHoleInList+1] = 0.0;
    FacWall->holelist[3*iHoleInList+2] = 0.5*(pItem->min_[2]+pItem->max_[2]);
    iHoleInList++;
    tetgen_in.facetmarkerlist[iFacInList] = iFacMakerUsed++;
    SetFacette(tetgen_in.facetlist+(iFacInList++), iFirstItemPt, 0, 1, 3, 2);
  }

  // x1==_RoomSize[1]
  iWall++;
  iHoleInList = 0;
  tetgen_in.facetmarkerlist[iFacInList] = 0;
  FacWall = tetgen_in.facetlist+(iFacInList++);
  SetFacette(FacWall, 0, 2, 3, 7, 6, NumOfWallItems[iWall]);
  ppItem=_ItemAtWall[iWall].begin();
  for (;ppItem!=_ItemAtWall[iWall].end();ppItem++) {
    pItem = *ppItem;
    iFirstItemPt = iPtInList;
    SetPoint(tetgen_in.pointlist,iPtInList++,pItem->min_[0],pItem->min_[1],pItem->min_[2]);
    SetPoint(tetgen_in.pointlist,iPtInList++,pItem->min_[0],pItem->min_[1],pItem->max_[2]);
    SetPoint(tetgen_in.pointlist,iPtInList++,pItem->max_[0],pItem->min_[1],pItem->min_[2]);
    SetPoint(tetgen_in.pointlist,iPtInList++,pItem->max_[0],pItem->min_[1],pItem->max_[2]);
    p = FacWall->polygonlist+1+iHoleInList;
    p->numberofvertices = 4;
    p->vertexlist = new int[p->numberofvertices];
    p->vertexlist[0]=iFirstItemPt+0;
    p->vertexlist[1]=iFirstItemPt+1;
    p->vertexlist[2]=iFirstItemPt+3;
    p->vertexlist[3]=iFirstItemPt+2;
    FacWall->holelist[3*iHoleInList+0] = 0.5*(pItem->min_[0]+pItem->max_[0]);
    FacWall->holelist[3*iHoleInList+1] = _RoomSize[1];
    FacWall->holelist[3*iHoleInList+2] = 0.5*(pItem->min_[2]+pItem->max_[2]);
    iHoleInList++;
    tetgen_in.facetmarkerlist[iFacInList] = iFacMakerUsed++;
    SetFacette(tetgen_in.facetlist+(iFacInList++), iFirstItemPt, 0, 1, 3, 2);
  }

  tetgenio tetgen_out;

  bool bCallExternal=true;
  if (bCallExternal) {
    int ProcID = GetProcID();
    if(ProcID==0) {
      std::ofstream os("Mesh3DRaw.poly",std::ios::out);
      PrintTetgenMesh(tetgen_in, os);
      char strBuf[256];
      sprintf(strBuf,"-nQpqa%e",_h*_h*_h);
      //sprintf(strBuf,"zpQqa%e",h*h*h);  // tetgen in silent mode
      char Buf[2014];
      sprintf(Buf,"tetgen %s Mesh3DRaw.poly",strBuf);
      std::cout << "calling \"" << Buf << "\"...";
      std::cout.flush();
      system(Buf);
      std::cout << " ... finished" << std::endl;
   
      ReadNodeFile("Mesh3DRaw.1.node",&tetgen_out);
      ReadEleFile("Mesh3DRaw.1.ele",&tetgen_out);
      ReadNeighFile("Mesh3DRaw.1.neigh",&tetgen_out);
    }
    MPI_Bcast(&(tetgen_out.numberofpoints), 1, MPI_INT, 0, MPI_COMM_WORLD);
    //MY_DBG_PRINT("tetgen_out.numberofpoints:" << sizeof(tetgen_out.numberofpoints) << ", " << tetgen_out.numberofpoints)

    if(ProcID!=0)
      tetgen_out.pointlist = new double[3*tetgen_out.numberofpoints];
    MPI_Bcast(tetgen_out.pointlist,3*tetgen_out.numberofpoints,MPI_DOUBLE,0,MPI_COMM_WORLD);

    MPI_Bcast(&(tetgen_out.numberoftetrahedra), 1, MPI_INT, 0, MPI_COMM_WORLD);
    if(ProcID!=0)
      tetgen_out.tetrahedronlist = new int[4*tetgen_out.numberoftetrahedra];
    MPI_Bcast(tetgen_out.tetrahedronlist,4*tetgen_out.numberoftetrahedra,MPI_INT,0,MPI_COMM_WORLD);

    if(ProcID!=0)
      tetgen_out.neighborlist = new int[4*tetgen_out.numberoftetrahedra];
    MPI_Bcast(tetgen_out.neighborlist,4*tetgen_out.numberoftetrahedra,MPI_INT,0,MPI_COMM_WORLD);

  }
  else {
    bool PrintMeshingData=true;
    if (PrintMeshingData) {
      std::ofstream os("Mesh3DRaw.poly",std::ios::out);
      PrintTetgenMesh(tetgen_in, os);
    }
    tetgenbehavior  tetgen_beh;
    char strBuf[256];
    sprintf(strBuf,"npqa%f",_h*_h*_h);
    //sprintf(strBuf,"-n");
    //sprintf(strBuf,"%s","pq1.414a0.1");
    //sprintf(strBuf,"zpQqa%f",h*h*h);  // tetgen in silent mode
    //tetgen_beh.parse_commandline(strBuf);
    //tetrahedralize(&tetgen_beh, &tetgen_in, &tetgen_out);
    printf("%s\n",strBuf);
    //tetrahedralize(&tetgen_beh, &tetgen_in, &tetgen_out);
    //tetrahedralize(&tetgen_beh, &tetgen_in, &tetgen_out);
    tetrahedralize(strBuf, &tetgen_in, &tetgen_out);
  }
  Tetgen2Mesh(tetgen_out, p_mesh);
  p_mesh->prepare_for_use();
// TODO: Clearify what's going on here, why program termination?
//  tetgen_in.deinitialize();
//  tetgen_out.deinitialize();
}

void ProblemGeometry::CreateMesh2D(libMesh::UnstructuredMesh* p_mesh)
{
  assert(GetDim()==2);
  int dim = GetDim();
  libMesh::Triangle::triangulateio tri_in, tri_out;
  libMesh::Triangle::init(tri_in);
  libMesh::Triangle::init(tri_out);

  tri_in.numberofpoints = 4+4*_Equip.size()+ 2*(_AC.size()+_Exh.size());
  tri_in.pointlist = static_cast<REAL*>(std::malloc(2*tri_in.numberofpoints*sizeof(REAL)));
  int iPtInList = 0;
  tri_in.numberofsegments = 4+4*_Equip.size()+2*(_AC.size()+_Exh.size());
  tri_in.segmentlist = static_cast<int*> (std::malloc(tri_in.numberofsegments*2*sizeof(int)));
  int iSegInList = 0;

  ////////////////////////
  // Outer Boundary (including AC's, Exhausts,
  // Roundtrip: x1=0 -> x2=1 -> x1=1 -> x2=0

  double epsilon = 1e-8;

  std::vector<double> BdPts;
  BdPts.reserve(_AC.size()+_Exh.size());
  int iItem;
  int iPt;
  std::vector<CompareItem> BdItems;
  std::vector<CompareItem>::iterator pSortItem;
  ValidateItemAtWall();

  BdItems.reserve(_AC.size()+_Exh.size());
  std::list<Item*>::iterator ppItem;

  int iWall(0);
  CompareItem item;
  double start, end;

  tri_in.pointlist[2*iPtInList] = 0.0;
  tri_in.pointlist[2*iPtInList+1] = 0.0;
  ++iPtInList;

  // Outer Boundary: x1=0, (0,0)->(0,Max);
  iWall=0;
  start=0.0;
  end=_RoomSize[1];
  if (_ItemAtWall[iWall].empty()) {
    tri_in.pointlist[2*iPtInList] = 0.0;
    tri_in.pointlist[2*iPtInList+1] = end;
    tri_in.segmentlist[2*iSegInList] = iPtInList-1;
    tri_in.segmentlist[2*iSegInList+1] = iPtInList;
    ++iPtInList;
    ++iSegInList;
  }
  else {
    BdItems.clear();
    ppItem=_ItemAtWall[iWall].begin();
    for (;ppItem!=_ItemAtWall[iWall].end();++ppItem) {
      item.pItem = *ppItem;
      item.SortValue = (*ppItem)->min_[1];
      BdItems.push_back(item);
    }
    sort(BdItems.begin(), BdItems.end(),CompareItems);
    for (iItem=0;iItem<BdItems.size();++iItem) {
      tri_in.pointlist[2*iPtInList] = 0.0;
      tri_in.pointlist[2*iPtInList+1] = BdItems[iItem].pItem->min_[1];
      tri_in.segmentlist[2*iSegInList] = iPtInList-1;
      tri_in.segmentlist[2*iSegInList+1] = iPtInList;
      ++iPtInList;
      ++iSegInList;
      tri_in.pointlist[2*iPtInList] = 0.0;
      tri_in.pointlist[2*iPtInList+1] = BdItems[iItem].pItem->max_[1];
      tri_in.segmentlist[2*iSegInList] = iPtInList-1;
      tri_in.segmentlist[2*iSegInList+1] = iPtInList;
      ++iPtInList;
      ++iSegInList;
    }
    if (fabs(BdItems.back().pItem->max_[1]-end)>epsilon) {
      tri_in.pointlist[2*iPtInList] = 0.0;
      tri_in.pointlist[2*iPtInList+1] = end;
      tri_in.segmentlist[2*iSegInList] = iPtInList-1;
      tri_in.segmentlist[2*iSegInList+1] = iPtInList;
      ++iPtInList;
      ++iSegInList;
    }
  }

  // Outer Boundary: x2=Max, (0,Max)->(Max,Max);
  iWall = 3;
  start=0.0;
  end=_RoomSize[0];
  if (_ItemAtWall[iWall].empty()) {
    tri_in.pointlist[2*iPtInList] = end;
    tri_in.pointlist[2*iPtInList+1] = _RoomSize[1];
    tri_in.segmentlist[2*iSegInList] = iPtInList-1;
    tri_in.segmentlist[2*iSegInList+1] = iPtInList;
    ++iPtInList;
    ++iSegInList;
  }
  else {
    BdItems.clear();
    ppItem=_ItemAtWall[iWall].begin();
    for (;ppItem!=_ItemAtWall[iWall].end();++ppItem) {
      item.pItem = *ppItem;
      item.SortValue = (*ppItem)->min_[0];
      BdItems.push_back(item);
    }
    sort(BdItems.begin(), BdItems.end(),CompareItems);
    for (iItem=0;iItem<BdItems.size();++iItem) {
      tri_in.pointlist[2*iPtInList] = BdItems[iItem].pItem->min_[0];
      tri_in.pointlist[2*iPtInList+1] = _RoomSize[1];
      tri_in.segmentlist[2*iSegInList] = iPtInList-1;
      tri_in.segmentlist[2*iSegInList+1] = iPtInList;
      ++iPtInList;
      ++iSegInList;
      tri_in.pointlist[2*iPtInList] = BdItems[iItem].pItem->max_[0];
      tri_in.pointlist[2*iPtInList+1] = _RoomSize[1];
      tri_in.segmentlist[2*iSegInList] = iPtInList-1;
      tri_in.segmentlist[2*iSegInList+1] = iPtInList;
      ++iPtInList;
      ++iSegInList;
    }
    if (fabs(BdItems.back().pItem->max_[0]-end)>epsilon) {
      tri_in.pointlist[2*iPtInList] = end;
      tri_in.pointlist[2*iPtInList+1] = _RoomSize[1];
      tri_in.segmentlist[2*iSegInList] = iPtInList-1;
      tri_in.segmentlist[2*iSegInList+1] = iPtInList;
      ++iPtInList;
      ++iSegInList;
    }
  }
  // Outer Boundary: x1=Max, (Max,Max)->(Max,0);
  iWall = 1;
  start=_RoomSize[1];
  end=0.0;
  if (_ItemAtWall[iWall].empty()) {
    tri_in.pointlist[2*iPtInList] = _RoomSize[0];
    tri_in.pointlist[2*iPtInList+1] = end;
    tri_in.segmentlist[2*iSegInList] = iPtInList-1;
    tri_in.segmentlist[2*iSegInList+1] = iPtInList;
    ++iPtInList;
    ++iSegInList;
  }
  else {
    BdItems.clear();
    ppItem=_ItemAtWall[iWall].begin();
    for (;ppItem!=_ItemAtWall[iWall].end();++ppItem) {
      item.pItem = *ppItem;
      item.SortValue = (*ppItem)->min_[1];
      BdItems.push_back(item);
    }
    sort(BdItems.begin(), BdItems.end(),CompareItems);
    for (iItem=BdItems.size()-1;iItem>=0;--iItem) {
      tri_in.pointlist[2*iPtInList] = _RoomSize[0];
      tri_in.pointlist[2*iPtInList+1] = BdItems[iItem].pItem->max_[1];
      tri_in.segmentlist[2*iSegInList] = iPtInList-1;
      tri_in.segmentlist[2*iSegInList+1] = iPtInList;
      ++iPtInList;
      ++iSegInList;
      tri_in.pointlist[2*iPtInList] = _RoomSize[0];
      tri_in.pointlist[2*iPtInList+1] = BdItems[iItem].pItem->min_[1];
      tri_in.segmentlist[2*iSegInList] = iPtInList-1;
      tri_in.segmentlist[2*iSegInList+1] = iPtInList;
      ++iPtInList;
      ++iSegInList;
    }
    if (fabs(BdItems.back().pItem->max_[1]-end)>epsilon) {
      tri_in.pointlist[2*iPtInList] = _RoomSize[0];
      tri_in.pointlist[2*iPtInList+1] = end;
      tri_in.segmentlist[2*iSegInList] = iPtInList-1;
      tri_in.segmentlist[2*iSegInList+1] = iPtInList;
      ++iPtInList;
      ++iSegInList;
    }
  }

  // Outer Boundary: x2=0, (Max,0)->(0,0);
  iWall = 2;
  start=_RoomSize[0];
  end=0.0;
  if (_ItemAtWall[iWall].empty()) {
    tri_in.segmentlist[2*iSegInList] = iPtInList-1;
    tri_in.segmentlist[2*iSegInList+1] = 0;
    ++iSegInList;
  }
  else {
    BdItems.clear();
    ppItem=_ItemAtWall[iWall].begin();
    for (;ppItem!=_ItemAtWall[iWall].end();++ppItem) {
      item.pItem = *ppItem;
      item.SortValue = (*ppItem)->min_[0];
      BdItems.push_back(item);
    }
    sort(BdItems.begin(), BdItems.end(),CompareItems);
    for (iItem=BdItems.size()-1;iItem>=0;--iItem) {
      tri_in.pointlist[2*iPtInList] = BdItems[iItem].pItem->max_[0];
      tri_in.pointlist[2*iPtInList+1] = 0.0;
      tri_in.segmentlist[2*iSegInList] = iPtInList-1;
      tri_in.segmentlist[2*iSegInList+1] = iPtInList;
      ++iPtInList;
      ++iSegInList;
      tri_in.pointlist[2*iPtInList] = BdItems[iItem].pItem->min_[0];
      tri_in.pointlist[2*iPtInList+1] = 0.0;
      tri_in.segmentlist[2*iSegInList] = iPtInList-1;
      tri_in.segmentlist[2*iSegInList+1] = iPtInList;
      ++iPtInList;
      ++iSegInList;
    }
    if (fabs(BdItems.back().pItem->max_[0]-end)>epsilon) {
      tri_in.segmentlist[2*iSegInList] = iPtInList-1;
      tri_in.segmentlist[2*iSegInList+1] = 0;
      ++iSegInList;
    }
    else {
      tri_in.segmentlist[2*iSegInList+1] = 0;
    }
  }
  ///////////////////////////////////////////
  // Add Holes
  tri_in.numberofholes = _Equip.size();
  tri_in.holelist = static_cast<REAL*>(std::malloc(2*tri_in.numberofholes*sizeof(REAL)));
  int iHoleInList = 0;
  Point Center;
  for (iItem=0; iItem<_Equip.size();iItem++) {
    // Add hole boundary points and segments
    tri_in.pointlist[2*iPtInList+0] = _Equip[iItem].min_[0];
    tri_in.pointlist[2*iPtInList+1] = _Equip[iItem].min_[1];
    tri_in.pointlist[2*iPtInList+2] = _Equip[iItem].min_[0];
    tri_in.pointlist[2*iPtInList+3] = _Equip[iItem].max_[1];
    tri_in.pointlist[2*iPtInList+4] = _Equip[iItem].max_[0];
    tri_in.pointlist[2*iPtInList+5] = _Equip[iItem].min_[1];
    tri_in.pointlist[2*iPtInList+6] = _Equip[iItem].max_[0];
    tri_in.pointlist[2*iPtInList+7] = _Equip[iItem].max_[1];
    tri_in.segmentlist[2*iSegInList+0] = iPtInList+0;
    tri_in.segmentlist[2*iSegInList+1] = iPtInList+1;
    tri_in.segmentlist[2*iSegInList+2] = iPtInList+1;
    tri_in.segmentlist[2*iSegInList+3] = iPtInList+3;
    tri_in.segmentlist[2*iSegInList+4] = iPtInList+3;
    tri_in.segmentlist[2*iSegInList+5] = iPtInList+2;
    tri_in.segmentlist[2*iSegInList+6] = iPtInList+2;
    tri_in.segmentlist[2*iSegInList+7] = iPtInList+0;
    iPtInList+=4;
    iSegInList+=4;
    tri_in.holelist[2*iHoleInList+0] = 0.5*(_Equip[iItem].min_[0]+_Equip[iItem].max_[0]);
    tri_in.holelist[2*iHoleInList+1] = 0.5*(_Equip[iItem].min_[1]+_Equip[iItem].max_[1]);
    iHoleInList++;
  }

  bool PrintMeshingData=false;
  if (PrintMeshingData) {
    std::ofstream file("Mesh2DIn.poly");
    PrintTriangleMesh(tri_in, file);
  }
  char strBuf[256];
  sprintf(strBuf,"znqpa%f",_h*_h);
  triangulate(strBuf, &tri_in, &tri_out, NULL);
  //libMesh::Triangle::copy_tri_to_mesh(tri_out,*p_mesh,TRI3);

  Triangle2Mesh(tri_out,p_mesh);
  p_mesh->prepare_for_use();
  libMesh::Triangle::destroy(tri_in,Triangle::INPUT);   // deletes also the memory
  libMesh::Triangle::destroy(tri_out,Triangle::OUTPUT);
}

void ProblemGeometry::ReadNodeFile(std::string str, tetgenio* tet)
{
  char Buf[1024];
  double Vals[4];
  int n_nodes, dim;
  int read;

  std::ifstream f(str.c_str(),std::ios::in);
  f.getline(Buf,1024);

  read = sscanf(Buf,"%d %d %lf %lf", &n_nodes,&dim,Vals,Vals+1);
  if ( (read<2) || (read>4) ) {
    std::string str("Can't read header line of node file:");
    str += Buf;
    throw std::runtime_error(str);
  }
  if (dim!=3) {
    std::string str("wrong dimension in node file:");
    str += Buf;
    throw std::runtime_error(str);
  }

  tet->numberofpoints = n_nodes;
  tet->pointlist = new double[3*tet->numberofpoints];
  int i_node;
  for (int iPt=0;iPt<tet->numberofpoints;iPt++) {
    f.getline(Buf,1024);
    read = sscanf(Buf,"%d %lf %lf %lf", &i_node,Vals,Vals+1,Vals+2);
    if ( read!=4 ) {
      std::string str("Can't read node from line:");
      str += Buf;
      throw std::runtime_error(str);
    }
    else {
      assert(iPt==i_node-1);
      tet->pointlist[3*iPt+0] = Vals[0];
      tet->pointlist[3*iPt+1] = Vals[1];
      tet->pointlist[3*iPt+2] = Vals[2];
    }
  }
}

void ProblemGeometry::ReadEleFile(std::string str, tetgenio* tet)
{
  char Buf[1024];
  int Vals[5];
  int n_elems, n_PtsPerTet;
  int read;

  std::ifstream f(str.c_str(),std::ios::in);
  f.getline(Buf,1024);

  read = sscanf(Buf,"%d %d %d %d", &n_elems,&n_PtsPerTet,Vals,Vals+1);
  if ( (read<2) || (read>4) ) {
    std::string str("Can't read header line of ele file:");
    str += Buf;
    throw std::runtime_error(str);
  }
  if (n_PtsPerTet!=4) {
    std::string str("Wrong number of points per tetrahedron:");
    str += Buf;
    throw std::runtime_error(str);
  }

  tet->numberoftetrahedra = n_elems;
  tet->tetrahedronlist = new int[4*tet->numberoftetrahedra];
  int i_elem;
  for (int i_elem=0;i_elem<tet->numberoftetrahedra;i_elem++) {
    f.getline(Buf,1024);
    read = sscanf(Buf,"%d %d %d %d %d", Vals,Vals+1,Vals+2,Vals+3,Vals+4);
    if ( read!=5 ) {
      std::string str("Can't read tetrahedron from line:");
      str += Buf;
      throw std::runtime_error(str);
    }
    else {
      assert(i_elem==Vals[0]-1);
      tet->tetrahedronlist[4*i_elem+0] = Vals[1]-1;
      tet->tetrahedronlist[4*i_elem+1] = Vals[2]-1;
      tet->tetrahedronlist[4*i_elem+2] = Vals[3]-1;
      tet->tetrahedronlist[4*i_elem+3] = Vals[4]-1;
    }
  }
}

void ProblemGeometry::ReadNeighFile(std::string str, tetgenio* tet)
{
  char Buf[1024];
  int Vals[5];
  int n_elems, n_NeighPerTet;
  int read;

  std::ifstream f(str.c_str(),std::ios::in);
  f.getline(Buf,1024);

  read = sscanf(Buf,"%d %d", &n_elems,&n_NeighPerTet);
  if ( (read<2) || (read>4) ) {
    std::string str("Can't read header line of neigh file:");
    str += Buf;
    throw std::runtime_error(str);
  }
  if (n_NeighPerTet!=4) {
    std::string str("Wrong number of neighbors per tetrahedron:");
    str += Buf;
    throw std::runtime_error(str);
  }

  assert(n_elems==tet->numberoftetrahedra);
  tet->neighborlist = new int[4*tet->numberoftetrahedra];
  int i_elem;
  for (int i_elem=0;i_elem<tet->numberoftetrahedra;i_elem++) {
    f.getline(Buf,1024);
    read = sscanf(Buf,"%d %d %d %d %d", Vals,Vals+1,Vals+2,Vals+3,Vals+4);
    if ( read!=5 ) {
      std::string str("Can't read neighbors from line:");
      str += Buf;
      throw std::runtime_error(str);
    }
    else {
      assert(i_elem==Vals[0]-1);
      tet->neighborlist[4*i_elem+0] = Vals[1]<0 ? -1 : (Vals[1]-1);
      tet->neighborlist[4*i_elem+1] = Vals[2]<0 ? -1 : (Vals[2]-1);
      tet->neighborlist[4*i_elem+2] = Vals[3]<0 ? -1 : (Vals[3]-1);
      tet->neighborlist[4*i_elem+3] = Vals[4]<0 ? -1 : (Vals[4]-1);
    }
  }
}

// for debugging:
void PrintTetgenMesh(const tetgenio& tet, std::ostream& os)
{
  os.flush();
  os << tet.numberofpoints << " " << tet.mesh_dim << std::endl;
  os << "# Nodes" << std::endl;
  for (int iPt=0;iPt<tet.numberofpoints;iPt++)
    os << iPt+1 << " " << tet.pointlist[3*iPt+0] << " " << tet.pointlist[3*iPt+1] << " " << tet.pointlist[3*iPt+2] << std::endl;
  os << "# Facets" << std::endl;
  tetgenio::facet *f;
  tetgenio::polygon *p;

  os << tet.numberoffacets << " 1" << std::endl;
  for (int iFac=0;iFac<tet.numberoffacets;iFac++) {
    f = tet.facetlist+iFac;
    os << f->numberofpolygons << " " << f->numberofholes << " " << iFac << std::endl;
    for (int iPoly=0; iPoly<f->numberofpolygons; iPoly++) {
      p = f->polygonlist+iPoly;
      os << p->numberofvertices << " ";
      for ( int iPt=0;iPt<p->numberofvertices; iPt++)
        os << p->vertexlist[iPt]+1 << " ";
      os << std::endl;
    }
    for (int iPt=0;iPt<f->numberofholes;iPt++) {
      os << iPt+1 << " ";
      for (int iDim=0;iDim<tet.mesh_dim;iDim++)
        os << f->holelist[iPt*tet.mesh_dim+iDim] << " ";
      os << std::endl;
    }
  }

  os << "# volume holes" << std::endl << "0" << std::endl << std::endl;
  os << "# regions" << std::endl << "0" << std::endl << std::endl;
  /*  os << tet.numberoftetrahedra << std::endl;
    for(int iEl=0;iEl<tet.numberoftetrahedra;iEl++)
      os << iEl << " " << tet.tetrahedronlist[4*iEl+0]+1 << " " << tet.tetrahedronlist[4*iEl+1]+1 << " " << tet.tetrahedronlist[4*iEl+2]+1 << " " << tet.tetrahedronlist[4*iEl+3]+1 << std::endl;
    os.flush();
  */
}

//////////////////////////////////////////////////////////////////////////////////////////
// 2D
void ProblemGeometry::Triangle2Mesh(const libMesh::Triangle::triangulateio& tri, libMesh::UnstructuredMesh* p_mesh)
{
  const double eps=1e-8;
  assert(GetDim()==2);
  int dim = GetDim();
  p_mesh->set_mesh_dimension(GetDim());
  // Add Nodes:
  unsigned long int iPt;
  p_mesh->reserve_elem(tri.numberoftriangles);
  p_mesh->reserve_nodes(tri.numberofpoints);
  for (iPt=0;iPt<tri.numberofpoints;iPt++)
    p_mesh->add_point(libMesh::Point(tri.pointlist[dim*iPt+0],tri.pointlist[dim*iPt+1],0),iPt);
  // elements
  libMesh::Tri3 *pElem;
  for (unsigned int iEl=0; iEl<tri.numberoftriangles; iEl++) {
    pElem = new libMesh::Tri3;
    for (iPt=0;iPt<dim+1;iPt++)
      pElem->set_node(iPt) = p_mesh->node_ptr(tri.trianglelist[(dim+1)*iEl+iPt]);
    // Finally, add this element to the mesh
    p_mesh->add_elem(pElem);
  }

  if (tri.neighborlist) {
    for (unsigned int iEl=0; iEl<tri.numberoftriangles; iEl++) {
      for (int iNeig=0;iNeig<dim+1;++iNeig) {
        if (tri.neighborlist[3*iEl+iNeig]==-1) {
          pElem = static_cast<libMesh::Tri3*>(p_mesh->elem(iEl));
          int iBdNode1 = tri.trianglelist[(dim+1)*iEl+((iNeig+1)%(dim+1))];
          int iBdNode2 = tri.trianglelist[(dim+1)*iEl+((iNeig+2)%(dim+1))];
          int sum = iBdNode1+iBdNode2;
          for (int iSide=0;iSide<GetDim()+1;iSide++) {
            // sum is correct -> sid is correct (holds only, if side has one node less than elem)
            if ( pElem->node(pElem->side_nodes_map[iSide][0])+pElem->node(pElem->side_nodes_map[iSide][1]) == sum ) {
              std::vector<double> Center;
              Center.resize(dim);
              for (int iDim=0;iDim<dim;iDim++)
                Center[iDim] = 0.5*(tri.pointlist[dim*iBdNode1+iDim] + tri.pointlist[dim*iBdNode2+iDim]);
              int BoundaryMarker = GetBoundaryMarker(Center);
              assert(BoundaryMarker!=-1);
              p_mesh->boundary_info->add_side(iEl,iSide,BoundaryMarker);
              if ( _BoundCond[BoundaryMarker]->IsPhiDirichlet() || _BoundCond[BoundaryMarker]->IsTDirichlet() ) {
                p_mesh->boundary_info->add_node(iBdNode1,BoundaryMarker);
                p_mesh->boundary_info->add_node(iBdNode2,BoundaryMarker);
              }
            }
          }
        }
      }
    }
  }
  p_mesh->prepare_for_use(false);
/*
  int dim = GetDim();
  libMesh::MeshTools::Generation::build_square(*p_mesh, 2, 2,0.0,1.0,0.0,1.0,TRI3);
  libMesh::MeshBase::const_element_iterator itCurEl = p_mesh->active_local_elements_begin();
  const MeshBase::const_element_iterator itEndEl = p_mesh->active_local_elements_end();
  for ( ; itCurEl != itEndEl; ++itCurEl) {
    const libMesh::Elem* CurElem = *itCurEl;
    for (unsigned int side=0; side<CurElem->n_sides(); side++) {
      if (CurElem->neighbor(side) == NULL) {
        AutoPtr<Elem> CurSide = CurElem->build_side(side);
        std::vector<double> Center;
        Center.resize(dim);
        libMesh::Point CenterPt = CurSide->centroid();
        for (int iDim=0;iDim<dim;iDim++)
          Center[iDim] = CenterPt(iDim);
        int BoundaryMarker = GetBoundaryMarker(Center);
        assert(BoundaryMarker!=-1);
        std::cout << "Adding boundary " << BoundaryMarker << "at point " << CenterPt << std::endl;
        p_mesh->boundary_info->add_side(CurElem,side,BoundaryMarker);
        if(_BoundCond[BoundaryMarker]->IsPhiDirichlet() || _BoundCond[BoundaryMarker]->IsTDirichlet() ) {
          p_mesh->boundary_info->add_node(CurSide->get_node(0),BoundaryMarker);
          p_mesh->boundary_info->add_node(CurSide->get_node(1),BoundaryMarker);
        }
      }
    }
  }
  p_mesh->prepare_for_use(false);
*/
}

void PrintTriangleMesh(const libMesh::Triangle::triangulateio& tri, std::ostream& os)
{
  os.flush();
  os << "# Nodes:" << std::endl;
  os << tri.numberofpoints << " 2 0 0" << std::endl;
  for (int iPt=0;iPt<tri.numberofpoints; ++iPt)
    os << iPt << " " << tri.pointlist[2*iPt] << " " << tri.pointlist[2*iPt+1] << std::endl;
  os << std::endl;
  os << "# Edges:" << std::endl;
  os << tri.numberofsegments << " " << (int)(tri.segmentmarkerlist!=NULL) << std::endl;
  for (int iEdge=0;iEdge<tri.numberofsegments;iEdge++) {
    os << iEdge << " " << tri.segmentlist[2*iEdge] << " " << tri.segmentlist[2*iEdge+1];
    if (tri.segmentmarkerlist!=NULL)
      os << " " << tri.segmentmarkerlist[iEdge];
    os << std::endl;
  }
  os << std::endl;
  os << "# Holes:" << std::endl;
  os << tri.numberofholes << " 0" << std::endl;
  for (int iHole=0;iHole<tri.numberofholes;iHole++)
    os << iHole << " " << tri.holelist[2*iHole] << " " << tri.holelist[2*iHole+1] << std::endl;
  os << std::endl;
  os << "# Region attrib., area constr.:" << std::endl;
  os << "0";
}


/*
std::ostream& operator << (std::ostream& os, const BoundaryCondition& BC)
{
  os << "BC: " << BC.PhiDiricheltCoef << ", " << BC.PhiNeumannCoef << ", " << BC.PhiRhs << ", " << BC.TDiricheltCoef << ", " << BC.TNeumannCoef << ", " << BC.TRhs << std::endl;
  return os;
}


std::ostream& operator << (std::ostream& os, const std::vector<BoundaryCondition>& BCs)
{
  for (int i=0;i<BCs.size(); i++) os << BCs[i];
  return os;
}*/

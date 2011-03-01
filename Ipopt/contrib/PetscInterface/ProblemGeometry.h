// Copyright (C) 2010 International Business Machines.
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// $Id$
//
// Authors:  Johannes Huber     IBM        2010-09-03

#ifndef PROBLEM_GEOMETRY_INCLUDED
#define PROBLEM_GEOMETRY_INCLUDED
#include "mesh_triangle_support.h"
#include "mesh_tetgen_support.h"
#include "mesh.h"
#include <boundary_info.h>

#include <vector>
#include <list>

#define EXHAUST_AS_CONTROL

void PrintTetgenMesh(const tetgenio& tet, std::ostream& os);
void PrintTriangleMesh(const libMesh::Triangle::triangulateio& tri, std::ostream& os);

class BoundaryConditionBase
{
  //Boundary Conditions:  \alpha\del \Phi/\del n + \beta \Phi = \gamma
public:
  virtual double PhiDirichletCoef(const libMesh::Point& x) = 0; // \beta
  virtual double PhiNeumannCoef(const libMesh::Point& x) = 0;  // \alpha
  virtual double PhiRhs(const libMesh::Point& x) = 0; // \gamma
  virtual double TDirichletCoef(const libMesh::Point& x) = 0;  // these three will be for temparature
  virtual double TNeumannCoef(const libMesh::Point& x) = 0;
  virtual double TRhs(const libMesh::Point& x) = 0;
  virtual bool IsPhiDirichlet() = 0;
  virtual bool IsTDirichlet() = 0;
  BoundaryConditionBase() {}
  virtual ~BoundaryConditionBase() {}
  virtual void print() = 0;
};

class BoundaryConditionConstValues : public BoundaryConditionBase
{
  //Boundary Conditions:  \alpha\del \Phi/\del n + \beta \Phi = \gamma
public:
  double _PhiDirichletCoef; // \beta
  double _PhiNeumannCoef;  // \alpha
  double _PhiRhs; // \gamma
  double _TDirichletCoef;  // these three will be for temparature
  double _TNeumannCoef;
  double _TRhs;

  virtual double PhiDirichletCoef(const libMesh::Point& x) {return _PhiDirichletCoef;} // \beta
  virtual double PhiNeumannCoef(const libMesh::Point& x) {return _PhiNeumannCoef;}  // \alpha
  virtual double PhiRhs(const libMesh::Point& x) { return _PhiRhs;} // \gamma
  virtual double TDirichletCoef(const libMesh::Point& x) {return _TDirichletCoef;}  // these three will be for temparature
  virtual double TNeumannCoef(const libMesh::Point& x) {return _TNeumannCoef;}
  virtual double TRhs(const libMesh::Point& x)  {return _TRhs;}
  virtual bool IsPhiDirichlet() { if (fabs(_PhiDirichletCoef)<1e-16) return false; else return (fabs(_PhiNeumannCoef/_PhiDirichletCoef)<1e-16); }
  virtual bool IsTDirichlet() { if (fabs(_TDirichletCoef)<1e-16) return false; else return (fabs(_TNeumannCoef/_TDirichletCoef)<1e-16); }
  virtual void print() { std::cout << _PhiDirichletCoef << " " << _PhiNeumannCoef << " " << _PhiRhs << " " << _TDirichletCoef << " " << _TNeumannCoef << " " << _TRhs << std::endl; }

  BoundaryConditionConstValues() : _PhiDirichletCoef(0.0), _PhiNeumannCoef(0.0), _PhiRhs(0.0),
      _TDirichletCoef(0.0), _TNeumannCoef(0.0), _TRhs(0.0)
  {}
  BoundaryConditionConstValues(double PhiDirC, double PhiNeumC, double PhiRhsVal, double TDirC, double TNeumC, double TRhsVal) : _PhiDirichletCoef(PhiDirC), _PhiNeumannCoef(PhiNeumC), _PhiRhs(PhiRhsVal), _TDirichletCoef(TDirC), _TNeumannCoef(TNeumC), _TRhs(TRhsVal)
  {}
};

class BoundaryConditionSquarePhiRhs : public BoundaryConditionBase
{
  //Boundary Conditions:  \alpha\del \Phi/\del n + \beta \Phi = \gamma
public:
  double _PhiDirichletCoef; // \beta
  double _PhiNeumannCoef;  // \alpha
  const std::vector<double> _min, _max;
  double _PhiRhsScale; // \gamma
  double _TDirichletCoef;  // these three will be for temparature
  double _TNeumannCoef;
  double _TRhs;

  virtual double PhiDirichletCoef(const libMesh::Point& x) {return _PhiDirichletCoef;} // \beta
  virtual double PhiNeumannCoef(const libMesh::Point& x) {return _PhiNeumannCoef;}  // \alpha
  virtual double PhiRhs(const libMesh::Point& x)
  {
    double val = 0.0;
    double tmp;
    for(int idim=0; idim<_min.size();idim++) {    // in LibMesh
      if(_max[idim]-_min[idim]>1e-12) {
        tmp = (_min[idim]+_max[idim])/2.0;
        val += (x(idim)-_min[idim]) * (_max[idim]-x(idim)) / ((tmp-_min[idim]) * (_max[idim]-tmp));
        // std::cout << _min[idim] << ", " << _max[idim] << ", " << tmp << std::endl;
      }
    }
    return val*_PhiRhsScale;
  }
  virtual double TDirichletCoef(const libMesh::Point& x) {return _TDirichletCoef;}  // these three will be for temparature
  virtual double TNeumannCoef(const libMesh::Point& x) {return _TNeumannCoef;}
  virtual double TRhs(const libMesh::Point& x)  {return _TRhs;}
  virtual bool IsPhiDirichlet() { if (fabs(_PhiDirichletCoef)<1e-16) return false; else return (fabs(_PhiNeumannCoef/_PhiDirichletCoef)<1e-16); }
  virtual bool IsTDirichlet() { if (fabs(_TDirichletCoef)<1e-16) return false; else return (fabs(_TNeumannCoef/_TDirichletCoef)<1e-16); }

  BoundaryConditionSquarePhiRhs() : _PhiDirichletCoef(0.0), _PhiNeumannCoef(0.0), _PhiRhsScale(0.0),
      _TDirichletCoef(0.0), _TNeumannCoef(0.0), _TRhs(0.0)
  {}
  BoundaryConditionSquarePhiRhs(double PhiDirC, double PhiNeumC, double PhiRhsVal,
				double TDirC, double TNeumC, double TRhsVal,
				const std::vector<double>& min, const std::vector<double>& max)
			 : _PhiDirichletCoef(PhiDirC), _PhiNeumannCoef(PhiNeumC), _PhiRhsScale(PhiRhsVal),
			   _TDirichletCoef(TDirC), _TNeumannCoef(TNeumC), _TRhs(TRhsVal), _min(min), _max(max)
  {}
  virtual void print() { std::cout << _PhiDirichletCoef << " " << _PhiNeumannCoef << " " << _min[0] << " " << _max[0] << " " << _min[1] << " " << _max[1] << " " << _min[2] << " " << _max[2] << " " <<  _PhiRhsScale << " " << _TDirichletCoef << " " << _TNeumannCoef << " " << _TRhs << std::endl; }
};

//std::ostream& operator << (std::ostream& os, const BoundaryCondition& BC);
//std::ostream& operator << (std::ostream& os, const std::vector<BoundaryCondition>& BCs);

template<class T>
std::ostream& operator << (std::ostream& os, const std::vector<T>& vec)
{
  for (int i=0;i<vec.size(); i++) os << vec[i] << " ";
  return os;
}


class ProblemGeometry
{
private:
  // class holding data for an item (AC, exhaust, equipment), modelled as cube, perhaps with thickness 0
  class Item
  {
  public:
    std::vector<double> min_, max_; // extreme corners of box (dim: 2 or 3 for 2D/3D)
    std::vector<int> BoundaryMarker;    // [0]:+x0 (from midpoint in positive x0 direction), [1]:-x0, 2:+x1, 3:-x1, 4:+x2, 5:-x2
    // dimension: 4 or 6 (number of sides of box)
    // for AC: depending on which wall
    Item(const std::vector<double>& min, const std::vector<double>& max);

    // Check, if pt is within epsilon of Item
    bool IsInArea(const std::vector<double>& pt, double epsilon=1e-8);
    // which item boundary is clodest to pt?
    // return value: 0:+x0, 1:-x0, 2:+x1, 3:-x1, 4:+x2, 5:-x2 (same as BoundaryMaker index)
    int GetClosestBoundary(const std::vector<double>& pt); //return: which component of member BounaryMarker is closest

    int GetBoundaryMarker(const std::vector<double>& pt); // same as before but returns boundary marker, and -1 is not inside Item

  private:
    // ensures Min<=Max
    void ValidateMinMax();
  };

public:
  // class to represent one control parameter of the problem, which is defined by
  // a bounday marker (where holds the BC depending on this control variable) and
  // a BCParameter (which coefficient is controled)
  class ControlParameter
  {
  public:
    int BoundaryMarker; // same as boundary marker in tetview
    int BCParameter;  // 0:PhiDiriCoef, 1: PhiNeumCoef, 2: PhiRhs, 3:TDiriCoef, 4: TNeumCoef, 5: TiRhs
    ControlParameter(int BM, int BCP) : BoundaryMarker(BM), BCParameter(BCP)
    {}
  };

private:
  std::vector<double> _RoomSize;    // always starts at (0,0)
  std::vector<Item> _Equip;
  std::vector<Item> _AC;
  std::vector<Item> _Exh;
  std::vector<BoundaryConditionBase*> _BoundCond;  // mapping boundary marker -> boundary condition
  std::vector<ControlParameter> _ParamIdx2BCParam;  // mapping control parameter index in optimization problem numbering -> control parameter
  std::vector< std::list<Item*> > _ItemAtWall; // Array of List of Item on each wall (dimension always 4: +x0 -x0 +x1 -x1)
  int NextFreeBoundaryMarker;                 // used to increment boundary markers when building the model (Add...)
  double _h;                                   // initial discretization width (refinement not tracked)
  unsigned short _FE_degree;                  // degree of finite elements (linear=1, quadratic=2)
public:
  ProblemGeometry();
  ~ProblemGeometry();
  inline int GetDim() const
  {
    return _RoomSize.size();
  }
  inline unsigned short GetFE_Degree() const
  {
    return _FE_degree;
  }
  void AddEquipment(const std::vector<double>& min, const std::vector<double>& max, double TempEquip, double kappa); // TEmpEquip and kappa for T variable (later)
  void AddAC(const std::vector<double>& min, const std::vector<double>& max, double vAc, double TempAc); // This vAc will be changed during optimization
  void AddExhaust(const std::vector<double>& min, const std::vector<double>& max);
  void GetHeatExchangeBoundaryMarkers(int iEquip, std::set<int>* pVals);
  void CreateMesh(libMesh::UnstructuredMesh* p_mesh, int order);
  const std::vector<BoundaryConditionBase*>& GetBoundaryConditions() const
  {
    return _BoundCond;
  }
  void ReadFromStream(std::istream& is); // read Problem.dat input file
private:
  void SetBoundaryInfo( libMesh::Mesh* p_mesh );
  void Tetgen2Mesh(const tetgenio& tet, libMesh::UnstructuredMesh* p_mesh, int order);
  int GetWallItemWall(const Item& item);
  void SetFacette(tetgenio::facet *f, int Base, int Offset0, int Offset1, int Offset2, int Offset3, int nHoles=0);
  void SetPoint(REAL* pointlist, int idx, REAL x, REAL y, REAL z);
  void ValidateItemAtWall();
  void Triangle2Mesh(const libMesh::Triangle::triangulateio& tri, libMesh::UnstructuredMesh* p_mesh, int order);
  int GetBoundaryMarker(const std::vector<double>& pt);
  int GetWall(const std::vector<double>& pt);
  void CreateMesh3D(libMesh::UnstructuredMesh* p_mesh, int order);
  void CreateMesh2D(libMesh::UnstructuredMesh* p_mesh, int order);
  void ReadNodeFile(std::string str, tetgenio* tet);
  void ReadEleFile(std::string str, tetgenio* tet);
  void ReadNeighFile(std::string str, tetgenio* tet);

private:
  class CompareItem // Helper class, used to sort items according to their possition on the wall (only required for 2D)
  {
  public:
    Item* pItem;
    double SortValue;
  };
  static bool CompareItems(CompareItem a, CompareItem b)
  {
    return a.SortValue < b.SortValue;
  }

  friend class LibMeshPDEBase;
};
#endif //PROBLEM_GEOMETRY_INCLUDED

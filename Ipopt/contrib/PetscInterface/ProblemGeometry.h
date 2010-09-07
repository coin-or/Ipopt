// Copyright (C) 2010 International Business Machines.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id$
//
// Authors:  Johannes Huber     IBM        2010-09-03

#ifndef PROBLEM_GEOMETRY_INCLUDED
#define PROBLEM_GEOMETRY_INCLUDED
#include "mesh_triangle_support.h"
#include "mesh_tetgen_support.h"
#include "mesh.h"
#include <cell_tet4.h>
#include <face_tri3.h>
#include <boundary_info.h>

#include <vector>
#include <list>

void PrintTetgenMesh(const tetgenio& tet, std::ostream& os);
void PrintTriangleMesh(const libMesh::Triangle::triangulateio& tri, std::ostream& os);

class BoundaryCondition
{
	public:
	double PhiDiricheltCoef;
	double PhiNeumannCoef;
	double PhiRhs;
	double TDiricheltCoef;
	double TNeumannCoef;
	double TRhs;
	BoundaryCondition() : PhiDiricheltCoef(0.0), PhiNeumannCoef(0.0), PhiRhs(0.0),
												TDiricheltCoef(0.0), TNeumannCoef(0.0), TRhs(0.0) {}
	BoundaryCondition(double PhiDirC, double PhiNeumC, double PhiRhsVal, double TDirC, double TNeumC, double TRhsVal) : PhiDiricheltCoef(PhiDirC), PhiNeumannCoef(PhiNeumC), PhiRhs(PhiRhsVal), TDiricheltCoef(TDirC), TNeumannCoef(TNeumC), TRhs(TRhsVal) {}
};

std::ostream& operator << (std::ostream& os, const BoundaryCondition& BC);
//std::ostream& operator << (std::ostream& os, const std::vector<BoundaryCondition>& BCs);

template<class T>
std::ostream& operator << (std::ostream& os, const std::vector<T>& vec)
{ for(int i=0;i<vec.size(); i++) os << vec[i] << " "; return os; }


class ProblemGeometry
{
  private:
    // class holding data for an item (AC, exhaust, equipment), modelled as cube, perhaps with thickness 0
    class Item
    {
      public:
        std::vector<double> min_, max_;
        std::vector<int> BoundaryMarker;    // 0:+x0, 1:-x0, 2:+x1, 3:-x1, 4:+x2, 5:-x2
        Item(const std::vector<double>& min, const std::vector<double>& max);

        // Check, if pt is within epsilon of Item
        bool IsInArea(const std::vector<double>& pt, double epsilon=1e-8);
        // which item boundary is clodest to pt?
        // return value: 0:+x0, 1:-x0, 2:+x1, 3:-x1, 4:+x2, 5:-x2 (same as BoundaryMaker index)
        int GetClosestBoundary(const std::vector<double>& pt);

        int GetBoundaryMarker(const std::vector<double>& pt);

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
      int BoundaryMarker;
      int BCParameter;  // 0:PhiDiriCoef, 1: PhiNeumCoef, 2: PhiRhs, 3:TDiriCoef, 4: TNeumCoef, 5: TiRhs
      ControlParameter(int BM, int BCP) : BoundaryMarker(BM), BCParameter(BCP) {}
    };

  private:
    std::vector<double> _RoomSize;    // always starts at (0,0)
    std::vector<Item> _Equip;
    std::vector<Item> _AC;
    std::vector<Item> _Exh;
    std::vector<BoundaryCondition> _BoundCond;  // mapping boundary marker -> boundary condition
    std::vector<ControlParameter> _ParamIdx2BCParam;  // mapping control parameter index -> control parameter
    std::vector< std::list<Item*> > _ItemAtWall;	// Array of List of Item on each wall
    int NextFreeBoundaryMarker;                 // used to increment boundary markers when building the model (Add...)
    double _h;                                   // initial discretization width (refinement not tracked) 
  public:
    ProblemGeometry();
    inline int GetDim() const {return _RoomSize.size();}
    void AddEquipment(const std::vector<double>& min, const std::vector<double>& max, double TempEquip, double kappa);
    void AddAC(const std::vector<double>& min, const std::vector<double>& max, double vAc, double TempAc);
    void AddExhaust(const std::vector<double>& min, const std::vector<double>& max);

	  void CreateMesh(libMesh::UnstructuredMesh* p_mesh);
    const std::vector<BoundaryCondition>& GetBoundaryConditions() const { return _BoundCond; }
    void ReadFromStream(std::istream& is);
  private:
    void SetBoundaryInfo( libMesh::Mesh* p_mesh );
	  void Tetgen2Mesh(const tetgenio& tet, libMesh::UnstructuredMesh* p_mesh);
    int GetWallItemWall(const Item& item);
	  void SetFacette(tetgenio::facet *f, int Base, int Offset0, int Offset1, int Offset2, int Offset3, int nHoles=0);
	  void SetPoint(REAL* pointlist, int idx, REAL x, REAL y, REAL z);
	  void ValidateItemAtWall();
	  void Triangle2Mesh(const libMesh::Triangle::triangulateio& tri, libMesh::UnstructuredMesh* p_mesh);
	  int GetBoundaryMarker(const std::vector<double>& pt);
    int GetWall(const std::vector<double>& pt);
    void CreateMesh3D(libMesh::UnstructuredMesh* p_mesh);
    void CreateMesh2D(libMesh::UnstructuredMesh* p_mesh);
    void ReadNodeFile(std::string str, tetgenio* tet);
    void ReadEleFile(std::string str, tetgenio* tet);
    void ReadNeighFile(std::string str, tetgenio* tet);

  private:
    class CompareItem // Helper class, used to sort items according to their possition on the wall
    {
      public:
      Item* pItem;
      double SortValue;
    };
    static bool CompareItems(CompareItem a, CompareItem b) { return a.SortValue < b.SortValue; }

  friend class LibMeshPDEBase;
};
#endif //PROBLEM_GEOMETRY_INCLUDED

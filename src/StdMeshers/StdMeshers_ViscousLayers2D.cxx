// Copyright (C) 2007-2026  CEA, EDF, OPEN CASCADE
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
//
// See http://www.salome-platform.org/ or email : webmaster.salome@opencascade.com
//

// File      : StdMeshers_ViscousLayers2D.cxx
// Created   : 23 Jul 2012
// Author    : Edward AGAPOV (eap)

#include "StdMeshers_ViscousLayers2D.hxx"

#include "SMDS_EdgePosition.hxx"
#include "SMDS_FaceOfNodes.hxx"
#include "SMDS_FacePosition.hxx"
#include "SMDS_MeshNode.hxx"
#include "SMDS_SetIterator.hxx"
#include "SMESHDS_Group.hxx"
#include "SMESHDS_Hypothesis.hxx"
#include "SMESHDS_Mesh.hxx"
#include "SMESH_Algo.hxx"
#include "SMESH_ComputeError.hxx"
#include "SMESH_ControlsDef.hxx"
#include "SMESH_Gen.hxx"
#include "SMESH_Group.hxx"
#include "SMESH_HypoFilter.hxx"
#include "SMESH_Mesh.hxx"
#include "SMESH_MeshEditor.hxx"
#include "SMESH_MesherHelper.hxx"
#include "SMESH_ProxyMesh.hxx"
#include "SMESH_Quadtree.hxx"
#include "SMESH_subMesh.hxx"
#include "SMESH_subMeshEventListener.hxx"
#include "StdMeshers_FaceSide.hxx"

#include "utilities.h"

#include <BRepAdaptor_Curve.hxx>
#include <BRepAdaptor_Curve2d.hxx>
#include <BRepTools.hxx>
#include <BRep_Tool.hxx>
#include <Bnd_B2d.hxx>
#include <Bnd_B3d.hxx>
#include <ElCLib.hxx>
#include <GCPnts_AbscissaPoint.hxx>
#include <Geom2dAdaptor_Curve.hxx>
#include <Geom2dInt_GInter.hxx>
#include <Geom2d_Circle.hxx>
#include <Geom2d_Line.hxx>
#include <Geom2d_TrimmedCurve.hxx>
#include <GeomAdaptor_Curve.hxx>
#include <Geom_Circle.hxx>
#include <Geom_Curve.hxx>
#include <Geom_Line.hxx>
#include <Geom_Surface.hxx>
#include <Geom_TrimmedCurve.hxx>
#include <IntRes2d_IntersectionPoint.hxx>
#include <Precision.hxx>
#include <Standard_ErrorHandler.hxx>
#include <TColStd_Array1OfReal.hxx>
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <TopTools_IndexedDataMapOfShapeListOfShape.hxx>
#include <TopTools_IndexedMapOfShape.hxx>
#include <TopTools_ListIteratorOfListOfShape.hxx>
#include <TopTools_ListOfShape.hxx>
#include <TopTools_MapOfShape.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Vertex.hxx>
#include <gp_Ax1.hxx>
#include <gp_Dir.hxx>
#include <gp_Pnt.hxx>
#include <gp_Pnt2d.hxx>
#include <gp_Trsf.hxx>
#include <gp_Vec.hxx>
#include <gp_XY.hxx>
#include <smIdType.hxx>

#include <list>
#include <string>
#include <cmath>
#include <limits>

#ifdef _DEBUG_
//#define __myDEBUG
#endif

using namespace std;

//================================================================================
namespace VISCOUS_2D
{
  typedef int TGeomID;

  //--------------------------------------------------------------------------------
  /*!
   * \brief Proxy Mesh of FACE with viscous layers. It's needed only to 
   *        redefine newSubmesh().
   */
  struct _ProxyMeshOfFace : public SMESH_ProxyMesh
  {
    //---------------------------------------------------
    // Proxy sub-mesh of an EDGE. It contains nodes in _uvPtStructVec.
    struct _EdgeSubMesh : public SMESH_ProxyMesh::SubMesh
    {
      _EdgeSubMesh(const SMDS_Mesh* mesh, int index=0): SubMesh(mesh,index) {}
      //virtual int NbElements() const { return _elements.size()+1; }
      virtual smIdType NbNodes() const { return Max( 0, _uvPtStructVec.size()-2 ); }
      void SetUVPtStructVec(UVPtStructVec& vec) { _uvPtStructVec.swap( vec ); }
      UVPtStructVec& GetUVPtStructVec() { return _uvPtStructVec; }
    };
    _ProxyMeshOfFace(const SMESH_Mesh& mesh): SMESH_ProxyMesh(mesh) {}
    _EdgeSubMesh* GetEdgeSubMesh(int ID) { return (_EdgeSubMesh*) getProxySubMesh(ID); }
    virtual SubMesh* newSubmesh(int index=0) const { return new _EdgeSubMesh( GetMeshDS(), index); }
  };
  //--------------------------------------------------------------------------------
  /*!
   * \brief SMESH_subMeshEventListener used to store _ProxyMeshOfFace, computed
   *        by _ViscousBuilder2D, in a SMESH_subMesh of the FACE.
   *        This is to delete _ProxyMeshOfFace when StdMeshers_ViscousLayers2D
   *        hypothesis is modified
   */
  struct _ProxyMeshHolder : public SMESH_subMeshEventListener
  {
    _ProxyMeshHolder( const TopoDS_Face&    face,
                      SMESH_ProxyMesh::Ptr& mesh)
      : SMESH_subMeshEventListener( /*deletable=*/true, Name() )
    {
      SMESH_subMesh* faceSM = mesh->GetMesh()->GetSubMesh( face );
      faceSM->SetEventListener( this, new _Data( mesh ), faceSM );
    }
    // Finds a proxy mesh of face
    static SMESH_ProxyMesh::Ptr FindProxyMeshOfFace( const TopoDS_Shape& face,
                                                     SMESH_Mesh&         mesh )
    {
      SMESH_ProxyMesh::Ptr proxy;
      SMESH_subMesh* faceSM = mesh.GetSubMesh( face );
      if ( EventListenerData* ld = faceSM->GetEventListenerData( Name() ))
        proxy = static_cast< _Data* >( ld )->_mesh;
      return proxy;
    }
    // Treat events
    void ProcessEvent(const int          event,
                      const int          eventType,
                      SMESH_subMesh*     /*subMesh*/,
                      EventListenerData* data,
                      const SMESH_Hypothesis*  /*hyp*/)
    {
      if ( event == SMESH_subMesh::CLEAN && eventType == SMESH_subMesh::COMPUTE_EVENT)
        ((_Data*) data)->_mesh.reset();
    }
  private:
    // holder of a proxy mesh
    struct _Data : public SMESH_subMeshEventListenerData
    {
      SMESH_ProxyMesh::Ptr _mesh;
      _Data( SMESH_ProxyMesh::Ptr& mesh )
        :SMESH_subMeshEventListenerData( /*isDeletable=*/true), _mesh( mesh )
      {}
    };
    // Returns identifier string
    static const char* Name() { return "VISCOUS_2D::_ProxyMeshHolder"; }
  };
  //--------------------------------------------------------------------------------
  /*!
   * \brief Data shared between the two "arm" FACEs (each carrying viscous layers
   *        on an EDGE meeting at a common VERTEX through an intervening bare
   *        FACE) needed to bevel that VERTEX along their geometric bisector B13
   *        instead of both shrinking all the way to the VERTEX. Created once per
   *        VERTEX by whichever arm's _ViscousBuilder2D::shrink() reaches it
   *        first; the other arm and the bare FACE's Compute() reuse the same
   *        cached P1Node/P2Node so the mesh stays conformal across all three
   *        FACEs. See StdMeshers_ViscousLayers2D::ComputeCornerBisector().
   */
  struct _CornerData
  {
    bool                  _isValid;    // true once B13/P1/P2 are known
    gp_Dir                _B13;
    gp_Pnt                _p1Point;    // target position owned by _armFaceID1
    gp_Pnt                _p2Point;    // target position owned by _armFaceID3
    TGeomID               _armFaceID1; // ShapeToIndex of the FACE that computed
                                       // _CornerData and owns _p1Point/_p1Node
    TGeomID               _armFaceID3; // ShapeToIndex of the other arm FACE,
                                       // owning _p2Point/_p2Node
    const SMDS_MeshNode*  _p1Node;
    const SMDS_MeshNode*  _p2Node;
    smIdType              _p1NodeID;   // GetID() of _p1Node/_p2Node at the moment
    smIdType              _p2NodeID;   // they were cached -- see FindCornerData()'s
                                       // staleness check: a Mesh.Clear() deletes all
                                       // nodes without necessarily notifying this
                                       // VERTEX-keyed cache, which would otherwise
                                       // hand back dangling pointers on a later re-Compute
    bool                  _wedgeBuilt; // true once the (O,P1,P2) closing triangle
                                       // has been created by its designated owner

    // Per-arm chain of one bend node per layer (index iF = layer number, same
    // order/count as StdMeshers_ViscousLayers2D::GetNumberLayers()), populated
    // once by buildCornerBevel() purely so the OTHER arm (running second, when
    // thicknesses/layer counts match and bend points coincide spatially) can
    // reuse the same node instead of creating a spatially-duplicate one at the
    // same position -- lookup only, never fed into refine()'s own per-layer
    // bands (confirmed unsafe, see buildCornerBevel()'s own comments).
    vector< const SMDS_MeshNode* > _bendChain1; // owned by _armFaceID1
    vector< const SMDS_MeshNode* > _bendChain2; // owned by _armFaceID3

    // Per-arm chain of one "near" node per layer (Qi = V + height*N, N = that
    // arm's OWN inward growth normal at V -- see buildCornerBevel()), i.e. the
    // arm's own official per-layer band end AT V. Unlike _bendChain1/2 above,
    // THIS chain IS fed into refine()'s own per-layer bands (isCornerBevelL/R),
    // which is only safe together with splicing the same nodes into this arm's
    // OWN proxy submesh for its bare-adjacent EDGE (see the
    // spliceCornerIntoEdgeProxy call in shrink()) -- without that splice, the
    // underlying 2D algorithm's core fill still reaches V and any real element
    // built on this chain would overlap it (see the isCornerBevelL/R comment
    // in refine() for the full story of why a naive version of this failed).
    vector< const SMDS_MeshNode* > _nearChain1; // owned by _armFaceID1
    vector< const SMDS_MeshNode* > _nearChain2; // owned by _armFaceID3

    _CornerData()
      : _isValid(false), _armFaceID1(0), _armFaceID3(0),
        _p1Node(0), _p2Node(0), _p1NodeID(-1), _p2NodeID(-1), _wedgeBuilt(false) {}
  };
  //--------------------------------------------------------------------------------
  /*!
   * \brief SMESH_subMeshEventListener caching _CornerData for a VERTEX shared by
   *        two independently viscous-layered FACEs through an intervening bare
   *        FACE. Modeled on _ProxyMeshHolder but keyed off the VERTEX's own
   *        sub-mesh rather than a FACE's, since B13/P1/P2 belong to neither arm
   *        FACE alone. Invalidated on CLEAN of the VERTEX's sub-mesh; callers
   *        are responsible for making sure that sub-mesh's CLEAN is cascaded
   *        from either arm FACE (see the VISCOUS_3D::ToClearSubWithMain call in
   *        _ViscousBuilder2D::shrink() and in RestoreListeners()).
   */
  struct _CornerDataHolder : public SMESH_subMeshEventListener
  {
    _CornerDataHolder( const TopoDS_Vertex& vertex, SMESH_Mesh& mesh )
      : SMESH_subMeshEventListener( /*deletable=*/true, Name() )
    {
      SMESH_subMesh* vertexSM = mesh.GetSubMesh( vertex );
      vertexSM->SetEventListener( this, new _Data(), vertexSM );
    }
    // Finds (without creating) the _CornerData cached for a VERTEX, or NULL.
    // Self-heals a stale cache: a Mesh.Clear() (e.g. a user-triggered
    // Clear()+Compute() cycle) deletes every SMDS_MeshNode, but does not
    // reliably fire a CLEAN event on THIS VERTEX's own sub-mesh (unlike a
    // FACE's sub-mesh, whose CLEAN event _ProxyMeshHolder above relies on) --
    // confirmed by a SIGSEGV reproduction: a second Compute() dereferenced
    // _p1Node/_bendChain pointers left over from the first one, now dangling.
    // Since dereferencing a possibly-freed pointer to check it is itself
    // unsafe, validate by identity through the mesh's own node registry
    // instead: if the node currently registered under the ID we cached no
    // longer IS that same pointer, the mesh was cleared/rebuilt since, and
    // the whole cache (not just that one node) is discarded and rebuilt.
    static _CornerData* FindCornerData( const TopoDS_Vertex& vertex, SMESH_Mesh& mesh )
    {
      SMESH_subMesh* vertexSM = mesh.GetSubMesh( vertex );
      if ( EventListenerData* ld = vertexSM->GetEventListenerData( Name() ))
      {
        _CornerData* cd = & static_cast< _Data* >( ld )->_corner;
        if ( cd->_isValid &&
             (( cd->_p1Node && mesh.GetMeshDS()->FindNode( cd->_p1NodeID ) != cd->_p1Node ) ||
              ( cd->_p2Node && mesh.GetMeshDS()->FindNode( cd->_p2NodeID ) != cd->_p2Node )))
          *cd = _CornerData();
        return cd;
      }
      return 0;
    }
    // Finds the _CornerData cached for a VERTEX, creating an empty (invalid)
    // one and registering this listener if none exists yet.
    static _CornerData* GetCornerData( const TopoDS_Vertex& vertex, SMESH_Mesh& mesh )
    {
      if ( _CornerData* cd = FindCornerData( vertex, mesh ))
        return cd;
      new _CornerDataHolder( vertex, mesh );
      return FindCornerData( vertex, mesh );
    }
    // Treat events
    void ProcessEvent(const int          event,
                      const int          eventType,
                      SMESH_subMesh*     /*subMesh*/,
                      EventListenerData* data,
                      const SMESH_Hypothesis*  /*hyp*/)
    {
      if ( event == SMESH_subMesh::CLEAN && eventType == SMESH_subMesh::COMPUTE_EVENT)
        static_cast< _Data* >( data )->_corner = _CornerData();
    }
  private:
    // holder of a _CornerData
    struct _Data : public SMESH_subMeshEventListenerData
    {
      _CornerData _corner;
      _Data(): SMESH_subMeshEventListenerData( /*isDeletable=*/true ) {}
    };
    // Returns identifier string
    static const char* Name() { return "VISCOUS_2D::_CornerDataHolder"; }
  };

  struct _PolyLine;
  //--------------------------------------------------------------------------------
  /*!
   * \brief Segment connecting inner ends of two _LayerEdge's.
   */
  struct _Segment
  {
    const gp_XY* _uv[2];       // pointer to _LayerEdge::_uvIn
    int          _indexInLine; // position in _PolyLine

    _Segment() {}
    _Segment(const gp_XY& p1, const gp_XY& p2):_indexInLine(-1) { _uv[0] = &p1; _uv[1] = &p2; }
    const gp_XY& p1() const { return *_uv[0]; }
    const gp_XY& p2() const { return *_uv[1]; }
  };
  //--------------------------------------------------------------------------------
  /*!
   * \brief Tree of _Segment's used for a faster search of _Segment's.
   */
  struct _SegmentTree : public SMESH_Quadtree
  {
    typedef boost::shared_ptr< _SegmentTree > Ptr;

    _SegmentTree( const vector< _Segment >& segments );
    void GetSegmentsNear( const _Segment& seg, vector< const _Segment* >& found );
    void GetSegmentsNear( const gp_Ax2d& ray, vector< const _Segment* >& found );
  protected:
    _SegmentTree() {}
    _SegmentTree* newChild() const { return new _SegmentTree; }
    void          buildChildrenData();
    Bnd_B2d*      buildRootBox();
  private:
    static int    maxNbSegInLeaf() { return 5; }
    struct _SegBox
    {
      const _Segment* _seg;
      bool            _iMin[2];
      void Set( const _Segment& seg )
      {
        _seg = &seg;
        _iMin[0] = ( seg._uv[1]->X() < seg._uv[0]->X() );
        _iMin[1] = ( seg._uv[1]->Y() < seg._uv[0]->Y() );
      }
      bool IsOut( const _Segment& seg ) const;
      bool IsOut( const gp_Ax2d& ray ) const;
    };
    vector< _SegBox > _segments;
  };
  //--------------------------------------------------------------------------------
  /*!
   * \brief Edge normal to FACE boundary, connecting a point on EDGE (_uvOut)
   * and a point of a layer internal boundary (_uvIn)
   */
  struct _LayerEdge
  {
    gp_XY         _uvOut;    // UV on the FACE boundary
    gp_XY         _uvIn;     // UV inside the FACE
    double        _length2D; // distance between _uvOut and _uvIn

    bool          _isBlocked;// is more inflation possible or not

    gp_XY         _normal2D; // to curve
    double        _len2dTo3dRatio; // to pass 2D <--> 3D
    gp_Ax2d       _ray;      // a ray starting at _uvOut

    vector<gp_XY> _uvRefined; // divisions by layers

    bool SetNewLength( const double length );

    int           _ID; // debug
  };
  //--------------------------------------------------------------------------------
  /*!
   * \brief Poly line composed of _Segment's of one EDGE.
   *        It's used to detect intersection of inflated layers by intersecting
   *        _Segment's in 2D.
   */
  struct _PolyLine
  {
    StdMeshers_FaceSide* _wire;
    int                  _edgeInd;     // index of my EDGE in _wire
    bool                 _advancable;  // true if there is a viscous layer on my EDGE
    bool                 _isStraight2D;// pcurve type
    _PolyLine*           _leftLine;    // lines of neighbour EDGE's
    _PolyLine*           _rightLine;
    int                  _firstPntInd; // index in vector<UVPtStruct> of _wire
    int                  _lastPntInd;
    int                  _index;       // index in _ViscousBuilder2D::_polyLineVec

    vector< _LayerEdge > _lEdges;      /* _lEdges[0] is usually is not treated
                                          as it is equal to the last one of the _leftLine */
    vector< _Segment >   _segments;    // segments connecting _uvIn's of _lEdges
    _SegmentTree::Ptr    _segTree;

    vector< _PolyLine* > _reachableLines;       // lines able to interfere with my layer

    vector< const SMDS_MeshNode* > _leftNodes;  // nodes built from a left VERTEX
    vector< const SMDS_MeshNode* > _rightNodes; // nodes built from a right VERTEX

    typedef vector< _Segment >::iterator   TSegIterator;
    typedef vector< _LayerEdge >::iterator TEdgeIterator;

    TIDSortedElemSet     _newFaces; // faces generated from this line

    bool IsCommonEdgeShared( const _PolyLine& other );
    size_t FirstLEdge() const
    {
      return ( _leftLine->_advancable && _lEdges.size() > 2 ) ? 1 : 0;
    }
    bool IsAdjacent( const _Segment& seg, const _LayerEdge* LE=0 ) const
    {
      if ( LE /*&& seg._indexInLine < _lEdges.size()*/ )
        return ( seg._uv[0] == & LE->_uvIn ||
                 seg._uv[1] == & LE->_uvIn );
      return ( & seg == &_leftLine->_segments.back() ||
               & seg == &_rightLine->_segments[0] );
    }
    bool IsConcave() const;
  };
  //--------------------------------------------------------------------------------
  /*!
   * \brief Intersector of _Segment's
   */
  struct _SegmentIntersection
  {
    gp_XY    _vec1, _vec2;     // Vec( _seg.p1(), _seg.p2() )
    gp_XY    _vec21;           // Vec( _seg2.p1(), _seg1.p1() )
    double   _D;               // _vec1.Crossed( _vec2 )
    double   _param1, _param2; // intersection param on _seg1 and _seg2

    _SegmentIntersection(): _D(0), _param1(0), _param2(0) {}

    bool Compute(const _Segment& seg1, const _Segment& seg2, bool seg2IsRay = false )
    {
      // !!! If seg2IsRay, returns true at any _param2 !!!
      const double eps = 1e-10;
      _vec1  = seg1.p2() - seg1.p1(); 
      _vec2  = seg2.p2() - seg2.p1(); 
      _vec21 = seg1.p1() - seg2.p1(); 
      _D = _vec1.Crossed(_vec2);
      if ( fabs(_D) < std::numeric_limits<double>::min())
        return false;
      _param1 = _vec2.Crossed(_vec21) / _D; 
      if (_param1 < -eps || _param1 > 1 + eps )
        return false;
      _param2 = _vec1.Crossed(_vec21) / _D;
      return seg2IsRay || ( _param2 > -eps && _param2 < 1 + eps );
    }
    bool Compute( const _Segment& seg1, const gp_Ax2d& ray )
    {
      gp_XY segEnd = ray.Location().XY() + ray.Direction().XY();
      _Segment seg2( ray.Location().XY(), segEnd );
      return Compute( seg1, seg2, true );
    }
    //gp_XY GetPoint() { return _seg1.p1() + _param1 * _vec1; }
  };
  //--------------------------------------------------------------------------------

  typedef map< const SMDS_MeshNode*, _LayerEdge*, TIDCompare > TNode2Edge;
  typedef StdMeshers_ViscousLayers2D                           THypVL;
  
  //--------------------------------------------------------------------------------
  /*!
   * \brief Builder of viscous layers
   */
  class _ViscousBuilder2D
  {
  public:
    _ViscousBuilder2D(SMESH_Mesh&                       theMesh,
                      const TopoDS_Face&                theFace,
                      vector< const THypVL* > &         theHyp,
                      vector< TopoDS_Shape > &          theHypShapes);
    SMESH_ComputeErrorPtr GetError() const { return _error; }
    // does it's job
    SMESH_ProxyMesh::Ptr  Compute();

  private:

    friend class ::StdMeshers_ViscousLayers2D;

    bool findEdgesWithLayers();
    bool makePolyLines();
    bool inflate();
    bool fixCollisions();
    bool refine();
    bool shrink();
    bool improve();
    bool toShrinkForAdjacent( const TopoDS_Face& adjFace,
                              const TopoDS_Edge& E,
                              const TopoDS_Vertex& V);
    bool findBareCornerPattern( const TopoDS_Vertex& V,
                               const TopoDS_Edge&   myArmEdge,
                               const TopoDS_Edge&   myBareEdge,
                               const TopoDS_Face&   bareFace,
                               TopoDS_Edge&         otherArmEdge,
                               TopoDS_Face&         otherArmFace,
                               double&              otherArmThickness,
                               int&                 otherArmNbLayers );
    const SMDS_MeshNode* buildCornerBevel( const TopoDS_Vertex&  V,
                           const TopoDS_Edge&                 myArmEdge,
                           const TopoDS_Edge&                 myBareEdge,
                           const TopoDS_Face&                 bareFace,
                           const SMDS_MeshNode*                vertexNode,
                           const THypVL*                       myHyp,
                           const SMDS_MeshNode*&                outNearCornerNode );
    bool tryBuildBareCornerProxy();
    bool findArmAcrossEdge( const TopoDS_Vertex& V,
                            const TopoDS_Edge&   bareEdge,
                            const TopoDS_Face&   excludeFace,
                            TopoDS_Edge&         armEdge,
                            TopoDS_Face&         armFace,
                            int&                 armNbLayers );
    void spliceCornerIntoEdgeProxy( StdMeshers_FaceSidePtr     wire,
                                    int                        iE,
                                    const vector<UVPtStruct>&  points,
                                    bool                       cornerAtEnd,
                                    const SMDS_MeshNode*       cornerNode,
                                    int                        numWallNodesToExclude );
    void setLenRatio( _LayerEdge& LE, const gp_Pnt& pOut );
    void setLayerEdgeData( _LayerEdge&                 lEdge,
                           const double                u,
                           Handle(Geom2d_Curve)&       pcurve,
                           Handle(Geom_Curve)&         curve,
                           const gp_Pnt                pOut,
                           const bool                  reverse,
                           GeomAPI_ProjectPointOnSurf* faceProj);
    void adjustCommonEdge( _PolyLine& LL, _PolyLine& LR );
    void calcLayersHeight(const double    totalThick,
                          vector<double>& heights,
                          const THypVL*   hyp);
    bool removeMeshFaces(const TopoDS_Shape& face);

    const THypVL*     getLineHypothesis(int iPL);
    double            getLineThickness (int iPL);

    bool              error( const string& text );
    SMESHDS_Mesh*     getMeshDS() { return _mesh->GetMeshDS(); }
    _ProxyMeshOfFace* getProxyMesh();

    // debug
    //void makeGroupOfLE();

  private:

    // input data
    SMESH_Mesh*                 _mesh;
    TopoDS_Face                 _face;
    vector< const THypVL* >     _hyps;
    vector< TopoDS_Shape >      _hypShapes;

    // result data
    SMESH_ProxyMesh::Ptr        _proxyMesh;
    SMESH_ComputeErrorPtr       _error;

    // working data
    Handle(Geom_Surface)        _surface;
    SMESH_MesherHelper          _helper;
    TSideVector                 _faceSideVec; // wires (StdMeshers_FaceSide) of _face
    vector<_PolyLine>           _polyLineVec; // fronts to advance
    vector< const THypVL* >     _hypOfEdge; // a hyp per an EDGE of _faceSideVec
    bool                        _is2DIsotropic; // is same U and V resoulution of _face
    vector<TopoDS_Face>         _clearedFaces; // FACEs whose mesh was removed by shrink()

    //double                      _fPowN; // to compute thickness of layers
    double                      _maxThickness; // max possible layers thickness

    // sub-shapes of _face 
    set<TGeomID>                _ignoreShapeIds; // ids of EDGEs w/o layers
    set<TGeomID>                _noShrinkVert;   // ids of VERTEXes that are extremities
    // of EDGEs along which _LayerEdge can't be inflated because no viscous layers
    // defined on neighbour FACEs sharing an EDGE. Nonetheless _LayerEdge's
    // are inflated along such EDGEs but then such _LayerEdge's are turned into
    // a node on VERTEX, i.e. all nodes on a _LayerEdge are melded into one node.


    int                         _nbLE; // for DEBUG
  };

  //================================================================================
  /*!
   * \brief Returns StdMeshers_ViscousLayers2D for the FACE
   */
  bool findHyps(SMESH_Mesh&                                   theMesh,
                const TopoDS_Face&                            theFace,
                vector< const StdMeshers_ViscousLayers2D* > & theHyps,
                vector< TopoDS_Shape > &                      theAssignedTo)
  {
    theHyps.clear();
    theAssignedTo.clear();
    SMESH_HypoFilter hypFilter
      ( SMESH_HypoFilter::HasName( StdMeshers_ViscousLayers2D::GetHypType() ));
    list< const SMESHDS_Hypothesis * > hypList;
    list< TopoDS_Shape >               hypShapes;
    int nbHyps = theMesh.GetHypotheses
      ( theFace, hypFilter, hypList, /*ancestors=*/true, &hypShapes );
    if ( nbHyps )
    {
      theHyps.reserve( nbHyps );
      theAssignedTo.reserve( nbHyps );
      list< const SMESHDS_Hypothesis * >::iterator hyp = hypList.begin();
      list< TopoDS_Shape >::iterator               shape = hypShapes.begin();
      for ( ; hyp != hypList.end(); ++hyp, ++shape )
      {
        theHyps.push_back( static_cast< const StdMeshers_ViscousLayers2D* > ( *hyp ));
        theAssignedTo.push_back( *shape );
      }
    }
    return nbHyps;
  }

  //================================================================================
  /*!
   * \brief Returns ids of EDGEs not to create Viscous Layers on
   *  \param [in] theHyp - the hypothesis, holding edges either to ignore or not to.
   *  \param [in] theFace - the FACE whose EDGEs are checked.
   *  \param [in] theMesh - the mesh.
   *  \param [in,out] theEdgeIds - container returning EDGEs to ignore.
   *  \return int - number of found EDGEs of the FACE.
   */
  //================================================================================

  int getEdgesToIgnore( const StdMeshers_ViscousLayers2D* theHyp,
                        const TopoDS_Shape&               theFace,
                        const SMESHDS_Mesh*               theMesh,
                        set< int > &                      theEdgeIds)
  {
    int nbEdgesToIgnore = 0;
    vector<TGeomID> ids = theHyp->GetBndShapes();
    if ( theHyp->IsToIgnoreShapes() ) // EDGEs to ignore are given
    {
      for ( size_t i = 0; i < ids.size(); ++i )
      {
        const TopoDS_Shape& E = theMesh->IndexToShape( ids[i] );
        if ( !E.IsNull() &&
             E.ShapeType() == TopAbs_EDGE &&
             SMESH_MesherHelper::IsSubShape( E, theFace ))
        {
          theEdgeIds.insert( ids[i] );
          ++nbEdgesToIgnore;
        }
      }
    }
    else // EDGEs to make the Viscous Layers on are given
    {
      TopExp_Explorer E( theFace, TopAbs_EDGE );
      for ( ; E.More(); E.Next(), ++nbEdgesToIgnore )
        theEdgeIds.insert( theMesh->ShapeToIndex( E.Current() ));

      for ( size_t i = 0; i < ids.size(); ++i )
        nbEdgesToIgnore -= theEdgeIds.erase( ids[i] );
    }
    return nbEdgesToIgnore;
  }

  //================================================================================
  /*!
   * \brief Returns true if theFace has at least one of its OWN EDGEs actually
   *        getting viscous layers. Unlike a plain findHyps() truthy check,
   *        this is correct even when the ViscousLayers2D hypothesis is
   *        assigned once to an ANCESTOR shape (e.g. the whole SHELL/mesh main
   *        shape) rather than per-FACE: findHyps() (ancestors=true) then
   *        finds the SAME hypothesis for every FACE in the shell, including
   *        ones none of whose EDGEs are actually in the hypothesis's edge
   *        list -- those must still be treated as "bare" for the reentrant-
   *        corner bevel pattern.
   */
  //================================================================================

  bool hasOwnViscousLayers( SMESH_Mesh& theMesh, const TopoDS_Face& theFace )
  {
    vector< const StdMeshers_ViscousLayers2D* > hyps;
    vector< TopoDS_Shape >                      hypShapes;
    if ( !findHyps( theMesh, theFace, hyps, hypShapes ))
      return false;

    int totalNbEdges = 0;
    for ( TopExp_Explorer e( theFace, TopAbs_EDGE ); e.More(); e.Next() )
      ++totalNbEdges;

    set<int> ignoreIds;
    for ( size_t i = 0; i < hyps.size(); ++i )
      getEdgesToIgnore( hyps[i], theFace, theMesh.GetMeshDS(), ignoreIds );

    return (int) ignoreIds.size() < totalNbEdges;
  }

  //================================================================================
  /*!
   * \brief Compute, at VERTEX V, the 3D tangent along EDGE E and the inward
   *        growth-normal within FACE F (perpendicular to E, in F's local tangent
   *        plane, pointing into F's interior) -- the 3D analog of the 2D
   *        pcurve-tangent-rotated-by-90-degrees used elsewhere in this file
   *        (see setLayerEdgeData()), but not tied to any single FACE's own UV
   *        parametrization, so it can be combined across two different FACEs.
   *        The inward orientation is verified against the FACE's own parametric
   *        midpoint rather than relying on TopAbs_Orientation bookkeeping, to
   *        stay correct regardless of how E/F's orientations were obtained.
   */
  //================================================================================

  bool edgeTangentAndNormalAt( const TopoDS_Vertex& V,
                               const TopoDS_Edge&   E,
                               const TopoDS_Face&   F,
                               gp_Pnt&              P,
                               gp_Vec&              tangent,
                               gp_Vec&              normal,
                               gp_Vec*              faceNormOut = 0 )
  {
    TopLoc_Location eLoc;
    double f,l;
    Handle(Geom_Curve) curve = BRep_Tool::Curve( E, eLoc, f, l );
    if ( curve.IsNull() )
      return false;
    const double u = BRep_Tool::Parameter( V, E );

    gp_Vec d1;
    curve->D1( u, P, d1 );
    if ( !eLoc.IsIdentity() )
    {
      const gp_Trsf& trsf = eLoc.Transformation();
      P.Transform( trsf );
      d1.Transform( trsf );
    }
    if ( d1.Magnitude() < gp::Resolution() )
      return false;
    tangent = d1.Normalized();

    TopLoc_Location fLoc;
    Handle(Geom_Surface) surf = BRep_Tool::Surface( F, fLoc );
    if ( surf.IsNull() )
      return false;
    const gp_Pnt2d uv = BRep_Tool::Parameters( V, F );

    gp_Pnt ps; gp_Vec du, dv;
    surf->D1( uv.X(), uv.Y(), ps, du, dv );
    if ( !fLoc.IsIdentity() )
    {
      const gp_Trsf& trsf = fLoc.Transformation();
      du.Transform( trsf );
      dv.Transform( trsf );
    }
    const gp_Vec faceNorm = du ^ dv;
    normal = faceNorm ^ tangent;
    if ( normal.Magnitude() < gp::Resolution() )
      return false;
    normal.Normalize();

    // Orient normal to point toward the FACE interior. This is checked
    // explicitly (rather than derived from TopAbs_Orientation composition)
    // to stay robust regardless of how E/F were obtained by the caller.
    double umin, umax, vmin, vmax;
    BRepTools::UVBounds( F, umin, umax, vmin, vmax );
    gp_Pnt center = surf->Value( 0.5 * ( umin + umax ), 0.5 * ( vmin + vmax ));
    if ( !fLoc.IsIdentity() )
      center.Transform( fLoc.Transformation() );
    if ( normal.Dot( gp_Vec( P, center )) < 0 )
      normal.Reverse();

    if ( faceNormOut )
      *faceNormOut = faceNorm;

    return true;
  }

  //================================================================================
  /*!
   * \brief Find the point on the ray { s*B13 : s >= 0 } from O where it crosses
   *        the offset line { pv + t*T + height*N : t in R } (parallel to an
   *        EDGE at pv, offset by \a height along its inward growth normal N).
   *        Used both for the overall (total-thickness) bevel point and, per
   *        layer, for each individual layer's own bend point on B13 (see
   *        buildCornerBevel(): each layer of a viscous-layer stack bends at
   *        its OWN point on the bisector, farther from O the deeper the
   *        layer, not all layers piling onto a single shared point).
   */
  //================================================================================

  bool pointOnBisectorAtHeight( const gp_Pnt& pv, const gp_Vec& T, const gp_Vec& N,
                                const gp_Dir& B13, double height, gp_Pnt& P )
  {
    const gp_Vec b13vec( B13 );
    const double bx = b13vec.Dot( T );
    const double by = b13vec.Dot( N );
    if ( Abs( by ) < 1e-12 )
      return false; // B13 (near-)parallel to the EDGE
    const double s = height / by;
    if ( s <= 0 )
      return false; // intersection behind O: pattern not applicable
    const double t = s * bx;
    P = pv.XYZ() + t * T.XYZ() + height * N.XYZ();
    return true;
  }

  //================================================================================
  /*!
   * \brief Compute B13 (bisector of E1,E3 at shared VERTEX O) and P1,P2 (where
   *        each EDGE's own inward VL offset line, extended past O, crosses B13,
   *        at the EDGE's full thickness -- i.e. its deepest/last layer).
   *        See StdMeshers_ViscousLayers2D::ComputeCornerBisector() for details;
   *        this is the actual implementation, kept in VISCOUS_2D like the rest
   *        of _ViscousBuilder2D's internals.
   */
  //================================================================================

  bool computeCornerBisector( const TopoDS_Vertex& O,
                              const TopoDS_Edge&   E1,
                              const TopoDS_Face&   F1,
                              double               thick1,
                              const TopoDS_Edge&   E3,
                              const TopoDS_Face&   F3,
                              double               thick3,
                              gp_Dir&              B13,
                              gp_Pnt&              P1,
                              gp_Pnt&              P2 )
  {
    gp_Pnt pv1, pv3;
    gp_Vec t1, n1, t3, n3;
    if ( !edgeTangentAndNormalAt( O, E1, F1, pv1, t1, n1 ) ||
         !edgeTangentAndNormalAt( O, E3, F3, pv3, t3, n3 ))
      return false;

    const gp_Vec sumN = n1 + n3;
    if ( sumN.Magnitude() < gp::Resolution() )
      return false; // E1/E3 growth normals cancel out: no well-defined bisector

    B13 = gp_Dir( sumN );

    if ( !pointOnBisectorAtHeight( pv1, t1, n1, B13, thick1, P1 ) ||
         !pointOnBisectorAtHeight( pv3, t3, n3, B13, thick3, P2 ))
      return false;

    return true;
  }

} // namespace VISCOUS_2D

//================================================================================
// StdMeshers_ViscousLayers hypothesis
//
StdMeshers_ViscousLayers2D::StdMeshers_ViscousLayers2D(int hypId, SMESH_Gen* gen)
  :StdMeshers_ViscousLayers(hypId, gen)
{
  _name = StdMeshers_ViscousLayers2D::GetHypType();
  _param_algo_dim = -2; // auxiliary hyp used by 2D algos
}
// --------------------------------------------------------------------------------
bool StdMeshers_ViscousLayers2D::SetParametersByMesh(const SMESH_Mesh*   /*theMesh*/,
                                                     const TopoDS_Shape& /*theShape*/)
{
  // TODO ???
  return false;
}
// --------------------------------------------------------------------------------
SMESH_ProxyMesh::Ptr
StdMeshers_ViscousLayers2D::Compute(SMESH_Mesh&        theMesh,
                                    const TopoDS_Face& theFace)
{
  using namespace VISCOUS_2D;
  vector< const StdMeshers_ViscousLayers2D* > hyps;
  vector< TopoDS_Shape >                      hypShapes;

  SMESH_ProxyMesh::Ptr pm = _ProxyMeshHolder::FindProxyMeshOfFace( theFace, theMesh );
  if ( !pm )
  {
    // Whether theFace itself ends up with real viscous layers on at least
    // one of its own EDGEs. NB this can differ from "findHyps() found
    // something": a ViscousLayers2D hypothesis assigned once to an ancestor
    // shape (e.g. a whole SHELL, restricted to specific EDGEs elsewhere in
    // it) is found by findHyps(ancestors=true) for EVERY FACE of that shell,
    // including ones none of whose own EDGEs are actually covered -- for
    // those, _ViscousBuilder2D::Compute() (via findEdgesWithLayers()) simply
    // returns a null proxy, same as if no hypothesis had been found at all.
    bool gotRealLayers = false;
    if ( findHyps( theMesh, theFace, hyps, hypShapes ))
    {
      VISCOUS_2D::_ViscousBuilder2D builder( theMesh, theFace, hyps, hypShapes );
      pm = builder.Compute();
      SMESH_ComputeErrorPtr error = builder.GetError();
      if ( error && !error->IsOK() )
      {
        theMesh.GetSubMesh( theFace )->GetComputeError() = error;
        gotRealLayers = true; // a real error on a real attempt
      }
      else if ( pm )
      {
        gotRealLayers = true;
      }
      if ( getenv("__ONLY__VL2D__"))
      {
        pm.reset();
        gotRealLayers = true; // preserve old behaviour: always bail out here
      }
    }
    if ( !gotRealLayers )
    {
      // theFace may still be the "bare" FACE of a reentrant corner shared by
      // two other, independently viscous-layered FACEs (see StdMeshers_
      // ViscousLayers2D document / "corner bevel" feature) -- check for that
      // pattern and, if found, splice the bevel into theFace's own proxy.
      VISCOUS_2D::_ViscousBuilder2D bareBuilder( theMesh, theFace, hyps, hypShapes );
      bareBuilder.tryBuildBareCornerProxy();
      pm = bareBuilder._proxyMesh; // friend access; stays null if no corner matched
      if ( !pm )
        pm.reset( new SMESH_ProxyMesh( theMesh ));
    }
  }
  return pm;
}
// --------------------------------------------------------------------------------
void StdMeshers_ViscousLayers2D::SetProxyMeshOfEdge( const StdMeshers_FaceSide& edgeNodes )
{
  using namespace VISCOUS_2D;
  SMESH_ProxyMesh::Ptr pm =
    _ProxyMeshHolder::FindProxyMeshOfFace( edgeNodes.Face(), *edgeNodes.GetMesh() );
  if ( !pm ) {
    _ProxyMeshOfFace* proxyMeshOfFace = new _ProxyMeshOfFace( *edgeNodes.GetMesh() );
    pm.reset( proxyMeshOfFace );
    new _ProxyMeshHolder( edgeNodes.Face(), pm );
  }
  _ProxyMeshOfFace*  proxyMeshOfFace = static_cast<_ProxyMeshOfFace*>( pm.get() );
  _ProxyMeshOfFace::_EdgeSubMesh* sm = proxyMeshOfFace->GetEdgeSubMesh( edgeNodes.EdgeID(0) );
  sm->GetUVPtStructVec() = edgeNodes.GetUVPtStruct();
}
// --------------------------------------------------------------------------------
bool StdMeshers_ViscousLayers2D::HasProxyMesh( const TopoDS_Face& face, SMESH_Mesh& mesh )
{
  return VISCOUS_2D::_ProxyMeshHolder::FindProxyMeshOfFace( face, mesh ).get();
}
// --------------------------------------------------------------------------------
bool StdMeshers_ViscousLayers2D::ComputeCornerBisector( const TopoDS_Vertex& O,
                                                        const TopoDS_Edge&   E1,
                                                        const TopoDS_Face&   F1,
                                                        double               thickness1,
                                                        const TopoDS_Edge&   E3,
                                                        const TopoDS_Face&   F3,
                                                        double               thickness3,
                                                        gp_Dir&              B13,
                                                        gp_Pnt&              P1,
                                                        gp_Pnt&              P2 )
{
  return VISCOUS_2D::computeCornerBisector( O, E1, F1, thickness1,
                                            E3, F3, thickness3, B13, P1, P2 );
}
// --------------------------------------------------------------------------------
SMESH_ComputeErrorPtr
StdMeshers_ViscousLayers2D::CheckHypothesis(SMESH_Mesh&                          theMesh,
                                            const TopoDS_Shape&                  theShape,
                                            SMESH_Hypothesis::Hypothesis_Status& theStatus)
{
  SMESH_ComputeErrorPtr error = SMESH_ComputeError::New(COMPERR_OK);
  theStatus = SMESH_Hypothesis::HYP_OK;

  TopExp_Explorer exp( theShape, TopAbs_FACE );
  for ( ; exp.More() && theStatus == SMESH_Hypothesis::HYP_OK; exp.Next() )
  {
    const TopoDS_Face& face = TopoDS::Face( exp.Current() );
    vector< const StdMeshers_ViscousLayers2D* > hyps;
    vector< TopoDS_Shape >                      hypShapes;
    if ( VISCOUS_2D::findHyps( theMesh, face, hyps, hypShapes ))
    {
      VISCOUS_2D::_ViscousBuilder2D builder( theMesh, face, hyps, hypShapes );
      builder._faceSideVec =
        StdMeshers_FaceSide::GetFaceWires( face, theMesh, true, error,
                                           NULL, SMESH_ProxyMesh::Ptr(),
                                           /*theCheckVertexNodes=*/false);
      if ( error->IsOK() && !builder.findEdgesWithLayers())
      {
        error = builder.GetError();
        if ( error && !error->IsOK() )
          theStatus = SMESH_Hypothesis::HYP_INCOMPAT_HYPS;
      }
    }
  }
  return error;
}
// --------------------------------------------------------------------------------
void StdMeshers_ViscousLayers2D::RestoreListeners() const
{
  StudyContextStruct* sc = _gen->GetStudyContext();
  std::map < int, SMESH_Mesh * >::iterator i_smesh = sc->mapMesh.begin();
  for ( ; i_smesh != sc->mapMesh.end(); ++i_smesh )
  {
    SMESH_Mesh* smesh = i_smesh->second;
    if ( !smesh ||
         !smesh->HasShapeToMesh() ||
         !smesh->GetMeshDS() ||
         !smesh->GetMeshDS()->IsUsedHypothesis( this ))
      continue;

    // set event listeners to EDGE's of FACE where this hyp is used
    TopoDS_Shape shape = i_smesh->second->GetShapeToMesh();
    for ( TopExp_Explorer face( shape, TopAbs_FACE); face.More(); face.Next() )
      if ( SMESH_Algo* algo = _gen->GetAlgo( *smesh, face.Current() ))
      {
        const std::list <const SMESHDS_Hypothesis *> & usedHyps =
          algo->GetUsedHypothesis( *smesh, face.Current(), /*ignoreAuxiliary=*/false );
        if ( std::find( usedHyps.begin(), usedHyps.end(), this ) != usedHyps.end() )
        {
          for ( TopExp_Explorer edge( face.Current(), TopAbs_EDGE); edge.More(); edge.Next() )
            VISCOUS_3D::ToClearSubWithMain( smesh->GetSubMesh( edge.Current() ), face.Current() );
          // also VERTEXes: a VERTEX may hold cached VISCOUS_2D::_CornerData
          // (reentrant-corner bevel, shared with another FACE through an
          // intervening bare FACE) that must likewise be invalidated when
          // this FACE's viscous-layers hypothesis is edited after study load.
          for ( TopExp_Explorer vertex( face.Current(), TopAbs_VERTEX); vertex.More(); vertex.Next() )
            VISCOUS_3D::ToClearSubWithMain( smesh->GetSubMesh( vertex.Current() ), face.Current() );
        }
      }
  }
}
// END StdMeshers_ViscousLayers2D hypothesis
//================================================================================

using namespace VISCOUS_2D;

//================================================================================
/*!
 * \brief Constructor of _ViscousBuilder2D
 */
//================================================================================

_ViscousBuilder2D::_ViscousBuilder2D(SMESH_Mesh&               theMesh,
                                     const TopoDS_Face&        theFace,
                                     vector< const THypVL* > & theHyps,
                                     vector< TopoDS_Shape > &  theAssignedTo):
  _mesh( &theMesh ), _face( theFace ), _helper( theMesh )
{
  _hyps.swap( theHyps );
  _hypShapes.swap( theAssignedTo );

  _helper.SetSubShape( _face );
  _helper.SetElementsOnShape( true );

  _face.Orientation( TopAbs_FORWARD ); // 2D logic works only in this case
  _surface = BRep_Tool::Surface( _face );

  _error = SMESH_ComputeError::New(COMPERR_OK);

  _nbLE = 0;
}

//================================================================================
/*!
 * \brief Stores error description and returns false
 */
//================================================================================

bool _ViscousBuilder2D::error(const string& text )
{
  _error->myName    = COMPERR_ALGO_FAILED;
  _error->myComment = string("Viscous layers builder 2D: ") + text;
  if ( SMESH_subMesh* sm = _mesh->GetSubMesh( _face ) )
  {
    SMESH_ComputeErrorPtr& smError = sm->GetComputeError();
    if ( smError && smError->myAlgo )
      _error->myAlgo = smError->myAlgo;
    smError = _error;
  }

  if (SALOME::VerbosityActivated())
    cout << "_ViscousBuilder2D::error " << text << endl;

  return false;
}

//================================================================================
/*!
 * \brief Does its job
 */
//================================================================================

SMESH_ProxyMesh::Ptr _ViscousBuilder2D::Compute()
{
  _faceSideVec = StdMeshers_FaceSide::GetFaceWires( _face, *_mesh, true, _error, &_helper );

  if ( !_error->IsOK() )
    return _proxyMesh;

  if ( !findEdgesWithLayers() ) // analysis of a shape
    return _proxyMesh;

  if ( ! makePolyLines() ) // creation of fronts
    return _proxyMesh;
    
  if ( ! inflate() ) // advance fronts
    return _proxyMesh;

  // remove elements and nodes from _face
  removeMeshFaces( _face );

  if ( !shrink() ) // shrink segments on edges w/o layers
    return _proxyMesh;

  if ( ! refine() ) // make faces
    return _proxyMesh;

  //improve();

  return _proxyMesh;
}

//================================================================================
/*!
 * \brief Finds EDGE's to make viscous layers on.
 */
//================================================================================

bool _ViscousBuilder2D::findEdgesWithLayers()
{
  // collect all EDGEs to ignore defined by _hyps
  typedef std::pair< set<TGeomID>, const THypVL* > TEdgesOfHyp;
  vector< TEdgesOfHyp > ignoreEdgesOfHyp( _hyps.size() );
  for ( size_t i = 0; i < _hyps.size(); ++i )
  {
    ignoreEdgesOfHyp[i].second = _hyps[i];
    getEdgesToIgnore( _hyps[i], _face, getMeshDS(), ignoreEdgesOfHyp[i].first );
  }

  // get all shared EDGEs
  TopTools_MapOfShape sharedEdges;
  TopTools_IndexedMapOfShape hypFaces; // faces with VL hyps
  for ( size_t i = 0; i < _hypShapes.size(); ++i )
    TopExp::MapShapes( _hypShapes[i], TopAbs_FACE, hypFaces );
  TopTools_IndexedDataMapOfShapeListOfShape facesOfEdgeMap;
  for ( int iF = 1; iF <= hypFaces.Extent(); ++iF )
    TopExp::MapShapesAndAncestors( hypFaces(iF), TopAbs_EDGE, TopAbs_FACE, facesOfEdgeMap);
  for ( int iE = 1; iE <= facesOfEdgeMap.Extent(); ++iE )
    if ( facesOfEdgeMap( iE ).Extent() > 1 )
      sharedEdges.Add( facesOfEdgeMap.FindKey( iE ));

  // fill _hypOfEdge
  if ( _hyps.size() > 1 )
  {
    // check if two hypotheses define different parameters for the same EDGE
    for ( size_t iWire = 0; iWire < _faceSideVec.size(); ++iWire )
    {
      StdMeshers_FaceSidePtr wire = _faceSideVec[ iWire ];
      for ( int iE = 0; iE < wire->NbEdges(); ++iE )
      {
        const THypVL* hyp = 0;
        const TGeomID edgeID = wire->EdgeID( iE );
        if ( !sharedEdges.Contains( wire->Edge( iE )))
        {
          for ( size_t i = 0; i < ignoreEdgesOfHyp.size(); ++i )
            if ( ! ignoreEdgesOfHyp[i].first.count( edgeID ))
            {
              if ( hyp )
                return error(SMESH_Comment("Several hypotheses define "
                                           "Viscous Layers on the edge #") << edgeID );
              hyp = ignoreEdgesOfHyp[i].second;
            }
        }
        _hypOfEdge.push_back( hyp );
        if ( !hyp )
          _ignoreShapeIds.insert( edgeID );
      }
      // check if two hypotheses define different number of viscous layers for
      // adjacent EDGEs
      const THypVL *hyp, *prevHyp = _hypOfEdge.back();
      size_t iH = _hypOfEdge.size() - wire->NbEdges();
      for ( ; iH < _hypOfEdge.size(); ++iH )
      {
        hyp = _hypOfEdge[ iH ];
        if ( hyp && prevHyp &&
             hyp->GetNumberLayers() != prevHyp->GetNumberLayers() )
        {
          return error("Two hypotheses define different number of "
                       "viscous layers on adjacent edges");
        }
        prevHyp = hyp;
      }
    }
  }
  else if ( _hyps.size() == 1 )
  {
    _ignoreShapeIds.swap( ignoreEdgesOfHyp[0].first );
  }

  // check all EDGEs of the _face to fill _ignoreShapeIds and _noShrinkVert

  int totalNbEdges = 0;
  for ( size_t iWire = 0; iWire < _faceSideVec.size(); ++iWire )
  {
    StdMeshers_FaceSidePtr wire = _faceSideVec[ iWire ];
    totalNbEdges += wire->NbEdges();
    for ( int iE = 0; iE < wire->NbEdges(); ++iE )
    {
      if ( sharedEdges.Contains( wire->Edge( iE )))
      {
        // ignore internal EDGEs (shared by several FACEs)
        const TGeomID edgeID = wire->EdgeID( iE );
        _ignoreShapeIds.insert( edgeID );

        // check if ends of an EDGE are to be added to _noShrinkVert
        const TopTools_ListOfShape& faceList = facesOfEdgeMap.FindFromKey( wire->Edge( iE ));
        TopTools_ListIteratorOfListOfShape faceIt( faceList );
        for ( ; faceIt.More(); faceIt.Next() )
        {
          const TopoDS_Shape& neighbourFace = faceIt.Value();
          if ( neighbourFace.IsSame( _face )) continue;
          SMESH_Algo* algo = _mesh->GetGen()->GetAlgo( *_mesh, neighbourFace );
          if ( !algo ) continue;

          const StdMeshers_ViscousLayers2D* viscHyp = 0;
          const list <const SMESHDS_Hypothesis *> & allHyps =
            algo->GetUsedHypothesis(*_mesh, neighbourFace, /*noAuxiliary=*/false);
          list< const SMESHDS_Hypothesis *>::const_iterator hyp = allHyps.begin();
          for ( ; hyp != allHyps.end() && !viscHyp; ++hyp )
            viscHyp = dynamic_cast<const StdMeshers_ViscousLayers2D*>( *hyp );

          // set<TGeomID> neighbourIgnoreEdges;
          // if (viscHyp)
          //   getEdgesToIgnore( viscHyp, neighbourFace, getMeshDS(), neighbourIgnoreEdges );

          for ( int iV = 0; iV < 2; ++iV )
          {
            TopoDS_Vertex vertex = iV ? wire->LastVertex(iE) : wire->FirstVertex(iE);

            // Exempt from _noShrinkVert when the reentrant-corner bevel
            // pattern applies here: wire->Edge(iE) is the bare EDGE between
            // two independently viscous-layered FACEs meeting only at
            // `vertex` through neighbourFace. buildCornerBevel() handles
            // this case properly instead of the plain "meld all layer nodes
            // into the VERTEX" fallback _noShrinkVert exists for -- but only
            // gets a chance to if `vertex` isn't blocked here first. This is
            // the same 2-hop-vs-3-hop gap findBareCornerPattern() itself had
            // to be fixed for; myArmEdge is passed null since it's only used
            // here to exclude a specific EDGE from consideration, which is
            // unnecessary in this generic context.
            TopoDS_Edge tmpArmEdge; TopoDS_Face tmpArmFace;
            double      tmpThick = 0; int tmpNbLayers = 0;
            if ( findBareCornerPattern( vertex, TopoDS_Edge(), wire->Edge( iE ),
                                        TopoDS::Face( neighbourFace ),
                                        tmpArmEdge, tmpArmFace, tmpThick, tmpNbLayers ))
              continue;

            if ( !viscHyp )
              _noShrinkVert.insert( getMeshDS()->ShapeToIndex( vertex ));
            else
            {
              PShapeIteratorPtr edgeIt = _helper.GetAncestors( vertex, *_mesh, TopAbs_EDGE );
              while ( const TopoDS_Shape* edge = edgeIt->next() )
                if ( !edge->IsSame( wire->Edge( iE )) &&
                     _helper.IsSubShape( *edge, neighbourFace ))
                {
                  const TGeomID neighbourID = getMeshDS()->ShapeToIndex( *edge );
                  bool hasVL = !sharedEdges.Contains( *edge );
                  if ( hasVL )
                  {
                    hasVL = false;
                    for ( hyp = allHyps.begin(); hyp != allHyps.end() && !hasVL; ++hyp )
                      if (( viscHyp = dynamic_cast<const THypVL*>( *hyp )))
                        hasVL = viscHyp->IsShapeWithLayers( neighbourID );
                  }
                  if ( !hasVL )
                  {
                    _noShrinkVert.insert( getMeshDS()->ShapeToIndex( vertex ));
                    break;
                  }
                }
            }
          }
        }
      }
    }
  }

  int nbMyEdgesIgnored = _ignoreShapeIds.size();

  // add VERTEXes w/o layers to _ignoreShapeIds (this is used by toShrinkForAdjacent())
  // for ( size_t iWire = 0; iWire < _faceSideVec.size(); ++iWire )
  // {
  //   StdMeshers_FaceSidePtr wire = _faceSideVec[ iWire ];
  //   for ( int iE = 0; iE < wire->NbEdges(); ++iE )
  //   {
  //     TGeomID edge1 = wire->EdgeID( iE );
  //     TGeomID edge2 = wire->EdgeID( iE+1 );
  //     if ( _ignoreShapeIds.count( edge1 ) && _ignoreShapeIds.count( edge2 ))
  //       _ignoreShapeIds.insert( getMeshDS()->ShapeToIndex( wire->LastVertex( iE )));
  //   }
  // }

  return ( nbMyEdgesIgnored < totalNbEdges );
}

//================================================================================
/*!
 * \brief Create the inner front of the viscous layers and prepare data for inflation
 */
//================================================================================

bool _ViscousBuilder2D::makePolyLines()
{
  // Create _PolyLines and _LayerEdge's

  // count total nb of EDGEs to allocate _polyLineVec
  int nbEdges = 0;
  for ( size_t iWire = 0; iWire < _faceSideVec.size(); ++iWire )
  {
    StdMeshers_FaceSidePtr wire = _faceSideVec[ iWire ];
    nbEdges += wire->NbEdges();
    if ( wire->GetUVPtStruct().empty() && wire->NbPoints() > 0 )
      return error("Invalid node parameters on some EDGE");
  }
  _polyLineVec.resize( nbEdges );

  // check if 2D normal should be computed by 3D one by means of projection
  GeomAPI_ProjectPointOnSurf* faceProj = 0;
  TopLoc_Location loc;
  {
    _LayerEdge tmpLE;
    const UVPtStruct& uv = _faceSideVec[0]->GetUVPtStruct()[0];
    gp_Pnt p = SMESH_TNodeXYZ( uv.node );
    tmpLE._uvOut.SetCoord( uv.u, uv.v );
    tmpLE._normal2D.SetCoord( 1., 0. );
    setLenRatio( tmpLE, p );
    const double r1 = tmpLE._len2dTo3dRatio;
    tmpLE._normal2D.SetCoord( 0., 1. );
    setLenRatio( tmpLE, p );
    const double r2 = tmpLE._len2dTo3dRatio;
    // projection is needed if two _len2dTo3dRatio's differ too much
    const double maxR = Max( r2, r1 );
    if ( Abs( r2-r1 )/maxR > 0.2*maxR )
      faceProj = & _helper.GetProjector( _face, loc );
  }
  _is2DIsotropic = !faceProj;

  // Assign data to _PolyLine's
  // ---------------------------

  size_t iPoLine = 0;
  for ( size_t iWire = 0; iWire < _faceSideVec.size(); ++iWire )
  {
    StdMeshers_FaceSidePtr      wire = _faceSideVec[ iWire ];
    const vector<UVPtStruct>& points = wire->GetUVPtStruct();
    int iPnt = 0;
    for ( int iE = 0; iE < wire->NbEdges(); ++iE )
    {
      _PolyLine& L  = _polyLineVec[ iPoLine++ ];
      L._index      = iPoLine-1;
      L._wire       = wire.get();
      L._edgeInd    = iE;
      L._advancable = !_ignoreShapeIds.count( wire->EdgeID( iE ));

      int iRight    = iPoLine - (( iE+1 < wire->NbEdges() ) ? 0 : wire->NbEdges() );
      L._rightLine  = &_polyLineVec[ iRight ];
      _polyLineVec[ iRight ]._leftLine = &L;

      L._firstPntInd = iPnt;
      double lastNormPar = wire->LastParameter( iE ) - 1e-10;
      while ( points[ iPnt ].normParam < lastNormPar )
        ++iPnt;
      L._lastPntInd = iPnt;
      L._lEdges.resize( Max( 3, L._lastPntInd - L._firstPntInd + 1 )); // 3 edges minimum

      // TODO: add more _LayerEdge's to strongly curved EDGEs
      // in order not to miss collisions

      double u; gp_Pnt p;
      Handle(Geom_Curve)   curve  = BRep_Tool::Curve( L._wire->Edge( iE ), loc, u, u );
      Handle(Geom2d_Curve) pcurve = L._wire->Curve2d( L._edgeInd );
      const bool reverse = (( L._wire->Edge( iE ).Orientation() == TopAbs_REVERSED ) ^
                            (_face.Orientation()                == TopAbs_REVERSED ));
      for ( int i = L._firstPntInd; i <= L._lastPntInd; ++i )
      {
        _LayerEdge& lEdge = L._lEdges[ i - L._firstPntInd ];
        u = ( i == L._firstPntInd ? wire->FirstU(iE) : points[ i ].param );
        p = SMESH_TNodeXYZ( points[ i ].node );
        setLayerEdgeData( lEdge, u, pcurve, curve, p, reverse, faceProj );
        setLenRatio( lEdge, p );
      }
      if ( L._lastPntInd - L._firstPntInd + 1 < 3 ) // add 3-d _LayerEdge in the middle
      {
        L._lEdges[2] = L._lEdges[1];
        u = 0.5 * ( wire->FirstU(iE) + wire->LastU(iE) );
        if ( !curve.IsNull() )
          p = curve->Value( u );
        else
          p = 0.5 * ( SMESH_TNodeXYZ( points[ L._firstPntInd ].node ) +
                      SMESH_TNodeXYZ( points[ L._lastPntInd ].node ));
        setLayerEdgeData( L._lEdges[1], u, pcurve, curve, p, reverse, faceProj );
        setLenRatio( L._lEdges[1], p );
      }
    }
  }

  // Fill _PolyLine's with _segments
  // --------------------------------

  double maxLen2dTo3dRatio = 0;
  for ( iPoLine = 0; iPoLine < _polyLineVec.size(); ++iPoLine )
  {
    _PolyLine& L = _polyLineVec[ iPoLine ];
    L._segments.resize( L._lEdges.size() - 1 );
    for ( size_t i = 1; i < L._lEdges.size(); ++i )
    {
      _Segment & S   = L._segments[i-1];
      S._uv[0]       = & L._lEdges[i-1]._uvIn;
      S._uv[1]       = & L._lEdges[i  ]._uvIn;
      S._indexInLine = i-1;
      if ( maxLen2dTo3dRatio < L._lEdges[i]._len2dTo3dRatio )
        maxLen2dTo3dRatio = L._lEdges[i]._len2dTo3dRatio;
    }
    // // connect _PolyLine's with segments, the 1st _LayerEdge of every _PolyLine
    // // becomes not connected to any segment
    // if ( L._leftLine->_advancable )
    //   L._segments[0]._uv[0] = & L._leftLine->_lEdges.back()._uvIn;

    L._segTree.reset( new _SegmentTree( L._segments ));
  }

  // Evaluate max possible _thickness if required layers thickness seems too high
  // ----------------------------------------------------------------------------

  _maxThickness = _hyps[0]->GetTotalThickness();
  for ( size_t iH = 1; iH < _hyps.size(); ++iH )
    _maxThickness = Max( _maxThickness, _hyps[iH]->GetTotalThickness() );

  _SegmentTree::box_type faceBndBox2D;
  for ( iPoLine = 0; iPoLine < _polyLineVec.size(); ++iPoLine )
    faceBndBox2D.Add( *_polyLineVec[ iPoLine]._segTree->getBox() );
  const double boxTol = 1e-3 * sqrt( faceBndBox2D.SquareExtent() );

  if ( _maxThickness * maxLen2dTo3dRatio > sqrt( faceBndBox2D.SquareExtent() ) / 2 )
  {
    vector< const _Segment* > foundSegs;
    double maxPossibleThick = 0;
    _SegmentIntersection intersection;
    for ( size_t iL1 = 0; iL1 < _polyLineVec.size(); ++iL1 )
    {
      _PolyLine& L1 = _polyLineVec[ iL1 ];
      _SegmentTree::box_type boxL1 = * L1._segTree->getBox();
      boxL1.Enlarge( boxTol );
      // consider case of a circle as well!
      for ( size_t iL2 = iL1; iL2 < _polyLineVec.size(); ++iL2 )
      {
        _PolyLine& L2 = _polyLineVec[ iL2 ];
        _SegmentTree::box_type boxL2 = * L2._segTree->getBox();
        boxL2.Enlarge( boxTol );
        if ( boxL1.IsOut( boxL2 ))
          continue;
        for ( size_t iLE = 1; iLE < L1._lEdges.size(); ++iLE )
        {
          foundSegs.clear();
          L2._segTree->GetSegmentsNear( L1._lEdges[iLE]._ray, foundSegs );
          for ( size_t i = 0; i < foundSegs.size(); ++i )
            // In this block periodically finding strange intersection with parameter close to 0
            // to investigate
            if ( intersection.Compute( *foundSegs[i], L1._lEdges[iLE]._ray ))
            {
              double  distToL2 = intersection._param2 / L1._lEdges[iLE]._len2dTo3dRatio;
              double psblThick = distToL2 / ( 1 + L1._advancable + L2._advancable );
              maxPossibleThick = Max( psblThick, maxPossibleThick );
            }
        }
      }
    }
    if ( maxPossibleThick > 0. )
      _maxThickness = Min( _maxThickness, maxPossibleThick );
  }

  // Adjust _LayerEdge's at _PolyLine's extremities
  // -----------------------------------------------

  for ( iPoLine = 0; iPoLine < _polyLineVec.size(); ++iPoLine )
  {
    _PolyLine& LL = _polyLineVec[ iPoLine ];
    _PolyLine& LR = *LL._rightLine;
    adjustCommonEdge( LL, LR );
  }
  // recreate _segments if some _LayerEdge's have been removed by adjustCommonEdge()
  for ( iPoLine = 0; iPoLine < _polyLineVec.size(); ++iPoLine )
  {
    _PolyLine& L = _polyLineVec[ iPoLine ];
    // if ( L._segments.size() ==  L._lEdges.size() - 1 )
    //   continue;
    L._segments.resize( L._lEdges.size() - 1 );
    for ( size_t i = 1; i < L._lEdges.size(); ++i )
    {
      _Segment & S   = L._segments[i-1];
      S._uv[0]       = & L._lEdges[i-1]._uvIn;
      S._uv[1]       = & L._lEdges[i  ]._uvIn;
      S._indexInLine = i-1;
    }
    L._segTree.reset( new _SegmentTree( L._segments ));
  }
  // connect _PolyLine's with segments, the 1st _LayerEdge of every _PolyLine
  // becomes not connected to any segment
  for ( iPoLine = 0; iPoLine < _polyLineVec.size(); ++iPoLine )
  {
    _PolyLine& L = _polyLineVec[ iPoLine ];
    if ( L._leftLine->_advancable )
      L._segments[0]._uv[0] = & L._leftLine->_lEdges.back()._uvIn;
  }

  // Fill _reachableLines.
  // ----------------------

  // compute bnd boxes taking into account the layers total thickness
  vector< _SegmentTree::box_type > lineBoxes( _polyLineVec.size() );
  for ( iPoLine = 0; iPoLine < _polyLineVec.size(); ++iPoLine )
  {
    lineBoxes[ iPoLine ] = *_polyLineVec[ iPoLine ]._segTree->getBox();
    lineBoxes[ iPoLine ].Enlarge( maxLen2dTo3dRatio * getLineThickness( iPoLine ) * 
                                  ( _polyLineVec[ iPoLine ]._advancable ? 2. : 1.2 ));
  }
  // _reachableLines
  for ( iPoLine = 0; iPoLine < _polyLineVec.size(); ++iPoLine )
  {
    _PolyLine& L1 = _polyLineVec[ iPoLine ];
    const double thick1 = getLineThickness( iPoLine );
    for ( size_t iL2 = 0; iL2 < _polyLineVec.size(); ++iL2 )
    {
      _PolyLine& L2 = _polyLineVec[ iL2 ];
      if ( iPoLine == iL2 || lineBoxes[ iPoLine ].IsOut( lineBoxes[ iL2 ]))
        continue;
      if ( !L1._advancable && ( L1._leftLine == &L2 || L1._rightLine == &L2 ))
        continue;
      // check reachability by _LayerEdge's
      int iDelta = 1; //Max( 1, L1._lEdges.size() / 100 );
      for ( size_t iLE = 1; iLE < L1._lEdges.size(); iLE += iDelta )
      {
        _LayerEdge& LE = L1._lEdges[iLE];
        if ( !lineBoxes[ iL2 ].IsOut ( LE._uvOut,
                                       LE._uvOut + LE._normal2D * thick1 * LE._len2dTo3dRatio ))
        {
          L1._reachableLines.push_back( & L2 );
          break;
        }
      }
    }
    // add self to _reachableLines
    Geom2dAdaptor_Curve pcurve( L1._wire->Curve2d( L1._edgeInd ));
    L1._isStraight2D = ( pcurve.GetType() == GeomAbs_Line );
    if ( !L1._isStraight2D )
    {
      // TODO: check carefully
      L1._reachableLines.push_back( & L1 );
    }
  }

  return true;
}

//================================================================================
/*!
 * \brief adjust common _LayerEdge of two adjacent _PolyLine's
 *  \param LL - left _PolyLine
 *  \param LR - right _PolyLine
 */
//================================================================================

void _ViscousBuilder2D::adjustCommonEdge( _PolyLine& LL, _PolyLine& LR )
{
  int nbAdvancableL = LL._advancable + LR._advancable;
  if ( nbAdvancableL == 0 )
    return;

  _LayerEdge& EL = LL._lEdges.back();
  _LayerEdge& ER = LR._lEdges.front();
  gp_XY normL    = EL._normal2D;
  gp_XY normR    = ER._normal2D;
  gp_XY tangL ( normL.Y(), -normL.X() );

  // set common direction to a VERTEX _LayerEdge shared by two _PolyLine's
  gp_XY normCommon = ( normL * int( LL._advancable ) +
                       normR * int( LR._advancable )).Normalized();
  EL._normal2D = normCommon;
  EL._ray.SetLocation ( EL._uvOut );
  EL._ray.SetDirection( EL._normal2D );
  if ( nbAdvancableL == 1 ) { // _normal2D is true normal (not average)
    EL._isBlocked = true; // prevent intersecting with _Segments of _advancable line
    EL._length2D  = 0;
  }
  // update _LayerEdge::_len2dTo3dRatio according to a new direction
  const vector<UVPtStruct>& points = LL._wire->GetUVPtStruct();
  setLenRatio( EL, SMESH_TNodeXYZ( points[ LL._lastPntInd ].node ));

  ER = EL;

  const double dotNormTang = normR * tangL;
  const bool    largeAngle = Abs( dotNormTang ) > 0.2;
  if ( largeAngle ) // not 180 degrees
  {
    // recompute _len2dTo3dRatio to take into account angle between EDGEs
    gp_Vec2d oldNorm( LL._advancable ? normL : normR );
    double angleFactor  = 1. / Max( 0.3, Cos( oldNorm.Angle( normCommon )));
    EL._len2dTo3dRatio *= angleFactor;
    ER._len2dTo3dRatio  = EL._len2dTo3dRatio;

    gp_XY normAvg = ( normL + normR ).Normalized(); // average normal at VERTEX

    if ( dotNormTang < 0. ) // ---------------------------- CONVEX ANGLE
    {
      // Remove _LayerEdge's intersecting the normAvg to avoid collisions
      // during inflate().
      //
      // find max length of the VERTEX-based _LayerEdge whose direction is normAvg
      double       maxLen2D  = _maxThickness * EL._len2dTo3dRatio;
      const gp_XY& pCommOut  = ER._uvOut;
      gp_XY        pCommIn   = pCommOut + normAvg * maxLen2D;
      _Segment segCommon( pCommOut, pCommIn );
      _SegmentIntersection intersection;
      vector< const _Segment* > foundSegs;
      for ( size_t iL1 = 0; iL1 < _polyLineVec.size(); ++iL1 )
      {
        _PolyLine& L1 = _polyLineVec[ iL1 ];
        const _SegmentTree::box_type* boxL1 = L1._segTree->getBox();
        if ( boxL1->IsOut ( pCommOut, pCommIn ))
          continue;
        for ( size_t iLE = 1; iLE < L1._lEdges.size(); ++iLE )
        {
          foundSegs.clear();
          L1._segTree->GetSegmentsNear( segCommon, foundSegs );
          for ( size_t i = 0; i < foundSegs.size(); ++i )
            if ( intersection.Compute( *foundSegs[i], segCommon ) &&
                 intersection._param2 > 1e-10 )
            {
              double len2D = intersection._param2 * maxLen2D / ( 2 + L1._advancable );
              if ( len2D < maxLen2D ) {
                maxLen2D = len2D;
                pCommIn  = pCommOut + normAvg * maxLen2D; // here length of segCommon changes
              }
            }
        }
      }

      // remove _LayerEdge's intersecting segCommon
      for ( int isR = 0; isR < 2; ++isR ) // loop on [ LL, LR ]
      {
        _PolyLine&                 L = isR ? LR : LL;
        _PolyLine::TEdgeIterator eIt = isR ? L._lEdges.begin()+1 : L._lEdges.end()-2;
        int                      dIt = isR ? +1 : -1;
        if ( nbAdvancableL == 1 && L._advancable && normL * normR > -0.01 )
          continue;  // obtuse internal angle
        // at least 3 _LayerEdge's should remain in a _PolyLine
        if ( L._lEdges.size() < 4 ) continue;
        size_t iLE = 1;
        _SegmentIntersection lastIntersection;
        for ( ; iLE < L._lEdges.size(); ++iLE, eIt += dIt )
        {
          gp_XY uvIn = eIt->_uvOut + eIt->_normal2D * _maxThickness * eIt->_len2dTo3dRatio;
          _Segment segOfEdge( eIt->_uvOut, uvIn );
          if ( !intersection.Compute( segCommon, segOfEdge ))
            break;
          lastIntersection._param1 = intersection._param1;
          lastIntersection._param2 = intersection._param2;
        }
        if ( iLE >= L._lEdges.size() - 1 )
        {
          // all _LayerEdge's intersect the segCommon, limit inflation
          // of remaining 3 _LayerEdge's
          vector< _LayerEdge > newEdgeVec( Min( 3, L._lEdges.size() ));
          newEdgeVec.front() = L._lEdges.front();
          newEdgeVec.back()  = L._lEdges.back();
          if ( newEdgeVec.size() == 3 )
          {
            newEdgeVec[1] = L._lEdges[ isR ? (L._lEdges.size() - 2) : 1 ];
            newEdgeVec[1]._len2dTo3dRatio *= lastIntersection._param2;
          }
          L._lEdges.swap( newEdgeVec );
          if ( !isR ) std::swap( lastIntersection._param1 , lastIntersection._param2 );
          L._lEdges.front()._len2dTo3dRatio *= lastIntersection._param1; // ??
          L._lEdges.back ()._len2dTo3dRatio *= lastIntersection._param2;
        }
        else if ( iLE != 1 )
        {
          // eIt points to the _LayerEdge not intersecting with segCommon
          if ( isR )
            LR._lEdges.erase( LR._lEdges.begin()+1, eIt );
          else
            LL._lEdges.erase( eIt+1, --LL._lEdges.end() );
          // eIt = isR ? L._lEdges.begin()+1 : L._lEdges.end()-2;
          // for ( size_t i = 1; i < iLE; ++i, eIt += dIt )
          //   eIt->_isBlocked = true;
        }
      }
    }
    else // ------------------------------------------ CONCAVE ANGLE
    {
      if ( nbAdvancableL == 1 )
      {
        // make that the _LayerEdge at VERTEX is not shared by LL and LR:
        // different normals is a sign that they are not shared
        _LayerEdge& notSharedEdge = LL._advancable ? LR._lEdges[0] : LL._lEdges.back();
        _LayerEdge&    sharedEdge = LR._advancable ? LR._lEdges[0] : LL._lEdges.back();

        notSharedEdge._normal2D.SetCoord( 0.,0. );
        sharedEdge._normal2D     = normAvg;
        sharedEdge._isBlocked    = false;
        notSharedEdge._isBlocked = true;
      }
    }
  }
}

//================================================================================
/*!
 * \brief initialize data of a _LayerEdge
 */
//================================================================================

void _ViscousBuilder2D::setLayerEdgeData( _LayerEdge&                 lEdge,
                                          const double                u,
                                          Handle(Geom2d_Curve)&       pcurve,
                                          Handle(Geom_Curve)&         curve,
                                          const gp_Pnt                pOut,
                                          const bool                  reverse,
                                          GeomAPI_ProjectPointOnSurf* faceProj)
{
  gp_Pnt2d uv;
  if ( faceProj && !curve.IsNull() )
  {
    uv = pcurve->Value( u );
    gp_Vec tangent; gp_Pnt p; gp_Vec du, dv;
    curve->D1( u, p, tangent );
    if ( reverse )
      tangent.Reverse();
    _surface->D1( uv.X(), uv.Y(), p, du, dv );
    gp_Vec faceNorm = du ^ dv;
    gp_Vec normal   = faceNorm ^ tangent;
    normal.Normalize();
    p = pOut.XYZ() + normal.XYZ() * /*1e-2 * */_hyps[0]->GetTotalThickness() / _hyps[0]->GetNumberLayers();
    faceProj->Perform( p );
    if ( !faceProj->IsDone() || faceProj->NbPoints() < 1 )
      return setLayerEdgeData( lEdge, u, pcurve, curve, p, reverse, NULL );
    Standard_Real U,V;
    faceProj->LowerDistanceParameters(U,V);
    lEdge._normal2D.SetCoord( U - uv.X(), V - uv.Y() );
    lEdge._normal2D.Normalize();
  }
  else
  {
    gp_Vec2d tangent;
    pcurve->D1( u, uv, tangent );
    tangent.Normalize();
    if ( reverse )
      tangent.Reverse();
    lEdge._normal2D.SetCoord( -tangent.Y(), tangent.X() );
  }
  lEdge._uvOut = lEdge._uvIn = uv.XY();
  lEdge._ray.SetLocation ( lEdge._uvOut );
  lEdge._ray.SetDirection( lEdge._normal2D );
  lEdge._isBlocked = false;
  lEdge._length2D  = 0;

  if (SALOME::VerbosityActivated())
    lEdge._ID        = _nbLE++;
}

//================================================================================
/*!
 * \brief Compute and set _LayerEdge::_len2dTo3dRatio
 */
//================================================================================

void _ViscousBuilder2D::setLenRatio( _LayerEdge& LE, const gp_Pnt& pOut )
{
  const double probeLen2d = 1e-3;

  gp_Pnt2d p2d = LE._uvOut + LE._normal2D * probeLen2d;
  gp_Pnt   p3d = _surface->Value( p2d.X(), p2d.Y() );
  double len3d = p3d.Distance( pOut );
  if ( len3d < std::numeric_limits<double>::min() )
    LE._len2dTo3dRatio = std::numeric_limits<double>::min();
  else
    LE._len2dTo3dRatio = probeLen2d / len3d;
}

//================================================================================
/*!
 * \brief Increase length of _LayerEdge's to reach the required thickness of layers
 */
//================================================================================

bool _ViscousBuilder2D::inflate()
{
  // Limit size of inflation step by geometry size found by
  // itersecting _LayerEdge's with _Segment's
  double minSize = _maxThickness, maxSize = 0;
  vector< const _Segment* > foundSegs;
  _SegmentIntersection intersection;
  for ( size_t iL1 = 0; iL1 < _polyLineVec.size(); ++iL1 )
  {
    _PolyLine& L1 = _polyLineVec[ iL1 ];
    for ( size_t iL2 = 0; iL2 < L1._reachableLines.size(); ++iL2 )
    {
      _PolyLine& L2 = * L1._reachableLines[ iL2 ];
      for ( size_t iLE = 1; iLE < L1._lEdges.size(); ++iLE )
      {
        foundSegs.clear();
        L2._segTree->GetSegmentsNear( L1._lEdges[iLE]._ray, foundSegs );
        for ( size_t i = 0; i < foundSegs.size(); ++i )
          if ( ! L1.IsAdjacent( *foundSegs[i], & L1._lEdges[iLE] ) &&
               intersection.Compute( *foundSegs[i], L1._lEdges[iLE]._ray ))
          {
            double distToL2 = intersection._param2 / L1._lEdges[iLE]._len2dTo3dRatio;
            double     size = distToL2 / ( 1 + L1._advancable + L2._advancable );
            if ( 1e-10 < size && size < minSize )
              minSize = size;
            if ( size > maxSize )
              maxSize = size;
          }
      }
    }
  }
  if ( minSize > maxSize ) // no collisions possible
    maxSize = _maxThickness;
#ifdef __myDEBUG
  cout << "-- minSize = " << minSize << ", maxSize = " << maxSize << endl;
#endif

  double curThick = 0, stepSize = minSize;
  int nbSteps = 0;
  if ( maxSize > _maxThickness )
    maxSize = _maxThickness;
  while ( curThick < maxSize )
  {
    curThick += stepSize * 1.25;
    if ( curThick > _maxThickness )
      curThick = _maxThickness;

    // Elongate _LayerEdge's
    for ( size_t iL = 0; iL < _polyLineVec.size(); ++iL )
    {
      _PolyLine& L = _polyLineVec[ iL ];
      if ( !L._advancable ) continue;
      const double lineThick = Min( curThick, getLineThickness( iL ));
      bool lenChange = false;
      for ( size_t iLE = L.FirstLEdge(); iLE < L._lEdges.size(); ++iLE )
        lenChange |= L._lEdges[iLE].SetNewLength( lineThick );
      // for ( int k=0; k<L._segments.size(); ++k)
      //   cout << "( " << L._segments[k].p1().X() << ", " <<L._segments[k].p1().Y() << " ) "
      //        << "( " << L._segments[k].p2().X() << ", " <<L._segments[k].p2().Y() << " ) "
      //        << endl;
      if ( lenChange )
        L._segTree.reset( new _SegmentTree( L._segments ));
    }

    // Avoid intersection of _Segment's
    bool allBlocked = fixCollisions();
    if ( allBlocked )
    {
      break; // no more inflating possible
    }
    stepSize = Max( stepSize , _maxThickness / 10. );
    nbSteps++;
  }

  // if (nbSteps == 0 )
  //   return error("failed at the very first inflation step");


  // remove _LayerEdge's of one line intersecting with each other
  for ( size_t iL = 0; iL < _polyLineVec.size(); ++iL )
  {
    _PolyLine& L = _polyLineVec[ iL ];
    if ( !L._advancable ) continue;

    // replace an inactive (1st) _LayerEdge with an active one of a neighbour _PolyLine
    if ( /*!L._leftLine->_advancable &&*/ L.IsCommonEdgeShared( *L._leftLine ) ) {
      L._lEdges[0] = L._leftLine->_lEdges.back();
    }
    if ( !L._rightLine->_advancable && L.IsCommonEdgeShared( *L._rightLine ) ) {
      L._lEdges.back() = L._rightLine->_lEdges[0];
    }

    _SegmentIntersection intersection;
    for ( int isR = 0; ( isR < 2 && L._lEdges.size() > 2 ); ++isR )
    {
      int nbRemove = 0, deltaIt = isR ? -1 : +1;
      _PolyLine::TEdgeIterator eIt = isR ? L._lEdges.end()-1 : L._lEdges.begin();
      if ( eIt->_length2D == 0 ) continue;
      _Segment seg1( eIt->_uvOut, eIt->_uvIn );
      for ( eIt += deltaIt; nbRemove < (int)L._lEdges.size()-1; eIt += deltaIt )
      {
        _Segment seg2( eIt->_uvOut, eIt->_uvIn );
        if ( !intersection.Compute( seg1, seg2 ))
          break;
        ++nbRemove;
      }
      if ( nbRemove > 0 ) {
        if ( nbRemove == (int)L._lEdges.size()-1 ) // 1st and last _LayerEdge's intersect
        {
          --nbRemove;
          _LayerEdge& L0 = L._lEdges.front();
          _LayerEdge& L1 = L._lEdges.back();
          L0._length2D *= intersection._param1 * 0.5;
          L1._length2D *= intersection._param2 * 0.5;
          L0._uvIn = L0._uvOut + L0._normal2D * L0._length2D;
          L1._uvIn = L1._uvOut + L1._normal2D * L1._length2D;
          if ( L.IsCommonEdgeShared( *L._leftLine ))
            L._leftLine->_lEdges.back() = L0;
        }
        if ( isR )
          L._lEdges.erase( L._lEdges.end()-nbRemove-1,
                           L._lEdges.end()-nbRemove );
        else
          L._lEdges.erase( L._lEdges.begin()+1,
                           L._lEdges.begin()+1+nbRemove );
      }
    }
  }
  return true;
}

//================================================================================
/*!
 * \brief Remove intersection of _PolyLine's
 */
//================================================================================

bool _ViscousBuilder2D::fixCollisions()
{
  // look for intersections of _Segment's by intersecting _LayerEdge's with
  // _Segment's
  vector< const _Segment* > foundSegs;
  _SegmentIntersection intersection;

  list< pair< _LayerEdge*, double > > edgeLenLimitList;
  list< _LayerEdge* >                 blockedEdgesList;

  for ( size_t iL1 = 0; iL1 < _polyLineVec.size(); ++iL1 )
  {
    _PolyLine& L1 = _polyLineVec[ iL1 ];
    //if ( !L1._advancable ) continue;
    for ( size_t iL2 = 0; iL2 < L1._reachableLines.size(); ++iL2 )
    {
      _PolyLine& L2 = * L1._reachableLines[ iL2 ];
      for ( size_t iLE = L1.FirstLEdge(); iLE < L1._lEdges.size(); ++iLE )
      {
        _LayerEdge& LE1 = L1._lEdges[iLE];
        if ( LE1._isBlocked ) continue;
        foundSegs.clear();
        L2._segTree->GetSegmentsNear( LE1._ray, foundSegs );
        for ( size_t i = 0; i < foundSegs.size(); ++i )
        {
          if ( ! L1.IsAdjacent( *foundSegs[i], &LE1 ) &&
               intersection.Compute( *foundSegs[i], LE1._ray ))
          {
            const double dist2DToL2 = intersection._param2;
            double         newLen2D = dist2DToL2 / 2;
            if ( newLen2D < 1.1 * LE1._length2D ) // collision!
            {
              if ( newLen2D > 0 || !L1._advancable )
              {
                blockedEdgesList.push_back( &LE1 );
                if ( L1._advancable && newLen2D > 0 )
                {
                  edgeLenLimitList.push_back( make_pair( &LE1, newLen2D ));
                  blockedEdgesList.push_back( &L2._lEdges[ foundSegs[i]->_indexInLine     ]);
                  blockedEdgesList.push_back( &L2._lEdges[ foundSegs[i]->_indexInLine + 1 ]);
                }
                else // here dist2DToL2 < 0 and LE1._length2D == 0
                {
                  _LayerEdge* LE2[2] = { & L2._lEdges[ foundSegs[i]->_indexInLine     ],
                                         & L2._lEdges[ foundSegs[i]->_indexInLine + 1 ] };
                  _Segment outSeg2( LE2[0]->_uvOut, LE2[1]->_uvOut );
                  intersection.Compute( outSeg2, LE1._ray );
                  newLen2D = intersection._param2 / 2;
                  if ( newLen2D > 0 )
                  {
                    edgeLenLimitList.push_back( make_pair( LE2[0], newLen2D ));
                    edgeLenLimitList.push_back( make_pair( LE2[1], newLen2D ));
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  // limit length of _LayerEdge's that are extrema of _PolyLine's
  // to avoid intersection of these _LayerEdge's
  for ( size_t iL1 = 0; iL1 < _polyLineVec.size(); ++iL1 )
  {
    _PolyLine& L = _polyLineVec[ iL1 ];
    if ( L._lEdges.size() < 4 ) // all intermediate _LayerEdge's intersect with extremum ones
    {
      _LayerEdge& LEL = L._leftLine->_lEdges.back();
      _LayerEdge& LER = L._lEdges.back();
      _Segment segL( LEL._uvOut, LEL._uvIn );
      _Segment segR( LER._uvOut, LER._uvIn );
      double newLen2DL, newLen2DR;
      if ( intersection.Compute( segL, LER._ray ))
      {
        newLen2DR = intersection._param2 / 2;
        newLen2DL = LEL._length2D * intersection._param1 / 2;
      }
      else if ( intersection.Compute( segR, LEL._ray ))
      {
        newLen2DL = intersection._param2 / 2;
        newLen2DR = LER._length2D * intersection._param1 / 2;
      }
      else
      {
        continue;
      }
      if ( newLen2DL > 0 && newLen2DR > 0 )
      {
        if ( newLen2DL < 1.1 * LEL._length2D )
          edgeLenLimitList.push_back( make_pair( &LEL, newLen2DL ));
        if ( newLen2DR < 1.1 * LER._length2D )
          edgeLenLimitList.push_back( make_pair( &LER, newLen2DR ));
      }
    }
  }

  // set limited length to _LayerEdge's
  list< pair< _LayerEdge*, double > >::iterator edge2Len = edgeLenLimitList.begin();
  for ( ; edge2Len != edgeLenLimitList.end(); ++edge2Len )
  {
    _LayerEdge* LE = edge2Len->first;
    if ( LE->_length2D > edge2Len->second )
    {
      LE->_isBlocked = false;
      LE->SetNewLength( edge2Len->second / LE->_len2dTo3dRatio );
    }
    LE->_isBlocked = true;
  }

  // block inflation of _LayerEdge's
  list< _LayerEdge* >::iterator edge = blockedEdgesList.begin();
  for ( ; edge != blockedEdgesList.end(); ++edge )
    (*edge)->_isBlocked = true;

  // find a not blocked _LayerEdge
  for ( size_t iL = 0; iL < _polyLineVec.size(); ++iL )
  {
    _PolyLine& L = _polyLineVec[ iL ];
    if ( !L._advancable ) continue;
    for ( size_t iLE = L.FirstLEdge(); iLE < L._lEdges.size(); ++iLE )
      if ( !L._lEdges[ iLE ]._isBlocked )
        return false;
  }

  return true;
}

//================================================================================
/*!
 * \brief Create new edges and shrink edges existing on a non-advancable _PolyLine
 *        adjacent to an advancable one.
 */
//================================================================================

bool _ViscousBuilder2D::shrink()
{
  gp_Pnt2d uv; //gp_Vec2d tangent;
  _SegmentIntersection intersection;
  double sign;

  for ( size_t iL1 = 0; iL1 < _polyLineVec.size(); ++iL1 )
  {
    _PolyLine& L = _polyLineVec[ iL1 ]; // line with no layers
    if ( L._advancable )
      continue;
    const int nbAdvancable = ( L._rightLine->_advancable + L._leftLine->_advancable );
    if ( nbAdvancable == 0 )
      continue;

    const TopoDS_Vertex&  V1 = L._wire->FirstVertex( L._edgeInd );
    const TopoDS_Vertex&  V2 = L._wire->LastVertex ( L._edgeInd );
    const int           v1ID = getMeshDS()->ShapeToIndex( V1 );
    const int           v2ID = getMeshDS()->ShapeToIndex( V2 );
    const bool isShrinkableL = ! _noShrinkVert.count( v1ID ) && L._leftLine->_advancable;
    const bool isShrinkableR = ! _noShrinkVert.count( v2ID ) && L._rightLine->_advancable;
    if ( !isShrinkableL && !isShrinkableR )
      continue;

    const TopoDS_Edge&        E = L._wire->Edge       ( L._edgeInd );
    const int            edgeID = L._wire->EdgeID     ( L._edgeInd );
    const double        edgeLen = L._wire->EdgeLength ( L._edgeInd );
    Handle(Geom2d_Curve) pcurve = L._wire->Curve2d    ( L._edgeInd );
    const bool     edgeReversed = ( E.Orientation() == TopAbs_REVERSED );

    SMESH_MesherHelper helper( *_mesh ); // to create nodes and edges on E
    helper.SetSubShape( E );
    helper.SetElementsOnShape( true );

    // Check a FACE adjacent to _face by E
    bool existingNodesFound = false;
    TopoDS_Face adjFace;
    PShapeIteratorPtr faceIt = _helper.GetAncestors( E, *_mesh, TopAbs_FACE );
    while ( const TopoDS_Shape* f = faceIt->next() )
      if ( !_face.IsSame( *f ))
      {
        adjFace = TopoDS::Face( *f );
        SMESH_ProxyMesh::Ptr pm = _ProxyMeshHolder::FindProxyMeshOfFace( adjFace, *_mesh );
        if ( !pm || pm->NbProxySubMeshes() == 0 /*|| !pm->GetProxySubMesh( E )*/)
        {
          // There are no viscous layers on an adjacent FACE, clear it's 2D mesh
          removeMeshFaces( adjFace );
          // if ( removeMeshFaces( adjFace ))
          //   _clearedFaces.push_back( adjFace ); // to re-compute after all
        }
        else
        {
          // There are viscous layers on the adjacent FACE; shrink must be already done;
          //
          // copy layer nodes
          //
          const vector<UVPtStruct>& points = L._wire->GetUVPtStruct();
          int iPFrom = L._firstPntInd, iPTo = L._lastPntInd;
          if ( isShrinkableL )
          {
            const THypVL* hyp = getLineHypothesis( L._leftLine->_index );
            vector<gp_XY>& uvVec = L._lEdges.front()._uvRefined;
            for ( int i = 0; i < hyp->GetNumberLayers(); ++i ) {
              const UVPtStruct& uvPt = points[ iPFrom + i + 1 ];
              L._leftNodes.push_back( uvPt.node );
              uvVec.push_back ( pcurve->Value( uvPt.param ).XY() );
            }
            iPFrom += hyp->GetNumberLayers();
          }
          if ( isShrinkableR )
          {
            const THypVL* hyp = getLineHypothesis( L._rightLine->_index );
            vector<gp_XY>& uvVec = L._lEdges.back()._uvRefined;
            for ( int i = 0; i < hyp->GetNumberLayers(); ++i ) {
              const UVPtStruct& uvPt = points[ iPTo - i - 1 ];
              L._rightNodes.push_back( uvPt.node );
              uvVec.push_back ( pcurve->Value( uvPt.param ).XY() );
            }
            iPTo -= hyp->GetNumberLayers();
          }
          // make proxy sub-mesh data of present nodes
          //
          UVPtStructVec nodeDataVec( & points[ iPFrom ], & points[ iPTo + 1 ]);

          double normSize = nodeDataVec.back().normParam - nodeDataVec.front().normParam;
          for ( int iP = nodeDataVec.size()-1; iP >= 0 ; --iP )
            nodeDataVec[iP].normParam =
              ( nodeDataVec[iP].normParam - nodeDataVec[0].normParam ) / normSize;

          const SMDS_MeshNode* n = nodeDataVec.front().node;
          if ( n->GetPosition()->GetTypeOfPosition() == SMDS_TOP_VERTEX )
            nodeDataVec.front().param = L._wire->FirstU( L._edgeInd );
          n = nodeDataVec.back().node;
          if ( n->GetPosition()->GetTypeOfPosition() == SMDS_TOP_VERTEX )
            nodeDataVec.back().param = L._wire->LastU( L._edgeInd );

          _ProxyMeshOfFace::_EdgeSubMesh* myEdgeSM = getProxyMesh()->GetEdgeSubMesh( edgeID );
          myEdgeSM->SetUVPtStructVec( nodeDataVec );

          existingNodesFound = true;
          break;
        }
      } // loop on FACEs sharing E

    // Check if L is an already shrinked seam
    if ( adjFace.IsNull() && _helper.IsRealSeam( edgeID ))
      if ( L._wire->Edge( L._edgeInd ).Orientation() == TopAbs_FORWARD )
        continue;
    // Commented as a case with a seam EDGE (issue 0052461) is hard to support
    // because SMESH_ProxyMesh can't hold different sub-meshes for two
    // 2D representations of the seam. But such a case is not a real practice one.
    // {
    //   for ( int iL2 = iL1-1; iL2 > -1; --iL2 )
    //   {
    //     _PolyLine& L2 = _polyLineVec[ iL2 ];
    //     if ( edgeID == L2._wire->EdgeID( L2._edgeInd ))
    //     {
    //       // copy layer nodes
    //       const int seamPar = _helper.GetPeriodicIndex();
    //       vector<gp_XY>& uvVec = L._lEdges.front()._uvRefined;
    //       if ( isShrinkableL )
    //       {
    //         L._leftNodes = L2._rightNodes;
    //         uvVec = L2._lEdges.back()._uvRefined;
    //       }
    //       if ( isShrinkableR )
    //       {
    //         L._rightNodes = L2._leftNodes;
    //         uvVec = L2._lEdges.front()._uvRefined;
    //       }
    //       for ( size_t i = 0; i < uvVec.size(); ++i )
    //       {
    //         gp_XY & uv = uvVec[i];
    //         uv.SetCoord( seamPar, _helper.GetOtherParam( uv.Coord( seamPar )));
    //       }

    //       existingNodesFound = true;
    //       break;
    //     }
    //   }
    // }

    if ( existingNodesFound )
      continue; // nothing more to do in this case

    double u1 = L._wire->FirstU( L._edgeInd ), uf = u1;
    double u2 = L._wire->LastU ( L._edgeInd ), ul = u2;

    // a ratio to pass 2D <--> 1D
    const double len1D = 1e-3;
    const double len2D = pcurve->Value(uf).Distance( pcurve->Value(uf+len1D));
    double len1dTo2dRatio = len1D / len2D;

    // create a vector of proxy nodes
    const vector<UVPtStruct>& points = L._wire->GetUVPtStruct();
    UVPtStructVec nodeDataVec( & points[ L._firstPntInd ],
                               & points[ L._lastPntInd + 1 ]);
    nodeDataVec.front().param = u1; // U on vertex is correct on only one of shared edges
    nodeDataVec.back ().param = u2;
    nodeDataVec.front().normParam = 0;
    nodeDataVec.back ().normParam = 1;

    // Get length of existing segments (from an edge start to a node) and their nodes
    vector< double > segLengths( nodeDataVec.size() - 1 );
    BRepAdaptor_Curve curve( E );
    for ( size_t iP = 1; iP < nodeDataVec.size(); ++iP )
    {
      const double len = GCPnts_AbscissaPoint::Length( curve, uf, nodeDataVec[iP].param );
      segLengths[ iP-1 ] = len;
    }

    // Move first and last parameters on EDGE (U of n1) according to layers' thickness
    // and create nodes of layers on EDGE ( -x-x-x )

    // Before
    //  n1    n2    n3    n4
    //  x-----x-----x-----x-----
    //  |  e1    e2    e3    e4

    // After
    //  n1          n2    n3
    //  x-x-x-x-----x-----x----
    //  | | | |  e1    e2    e3

    int isRShrinkedForAdjacent = 0;
    UVPtStructVec nodeDataForAdjacent;
    for ( int isR = 0; isR < 2; ++isR )
    {
      _PolyLine* L2 = isR ? L._rightLine : L._leftLine; // line with layers
      if ( !L2->_advancable &&
           !toShrinkForAdjacent( adjFace, E, L._wire->FirstVertex( L._edgeInd + isR )))
        continue;
      if ( isR ? !isShrinkableR : !isShrinkableL )
        continue;

      // If L2 (my own EDGE with viscous layers) meets, at this VERTEX, another
      // EDGE with viscous layers on a DIFFERENT FACE reachable only through a
      // bare adjFace, bevel the corner along their bisector instead of the
      // usual wall-shrink below: buildCornerBevel() builds its own explicit
      // 2D elements covering the EXTRA wedge beyond V, out to the bevel
      // point. This FACE's own regular (non-corner) fill still ends at the
      // plain VERTEX V exactly as it always did (nodeDataVec/E's own
      // discretization near V is left untouched on purpose, unlike the
      // wall-shrink below) -- V's own natural interior angle and the extra
      // bevel wedge beyond it are two different, non-overlapping angular
      // ranges from V, so leaving V alone here and letting buildCornerBevel's
      // quads cover the rest, independently, is what makes the two pieces
      // meet with no gap and no double coverage (confirmed by an area check).
      if ( L2->_advancable )
      {
        // NB: NOT L._wire->FirstVertex(L._edgeInd+isR) -- that doesn't wrap
        // when L is the wire's last EDGE (unlike EdgeID(), FirstVertex()/
        // LastVertex() don't take a modulo), silently returning a null
        // VERTEX. V1/V2 (computed once above from L._edgeInd alone) are the
        // safe way to get either end of L.
        const TopoDS_Vertex& vAtEnd    = isR ? V2 : V1;
        const TopoDS_Edge&   armE      = L2->_wire->Edge( L2->_edgeInd );
        const int            iPEndHere = isR ? nodeDataVec.size() - 1 : 0;
        const SMDS_MeshNode* vNode     = nodeDataVec[ iPEndHere ].node;
        const SMDS_MeshNode* nearCornerNode = 0;
        if ( buildCornerBevel( vAtEnd, armE, E, adjFace, vNode,
                               getLineHypothesis( L2->_index ), nearCornerNode ))
        {
          // Splice THIS FACE's own proxy for E (myBareEdge, from
          // buildCornerBevel()'s point of view) so it, too, retreats from V
          // to the same near-corner node used for this arm's own official
          // per-layer band at V (isCornerBevelL/R in refine()) -- otherwise
          // the underlying 2D algorithm's core fill for _face still reaches
          // V and would overlap whatever refine() builds on that chain (see
          // the long isCornerBevelL/R comment in refine() for why a naive
          // version of that chain, without this splice, measurably
          // overlapped in an earlier attempt).
          if ( nearCornerNode && nearCornerNode != vNode && nodeDataVec.size() > 1 )
          {
            const int neighbourIdx = isR ? iPEndHere - 1 : iPEndHere + 1;
            const gp_Pnt2d seedUV( nodeDataVec[ neighbourIdx ].u, nodeDataVec[ neighbourIdx ].v );
            const gp_Pnt   cornerXYZ = SMESH_TNodeXYZ( nearCornerNode );
            const gp_Pnt2d cornerUV = _helper.GetSurface( _face )->NextValueOfUV
              ( seedUV, cornerXYZ, BRep_Tool::Tolerance( _face ));
            nodeDataVec[ iPEndHere ].node  = nearCornerNode;
            nodeDataVec[ iPEndHere ].u     = cornerUV.X();
            nodeDataVec[ iPEndHere ].v     = cornerUV.Y();
            nodeDataVec[ iPEndHere ].param = nodeDataVec[ neighbourIdx ].param;
          }
          continue; // no wall-shrink node created on E for this VERTEX end
        }
      }

      double & u = isR ? u2 : u1; // param to move
      double  u0 = isR ? ul : uf; // init value of the param to move
      int  iPEnd = isR ? nodeDataVec.size() - 1 : 0;

      _LayerEdge& nearLE = isR ? L._lEdges.back() : L._lEdges.front();
      _LayerEdge&  farLE = isR ? L._lEdges.front() : L._lEdges.back();

      // try to find length of advancement along L by intersecting L with
      // an adjacent _Segment of L2

      double& length2D = nearLE._length2D;
      double  length1D = 0;
      sign = ( isR ^ edgeReversed ) ? -1. : 1.;

      bool isConvex = false;
      if ( L2->_advancable )
      {
        const uvPtStruct& tang2P1 = points[ isR ? L2->_firstPntInd   : L2->_lastPntInd ];
        const uvPtStruct& tang2P2 = points[ isR ? L2->_firstPntInd+1 : L2->_lastPntInd-1 ];
        gp_XY seg2Dir( tang2P2.u - tang2P1.u,
                       tang2P2.v - tang2P1.v );
        int iFSeg2 = isR ? 0 : L2->_segments.size() - 1;
        int iLSeg2 = isR ? 1 : L2->_segments.size() - 2;
        gp_XY uvLSeg2In  = L2->_lEdges[ iLSeg2 ]._uvIn;
        Handle(Geom2d_Line) seg2Line = new Geom2d_Line( uvLSeg2In, seg2Dir );

        Geom2dAdaptor_Curve edgeCurve( pcurve, Min( uf, ul ), Max( uf, ul ));
        Geom2dAdaptor_Curve seg2Curve( seg2Line );
        Geom2dInt_GInter     curveInt( edgeCurve, seg2Curve, 1e-7, 1e-7 );

        // In the older version length2D was set to this value only inside the !convex if block
        // But it seems that length2D can be set here anyway, because if not set valid value of length2D here,
        // it will be calculated later using length1D, and it can be not valid in cases if length1D is too large or too small.
        length2D = L2->_lEdges[iFSeg2]._length2D;

        /*  convex VERTEX
         *                     L  seg2
         *                    |  o---o---
         *                    | /    |
         *                    |/     |  L2
         *                    x------x---      */
        /* concave VERTEX
         *                     o-----o---
         *                      \    |
         *                       \   |  L2
         *                        x--x---
         *                       /
         *                    L /               */
        isConvex = ( curveInt.IsDone() && !curveInt.IsEmpty() );
        if ( isConvex )
        {
          length1D = Abs( u - curveInt.Point( 1 ).ParamOnFirst() );
          length2D = Max(length2D, length1D / len1dTo2dRatio);
        }
      }
      else // L2 is advancable but in the face adjacent by L
      {
        length2D = farLE._length2D;
        if ( length2D == 0 ) {
          _LayerEdge& neighborLE =
            ( isR ? L._leftLine->_lEdges.back() : L._rightLine->_lEdges.front() );
          length2D = neighborLE._length2D;
          if ( length2D == 0 )
            length2D = _maxThickness * nearLE._len2dTo3dRatio;
        }
      }

      // move u to the internal boundary of layers
      //  u --> u
      //  x-x-x-x-----x-----x----
      double maxLen3D = Min( _maxThickness, edgeLen / ( 1 + nbAdvancable ));
      double maxLen2D = maxLen3D * nearLE._len2dTo3dRatio;
      if ( !length2D ) length2D = length1D / len1dTo2dRatio;
      if ( Abs( length2D ) > maxLen2D )
        length2D = maxLen2D;
      nearLE._uvIn = nearLE._uvOut + nearLE._normal2D * length2D;

      u += length2D * len1dTo2dRatio * sign;
      nodeDataVec[ iPEnd ].param = u;

      gp_Pnt2d newUV = pcurve->Value( u );
      nodeDataVec[ iPEnd ].u = newUV.X();
      nodeDataVec[ iPEnd ].v = newUV.Y();

      // compute params of layers on L
      vector<double> heights;
      const THypVL* hyp = getLineHypothesis( L2->_index );
      calcLayersHeight( u - u0, heights, hyp );
      //
      vector< double > params( heights.size() );
      for ( size_t i = 0; i < params.size(); ++i )
        params[ i ] = u0 + heights[ i ];

      // create nodes of layers and edges between them
      //  x-x-x-x---
      vector< const SMDS_MeshNode* >& layersNode = isR ? L._rightNodes : L._leftNodes;
      vector<gp_XY>& nodeUV = ( isR ? L._lEdges.back() : L._lEdges[0] )._uvRefined;
      nodeUV.resize    ( hyp->GetNumberLayers() );
      layersNode.resize( hyp->GetNumberLayers() );
      const SMDS_MeshNode* vertexNode = nodeDataVec[ iPEnd ].node;
      const SMDS_MeshNode *  prevNode = vertexNode;
      for ( size_t i = 0; i < params.size(); ++i )
      {
        const gp_Pnt p  = curve.Value( params[i] );
        layersNode[ i ] = helper.AddNode( p.X(), p.Y(), p.Z(), /*id=*/0, params[i] );
        nodeUV    [ i ] = pcurve->Value( params[i] ).XY();
        helper.AddEdge( prevNode, layersNode[ i ] );
        prevNode = layersNode[ i ];
      }

      // store data of layer nodes made for adjacent FACE
      if ( !L2->_advancable )
      {
        isRShrinkedForAdjacent = isR;
        nodeDataForAdjacent.resize( hyp->GetNumberLayers() );

        size_t iFrw = 0, iRev = nodeDataForAdjacent.size()-1, *i = isR ? &iRev : &iFrw;
        nodeDataForAdjacent[ *i ] = points[ isR ? L._lastPntInd : L._firstPntInd ];
        nodeDataForAdjacent[ *i ].param     = u0;
        nodeDataForAdjacent[ *i ].normParam = isR;
        for ( ++iFrw, --iRev; iFrw < layersNode.size(); ++iFrw, --iRev )
        {
          nodeDataForAdjacent[ *i ].node  = layersNode[ iFrw - 1 ];
          nodeDataForAdjacent[ *i ].u     = nodeUV    [ iFrw - 1 ].X();
          nodeDataForAdjacent[ *i ].v     = nodeUV    [ iFrw - 1 ].Y();
          nodeDataForAdjacent[ *i ].param = params    [ iFrw - 1 ];
        }
      }
      // replace a node on vertex by a node of last (most internal) layer
      // in a segment on E
      SMDS_ElemIteratorPtr segIt = vertexNode->GetInverseElementIterator( SMDSAbs_Edge );
      const SMDS_MeshNode* segNodes[3];
      while ( segIt->more() )
      {
        const SMDS_MeshElement* segment = segIt->next();
        if ( segment->getshapeId() != edgeID ) continue;

        const int nbNodes = segment->NbNodes();
        for ( int i = 0; i < nbNodes; ++i )
        {
          const SMDS_MeshNode* n = segment->GetNode( i );
          segNodes[ i ] = ( n == vertexNode ? layersNode.back() : n );
        }
        getMeshDS()->ChangeElementNodes( segment, segNodes, nbNodes );
        break;
      }
      nodeDataVec[ iPEnd ].node = layersNode.back();

    } // loop on the extremities of L

    // Shrink edges to fit in between the layers at EDGE ends

    double newLength = GCPnts_AbscissaPoint::Length( curve, u1, u2 );
    double lenRatio  = newLength / edgeLen * ( edgeReversed ? -1. : 1. );
    for ( size_t iP = 1; iP < nodeDataVec.size()-1; ++iP )
    {
      const SMDS_MeshNode* oldNode = nodeDataVec[iP].node;

      GCPnts_AbscissaPoint discret( curve, segLengths[iP-1] * lenRatio, u1 );
      if ( !discret.IsDone() )
        throw SALOME_Exception(LOCALIZED("GCPnts_AbscissaPoint failed"));

      nodeDataVec[iP].param = discret.Parameter();
      if ( oldNode->GetPosition()->GetTypeOfPosition() != SMDS_TOP_EDGE )
        throw SALOME_Exception(SMESH_Comment("ViscousBuilder2D: not SMDS_TOP_EDGE node position: ")
                               << oldNode->GetPosition()->GetTypeOfPosition()
                               << " of node " << oldNode->GetID());
      SMDS_EdgePositionPtr pos = oldNode->GetPosition();
      pos->SetUParameter( nodeDataVec[iP].param );

      gp_Pnt newP = curve.Value( nodeDataVec[iP].param );
      getMeshDS()->MoveNode( oldNode, newP.X(), newP.Y(), newP.Z() );

      gp_Pnt2d newUV = pcurve->Value( nodeDataVec[iP].param ).XY();
      nodeDataVec[iP].u         = newUV.X();
      nodeDataVec[iP].v         = newUV.Y();
      nodeDataVec[iP].normParam = segLengths[iP-1] / edgeLen;
      // nodeDataVec[iP].x         = segLengths[iP-1] / edgeLen;
      // nodeDataVec[iP].y         = segLengths[iP-1] / edgeLen;
    }

    // Add nodeDataForAdjacent to nodeDataVec

    if ( !nodeDataForAdjacent.empty() )
    {
      const double par1      = isRShrinkedForAdjacent ? u2 : uf;
      const double par2      = isRShrinkedForAdjacent ? ul : u1;
      const double shrinkLen = GCPnts_AbscissaPoint::Length( curve, par1, par2 );

      // compute new normParam for nodeDataVec
      for ( size_t iP = 0; iP < nodeDataVec.size()-1; ++iP )
        nodeDataVec[iP+1].normParam = segLengths[iP] / ( edgeLen + shrinkLen );
      double normDelta = 1 - nodeDataVec.back().normParam;
      if ( !isRShrinkedForAdjacent )
        for ( size_t iP = 0; iP < nodeDataVec.size(); ++iP )
          nodeDataVec[iP].normParam += normDelta;

      // compute new normParam for nodeDataForAdjacent
      const double deltaR = isRShrinkedForAdjacent ? nodeDataVec.back().normParam : 0;
      for ( size_t iP = !isRShrinkedForAdjacent; iP < nodeDataForAdjacent.size(); ++iP )
      {
        double lenFromPar1 =
          GCPnts_AbscissaPoint::Length( curve, par1, nodeDataForAdjacent[iP].param );
        nodeDataForAdjacent[iP].normParam = deltaR + normDelta * lenFromPar1 / shrinkLen;
      }
      // concatenate nodeDataVec and nodeDataForAdjacent
      nodeDataVec.insert(( isRShrinkedForAdjacent ? nodeDataVec.end() : nodeDataVec.begin() ),
                         nodeDataForAdjacent.begin(), nodeDataForAdjacent.end() );
    }

    // Extend nodeDataVec by a node located at the end of not shared _LayerEdge
    /*      n - to add to nodeDataVec
     *      o-----o---
     *      |\    |
     *      | o---o---
     *      | |x--x--- L2
     *      | /
     *      |/ L
     *      x
     *     /    */
    for ( int isR = 0; isR < 2; ++isR )
    {
      _PolyLine& L2 = *( isR ? L._rightLine : L._leftLine ); // line with layers
      if ( ! L2._advancable || L.IsCommonEdgeShared( L2 ) )
        continue;
      vector< const SMDS_MeshNode* >& layerNodes2 = isR ? L2._leftNodes : L2._rightNodes;
      _LayerEdge& LE2 = isR ? L2._lEdges.front() : L2._lEdges.back();
      if ( layerNodes2.empty() )
      {
        // refine the not shared _LayerEdge
        vector<double> layersHeight;
        calcLayersHeight( LE2._length2D, layersHeight, getLineHypothesis( L2._index ));

        vector<gp_XY>& nodeUV2 = LE2._uvRefined;
        nodeUV2.resize    ( layersHeight.size() );
        layerNodes2.resize( layersHeight.size() );
        for ( size_t i = 0; i < layersHeight.size(); ++i )
        {
          gp_XY uv = LE2._uvOut + LE2._normal2D * layersHeight[i];
          gp_Pnt p = _surface->Value( uv.X(), uv.Y() );
          nodeUV2    [ i ] = uv;
          layerNodes2[ i ] = _helper.AddNode( p.X(), p.Y(), p.Z(), /*id=*/0, uv.X(), uv.Y() );
        }
      }
      UVPtStruct ptOfNode;
      ptOfNode.u         = LE2._uvRefined.back().X();
      ptOfNode.v         = LE2._uvRefined.back().Y();
      ptOfNode.node      = layerNodes2.back();
      ptOfNode.param     = isR ? ul : uf;
      ptOfNode.normParam = isR ? 1 : 0;

      nodeDataVec.insert(( isR ? nodeDataVec.end() : nodeDataVec.begin() ), ptOfNode );

      // recompute normParam of nodes in nodeDataVec
      newLength = GCPnts_AbscissaPoint::Length( curve,
                                                nodeDataVec.front().param,
                                                nodeDataVec.back().param);
      for ( size_t iP = 1; iP < nodeDataVec.size(); ++iP )
      {
        const double len = GCPnts_AbscissaPoint::Length( curve,
                                                         nodeDataVec.front().param,
                                                         nodeDataVec[iP].param );
        nodeDataVec[iP].normParam = len / newLength;
      }
    }

    // create a proxy sub-mesh containing the moved nodes
    _ProxyMeshOfFace::_EdgeSubMesh* edgeSM = getProxyMesh()->GetEdgeSubMesh( edgeID );
    edgeSM->SetUVPtStructVec( nodeDataVec );

    // set a sub-mesh event listener to remove just created edges when
    // "ViscousLayers2D" hypothesis is modified
    VISCOUS_3D::ToClearSubWithMain( _mesh->GetSubMesh( E ), _face );

  } // loop on _polyLineVec

  return true;
}

//================================================================================
/*!
 * \brief Returns true if there will be a shrinked mesh on EDGE E of FACE adjFace
 *        near VERTEX V
 */
//================================================================================

bool _ViscousBuilder2D::toShrinkForAdjacent( const TopoDS_Face&   adjFace,
                                             const TopoDS_Edge&   E,
                                             const TopoDS_Vertex& V)
{
  if ( _noShrinkVert.count( getMeshDS()->ShapeToIndex( V )) || adjFace.IsNull() )
    return false;

  vector< const StdMeshers_ViscousLayers2D* > hyps;
  vector< TopoDS_Shape >                      hypShapes;
  if ( VISCOUS_2D::findHyps( *_mesh, adjFace, hyps, hypShapes ))
  {
    VISCOUS_2D::_ViscousBuilder2D builder( *_mesh, adjFace, hyps, hypShapes );
    builder._faceSideVec = StdMeshers_FaceSide::GetFaceWires( adjFace, *_mesh, true, _error );
    builder.findEdgesWithLayers();

    PShapeIteratorPtr edgeIt = _helper.GetAncestors( V, *_mesh, TopAbs_EDGE );
    while ( const TopoDS_Shape* edgeAtV = edgeIt->next() )
    {
      if ( !edgeAtV->IsSame( E ) &&
           _helper.IsSubShape( *edgeAtV, adjFace ) &&
           !builder._ignoreShapeIds.count( getMeshDS()->ShapeToIndex( *edgeAtV )))
      {
        return true;
      }
    }
  }
  return false;
}

//================================================================================
/*!
 * \brief Detect the "bare-face reentrant corner" pattern at VERTEX V: myArmEdge
 *        (advancable, on _face) and myBareEdge (non-advancable, on _face) meet
 *        at V; bareFace is the FACE on the other side of myBareEdge. The
 *        pattern matches when bareFace itself carries no viscous layers, and
 *        exactly one of bareFace's own EDGEs at V (other than myBareEdge)
 *        carries viscous layers on ITS other incident FACE -- that EDGE/FACE
 *        become otherArmEdge/otherArmFace. Returns false for anything else
 *        (no matching other arm, or more than one candidate -- more than 2
 *        viscous-layered arms meeting at one VERTEX is out of scope for the
 *        corner bevel and falls back to today's plain-VERTEX behaviour).
 */
//================================================================================

bool _ViscousBuilder2D::findBareCornerPattern( const TopoDS_Vertex& V,
                                               const TopoDS_Edge&   myArmEdge,
                                               const TopoDS_Edge&   myBareEdge,
                                               const TopoDS_Face&   bareFace,
                                               TopoDS_Edge&         otherArmEdge,
                                               TopoDS_Face&         otherArmFace,
                                               double&              otherArmThickness,
                                               int&                 otherArmNbLayers )
{
  if ( bareFace.IsNull() )
    return false;

  // bareFace must itself have no EDGE of its own actually getting viscous
  // layers (not just "no hypothesis findable" -- a hyp assigned to an
  // ancestor shape would be found for bareFace too even with zero of its own
  // EDGEs covered, see hasOwnViscousLayers()).
  if ( VISCOUS_2D::hasOwnViscousLayers( *_mesh, bareFace ))
    return false;

  bool found = false;

  // Walk bareFace's own boundary EDGEs at V (besides myBareEdge); for each,
  // hop across it (via findArmAcrossEdge, excluding bareFace itself) to find
  // a FACE that has, at V, ANOTHER edge (besides the one just crossed)
  // carrying viscous layers on itself -- e.g. from F1: myBareEdge=E12,
  // candidate=E23 (bareFace=F2's other edge at V), hop across E23 into F3,
  // and within F3's own wire find E3 (F3's other edge at V) carrying VL.
  PShapeIteratorPtr edgeIt = _helper.GetAncestors( V, *_mesh, TopAbs_EDGE );
  while ( const TopoDS_Shape* edgeAtV = edgeIt->next() )
  {
    const TopoDS_Edge& E = TopoDS::Edge( *edgeAtV );
    if ( E.IsSame( myBareEdge ) || E.IsSame( myArmEdge ))
      continue;
    if ( !_helper.IsSubShape( E, bareFace ))
      continue; // consider only bareFace's own boundary EDGEs at V

    TopoDS_Edge candidateArmEdge;
    TopoDS_Face candidateArmFace;
    int         candidateArmNbLayers = 0;
    if ( !findArmAcrossEdge( V, E, bareFace, candidateArmEdge, candidateArmFace, candidateArmNbLayers ))
      continue;

    if ( found )
      return false; // a second candidate arm: >2-arm corner, out of scope

    // findArmAcrossEdge() already confirmed candidateArmEdge carries VL on
    // candidateArmFace; look the hypothesis back up just to read its
    // thickness/layer count.
    vector< const StdMeshers_ViscousLayers2D* > hyps;
    vector< TopoDS_Shape >                      hypShapes;
    VISCOUS_2D::findHyps( *_mesh, candidateArmFace, hyps, hypShapes );
    for ( size_t i = 0; i < hyps.size(); ++i )
    {
      set<int> ignoreIds;
      VISCOUS_2D::getEdgesToIgnore( hyps[i], candidateArmFace, getMeshDS(), ignoreIds );
      if ( !ignoreIds.count( getMeshDS()->ShapeToIndex( candidateArmEdge )))
      {
        otherArmEdge      = candidateArmEdge;
        otherArmFace      = candidateArmFace;
        otherArmThickness = hyps[i]->GetTotalThickness();
        otherArmNbLayers  = hyps[i]->GetNumberLayers();
        found = true;
        break;
      }
    }
  }
  return found;
}

//================================================================================
/*!
 * \brief Bevel a reentrant corner shared by two independently viscous-layered
 *        FACEs (this FACE, "the arm", and another one) through an intervening
 *        bare FACE, instead of letting both sides' shrunk boundaries converge
 *        on the plain VERTEX V. If V matches the pattern (see
 *        findBareCornerPattern()), this ADDS explicit 2D elements covering
 *        the EXTRA wedge beyond V (between myBareEdge's own direction and the
 *        bisector B13, out to the bevel point) -- it never touches V's own
 *        natural interior angle, which the caller (shrink()) leaves for this
 *        FACE's regular (non-corner) fill to cover exactly as it always did
 *        (nodeDataVec/myBareEdge's own discretization near V is left
 *        untouched on purpose). Those two pieces -- V's natural angle and
 *        this extra wedge -- are different, non-overlapping angular ranges
 *        from V, so leaving one alone and adding the other, independently,
 *        is what makes them meet with no gap and no double coverage.
 *        Specifically, this:
 *         - computes (once, cached in a per-VERTEX _CornerData shared with the
 *           other arm) B13 and the two bevel points P1/P2;
 *         - creates *this* FACE's own bevel node at its point (P1 or P2,
 *           whichever this FACE owns) if not already created by an earlier
 *           call for the same VERTEX;
 *         - tiles the wedge between V and the bevel point with one quad per
 *           layer (a triangle for the first layer);
 *         - closes the remaining strip between the deepest quad and
 *           myBareEdge's REAL nearest node to V (wherever the plain 1D
 *           discretization happens to have put it, regardless of layer
 *           thickness) with one more triangle -- see the NB at that point in
 *           the body for why this strip exists and must be closed here.
 *        No NEW node is ever created on myArmEdge/myBareEdge itself: the
 *        "near" corners are free points, geometrically off the wall (at
 *        V + height*N, N = myArmEdge's own inward growth normal at V), never
 *        added to myBareEdge's own mesh/proxy, so the bare FACE's own
 *        boundary never sees them.
 *        IMPORTANT: each bend point (on B13, at cumulative layer thickness)
 *        is, by construction of pointOnBisectorAtHeight(), always exactly
 *        colinear with V (it is defined as the point on the ray from V along
 *        B13). A "fan" of triangles from V straight through consecutive bend
 *        points is therefore ALWAYS degenerate (zero area, verified
 *        numerically) and must never be used -- it silently leaves the whole
 *        wedge unmeshed; that is why every triangle/quad built here has its
 *        "near" corner(s) off that line instead.
 *        Does nothing (falls back to today's plain-VERTEX behaviour) if V
 *        doesn't match the pattern, or the corner is geometrically degenerate.
 *        The final (V,P1,P2) triangle stitching both arms' wedges together is
 *        built separately, once both arms are known to be done, by
 *        StdMeshers_ViscousLayers2D::Compute()'s bare-FACE branch.
 *  \return non-NULL if V matched the pattern and the wedge was built (caller
 *          then skips its own wall-shrink node creation on myBareEdge for
 *          this VERTEX end), NULL if V doesn't match the pattern or the
 *          corner is degenerate (caller falls back to plain-VERTEX
 *          behaviour).
 */
//================================================================================

const SMDS_MeshNode* _ViscousBuilder2D::buildCornerBevel( const TopoDS_Vertex& V,
                                          const TopoDS_Edge&                     myArmEdge,
                                          const TopoDS_Edge&                     myBareEdge,
                                          const TopoDS_Face&                     bareFace,
                                          const SMDS_MeshNode*                   vertexNode,
                                          const THypVL*                          myHyp,
                                          const SMDS_MeshNode*&                  outNearCornerNode )
{
  outNearCornerNode = 0;
  if ( !myHyp || myHyp->GetNumberLayers() <= 0 )
    return 0;

  TopoDS_Edge otherArmEdge;
  TopoDS_Face otherArmFace;
  double      otherArmThickness = 0;
  int         otherArmNbLayers  = 0;
  if ( !findBareCornerPattern( V, myArmEdge, myBareEdge, bareFace,
                               otherArmEdge, otherArmFace,
                               otherArmThickness, otherArmNbLayers ))
    return 0;

  const double myThickness = myHyp->GetTotalThickness();

  gp_Pnt pv; gp_Vec T, N, faceNorm;
  if ( !VISCOUS_2D::edgeTangentAndNormalAt( V, myArmEdge, _face, pv, T, N, &faceNorm ))
    return 0; // can't build a local frame: fall back to today's plain-VERTEX behaviour

  _CornerData* cd = VISCOUS_2D::_CornerDataHolder::GetCornerData( V, *_mesh );
  const int myFaceID = getMeshDS()->ShapeToIndex( _face );

  if ( !cd->_isValid )
  {
    gp_Dir B13; gp_Pnt p1, p2;
    if ( !StdMeshers_ViscousLayers2D::ComputeCornerBisector(
           V, myArmEdge, _face, myThickness,
           otherArmEdge, otherArmFace, otherArmThickness, B13, p1, p2 ))
      return 0; // degenerate angle: fall back to today's plain-VERTEX behaviour

    cd->_isValid    = true;
    cd->_B13        = B13;
    cd->_p1Point    = p1;
    cd->_p2Point    = p2;
    cd->_armFaceID1 = myFaceID;                                  // I always call
    cd->_armFaceID3 = getMeshDS()->ShapeToIndex( otherArmFace );  // myself as "arm 1"

    // Make sure this VERTEX's sub-mesh is cleared (invalidating this cached
    // _CornerData) whenever either arm FACE is recomputed, mirroring what
    // RestoreListeners() already does for the shared EDGEs.
    VISCOUS_3D::ToClearSubWithMain( _mesh->GetSubMesh( V ), _face );
    VISCOUS_3D::ToClearSubWithMain( _mesh->GetSubMesh( V ), otherArmFace );
  }

  const bool             iOwnP1 = ( myFaceID == cd->_armFaceID1 );
  const SMDS_MeshNode*&  myNode = iOwnP1 ? cd->_p1Node : cd->_p2Node;
  const gp_Pnt&         myPoint = iOwnP1 ? cd->_p1Point : cd->_p2Point;
  vector< const SMDS_MeshNode* >&       myBendChain = iOwnP1 ? cd->_bendChain1 : cd->_bendChain2;
  const vector< const SMDS_MeshNode* >& otherBendChain = iOwnP1 ? cd->_bendChain2 : cd->_bendChain1;

  // When both arms share equal thickness/layer counts (the common symmetric
  // case), every bend point -- including the deepest, P1/P2 -- coincides
  // spatially with the OTHER arm's corresponding point. Reuse the other arm's
  // already-created node instead of creating a spatially-duplicate one: two
  // distinct node objects at the same position would otherwise (a) show up as
  // a "double node" defect mesh-wide, and (b) make tryBuildBareCornerProxy()'s
  // (V,P1,P2) closing triangle degenerate (zero area) without collapsing it.
  auto reuseIfCoincident = [&]( const gp_Pnt& p )->const SMDS_MeshNode*
  {
    for ( size_t i = 0; i < otherBendChain.size(); ++i )
      if ( otherBendChain[i] && ( SMESH_TNodeXYZ( otherBendChain[i] ) - p.XYZ() ).Modulus() <= Precision::Confusion() )
        return otherBendChain[i];
    return 0;
  };

  if ( !myNode )
  {
    myNode = reuseIfCoincident( myPoint );
    if ( !myNode )
      myNode = _helper.AddNode( myPoint.X(), myPoint.Y(), myPoint.Z() );
  }
  ( iOwnP1 ? cd->_p1NodeID : cd->_p2NodeID ) = myNode->GetID();

  // Each individual layer bends at its OWN point on B13, farther from V the
  // deeper the layer, at the SAME cumulative thickness this layer would reach
  // on a plain (non-corner) wall -- computed analytically by calcLayersHeight,
  // never by measuring a node actually placed on the wall.
  vector<double> heights;
  calcLayersHeight( myThickness, heights, myHyp );

  vector< const SMDS_MeshNode* > bendNodes( heights.size() );
  bendNodes.back() = myNode; // deepest layer's bend point is the shared node

  for ( size_t i = 0; i + 1 < heights.size(); ++i )
  {
    gp_Pnt Pi;
    if ( VISCOUS_2D::pointOnBisectorAtHeight( pv, T, N, cd->_B13, heights[i], Pi ))
    {
      bendNodes[i] = reuseIfCoincident( Pi );
      if ( !bendNodes[i] )
        bendNodes[i] = _helper.AddNode( Pi.X(), Pi.Y(), Pi.Z() );
    }
    else
      bendNodes[i] = myNode; // degenerate fallback: collapse to the deepest point
  }
  myBendChain = bendNodes; // let the OTHER arm (if it runs after us) reuse these

  // "Near" corners, one per layer boundary, off B13: Qi = V + height*N
  // (N = myArmEdge's own inward growth normal at V, into THIS FACE) -- a
  // purely analytic, free point (never added to any EDGE's mesh/proxy), off
  // the wall entirely, strictly inside this FACE's own territory. These
  // close the wedge between myArmEdge's own regular (collapsed-onto-V)
  // per-layer bands and B13, which those bands (built by refine(), always
  // collapsed onto the plain VERTEX for a corner-bevel end, see
  // isCornerBevelL/R there) do NOT cover on their own -- confirmed: without
  // this wedge, a real hole reopens between the wall and B13 (visible as a
  // triangular gap in a filled-polygon plot, even though the WHOLE-mesh area
  // check can net out close to zero if this
  // wedge's absence happens to be offset by unrelated overlap elsewhere).
  vector< const SMDS_MeshNode* > nearNodes( heights.size() );
  for ( size_t i = 0; i < heights.size(); ++i )
  {
    const gp_Pnt Qi = pv.XYZ() + heights[i] * N.XYZ();
    nearNodes[i] = _helper.AddNode( Qi.X(), Qi.Y(), Qi.Z() );
  }
  ( iOwnP1 ? cd->_nearChain1 : cd->_nearChain2 ) = nearNodes;
  outNearCornerNode = nearNodes.back();

  const bool isReverse =
    ( _helper.GetSubShapeOri( _mesh->GetShapeToMesh(), _face ) == TopAbs_REVERSED );
  const gp_Vec reqNormal = isReverse ? -faceNorm : faceNorm;

  auto addOriented3 = [&]( const SMDS_MeshNode* a,
                           const SMDS_MeshNode* b,
                           const SMDS_MeshNode* c )
  {
    const gp_Pnt pa = SMESH_TNodeXYZ( a ), pb = SMESH_TNodeXYZ( b ), pc = SMESH_TNodeXYZ( c );
    const gp_Vec actual = gp_Vec( pa, pb ) ^ gp_Vec( pa, pc );
    if ( actual.Dot( reqNormal ) >= 0 )
      _helper.AddFace( a, b, c );
    else
      _helper.AddFace( a, c, b );
  };
  auto addOriented4 = [&]( const SMDS_MeshNode* a, const SMDS_MeshNode* b,
                           const SMDS_MeshNode* c, const SMDS_MeshNode* d )
  {
    const gp_Pnt pa = SMESH_TNodeXYZ( a ), pb = SMESH_TNodeXYZ( b ), pc = SMESH_TNodeXYZ( c );
    const gp_Vec actual = gp_Vec( pa, pb ) ^ gp_Vec( pa, pc );
    if ( actual.Dot( reqNormal ) >= 0 )
      _helper.AddFace( a, b, c, d );
    else
      _helper.AddFace( a, d, c, b );
  };

  // Layer 0's band is the triangle (V, nearNodes[0], bendNodes[0]) -- V has
  // no "previous" near/bend point. Band i>0 (between layer i-1 and layer i)
  // is the quad (nearNodes[i-1], nearNodes[i], bendNodes[i], bendNodes[i-1]).
  addOriented3( vertexNode, nearNodes[0], bendNodes[0] );
  for ( size_t i = 1; i < bendNodes.size(); ++i )
    addOriented4( nearNodes[i-1], nearNodes[i], bendNodes[i], bendNodes[i-1] );

  // Close the strip between nearNodes.back() (cumulative layer thickness from
  // V) and myBareEdge's REAL nearest node to V (untouched by this whole
  // feature, still wherever the plain 1D discretization put it -- possibly
  // much farther from V than the total layer thickness) with one more
  // triangle -- see the equivalent NB in the wall-shrink caller (shrink()).
  const SMDS_MeshNode* wallNeighbor = 0;
  {
    const int bareEdgeID = getMeshDS()->ShapeToIndex( myBareEdge );
    SMDS_ElemIteratorPtr segIt = vertexNode->GetInverseElementIterator( SMDSAbs_Edge );
    while ( segIt->more() )
    {
      const SMDS_MeshElement* seg = segIt->next();
      if ( seg->getshapeId() != bareEdgeID )
        continue;
      for ( int i = 0; i < seg->NbNodes(); ++i )
        if ( seg->GetNode( i ) != vertexNode )
        {
          wallNeighbor = seg->GetNode( i );
          break;
        }
      break;
    }
  }
  if ( wallNeighbor )
  {
    addOriented3( nearNodes.back(), wallNeighbor, myNode );
    // Close the OTHER remaining sliver, between V itself, wallNeighbor and
    // nearNodes.back(): now that this arm's OWN proxy submesh for myBareEdge
    // is spliced (by the caller, shrink()) to jump straight from wallNeighbor
    // to nearNodes.back() instead of passing through V, the underlying 2D
    // algorithm's core fill for THIS FACE no longer reaches V at all near
    // this corner -- V, nearNodes[0..] are exactly colinear (all lie on the
    // ray V+height*N), so V-nearNodes.back() is already a real chain of mesh
    // edges (the wedge's own triangle/quads above), making (V, wallNeighbor,
    // nearNodes.back()) the exact remaining triangle, degenerate (zero area,
    // harmless) only in the special case where myBareEdge happens to be
    // colinear with N itself (e.g. a plain right-angle corner).
    addOriented3( vertexNode, wallNeighbor, nearNodes.back() );
  }

  // The caller only checks this return value for null/non-null (whether the
  // pattern matched, i.e. whether to skip its own wall-shrink node creation
  // on myBareEdge for this VERTEX end) -- any non-null node works.
  return myNode;
}

//================================================================================
/*!
 * \brief Look, among EDGEs of _face (a bare FACE with no viscous layers of its
 *        own) incident to VERTEX V, across \a bareEdge (also incident to V, one
 *        of _face's own boundary EDGEs) for the FACE on its other side; if that
 *        FACE has, at V, another EDGE (besides bareEdge) carrying viscous
 *        layers on itself, return that EDGE/FACE as the "arm" reached via
 *        bareEdge. This is the mirror of findBareCornerPattern(), starting
 *        from the bare FACE's own side instead of from an arm.
 */
//================================================================================

bool _ViscousBuilder2D::findArmAcrossEdge( const TopoDS_Vertex& V,
                                           const TopoDS_Edge&   bareEdge,
                                           const TopoDS_Face&   excludeFace,
                                           TopoDS_Edge&         armEdge,
                                           TopoDS_Face&         armFace,
                                           int&                 armNbLayers )
{
  PShapeIteratorPtr faceIt = _helper.GetAncestors( bareEdge, *_mesh, TopAbs_FACE );
  while ( const TopoDS_Shape* f = faceIt->next() )
  {
    const TopoDS_Face& F = TopoDS::Face( *f );
    if ( F.IsSame( excludeFace ) || F.IsSame( _face ))
      continue;

    PShapeIteratorPtr edgeIt = _helper.GetAncestors( V, *_mesh, TopAbs_EDGE );
    while ( const TopoDS_Shape* e2 = edgeIt->next() )
    {
      const TopoDS_Edge& E = TopoDS::Edge( *e2 );
      if ( E.IsSame( bareEdge ) || !_helper.IsSubShape( E, F ))
        continue;

      vector< const StdMeshers_ViscousLayers2D* > hyps;
      vector< TopoDS_Shape >                      hypShapes;
      if ( !VISCOUS_2D::findHyps( *_mesh, F, hyps, hypShapes ))
        continue;

      for ( size_t i = 0; i < hyps.size(); ++i )
      {
        set<int> ignoreIds;
        VISCOUS_2D::getEdgesToIgnore( hyps[i], F, getMeshDS(), ignoreIds );
        if ( !ignoreIds.count( getMeshDS()->ShapeToIndex( E )))
        {
          armEdge     = E;
          armFace     = F;
          armNbLayers = hyps[i]->GetNumberLayers();
          return true;
        }
      }
    }
  }
  return false;
}

//================================================================================
/*!
 * \brief Rebuild EDGE iE's own proxy sub-mesh (of _face's wire \a wire) so it
 *        ends (cornerAtEnd) or starts (!cornerAtEnd) at \a cornerNode instead
 *        of at the plain shared VERTEX, splicing in the corner bevel computed
 *        by the two arm FACEs (see buildCornerBevel()). \a points is the
 *        already-computed whole-wire StdMeshers_FaceSide::GetUVPtStruct(),
 *        used to read back EDGE iE's own real point range before splicing.
 */
//================================================================================

void _ViscousBuilder2D::spliceCornerIntoEdgeProxy( StdMeshers_FaceSidePtr     wire,
                                                   int                        iE,
                                                   const vector<UVPtStruct>&  points,
                                                   bool                       cornerAtEnd,
                                                   const SMDS_MeshNode*       cornerNode,
                                                   int                        numWallNodesToExclude )
{
  if ( !cornerNode )
    return;

  const double parF = wire->FirstParameter( iE );
  const double parL = wire->LastParameter ( iE );
  if ( parL - parF <= 0 )
    return;

  int iFirst = 0;
  while ( iFirst < (int)points.size() && points[ iFirst ].normParam < parF - 1e-10 )
    ++iFirst;
  int iLast = iFirst;
  while ( iLast+1 < (int)points.size() && points[ iLast ].normParam < parL - 1e-10 )
    ++iLast;
  if ( iLast <= iFirst )
    return;

  UVPtStructVec edgePoints( points.begin()+iFirst, points.begin()+iLast+1 );

  // Drop the VERTEX itself and replace it with the bevel node. The arm's own
  // near-corner layer nodes live in *that arm's own* proxy submesh for this
  // EDGE, a per-FACE override this bare FACE's wire never reads (confirmed
  // by inspection: this FACE's own point list for this EDGE is the plain,
  // untouched regular discretization, with no trace of the arm's layer
  // nodes) -- so there is nothing of the arm's to skip past here beyond V
  // itself.
  const int nExcl = Min( numWallNodesToExclude, (int)edgePoints.size() - 1 );
  if ( nExcl > 0 )
  {
    if ( cornerAtEnd )
      edgePoints.erase( edgePoints.end() - nExcl, edgePoints.end() );
    else
      edgePoints.erase( edgePoints.begin(), edgePoints.begin() + nExcl );
  }

  UVPtStruct cornerPt;
  cornerPt.node = cornerNode;
  // cornerNode is a free node (added via _helper.AddNode(x,y,z), not
  // positioned on any sub-shape), so plain GetNodeUV() has no EDGE/VERTEX/
  // FACE position to read and falls back to a cold-start surface
  // projection that can fail to converge and silently return (0,0) -- a
  // wildly wrong point in this FACE's parametric space that corrupts
  // NETGEN's 2D input (confirmed: this is what produced the long sliver
  // triangles/gaps seen in practice). Seed the projection with the
  // immediately adjacent, already-correctly-positioned wire point's own UV
  // instead of starting cold.
  const UVPtStruct& neighbourPt =
    cornerAtEnd ? edgePoints.back() : edgePoints.front();
  const gp_Pnt2d seedUV( neighbourPt.u, neighbourPt.v );
  const gp_Pnt cornerXYZ = SMESH_TNodeXYZ( cornerNode );
  const gp_Pnt2d cornerUV = _helper.GetSurface( _face )->NextValueOfUV
    ( seedUV, cornerXYZ, BRep_Tool::Tolerance( _face ));
  cornerPt.u     = cornerUV.X();
  cornerPt.v     = cornerUV.Y();
  cornerPt.param = cornerAtEnd ? edgePoints.back().param : edgePoints.front().param;
  cornerPt.normParam = 0; // recomputed below

  if ( cornerAtEnd )
    edgePoints.push_back( cornerPt );
  else
    edgePoints.insert( edgePoints.begin(), cornerPt );

  // Renormalize normParam in [0,1] over the (now corner-extended) point list,
  // by cumulative 3D chord length -- the same length-based normalization
  // convention used elsewhere in this file (e.g. refine()'s segLen).
  vector<double> len( edgePoints.size(), 0.0 );
  for ( size_t i = 1; i < edgePoints.size(); ++i )
  {
    gp_Pnt p0 = SMESH_TNodeXYZ( edgePoints[i-1].node );
    gp_Pnt p1 = SMESH_TNodeXYZ( edgePoints[i  ].node );
    len[i] = len[i-1] + p0.Distance( p1 );
  }
  const double totalLen = len.back();
  if ( totalLen > gp::Resolution() )
    for ( size_t i = 0; i < edgePoints.size(); ++i )
      edgePoints[i].normParam = len[i] / totalLen;

  getProxyMesh()->GetEdgeSubMesh( wire->EdgeID( iE ))->SetUVPtStructVec( edgePoints );
}

//================================================================================
/*!
 * \brief For this bare FACE (no viscous layers of its own), find any VERTEX of
 *        its own wire matching the 2-arm/1-bare reentrant-corner pattern (see
 *        findArmAcrossEdge()) and splice the corner bevel -- computed by the
 *        two arm FACEs' own buildCornerBevel(), forced to run here if they
 *        haven't yet -- into this FACE's own proxy mesh, so its own 2D
 *        algorithm routes through P1/P2 instead of the plain VERTEX. Also
 *        creates the (VERTEX,P1,P2) triangle stitching both arms' fans
 *        together, once, since this bare FACE is the one point guaranteed to
 *        run only after both arms are done.
 *  \return true if a proxy mesh was actually built for at least one VERTEX
 */
//================================================================================

bool _ViscousBuilder2D::tryBuildBareCornerProxy()
{
  _faceSideVec = StdMeshers_FaceSide::GetFaceWires( _face, *_mesh, true, _error, &_helper );
  if ( !_error->IsOK() )
    return false;

  bool builtAny = false;
  for ( size_t iWire = 0; iWire < _faceSideVec.size(); ++iWire )
  {
    StdMeshers_FaceSidePtr wire = _faceSideVec[ iWire ];
    const vector<UVPtStruct>& points = wire->GetUVPtStruct();
    if ( points.empty() )
      continue;
    const int nbE = wire->NbEdges();
    for ( int iE = 0; iE < nbE; ++iE )
    {
      const int             iE2        = ( iE + 1 ) % nbE;
      const TopoDS_Edge&    edgeBefore = wire->Edge( iE );
      const TopoDS_Edge&    edgeAfter  = wire->Edge( iE2 );
      const TopoDS_Vertex&  V          = wire->LastVertex( iE );

      TopoDS_Edge armEdgeBefore, armEdgeAfter;
      TopoDS_Face armFaceBefore, armFaceAfter;
      int         armNbLayersBefore = 0, armNbLayersAfter = 0;
      if ( !findArmAcrossEdge( V, edgeBefore, _face, armEdgeBefore, armFaceBefore, armNbLayersBefore ) ||
           !findArmAcrossEdge( V, edgeAfter,  _face, armEdgeAfter,  armFaceAfter,  armNbLayersAfter ))
        continue; // not the 2-arm/1-bare pattern at this VERTEX

      _CornerData* cd = VISCOUS_2D::_CornerDataHolder::GetCornerData( V, *_mesh );
      if ( !cd->_p1Node || !cd->_p2Node )
      {
        // Force the two arms' shrink() to run (idempotent: StdMeshers_
        // ViscousLayers2D::Compute() is itself cached) so _CornerData gets
        // populated. NB gated on the NODEs, not cd->_isValid: one arm may
        // have already run (independently of this bare FACE's own Compute())
        // and made _isValid true while only creating its OWN node -- the
        // other arm still needs forcing.
        StdMeshers_ViscousLayers2D::Compute( *_mesh, armFaceBefore );
        StdMeshers_ViscousLayers2D::Compute( *_mesh, armFaceAfter );
      }
      if ( !cd->_isValid || !cd->_p1Node || !cd->_p2Node )
        continue; // degenerate angle, or an arm failed: leave plain VERTEX

      const int idBefore = getMeshDS()->ShapeToIndex( armFaceBefore );
      const SMDS_MeshNode* nodeBefore =
        ( idBefore == cd->_armFaceID1 ) ? cd->_p1Node : cd->_p2Node;
      const int idAfter = getMeshDS()->ShapeToIndex( armFaceAfter );
      const SMDS_MeshNode* nodeAfter =
        ( idAfter == cd->_armFaceID1 ) ? cd->_p1Node : cd->_p2Node;

      // Exclude only the VERTEX itself and replace it with the bevel node.
      // The arm's own near-corner layer nodes (created by that arm's own
      // shrink()) live in *that arm's own* proxy submesh for this EDGE, a
      // per-FACE override this bare FACE's wire never reads (proxies are
      // per-consuming-FACE, confirmed by inspection: this FACE's own
      // GetUVPtStruct() for this EDGE returns the plain, untouched regular
      // discretization, with no trace of the arm's layer nodes) -- so there
      // is nothing of the arm's to skip past here beyond V itself. Excluding
      // more than that (as done previously) cuts into this FACE's own
      // legitimate regular boundary points instead, leaving a large gap.
      spliceCornerIntoEdgeProxy( wire, iE,  points, /*cornerAtEnd=*/true,  nodeBefore, 1 );
      spliceCornerIntoEdgeProxy( wire, iE2, points, /*cornerAtEnd=*/false, nodeAfter,  1 );

      if ( !cd->_wedgeBuilt && nodeBefore == nodeAfter )
      {
        // Both arms' deepest bend node coincide (equal thickness/layers on
        // both arms, the common symmetric case) -- the corner is already
        // fully closed by the two arms' own wedges meeting at that single
        // shared node; a (V,P1,P2) triangle here would be degenerate
        // (zero area) since P1==P2, so there is nothing left to build.
        cd->_wedgeBuilt = true;
      }
      if ( !cd->_wedgeBuilt )
      {
        const SMDS_MeshNode* vNode = 0;
        if ( SMESH_subMesh* vSM = _mesh->GetSubMesh( V ))
          if ( SMESHDS_SubMesh* smDS = vSM->GetSubMeshDS() )
          {
            SMDS_NodeIteratorPtr nIt = smDS->GetNodes();
            if ( nIt->more() )
              vNode = nIt->next();
          }
        if ( vNode )
        {
          // Orient geometrically rather than assuming nodeBefore/nodeAfter's
          // wire order maps directly to a consistent winding (buildCornerBevel
          // found the analogous fixed-order assumption to be correct for one
          // arm of a corner and reversed for the other).
          gp_Pnt pv; gp_Vec T, N, faceNorm;
          const bool haveFrame =
            VISCOUS_2D::edgeTangentAndNormalAt( V, edgeBefore, _face, pv, T, N, &faceNorm );
          const bool isReverse =
            ( _helper.GetSubShapeOri( _mesh->GetShapeToMesh(), _face ) == TopAbs_REVERSED );
          const gp_Vec reqNormal = isReverse ? -faceNorm : faceNorm;

          const gp_Pnt pV = SMESH_TNodeXYZ( vNode );
          const gp_Pnt pB = SMESH_TNodeXYZ( nodeBefore );
          const gp_Pnt pA = SMESH_TNodeXYZ( nodeAfter );
          const gp_Vec actual = gp_Vec( pV, pB ) ^ gp_Vec( pV, pA );
          if ( !haveFrame || actual.Dot( reqNormal ) >= 0 )
            _helper.AddFace( vNode, nodeBefore, nodeAfter );
          else
            _helper.AddFace( vNode, nodeAfter, nodeBefore );
          cd->_wedgeBuilt = true;
        }
      }
      builtAny = true;
    }
  }
  return builtAny;
}

//================================================================================
/*!
 * \brief Make faces
 */
//================================================================================

bool _ViscousBuilder2D::refine()
{
  // find out orientation of faces to create
  bool isReverse = 
    ( _helper.GetSubShapeOri( _mesh->GetShapeToMesh(), _face ) == TopAbs_REVERSED );

  // store a proxyMesh in a sub-mesh
  // make faces on each _PolyLine
  vector< double > layersHeight;
  //double prevLen2D = -1;
  for ( size_t iL = 0; iL < _polyLineVec.size(); ++iL )
  {
    _PolyLine& L = _polyLineVec[ iL ];
    if ( !L._advancable ) continue;

    // replace an inactive (1st) _LayerEdge with an active one of a neighbour _PolyLine
    //size_t iLE = 0, nbLE = L._lEdges.size();
    const bool leftEdgeShared  = L.IsCommonEdgeShared( *L._leftLine );
    const bool rightEdgeShared = L.IsCommonEdgeShared( *L._rightLine );
    if ( /*!L._leftLine->_advancable &&*/ leftEdgeShared )
    {
      L._lEdges[0] = L._leftLine->_lEdges.back();
      //iLE += int( !L._leftLine->_advancable );
    }
    if ( !L._rightLine->_advancable && rightEdgeShared )
    {
      L._lEdges.back() = L._rightLine->_lEdges[0];
      //--nbLE;
    }

    // limit length of neighbour _LayerEdge's to avoid sharp change of layers thickness

    vector< double > segLen( L._lEdges.size() );
    segLen[0] = 0.0;

    // check if length modification is useful: look for _LayerEdge's
    // with length limited due to collisions
    bool lenLimited = false;
    for ( size_t iLE = 1; ( iLE < L._lEdges.size()-1 && !lenLimited ); ++iLE )
      lenLimited = L._lEdges[ iLE ]._isBlocked;

    if ( lenLimited )
    {
      for ( size_t i = 1; i < segLen.size(); ++i )
      {
        // accumulate length of segments
        double sLen = (L._lEdges[i-1]._uvOut - L._lEdges[i]._uvOut ).Modulus();
        segLen[i] = segLen[i-1] + sLen;
      }
      const double totSegLen = segLen.back();
      // normalize the accumulated length
      for ( size_t iS = 1; iS < segLen.size(); ++iS )
        segLen[iS] /= totSegLen;

      for ( int isR = 0; isR < 2; ++isR )
      {
        size_t iF = 0, iL = L._lEdges.size()-1;
        size_t *i = isR ? &iL : &iF;
        _LayerEdge* prevLE = & L._lEdges[ *i ];
        double weight = 0;
        for ( ++iF, --iL; iF < L._lEdges.size()-1; ++iF, --iL )
        {
          _LayerEdge& LE = L._lEdges[*i];
          if ( prevLE->_length2D > 0 )
          {
            gp_XY tangent ( LE._normal2D.Y(), -LE._normal2D.X() );
            weight += Abs( tangent * ( prevLE->_uvIn - LE._uvIn )) / totSegLen;
            // gp_XY prevTang( LE._uvOut - prevLE->_uvOut );
            // gp_XY prevNorm( -prevTang.Y(), prevTang.X() );
            gp_XY prevNorm = LE._normal2D;
            double prevProj = prevNorm * ( prevLE->_uvIn - prevLE->_uvOut );
            if ( prevProj > 0 ) {
              prevProj /= prevNorm.Modulus();
              if ( LE._length2D < prevProj )
                weight += 0.75 * ( 1 - weight ); // length decrease is more preferable
              LE._length2D  = weight * LE._length2D + ( 1 - weight ) * prevProj;
              LE._uvIn = LE._uvOut + LE._normal2D * LE._length2D;
            }
          }
          prevLE = & LE;
        }
      }
    }
    // DEBUG:  to see _uvRefined. cout can be redirected to hide NETGEN output
    // cerr << "import smesh" << endl << "mesh = smesh.Mesh()"<< endl;

    const vector<UVPtStruct>& points = L._wire->GetUVPtStruct();

    // analyse extremities of the _PolyLine to find existing nodes
    const TopoDS_Vertex&  V1 = L._wire->FirstVertex( L._edgeInd );
    const TopoDS_Vertex&  V2 = L._wire->LastVertex ( L._edgeInd );
    const int           v1ID = getMeshDS()->ShapeToIndex( V1 );
    const int           v2ID = getMeshDS()->ShapeToIndex( V2 );
    const bool isShrinkableL = ! _noShrinkVert.count( v1ID );
    const bool isShrinkableR = ! _noShrinkVert.count( v2ID );

    // Reentrant-corner bevel: V1/V2 may be a bevel VERTEX exempted from
    // _noShrinkVert purely so shrink()'s own buildCornerBevel() branch can
    // run (see the exemption in findEdgesWithLayers()). This ARM's own
    // official per-layer band reaches, at such a VERTEX, the SAME near-chain
    // node buildCornerBevel() computed (Qi = V + height*N along this ARM's
    // own inward normal at V) instead of collapsing onto the plain VERTEX --
    // safe (unlike two earlier, reverted attempts at this exact substitution,
    // see git history / project memory) ONLY because shrink()'s call site
    // that invokes buildCornerBevel() now ALSO splices this same near-chain
    // into THIS FACE's own proxy sub-mesh for its bare-adjacent EDGE, so the
    // underlying 2D algorithm's core fill at V no longer reaches V either --
    // both retreat to the same node, so there is nothing left to double-cover.
    const int myFaceID = getMeshDS()->ShapeToIndex( _face );
    bool isCornerBevelL = false, isCornerBevelR = false;
    const vector< const SMDS_MeshNode* >* nearChainL = 0;
    const vector< const SMDS_MeshNode* >* nearChainR = 0;
    if ( _CornerData* cd1 = VISCOUS_2D::_CornerDataHolder::FindCornerData( V1, *_mesh ))
      if (( isCornerBevelL = cd1->_isValid ))
        nearChainL = ( myFaceID == cd1->_armFaceID1 ) ? &cd1->_nearChain1 : &cd1->_nearChain2;
    if ( V2.IsSame( V1 ))
    {
      isCornerBevelR = isCornerBevelL;
      nearChainR     = nearChainL;
    }
    else if ( _CornerData* cd2 = VISCOUS_2D::_CornerDataHolder::FindCornerData( V2, *_mesh ))
      if (( isCornerBevelR = cd2->_isValid ))
        nearChainR = ( myFaceID == cd2->_armFaceID1 ) ? &cd2->_nearChain1 : &cd2->_nearChain2;
    if ( nearChainL && nearChainL->empty() ) { isCornerBevelL = false; nearChainL = 0; }
    if ( nearChainR && nearChainR->empty() ) { isCornerBevelR = false; nearChainR = 0; }

    bool hasLeftNode     = ( !L._leftLine->_rightNodes.empty() && leftEdgeShared  );
    bool hasRightNode    = ( !L._rightLine->_leftNodes.empty() && rightEdgeShared );
    bool hasOwnLeftNode  = ( !L._leftNodes.empty() );
    bool hasOwnRightNode = ( !L._rightNodes.empty() );
    bool isClosedEdge    = ( points[ L._firstPntInd ].node == points[ L._lastPntInd ].node );
    const size_t
      nbN = L._lastPntInd - L._firstPntInd + 1,
      iN0 = ( hasLeftNode || hasOwnLeftNode || isClosedEdge || !isShrinkableL || isCornerBevelL ),
      iNE = nbN - ( hasRightNode || hasOwnRightNode || !isShrinkableR || isCornerBevelR );

    // update _uvIn of end _LayerEdge's by existing nodes
    const SMDS_MeshNode *nL = 0, *nR = 0;
    if ( hasOwnLeftNode )    nL = L._leftNodes.back();
    else if ( hasLeftNode )  nL = L._leftLine->_rightNodes.back();
    if ( hasOwnRightNode )   nR = L._rightNodes.back();
    else if ( hasRightNode ) nR = L._rightLine->_leftNodes.back();
    if ( nL )
      L._lEdges[0]._uvIn = _helper.GetNodeUV( _face, nL, points[ L._firstPntInd + 1 ].node );
    if ( nR )
      L._lEdges.back()._uvIn = _helper.GetNodeUV( _face, nR, points[ L._lastPntInd - 1 ].node );

    // compute normalized [0;1] node parameters of nodes on a _PolyLine
    vector< double > normPar( nbN );
    const double
      normF    = L._wire->FirstParameter( L._edgeInd ),
      normL    = L._wire->LastParameter ( L._edgeInd ),
      normDist = normL - normF;
    for ( int i = L._firstPntInd; i <= L._lastPntInd; ++i )
      normPar[ i - L._firstPntInd ] = ( points[i].normParam - normF ) / normDist;

    // Calculate UV of most inner nodes

    vector< gp_XY > innerUV( nbN );

    // check if innerUV should be interpolated between _LayerEdge::_uvIn's
    const size_t nbLE = L._lEdges.size();
    bool needInterpol = ( nbN != nbLE );
    if ( !needInterpol )
    {
      // more check: compare length of inner and outer end segments
      double lenIn, lenOut;
      for ( int isR = 0; isR < 2 && !needInterpol; ++isR )
      {
        const _Segment& segIn = isR ? L._segments.back() : L._segments[0];
        const gp_XY&   uvIn1  = segIn.p1();
        const gp_XY&   uvIn2  = segIn.p2();
        const gp_XY&   uvOut1 = L._lEdges[ isR ? nbLE-1 : 0 ]._uvOut;
        const gp_XY&   uvOut2 = L._lEdges[ isR ? nbLE-2 : 1 ]._uvOut;
        if ( _is2DIsotropic )
        {
          lenIn  = ( uvIn1  - uvIn2  ).Modulus();
          lenOut = ( uvOut1 - uvOut2 ).Modulus();
        }
        else
        {
          lenIn  = _surface->Value( uvIn1.X(), uvIn1.Y() )
            .Distance( _surface->Value( uvIn2.X(), uvIn2.Y() ));
          lenOut  = _surface->Value( uvOut1.X(), uvOut1.Y() )
            .Distance( _surface->Value( uvOut2.X(), uvOut2.Y() ));
        }
        needInterpol = ( lenIn < 0.66 * lenOut );
      }
    }

    if ( needInterpol )
    {
      // compute normalized accumulated length of inner segments
      size_t iS;
      if ( _is2DIsotropic )
        for ( iS = 1; iS < segLen.size(); ++iS )
        {
          double sLen = ( L._lEdges[iS-1]._uvIn - L._lEdges[iS]._uvIn ).Modulus();
          segLen[iS] = segLen[iS-1] + sLen;
        }
      else
        for ( iS = 1; iS < segLen.size(); ++iS )
        {
          const gp_XY& uv1 = L._lEdges[iS-1]._uvIn;
          const gp_XY& uv2 = L._lEdges[iS  ]._uvIn;
          gp_Pnt p1 = _surface->Value( uv1.X(), uv1.Y() );
          gp_Pnt p2 = _surface->Value( uv2.X(), uv2.Y() );
          double sLen = p1.Distance( p2 );
          segLen[iS] = segLen[iS-1] + sLen;
        }
      // normalize the accumulated length
      for ( iS = 1; iS < segLen.size(); ++iS )
        segLen[iS] /= segLen.back();

      // calculate UV of most inner nodes according to the normalized node parameters
      iS = 0;
      for ( size_t i = 0; i < innerUV.size(); ++i )
      {
        while ( normPar[i] > segLen[iS+1] )
          ++iS;
        double r = ( normPar[i] - segLen[iS] ) / ( segLen[iS+1] - segLen[iS] );
        innerUV[ i ] = r * L._lEdges[iS+1]._uvIn + (1-r) * L._lEdges[iS]._uvIn;
      }
    }
    else // ! needInterpol
    {
      for ( size_t i = 0; i < nbLE; ++i )
        innerUV[ i ] = L._lEdges[i]._uvIn;
    }

    // normalized height of layers
    const THypVL* hyp = getLineHypothesis( iL );
    calcLayersHeight( 1., layersHeight, hyp);

    // Create layers of faces

    // nodes to create 1 layer of faces
    vector< const SMDS_MeshNode* > outerNodes( nbN );
    vector< const SMDS_MeshNode* > innerNodes( nbN );

    // initialize outerNodes by nodes of the L._wire
    for ( int i = L._firstPntInd; i <= L._lastPntInd; ++i )
      outerNodes[ i-L._firstPntInd ] = points[i].node;

    L._leftNodes .reserve( hyp->GetNumberLayers() );
    L._rightNodes.reserve( hyp->GetNumberLayers() );
    int cur = 0, prev = -1; // to take into account orientation of _face
    if ( isReverse ) std::swap( cur, prev );
    for ( int iF = 0; iF < hyp->GetNumberLayers(); ++iF ) // loop on layers of faces
    {
      // create innerNodes of a current layer
      for ( size_t i = iN0; i < iNE; ++i )
      {
        gp_XY uvOut = points[ L._firstPntInd + i ].UV();
        gp_XY& uvIn = innerUV[ i ];
        gp_XY    uv = layersHeight[ iF ] * uvIn + ( 1.-layersHeight[ iF ]) * uvOut;
        gp_Pnt    p = _surface->Value( uv.X(), uv.Y() );
        innerNodes[i] = _helper.AddNode( p.X(), p.Y(), p.Z(), /*id=*/0, uv.X(), uv.Y() );
      }
      // use nodes created for adjacent _PolyLine's
      if ( hasOwnLeftNode )    innerNodes.front() = L._leftNodes [ iF ];
      else if ( hasLeftNode )  innerNodes.front() = L._leftLine->_rightNodes[ iF ];
      if ( hasOwnRightNode )   innerNodes.back()  = L._rightNodes[ iF ];
      else if ( hasRightNode ) innerNodes.back()  = L._rightLine->_leftNodes[ iF ];
      if ( isClosedEdge )      innerNodes.front() = innerNodes.back(); // circle
      if ( !isShrinkableL ) innerNodes.front() = outerNodes.front();
      else if ( isCornerBevelL && nearChainL && (size_t)iF < nearChainL->size() )
        innerNodes.front() = (*nearChainL)[ iF ];
      if ( !isShrinkableR ) innerNodes.back()  = outerNodes.back();
      else if ( isCornerBevelR && nearChainR && (size_t)iF < nearChainR->size() )
        innerNodes.back() = (*nearChainR)[ iF ];
      if ( !hasOwnLeftNode )   L._leftNodes.push_back( innerNodes.front() );
      if ( !hasOwnRightNode )  L._rightNodes.push_back( innerNodes.back() );

      // create faces
      for ( size_t i = 1; i < innerNodes.size(); ++i )
      {
        const SMDS_MeshNode* n1 = outerNodes[ i+prev ];
        const SMDS_MeshNode* n2 = outerNodes[ i+cur  ];
        const SMDS_MeshNode* n3 = innerNodes[ i+cur  ];
        const SMDS_MeshNode* n4 = innerNodes[ i+prev ];
        // Reentrant-corner bevel (V1/V2 == isCornerBevelL/R): the
        // corner-adjacent element at that end reaches the near-chain node
        // (see innerNodes.front()/.back() override above) instead of
        // collapsing onto the plain VERTEX -- this element is this ARM's own
        // regular band, connecting its own near-chain to the far end of the
        // SAME EDGE (e.g. the requested "70-76" style connection), and is
        // safe from overlapping the underlying 2D algorithm's core fill only
        // because the near-chain is ALSO spliced into this FACE's own proxy
        // for its bare-adjacent EDGE (see shrink()'s buildCornerBevel() call
        // site). The bisector-ward bend itself is still handled entirely by
        // buildCornerBevel()'s own separate, small, explicitly-built wedge,
        // which shares the near-chain nodes with this element (not just the
        // plain VERTEX any more).
        if ( SMDS_MeshElement* f = _helper.AddFace( n1, n2, n3, n4 ))
          L._newFaces.insert( L._newFaces.end(), f );
      }

      outerNodes.swap( innerNodes );
    }

    // Add faces to a group
    SMDS_MeshGroup* group = StdMeshers_ViscousLayers::CreateGroup( hyp->GetGroupName(),
                                                                   *_helper.GetMesh(),
                                                                   SMDSAbs_Face );
    if ( group )
    {
      TIDSortedElemSet::iterator fIt = L._newFaces.begin();
      for ( ; fIt != L._newFaces.end(); ++fIt )
        group->Add( *fIt );
    }

    // faces between not shared _LayerEdge's (at concave VERTEX)
    for ( int isR = 0; isR < 2; ++isR )
    {
      if ( isR ? rightEdgeShared : leftEdgeShared )
        continue;
      vector< const SMDS_MeshNode* > &
        lNodes = (isR ? L._rightNodes : L._leftLine->_rightNodes ),
        rNodes = (isR ? L._rightLine->_leftNodes : L._leftNodes );
      if ( lNodes.empty() || rNodes.empty() || lNodes.size() != rNodes.size() )
        continue;

      const SMDS_MeshElement* face = 0;
      for ( size_t i = 1; i < lNodes.size(); ++i )
      {
        face = _helper.AddFace( lNodes[ i+prev ], rNodes[ i+prev ],
                                rNodes[ i+cur ],  lNodes[ i+cur ]);
        if ( group )
          group->Add( face );
      }

      const UVPtStruct& ptOnVertex = points[ isR ? L._lastPntInd : L._firstPntInd ];
      if ( isReverse )
        face = _helper.AddFace( ptOnVertex.node, lNodes[ 0 ], rNodes[ 0 ]);
      else
        face = _helper.AddFace( ptOnVertex.node, rNodes[ 0 ], lNodes[ 0 ]);
      if ( group )
        group->Add( face );
    }

    // Fill the _ProxyMeshOfFace

    UVPtStructVec nodeDataVec( outerNodes.size() ); // outerNodes swapped with innerNodes
    for ( size_t i = 0; i < outerNodes.size(); ++i )
    {
      gp_XY uv = _helper.GetNodeUV( _face, outerNodes[i] );
      nodeDataVec[i].u         = uv.X();
      nodeDataVec[i].v         = uv.Y();
      nodeDataVec[i].node      = outerNodes[i];
      nodeDataVec[i].param     = points [i + L._firstPntInd].param;
      nodeDataVec[i].normParam = normPar[i];
      nodeDataVec[i].x         = normPar[i];
      nodeDataVec[i].y         = normPar[i];
    }
    nodeDataVec.front().param = L._wire->FirstU( L._edgeInd );
    nodeDataVec.back() .param = L._wire->LastU ( L._edgeInd );

    if (( nodeDataVec[0].node == nodeDataVec.back().node ) &&
        ( _helper.GetPeriodicIndex() == 1 || _helper.GetPeriodicIndex() == 2 )) // closed EDGE
    {
      const int iCoord = _helper.GetPeriodicIndex();
      gp_XY uv = nodeDataVec[0].UV();
      uv.SetCoord( iCoord, L._lEdges[0]._uvOut.Coord( iCoord ));
      nodeDataVec[0].SetUV( uv );

      uv = nodeDataVec.back().UV();
      uv.SetCoord( iCoord, L._lEdges.back()._uvOut.Coord( iCoord ));
      nodeDataVec.back().SetUV( uv );
    }

    _ProxyMeshOfFace::_EdgeSubMesh* edgeSM
      = getProxyMesh()->GetEdgeSubMesh( L._wire->EdgeID( L._edgeInd ));
    edgeSM->SetUVPtStructVec( nodeDataVec );

  } // loop on _PolyLine's

  // re-compute FACEs whose mesh was removed by shrink()
  for ( size_t i = 0; i < _clearedFaces.size(); ++i )
  {
    SMESH_subMesh* sm = _mesh->GetSubMesh( _clearedFaces[i] );
    if ( sm->GetComputeState() == SMESH_subMesh::READY_TO_COMPUTE )
      sm->ComputeStateEngine( SMESH_subMesh::COMPUTE );
  }

  return true;
}

//================================================================================
/*!
 * \brief Improve quality of the created mesh elements
 */
//================================================================================

bool _ViscousBuilder2D::improve()
{
  if ( !_proxyMesh )
    return false;

  // fixed nodes on EDGE's
  std::set<const SMDS_MeshNode*> fixedNodes;
  for ( size_t iWire = 0; iWire < _faceSideVec.size(); ++iWire )
  {
    StdMeshers_FaceSidePtr      wire = _faceSideVec[ iWire ];
    const vector<UVPtStruct>& points = wire->GetUVPtStruct();
    for ( size_t i = 0; i < points.size(); ++i )
      fixedNodes.insert( fixedNodes.end(), points[i].node );
  }
  // fixed proxy nodes
  for ( size_t iL = 0; iL < _polyLineVec.size(); ++iL )
  {
    _PolyLine&         L = _polyLineVec[ iL ];
    const TopoDS_Edge& E = L._wire->Edge( L._edgeInd );
    if ( const SMESH_ProxyMesh::SubMesh* sm = _proxyMesh->GetProxySubMesh( E ))
    {
      const UVPtStructVec& points = sm->GetUVPtStructVec();
      for ( size_t i = 0; i < points.size(); ++i )
        fixedNodes.insert( fixedNodes.end(), points[i].node );
    }
    for ( size_t i = 0; i < L._rightNodes.size(); ++i )
      fixedNodes.insert( fixedNodes.end(), L._rightNodes[i] );
  }

  // smoothing
  SMESH_MeshEditor editor( _mesh );
  for ( size_t iL = 0; iL < _polyLineVec.size(); ++iL )
  {
    _PolyLine& L = _polyLineVec[ iL ];
    if ( L._isStraight2D ) continue;
    // SMESH_MeshEditor::SmoothMethod how =
    //   L._isStraight2D ? SMESH_MeshEditor::LAPLACIAN : SMESH_MeshEditor::CENTROIDAL;
    //editor.Smooth( L._newFaces, fixedNodes, how, /*nbIt = */3 );
    //editor.Smooth( L._newFaces, fixedNodes, SMESH_MeshEditor::LAPLACIAN, /*nbIt = */1 );
    editor.Smooth( L._newFaces, fixedNodes, SMESH_MeshEditor::CENTROIDAL, /*nbIt = */3 );
  }
  return true;
}

//================================================================================
/*!
 * \brief Remove elements and nodes from a face
 */
//================================================================================

bool _ViscousBuilder2D::removeMeshFaces(const TopoDS_Shape& face)
{
  // we don't use SMESH_subMesh::ComputeStateEngine() because of a listener
  // which clears EDGEs together with _face.
  bool thereWereElems = false;
  SMESH_subMesh* sm = _mesh->GetSubMesh( face );
  if ( SMESHDS_SubMesh* smDS = sm->GetSubMeshDS() )
  {
    SMDS_ElemIteratorPtr eIt = smDS->GetElements();
    thereWereElems = eIt->more();
    while ( eIt->more() ) getMeshDS()->RemoveFreeElement( eIt->next(), smDS );
    SMDS_NodeIteratorPtr nIt = smDS->GetNodes();
    while ( nIt->more() ) getMeshDS()->RemoveFreeNode( nIt->next(), smDS );
  }
  sm->ComputeStateEngine( SMESH_subMesh::CHECK_COMPUTE_STATE );

  return thereWereElems;
}

//================================================================================
/*!
 * \brief Returns a hypothesis for a _PolyLine
 */
//================================================================================

const StdMeshers_ViscousLayers2D* _ViscousBuilder2D::getLineHypothesis(int iPL)
{
  return iPL < (int)_hypOfEdge.size() ? _hypOfEdge[ iPL ] : _hyps[0];
}

//================================================================================
/*!
 * \brief Returns a layers thickness for a _PolyLine
 */
//================================================================================

double _ViscousBuilder2D::getLineThickness(int iPL)
{
  if ( const StdMeshers_ViscousLayers2D* h = getLineHypothesis( iPL ))
    return Min( _maxThickness, h->GetTotalThickness() );
  return _maxThickness;
}

//================================================================================
/*!
 * \brief Creates a _ProxyMeshOfFace and store it in a sub-mesh of FACE
 */
//================================================================================

_ProxyMeshOfFace* _ViscousBuilder2D::getProxyMesh()
{
  if ( _proxyMesh.get() )
    return (_ProxyMeshOfFace*) _proxyMesh.get();

  _ProxyMeshOfFace* proxyMeshOfFace = new _ProxyMeshOfFace( *_mesh );
  _proxyMesh.reset( proxyMeshOfFace );
  new _ProxyMeshHolder( _face, _proxyMesh );

  return proxyMeshOfFace;
}

//================================================================================
/*!
 * \brief Calculate height of layers for the given thickness. Height is measured
 *        from the outer boundary
 */
//================================================================================

void _ViscousBuilder2D::calcLayersHeight(const double    totalThick,
                                         vector<double>& heights,
                                         const THypVL*   hyp)
{
  const double fPowN = pow( hyp->GetStretchFactor(), hyp->GetNumberLayers() );
  heights.resize( hyp->GetNumberLayers() );
  double h0;
  if ( fPowN - 1 <= numeric_limits<double>::min() )
    h0 = totalThick / hyp->GetNumberLayers();
  else
    h0 = totalThick * ( hyp->GetStretchFactor() - 1 )/( fPowN - 1 );

  double hSum = 0, hi = h0;
  for ( int i = 0; i < hyp->GetNumberLayers(); ++i )
  {
    hSum += hi;
    heights[ i ] = hSum;
    hi *= hyp->GetStretchFactor();
  }
}

//================================================================================
/*!
 * \brief Elongate this _LayerEdge
 */
//================================================================================

bool _LayerEdge::SetNewLength( const double length3D )
{
  if ( _isBlocked ) return false;

  //_uvInPrev = _uvIn;
  _length2D = length3D * _len2dTo3dRatio;
  _uvIn     = _uvOut + _normal2D * _length2D;
  return true;
}

//================================================================================
/*!
 * \brief Return true if _LayerEdge at a common VERTEX between EDGEs with
 *  and w/o layer is common to the both _PolyLine's. If this is true, nodes
 *  of this _LayerEdge are inflated along a _PolyLine w/o layer, else the nodes
 *  are inflated along _normal2D of _LayerEdge of EDGE with layer
 */
//================================================================================

bool _PolyLine::IsCommonEdgeShared( const _PolyLine& other )
{
  const double tol = 1e-30;

  if ( & other == _leftLine )
    return _lEdges[0]._normal2D.IsEqual( _leftLine->_lEdges.back()._normal2D, tol );

  if ( & other == _rightLine )
    return _lEdges.back()._normal2D.IsEqual( _rightLine->_lEdges[0]._normal2D, tol );

  return false;
}

//================================================================================
/*!
 * \brief Return \c true if the EDGE of this _PolyLine is concave
 */
//================================================================================

bool _PolyLine::IsConcave() const
{
  if ( _lEdges.size() < 2 )
    return false;

  gp_Vec2d v1( _lEdges[0]._uvOut, _lEdges[1]._uvOut );
  gp_Vec2d v2( _lEdges[0]._uvOut, _lEdges[2]._uvOut );
  const double size2 = v2.Magnitude();

  return ( v1 ^ v2 ) / size2 < -1e-3 * size2;
}

//================================================================================
/*!
 * \brief Constructor of SegmentTree
 */
//================================================================================

_SegmentTree::_SegmentTree( const vector< _Segment >& segments ):
  SMESH_Quadtree()
{
  _segments.resize( segments.size() );
  for ( size_t i = 0; i < segments.size(); ++i )
    _segments[i].Set( segments[i] );

  compute();
}

//================================================================================
/*!
 * \brief Return the maximal bnd box
 */
//================================================================================

_SegmentTree::box_type* _SegmentTree::buildRootBox()
{
  _SegmentTree::box_type* box = new _SegmentTree::box_type;
  for ( size_t i = 0; i < _segments.size(); ++i )
  {
    box->Add( *_segments[i]._seg->_uv[0] );
    box->Add( *_segments[i]._seg->_uv[1] );
  }
  return box;
}

//================================================================================
/*!
 * \brief Redistribute _segments among children
 */
//================================================================================

void _SegmentTree::buildChildrenData()
{
  for ( size_t i = 0; i < _segments.size(); ++i )
    for (int j = 0; j < nbChildren(); j++)
      if ( !myChildren[j]->getBox()->IsOut( *_segments[i]._seg->_uv[0],
                                            *_segments[i]._seg->_uv[1] ))
        ((_SegmentTree*)myChildren[j])->_segments.push_back( _segments[i]);

  SMESHUtils::FreeVector( _segments ); // = _elements.clear() + free memory

  for (int j = 0; j < nbChildren(); j++)
  {
    _SegmentTree* child = static_cast<_SegmentTree*>( myChildren[j]);
    child->myIsLeaf = ((int) child->_segments.size() <= maxNbSegInLeaf() );
  }
}

//================================================================================
/*!
 * \brief Return elements which can include the point
 */
//================================================================================

void _SegmentTree::GetSegmentsNear( const _Segment&            seg,
                                    vector< const _Segment* >& found )
{
  if ( getBox()->IsOut( *seg._uv[0], *seg._uv[1] ))
    return;

  if ( isLeaf() )
  {
    for ( size_t i = 0; i < _segments.size(); ++i )
      if ( !_segments[i].IsOut( seg ))
        found.push_back( _segments[i]._seg );
  }
  else
  {
    for (int i = 0; i < nbChildren(); i++)
      ((_SegmentTree*) myChildren[i])->GetSegmentsNear( seg, found );
  }
}


//================================================================================
/*!
 * \brief Return segments intersecting a ray
 */
//================================================================================

void _SegmentTree::GetSegmentsNear( const gp_Ax2d&             ray,
                                    vector< const _Segment* >& found )
{
  if ( getBox()->IsOut( ray ))
    return;

  if ( isLeaf() )
  {
    for ( size_t i = 0; i < _segments.size(); ++i )
      if ( !_segments[i].IsOut( ray ))
        found.push_back( _segments[i]._seg );
  }
  else
  {
    for (int i = 0; i < nbChildren(); i++)
      ((_SegmentTree*) myChildren[i])->GetSegmentsNear( ray, found );
  }
}

//================================================================================
/*!
 * \brief Classify a _Segment
 */
//================================================================================

bool _SegmentTree::_SegBox::IsOut( const _Segment& seg ) const
{
  const double eps = std::numeric_limits<double>::min();
  for ( int iC = 0; iC < 2; ++iC )
  {
    if ( seg._uv[0]->Coord(iC+1) < _seg->_uv[ _iMin[iC]]->Coord(iC+1)+eps &&
         seg._uv[1]->Coord(iC+1) < _seg->_uv[ _iMin[iC]]->Coord(iC+1)+eps )
      return true;
    if ( seg._uv[0]->Coord(iC+1) > _seg->_uv[ 1-_iMin[iC]]->Coord(iC+1)-eps &&
         seg._uv[1]->Coord(iC+1) > _seg->_uv[ 1-_iMin[iC]]->Coord(iC+1)-eps )
      return true;
  }
  return false;
}

//================================================================================
/*!
 * \brief Classify a ray
 */
//================================================================================

bool _SegmentTree::_SegBox::IsOut( const gp_Ax2d& ray ) const
{
  double distBoxCenter2Ray =
    ray.Direction().XY() ^ ( ray.Location().XY() - 0.5 * (*_seg->_uv[0] + *_seg->_uv[1]));

  double boxSectionDiam =
    Abs( ray.Direction().X() ) * ( _seg->_uv[1-_iMin[1]]->Y() - _seg->_uv[_iMin[1]]->Y() ) +
    Abs( ray.Direction().Y() ) * ( _seg->_uv[1-_iMin[0]]->X() - _seg->_uv[_iMin[0]]->X() );

  return Abs( distBoxCenter2Ray ) > 0.5 * boxSectionDiam;
}

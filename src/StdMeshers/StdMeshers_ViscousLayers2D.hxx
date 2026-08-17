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

// File      : StdMeshers_ViscousLayers2D.hxx
// Created   : 23 Jul 2012
// Author    : Edward AGAPOV (eap)

#ifndef __StdMeshers_ViscousLayers2D_HXX__
#define __StdMeshers_ViscousLayers2D_HXX__

#include "StdMeshers_ViscousLayers.hxx"

#include <gp_Dir.hxx>
#include <gp_Pnt.hxx>

class TopoDS_Edge;
class TopoDS_Face;
class TopoDS_Vertex;
class StdMeshers_FaceSide;

/*!
 * \brief Hypothesis defining parameters of viscous layers
 */
class STDMESHERS_EXPORT StdMeshers_ViscousLayers2D : public StdMeshers_ViscousLayers
{
public:
  StdMeshers_ViscousLayers2D(int hypId, SMESH_Gen* gen);
  /*!
   * \brief Computes temporary 2D mesh to be used by 2D algorithm.
   *        Return SMESH_ProxyMesh for the given FACE, or NULL in case of error
   */
  static SMESH_ProxyMesh::Ptr Compute(SMESH_Mesh&        theMesh,
                                      const TopoDS_Face& theShape);
  /*!
   * \brief At study restoration, restore event listeners used to clear an inferior
   *  dim sub-mesh modified by viscous layers
   */
  void RestoreListeners() const;

  /*!
   * \brief Checks compatibility of assigned StdMeshers_ViscousLayers2D hypotheses
   */
  static SMESH_ComputeErrorPtr CheckHypothesis(SMESH_Mesh&         aMesh,
                                               const TopoDS_Shape& aShape,
                                               Hypothesis_Status&  aStatus);
  /*!
   * \brief Initialize my parameter values by the mesh built on the geometry
   * \param theMesh - the built mesh
   * \param theShape - the geometry of interest
   * \retval bool - true if parameter values have been successfully defined
   *
   * Just return false as this hypothesis does not have parameters values
   */
  virtual bool SetParametersByMesh(const SMESH_Mesh* theMesh, const TopoDS_Shape& theShape);

  /*!
   * \brief Initialize my parameter values by default parameters.
   *  \retval bool - true if parameter values have been successfully defined
   */
  virtual bool SetParametersByDefaults(const TDefaults& /*dflts*/, const SMESH_Mesh* /*theMesh*/=0)
  { return false; }

  static const char* GetHypType() { return "ViscousLayers2D"; }

  static void SetProxyMeshOfEdge( const StdMeshers_FaceSide& edgeNodes );
  static bool HasProxyMesh( const TopoDS_Face& face, SMESH_Mesh& theMesh );

  /*!
   * \brief Compute the geometric bisector B13 of two EDGEs E1, E3 sharing a
   *        VERTEX O (each carrying viscous layers on a different FACE), together
   *        with the points P1, P2 where each EDGE's own inward VL offset line
   *        (parallel to the EDGE, offset by its total layers thickness, extended
   *        past O) crosses B13. Pure geometry, no mesh dependency: used to bevel
   *        a reentrant corner shared by two independently viscous-layered FACEs
   *        through an intervening FACE with no viscous layers of its own.
   *  \return false if no well-defined bisector/intersection exists (degenerate
   *          angle between E1/E3, or the intersection lies behind O).
   */
  static bool ComputeCornerBisector( const TopoDS_Vertex& O,
                                     const TopoDS_Edge&   E1,
                                     const TopoDS_Face&   F1,
                                     double               thickness1,
                                     const TopoDS_Edge&   E3,
                                     const TopoDS_Face&   F3,
                                     double               thickness3,
                                     gp_Dir&              B13,
                                     gp_Pnt&              P1,
                                     gp_Pnt&              P2 );

 private:
};

#endif

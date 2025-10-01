#pragma once

enum eGUIGroupType{
    kBCMeshGroup = 0,
    kContactSlaveMeshGroup = 1,
    kContactMasterMeshGroup = 2,
    kMaterialMeshGroup = 3,
    kBCAmplitudeGroup = 5,
    kCESlaveMeshGroup = 6,
    kCEMasterMeshGroup = 7,
    kBodyCMeshGroup = 8,
    kBCMeshElemGroup = 9, // for shell pressure bc
    kBCMeshRigidBodyGroup = 10,
    kContactSelfGroup = 11
};

///////////////////////////////////////////////////////////////////////////////
/// Type (node, edge, face or volume) of elements
///////////////////////////////////////////////////////////////////////////////
enum WMDataElemType
  {
    WMDataAbs_All,
    WMDataAbs_Node,
    WMDataAbs_Edge,
    WMDataAbs_Face,
    WMDataAbs_Volume,
    WMDataAbs_0DElement,
    WMDataAbs_Ball,
    WMDataAbs_Vertex,
    WMDataAbs_ShapeBody,
    WMDataAbs_NbElementTypes
  };

/*! enumeration for element geometry type */
enum WMDataGeometryType
  {
    // 0D element
    WMDataGeom_POINT,
    // 1D element
    WMDataGeom_EDGE,
    // 2D element
    WMDataGeom_TRIANGLE,
    WMDataGeom_QUADRANGLE,
    WMDataGeom_POLYGON,
    // 3D element
    WMDataGeom_TETRA,
    WMDataGeom_PYRAMID,
    WMDataGeom_HEXA,
    WMDataGeom_PENTA,
    WMDataGeom_HEXAGONAL_PRISM,
    WMDataGeom_POLYHEDRA,
    // Discrete elements
    WMDataGeom_BALL,
    //
    WMDataGeom_NONE
  };


enum WMDataElementOrder {
  ORDER_ANY,          /*! entities of any order */
  ORDER_LINEAR,       /*! entities of 1st order */
  ORDER_QUADRATIC     /*! entities of 2nd order */
};

/*!
 * Enumeration of entity type used in mesh info array
 */
enum WMDataEntityType {
  WMDataEntity_Node,
  WMDataEntity_0D,
  WMDataEntity_Edge,
  WMDataEntity_Quad_Edge,
  WMDataEntity_Triangle,
  WMDataEntity_Quad_Triangle,
  WMDataEntity_BiQuad_Triangle,
  WMDataEntity_Quadrangle,
  WMDataEntity_Quad_Quadrangle,
  WMDataEntity_BiQuad_Quadrangle,
  WMDataEntity_Polygon,
  WMDataEntity_Quad_Polygon,
  WMDataEntity_Tetra,
  WMDataEntity_Quad_Tetra,
  WMDataEntity_Pyramid,
  WMDataEntity_Quad_Pyramid,
  WMDataEntity_Hexa,
  WMDataEntity_Quad_Hexa,
  WMDataEntity_TriQuad_Hexa,
  WMDataEntity_Penta,
  WMDataEntity_Quad_Penta,
  WMDataEntity_Hexagonal_Prism,
  WMDataEntity_Polyhedra,
  WMDataEntity_Quad_Polyhedra,
  WMDataEntity_Ball,
  WMDataEntity_Unknown,
  WMDataEntity_Last
};

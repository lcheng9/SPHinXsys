#pragma once


// Name:   Coordinate Systems
// -----------------------------------------------------------------------

// Record 1:        FORMAT (1I10)
//                  Field 1       -- Part UID

// Record 2:        FORMAT (40A2)
//                  Field 1       -- Part Name

// Record 3:        FORMAT (3I10)
//                  Field 1       -- Coordinate System Label
//                  Field 2       -- Coordinate System Type
//                                   = 0, Cartesian
//                                   = 1, Cylindrical
//                                   = 2, Spherical
//                  Field 3       -- Coordinate System Color

// Record 4:        FORMAT (40A2)
//                  Field 1       -- Coordinate System Name

// Record 5:        FORMAT (1P3D25.16)
//                  Field 1-3     -- Transformation Matrix Row 1

// Record 6:        FORMAT (1P3D25.16)
//                  Field 1-3     -- Transformation Matrix Row 2

// Record 7:        FORMAT (1P3D25.16)
//                  Field 1-3     -- Transformation Matrix Row 3

// Record 8:        FORMAT (1P3D25.16)
//                  Field 1-3     -- Transformation Matrix Row 4

// Records 3 thru 8 are repeated for each Coordinate System in the Part.

// Example:
//     -1
//   2420
//        100
// Untitled
//          6         1        15
// FEMAP Global Cylindrical (6)
//     1.0000000000000000E+0    0.0000000000000000E+0    0.0000000000000000E+0
//     0.0000000000000000E+0    1.0000000000000000E+0    0.0000000000000000E+0
//     0.0000000000000000E+0    0.0000000000000000E+0    1.0000000000000000E+0
//     0.0000000000000000E+0    0.0000000000000000E+0    0.0000000000000000E+0
//          7         2        15
// Coordinate System 4
//     1.0000000000000000E+0    0.0000000000000000E+0    0.0000000000000000E+0
//     0.0000000000000000E+0    1.0000000000000000E+0    0.0000000000000000E+0
//     0.0000000000000000E+0    0.0000000000000000E+0    1.0000000000000000E+0
//     0.0000000000000000E+0    0.0000000000000000E+0    0.0000000000000000E+0
//     -1

//#include "WMESH_TransferUNV.h"

#include <string>
#include <vector>

namespace unvMesh2420
{
  enum { Cartesian=0, Cylindrical, Spherical };

  typedef int TCSLabel; // type of coord system label

  struct UNVRecordData
  {
    TCSLabel    coord_sys_label; 
    int         coord_sys_type;  // { Cartesian=0, Cylindrical, Spherical }
    int         coord_sys_color;
    std::string coord_sys_name;
    double      matrix[4][3];

    bool        is_identity_tensor() const;
    void        modify_by_tensor      ( double* coords ) const;
    static void convert_by_CCS( double* coords );
    static void convert_by_SCS  ( double* coords );
  };
  
  typedef std::vector<UNVRecordData> TDataSet;

  void
  read_stream(std::ifstream& in_stream,
       std::string&   part_name, // can re-store a mesh name
       TDataSet&      theDataSet);

  void write_stream(std::ofstream& out_stream, const std::string& part_name); // can store a mesh name
  void write_stream(std::ostream&  out_stream, const std::string& part_name); // can store a mesh name

};


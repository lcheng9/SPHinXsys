#pragma once

//#include "../MeshDB_PreCompiled.h"

// Universal Dataset Number: 164
// Name:   Units
// Status: Current
// Owner:  General
// Revision Date: 19-AUG-1987
// -----------------------------------------------------------------------

// Record 1:       FORMAT(I10,20A1,I10)
//                 Field 1      -- units code
//                                 = 1 - SI: Meter (newton)
//                                 = 2 - BG: Foot (pound f)
//                                 = 3 - MG: Meter (kilogram f)
//                                 = 4 - BA: Foot (poundal)
//                                 = 5 - MM: mm (milli newton)
//                                 = 6 - CM: cm (centi newton)
//                                 = 7 - IN: Inch (pound f)
//                                 = 8 - GM: mm (kilogram f)
//                                 = 9 - US: USER_DEFINED
//                                 = 10- MN: mm (newton)
//                 Field 2      -- units description (used for
//                                 documentation only)
//                 Field 3      -- temperature mode
//                                 = 1 - absolute
//                                 = 2 - relative
// Record 2:       FORMAT(3D25.17)
//                 Unit factors for converting universal file units to SI.
//                 To convert from universal file units to SI divide by
//                 the appropriate factor listed below.
//                 Field 1      -- length
//                 Field 2      -- force
//                 Field 3      -- temperature
//                 Field 4      -- temperature offset

// Example:

//     -1
//    164
//          2Foot (pound f)               2
//   3.28083989501312334D+00  2.24808943099710480D-01  1.79999999999999999D+00
//   4.59670000000000002D+02
//     -1

//#include "WMESH_TransferUNV.h"

#include <string>

namespace unvMesh164
{
  enum { LENGTH_FACTOR, FORCE_FACTOR, TEMP_FACTOR, TEMP_OFFSET };

  struct UNVRecordData
  {
    int         units_code;
    std::string units_description;
    int         temp_mode;
    double      factors[4];
    UNVRecordData();
    void set_length_unit(const std::string& s);
  };
  
  void
  read_stream(std::ifstream& in_stream, UNVRecordData& theUnitsRecord);

  void write_stream(std::ofstream& out_stream, const UNVRecordData& theDataSet);
  void write_stream(std::ostream& out_stream, const UNVRecordData& theDataSet);
};


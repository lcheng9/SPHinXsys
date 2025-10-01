#pragma once

//#include "WMESH_TransferUNV.h"

#include <vector>
#include <fstream>      

//#include "../MeshDB_PreCompiled.h"

namespace unvMesh2411{

  struct UNVRecordData{
    UNVRecordData();
    unsigned long label;
    int exp_coord_sys_num;  // export coordinate system number
    int disp_coord_sys_num;  // displacement coordinate system number
    int color;  // color                                
    double coord[3];  // node coordinates in the part coordinate system
  };
  
  typedef std::vector<UNVRecordData> TDataSet;

  void read_stream(std::ifstream& in_stream, TDataSet& theDataSet);

  void write_stream(std::ofstream& out_stream, const TDataSet& theDataSet);
  void write_stream(std::ostream& out_stream, const TDataSet& theDataSet);
};


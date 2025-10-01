#pragma once

//#include "WMESH_TransferUNV.h"

#include "../Transfer/TransferWMDataMesh.h"
#include "../WMESHData/WMESHDataGroupBase.h"
#include <list>


typedef std::list<WMESHDataGroupBase*> TGroupList;

class TransferUNVWMDataWriteMesh: public TransferWMDataMesh
{
 public:
  ~TransferUNVWMDataWriteMesh();
  virtual Status Execute(bool binary);

  void add_mesh_bc_group(WMESHDataGroupBase* theGroup);
  void add_mesh_mat_group(WMESHDataGroupBase* theGroup);



private:
    Status executeBinary();
    Status executeAscii();

    void clearMatGroups();
    void clearBCGroups();



 private:
  TGroupList myMatGroups;
  TGroupList myBCGroups;
};


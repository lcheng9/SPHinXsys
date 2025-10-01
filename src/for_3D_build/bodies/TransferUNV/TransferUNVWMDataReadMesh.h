#pragma once

//#include "WMESH_TransferUNV.h"

#include "../Transfer/TransferWMDataMesh.h"
#include <map>
#include <string>


class WMDataMeshManager;
class WMDataMeshElemGroup;


typedef std::map<WMDataMeshElemGroup*, std::string> TGroupNamesMap;
typedef std::map<WMDataMeshElemGroup*, int> TGroupIdMap;

typedef std::map<WMDataMeshElemGroup*, std::string> TGroupNamesMap;
typedef std::map<WMDataMeshElemGroup*, int> TGroupIdMap;

class TransferUNVWMDataReadMesh: public TransferWMDataMesh
{
 public:
  TransferUNVWMDataReadMesh():TransferWMDataMesh(),myGroup(0) {};
  ~TransferUNVWMDataReadMesh();
 
  virtual Status Execute(bool binary);

  const WMDataMeshElemGroup* get_mesh_group() const { return myGroup;}
  const TGroupNamesMap& get_mesh_group_names() const { return myGroupNames; }
  const TGroupIdMap&    get_mesh_group_ids() const { return myGroupId; }

private:
    TransferFemMesh::Status executeAscii();


 private:
  WMDataMeshElemGroup* myGroup;
  TGroupNamesMap myGroupNames;
  TGroupIdMap    myGroupId;
};


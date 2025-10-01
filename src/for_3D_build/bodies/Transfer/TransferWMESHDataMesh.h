#ifndef _TRANSFER_WMESHData_MESH_H_
#define _TRANSFER_WMESHData_MESH_H_

#include "TransferFemMesh.h"
#include "../MeshDB_PreCompiled.h"

class WMESHDataMesh;

class MESHTRANSFER_EXPORT TransferWMESHDataMesh: public TransferFemMesh
{
 public:
  TransferWMESHDataMesh();
  void set_mesh(WMESHDataMesh *theMesh);
  
 protected:
  WMESHDataMesh *myMesh;
};

#endif

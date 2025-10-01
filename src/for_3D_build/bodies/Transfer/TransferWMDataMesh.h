#pragma once

#include "TransferFemMesh.h"
#include "../MeshDB_PreCompiled.h"

class WMESHDataMesh;

class MESHTRANSFER_EXPORT TransferWMDataMesh : public TransferFemMesh
{
public:
    TransferWMDataMesh();
    void set_mesh(WMESHDataMesh *theMesh);

protected:
    WMESHDataMesh *myMesh;

};


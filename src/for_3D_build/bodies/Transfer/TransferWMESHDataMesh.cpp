#include "TransferWMESHDataMesh.h"

using namespace std;

TransferWMESHDataMesh::TransferWMESHDataMesh():
  myMesh(NULL)
{}

void TransferWMESHDataMesh::set_mesh(WMESHDataMesh *theMesh)
{
  myMesh = theMesh;
}

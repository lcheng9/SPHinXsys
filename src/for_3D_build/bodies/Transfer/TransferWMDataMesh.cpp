#include "TransferWMDataMesh.h"

using namespace std;

TransferWMDataMesh::TransferWMDataMesh() :
    myMesh(NULL)
{

}

void TransferWMDataMesh::set_mesh(WMESHDataMesh *theMesh)
{
    myMesh = theMesh;
}

#include "TransferUNVTools.h"

using namespace std;
using namespace unvMesh;

int unvMesh::UNVCounter::myCounter = 0;

unvMesh::UNVCounter::UNVCounter()
{
    myCounter++;
}

unvMesh::UNVCounter::~UNVCounter()
{
    myCounter--;
}

string unvMesh::UNVCounter::GetUNVCountNumber()
{
    if (myCounter)
        return string(myCounter * 2, ' ');
    return "";
}

namespace unvMesh {

void writeFoamFile(const std::string& sClass, const std::string& sObject, std::ofstream& out_stream)
{
    //
    out_stream << "FoamFile" << std::endl;
    out_stream << "{" << std::endl;
    out_stream << "    version     2.0;" << std::endl;
    out_stream << "    format      ascii;" << std::endl;
    out_stream << "    class       " << sClass << ";" << std::endl;
    out_stream << "    object      " << sObject << ";" << std::endl;
    out_stream << "}" << std::endl;
}


int getElmerElemType(WMDataElemType elemType, int nNodes)
{
    int iVal = 0;
    switch (elemType)
    {
    case WMDataAbs_Edge:
    {
        if (nNodes == 2) {
            iVal = 202;
        }
        else if (nNodes == 3) {
            iVal = 203;
        }
        else if (nNodes == 4) {
            iVal = 204;
        }
        else {
            assert(0);
        }
    }
    break;
    case WMDataAbs_Face:
    {
        if (nNodes == 3) {
            iVal = 303;
        }
        else if (nNodes == 4) {
            iVal = 404;
        }
        else if (nNodes == 6) {
            iVal = 306;
        }
        else if (nNodes == 8) {
            iVal = 408;
        }
        else if (nNodes == 9) {
            iVal = 409;
        }
        else {
            assert(0);
        }
    }
    break;
    case WMDataAbs_Volume:
    {
        if (nNodes == 4) {
            iVal = 504;
        }
        else if (nNodes == 10) {
            iVal = 510;
        }
        else if (nNodes == 8) {
            iVal = 808;
        }
        else if (nNodes == 20) {
            iVal = 820;
        }
        else if (nNodes == 27) {
            iVal = 827;
        }
        else {
            assert(0);
        }
    }
    break;
    case WMDataAbs_Vertex:
        iVal = 101;
        break;
    default:
        assert(0);
        break;
    }

    return iVal;
}

}
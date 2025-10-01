#include "TransferUNVWMDataReadMesh.h"

#include "../WMData/WMDataMeshElemGroup.h"
#include "../WMESHData/WMESHDataMesh.h"

#include "../WELSIMLocalTrace/TraceTools.h"

#include "TransferUNV164DS.h"
#include "TransferUNV2411DS.h"
#include "TransferUNV2412DS.h"
#include "TransferUNV2417DS.h"
#include "TransferUNV2420DS.h"
#include "TransferUNVTools.h"

#include <gzstream.h> // For reading/writing compressed streams

#include "../Basics/Fundamental_Tools.h"

using namespace std;


#if defined(_DEBUG_) || defined(_DEBUG)
static int MYDEBUG = 1;
#else
static int MYDEBUG = 0;
#endif

namespace
{
/*!
 * \brief Move node coordinates to the global Cartesian CS
 */
void transformNodes(unvMesh2411::TDataSet::const_iterator fromNode,
    unvMesh2411::TDataSet::const_iterator endNode,
    const unvMesh2420::UNVRecordData &          csRecord)
{
    const int csLabel = fromNode->exp_coord_sys_num;

    unvMesh2411::TDataSet::const_iterator nodeIt;

    // apply Transformation Matrix
    if (!csRecord.is_identity_tensor())
    {
        for (nodeIt = fromNode; nodeIt != endNode; ++nodeIt)
        {
            const unvMesh2411::UNVRecordData& nodeRec = *nodeIt;
            if (nodeRec.exp_coord_sys_num == csLabel)
                csRecord.modify_by_tensor((double*)nodeRec.coord);
        }
    }

    // transform from Cylindrical CS
    if (csRecord.coord_sys_type == unvMesh2420::Cylindrical)
    {
        for (nodeIt = fromNode; nodeIt != endNode; ++nodeIt)
        {
            const unvMesh2411::UNVRecordData& nodeRec = *nodeIt;
            if (nodeRec.exp_coord_sys_num == csLabel)
                csRecord.convert_by_CCS((double*)nodeRec.coord);
        }
    }
    // transform from Spherical CS
    else if (csRecord.coord_sys_type == unvMesh2420::Spherical)
    {
        for (nodeIt = fromNode; nodeIt != endNode; ++nodeIt)
        {
            const unvMesh2411::UNVRecordData& nodeRec = *nodeIt;
            if (nodeRec.exp_coord_sys_num == csLabel)
                csRecord.convert_by_SCS((double*)nodeRec.coord);
        }
    }
}
}

TransferUNVWMDataReadMesh::~TransferUNVWMDataReadMesh()
{
    if (myGroup != 0)
        delete myGroup;
}

TransferFemMesh::Status TransferUNVWMDataReadMesh::Execute(bool binary)
{
    if (binary) {
        assert(0);
        return TransferFemMesh::DRS_EMPTY;
    }

    return executeAscii();
}


TransferFemMesh::Status TransferUNVWMDataReadMesh::executeAscii()
{
    Kernel_Utils::CLocalizer loc;
    Status aResult = DRS_OK;
    std::ifstream in_stream(myFile.c_str());
    try
    {
        {
            // Read Units
            unvMesh164::UNVRecordData aUnitsRecord;
            unvMesh164::read_stream(in_stream, aUnitsRecord);

            // Read Coordinate systems
            unvMesh2420::TDataSet aCoordSysDataSet;
            unvMesh2420::read_stream(in_stream, myMeshName, aCoordSysDataSet);

            // Read nodes
            using namespace unvMesh2411;
            TDataSet aDataSet2411;
            unvMesh2411::read_stream(in_stream, aDataSet2411);
            //if(MYDEBUG) MESSAGE("Perform - aDataSet2411.size() = "<<aDataSet2411.size());

            // Move nodes in a global CS
            if (!aCoordSysDataSet.empty())
            {
                unvMesh2420::TDataSet::const_iterator csIter = aCoordSysDataSet.begin();
                for (; csIter != aCoordSysDataSet.end(); ++csIter)
                {
                    // find any node in this CS
                    TDataSet::const_iterator nodeIter = aDataSet2411.begin();
                    for (; nodeIter != aDataSet2411.end(); nodeIter++)
                        if (nodeIter->exp_coord_sys_num == csIter->coord_sys_label)
                        {
                            transformNodes(nodeIter, aDataSet2411.end(), *csIter);
                            break;
                        }
                }
            }
            // Move nodes to SI unit system
            const double lenFactor = aUnitsRecord.factors[unvMesh164::LENGTH_FACTOR];
            if (lenFactor != 1.)
            {
                TDataSet::iterator nodeIter = aDataSet2411.begin(), nodeEnd;
                for (nodeEnd = aDataSet2411.end(); nodeIter != nodeEnd; nodeIter++)
                {
                    unvMesh2411::UNVRecordData& nodeRec = *nodeIter;
                    nodeRec.coord[0] *= lenFactor;
                    nodeRec.coord[1] *= lenFactor;
                    nodeRec.coord[2] *= lenFactor;
                }
            }

            // Create nodes in the mesh
            TDataSet::const_iterator anIter = aDataSet2411.begin();
            for (; anIter != aDataSet2411.end(); anIter++)
            {
                const UNVRecordData& aRec = *anIter;
                myMesh->add_node_and_id(aRec.coord[0], aRec.coord[1], aRec.coord[2], aRec.label);
            }
        }
        {
            using namespace unvMesh2412;
            TDataSet aDataSet2412;
            unvMesh2412::read_stream(in_stream, aDataSet2412);
            TDataSet::const_iterator anIter = aDataSet2412.begin();
            //if(MYDEBUG) MESSAGE("Perform - aDataSet2412.size() = "<<aDataSet2412.size());
            for (; anIter != aDataSet2412.end(); anIter++)
            {
                WMDataMeshElem* anElement = NULL;
                const UNVRecordData& aRec = *anIter;
                if (is_beam_element(aRec.fe_descriptor_id)) {
                    switch (aRec.node_labels.size()) {
                    case 2: // edge with two nodes
                        MESSAGE("add edge " << aRec.label << " " << aRec.node_labels[0] << " " << aRec.node_labels[1]);
                        anElement = myMesh->add_edge_and_id(aRec.node_labels[0],
                            aRec.node_labels[1],
                            aRec.label);
                        break;
                    case 3: // quadratic edge (with 3 nodes)
                        MESSAGE("add edge " << aRec.label << " " << aRec.node_labels[0] << " " << aRec.node_labels[1] << " " << aRec.node_labels[2]);
                        anElement = myMesh->add_edge_and_id(aRec.node_labels[0],
                            aRec.node_labels[2],
                            aRec.node_labels[1],
                            aRec.label);
                    }
                }
                else if (is_face_element(aRec.fe_descriptor_id)) {
                    MESSAGE("add face " << aRec.label);
                    switch (aRec.fe_descriptor_id) {
                    case 41: // Plane Stress Linear Triangle
                    case 51: // Plane Strain Linear Triangle
                    case 61: // Plate Linear Triangle
                    case 74: // Membrane Linear Triangle
                    case 81: // Axisymetric Solid Linear Triangle
                    case 91: // Thin Shell Linear Triangle
                        anElement = myMesh->add_face_and_id(aRec.node_labels[0],
                            aRec.node_labels[1],
                            aRec.node_labels[2],
                            aRec.label);
                        break;

                    case 42: //  Plane Stress Parabolic Triangle
                    case 52: //  Plane Strain Parabolic Triangle
                    case 62: //  Plate Parabolic Triangle
                    case 72: //  Membrane Parabolic Triangle
                    case 82: //  Axisymetric Solid Parabolic Triangle
                    case 92: //  Thin Shell Parabolic Triangle
                        anElement = myMesh->add_face_and_id(aRec.node_labels[0],
                            aRec.node_labels[2],
                            aRec.node_labels[4],
                            aRec.node_labels[1],
                            aRec.node_labels[3],
                            aRec.node_labels[5],
                            aRec.label);
                        break;

                    case 44: // Plane Stress Linear Quadrilateral
                    case 54: // Plane Strain Linear Quadrilateral
                    case 64: // Plate Linear Quadrilateral
                    case 71: // Membrane Linear Quadrilateral
                    case 84: // Axisymetric Solid Linear Quadrilateral
                    case 94: // Thin Shell Linear Quadrilateral
                        anElement = myMesh->add_face_and_id(aRec.node_labels[0],
                            aRec.node_labels[1],
                            aRec.node_labels[2],
                            aRec.node_labels[3],
                            aRec.label);
                        break;

                    case 45: // Plane Stress Parabolic Quadrilateral
                    case 55: // Plane Strain Parabolic Quadrilateral
                    case 65: // Plate Parabolic Quadrilateral
                    case 75: // Membrane Parabolic Quadrilateral
                    case 85: // Axisymetric Solid Parabolic Quadrilateral
                    case 95: // Thin Shell Parabolic Quadrilateral
                        if (aRec.node_labels.size() == 9) {
                            assert(0);
                            //anElement = myMesh->add_face_and_id(aRec.node_labels[0],
                            //    aRec.node_labels[2],
                            //    aRec.node_labels[4],
                            //    aRec.node_labels[6],
                            //    aRec.node_labels[1],
                            //    aRec.node_labels[3],
                            //    aRec.node_labels[5],
                            //    aRec.node_labels[7],
                            //    aRec.node_labels[8],
                            //    aRec.label);
                        }
                        else
                            anElement = myMesh->add_face_and_id(aRec.node_labels[0],
                                aRec.node_labels[2],
                                aRec.node_labels[4],
                                aRec.node_labels[6],
                                aRec.node_labels[1],
                                aRec.node_labels[3],
                                aRec.node_labels[5],
                                aRec.node_labels[7],
                                aRec.label);
                        break;
                    }
                }
                else if (is_body_element(aRec.fe_descriptor_id)) {
                    MESSAGE("add volume " << aRec.label);
                    switch (aRec.fe_descriptor_id) {

                    case 111: // Solid Linear Tetrahedron - TET4
                        anElement = myMesh->add_volume_and_id(aRec.node_labels[0],
                            aRec.node_labels[1],
                            aRec.node_labels[2],
                            aRec.node_labels[3],
                            aRec.label);
                        break;

                    case 118: // Solid Quadratic Tetrahedron - TET10
                        anElement = myMesh->add_volume_and_id(aRec.node_labels[0],
                            aRec.node_labels[2],
                            aRec.node_labels[4],

                            aRec.node_labels[9],

                            aRec.node_labels[1],
                            aRec.node_labels[3],
                            aRec.node_labels[5],

                            aRec.node_labels[6],
                            aRec.node_labels[7],
                            aRec.node_labels[8],
                            aRec.label);
                        break;

                    case 112: // Solid Linear Prism - PRISM6
                        anElement = myMesh->add_volume_and_id(aRec.node_labels[0],
                            aRec.node_labels[2],
                            aRec.node_labels[1],
                            aRec.node_labels[3],
                            aRec.node_labels[5],
                            aRec.node_labels[4],
                            aRec.label);
                        break;

                    case 113: // Solid Quadratic Prism - PRISM15
                        anElement = myMesh->add_volume_and_id(aRec.node_labels[0],
                            aRec.node_labels[4],
                            aRec.node_labels[2],

                            aRec.node_labels[9],
                            aRec.node_labels[13],
                            aRec.node_labels[11],

                            aRec.node_labels[5],
                            aRec.node_labels[3],
                            aRec.node_labels[1],

                            aRec.node_labels[14],
                            aRec.node_labels[12],
                            aRec.node_labels[10],

                            aRec.node_labels[6],
                            aRec.node_labels[8],
                            aRec.node_labels[7],
                            aRec.label);
                        break;

                    case 115: // Solid Linear Brick - HEX8
                        anElement = myMesh->add_volume_and_id(aRec.node_labels[0],
                            aRec.node_labels[3],
                            aRec.node_labels[2],
                            aRec.node_labels[1],
                            aRec.node_labels[4],
                            aRec.node_labels[7],
                            aRec.node_labels[6],
                            aRec.node_labels[5],
                            aRec.label);
                        break;

                    case 116: // Solid Quadratic Brick - HEX20
                        anElement = myMesh->add_volume_and_id(aRec.node_labels[0],
                            aRec.node_labels[6],
                            aRec.node_labels[4],
                            aRec.node_labels[2],

                            aRec.node_labels[12],
                            aRec.node_labels[18],
                            aRec.node_labels[16],
                            aRec.node_labels[14],

                            aRec.node_labels[7],
                            aRec.node_labels[5],
                            aRec.node_labels[3],
                            aRec.node_labels[1],

                            aRec.node_labels[19],
                            aRec.node_labels[17],
                            aRec.node_labels[15],
                            aRec.node_labels[13],

                            aRec.node_labels[8],
                            aRec.node_labels[11],
                            aRec.node_labels[10],
                            aRec.node_labels[9],
                            aRec.label);
                        break;

                    case 114: // pyramid of 13 nodes (quadratic) - PIRA13
                        anElement = myMesh->add_volume_and_id(aRec.node_labels[0],
                            aRec.node_labels[6],
                            aRec.node_labels[4],
                            aRec.node_labels[2],
                            aRec.node_labels[7],
                            aRec.node_labels[5],
                            aRec.node_labels[3],
                            aRec.node_labels[1],

                            aRec.node_labels[8],
                            aRec.node_labels[11],
                            aRec.node_labels[10],
                            aRec.node_labels[9],
                            aRec.node_labels[12],
                            aRec.label);
                        break;

                    }
                }
                if (!anElement) {
                    MESSAGE("TransferUNVWMDataReadMesh::Perform - can not add element with ID = " << aRec.label << " and type = " << aRec.fe_descriptor_id);
                }
            }
        }


        {
            using namespace unvMesh2417;
            TDataSet aDataSet2417;
            unvMesh2417::read_stream(in_stream, aDataSet2417);
            if (MYDEBUG) MESSAGE("Perform - aDataSet2417.size() = " << aDataSet2417.size());
            if (aDataSet2417.size() > 0) {
                myGroup = new WMDataMeshElemGroup(myMesh);
                TDataSet::const_iterator anIter = aDataSet2417.begin();
                for (; anIter != aDataSet2417.end(); anIter++) {
                    const TGroupId& aLabel = anIter->first;
                    const UNVRecordData& aRec = anIter->second;

                    int aNodesNb = aRec.NodeList.size();
                    int aElementsNb = aRec.ElementList.size();

                    bool useSuffix = ((aNodesNb > 0) && (aElementsNb > 0));
                    int i;
                    if (aNodesNb > 0) {
                        WMDataMeshElemGroup* aNodesGroup = (WMDataMeshElemGroup*)myGroup->add_subgroup(WMDataAbs_Node);
                        std::string aGrName = (useSuffix) ? aRec.GroupName + "_Nodes" : aRec.GroupName;
                        int i = aGrName.find("\r");
                        if (i > 0)
                            aGrName.erase(i, 2);
                        myGroupNames.insert(TGroupNamesMap::value_type(aNodesGroup, aGrName));
                        myGroupId.insert(TGroupIdMap::value_type(aNodesGroup, aLabel));

                        for (i = 0; i < aNodesNb; i++) {
                            const WMDataMeshNode* aNode = myMesh->look_for_node(aRec.NodeList[i]);
                            if (aNode)
                                aNodesGroup->add_elem(aNode);
                        }
                    }
                    if (aElementsNb > 0) {
                        WMDataMeshElemGroup* aEdgesGroup = 0;
                        WMDataMeshElemGroup* aFacesGroup = 0;
                        WMDataMeshElemGroup* aVolumeGroup = 0;
                        bool createdGroup = false;

                        for (i = 0; i < aElementsNb; i++) {
                            const WMDataMeshElem* aElement = myMesh->look_for_element(aRec.ElementList[i]);
                            if (aElement) {
                                switch (aElement->get_element_type()) {
                                case WMDataAbs_Edge:
                                    if (!aEdgesGroup) {
                                        aEdgesGroup = (WMDataMeshElemGroup*)myGroup->add_subgroup(WMDataAbs_Edge);
                                        if (!useSuffix && createdGroup) useSuffix = true;
                                        std::string aEdgesGrName = (useSuffix) ? aRec.GroupName + "_Edges" : aRec.GroupName;
                                        int i = aEdgesGrName.find("\r");
                                        if (i > 0)
                                            aEdgesGrName.erase(i, 2);
                                        myGroupNames.insert(TGroupNamesMap::value_type(aEdgesGroup, aEdgesGrName));
                                        myGroupId.insert(TGroupIdMap::value_type(aEdgesGroup, aLabel));
                                        createdGroup = true;
                                    }
                                    aEdgesGroup->add_elem(aElement);
                                    break;
                                case WMDataAbs_Face:
                                    if (!aFacesGroup) {
                                        aFacesGroup = (WMDataMeshElemGroup*)myGroup->add_subgroup(WMDataAbs_Face);
                                        if (!useSuffix && createdGroup) useSuffix = true;
                                        std::string aFacesGrName = (useSuffix) ? aRec.GroupName + "_Faces" : aRec.GroupName;
                                        int i = aFacesGrName.find("\r");
                                        if (i > 0)
                                            aFacesGrName.erase(i, 2);
                                        myGroupNames.insert(TGroupNamesMap::value_type(aFacesGroup, aFacesGrName));
                                        myGroupId.insert(TGroupIdMap::value_type(aFacesGroup, aLabel));
                                        createdGroup = true;
                                    }
                                    aFacesGroup->add_elem(aElement);
                                    break;
                                case WMDataAbs_Volume:
                                    if (!aVolumeGroup) {
                                        aVolumeGroup = (WMDataMeshElemGroup*)myGroup->add_subgroup(WMDataAbs_Volume);
                                        if (!useSuffix && createdGroup) useSuffix = true;
                                        std::string aVolumeGrName = (useSuffix) ? aRec.GroupName + "_Volumes" : aRec.GroupName;
                                        int i = aVolumeGrName.find("\r");
                                        if (i > 0)
                                            aVolumeGrName.erase(i, 2);
                                        myGroupNames.insert(TGroupNamesMap::value_type(aVolumeGroup, aVolumeGrName));
                                        myGroupId.insert(TGroupIdMap::value_type(aVolumeGroup, aLabel));
                                        createdGroup = true;
                                    }
                                    aVolumeGroup->add_elem(aElement);
                                    break;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    catch (const std::exception& exc) {
        INFOS("Exception:\n\t" << exc.what());
    }
    catch (...) {
        INFOS("Unknown exception");
    }

    return aResult;
}

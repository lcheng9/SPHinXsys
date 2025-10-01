#include <algorithm>
#include <iostream>
#include <string>


#include "TransferUNVWMDataWriteMesh.h"

#include "../WMData/WMDataEdgeQuadratic.h"
#include "../WMData/WMDataFaceQuadByNodes.h"
#include "../WMData/WMDataVolumePolyByNodes.h"
#include "../WMESHData/WMESHDataGroupBase.h"
#include "../WMESHData/WMESHDataMesh.h"

#include "../WELSIMLocalTrace/TraceTools.h"

#include "TransferUNV164DS.h"
#include "TransferUNV2411DS.h"
#include "TransferUNV2412DS.h"
#include "TransferUNV2417DS.h"
#include "TransferUNV2420DS.h"
#include "TransferUNVTools.h"

#include "../Basics/Fundamental_Tools.h"

#include <gzstream.h> // For reading/writing compressed streams

using namespace std;
using namespace unvMesh;

namespace {
typedef std::vector<size_t> TConnect;

int GetConnect(const WMDataElemIteratorPtr& theNodesIter,
    TConnect& theConnect)
{
    theConnect.clear();
    for (; theNodesIter->more();) {
        const WMDataMeshElem* anElem = theNodesIter->next();
        theConnect.push_back(anElem->get_elem_ID());
    }
    return theConnect.size();
}

}


TransferUNVWMDataWriteMesh::~TransferUNVWMDataWriteMesh()
{
    clearMatGroups();
    clearBCGroups();
}

void TransferUNVWMDataWriteMesh::add_mesh_bc_group(WMESHDataGroupBase* theGroup)
{
    myBCGroups.push_back(theGroup);
}

void TransferUNVWMDataWriteMesh::add_mesh_mat_group(WMESHDataGroupBase* theGroup)
{
    myMatGroups.push_back(theGroup);
}


TransferFemMesh::Status TransferUNVWMDataWriteMesh::Execute(bool binary)
{
    if (binary) {
        return executeBinary();
    }
    else {
        return executeAscii();
    }
}

TransferFemMesh::Status TransferUNVWMDataWriteMesh::executeBinary()
{
    //#ifdef WELSIM_MESH_BZIP2
    //  std::ofstream out_stream(myFile.c_str(), std::ios::binary);
    //  boost::iostreams::filtering_ostreambuf zdat;
    //  zdat.push(boost::iostreams::bzip2_compressor());  // your compressor here
    //  zdat.push(boost::iostreams::file_sink(myFile.c_str()));
    //#endif

    Kernel_Utils::CLocalizer loc;
    Status aResult = DRS_OK;

    ogzstream out_stream(myFile.c_str());

    try {
        // Units
        {
            using namespace unvMesh164;
            UNVRecordData aDataSet164;
            //aDataSet164.set_length_unit(myLengthUnit);
            unvMesh164::write_stream(out_stream, aDataSet164);
        }

        // Coordinate system
        unvMesh2420::write_stream(out_stream, myMeshName);

        // Nodes
        {
            using namespace unvMesh2411;
            TDataSet aDataSet2411;
            // Storing WMDS nodes to the UNV file
            //-----------------------------------
            MESSAGE("Perform - myMesh->get_number_of_nodes() = " << myMesh->get_number_of_nodes());
            WMDataNodeIteratorPtr aNodesIter = myMesh->get_nodes_iterator();
            UNVRecordData aRec;
            while (aNodesIter->more())
            {
                const WMDataMeshNode* aNode = aNodesIter->next();
                aRec.label = aNode->get_elem_ID();
                aRec.coord[0] = aNode->X();
                aRec.coord[1] = aNode->Y();
                aRec.coord[2] = aNode->Z();
                aDataSet2411.push_back(aRec);
            }
            MESSAGE("Perform - aDataSet2411.size() = " << aDataSet2411.size());
            unvMesh2411::write_stream(out_stream, aDataSet2411);
        }

        // Elements
        {
            using namespace unvMesh2412;
            TDataSet aDataSet2412;

            // Storing WMDS Edges
            MESSAGE("Perform - myMesh->get_number_of_edges() = " << myMesh->get_number_of_edges());
            if (myMesh->get_number_of_edges()) {
                WMDataEdgeIteratorPtr anIter = myMesh->get_edges_iterator();
                while (anIter->more())
                {
                    const WMDataEdge* anElem = anIter->next();
                    int aNbNodes = anElem->get_number_of_nodes();
                    UNVRecordData aRec;
                    aRec.label = anElem->get_elem_ID();
                    aRec.mat_prop_tab_num = anElem->getComponentId();
                    aRec.node_labels.reserve(aNbNodes);
                    WMDataNodeIteratorPtr aNodesIter = anElem->get_nodes_iter();
                    if (anElem->is_quadratic_element()) {
                        aRec.fe_descriptor_id = 22;
                        const WMDataMeshNode* aNode0 = aNodesIter->next();
                        const WMDataMeshNode* aNode1 = aNodesIter->next();
                        const WMDataMeshNode* aNode2 = aNodesIter->next();
                        aRec.node_labels.push_back(aNode0->get_elem_ID());
                        aRec.node_labels.push_back(aNode2->get_elem_ID());
                        aRec.node_labels.push_back(aNode1->get_elem_ID());
                    }
                    else {
                        aRec.fe_descriptor_id = 11;
                        const WMDataMeshNode* aNode0 = aNodesIter->next();
                        const WMDataMeshNode* aNode1 = aNodesIter->next();
                        aRec.node_labels.push_back(aNode0->get_elem_ID());
                        aRec.node_labels.push_back(aNode1->get_elem_ID());
                    }
                    if (aRec.fe_descriptor_id > 0) {
                        aDataSet2412.push_back(aRec);
                    }
                }
                MESSAGE("Perform - aDataSet2412.size() = " << aDataSet2412.size());
            }

            MESSAGE("Perform - myMesh->get_number_of_faces() = " << myMesh->get_number_of_faces());
            if (myMesh->get_number_of_faces())
            {
                WMDataFaceIteratorPtr anIter = myMesh->get_faces_iterator();
                while (anIter->more())
                {
                    const WMDataFace* anElem = anIter->next();
                    if (anElem->is_polyhedron()) continue;
                    int aNbNodes = anElem->get_number_of_nodes();
                    UNVRecordData aRec;
                    aRec.label = anElem->get_elem_ID();
                    aRec.mat_prop_tab_num = anElem->getComponentId();
                    aRec.node_labels.reserve(aNbNodes);
                    WMDataNodeIteratorPtr aNodesIter = anElem->get_nodes_iter();
                    switch (aNbNodes) {
                    case 3:
                    {
                        aRec.fe_descriptor_id = 41;
                        const WMDataMeshNode* aNode0 = aNodesIter->next();
                        const WMDataMeshNode* aNode1 = aNodesIter->next();
                        const WMDataMeshNode* aNode2 = aNodesIter->next();
                        aRec.node_labels.push_back(aNode0->get_elem_ID());
                        aRec.node_labels.push_back(aNode1->get_elem_ID());
                        aRec.node_labels.push_back(aNode2->get_elem_ID());
                    }
                    break;
                    case 4:
                    {
                        aRec.fe_descriptor_id = 44;
                        assert(0);
                        while (aNodesIter->more())
                        {
                            const WMDataMeshNode* aNode = aNodesIter->next();
                            aRec.node_labels.push_back(aNode->get_elem_ID());
                        }
                    }
                    break;
                    case 6:
                    {
                        aRec.fe_descriptor_id = 42;
                        const WMDataMeshNode* aNode0 = aNodesIter->next();
                        const WMDataMeshNode* aNode1 = aNodesIter->next();
                        const WMDataMeshNode* aNode2 = aNodesIter->next();
                        const WMDataMeshNode* aNode3 = aNodesIter->next();
                        const WMDataMeshNode* aNode4 = aNodesIter->next();
                        const WMDataMeshNode* aNode5 = aNodesIter->next();
                        aRec.node_labels.push_back(aNode0->get_elem_ID());
                        aRec.node_labels.push_back(aNode3->get_elem_ID());
                        aRec.node_labels.push_back(aNode1->get_elem_ID());
                        aRec.node_labels.push_back(aNode4->get_elem_ID());
                        aRec.node_labels.push_back(aNode2->get_elem_ID());
                        aRec.node_labels.push_back(aNode5->get_elem_ID());
                    }
                    break;
                    case 7:
                    {
                        aRec.fe_descriptor_id = 42;
                        assert(0);
                        while (aNodesIter->more())
                        {
                            const WMDataMeshNode* aNode = aNodesIter->next();
                            aRec.node_labels.push_back(aNode->get_elem_ID());
                        }
                    }
                    break;
                    case 8:
                    {
                        aRec.fe_descriptor_id = 45;
                        assert(0);
                        while (aNodesIter->more())
                        {
                            const WMDataMeshNode* aNode = aNodesIter->next();
                            aRec.node_labels.push_back(aNode->get_elem_ID());
                        }
                    }
                    break;
                    case 9:
                    {
                        aRec.fe_descriptor_id = 45;
                        aRec.node_labels.resize(8);
                        while (aNodesIter->more())
                        {
                            const WMDataMeshNode* aNode = aNodesIter->next();
                            aRec.node_labels.push_back(aNode->get_elem_ID());
                        }
                        assert(0);
                    }
                    break;
                    default:
                        assert(0);
                        continue;
                    }
                    if (aRec.fe_descriptor_id > 0) {
                        aDataSet2412.push_back(aRec);
                    }
                }
                MESSAGE("Perform - aDataSet2412.size() = " << aDataSet2412.size());
            }

            MESSAGE("Perform - myMesh->get_volumes_count() = " << myMesh->get_volumes_count());
            if (myMesh->get_volumes_count() > 0)
            {
                WMDataVolumeIteratorPtr anIter = myMesh->get_volumes_iterator();
                while (anIter->more())
                {
                    const WMDataVolume* anElem = anIter->next();
                    if (anElem->is_polyhedron())
                        continue;
                    UNVRecordData aRec;
                    aRec.label = anElem->get_elem_ID();
                    aRec.mat_prop_tab_num = anElem->getComponentId();
                    int aNbNodes = anElem->get_number_of_nodes();
                    aRec.node_labels.reserve(aNbNodes);
                    WMDataNodeIteratorPtr aNodesIter = anElem->get_nodes_iter();
                    switch (aNbNodes) {
                    case 4:
                    {
                        aRec.fe_descriptor_id = 111;
                        const WMDataMeshNode* aNode0 = aNodesIter->next();
                        const WMDataMeshNode* aNode1 = aNodesIter->next();
                        const WMDataMeshNode* aNode2 = aNodesIter->next();
                        const WMDataMeshNode* aNode3 = aNodesIter->next();
                        aRec.node_labels.push_back(aNode0->get_elem_ID());
                        aRec.node_labels.push_back(aNode1->get_elem_ID());
                        aRec.node_labels.push_back(aNode2->get_elem_ID());
                        aRec.node_labels.push_back(aNode3->get_elem_ID());
                    }
                    break;
                    case 6:
                    {
                        aRec.fe_descriptor_id = 112;
                        assert(0);
                        while (aNodesIter->more() && aRec.node_labels.size() < aNbNodes)
                        {
                            const WMDataMeshElem* aNode = aNodesIter->next();
                            aRec.node_labels.push_back(aNode->get_elem_ID());
                        }
                    }
                    break;
                    case 8:
                    {
                        aRec.fe_descriptor_id = 115;
                        assert(0);
                        while (aNodesIter->more() && aRec.node_labels.size() < aNbNodes)
                        {
                            const WMDataMeshElem* aNode = aNodesIter->next();
                            aRec.node_labels.push_back(aNode->get_elem_ID());
                        }
                    }
                    break;
                    case 10:
                    {
                        aRec.fe_descriptor_id = 118;
                        const WMDataMeshNode* aNode0 = aNodesIter->next();
                        const WMDataMeshNode* aNode1 = aNodesIter->next();
                        const WMDataMeshNode* aNode2 = aNodesIter->next();
                        const WMDataMeshNode* aNode3 = aNodesIter->next();
                        const WMDataMeshNode* aNode4 = aNodesIter->next();
                        const WMDataMeshNode* aNode5 = aNodesIter->next();
                        const WMDataMeshNode* aNode6 = aNodesIter->next();
                        const WMDataMeshNode* aNode7 = aNodesIter->next();
                        const WMDataMeshNode* aNode8 = aNodesIter->next();
                        const WMDataMeshNode* aNode9 = aNodesIter->next();
                        aRec.node_labels.push_back(aNode0->get_elem_ID());
                        aRec.node_labels.push_back(aNode4->get_elem_ID());
                        aRec.node_labels.push_back(aNode1->get_elem_ID());
                        aRec.node_labels.push_back(aNode5->get_elem_ID());
                        aRec.node_labels.push_back(aNode2->get_elem_ID());
                        aRec.node_labels.push_back(aNode6->get_elem_ID());

                        aRec.node_labels.push_back(aNode7->get_elem_ID());
                        aRec.node_labels.push_back(aNode8->get_elem_ID());
                        aRec.node_labels.push_back(aNode9->get_elem_ID());
                        aRec.node_labels.push_back(aNode3->get_elem_ID());
                    }
                    break;
                    case 13:
                    {
                        aRec.fe_descriptor_id = 114;
                        assert(0);
                        while (aNodesIter->more() && aRec.node_labels.size() < aNbNodes)
                        {
                            const WMDataMeshElem* aNode = aNodesIter->next();
                            aRec.node_labels.push_back(aNode->get_elem_ID());
                        }
                    }
                    break;
                    case 15:
                    {
                        aRec.fe_descriptor_id = 113;
                        assert(0);
                        while (aNodesIter->more() && aRec.node_labels.size() < aNbNodes)
                        {
                            const WMDataMeshElem* aNode = aNodesIter->next();
                            aRec.node_labels.push_back(aNode->get_elem_ID());
                        }
                    }
                    break;
                    case 20:
                    case 27:
                    {
                        aRec.fe_descriptor_id = 116;
                        aNbNodes = 20;
                        assert(0);
                        while (aNodesIter->more() && aRec.node_labels.size() < aNbNodes)
                        {
                            const WMDataMeshElem* aNode = aNodesIter->next();
                            aRec.node_labels.push_back(aNode->get_elem_ID());
                        }
                    }
                    break;
                    default:
                        assert(0);
                        continue;
                    }
                    if (aRec.fe_descriptor_id > 0) {
                        aDataSet2412.push_back(aRec);
                    }
                }
                MESSAGE("Perform - aDataSet2412.size() = " << aDataSet2412.size());
            }
            unvMesh2412::write_stream(out_stream, aDataSet2412);
        }

        // Material group (volume/face/edge)
        {
            using namespace unvMesh2417;
            if (myMatGroups.size() > 0)
            {
                TDataSet aDataSet2417;
                TGroupList::const_iterator aIter = myMatGroups.begin();
                for (; aIter != myMatGroups.end(); aIter++)
                {
                    WMESHDataGroupBase* aGroupDS = *aIter;
                    if (aGroupDS->get_group_type() != kMaterialMeshGroup)
                        continue;

                    UNVRecordData aRec;
                    aRec.GroupName = aGroupDS->get_internal_name();

                    int i;
                    WMDataElemIteratorPtr aIter = aGroupDS->get_elements();
                    if (aGroupDS->get_element_type() == WMDataAbs_Node) {
                        assert(0);
                        continue;
                    }

                    aRec.ElementList.resize(aGroupDS->get_size());
                    i = 0;
                    while (aIter->more())
                    {
                        const WMDataMeshElem* aElem = aIter->next();
                        aRec.ElementList[i] = aElem->get_elem_ID();
                        i++;
                    }

                    aDataSet2417.insert(TDataSet::value_type(aGroupDS->get_elem_ID(), aRec));
                }
                unvMesh2417::write_stream(out_stream, aDataSet2417);
                this->clearMatGroups();
            }
        }

        // Boundary condition group
        {
            using namespace unvMesh2417;
            if (myBCGroups.size() > 0)
            {
                TDataSet aDataSet2417;
                TGroupList::const_iterator aIter = myBCGroups.begin();
                for (; aIter != myBCGroups.end(); aIter++)
                {
                    WMESHDataGroupBase* aGroupDS = *aIter;
                    UNVRecordData aRec;
                    aRec.GroupName = aGroupDS->get_internal_name();

                    int i;
                    WMDataElemIteratorPtr aIter = aGroupDS->get_elements();
                    if (aGroupDS->get_element_type() == WMDataAbs_Node)
                    {
                        aRec.NodeList.resize(aGroupDS->get_size());
                        i = 0;
                        while (aIter->more())
                        {
                            const WMDataMeshElem* aElem = aIter->next();
                            aRec.NodeList[i] = aElem->get_elem_ID();
                            i++;
                        }
                    }
                    else
                    {
                        aRec.ElementList.resize(aGroupDS->get_size());
                        i = 0;
                        while (aIter->more())
                        {
                            const WMDataMeshElem* aElem = aIter->next();
                            aRec.ElementList[i] = aElem->get_elem_ID();
                            i++;
                        }
                    }
                    aDataSet2417.insert(TDataSet::value_type(aGroupDS->get_elem_ID(), aRec));
                }
                unvMesh2417::write_stream(out_stream, aDataSet2417);
                myBCGroups.clear();
            }
        }


        out_stream.flush();
        out_stream.close();
    }
    catch (const std::exception& exc) {
        INFOS("Exception:\n\t" << exc.what());
        throw;
    }
    catch (...) {
        INFOS("Unknown exception");
        throw;
    }
    return aResult;
}

TransferFemMesh::Status TransferUNVWMDataWriteMesh::executeAscii()
{
    Kernel_Utils::CLocalizer loc;
    Status aResult = DRS_OK;

    std::ofstream out_stream(myFile.c_str());

    try {
        // Unit System
        {
            using namespace unvMesh164;
            UNVRecordData aDataSet164;
            //aDataSet164.set_length_unit(myLengthUnit);
            unvMesh164::write_stream(out_stream, aDataSet164);
        }

        // Coordinate system
        unvMesh2420::write_stream(out_stream, myMeshName);

        // Nodes
        {
            using namespace unvMesh2411;
            TDataSet aDataSet2411;
            MESSAGE("Perform - myMesh->get_number_of_nodes() = " << myMesh->get_number_of_nodes());
            WMDataNodeIteratorPtr aNodesIter = myMesh->get_nodes_iterator();
            UNVRecordData aRec;
            while (aNodesIter->more())
            {
                const WMDataMeshNode* aNode = aNodesIter->next();
                aRec.label = aNode->get_elem_ID();
                aRec.coord[0] = aNode->X();
                aRec.coord[1] = aNode->Y();
                aRec.coord[2] = aNode->Z();
                aDataSet2411.push_back(aRec);
            }
            MESSAGE("Perform - aDataSet2411.size() = " << aDataSet2411.size());
            unvMesh2411::write_stream(out_stream, aDataSet2411);
        }

        // Elements
        {
            using namespace unvMesh2412;
            TDataSet aDataSet2412;

            // Storing WMDS Edges
            MESSAGE("Perform - myMesh->get_number_of_edges() = " << myMesh->get_number_of_edges());
            if (myMesh->get_number_of_edges()) {
                WMDataEdgeIteratorPtr anIter = myMesh->get_edges_iterator();
                while (anIter->more())
                {
                    const WMDataEdge* anElem = anIter->next();
                    int aNbNodes = anElem->get_number_of_nodes();
                    UNVRecordData aRec;
                    aRec.label = anElem->get_elem_ID();
                    aRec.mat_prop_tab_num = anElem->getComponentId();
                    aRec.node_labels.reserve(aNbNodes);
                    WMDataNodeIteratorPtr aNodesIter = anElem->get_nodes_iter();
                    if (anElem->is_quadratic_element()) {
                        aRec.fe_descriptor_id = 22;
                        const WMDataMeshNode* aNode0 = aNodesIter->next();
                        const WMDataMeshNode* aNode1 = aNodesIter->next();
                        const WMDataMeshNode* aNode2 = aNodesIter->next();
                        aRec.node_labels.push_back(aNode0->get_elem_ID());
                        aRec.node_labels.push_back(aNode2->get_elem_ID());
                        aRec.node_labels.push_back(aNode1->get_elem_ID());
                    }
                    else {
                        aRec.fe_descriptor_id = 11;
                        const WMDataMeshNode* aNode0 = aNodesIter->next();
                        const WMDataMeshNode* aNode1 = aNodesIter->next();
                        aRec.node_labels.push_back(aNode0->get_elem_ID());
                        aRec.node_labels.push_back(aNode1->get_elem_ID());
                    }
                    aDataSet2412.push_back(aRec);
                }
                MESSAGE("Perform - aDataSet2412.size() = " << aDataSet2412.size());
            }

            MESSAGE("Perform - myMesh->get_number_of_faces() = " << myMesh->get_number_of_faces());
            if (myMesh->get_number_of_faces())
            {
                WMDataFaceIteratorPtr anIter = myMesh->get_faces_iterator();
                while (anIter->more())
                {
                    const WMDataFace* anElem = anIter->next();
                    if (anElem->is_polyhedron()) continue;
                    int aNbNodes = anElem->get_number_of_nodes();
                    UNVRecordData aRec;
                    aRec.label = anElem->get_elem_ID();
                    aRec.mat_prop_tab_num = anElem->getComponentId();
                    aRec.node_labels.reserve(aNbNodes);
                    WMDataNodeIteratorPtr aNodesIter = anElem->get_nodes_iter();
                    switch (aNbNodes) {
                    case 3:
                    {
                        aRec.fe_descriptor_id = 41;
                        const WMDataMeshNode* aNode0 = aNodesIter->next();
                        const WMDataMeshNode* aNode1 = aNodesIter->next();
                        const WMDataMeshNode* aNode2 = aNodesIter->next();
                        aRec.node_labels.push_back(aNode0->get_elem_ID());
                        aRec.node_labels.push_back(aNode1->get_elem_ID());
                        aRec.node_labels.push_back(aNode2->get_elem_ID());
                    }
                    break;
                    case 4:
                    {
                        aRec.fe_descriptor_id = 44;
                        assert(0);
                        while (aNodesIter->more())
                        {
                            const WMDataMeshNode* aNode = aNodesIter->next();
                            aRec.node_labels.push_back(aNode->get_elem_ID());
                        }
                    }
                    break;
                    case 6:
                    {
                        aRec.fe_descriptor_id = 42;
                        const WMDataMeshNode* aNode0 = aNodesIter->next();
                        const WMDataMeshNode* aNode1 = aNodesIter->next();
                        const WMDataMeshNode* aNode2 = aNodesIter->next();
                        const WMDataMeshNode* aNode3 = aNodesIter->next();
                        const WMDataMeshNode* aNode4 = aNodesIter->next();
                        const WMDataMeshNode* aNode5 = aNodesIter->next();
                        aRec.node_labels.push_back(aNode0->get_elem_ID());
                        aRec.node_labels.push_back(aNode3->get_elem_ID());
                        aRec.node_labels.push_back(aNode1->get_elem_ID());
                        aRec.node_labels.push_back(aNode4->get_elem_ID());
                        aRec.node_labels.push_back(aNode2->get_elem_ID());
                        aRec.node_labels.push_back(aNode5->get_elem_ID());
                    }
                    break;
                    case 7:
                    {
                        aRec.fe_descriptor_id = 42;
                        assert(0);
                        while (aNodesIter->more())
                        {
                            const WMDataMeshNode* aNode = aNodesIter->next();
                            aRec.node_labels.push_back(aNode->get_elem_ID());
                        }
                    }
                    break;
                    case 8:
                    {
                        aRec.fe_descriptor_id = 45;
                        assert(0);
                        while (aNodesIter->more())
                        {
                            const WMDataMeshNode* aNode = aNodesIter->next();
                            aRec.node_labels.push_back(aNode->get_elem_ID());
                        }
                    }
                    break;
                    case 9:
                    {
                        aRec.fe_descriptor_id = 45;
                        aRec.node_labels.resize(8);
                        while (aNodesIter->more())
                        {
                            const WMDataMeshNode* aNode = aNodesIter->next();
                            aRec.node_labels.push_back(aNode->get_elem_ID());
                        }
                        assert(0);
                    }
                    break;
                    default:
                        assert(0);
                        continue;
                    }
                    if (aRec.fe_descriptor_id > 0) {
                        aDataSet2412.push_back(aRec);
                    }
                }
                MESSAGE("Perform - aDataSet2412.size() = " << aDataSet2412.size());
            }

            MESSAGE("Perform - myMesh->get_volumes_count() = " << myMesh->get_volumes_count());
            if (myMesh->get_volumes_count() > 0)
            {
                WMDataVolumeIteratorPtr anIter = myMesh->get_volumes_iterator();
                while (anIter->more())
                {
                    const WMDataVolume* anElem = anIter->next();
                    if (anElem->is_polyhedron())
                        continue;
                    UNVRecordData aRec;
                    aRec.label = anElem->get_elem_ID();
                    aRec.mat_prop_tab_num = anElem->getComponentId();
                    int aNbNodes = anElem->get_number_of_nodes();
                    aRec.node_labels.reserve(aNbNodes);
                    WMDataNodeIteratorPtr aNodesIter = anElem->get_nodes_iter();
                    switch (aNbNodes) {
                    case 4:
                    {
                        aRec.fe_descriptor_id = 111;
                        const WMDataMeshNode* aNode0 = aNodesIter->next();
                        const WMDataMeshNode* aNode1 = aNodesIter->next();
                        const WMDataMeshNode* aNode2 = aNodesIter->next();
                        const WMDataMeshNode* aNode3 = aNodesIter->next();
                        aRec.node_labels.push_back(aNode0->get_elem_ID());
                        aRec.node_labels.push_back(aNode1->get_elem_ID());
                        aRec.node_labels.push_back(aNode2->get_elem_ID());
                        aRec.node_labels.push_back(aNode3->get_elem_ID());
                    }
                    break;
                    case 6:
                    {
                        aRec.fe_descriptor_id = 112;
                        assert(0);
                        while (aNodesIter->more() && aRec.node_labels.size() < aNbNodes)
                        {
                            const WMDataMeshElem* aNode = aNodesIter->next();
                            aRec.node_labels.push_back(aNode->get_elem_ID());
                        }
                    }
                    break;
                    case 8:
                    {
                        aRec.fe_descriptor_id = 115;
                        assert(0);
                        while (aNodesIter->more() && aRec.node_labels.size() < aNbNodes)
                        {
                            const WMDataMeshElem* aNode = aNodesIter->next();
                            aRec.node_labels.push_back(aNode->get_elem_ID());
                        }
                    }
                    break;
                    case 10:
                    {
                        aRec.fe_descriptor_id = 118;
                        const WMDataMeshNode* aNode0 = aNodesIter->next();
                        const WMDataMeshNode* aNode1 = aNodesIter->next();
                        const WMDataMeshNode* aNode2 = aNodesIter->next();
                        const WMDataMeshNode* aNode3 = aNodesIter->next();
                        const WMDataMeshNode* aNode4 = aNodesIter->next();
                        const WMDataMeshNode* aNode5 = aNodesIter->next();
                        const WMDataMeshNode* aNode6 = aNodesIter->next();
                        const WMDataMeshNode* aNode7 = aNodesIter->next();
                        const WMDataMeshNode* aNode8 = aNodesIter->next();
                        const WMDataMeshNode* aNode9 = aNodesIter->next();
                        aRec.node_labels.push_back(aNode0->get_elem_ID());
                        aRec.node_labels.push_back(aNode4->get_elem_ID());
                        aRec.node_labels.push_back(aNode1->get_elem_ID());
                        aRec.node_labels.push_back(aNode5->get_elem_ID());
                        aRec.node_labels.push_back(aNode2->get_elem_ID());
                        aRec.node_labels.push_back(aNode6->get_elem_ID());

                        aRec.node_labels.push_back(aNode7->get_elem_ID());
                        aRec.node_labels.push_back(aNode8->get_elem_ID());
                        aRec.node_labels.push_back(aNode9->get_elem_ID());
                        aRec.node_labels.push_back(aNode3->get_elem_ID());
                    }
                    break;
                    case 13:
                    {
                        aRec.fe_descriptor_id = 114;
                        assert(0);
                        while (aNodesIter->more() && aRec.node_labels.size() < aNbNodes)
                        {
                            const WMDataMeshElem* aNode = aNodesIter->next();
                            aRec.node_labels.push_back(aNode->get_elem_ID());
                        }
                    }
                    break;
                    case 15:
                    {
                        aRec.fe_descriptor_id = 113;
                        assert(0);
                        while (aNodesIter->more() && aRec.node_labels.size() < aNbNodes)
                        {
                            const WMDataMeshElem* aNode = aNodesIter->next();
                            aRec.node_labels.push_back(aNode->get_elem_ID());
                        }
                    }
                    break;
                    case 20:
                    case 27:
                    {
                        aRec.fe_descriptor_id = 116;
                        aNbNodes = 20;
                        assert(0);
                        while (aNodesIter->more() && aRec.node_labels.size() < aNbNodes)
                        {
                            const WMDataMeshElem* aNode = aNodesIter->next();
                            aRec.node_labels.push_back(aNode->get_elem_ID());
                        }
                    }
                    break;
                    default:
                        assert(0);
                        continue;
                    }
                    if (aRec.fe_descriptor_id > 0) {
                        aDataSet2412.push_back(aRec);
                    }
                }
                MESSAGE("Perform - aDataSet2412.size() = " << aDataSet2412.size());
            }
            unvMesh2412::write_stream(out_stream, aDataSet2412);
        }

        // Material group (volume/face/edge)
        {
            using namespace unvMesh2417;
            if (myMatGroups.size() > 0)
            {
                TDataSet aDataSet2417;
                TGroupList::const_iterator aIter = myMatGroups.begin();
                for (; aIter != myMatGroups.end(); aIter++)
                {
                    WMESHDataGroupBase* aGroupDS = *aIter;
                    if (aGroupDS->get_group_type() != kMaterialMeshGroup)
                        continue;

                    UNVRecordData aRec;
                    aRec.GroupName = aGroupDS->get_internal_name();

                    int i;
                    WMDataElemIteratorPtr aIter = aGroupDS->get_elements();
                    if (aGroupDS->get_element_type() == WMDataAbs_Node) {
                        assert(0);
                        continue;
                    }

                    aRec.ElementList.resize(aGroupDS->get_size());
                    i = 0;
                    while (aIter->more())
                    {
                        const WMDataMeshElem* aElem = aIter->next();
                        aRec.ElementList[i] = aElem->get_elem_ID();
                        i++;
                    }

                    aDataSet2417.insert(TDataSet::value_type(aGroupDS->get_elem_ID(), aRec));
                }
                unvMesh2417::write_stream(out_stream, aDataSet2417);
                this->clearMatGroups();
            }
        }

        // bc group
        {
            using namespace unvMesh2417;
            if (myBCGroups.size() > 0)
            {
                TDataSet aDataSet2417;
                TGroupList::const_iterator aIter = myBCGroups.begin();
                for (; aIter != myBCGroups.end(); aIter++)
                {
                    WMESHDataGroupBase* aGroupDS = *aIter;
                    UNVRecordData aRec;
                    aRec.GroupName = aGroupDS->get_internal_name();

                    int i;
                    WMDataElemIteratorPtr aIter = aGroupDS->get_elements();
                    if (aGroupDS->get_element_type() == WMDataAbs_Node)
                    {
                        aRec.NodeList.resize(aGroupDS->get_size());
                        i = 0;
                        while (aIter->more())
                        {
                            const WMDataMeshElem* aElem = aIter->next();
                            aRec.NodeList[i] = aElem->get_elem_ID();
                            i++;
                        }
                    }
                    else
                    {
                        aRec.ElementList.resize(aGroupDS->get_size());
                        i = 0;
                        while (aIter->more())
                        {
                            const WMDataMeshElem* aElem = aIter->next();
                            aRec.ElementList[i] = aElem->get_elem_ID();
                            i++;
                        }
                    }
                    aDataSet2417.insert(TDataSet::value_type(aGroupDS->get_elem_ID(), aRec));
                }
                unvMesh2417::write_stream(out_stream, aDataSet2417);
                myBCGroups.clear();
            }
        }

        out_stream.flush();
        out_stream.close();

        if (!check_file(myFile)) {
            EXCEPTION(runtime_error, "Error: bad UNV output file at the final check.");
        }
    }
    catch (const std::exception& exc) {
        INFOS("Exception:\n\t" << exc.what());
        throw;
    }
    catch (...) {
        INFOS("Unknown exception");
        throw;
    }
    return aResult;
}


void TransferUNVWMDataWriteMesh::clearMatGroups()
{
    //const TGroupList::iterator it_end = myMatGroups.end();
    //for (TGroupList::iterator it = myMatGroups.begin(); it != it_end; ++it)
    //{
    //    delete *it;
    //}
    myMatGroups.clear();
}
void TransferUNVWMDataWriteMesh::clearBCGroups()
{
    //const TGroupList::iterator it_end = myBCGroups.end();
    //for (TGroupList::iterator it = myBCGroups.begin(); it != it_end; ++it)
    //{
    //    delete *it;
    //}
    myBCGroups.clear();
}


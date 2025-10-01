#include "mesh_helper.h"
#include "unstructured_mesh2.h"

#include "base_particle_dynamics.h"

#include <filesystem>
#include <algorithm>
#include "TransferUNV/TransferUNV164DS.h"
#include "TransferUNV/TransferUNV2411DS.h"
#include "TransferUNV/TransferUNV2412DS.h"
#include "TransferUNV/TransferUNV2417DS.h"
#include "TransferUNV/TransferUNV2420DS.h"
#include "TransferUNV/TransferUNVTools.h"


namespace fs = std::filesystem;

namespace SPH
{
template <typename T>
bool IsSubset(std::vector<T>& A, std::vector<T>& B)
{
    std::sort(A.begin(), A.end());
    std::sort(B.begin(), B.end());
    return std::includes(A.begin(), A.end(), B.begin(), B.end());
}

int IsTetSubset(std::vector<size_t> largeVec, std::vector<size_t> smallVec)
{
    std::sort(smallVec.begin(), smallVec.end());

    std::vector<size_t> newVec = {largeVec[0], largeVec[1], largeVec[2]};
    std::sort(newVec.begin(), newVec.end());
    if (newVec == smallVec) { // std::includes(A.begin(), A.end(), B.begin(), B.end());
        return 0;
    }

    newVec = {largeVec[0], largeVec[1], largeVec[3]};
    std::sort(newVec.begin(), newVec.end());
    if (newVec == smallVec) {
        return 1;
    }

    newVec = {largeVec[1], largeVec[2], largeVec[3]};
    std::sort(newVec.begin(), newVec.end());
    if (newVec == smallVec)
    {
        return 2;
    }

    newVec = {largeVec[0], largeVec[2], largeVec[3]};
    std::sort(newVec.begin(), newVec.end());
    if (newVec == smallVec)
    {
        return 3;
    }

    return -1;
}

bool containsAllElements(const std::vector<size_t> &smaller_vec, const std::vector<size_t> &larger_vec)
{
    // If the smaller vector is larger than the larger vector, it cannot be included.
    if (smaller_vec.size() > larger_vec.size())
    {
        return false;
    }

    // Check if all elements of smaller_vec are present in larger_vec
    return std::all_of(smaller_vec.begin(), smaller_vec.end(),
                       [&](int element)
                       {
                           return std::find(larger_vec.begin(), larger_vec.end(), element) != larger_vec.end();
                       });
}


/*!
 * \brief Move node coordinates to the global Cartesian CS
 */
void transformNodes(unvMesh2411::TDataSet::const_iterator fromNode,
                    unvMesh2411::TDataSet::const_iterator endNode,
                    const unvMesh2420::UNVRecordData &csRecord)
{
    const int csLabel = fromNode->exp_coord_sys_num;

    unvMesh2411::TDataSet::const_iterator nodeIt;

    // apply Transformation Matrix
    if (!csRecord.is_identity_tensor())
    {
        for (nodeIt = fromNode; nodeIt != endNode; ++nodeIt)
        {
            const unvMesh2411::UNVRecordData &nodeRec = *nodeIt;
            if (nodeRec.exp_coord_sys_num == csLabel)
                csRecord.modify_by_tensor((double *)nodeRec.coord);
        }
    }

    // transform from Cylindrical CS
    if (csRecord.coord_sys_type == unvMesh2420::Cylindrical)
    {
        for (nodeIt = fromNode; nodeIt != endNode; ++nodeIt)
        {
            const unvMesh2411::UNVRecordData &nodeRec = *nodeIt;
            if (nodeRec.exp_coord_sys_num == csLabel)
                csRecord.convert_by_CCS((double *)nodeRec.coord);
        }
    }
    // transform from Spherical CS
    else if (csRecord.coord_sys_type == unvMesh2420::Spherical)
    {
        for (nodeIt = fromNode; nodeIt != endNode; ++nodeIt)
        {
            const unvMesh2411::UNVRecordData &nodeRec = *nodeIt;
            if (nodeRec.exp_coord_sys_num == csLabel)
                csRecord.convert_by_SCS((double *)nodeRec.coord);
        }
    }
}

ANSYSMesh2::ANSYSMesh2()
{
}



//=================================================================================================//
ANSYSMesh2::ANSYSMesh2(const std::string &full_path)
{
    getDataFromMeshFile(full_path);
    getElementCenterCoordinates();
    getMinimumDistanceBetweenNodes();
}
//=================================================================================================//
void ANSYSMesh2::getDataFromMeshFile(const std::string &full_path)
{
    Real ICEM = 0;
    std::ifstream mesh_file; /*!< \brief File object for the Ansys ASCII mesh file. */
    mesh_file.open(full_path);
    if (mesh_file.fail())
    {
        std::cout << "Error:Check if the file exists." << std::endl;
    }

    std::string suffix = fs::path(full_path).filename().extension().string();
    transform(suffix.begin(), suffix.end(), suffix.begin(), ::tolower);
    if (suffix == ".unv" )
    {
        getDataFromUnvFile(full_path);
        return;
    }


    std::string text_line;
    /*Check mesh file generation software*/
    (getline(mesh_file, text_line));
    text_line.erase(0, 4);
    std::istringstream value(text_line);
    if (text_line.find("Created by", 0) != std::string::npos)
    {
        ICEM = 1;
        std::cout << "Reading 3D ICEM mesh." << std::endl;
    }

    if (ICEM == 1)
    {

        /*--- Read the dimension of the mesh ---*/
        size_t dimension(0);
        MeshFileHelpers::meshDimension(mesh_file, dimension, text_line);

        /*--- Read the node data (index is starting from zero) ---*/
        size_t number_of_points(0);
        MeshFileHelpers::numberOfNodes(mesh_file, number_of_points, text_line);
        MeshFileHelpers::nodeCoordinates(mesh_file, node_coordinates_, text_line, dimension);

        size_t boundary_type(0);
        size_t number_of_elements(0);
        size_t mesh_type = 4;
        MeshFileHelpers::numberOfElements(mesh_file, number_of_elements, text_line);

        /*Preparing and initializing the data structure of mesh topology and element node connection*/
        MeshFileHelpers::dataStruct(mesh_topology_, elements_nodes_connection_, number_of_elements, mesh_type, dimension);

        /*--- find the elements lines ---*/
        while (getline(mesh_file, text_line))
        {
            if (text_line.find("(13", 0) != std::string::npos && text_line.find(")(", 0) != std::string::npos)
            {
                boundary_type = MeshFileHelpers::findBoundaryType(text_line, boundary_type);
                types_of_boundary_condition_.push_back(boundary_type);
                while (getline(mesh_file, text_line))
                {
                    if (text_line.find(")", 0) == std::string::npos)
                    {
                        Vecd nodes = MeshFileHelpers::nodeIndex(text_line);
                        Vec2d cells = MeshFileHelpers::cellIndex(text_line);

                        /*--- build up all topology---*/
                        bool check_neighbor_cell1 = 1;
                        bool check_neighbor_cell2 = 0;
                        for (int cell1_cell2 = 0; cell1_cell2 != cells.size(); ++cell1_cell2)
                        {
                            if (mesh_type == 4)
                            {
                                if (cells[check_neighbor_cell2] != 0)
                                {
                                    MeshFileHelpers::updateElementsNodesConnection(elements_nodes_connection_, nodes, cells, check_neighbor_cell1, check_neighbor_cell2);
                                    MeshFileHelpers::updateCellLists(mesh_topology_, elements_nodes_connection_, nodes, cells, check_neighbor_cell1, check_neighbor_cell2, boundary_type);
                                    MeshFileHelpers::updateBoundaryCellLists(mesh_topology_, elements_nodes_connection_, nodes, cells, check_neighbor_cell1, check_neighbor_cell2, boundary_type);
                                }
                                if (cells[check_neighbor_cell2] == 0)
                                {
                                    break;
                                }
                            }
                        }
                    }
                    else
                        break;
                }
            }
            if (text_line.find("Zone Sections", 0) != std::string::npos)
                break;
            if (text_line.find(")") != std::string::npos)
                continue;
        }
        mesh_topology_.erase(mesh_topology_.begin());
    }
    else /*This section is for mesh files created from fluent*/
    {
        std::cout << "Reading 3D Fluent mesh." << std::endl;
        size_t dimension(0);
        MeshFileHelpers::meshDimension(mesh_file, dimension, text_line);

        /*--- Read the node data (index is starting from zero) ---*/
        size_t number_of_points(0);
        MeshFileHelpers::numberOfNodesFluent(mesh_file, number_of_points, text_line);

        size_t boundary_type(0);
        size_t number_of_elements(0);
        size_t mesh_type = 4;

        MeshFileHelpers::numberOfElementsFluent(mesh_file, number_of_elements, text_line);
        MeshFileHelpers::dataStruct(mesh_topology_, elements_nodes_connection_, number_of_elements, mesh_type, dimension);
        MeshFileHelpers::nodeCoordinatesFluent(mesh_file, node_coordinates_, text_line, dimension);

        while (getline(mesh_file, text_line))
        {

            if (text_line.find("(13", 0) != std::string::npos && text_line.find(") (", 0) != std::string::npos)
            {
                boundary_type = MeshFileHelpers::findBoundaryType(text_line, boundary_type);
                types_of_boundary_condition_.push_back(boundary_type);
                while (getline(mesh_file, text_line))
                {

                    if (text_line.find(")", 0) == std::string::npos)
                    {
                        Vecd nodes = MeshFileHelpers::nodeIndex(text_line);
                        Vec2d cells = MeshFileHelpers::cellIndex(text_line);
                        /*--- build up all topology---*/
                        bool check_neighbor_cell1 = 1;
                        bool check_neighbor_cell2 = 0;
                        for (int cell1_cell2 = 0; cell1_cell2 != cells.size(); ++cell1_cell2)
                        {
                            if (mesh_type == 4)
                            {
                                if (cells[check_neighbor_cell2] != 0)
                                {
                                    MeshFileHelpers::updateElementsNodesConnection(elements_nodes_connection_, nodes, cells, check_neighbor_cell1, check_neighbor_cell2);
                                    MeshFileHelpers::updateCellLists(mesh_topology_, elements_nodes_connection_, nodes, cells, check_neighbor_cell1,
                                                                     check_neighbor_cell2, boundary_type);
                                    MeshFileHelpers::updateBoundaryCellListsFluent(mesh_topology_, elements_nodes_connection_, nodes, cells,
                                                                                   check_neighbor_cell1, check_neighbor_cell2, boundary_type);
                                }
                                if (cells[check_neighbor_cell2] == 0)
                                {
                                    break;
                                }
                            }
                        }
                    }
                    else
                        break;
                }
            }
            if (text_line.find("(12", 0) != std::string::npos && text_line.find("))", 0) != std::string::npos)
                break;
            if (text_line.find(")") != std::string::npos)
                continue;
        }
        mesh_topology_.erase(mesh_topology_.begin());
    }
}
//=================================================================================================//

void ANSYSMesh2::getDataFromUnvFile(const std::string &full_path)
{
    Real ICEM = 0;
    std::ifstream in_stream; /*!< \brief File object for the Ansys ASCII mesh file. */
    in_stream.open(full_path);
    if (in_stream.fail())
    {
        std::cout << "Error:Check if the file exists." << std::endl;
    }

    try
    {
        size_t dimension = 3;
        size_t mesh_type = 4; // unstructured, tet4

        {
            // Read Units
            unvMesh164::UNVRecordData aUnitsRecord;
            unvMesh164::read_stream(in_stream, aUnitsRecord);

            // Read Coordinate systems
            std::string myMeshName;
            unvMesh2420::TDataSet aCoordSysDataSet;
            unvMesh2420::read_stream(in_stream, myMeshName, aCoordSysDataSet);

            // Read nodes
            using namespace unvMesh2411;
            TDataSet aDataSet2411;
            unvMesh2411::read_stream(in_stream, aDataSet2411);

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
                    unvMesh2411::UNVRecordData &nodeRec = *nodeIter;
                    nodeRec.coord[0] *= lenFactor;
                    nodeRec.coord[1] *= lenFactor;
                    nodeRec.coord[2] *= lenFactor;
                }
            }

            // Create nodes in the mesh
            Vecd Coords = {0., 0., 0.};
            node_coordinates_.push_back(Coords); // dummy for 1-based vector
            TDataSet::const_iterator anIter = aDataSet2411.begin();
            for (; anIter != aDataSet2411.end(); anIter++)
            {
                const UNVRecordData &aRec = *anIter;
                // myMesh->add_node_and_id(aRec.coord[0], aRec.coord[1], aRec.coord[2], aRec.label);
                Vecd Coords = Vecd::Zero();
                for (unsigned int i = 0; i < 3; ++i)
                    Coords[i] = aRec.coord[i];
                node_coordinates_.push_back(Coords);
            }
        }

       

        {
            using namespace unvMesh2412;
            TDataSet aDataSet2412;
            unvMesh2412::read_stream(in_stream, aDataSet2412);
            unsigned int number_of_elements = 0;
            for (TDataSet::const_iterator anIter = aDataSet2412.begin(); anIter != aDataSet2412.end(); anIter++)
            {
                const UNVRecordData &aRec = *anIter;
                if (is_body_element(aRec.fe_descriptor_id))
                    number_of_elements++;
            }
            /*Preparing and initializing the data structure of mesh topology and element node connection*/
            MeshFileHelpers::dataStruct(mesh_topology_, elements_nodes_connection_, number_of_elements, mesh_type, dimension);

            unsigned long long count = 1; // 1-based
            for (TDataSet::const_iterator anIter = aDataSet2412.begin(); anIter != aDataSet2412.end(); anIter++)
            {
                //WMDataMeshElem *anElement = NULL;
                const UNVRecordData &aRec = *anIter;
                if (is_beam_element(aRec.fe_descriptor_id))
                {
                    continue;
                    switch (aRec.node_labels.size())
                    {
                    case 2:
                        elements_nodes_connection_.push_back({aRec.node_labels[0], aRec.node_labels[1]});
                        break;
                    case 3:
                        elements_nodes_connection_.push_back({aRec.node_labels[0], aRec.node_labels[1], aRec.node_labels[2]});
                        break;
                    default:
                        break;
                    }
                }
                else if (is_face_element(aRec.fe_descriptor_id))
                {
                    continue;
                    switch (aRec.fe_descriptor_id)
                    {
                    case 41: // Plane Stress Linear Triangle
                    case 51: // Plane Strain Linear Triangle
                    case 61: // Plate Linear Triangle
                    case 74: // Membrane Linear Triangle
                    case 81: // Axisymetric Solid Linear Triangle
                    case 91: // Thin Shell Linear Triangle
                        elements_nodes_connection_.push_back({aRec.node_labels[0], aRec.node_labels[1], aRec.node_labels[2]});
                        break;

                    case 42: //  Plane Stress Parabolic Triangle
                    case 52: //  Plane Strain Parabolic Triangle
                    case 62: //  Plate Parabolic Triangle
                    case 72: //  Membrane Parabolic Triangle
                    case 82: //  Axisymetric Solid Parabolic Triangle
                    case 92: //  Thin Shell Parabolic Triangle
                        elements_nodes_connection_.push_back({ aRec.node_labels[0], aRec.node_labels[2], aRec.node_labels[4],
                                                               aRec.node_labels[1], aRec.node_labels[3], aRec.node_labels[5]});
                        break;

                    case 44: // Plane Stress Linear Quadrilateral
                    case 54: // Plane Strain Linear Quadrilateral
                    case 64: // Plate Linear Quadrilateral
                    case 71: // Membrane Linear Quadrilateral
                    case 84: // Axisymetric Solid Linear Quadrilateral
                    case 94: // Thin Shell Linear Quadrilateral
                        elements_nodes_connection_.push_back({aRec.node_labels[0], aRec.node_labels[1], aRec.node_labels[2], aRec.node_labels[3]});
                        break;

                    case 45: // Plane Stress Parabolic Quadrilateral
                    case 55: // Plane Strain Parabolic Quadrilateral
                    case 65: // Plate Parabolic Quadrilateral
                    case 75: // Membrane Parabolic Quadrilateral
                    case 85: // Axisymetric Solid Parabolic Quadrilateral
                    case 95: // Thin Shell Parabolic Quadrilateral
                        if (aRec.node_labels.size() == 9)
                        {
                            assert(0);
                            // anElement = myMesh->add_face_and_id(aRec.node_labels[0],
                            //     aRec.node_labels[2],
                            //     aRec.node_labels[4],
                            //     aRec.node_labels[6],
                            //     aRec.node_labels[1],
                            //     aRec.node_labels[3],
                            //     aRec.node_labels[5],
                            //     aRec.node_labels[7],
                            //     aRec.node_labels[8],
                            //     aRec.label);
                        }
                        else
                        {
                            elements_nodes_connection_.push_back({aRec.node_labels[0], aRec.node_labels[2], aRec.node_labels[4], aRec.node_labels[6],
                                                                  aRec.node_labels[1], aRec.node_labels[3], aRec.node_labels[5], aRec.node_labels[7]});
                        }
                        break;
                    }
                }
                else if (is_body_element(aRec.fe_descriptor_id))
                {
                    switch (aRec.fe_descriptor_id)
                    {
                    case 111: // Solid Linear Tetrahedron - TET4
                    {
                        std::vector<size_t> elemNodes = {aRec.node_labels[0], aRec.node_labels[1], aRec.node_labels[2], aRec.node_labels[3]};
                        elements_nodes_connection_[count] = elemNodes;

                        // fill the first cell data
                        std::vector<std::vector<size_t>> faces;
                        faces.push_back({elemNodes[0], elemNodes[1], elemNodes[2]});
                        faces.push_back({elemNodes[0], elemNodes[1], elemNodes[3]});
                        faces.push_back({elemNodes[1], elemNodes[2], elemNodes[3]});
                        faces.push_back({elemNodes[0], elemNodes[2], elemNodes[3]});
                        //if (count == 1)
                        {
                            mesh_topology_[1][0] = {0, 3, faces[0][0], faces[0][1], faces[0][2]};
                            mesh_topology_[1][1] = {0, 3, faces[1][0], faces[1][1], faces[1][2]};
                            mesh_topology_[1][2] = {0, 3, faces[2][0], faces[2][1], faces[2][2]};
                            mesh_topology_[1][3] = {0, 3, faces[3][0], faces[3][1], faces[3][2]};
                        }

                        //parallel_for(IndexRange(1, count), [&](const IndexRange &r)
                        //{
                        //    for (unsigned int i = r.begin(); i != r.end(); ++i) // all processed elements
                        //    {
                        //        const std::vector<size_t> &curElemNodes = elements_nodes_connection_[i];
                        //        for (unsigned int j = 0; j < mesh_topology_[count].size(); ++j) // 4 faces for tet
                        //        {
                        //            if (mesh_topology_[i][j][1] == 2) // already found the neighber element
                        //                continue;

                        //            int idxSub = IsTetSubset(curElemNodes, faces[j]);
                        //            if (idxSub == -1)
                        //            {
                        //                mesh_topology_[count][j] = {0, 3, faces[j][0], faces[j][1], faces[j][2]};
                        //            }
                        //            else
                        //            {
                        //                mesh_topology_[count][j] = {i, 2, faces[j][0], faces[j][1], faces[j][2]};
                        //                mesh_topology_[i][idxSub][0] = count;
                        //                mesh_topology_[i][idxSub][1] = 2;
                        //            }
                        //        }
                        //    } 
                        //},
                        //ap);
                    }
                    break;
                    case 118: // Solid Quadratic Tetrahedron - TET10
                    {
                        assert(0);
                        elements_nodes_connection_[count] = {aRec.node_labels[0], aRec.node_labels[2], aRec.node_labels[4],
                                                              aRec.node_labels[9],
                                                              aRec.node_labels[1], aRec.node_labels[3], aRec.node_labels[5],
                                                              aRec.node_labels[6], aRec.node_labels[7], aRec.node_labels[8]};
                        mesh_topology_[count][0] = {0, 2, aRec.node_labels[0], aRec.node_labels[1], aRec.node_labels[2]};
                        mesh_topology_[count][1] = {0, 2, aRec.node_labels[0], aRec.node_labels[1], aRec.node_labels[3]};
                        mesh_topology_[count][2] = {0, 2, aRec.node_labels[1], aRec.node_labels[2], aRec.node_labels[3]};
                        mesh_topology_[count][3] = {0, 2, aRec.node_labels[0], aRec.node_labels[2], aRec.node_labels[3]};
                    }
                    break;
                    case 112: // Solid Linear Prism - PRISM6
                    {
                        assert(0);
                        elements_nodes_connection_[count] = {aRec.node_labels[0], aRec.node_labels[2], aRec.node_labels[1],
                                                              aRec.node_labels[3], aRec.node_labels[5], aRec.node_labels[4]};
                        //mesh_topology_[count][0] = {0, 2, aRec.node_labels[0], aRec.node_labels[1], aRec.node_labels[2], aRec.node_labels[3]};
                        //mesh_topology_[count][1] = {0, 2, aRec.node_labels[0], aRec.node_labels[1], aRec.node_labels[4], aRec.node_labels[5]};
                        //mesh_topology_[count][2] = {0, 2, aRec.node_labels[1], aRec.node_labels[2], aRec.node_labels[5], aRec.node_labels[6]};
                        //mesh_topology_[count][3] = {0, 2, aRec.node_labels[2], aRec.node_labels[3], aRec.node_labels[6], aRec.node_labels[7]};
                        //mesh_topology_[count][4] = {0, 2, aRec.node_labels[0], aRec.node_labels[3], aRec.node_labels[4], aRec.node_labels[7]};
                        //mesh_topology_[count][5] = {0, 2, aRec.node_labels[4], aRec.node_labels[5], aRec.node_labels[6], aRec.node_labels[7]};
                    }
                    break;
                    case 113: // Solid Quadratic Prism - PRISM15
                    {
                        assert(0);
                        elements_nodes_connection_[count] = {aRec.node_labels[0], aRec.node_labels[4], aRec.node_labels[2],
                                                              aRec.node_labels[9], aRec.node_labels[13], aRec.node_labels[11],
                                                              aRec.node_labels[5], aRec.node_labels[3], aRec.node_labels[1],
                                                              aRec.node_labels[14], aRec.node_labels[12], aRec.node_labels[10],
                                                              aRec.node_labels[6], aRec.node_labels[8], aRec.node_labels[7]};
                        //mesh_topology_[count] = {};
                    }
                    break;
                    case 115: // Solid Linear Brick - HEX8
                    {
                        elements_nodes_connection_[count] = {aRec.node_labels[0], aRec.node_labels[3], aRec.node_labels[2],
                                                              aRec.node_labels[1], aRec.node_labels[4], aRec.node_labels[7],
                                                              aRec.node_labels[6], aRec.node_labels[5]};
                        mesh_topology_[count][0] = {0, 2, aRec.node_labels[0], aRec.node_labels[1], aRec.node_labels[2], aRec.node_labels[3]};
                        mesh_topology_[count][1] = {0, 2, aRec.node_labels[0], aRec.node_labels[1], aRec.node_labels[4], aRec.node_labels[5]};
                        mesh_topology_[count][2] = {0, 2, aRec.node_labels[1], aRec.node_labels[2], aRec.node_labels[5], aRec.node_labels[6]};
                        mesh_topology_[count][3] = {0, 2, aRec.node_labels[2], aRec.node_labels[3], aRec.node_labels[6], aRec.node_labels[7]};
                        mesh_topology_[count][4] = {0, 2, aRec.node_labels[0], aRec.node_labels[3], aRec.node_labels[4], aRec.node_labels[7]};
                        mesh_topology_[count][5] = {0, 2, aRec.node_labels[4], aRec.node_labels[5], aRec.node_labels[6], aRec.node_labels[7]};
                    }
                    break;
                    case 116: // Solid Quadratic Brick - HEX20
                    {
                        elements_nodes_connection_[count] = {aRec.node_labels[0], aRec.node_labels[6], aRec.node_labels[4], aRec.node_labels[2],
                                                              aRec.node_labels[12], aRec.node_labels[18], aRec.node_labels[16], aRec.node_labels[14],
                                                              aRec.node_labels[7], aRec.node_labels[5], aRec.node_labels[3], aRec.node_labels[1],
                                                              aRec.node_labels[19], aRec.node_labels[17], aRec.node_labels[15], aRec.node_labels[13],
                                                              aRec.node_labels[8], aRec.node_labels[11], aRec.node_labels[10], aRec.node_labels[9]};
                        mesh_topology_[count][0] = {0, 2, aRec.node_labels[0], aRec.node_labels[1], aRec.node_labels[2], aRec.node_labels[3]};
                        mesh_topology_[count][1] = {0, 2, aRec.node_labels[0], aRec.node_labels[1], aRec.node_labels[4], aRec.node_labels[5]};
                        mesh_topology_[count][2] = {0, 2, aRec.node_labels[1], aRec.node_labels[2], aRec.node_labels[5], aRec.node_labels[6]};
                        mesh_topology_[count][3] = {0, 2, aRec.node_labels[2], aRec.node_labels[3], aRec.node_labels[6], aRec.node_labels[7]};
                        mesh_topology_[count][4] = {0, 2, aRec.node_labels[0], aRec.node_labels[3], aRec.node_labels[4], aRec.node_labels[7]};
                        mesh_topology_[count][5] = {0, 2, aRec.node_labels[4], aRec.node_labels[5], aRec.node_labels[6], aRec.node_labels[7]};
                    }
                    break;
                    case 114: // pyramid of 13 nodes (quadratic) - PIRA13
                    {
                        assert(0);
                        elements_nodes_connection_[count] = {aRec.node_labels[0], aRec.node_labels[6], aRec.node_labels[4], aRec.node_labels[2],
                                                              aRec.node_labels[7], aRec.node_labels[5], aRec.node_labels[3], aRec.node_labels[1],
                                                              aRec.node_labels[8], aRec.node_labels[11], aRec.node_labels[10], aRec.node_labels[9], aRec.node_labels[12]};
                    }
                    break;
                    default:
                        break;
                    }
                    count++;
                }
                //if (!anElement)
                //{
                    //assert(0);
                    //MESSAGE("TransferUNVWMDataReadMesh::Perform - can not add element with ID = " << aRec.label << " and type = " << aRec.fe_descriptor_id);
                //}
            }
        }
    }
    catch (std::exception &e)
    {
        std::cerr << "error: " << e.what() << "\n";
        exit(1);
    }
    catch (...)
    {
        std::cerr << "Exception of unknown type!\n";
        mesh_topology_.erase(mesh_topology_.begin());
    }

    return;
}


void ANSYSMesh2::updateDataForParticleGen()
{
    this->getElementCenterCoordinates();
    this->getMinimumDistanceBetweenNodes();
}


void ANSYSMesh2::getElementCenterCoordinates()
{
    elements_centroids_.resize(elements_nodes_connection_.size());
    elements_volumes_.resize(elements_nodes_connection_.size());
    for (std::size_t element = 1; element != elements_nodes_connection_.size(); ++element)
    {
        Vecd center_coordinate = Vecd::Zero();
        MeshFileHelpers::cellCenterCoordinates(elements_nodes_connection_, element, node_coordinates_, elements_centroids_, center_coordinate);
        MeshFileHelpers::elementVolume(elements_nodes_connection_, element, node_coordinates_, elements_volumes_);
    }
    elements_volumes_.erase(elements_volumes_.begin());
    elements_centroids_.erase(elements_centroids_.begin());
    elements_nodes_connection_.erase(elements_nodes_connection_.begin());
}
//=================================================================================================//

void ANSYSMesh2::getMinimumDistanceBetweenNodes()
{
    StdVec<Real> all_data_of_distance_between_nodes;
    all_data_of_distance_between_nodes.resize(0);
    MeshFileHelpers::minimumDistance(all_data_of_distance_between_nodes, elements_volumes_, mesh_topology_, node_coordinates_);
    auto min_distance_iter = std::min_element(all_data_of_distance_between_nodes.begin(), all_data_of_distance_between_nodes.end());
    if (min_distance_iter != all_data_of_distance_between_nodes.end())
    {
        min_distance_between_nodes_ = *min_distance_iter;
    }
    else
    {
        std::cout << "The array of all distance between nodes is empty " << std::endl;
    }
}
//=================================================================================================//
void BaseInnerRelationInFVM2::resetNeighborhoodCurrentSize()
{
    parallel_for(
        IndexRange(0, base_particles_.TotalRealParticles()),
        [&](const IndexRange &r)
        {
            for (size_t num = r.begin(); num != r.end(); ++num)
            {
                inner_configuration_[num].current_size_ = 0;
            }
        },
        ap);
}

//=================================================================================================//
InnerRelationInFVM2::InnerRelationInFVM2(RealBody &real_body, ANSYSMesh2 &ansys_mesh)
    : BaseInnerRelationInFVM2(real_body, ansys_mesh), get_inner_neighbor_(&real_body){};
//=================================================================================================//
template <typename GetParticleIndex, typename GetNeighborRelation>
void InnerRelationInFVM2::searchNeighborsByParticles(size_t total_particles, BaseParticles &source_particles,
                                                    ParticleConfiguration &particle_configuration,
                                                    GetParticleIndex &get_particle_index, GetNeighborRelation &get_neighbor_relation)
{
    parallel_for(
        IndexRange(0, base_particles_.TotalRealParticles()),
        [&](const IndexRange &r)
        {
            for (size_t num = r.begin(); num != r.end(); ++num)
            {
                size_t index_i = get_particle_index(num);
                Vecd &particle_position = pos_[index_i];

                Neighborhood &neighborhood = particle_configuration[index_i];
                for (std::vector<std::vector<long unsigned int>>::size_type neighbor = 0; neighbor != mesh_topology_[index_i].size(); ++neighbor)
                {
                    size_t index_j = mesh_topology_[index_i][neighbor][0] - 1;
                    size_t boundary_type = mesh_topology_[index_i][neighbor][1];
                    size_t interface_node1_index = mesh_topology_[index_i][neighbor][2];
                    size_t interface_node2_index = mesh_topology_[index_i][neighbor][3];
                    size_t interface_node3_index = mesh_topology_[index_i][neighbor][4];
                    Vecd node1_position = Vecd(node_coordinates_[interface_node1_index][0], node_coordinates_[interface_node1_index][1], node_coordinates_[interface_node1_index][2]);
                    Vecd node2_position = Vecd(node_coordinates_[interface_node2_index][0], node_coordinates_[interface_node2_index][1], node_coordinates_[interface_node2_index][2]);
                    Vecd node3_position = Vecd(node_coordinates_[interface_node3_index][0], node_coordinates_[interface_node3_index][1], node_coordinates_[interface_node3_index][2]);
                    Vecd interface_area_vector1 = node1_position - node2_position;
                    Vecd interface_area_vector2 = node1_position - node3_position;
                    Vecd normal_vector = interface_area_vector1.cross(interface_area_vector2);
                    Real magnitude = normal_vector.norm();
                    Vecd normalized_normal_vector = normal_vector / magnitude;
                    Vecd node1_to_center_direction = particle_position - node1_position;
                    if (node1_to_center_direction.dot(normalized_normal_vector) < 0)
                    {
                        normalized_normal_vector = -normalized_normal_vector;
                    };
                    Real r_ij = 0; // we need r_ij to calculate the viscous force
                    // boundary_type == 2 means both of them are inside of fluid
                    if (boundary_type == 2)
                    {
                        r_ij = (particle_position - pos_[index_j]).dot(normalized_normal_vector);
                    }
                    // this refer to the different types of wall boundary conditions
                    if ((boundary_type == 3) | (boundary_type == 4) | (boundary_type == 5) | (boundary_type == 7) | (boundary_type == 9) |
                        (boundary_type == 10) | (boundary_type == 36))
                    {
                        r_ij = node1_to_center_direction.dot(normalized_normal_vector) * 2.0;
                    }
                    Real dW_ij = (-0.5 * magnitude) / (2.0 * Vol_[index_i] * Vol_[index_j]);
                    get_neighbor_relation(neighborhood, r_ij, dW_ij, normalized_normal_vector, index_j);
                }
            }
        },
        ap);
}
//=================================================================================================//
void InnerRelationInFVM2::updateConfiguration()
{
    resetNeighborhoodCurrentSize();
    searchNeighborsByParticles(base_particles_.TotalRealParticles(),
                               base_particles_, inner_configuration_,
                               get_particle_index_, get_inner_neighbor_);
}

} // namespace SPH

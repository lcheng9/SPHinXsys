#pragma once

#include <fstream>
#include <vector>

namespace unvMesh2412
{

typedef std::vector<unsigned int> TNodeLabels; // Nodal connectivities

struct UNVRecordData
{
    UNVRecordData();

    unsigned long long label;
    int fe_descriptor_id = -1; // FE descriptor id
    int phys_prop_tab_num;     // physical property table number
    int mat_prop_tab_num;      // material property table number
    int color;                 // color
    TNodeLabels node_labels;   // node labels defining element

    // FOR BEAM ELEMENTS ONLY
    int beam_orientation; // beam orientation node number
    int beam_fore_end;    // beam fore-end cross section number
    int beam_aft_end;     // beam  aft-end cross section number
};

typedef std::vector<UNVRecordData> TDataSet;

void read_stream(std::ifstream &in_stream, TDataSet &theDataSet);

void write_stream(std::ofstream &out_stream, const TDataSet &theDataSet);
void write_stream(std::ostream &out_stream, const TDataSet &theDataSet);

bool is_beam_element(int theFeDescriptorId);
bool is_face_element(int theFeDescriptorId);
bool is_body_element(int theFeDescriptorId);

}; // namespace unvMesh2412

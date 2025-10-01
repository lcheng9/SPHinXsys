#include "TransferUNV164DS.h"
#include "TransferUNVTools.h"

#include <fstream>
#include <cstdio>
#include <cmath>

using namespace std;
using namespace unvMesh;
using namespace unvMesh164;

static string _label_dataset = "164";

void unvMesh164::read_stream(std::ifstream& in_stream, UNVRecordData& theUnitsRecord )
{
  if(!in_stream.good())
    EXCEPTION(runtime_error,"Error: bad input UNV file.");

  if(!beginning_of_dataset(in_stream,_label_dataset))
    return;

  string num_buf;
  char line[theMaxLineLen] = "";

  in_stream >> theUnitsRecord.units_code;
  in_stream.readsome( line, 20 );
  theUnitsRecord.units_description = line;
  in_stream >> theUnitsRecord.temp_mode;

  for ( int i = 0; i < 4; i++ )
  {
    in_stream >> num_buf;
    theUnitsRecord.factors[i] = D_to_e(num_buf);
  }
}

void unvMesh164::write_stream(std::ofstream& out_stream, const UNVRecordData& theDataSet)
{
  if(!out_stream.good())
    EXCEPTION(runtime_error,"Error: bad UNV output file for header data.");
  
  out_stream << "    -1" << endl;
  out_stream << "   " << _label_dataset << endl;
  out_stream <<"         1  SI: Meter          2" << endl;
  //out_stream << "         " << theDataSet.units_code << "  " << theDataSet.units_description << "          " << theDataSet.temp_mode << endl;
  out_stream << "    1.0000E+0    1.0000E+0    1.0000E+0" << endl;
  out_stream << "    2.7315E+2" << endl;

  out_stream << "    -1" << endl;
}

void unvMesh164::write_stream(std::ostream& out_stream, const UNVRecordData& theDataSet)
{
    if (!out_stream.good())
        EXCEPTION(runtime_error, "Error: bad UNV output file for header data.");

    out_stream << "    -1" << endl;
    out_stream << "   " << _label_dataset << endl;
    out_stream <<"         1  SI: Meter          2" << endl;
    //out_stream << "         " << theDataSet.units_code << "  " << theDataSet.units_description << "          " << theDataSet.temp_mode << endl;
    out_stream << "    1.0000E+0    1.0000E+0    1.0000E+0" << endl;
    out_stream << "    2.7315E+2" << endl;

    out_stream << "    -1" << endl;
}

unvMesh164::UNVRecordData::UNVRecordData()
{
  units_code        = 1;
  units_description = "SI: Meter";
  temp_mode         = 2;
  factors[0]        = 1.0;
  factors[1]        = 1.0;
  factors[2]        = 1.0;
  factors[3]        = 273.15;
}

void UNVRecordData::set_length_unit(const std::string& s)
{
    if (s == "M") {
        units_code = 1;
        units_description = "SI: Meter";
    }
    else if (s == "FT") {
        units_code = 2;
        units_description = "BG: Foot";
    }
    else if (s == "MM") {
        units_code = 5;
        units_description = "MM: mm";
    }
    else if (s == "CM") {
        units_code = 6;
        units_description = "CM: cm";
    }
    else if (s == "IN") {
        units_code = 7;
        units_description = "IN: Inch";
    }
    else
    {
        assert(0);
        units_code = 0;
        units_description = "Unknown Unit";
    }
}

#include "TransferUNV2420DS.h"
#include "TransferUNVTools.h"

#include <fstream>
#include <cstdio>
#include <cmath>

using namespace std;
using namespace unvMesh;
using namespace unvMesh2420;

static string _label_dataset = "2420";

void unvMesh2420::read_stream(std::ifstream& in_stream,
                   std::string&   part_name, // can re-store a mesh name
                   TDataSet&      theDataSet)
{
  if(!in_stream.good())
    EXCEPTION(runtime_error,"-=*Error: bad input file");

  /*
   * adjust the \p istream to our
   * position
   */
  if(!beginning_of_dataset(in_stream,_label_dataset))
    return;

  string num_buf;
  int part_uid;

  in_stream >> part_uid; // Record 1
  part_name = read_line( in_stream );  // Record 2

  while ( !in_stream.eof() )
  {
    UNVRecordData aRec;

    // Record 3
    in_stream >> aRec.coord_sys_label;
    if ( aRec.coord_sys_label == -1 ) // end of dataset is reached
      break;
    in_stream >> aRec.coord_sys_type;
    in_stream >> aRec.coord_sys_color;

    aRec.coord_sys_name = read_line( in_stream ); // Record 4

    // Records 5-8: rows of Transformation Matrix
    for ( int row = 0; row < 4; ++row )
      for ( int i = 0; i < 3; i++ )
      {
        in_stream >> num_buf;
        aRec.matrix[row][i] = D_to_e(num_buf);
      }
    // Store a CS data only if it requires conversion into the global Cartesian CS
    if ( aRec.coord_sys_type != 0 || !aRec.is_identity_tensor() ) // 0 - Cartesian CS
      theDataSet.push_back( aRec );
  }
}


void unvMesh2420::write_stream(std::ofstream&     out_stream,
                    const std::string& part_name)
//                    const TDataSet& theDataSet)
{
  if(!out_stream.good())
    EXCEPTION(runtime_error,"Error: bad UNV output file.");
  
  out_stream<<"    -1"  << endl;
  out_stream<<"  "<<_label_dataset << endl;

  out_stream<<"         1"                     << endl; // R1: Part UID
  if ( part_name.empty() )
    out_stream<<"Welsim_Mesh"                   << endl; // R2: Part Name
  else
    out_stream<< part_name                     << endl;
  out_stream<<"         1         0         0" << endl; // R3: Label, Type, Color

  out_stream<<"Global coordinate system" << endl; // R4: Name
  out_stream<<"    1.0000E+0    0.0000E+0    0.0000E+0" << endl;
  out_stream<<"    0.0000E+0    1.0000E+0    0.0000E+0" << endl;
  out_stream<<"    0.0000E+0    0.0000E+0    1.0000E+0" << endl;
  out_stream<<"    0.0000E+0    0.0000E+0    0.0000E+0" << endl;

  out_stream<<"    -1"  << endl;
}


void unvMesh2420::write_stream(std::ostream& out_stream, const std::string& part_name)
{
    if (!out_stream.good())
        EXCEPTION(runtime_error, "Error: bad UNV output file.");

    out_stream << "    -1" << endl;
    out_stream << "  " << _label_dataset << endl;

    out_stream << "         1" << endl; // R1: Part UID
    if (part_name.empty())
        out_stream << "Welsim_Mesh" << endl; // R2: Part Name
    else
        out_stream << part_name << endl;
    out_stream << "         1         0         0" << endl; // R3: Label, Type, Color

    out_stream << "Global coordinate system" << endl; // R4: Name
    out_stream << "    1.0000E+0    0.0000E+0    0.0000E+0" << endl;
    out_stream << "    0.0000E+0    1.0000E+0    0.0000E+0" << endl;
    out_stream << "    0.0000E+0    0.0000E+0    1.0000E+0" << endl;
    out_stream << "    0.0000E+0    0.0000E+0    0.0000E+0" << endl;

    out_stream << "    -1" << endl;
}

bool unvMesh2420::UNVRecordData::is_identity_tensor() const
{
  bool isIdentity = true;
  for ( int row = 0; row < 4 && isIdentity; ++row )
    for ( int i = 0; i < 3; i++ )
    {
      if ( matrix[row][i] != ( row==i ? 1. : 0. ))
      {
        isIdentity = false;
        break;
      }
    }
  return isIdentity;
}

void unvMesh2420::UNVRecordData::modify_by_tensor( double* c ) const
{
  const double x = matrix[0][0] * c[0] + matrix[0][1] * c[1] + matrix[0][2] * c[2];
  const double y = matrix[1][0] * c[0] + matrix[1][1] * c[1] + matrix[1][2] * c[2];
  const double z = matrix[2][0] * c[0] + matrix[2][1] * c[1] + matrix[2][2] * c[2];
  c[0] = x + matrix[3][0];
  c[1] = y + matrix[3][1];
  c[2] = z + matrix[3][2];
}

void unvMesh2420::UNVRecordData::convert_by_CCS( double* coords )
{
  const double x = coords[0] * cos( coords[1] );
  const double y = coords[0] * sin( coords[1] );
  coords[0] = x;
  coords[1] = y;
}

void unvMesh2420::UNVRecordData::convert_by_SCS  ( double* coords )
{
  const double sin2 = sin( coords[2] );
  const double x = coords[0] * cos( coords[1] ) * sin2;
  const double y = coords[0] * sin( coords[1] ) * sin2;
  const double z = coords[0] * cos( coords[2] );
  coords[0] = x;
  coords[1] = y;
  coords[2] = z;
}

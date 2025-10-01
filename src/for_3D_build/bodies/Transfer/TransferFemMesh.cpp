#include "TransferFemMesh.h"

#include "../WMESHTools/WMESHTextLabel.h"

#include "../WELSIMLocalTrace/TraceTools.h"

using namespace std;

TransferFemMesh::TransferFemMesh():
  myFile(""),
  myFolder(""),
  myMeshId(-1),
  myStatus( DRS_OK ),
  myPrimaryTopo(wsMsh::kMeshTopoUnknown), 
  myLengthUnit("")
{}


void TransferFemMesh::set_mesh_id(int theMeshId)
{
  myMeshId = theMeshId;
}

void TransferFemMesh::set_mesh_name(const std::string& theMeshName)
{
  myMeshName = theMeshName;
}
std::string TransferFemMesh::get_mesh_name() const
{
  return myMeshName;
}

void TransferFemMesh::set_file_name(const std::string& theFileName)
{
  myFile = theFileName;
}
void TransferFemMesh::set_folder(const QString& path)
{
    myFolder = path;
}
void TransferFemMesh::set_length_unit(const std::string& s)
{
    myLengthUnit = s;
}


//================================================================================
/*!
 * \brief Stores an error message
 *
 * We consider an error fatal if none mesh can be read
 */
//================================================================================

TransferFemMesh::Status TransferFemMesh::add_mesh_message(const std::string& msg,
                                            const bool         isFatal/*=false*/)
{
  if ( isFatal )
    myErrorMessages.clear(); // warnings are useless if a fatal error encounters

  myErrorMessages.push_back( msg );

  MESSAGE(msg);
#ifdef _DEBUG
  cout << msg << endl;
#endif
  return ( myStatus = isFatal ? DRS_FAIL : DRS_WARN_SKIP_ELEM );
}

//================================================================================
/*!
 * \brief Return a structure containing description of errors
 */
//================================================================================

WMESHDoMeshingErrorPtr TransferFemMesh::get_mesh_error()
{
  WMESHTextLabel msg;
  for ( size_t i = 0; i < myErrorMessages.size(); ++i )
  {
    msg << myErrorMessages[i];
    if ( i+1 < myErrorMessages.size() )
      msg << "\n";
  }
  return WMESHDoMeshingError::New( myStatus == DRS_OK ? int(COMPERR_OK) : int(myStatus), msg );
}

#ifndef _TRANSFER_FEM_MESH_H_
#define _TRANSFER_FEM_MESH_H_

//#include <Mesher/MeshRep/mesh_enums.h>

//#include "../WMESHTools/WMESHDoMeshingError.h"

#include <string>
#include <vector>
//#include <QString>

#include "../MeshDB_PreCompiled.h"


class MESHTRANSFER_EXPORT TransferFemMesh
{
 public:
  TransferFemMesh();
  virtual ~TransferFemMesh(){}

  enum Status {
    DRS_OK,
    DRS_EMPTY,          // a file contains no mesh with the given name
    DRS_WARN_RENUMBER,  // a file has overlapped ranges of element numbers,
                        // so the numbers from the file are ignored
    DRS_WARN_SKIP_ELEM, // some elements were skipped due to incorrect file data
    DRS_WARN_DESCENDING, // some elements were skipped due to descending connectivity
    DRS_FAIL            // general failure (exception etc.)
  };

  void                set_mesh_id(int theMeshId);
  void                set_length_unit(const std::string& s);
  virtual void        set_file_name(const std::string& theFileName);
  virtual void        set_mesh_name(const std::string& theMeshName);
  virtual std::string get_mesh_name() const;
  void set_folder(const QString& path);

  virtual void        set_option_value(const std::string& optionName,
                                const std::string& optionValue) {}

  virtual Status Execute(bool binary) = 0;

  virtual WMESHDoMeshingErrorPtr get_mesh_error();

 protected:
  std::string myFile;
  QString myFolder;
  std::string myMeshName;
  int         myMeshId;
  wsMsh::MeshTopoType myPrimaryTopo;
  std::string myLengthUnit;

  Status add_mesh_message(const std::string& msg, const bool isFatal=false);
  std::vector< std::string > myErrorMessages;
  Status                     myStatus;
};

#endif

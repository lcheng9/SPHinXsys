#ifndef _TRANSFER_FILE_H
#define _TRANSFER_FILE_H

#include "TransferFemMesh.h"

#include <string>
#include "../MeshDB_PreCompiled.h"

class WMESHDataProject;

class MESHTRANSFER_EXPORT TransferFile
{
 public:
  TransferFile();
  virtual ~TransferFile(){}

  virtual void execute() = 0;
  void set_file_name(const std::string& theFileName);
  void set_meshdata_file(WMESHDataProject *theDocument);

 protected:
  WMESHDataProject * myDocument;
  std::string myFile;

};


#endif

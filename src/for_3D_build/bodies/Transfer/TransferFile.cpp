#include "TransferFile.h"

TransferFile::TransferFile():
  myDocument(NULL)
{}


void TransferFile::set_file_name(const std::string& theFileName)
{
  myFile = theFileName;
}


void TransferFile::set_meshdata_file(WMESHDataProject * theDocument)
{
  myDocument = theDocument;
}

#pragma once

#include <map>
#include <vector>
#include <fstream>      
#include <string>       


namespace unvMesh2417{

  typedef std::vector<unsigned long> TListOfId; // Nodal connectivitiesList of Id

  struct UNVRecordData{
    std::string GroupName;
    TListOfId NodeList;
    TListOfId ElementList;
  };

  typedef int TGroupId; // type of element label
  typedef std::map<TGroupId, UNVRecordData> TDataSet;

  void read_stream(std::ifstream& in_stream, TDataSet& theDataSet);
  void read_stream_group(const std::string& myGroupLabel, std::ifstream& in_stream, TDataSet& theDataSet);

  void write_stream(std::ofstream& out_stream, const TDataSet& theDataSet);
  void write_stream(std::ostream& out_stream, const TDataSet& theDataSet);
};


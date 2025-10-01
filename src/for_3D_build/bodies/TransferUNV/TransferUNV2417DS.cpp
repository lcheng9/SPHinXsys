#include "TransferUNV2417DS.h"
#include "TransferUNVTools.h"

#include <fstream>      
#include <iomanip>

using namespace std;
using namespace unvMesh;
using namespace unvMesh2417;

static string _group_labels[] = {"2417", "2429", "2430", "2432",
                                 "2435", "2452", "2467", "2477"};
#define NBGROUP 8

static string _label_dataset = "2467";

void unvMesh2417::read_stream(std::ifstream& in_stream, TDataSet& theDataSet)
{
  if(!in_stream.good())
    EXCEPTION(runtime_error,"-=*Error: bad input file");

  std::string olds, news;
  
  while(true){
    in_stream >> olds >> news;
    /*
     * a "-1" followed by a number means the beginning of a dataset
     * stop combing at the end of the file
     */
    while( ((olds != "-1") || (news == "-1") ) && !in_stream.eof() ){     
      olds = news;
      in_stream >> news;
    }
    if(in_stream.eof())
      return;
    for (int i = 0; i < NBGROUP; i++) {
      if (news == _group_labels[i]) {
        read_stream_group(news, in_stream, theDataSet);
      }
    }
  }
}



void unvMesh2417::read_stream_group(const std::string& myGroupLabel, std::ifstream& in_stream, TDataSet& theDataSet)
{
  TGroupId aId;
  for(; !in_stream.eof();){
    in_stream >> aId ;
    if(aId == -1){
      // end of dataset is reached
      break;
    }

    int n_nodes;
    UNVRecordData aRec;
    int aTmp;
    in_stream>>aTmp; // miss not necessary values
    in_stream>>aTmp;
    in_stream>>aTmp;
    in_stream>>aTmp;
    in_stream>>aTmp;
    in_stream>>aTmp;
    in_stream>>n_nodes;

    std::getline(in_stream, aRec.GroupName, '\n'); // Finalise previous reading
    std::getline(in_stream, aRec.GroupName, '\n');

    int aElType;
    int aElId;
    int aNum;
    for(int j=0; j < n_nodes; j++){
      in_stream>>aElType;
      in_stream>>aElId;
      if ((myGroupLabel.compare("2435") == 0) ||
          (myGroupLabel.compare("2452") == 0) ||
          (myGroupLabel.compare("2467") == 0) ||
          (myGroupLabel.compare("2477") == 0)) {
        in_stream>>aTmp;
        in_stream>>aTmp;
      }
      switch (aElType) {
      case 7: // Nodes
        aNum = aRec.NodeList.size();
        aRec.NodeList.resize(aNum + 1);
        aRec.NodeList[aNum] = aElId;
        break;
      case 8: // Elements
        aNum = aRec.ElementList.size();
        aRec.ElementList.resize(aNum + 1);
        aRec.ElementList[aNum] = aElId;
        break;
      }
    }
    theDataSet.insert(TDataSet::value_type(aId,aRec));
  }

}


//Record 1: FORMAT(8I10)
//Field 1 --group number
//Field 2 --active constraint set no. for group
//Field 3 --active restraint set no. for group
//Field 4 --active load set no. for group
//Field 5 --active dof set no. for group
//Field 6 --active temperature set no. for group
//Field 7 --active contact set no. for group
//Field 8 --number of entities in group
void unvMesh2417::write_stream(std::ofstream& out_stream, const TDataSet& theDataSet)
{
  if(!out_stream.good())
    EXCEPTION(runtime_error,"-=*Error: bad output file");
  
  /*
   * Write beginning of dataset
   */
  out_stream<<"    -1\n";
  out_stream<<"  "<<_label_dataset<<"\n";

  TDataSet::const_iterator anIter = theDataSet.begin();
  for(; anIter != theDataSet.end(); anIter++){
    const TGroupId& aLabel = anIter->first;
    const UNVRecordData& aRec = anIter->second;
    int aNbNodes = aRec.NodeList.size();
    int aNbElements = aRec.ElementList.size();
    int aNbRecords = aNbNodes + aNbElements;

    out_stream<<std::setw(6)<<aLabel;  /* group ID */
    out_stream<<std::setw(6)<<0;  
    out_stream<<std::setw(6)<<0;
    out_stream<<std::setw(6)<<0;
    out_stream<<std::setw(6)<<0;
    out_stream<<std::setw(6)<<0;
    out_stream<<std::setw(6)<<0;
    out_stream<<std::setw(10)<<aNbRecords<<std::endl; 

    // GroupName must be one single name
    out_stream<<aRec.GroupName<<std::endl;
    for (int i = 0; i < aRec.GroupName.size(); ++i) {
        if (isspace(aRec.GroupName[i])) {
            assert(0);
            break;
        }
    }


    int aRow = 0;
    int i;
    for (i = 0; i < aNbNodes; i++) {
      if (aRow == 2) {
        out_stream<<std::endl; 
        aRow = 0;
      }
      out_stream<<std::setw(6)<<7;
      out_stream<<std::setw(10)<<aRec.NodeList[i];
      out_stream<<std::setw(6)<<0;
      out_stream<<std::setw(6)<<0;
      aRow++;
    }
    for (i = 0; i < aNbElements; i++) {
      if (aRow == 2) {
        out_stream<<std::endl; 
        aRow = 0;
      }
      out_stream<<std::setw(6)<<8;
      out_stream<<std::setw(10)<<aRec.ElementList[i];
      out_stream<<std::setw(6)<<0;
      out_stream<<std::setw(6)<<0;
      aRow++;
    }
    out_stream<<std::endl; 
  }

  /*
   * Write end of dataset
   */
  out_stream<<"    -1\n";
}

void unvMesh2417::write_stream(std::ostream& out_stream, const TDataSet& theDataSet)
{
    if (!out_stream.good())
        EXCEPTION(runtime_error, "-=*Error: bad output file");

    /*
    * Write beginning of dataset
    */
    out_stream << "    -1\n";
    out_stream << "  " << _label_dataset << "\n";

    TDataSet::const_iterator anIter = theDataSet.begin();
    for (; anIter != theDataSet.end(); anIter++) {
        const TGroupId& aLabel = anIter->first;
        const UNVRecordData& aRec = anIter->second;
        int aNbNodes = aRec.NodeList.size();
        int aNbElements = aRec.ElementList.size();
        int aNbRecords = aNbNodes + aNbElements;

        out_stream << std::setw(9) << aLabel;  /* group ID */
        out_stream << std::setw(6) << 0;
        out_stream << std::setw(6) << 0;
        out_stream << std::setw(6) << 0;
        out_stream << std::setw(6) << 0;
        out_stream << std::setw(6) << 0;
        out_stream << std::setw(6) << 0;
        out_stream << std::setw(10) << aNbRecords << std::endl;

        // GroupName must be one single name
        out_stream << aRec.GroupName << std::endl;
        for (int i = 0; i < aRec.GroupName.size(); ++i) {
            if (isspace(aRec.GroupName[i])) {
                assert(0);
                break;
            }
        }

        int aRow = 0;
        int i;
        for (i = 0; i < aNbNodes; i++) {
            if (aRow == 2) {
                out_stream << std::endl;
                aRow = 0;
            }
            out_stream << std::setw(6) << 7;
            out_stream << std::setw(10) << aRec.NodeList[i];
            out_stream << std::setw(6) << 0;
            out_stream << std::setw(6) << 0;
            aRow++;
        }
        for (i = 0; i < aNbElements; i++) {
            if (aRow == 2) {
                out_stream << std::endl;
                aRow = 0;
            }
            out_stream << std::setw(6) << 8;
            out_stream << std::setw(10) << aRec.ElementList[i];
            out_stream << std::setw(6) << 0;
            out_stream << std::setw(6) << 0;
            aRow++;
        }
        out_stream << std::endl;
    }

    /*
    * Write end of dataset
    */
    out_stream << "    -1\n";
}

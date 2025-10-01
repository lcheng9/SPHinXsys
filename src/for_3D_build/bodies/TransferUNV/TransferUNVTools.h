#pragma once

#include <iostream>     
#include <sstream>      
#include <fstream>
#include <string>
#include <stdexcept>
#include <cassert>
#include <cstdlib>

#include "WMDataEnums.h"

namespace unvMesh{
  using namespace std;

  const size_t theMaxLineLen = 82; // 80 for text + 2 for "\r\n"

  class UNVCounter{
    static int myCounter;
  public:
    UNVCounter();
    ~UNVCounter();

    static string GetUNVCountNumber();
  };

  /**
   * @returns \p false when error occured, \p true otherwise.
   * Adjusts the \p in_stream to the beginning of the
   * dataset \p ds_name.
   */
  inline bool beginning_of_dataset(std::istream& in_file, const std::string& ds_name)
  {
    assert (in_file.good());
    assert (!ds_name.empty());
    
    std::string olds, news;
    
    in_file.seekg(0);
    while(true){
      in_file >> olds >> news;
      /*
       * a "-1" followed by a number means the beginning of a dataset
       * stop combing at the end of the file
       */
      while( ((olds != "-1") || (news == "-1") ) && !in_file.eof() ){     
        olds = news;
        in_file >> news;
      }
      if(in_file.eof())
      {
        in_file.clear();
        return false;
      }
      if (news == ds_name)
        return true;
    }
    // should never end up here
    return false;
  }

  /**
   * Method for converting exponential notation
   * from "D" to "e", for example
   * \p 3.141592654D+00 \p --> \p 3.141592654e+00
   * in order to make it readable for C++.
   */
  inline double D_to_e(std::string& number)
  {
    /* find "D" in string, start looking at 
     * 6th element, to improve speed.
     * We dont expect a "D" earlier
     */
    const int position = number.find("D",6);
    if(position != std::string::npos){
      number.replace(position, 1, "e"); 
    }
    return atof (number.c_str());
  }
  
  /**
   * @returns \p false when file is incorrect, \p true otherwise.
   * Check file with name \p theFileName for correct terminate
   * string, i.e. the next to the last line is equal to "    -1",
   */
  inline bool check_file(const std::string theFileName)
  {
    std::ifstream in_stream(theFileName.c_str());
    if (!in_stream)
      return false;
    std::string olds, news;
    while (!in_stream.eof()){
      olds = news;
      std::getline(in_stream, news, '\n');
    }
    return (olds == "    -1");
  }

  /*!
   * \brief reads a whole line
   *  \param in_stream - source stream
   *  \param next - if true, first reads the current line up to the end
   *  which is necessary after reading using >> operator
   *  \retval std::string - the result line
   */
  inline std::string read_line(std::ifstream& in_stream, const bool next=true)
  {
    std::string resLine;
    std::getline( in_stream, resLine );
    if ( next )
      std::getline( in_stream, resLine );

    if ( resLine.size() > 0 && resLine[ resLine.size()-1 ] == '\r' )
      resLine.resize( resLine.size()-1 );
    return resLine;
  }

  void writeFoamFile(const std::string& sClass, const std::string& sObject, std::ofstream& out_stream);
  int getElmerElemType(WMDataElemType elemType, int nNodes);
};


#ifndef MESSAGE

#define MESSAGE(msg) std::cout<<__FILE__<<"["<<__LINE__<<"]::"<<msg<<endl;

#define BEGMSG(msg) std::cout<<UNV::UNVCounter::GetUNVCountNumber()<<msg

#define ADDMSG(msg) std::cout<<msg

#endif


#ifndef EXCEPTION

#define EXCEPTION(TYPE, MSG) {\
  std::ostringstream aStream;\
  aStream<<__FILE__<<"["<<__LINE__<<"]::"<<MSG;\
  throw TYPE(aStream.str());\
}

#endif

#ifndef TOOL_H
#define TOOL_H
#include "TH1F.h"
#include <string>
#include <map>
class tool{
  private:

  public:
    tool(){
      frameSide ="NA";
      frameName ="";
      frameTliquid = -999.;
    }
    ~tool(){}

    std::string frameSide;
    std::string frameName;
    float frameTliquid;
    bool read_config(std::string sconfig="config");
    bool find_pipes(std::string sfile_input, std::string sfolder_input = "../roo", std::string sname_output = "fitPipe", std::string sfile_output = "a1_pipes.root");

};

#endif

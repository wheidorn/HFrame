#ifndef TOOL_H
#define TOOL_H
#include "HFrame.h"
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
    bool read_frame(std::string sfile, std::string sfolder = "roo");

};

#endif

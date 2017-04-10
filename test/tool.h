#ifndef TOOL_H
#define TOOL_H
#include "HFrame.h"
#include "TH1F.h"
#include <string>
#include <map>
class tool{
  private:
    std::map<std::string, TH1F*> pipeMap;

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
    HFrame * read_frame(std::string sfile, std::string sfolder = "roo");
    bool read_pipe(std::string sprefix, std::string sfile, std::string sfolder = "roo");

    TH1F* getPipe(std::string hname){ return pipeMap[hname];}
};

#endif

#include <fstream>
#include <sstream>
#include <string>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include "tool.h"
#include "HPlot.h"
#include "HTools.h"
using namespace std;

template <typename T>
string NumberToString ( T Number )
{
   ostringstream ss; 
   ss << Number;
   return ss.str();
}

bool tool::read_config(string sconfig){

  ifstream inconf(sconfig.c_str());
  string line;
  while (getline(inconf, line))
  {
    size_t k=line.find(":");
    if (line.size() <=0 || k ==string::npos) continue;
    if (line[0]=='#' || line[0]=='*') continue;

    if      (line.substr(0,k) =="Name" ) frameName = line.substr(k+1);
    else if (line.substr(0,k) =="Side") frameSide = line.substr(k+1);
    else if (line.substr(0,k) =="Tliquid") frameTliquid = stoi( string(line.substr(k+1)) );
    else cout<<"key "<<line.substr(0,k)<<"not known. Use: Name, Side, Tliquid. "<<endl;
  }
  if (frameTliquid <-998.) return false;
  if (frameSide != "J" && frameSide != "L") return false;
  return true;
}

bool tool::find_pipes(string sfile_input, string sfolder_input, string sname_output, string sfile_output){

  bool negative = false;
  if(frameTliquid < 10.) negative = true;
  else if(frameTliquid <29.) return false;
  // read hists
  string m_sRawFile = sfolder_input + "/Pipes/" + sfile_input + ".root";
  TFile * ff0 = new TFile(m_sRawFile.c_str(), "read");
  HTools * pHT = new HTools();
  vector<string> snames { "fitPipe_temperature_top", "fitPipe_temperature_bot", "fitPipe_temperature"};
  bool ok = true;
  for(auto aname : snames){
    TH1F* h_pipe = (TH1F *) ff0 -> Get(aname.c_str());
    ok = ok &  pHT -> getPeaks(h_pipe, sname_output, negative);
  }

  delete pHT;
  return ok;
}


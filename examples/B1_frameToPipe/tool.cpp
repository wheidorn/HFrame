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

// http://www.cplusplus.com/forum/articles/9645/
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

  // read hists
  string m_sRawFile = sfolder_input + "/combined/" + sfile_input + ".root";
  TFile * ff0 = new TFile(m_sRawFile.c_str(), "read");
  string stitle = frameSide + Form(" side, T_{liquid} = %3.0f #circ C", frameTliquid);

  // get the stave pipe 2D frame
  TH2F* h_pipe = (TH2F *) ff0 -> Get("T_stave_pipe");
  // add error on to h_pipe, should be done at production level.
  for(int ix=1; ix<=h_pipe->GetNbinsX(); ix++){
    for(int iy=1; iy<=h_pipe->GetNbinsY(); iy++){
      h_pipe -> SetBinError(ix, iy, h_pipe -> GetBinContent(ix, iy) * 0.005 );
    }
  }
  if(h_pipe){
    h_pipe -> SetName("pipe"); // output png name
    h_pipe -> GetXaxis() -> SetTitle( "X direction (cm)" );
    h_pipe -> GetYaxis() -> SetTitle( "Y direction (cm)" );
    h_pipe -> SetTitle( stitle.c_str() );
  }

  HTools * pHT = new HTools();
  bool ok = pHT -> fitFramePipesGaus(h_pipe, sname_output, frameTliquid, frameSide, sfile_output);

  delete pHT;
  return ok;
}


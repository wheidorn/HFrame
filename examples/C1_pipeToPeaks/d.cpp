#include <fstream>
#include <sstream>
#include <string>
#include <iostream>
#include <stdlib.h>
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include <TMinuit.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TLegend.h>
#include <TLatex.h>
#include "TFile.h"
#include "TTree.h"
#include "TGaxis.h"

#include "tool.h"

using namespace std;
int main(int argc, char* argv[]){


  // read the file name
  //
  std::string sfolder = "../../roo/";
  std::string sfile = "";
  for (int i1 = 1; i1<argc; ++i1) {
    if (strcmp(argv[i1],"-name")==0) { i1++; sfile = string(argv[i1]); continue;}
    if (argv[i1][0] =='-'){cout<<"arg: "<<argv[i1]<<" not found "<<endl; continue;}
  }
  if (sfile==""){return 0;}
  string sSide = "L";
  if(sfile.find("_sJ") != string::npos || sfile.find("_J") != string::npos) sSide = "J";


  tool * pTool = new tool(); 
  // assuming config file name: "config"
  if(!pTool -> read_config()) return 1;

  if( ! pTool -> find_pipes(sfile, sfolder) ) 
  {
    cout<<"ERROR: read frame failed!"<<endl;
  }
  return 3;
}

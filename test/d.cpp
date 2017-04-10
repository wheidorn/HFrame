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

#include "HPeaks.h"
#include "HFrame.h"
#include "HTools.h"
#include "HPlot.h"
#include "tool.h"

using namespace std;
int main(int argc, char* argv[]){


  // read the file name
  std::string sfolder = "roo";
  std::string sfile = "";
  for (int i1 = 1; i1<argc; ++i1) {
    if (strcmp(argv[i1],"-name")==0) { i1++; sfile = string(argv[i1]); continue;}
    if (argv[i1][0] =='-'){cout<<"arg: "<<argv[i1]<<" not found "<<endl; continue;}
  }
  if (sfile==""){return 0;}

  tool * pTool = new tool(); 
  HTools * ftol = new HTools();
  HPlot  * fplt = new HPlot();
  HPeaks  * fpek = new HPeaks();

  // assuming config file name: "config"
  if(!pTool -> read_config()) return 1;


  string pre_pipename = "fitPipe"; //"fitPipe" + "_temperature_top";

  bool is_pipe_ok = pTool -> read_pipe(pre_pipename, sfile, sfolder);
  // --------------------------------------------
  // if pipe information are not ready, make them
  // --------------------------------------------
  if(!is_pipe_ok){
   HFrame * aframe = pTool -> read_frame(sfile, sfolder);
    if(!aframe) return 2;
    if(pTool -> frameSide=="L") aframe -> setSideL(true);
    else aframe -> setSideL(false);

    float xsize_cm = 128., ysize_cm = 11.6;
    cout<<"result of raw to crop: "<< aframe -> cropRawToStave(xsize_cm, ysize_cm)<<endl;
    float Tliquid = aframe -> getTLiquid();
    cout<<"T liquid: "<< Tliquid <<endl;
    float Tmax = Tliquid + 35, Tmin = Tliquid + 5.;
    if(! aframe -> isNegative()){Tmax = Tliquid; Tmin = Tliquid - 15;}

    fplt -> drawFrame2D(aframe, Tmin, Tmax);
    if(!aframe -> cropStaveToPipe()   ) cout<<"ERROR: crop stave to pipe failed. "<<endl;

    fplt -> setCanvasSize(1000., 300.);
    fplt -> setOffset(1.25, 0.5, 0.6);
    fplt -> drawFramePipe2D(aframe, Tmin, Tmax);
 
    int method = 0;
    if(!ftol -> findFramePipes(aframe, pre_pipename, method))cout<<"ERROR: result of finding pipes failed. "<<endl;

    string m_sCropFile = sfolder + "/cropped/" + sfile + ".root";
    aframe -> writeToRootFile(m_sCropFile);
  }

  double sigma = 2., threshold=0.2; int niters=15; //20
  is_pipe_ok = pTool -> read_pipe(pre_pipename, sfile, sfolder);
  if(!is_pipe_ok){
    cout<<"ERROR: cannot read pipe information! "<<endl;
    return 3;
  }

  TH1F* h_temp_top = pTool -> getPipe(pre_pipename+"_temperature_top");
  TH1F* h_temp_bot = pTool -> getPipe(pre_pipename+"_temperature_bot");
  TH1F* h_temp = pTool -> getPipe(pre_pipename+"_temperature");
  bool negative = false;
  if (pTool -> frameTliquid < 20.) negative = true;
  if (h_temp_top) fpek -> getPeaks(h_temp_top, pTool -> frameName, negative, sigma, threshold, niters);
  if (h_temp_bot) fpek -> getPeaks(h_temp_bot, pTool -> frameName, negative, sigma, threshold, niters);
  if (h_temp)     fpek -> getPeaks(h_temp,     pTool -> frameName, negative, sigma, threshold, niters);

  TH1F* h_pos_top = pTool -> getPipe(pre_pipename+"_peakposition_top");
  TH1F* h_pos_bot = pTool -> getPipe(pre_pipename+"_peakposition_bot");
  TH1F* h_wid_top = pTool -> getPipe(pre_pipename+"_peakwidth_top");
  TH1F* h_wid_bot = pTool -> getPipe(pre_pipename+"_peakwidth_bot");
  TH1F* h_chi_top = pTool -> getPipe(pre_pipename+"_chi2overndf_top");
  TH1F* h_chi_bot = pTool -> getPipe(pre_pipename+"_chi2overndf_bot");

  if( !fplt -> drawTwoPipes(h_pos_top, h_pos_bot) || 
      !fplt -> drawTwoPipes(h_wid_top, h_wid_bot) || 
      !fplt -> drawTwoPipes(h_chi_top, h_chi_bot) ){
    cout<<"ERROR: cannot draw pipe inform plots. "<<endl;
    return 4;
  }

  return 3;
}

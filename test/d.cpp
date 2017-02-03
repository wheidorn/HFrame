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

#include "HFrame.h"
#include "HTools.h"
#include "HPlot.h"

using namespace std;
int main(int argc, char* argv[]){


  HTools * ftol = new HTools();
  HPlot  * fplt = new HPlot();
  int nfile = 1;
  HFrame * fbox = new HFrame();
  fbox -> fillboxHist("m55.root");
  fbox -> crop(42, 215,  590, 221,  41, 270);
  fplt -> drawFrame2D(fbox);

  int nbinx= fbox -> getFrame() -> GetNbinsX();
  int nbiny= fbox -> getFrame() -> GetNbinsY();
  int x_min = 5, x_max = 540;
  int y_low =  11, y_hig = 44;
  bool negative = true;
  ftol -> fitPipeTwoGaus(fbox, "fitPipe", x_min, x_max, y_low, y_hig, negative);

  ftol -> averagePipe(fbox, "avgPipe_top", x_min, x_max, 34, 40, negative); // average pixels around 35 -- 37
  ftol -> averagePipe(fbox, "avgPipe_bot", x_min, x_max, 13, 19, negative); // a few Y pixels along X

  fbox -> writeToRootFile("a1.root");

  return 3;
}

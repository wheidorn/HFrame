#ifndef HTOOLS_H
#define HTOOLS_H

#include <sstream>
#include <string>
#include <iostream>
#include <map>
#include <algorithm>
#include <utility>
#include <vector>
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TMath.h"

class HTools{

  private:
    // providing 2D histogram
    // get 1D in x=ibin, Y=1 to NYBins if isY = true
    TH1F* getHist1D(TH2F* h2, int ibin, bool isY = true);

  public:
    HTools(){}
    ~HTools(){}

    // fitting stave cooling pipe temperature profile
    // at each x in the range [x_min, x_max]
    // get TH1F for [1, nybins]
    // fit two gaus [y1-rangefity, y1+rangefity] and [y2-rangefity, y2+rangefity] 
    // store the information of the fitted results
    bool fitFramePipesGaus(TH2F* hframe, std::string s_out, float Tliquid, std::string s_side = "L", std::string sfile = "a1_pipes.root");
    bool getPeaks(TH1F* hp, std::string s_out, bool negative, double sigma = 2., double threshold=0.15, int niters=15); //20

};
#endif

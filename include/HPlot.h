#ifndef HPLOT_H
#define HPLOT_H

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

class HPlot{

  /* class HPlot aim to
   *
   * -- draw nice plots for the stave T profile by resetting the margins, etc.
   * -- add necessary information into the plot.
   */
  
  private:
    TCanvas * fCurCanvas2D;
    TCanvas * fCurCanvas1D;
    float fSizeX;
    float fSizeY;
    float fMarginL;
    float fMarginR;
    float fMarginT;
    float fMarginB;
    float fOffsetX;
    float fOffsetY;
    float fOffsetZ;
    float fLabelSizeX;
    float fLabelSizeY;
    float fLabelSizeZ;
 
  public:
    HPlot(){ 
      fCurCanvas1D = NULL; fCurCanvas2D = NULL;
      fSizeX = 1000.;
      fSizeY = 450.;
      fMarginL = 0.06;
      fMarginR = 0.11;
      fMarginT = 0.07;
      fMarginB = 0.14;
      fOffsetX = 1.25;
      fOffsetY = 0.65;
      fOffsetZ = 0.75;
      fLabelSizeX = 1.1;
      fLabelSizeY = 1.1;
      fLabelSizeZ = 1.1;
    }

    ~HPlot(){ 
      if(!fCurCanvas1D) delete fCurCanvas1D; 
      if(!fCurCanvas2D) delete fCurCanvas2D; 
    }
    void setCanvasSize(float sizex, float sizey){
      fSizeX = sizex;
      fSizeY = sizey;
    }
    void setMargins(float left, float right, float bottom, float top){
      fMarginL = left;
      fMarginR = right;
      fMarginT = top;
      fMarginB = bottom;
    }

    void setLeftMargin(float left){
      fMarginL = left;
    }

    void setRightMargin(float right){
      fMarginR = right;
    }

    void setTopMargin(float top){
      fMarginT = top;
    }

    void setBottomMargin(float bottom){
      fMarginB = bottom;
    }

  
    void setOffset(float xaxis, float yaxis, float zaxis){
      fOffsetX = xaxis;
      fOffsetY = yaxis;
      fOffsetZ = zaxis;
    }
    void setLabelSize(float xaxis, float yaxis, float zaxis){
      fLabelSizeX = xaxis;
      fLabelSizeY = yaxis;
      fLabelSizeZ = zaxis;
    }
 
    /*
     * plot mode (pMode) = 
     * <=0: png only
     * 1: pdf only
     * 2: eps only
     * 3: png + pdf
     * 4: png + eps
     * >=5: png + pdf + eps
     */
    bool drawHist2D(TH2F * h2, int pMode = 0, float ValMin = 999., float ValMax = -999.);
    bool drawHistsUpDown(std::vector<TH1F*> h_ups, std::vector<TH1F*> h_downs);
};
#endif

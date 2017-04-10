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

#include "HFrame.h"

class HPlot{
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
      fSizeY = 400.;
      fMarginL = 0.06;
      fMarginR = 0.11;
      fMarginT = 0.05;
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
 
    bool drawFrame2D(HFrame * m_fbox, float ValMin = 999., float ValMax = -999.);
    bool drawFramePipe2D(HFrame * m_fbox, float ValMin = 999., float ValMax = -999.);
    bool drawFrame2D(TH2F * h2, float ValMin = 999., float ValMax = -999.);

    // draw all the lines, having word "lname" into one plot
    bool drawOneLine(HFrame * m_fbox, std::string lname);

    // draw two lines, up, down
    bool drawTwoPipes(HFrame * m_fbox, std::string lname);
    bool drawTwoPipes(HFrame * m_fbox, std::string lname_up, std::string lname_down);
    bool drawTwoPipes(TH1F* h_up, TH1F* h_down);

    //one for each line
    bool drawAllLine(HFrame * m_fbox);
};
#endif

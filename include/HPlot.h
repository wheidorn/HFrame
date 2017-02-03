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

  public:
    HPlot(){ fCurCanvas1D = NULL; fCurCanvas2D = NULL;}
    ~HPlot(){ 
      if(!fCurCanvas1D) delete fCurCanvas1D; 
      if(!fCurCanvas2D) delete fCurCanvas2D; 
    }

    bool drawFrame2D(HFrame * m_fbox, float ValMin = 999., float ValMax = -999.);

    // draw all the lines, having word "lname" into one plot
    bool drawOneLine(HFrame * m_fbox, std::string lname);
    //
    //one for each line
    bool drawAllLine(HFrame * m_fbox);
    //TCanvas * getCanvas2D(){ return fCurCanvas2D; }
    //TCanvas * createCanvas();
    void getCanvasMarginOffset(float frac, float &xcan, float &margin_left, float &margin_right, float &margin_bottom, float &margin_top, float &xoff, float &yoff, float &zoff, bool is_2D = true);

};
#endif

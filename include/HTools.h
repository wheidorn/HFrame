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

#include "HFrame.h"

class HTools{

  private:

  public:
    HTools(){}
    ~HTools(){}

    // fitting stave cooling pipe temperature profile
    // at each x in the range [x_min, x_max]
    // get TH1F for [1, nybins]
    // fit two gaus [y1-rangefity, y1+rangefity] and [y2-rangefity, y2+rangefity] 
    // store the information of the fitted results
    //
    // if T is negative, instead of peak, we are looking for dips, so need to first covert data into positive!!!
    //bool fitPipeTwoGaus(HFrame * fbox, std::string s_out, int x_min, int x_max, int y1, int y2, int range_fity, bool negative);

    bool fitPipeTwoGaus(HFrame * fbox, std::string s_out, int x_min, int x_max, int y_min, int y_max, bool negative, bool convert_to_cm = true);

      // calculate the average T over range [y1, y2], return the TH1F across x range.
    TH1F* averagePipe(HFrame * fbox, std::string s_out, int x_min, int x_max, int y1, int y2, bool negative, bool convert_to_cm = true);
    HFrame* averageFrames(std::vector<HFrame *> fboxes, std::string shout = "T_average");
};
#endif

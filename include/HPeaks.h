#ifndef HPEAKS_H
#define HPEAKS_H

#include <sstream>
#include <string>
#include <iostream>
#include <map>
#include <algorithm>
#include <utility>
#include <vector>
#include "TH1F.h"
#include "TMath.h"

#include "HFrame.h"

class HPeaks{

  private:

  public:
    HPeaks(){}
    ~HPeaks(){}

    // starting from a pointer to the frame and knowing the histogram name of the peaks you would like to search for.
    //
    bool getPeaks(HFrame * fbox, std::string hname, std::string s_out, bool negative, double sigma = 2., double threshold=0.15, int niters=15); //20

    // more general, if you have a histogram and want to find its peaks, do the following.
    bool getPeaks(TH1F* hp, std::string s_out, bool negative, double sigma = 2., double threshold=0.15, int niters=15); //20
};
#endif

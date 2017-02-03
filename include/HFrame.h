#ifndef HFRAME_H
#define HFRAME_H

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

//double cmperpix = 128./550.; // 1.28 meter corresponds to 550 pixels, aproximately
const double cmperpix = 150./640.; // 1.28 meter corresponds to 550 pixels, aproximately
const int NXPIX = 640, NYPIX = 480;

class HFrame{

  private:
    TH2F* h_frame;
    std::vector<TH1F*> h_lines;
    std::string m_name;
    std::string m_title;
    std::string m_xtitle;
    std::string m_ytitle;
    TFile* fs; // source file

  public:
    HFrame(){
      h_frame = NULL;
    }

    HFrame(TFile* f, TH2F* h) : 
      fs( f ),
      h_frame( h ){ }
 
    ~HFrame(){
      if(h_frame) delete h_frame;
      if(fs){ fs->Close(); delete fs; }
      for(auto h1 : h_lines) delete h1;
    }

    enum angle {
      C090 =  90,
      C180 = 180,
      C270 = 270
    };
    void setTitle(std::string title){ m_title = title; }
    std::string getTitle(){ return m_title; }

    TH2F* getFrame(){ return h_frame; }
    bool setFrame(TH2F * h2, std::string shout = "T_average");

    //any line at any direction
    //using larger distance in X or Y as number of bins
    // ---- x1, y1 --------------------------------
    // ----------------------- x2, y2 -------------
    //
    // ----------------------- x2, y2 -------------
    // ---- x1, y1 --------------------------------
    TH1F* getLine(int x1, int y1, int x2, int y2, bool m_record = false); // to record line into HFrame
    int getNbinX(){if(!h_frame) return 0; else return h_frame->GetNbinsX();}
    int getNbinY(){if(!h_frame) return 0; else return h_frame->GetNbinsY();}
    void recordLine(TH1F* h0){ h_lines.push_back(h0);}

    // file the frame box with TTree
    bool fillboxTrees(std::vector<std::string> sfile, std::string shout = "T_average", std::string stree="atree", 
        std::string stemp="temperature", std::string sxbrh="xpos", std::string sybrh="ypos");

    // fill the frame box with a 2D histogram
    bool fillboxHist(std::string sfile, std::string sname = "T_average", bool convert_to_cm = true);
    bool rotate(HFrame::angle ang_clock); //support clockwise: 90, 180, 270
    bool mirror( bool byX=true); // convert to a mirror image either by Xaxis or Yaxis

    // crop to:
    // --------------------------
    // --- x1,y1 ------ x2,y2 ---
    // --------------------------
    // ----- x3,y3 ------ nanan -
    // --------------------------
    bool crop( int x1, int y1, int x2, int y2, int x3, int y3, bool convert_to_cm = true);
    bool shift( float T_shift); 
    // set X Y starting from 1.
    //bool resetXY(); 
    bool writeToRootFile(std::string sfile, std::string sname = "T_average");
};
#endif

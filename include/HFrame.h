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

//const double cmperpix = 150./640.; // 1.28 meter corresponds to 550 pixels, aproximately
const double CMPERPIXELPRE = 128./550.; // 1.28 meter corresponds to 550 pixels, aproximately
const int NXPIXPRE = 640, NYPIXPRE = 480;

class HFrame{

  private:
    TH2F* h_frame;
    TH2F* h_fpipe; // frame including cooling pipes; sub of h_frame
    std::map<std::string, TH1F*> m_pipemap;
    float T_liquid; // flowing liquid temperature
    bool conversionToCM;
    int NXPIX;
    int NYPIX;
    bool side_L; // "L" or "J"
    std::pair<int, int> crop_point1; // point 1 pixel number of X, Y in raw frame after crop
    std::pair<int, int> crop_point2;
    std::pair<int, int> pipe_point1;
    std::pair<int, int> pipe_point2;
    double cmperpix;
    int lens;
    float guessTLiquid();

  public:
    HFrame(){
      h_frame = NULL;
      h_fpipe = NULL;
      T_liquid = -999.;
      conversionToCM = true;
      NXPIX = NXPIXPRE;
      NYPIX = NYPIXPRE;
      side_L = true;
      crop_point1 = std::make_pair(-1, -1);
      crop_point2 = std::make_pair(-1, -1);
      pipe_point1 = std::make_pair(-1, -1);
      pipe_point2 = std::make_pair(-1, -1);
      cmperpix = CMPERPIXELPRE;
      lens = 80;
    }

    HFrame(TH2F* kframe) : 
      h_frame( kframe ){ }
    // file the frame with TTree from a list of input files
    HFrame(std::vector<std::string> sfile, int n_xpixels = 640, int n_ypixels = 480, std::string shout = "T_average", std::string stree="atree", 
        std::string stemp="temperature", std::string sxbrh="xpos", std::string sybrh="ypos", bool convert_to_cm = true);

    // fill the frame box with a 2D histogram
    HFrame(std::string sfile, std::string sname = "T_average", bool convert_to_cm = true);
  
    ~HFrame(){
      if(h_frame) delete h_frame;
      if(h_fpipe) delete h_fpipe;
      for(std::map<std::string, TH1F*>::iterator mpip = m_pipemap.begin(); mpip != m_pipemap.end(); mpip++){
        if(!mpip->second) delete mpip->second;
      }
    }

    void Print(){
      std::cout<<"---------------------------------------------" <<std::endl;
      std::cout<<"---- Information about the stave frame   ----" <<std::endl;
      std::cout<<"---------------------------------------------" <<std::endl;
      std::cout<<"-- Liquid Temperature: "<< T_liquid <<std::endl;
      std::cout<<"-- Convert to cm: "<< (conversionToCM ? "Yes" : "No" ) <<std::endl;
      std::cout<<"-- frame name: "<< h_frame -> GetName() <<std::endl;
      std::cout<<"-- frame side: "<< (side_L ? "L" : "J") <<std::endl;
      std::cout<<"-- frame number of pix, X: "<< NXPIX << ", Y: "<< NYPIX <<std::endl;
      std::cout<<"-- frame from raw, X1: "<<crop_point1.first<<" Y1: "<<crop_point1.second<<"; X2: "<<crop_point2.first <<" Y2: "<< crop_point2.second <<std::endl;
      std::cout<<"-- pipe region frame name: "<< h_fpipe -> GetName() <<std::endl;
      std::cout<<"-- pipe region from cropped, X1: "<<pipe_point1.first<<" Y1: "<<pipe_point1.second<<"; X2: "<<pipe_point2.first <<" Y2: "<< pipe_point2.second <<std::endl;
      for (auto pone : m_pipemap){
        std::cout<<"---- pipe results: "<<pone.first<<std::endl;
      }
      std::cout<<"---------------------------------------------" <<std::endl;
    }

    enum angle {
      C090 =  90,
      C180 = 180,
      C270 = 270
    };

    void setLens(int idegree){ lens = idegree;}
    int  getLens(){ return lens;}
    void setCMperPixel(double icmperpix){ cmperpix = icmperpix;}
    double getCMperPixel(){ return cmperpix;}
    void setConversionCM(bool isCM){conversionToCM = isCM;}
    bool isConversionCM(){return conversionToCM;}
    void setTLiquid(float kTliquid){ T_liquid = kTliquid;}
    float getTLiquid(){ return T_liquid;}
    bool isNegative(); 

    TH2F* getFrame(){ return h_frame; }
    void  setFrame(TH2F* kframe){ if(!h_frame) delete h_frame; h_frame = kframe; }
    TH2F* getFramePipe(){ return h_fpipe; }
    void  setSideL(bool sside){ side_L = sside; }
    bool  isSideL(){ return side_L; }
    std::pair<int, int> getCropPointOne(){return crop_point1;}
    std::pair<int, int> getCropPointTwo(){return crop_point2;}
    std::pair<int, int> getPipePointOne(){return pipe_point1;}
    std::pair<int, int> getPipePointTwo(){return pipe_point2;}

    //any line at any direction
    //using larger distance in X or Y as number of bins
    // ---- x1, y1 --------------------------------
    // ----------------------- x2, y2 -------------
    //
    // ----------------------- x2, y2 -------------
    // ---- x1, y1 --------------------------------
    TH1F* getLine(int x1, int y1, int x2, int y2, bool is_fpipe); 
    void recordPipeInfo(std::string hname, TH1F* h0){ m_pipemap.insert(std::make_pair(hname,h0));}
    TH1F* getPipeInfo(std::string hname){if(m_pipemap.find(hname) == m_pipemap.end()) return NULL; else return m_pipemap[hname];}

    bool rotate(HFrame::angle ang_clock); //support clockwise: 90, 180, 270
    bool mirror( bool byX=true); // convert to a mirror image either by Xaxis or Yaxis
    bool shift( float T_shift); 

    // crop to:
    // --------------------------
    // --- x1,y1 ------ x2,y2 ---
    // --------------------------
    // --- x3,y3 ------ x4,y4 ---
    // --------------------------
    bool crop( int x1, int y1, int x4, int y4, bool to_fpipe = false);

    // crop to:
    // --------------------------
    // --- x1,y1 ------ x2,y2 ---
    // --------------------------
    // ----- x3,y3 ------ nanan -
    // --------------------------
    //
    // if to_fpipe = true, h_frame is unchanged, h_fpipe is the result of the cropped frame
    // otherwise, h_frame is cropped.
    bool crop( int x1, int y1, int x2, int y2, int x3, int y3, bool to_fpipe = false);

    // find the sub frame of xpix * ypix having max/min average temperature
    // depending on high / low T experiment.
    // for 85 degree lens, xpix = 550 out of 640, ypix = 50 out of 480 the central region
    // should add the end of the stave card on top of it.
    //
    // |--------------------------------------------------------------|  
    // |                                                              |
    // |                                                              |
    // |  |--|                                                        |
    // |  |--------------------------------------------------------|  |
    // |  |                                                        |  |
    // |  |--------------------------------------------------------|  |
    // |--------------------------------------------------------------|
    //
    bool cropRawToStave( float xsize_cm = 128., float ysize_cm = 11.6);
    //
    // provide the pipe size in X and distance between two pipes in Y
    // and the pipe to edge size that you would like to keep in the frame
    bool cropStaveToPipe( float xsize_cm = 120., float ysize_cm = 5.0, float yedge_cm = 1.2); //was 124.
    // set X Y starting from 1.
    //bool resetXY(); 
    bool writeToRootFile(std::string sfile);
};
#endif

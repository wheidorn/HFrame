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

class HFrame
{
  /* class HFrame aim to do
   * - convert a raw image from IR camera to a frame, which fits the stave 
   * --- see function findCropPositions(...) by feeding x size (cm) and y size of the frame.
   * --- the coordinates is from 0 (left) of the frame.
   *
   *       |-|
   *       | --------|
   * (0,0) |---------| (Xmax,0)
   *
   *
   * - zoom a frame image to a smaller area including only two pipes.
   * --- see function findPipePositions(...) by feeding x size (cm) and y size of the pipe zoom area.
   *
   * - rotate the frame by 90, 180, 270 degrees
   * --- see function rotate(...)
   *
   * - mirror reflection of the frame
   * --- see function mirror(...)
   *
   * - store the results
   * --- see functions writeToRootFile(...) and print(...)
   *
   */
  private:
    /*
     * temperature data in 2D
     */
    std::vector<std::vector<float>> m_data;

    /* 
     * frame name
     */
    std::string m_name;

    /* 
     * output histogram names
     */
    std::string m_name_hist;
    std::string m_name_crop;
    std::string m_name_pipe;

    /* 
     * frame side: L or J
     */
    char  m_side; 

    /* 
     * lens: 80, 45, 25
     */
    int   m_lens;

    /* 
     * temperature of liquid
     */
    float m_T_liquid; 

    bool m_conversionToCM;
    const int m_NXPixel;
    const int m_NYPixel;
    const float m_CMperPixel;

    /* 
     * pixel position of the cropped frame @ original frame
     */
    int m_crop_x1; 
    int m_crop_y1; 
    int m_crop_x2;
    int m_crop_y2;

    /* 
     * pixel position of the cropped frame (cooling pipe) @ original frame
     */
    int m_pipe_x1; 
    int m_pipe_y1; 
    int m_pipe_x2;
    int m_pipe_y2;

  public:
    /*
     * constructor, providing:
     * - number of X pixels (default: 640)
     * - number of Y pixels (default: 480)
     * - cm per pixel (default: 128./550.) or 150./640.?
     */
    HFrame( std::string name = "a_frame", int NXPixel = 640, int NYPixel = 480, float CMperPixel = 0.23272727) : 
      m_NXPixel(NXPixel), 
      m_NYPixel(NYPixel), 
      m_CMperPixel(CMperPixel){

        m_name = name;
        m_data = std::vector<std::vector<float>> (NXPixel, std::vector<float>(NYPixel, 0.));
        m_T_liquid = -999.;
        m_conversionToCM = true;
        m_side = '0';
        m_lens = 80;
        m_crop_x1 = -1;
        m_crop_y1 = -1;
        m_crop_x2 = -1;
        m_crop_y2 = -1;
        m_pipe_x1 = -1;
        m_pipe_y1 = -1;
        m_pipe_x2 = -1;
        m_pipe_y2 = -1;
        m_name_hist = "T_average";
        m_name_crop = "T_stave";
        m_name_pipe = "T_stave_pipe";
      }
 
    ~HFrame(){
      m_data.clear();
    }

    /*
     * read frame information from files holding TTree 
     * should provide filenames, tree name, branch names and output 
     */
    bool read(std::string sfile_one, std::string stree="atree", 
        std::string stemp="temperature", std::string sxbrh="xpos", std::string sybrh="ypos", float scale = 1.0);
    bool read(std::vector<std::string> sfile, std::string stree="atree", 
        std::string stemp="temperature", std::string sxbrh="xpos", std::string sybrh="ypos");


    /*
     * fill the frame box with a 2D histogram
     */
    bool read(std::string sfile, std::string sname = "T_average");


    enum angle {
      C090 =  90,
      C180 = 180,
      C270 = 270
    };

    void setLens(int idegree){ m_lens = idegree;}
    int  getLens(){ return m_lens;}

    float getCMperPixel(){ return m_CMperPixel;}

    void setConversionCM(bool isCM){m_conversionToCM = isCM;}
    bool isConversionCM(){return m_conversionToCM;}

    void setTLiquid(float kTliquid){ m_T_liquid = kTliquid;}
    float getTLiquid(){ return m_T_liquid;}

    bool isNegative(); 
    void clearData(){ m_data.clear();}

    void  setSide(char side){ m_side = side; }
    bool  isSideL(){ if(m_side == 'L') return true; else return false; }
    void  getCropPoints(int& crop_x1, int& crop_y1, int& crop_x2, int& crop_y2){
      crop_x1 = m_crop_x1;
      crop_x2 = m_crop_x2;
      crop_y1 = m_crop_y1;
      crop_y2 = m_crop_y2;
    }
    void  getPipePoints(int& pipe_x1, int& pipe_y1, int& pipe_x2, int& pipe_y2){
      pipe_x1 = m_pipe_x1;
      pipe_x2 = m_pipe_x2;
      pipe_y1 = m_pipe_y1;
      pipe_y2 = m_pipe_y2;
    }


    void setCropPoints(int crop_x1, int crop_y1, int crop_x2, int crop_y2){ 
      m_crop_x1 = crop_x1;
      m_crop_x2 = crop_x2;
      m_crop_y1 = crop_y1;
      m_crop_y2 = crop_y2;
    }
    void setPipePoints(int pipe_x1, int pipe_y1, int pipe_x2, int pipe_y2){ 
      m_pipe_x1 = pipe_x1;
      m_pipe_x2 = pipe_x2;
      m_pipe_y1 = pipe_y1;
      m_pipe_y2 = pipe_y2;
    }
    void setOutHistNames(std::string name_hist, std::string name_crop, std::string name_pipe){
      m_name_hist = name_hist;
      m_name_crop = name_crop;
      m_name_pipe = name_pipe;
    }

    /*
     * frame temperature average
     * from point x1, y1 to point x2, y2
     */
    float average(int x1 = -1, int y1 = -1, int x4 = -1, int y4 = -1); 


    /*
     * frame rotation,
     * support clockwise: 90, 180, 270
     */
    bool rotate(HFrame::angle ang_clock); 

    /* 
     * frame mirror image
     * either by X or Y
     */
    bool mirror( bool byX=true); 

    /* 
     * frame shift a constant temeprature
     */
    bool shift( float T_shift); 
    bool crop( int x1, int y1, int x4, int y4, bool to_fpipe = false);

    /*
     * find the sub frame of xpix * ypix having max/min average temperature
     * depending on high / low T experiment.
     * for 85 degree lens, xpix = 550 out of 640, ypix = 50 out of 480 the central region
     * should add the end of the stave card on top of it.
     *
     * |--------------------------------------------------------------|  
     * |                                                              |
     * |                                                              |
     * |  |--|                                                        |
     * |  |--------------------------------------------------------|  |
     * |  |                                                        |  |
     * |  |--------------------------------------------------------|  |
     * |--------------------------------------------------------------|
     */

    /*
     * crop to:
     * --------------------------
     * ---------------- x2,y2 ---
     * --------------------------
     * --- x1,y1 ----------------
     * --------------------------
     */

    bool findCropPositions( float xsize_cm = 126., float ysize_cm = 11.6);

    /*
     * provide the pipe size in X and distance between two pipes in Y
     * and the pipe to edge size that you would like to keep in the frame
     */
    bool findPipePositions( float xsize_cm = 120., float ysize_cm = 5.0, float yedge_cm = 1.2); //was 124.

    /*
     * provide the starting and end points w.r.t the cropped stave
     */
    void setPipePositions(float x1_cm, float y1_cm, float x2_cm, float y2_cm);

    /*
     * write the frames into root
     */
    bool writeToRootFile(std::string sfile);
 
    /* 
     * write frame information into text
     */
    void print(std::string sout = "out.txt");
};

#endif

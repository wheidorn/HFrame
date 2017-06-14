#include <fstream>
#include <sstream>
#include <string>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include "tool.h"
using namespace std;

template <typename T>
string NumberToString ( T Number )
{
   ostringstream ss; 
   ss << Number;
   return ss.str();
}

bool tool::read_config(string sconfig){

  ifstream inconf(sconfig.c_str());
  string line;
  while (getline(inconf, line))
  {
    size_t k=line.find(":");
    if (line.size() <=0 || k ==string::npos) continue;
    if (line[0]=='#' || line[0]=='*') continue;

    if      (line.substr(0,k) =="Name" ) frameName = line.substr(k+1);
    else if (line.substr(0,k) =="Side") frameSide = line.substr(k+1);
    else if (line.substr(0,k) =="Tliquid") frameTliquid = stoi( string(line.substr(k+1)) );
    else cout<<"key "<<line.substr(0,k)<<"not known. Use: Name, Side, Tliquid. "<<endl;
  }
  if (frameTliquid <-998.) return false;
  if (frameSide != "J" && frameSide != "L") return false;

  return true;
}

bool tool::read_frame(string sfile, string sfolder){
  // check the existence of the raw averaged frame
  // if so, use them, if not, then make them.

  string m_sRawFile = sfolder + "/combined/" + sfile + ".root";

  // get the frame of averaged Temperature over time.
  HFrame * aframe = new HFrame();
  const unsigned int nframes = 200;
  vector<string> filenames;
  for (unsigned int jf=1; jf<=nframes; jf++){
    string sfone = sfolder + "/" + sfile + "/seq_n" + NumberToString(jf) + ".root";
    filenames.push_back(sfone);
  }
  aframe -> read( filenames );
  
  // make a mirror of the frame, 
  // (Y starts on top by default in the frame software, here
  // we want it start from bottom)
  if(frameSide == "J") aframe -> mirror();
  // no Mirror image for "L" side, trying to have the same image as J size

  aframe -> setSide( frameSide[0] );
  aframe -> setTLiquid( frameTliquid );

  //default: float xsize_cm = 128., ysize_cm = 11.6;
  aframe -> findCropPositions(126.);
  //default: float xsize_cm = 120., float ysize_cm = 5.0, float yedge_cm = 1.2
  //aframe -> findPipePositions(124., 5.4, 0.8);
  //x: 0 -- 122, y: 2 -- 10 (cm)
  aframe -> setPipePositions(0., 2., 122., 10.);

  // write out the raw frame
  aframe -> writeToRootFile(m_sRawFile);
  aframe -> print();

  return true;
}


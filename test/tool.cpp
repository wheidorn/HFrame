#include <fstream>
#include <sstream>
#include <string>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include "tool.h"
using namespace std;

// http://www.cplusplus.com/forum/articles/9645/
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

HFrame * tool::read_frame(string sfile, string sfolder){
  // check the existence of the raw averaged frame
  // if so, use them, if not, then make them.

  string m_sRawFile = sfolder + "/combined/" + sfile + ".root";
  bool m_rawExist = false;
  if (FILE *file = fopen(m_sRawFile.c_str(), "r")) {
    fclose(file);
    m_rawExist = true;
  }  

  // get the frame of averaged Temperature over time.
  HFrame * aframe = NULL;
  if(!m_rawExist){
    const unsigned int nframes = 200;
    vector<string> filenames;
    for (unsigned int jf=1; jf<=nframes; jf++){
      string sfone = sfolder + "/" + sfile + "/seq_n" + NumberToString(jf) + ".root";
      filenames.push_back(sfone);
    }
    aframe = new HFrame(filenames);

    // make a mirror of the frame, (Y starts on top by default in the frame software, here
    // we want it start from bottom)
    aframe -> mirror();

    // write out the raw frame
    aframe -> writeToRootFile(m_sRawFile);
  }else{
    aframe = new HFrame(m_sRawFile);
    // Tmp: Will's old T not flipped
    // do it now !!!!!!!!!!!!!!!!!!!!
    if(m_sRawFile.find("Will") != string::npos){
      aframe -> mirror();
    }
    if(m_sRawFile.find("QMUL") != string::npos){
      aframe -> rotate( HFrame::C180 );
    }
  }


  //aframe -> setConversionCM(true); // convert pixels to CM
  if(frameTliquid< -998.) { cout<<"WARNING: liquid temperature is not set, using "<<frameTliquid<<endl; }
  if(frameSide != "L" && frameSide != "J") { cout<<"WARNING: frame side not set to L or J, but: "<< frameSide <<endl; }

  aframe -> setTLiquid( frameTliquid );
  if (frameSide == "L") aframe -> setSideL(true);
  else aframe -> setSideL(false);
  //if(frameName != "") aframe -> getFrame() -> SetName(frameName.c_str());

  return aframe;
}

bool tool::read_pipe(string sprefix, string sfile, string sfolder){
  // check the existence of the pipe information file already 
  // if so, use them, if not return false.

  string m_sPipFile = sfolder + "/cropped/" + sfile + ".root";
  if (FILE *file = fopen(m_sPipFile.c_str(), "r")) {
    fclose(file);
  }else return false;
  TFile * f0 = new TFile(m_sPipFile.c_str(), "read");
  if (!f0){
    return false;
  }

  TH1F* h_temp_top = (TH1F*) f0 -> Get((sprefix+"_temperature_top").c_str()) -> Clone();
  TH1F* h_posi_top = (TH1F*) f0 -> Get((sprefix+"_peakposition_top").c_str()) -> Clone();
  TH1F* h_sigm_top = (TH1F*) f0 -> Get((sprefix+"_peakwidth_top").c_str()) -> Clone();
  TH1F* h_chi2_top = (TH1F*) f0 -> Get((sprefix+"_chi2overndf_top").c_str()) -> Clone();
  TH1F* h_temp_bot = (TH1F*) f0 -> Get((sprefix+"_temperature_bot").c_str()) -> Clone();
  TH1F* h_posi_bot = (TH1F*) f0 -> Get((sprefix+"_peakposition_bot").c_str()) -> Clone();
  TH1F* h_sigm_bot = (TH1F*) f0 -> Get((sprefix+"_peakwidth_bot").c_str()) -> Clone();
  TH1F* h_chi2_bot = (TH1F*) f0 -> Get((sprefix+"_chi2overndf_bot").c_str()) -> Clone();
  TH1F* h_temp = (TH1F*) f0 -> Get((sprefix+"_temperature").c_str()) -> Clone();
  TH1F* h_posi = (TH1F*) f0 -> Get((sprefix+"_peakposition").c_str()) -> Clone();
  TH1F* h_sigm = (TH1F*) f0 -> Get((sprefix+"_peakwidth").c_str()) -> Clone();
  TH1F* h_chi2 = (TH1F*) f0 -> Get((sprefix+"_chi2overndf").c_str()) -> Clone();
  if (!h_temp_top || !h_temp_bot || !h_temp) return false;
  h_temp_top -> SetDirectory(NULL);
  h_posi_top -> SetDirectory(NULL);
  h_sigm_top -> SetDirectory(NULL);
  h_chi2_top -> SetDirectory(NULL);
  h_temp_bot -> SetDirectory(NULL);
  h_posi_bot -> SetDirectory(NULL);
  h_sigm_bot -> SetDirectory(NULL);
  h_chi2_bot -> SetDirectory(NULL);
  h_temp -> SetDirectory(NULL);
  h_posi -> SetDirectory(NULL);
  h_sigm -> SetDirectory(NULL);
  h_chi2 -> SetDirectory(NULL);

  pipeMap.insert(make_pair(sprefix+"_temperature_top", h_temp_top));
  pipeMap.insert(make_pair(sprefix+"_peakposition_top", h_posi_top));
  pipeMap.insert(make_pair(sprefix+"_peakwidth_top", h_sigm_top));
  pipeMap.insert(make_pair(sprefix+"_chi2overndf_top", h_chi2_top));
  pipeMap.insert(make_pair(sprefix+"_temperature_bot", h_temp_bot));
  pipeMap.insert(make_pair(sprefix+"_peakposition_bot", h_posi_bot));
  pipeMap.insert(make_pair(sprefix+"_peakwidth_bot", h_sigm_bot));
  pipeMap.insert(make_pair(sprefix+"_chi2overndf_bot", h_chi2_bot));
  pipeMap.insert(make_pair(sprefix+"_temperature", h_temp));
  pipeMap.insert(make_pair(sprefix+"_peakposition", h_posi));
  pipeMap.insert(make_pair(sprefix+"_peakwidth", h_sigm));
  pipeMap.insert(make_pair(sprefix+"_chi2overndf", h_chi2));
 
  f0 ->Close();
  delete f0;
  return true;
}

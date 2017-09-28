#include <fstream>
#include <sstream>
#include <string>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include "tool.h"
#include "HPlot.h"
using namespace std;

// http://www.cplusplus.com/forum/articles/9645/
template <typename T>
string NumberToString ( T Number )
{
   ostringstream ss; 
   ss << Number;
   return ss.str();
}


class ThermalImpedance{
  /*
   Z = (T_object - T_ambient) / P
   P = Power = Energy / time = rho * ( V / t ) * c * deltaT = const * deltaT 
   where,
   rho = liquid desity = 1.52 g / ml = 1.52 Kg / L
   c = specific heat capacity = 1183 J / (kg * K)
   Q = liquid flow rate = V/t = 1.25 L / min
   deltaT = T_liquid_out - T_liquid_in

   const = 1.52 * 1183 * 1.25 / 60  ( J / (K * s) ) = 37.462 J / (K * s) = ~ watt / K
   */
  public: 
  ThermalImpedance(float par_liquid = 37.462) : m_par_liquid(par_liquid)
  {
  }
  float getImpedance(float T_obj, float T_amb, float deltaT)
  {
    float imped = 1e10;
    if(deltaT==0.) return imped;
    // ERROR: imped = (T_obj - T_amb) / deltaT;
    imped = (T_obj - T_amb) / ( m_par_liquid * deltaT);
    return imped;
  }
  float getTemperatureBox(float T_chiller){
   return ( 18.1 + 0.133 * T_chiller);
  }
  float getTemperatureInlet(float T_chiller){
   return ( 5.78 + 0.850 * T_chiller);
  }
  float getTemperatureOutlet(float T_chiller){
   return ( 7.97 + 0.802 * T_chiller);
  }



  private:
    const float m_par_liquid;
};


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

bool tool::draw_frame(string sfile, string sfolder){

  HPlot * pHP = new HPlot();
  // read 2D hists
  string m_sRawFile = sfolder + "/combined/" + sfile + ".root";
  TFile * ff0 = new TFile(m_sRawFile.c_str(), "read");
  string stitle = frameSide + Form(" side, T_{liquid} = %3.0f #circ C", frameTliquid);

  TH2F* h_stave = (TH2F *) ff0 -> Get("T_stave");
  h_stave -> SetName("stave"); // output png name
  h_stave -> GetXaxis() -> SetTitle( "X direction (cm)" );
  h_stave -> GetYaxis() -> SetTitle( "Y direction (cm)" );
  h_stave -> SetTitle( stitle.c_str() );
  bool ok = pHP -> drawHist2D(h_stave);


  // new plot setup
  HPlot * pHQ = new HPlot();
  pHQ -> setLeftMargin(0.08);
  pHQ -> setRightMargin(0.14);
  TH2F* h_pipe = (TH2F *) ff0 -> Get("T_stave_pipe");
  if(h_pipe){
  h_pipe -> SetName("pipe"); // output png name
  h_pipe -> GetXaxis() -> SetTitle( "X direction (cm)" );
  h_pipe -> GetYaxis() -> SetTitle( "Y direction (cm)" );
  h_pipe -> SetTitle( stitle.c_str() );
  }else{
    ok = false;
  }
  ok = ok & pHQ -> drawHist2D(h_pipe);

  // calculate the thermal impedance with pipe frame
  ThermalImpedance * pIP = new ThermalImpedance();
  float Tbox = pIP -> getTemperatureBox(frameTliquid);
  float Tinl = pIP -> getTemperatureInlet(frameTliquid);
  float Tout = pIP -> getTemperatureOutlet(frameTliquid);
  float Tdelta = Tinl - Tout;


  TH2F* h_impd = (TH2F*) h_pipe -> Clone( "impedance" );
  int nbinx = h_pipe -> GetNbinsX(), nbiny = h_pipe -> GetNbinsY();
  float ImpMean = 0., PowMean;
  for(int ix=1; ix<=nbinx; ix++){
    for(int iy=1; iy<=nbiny; iy++){
      float Tobject = h_pipe->GetBinContent(ix, iy);
      float Zimpd = pIP -> getImpedance( Tobject, Tbox, Tdelta);
      h_impd -> SetBinContent(ix, iy, Zimpd);
      ImpMean += Zimpd / (nbinx * nbiny);

      float Zpow = (Tobject - Tbox) / Zimpd;
      PowMean += Zpow / (nbinx * nbiny);
    }
  }
  ok = ok & pHQ -> drawHist2D(h_impd);
  cout<<"T liquid set up = "<<frameTliquid<<" T box = "<<Tbox<<" T inlet = "<<Tinl<<" T outlet = "<<Tout<<" avg InOut "<<(Tinl + Tout)/2.<< " avg power "<<PowMean <<", average impedance: "<<ImpMean<<endl;
 
  delete pHP;
  delete pHQ;
  return ok;
}


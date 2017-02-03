#define HPLOT_CPP
#include "HPlot.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TF1.h"
#include <math.h>

void HPlot::getCanvasMarginOffset(float frac, float &xcan, float &margin_left, float &margin_right, float &margin_bottom, float &margin_top, float &xoff, float &yoff, float &zoff, bool is_2D){
  if( frac >= 1.2 && frac < 1.5){
    xcan = 700; margin_left = 0.09; margin_right = 0.14; margin_top = 0.05; margin_bottom = 0.12;
    xoff = 1.45; yoff = 0.8; zoff = 1.2;
  }
  else if( frac >= 1.5 && frac < 2){
    xcan = 800; margin_left = 0.08; margin_right = 0.135; margin_top = 0.05; margin_bottom = 0.12;
    xoff = 1.4; yoff = 0.7; zoff = 1.1;
  }
  else if( frac >= 2 && frac < 3){
    xcan = 1000; margin_left = 0.06; margin_right = 0.12; margin_top = 0.05; margin_bottom = 0.12;
    xoff = 1.45; yoff = 0.6; zoff = 0.9;
  }
  else if( frac >= 3 && frac < 4){
    xcan = 1200; margin_left = 0.05; margin_right = 0.11; margin_top = 0.05; margin_bottom = 0.12;
    xoff = 1.35; yoff = 0.5; zoff = 0.85;
  }
  else if( frac >= 4 && frac < 6){
    xcan = 1400; margin_left = 0.04; margin_right = 0.1; margin_top = 0.05; margin_bottom = 0.12;
    xoff = 1.2; yoff = 0.4; zoff = 0.75;
  }
  else if( frac >= 6){
    xcan = 1600; margin_left = 0.035; margin_right = 0.1; margin_top = 0.05; margin_bottom = 0.12;
    xoff = 1.1; yoff = 0.35; zoff = 0.7;
  }
  if(!is_2D) margin_right = margin_left / 2.;
}


bool HPlot::drawFrame2D(HFrame * fbox, float ValMin, float ValMax){
  if(!fbox) return false;

  TH2F * h2 = fbox -> getFrame();
  if(!h2) return false;

  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);

  int nbinx = h2 -> GetNbinsX();
  int nbiny = h2 -> GetNbinsY();
  std::string fName = h2->GetName();
  size_t itt = fName.find("_copy");
  fName = fName.substr(0, itt);

  if (nbinx < 1 || nbiny < 1) return false;

  float xcan = 600, ycan = 600, margin_left = 0.1, margin_right = 0.15, margin_top = 0.05, margin_bottom = 0.12;
  float frac = float(nbinx) / nbiny;
  float xoff = 1.5, yoff = 0.85, zoff = 1.3;
  bool is_2D = true;
  this->getCanvasMarginOffset(frac, xcan, margin_left, margin_right, margin_bottom, margin_top, xoff, yoff, zoff, is_2D);

  std::string cname = std::string("can2D_") + fName;
  TCanvas * c0 = new TCanvas(cname.c_str(),"",xcan, ycan);
  c0->SetMargin( margin_left, margin_right, margin_bottom, margin_top);


  if (nbinx < 1 || nbiny < 1) return false;

  if (ValMin < ValMax){
    h2->SetMaximum(ValMax);
    h2->SetMinimum(ValMin);
  }
  h2->SetZTitle("Temperature (#circ C)");
  h2->GetZaxis()->SetTitle("Temperature (#circ C)");
  h2->GetXaxis()->SetTitleOffset( xoff );
  h2->GetYaxis()->SetTitleOffset( yoff );
  h2->GetZaxis()->SetTitleOffset( zoff );
  h2->Draw("SURF2");
  std::string sout="frame2D_";
  c0->Print( (sout+fName+".png").c_str());

  h2->SetTitle("");
  h2->Draw("COLZ");
  TLatex Tl; Tl.SetTextSize(20); Tl.SetTextFont(43);
  Tl.DrawLatex(h2->GetXaxis()->GetXmax()*0.7, h2->GetYaxis()->GetXmax() * 1.01 , h2->GetTitle());
  sout += "COLZ_";
  c0->Print( (sout+fName+".png").c_str());

  delete c0;
  return true;
}



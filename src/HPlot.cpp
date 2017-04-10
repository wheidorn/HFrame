#define HPLOT_CPP
#include "HPlot.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TF1.h"
#include <math.h>
#include "AtlasLabels.h"
#include "AtlasStyle.h"
#include "AtlasUtils.h"


bool HPlot::drawOneLine(HFrame * m_fbox, std::string lname){

  return true;
}
bool HPlot::drawTwoPipes(HFrame * fbox, std::string lname){
  if(!fbox) return false;

  TH1F * h_up = fbox -> getPipeInfo(lname+"_top");
  TH1F * h_down = fbox -> getPipeInfo(lname+"_bot");
  if(!h_up || !h_down) return false;
  return this-> drawTwoPipes(h_up, h_down);
}

bool HPlot::drawTwoPipes(HFrame * fbox, std::string lname_up, std::string lname_down){
  if(!fbox) return false;

  TH1F * h_up = fbox -> getPipeInfo(lname_up);
  TH1F * h_down = fbox -> getPipeInfo(lname_down);
  if(!h_up || !h_down) return false;
  return this-> drawTwoPipes(h_up, h_down);
}

bool HPlot::drawTwoPipes(TH1F* h_up, TH1F* h_down){

  if(!h_up || !h_down) return false;

  TCanvas *c0 = NULL;
  TPad * fPad0 = NULL;
  TPad * fPad1 = NULL;
  c0 = new TCanvas("c0","c0",1000,400);
  fPad0 = new TPad("pad0","pad0",0,0.5,1,1,   0,0,0);
  fPad1 = new TPad("pad1","pad1",0,0,  1,0.5,0,0,0);
  fPad0 -> SetMargin(0.07, 0.01, 0.0,  0.04); 
  fPad1 -> SetMargin(0.07, 0.01, 0.25, 0.0); 
  h_up -> GetYaxis()->SetTitleSize(0.08);
  h_up -> GetYaxis()->SetLabelSize(0.08);
  h_up -> GetYaxis()->SetTitleOffset(0.285);
  h_down -> GetYaxis()->SetLabelSize(0.08);
  // one for normal distribution, the other for background subtracted
  fPad0 -> SetGrid();
  fPad1 -> SetGrid();
  c0 -> cd();
  fPad0 -> Draw();
  fPad1 -> Draw();
  fPad0 -> cd();

  std::string hname = h_up -> GetName();
  std::string htitleY = "Arbitrary";
  if (hname.find("peakposition") != std::string::npos) htitleY = "Peak Position (pix index)";
  else if(hname.find("peakwidth") != std::string::npos) htitleY = "Peak Width (N pixel)";
  else if(hname.find("chi2overndf") != std::string::npos) htitleY = "#chi^{2} / NDF";

  TH1F * hsm_up = (TH1F*) h_up -> Clone( (hname+"_smooth").c_str());
  hsm_up -> Smooth(3);
  hsm_up -> GetYaxis() -> SetNdivisions(505);
  hsm_up -> GetXaxis()->SetTitleSize(0.);
  hsm_up -> GetXaxis()->SetLabelSize(0.);
  hsm_up -> SetLineStyle(1);
  hsm_up -> SetLineWidth(2);
  hsm_up -> SetLineColor(kBlue);
  hsm_up -> GetYaxis()->SetTitle( htitleY.c_str());
  hsm_up -> GetYaxis()->SetTitleOffset(0.375);
  h_up -> SetLineStyle(3);
  h_up -> SetLineWidth(2);
  h_up -> SetLineColor(kBlue);
  hsm_up -> Draw();
  h_up -> Draw("same");
  TLatex fTxt; 
  fTxt.SetTextSize(0.12);
  fTxt.SetNDC();
  fTxt.DrawLatex(0.35, 0.7, "pipe inlet");
 

  fPad1 -> cd();

  TH1F * hsm_down = (TH1F*) h_down -> Clone( (hname+"_smooth").c_str());
  hsm_down -> Smooth(3);
  hsm_down -> GetYaxis() -> SetNdivisions(505);
  hsm_down -> GetXaxis()->SetTitleSize(0.1);
  hsm_down -> GetXaxis()->SetLabelSize(0.12);
  hsm_down -> GetXaxis()->SetTitle("Stave X direction (cm)");
  hsm_down -> GetXaxis()->SetTitleOffset(1.05);
  hsm_down -> SetLineStyle(1);
  hsm_down -> SetLineWidth(2);
  hsm_down -> SetLineColor(kBlue);
  h_down -> SetLineStyle(3);
  h_down -> SetLineWidth(2);
  h_down -> SetLineColor(kBlue);
  hsm_down -> Draw();
  h_down -> Draw("same");
  fTxt.DrawLatex(0.35, 0.7, "pipe outlet");
 

  std::string pname = "peaks_";
  size_t t0 = hname.find("_top");
  size_t t1 = hname.find("_bot");
  if     (t0 != std::string::npos) pname += hname.substr(0,t0);
  else if(t1 != std::string::npos) pname += hname.substr(0,t1);
  else   pname += hname;
  c0->Print( (pname+".png").c_str());
  c0->Print( (pname+".eps").c_str());

  delete fPad0;
  delete fPad1;
  delete c0;

  delete hsm_up;
  delete hsm_down;
  return true; 

}

bool HPlot::drawFrame2D(HFrame * fbox, float ValMin, float ValMax){
  if(!fbox) return false;

  TH2F * h2 = fbox -> getFrame();
  if(!h2) return false;
  return this -> drawFrame2D(h2, ValMin, ValMax);
}
bool HPlot::drawFramePipe2D(HFrame * fbox, float ValMin, float ValMax){
  if(!fbox) return false;

  TH2F * hp2 = fbox -> getFramePipe();
  if(!hp2) return false;
  return this -> drawFrame2D(hp2, ValMin, ValMax);

}

bool HPlot::drawFrame2D(TH2F* h2, float ValMin, float ValMax){
  if(!h2) return false;

  //gROOT->SetStyle("Plain");
  SetAtlasStyle();
  gStyle->SetOptStat(0);

  int nbinx = h2 -> GetNbinsX();
  int nbiny = h2 -> GetNbinsY();
  std::string fName = h2->GetName();
  if (nbinx < 1 || nbiny < 1) return false;
  std::string cname = std::string("can2D_") + fName;
  TCanvas * c0 = new TCanvas(cname.c_str(),"",fSizeX, fSizeY);
  c0->SetMargin( fMarginL, fMarginR, fMarginB, fMarginT);


  if (nbinx < 1 || nbiny < 1) return false;

  if (ValMin < ValMax){
    h2->SetMaximum(ValMax);
    h2->SetMinimum(ValMin);
  }
  h2->SetZTitle("Temperature (#circ C)");
  h2->GetZaxis()->SetTitle("Temperature (#circ C)");
  h2->GetXaxis()->SetTitleOffset( fOffsetX );
  h2->GetYaxis()->SetTitleOffset( fOffsetY );
  h2->GetZaxis()->SetTitleOffset( fOffsetZ );

  h2->GetXaxis()->SetLabelSize( h2->GetXaxis()->GetLabelSize() * fLabelSizeX );
  h2->GetYaxis()->SetLabelSize( h2->GetYaxis()->GetLabelSize() * fLabelSizeY );
  h2->GetZaxis()->SetLabelSize( h2->GetZaxis()->GetLabelSize() * fLabelSizeZ );
 
  h2->Draw("SURF2");
  std::string sout="frame2D_";
  c0->Print( (sout+fName+".png").c_str());

  h2->SetTitle("");
  h2->Draw("COLZ");
  TLatex fTxt; //Tl.SetTextSize(20); Tl.SetTextFont(43);
  fTxt.SetTextSize(0.04);
  fTxt.SetNDC();
  fTxt.DrawLatex(0.7, 0.98, h2->GetTitle());
  sout += "COLZ_";
  c0->Print( (sout+fName+".png").c_str());

  delete c0;
  return true;
}


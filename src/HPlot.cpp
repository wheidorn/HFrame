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

bool HPlot::drawHist2D(TH2F* h2, int pMode, float ValMin, float ValMax){
  if(!h2) return false;

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
  std::string sout="fSURF_";
  c0->Print( (sout+fName+".png").c_str());
  c0->Print( (sout+fName+".eps").c_str());

  //h2->SetTitle("");
  h2->Draw("COLZ");
  TLatex fTxt; //Tl.SetTextSize(20); Tl.SetTextFont(43);
  fTxt.SetTextSize(0.04);
  fTxt.SetNDC();
  fTxt.DrawLatex(0.7, 0.95, h2->GetTitle()); 
  sout = "fCOLZ_";

  if(pMode <=0     || pMode >=4 ) c0->Print( (sout+fName+".png").c_str());
  if(pMode %2 == 1 || pMode >=5 ) c0->Print( (sout+fName+".pdf").c_str());
  if(pMode ==2     || pMode >=4 ) c0->Print( (sout+fName+".eps").c_str());

  delete c0;
  return true;
}

bool HPlot::drawHistsUpDown(std::vector<TH1F*> h_ups, std::vector<TH1F*> h_downs)
{

  if(!h_ups.size() || !h_downs.size() || h_ups.size() != h_downs.size()) return false;

  int nhist = h_ups.size();

  TCanvas *c0 = NULL;
  TPad * fPad0 = NULL;
  TPad * fPad1 = NULL;
  c0 = new TCanvas("c0","c0",1000,400);
  fPad0 = new TPad("pad0","pad0",0,0.55,1,1,   0,0,0);
  fPad1 = new TPad("pad1","pad1",0,0,   1,0.5,0,0,0);
  fPad0 -> SetMargin(0.07, 0.01, 0.0,  0.04); 
  fPad1 -> SetMargin(0.07, 0.01, 0.25, 0.0); 
  h_ups.at(0) -> GetYaxis()->SetTitleSize(0.08);
  h_ups.at(0) -> GetYaxis()->SetLabelSize(0.08);
  h_ups.at(0) -> GetYaxis()->SetTitleOffset(0.285);
  h_downs.at(0) -> GetYaxis()->SetLabelSize(0.08);

  fPad0 -> SetGrid();
  fPad1 -> SetGrid();
  c0 -> cd();
  fPad0 -> Draw();
  fPad1 -> Draw();

//  TLegend * leg = new TLegend(0.55, 0.75, 0.8, 0.92);

  for(int ih=0; ih<nhist; ih++){
    TH1F* h_up = h_ups.at(ih);
    TH1F* h_down = h_downs.at(ih);

    std::string hname = h_up -> GetName();
    fPad0 -> cd();
    h_up -> Smooth(3);
    h_up -> GetYaxis() -> SetNdivisions(505);
    h_up -> GetXaxis()->SetTitleSize(0.);
    h_up -> GetXaxis()->SetLabelSize(0.);
    h_up -> SetLineStyle(1);
    h_up -> SetLineWidth(2);
    h_up -> SetLineColor(kBlue);
    //h_up -> GetYaxis()->SetTitle( htitleY.c_str());
    h_up -> GetYaxis()->SetTitleOffset(0.375);
    h_up -> SetLineStyle(3);
    h_up -> SetLineWidth(2);
    h_up -> SetLineColor(kBlue);

    TLatex fTxt; 
    if(!ih){ 
      h_up -> Draw();
      //h_up -> Draw("same");

      fTxt.SetTextSize(0.12);
      fTxt.SetNDC();
      fTxt.DrawLatex(0.12, 0.8, h_up -> GetTitle());

    }else h_up -> Draw("same");
    //leg -> AddEntry(hsm_up, h_up -> GetTitle(), "L");


    fPad1 -> cd();

    h_down -> Smooth(3);
    h_down -> GetYaxis() -> SetNdivisions(505);
    h_down -> GetXaxis()->SetTitleSize(0.1);
    h_down -> GetXaxis()->SetLabelSize(0.12);
    h_down -> GetXaxis()->SetTitle("Stave X direction (cm)");
    h_down -> GetXaxis()->SetTitleOffset(1.05);
    h_down -> SetLineStyle(1);
    h_down -> SetLineWidth(2);
    h_down -> SetLineColor(kBlue);
    h_down -> SetLineStyle(3);
    h_down -> SetLineWidth(2);
    h_down -> SetLineColor(kBlue);
    if(!ih){ 
    h_down -> Draw();
    //fTxt.DrawLatex(0.12, 0.8, "pipe at bottom");
    }else{
    h_down -> Draw("same");
    }

  }
//fPad0 -> cd();
//leg -> Draw();

    std::string pname = "peaks_";
    std::string hname = h_ups.at(0) -> GetName();
    size_t t0 = hname.find("_top");
    size_t t1 = hname.find("_bot");
    if     (t0 != std::string::npos) pname += hname.substr(0,t0);
    else if(t1 != std::string::npos) pname += hname.substr(0,t1);
    else   pname += hname;
  std::cout<<10<<" "<<c0<<std::endl;
  c0 -> cd();
    c0->Print( (pname+".png").c_str());
    c0->Print( (pname+".eps").c_str());
  std::cout<<11<<std::endl;

//    delete leg;
    delete fPad0;
    delete fPad1;
    delete c0;

  return true; 
}

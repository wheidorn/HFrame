#define HPEAKS_CPP
#include "HPeaks.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TROOT.h"
#include "TStyle.h"
#include <math.h>
#include "TGraph.h"

#include "TF1.h"
#include <stdio.h> // fprintf
#include "TSpectrum.h"
#include "TVirtualFitter.h"
#include "TPolyMarker.h"
#include "TCanvas.h"
#include "TPad.h"

int npeaks;
// used for peak finder!
double fpeaks(double *x, double *pars) {
  double result = pars[0] + pars[1]*x[0];
  for (int p=0;p<npeaks;p++) {
    double norm  = pars[3*p+2];
    double mean  = pars[3*p+3];
    double sigma = pars[3*p+4];
    result += norm*TMath::Gaus(x[0],mean,sigma);
  }
  return result;
}


bool HPeaks::getPeaks(HFrame * fbox, std::string hname, std::string s_out, bool negative, double sigma, double threshold, int niters){

  TH1F * hp = fbox -> getPipeInfo(hname);
  if(!hp){std::cout<<"ERROR: no pipe information: "<<hname<<" is found. Return false!"<<std::endl; return false;}

  return getPeaks(hp, s_out, negative, sigma, threshold, niters);
}

bool HPeaks::getPeaks(TH1F* hp, std::string s_out, bool negative, double sigma, double threshold, int niters){
  std::string hname = std::string(hp -> GetName());

  if(!hp){std::cout<<"ERROR: no pipe information: "<<hname<<" is found. Return false!"<<std::endl; return false;}
  gStyle->SetOptStat(0);
  FILE * fp = fopen((s_out+".txt").c_str(), "w");

  TH1F * hp_sm = (TH1F*) hp -> Clone( (hname+"_smooth").c_str());
  hp_sm -> Smooth(3);
  TH1F * hp2 = (TH1F*) hp_sm -> Clone( (hname+"_copy").c_str());
  if(!negative){
    hp2->Scale(-1.); // for non-negative, looking for drops, code is looking for peaks!!!
  }

  TCanvas *c0 = NULL;
  TPad * fPad0 = NULL;
  TPad * fPad1 = NULL;
  if(hname.find("_top") != std::string::npos || hname.find("_bot") != std::string::npos){
    c0 = new TCanvas("c0","c0",1000,400);
    //c0 -> SetMargin(0.06, 0.01, 0.25, 0.04); // left, right, bottom, top
    fPad0 = new TPad("pad0","pad0",0,0.5,1,1,   0,0,0);
    fPad1 = new TPad("pad1","pad1",0,0,  1,0.5,0,0,0);
    fPad0 -> SetMargin(0.07, 0.01, 0.0,  0.04); 
    fPad1 -> SetMargin(0.07, 0.01, 0.25, 0.0); 
    hp_sm -> GetYaxis()->SetTitleSize(0.08);
    hp_sm -> GetYaxis()->SetLabelSize(0.08);
    hp_sm -> GetYaxis()->SetTitleOffset(0.285);
  }else{
    c0 = new TCanvas("c0","c0",2000,400);
    //c0 -> SetMargin(0.03, 0.005, 0.25, 0.04); // left, right, bottom, top
    fPad0 = new TPad("pad0","pad0",0,0.5,1,1,   0,0,0);
    fPad1 = new TPad("pad1","pad1",0,0,  1,0.5,0,0,0);
    fPad0 -> SetMargin(0.035, 0.005, 0.0,  0.04); 
    fPad1 -> SetMargin(0.035, 0.005, 0.25, 0.0); 
    hp_sm -> GetYaxis()->SetTitleSize(0.1);
    hp_sm -> GetYaxis()->SetLabelSize(0.1);
    hp_sm -> GetYaxis()->SetTitleOffset(0.15);
  }
  // one for normal distribution, the other for background subtracted
  fPad0 -> SetGrid();
  fPad1 -> SetGrid();
  c0 -> cd();
  fPad0 -> Draw();
  fPad1 -> Draw();
  fPad0 -> cd();

  hp_sm -> GetYaxis() -> SetNdivisions(505);
  hp_sm -> GetXaxis()->SetTitleSize(0.);
  hp_sm -> GetXaxis()->SetLabelSize(0.);
  hp_sm -> GetXaxis()->SetTitleOffset(1.);
  hp_sm -> SetLineStyle(1);
  hp_sm -> SetLineWidth(2);
  hp_sm -> SetLineColor(kBlue);
  hp_sm -> GetYaxis()->SetTitle("Temperature #circ C");
  hp -> SetLineStyle(3);
  hp -> SetLineWidth(2);
  hp -> SetLineColor(kBlue);
  double Tmax = hp_sm -> GetMaximum();
  Tmax = int(Tmax) + (Tmax > 0 ? 3. : 2.);
  double Tmin = hp_sm -> GetMinimum();
  Tmin = int(Tmin) - (Tmin < 0 ? 1. : 0.);
  hp_sm -> SetMaximum( Tmax );
  hp_sm -> SetMinimum( Tmin );
  hp_sm -> Draw();
  hp -> Draw("same");


  double xmin = hp -> GetXaxis()->GetXmin();
  double xmax = hp -> GetXaxis()->GetXmax();
  npeaks = 20; // guess how many peaks on the curve.
  //Use TSpectrum to find the peak candidates
  TSpectrum *spkm = new TSpectrum(2*npeaks);
  //double sigma = 1.75, threshold = 0.2;
  // searching for the peaks!
  // using hp2 the cloned version to be avoid of drawing the points found
  int nfound = spkm->Search(hp2,sigma,"nodraw",threshold); //nobackgroundnomarkov

  double * xpeaks = spkm->GetPositionX();
  double * ypeaks = spkm->GetPositionY();
  if(!negative){
    for(int ip=0; ip<nfound; ip++) ypeaks[ip] *= -1.; // turn them back to positive
  }
  TPolyMarker * pm = new TPolyMarker(nfound, xpeaks, ypeaks);
  pm->SetMarkerStyle(23);
  pm->SetMarkerColor(kGray + 1);
  pm->SetMarkerSize(1.35);
  pm->Draw("same");


  //int niters = 20;
  //
  // looking for the background.
  TH1F *hpbkg = (TH1F*) spkm->Background(hp2,niters,""); //if "same", the background is plotted 
  if(!negative){
    hpbkg -> Scale(-1.);
  }
  hpbkg -> SetLineWidth(2);
  hpbkg -> SetLineStyle(3);
  hpbkg -> SetLineColor(kGray+2);
  hpbkg -> Draw("same");
  if (hpbkg) c0->Update();

  // only keep the peaks over background level!!!
  std::vector<double> xpkeep, ypkeep;
  for (int p=0;p<nfound;p++) {
    double xp = xpeaks[p];
    double yp = ypeaks[p];
    int ixbin = hpbkg -> GetXaxis() -> FindBin(xp);
    double yp_bkg = hpbkg -> GetBinContent(ixbin);
    if (yp_bkg ==0.) continue;
    //if ( fabs(yp - yp_bkg) < threshold * fabs(yp_bkg) ) continue;
    //
    //
    //if ( fabs(yp - yp_bkg) < .25 ) continue;
    xpkeep.push_back(xp);
    ypkeep.push_back(yp);
  }
  int npkeep = xpkeep.size();
  double * xpkeep_a = new double[npkeep];
  double * ypkeep_a = new double[npkeep];
  for(int ip =0; ip<npkeep; ip++){ 
    xpkeep_a[ip] = xpkeep.at(ip); 
    ypkeep_a[ip] = ypkeep.at(ip); 
  }
  TPolyMarker * pmkeep = new TPolyMarker(npkeep, xpkeep_a, ypkeep_a);
  pmkeep->SetMarkerStyle(23);
  pmkeep->SetMarkerColor(kRed);
  pmkeep->SetMarkerSize(1.45);
  pmkeep->Draw("same");

  TLatex fTxt; 
  fTxt.SetTextSize(0.12); //0.04
  fTxt.SetNDC();
  fTxt.DrawLatex(0.3, 0.85, Form("Bkg NIter=%d, peak sigma=%4.2f, threshold=%4.2f",niters,sigma,threshold));
 


  const int npar_max = 1000;
  double pars[npar_max];

  TF1 *fline = new TF1("fline","pol1",xmin, xmax);
  hpbkg->Fit("fline","qn");
  pars[0] = fline->GetParameter(0);
  pars[1] = fline->GetParameter(1);
  npeaks = 0; // reset all peaks
  //Double_t *xpeaks = spkm->GetPositionX();
  for (int p=0;p<nfound;p++) {
    Double_t xp = xpeaks[p];
    Int_t bin = hp_sm->GetXaxis()->FindBin(xp);
    Double_t yp = hp_sm->GetBinContent(bin);
    if (yp-TMath::Sqrt( fabs(yp) ) < fline->Eval(xp)) continue;
    pars[3*npeaks+2] = yp;
    pars[3*npeaks+3] = xp;
    pars[3*npeaks+4] = 1.; // guess the sigma
    npeaks++;
  }


  ///////   TF1 *fit = new TF1("fit_sb", fpeaks,xmin,xmax,2+3*npeaks);
  ///////   TVirtualFitter::Fitter(hp2,10+3*npeaks);
  ///////   fit->SetParameters(par);
  ///////   fit->SetNpx(1000);
  ///////   hp2->Fit("fit_sb");


  // background subtraction
  fPad1 -> cd();
  TH1F* hp_bsub = (TH1F*)hp_sm -> Clone((hname+"_sm_noBkg").c_str());
  hp_bsub -> GetXaxis()->SetTitleSize(0.1);
  hp_bsub -> GetXaxis()->SetLabelSize(0.12);
  hp_bsub -> GetYaxis()->SetTitle("T subtracted (#circ C)");
  hp_bsub -> GetXaxis()->SetTitle("Stave X direction (cm)");

  int nbinx = hp_bsub -> GetNbinsX();
  if(nbinx != hpbkg -> GetNbinsX()){
    std::cout<<"ERROR: histogram background has different number of bins than norminal: "<<nbinx<<" != "<<hpbkg -> GetNbinsX()<<std::endl;
    return false;
  }
  hp_bsub -> Add(hpbkg, -1.);

  double SBmax = hp_bsub -> GetMaximum();
  SBmax = int(SBmax) + (SBmax > 0 ? 2. : 1.);
  double SBmin = hp_bsub -> GetMinimum();
  SBmin = int(SBmin) - (SBmin < 0 ? 1. : 0.);
  hp_bsub -> SetMaximum( SBmax );
  hp_bsub -> SetMinimum( SBmin );
  hp_bsub -> Draw(); 


//  TH1F* hp_bsub2 = (TH1F*) hp_bsub -> Clone((hname+"_sm_noBkg2").c_str());
//  if (!negative) hp_bsub2 -> Scale(-1.);
//  TH1F *hpbkg_sb = (TH1F*) spkm->Background(hp_bsub2,niters,"");
//  if(!negative) hpbkg_sb -> Scale(-1.);
//  hpbkg_sb -> SetLineWidth(2);
//  hpbkg_sb -> SetLineStyle(2);
//  hpbkg_sb -> SetLineColor(kGreen+8);
//  hpbkg_sb -> Draw("same");
  
  std::string pname = "peaks_" + s_out + "_";
  pname += hname;
  c0->Print( (pname+".png").c_str());
  c0->Print( (pname+".eps").c_str());

  delete fPad0;
  delete fPad1;
  delete c0;

  fclose(fp);
  delete hp2;
  delete hpbkg;
  delete hp_bsub;
  delete hp_sm;
  return true;
}

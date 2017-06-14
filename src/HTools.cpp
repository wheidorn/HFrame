#define HTOOLS_CPP
#include "HTools.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TROOT.h"
#include "TStyle.h"
#include <math.h>
#include "TGraph.h"
#include <stdio.h> // fprintf
#include "TF1.h"
#include <iostream>
#include <string>
#include "TSpectrum.h"
#include "TPolyMarker.h"

using namespace std;
TH1F* HTools::getHist1D(TH2F* h2, int ibin, bool isY){
  TH1F * h1 = NULL;

  if(!h2) return h1;

  int nxbin = h2 -> GetNbinsX();
  int nybin = h2 -> GetNbinsY();

  if(nxbin<=0 || nybin<=0) return h1;
  
  string name = "h1_";
  name += to_string(ibin);
  if (isY){
    double bmin = h2 -> GetYaxis() -> GetXmin();
    double bmax = h2 -> GetYaxis() -> GetXmax();
    h1 = new TH1F(name.c_str(), "", nybin, bmin, bmax);
    for(int ib=1; ib<=nybin; ib++){
      h1 -> SetBinContent( ib, h2 -> GetBinContent(ibin, ib));
      h1 -> SetBinError( ib, h2 -> GetBinError(ibin, ib));
    }
  }else{
    double bmin = h2 -> GetXaxis() -> GetXmin();
    double bmax = h2 -> GetXaxis() -> GetXmax();
    h1 = new TH1F(name.c_str(), "", nxbin, bmin, bmax);
    for(int ib=1; ib<=nxbin; ib++){
      h1 -> SetBinContent( ib, h2 -> GetBinContent(ib, ibin));
      h1 -> SetBinError( ib, h2 -> GetBinError(ib, ibin));
    }
  }
  return h1;
}

bool HTools::fitFramePipesGaus(TH2F * hframe, std::string s_out,  float Tliquid, std::string s_side/* = "L"*/, string sfile /* = "a1_pipe.root" */)
{
  /*
  // positive: Tliquid > 30
  // negative: Tliquid < 15
  // in between: do not know if it is negative or positive
  */
  bool negative = true;
  if(Tliquid > 30) negative = false;
  else if (Tliquid >15) return false; 

  /*
  // Template for output histograms
  */
  gStyle->SetOptStat(0);
//(TH1F*) hframe -> ProjectionX() -> Clone("tmp");
  //h_tmp -> Reset();
  //
  TH1F * h_tmp = new TH1F("tmp","", hframe->GetNbinsX(), hframe->GetXaxis()->GetXmin(), hframe->GetXaxis()->GetXmax());

  /*
     1) Temperature peak value
     2) Peak pixel position
     3) Sigma of Gaussian out of fit
     4) chi2 / NDF out of fit
     */
  TFile * fout = new TFile(sfile.c_str(), "recreate");
  TH1F* h_temp_top = (TH1F*) h_tmp -> Clone((s_out+"_temperature_top").c_str());
  TH1F* h_posi_top = (TH1F*) h_tmp -> Clone((s_out+"_peakposition_top").c_str());
  TH1F* h_sigm_top = (TH1F*) h_tmp -> Clone((s_out+"_peakwidth_top").c_str());
  TH1F* h_chi2_top = (TH1F*) h_tmp -> Clone((s_out+"_chi2overndf_top").c_str());
  TH1F* h_temp_bot = (TH1F*) h_tmp -> Clone((s_out+"_temperature_bot").c_str());
  TH1F* h_posi_bot = (TH1F*) h_tmp -> Clone((s_out+"_peakposition_bot").c_str());
  TH1F* h_sigm_bot = (TH1F*) h_tmp -> Clone((s_out+"_peakwidth_bot").c_str());
  TH1F* h_chi2_bot = (TH1F*) h_tmp -> Clone((s_out+"_chi2overndf_bot").c_str());

  int nbinsx = h_tmp -> GetNbinsX();
  float hxmin = h_tmp -> GetXaxis() -> GetXmin(); 
  float hxmax = h_tmp -> GetXaxis() -> GetXmax(); 
  delete h_tmp;

  /*
  // fit two gaussian through line in Y at each X point
  // gaus: gaussian with 3 parameters: f(x) = p0*exp(-0.5*((x-p1)/p2)^2)).
  */
  const int npars=6;
  double pars[npars];
  for(int ix=1; ix<=nbinsx; ix++){

    TH1F* h1 = this -> getHist1D(hframe, ix); 
    // h1 -> Smooth(); // do smoothing before fitting??

    TAxis *h1xax = h1 -> GetXaxis();
    float binWidth = h1xax -> GetBinWidth(1);
    if (!h1){ std::cout<<"error finding the line in fit "<<std::endl; return false;}

    TCanvas *c1 = new TCanvas("c1","c1",600,600);
    //SetMargin (Float_t left, Float_t right, Float_t bottom, Float_t top)
    c1 -> SetMargin(0.12,0.03,0.12,0.03);
    double Maxi = int(h1->GetMaximum() ) + 1.; 
    double Mini = int(h1->GetMinimum() ) - 2.; 
    h1->SetMaximum( Maxi+6.);
    h1->SetMinimum( Mini);
    h1->SetLineWidth(3);
    h1->SetLineColor(1);
    h1->SetMarkerStyle(23);
    h1->SetMarkerColor(1);
    h1->GetXaxis() -> SetTitle("Y (cm)");
    h1->GetYaxis() -> SetTitle("Temperature (#circ C)");
    h1->GetXaxis() -> SetTitleOffset(1.2* h1->GetXaxis() -> GetTitleOffset());
    h1->GetYaxis() -> SetTitleOffset(1.3* h1->GetYaxis() -> GetTitleOffset());
    h1->Draw();

    /* 
    // find the maximum (negative temperature) or minimum (positive temperature) position of Y
    // use as the gaussian fit center to start with.
    // then set up the fit range around the min or max.
    */
    int ib_mmy_1 = 0., ib_mmy_2 = 0.; // mm = max or min
    double mmy_1 = -99., mmy_2 = -99.;  if(negative){mmy_1 = 9999.; mmy_2 = 9999.;}
    int y_nbin = h1->GetNbinsX();
    int y_thr = h1->GetNbinsX() / 3;

    // find minimum by averaging of 3 pixels
    for(int ib=2; ib<y_thr; ib++){
      int irb = y_nbin - ib + 1; // reversed order bin, starting from (y_nbin - 1)
      double avg_1 = h1->Integral(ib-1,   ib+1) / 3.;
      double avg_2 = h1->Integral(irb-1, irb+1) / 3.;
      if( ( negative && mmy_1 > avg_1) || (!negative && mmy_1 < avg_1) ){ mmy_1 = avg_1; ib_mmy_1 = ib;  }
      if( ( negative && mmy_2 > avg_2) || (!negative && mmy_2 < avg_2) ){ mmy_2 = avg_2; ib_mmy_2 = irb; }
    }



    /*
    // using the several bins around the maximum, minimum found!
    // --- |-|-|x|x|x|o|x|x|x|-|-|
    // o: found center
    // x: used for fitting.
    // use 1/10 nbinsy as half number of bins for fit
    */
    int np_side = y_nbin / 10;
    int np_used = 2 * np_side + 1;
    double *pos1 = new double[np_used], *pos2 = new double[np_used];
    double *temp1 = new double[np_used], *temp2 = new double[np_used];
    int npfill1 = 0, npfill2 =0;
    double fit_min1 = 99999., fit_max1 = -9999., fit_min2 = 99999., fit_max2 = -9999.;
    for(int ip=0; ip<np_used; ip++){
      int ip1 = ib_mmy_1 - np_side + ip;
      int ip2 = ib_mmy_2 - np_side + ip;
      if(ip1>=1){
        pos1[npfill1] = h1->GetXaxis()->GetBinCenter(ip1);
        temp1[npfill1] = h1->GetBinContent(ip1);
        if     ( pos1[npfill1] < fit_min1) fit_min1 = pos1[npfill1];
        else if( pos1[npfill1] > fit_max1) fit_max1 = pos1[npfill1];
        npfill1++;
      }
      if(ip2<=y_nbin){
        pos2[npfill2] = h1->GetXaxis()->GetBinCenter(ip2);
        temp2[npfill2] = h1->GetBinContent(ip2);
        if     ( pos2[npfill2] < fit_min2) fit_min2 = pos2[npfill2];
        else if( pos2[npfill2] > fit_max2) fit_max2 = pos2[npfill2];
        npfill2++;
      }
    }
    TGraph *gr1 = new TGraph(npfill1, pos1, temp1);
    TGraph *gr2 = new TGraph(npfill2, pos2, temp2);
    TF1 *tfs1 = NULL;
    TF1 *tfs2 = NULL;

    // results initialized with average method
    double res_mean1 = ib_mmy_1, res_width1 = -1., res_chisq1 = -1., res_temp1 = mmy_1;
    double res_mean2 = ib_mmy_2, res_widthframe = -1., res_chisq2 = -1., res_temp2 = mmy_2;
    int res_ndf1 = 2, res_ndf2 = 2;
    float frange_min1 = fit_min1 - binWidth;
    float frange_max1 = fit_max1 + binWidth;
    float frange_min2 = fit_min2 - binWidth;
    float frange_max2 = fit_max2 + binWidth;
    tfs1 = new TF1("SigFuncOne","gaus", frange_min1, frange_max1);
    tfs2 = new TF1("SigFuncTwo","gaus", frange_min2, frange_max2);

    if(negative){
      // p0: temperature, p0 < 0. negative temperature.
      // tfs1 -> SetParLimits(0, -100., -0.00001); 
      // tfs2 -> SetParLimits(0, -100., -0.00001);
      tfs1 -> SetParLimits(0, -100., 10.); 
      tfs2 -> SetParLimits(0, -100., 10.);
    }else{
      // p0 > 0. positive temperature.
      tfs1 -> SetParLimits(0, 0.00001, 100.);
      tfs2 -> SetParLimits(0, 0.00001, 100.);
      // tfs1 -> SetParLimits(0, 10., 100.);
      // tfs2 -> SetParLimits(0, 10., 100.);
    }

    // constraints on sigma
    tfs1 -> SetParameter(2, 2* np_side * binWidth); 
    tfs2 -> SetParameter(2, 2* np_side * binWidth); 
    tfs1 -> SetParLimits(2, np_side * binWidth, (np_side * 4) * binWidth ); 
    tfs2 -> SetParLimits(2, np_side * binWidth, (np_side * 4) * binWidth ); 

    // constraints on peak position
    tfs1 -> SetParLimits(1, frange_min1, frange_max1);
    tfs2 -> SetParLimits(1, frange_min2, frange_max2);
    tfs1 -> SetParameter(1, h1->GetXaxis()->GetBinCenter(ib_mmy_1) );
    tfs2 -> SetParameter(1, h1->GetXaxis()->GetBinCenter(ib_mmy_2) );
    //"0" without drawing
    gr1->Fit(tfs1,"0");
    gr2->Fit(tfs2,"0");
    tfs1->GetParameters(&pars[0]);
    tfs2->GetParameters(&pars[3]);

    // how large the chi2 needs to be checked?
    // here chi2 > 1.0
    if(tfs1->GetChisquare() > 1. && ix > 1){
      // check if width is larger than usual
      double sigma_avg = 0.;

      // average over the previous 10 bins, if the current one is over ix=10.
      // otherwise, just average all previous bins.
      int ixbins = 10; 
      if (ix <= ixbins) ixbins = ix - 1;
      for (int jx=(ix - ixbins); jx<=ix-1; jx++) sigma_avg += h_sigm_top -> GetBinContent(jx) / ixbins;
      // can do an iteration here.
      // if sigma is 40% over average
      if( pars[2] > sigma_avg * 1.4){
        tfs1 -> SetParLimits(2, np_side * binWidth, sigma_avg * 1.4);
        gr1->Fit(tfs1, "0");
        tfs1->GetParameters(&pars[0]);
      }
    }
    if(tfs2->GetChisquare() > 1. && ix > 1){
      double sigma_avg = 0.;
      int ixbins = 10; 
      if (ix <= ixbins) ixbins = ix - 1;
      for (int jx=(ix - ixbins); jx<=ix-1; jx++) sigma_avg += h_sigm_bot -> GetBinContent(jx) / ixbins;
      if( pars[5] > sigma_avg * 1.4){
        tfs2 -> SetParLimits(2, np_side * binWidth, sigma_avg * 1.4);
        gr2->Fit(tfs2, "0");
        tfs2->GetParameters(&pars[3]);
      }
    }

    res_temp1 = pars[0];
    res_mean1 = pars[1];
    res_width1 = pars[2];
    res_chisq1 = tfs1->GetChisquare();
    res_ndf1 = tfs1->GetNDF();

    res_temp2 = pars[3];
    res_mean2 = pars[4];
    res_widthframe = pars[5];
    res_chisq2 = tfs2->GetChisquare();
    res_ndf2 = tfs2->GetNDF();

    tfs1->SetLineWidth(2);
    tfs2->SetLineWidth(2);
    tfs1->SetLineColor(kBlue);
    tfs2->SetLineColor(kRed);
    tfs1->Draw("same");
    tfs2->Draw("same");

    //c1 -> Update();

    gr1 ->SetLineColor(4);
    gr1 ->SetMarkerColor(4);
    gr1 ->SetMarkerStyle(20);
    gr2 ->SetLineColor(2);
    gr2 ->SetMarkerColor(2);
    gr2 ->SetMarkerStyle(21);
    gr1 -> Draw("samePC");
    gr2 -> Draw("samePC");

    // color: 2 = red, 4 = blue
    //Draw in NDC (Normalized Device Coordinates) [0,1]
    TLatex Tl; Tl.SetTextSize(20); Tl.SetTextFont(43);
    Tl.SetNDC(); 
    std::string ss = Form("#color[4]{Left: #chi^{2}/NDF = %.1f / %d}", res_chisq1, res_ndf1);
    double x_lef = 0.18;
    double x_rig = 0.58;
    double y_txt = 0.88, y_gap = 0.04;
    Tl.DrawLatex(x_lef, y_txt , ss.c_str());
    ss = Form("#color[4]{mean: %5.3f}", res_mean1);
    Tl.DrawLatex(x_lef, y_txt - y_gap, ss.c_str());
    ss = Form("#color[4]{width: %5.3f}", res_width1);
    Tl.DrawLatex(x_lef, y_txt - y_gap * 2, ss.c_str());
    ss = Form("#color[4]{T: %3.1f #circ C}", res_temp1);
    Tl.DrawLatex(x_lef, y_txt - y_gap * 3, ss.c_str());

    ss = Form("#color[2]{Right: #chi^{2}/NDF = %.1f / %d}", res_chisq2, res_ndf2);
    Tl.DrawLatex(x_rig, y_txt, ss.c_str());
    ss = Form("#color[2]{mean: %5.3f}", res_mean2);
    Tl.DrawLatex(x_rig, y_txt - y_gap, ss.c_str());
    ss = Form("#color[2]{width: %5.3f}", res_widthframe);
    Tl.DrawLatex(x_rig, y_txt - y_gap * 2, ss.c_str());
    ss = Form("#color[2]{T: %3.1f #circ C}", res_temp2);
    Tl.DrawLatex(x_rig, y_txt - y_gap * 3, ss.c_str());


    ss = Form("#color[1]{side:%s, T_{set}=%3.1f}", s_side.c_str(), Tliquid);
    Tl.DrawLatex(x_lef, y_txt - y_gap * 4.5, ss.c_str());
    double cmXpos = hframe -> GetXaxis() -> GetBinCenter(ix);
    ss = Form("#color[1]{X pixel=%3d (@%3.1fcm)}", ix, cmXpos);
    Tl.DrawLatex(x_lef, y_txt - y_gap * 5.8, ss.c_str());


    // fit results
    h_temp_bot -> SetBinContent(ix, res_temp1);
    h_posi_bot -> SetBinContent(ix, res_mean1);
    h_sigm_bot -> SetBinContent(ix, res_width1);
    h_chi2_bot -> SetBinContent(ix, res_chisq1 / res_ndf1);
    h_temp_top -> SetBinContent(ix, res_temp2);
    h_posi_top -> SetBinContent(ix, res_mean2);
    h_sigm_top -> SetBinContent(ix, res_widthframe);
    h_chi2_top -> SetBinContent(ix, res_chisq2 / res_ndf2);
    std::string hname = "fit1D_";
    hname += Form("x%d",ix);
    c1->Print((hname+".png").c_str());

    delete h1; 
    delete c1;
    delete pos1;
    delete pos2;
    delete temp1; 
    delete temp2;
    delete gr1;
    delete gr2;
    if(tfs1) delete tfs1;
    if(tfs2) delete tfs2;
  }

  h_temp_top->Write();
  h_posi_top->Write();
  h_sigm_top->Write();
  h_chi2_top->Write();
  h_temp_bot->Write();
  h_posi_bot->Write();
  h_sigm_bot->Write();
  h_chi2_bot->Write();


  // combining _top and _bot into one single line, end of top should connect to end of bottom
  // and bottom head is the new end of the combined line.
  TH1F* h_temp = new TH1F((s_out+"_temperature").c_str(), "", 2 * nbinsx, hxmin, hxmax * 2 - hxmin); 
  TH1F* h_posi = new TH1F((s_out+"_peakposition").c_str(), "", 2 * nbinsx, hxmin, hxmax * 2 - hxmin);
  TH1F* h_sigm = new TH1F((s_out+"_peakwidth").c_str(), "", 2 * nbinsx, hxmin, hxmax * 2 - hxmin); 
  TH1F* h_chi2 = new TH1F((s_out+"_chi2overndf").c_str(), "", 2 * nbinsx, hxmin, hxmax * 2 - hxmin);
  for(int ib=1; ib<=nbinsx; ib++){
    h_temp -> SetBinContent(ib, h_temp_top -> GetBinContent(ib) );
    h_posi -> SetBinContent(ib, h_posi_top -> GetBinContent(ib) );
    h_sigm -> SetBinContent(ib, h_sigm_top -> GetBinContent(ib) );
    h_chi2 -> SetBinContent(ib, h_chi2_top -> GetBinContent(ib) );

    h_temp -> SetBinContent(2 * nbinsx + 1 - ib, h_temp_bot -> GetBinContent(ib) );
    h_posi -> SetBinContent(2 * nbinsx + 1 - ib, h_posi_bot -> GetBinContent(ib) );
    h_sigm -> SetBinContent(2 * nbinsx + 1 - ib, h_sigm_bot -> GetBinContent(ib) );
    h_chi2 -> SetBinContent(2 * nbinsx + 1 - ib, h_chi2_bot -> GetBinContent(ib) );
    // artificially reduce the difference in the connection part to 1/3 of the initial difference !!!
    if(ib==nbinsx){
      double Tdif = h_temp_top -> GetBinContent(ib) - h_temp_bot -> GetBinContent(ib);
      h_temp -> SetBinContent(ib, h_temp_top -> GetBinContent(ib) - Tdif / 3. );
      h_temp -> SetBinContent(2 * nbinsx + 1 - ib, h_temp_bot -> GetBinContent(ib) + Tdif / 3.);
    }

  }

  h_temp->Write();
  h_posi->Write();
  h_sigm->Write();
  h_chi2->Write();


  delete h_temp_top;
  delete h_posi_top;
  delete h_sigm_top;
  delete h_chi2_top;
  delete h_temp_bot;
  delete h_posi_bot;
  delete h_sigm_bot;
  delete h_chi2_bot;
  delete h_temp;
  delete h_posi;
  delete h_sigm;
  delete h_chi2;

  fout->Close();
  if(fout) delete fout;

  return true;
}

bool HTools::getPeaks(TH1F* hp, std::string s_out, bool negative, double sigma, double threshold, int niters){
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
  /*
  fPad0 -> SetGrid();
  fPad1 -> SetGrid();
  */
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
  hp_sm -> SetMaximum( Tmax + 2. );
  hp_sm -> SetMinimum( Tmin );
  hp_sm -> Draw();
  hp -> Draw("same");


  double xmin = hp -> GetXaxis()->GetXmin();
  double xmax = hp -> GetXaxis()->GetXmax();
  int npeaks = 20; // guess how many peaks on the curve.
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
  hpbkg -> SetLineColor(kGray+3);
  hpbkg -> Draw("same");
  hpbkg -> Fit("pol1");
  TF1 * fl1 = hpbkg -> GetFunction("pol1");
  double hpb0 = fl1 -> GetParameter(0);
  double hpb1 = fl1 -> GetParameter(1);
  double hpb0_err = fl1 -> GetParError(0);
  double hpb1_err = fl1 -> GetParError(1);


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
  fTxt.SetTextSize(0.15); //0.04
  fTxt.SetNDC();
  fTxt.DrawLatex(0.3, 0.85, Form("Bkg NIter=%d, peak sigma=%4.2f, threshold=%4.2f",niters,sigma,threshold));
  fTxt.DrawLatex(0.3, 0.62, Form("#color[9]{Bkg pol1, p0 = %6.4f #pm %5.4f, p1 = %6.4f #pm %5.4f.}", hpb0, hpb0_err, hpb1, hpb1_err));

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
  
  std::string pname = "peaks_"; // + s_out + "_";
  pname += hname;
  c0->Print( (pname+".png").c_str());
  c0->Print( (pname+".pdf").c_str());

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

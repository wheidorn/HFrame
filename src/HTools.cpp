#define HTOOLSS_CPP
#include "HTools.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TROOT.h"
#include "TStyle.h"
#include <math.h>

#include "TF1.h"
#include "func.h"


bool HTools::fitPipeTwoGaus(HFrame * fbox, std::string s_out, int x_min, int x_max, int y_min, int y_max, bool negative, bool convert_to_cm){
  if(!fbox) return false;
  TH2F* h2 = fbox->getFrame();
  if(!h2) return false;

  gStyle->SetOptStat(0);

  int nbinx = h2->GetNbinsX();
  int nbiny = h2->GetNbinsY();
  if(x_min <=0 || x_max >nbinx) return false;
  if(y_min <=0 || y_max >nbiny) return false;

  int nbinx_fit = x_max - x_min + 1;
  float hxmin = x_min, hxmax = x_max + 1;
  if(convert_to_cm){
    hxmin = x_min * cmperpix;
    hxmax = (x_max + 1) * cmperpix;
  }
  TH1F* h_temp_top = new TH1F((s_out+"_temperature_top").c_str(), "", nbinx_fit, hxmin, hxmax); // Temperature peak value
  TH1F* h_posi_top = new TH1F((s_out+"_peakposition_top").c_str(), "", nbinx_fit, hxmin, hxmax); // peak pixel position
  TH1F* h_sigm_top = new TH1F((s_out+"_gaussigma_top").c_str(), "", nbinx_fit, hxmin, hxmax); // sigma out of fit, gaussian
  TH1F* h_chi2_top = new TH1F((s_out+"_chi2overndf_top").c_str(), "", nbinx_fit, hxmin, hxmax); // chi2 / NDF out of fit
  TH1F* h_temp_bot = new TH1F((s_out+"_temperature_bot").c_str(), "", nbinx_fit, hxmin, hxmax);
  TH1F* h_posi_bot = new TH1F((s_out+"_peakposition_bot").c_str(), "", nbinx_fit, hxmin, hxmax);
  TH1F* h_sigm_bot = new TH1F((s_out+"_gaussigma_bot").c_str(), "", nbinx_fit, hxmin, hxmax);
  TH1F* h_chi2_bot = new TH1F((s_out+"_chi2overndf_bot").c_str(), "", nbinx_fit, hxmin, hxmax); 

  const int npars=6;
  double pars[npars];
  for(int ix=x_min; ix<=x_max; ix++){


    // h0 histogram with pixels.
    //TH1F* h0 = fbox->getLine(ix, 1, ix, nbiny);
    TH1F* h0 = fbox->getLine(ix, y_min, ix, y_max);
    TH1F* h1 = (TH1F*) h0 -> Clone( h0->GetName() + char('X') );
    if (!h1){ std::cout<<"error finding the line in fit "<<std::endl; return false;}

     // find the x position of Y min or max
     // then set up the fit range around the min or max
    int ib_mmy_1 = 0., ib_mmy_2 = 0.; // mm = max or min
    double mmy_1 = -99., mmy_2 = -99.;
    if(negative){mmy_1 = 9999.; mmy_2 = 9999.;}
    int y_thr = h1->GetNbinsX() / 3;

    std::cout<<" bins y for fit "<<h1->GetNbinsX()<<std::endl;

    for(int ib=1; ib<=h1->GetNbinsX(); ib++){
      if (ib<= y_thr){
        if      ( negative && mmy_1 > h1->GetBinContent(ib)){ mmy_1 = h1->GetBinContent(ib); ib_mmy_1 = ib; }
        else if (!negative && mmy_1 < h1->GetBinContent(ib)){ mmy_1 = h1->GetBinContent(ib); ib_mmy_1 = ib; }
      }else if (ib >= h1->GetNbinsX() - y_thr){
        if      ( negative && mmy_2 > h1->GetBinContent(ib)){ mmy_2 = h1->GetBinContent(ib); ib_mmy_2 = ib; }
        else if (!negative && mmy_2 < h1->GetBinContent(ib)){ mmy_2 = h1->GetBinContent(ib); ib_mmy_2 = ib; }
      }
    }


    // using the 7 bins around the maximum, minimum found!
    // --- |-|-|x|x|x|o|x|x|x|-|-|
    // o: found center
    // x: used for fitting.
    TF1 *tfs1 = new TF1("SigFuncOne","gaus", (ib_mmy_1 >= 4 ? float(ib_mmy_1 - 4) : 0.), float( ib_mmy_1 + 3));
    TF1 *tfs2 = new TF1("SigFuncTwo","gaus", float( ib_mmy_2 - 4.), (ib_mmy_2 + 3 <= nbiny ? float(ib_mmy_2 + 3) : float(nbiny)) );

    TCanvas *c1 = new TCanvas("c1","c1",10,10,600,600);
    h1->SetMaximum( h1->GetMaximum() * 1.2);
    if(negative) h1->SetMaximum( h1->GetMaximum() * 0.75);
    //h1->SetMinimum( h0->GetMinimum());
    h1->SetLineWidth(3);
    h1->SetLineColor(1);
    h1->SetMarkerStyle(23);
    h1->SetMarkerColor(1);
    h1->Draw();

    if(negative){
      tfs1 -> SetParLimits(0, -100., -0.00001); // p0 < 0. negative temperature.
      tfs2 -> SetParLimits(0, -100., -0.00001);
      tfs1 -> SetParLimits(2, 0.1, 100.); // sigma
      tfs2 -> SetParLimits(2, 0.1, 100.);
      // peak position
      //tfs1 -> SetParLimits(1, y_range / 4., y_range * 3. / 4.);
      //tfs2 -> SetParLimits(1, nbiny - y_range * 3. / 4., nbiny - y_range / 4.);
    }
    else{
      tfs1 -> SetParLimits(0, 0.00001, 100.); // p0 < 0. negative temperature.
      tfs2 -> SetParLimits(0, 0.00001, 100.);
      tfs1 -> SetParLimits(2, 0.1, 100.); // sigma
      tfs2 -> SetParLimits(2, 0.1, 100.);
    }


    h1->Fit(tfs1,"R+");
    h1->Fit(tfs2,"R+");
    tfs1->GetParameters(&pars[0]);
    tfs2->GetParameters(&pars[3]);

    // improve the picture:
    tfs1->SetLineWidth(2);
    tfs2->SetLineWidth(2);
    tfs1->SetLineColor(kBlue);
    tfs2->SetLineColor(kRed);

    // writes the fit results into the par array
    tfs1->SetParameters(pars);
    tfs1->Draw("same");
    tfs2->SetParameters(&pars[3]);
    tfs2->Draw("same");


    // color: 2 = red, 4 = blue
    TLatex Tl; Tl.SetTextSize(22); Tl.SetTextFont(43);
    std::string ss = Form("#color[4]{Left: #chi^{2}/NDF = %.1f / %d}", tfs1->GetChisquare(), tfs1->GetNDF());
    double x_txt = 0.6*h1->GetXaxis()->GetXmin() + 0.4*h1->GetXaxis()->GetXmax();
    double y_txt = 0.2*h1->GetMinimum() + 0.8*h1->GetMaximum();
    Tl.DrawLatex(x_txt, y_txt , ss.c_str());
    ss = Form("#color[2]{Right: #chi^{2}/NDF = %.1f / %d}", tfs2->GetChisquare(), tfs2->GetNDF());
    y_txt = 0.3*h1->GetMinimum() + 0.7*h1->GetMaximum();
    Tl.DrawLatex(x_txt, y_txt , ss.c_str());
      
    // fit results
    // check the definition of the TH1 to decide which bin to set
    h_temp_bot -> SetBinContent(ix - x_min+1, pars[0]);
    h_posi_bot -> SetBinContent(ix - x_min+1, pars[1]);
    h_sigm_bot -> SetBinContent(ix - x_min+1, pars[2]);
    h_chi2_bot -> SetBinContent(ix - x_min+1, tfs1->GetChisquare() / tfs1->GetNDF());
    h_temp_top -> SetBinContent(ix - x_min+1, pars[3]);
    h_posi_top -> SetBinContent(ix - x_min+1, pars[4]);
    h_sigm_top -> SetBinContent(ix - x_min+1, pars[5]);
    h_chi2_top -> SetBinContent(ix - x_min+1, tfs2->GetChisquare() / tfs2->GetNDF());
    std::string hname = "fit1D_";
    hname += Form("x%d",ix);
    c1->Print((hname+".png").c_str());
    delete c1;
    delete tfs1;
    delete tfs2;
  }
  fbox -> recordLine(h_temp_top);
  fbox -> recordLine(h_posi_top);
  fbox -> recordLine(h_sigm_top);
  fbox -> recordLine(h_chi2_top);
  fbox -> recordLine(h_temp_bot);
  fbox -> recordLine(h_posi_bot);
  fbox -> recordLine(h_sigm_bot);
  fbox -> recordLine(h_chi2_bot);
  return true;
}


TH1F* HTools::averagePipe(HFrame * fbox, std::string s_out, int x_min, int x_max, int y1, int y2, bool negative, bool convert_to_cm){
  if(!fbox) return NULL;
  TH2F* h2 = fbox->getFrame();
  if(!h2) return NULL;
  int nbinx = h2->GetNbinsX();
  int nbiny = h2->GetNbinsY();
  if (y1 > y2 || y2 > nbiny) return NULL;
  if (x_min > x_max || x_max > nbinx) return NULL;

  int navgbin = y2 - y1 + 1, nbinx_out = x_max - x_min + 1;

  float lmin = x_min, lmax = x_max+1;
  if(convert_to_cm){lmin = x_min * cmperpix; lmax = (x_max + 1) * cmperpix;}

  TH1F* h0 = new TH1F(s_out.c_str(), "", nbinx_out, lmin, lmax);
  for(int ix=x_min; ix<=x_max; ix++){
    //for (int iy=y1; iy<=y2; iy++){
    //  avg += h2->GetBinContent(ix,iy) / navgbin;
    //}
    //
    // NOT average, but min or max
    double avg = -999.; if (negative) avg = 999.;
    for (int iy=y1; iy<=y2; iy++){
      double cur = h2->GetBinContent(ix,iy);
      if( negative && avg > cur) avg = cur;
      if(!negative && avg < cur) avg = cur;
    }

    h0->SetBinContent(ix - x_min + 1, avg);
    //delete h1;
  }
  //std::string sfilename = "b1_fit.root";
  //fbox->writeToRootFile(sfilename);
  fbox -> recordLine(h0);

  gStyle->SetOptStat(0);
  TCanvas *c1 = new TCanvas("c1","c1",10,10,600,600);
  h0->SetMaximum( h0->GetMaximum() * 1.1);
  if( h0->GetMaximum() < 0.) h0->SetMaximum( h0->GetMaximum() * 0.8);
  h0->SetLineWidth(3);
  h0->SetLineColor(1);
  h0->SetMarkerStyle(23);
  h0->SetMarkerColor(1);
  h0->SetYTitle("Temperature (#circ C)");
  if(convert_to_cm)
  h0->SetXTitle("X direction (cm)");
  else 
  h0->SetXTitle("X direction pixel");

  h0->Draw();
  c1->Print((s_out+".png").c_str());


  delete c1;
  return h0;

}


HFrame* HTools::averageFrames(std::vector<HFrame *> fboxes, std::string shout){
  int nbox = fboxes.size();
  if(nbox<=0) return NULL;
  else if(nbox==1) return fboxes.at(0);

  TH2F* h2 = NULL;
  for ( auto fbox : fboxes){
    TH2F * hbox = fbox->getFrame();
    float scl = 1. / nbox;
    if (!h2){ h2 = (TH2F*) hbox -> Clone(shout.c_str()); h2->Scale( scl );}
    else h2->Add( hbox, scl ); 
  }
  HFrame * fboxout = new HFrame();
  fboxout -> setFrame(h2);
  std::cout<<" average "<<h2<<" frame "<<fboxout<<" getframe "<<fboxout->getFrame()<<std::endl;
  return fboxout;
}

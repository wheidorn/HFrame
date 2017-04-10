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

bool HTools::findFramePipes(HFrame * fbox, std::string s_out, int method_id){
  if (method_id <0 || method_id >2){
    std::cout<<"method_id = 0: fit gaussian about the T peaks!"<<std::endl;
    std::cout<<"method_id = 1: smooth the input Histogram to find T peaks!"<<std::endl;
    std::cout<<"method_id = 2: average over the T peaks!"<<std::endl;
    return false;
  }

  if(!fbox) return false;
  TH2F* h2 = fbox->getFramePipe();
  if(!h2){
    std::cout<<"ERROR: frame for pipe lines NOT found! "<<std::endl;
    return false;
  }
  bool negative = fbox->isNegative();

  int nbinsx = h2 -> GetNbinsX();
  int nbinsy = h2 -> GetNbinsY();
  float hxmin = h2 -> GetXaxis() -> GetXmin();
  float hxmax = h2 -> GetXaxis() -> GetXmax();

  gStyle->SetOptStat(0);
  TH1F * h_tmp = (TH1F*) h2 -> ProjectionX() -> Clone("tmp"); h_tmp -> Reset();
  // Temperature peak value
  // peak pixel position
  // sigma out of fit, gaussian
  // chi2 / NDF out of fit
  TH1F* h_temp_top = (TH1F*) h_tmp -> Clone((s_out+"_temperature_top").c_str());
  TH1F* h_posi_top = (TH1F*) h_tmp -> Clone((s_out+"_peakposition_top").c_str());
  TH1F* h_sigm_top = (TH1F*) h_tmp -> Clone((s_out+"_peakwidth_top").c_str());
  TH1F* h_chi2_top = (TH1F*) h_tmp -> Clone((s_out+"_chi2overndf_top").c_str());
  TH1F* h_temp_bot = (TH1F*) h_tmp -> Clone((s_out+"_temperature_bot").c_str());
  TH1F* h_posi_bot = (TH1F*) h_tmp -> Clone((s_out+"_peakposition_bot").c_str());
  TH1F* h_sigm_bot = (TH1F*) h_tmp -> Clone((s_out+"_peakwidth_bot").c_str());
  TH1F* h_chi2_bot = (TH1F*) h_tmp -> Clone((s_out+"_chi2overndf_bot").c_str());
  delete h_tmp;

  const int npars=6;
  double pars[npars];
  for(int ix=1; ix<=nbinsx; ix++){

    TH1F* h1 = fbox->getLine(ix, 1, ix, nbinsy, true); // true: read from frame pipe
    TH1F* h1_sm = (TH1F*) h1 -> Clone( h1 -> GetName() + char('S') ); // smoothing
    h1_sm -> Smooth();
    if (!h1){ std::cout<<"error finding the line in fit "<<std::endl; return false;}

    TCanvas *c1 = new TCanvas("c1","c1",600,600);
    //SetMargin (Float_t left, Float_t right, Float_t bottom, Float_t top)
    c1 -> SetMargin(0.12,0.03,0.12,0.03);
    double Maxi = int(h1->GetMaximum() ) + 1.; 
    double Mini = int(h1->GetMinimum() ) - 1.; 
    h1->SetMaximum( Maxi+6.);
    h1->SetMinimum( Mini);
    h1->SetLineWidth(3);
    h1->SetLineColor(1);
    h1->SetMarkerStyle(23);
    h1->SetMarkerColor(1);
    h1->GetXaxis() -> SetTitle("Y pixel index");
    h1->GetYaxis() -> SetTitle("Temperature (#circ C)");
    h1->GetXaxis() -> SetTitleOffset(1.2* h1->GetXaxis() -> GetTitleOffset());
    h1->GetYaxis() -> SetTitleOffset(1.3* h1->GetYaxis() -> GetTitleOffset());
    h1->Draw();


    // find the x position of Y min or max
    // then set up the fit range around the min or max
    int ib_mmy_1 = 0., ib_mmy_2 = 0.; // mm = max or min
    double mmy_1 = -99., mmy_2 = -99.;
    if(negative){mmy_1 = 9999.; mmy_2 = 9999.;}
    int y_nbin = h1->GetNbinsX();
    int y_thr = h1->GetNbinsX() / 3;

    // find minimum by averaging of 3 pixels
    for(int ib=2; ib<h1->GetNbinsX(); ib++){
      double avg3 = (h1->GetBinContent(ib-1) + h1->GetBinContent(ib) + h1->GetBinContent(ib+1) )/ 3.;
      if (ib<= y_thr){
        if      ( negative && mmy_1 > avg3){ mmy_1 = avg3; ib_mmy_1 = ib; }
        else if (!negative && mmy_1 < avg3){ mmy_1 = avg3; ib_mmy_1 = ib; }
      }else if (ib >= h1->GetNbinsX() - y_thr){
        if      ( negative && mmy_2 > avg3){ mmy_2 = avg3; ib_mmy_2 = ib; }
        else if (!negative && mmy_2 < avg3){ mmy_2 = avg3; ib_mmy_2 = ib; }
      }
    }


    // using the several bins around the maximum, minimum found!
    // --- |-|-|x|x|x|o|x|x|x|-|-|
    // o: found center
    // x: used for fitting.
    int np_side = y_nbin / 10; // use 1/10 y_nbin as half number of bins for fit
    int np_used = 2 * np_side + 1;
    double *pos1 = new double[np_used], *pos2 = new double[np_used];
    double *temp1 = new double[np_used], *temp2 = new double[np_used];
    int npfill1 = 0, npfill2 =0;
    for(int ip=0; ip<np_used; ip++){
      int ip1 = ib_mmy_1 - 3 + ip;
      int ip2 = ib_mmy_2 - 3 + ip;
      if(ip1>=1){
        if (method_id == 1){
          pos1[npfill1] = h1_sm->GetXaxis()->GetBinCenter(ip1);
          temp1[npfill1] = h1_sm->GetBinContent(ip1);
        }else{
          pos1[npfill1] = h1->GetXaxis()->GetBinCenter(ip1);
          temp1[npfill1] = h1->GetBinContent(ip1);
        }
        npfill1++;
      }
      if(ip2<=y_nbin){
        if (method_id == 1){
          pos2[npfill2] = h1_sm->GetXaxis()->GetBinCenter(ip2);
          temp2[npfill2] = h1_sm->GetBinContent(ip2);
        }else{
          pos2[npfill2] = h1->GetXaxis()->GetBinCenter(ip2);
          temp2[npfill2] = h1->GetBinContent(ip2);
        }
        npfill2++;
      }
    }
    TGraph *gr1 = new TGraph(npfill1, pos1, temp1);
    TGraph *gr2 = new TGraph(npfill2, pos2, temp2);
    TF1 *tfs1 = NULL;
    TF1 *tfs2 = NULL;



    // results initialized with average method
    double res_mean1 = ib_mmy_1, res_width1 = -1., res_chisq1 = -1., res_temp1 = mmy_1;
    double res_mean2 = ib_mmy_2, res_width2 = -1., res_chisq2 = -1., res_temp2 = mmy_2;
    int res_ndf1 = 2, res_ndf2 = 2;


    if(method_id ==0){
      // fitting gaussian to get peak and tempearatures
      float frange_min1 = (ib_mmy_1 > np_side ? ib_mmy_1 - np_side - 0.5 : 0.5);
      float frange_max1 = ib_mmy_1 + np_side + 0.5;
      float frange_min2 = ib_mmy_2 - np_side - 0.5;
      float frange_max2 = (ib_mmy_2 + np_side < nbinsy ? ib_mmy_2 + np_side + 0.5 : nbinsy - 0.5);
      tfs1 = new TF1("SigFuncOne","gaus", frange_min1, frange_max1);
      tfs2 = new TF1("SigFuncTwo","gaus", frange_min2, frange_max2);

      if(negative){
        tfs1 -> SetParLimits(0, -100., -0.00001); // p0 < 0. negative temperature.
        tfs2 -> SetParLimits(0, -100., -0.00001);
      }
      else{
        tfs1 -> SetParLimits(0, 0.00001, 100.); // p0 > 0. positive temperature.
        tfs2 -> SetParLimits(0, 0.00001, 100.);
      }
      tfs1 -> SetParLimits(2, y_thr / 4., y_thr * 3.); // constraints on sigma
      tfs2 -> SetParLimits(2, y_thr / 4., y_thr * 3.);
      tfs1 -> SetParLimits(1, frange_min1, frange_max1); // constraints on sigma
      tfs2 -> SetParLimits(1, frange_min2, frange_max2);
      tfs1 -> SetParameter(1, h1->GetXaxis()->GetBinCenter(ib_mmy_1) );
      tfs2 -> SetParameter(1, h1->GetXaxis()->GetBinCenter(ib_mmy_2) );
      gr1->Fit(tfs1,"0");
      gr2->Fit(tfs2,"0"); //without drawing
      tfs1->GetParameters(&pars[0]);
      tfs2->GetParameters(&pars[3]);
      // how large the chi2 needs to be checked?
      // here chi2 > 1.0
      if(tfs1->GetChisquare() > 1.){
        // check if width is larger than usual
        double sigma_avg = 0.;
        for (int jx=1; jx<=ix-1; jx++) sigma_avg += h_sigm_top -> GetBinContent(jx) / (ix-1);
        // can do an iteration here.
        // if sigma is 40% over average
        if( pars[2] > sigma_avg * 1.4){
          tfs1 -> SetParLimits(2, y_thr / 4., sigma_avg * 1.4);
          gr1->Fit(tfs1, "0");
          tfs1->GetParameters(&pars[0]);
        }
      }
      if(tfs2->GetChisquare() > 1.){
        double sigma_avg = 0.;
        for (int jx=1; jx<=ix-1; jx++) sigma_avg += h_sigm_bot -> GetBinContent(jx) / (ix-1);
        if( pars[5] > sigma_avg * 1.4){
          tfs2 -> SetParLimits(2, y_thr / 4., sigma_avg * 1.4);
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
      res_width2 = pars[5];
      res_chisq2 = tfs2->GetChisquare();
      res_ndf2 = tfs2->GetNDF();

      tfs1->SetLineWidth(2);
      tfs2->SetLineWidth(2);
      tfs1->SetLineColor(kBlue);
      tfs2->SetLineColor(kRed);
      tfs1->Draw("same");
      tfs2->Draw("same");
 
      //c1 -> Update();
    }else if (method_id ==1 ){
      res_chisq1 = 0.; 
      res_chisq2 = 0.;
      // searching for the point X where the T value is half of the max /min;
      double half_a1 = 999, half_b1 = 999., half_a2 = 999., half_b2 = 999.;
      double Tdif1 = (temp1[0] + temp1[npfill1 - 1]) / 2 - h1_sm -> GetBinContent(ib_mmy_1);
      double Tdif2 = (temp2[0] + temp2[npfill2 - 1]) / 2 - h1_sm -> GetBinContent(ib_mmy_2);
      int ix_a1 = -1, ix_b1 = -1, ix_a2 = -1, ix_b2 = -1;
      for(int ip=0; ip<npfill1; ip++){
        float posx = pos1[ip]; 
        int ipos = h1 -> GetXaxis() -> FindBin( posx );
        float val = h1 -> GetBinContent(ipos); // raw temparature before smoothing
        res_chisq1 += pow((val - temp1[ip]), 2) / fabs(temp1[ip]);

        if      ( ipos <= ib_mmy_1 && fabs(temp1[ip] - Tdif1) < half_a1 ){ix_a1 = ipos; half_a1 = fabs(temp1[ip] - Tdif1); }
        else if ( ipos >= ib_mmy_1 && fabs(temp1[ip] - Tdif1) < half_b1 ){ix_b1 = ipos; half_b1 = fabs(temp1[ip] - Tdif1); }
      }
      for(int ip=0; ip<npfill2; ip++){
        float posx = pos2[ip];
        int ipos = h1 -> GetXaxis() -> FindBin( posx );
        float val = h1 -> GetBinContent(ipos); // raw temparature before smoothing
        res_chisq2 += pow((val - temp2[ip]), 2) / fabs(temp2[ip]);

        if      ( ipos <= ib_mmy_2 && fabs(temp2[ip] - Tdif2) < half_a2 ){ix_a2 = ipos; half_a2 = fabs(temp2[ip] - Tdif2); }
        else if ( ipos >= ib_mmy_2 && fabs(temp2[ip] - Tdif2) < half_b2 ){ix_b2 = ipos; half_b2 = fabs(temp2[ip] - Tdif2); }
      }

      res_temp1 = mmy_1;
      res_mean1 = h1 -> GetXaxis() -> GetBinCenter(ib_mmy_1);
      res_width1 = ix_b1 - ix_a1;
      res_ndf1 = npfill1;

      res_temp2 = mmy_2;
      res_mean2 = h1 -> GetXaxis() -> GetBinCenter(ib_mmy_2);
      res_width2 = ix_b2 - ix_a2;
      res_ndf2 = npfill2;
    }

    gr1 ->SetLineColor(4);
    gr1 ->SetMarkerColor(4);
    gr1 ->SetMarkerStyle(20);
    gr2 ->SetLineColor(2);
    gr2 ->SetMarkerColor(2);
    gr2 ->SetMarkerStyle(21);
    gr1 -> Draw("samePC");
    gr2 -> Draw("samePC");
    //if(method_id ==0){
    //  // improve the picture:
    //  tfs1->SetLineWidth(2);
    //  tfs2->SetLineWidth(2);
    //  tfs1->SetLineColor(kBlue);
    //  tfs2->SetLineColor(kRed);
    //  tfs1->Draw("same");
    //  tfs2->Draw("same");
    //}


    // color: 2 = red, 4 = blue
    TLatex Tl; Tl.SetTextSize(20); Tl.SetTextFont(43);
    Tl.SetNDC(); //Draw in NDC (Normalized Device Coordinates) [0,1]
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
    ss = Form("#color[2]{width: %5.3f}", res_width2);
    Tl.DrawLatex(x_rig, y_txt - y_gap * 2, ss.c_str());
    ss = Form("#color[2]{T: %3.1f #circ C}", res_temp2);
    Tl.DrawLatex(x_rig, y_txt - y_gap * 3, ss.c_str());


    ss = Form("#color[1]{side:%s, T_{set}=%3.1f}", (fbox -> isSideL() ? "L" : "J"), fbox -> getTLiquid());
    Tl.DrawLatex(x_lef, y_txt - y_gap * 4.5, ss.c_str());
    double cmperpix = fbox -> getCMperPixel();
    ss = Form("#color[1]{X pixel=%3d (@%3.1fcm)}", ix, (ix - 1) * cmperpix);
    Tl.DrawLatex(x_lef, y_txt - y_gap * 5.8, ss.c_str());


    // fit results
    h_temp_bot -> SetBinContent(ix, res_temp1);
    h_posi_bot -> SetBinContent(ix, res_mean1);
    h_sigm_bot -> SetBinContent(ix, res_width1);
    h_chi2_bot -> SetBinContent(ix, res_chisq1 / res_ndf1);
    h_temp_top -> SetBinContent(ix, res_temp2);
    h_posi_top -> SetBinContent(ix, res_mean2);
    h_sigm_top -> SetBinContent(ix, res_width2);
    h_chi2_top -> SetBinContent(ix, res_chisq2 / res_ndf2);
    std::string hname = "hist1D_";
    if(method_id ==0 ) hname = "fit1D_";
    hname += Form("x%d",ix);
    c1->Print((hname+".png").c_str());

    delete h1; 
    delete h1_sm; 
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

  fbox -> recordPipeInfo(std::string(h_temp_top->GetName()), h_temp_top);
  fbox -> recordPipeInfo(std::string(h_posi_top->GetName()), h_posi_top);
  fbox -> recordPipeInfo(std::string(h_sigm_top->GetName()), h_sigm_top);
  fbox -> recordPipeInfo(std::string(h_chi2_top->GetName()), h_chi2_top);
  fbox -> recordPipeInfo(std::string(h_temp_bot->GetName()), h_temp_bot);
  fbox -> recordPipeInfo(std::string(h_posi_bot->GetName()), h_posi_bot);
  fbox -> recordPipeInfo(std::string(h_sigm_bot->GetName()), h_sigm_bot);
  fbox -> recordPipeInfo(std::string(h_chi2_bot->GetName()), h_chi2_bot);


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

  fbox -> recordPipeInfo(std::string(h_temp->GetName()), h_temp);
  fbox -> recordPipeInfo(std::string(h_posi->GetName()), h_posi);
  fbox -> recordPipeInfo(std::string(h_sigm->GetName()), h_sigm);
  fbox -> recordPipeInfo(std::string(h_chi2->GetName()), h_chi2);
  return true;
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



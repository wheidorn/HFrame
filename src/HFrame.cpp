#define HFRAME_CPP
#include "HFrame.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TTree.h"
#include <math.h>       /* std::abs */


HFrame::HFrame(std::string sfile, std::string sname, bool convert_to_cm){
  crop_point1 = std::make_pair(-1, -1);
  crop_point2 = std::make_pair(-1, -1);
  pipe_point1 = std::make_pair(-1, -1);
  pipe_point2 = std::make_pair(-1, -1);

  conversionToCM = convert_to_cm;
  cmperpix = CMPERPIXELPRE;
  lens = 80;
  if(sfile.find("sJ") != std::string::npos || sfile.find("_J") != std::string::npos){
    side_L = false;
  }else if(sfile.find("sL") != std::string::npos || sfile.find("_L") != std::string::npos){
    side_L = true;
  }else{
    std::cout<<"WARNING: Frame side information not available in file name: "<<sfile<<". Looking for sJ or _J or sL or _L. Set to side L."<<std::endl;
    std::cout<<"Use function p_to_frame -> setSideL( ... ) to modify it."<<std::endl;
    side_L = true;
  }

  TFile * fs = new TFile( sfile.c_str(), "read");
  if (!fs) return;

  if(convert_to_cm){
    TH2F* htmp = (TH2F*) fs -> Get(sname.c_str());

    std::string title = std::string("pixels to cm: ") + Form("%.2f",cmperpix) + " cm/pix";
    int nbinx = htmp -> GetNbinsX(), nbiny = htmp->GetNbinsY();
    h_frame = new TH2F((sname).c_str(), title.c_str(), nbinx, 0., nbinx*cmperpix, nbiny, 0., nbiny*cmperpix);
    h_frame -> GetXaxis()->SetTitle("X direction (cm)");
    h_frame -> GetYaxis()->SetTitle("Y direction (cm)");
    for (int ix=1; ix<nbinx; ix++){
      for (int iy=1; iy<nbiny; iy++){
        h_frame -> SetBinContent(ix, iy, htmp->GetBinContent(ix, iy));
        h_frame -> SetBinError(ix, iy, htmp->GetBinError(ix, iy));
      }
    }

    delete htmp;
  }
  else{
    h_frame = (TH2F*) fs -> Get(sname.c_str()) ->Clone((sname).c_str());
    h_frame -> GetXaxis()->SetTitle("X pixels");
    h_frame -> GetYaxis()->SetTitle("Y pixels");
  }
  if (!h_frame) return;

  // if h_fpipe exist, then, get that as well
  // no need to convert to cm for it since it should have
  // been converted to cm in the previous step
  fs -> cd();
  std::string snamepipe = sname+"_framepipe";
  TH2F* hptmp = (TH2F*) fs -> Get(snamepipe.c_str());
  if(hptmp){
    h_fpipe = (TH2F*) hptmp -> Clone();
    h_fpipe -> SetDirectory(NULL);
    delete hptmp;
    std::cout<<"INFO: got frame pipe at "<<h_fpipe<<std::endl;
  }else{
    std::cout<<"INFO: frame pipe NOT found. using name: "<< sname<<"_framepipe."<<std::endl;
  }


  T_liquid = this-> guessTLiquid();
  NXPIX = h_frame -> GetNbinsX();
  NYPIX = h_frame -> GetNbinsY();
  h_frame -> SetDirectory(NULL);

  TTree *t0 = (TTree *) fs->Get("info");
  if(t0){
    t0 -> SetBranchAddress("T_liquid", &T_liquid);
    t0 -> SetBranchAddress("conversionToCM", &conversionToCM);
    t0 -> SetBranchAddress("NXPIX", &NXPIX);
    t0 -> SetBranchAddress("NYPIX", &NYPIX);
    t0 -> SetBranchAddress("side_L", &side_L);
    t0 -> SetBranchAddress("cmperpix", &cmperpix);
    t0 -> SetBranchAddress("lens", &lens);
    t0 -> SetBranchAddress("crop_point1_x", &(crop_point1.first) );
    t0 -> SetBranchAddress("crop_point1_y", &(crop_point1.second));
    t0 -> SetBranchAddress("pipe_point1_x", &(pipe_point1.first) );
    t0 -> SetBranchAddress("pipe_point1_y", &(pipe_point1.second));
    t0 -> GetEntry(0);
    delete t0;
  }
  fs -> Close();
  delete fs;
}

HFrame::HFrame(std::vector<std::string> sfiles, int n_xpixels, int n_ypixels, std::string shout, std::string stree, std::string stemp, std::string sxbrh, std::string sybrh, bool convert_to_cm){
  crop_point1 = std::make_pair(-1, -1);
  crop_point2 = std::make_pair(-1, -1);
  pipe_point1 = std::make_pair(-1, -1);
  pipe_point2 = std::make_pair(-1, -1);

  conversionToCM = convert_to_cm;
  cmperpix = CMPERPIXELPRE;
  lens = 80;
  T_liquid = -999.;
  NXPIX = n_xpixels;
  NYPIX = n_ypixels;
  if (!h_frame){
    if(convert_to_cm) h_frame = new TH2F( shout.c_str(), "", NXPIX, 0., cmperpix * NXPIX, NYPIX, 0., cmperpix * NYPIX);
    else              h_frame = new TH2F( shout.c_str(), "", NXPIX, 0., NXPIX, NYPIX, 0., NYPIX);
  }
  h_frame -> SetDirectory(NULL);

  int nfile = sfiles.size();
  if (nfile<=0){
    std::cout<<"ERROR: number of files "<<nfile<<std::endl;
    return;
  }

  float scl = 1. / nfile;
  for (int jf =0; jf<nfile; jf++){

    std::cout<<jf<<"\r"; std::cout.flush();

    std::string sfile = sfiles.at(jf);
    TFile * f1 = new TFile( sfile.c_str(), "read");
    if (!f1){
      std::cout<<"ERROR: file "<<sfile<<" NOT found!!!"<<std::endl;
      return;
    }

    TTree * t1 = (TTree*)f1->Get(stree.c_str());
    if(!t1){ delete f1; return; }

    int xpos, ypos;
    float temperature;
    t1->SetBranchAddress(sxbrh.c_str(), &xpos);
    t1->SetBranchAddress(sybrh.c_str(), &ypos);
    t1->SetBranchAddress(stemp.c_str(), &temperature);

    long int nentry = t1->GetEntries();
    if (nentry <=0) return;

    int xmax_tree = 0, ymax_tree = 0;
    for (int je=0;je<nentry;je++){
      t1->GetEntry(je);

      // in tree file, position starting from 0, why ???????
      xpos += 1;
      ypos += 1;

      if(jf==0){
        h_frame -> SetBinContent(xpos, ypos, temperature * scl);
        h_frame -> SetBinError(xpos, ypos, 0.1);
      }else{
        float tmp = temperature * scl + h_frame -> GetBinContent(xpos, ypos);
        h_frame -> SetBinContent(xpos, ypos, tmp);
      }

      if (xmax_tree < xpos) xmax_tree = xpos;
      if (ymax_tree < ypos) ymax_tree = ypos;

    }
    delete t1;
    f1->Close();
    delete f1;
    if(xmax_tree !=NXPIX || ymax_tree !=NYPIX) std::cout<<sfile<<": found frame, xmax "<<xmax_tree<<", ymax "<<ymax_tree<<std::endl;
  }

  T_liquid = this-> guessTLiquid();
}


bool HFrame::writeToRootFile(std::string sfile){
  if( !h_frame ) return false;
  TFile * f0 = new TFile( sfile.c_str(), "recreate");
  if (!f0) return false;

  if(h_frame){
    TH2F* h1 = (TH2F*) h_frame->Clone(); // h_frame -> GetName() 
    h1 ->SetDirectory(f0); 
    h1->Write();
    delete h1;
  }

  if(h_fpipe){
    TH2F* h2 = (TH2F*) h_fpipe->Clone(); 
    h2 ->SetDirectory(f0); 
    h2->Write();
    delete h2;
  }
  for (auto pone : m_pipemap){
    TH1F * h1 = (TH1F*) pone.second -> Clone();
    h1 ->SetDirectory(f0); 
    h1 -> Write();
    delete h1;
  }
  // write out basic information only when T_liquid is set!
  // otherwise, consider the basic information is not provided.
  if(T_liquid > -998.){
    TTree *t0 = new TTree("info","information about the frame");
    t0 -> Branch("T_liquid", &T_liquid);
    t0 -> Branch("conversionToCM", &conversionToCM);
    t0 -> Branch("NXPIX", &NXPIX);
    t0 -> Branch("NYPIX", &NYPIX);
    t0 -> Branch("side_L", &side_L, "side_L/O");
    t0 -> Branch("cmperpix", &cmperpix);
    t0 -> Branch("lens", &lens);
    t0 -> Branch("crop_point1_x", &(crop_point1.first) );
    t0 -> Branch("crop_point1_y", &(crop_point1.second));
    t0 -> Branch("pipe_point1_x", &(pipe_point1.first) );
    t0 -> Branch("pipe_point1_y", &(pipe_point1.second));
    t0 -> Fill(); // fill the current values, once!
    t0-> Write();
    delete t0;
  }
  f0->Close();
  delete f0;
  return true;
}

bool HFrame::rotate(HFrame::angle ang_clock){

  if (!h_frame) return false;
  if (ang_clock != C090 && ang_clock != C180 && ang_clock != C270) return false;

  TH2F* h_tmp = (TH2F*) h_frame->Clone("h_tmp");
  int nbinx = h_tmp->GetNbinsX();
  int nbiny = h_tmp->GetNbinsY();

  for (int ix=1; ix<=nbinx; ix++){
    for (int iy=1; iy<=nbiny; iy++){

      int tmp_ix = iy, tmp_iy = nbinx - ix + 1; // C090
      if     (ang_clock == C180){tmp_ix = nbinx - ix + 1; tmp_iy = nbiny - iy + 1;}
      else if(ang_clock == C270){tmp_ix = nbiny - iy + 1; tmp_iy = ix;}

      float temperature = h_tmp->GetBinContent(ix, iy);
      float temp_error  = h_tmp->GetBinError(ix, iy);
      h_frame -> SetBinContent(tmp_ix, tmp_iy, temperature);
      h_frame -> SetBinError(tmp_ix, tmp_iy, temp_error);
    }
  }
  delete h_tmp;
  return true;
}

bool HFrame::mirror( bool byX){

  if (!h_frame) return false;

  TH2F* h_tmp = (TH2F*) h_frame->Clone("h_tmp");
  int nbinx = h_tmp->GetNbinsX();
  int nbiny = h_tmp->GetNbinsY();


  for (int ix=1; ix<=nbinx; ix++){
    for (int iy=1; iy<=nbiny; iy++){

      unsigned int tmp_ix = ix, tmp_iy = nbiny - iy + 1; // by X, change Y
      if(!byX){ tmp_ix = nbinx - ix + 1; tmp_iy = iy; }

      float temperature = h_tmp->GetBinContent(ix, iy);
      float temp_error  = h_tmp->GetBinError(ix, iy);
      h_frame -> SetBinContent(tmp_ix, tmp_iy, temperature);
      h_frame -> SetBinError(tmp_ix, tmp_iy, temp_error);

    }
  }
  delete h_tmp;
  return true;
}


TH1F* HFrame::getLine(int x1, int y1, int x2, int y2, bool is_fpipe){
  if ( (!is_fpipe && !h_frame) || (is_fpipe && !h_fpipe) ) return NULL;
  TH2F * h_tmp = NULL;
  if (is_fpipe) h_tmp = h_fpipe;
  else          h_tmp = h_frame;

  bool Xdir = true;
  int nbin = std::abs(x1 - x2) + 1;
  if( nbin < std::abs(y1 - y2) + 1){
    nbin = std::abs(y1 - y2) + 1;
    Xdir = false;
  }
  if(nbin<=1) return NULL;

  std::string sname = "frame1D_";//std::string(h_frame->GetName());
  sname += Form("x%dy%d_x%dy%d",x1,y1,x2,y2);
  TH1F* h1 = new TH1F(sname.c_str(),"",nbin,0.,float(nbin));
  for (int ib = 1; ib <= nbin; ib++){
    int tmp_ix = x1 - 1 + ib * (x1 < x2 ? 1 : -1);
    int tmp_iy = y1     + int( float(y2 - y1) * float(ib - 1) / std::abs(x2 - x1) );
    if(!Xdir){
      tmp_ix = x1     + int( float(x2 - x1) * float(ib - 1) / std::abs(y2 - y1) );
      tmp_iy = y1 - 1 + ib * (y1 < y2 ? 1 : -1);
    }
    float temperature = h_tmp->GetBinContent(tmp_ix, tmp_iy);
    float temp_error  = h_tmp->GetBinError(tmp_ix, tmp_iy);
    // if no error, fitting 1 D histogram wouldn't work.
    if (temp_error <=0.) temp_error = 0.1;
    h1 -> SetBinContent(ib, temperature);
    h1 -> SetBinError(ib, temp_error);
  }
  h_tmp = NULL;
  return h1;
}

bool HFrame::crop( int x1, int y1, int x4, int y4, bool to_fpipe){
  if (!h_frame) return false;
  int nbinx = x4 - x1 + 1;
  int nbiny = y4 - y1 + 1;
  if (nbinx <= 1 || nbiny <= 1) return false;

  if(crop_point1.first <=0){
    if(to_fpipe){
      pipe_point1 = std::make_pair(x1, y1); 
      pipe_point2 = std::make_pair(x4, y4); 
    }else{
      crop_point1 = std::make_pair(x1, y1); 
      crop_point2 = std::make_pair(x4, y4); 
    }
  }
  return this->crop(x1, y1, x4,y1, x1,y4, to_fpipe);

}

bool HFrame::crop( int x1, int y1, int x2, int y2, int x3, int y3, bool to_fpipe){

  if (!h_frame) return false;
  int nbinx = x2 - x1 + 1;
  int nbiny = y3 - y1 + 1;
  if (nbinx <= 1 || nbiny <= 1) return false;

  TH2F* h_tmp = (TH2F*) h_frame->Clone("h_tmp");
  std::string sname = static_cast<std::string>(h_frame -> GetName());
  float xmax = float(nbinx), ymax = float(nbiny);
  if(conversionToCM){
    xmax = cmperpix * nbinx;
    ymax = cmperpix * nbiny;
  }

  if(to_fpipe) sname += "_framepipe";
  TH2F * h_tmp2 = new TH2F(sname.c_str(), "", nbinx, 0., xmax, nbiny, 0., ymax);
  h_tmp2 -> GetXaxis() -> SetTitle( h_tmp -> GetXaxis() -> GetTitle() );
  h_tmp2 -> GetYaxis() -> SetTitle( h_tmp -> GetYaxis() -> GetTitle() );
  //  std::cout<<"in crop, x title "<<h_tmp -> GetXaxis() -> GetTitle()<<std::endl;

  for (int ix=1; ix<=nbinx; ix++){
    for (int iy=1; iy<=nbiny; iy++){

      int tmp_ix = x1 - 1 + ix +  int( float(x3 - x1) * float(iy - 1) / (y3 - y1) );
      int tmp_iy = y1 - 1 + iy +  int( float(y2 - y1) * float(ix - 1) / (x2 - x1) );
      float temperature = h_tmp->GetBinContent(tmp_ix, tmp_iy);
      float temp_error  = h_tmp->GetBinError(tmp_ix, tmp_iy);

      h_tmp2 -> SetBinContent(ix, iy, temperature);
      h_tmp2 -> SetBinError(ix, iy, temp_error);
    }
  }
  if(to_fpipe){
    h_fpipe = h_tmp2;
  }else{
    if(h_frame) delete h_frame;
    h_frame = h_tmp2;
  }

  delete h_tmp;
  return true;

}

bool HFrame::cropRawToStave( float xsize_cm, float ysize_cm){
  int xpix = static_cast<int>(xsize_cm / cmperpix), ypix = static_cast<int>(ysize_cm / cmperpix);
  if (!h_frame || xpix <=1 || ypix <= 1) return false;

  const float Tambient = 20.;
  double Tmin = h_frame -> GetMinimum();
  double Tmax = h_frame -> GetMaximum();
  bool is_negative = this->isNegative();

  int ix1 = -1, iy1 = -1;
  double Tmm = -999.;
  if(is_negative) Tmm = 999.;
  for (int ix=1; ix<=NXPIX-xpix; ix++){
    for (int iy=1; iy<=NYPIX-ypix; iy++){

      double Tavg = h_frame -> Integral(ix, ix+xpix-1, iy, iy+ypix-1);
      Tavg /= (xpix * ypix);
      if ((is_negative && Tavg < Tmm ) || (!is_negative && Tavg > Tmm) ){
        Tmm = Tavg;
        ix1 = ix;
        iy1 = iy;
      }
    }
  }
  std::cout<<"crop auto, before end-of-stave: p1 "<<ix1<<" "<< iy1<<" p2 "<<  ix1+xpix-1<<" "<< iy1<<" p3 "<<  ix1<<" "<< iy1+ypix-1<<std::endl;
  // include the end of stave card now. The size of the end of stave card is about
  // 7.5% * xpix and 35% * ypix
  int xpix_ec = xpix / 13, ypix_ec = ypix / 3;
  if (iy1 - ypix_ec + 1 <= 0){
    std::cout<<"crop_auto failed, found x1 "<<ix1<<" y1 "<<iy1<<" too low in Y"<<std::endl;
    return false;
  }else if (iy1 + ypix + ypix_ec - 2 > NYPIX){
    std::cout<<"crop_auto failed, found x1 "<<ix1<<" y1 "<<iy1<<" too high in Y"<<std::endl;
    return false;
  }
  double Tsum_bl = h_frame -> Integral(ix1, ix1+xpix_ec-1, iy1 - ypix_ec + 1, iy1); // bottom left
  double Tsum_tl = h_frame -> Integral(ix1, ix1+xpix_ec-1, iy1 + ypix - 1, iy1 + ypix + ypix_ec - 2); // top left
  double Tsum_br = h_frame -> Integral(ix1+xpix-xpix_ec, ix1+xpix-1, iy1 - ypix_ec + 1, iy1); // bottom right
  double Tsum_tr = h_frame -> Integral(ix1+xpix-xpix_ec, ix1+xpix-1, iy1 + ypix - 1, iy1 + ypix + ypix_ec - 2); // top right
  if(is_negative){
    if((Tsum_bl < Tsum_tl && Tsum_bl < Tsum_tr) || (Tsum_br < Tsum_tl && Tsum_br < Tsum_tr)){ iy1 -= ypix_ec; ypix += ypix_ec;}
    else  ypix += ypix_ec;
  }else{
    if((Tsum_bl > Tsum_tl && Tsum_bl > Tsum_tr) || (Tsum_br > Tsum_tl && Tsum_br > Tsum_tr)){ iy1 -= ypix_ec; ypix += ypix_ec;}
    else  ypix += ypix_ec;
  }

  std::cout<<"crop auto: p1 "<<ix1<<" "<< iy1<<" p2 "<<  ix1+xpix-1<<" "<< iy1<<" p3 "<<  ix1<<" "<< iy1+ypix-1<<std::endl;

  // found points one(X1, Y1) and two (X2,Y2)
  crop_point1.first  = ix1;
  crop_point1.second = iy1;
  crop_point2.first  = ix1+xpix-1;
  crop_point2.second = iy1+ypix-1;

  bool is_to_pipe = false;
  return this->crop(ix1, iy1,  ix1+xpix-1, iy1,  ix1, iy1+ypix-1, is_to_pipe);
}

bool HFrame::cropStaveToPipe( float xsize_cm, float ysize_cm, float yedge_cm){

  int xpix = static_cast<int>(xsize_cm / cmperpix), ypix = static_cast<int>(ysize_cm / cmperpix);
  int ypix_edge = static_cast<int>(yedge_cm / cmperpix);

  if (!h_frame || xpix <=1 || ypix <= 1) return false;
  int NXBin = h_frame -> GetNbinsX();
  int NYBin = h_frame -> GetNbinsY();
  if(xpix > NXBin || ypix > NYBin) return false;


  const float Tambient = 20.;
  double Tmin = h_frame -> GetMinimum();
  double Tmax = h_frame -> GetMaximum();
  bool is_negative = this->isNegative();

  int ix1 = -1, iy1 = -1;
  double Tmm = -999.;
  if(is_negative) Tmm = 999.;
  for (int ix=1; ix<=NXBin-xpix; ix++){
    for (int iy=1; iy<NYBin-ypix; iy++){

      double Tavg  = h_frame -> Integral(ix, ix+xpix-1, iy, iy);
      Tavg += h_frame -> Integral(ix, ix+xpix-1, iy+ypix, iy+ypix);
      Tavg /= (xpix * 2);
      if ((is_negative && Tavg < Tmm ) || (!is_negative && Tavg > Tmm) ){
        Tmm = Tavg;
        ix1 = ix;
        iy1 = iy;
      }
    }
  }
  iy1 -= ypix_edge;
  ypix += 2 * ypix_edge;
  if (iy1 < 1 || iy1 + ypix > NYBin ) return false;

  // found points one(X1, Y1) and two (X2,Y2)
  pipe_point1.first  = ix1;
  pipe_point1.second = iy1;
  pipe_point2.first  = ix1+xpix-1;
  pipe_point2.second = iy1+ypix-1;

  bool is_to_pipe = true;
  return this->crop(ix1, iy1,  ix1+xpix-1, iy1,  ix1, iy1+ypix-1, is_to_pipe);
}

float HFrame::guessTLiquid(){
  float TLguess = -999.;
  if (!h_frame) return TLguess;
  // if T_liquid is not updated
  // use the closest 5 C to the Min/Max temperature.
  int TempLow = int(h_frame -> GetMinimum());
  int TempHig = int(h_frame -> GetMaximum());
  if(TempLow < 10){
    TempLow = int(TempLow / 10) * 10 - 5 - (TempLow % 10 >= 5 ? 5 : 0);
    TLguess = float(TempLow);
  }else if(TempHig>30){
    TempHig = int(TempHig / 10) * 10 + 5 + (TempHig % 10 >= 5 ? 5 : 0);
    TLguess = float(TempHig);
  }else{
    TLguess = float( (TempLow + TempHig) / 2 );
  }
  return TLguess;
} 

bool HFrame::isNegative(){
  if(T_liquid < -998.){
    std::cout<<"WARNING: liquid temperature is not set! Guess: ";
    T_liquid = this -> guessTLiquid();
    std::cout<< T_liquid <<std::endl;
  }

  if(T_liquid < 20.){
    return true;
  }else if(T_liquid > 26.){
    return false;
  }else{
    std::cout<<"ERROR: liquid temperature is set to [20, 26], hard to decide the sign! Return false. "<<std::endl;
    return false;
  }
}


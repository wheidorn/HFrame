#define HFRAME_CPP
#include "HFrame.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TTree.h"
#include <math.h>

bool HFrame::writeToRootFile(std::string sfile, std::string sname){
  if( !h_frame ) return false;
  int nbinx = h_frame->GetNbinsX();
  int nbiny = h_frame->GetNbinsY();

  TFile * f0 = new TFile( sfile.c_str(), "recreate");

  TH2F* h2 = (TH2F*) h_frame->Clone((sname).c_str());
  h2 ->SetDirectory(f0); 
  h2->Write();
  delete h2;
  for (auto h1 : h_lines){
    std::string s1 = std::string(h1->GetName()); 
    TH1F* h1x = (TH1F*) h1->Clone(s1.c_str());
    h1x -> Write();
    delete h1x;
  }
  f0->Close();
  delete f0;
  return true;
}

bool HFrame::fillboxTrees(std::vector<std::string> sfiles, std::string shout, std::string stree, std::string stemp, std::string sxbrh, std::string sybrh){
  // create a 2D Histogram, if it doesn't exist.
  if (!h_frame) h_frame = new TH2F( shout.c_str(), "", NXPIX, 0., cmperpix * NXPIX, NYPIX, 0., cmperpix * NYPIX);
  h_frame -> SetDirectory(NULL);

  int nfile = sfiles.size();
  std::cout<<"filling frame: number of files "<<nfile<<std::endl;
  if (nfile<=0) return false;

  float scl = 1. / nfile;
  for (int jf =0; jf<nfile; jf++){

    std::cout<<jf<<"\r"; std::cout.flush();

    std::string sfile = sfiles.at(jf);
    TFile * f1 = new TFile( sfile.c_str(), "read");
    if (!f1) return false;

    TTree * t1 = (TTree*)f1->Get(stree.c_str());
    if(!t1){ delete f1; return false; }

    int xpos, ypos;
    float temperature;
    t1->SetBranchAddress(sxbrh.c_str(), &xpos);
    t1->SetBranchAddress(sybrh.c_str(), &ypos);
    t1->SetBranchAddress(stemp.c_str(), &temperature);

    long int nentry = t1->GetEntries();
    if (nentry <=0) return false;

    int xmax_tree = 0, ymax_tree = 0;
    for (int je=0;je<nentry;je++){
      t1->GetEntry(je);

      // in tree file, position starting from 0, why ??!!
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

  return true;
}


bool HFrame::setFrame(TH2F * h2, std::string shout){ 
  if( !h2 ) return false; 
  else{ 
    h_frame = h2; 
    h_frame -> SetName(shout.c_str()); 
  } 
  return true;
}


bool HFrame::fillboxHist(std::string sfile, std::string sname, bool convert_to_cm){
  fs = new TFile( sfile.c_str(), "read");
  if (!fs) return false;

  if(convert_to_cm){
    TH2F* htmp = (TH2F*) fs -> Get(sname.c_str());

    std::string title = std::string("pixels to cm: ") + Form("%.2f",cmperpix) + " cm/pix";
    int nbinx = htmp -> GetNbinsX(), nbiny = htmp->GetNbinsY();
    h_frame = new TH2F((sname+"_copy").c_str(), title.c_str(), nbinx, 0., nbinx*cmperpix, nbiny, 0., nbiny*cmperpix);
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
    h_frame = (TH2F*) fs -> Get(sname.c_str()) ->Clone((sname+"_copy").c_str());
    h_frame -> GetXaxis()->SetTitle("X pixels");
    h_frame -> GetYaxis()->SetTitle("Y pixels");
  }
  if (!h_frame) return false;

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

      // mirror image by X, then change Y
      unsigned int tmp_ix = ix, tmp_iy = nbiny - iy + 1;
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

bool HFrame::crop( int x1, int y1, int x2, int y2, int x3, int y3, bool convert_to_cm){

  if (!h_frame) return false;
  int nbinx = x2 - x1 + 1;
  int nbiny = y3 - y1 + 1;
  if (nbinx <= 1 || nbiny <= 1) return false;

  TH2F* h_tmp = (TH2F*) h_frame->Clone("h_tmp");
  std::string sname = static_cast<std::string>(h_frame -> GetName());
  delete h_frame;

  float xmax = float(nbinx), ymax = float(nbiny);
  if(convert_to_cm){
    xmax = cmperpix * nbinx;
    ymax = cmperpix * nbiny;
  }

  // Histogram starts from 0. always --> can be changed to point x1,y1.
  h_frame = new TH2F(sname.c_str(), "", nbinx, 0., xmax, nbiny, 0., ymax);
  h_frame -> GetXaxis() -> SetTitle( h_tmp -> GetXaxis() -> GetTitle() );
  h_frame -> GetYaxis() -> SetTitle( h_tmp -> GetYaxis() -> GetTitle() );

  for (int ix=1; ix<=nbinx; ix++){
    for (int iy=1; iy<=nbiny; iy++){

      int tmp_ix = x1 - 1 + ix +  int( float(x3 - x1) * float(iy - 1) / (y3 - y1) );
      int tmp_iy = y1 - 1 + iy +  int( float(y2 - y1) * float(ix - 1) / (x2 - x1) );
      float temperature = h_tmp->GetBinContent(tmp_ix, tmp_iy);
      float temp_error  = h_tmp->GetBinError(tmp_ix, tmp_iy);

      h_frame -> SetBinContent(ix, iy, temperature);
      h_frame -> SetBinError(ix, iy, temp_error);
    }
  }
  delete h_tmp;
  return true;

}

TH1F* HFrame::getLine(int x1, int y1, int x2, int y2, bool m_record){
  if (!h_frame) return NULL;
  bool Xdir = true;
  int nbin = std::abs(x1 - x2) + 1;

  if( nbin < std::abs(y1 - y2) + 1){
    nbin = std::abs(y1 - y2) + 1;
    Xdir = false;
  }
  if(nbin<=1) return NULL;

  std::string sname = "line_";
  sname += Form("x%dy%d_x%dy%d",x1,y1,x2,y2);
  TH1F* h1 = new TH1F(sname.c_str(),"",nbin,0.,float(nbin));

  for (int ib = 1; ib <= nbin; ib++){
    int tmp_ix = x1 - 1 + ib * (x1 < x2 ? 1 : -1);
    int tmp_iy = y1     + int( float(y2 - y1) * float(ib - 1) / std::abs(x2 - x1) );
    if(!Xdir){
      tmp_ix = x1     + int( float(x2 - x1) * float(ib - 1) / std::abs(y2 - y1) );
      tmp_iy = y1 - 1 + ib * (y1 < y2 ? 1 : -1);
    }


      float temperature = h_frame->GetBinContent(tmp_ix, tmp_iy);
      float temp_error  = h_frame->GetBinError(tmp_ix, tmp_iy);

      // if a 1D Histogram has no error, fitting wouldn't work!!
      // set the error of 0.1 C!
      if (temp_error <=0.) temp_error = 0.1;

      h1 -> SetBinContent(ib, temperature);
      h1 -> SetBinError(ib, temp_error);
  }
  if(m_record) h_lines.push_back( h1 );
  return h1;
}

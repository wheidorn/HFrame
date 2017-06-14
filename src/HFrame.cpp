#define HFRAME_CPP
#include "HFrame.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TTree.h"
#include "TFile.h"
#include <math.h>       /* std::abs */


using namespace std;
bool HFrame::read(std::string sfile, std::string sname){
  TFile * fs = new TFile( sfile.c_str(), "read");
  if (!fs) return false;

  TH2F* h_frame = (TH2F*) fs -> Get(sname.c_str());
  if(!h_frame) return false;
  int nbinx = h_frame -> GetNbinsX();
  int nbiny = h_frame -> GetNbinsY();
  if(nbinx > m_NXPixel || nbiny > m_NYPixel){
    cout<<"too many bins in input histogram: X "<<nbinx<<" > "<<m_NXPixel<<" ? or Y "<<nbiny<<" > "<<m_NYPixel<<endl;
    nbinx = m_NXPixel;
    nbiny = m_NYPixel;
  }

  for (int ix=0; ix<nbinx; ix++){
    for (int iy=0; iy<nbiny; iy++){
      m_data[ix][iy] = h_frame -> GetBinContent(ix+1, iy+1);
    }
  }
  return true;
}

bool HFrame::read(std::vector<std::string> sfiles, std::string stree, std::string stemp, std::string sxbrh, std::string sybrh){
  int nfile = sfiles.size();

  if (nfile<=0){
    std::cout<<"ERROR: number of files "<<nfile<<std::endl;
    return false;
  }

  float scl = 1. / nfile;
  for (int jf =0; jf<nfile; jf++){
    std::cout<<jf<<"\r"; std::cout.flush();
    std::string sfile = sfiles.at(jf);
    bool status = this -> read( sfile, stree, stemp, sxbrh, sybrh, scl);

    // if found one failed reading, return false;
    if (!status) return false;
  }

  return true;
}

bool HFrame::read(string sfile, std::string stree, std::string stemp, std::string sxbrh, std::string sybrh, float scale){

  TFile * f1 = new TFile( sfile.c_str(), "read");
  if (!f1){
    std::cout<<"ERROR: file "<<sfile<<" NOT found!!!"<<std::endl;
    return false;
  }

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

    // in tree file, position starting from 0, need to +1 if filling a Histogram!!!
    if(xpos >= m_NXPixel || ypos >= m_NYPixel){
      cout<<" over flow X "<<xpos << " >= "<< m_NXPixel<<" ? or Y "<<ypos <<" >= " << m_NYPixel<<" skip!!! "<<endl;
      continue;
    }

    // update m_data
    m_data[xpos][ypos] += temperature * scale;
  }

  if(t1) delete t1;
  f1->Close();
  if(f1) delete f1;
  return true;
}


bool HFrame::writeToRootFile(std::string sfile){

  if( sfile == "" ) return false;
  
  TFile * f0 = new TFile( sfile.c_str(), "recreate");
  if (!f0) return false;
  TH2F* h_hist = new TH2F(m_name_hist.c_str(), "", m_NXPixel, 0., m_NXPixel * m_CMperPixel, m_NYPixel, 0., m_NYPixel * m_CMperPixel);

  int NXPixel_crop = m_crop_x2 - m_crop_x1 + 1, NYPixel_crop = m_crop_y2 - m_crop_y1 + 1;
  TH2F* h_crop = new TH2F(m_name_crop.c_str(), "", NXPixel_crop, 0., NXPixel_crop * m_CMperPixel, NYPixel_crop, 0., NYPixel_crop * m_CMperPixel);

  // pipe frame is a sub-figure of the cropped frame
  int NXPixel_pipe = m_pipe_x2 - m_pipe_x1 + 1, NYPixel_pipe = m_pipe_y2 - m_pipe_y1 + 1;
  float pipe_xmin = (m_pipe_x1 - m_crop_x1) * m_CMperPixel;
  float pipe_xmax = NXPixel_pipe * m_CMperPixel + pipe_xmin;
  float pipe_ymin = (m_pipe_y1 - m_crop_y1) * m_CMperPixel;
  float pipe_ymax = NYPixel_pipe * m_CMperPixel + pipe_ymin;
  TH2F* h_pipe = new TH2F(m_name_pipe.c_str(), "", NXPixel_pipe, pipe_xmin, pipe_xmax, NYPixel_pipe, pipe_ymin, pipe_ymax);


  // add 0.5% error, since otherwise the fit of pipe peak wouldn't work.
  const float rel_error = 0.005; 
  for (int ix=0; ix< m_NXPixel; ix++){
    for (int iy=0; iy< m_NYPixel; iy++){
      h_hist -> SetBinContent(ix+1, iy+1, m_data[ix][iy]);
      h_hist -> SetBinError(ix+1, iy+1, m_data[ix][iy] * rel_error);

      int ix_crop = ix - m_crop_x1, iy_crop = iy - m_crop_y1;
      if(ix_crop >=0 && iy_crop >=0 && ix_crop<NXPixel_crop && iy_crop < NYPixel_crop){
        h_crop -> SetBinContent(ix_crop+1, iy_crop+1, m_data[ix][iy]);
        h_crop -> SetBinError(ix_crop+1, iy_crop+1, m_data[ix][iy] * rel_error);
      }

      int ix_pipe = ix - m_pipe_x1, iy_pipe = iy - m_pipe_y1;
      if(ix_pipe >=0 && iy_pipe >=0 && ix_pipe<NXPixel_pipe && iy_pipe < NYPixel_pipe){
        h_pipe -> SetBinContent(ix_pipe+1, iy_pipe+1, m_data[ix][iy]);
        h_pipe -> SetBinError(ix_pipe+1, iy_pipe+1, m_data[ix][iy] * rel_error);
      }
    }
  }
  h_hist ->Write();
  h_crop ->Write();
  h_pipe ->Write();

  TTree *t0 = new TTree("info","information about the frame");
  int NXPixel = m_NXPixel;
  int NYPixel = m_NYPixel;
  float CMperPixel = m_CMperPixel;
  t0 -> Branch("NXPixel", &NXPixel);
  t0 -> Branch("NYPixel", &NYPixel);
  t0 -> Branch("CMperPixel", &CMperPixel);

  t0 -> Branch("T_liquid", &m_T_liquid);
  t0 -> Branch("side", &m_side, "m_side/B");
  t0 -> Branch("lens", &m_lens);
  t0 -> Branch("crop_x1", &m_crop_x1);
  t0 -> Branch("crop_y1", &m_crop_y1);
  t0 -> Branch("crop_x2", &m_crop_x2);
  t0 -> Branch("crop_y2", &m_crop_y2);
  t0 -> Branch("pipe_x1", &m_pipe_x1);
  t0 -> Branch("pipe_y1", &m_pipe_y1);
  t0 -> Branch("pipe_x2", &m_pipe_x2);
  t0 -> Branch("pipe_y2", &m_pipe_y2);
  t0 -> Fill(); // fill the current values, once!
  t0-> Write();
  if(t0) delete t0;
  f0->Close();
  if(f0) delete f0;
  return true;
}
    
void HFrame::print(string sout){

  FILE* fout=fopen( sout.c_str(), "w");

  fprintf(fout, "lens:                   %3d \n",   m_lens);
  fprintf(fout, "side:                   %c \n",    m_side);
  fprintf(fout, "cm per pixel:           %5.3f \n", m_CMperPixel);
  fprintf(fout, "Temperature liquid (C): %4.1f \n", m_T_liquid);
  fprintf(fout, "Number of pixels X:     %4d \n",   m_NXPixel);
  fprintf(fout, "Number of pixels Y:     %4d \n",   m_NYPixel);
  fprintf(fout, "crop (x1,y1,x2,y2):     %3d,%3d,%3d,%3d \n", m_crop_x1, m_crop_y1, m_crop_x2, m_crop_y2);
  fprintf(fout, "pipe (x1,y1,x2,y2):     %3d,%3d,%3d,%3d \n", m_pipe_x1, m_pipe_y1, m_pipe_x2, m_pipe_y2);
  fclose(fout);
  //if(fout) delete fout;
}

bool HFrame::rotate(HFrame::angle ang_clock){

  if (ang_clock != C090 && ang_clock != C180 && ang_clock != C270) return false;

  vector<vector<float>> tmpdata = m_data;

  for (int ix=0; ix< m_NXPixel; ix++){
    for (int iy=0; iy< m_NYPixel; iy++){

      int to_ix = iy, to_iy = m_NXPixel - ix - 1; // C090
      if     (ang_clock == C180){to_ix = m_NXPixel - ix - 1; to_iy = m_NYPixel - iy - 1;}
      else if(ang_clock == C270){to_ix = m_NYPixel - iy - 1; to_iy = ix;}

      m_data[to_ix][to_iy] = tmpdata[ix][iy];

    }
  }
  tmpdata.clear();
  return true;
}

float HFrame::average(int x1, int y1, int x2, int y2){
  float Tavg = -999.;
  if(m_data.size() ==0){
    cout<<"WARNING: frame is empty now! No average found!"<<endl;
    return Tavg;
  }
  int ix1 = x1, iy1 = y1, ix2 = x2, iy2 = y2;
  if(x1 < 0 || x1 >= m_NXPixel) ix1  = 0;
  if(x2 < 0 || x2 >= m_NXPixel) ix2  = m_NXPixel - 1;
  if(y1 < 0 || y1 >= m_NYPixel) iy1  = 0;
  if(y2 < 0 || y2 >= m_NYPixel) iy2  = m_NYPixel - 1;
  if(ix1 > ix2){
    cout<<"WARNING: x1 "<<ix1<<" > x2 "<<x2<<" swap!"<<endl;
    int itmp = ix1; ix1 = ix2; ix2 = itmp;
  } 
  if(iy1 > iy2){
    cout<<"WARNING: y1 "<<iy1<<" > y2 "<<iy2<<" swap!"<<endl;
    int itmp = iy1; iy1 = iy2; iy2 = itmp;
  } 


  Tavg = 0.; 
  int NPixelAllAvg = (x2 - x1 + 1) * (y2 - y1 + 1);
  for (int jx=ix1; jx<=ix2; jx++){
    for (int jy=iy1; jy<=iy2; jy++){
      Tavg += m_data[jx][jy] / NPixelAllAvg;
    }
  }
  return Tavg;
}

bool HFrame::mirror( bool byX){

  vector<vector<float>> tmpdata = m_data;
  for (int ix=0; ix< m_NXPixel; ix++){
    for (int iy=0; iy< m_NYPixel; iy++){

      // by X, change Y
      int to_ix = ix, to_iy = m_NYPixel - iy - 1;
      if(!byX){ to_ix = m_NXPixel - ix - 1; to_iy = iy; }

      m_data[to_ix][to_iy] = tmpdata[ix][iy];
    }
  }
  return true;
}


bool HFrame::findCropPositions( float xsize_cm, float ysize_cm){
  int xpix = static_cast<int>(xsize_cm / m_CMperPixel);
  int ypix = static_cast<int>(ysize_cm / m_CMperPixel);
  if (m_data.size() == 0 || xpix <=1 || ypix <= 1) return false;
  int NPixCrop = xpix * ypix;

  bool is_negative = this->isNegative();

  int ix1 = -1, iy1 = -1;
  float Tmm = -999.;
  if(is_negative) Tmm = 999;
  for (int ix=0; ix<m_NXPixel-xpix; ix++){
    for (int iy=0; iy<m_NYPixel-ypix; iy++){

      float Tavg = this -> average(ix, iy, ix+xpix-1, ypix - 1 + iy);
      if ((is_negative && Tavg < Tmm ) || (!is_negative && Tavg > Tmm) ){
        Tmm = Tavg;
        ix1 = ix;
        iy1 = iy;
      }
    }
  }

  // include the end of stave card now. The size of the end of stave card is about
  // 7.5% * xpix and 35% * ypix
  int xpix_ec = xpix / 13, ypix_ec = ypix / 3;
  while (iy1 - ypix_ec < 0 || iy1 + ypix + ypix_ec >= m_NYPixel){
    cout<<"WARNING: not enough place to put End-of-Stave card. reduce size!"<<endl;
    ypix_ec--;
  }
  float Tavg_bl = this -> average(ix1, iy1 - ypix_ec, ix1 + xpix_ec - 1, iy1 - 1                 ); // bottom left
  float Tavg_tl = this -> average(ix1, iy1 + ypix,    ix1 + xpix_ec - 1, iy1 + ypix + ypix_ec - 1); // top left
  float Tavg_br = this -> average(ix1 + xpix - xpix_ec, iy1 - ypix_ec, ix1 + xpix - 1, iy1 - 1                 ); // bottom right
  float Tavg_tr = this -> average(ix1 + xpix - xpix_ec, iy1 + ypix,    ix1 + xpix - 1, iy1 + ypix + ypix_ec - 1); // top right

  float Tavg_bmin = Tavg_bl, Tavg_bmax = Tavg_br;
  if(Tavg_bl > Tavg_br){Tavg_bmin = Tavg_br; Tavg_bmax = Tavg_bl;}
  float Tavg_tmin = Tavg_tl, Tavg_tmax = Tavg_tr;
  if(Tavg_tl > Tavg_tr){Tavg_tmin = Tavg_tr; Tavg_tmax = Tavg_tl;}

  if((is_negative && Tavg_bmin < Tavg_tmin) || (!is_negative && Tavg_bmin > Tavg_tmin)) iy1 -= ypix_ec; 
  ypix += ypix_ec;

  // found points one(X1, Y1) and two (X2,Y2)
  m_crop_x1 = ix1;
  m_crop_y1 = iy1;
  m_crop_x2 = ix1+xpix-1;
  m_crop_y2 = iy1+ypix-1;

  return true;
}




bool HFrame::findPipePositions( float xsize_cm, float ysize_cm, float yedge_cm){

  int xpix = static_cast<int>(xsize_cm / m_CMperPixel), ypix = static_cast<int>(ysize_cm / m_CMperPixel);
  int ypix_edge = static_cast<int>(yedge_cm / m_CMperPixel);

  if (m_data.size() ==0 || xpix <=1 || ypix <= 1) return false;

  bool is_negative = this->isNegative();

  int ix_a = 0, iy_a = 0; //starting point
  int ix_z = m_NXPixel-xpix, iy_z = m_NYPixel-ypix;
  if(m_crop_x1 >=0) ix_a = m_crop_x1;
  if(m_crop_y1 >=0) iy_a = m_crop_y1;
  if(m_crop_x2 >=0) ix_z = m_crop_x2 - xpix;
  if(m_crop_y2 >=0) iy_z = m_crop_y2 - ypix;

  int ix1 = -1, iy1 = -1;
  float Tmm = -999.;
  if(is_negative) Tmm = 999.;
  for (int ix=ix_a; ix<ix_z; ix++){
    for (int iy=iy_a; iy<iy_z; iy++){

      float Tavg = this -> average(ix, iy, ix + xpix - 1, iy + ypix - 1);
      if ((is_negative && Tavg < Tmm ) || (!is_negative && Tavg > Tmm) ){
        Tmm = Tavg;
        ix1 = ix;
        iy1 = iy;
      }
    }
  }

  iy1 -= ypix_edge;
  ypix += 2 * ypix_edge;
  if (iy1 < 1 || iy1 + ypix > m_NYPixel ) return false;

  // found points one(X1, Y1) and two (X2,Y2)
  m_pipe_x1 = ix1;
  m_pipe_y1 = iy1;
  m_pipe_x2 = ix1+xpix-1;
  m_pipe_y2 = iy1+ypix-1;
  return true;
}

void HFrame::setPipePositions(float x1_cm, float y1_cm, float x2_cm, float y2_cm){

  if(x1_cm > x2_cm || y1_cm > y2_cm){
    std::cout<<"ERROR: setting pipe positions with x1 "<<x1_cm<<" > y2 "<<x2_cm<<" OR y1 "<<y1_cm<<" > y2 "<<y2_cm<<std::endl;
    return;
  }
  int ix1 = static_cast<int>(x1_cm / m_CMperPixel);
  int ix2 = static_cast<int>(x2_cm / m_CMperPixel);
  int iy1 = static_cast<int>(y1_cm / m_CMperPixel);
  int iy2 = static_cast<int>(y2_cm / m_CMperPixel);

  if(m_crop_x1 >=0 && m_crop_x1 + ix1 <m_NXPixel ) m_pipe_x1 = m_crop_x1 + ix1;
  if(m_crop_y1 >=0 && m_crop_y1 + iy1 <m_NYPixel ) m_pipe_y1 = m_crop_y1 + iy1;
  if(m_crop_x2 >=0 && m_crop_x1 + ix2 <m_NXPixel ) m_pipe_x2 = m_crop_x1 + ix2;
  if(m_crop_y2 >=0 && m_crop_y1 + iy2 <m_NYPixel ) m_pipe_y2 = m_crop_y1 + iy2;

  return;
}


bool HFrame::isNegative(){
  if(m_T_liquid < -998.){
    std::cout<<"ERROR: liquid temperature is not set! RETURN false! ";
    return false;
  }

  if(m_T_liquid < 15.1){
    return true;
  }else if(m_T_liquid > 29.9){
    return false;
  }else{
    std::cout<<"ERROR: liquid temperature is set to [15, 30], hard to decide the sign! RETURN false! "<<std::endl;
    return false;
  }
}


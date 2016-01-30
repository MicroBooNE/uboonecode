
////////////////////////////////////////////////////////////////////////
//
//    trivia : Class to hold the Muon Counter System (MuCS) reco fields 
//    author : Matt Bass
//    e-mail : Matthew.Bass@physics.ox.ac.uk
//
////////////////////////////////////////////////////////////////////////

#include "MuCSRecoData.h"

#include "cetlib/exception.h"

namespace MuCS {

MuCSRecoData::MuCSRecoData()
{}

MuCSRecoData::~MuCSRecoData()
{}
  
MuCSRecoData::MuCSRecoData( Float_t theta_xy, Float_t theta_xy_rms, Float_t x,
                    Float_t x_rms, Float_t theta_yz, Float_t theta_yz_rms,
                    Float_t z, Float_t z_rms, Float_t y, Int_t xmatches, Int_t zmatches ) 
{
  ftheta_xy=theta_xy;
  ftheta_xy_rms=theta_xy_rms;
  fx=x;
  fx_rms=x_rms;
  ftheta_yz=theta_yz;
  ftheta_yz_rms=theta_yz_rms;
  fz=z;
  fz_rms=z_rms;
  fy=y;
  fxmatches=xmatches;
  fzmatches=zmatches;
}

//compute from underlying variables
Float_t MuCSRecoData::theta() const{
  return -M_PI/2+acos((pow(1+pow(tan(ftheta_xy),-2)+pow(tan(ftheta_yz),-2),-0.5)));
}
Float_t MuCSRecoData::phi() const{
  Float_t ltheta = theta();
  return atan2(sin(ltheta)/tan(ftheta_yz),sin(ltheta)/tan(ftheta_xy));
}

Float_t MuCSRecoData::theta_xy() const{
  return ftheta_xy;
}
Float_t MuCSRecoData::theta_xy_rms() const{
  return ftheta_xy_rms;
}
Float_t MuCSRecoData::x() const{
  return fx;
}
Float_t MuCSRecoData::x_rms() const{
  return fx_rms;
}

Float_t MuCSRecoData::theta_yz() const{
  return ftheta_yz;
}
Float_t MuCSRecoData::theta_yz_rms() const{
  return ftheta_yz_rms;
}
Float_t MuCSRecoData::z() const{
  return fz;
}
Float_t MuCSRecoData::z_rms() const{
  return fz_rms;
}

Float_t MuCSRecoData::y() const{
  return fy;
}
Int_t MuCSRecoData::xmatches() const{
  return fxmatches;
}
Int_t MuCSRecoData::zmatches() const{
  return fzmatches;
}
}

////////////////////////////////////////////////////////////////////////
//
//    The end !
//
////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////
//
//    trivia : Class to hold the Muon Counter System (MuCS) reco fields 
//    author : Matt Bass
//    e-mail : Matthew.Bass@physics.ox.ac.uk
//
////////////////////////////////////////////////////////////////////////

#ifndef MuCSRecoData_H
#define MuCSRecoData_H

#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <iosfwd>

#include "TMath.h"

namespace MuCS {
 
class MuCSRecoData 
{

 public:
  MuCSRecoData();
  virtual ~MuCSRecoData();  
  MuCSRecoData( Float_t theta_xy, Float_t theta_xy_rms, Float_t x,
                    Float_t x_rms, Float_t theta_yz, Float_t theta_yz_rms,
                    Float_t z, Float_t z_rms, Float_t y, Int_t xmatches, Int_t zmatches ); 
                    
  Float_t theta_xy() const;
  Float_t theta_xy_rms() const;
  Float_t x() const;
  Float_t x_rms() const;
  Float_t theta_yz() const;
  Float_t theta_yz_rms() const;
  Float_t z() const;
  Float_t z_rms() const;
  Float_t y() const;
  Int_t xmatches() const;
  Int_t zmatches() const;
  
  Float_t theta() const;
  Float_t phi() const;
  
 private:
  Float_t ftheta_xy;
  Float_t ftheta_xy_rms;
  Float_t fx;
  Float_t fx_rms;
  Float_t ftheta_yz;
  Float_t ftheta_yz_rms;
  Float_t fz;
  Float_t fz_rms;  
  Float_t fy;
  Int_t fxmatches;
  Int_t fzmatches;
};
}
#endif 


#ifndef _GETFOM_H
#define _GETFOM_H

#include "datatypes/ub_BeamHeader.h"
#include "datatypes/ub_BeamData.h"
#include <string>
#include "TMath.h"
#include "TH1.h"

namespace bmd
{
  float getFOM(std::string beam, const  gov::fnal::uboone::datatypes::ub_BeamHeader& bh, const std::vector<gov::fnal::uboone::datatypes::ub_BeamData>& bd);
  
  
  double fractionmisstarget( double off_x, double off_y, double sig_x, double sig_y);
  double fpeaks(Double_t *x, Double_t *par);
  double profilewidth( std::vector<double>& y );
  Int_t profilesigmaROOT( TH1F* myhist, Int_t* npeaks, Double_t* par );
  int profilesigma( std::vector<double>& mwire, int* npeaks, double* par );
  double findedge( double offset, double width );
}

#endif

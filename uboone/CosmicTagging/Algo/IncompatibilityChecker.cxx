#ifndef INCOMPATIBILITYCHECKER_CXX
#define INCOMPATIBILITYCHECKER_CXX

#include "IncompatibilityChecker.h"
#include <iostream>

namespace flashana {

  IncompatibilityChecker::IncompatibilityChecker()
  {
    _sigmaThreshold = 5;
    _nBinsRequirement = 1;
    _useFlashPosition = false;
  }

  void IncompatibilityChecker::Configure(fhicl::ParameterSet const& pset)
  {
    _sigmaThreshold   = pset.get< double > ( "SigmaThreshold"   );
    _nBinsRequirement = pset.get< int    > ( "NBinsRequirement" );
    _useFlashPosition = pset.get< bool   > ( "UseFlashPosition" );
  }

  void IncompatibilityChecker::PrintConfig() {

    std::cout << "--- IncompatibilityChecker configuration:" << std::endl;
    std::cout << "---   _sigmaThreshold   = " << _sigmaThreshold << std::endl;
    std::cout << "---   _nBinsRequirement = " << _nBinsRequirement << std::endl;
    std::cout << "---   _useFlashPosition = " << _useFlashPosition << std::endl;

  }

  bool IncompatibilityChecker::CheckIncompatibility(const Flash_t &flash, const Flash_t &flash_hypo) {

    //std::cout << "IncompatibilityChecker::CheckIncompatibility starts" << std::endl;
    // Now we have two spectra: the 'true' one (flash), and the hypo one (flash_hypo)
    // We want to understand is they are incompatible
    //std::cout << "pe for pmt 0 from hypo " << flash_hypo.pe_v[0] << std::endl;

    if (flash.pe_v.size() != flash_hypo.pe_v.size()) {
     throw cet::exception("Flash and hypo flash pe vector size mismatch."); 
    }

    bool areIncompatibleByBin    = false;
    bool areIncompatibleByTotPe  = false;
    bool areIncompatibleByFlsPos = false;

    double totalPE_flash = 0.;
    double totalPE_hypo  = 0.;

    // Check bin by bin incompatibility
    int nIncompBins = 0;
    for (unsigned int pmt = 0; pmt < flash.pe_v.size(); pmt++) {
      double error  = std::sqrt(flash_hypo.pe_v[pmt]);
      double nsigma = (flash_hypo.pe_v[pmt] - flash.pe_v[pmt]) / error;
      if (nsigma > _sigmaThreshold) {
        nIncompBins ++;
        if (nIncompBins >= _nBinsRequirement){
          areIncompatibleByBin = true;
        }
      }
      totalPE_flash += flash.pe_v[pmt];
      totalPE_hypo  += flash_hypo.pe_v[pmt];
    }

    // Check overall incompatibility
    double error  = std::sqrt(totalPE_hypo);
    double nsigma = (totalPE_hypo - totalPE_flash) / error;
    if (nsigma > _sigmaThreshold) {
      areIncompatibleByTotPe = true;
    }

    if(_useFlashPosition) {
      double beamFlashZmin = flash.z - flash.z_err;
      double beamFlashZmax = flash.z + flash.z_err;
      //std::cout << "beamFlashZmin " << beamFlashZmin << std::endl;
      //std::cout << "beamFlashZmax " << beamFlashZmax << std::endl;
      //std::cout << "flash_hypo.z  " << flash_hypo.z  << std::endl;
      if (flash_hypo.z > beamFlashZmax || flash_hypo.z < beamFlashZmin) {
        areIncompatibleByFlsPos = true;
      }
    }

    if(!_useFlashPosition) {
      if( (areIncompatibleByBin || areIncompatibleByTotPe) ) {
      }
      if( areIncompatibleByBin ) {
        //std::cout << "Flashes are incompatible by bin or total pe." << std::endl;
        return true;
      }
    } else {
      if( areIncompatibleByBin && areIncompatibleByFlsPos) {
        //std::cout << "Flashes are incompatible by (bin or total pe) and flash position." << std::endl;
        return true;
      }
    }

    return false;
  
  }
}


#endif

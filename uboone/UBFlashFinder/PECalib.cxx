#ifndef PECALIB_CXX
#define PECALIB_CXX

#include "PECalib.h"

namespace pmtana {

  PECalib::PECalib()
  {}

  void PECalib::Configure(const Config_t &pset)
  {
    
    _spe_area_gain_v = pset.get<std::vector<double> >("SPEAreaGain");
    if(_spe_area_gain_v.size() != NOpDets()) {
      std::cerr << "SPEAreaGain array size (" << _spe_area_gain_v.size()
		<< ") != NOpDets (" << NOpDets() << ")..." << std::endl;
      throw std::exception();
    }

    _spe_amp_gain_v = pset.get<std::vector<double> >("SPEAmpGain");
    if(_spe_amp_gain_v.size() != NOpDets()) {
      std::cerr << "SPEAmpGain array size (" << _spe_amp_gain_v.size()
		<< ") != NOpDets (" << NOpDets() << ")..." << std::endl;
      throw std::exception();
    }

    _cosmic_ophit_correction_v = pset.get<std::vector<double> >("RecoPECorrection");
    _cosmic_ophit_ratio_mean_v = pset.get<std::vector<double> >("RecoPERatioMean");
    _cosmic_ophit_ratio_std_v  = pset.get<std::vector<double> >("RecoPERatioStd");
    _cosmic_ophit_amppe_cut_v  = pset.get<std::vector<double> >("RecoAmpPECut");
    _cosmic_ophit_areape_cut_v  = pset.get<std::vector<double> >("RecoAreaPECut");

    if(_cosmic_ophit_correction_v.size() != NOpDets()) {
      std::cerr << "RecoPECorrection array size (" << _cosmic_ophit_correction_v.size()
		<< ") != NOpDets (" << NOpDets() << ")..." << std::endl;
      throw std::exception();
    }
    if(_cosmic_ophit_ratio_mean_v.size() != NOpDets()) {
      std::cerr << "RecoPERatioMean array size (" << _cosmic_ophit_ratio_mean_v.size()
		<< ") != NOpDets (" << NOpDets() << ")..." << std::endl;
      throw std::exception();
    }
    if(_cosmic_ophit_ratio_std_v.size() != NOpDets()) {
      std::cerr << "RecoPERatioStd array size (" << _cosmic_ophit_ratio_std_v.size()
		<< ") != NOpDets (" << NOpDets() << ")..." << std::endl;
      throw std::exception();
    }
    if(_cosmic_ophit_amppe_cut_v.size() != NOpDets()) {
      std::cerr << "RecoAmpPECut array size (" << _cosmic_ophit_amppe_cut_v.size()
		<< ") != NOpDets (" << NOpDets() << ")..." << std::endl;
      throw std::exception();
    }
    if(_cosmic_ophit_areape_cut_v.size() != NOpDets()) {
      std::cerr << "RecoAreaPECut array size (" << _cosmic_ophit_areape_cut_v.size()
		<< ") != NOpDets (" << NOpDets() << ")..." << std::endl;
      throw std::exception();
    }

    _relative_qe_v = pset.get<std::vector<double> >("RelativeQE");
    if(_relative_qe_v.size() != NOpDets()) {
      std::cerr << "RelativeQE array size (" << _relative_qe_v.size()
		<< ") != NOpDets (" << NOpDets() << ")..." << std::endl;
      throw std::exception();
    }
  }

  double PECalib::CosmicPE(const size_t opdet, const double area, const double amp) const
  {
    if( opdet > NOpDets() ) {
      std::cerr << "OpDet ID " << opdet << " exceeding max # of OpDet (" << NOpDets() << ")" << std::endl;
      throw std::exception();
    }

    double area_pe = area / _spe_area_gain_v[opdet];
    double amp_pe = amp / _spe_amp_gain_v[opdet];

    double ratio_thresh = _cosmic_ophit_ratio_mean_v[opdet] - 2 * _cosmic_ophit_ratio_std_v[opdet] ;
    if( (area_pe > _cosmic_ophit_areape_cut_v[opdet])
	||
	(amp_pe  > _cosmic_ophit_amppe_cut_v[opdet] && (area_pe / amp_pe) > ratio_thresh) )

      area_pe *= _cosmic_ophit_correction_v[opdet];

    area_pe *= _relative_qe_v[opdet];

    return area_pe;
  }
  
  double PECalib::BeamPE(const size_t opdet, const double area, const double amp) const
  {
    if( opdet > NOpDets() ) {
      std::cerr << "OpDet ID " << opdet << " exceeding max # of OpDet (" << NOpDets() << ")" << std::endl;
      throw std::exception();
    }

    double area_pe = area / _spe_area_gain_v[opdet] * _relative_qe_v[opdet];

    return area_pe;

  }

}


#endif

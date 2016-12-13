#ifndef PHOTONLIBHYPOTHESIS_CXX
#define PHOTONLIBHYPOTHESIS_CXX

#include "PhotonLibHypothesis.h"
#include "larsim/PhotonPropagation/PhotonVisibilityService.h"
//#include "OpT0Finder/PhotonLibrary/PhotonVisibilityService.h"

#include "larcore/Geometry/Geometry.h"
#include "larcore/Geometry/OpDetGeo.h"
#include "uboone/Geometry/UBOpReadoutMap.h"

namespace flashana {

  static PhotonLibHypothesisFactory __global_PhotonLibHypothesisFactory__;

  PhotonLibHypothesis::PhotonLibHypothesis(const std::string name)
    : BaseFlashHypothesis(name)
  {}

  void PhotonLibHypothesis::_Configure_(const Config_t &pset)
  {
    _global_qe = pset.get<double>("GlobalQE");
    _qe_v      = pset.get<std::vector<double> >("CCVCorrection");

    ::art::ServiceHandle<geo::Geometry> geo;
    ::art::ServiceHandle<geo::UBOpReadoutMap> ub_geo;

    if(geo->NOpDets() != _qe_v.size()) {
      std::cout << "CCVCorrection factor array has size " << _qe_v.size() << " != # OpDet " << geo->NOpDets() << std::endl;
      throw std::exception();
    }
  }
  
  void PhotonLibHypothesis::FillEstimate(const QCluster_t& trk,
					 Flash_t &flash) const
  {
    art::ServiceHandle<phot::PhotonVisibilityService> vis;
    static double xyz[3] = {0.};

    size_t n_pmt = BaseAlgorithm::NOpDets();//n_pmt returns 0 now, needs to be fixed
    
    for ( auto& v : flash.pe_v ) v = 0;
    
    for ( size_t ipmt = 0; ipmt < n_pmt; ++ipmt) {

      for ( size_t ipt = 0; ipt < trk.size(); ++ipt) {
	
        auto const& pt = trk[ipt];
	
        double q = pt.q;
	
        //q *= ::phot::PhotonVisibilityService::GetME().GetVisibility( pt.x, pt.y, pt.z, ipmt)*0.0093;
        xyz[0] = pt.x;
	xyz[1] = pt.y;
	xyz[2] = pt.z;
	q *= vis->GetVisibility(xyz,ipmt) * 0.0093 * _global_qe / _qe_v[ipmt];
        flash.pe_v[ipmt] += q;
	//std::cout << "PMT : " << ipmt << " [x,y,z] -> [q] : [" << pt.x << ", " << pt.y << ", " << pt.z << "] -> [" << q << std::endl;
	
      }
    }

    return;
  }
}
#endif

#ifndef SCANNERALGO_CXX
#define SCANNERALGO_CXX
#include "ScannerAlgo.h"

namespace larlite {

  size_t ScannerAlgo::NameIndex(::larlite::data::DataType_t const data_type,
				std::string const& name) const
  {
    for(size_t i=0; i<fModuleLabel_v[(size_t)data_type].size(); ++i)

      if(fModuleLabel_v[(size_t)data_type][i] == name)

	return i;

    throw cet::exception(__PRETTY_FUNCTION__) << "Invalid producer name \""
					      << name.c_str()
					      << "\" for data type "
					      << ::larlite::data::kDATA_TREE_NAME[(size_t)data_type].c_str();

    return std::numeric_limits<size_t>::max();
  }

  void ScannerAlgo::EventClear()
  {
    for(auto& name_bool_v : fDataReadFlag_v) {
      for(auto& name_bool : name_bool_v)
	name_bool.second = false;
    }
    fPtrIndex_mctruth.clear();
    fPtrIndex_gtruth.clear();
    fPtrIndex_mcflux.clear();
    fPtrIndex_mcpart.clear();
    fPtrIndex_mcshower.clear();
    fPtrIndex_mctrack.clear();
    fPtrIndex_simch.clear();
    fPtrIndex_rawdigit.clear();
    fPtrIndex_wire.clear();
    fPtrIndex_hit.clear();
    fPtrIndex_ophit.clear();
    fPtrIndex_opflash.clear();
    fPtrIndex_trigger.clear();
    fPtrIndex_cluster.clear();
    fPtrIndex_shower.clear();
    fPtrIndex_vertex.clear();
    fPtrIndex_cosmictag.clear();
    fPtrIndex_track.clear();
    fPtrIndex_calo.clear();
    fPtrIndex_sps.clear();
    fPtrIndex_seed.clear();
    fPtrIndex_end2d.clear();
    fPtrIndex_partid.clear();
    fPtrIndex_pfpart.clear();
    fPtrIndex_pcaxis.clear();
  }
}

#endif

////////////////////////////////////////////////////////////////////////
/// \file  LArRawInputDriverUBooNE.h
/// \brief Source to convert raw binary files to root files for MicroBooNE
///
/// \Original Version from:
/// \version $Id: T962ConvertBinaryToROOT.h,v 1.7 2010/01/14 19:20:33 brebel Exp $
/// \author  brebel@fnal.gov, soderber@fnal.gov
/// \MicroBooNE author: jasaadi@fnal.gov (with much help from Eric and Wes)
////////////////////////////////////////////////////////////////////////

#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Core/ProductRegistryHelper.h"
#include "art/Framework/Core/PrincipalMaker.h"
#include "art/Framework/Core/FileBlock.h"
#include "art/Framework/Principal/RunPrincipal.h"
#include "art/Framework/Principal/SubRunPrincipal.h"
#include "art/Framework/Principal/EventPrincipal.h"
#include "art/Persistency/Provenance/SubRunID.h"

#include <boost/serialization/string.hpp>
#include <boost/serialization/version.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/list.hpp>
#include <boost/serialization/split_member.hpp>
#include <boost/serialization/binary_object.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include "boost/date_time/gregorian/gregorian.hpp"
#include "boost/date_time/posix_time/posix_time.hpp"
#include "boost/date_time/gregorian/greg_serialize.hpp"
#include "boost/date_time/posix_time/time_serialize.hpp"
#include "boost/serialization/vector.hpp"

#include "RawData/uboone_datatypes/share/boonetypes.h"

#include <fstream>
#include <string>
#include <vector>

namespace raw {
class RawDigit;
class OpDetPulse;
class BeamInfo;
class DAQHeader;
}

///Conversion of binary data to root files
namespace lris {

  struct daqid_t {
    daqid_t():crate(-1),card(-1),channel(-1) {};
    daqid_t(int crate_id, int card_id, int channel_id):
      crate(crate_id),card(card_id),channel(channel_id) {};

    int crate;
    int card;
    int channel;
  };

  bool operator<(daqid_t const& lhs,daqid_t const& rhs) {
    bool is_less=false;
    if (lhs.crate   == rhs.crate && 
	lhs.card    == rhs.card  && 
	lhs.channel <  rhs.channel) is_less=true;
    else if (lhs.crate == rhs.crate && 
	     lhs.card  <  rhs.card) is_less=true;
    else if (lhs.crate < rhs.crate) is_less=true;
    return is_less;
  }

  class LArRawInputDriverUBooNE {
    /// Class to fill the constraints on a template argument to the class,
    /// FileReaderSource
  public:
    // Required constructor
    LArRawInputDriverUBooNE(fhicl::ParameterSet const &pset,
			    art::ProductRegistryHelper &helper,
			    art::PrincipalMaker const &pm);
    
    // Required by FileReaderSource:
    void closeCurrentFile();
    void readFile(std::string const &name,
		  art::FileBlock* &fb);
    bool readNext(art::RunPrincipal* const &inR,
		  art::SubRunPrincipal* const &inSR,
		  art::RunPrincipal* &outR,
		  art::SubRunPrincipal* &outSR,
		  art::EventPrincipal* &outE);

    //Other functions
    void initChannelMap();
    bool processNextEvent(std::vector<raw::RawDigit>& digitList,
			  std::vector<raw::OpDetPulse>& pmtDigitList,
			  raw::DAQHeader& daqHeader,
			  raw::BeamInfo& beamInfo);

  private:
    
    art::PrincipalMaker            fPrincipalMaker;
    art::SubRunID                  fCurrentSubRunID;

    std::ifstream fInputStream;
    std::vector<std::streampos> fEventLocation;

    uint32_t fNumberOfEvents;
    uint32_t fEventCounter; //event counter
    bool fHuffmanDecode;

    std::map<daqid_t, int> fChannelMap;
    
  };  // LArRawInputDriverUBooNE

}

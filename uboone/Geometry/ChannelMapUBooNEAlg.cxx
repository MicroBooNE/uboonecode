////////////////////////////////////////////////////////////////////////
/// \file  ChannelMapUBooNEAlg.cxx
/// \brief MicroBooNE optical detector channel mapping
///   
/// Terminology:
///  OpDet: the physical PMT or light guide paddle
///  HardwareChannel: The gain copy number. 0=Low A, 1 = Low B, 3=High A, 3 = High B  
///  OpChannel: The readout channel number
///
/// \version $Id:  $
/// \author  taritree@mit.edu
////////////////////////////////////////////////////////////////////////

#include "uboonecode/uboone/Geometry/ChannelMapUBooNEAlg.h"
#include "Geometry/AuxDetGeo.h"
#include "Geometry/CryostatGeo.h"
#include "Geometry/TPCGeo.h"
#include "Geometry/PlaneGeo.h"
#include "Geometry/WireGeo.h"

#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "messagefacility/MessageLogger/MessageLogger.h" 

namespace geo {

  //----------------------------------------------------------------------------
  ChannelMapUBooNEAlg::ChannelMapUBooNEAlg( fhicl::ParameterSet const& pvals, fhicl::ParameterSet const& sortingParameters )
    : ChannelMapStandardAlg( sortingParameters )
  {
    // parameter set will come from UBooNEGomeotryHelper service
    LoadOpticalMapData( pvals );
  }

  //----------------------------------------------------------------------------
  ChannelMapUBooNEAlg::~ChannelMapUBooNEAlg()
  {
  }

  //----------------------------------------------------------------------------
  // OPTICAL CHANNELS
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  unsigned int ChannelMapUBooNEAlg::NOpChannels(unsigned int NOpDets) const
  {
    return fNReadoutChannels;
  }

  //----------------------------------------------------------------------------
  unsigned int ChannelMapUBooNEAlg::NOpHardwareChannels(unsigned int opDet) const
  {
    auto it = fPMT2channels.find( opDet );
    if ( it!=fPMT2channels.end() )
      return (*it).second.size();
    else
      throw std::runtime_error( "ChannelMapUBooNEAlg::NOpHardwareChannels : Invalid opdet value" );
  }

  //----------------------------------------------------------------------------
  unsigned int ChannelMapUBooNEAlg::OpChannel(unsigned int pmtID, unsigned int copynum) const
  {
    // converts OpDet and Gain Channel
    unsigned int uniqueChannel = fPMT2channels.at(pmtID).at(copynum);
    return uniqueChannel;
  }

  //----------------------------------------------------------------------------
  unsigned int ChannelMapUBooNEAlg::OpDetFromOpChannel(unsigned int opChannel) const
  {
    unsigned int pmtID = fChannel2pmt.at( opChannel );
    return pmtID;
  }

  //----------------------------------------------------------------------------
  unsigned int ChannelMapUBooNEAlg::HardwareChannelFromOpChannel(unsigned int opChannel) const
  {
    auto it = fPMT2channels.find( OpDetFromOpChannel(opChannel)  );
    if ( it!=fPMT2channels.end() ) {
      auto readoutChList = it->second;
      for ( unsigned int hwch=0; hwch<readoutChList.size(); hwch++ )
	if ( opChannel==readoutChList.at(hwch) )
	  return hwch;
    }

    throw std::runtime_error( "Could not find copy index of readout channel" );
  }

  //----------------------------------------------------------------------------
  bool ChannelMapUBooNEAlg::IsValidOpChannel(unsigned int opChannel, unsigned int NOpDets) const {
    auto it=fChannel2pmt.find( opChannel );
    if ( it!=fChannel2pmt.end() )
      return true;
    return false;
  }
  
  //----------------------------------------------------------------------------
  void ChannelMapUBooNEAlg::LoadOpticalMapData( fhicl::ParameterSet const& pset ) {
    fNOpDets = pset.get< unsigned int >("numberOfDetectors");
    fNReadoutChannels = 0;    
    // ----------------------------------------------------------------------
    // map between opdet ID and Readout Channel Number
    for (unsigned int iop=0; iop<fNOpDets; iop++) {
      char entryname[50];
      sprintf(entryname,"OpDet%d_channels",iop);
      std::vector< unsigned int > chinput =  pset.get< std::vector<unsigned int> >( entryname );
      fPMT2channels[ iop ] = chinput;

      //std::cout << entryname << ": [";
      for (std::vector<unsigned int>::iterator it_ch=chinput.begin(); it_ch!=chinput.end(); it_ch++) {
	fChannel2pmt[ *it_ch ] = iop;
	//std::cout << *it_ch << ",";
	fNReadoutChannels++;
      }
      //std::cout << "]" << std::endl;
    }
    
  }
  
} // namespace

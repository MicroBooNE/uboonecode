////////////////////////////////////////////////////////////////////////
/// \file  UBOpReadoutMap.h
/// \brief MicroBooNE optical detector mapping between reaodout channel number and optical channel categories
///
/// \version $Id:  $
/// \author taritree@mit.edu
////////////////////////////////////////////////////////////////////////
#ifndef GEO_CHANNELMAP_UBOONE_ALG_H
#define GEO_CHANNELMAP_UBOONE_ALG_H

#include <vector>
#include <set>
#include <iostream>
#include <map>

#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Registry/ServiceMacros.h"

#include "uboone/Geometry/UBOpChannelTypes.h" // uboonecode

namespace geo {

  class UBOpReadoutMap {

  public:

    UBOpReadoutMap(fhicl::ParameterSet const& pset, art::ActivityRegistry& reg );
    ~UBOpReadoutMap();
    

    // ----------------------------------------------------------------------------
    // OPDET MAPPING
  public:
    /* below is covered by geometry service */
    /* unsigned int NOpChannels(unsigned int NOpDets)                        const; */
    /* unsigned int NOpHardwareChannels(unsigned int opDet)                  const; */
    /* unsigned int OpChannel(unsigned int detNum, unsigned int channel = 0) const; */
    /* unsigned int OpDetFromOpChannel(unsigned int opChannel)               const; */
    /* unsigned int HardwareChannelFromOpChannel(unsigned int opChannel)     const; */
    /* bool         IsValidOpChannel(unsigned int opChannel, unsigned int NOpDets) const; */
    unsigned int NOpLogicChannels() const;
    void         GetLogicChannelList( std::vector< unsigned int >& channels ) const;
    opdet::UBOpticalChannelType_t     GetChannelType( unsigned int opchannel ) const;
    opdet::UBOpticalChannelCategory_t GetChannelCategory( unsigned int opChannel ) const;
    unsigned int GetNumberOfChannelsInCategory( opdet::UBOpticalChannelCategory_t category ) const;
    unsigned int GetChannelNumberFromCrateSlotFEMCh( unsigned int crate, unsigned int slot, unsigned int femch ) const;
    void GetCrateSlotFEMChFromReadoutChannel( unsigned int readoutch, unsigned int& crate, unsigned int& slot, unsigned int& femch ) const;

  private:
    void LoadOpticalReadoutMapData( fhicl::ParameterSet const& p);
    unsigned int fNReadoutChannels;
    std::set< unsigned int > fReadoutChannelSet;
    std::map< opdet::UBOpticalChannelCategory_t, std::set< unsigned int > > fCategoryChannels; // list of channels assigned to each category
    std::map< opdet::UBOpticalChannelType_t, std::set< unsigned int> > fTypeChannels; // list of chanels assignd to channel type
    std::map< unsigned int, opdet::UBOpticalChannelCategory_t > fChannelCategory;
    std::map< unsigned int, opdet::UBOpticalChannelType_t > fChannelType;
    std::set< unsigned int > fLogicChannels;

    class CrateSlotFEMCh {
    public:
      CrateSlotFEMCh( unsigned int _crate, unsigned int _slot, unsigned int _femch ) {
	crate = _crate;
	slot  = _slot;
	femch = _femch;
      }
      ~CrateSlotFEMCh() {};
      
      unsigned int crate;
      unsigned int slot;
      unsigned int femch;

      // comparison
      bool operator<( CrateSlotFEMCh other ) const {
	bool is_less = false;
	if ( this->crate == other.crate &&
	     this->slot  == other.slot  &&
	     this->femch <  other.femch) is_less=true;
	else if (this->crate == other.crate &&
		 this->slot  <  other.slot) is_less=true;
	else if (this->crate <  other.crate) is_less=true;
	return is_less;
      };
      
      // assignment
      CrateSlotFEMCh& operator=(CrateSlotFEMCh other) {
	std::swap(crate, other.crate);
	std::swap(slot,  other.slot);
	std::swap(femch,  other.femch);
        return *this;
      };
    };
    
    std::map< unsigned int, CrateSlotFEMCh > fReadout2CSF;
    std::map< CrateSlotFEMCh, unsigned int> fCSF2Readout;
    
    // ----------------------------------------------------------------
    // Prevent multiple loading
    static unsigned short __fInstances__;
    
  };


}

DECLARE_ART_SERVICE(geo::UBOpReadoutMap, LEGACY)

#endif // GEO_CHANNELMAPSTANDARDALG_H


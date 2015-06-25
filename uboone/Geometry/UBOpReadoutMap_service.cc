////////////////////////////////////////////////////////////////////////
/// \file  UBOpReadoutMap.cxx
/// \brief MicroBooNE optical detector channel mapping
///
/// \version $Id:  $
/// \author  taritree@mit.edu
////////////////////////////////////////////////////////////////////////

#include "UBOpReadoutMap.h"
#include <sstream>

namespace geo {

  unsigned short UBOpReadoutMap::__fInstances__ = 0;

  //----------------------------------------------------------------------------
  UBOpReadoutMap::UBOpReadoutMap(fhicl::ParameterSet const& pset, art::ActivityRegistry& reg )
  {
    // parameter set will come from UBooNEGomeotryHelper service
    LoadOpticalReadoutMapData( pset );
    if ( __fInstances__ > 0 ) {
      throw std::runtime_error("Creating second copy of this class.  Should only be one.");
    }
    __fInstances__++;
  }

  //----------------------------------------------------------------------------
  UBOpReadoutMap::~UBOpReadoutMap()
  {
  }

  //----------------------------------------------------------------------------
  // OPTICAL CHANNELS
  //----------------------------------------------------------------------------

  //----------------------------------------------------------------------------
  opdet::UBOpticalChannelCategory_t UBOpReadoutMap::GetChannelCategory( unsigned int opChannel ) const {
    auto it = fChannelCategory.find( opChannel );
    if ( it!=fChannelCategory.end() )
      return it->second;
    return opdet::Uncategorized;
  }

  //----------------------------------------------------------------------------
  opdet::UBOpticalChannelType_t UBOpReadoutMap::GetChannelType( unsigned int opChannel ) const {
    auto it = fChannelType.find( opChannel );
    if ( it!=fChannelType.end() )
      return it->second;
    throw std::runtime_error( "UBOpReadoutMap::GetChannelType : Could not find channel type (low,high,logic)" );
  }

  //----------------------------------------------------------------------------
  void  UBOpReadoutMap::GetLogicChannelList( std::vector< unsigned int >& channels ) const {
    channels.resize( fLogicChannels.size() );
    std::copy( fLogicChannels.begin(), fLogicChannels.end(), channels.begin() );
  }

  //----------------------------------------------------------------------------
  unsigned int UBOpReadoutMap::GetNumberOfChannelsInCategory( opdet::UBOpticalChannelCategory_t category ) const {
    auto it = fCategoryChannels.find( category );
    if ( it==fCategoryChannels.end() )
      return 0;
    return (*it).second.size();
  }

  //----------------------------------------------------------------------------
  unsigned int UBOpReadoutMap::GetChannelNumberFromCrateSlotFEMCh( unsigned int crate, unsigned int slot, unsigned int femch ) const {
    CrateSlotFEMCh csf( crate, slot, femch );
    auto it = fCSF2Readout.find( csf );
    if ( it==fCSF2Readout.end() ) {
      std::stringstream ss;
      ss << "(crate,slot,femch)=(" << crate << ", " << slot << ", " << femch << " ) entry not found in map to readout channel number.";
      throw std::runtime_error( ss.str() );
    }
    return (*it).second;
  }

  void UBOpReadoutMap::GetCrateSlotFEMChFromReadoutChannel( unsigned int readoutch, unsigned int& crate, unsigned int& slot, unsigned int& femch ) const {
    auto it = fReadout2CSF.find( readoutch );
    if ( it==fReadout2CSF.end() )
      throw std::runtime_error( "readout channel number entry not found in map to (crate,slot,femch)" );
    crate = (*it).second.crate;
    slot = (*it).second.slot;
    femch = (*it).second.femch;
  }
  
  //----------------------------------------------------------------------------
  void UBOpReadoutMap::LoadOpticalReadoutMapData( fhicl::ParameterSet const& pset ) {
    
    // ----------------------------------------------------------------------
    // read in channel types
    fNReadoutChannels = 0;

    for ( unsigned int icat=0; icat<(unsigned int)opdet::NumUBOpticalChannelCategories; icat++ ) {
      std::vector< unsigned int > cat_channels;
      try {
	cat_channels = pset.get< std::vector<unsigned int> >( opdet::UBOpChannelEnumName( (opdet::UBOpticalChannelCategory_t)icat ) );
      }
      catch (...) {
	continue;
      }

      fNReadoutChannels += cat_channels.size();

      // save category channels
      fCategoryChannels[ (opdet::UBOpticalChannelCategory_t)icat ] = std::set< unsigned int >( cat_channels.begin(), cat_channels.end() );
      // save gain channel 
      opdet::UBOpticalChannelType_t chtype = opdet::GetUBTypeFromCategory( (opdet::UBOpticalChannelCategory_t)icat );
      if ( fTypeChannels.find( chtype )==fTypeChannels.end() )
	fTypeChannels[ chtype ] = std::set< unsigned int >( cat_channels.begin(), cat_channels.end() ); // new set
      else {
	for ( auto v : cat_channels )
	  fTypeChannels[ chtype ].insert( v ); // append to set
      }

      for ( auto v : cat_channels ) {
	fChannelCategory[ v ] = (opdet::UBOpticalChannelCategory_t)icat;
	fChannelType[ v ] = chtype;
	if (chtype==opdet::LogicChannel)
	  fLogicChannels.insert( v );
	fReadoutChannelSet.insert( v );
      }
    }//end loop over categories

    //std::cout << "Number of defined readout channels: " << fNReadoutChannels << std::endl;

    // ----------------------------------------------------------------------
    // Read in Crate, Slot, FEMCh
    char readoutname[100];
    for ( auto v : fReadoutChannelSet ) {
      sprintf(readoutname,"ReadoutChannel%d",v);
      std::vector< unsigned int > fichl_csf = pset.get< std::vector<unsigned int> >( readoutname );
      if ( fichl_csf.size()!=3 ) {
	throw std::runtime_error( "Need to have 3 entries for Crate, Slot, FEMCh map." );
      }
      //std::cout << "reading in " << readoutname << " : " << fichl_csf.at(0) << ", " << fichl_csf.at(1) << ", " << fichl_csf.at(2) << std::endl;
      fReadout2CSF.insert( std::make_pair( v, CrateSlotFEMCh( fichl_csf.at(0), fichl_csf.at(1), fichl_csf.at(2) ) ) );
      fCSF2Readout.insert( std::make_pair( CrateSlotFEMCh( fichl_csf.at(0), fichl_csf.at(1), fichl_csf.at(2) ), v ) );
    }
    //std::cout << "size of csf2readout: " << fCSF2Readout.size() << std::endl;
  }

  DEFINE_ART_SERVICE(UBOpReadoutMap)  
} // namespace

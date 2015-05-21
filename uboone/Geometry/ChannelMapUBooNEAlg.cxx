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

  unsigned short ChannelMapUBooNEAlg::__fInstances__ = 0;

  //----------------------------------------------------------------------------
  ChannelMapUBooNEAlg::ChannelMapUBooNEAlg(fhicl::ParameterSet const& sort_pars, fhicl::ParameterSet const& helper_pset )
    : ChannelMapAlg(), fSorter(geo::GeoObjectSorterStandard(sort_pars))
  {
    // parameter set will come from UBooNEGomeotryHelper service
    LoadOpticalMapData( helper_pset );
    if ( __fInstances__ > 0 ) {
      throw std::runtime_error("Creating second copy of this class.  Should only be one.  Get it via Geometry service using GetChannelMapAlg() (and cast to ChannelMapUBooNEAlg)");
    }
    __fInstances__++;
  }

  //----------------------------------------------------------------------------
  ChannelMapUBooNEAlg::~ChannelMapUBooNEAlg()
  {
  }

  //----------------------------------------------------------------------------
  void ChannelMapUBooNEAlg::Initialize( std::vector<geo::CryostatGeo*> & cgeo, 
					  std::vector<geo::AuxDetGeo*>   & adgeo )
  {
    fNcryostat = cgeo.size();
    
    mf::LogInfo("ChannelMapUBooNEAlg") << "Initializing UBooNE ChannelMap...";

    fSorter.SortCryostats(cgeo);
    fSorter.SortAuxDets(adgeo);
    for(size_t c = 0; c < cgeo.size(); ++c) 
      cgeo[c]->SortSubVolumes(fSorter);
    
    fNTPC.resize(fNcryostat);
    fWireCounts.resize(fNcryostat);
    fNPlanes.resize(fNcryostat);
    fFirstWireProj.resize(fNcryostat);
    fOrthVectorsY.resize(fNcryostat);
    fOrthVectorsZ.resize(fNcryostat);
    fPlaneBaselines.resize(fNcryostat);
    fWiresPerPlane.resize(fNcryostat);
    fFirstChannelInNextPlane.resize(fNcryostat);
    fFirstChannelInThisPlane.resize(fNcryostat);
    fViews.clear();
    fPlaneIDs.clear();
    fTopChannel = 0;

    int RunningTotal = 0;

    for(unsigned int cs = 0; cs != fNcryostat; ++cs){
      
      fNTPC[cs] = cgeo[cs]->NTPC();
      
      // Size up all the vectors 
      fWireCounts[cs]             .resize(fNTPC[cs]);
      fFirstWireProj[cs] 	  .resize(fNTPC[cs]);
      fOrthVectorsY[cs]  	  .resize(fNTPC[cs]);
      fOrthVectorsZ[cs]  	  .resize(fNTPC[cs]);
      fPlaneBaselines[cs]	  .resize(fNTPC[cs]);
      fWiresPerPlane[cs] 	  .resize(fNTPC[cs]);
      fNPlanes[cs] 	  	  .resize(fNTPC[cs]);
      fFirstChannelInThisPlane[cs].resize(fNTPC[cs]);
      fFirstChannelInNextPlane[cs].resize(fNTPC[cs]);

      for(unsigned int TPCCount = 0; TPCCount != fNTPC[cs]; ++TPCCount){
	unsigned int PlanesThisTPC = cgeo[cs]->TPC(TPCCount).Nplanes();
	fWireCounts[cs][TPCCount]   .resize(PlanesThisTPC);
	fFirstWireProj[cs][TPCCount].resize(PlanesThisTPC);
	fOrthVectorsY[cs][TPCCount] .resize(PlanesThisTPC);
	fOrthVectorsZ[cs][TPCCount] .resize(PlanesThisTPC);
	fNPlanes[cs][TPCCount]=PlanesThisTPC;
	for(unsigned int PlaneCount = 0; PlaneCount != PlanesThisTPC; ++PlaneCount){

	  fViews.emplace(cgeo[cs]->TPC(TPCCount).Plane(PlaneCount).View());
	  fPlaneIDs.emplace(PlaneID(cs, TPCCount, PlaneCount));
	  double ThisWirePitch = cgeo[cs]->TPC(TPCCount).WirePitch(0, 1, PlaneCount);
	  fWireCounts[cs][TPCCount][PlaneCount] = cgeo[cs]->TPC(TPCCount).Plane(PlaneCount).Nwires();
	  
	  double  WireCentre1[3] = {0.,0.,0.};
	  double  WireCentre2[3] = {0.,0.,0.};
	  
	  const geo::WireGeo& firstWire = cgeo[cs]->TPC(TPCCount).Plane(PlaneCount).Wire(0);
	  const double sth = firstWire.SinThetaZ(), cth = firstWire.CosThetaZ();
	  
	  firstWire.GetCenter(WireCentre1,0);
	  cgeo[cs]->TPC(TPCCount).Plane(PlaneCount).Wire(1).GetCenter(WireCentre2,0);
	  
	  // figure out if we need to flip the orthogonal vector 
	  // (should point from wire n -> n+1)
	  double OrthY = cth, OrthZ = -sth;
	  if(((WireCentre2[1] - WireCentre1[1])*OrthY 
	      + (WireCentre2[2] - WireCentre1[2])*OrthZ) < 0){
	    OrthZ *= -1;
	    OrthY *= -1;
	  }
	  
	  // Overall we are trying to build an expression that looks like
	  //  int NearestWireNumber = round((worldPos.OrthVector - FirstWire.OrthVector)/WirePitch);     
	  // That runs as fast as humanly possible.
	  //
	  // Casting to an int is much faster than c rounding commands like floor().  We have to add 0.5
	  // to account for rounding up as well as down.  floor(A) ~ (int)(A+0.5).  We account for the
	  // 0.5 in the first wire constant to avoid adding it every time.
	  //
	  // We can also predivide everything by the wire pitch so we don't do this in the loop
	  //
	  // Putting this together into the useful constants we will use later per plane and tpc:
	  fOrthVectorsY[cs][TPCCount][PlaneCount] = OrthY / ThisWirePitch;
	  fOrthVectorsZ[cs][TPCCount][PlaneCount] = OrthZ / ThisWirePitch;
	  
	  fFirstWireProj[cs][TPCCount][PlaneCount]  = WireCentre1[1]*OrthY + WireCentre1[2]*OrthZ;
	  fFirstWireProj[cs][TPCCount][PlaneCount] /= ThisWirePitch;
	  
	  // now to count up wires in each plane and get first channel in each plane
	  int WiresThisPlane = cgeo[cs]->TPC(TPCCount).Plane(PlaneCount).Nwires();	
	  fWiresPerPlane[cs] .at(TPCCount).push_back(WiresThisPlane);
	  fPlaneBaselines[cs].at(TPCCount).push_back(RunningTotal);
	  
	  RunningTotal += WiresThisPlane;

	  fFirstChannelInThisPlane[cs].at(TPCCount).push_back(fTopChannel);
	  fTopChannel += WiresThisPlane;
	  fFirstChannelInNextPlane[cs].at(TPCCount).push_back(fTopChannel);

	}// end loop over planes
      }// end loop over TPCs
    }// end loop over cryostats

    // calculate the total number of channels in the detector
    fNchannels = fTopChannel;

    LOG_DEBUG("ChannelMapUBooNE") << "# of channels is " << fNchannels;

    
    

    return;
  }
   
  //----------------------------------------------------------------------------
  void ChannelMapUBooNEAlg::Uninitialize()
  {
  }

  //----------------------------------------------------------------------------
  std::vector<geo::WireID> ChannelMapUBooNEAlg::ChannelToWire(raw::ChannelID_t channel)  const
  {
    std::vector< geo::WireID > AllSegments; 
    unsigned int cstat = 0;
    unsigned int tpc   = 0;
    unsigned int plane = 0;
    unsigned int wire  = 0;
    
    // first check if this channel ID is legal
    if(channel > fTopChannel)
      throw cet::exception("Geometry") << "ILLEGAL CHANNEL ID for channel " << channel << "\n";
    
    // then go find which plane, tpc and cryostat it is in from the information we stored earlier
    bool foundWid(false);
    for(unsigned int csloop = 0; csloop != fNcryostat; ++csloop){
      for(unsigned int tpcloop = 0; tpcloop != fNTPC[csloop]; ++tpcloop){
	for(unsigned int planeloop = 0; planeloop != fFirstChannelInNextPlane[csloop][tpcloop].size(); ++planeloop){
	  if(channel < fFirstChannelInNextPlane[csloop][tpcloop][planeloop]){
	    cstat = csloop;
	    tpc   = tpcloop;
	    plane = planeloop;
	    wire  = channel - fFirstChannelInThisPlane[cstat][tpcloop][planeloop];
	    foundWid = true;
	    break;
	  }	    
	  if(foundWid) break;
	}// end plane loop
	if(foundWid) break;
      }// end tpc loop
      if(foundWid) break;
    }// end cryostat loop

    geo::WireID CodeWire(cstat, tpc, plane, wire);
    
    AllSegments.push_back(CodeWire);

    return AllSegments;
  }

  //----------------------------------------------------------------------------
  raw::ChannelID_t ChannelMapUBooNEAlg::Nchannels() const
  {
    return fNchannels;
  }

  //----------------------------------------------------------------------------
  double ChannelMapUBooNEAlg::WireCoordinate(double YPos, double ZPos,
                                               unsigned int PlaneNo,
                                               unsigned int TPCNo,
                                               unsigned int cstat) const
  {
    // Returns the wire number corresponding to a (Y,Z) position in PlaneNo 
    // with float precision.
    // B. Baller August 2014
    return YPos*fOrthVectorsY[cstat][TPCNo][PlaneNo] 
	 + ZPos*fOrthVectorsZ[cstat][TPCNo][PlaneNo]      
	 - fFirstWireProj[cstat][TPCNo][PlaneNo];
  }

  
  //----------------------------------------------------------------------------
  WireID ChannelMapUBooNEAlg::NearestWireID(const TVector3& worldPos,
					      unsigned int    PlaneNo,
					      unsigned int    TPCNo,
					      unsigned int    cstat) const
  {

    // This part is the actual calculation of the nearest wire number, where we assume
    //  uniform wire pitch and angle within a wireplane
    
    // add 0.5 to have the correct rounding
    int NearestWireNumber = int
      (0.5 + WireCoordinate(worldPos.Y(), worldPos.Z(), PlaneNo, TPCNo, cstat));
    
    // If we are outside of the wireplane range, throw an exception
    // (this response maintains consistency with the previous
    // implementation based on geometry lookup)
    if(NearestWireNumber < 0 ||
       NearestWireNumber >= fWireCounts[cstat][TPCNo][PlaneNo])
    {
      int wireNumber = NearestWireNumber; // save for the output
      
      if(NearestWireNumber < 0 ) NearestWireNumber = 0;
      else                       NearestWireNumber = fWireCounts[cstat][TPCNo][PlaneNo] - 1;
      
      throw InvalidWireIDError("Geometry", wireNumber, NearestWireNumber)
        << "Can't Find Nearest Wire for position (" 
        << worldPos[0] << "," << worldPos[1] << "," << worldPos[2] << ")"
        << " approx wire number # " << wireNumber
        << " (capped from " << NearestWireNumber << ")\n";
    }

    WireID wid(cstat, PlaneNo, TPCNo, (unsigned int)NearestWireNumber);
    return wid;

  }
  
  //----------------------------------------------------------------------------
  // This method returns the channel number, assuming the numbering scheme
  // is heirachical - that is, channel numbers run in order, for example:
  //                                             (Ben J Oct 2011)                   
  //                    Wire1     | 0
  //           Plane1 { Wire2     | 1
  //    TPC1 {          Wire3     | 2
  //           Plane2 { Wire1     | 3   increasing channel number
  //                    Wire2     | 4     (with no gaps)
  //    TPC2 { Plane1 { Wire1     | 5
  //           Plane2 { Wire1     | 6
  //                    Wire2     v 7
  //
  raw::ChannelID_t ChannelMapUBooNEAlg::PlaneWireToChannel(unsigned int plane,
						     unsigned int wire,
						     unsigned int tpc,
						     unsigned int cstat) const
  {
    // This is the actual lookup part - first make sure coordinates are legal
    if(tpc   < fNTPC[cstat] &&
       plane < fWiresPerPlane[cstat][tpc].size() &&
       wire  < fWiresPerPlane[cstat][tpc][plane] ){
      // if the channel has legal coordinates, its ID is given by the wire
      // number above the number of wires in lower planes, tpcs and cryostats
      return fPlaneBaselines[cstat][tpc][plane] + wire;
    }
    else{  
      // if the coordinates were bad, throw an exception
      throw cet::exception("ChannelMapUBooNEAlg") << "NO CHANNEL FOUND for cryostat,tpc,plane,wire: "
						    << cstat << "," << tpc << "," << plane << "," << wire;
    }
    
    // made it here, that shouldn't happen, return raw::InvalidChannelID
    mf::LogWarning("ChannelMapUBooNEAlg") << "should not be at the point in the function, returning "
					    << "invalid channel";
    return raw::InvalidChannelID;

  }


  //----------------------------------------------------------------------------
  SigType_t ChannelMapUBooNEAlg::SignalType(raw::ChannelID_t const channel) const
  {

    // still assume one cryostat for now -- faster
    unsigned int nChanPerTPC = fNchannels/fNTPC[0];
    // casting wil trunc towards 0 -- faster than floor
    unsigned int tpc = channel / nChanPerTPC;  
    //need number of planes to know Collection 
    unsigned int PlanesThisTPC = fNPlanes[0][tpc];
    
      
    
    SigType_t sigt = geo::kMysteryType;
    if(      (channel >= fFirstChannelInThisPlane[0][tpc][0]) &&
             (channel <  fFirstChannelInNextPlane[0][tpc][PlanesThisTPC-2])    ){ sigt = geo::kInduction; }
    else if( (channel >= fFirstChannelInThisPlane[0][tpc][PlanesThisTPC-1]) &&
             (channel <  fFirstChannelInNextPlane[0][tpc][PlanesThisTPC-1])    ){ sigt = geo::kCollection; }
    else
      mf::LogWarning("BadChannelSignalType") << "Channel " << channel
					     << " not given signal type." << std::endl;
            

    return sigt;
  }


  //----------------------------------------------------------------------------
  View_t ChannelMapUBooNEAlg::View(raw::ChannelID_t const channel) const
  {

    // still assume one cryostat for now -- faster
    unsigned int nChanPerTPC = fNchannels/fNTPC[0];
    // casting wil trunc towards 0 -- faster than floor
    unsigned int tpc = channel / nChanPerTPC; 


    View_t view = geo::kUnknown; 

    if(      (channel >= fFirstChannelInThisPlane[0][tpc][0]) &&
             (channel <  fFirstChannelInNextPlane[0][tpc][0])    ){ view = geo::kU; }
    else if( (channel >= fFirstChannelInThisPlane[0][tpc][1]) &&
             (channel <  fFirstChannelInNextPlane[0][tpc][1])    ){ view = geo::kV; }
    else if( (channel >= fFirstChannelInThisPlane[0][tpc][2]) &&
             (channel <  fFirstChannelInNextPlane[0][tpc][2])    ){ view = geo::kZ; }
    else
      mf::LogWarning("BadChannelSignalType") << "Channel " << channel
					     << " not given view type.";

    return view;
  }  

  //----------------------------------------------------------------------------
  std::set<View_t> const& ChannelMapUBooNEAlg::Views() const
  {
    return fViews;
  }

  //----------------------------------------------------------------------------
  std::set<PlaneID> const& ChannelMapUBooNEAlg::PlaneIDs() const
  {
    return fPlaneIDs;
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
  opdet::UBOpticalChannelCategory_t ChannelMapUBooNEAlg::GetChannelType( unsigned int opChannel ) const {
    auto it = fChannelCategory.find( opChannel );
    if ( it!=fChannelCategory.end() )
      return it->second;
    return opdet::Uncategorized;
  }

  //----------------------------------------------------------------------------
  opdet::UBOpticalChannelGain_t ChannelMapUBooNEAlg::GetChannelGain( unsigned int opChannel ) const {
    auto it = fChannelGain.find( opChannel );
    if ( it!=fChannelGain.end() )
      return it->second;
    throw std::runtime_error( "ChannelMapUBooNEAlg::GetChannelGain : Could not find channel type (low,high,logic)" );
  }

  //----------------------------------------------------------------------------
  void  ChannelMapUBooNEAlg::GetLogicChannelList( std::vector< unsigned int >& channels ) const {
    channels.resize( fLogicChannels.size() );
    std::copy( fLogicChannels.begin(), fLogicChannels.end(), channels.begin() );
  }

  //----------------------------------------------------------------------------
  void ChannelMapUBooNEAlg::LoadOpticalMapData( fhicl::ParameterSet const& pset ) {
    fNOpDets = pset.get< unsigned int >("numberOfDetectors");

    fNReadoutChannels = 0;

    // read in opdet to channel map, store
    for (unsigned int iop=0; iop<fNOpDets; iop++) {
      char entryname[50];
      sprintf(entryname,"OpDet%d_channels",iop);
      std::vector< unsigned int > chinput =  pset.get< std::vector<unsigned int> >( entryname );
      fPMT2channels[ iop ] = chinput;
      fNReadoutChannels += chinput.size();
      std::cout << entryname << ": [";
      for (std::vector<unsigned int>::iterator it_ch=chinput.begin(); it_ch!=chinput.end(); it_ch++) {
	fChannel2pmt[ *it_ch ] = iop;
	std::cout << *it_ch << ",";
      }
      std::cout << "]" << std::endl;
    }
    
    // read in channel types
    for ( unsigned int icat=0; icat<(unsigned int)opdet::NumUBOpticalChannelCategories; icat++ ) {
      std::vector< unsigned int > cat_channels;
      try {
	cat_channels = pset.get< std::vector<unsigned int> >( opdet::UBOpChannelEnumName( (opdet::UBOpticalChannelCategory_t)icat ) );
      }
      catch (...) {
	continue;
      }

      // save category channels
      fCategoryChannels[ (opdet::UBOpticalChannelCategory_t)icat ] = std::set< unsigned int >( cat_channels.begin(), cat_channels.end() );
      // save gain channel 
      opdet::UBOpticalChannelGain_t chtype = opdet::GetUBTypeFromCategory( (opdet::UBOpticalChannelCategory_t)icat );
      if ( fGainChannels.find( chtype )==fGainChannels.end() )
	fGainChannels[ chtype ] = std::set< unsigned int >( cat_channels.begin(), cat_channels.end() ); // new set
      else {
	for ( auto v : cat_channels )
	  fGainChannels[ chtype ].insert( v ); // append to set
      }

      for ( auto v : cat_channels ) {
	fChannelGain[ v ] = chtype;
	if (chtype==opdet::LogicChannel)
	  fLogicChannels.insert( v );
      }
    }
  }
  
} // namespace

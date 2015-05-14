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
  ChannelMapUBooNEAlg::ChannelMapUBooNEAlg(fhicl::ParameterSet const& sort_pars, fhicl::ParameterSet const& helper_pset )
    : ChannelMapAlg(), fSorter(geo::GeoObjectSorterStandard(sort_pars))
  {
    // parameter set will come from UBooNEGomeotryHelper service
    LoadOpticalMapData( helper_pset );
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
    return fNSplits*NOpDets;
  }

  //----------------------------------------------------------------------------
  unsigned int ChannelMapUBooNEAlg::NOpHardwareChannels(unsigned int opDet) const
  {
    // four channels per optical detector: LoA, LoB, HiA, HiB
    return fNSplits;
  }

  //----------------------------------------------------------------------------
  unsigned int ChannelMapUBooNEAlg::OpChannel(unsigned int pmtID, unsigned int copynum) const
  {
    // converts OpDet and Gain Channel
    // 0=low A, 1= low B, 2 = high A, 3 = high B
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
    unsigned int copyid = opChannel % 100;
    return copyid;
  }
  std::string ChannelMapUBooNEAlg::GetSplitGain( unsigned int isplit ) const 
  {
    return fSplitGains.at( isplit );
  }

  //----------------------------------------------------------------------------
  bool ChannelMapUBooNEAlg::IsValidOpChannel(unsigned int opChannel, unsigned int NOpDets) const {
    unsigned int pmtid = OpDetFromOpChannel( opChannel );
    if ( pmtid>=36 )
      return false;
    if ( opChannel>=336 )
      return false;

    return true;
  }

  //----------------------------------------------------------------------------
  unsigned int ChannelMapUBooNEAlg::NOpLogicChannels() const {
    return fLogicChannelList.size();
  }
  
  //----------------------------------------------------------------------------
  void ChannelMapUBooNEAlg::GetLogicChannelList( std::vector< unsigned int >& channels ) const {
    channels = fLogicChannelList; // copy
  }

  //----------------------------------------------------------------------------
  void ChannelMapUBooNEAlg::LoadOpticalMapData( fhicl::ParameterSet const& pset ) {
    fNOpDets = pset.get< unsigned int >("numberOfDetectors");
    fNSplits = pset.get<  unsigned int >("numberOfSplits");
    fSplitGains = pset.get< std::vector< std::string > >( "splitGains" );
    if ( fSplitGains.size()!=fNSplits ) {
      throw cet::exception("numberOfSplits does not match the number of specified gains in splitGains");
    }
    for (unsigned int i=0; i<fNSplits; i++) {
      if ( fSplitGains.at(i)!="HighGain" && fSplitGains.at(i)!="LowGain" ) {
	throw cet::exception( __FUNCTION__ ) << " invalid splitGains option.  Either 'HighGain' or 'LowGain'." << std::endl;
      }
    }
    
    // read in low gain, high gain ranges
    std::vector< unsigned int > lorange_input = pset.get< std::vector<unsigned int> >( "ReadoutLowChannelRange" );
    std::vector< unsigned int > hirange_input = pset.get< std::vector<unsigned int> >( "ReadoutHighChannelRange" );
    for (unsigned int i=0; i<(unsigned int)lorange_input.size(); i+=2 ) {
      unsigned int lo[2] = { lorange_input.at(i), lorange_input.at(i+1) };
      std::vector< unsigned int > vlo( lo, lo+2 );
      fLowgain_channel_ranges.push_back( vlo );
    }
    for (unsigned int i=0; i<(unsigned int)hirange_input.size(); i+=2 ) {
      unsigned int hi[2] = { hirange_input.at(i), hirange_input.at(i+1) };
      std::vector< unsigned int > vhi( hi, hi+2 );
      fHighgain_channel_ranges.push_back( vhi );
    }

    // read in opdet to channel map, store
    for (unsigned int iop=0; iop<fNOpDets; iop++) {
      char entryname[50];
      sprintf(entryname,"OpDet%d_channels",iop);
      std::vector< unsigned int > chinput =  pset.get< std::vector<unsigned int> >( entryname );
      fPMT2channels[ iop ] = chinput;
      std::cout << entryname << ": [";
      for (std::vector<unsigned int>::iterator it_ch=chinput.begin(); it_ch!=chinput.end(); it_ch++) {
	fChannel2pmt[ *it_ch ] = iop;
	std::cout << *it_ch << ",";
      }
      std::cout << "]" << std::endl;
    }
    
    // read in active logic channel lists
    fLogicChannelList = pset.get< std::vector<unsigned int> >( "LogicChannelList" );
  }

} // namespace

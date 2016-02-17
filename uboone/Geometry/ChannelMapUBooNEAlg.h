////////////////////////////////////////////////////////////////////////
/// \file  ChannelMapUBooNEAlg.h
/// \brief MicroBooNE optical detector channel mapping
///
/// Terminology:
///  OpDet: the physical PMT or light guide paddle
///  HardwareChannel: The gain copy number. 0=Low A, 1 = High A, 2=Low B, 3 = High B  
///  OpChannel: The readout channel number
///
/// \version $Id:  $
/// \author taritree@mit.edu
////////////////////////////////////////////////////////////////////////
#ifndef GEO_CHANNELMAP_UBOONE_ALG_H
#define GEO_CHANNELMAP_UBOONE_ALG_H

#include <vector>
#include <set>
#include <iostream>

#include "fhiclcpp/ParameterSet.h"

#include "larcore/Geometry/ChannelMapStandardAlg.h" // larcore
#include "larcore/Geometry/GeoObjectSorterStandard.h" //  larcore
#include "uboone/Geometry/UBOpChannelTypes.h" // uboonecode

namespace geo {

  class ChannelMapUBooNEAlg : public ChannelMapStandardAlg {

  public:

    ChannelMapUBooNEAlg(fhicl::ParameterSet const& pset, fhicl::ParameterSet const& sortingParameters );
    ~ChannelMapUBooNEAlg();

    // Below inherited from ChannelMapStandardAlg
    /* ChannelMapStandardAlg(fhicl::ParameterSet const& p); */
    
    /* void                     Initialize( GeometryData_t& geodata ) override; */
    /* void                     Uninitialize(); */
    /* std::vector<WireID>      ChannelToWire(raw::ChannelID_t channel)     const; */
    /* unsigned int             Nchannels()                                 const; */

    /* //@{ */
    /* virtual double WireCoordinate */
    /*   (double YPos, double ZPos, geo::PlaneID const& planeID) const override; */
    /* virtual double WireCoordinate(double YPos, double ZPos, */
    /*                              unsigned int PlaneNo, */
    /*                              unsigned int TPCNo, */
    /*                              unsigned int cstat) const */
    /*   { return WireCoordinate(YPos, ZPos, geo::PlaneID(cstat, TPCNo, PlaneNo)); } */
    /* //@} */
    
    /* //@{ */
    /* virtual WireID NearestWireID */
    /*   (const TVector3& worldPos, geo::PlaneID const& planeID) const override; */
    /* virtual WireID NearestWireID(const TVector3& worldPos, */
    /*                              unsigned int    PlaneNo, */
    /*                              unsigned int    TPCNo, */
    /*                              unsigned int    cstat) const override */
    /*   { return NearestWireID(worldPos, geo::PlaneID(cstat, TPCNo, PlaneNo)); } */
    /* //@} */
    
    /* //@{ */
    /* virtual raw::ChannelID_t PlaneWireToChannel */
    /*   (geo::WireID const& wireID) const override; */
    /* virtual raw::ChannelID_t PlaneWireToChannel(unsigned int plane, */
    /*                                             unsigned int wire, */
    /*                                             unsigned int tpc, */
    /*                                             unsigned int cstat) const override */
    /*   { return PlaneWireToChannel(geo::WireID(cstat, tpc, plane, wire)); } */
    /* //@} */
    
    /* virtual View_t                   View( raw::ChannelID_t const channel )       const override; */
    /* virtual SigType_t                SignalType( raw::ChannelID_t const channel ) const override; */
    /* virtual std::set<View_t>  const& Views()                                      const override; */
    /* virtual std::set<PlaneID> const& PlaneIDs()                                   const override; */


    // ----------------------------------------------------------------------------
    // OPDET MAPPING: What we implement
  protected:
    unsigned int fNOpDets;
    unsigned int fNReadoutChannels;
    unsigned int fMaxOpChannel;
    std::set< unsigned int > fReadoutChannelSet;
    std::map< unsigned int, unsigned int > fChannel2pmt;  // readout channel to opdet(pmt) id
    std::map< unsigned int, std::vector< unsigned int > > fPMT2channels; // opdet(pmt) id to readout channel

  public:
    unsigned int NOpChannels(unsigned int NOpDets)                        const;
    unsigned int MaxOpChannel(unsigned int NOpDets)                        const;
    unsigned int NOpHardwareChannels(unsigned int opDet)                  const;
    unsigned int OpChannel(unsigned int detNum, unsigned int channel = 0) const;
    unsigned int OpDetFromOpChannel(unsigned int opChannel)               const;
    unsigned int HardwareChannelFromOpChannel(unsigned int opChannel)     const;
    bool         IsValidOpChannel(unsigned int opChannel, unsigned int NOpDets) const;
    void         LoadOpticalMapData( fhicl::ParameterSet const& pset );
    
  };


}
#endif // GEO_CHANNELMAPSTANDARDALG_H


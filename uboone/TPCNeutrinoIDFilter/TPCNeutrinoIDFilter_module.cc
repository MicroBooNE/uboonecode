/**
 * @file    TPCNeutrinoIDFilter_module.cc
 * @brief   Module to filter neutrino candidate events based on TPC topology
 * @authors aschu@fnal.gov
 * 
 ******************************************************************************/


// Framework includes
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/View.h"
#include "art/Persistency/Common/Ptr.h"
#include "art/Persistency/Common/PtrVector.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "art/Framework/Core/FindMany.h"
#include "art/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "Geometry/Geometry.h"
#include "SimulationBase/MCTruth.h"
#include "MCBase/MCShower.h"
#include "MCBase/MCStep.h"
#include "SimulationBase/MCFlux.h"
#include "Simulation/SimChannel.h"
#include "Simulation/AuxDetSimChannel.h"
#include "AnalysisBase/Calorimetry.h"
#include "AnalysisBase/ParticleID.h"
#include "RawData/RawDigit.h"
#include "RawData/BeamInfo.h"
#include "Utilities/LArProperties.h"
#include "Utilities/AssociationUtil.h"
#include "Utilities/DetectorProperties.h"
#include "SummaryData/POTSummary.h"
#include "MCCheater/BackTracker.h"
#include "RecoBase/Track.h"
#include "RecoBase/Shower.h"
#include "RecoBase/Cluster.h"
#include "RecoBase/Hit.h"
#include "RecoBase/EndPoint2D.h"
#include "RecoBase/Vertex.h"
#include "RecoBase/OpFlash.h"
#include "SimpleTypesAndConstants/geo_types.h"
#include "RecoObjects/BezierTrack.h"
#include "RecoAlg/TrackMomentumCalculator.h"
#include "AnalysisBase/CosmicTag.h"
#include "AnalysisBase/FlashMatch.h"
	

#include <cstddef> // std::ptrdiff_t
#include <cstring> // std::memcpy()
#include <vector>
#include <map>
#include <iterator> // std::begin(), std::end()
#include <string>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <functional> // std::mem_fun_ref
#include <typeinfo>
#include <memory> // std::unique_ptr<>

#include "TTree.h"
#include "TTimeStamp.h"


class TPCNeutrinoIDFilter;
class TPCNeutrinoIDFilter : public art::EDFilter {

public:

   explicit TPCNeutrinoIDFilter(fhicl::ParameterSet const& pset);
   virtual ~TPCNeutrinoIDFilter();

   bool filter(const art::Event& evt);

private:

   std::string fTrackModuleLabel;
   std::string fVertexModuleLabel;
   std::string fCosmicTaggerAssocLabel;

}; // class microboone::TPCNeutrinoIDFilter


TPCNeutrinoIDFilter::TPCNeutrinoIDFilter(fhicl::ParameterSet const& pset) :
   fTrackModuleLabel         (pset.get< std::string >("TrackModuleLabel")),
   fVertexModuleLabel        (pset.get< std::string >("VertexModuleLabel")),
   fCosmicTaggerAssocLabel   (pset.get< std::string >("CosmicTaggerAssocLabel"))
{
} // microboone::TPCNeutrinoIDFilter::TPCNeutrinoIDFilter()

bool TPCNeutrinoIDFilter::filter(const art::Event& evt)
{

  bool pass = false;

  double trkstartx = 0;
  double trkstarty = 0;
  double trkstartz = 0;
  double trkendx = 0;
  double trkendy = 0;
  double trkendz = 0;

  // * vertices
  art::Handle< std::vector<recob::Vertex> > vertexListHandle;
  std::vector< art::Ptr<recob::Vertex> > vertexlist;
  if (evt.getByLabel(fVertexModuleLabel,vertexListHandle))
    art::fill_ptr_vector(vertexlist, vertexListHandle);
  
  // * tracks
  art::Handle< std::vector<recob::Track> > trackListHandle;
  std::vector< art::Ptr<recob::Track> > tracklist;
  if (evt.getByLabel(fTrackModuleLabel,trackListHandle))
     art::fill_ptr_vector(tracklist, trackListHandle);

  size_t NTracks = tracklist.size();
    
  for(size_t iTrk=0; iTrk < NTracks; ++iTrk) {//loop over tracks
      
     //Cosmic Tagger information
     //art::FindManyP<anab::CosmicTag> fmct(trackListHandle,evt,fCosmicTaggerAssocLabel);
     //if (fmct.isValid()) float cosmicscore = fmct.at(0)->CosmicScore();

     art::Ptr<recob::Track> ptrack(trackListHandle, iTrk);
     const recob::Track& track = *ptrack;

     TVector3 pos, end;

     int ntraj = 0;
     if(fTrackModuleLabel == "beziertracker") {
        trkf::BezierTrack btrack(*ptrack);
        ntraj = btrack.NSegments();
        if(ntraj > 0) {
           double xyz[3];
           btrack.GetTrackPoint(0,xyz);
           pos.SetXYZ(xyz[0],xyz[1],xyz[2]);
           btrack.GetTrackPoint(1,xyz);
           end.SetXYZ(xyz[0],xyz[1],xyz[2]);
        }
     }
     else {   //use the normal methods for other kinds of tracks
        ntraj = track.NumberTrajectoryPoints();
        if(ntraj > 0) {
           pos = track.Vertex();
           end = track.End();
        }
     }

     std::cout << "trajectory points " << ntraj << std::endl;  
     if(ntraj > 0) {
        trkstartx = pos.X();
        trkstarty = pos.Y();
        trkstartz = pos.Z();
        trkendx = end.X();
        trkendy = end.Y();
        trkendz = end.Z();
     }

     std::cout << trkstartx << "\t" << trkstarty << "\t" << trkstartz << std::endl;
     std::cout << trkendx << "\t" << trkendy << "\t" << trkendz << std::endl;
  }

  return pass;

} // microboone::TPCNeutrinoIDFilter::filter()


namespace microboone{

  DEFINE_ART_MODULE(TPCNeutrinoIDFilter)

}

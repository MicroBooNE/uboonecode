////////////////////////////////////////////////////////////////////////
// Class:       CosmicFlashTagger
// Plugin Type: producer (art v2_05_00)
// File:        CosmicFlashTagger_module.cc
//
// Generated at Wed Nov 30 09:45:20 2016 by Marco Del Tutto using cetskelgen
// from cetlib version v1_21_00.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "canvas/Persistency/Common/FindManyP.h"

#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/AnalysisBase/CosmicTag.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larcore/Geometry/Geometry.h"
#include "larcore/Geometry/CryostatGeo.h"
#include "larcore/Geometry/PlaneGeo.h"
#include "larcore/Geometry/OpDetGeo.h"
#include "uboone/Geometry/UBOpReadoutMap.h"

#include "lardataobj/RecoBase/PFParticle.h"
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"

#include "uboone/LLSelectionTool/OpT0Finder/Base/OpT0FinderTypes.h"
#include "uboone/LLBasicTool/GeoAlgo/GeoTrajectory.h"
#include "uboone/LLSelectionTool/OpT0Finder/Base/FlashMatchManager.h"
#include "uboone/LLSelectionTool/OpT0Finder/Algorithms/LightPath.h"
#include "uboone/LLSelectionTool/OpT0Finder/Algorithms/PhotonLibHypothesis.h"

#include "uboone/CosmicTagging/Algo/IncompatibilityChecker.h"

#include "TString.h"
#include "TTree.h"

#include <memory>

class CosmicFlashTagger;

namespace cosmic {
      class BeamFlashTrackMatchTagger;
}

class CosmicFlashTagger : public art::EDProducer {
public:
  explicit CosmicFlashTagger(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  CosmicFlashTagger(CosmicFlashTagger const &) = delete;
  CosmicFlashTagger(CosmicFlashTagger &&) = delete;
  CosmicFlashTagger & operator = (CosmicFlashTagger const &) = delete;
  CosmicFlashTagger & operator = (CosmicFlashTagger &&) = delete;

  // Required functions.
  void produce(art::Event & e) override;

private:

  int  GetTrajectory(std::vector<art::Ptr<recob::Track>> track, double xoffset, ::geoalgo::Trajectory &);
  void AddFlashPosition(::flashana::Flash_t &);

  ::flashana::FlashMatchManager       _mgr;
  ::flashana::IncompatibilityChecker  _incompChecker;
  std::vector<::flashana::Flash_t>    beam_flashes;
  std::vector<art::Ptr<recob::Track>> track_v;

  const anab::CosmicTagID_t TAGID_BEAM_INCOMPATIBLE = anab::CosmicTagID_t::kFlash_BeamIncompatible;
  const anab::CosmicTagID_t TAGID_NOT_TAGGED        = anab::CosmicTagID_t::kNotTagged;

  const std::vector<float> endPt1 = {-9999., -9999., -9999.};
  const std::vector<float> endPt2 = {-9999., -9999., -9999.};

  std::string _pfp_producer;
  std::string _track_producer;
  std::string _opflash_producer_beam;
  std::string _opflash_producer_cosmic;
  std::string _trigger_producer;
  std::string _beam_window_start_BNB;
  std::string _beam_window_end_BNB;
  double _flash_trange_start;
  double _flash_trange_end;
  int    _min_trj_pts;
  double _min_track_length;
  bool _debug;
  double fDriftVelocity;

  TTree* _tree1;
  int _run, _subrun, _event, _matchid;
  int _n_beam_flashes, _n_pfp;
  std::vector<std::vector<double>> _beam_flash_spec, _pfp_hypo_spec;
  std::vector<double> _beam_flash_time;
  int _beam_flash_exists; // 0 == no;  1 == yes
};


CosmicFlashTagger::CosmicFlashTagger(fhicl::ParameterSet const & p)
{
  _pfp_producer            = p.get<std::string>("PFParticleProducer");
  _track_producer          = p.get<std::string>("TrackProducer");
  _opflash_producer_beam   = p.get<std::string>("BeamOpFlashProducer");
  _opflash_producer_cosmic = p.get<std::string>("CosmicOpFlashProducer");
  _trigger_producer        = p.get<std::string>("TriggerProducer");
  _flash_trange_start      = p.get<double>     ("FlashVetoTimeStart");
  _flash_trange_end        = p.get<double>     ("FlashVetoTimeEnd");
  _min_trj_pts             = p.get<int>        ("MinimumNumberOfTrajectoryPoints");
  _beam_window_start_BNB   = p.get<int>        ("BeamWindowStartBNB");
  _beam_window_end_BNB     = p.get<int>        ("BeamWindowEndBNB");
  _min_track_length        = p.get<double>     ("MinimumTrackLength");
  _debug                   = p.get<bool>       ("DebugMode");

  _mgr.Configure(p.get<flashana::Config_t>("FlashMatchConfig"));
  _incompChecker.Configure(p.get<fhicl::ParameterSet>("IncompCheckConfig"));
  if(_debug) _incompChecker.PrintConfig();

  // Use '_detp' to find 'efield' and 'temp'
  auto const* _detp = lar::providerFrom<detinfo::DetectorPropertiesService>();
  double efield = _detp -> Efield();
  double temp   = _detp -> Temperature();
  // Determine the drift velocity from 'efield' and 'temp'
  fDriftVelocity = _detp -> DriftVelocity(efield,temp); 
  if (_debug) std::cout << "Using drift velocity = " << fDriftVelocity << " cm/us, with E = " << efield << ", and T = " << temp << std::endl;

  if (_debug) {
    art::ServiceHandle<art::TFileService> fs;
    _tree1 = fs->make<TTree>("tree","");
    _tree1->Branch("run",&_run,"run/I");
    _tree1->Branch("subrun",&_subrun,"subrun/I");
    _tree1->Branch("event",&_event,"event/I");
    _tree1->Branch("n_beam_flashes",&_n_beam_flashes,"n_beam_flashes/I");
    _tree1->Branch("beam_flash_spec","std::vector<std::vector<double>>",&_beam_flash_spec);
    _tree1->Branch("beam_flash_time","std::vector<double>",&_beam_flash_time);
    _tree1->Branch("beam_flash_exists",&_beam_flash_exists,"beam_flash_exists/I");
    _tree1->Branch("n_pfp",&_n_pfp,"n_pfp/I");
    _tree1->Branch("pfp_hypo_spec","std::vector<std::vector<double>>",&_pfp_hypo_spec);
  }

  produces< std::vector<anab::CosmicTag>>();  
  produces< art::Assns<anab::CosmicTag,   recob::Track>>();  
  produces< art::Assns<recob::PFParticle, anab::CosmicTag>>();
}

void CosmicFlashTagger::produce(art::Event & e)
{
  if(_debug) {
    std::cout << "CosmicFlashTagger::produce starts." << std::endl;
    std::cout << "This is run | subrun | event: " << e.id().run() << " | " << e.id().subRun() << " | " << e.id().event() << std::endl;
    _run    = e.id().run();
    _subrun = e.id().subRun();
    _event  = e.id().event();
  }

  // Instantiate the output
  std::unique_ptr< std::vector< anab::CosmicTag>>                  cosmicTagTrackVector      (new std::vector<anab::CosmicTag>);
  std::unique_ptr< art::Assns<anab::CosmicTag, recob::Track>>      assnOutCosmicTagTrack     (new art::Assns<anab::CosmicTag, recob::Track>);
  std::unique_ptr< art::Assns<recob::PFParticle, anab::CosmicTag>> assnOutCosmicTagPFParticle(new art::Assns<recob::PFParticle, anab::CosmicTag>);

  // Reset the flash match manager
  _mgr.Reset();

  ::art::ServiceHandle<geo::Geometry> geo;

  // Get Beam Flashes from the ART event
  ::art::Handle<std::vector<recob::OpFlash> > beamflash_h;
  e.getByLabel(_opflash_producer_beam,beamflash_h);
  if( !beamflash_h.isValid() || beamflash_h->empty() ) {
    std::cerr << "Don't have good flashes." << std::endl;
    return;
  }

  // Get PFParticles and map to tracks from the ART event
  lar_pandora::TrackVector         trackVector; 
  lar_pandora::PFParticlesToTracks PFPtoTracks; 
  lar_pandora::LArPandoraHelper::CollectTracks(e, _pfp_producer, trackVector, PFPtoTracks);

  // Get Tracks from the ART event
  ::art::Handle<std::vector<recob::Track> > track_h;
  e.getByLabel(_track_producer,track_h);
  if( !track_h.isValid() || track_h->empty() )  {
    std::cerr << "Don't have tracks, or they are not valid." << std::endl;
    throw std::exception();
  }

  // Loop through beam flashes 
  _n_beam_flashes = 0;
  _beam_flash_exists = 0;
  beam_flashes.clear();
  for (size_t n = 0; n < beamflash_h->size(); n++) {

    auto const& flash = (*beamflash_h)[n];

    if(flash.Time() > _beam_window_start_BNB && flash.Time() < _beam_window_end_BNB){
      _beam_flash_exists = 1;
    }
    if(flash.Time() < _flash_trange_start || _flash_trange_end < flash.Time()) {
      if(_debug) std::cout << "Flash is in veto region (flash time is " << flash.Time() << "). Continue." << std::endl;
      continue;
    } 
   
    _n_beam_flashes++;
    if (_debug) {
      _beam_flash_spec.resize(_n_beam_flashes);
      _beam_flash_spec[_n_beam_flashes-1].resize(geo->NOpDets());
      for (unsigned int i = 0; i < geo->NOpDets(); i++) {
        unsigned int opdet = geo->OpDetFromOpChannel(i);
        _beam_flash_spec[_n_beam_flashes-1][opdet] = flash.PE(i);
      }
      _beam_flash_time.resize(_n_beam_flashes);
      _beam_flash_time[_n_beam_flashes-1] = flash.Time();
    }

    // Construct a Flash_t
    ::flashana::Flash_t f;
    f.x = f.x_err = 0;
    f.y = flash.YCenter();
    f.z = flash.ZCenter();
    f.y_err = flash.YWidth();
    f.z_err = flash.ZWidth();
    f.pe_v.resize(geo->NOpDets());
    f.pe_err_v.resize(geo->NOpDets());
    for (unsigned int i = 0; i < f.pe_v.size(); i++) {
      unsigned int opdet = geo->OpDetFromOpChannel(i);
      f.pe_v[opdet] = flash.PE(i);
      f.pe_err_v[opdet] = sqrt(flash.PE(i));
    }
    f.time = flash.Time();
    //f.idx = flash_id;
    beam_flashes.resize(_n_beam_flashes);
    beam_flashes[_n_beam_flashes-1] = f;

  } // end of flash loop

  if (_debug) std::cout << "Number of beam flashes in this event: " << beam_flashes.size() << std::endl;

  _n_pfp = 0;

  // --- Loop over PFParticles
  for (lar_pandora::PFParticlesToTracks::iterator it = PFPtoTracks.begin(); it != PFPtoTracks.end(); ++it) {

    bool beamIncompatible = false;
    art::Ptr<recob::PFParticle> pfParticle;
      
    // --- Loop over beam flashes ---
    for (unsigned int bf = 0; bf < beam_flashes.size(); bf++) {

      // Get the PFParticle
      pfParticle = it->first;
      if(_debug){
        std::cout << "This is PFP with ID " << pfParticle->Self() << std::endl;
        _n_pfp++;
        _pfp_hypo_spec.resize(_n_pfp);
      }

      // Get the tracks associated with this PFParticle
      lar_pandora::TrackVector track_v = it->second;

      // Get the beam flash
      ::flashana::Flash_t flashBeam = beam_flashes[bf];

      // Calculate x offset, assuming this track caused this beam flash
      double Xoffset = flashBeam.time * 0.1114359;
      if(_debug) std::cerr << "Xoffset is " << Xoffset << std::endl;

      // Get trajectory (1 trajectory for all these tracks)
      ::geoalgo::Trajectory trkTrj;
      int statusCode = this->GetTrajectory(track_v, Xoffset, trkTrj);
      if(statusCode != 0) break;

      // From the trajectory construct a QCluster
      auto qcluster = ((flashana::LightPath*)(_mgr.GetCustomAlgo("LightPath")))->FlashHypothesis(trkTrj);

      // From the QCluster get the flash hypothesis using registered FlashHypothesis algorithm
      flashana::Flash_t flashHypo;
      flashHypo.pe_v.resize(geo->NOpDets());
      ((flashana::PhotonLibHypothesis*)(_mgr.GetAlgo(flashana::kFlashHypothesis)))->FillEstimate(qcluster,flashHypo);
      
      if(_debug) {
        _pfp_hypo_spec[_n_pfp-1] = flashHypo.pe_v;
        std::cout << "***The beam flash has Z = " << flashBeam.z << " +- " << flashBeam.z_err << std::endl;
        this->AddFlashPosition(flashHypo);
        std::cout << "***The hypo flash has Z = " << flashHypo.z << " +- " << flashHypo.z_err << std::endl;
      }
      // CORE FUNCTION: Check if this beam flash and this flash hypothesis are incompatible
      bool areIncompatible = _incompChecker.CheckIncompatibility(flashBeam,flashHypo); 
      if(_debug) std::cout << "For this PFP: " << (areIncompatible ? "are INcompatible" : "are compatible") << std::endl;
      
      if (areIncompatible == false) break;
      else if (areIncompatible && bf == beam_flashes.size() - 1) {
        // This PFP is not compatible with any of the beam flashes
        beamIncompatible = true;

      } else if (!areIncompatible && bf == beam_flashes.size() - 1) {
        // Can't tell anything for this PFP
        beamIncompatible = false;
      }
      
    } // end of beam flash loop

    if (beamIncompatible) {
      float cosmicScore = 1;
      cosmicTagTrackVector->emplace_back(endPt1, endPt2, cosmicScore, TAGID_BEAM_INCOMPATIBLE);
      util::CreateAssn(*this, e, *cosmicTagTrackVector, track_v,    *assnOutCosmicTagTrack );
      util::CreateAssn(*this, e, *cosmicTagTrackVector, pfParticle, *assnOutCosmicTagPFParticle);
    } else {
      float cosmicScore = 0;
      cosmicTagTrackVector->emplace_back(endPt1, endPt2, cosmicScore, TAGID_NOT_TAGGED);
      util::CreateAssn(*this, e, *cosmicTagTrackVector, track_v,    *assnOutCosmicTagTrack );
      util::CreateAssn(*this, e, *cosmicTagTrackVector, pfParticle, *assnOutCosmicTagPFParticle);
    }

  } // end of PFP loop

  _tree1->Fill();

  e.put(std::move(cosmicTagTrackVector));
  e.put(std::move(assnOutCosmicTagTrack));
  e.put(std::move(assnOutCosmicTagPFParticle));

  if(_debug) std::cout << "CosmicFlashTagger::produce ends." << std::endl;
}


//______________________________________________________________________________________________________________________________________
int CosmicFlashTagger::GetTrajectory(std::vector<art::Ptr<recob::Track>> track_v, double Xoffset, ::geoalgo::Trajectory &track_geotrj) {

  if (_debug) std::cout << "Creating trajectory for " << track_v.size() << " tracks." << std::endl;
  int statusCode;

  int totalPoints = 0;
  for (unsigned int trk = 0; trk < track_v.size(); trk++) {
    art::Ptr<recob::Track> trk_ptr = track_v.at(trk);
    totalPoints += trk_ptr->NumberTrajectoryPoints();
  }

  if (_debug) std::cout << "Number of points for this trajectory: " << totalPoints << std::endl;

  if (totalPoints <= _min_trj_pts) {
    statusCode = 1;
    return statusCode;
  }
  track_geotrj.resize(totalPoints,::geoalgo::Vector(0.,0.,0.));
  int trj_pt = 0;

  for (unsigned int trk = 0; trk < track_v.size(); trk++) {
    art::Ptr<recob::Track> trk_ptr = track_v.at(trk);
    for (size_t pt_idx=0; pt_idx < trk_ptr->NumberTrajectoryPoints(); ++pt_idx) {
      auto const& pt = trk_ptr->LocationAtPoint(pt_idx);
      track_geotrj[trj_pt][0] = pt[0] - Xoffset;
      track_geotrj[trj_pt][1] = pt[1];
      track_geotrj[trj_pt][2] = pt[2];
      trj_pt ++;
    }
  }
 
  statusCode = 0;
  return statusCode;
}

//______________________________________________________________________________________________________________________________________
void CosmicFlashTagger::AddFlashPosition(::flashana::Flash_t & flash) {

  // Reset variables
  double Ycenter = 0.;
  double Zcenter = 0.;
  double Ywidth  = -999.;
  double Zwidth  = -999.;
  double totalPE = 0.;
  double sumy = 0., sumz = 0., sumy2 = 0., sumz2 = 0.;

  std::vector<double> pePerOpDetId = flash.pe_v;

  for (unsigned int opdet = 0; opdet < pePerOpDetId.size(); opdet++) {

    // Get physical detector location for this opDet
    double PMTxyz[3] = {0.,0.,0.};
    ::art::ServiceHandle<geo::Geometry> geo;
    geo->OpDetGeoFromOpDet(opdet).GetCenter(PMTxyz);
 
    // Add up the position, weighting with PEs
    sumy    += pePerOpDetId[opdet]*PMTxyz[1];
    sumy2   += pePerOpDetId[opdet]*PMTxyz[1]*PMTxyz[1];
    sumz    += pePerOpDetId[opdet]*PMTxyz[2];
    sumz2   += pePerOpDetId[opdet]*PMTxyz[2]*PMTxyz[2];

    totalPE += pePerOpDetId[opdet];
  }

  Ycenter = sumy/totalPE;
  Zcenter = sumz/totalPE;

  // This is just sqrt(<x^2> - <x>^2)
  if ( (sumy2*totalPE - sumy*sumy) > 0. )
    Ywidth = std::sqrt(sumy2*totalPE - sumy*sumy)/totalPE;

  if ( (sumz2*totalPE - sumz*sumz) > 0. )
    Zwidth = std::sqrt(sumz2*totalPE - sumz*sumz)/totalPE;

  flash.y = Ycenter;
  flash.z = Zcenter;

  flash.y_err = Ywidth;
  flash.z_err = Zwidth;

}


DEFINE_ART_MODULE(CosmicFlashTagger)

////////////////////////////////////////////////////////////////////////
// Class:       MuCSTrackFinder
// Plugin Type: producer (art v2_04_00)
// File:        MuCSTrackFinder_module.cc
//
// Generated at Wed Nov 16 02:58:12 2016 by Marco Del Tutto using cetskelgen
// from cetlib version v1_20_00.
////////////////////////////////////////////////////////////////////////

#include "canvas/Persistency/Common/FindManyP.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include <memory>

#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/AnalysisBase/CosmicTag.h"
#include "larcore/Geometry/Geometry.h"
#include "larcore/Geometry/CryostatGeo.h"
#include "larcore/Geometry/PlaneGeo.h"
#include "larcore/Geometry/OpDetGeo.h"
#include "uboone/Geometry/UBOpReadoutMap.h"

#include "uboone/LLBasicTool/GeoAlgo/GeoTrajectory.h"
#include "uboone/LLSelectionTool/OpT0Finder/Base/FlashMatchManager.h"
#include "uboone/LLSelectionTool/OpT0Finder/Algorithms/LightPath.h"
#include "TString.h"
#include "TTree.h"

class MuCSTrackFinder;


class MuCSTrackFinder : public art::EDProducer {
public:
  explicit MuCSTrackFinder(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.
  
  // Plugins should not be copied or assigned.
  MuCSTrackFinder(MuCSTrackFinder const &) = delete;
  MuCSTrackFinder(MuCSTrackFinder &&) = delete;
  MuCSTrackFinder & operator = (MuCSTrackFinder const &) = delete;
  MuCSTrackFinder & operator = (MuCSTrackFinder &&) = delete;
  
  // Required functions.
  void produce(art::Event & e) override;
  
private:
  
  ::flashana::FlashMatchManager _mgr;
  
  std::vector<flashana::FlashMatch_t> _result;
  
  std::string _config_file;
  std::string _track_producer;
  std::string _opflash_producer_beam;
  std::string _opflash_producer_cosmic;
  std::string _trigger_producer;
  double _flash_trange_start;
  double _flash_trange_end;
  size_t _num_tracks;
  std::vector<double> _gain_correction;
  
  TTree* _tree1;
  int _run, _subrun, _event, _matchid;
  int _tpcid, _flashid;
  double _tpc_xmin, _qll_xmin;
  double _t0, _score;
  double _hypo_pe, _flash_pe;
  int _mucs_flash; // 0: there is no mucs flash; 1: there is mucs flash
  double _mucs_track_startx, _mucs_track_starty, _mucs_track_startz;
  double _mucs_track_endx, _mucs_track_endy, _mucs_track_endz;
  
  TTree* _tree2;
  std::vector<double> _flash_spec;
  std::vector<double> _hypo_spec;
  std::vector<double> _mucs_flash_spec;
};


MuCSTrackFinder::MuCSTrackFinder(fhicl::ParameterSet const & p)
// :
// Initialize member data here.
{
  _track_producer          = p.get<std::string>("TrackProducer");
  _opflash_producer_beam   = p.get<std::string>("BeamOpFlashProducer");
  _opflash_producer_cosmic = p.get<std::string>("CosmicOpFlashProducer");
  _trigger_producer        = p.get<std::string>("TriggerProducer");
  _flash_trange_start      = p.get<double>("FlashVetoTimeStart");
  _flash_trange_end        = p.get<double>("FlashVetoTimeEnd");
  _num_tracks              = p.get<size_t>("MaxTrackCount");
  _gain_correction         = p.get<std::vector<double> >("GainCorrection");
  
  ::art::ServiceHandle<geo::Geometry> geo;
  ::art::ServiceHandle<geo::UBOpReadoutMap> ub_geo;
  
  if(geo->NOpDets() != _gain_correction.size()) {
    std::cout << "GainCorrection array size is " << _gain_correction.size() << " != # OpDet " << geo->NOpDets() << std::endl;
    throw std::exception();
  }
  
  _mgr.Configure(p.get<flashana::Config_t>("FlashMatchConfig"));
  
  _flash_spec.resize(geo->NOpDets(),0.);
  _hypo_spec.resize(geo->NOpDets(),0.);
  _mucs_flash_spec.resize(geo->NOpDets(),-999.);
  
  art::ServiceHandle<art::TFileService> fs;
  
  _tree1 = fs->make<TTree>("tree","");
  _tree1->Branch("run",&_run,"run/I");
  _tree1->Branch("subrun",&_subrun,"subrun/I");
  _tree1->Branch("event",&_event,"event/I");
  _tree1->Branch("matchid",&_matchid,"matchid/I");
  _tree1->Branch("tpcid",&_tpcid,"tpcid/I");
  _tree1->Branch("flashid",&_flashid,"flashid/I");
  _tree1->Branch("tpc_xmin",&_tpc_xmin,"tpc_xmin/D");
  _tree1->Branch("qll_xmin",&_qll_xmin,"qll_xmin/D");
  _tree1->Branch("t0",&_t0,"t0/D");
  _tree1->Branch("score",&_score,"score/D");
  _tree1->Branch("hypo_pe",&_hypo_pe,"hypo_pe/D");
  _tree1->Branch("flash_pe",&_flash_pe,"flash_pe/D");
  _tree1->Branch("mucs_flash",&_mucs_flash,"mucs_flash/I");
  _tree1->Branch("mucs_track_startx",&_mucs_track_startx,"mucs_track_startx/D");
  _tree1->Branch("mucs_track_starty",&_mucs_track_starty,"mucs_track_starty/D");
  _tree1->Branch("mucs_track_startz",&_mucs_track_startz,"mucs_track_startz/D");
  _tree1->Branch("mucs_track_endx",&_mucs_track_endx,"mucs_track_endx/D");
  _tree1->Branch("mucs_track_endy",&_mucs_track_endy,"mucs_track_endy/D");
  _tree1->Branch("mucs_track_endz",&_mucs_track_endz,"mucs_track_endz/D");
  
  
  _tree2 = fs->make<TTree>("spectree","");
  _tree2->Branch("run",&_run,"run/I");
  _tree2->Branch("subrun",&_subrun,"subrun/I");
  _tree2->Branch("event",&_event,"event/I");
  _tree2->Branch("matchid",&_matchid,"matchid/I");
  _tree2->Branch("flash_spec","std::vector<double>",&_flash_spec);
  _tree2->Branch("hypo_spec","std::vector<double>",&_hypo_spec);
  _tree2->Branch("mucs_flash_spec","std::vector<double>",&_mucs_flash_spec);
  
}

void MuCSTrackFinder::produce(art::Event & e)
{
  std::cout << "MuCSTrackFinder starts." << std::endl;
  
  _mgr.Reset();
  _result.clear();
  
  ::art::ServiceHandle<geo::Geometry> geo;
  ::art::ServiceHandle<geo::UBOpReadoutMap> ub_geo;
  
  ::art::Handle<std::vector<recob::OpFlash> > beamflash_h;
  e.getByLabel(_opflash_producer_beam,beamflash_h);
  
  ::art::Handle<std::vector<recob::OpFlash> > cosmicflash_h;
  e.getByLabel(_opflash_producer_cosmic,cosmicflash_h);
  
  if( (!beamflash_h.isValid() || beamflash_h->empty()) && (!cosmicflash_h.isValid() || cosmicflash_h->empty()) ) {
    std::cerr << "Don't have good flashes." << std::endl;
    return;
  }
  
  _mucs_flash = 0;
  
  size_t flash_id=0;
  
  if(beamflash_h.isValid()) {
    for (size_t n = 0; n < beamflash_h->size(); n++) {
      
      auto const& flash = (*beamflash_h)[n];
      
      if (flash.Time() < -1.5 || flash.Time() > -0.5) {
        continue;
      }
      
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
        f.pe_v[opdet] = flash.PE(i) / _gain_correction[i];
        f.pe_err_v[opdet] = sqrt(flash.PE(i) / _gain_correction[i]);
      }
      f.time = flash.Time();
      f.idx = flash_id;
      ++flash_id;
      _mgr.Emplace(std::move(f));
    }
  }
  
  if(cosmicflash_h.isValid()) {
    for (size_t n = 0; n < cosmicflash_h->size(); n++) {
      
      auto const& flash = (*cosmicflash_h)[n];
      
      if (flash.Time() < -1.5 || flash.Time() > -0.5) {
        continue;
      }
      
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
        if(flash.PE(i) == 0.) {
          std::cout << "op det " << opdet << "has 0 pe for this flash" << std::endl;
          f.pe_v[opdet]=-1.;
          f.pe_err_v[opdet]=-1.;
        }else{
          std::cout << "op det " << opdet << "has " << flash.PE(i) / _gain_correction[i] / 0.424 << " pe for this flash" << std::endl;
          f.pe_v[opdet] = flash.PE(i) / _gain_correction[i] / 0.424;
          f.pe_err_v[opdet] = sqrt(flash.PE(i) / _gain_correction[i]) / 0.424;
        }
      }
      
      if (flash.Time() > -1.5 && flash.Time() < -0.5) {    // This is mucs flash
        _mucs_flash = 1;
        for (unsigned int pmtid = 0; pmtid < f.pe_v.size(); pmtid++) _mucs_flash_spec[pmtid] = f.pe_v[pmtid];  // Save the mucs flash pe spectrum
      }
      
      f.time = flash.Time();
      f.idx = flash_id;
      ++flash_id;
      _mgr.Emplace(std::move(f));
    }
  }
  
  if(!flash_id) {
    std::cerr << "No relevant flash found... " << std::endl;
    return;
  }
  if (flash_id>1) {
    std::cerr << "Found 2 candidates for MuCS flash. Aborting. " << std::endl;
    return;
  }
  
  ::art::Handle<std::vector<recob::Track> > track_h;
  e.getByLabel(_track_producer,track_h);
  if(!track_h.isValid() || track_h->empty())  {
    std::cerr << "No tracks found with" << _track_producer << std::endl;
    throw std::exception();
  }
  

  // ********************
  // Creating trajectory for the _TPC tracks_ and passing it to the Manager
  // ********************
  ::geoalgo::Trajectory mucs_geotrj;

  for (size_t trk_idx=0; trk_idx<track_h->size(); trk_idx++) {
    
    const art::Ptr<recob::Track> track_ptr(track_h,trk_idx);
    
    if(_mgr.QClusterArray().size() >= _num_tracks) {
      std::cout << "Reached maximum number of tracks (" << _num_tracks << "), breaking." << std::endl;
      break;
    }
    
    mucs_geotrj.resize(track_ptr->NumberTrajectoryPoints(),::geoalgo::Vector(0.,0.,0.));
    for (size_t pt_idx=0; pt_idx < track_ptr->NumberTrajectoryPoints(); ++pt_idx) {
      auto const& pt = track_ptr->LocationAtPoint(pt_idx);
      mucs_geotrj[pt_idx][0] = pt[0];
      mucs_geotrj[pt_idx][1] = pt[1];
      mucs_geotrj[pt_idx][2] = pt[2];
    }
    auto qcluster = ((flashana::LightPath*)(_mgr.GetCustomAlgo("LightPath")))->FlashHypothesis(mucs_geotrj);
    _mgr.Emplace(std::move(qcluster));
  }
  
  
  // ********************
  // Run the matching
  // ********************
  _result = _mgr.Match();
  
  
  // ********************
  // Save to tree
  // ********************
  _run    = e.id().run();
  _subrun = e.id().subRun();
  _event  = e.id().event();
  _matchid = 0;
  for(_matchid=0; _matchid < (int)(_result.size()); ++_matchid) {
    
    auto const& match = _result[_matchid];
    _tpcid    = match.tpc_id;
    _flashid  = match.flash_id;
    _score    = match.score;
    _qll_xmin = match.tpc_point.x;
    
    _tpc_xmin = 1.e4;
    for(auto const& pt : _mgr.QClusterArray()[_tpcid]) {
      if(pt.x < _tpc_xmin) _tpc_xmin = pt.x;
    }
    
    // The tagged track
    const art::Ptr<recob::Track> track_ptr(track_h,match.tpc_id);
    const TVector3 & start = track_ptr->Vertex();
    const TVector3 & end   = track_ptr->End();

    _mucs_track_startx = start.X();
    _mucs_track_starty = start.Y();
    _mucs_track_startz = start.Z();
    _mucs_track_endx = end.X();
    _mucs_track_endy = end.Y();
    _mucs_track_endz = end.Z();
    
    auto const& flash = _mgr.FlashArray()[_flashid];
    
    _t0 = flash.time;
    
    if(_hypo_spec.size() != match.hypothesis.size()) {
      std::cout << "Hypothesis size mismatch!" << std::endl;
      throw std::exception();
    }
    
    for(size_t pmt=0; pmt<_hypo_spec.size(); ++pmt) _hypo_spec[pmt]  = match.hypothesis[pmt];
    for(size_t pmt=0; pmt<_hypo_spec.size(); ++pmt) _flash_spec[pmt] = flash.pe_v[pmt];
    
    _flash_pe = 0.;
    _hypo_pe  = 0.;
    for(auto const& v : _hypo_spec) _hypo_pe += v;
    for(auto const& v : _flash_spec) _flash_pe += v;
    
    std::cout << "MuCSTrackFinder ends." << std::endl;
    _tree1->Fill();
    _tree2->Fill();
  }
  
  
  
  
  
  
}

DEFINE_ART_MODULE(MuCSTrackFinder)

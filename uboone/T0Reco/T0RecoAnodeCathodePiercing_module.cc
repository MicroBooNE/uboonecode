////////////////////////////////////////////////////////////////////////
// Class:       T0RecoAnodeCathodePiercing
// Module Type: producer
// File:        T0RecoAnodeCathodePiercing_module.cc
//
// David Caratelli - davidc1@fnal.gov - July 13 2016
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

// services etc...
#include "larcore/Geometry/Geometry.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

// data-products
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "lardata/Utilities/AssociationUtil.h"

// ROOT
#include "TVector3.h"

// C++
#include <memory>
#include <iostream>
#include <utility>

class T0RecoAnodeCathodePiercing;

class T0RecoAnodeCathodePiercing : public art::EDProducer {
public:
  explicit T0RecoAnodeCathodePiercing(fhicl::ParameterSet const & p);
  // The destructor generated by the compiler is fine for classes
  // without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  T0RecoAnodeCathodePiercing(T0RecoAnodeCathodePiercing const &) = delete;
  T0RecoAnodeCathodePiercing(T0RecoAnodeCathodePiercing &&) = delete;
  T0RecoAnodeCathodePiercing & operator = (T0RecoAnodeCathodePiercing const &) = delete;
  T0RecoAnodeCathodePiercing & operator = (T0RecoAnodeCathodePiercing &&) = delete;

  // Required functions.
  void produce(art::Event & e) override;


private:

  // producer of 3D reconstructed track to be used
  std::string fTrackProducer;

  // producer of reconstructed optical flashes
  std::string fFlashProducer;

  // set "resolution". How far away from the detector bounds
  // do we want to be to make a claim.
  double fResolution; // [cm]

  // drift velocity // cm / us
  double fDriftVelocity;

  // tag which types of tracks to reconstruct
  bool top2side;
  bool side2bottom;

  // debug (verbose) mode?
  bool _debug;

  // define top, bottom, front and back boundaries of TPC
  double _TOP, _BOTTOM, _FRONT, _BACK;
  
  // vector to hold flash-times for the event
  std::vector<double> _flash_times;
  std::vector<size_t> _flash_idx_v;

  // detector width [drift-coord]
  double _det_width; // [cm]

  // time resolution for flashes to match [us] (separate values for Anode and Cathode)
  double fTimeResA, fTimeResC;
  
  // minimum PE threshold for flash to make it into match
  double fPEmin;

  // time-adjustment (us) to align reconstructed T0 for Anode and Cathode-crossing tracks to reconstructed flash-times
  double fRecoT0TimeOffsetA, fRecoT0TimeOffsetC;

  // functions to be used throughout module
  bool   TrackEntersTop     (const std::vector<TVector3>& sorted_trk);
  bool   TrackEntersAnode   (const std::vector<TVector3>& sorted_trk);
  bool   TrackEntersSide    (const std::vector<TVector3>& sorted_trk);
  bool   TrackExitsBottom   (const std::vector<TVector3>& sorted_trk);
  bool   TrackExitsAnode    (const std::vector<TVector3>& sorted_trk);
  bool   TrackExitsSide     (const std::vector<TVector3>& sorted_trk);

  // functions to be used for organization in the module
  void   SortTrackPoints      (const recob::Track& track,
			       std::vector<TVector3>& sorted_trk);
  double GetEnteringTimeCoord (const std::vector<TVector3>& sorted_trk);
  double GetExitingTimeCoord  (const std::vector<TVector3>& sorted_trk);

  // validate flash matching by requiring PMT flash
  std::pair<double,size_t> FlashMatch(const double reco_time);

};


T0RecoAnodeCathodePiercing::T0RecoAnodeCathodePiercing(fhicl::ParameterSet const & p)
// :
// Initialize member data here.
{

  produces< std::vector< anab::T0 > >();
  produces< art::Assns <recob::Track, anab::T0> >();
  produces< art::Assns <recob::Track, recob::OpFlash> >();

  fTrackProducer     = p.get<std::string>("TrackProducer");
  fFlashProducer     = p.get<std::string>("FlashProducer");
  fResolution        = p.get<double>     ("Resolution");
  fTimeResA          = p.get<double>     ("TimeResA");
  fTimeResC          = p.get<double>     ("TimeResC");
  fRecoT0TimeOffsetA = p.get<double>     ("RecoT0TimeOffsetA");
  fRecoT0TimeOffsetC = p.get<double>     ("RecoT0TimeOffsetC");
  fPEmin             = p.get<double>     ("PEmin");
  top2side           = p.get<bool>       ("top2side");
  side2bottom        = p.get<bool>       ("side2bottom");
  _debug             = p.get<bool>       ("debug");
  
  // get boundaries based on detector bounds
  auto const* geom = lar::providerFrom<geo::Geometry>();

  _TOP    =   geom->DetHalfHeight() - fResolution;
  _BOTTOM = - geom->DetHalfHeight() + fResolution;
  _FRONT  =   fResolution;
  _BACK   =   geom->DetLength() - fResolution;
  
  _det_width = geom->DetHalfWidth() * 2;

  // Use '_detp' to find 'efield' and 'temp'
  auto const* _detp = lar::providerFrom<detinfo::DetectorPropertiesService>();
  double efield = _detp -> Efield();
  double temp   = _detp -> Temperature();
  // Determine the drift velocity from 'efield' and 'temp'
  fDriftVelocity = _detp -> DriftVelocity(efield,temp);

}

void T0RecoAnodeCathodePiercing::produce(art::Event & e)
{

  if (_debug)
  std::cout << "NEW EVENT" << std::endl;

  _flash_times.clear();
  _flash_idx_v.clear();

  // produce OpFlash data-product to be filled within module
  std::unique_ptr< std::vector<anab::T0> > T0_v(new std::vector<anab::T0>);
  std::unique_ptr< art::Assns <recob::Track, anab::T0> >       trk_t0_assn_v   ( new art::Assns<recob::Track, anab::T0>       );
  std::unique_ptr< art::Assns <recob::Track, recob::OpFlash> > trk_flash_assn_v( new art::Assns<recob::Track, recob::OpFlash> );

  // load Flash
  art::Handle<std::vector<recob::OpFlash> > flash_h;
  e.getByLabel(fFlashProducer,flash_h);

  // make sure flash look good
  if(!flash_h.isValid()) {
    std::cerr<<"\033[93m[ERROR]\033[00m ... could not locate Flash!"<<std::endl;
    throw std::exception();
  }

  // load tracks previously created for which T0 reconstruction should occur
  art::Handle<std::vector<recob::Track> > track_h;
  e.getByLabel(fTrackProducer,track_h);

  // make sure tracks look good
  if(!track_h.isValid()) {
    std::cerr<<"\033[93m[ERROR]\033[00m ... could not locate Track!"<<std::endl;
    throw std::exception();
  }

  std::vector<art::Ptr<recob::Track> > TrkVec;
  art::fill_ptr_vector(TrkVec, track_h);

  // prepare a vector of optical flash times, if flash above some PE cut value

  size_t flash_ctr = 0;
  for (auto const& flash : *flash_h){
    if (flash.TotalPE() > fPEmin){
      _flash_times.push_back( flash.Time() );
      _flash_idx_v.push_back(flash_ctr);
    }
    flash_ctr += 1;
  }// for all flashes

  if (_debug)
  std::cout << "Selected a total of " << _flash_times.size() << " OpFlashes" << std::endl;

  // loop through reconstructed tracks
  size_t trk_ctr = 0;

  for (auto& track : TrkVec){

    trk_ctr += 1;

    if (_debug)
    std::cout << "Looping through reco track " << trk_ctr << std::endl;

    // get sorted points for the track object [assuming downwards going]
    std::vector<TVector3> sorted_trk;
    SortTrackPoints(*track,sorted_trk);

    // Declare the variable 'trkT' up here so that I can continue and not fill the t0 object if trkT is still equal to 0
    double trkT = 0.;

    // keep track of whether it goes thorugh the anode or cathode
    bool anode = 0;
    
    // check if the track is of good quality:                                                                                                                              
    // the first type: must exit through the bottom [indicates track reconstruction has gone 'till the end] and enter through the anode or cathode                         
    if ( (TrackExitsBottom(sorted_trk) == true) and (side2bottom == true) ) { // This begins the loop for the tracks that exit the bottom

      if (_debug)
      std::cout << "\t track exits bottom" << std::endl;

      // If the program reaches this loop and the boolean is 'false'
      // then this track enters either through the top, front or back.
      // We aren't interested in those types of tracks (we can't reconstruct the T0) so we can just continue
      if (TrackEntersSide(sorted_trk) == false)
	continue;

      if (_debug)
      std::cout << "\t track enters side" << std::endl;

      // made it this far -> the track is good to be used                                                                                                                    
      // figure out if it enters the anode or cathode                                                                                                              
      bool enters_anode = TrackEntersAnode(sorted_trk);
      
      // get the X coordinate of the point piercing the anode/cathode (upon ENTERING) 
      double trkX = GetEnteringTimeCoord(sorted_trk);
      
      // reconstruct track T0 w.r.t. trigger time                                                                                                                              

      // The 'trkX' enters on the anode, the side of the TPC with a lower x value than the cathode
      if (enters_anode){
	trkT = trkX / fDriftVelocity + fRecoT0TimeOffsetA;
	anode = 1;
      }

      // This will also give a small T0 value, because the cathode is a distance of _det_width from the anode
      else{
	trkT = (trkX - _det_width) / fDriftVelocity + fRecoT0TimeOffsetC; 
	anode = 0;
      }
      
    } // This can end the case in which the track exits through the bottom 
    
    // the second type: the track must enter through the top and exit through either the anode or cathode

    // This begins the loop for the tracks that enter the top
    if ( (TrackEntersTop(sorted_trk) == true) and (top2side == true) ) { 

      if (_debug)
      std::cout << "\t track enters top" << std::endl;
      
      // If the program reaches this loop and the boolean is 'false',
      // then this track exits either through the bottom, front or back.
      // We aren't interested in those types of tracks (we can't reconstruct the T0) so we can just continue
      if (TrackExitsSide(sorted_trk) == false)
	continue;

      if (_debug)
      std::cout << "\t track exits top" << std::endl;
      
      // Now, use a function to find if the track exits through the anode or the cathode 
      bool exits_anode = TrackExitsAnode(sorted_trk);

      // Get the X coordinate of the point piercing the anode/cathode (upon EXITING)
      double trkX = GetExitingTimeCoord(sorted_trk);

      // reconstruct track T0 w.r.t. trigger time
      // The 'trkX' exits on the anode, the side of the TPC with a lower x value than the cathode
      if (exits_anode){
	trkT = trkX / fDriftVelocity + fRecoT0TimeOffsetA; 
	anode = 1;
      }
      // This will also give a small T0 value, because the cathode is a distance of _det_width from the anode
      else{
	trkT = (trkX - _det_width) / fDriftVelocity + fRecoT0TimeOffsetC; 
	anode = 0;
      }
      
    } // This can end the case in which the track enters through the top 
    
    if (_debug)
    std::cout << "\t this track has a reconstructed time = " << trkT << std::endl;

    // if the time does not match one from optical flashes -> don't reconstruct
    auto const& flash_match_result = FlashMatch(trkT);
    if ( (flash_match_result.first > fTimeResA) && (anode == 1) )
      continue;
    if ( (flash_match_result.first > fTimeResC) && (anode == 0) )
      continue;
    
    // DON'T CREATE the t0 object unless the reconstructed t0 is some value other than 0

    if (trkT != 0.0) {
    // create T0 object with this information!
    anab::T0 t0(trkT, 0, 0);
    
    T0_v->emplace_back(t0);
    util::CreateAssn(*this, e, *T0_v, track, *trk_t0_assn_v);

    // get pointer to individual track
    // TMP const art::Ptr<recob::Track>   trk_ptr(track_h,trk_ctr-1);
    const art::Ptr<recob::OpFlash> flash_ptr(flash_h, flash_match_result.second );
    if (_debug)
    std::cout << "\t mathed to flash w/ index " << flash_match_result.second << " w/ PE " << flash_ptr->TotalPE() << " and time " << flash_ptr->Time() << " vs reco time " << trkT << std::endl;
    trk_flash_assn_v->addSingle( track, flash_ptr );
    //util::CreateAssn(*this, e, flash_ptr, track, *trk_flash_assn_v);

    }

  }// for all reconstructed tracks
  
  e.put(std::move(T0_v));
  e.put(std::move(trk_t0_assn_v));
  if (_debug)
    std::cout << "create track flash association " << std::endl;
  e.put(std::move(trk_flash_assn_v));

}

std::pair<double,size_t> T0RecoAnodeCathodePiercing::FlashMatch(const double reco_time){
  
  // loop through all reco'd flash times and see if one matches
  // the reco time from the track
  double dt_min = 4000.; // us
  size_t idx_min = _flash_times.size();

  for (size_t i=0; i < _flash_times.size(); i++){
    auto const& time = _flash_times[i];
    double dt = fabs(time - reco_time);
    if (dt < dt_min){
      dt_min  = dt;
      idx_min = _flash_idx_v[i];
    }
  }

  std::pair<double,size_t> ret(dt_min,idx_min);
  return ret;
}


bool   T0RecoAnodeCathodePiercing::TrackEntersTop(const std::vector<TVector3>& sorted_trk)
{
  // check that the first point in the track
  // pierces the top boundary of the TPC
  // This track either will pierce the top of the TPC or is just about to (the '_TOP' variable is just below the actual coordinate position of the top in Y)

  if (sorted_trk.at(0).Y() > _TOP)
    return true;

  return false;
}

bool   T0RecoAnodeCathodePiercing::TrackEntersAnode(const std::vector<TVector3>& sorted_trk)
{

  // we know the track enters either the                                                                                                                               
  // anode or cathode                                                                                                                                                        
  // at this point figure out                                                                                                                                                
  // if it ENTERS the ANODE or CATHODE                                                                                         
  // ANODE: top point must be at lower X-coord                                                                                                
  // than bottom point                                                                                              
  // CATHODE: top point must be at larger X-coord                                                                                                                              
  // than bottom point
  // assume track has already been sorted                                                                                                                                     
  // such that the 1st point is the most elevated in Y coord.                                                                                                                 
  // return TRUE if passes the ANODE                                                                                                                                                
  
  auto const& top    = sorted_trk.at(0);
  auto const& bottom = sorted_trk.at( sorted_trk.size() - 1 );

  if (top.X() < bottom.X())
    return true;
  
  return false;
}


bool   T0RecoAnodeCathodePiercing::TrackEntersSide(const std::vector<TVector3>& sorted_trk)
{
  
  // check that the top-most point                                                                                                                                            
  // is not on the top of the TPC                                                                                                                                              
  // nor on the front & back of the TPC                                                                                                                                           
  
  auto const& top_pt = sorted_trk.at(0);

  // if highest point above the TOP -> false                                                                                                                                   
  if (top_pt.Y() > _TOP)
    return false;

  // if highest point in Z close to front or back                                                                                                                              
  // -> FALSE                                                                                                                                                                 
  if ( (top_pt.Z() < _FRONT) or (top_pt.Z() > _BACK) )
    return false;


  // If the function makes it this far, then it will enter through one of the sides of the TPC
  return true;
}

bool   T0RecoAnodeCathodePiercing::TrackExitsBottom(const std::vector<TVector3>& sorted_trk)
{

  // check that the last point in the track                                                                                                                                    
  // pierces the bottom boundary of the TPC                                                                                                                                   
  if ( sorted_trk.at( sorted_trk.size() - 1).Y() < _BOTTOM )
    return true;

  return false;
}

bool   T0RecoAnodeCathodePiercing::TrackExitsAnode(const std::vector<TVector3>& sorted_trk)
{

  // Check, once it's known that the track doesn't exit out of the bottom, whether it's the anode or                                                                           
  // the cathode that it exits out of                                                                                                                                         
  // This can be done by direct analogy with the 'Anode' function (shown in this file as the 'TrackEntersAnode') function written by D. Caratelli                             
  // Define 'top' as the point at the start of the track, and 'bottom' as the point at the end of the track                                                                     

  auto const& top    = sorted_trk.at(0);
  auto const& bottom = sorted_trk.at(sorted_trk.size() - 1);

  // Check to see which point has a lower x coordinate                                                                                                                         
  // If the bottom does, then it exits out of the anode                                                                                                                       
  // If the top does, then it exits out of the cathode                                                                                                                          
  if (bottom.X() < top.X()) 
    return true;

  return false; // Otherwise, the top is less than the bottom, so the track ended closer to the cathode and exited there                                                          
}


bool   T0RecoAnodeCathodePiercing::TrackExitsSide(const std::vector<TVector3>& sorted_trk)
{

  // check that the bottom-most point                                                                                                                                           
  // is not on the bottom of the TPC                                                                                                                                            
  // nor on the front & back of the TPC                                                                                                                                              

  auto const& bottom_pt = sorted_trk.at(sorted_trk.size() - 1);

  // if lowest point below the BOTTOM -> false                                                                                                            
  // Within this resolution, this means that it's likely that the track exited out of the bottom (at a point earlier on in the process than the last point) OR is just about to

  if (bottom_pt.Y() <  _BOTTOM)
    return false;

  // if lowest point in Z close to front or back                                                                                                         
  // -> FALSE                                                                                                                                              
  // If the the bottom point is less than the front, then the track has already pierced the front of the TPC and exited that way OR is likely just about to
  // If the bottom point is greater than the back, then the track has already pierced the back of the TPC and exited that way OR is likely just about to
  if ( (bottom_pt.Z() < _FRONT) or (bottom_pt.Z() > _BACK) )
    return false;

  return true;
}

void   T0RecoAnodeCathodePiercing::SortTrackPoints(const recob::Track& track, std::vector<TVector3>& sorted_trk)
{

  // vector to store 3D coordinates of                                                                                                                                           
  // ordered track                              
  sorted_trk.clear();

  // take the reconstructed 3D track                                                                                                                                           
  // and assuming it is downwards                                                                                                    
  // going, sort points so that                                                                                                              
  // the track starts at the top                                                                                                     
  // which point is further up in Y coord?                                                                                                                  
  // start or end?                                                                                                                 
  auto const&N = track.NumberTrajectoryPoints();
  auto const&start = track.LocationAtPoint(0);
  auto const&end   = track.LocationAtPoint( N - 1 );

  // if points are ordered correctly                                                                                                                                       
  if (start.Y() > end.Y()){
    for (size_t i=0; i < N; i++)
      sorted_trk.push_back( track.LocationAtPoint(i) );
  }
  
  // otherwise flip order                                                                                                                                                 
  else {
    for (size_t i=0; i < N; i++)
      sorted_trk.push_back( track.LocationAtPoint( N - i - 1) );
  }
}


double T0RecoAnodeCathodePiercing::GetEnteringTimeCoord(const std::vector<TVector3>& sorted_trk)
{

  // get the drift-coordinate value                                                                                                             
  // associated with the point                                                                                                                                      
  // along the track piercing the anode / cathode                                                                                                 
  // ** WHEN the track enters the anode / cathode
  return sorted_trk.at(0).X();
}


double T0RecoAnodeCathodePiercing::GetExitingTimeCoord(const std::vector<TVector3>& sorted_trk) 
{
  // get the drift-coordinate value                                                                                                                                           
  // associated with the point                                                                                                                                                
  // along the track piercing the anode / cathode                                                                                                                            
  // ** WHEN the track exits the anode / cathode
  return sorted_trk.at(sorted_trk.size() - 1).X();
}

DEFINE_ART_MODULE(T0RecoAnodeCathodePiercing)

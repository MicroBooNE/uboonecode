#ifndef SCANNERALGO_TEMPLATE_H
#define SCANNERALGO_TEMPLATE_H

#include "DataFormat/event_ass.h"
#include "DataFormat/sparse_vector.h"
#include "DataFormat/opdetwaveform.h"
#include "DataFormat/simphotons.h"
#include "DataFormat/trigger.h"
#include "DataFormat/potsummary.h"
#include "DataFormat/hit.h"
#include "DataFormat/track.h"
#include "DataFormat/mctruth.h"
#include "DataFormat/mctree.h"
#include "DataFormat/user_info.h"
#include "DataFormat/spacepoint.h"
#include "DataFormat/rawdigit.h"
#include "DataFormat/wire.h"
#include "DataFormat/hit.h"
#include "DataFormat/cluster.h"
#include "DataFormat/shower.h"
#include "DataFormat/mcshower.h"
#include "DataFormat/mctrack.h"
#include "DataFormat/simch.h"
#include "DataFormat/auxsimch.h"
#include "DataFormat/calorimetry.h"
#include "DataFormat/vertex.h"
#include "DataFormat/endpoint2d.h"
#include "DataFormat/seed.h"
#include "DataFormat/cosmictag.h"
#include "DataFormat/opflash.h"
#include "DataFormat/ophit.h"
#include "DataFormat/mcflux.h"
#include "DataFormat/pfpart.h"
#include "DataFormat/partid.h"
#include "DataFormat/gtruth.h"
#include "DataFormat/minos.h"
#include "DataFormat/pcaxis.h"
#include "DataFormat/flashmatch.h"
#include "DataFormat/mucsdata.h"
#include "DataFormat/mucsreco.h"
#include <TStopwatch.h>
/*
  This file defines certain specilization of templated functions.
  In particular it implements:

  ScannerAlgo::ScanData
  ScannerAlgo::GetPtrMap
  ScannerAlgo::LiteDataType

  One has to implement the relevant specilization for introducing a new data product!

 */

namespace larlite {

  //
  // ScanData specialization
  //
  template <>
  void ScannerAlgo::ScanData(art::Handle<std::vector< ::simb::MCTruth> > const &dh,
			     ::larlite::event_base* lite_dh)
  { 
    fDataReadFlag_v[lite_dh->data_type()][lite_dh->name()] = true;
    //auto name_index = NameIndex(lite_dh->data_type(),lite_dh->name());
    auto lite_data = (::larlite::event_mctruth*)lite_dh;
    for(size_t i=0; i<dh->size(); ++i) {
      
      const art::Ptr<::simb::MCTruth> mct_ptr(dh,i);
      
      larlite::mctruth lite_mct;
      
      // Generator (Origin) ID
      lite_mct.SetOrigin((larlite::simb::Origin_t) mct_ptr->Origin() );

      // Particle Information
      for(size_t j=0; j<(size_t)(mct_ptr->NParticles()); ++j) {
	
	const ::simb::MCParticle lar_part(mct_ptr->GetParticle(j));
	
	larlite::mcpart lite_part(lar_part.TrackId(),
				  lar_part.PdgCode(),
				  lar_part.Process(), 
				  lar_part.Mother(),
				  lar_part.Mass(),
				  lar_part.StatusCode());
	
	lite_part.SetPolarization(lar_part.Polarization());
	lite_part.SetRescatter(lar_part.Rescatter());
	lite_part.SetWeight(lar_part.Weight());
	
	for(size_t k=0; k<(size_t)(lar_part.NumberDaughters()); ++k)
	  
	  lite_part.AddDaughter(lar_part.Daughter(k));
	
	// Add trajectory points
	larlite::mctrajectory lite_track;
	lite_track.reserve(lar_part.NumberTrajectoryPoints());
	
	bool   loopInFV=false;
	size_t start_FV=0;
	for(size_t l=0; l<(size_t)lar_part.NumberTrajectoryPoints(); ++l) {
	  
	  lite_track.push_back(lar_part.Position(l),lar_part.Momentum(l));
	  
	  // Record fiducial volume tracking
	  bool inFV = false;//IsFV(lar_part.Vx(l),lar_part.Vy(l),lar_part.Vz(l),lar_part.T(l));
	  
	  if(!loopInFV) {
	    
	    if(inFV) { loopInFV=true; start_FV=l; }
	    
	  }else if(!inFV) {
	    lite_part.AddFiducialTrack(start_FV,l-1);
	    loopInFV=false;
	  }
	}
	if(loopInFV){ lite_part.AddFiducialTrack(start_FV,lar_part.NumberTrajectoryPoints()-1); }
	
	lite_part.SetTrajectory(lite_track);
	
	lite_mct.Add(lite_part);
      }
      
      // Neutrino Information
      const ::simb::MCNeutrino lar_nu(mct_ptr->GetNeutrino());

      if (lite_mct.GetParticles().size() )      
        lite_mct.SetNeutrino( lar_nu.CCNC(),
			      lar_nu.Mode(),
			      lar_nu.InteractionType(),
			      lar_nu.Target(),
			      lar_nu.HitNuc(),
			      lar_nu.HitQuark(),
			      lar_nu.W(),
			      lar_nu.X(),
			      lar_nu.Y(),
			      lar_nu.QSqr() );

      // Store address map for downstream association
      //fPtrIndex_mctruth[mct_ptr] = std::make_pair(lite_data->size(),name_index);

      // Save
      lite_data->push_back(lite_mct);
    }
  }
  
  template <>
  void ScannerAlgo::ScanData(art::Handle<std::vector< ::simb::GTruth> > const &dh,
			     ::larlite::event_base* lite_dh)
  { 
    fDataReadFlag_v[lite_dh->data_type()][lite_dh->name()] = true;  
    //auto name_index = NameIndex(lite_dh->data_type(),lite_dh->name());
    auto lite_data = (::larlite::event_gtruth*)lite_dh;
    for(size_t i=0; i<dh->size(); ++i) {
      
      const art::Ptr<::simb::GTruth> gtruth_ptr(dh,i);
      
      larlite::gtruth lite_gtruth;
      
      lite_gtruth.fGint        = gtruth_ptr->fGint;
      lite_gtruth.fGscatter    = gtruth_ptr->fGscatter;
      lite_gtruth.fweight      = gtruth_ptr->fweight;
      lite_gtruth.fprobability = gtruth_ptr->fprobability;
      lite_gtruth.fXsec        = gtruth_ptr->fXsec;
      lite_gtruth.fDiffXsec    = gtruth_ptr->fDiffXsec;
      lite_gtruth.fNumPiPlus   = gtruth_ptr->fNumPiPlus;
      lite_gtruth.fNumPiMinus  = gtruth_ptr->fNumPiMinus;
      lite_gtruth.fNumPi0      = gtruth_ptr->fNumPi0;
      lite_gtruth.fNumProton   = gtruth_ptr->fNumProton;
      lite_gtruth.fNumNeutron  = gtruth_ptr->fNumNeutron;
      lite_gtruth.fIsCharm     = gtruth_ptr->fIsCharm;
      lite_gtruth.fResNum      = gtruth_ptr->fResNum;
      lite_gtruth.fgQ2         = gtruth_ptr->fgQ2;
      lite_gtruth.fgq2         = gtruth_ptr->fgq2;
      lite_gtruth.fgW          = gtruth_ptr->fgW;
      lite_gtruth.fgT          = gtruth_ptr->fgT;
      lite_gtruth.fgX          = gtruth_ptr->fgX;
      lite_gtruth.fgY          = gtruth_ptr->fgY;
      lite_gtruth.fFShadSystP4 = gtruth_ptr->fFShadSystP4;
      lite_gtruth.fIsSeaQuark  = gtruth_ptr->fIsSeaQuark;
      lite_gtruth.fHitNucP4    = gtruth_ptr->fHitNucP4;
      lite_gtruth.ftgtZ        = gtruth_ptr->ftgtZ;
      lite_gtruth.ftgtA        = gtruth_ptr->ftgtA;
      lite_gtruth.ftgtPDG      = gtruth_ptr->ftgtPDG;

      //fPtrIndex_gtruth[gtruth_ptr] = std::make_pair(lite_data->size(),name_index);
      
      lite_data->push_back(lite_gtruth);
      
    }
    
  }

  template <>
  void ScannerAlgo::ScanData(art::Handle<std::vector< ::simb::MCParticle> > const &dh,
			     ::larlite::event_base* lite_dh)
  { 
    fDataReadFlag_v[lite_dh->data_type()][lite_dh->name()] = true;  
    //auto name_index = NameIndex(lite_dh->data_type(),lite_dh->name());
    auto lite_data = (::larlite::event_mcpart*)lite_dh;
    for(size_t i=0; i<dh->size(); ++i) {
      
      art::Ptr<::simb::MCParticle> mcparticle_ptr(dh,i);
      
      larlite::mcpart mcpart_lite(mcparticle_ptr->TrackId(),
				  mcparticle_ptr->PdgCode(),
				  mcparticle_ptr->Process(), 
				  mcparticle_ptr->Mother(),
				  mcparticle_ptr->Mass(),
				  mcparticle_ptr->StatusCode());
      
      mcpart_lite.SetPolarization(mcparticle_ptr->Polarization());
      mcpart_lite.SetRescatter(mcparticle_ptr->Rescatter());
      mcpart_lite.SetWeight(mcparticle_ptr->Weight());
      
      for(size_t k=0; k<(size_t)(mcparticle_ptr->NumberDaughters()); ++k)
	mcpart_lite.AddDaughter(mcparticle_ptr->Daughter(k));
      
      // Add trajectory points
      larlite::mctrajectory lite_track;
      lite_track.reserve(mcparticle_ptr->NumberTrajectoryPoints());
      
      bool   loopInFV=false;
      size_t start_FV=0;
      for(size_t l=0; l<(size_t)mcparticle_ptr->NumberTrajectoryPoints(); ++l) {
	
	lite_track.push_back(mcparticle_ptr->Position(l),
			     mcparticle_ptr->Momentum(l));
	
	// Record fiducial volume tracking
	bool inFV = false;//IsFV(mcparticle_ptr->Vx(l), mcparticle_ptr->Vy(l),mcparticle_ptr->Vz(l),mcparticle_ptr->T(l));
	
	if(!loopInFV) {
	  
	  if(inFV) { loopInFV=true; start_FV=l; }
	  
	}else if(!inFV) {
	  mcpart_lite.AddFiducialTrack(start_FV,l-1);
	  loopInFV=false;
	}
      }
      if(loopInFV){ mcpart_lite.AddFiducialTrack(start_FV,mcparticle_ptr->NumberTrajectoryPoints()-1); }
      
      mcpart_lite.SetTrajectory(lite_track);
      
      // Store address map for downstream association
      //fPtrIndex_mcpart[mcparticle_ptr] = std::make_pair(lite_data->size(),name_index);
      
      lite_data->push_back(mcpart_lite);
    }
  }
  
  template <>
  void ScannerAlgo::ScanData(art::Handle<std::vector< ::simb::MCFlux> > const &dh,
			     ::larlite::event_base* lite_dh)
  { 
    fDataReadFlag_v[lite_dh->data_type()][lite_dh->name()] = true;  
    //auto name_index = NameIndex(lite_dh->data_type(),lite_dh->name());
    auto lite_data = (::larlite::event_mcflux*)lite_dh;
    for(size_t i=0; i<dh->size(); ++i) {
      
      const art::Ptr<::simb::MCFlux> mcflux_ptr(dh,i);
      
      larlite::mcflux lite_mcflux;

      lite_mcflux.frun = mcflux_ptr->frun;
      lite_mcflux.fevtno = mcflux_ptr->fevtno;
      lite_mcflux.fndxdz = mcflux_ptr->fndxdz;
      lite_mcflux.fndydz = mcflux_ptr->fndydz;
      lite_mcflux.fnpz = mcflux_ptr->fnpz;
      lite_mcflux.fnenergy = mcflux_ptr->fnenergy;
      lite_mcflux.fndxdznea = mcflux_ptr->fndxdznea;
      lite_mcflux.fndydznea = mcflux_ptr->fndydznea;
      lite_mcflux.fnenergyn = mcflux_ptr->fnenergyn;
      lite_mcflux.fnwtnear = mcflux_ptr->fnwtnear;
      lite_mcflux.fndxdzfar = mcflux_ptr->fndxdzfar;
      lite_mcflux.fndydzfar = mcflux_ptr->fndydzfar;
      lite_mcflux.fnenergyf = mcflux_ptr->fnenergyf;
      lite_mcflux.fnwtfar = mcflux_ptr->fnwtfar;
      lite_mcflux.fnorig = mcflux_ptr->fnorig;
      lite_mcflux.fndecay = mcflux_ptr->fndecay;
      lite_mcflux.fntype = mcflux_ptr->fntype;
      lite_mcflux.fvx = mcflux_ptr->fvx;
      lite_mcflux.fvy = mcflux_ptr->fvy;
      lite_mcflux.fvz = mcflux_ptr->fvz;
      lite_mcflux.fpdpx = mcflux_ptr->fpdpx;
      lite_mcflux.fpdpy = mcflux_ptr->fpdpy;
      lite_mcflux.fpdpz = mcflux_ptr->fpdpz;
      lite_mcflux.fppdxdz = mcflux_ptr->fppdxdz;
      lite_mcflux.fppdydz = mcflux_ptr->fppdydz;
      lite_mcflux.fpppz = mcflux_ptr->fpppz;
      lite_mcflux.fppenergy = mcflux_ptr->fppenergy;
      lite_mcflux.fppmedium = mcflux_ptr->fppmedium;
      lite_mcflux.fptype = mcflux_ptr->fptype;
      lite_mcflux.fppvx = mcflux_ptr->fppvx;
      lite_mcflux.fppvy = mcflux_ptr->fppvy;
      lite_mcflux.fppvz = mcflux_ptr->fppvz;
      lite_mcflux.fmuparpx = mcflux_ptr->fmuparpx;
      lite_mcflux.fmuparpy = mcflux_ptr->fmuparpy;
      lite_mcflux.fmuparpz = mcflux_ptr->fmuparpz;
      lite_mcflux.fmupare = mcflux_ptr->fmupare;
      lite_mcflux.fnecm = mcflux_ptr->fnecm;
      lite_mcflux.fnimpwt = mcflux_ptr->fnimpwt;
      lite_mcflux.fxpoint = mcflux_ptr->fxpoint;
      lite_mcflux.fypoint = mcflux_ptr->fypoint;
      lite_mcflux.fzpoint = mcflux_ptr->fzpoint;
      lite_mcflux.ftvx = mcflux_ptr->ftvx;
      lite_mcflux.ftvy = mcflux_ptr->ftvy;
      lite_mcflux.ftvz = mcflux_ptr->ftvz;
      lite_mcflux.ftpx = mcflux_ptr->ftpx;
      lite_mcflux.ftpy = mcflux_ptr->ftpy;
      lite_mcflux.ftpz = mcflux_ptr->ftpz;
      lite_mcflux.ftptype = mcflux_ptr->ftptype;

      lite_mcflux.ftgen = mcflux_ptr->ftgen;
      lite_mcflux.ftgptype = mcflux_ptr->ftgptype;

      lite_mcflux.ftgppx = mcflux_ptr->ftgppx;
      lite_mcflux.ftgppy = mcflux_ptr->ftgppy;
      lite_mcflux.ftgppz = mcflux_ptr->ftgppz;
      lite_mcflux.ftprivx = mcflux_ptr->ftprivx;
      lite_mcflux.ftprivy = mcflux_ptr->ftprivy;
      lite_mcflux.ftprivz = mcflux_ptr->ftprivz;
      lite_mcflux.fbeamx = mcflux_ptr->fbeamx;
      lite_mcflux.fbeamy = mcflux_ptr->fbeamy;
      lite_mcflux.fbeamz = mcflux_ptr->fbeamz;
      lite_mcflux.fbeampx = mcflux_ptr->fbeampx;
      lite_mcflux.fbeampy = mcflux_ptr->fbeampy;
      lite_mcflux.fbeampz = mcflux_ptr->fbeampz;

      lite_mcflux.fFluxType = (::larlite::simb::flux_code_)(mcflux_ptr->fFluxType);

      lite_mcflux.fgenx = mcflux_ptr->fgenx;

      lite_mcflux.fgeny = mcflux_ptr->fgeny;
      lite_mcflux.fgenz = mcflux_ptr->fgenz;
      lite_mcflux.fdk2gen = mcflux_ptr->fdk2gen;

      lite_mcflux.fgen2vtx = mcflux_ptr->fgen2vtx;

      lite_mcflux.SetFluxGen( mcflux_ptr->Flux( 12,::simb::kGenerator), 
			      mcflux_ptr->Flux(-12,::simb::kGenerator),
			      mcflux_ptr->Flux( 14,::simb::kGenerator),
			      mcflux_ptr->Flux(-14,::simb::kGenerator),
			      mcflux_ptr->Flux( 16,::simb::kGenerator),
			      mcflux_ptr->Flux(-16,::simb::kGenerator) );

      lite_mcflux.SetFluxPos( mcflux_ptr->Flux( 12,::simb::kHistPlusFocus), 
			      mcflux_ptr->Flux(-12,::simb::kHistPlusFocus),
			      mcflux_ptr->Flux( 14,::simb::kHistPlusFocus),
			      mcflux_ptr->Flux(-14,::simb::kHistPlusFocus),
			      mcflux_ptr->Flux( 16,::simb::kHistPlusFocus),
			      mcflux_ptr->Flux(-16,::simb::kHistPlusFocus) );

      lite_mcflux.SetFluxPos( mcflux_ptr->Flux( 12,::simb::kHistMinusFocus), 
			      mcflux_ptr->Flux(-12,::simb::kHistMinusFocus),
			      mcflux_ptr->Flux( 14,::simb::kHistMinusFocus),
			      mcflux_ptr->Flux(-14,::simb::kHistMinusFocus),
			      mcflux_ptr->Flux( 16,::simb::kHistMinusFocus),
			      mcflux_ptr->Flux(-16,::simb::kHistMinusFocus) );

      // Register pointer to association look-up map
      //fPtrIndex_mcflux[mcflux_ptr] = std::make_pair(lite_data->size(),name_index);

      lite_data->push_back(lite_mcflux);
    }
  }
  
  template <>
  void ScannerAlgo::ScanData(art::Handle<std::vector< ::sim::MCShower> > const &dh,
			     ::larlite::event_base* lite_dh)
  { 
    fDataReadFlag_v[lite_dh->data_type()][lite_dh->name()] = true;  
    //auto name_index = NameIndex(lite_dh->data_type(),lite_dh->name());
    auto lite_data = (::larlite::event_mcshower*)lite_dh;
    for(size_t i=0; i<dh->size(); ++i) {

      art::Ptr<sim::MCShower> mcs_ptr(dh,i);
      
      larlite::mcshower lite_mcs;

      lite_mcs.Origin  ( (::larlite::simb::Origin_t)(mcs_ptr->Origin())  );

      lite_mcs.PdgCode ( mcs_ptr->PdgCode() );
      lite_mcs.TrackID ( mcs_ptr->TrackID() );
      lite_mcs.Process ( mcs_ptr->Process() );
      lite_mcs.Start   ( ::larlite::mcstep( mcs_ptr->Start().Position(), mcs_ptr->Start().Momentum() )  );
      lite_mcs.End     ( ::larlite::mcstep( mcs_ptr->End().Position(),   mcs_ptr->End().Momentum()   )  );

      lite_mcs.MotherPdgCode ( mcs_ptr->MotherPdgCode() );
      lite_mcs.MotherTrackID ( mcs_ptr->MotherTrackID() );
      lite_mcs.MotherProcess ( mcs_ptr->MotherProcess() );
      lite_mcs.MotherStart   ( ::larlite::mcstep( mcs_ptr->MotherStart().Position(), mcs_ptr->MotherStart().Momentum() ) );
      lite_mcs.MotherEnd     ( ::larlite::mcstep( mcs_ptr->MotherEnd().Position(),   mcs_ptr->MotherEnd().Momentum()   ) );

      lite_mcs.AncestorPdgCode ( mcs_ptr->AncestorPdgCode() );
      lite_mcs.AncestorTrackID ( mcs_ptr->AncestorTrackID() );
      lite_mcs.AncestorProcess ( mcs_ptr->AncestorProcess() );
      lite_mcs.AncestorStart   ( ::larlite::mcstep( mcs_ptr->AncestorStart().Position(), mcs_ptr->AncestorStart().Momentum() ) );
      lite_mcs.AncestorEnd     ( ::larlite::mcstep( mcs_ptr->AncestorEnd().Position(),   mcs_ptr->AncestorEnd().Momentum()   ) );

      lite_mcs.Charge( mcs_ptr->Charge() );

      lite_mcs.DetProfile( ::larlite::mcstep( mcs_ptr->DetProfile().Position(), mcs_ptr->DetProfile().Momentum() ) );

      lite_mcs.DaughterTrackID( mcs_ptr->DaughterTrackID() );

      lite_mcs.StartDir( mcs_ptr->StartDir() );
      lite_mcs.dEdx( mcs_ptr->dEdx() );
      lite_mcs.dQdx( mcs_ptr->dQdx() );


      //fPtrIndex_mcshower[mcs_ptr] = std::make_pair(lite_data->size(),name_index);

      lite_data->push_back(lite_mcs);
    }
  }


  template <>
  void ScannerAlgo::ScanData(art::Handle<std::vector< ::sim::MCTrack> > const &dh,
			     ::larlite::event_base* lite_dh)
  { 
    fDataReadFlag_v[lite_dh->data_type()][lite_dh->name()] = true;  
    //auto name_index = NameIndex(lite_dh->data_type(),lite_dh->name());
    auto lite_data = (::larlite::event_mctrack*)lite_dh;
    for(size_t i=0; i<dh->size(); ++i) {

      art::Ptr<sim::MCTrack> mct_ptr(dh,i);
      
      larlite::mctrack lite_mct;

      lite_mct.Origin  ( (::larlite::simb::Origin_t)(mct_ptr->Origin())  );
     
      lite_mct.PdgCode ( mct_ptr->PdgCode() );
      lite_mct.TrackID ( mct_ptr->TrackID() );
      lite_mct.Process ( mct_ptr->Process() );
      lite_mct.Start   ( ::larlite::mcstep( mct_ptr->Start().Position(), mct_ptr->Start().Momentum() )  );
      lite_mct.End     ( ::larlite::mcstep( mct_ptr->End().Position(),   mct_ptr->End().Momentum()   )  );
      lite_mct.dEdx( mct_ptr->dEdx() );
      lite_mct.dQdx( mct_ptr->dQdx() );

      lite_mct.MotherPdgCode ( mct_ptr->MotherPdgCode() );
      lite_mct.MotherTrackID ( mct_ptr->MotherTrackID() );
      lite_mct.MotherProcess ( mct_ptr->MotherProcess() );
      lite_mct.MotherStart   ( ::larlite::mcstep( mct_ptr->MotherStart().Position(), mct_ptr->MotherStart().Momentum() ) );
      lite_mct.MotherEnd     ( ::larlite::mcstep( mct_ptr->MotherEnd().Position(),   mct_ptr->MotherEnd().Momentum()   ) );

      lite_mct.AncestorPdgCode ( mct_ptr->AncestorPdgCode() );
      lite_mct.AncestorTrackID ( mct_ptr->AncestorTrackID() );
      lite_mct.AncestorProcess ( mct_ptr->AncestorProcess() );
      lite_mct.AncestorStart   ( ::larlite::mcstep( mct_ptr->AncestorStart().Position(), mct_ptr->AncestorStart().Momentum() ) );
      lite_mct.AncestorEnd     ( ::larlite::mcstep( mct_ptr->AncestorEnd().Position(),   mct_ptr->AncestorEnd().Momentum()   ) );

      const ::sim::MCTrack& mct(*mct_ptr);
 
      lite_mct.reserve(mct.size());

      for(auto const& s : mct)
      
	lite_mct.push_back( ::larlite::mcstep(s.Position(), s.Momentum()) );

      //fPtrIndex_mctrack[mct_ptr] = std::make_pair(lite_data->size(),name_index);
      
      lite_data->push_back(lite_mct);
    }
  }

  template <>
  void ScannerAlgo::ScanData(art::Handle<std::vector< ::sim::SimChannel> > const &dh,
			     ::larlite::event_base* lite_dh)
  { 
    fDataReadFlag_v[lite_dh->data_type()][lite_dh->name()] = true;  
    //auto name_index = NameIndex(lite_dh->data_type(),lite_dh->name());
    auto lite_data = (::larlite::event_simch*)lite_dh;

    for(size_t i=0; i<dh->size(); ++i ) {
      
      const art::Ptr<::sim::SimChannel> sch_ptr(dh,i);
      
      std::map<unsigned short,std::vector<sim::IDE> > sch_map(sch_ptr->TDCIDEMap());
      
      larlite::simch lite_sch;
      lite_sch.set_channel(sch_ptr->Channel());
      
      for(auto sch_iter = sch_map.begin(); sch_iter!=sch_map.end(); ++sch_iter) {
	
	unsigned short tdc = (*sch_iter).first;
	
	for(auto const this_ide : (*sch_iter).second) {
	  
	  larlite::ide lite_ide;
	  lite_ide.trackID = this_ide.trackID;
	  lite_ide.numElectrons = this_ide.numElectrons;
	  lite_ide.energy = this_ide.energy;
	  lite_ide.x = this_ide.x;
	  lite_ide.y = this_ide.y;
	  lite_ide.z = this_ide.z;
	  
	  lite_sch.add_ide(tdc,lite_ide);
	}
      }

      //fPtrIndex_simch[sch_ptr] = std::make_pair(lite_data->size(),name_index);

      lite_data->push_back(lite_sch);
    }

  }

  template <>
  void ScannerAlgo::ScanData(art::Handle<std::vector< ::sim::AuxDetSimChannel> > const &dh,
			     ::larlite::event_base* lite_dh)
  { 
    fDataReadFlag_v[lite_dh->data_type()][lite_dh->name()] = true;  
    //auto name_index = NameIndex(lite_dh->data_type(),lite_dh->name());
    auto lite_data = (::larlite::event_auxsimch*)lite_dh;

    for(size_t i=0; i<dh->size(); ++i ) {
      
      const art::Ptr<::sim::AuxDetSimChannel> auxsch_ptr(dh,i);

      std::vector<larlite::auxide> lite_ide_v;

      for(auto const& ide : auxsch_ptr->AuxDetIDEs()) {

	::larlite::auxide lite_ide;
	lite_ide.trackID         = ide.trackID;
	lite_ide.energyDeposited = ide.energyDeposited;
	lite_ide.entryX          = ide.entryX;
	lite_ide.entryY          = ide.entryY;
	lite_ide.entryZ          = ide.entryZ;
	lite_ide.entryT          = ide.entryT;
	lite_ide.exitX           = ide.exitX;
	lite_ide.exitY           = ide.exitY;
	lite_ide.exitZ           = ide.exitZ;
	lite_ide.exitT           = ide.exitT;
	lite_ide.exitMomentumX   = ide.exitMomentumX;
	lite_ide.exitMomentumY   = ide.exitMomentumY;
	lite_ide.exitMomentumZ   = ide.exitMomentumZ;

	lite_ide_v.emplace_back(lite_ide);
      }

      ::larlite::auxsimch lite_auxsch(auxsch_ptr->AuxDetID(),
				      std::move(lite_ide_v),
				      auxsch_ptr->AuxDetSensitiveID());
      lite_data->emplace_back(lite_auxsch);
    }

  }

  template <>
  void ScannerAlgo::ScanData(art::Handle<std::vector< ::MuCS::MuCSData> > const &dh,
			     ::larlite::event_base* lite_dh)
  { 
    fDataReadFlag_v[lite_dh->data_type()][lite_dh->name()] = true;  
    auto lite_data = (::larlite::event_mucsdata*)lite_dh;

    for(size_t i=0; i<dh->size(); ++i ) {
      
      const art::Ptr<::MuCS::MuCSData> mucs_ptr(dh,i);
      
      auto adc1 = mucs_ptr->ADC1();
      auto adc2 = mucs_ptr->ADC2();
      auto adc3 = mucs_ptr->ADC3();
      auto adc7 = mucs_ptr->ADC7();
      auto hits1 = mucs_ptr->Hits1();
      auto hits2 = mucs_ptr->Hits2();
      auto hits3 = mucs_ptr->Hits3();
      auto hits7 = mucs_ptr->Hits7();
      
      larlite::mucsdata lite_mucs( mucs_ptr->T0(), 
				   &(adc1[0]), &(adc2[0]), &(adc3[0]), &(adc7[0]),
				   hits1, hits2, hits3, hits7);

      lite_data->push_back(lite_mucs);
    }

  }

  template <>
  void ScannerAlgo::ScanData(art::Handle<std::vector< ::MuCS::MuCSRecoData> > const &dh,
			     ::larlite::event_base* lite_dh)
  { 
    fDataReadFlag_v[lite_dh->data_type()][lite_dh->name()] = true;  
    auto lite_data = (::larlite::event_mucsreco*)lite_dh;

    for(size_t i=0; i<dh->size(); ++i ) {
      
      const art::Ptr<::MuCS::MuCSRecoData> mucs_ptr(dh,i);

      larlite::mucsreco lite_mucs(mucs_ptr->theta_xy(),
				  mucs_ptr->theta_xy_rms(),
				  mucs_ptr->x(),
				  mucs_ptr->x_rms(),
				  mucs_ptr->theta_yz(),
				  mucs_ptr->theta_yz_rms(),
				  mucs_ptr->z(),
				  mucs_ptr->z_rms(),
				  mucs_ptr->y(),
				  mucs_ptr->xmatches(),
				  mucs_ptr->zmatches());

      lite_data->push_back(lite_mucs);
    }

  }

  template <>
  void ScannerAlgo::ScanData(art::Handle<std::vector< ::raw::RawDigit> > const &dh,
			     ::larlite::event_base* lite_dh)
  { 
    fDataReadFlag_v[lite_dh->data_type()][lite_dh->name()] = true;  
    //auto name_index = NameIndex(lite_dh->data_type(),lite_dh->name());
    auto lite_data = (::larlite::event_rawdigit*)lite_dh;

    for(size_t i=0; i<dh->size(); i++){

      const art::Ptr<::raw::RawDigit> rawdigit_ptr(dh,i);
      std::vector<short> adclist(rawdigit_ptr->NADC(),0);
      for(size_t tick=0;tick<adclist.size();tick++)
	adclist[tick]=rawdigit_ptr->ADC(tick);
    
      larlite::rawdigit rawdigit_lite( rawdigit_ptr->Channel(),
				       rawdigit_ptr->NADC(),
				       adclist,
				       (::larlite::raw::Compress_t)(rawdigit_ptr->Compression()));
      //fPtrIndex_rawdigit[rawdigit_ptr] = std::make_pair(lite_data->size(),name_index);
      
      lite_data->push_back(rawdigit_lite);
    }  
  }

  template <>
  void ScannerAlgo::ScanData(art::Handle<std::vector< ::raw::OpDetWaveform> > const &dh,
			     ::larlite::event_base* lite_dh)
  { 
    fDataReadFlag_v[lite_dh->data_type()][lite_dh->name()] = true;  
    //auto name_index = NameIndex(lite_dh->data_type(),lite_dh->name());
    auto lite_data = (::larlite::event_opdetwaveform*)lite_dh;

    for(size_t i=0; i<dh->size(); i++){

      const art::Ptr<::raw::OpDetWaveform> opdigit_ptr(dh,i);
      ::larlite::opdetwaveform lite_opdigit;
      lite_opdigit.SetChannelNumber(opdigit_ptr->ChannelNumber());
      lite_opdigit.SetTimeStamp(opdigit_ptr->TimeStamp());

      lite_opdigit.reserve((*opdigit_ptr).size());
      for(auto const& adc : *opdigit_ptr)
	lite_opdigit.push_back(adc);
      
      lite_data->emplace_back(lite_opdigit);
    }  
  }


  template <>
    void ScannerAlgo::ScanData(art::Handle< std::vector<::raw::Trigger> > const &dh,
			       ::larlite::event_base* lite_dh)
    { 
      
      //fDataReadFlag_v[lite_dh->data_type()][lite_dh->name()] = true;  
      //auto name_index = NameIndex(lite_dh->data_type(),lite_dh->name());
      auto lite_data = (::larlite::trigger*)lite_dh;
      
      if (dh->size() == 0)
	return;
      const art::Ptr<::raw::Trigger> trigger_ptr(dh,0);
      lite_data->TriggerNumber(trigger_ptr->TriggerNumber());
      lite_data->TriggerTime(trigger_ptr->TriggerTime());
      lite_data->BeamGateTime(trigger_ptr->BeamGateTime());
      lite_data->TriggerBits(trigger_ptr->TriggerBits());
      
    }
  
  template <>
  void ScannerAlgo::ScanData(art::Handle<std::vector< ::recob::Wire> > const &dh,
			     ::larlite::event_base* lite_dh)
  { 
    fDataReadFlag_v[lite_dh->data_type()][lite_dh->name()] = true;  
    //auto name_index = NameIndex(lite_dh->data_type(),lite_dh->name());
    auto lite_data = (::larlite::event_wire*)lite_dh;

    for(size_t i=0; i<dh->size(); i++){

      const art::Ptr<::recob::Wire> wire_ptr(dh,i);

      ::larlite::sparse_vector<float> rois;

      auto const& signalROI = wire_ptr->SignalROI();

      for(const auto& range : signalROI.get_ranges())
	rois.add_range(range.begin_index(),range.data());
      
      larlite::wire wire_lite(rois,
			      wire_ptr->Channel(),
			      (::larlite::geo::View_t)(wire_ptr->View()));
      
      lite_data->push_back(wire_lite);
    }  
  }

  template <>
  void ScannerAlgo::ScanData(art::Handle<std::vector< ::recob::Hit> > const &dh,
			     ::larlite::event_base* lite_dh)
  { 
    fDataReadFlag_v[lite_dh->data_type()][lite_dh->name()] = true;
    //auto name_index = NameIndex(lite_dh->data_type(),lite_dh->name());    
    auto lite_data = (::larlite::event_hit*)lite_dh;
    art::ServiceHandle<::geo::Geometry> geo;

    for(size_t i=0; i<dh->size(); i++){
      
      art::Ptr<::recob::Hit> hit_ptr(dh,i);

      //std::cout<<i<<" "<<hit_ptr.id().productIndex()<<std::endl;

      larlite::hit lite_hit;
      lite_hit.set_rms(hit_ptr->RMS());
      lite_hit.set_time_range(hit_ptr->StartTick(),hit_ptr->EndTick());
      lite_hit.set_time_peak(hit_ptr->PeakTime(),hit_ptr->SigmaPeakTime());
      lite_hit.set_amplitude(hit_ptr->PeakAmplitude(),hit_ptr->SigmaPeakAmplitude());
      lite_hit.set_sumq(hit_ptr->SummedADC());
      lite_hit.set_integral(hit_ptr->Integral(),hit_ptr->SigmaIntegral());
      lite_hit.set_multiplicity(hit_ptr->Multiplicity());
      lite_hit.set_local_index(hit_ptr->LocalIndex());
      lite_hit.set_goodness(hit_ptr->GoodnessOfFit());
      lite_hit.set_ndf(hit_ptr->DegreesOfFreedom());
      lite_hit.set_channel(hit_ptr->Channel());
      lite_hit.set_view((::larlite::geo::View_t)(hit_ptr->View()));
      lite_hit.set_signal_type((::larlite::geo::SigType_t)(hit_ptr->SignalType()));

      auto const& wid = hit_ptr->WireID();
      lite_hit.set_wire(::larlite::geo::WireID(wid.Cryostat,wid.TPC,wid.Plane,wid.Wire));
      
      // Store address map for downstream association
      //fPtrIndex_hit[hit_ptr] = std::make_pair(lite_data->size(),name_index);
      
      lite_data->push_back(lite_hit);
    }
  }

  template <>
  void ScannerAlgo::ScanData(art::Handle<std::vector< ::recob::OpHit> > const &dh,
			     ::larlite::event_base* lite_dh)
  { 
    fDataReadFlag_v[lite_dh->data_type()][lite_dh->name()] = true;  

    auto lite_data = (::larlite::event_ophit*)lite_dh;

    for(size_t i=0; i<dh->size(); ++i) {
      
      art::Ptr<::recob::OpHit> hit_ptr(dh,i);
      
      ::larlite::ophit lite_hit( hit_ptr->OpChannel(),
				 hit_ptr->PeakTime(),
				 hit_ptr->PeakTimeAbs(),
				 hit_ptr->Frame(),
				 hit_ptr->Width(),
				 hit_ptr->Area(),
				 hit_ptr->Amplitude(),
				 hit_ptr->PE(),
				 hit_ptr->FastToTotal() );

      //fPtrIndex_ophit[hit_ptr] = std::make_pair(lite_data->size(),name_index);
      
      lite_data->push_back(lite_hit);
    }
    
  }

  template <>
  void ScannerAlgo::ScanData(art::Handle<std::vector< ::recob::OpFlash> > const &dh,
			     ::larlite::event_base* lite_dh)
  { 
    fDataReadFlag_v[lite_dh->data_type()][lite_dh->name()] = true;  
    //auto name_index = NameIndex(lite_dh->data_type(),lite_dh->name());
    auto lite_data = (::larlite::event_opflash*)lite_dh;
    art::ServiceHandle<::geo::Geometry> geo;

    for(size_t i=0; i<dh->size(); ++i) {

      art::Ptr<::recob::OpFlash> flash_ptr(dh,i);
      
      std::vector<double> pe_per_opdet;
      pe_per_opdet.reserve(geo->NOpChannels());
      for(size_t j=0; j<geo->NOpChannels(); ++j)
	pe_per_opdet.push_back(flash_ptr->PE(j));
      
      ::larlite::opflash lite_flash( flash_ptr->Time(),
				     flash_ptr->TimeWidth(),
				     flash_ptr->AbsTime(),
				     flash_ptr->Frame(),
				     pe_per_opdet,
				     flash_ptr->InBeamFrame(),
				     flash_ptr->OnBeamTime(),
				     flash_ptr->FastToTotal(),
				     flash_ptr->YCenter(),
				     flash_ptr->YWidth(),
				     flash_ptr->ZCenter(),
				     flash_ptr->ZWidth(),
				     flash_ptr->WireCenters(),
				     flash_ptr->WireWidths());
      
      //fPtrIndex_opflash[flash_ptr] = std::make_pair(lite_data->size(),name_index);

      lite_data->push_back(lite_flash);
    }
  }

  template <>
  void ScannerAlgo::ScanData(art::Handle<std::vector< ::anab::CosmicTag> > const &dh,
			     ::larlite::event_base* lite_dh)
  { 
    fDataReadFlag_v[lite_dh->data_type()][lite_dh->name()] = true;  
    //auto name_index = NameIndex(lite_dh->data_type(),lite_dh->name());
    auto lite_data = (::larlite::event_cosmictag*)lite_dh;
    for(size_t i=0; i<dh->size(); ++i) {
      
      art::Ptr<::anab::CosmicTag> tag_ptr(dh,i);
      
      larlite::cosmictag lite_tag(tag_ptr->EndPoint1(),
				  tag_ptr->EndPoint2(),
				  tag_ptr->CosmicScore(),
				  (::larlite::anab::CosmicTagID_t)((int)(tag_ptr->CosmicType())));
      
      // store product ptr for association
      //fPtrIndex_cosmictag[tag_ptr] = std::make_pair(lite_data->size(),name_index);
      
      lite_data->push_back(lite_tag);
      
    }
  }

  template <>
  void ScannerAlgo::ScanData(art::Handle<std::vector< ::recob::Cluster> > const &dh,
			     ::larlite::event_base* lite_dh)
  { 
    fDataReadFlag_v[lite_dh->data_type()][lite_dh->name()] = true;  
    //auto name_index = NameIndex(lite_dh->data_type(),lite_dh->name());
    auto lite_data = (::larlite::event_cluster*)lite_dh;
    for(size_t i=0; i<dh->size(); ++i) {
      
      art::Ptr<::recob::Cluster> cluster_ptr(dh,i);
      
      larlite::cluster lite_cluster;
      lite_cluster.set_n_hits(cluster_ptr->NHits());

      lite_cluster.set_start_wire(cluster_ptr->StartWire(), cluster_ptr->SigmaStartWire());
      lite_cluster.set_start_tick(cluster_ptr->StartTick(), cluster_ptr->SigmaStartTick());
      lite_cluster.set_start_angle(cluster_ptr->StartAngle());
      lite_cluster.set_start_charge(cluster_ptr->StartCharge());
      lite_cluster.set_start_opening(cluster_ptr->StartOpeningAngle());

      lite_cluster.set_end_wire(cluster_ptr->EndWire(), cluster_ptr->SigmaEndWire());
      lite_cluster.set_end_tick(cluster_ptr->EndTick(), cluster_ptr->SigmaEndTick());
      lite_cluster.set_end_angle(cluster_ptr->EndAngle());
      lite_cluster.set_end_charge(cluster_ptr->EndCharge());
      lite_cluster.set_end_opening(cluster_ptr->EndOpeningAngle());

      lite_cluster.set_integral(cluster_ptr->Integral(), 
				cluster_ptr->IntegralStdDev(), 
				cluster_ptr->IntegralAverage());
      
      lite_cluster.set_summedADC(cluster_ptr->SummedADC(),
				 cluster_ptr->SummedADCstdDev(),
				 cluster_ptr->SummedADCaverage());

      lite_cluster.set_multiple_hit_density(cluster_ptr->MultipleHitDensity());
      lite_cluster.set_width(cluster_ptr->Width());

      lite_cluster.set_id(cluster_ptr->ID());
      lite_cluster.set_view((::larlite::geo::View_t)cluster_ptr->View());
      
      auto const& pid = cluster_ptr->Plane();
      lite_cluster.set_planeID( ::larlite::geo::PlaneID( pid.Cryostat, pid.TPC, pid.Plane ) );
      
      // Store address map for downstream association
      //fPtrIndex_cluster[cluster_ptr] = std::make_pair(lite_data->size(),name_index);
      
      lite_data->push_back(lite_cluster);
    }
  }

  template <>
  void ScannerAlgo::ScanData(art::Handle<std::vector< ::recob::Seed> > const &dh,
			     ::larlite::event_base* lite_dh)
  { 
    fDataReadFlag_v[lite_dh->data_type()][lite_dh->name()] = true;  
    //auto name_index = NameIndex(lite_dh->data_type(),lite_dh->name());
    auto lite_data = (::larlite::event_seed*)lite_dh;

    double fSeedPoint[3];
    double fSeedDirection[3];
    double fSeedPointError[3];
    double fSeedDirectionError[3];

    for(size_t i=0; i<dh->size(); ++i) {
      
      art::Ptr<::recob::Seed> seed_ptr(dh,i);
      larlite::seed lite_seed;

      seed_ptr->GetPoint(fSeedPoint,fSeedPointError);
      seed_ptr->GetDirection(fSeedDirection,fSeedDirectionError);

      lite_seed.SetPoint(fSeedPoint,fSeedPointError);
      lite_seed.SetDirection(fSeedDirection,fSeedDirectionError);
      lite_seed.SetValidity(seed_ptr->IsValid());
      
      // Store address map for downstream association
      //fPtrIndex_seed[seed_ptr] = std::make_pair(lite_data->size(),name_index);
      
      lite_data->push_back(lite_seed);
    }
  }

  template <>
  void ScannerAlgo::ScanData(art::Handle<std::vector< ::recob::EndPoint2D> > const &dh,
			     ::larlite::event_base* lite_dh)
  { 
    fDataReadFlag_v[lite_dh->data_type()][lite_dh->name()] = true;  
    //auto name_index = NameIndex(lite_dh->data_type(),lite_dh->name());
    auto lite_data = (::larlite::event_endpoint2d*)lite_dh;
    for(size_t i=0; i < dh->size(); ++i) {
      
      const art::Ptr<::recob::EndPoint2D> end2d_ptr(dh,i);
      
      // get vertex info
      larlite::endpoint2d lite_end2d(end2d_ptr->DriftTime(),
				     end2d_ptr->WireID().Wire,
				     end2d_ptr->Strength(),
				     end2d_ptr->ID(),
				     (larlite::geo::View_t)(end2d_ptr->View()),
				     end2d_ptr->Charge());
      
      //fPtrIndex_end2d[end2d_ptr] = std::make_pair(lite_data->size(),name_index);
      
      // Store data
      lite_data->push_back(lite_end2d);
    }
  }

  template <>
  void ScannerAlgo::ScanData(art::Handle<std::vector< ::recob::SpacePoint> > const &dh,
			     ::larlite::event_base* lite_dh)
  { 
    fDataReadFlag_v[lite_dh->data_type()][lite_dh->name()] = true;  
    //auto name_index = NameIndex(lite_dh->data_type(),lite_dh->name());
    auto lite_data = (::larlite::event_spacepoint*)lite_dh;
    for(size_t i=0; i<dh->size(); ++i) {
      
      art::Ptr<::recob::SpacePoint> spacepoint_ptr(dh,i);
      
      larlite::spacepoint lite_spacepoint(spacepoint_ptr->ID(),
					  spacepoint_ptr->XYZ()[0],
					  spacepoint_ptr->XYZ()[1],
					  spacepoint_ptr->XYZ()[2],
					  spacepoint_ptr->ErrXYZ()[0],
					  spacepoint_ptr->ErrXYZ()[1],
					  spacepoint_ptr->ErrXYZ()[2],
					  spacepoint_ptr->Chisq() );
      
      // Store address map for downstream association
      //fPtrIndex_sps[spacepoint_ptr] = std::make_pair(lite_data->size(),name_index);
      
      lite_data->push_back(lite_spacepoint);
    }
  }
  
  template <>
  void ScannerAlgo::ScanData(art::Handle<std::vector< ::recob::Track> > const &dh,
			     ::larlite::event_base* lite_dh)
  { 
    fDataReadFlag_v[lite_dh->data_type()][lite_dh->name()] = true;  
    //auto name_index = NameIndex(lite_dh->data_type(),lite_dh->name());
    auto lite_data = (::larlite::event_track*)lite_dh;
    for(size_t i=0; i<dh->size(); ++i) {
      
      art::Ptr<::recob::Track> track_ptr(dh,i);
      
      larlite::track track_lite;
      track_lite.set_track_id ( track_ptr->ID()   );
      // Direction & points
      for(size_t i=0; i<track_ptr->NumberTrajectoryPoints(); i++) {
	track_lite.add_vertex     (track_ptr->LocationAtPoint(i));
	track_lite.add_direction  (track_ptr->DirectionAtPoint(i));
      }
      // Covariance
      for(size_t i=0; i<track_ptr->NumberCovariance(); i++)
	track_lite.add_covariance (track_ptr->CovarianceAtPoint(i));
      // Momentum
      for(size_t i=0; i<track_ptr->NumberFitMomentum(); i++)
	track_lite.add_momentum   (track_ptr->MomentumAtPoint(i));
      
      // Store address map for downstream association
      //fPtrIndex_track[track_ptr] = std::make_pair(lite_data->size(),name_index);
      
      lite_data->push_back(track_lite);
    }
  }

  template <>
  void ScannerAlgo::ScanData(art::Handle<std::vector< ::recob::Shower> > const &dh,
			     ::larlite::event_base* lite_dh)
  { 
    fDataReadFlag_v[lite_dh->data_type()][lite_dh->name()] = true;  
    //auto name_index = NameIndex(lite_dh->data_type(),lite_dh->name());
    auto lite_data = (::larlite::event_shower*)lite_dh;
    for(size_t i=0; i<dh->size(); ++i) {
      
      const art::Ptr<::recob::Shower> shower_ptr(dh,i);
      
      larlite::shower lite_shower;

      lite_shower.set_id(shower_ptr->ID());
      lite_shower.set_total_energy_v(shower_ptr->Energy());
      lite_shower.set_total_energy_err_v(shower_ptr->EnergyErr());
      lite_shower.set_total_MIPenergy_v(shower_ptr->MIPEnergy());
      lite_shower.set_total_MIPenergy_err_v(shower_ptr->MIPEnergyErr());
      lite_shower.set_total_best_plane(shower_ptr->best_plane());
      lite_shower.set_direction(shower_ptr->Direction());
      lite_shower.set_direction_err(shower_ptr->DirectionErr());
      lite_shower.set_start_point(shower_ptr->ShowerStart());
      lite_shower.set_start_point_err(shower_ptr->ShowerStartErr());
      lite_shower.set_dedx_v(shower_ptr->dEdx());
      lite_shower.set_dedx_err_v(shower_ptr->dEdxErr());
      lite_shower.set_length(shower_ptr->Length());
      lite_shower.set_direction(shower_ptr->Direction());
      lite_shower.set_direction_err(shower_ptr->DirectionErr());

      //fPtrIndex_shower[shower_ptr] = std::make_pair(lite_data->size(),name_index);
      
      lite_data->push_back(lite_shower);
    }
  }

  template <>
  void ScannerAlgo::ScanData(art::Handle<std::vector< ::recob::Vertex> > const &dh,
			     ::larlite::event_base* lite_dh)
  { 
    fDataReadFlag_v[lite_dh->data_type()][lite_dh->name()] = true;  
    //auto name_index = NameIndex(lite_dh->data_type(),lite_dh->name());
    auto lite_data = (::larlite::event_vertex*)lite_dh;
    Double_t xyz[3]={0};
    
    for(size_t i=0; i < dh->size(); ++i) {
      
      const art::Ptr<::recob::Vertex> vtx_ptr(dh,i);
      
      vtx_ptr->XYZ(xyz);
      
      // get vertex info
      larlite::vertex lite_vtx(xyz,
			       vtx_ptr->ID());
      
      //fPtrIndex_vertex[vtx_ptr] = std::make_pair(lite_data->size(),name_index);
      
      // Store data
      lite_data->push_back(lite_vtx);
      
    }
  }

  template <>
  void ScannerAlgo::ScanData(art::Handle<std::vector< ::anab::Calorimetry> > const &dh,
			     ::larlite::event_base* lite_dh)
  { 
    fDataReadFlag_v[lite_dh->data_type()][lite_dh->name()] = true;  
    //auto name_index = NameIndex(lite_dh->data_type(),lite_dh->name());
    auto lite_data = (::larlite::event_calorimetry*)lite_dh;
    for(size_t i=0; i < dh->size(); ++i) {
      
      const art::Ptr<::anab::Calorimetry> calo_ptr(dh,i);
      
      larlite::calorimetry lite_calo;
      
      lite_calo.set_dedx(calo_ptr->dEdx());
      lite_calo.set_dqdx(calo_ptr->dQdx());
      lite_calo.set_residual_range(calo_ptr->ResidualRange());
      lite_calo.set_deadwire_range(calo_ptr->DeadWireResRC());
      lite_calo.set_kinetic_energy(calo_ptr->KineticEnergy());
      lite_calo.set_range(calo_ptr->Range());
      lite_calo.set_track_pitch(calo_ptr->TrkPitchVec());
      lite_calo.set_plane_id( larlite::geo::PlaneID( calo_ptr->PlaneID().Cryostat,
						     calo_ptr->PlaneID().TPC,
						     calo_ptr->PlaneID().Plane ) );

      //fPtrIndex_calo[calo_ptr] = std::make_pair(lite_data->size(),name_index);
      
      lite_data->push_back(lite_calo);
      
    }
  }

  template <>
  void ScannerAlgo::ScanData(art::Handle<std::vector< ::anab::ParticleID> > const &dh,
			     ::larlite::event_base* lite_dh)
  { 
    fDataReadFlag_v[lite_dh->data_type()][lite_dh->name()] = true;  
    //auto name_index = NameIndex(lite_dh->data_type(),lite_dh->name());
    auto lite_data = (::larlite::event_partid*)lite_dh;
    
    for(size_t i=0; i<dh->size(); ++i) {
      
      art::Ptr<::anab::ParticleID> partid_ptr(dh,i);

      auto const& pid = partid_ptr->PlaneID();
      
      larlite::partid lite_partid( partid_ptr->Pdg(),
				   partid_ptr->Ndf(),
				   partid_ptr->MinChi2(),
				   partid_ptr->DeltaChi2(),
				   partid_ptr->Chi2Proton(),
				   partid_ptr->Chi2Kaon(),
				   partid_ptr->Chi2Pion(),
				   partid_ptr->Chi2Muon(),
				   partid_ptr->MissingE(),
				   partid_ptr->MissingEavg(),
				   partid_ptr->PIDA(),
				   ::larlite::geo::PlaneID(pid.Cryostat,pid.TPC,pid.Plane));
      
      //fPtrIndex_partid[partid_ptr] = std::make_pair(lite_data->size(),name_index);
      
      lite_data->push_back(lite_partid);
    }
  }

  template <>
  void ScannerAlgo::ScanData(art::Handle<std::vector< ::recob::PFParticle> > const &dh,
			     ::larlite::event_base* lite_dh)
  { 
    fDataReadFlag_v[lite_dh->data_type()][lite_dh->name()] = true;  
    //auto name_index = NameIndex(lite_dh->data_type(),lite_dh->name());
    auto lite_data = (::larlite::event_pfpart*)lite_dh;
    for(size_t i=0; i<dh->size(); ++i) {
      
      art::Ptr<::recob::PFParticle> pfpart_ptr(dh,i);
      
      larlite::pfpart lite_pfpart(pfpart_ptr->PdgCode(),
				  pfpart_ptr->Self(),
				  pfpart_ptr->Parent(),
				  pfpart_ptr->Daughters());
      
      //fPtrIndex_pfpart[pfpart_ptr] = std::make_pair(lite_data->size(),name_index);
      
      lite_data->push_back(lite_pfpart);
    }
  }

  template <>
  void ScannerAlgo::ScanData(art::Handle<std::vector< ::recob::PCAxis> > const &dh,
			     ::larlite::event_base* lite_dh)
  { 
    fDataReadFlag_v[lite_dh->data_type()][lite_dh->name()] = true;  
    //auto name_index = NameIndex(lite_dh->data_type(),lite_dh->name());
    auto lite_data = (::larlite::event_pcaxis*)lite_dh;
    for(size_t i=0; i<dh->size(); ++i) {
      
      art::Ptr<::recob::PCAxis> pcaxis_ptr(dh,i);

      larlite::pcaxis lite_pcaxis( pcaxis_ptr->getSvdOK(),
				   pcaxis_ptr->getNumHitsUsed(),
				   pcaxis_ptr->getEigenValues(),
				   pcaxis_ptr->getEigenVectors(),
				   pcaxis_ptr->getAvePosition(),
				   pcaxis_ptr->getAveHitDoca(),
				   pcaxis_ptr->getID() );
      
      lite_data->push_back(lite_pcaxis);
    }
  }

  template <>
  void ScannerAlgo::ScanData(art::Handle<std::vector< ::anab::FlashMatch> > const &dh,
			       ::larlite::event_base* lite_dh)
  { 
    fDataReadFlag_v[lite_dh->data_type()][lite_dh->name()] = true;  
    //auto name_index = NameIndex(lite_dh->data_type(),lite_dh->name());
    auto lite_data = (::larlite::event_flashmatch*)lite_dh;
    for(size_t i=0; i<dh->size(); ++i) {
      
      art::Ptr<::anab::FlashMatch> fmatch_ptr(dh,i);

      larlite::flashmatch lite_fmatch( fmatch_ptr->Chi2(),
				       fmatch_ptr->FlashID(),
				       fmatch_ptr->SubjectID(),
				       fmatch_ptr->InBeam() );
      
      lite_data->push_back(lite_fmatch);
    }
  }


  template <class T>
  void ScanData(art::Handle<std::vector<T> > const &dh,
		::larlite::event_base* lite_dh)
  { throw cet::exception(__PRETTY_FUNCTION__)<<"Not implemented!"; }


  //
  // PtrMap key generator
  //
  template <> void ScannerAlgo::ProducePtrMapKey(const art::Ptr< ::recob::Hit>& ptr, size_t& key1, size_t& key2)
  {
    key1 = (size_t)(ptr->WireID().Plane);
    key2 = (size_t)(ptr->WireID().Wire) / 200;
  }

  template <> void ScannerAlgo::ProducePtrMapKey(const art::Ptr< ::recob::OpHit>& ptr, size_t& key1, size_t& key2)
  {
    key1 = ptr->OpChannel();
    key2 = 0;
    //key2 = (size_t)(ptr->PeakTime() / 1000.);
  }

  template <> void ScannerAlgo::ProducePtrMapKey(const art::Ptr< ::recob::Cluster>& ptr, size_t& key1, size_t& key2)
  {
    key1 = 0;
    key2 = (size_t)(ptr->Plane().Plane);
  }

  template <class T>
  void ScannerAlgo::ProducePtrMapKey(const art::Ptr<T>& ptr,size_t& key1, size_t&key2)
  {
    if(ptr)
      key1 = key2 = 0;
    else
      throw std::exception();
  }

  //
  // Getter for associated data product pointer 
  //
  template <> std::map<art::Ptr< ::simb::MCTruth>,std::pair<size_t,size_t> >& ScannerAlgo::GetPtrMap(size_t key1, size_t key2)
  { if(fPtrIndex_mctruth.size()<=key1) fPtrIndex_mctruth.resize(key1+1);
    if(fPtrIndex_mctruth[key1].size()<=key2) fPtrIndex_mctruth[key1].resize(key2+1);
    return fPtrIndex_mctruth[key1][key2]; 
  }

  template <> std::map<art::Ptr< ::simb::GTruth>,std::pair<size_t,size_t> >& ScannerAlgo::GetPtrMap(size_t key1, size_t key2)
  { if(fPtrIndex_gtruth.size()<=key1) fPtrIndex_gtruth.resize(key1+1);
    if(fPtrIndex_gtruth[key1].size()<=key2) fPtrIndex_gtruth[key1].resize(key2+1);
    return fPtrIndex_gtruth[key1][key2]; 
  }

  template <> std::map<art::Ptr< ::simb::MCFlux>,std::pair<size_t,size_t> >& ScannerAlgo::GetPtrMap(size_t key1, size_t key2)
  { if(fPtrIndex_mcflux.size()<=key1) fPtrIndex_mcflux.resize(key1+1);
    if(fPtrIndex_mcflux[key1].size()<=key2) fPtrIndex_mcflux[key1].resize(key2+1);
    return fPtrIndex_mcflux[key1][key2]; 
  }

  template <> std::map<art::Ptr< ::simb::MCParticle>,std::pair<size_t,size_t> >& ScannerAlgo::GetPtrMap(size_t key1, size_t key2)
  { if(fPtrIndex_mcpart.size()<=key1) fPtrIndex_mcpart.resize(key1+1);
    if(fPtrIndex_mcpart[key1].size()<=key2) fPtrIndex_mcpart[key1].resize(key2+1);
    return fPtrIndex_mcpart[key1][key2]; 
  }

  template <> std::map<art::Ptr< ::sim::SimChannel>,std::pair<size_t,size_t> >& ScannerAlgo::GetPtrMap(size_t key1, size_t key2)
  { if(fPtrIndex_simch.size()<=key1) fPtrIndex_simch.resize(key1+1);
    if(fPtrIndex_simch[key1].size()<=key2) fPtrIndex_simch[key1].resize(key2+1);
    return fPtrIndex_simch[key1][key2]; 
  }

  template <> std::map<art::Ptr< ::sim::AuxDetSimChannel>,std::pair<size_t,size_t> >& ScannerAlgo::GetPtrMap(size_t key1, size_t key2)
  { if(fPtrIndex_auxsimch.size()<=key1) fPtrIndex_auxsimch.resize(key1+1);
    if(fPtrIndex_auxsimch[key1].size()<=key2) fPtrIndex_auxsimch[key1].resize(key2+1);
    return fPtrIndex_auxsimch[key1][key2]; 
  }

  template <> std::map<art::Ptr< ::sim::MCShower>,std::pair<size_t,size_t> >& ScannerAlgo::GetPtrMap(size_t key1, size_t key2)
  { if(fPtrIndex_mcshower.size()<=key1) fPtrIndex_mcshower.resize(key1+1);
    if(fPtrIndex_mcshower[key1].size()<=key2) fPtrIndex_mcshower[key1].resize(key2+1);
    return fPtrIndex_mcshower[key1][key2]; 
  }

  template <> std::map<art::Ptr< ::sim::MCTrack>,std::pair<size_t,size_t> >& ScannerAlgo::GetPtrMap(size_t key1, size_t key2)
  { if(fPtrIndex_mctrack.size()<=key1) fPtrIndex_mctrack.resize(key1+1);
    if(fPtrIndex_mctrack[key1].size()<=key2) fPtrIndex_mctrack[key1].resize(key2+1);
    return fPtrIndex_mctrack[key1][key2]; 
  }

  template <> std::map<art::Ptr< ::recob::OpHit>,std::pair<size_t,size_t> >& ScannerAlgo::GetPtrMap(size_t key1, size_t key2) 
  { if(fPtrIndex_ophit.size()<=key1) fPtrIndex_ophit.resize(key1+1);
    if(fPtrIndex_ophit[key1].size()<=key2) fPtrIndex_ophit[key1].resize(key2+1);
    return fPtrIndex_ophit[key1][key2]; 
  }

  template <> std::map<art::Ptr< ::recob::OpFlash>,std::pair<size_t,size_t> >& ScannerAlgo::GetPtrMap(size_t key1, size_t key2)
  { if(fPtrIndex_opflash.size()<=key1) fPtrIndex_opflash.resize(key1+1);
    if(fPtrIndex_opflash[key1].size()<=key2) fPtrIndex_opflash[key1].resize(key2+1);
    return fPtrIndex_opflash[key1][key2]; 
  }

  template <> std::map<art::Ptr< ::recob::Hit>,std::pair<size_t,size_t> >& ScannerAlgo::GetPtrMap(size_t key1, size_t key2)
  { if(fPtrIndex_hit.size()<=key1) fPtrIndex_hit.resize(key1+1);
    if(fPtrIndex_hit[key1].size()<=key2) fPtrIndex_hit[key1].resize(key2+1);
    return fPtrIndex_hit[key1][key2]; 
  }

  template <> std::map<art::Ptr< ::raw::RawDigit>,std::pair<size_t,size_t> >& ScannerAlgo::GetPtrMap(size_t key1, size_t key2)
  { if(fPtrIndex_rawdigit.size()<=key1) fPtrIndex_rawdigit.resize(key1+1);
    if(fPtrIndex_rawdigit[key1].size()<=key2) fPtrIndex_rawdigit[key1].resize(key2+1);
    return fPtrIndex_rawdigit[key1][key2]; 
  }

  template <> std::map<art::Ptr< ::raw::OpDetWaveform>,std::pair<size_t,size_t> >& ScannerAlgo::GetPtrMap(size_t key1, size_t key2)
  { if(fPtrIndex_opdigit.size()<=key1) fPtrIndex_opdigit.resize(key1+1);
    if(fPtrIndex_opdigit[key1].size()<=key2) fPtrIndex_opdigit[key1].resize(key2+1);
    return fPtrIndex_opdigit[key1][key2]; 
  }

  template <> std::map<art::Ptr< ::raw::Trigger>,std::pair<size_t,size_t> >& ScannerAlgo::GetPtrMap(size_t key1, size_t key2)
  { if(fPtrIndex_trigger.size()<=key1) fPtrIndex_trigger.resize(key1+1);
    if(fPtrIndex_trigger[key1].size()<=key2) fPtrIndex_trigger[key1].resize(key2+1);
    return fPtrIndex_trigger[key1][key2]; 
  }

  template <> std::map<art::Ptr< ::recob::Wire>,std::pair<size_t,size_t> >& ScannerAlgo::GetPtrMap(size_t key1, size_t key2)
  { if(fPtrIndex_wire.size()<=key1) fPtrIndex_wire.resize(key1+1);
    if(fPtrIndex_wire[key1].size()<=key2) fPtrIndex_wire[key1].resize(key2+1);
    return fPtrIndex_wire[key1][key2]; 
  }

  template <> std::map<art::Ptr< ::recob::Cluster>,std::pair<size_t,size_t> >& ScannerAlgo::GetPtrMap(size_t key1, size_t key2)
  { if(fPtrIndex_cluster.size()<=key1) fPtrIndex_cluster.resize(key1+1);
    if(fPtrIndex_cluster[key1].size()<=key2) fPtrIndex_cluster[key1].resize(key2+1);
    return fPtrIndex_cluster[key1][key2]; 
  }

  template <> std::map<art::Ptr< ::recob::Seed>,std::pair<size_t,size_t> >& ScannerAlgo::GetPtrMap(size_t key1, size_t key2)
  { if(fPtrIndex_seed.size()<=key1) fPtrIndex_seed.resize(key1+1);
    if(fPtrIndex_seed[key1].size()<=key2) fPtrIndex_seed[key1].resize(key2+1);
    return fPtrIndex_seed[key1][key2]; 
  }

  template <> std::map<art::Ptr< ::recob::Track>,std::pair<size_t,size_t> >& ScannerAlgo::GetPtrMap(size_t key1, size_t key2)
  { if(fPtrIndex_track.size()<=key1) fPtrIndex_track.resize(key1+1);
    if(fPtrIndex_track[key1].size()<=key2) fPtrIndex_track[key1].resize(key2+1);
    return fPtrIndex_track[key1][key2]; 
  }

  template <> std::map<art::Ptr< ::recob::Shower>,std::pair<size_t,size_t> >& ScannerAlgo::GetPtrMap(size_t key1, size_t key2)
  { if(fPtrIndex_shower.size()<=key1) fPtrIndex_shower.resize(key1+1);
    if(fPtrIndex_shower[key1].size()<=key2) fPtrIndex_shower[key1].resize(key2+1);
    return fPtrIndex_shower[key1][key2]; 
  }

  template <> std::map<art::Ptr< ::recob::Vertex>,std::pair<size_t,size_t> >& ScannerAlgo::GetPtrMap(size_t key1, size_t key2)
  { if(fPtrIndex_vertex.size()<=key1) fPtrIndex_vertex.resize(key1+1);
    if(fPtrIndex_vertex[key1].size()<=key2) fPtrIndex_vertex[key1].resize(key2+1);
    return fPtrIndex_vertex[key1][key2]; 
  }

  template <> std::map<art::Ptr< ::recob::SpacePoint>,std::pair<size_t,size_t> >& ScannerAlgo::GetPtrMap(size_t key1, size_t key2)
  { if(fPtrIndex_sps.size()<=key1) fPtrIndex_sps.resize(key1+1);
    if(fPtrIndex_sps[key1].size()<=key2) fPtrIndex_sps[key1].resize(key2+1);
    return fPtrIndex_sps[key1][key2]; 
  }

  template <> std::map<art::Ptr< ::recob::EndPoint2D>,std::pair<size_t,size_t> >& ScannerAlgo::GetPtrMap(size_t key1, size_t key2)
  { if(fPtrIndex_end2d.size()<=key1) fPtrIndex_end2d.resize(key1+1);
    if(fPtrIndex_end2d[key1].size()<=key2) fPtrIndex_end2d[key1].resize(key2+1);
    return fPtrIndex_end2d[key1][key2]; 
  }

  template <> std::map<art::Ptr< ::anab::CosmicTag>,std::pair<size_t,size_t> >& ScannerAlgo::GetPtrMap(size_t key1, size_t key2)
  { if(fPtrIndex_cosmictag.size()<=key1) fPtrIndex_cosmictag.resize(key1+1);
    if(fPtrIndex_cosmictag[key1].size()<=key2) fPtrIndex_cosmictag[key1].resize(key2+1);
    return fPtrIndex_cosmictag[key1][key2]; 
  }

  template <> std::map<art::Ptr< ::anab::Calorimetry>,std::pair<size_t,size_t> >& ScannerAlgo::GetPtrMap(size_t key1, size_t key2)
  { if(fPtrIndex_calo.size()<=key1) fPtrIndex_calo.resize(key1+1);
    if(fPtrIndex_calo[key1].size()<=key2) fPtrIndex_calo[key1].resize(key2+1);
    return fPtrIndex_calo[key1][key2]; 
  }

  template <> std::map<art::Ptr< ::anab::ParticleID>,std::pair<size_t,size_t> >& ScannerAlgo::GetPtrMap(size_t key1, size_t key2)
  { if(fPtrIndex_partid.size()<=key1) fPtrIndex_partid.resize(key1+1);
    if(fPtrIndex_partid[key1].size()<=key2) fPtrIndex_partid[key1].resize(key2+1);
    return fPtrIndex_partid[key1][key2]; 
  }

  template <> std::map<art::Ptr< ::recob::PFParticle>,std::pair<size_t,size_t> >& ScannerAlgo::GetPtrMap(size_t key1, size_t key2)
  { if(fPtrIndex_pfpart.size()<=key1) fPtrIndex_pfpart.resize(key1+1);
    if(fPtrIndex_pfpart[key1].size()<=key2) fPtrIndex_pfpart[key1].resize(key2+1);
    return fPtrIndex_pfpart[key1][key2]; 
  }

  template <> std::map<art::Ptr< ::recob::PCAxis>,std::pair<size_t,size_t> >& ScannerAlgo::GetPtrMap(size_t key1, size_t key2)
  { if(fPtrIndex_pcaxis.size()<=key1) fPtrIndex_pcaxis.resize(key1+1);
    if(fPtrIndex_pcaxis[key1].size()<=key2) fPtrIndex_pcaxis[key1].resize(key2+1);
    return fPtrIndex_pcaxis[key1][key2]; 
  }

  template <> std::map<art::Ptr< ::anab::FlashMatch>,std::pair<size_t,size_t> >& ScannerAlgo::GetPtrMap(size_t key1, size_t key2)
  { if(fPtrIndex_fmatch.size()<=key1) fPtrIndex_fmatch.resize(key1+1);
    if(fPtrIndex_fmatch[key1].size()<=key2) fPtrIndex_fmatch[key1].resize(key2+1);
    return fPtrIndex_fmatch[key1][key2]; 
  }

  template <class T>
  std::map<art::Ptr<T>,std::pair<size_t,size_t> >& ScannerAlgo::GetPtrMap(size_t key1, size_t key2)
  { throw cet::exception(__PRETTY_FUNCTION__)<<"Not implemented for a specified data product type..."; }

  //
  // Type identifier functions
  //
  // simb
  template <> const ::larlite::data::DataType_t ScannerAlgo::LiteDataType<::simb::GTruth> () const
  { return ::larlite::data::kGTruth; }
  template <> const ::larlite::data::DataType_t ScannerAlgo::LiteDataType<::simb::MCTruth> () const
  { return ::larlite::data::kMCTruth; }
  template <> const ::larlite::data::DataType_t ScannerAlgo::LiteDataType<::simb::MCParticle> () const
  { return ::larlite::data::kMCParticle; }
  template <> const ::larlite::data::DataType_t ScannerAlgo::LiteDataType<::simb::MCFlux> () const
  { return ::larlite::data::kMCFlux; }
  template <> const ::larlite::data::DataType_t ScannerAlgo::LiteDataType<::simb::MCTrajectory> () const
  { return ::larlite::data::kMCTrajectory; }
  template <> const ::larlite::data::DataType_t ScannerAlgo::LiteDataType<::simb::MCNeutrino> () const
  { return ::larlite::data::kMCNeutrino; }
  // sim
  template <> const ::larlite::data::DataType_t ScannerAlgo::LiteDataType<::sim::SimPhotonsCollection> () const
  { return ::larlite::data::kSimPhotons; }
  template <> const ::larlite::data::DataType_t ScannerAlgo::LiteDataType<::sim::SimChannel> () const
  { return ::larlite::data::kSimChannel; }
  template <> const ::larlite::data::DataType_t ScannerAlgo::LiteDataType<::sim::AuxDetSimChannel> () const
  { return ::larlite::data::kAuxDetSimChannel; }
  template <> const ::larlite::data::DataType_t ScannerAlgo::LiteDataType<::sim::MCShower> () const
  { return ::larlite::data::kMCShower; }
  template <> const ::larlite::data::DataType_t ScannerAlgo::LiteDataType<::sim::MCTrack> () const
  { return ::larlite::data::kMCTrack; }

  // raw
  template <> const ::larlite::data::DataType_t ScannerAlgo::LiteDataType<::raw::RawDigit> () const
  { return ::larlite::data::kRawDigit; }
  template <> const ::larlite::data::DataType_t ScannerAlgo::LiteDataType<::raw::OpDetWaveform> () const
  { return ::larlite::data::kOpDetWaveform; }
  template <> const ::larlite::data::DataType_t ScannerAlgo::LiteDataType<::raw::Trigger> () const
  { return ::larlite::data::kTrigger; }
  // recob
  template <> const ::larlite::data::DataType_t ScannerAlgo::LiteDataType<::recob::Wire> () const
  { return ::larlite::data::kWire; }
  template <> const ::larlite::data::DataType_t ScannerAlgo::LiteDataType<::recob::Hit> () const
  { return ::larlite::data::kHit; }
  template <> const ::larlite::data::DataType_t ScannerAlgo::LiteDataType<::recob::Cluster> () const
  { return ::larlite::data::kCluster; }
  template <> const ::larlite::data::DataType_t ScannerAlgo::LiteDataType<::recob::SpacePoint> () const
  { return ::larlite::data::kSpacePoint; }
  template <> const ::larlite::data::DataType_t ScannerAlgo::LiteDataType<::recob::OpHit> () const
  { return ::larlite::data::kOpHit; }
  template <> const ::larlite::data::DataType_t ScannerAlgo::LiteDataType<::recob::OpFlash> () const
  { return ::larlite::data::kOpFlash; }
  template <> const ::larlite::data::DataType_t ScannerAlgo::LiteDataType<::recob::Seed> () const
  { return ::larlite::data::kSeed; }
  template <> const ::larlite::data::DataType_t ScannerAlgo::LiteDataType<::recob::Track> () const
  { return ::larlite::data::kTrack; }
  template <> const ::larlite::data::DataType_t ScannerAlgo::LiteDataType<::recob::Shower> () const
  { return ::larlite::data::kShower; }
  template <> const ::larlite::data::DataType_t ScannerAlgo::LiteDataType<::recob::Vertex> () const
  { return ::larlite::data::kVertex; }
  template <> const ::larlite::data::DataType_t ScannerAlgo::LiteDataType<::recob::EndPoint2D> () const
  { return ::larlite::data::kEndPoint2D; }
  template <> const ::larlite::data::DataType_t ScannerAlgo::LiteDataType<::recob::PFParticle> () const
  { return ::larlite::data::kPFParticle; }
  template <> const ::larlite::data::DataType_t ScannerAlgo::LiteDataType<::recob::PCAxis> () const
  { return ::larlite::data::kPCAxis; }
  // anab
  template <> const ::larlite::data::DataType_t ScannerAlgo::LiteDataType<::anab::CosmicTag> () const
  { return ::larlite::data::kCosmicTag; }
  template <> const ::larlite::data::DataType_t ScannerAlgo::LiteDataType<::anab::Calorimetry> () const
  { return ::larlite::data::kCalorimetry; }
  template <> const ::larlite::data::DataType_t ScannerAlgo::LiteDataType<::anab::ParticleID> () const
  { return ::larlite::data::kParticleID; }
  template <> const ::larlite::data::DataType_t ScannerAlgo::LiteDataType<::anab::FlashMatch> () const
  { return ::larlite::data::kFlashMatch; }
  // MuCS
  template <> const ::larlite::data::DataType_t ScannerAlgo::LiteDataType<::MuCS::MuCSData> () const
  { return ::larlite::data::kMuCSData; }
  template <> const ::larlite::data::DataType_t ScannerAlgo::LiteDataType<::MuCS::MuCSRecoData> () const
  { return ::larlite::data::kMuCSReco; }
  //
  // LocateLiteProduct implementation
  //
  template <class T>
  bool ScannerAlgo::LocateLiteProduct(art::Ptr<T> const ptr,
				      std::pair<size_t,size_t> &loc)
  { 
    size_t key1, key2;
    ProducePtrMapKey(ptr,key1,key2);
    auto const& ptr_map = GetPtrMap<T>(key1,key2);
    auto id_iter = ptr_map.find(ptr);
    if(id_iter == ptr_map.end()) return false; 
    
    loc.first = (*id_iter).second.first;
    loc.second = (*id_iter).second.second;
    return true;
  }

  //
  // ScanAssociation implementation 
  //
  template <class T,class U>
  void ScannerAlgo::ScanAssociation(art::Event const& e,
				    art::Handle<std::vector<T> > &dh,
				    ::larlite::event_ass* lite_dh)
  { 
    art::FindManyP<U> ptr_coll_v(dh, e, lite_dh->name());
    auto ass_type_a = LiteDataType<T>();
    auto ass_type_b = LiteDataType<U>();

    larlite::product_id ass_id_a(ass_type_a,lite_dh->name());

    try{
      if(!ptr_coll_v.size()) return;
      const std::vector<art::Ptr<U> > ptr_coll = ptr_coll_v.at(0);
    }catch( art::Exception const& e){
      return;
    }
    // Instantiate association container. length = # of producers for associated data type
    std::vector< ::larlite::AssSet_t> ass_set_v(fAssModuleLabel_v[ass_type_b].size(),
						::larlite::AssSet_t());
    
    // Return if there's no data products stored for associated data type
    if(!(ass_set_v.size())) return;

    auto watch = TStopwatch();
    //std::cout<<"Scanning association start... " << GetPtrMap<U>().size() << " = pool size... " << std::endl;
    std::pair<size_t,size_t> lite_location;

    // Loop over origin data product vector, find associated objects and store association info
    for(size_t i=0; i<dh->size(); ++i) {
      watch.Start();
      //std::cout<<"Association for " << i << std::endl;
      const std::vector<art::Ptr<U> >& ptr_coll = ptr_coll_v.at(i);
      //std::cout<<"Retrieved pointer vector "<< watch.RealTime() << std::endl;
      watch.Start();
      // Association vector: one per associated data product producers
      std::vector< ::larlite::AssUnit_t> ass_unit_v(ass_set_v.size(),::larlite::AssUnit_t());

      for(auto& au : ass_unit_v) au.reserve(ptr_coll.size());
      //std::cout<<"Start loop for location... " << watch.RealTime() << std::endl;
      watch.Start();
      for(auto& art_ptr : ptr_coll) {

	if(!LocateLiteProduct(art_ptr,lite_location)) continue;

	ass_unit_v[lite_location.second].push_back(lite_location.first);
      }
      //std::cout<<"Located "<<ptr_coll.size()<<" "<<watch.RealTime()<<std::endl;

      // Alternative
      /*
      watch.Start();
      if(ptr_coll.size()) {
	
	auto const& first_ptr = ptr_coll.front();
	auto const& pid = first_ptr.id().productID();

	art::Handle< std::vector<U> > u_handle;
	e.get(pid,u_handle);

	std::vector<deque> 
      }
      */
      for(size_t i=0; i<ass_set_v.size(); ++i)

	ass_set_v[i].emplace_back(ass_unit_v[i]);
      
    } // end looping over origin data products
    //std::cout<<"Located a=>b "<<watch.RealTime()<<std::endl;
    
    // Store associations in larlite data products
    for(size_t i=0; i<ass_set_v.size(); ++i) {

      // Loop over in one association set, store if there's any
      for(auto const& ass_unit : ass_set_v[i]) {

	if(ass_unit.size()) {

	  auto const& ass_name = fAssModuleLabel_v[(size_t)ass_type_b][i];

	  larlite::product_id ass_id_b(ass_type_b,ass_name);
	  
	  lite_dh->set_association(ass_id_a,ass_id_b,ass_set_v[i]);

	  //std::cout<<"Listing association: "<<ass_id_a.second.c_str()<<" => "<<ass_id_b.second.c_str()<<std::endl;
	  //lite_dh->list_association();
	  break;
	}
      }// end looping over association set
    }// end looping over a vector of association set

  }


  template <> void ScannerAlgo::ScanAssociation <::recob::Cluster,::recob::Cluster>(art::Event const& e,
										    art::Handle<std::vector<::recob::Cluster> > &dh,
										    ::larlite::event_ass* lite_dh)
  { throw cet::exception(__PRETTY_FUNCTION__) << " not implemented!"; }

  template <> void ScannerAlgo::ScanAssociation <::recob::EndPoint2D,::recob::EndPoint2D>(art::Event const& e,
											  art::Handle<std::vector<::recob::EndPoint2D> > &dh,
											  ::larlite::event_ass* lite_dh)
  { throw cet::exception(__PRETTY_FUNCTION__) << " not implemented!"; }
  
  template <> void ScannerAlgo::ScanAssociation <::recob::Vertex,::recob::Vertex>(art::Event const& e,
										  art::Handle<std::vector<::recob::Vertex> > &dh,
										  ::larlite::event_ass* lite_dh)
  { throw cet::exception(__PRETTY_FUNCTION__) << " not implemented!"; }

  template <> void ScannerAlgo::ScanAssociation <::anab::CosmicTag,::anab::CosmicTag>(art::Event const& e,
										      art::Handle<std::vector<::anab::CosmicTag> > &dh,
										      ::larlite::event_ass* lite_dh)
  { throw cet::exception(__PRETTY_FUNCTION__) << " not implemented!"; }

  template <> void ScannerAlgo::ScanAssociation <::recob::SpacePoint,::recob::SpacePoint>(art::Event const& e,
											  art::Handle<std::vector<::recob::SpacePoint> > &dh,
											  ::larlite::event_ass* lite_dh)
  { throw cet::exception(__PRETTY_FUNCTION__) << " not implemented!"; }

  template <> void ScannerAlgo::ScanAssociation <::recob::Track,::recob::Track>(art::Event const& e,
										art::Handle<std::vector<::recob::Track> > &dh,
										::larlite::event_ass* lite_dh)
  { throw cet::exception(__PRETTY_FUNCTION__) << " not implemented!"; }

  template <> void ScannerAlgo::ScanAssociation <::recob::Shower,::recob::Shower>(art::Event const& e,
										  art::Handle<std::vector<::recob::Shower> > &dh,
										  ::larlite::event_ass* lite_dh)
  { throw cet::exception(__PRETTY_FUNCTION__) << " not implemented!"; }

  template <> void ScannerAlgo::ScanAssociation <::anab::Calorimetry,::anab::Calorimetry>(art::Event const& e,
											  art::Handle<std::vector<::anab::Calorimetry> > &dh,
											  ::larlite::event_ass* lite_dh)
  { throw cet::exception(__PRETTY_FUNCTION__) << " not implemented!"; }

  template <> void ScannerAlgo::ScanAssociation <::anab::ParticleID,::anab::ParticleID>(art::Event const& e,
											art::Handle<std::vector<::anab::ParticleID> > &dh,
											::larlite::event_ass* lite_dh)
  { throw cet::exception(__PRETTY_FUNCTION__) << " not implemented!"; }

  template <> void ScannerAlgo::ScanAssociation <::simb::MCTruth,::simb::MCTruth>(art::Event const& e,
										  art::Handle<std::vector<::simb::MCTruth> > &dh,
										  ::larlite::event_ass* lite_dh)
  { throw cet::exception(__PRETTY_FUNCTION__) << " not implemented!"; }
  
  template <> void ScannerAlgo::ScanAssociation <::simb::MCParticle,::simb::MCParticle>(art::Event const& e,
											art::Handle<std::vector<::simb::MCParticle> > &dh,
											::larlite::event_ass* lite_dh)
  { throw cet::exception(__PRETTY_FUNCTION__) << " not implemented!"; }

  template <> void ScannerAlgo::ScanAssociation <::recob::PFParticle,::recob::PFParticle>(art::Event const& e,
											  art::Handle<std::vector<::recob::PFParticle> > &dh,
											  ::larlite::event_ass* lite_dh)
  { throw cet::exception(__PRETTY_FUNCTION__) << " not implemented!"; }

  template <> void ScannerAlgo::ScanAssociation <::recob::PCAxis,::recob::PCAxis>(art::Event const& e,
										  art::Handle<std::vector<::recob::PCAxis> > &dh,
										  ::larlite::event_ass* lite_dh)
  { throw cet::exception(__PRETTY_FUNCTION__) << " not implemented!"; }


  //
  // LiteDataType
  //
  template <class T>
  const ::larlite::data::DataType_t ScannerAlgo::LiteDataType() const
  { throw cet::exception(__PRETTY_FUNCTION__)<<"Unsupported data type conversion attempted!";
    return ::larlite::data::kUndefined;
  }

  //
  // ProductID
  //
  template <class T>
  const ::larlite::product_id ScannerAlgo::ProductID(size_t name_index) const
  { auto data_type = LiteDataType<T>();
    if(fModuleLabel_v[(size_t)data_type].size() <= name_index)
      throw cet::exception(__PRETTY_FUNCTION__)
	<< "Length of registered products for data type " << ::larlite::data::kDATA_TREE_NAME[data_type].c_str()
	<< " is " << fModuleLabel_v[(size_t)data_type].size()
	<< " while you requested " << name_index;
    return ::larlite::product_id(data_type,fModuleLabel_v[(size_t)(data_type)][name_index]);
  }

  //
  // Associated ProductID
  //
  template <class T>
  const ::larlite::product_id ScannerAlgo::AssProductID(size_t name_index) const
  { auto data_type = LiteDataType<T>();
    if(fAssModuleLabel_v[(size_t)data_type].size() <= name_index)
      throw cet::exception(__PRETTY_FUNCTION__)
	<< "Length of registered products for data type " << ::larlite::data::kDATA_TREE_NAME[data_type].c_str()
	<< " is " << fAssModuleLabel_v[(size_t)data_type].size()
	<< " while you requested " << name_index;
    return ::larlite::product_id(data_type,fAssModuleLabel_v[(size_t)(data_type)][name_index]);
  }
}
#endif

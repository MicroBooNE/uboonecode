////////////////////////////////////////////////////////////////////////
// Class:       CosmicFlashTaggerAna
// Plugin Type: analyzer (art v2_05_00)
// File:        CosmicFlashTaggerAna_module.cc
//
// Generated at Fri Dec  9 09:44:39 2016 by Marco Del Tutto using cetskelgen
// from cetlib version v1_21_00.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
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

#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "larsim/MCCheater/BackTracker.h"
#include "lardataobj/AnalysisBase/CosmicTag.h"
#include "lardata/Utilities/AssociationUtil.h"

#include "lardataobj/RecoBase/PFParticle.h"
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"

#include "TString.h"
#include "TTree.h"

class CosmicFlashTaggerAna;


class CosmicFlashTaggerAna : public art::EDAnalyzer {
public:
  explicit CosmicFlashTaggerAna(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  CosmicFlashTaggerAna(CosmicFlashTaggerAna const &) = delete;
  CosmicFlashTaggerAna(CosmicFlashTaggerAna &&) = delete;
  CosmicFlashTaggerAna & operator = (CosmicFlashTaggerAna const &) = delete;
  CosmicFlashTaggerAna & operator = (CosmicFlashTaggerAna &&) = delete;

  // Required functions.
  void analyze(art::Event const & e) override;

private:

  std::string _hitfinderLabel;
  std::string _pfp_producer;
  std::string _geantModuleLabel;
  std::string _spacepointLabel;
  std::string _cosmic_tag_producer;
  bool _recursiveMatching = false;
  bool _debug = true;

  const simb::Origin_t NEUTRINO_ORIGIN = simb::kBeamNeutrino;

  /// Maps used for PFParticle truth matching
  typedef std::map< art::Ptr<recob::PFParticle>, unsigned int > RecoParticleToNMatchedHits;
  typedef std::map< art::Ptr<simb::MCParticle>,  RecoParticleToNMatchedHits > ParticleMatchingMap;
  typedef std::set< art::Ptr<recob::PFParticle> > PFParticleSet;
  typedef std::set< art::Ptr<simb::MCParticle> >  MCParticleSet;

  /**
   *  @brief Perform matching between true and reconstructed particles
   *
   *  @param recoParticlesToHits the mapping from reconstructed particles to hits
   *  @param trueHitsToParticles the mapping from hits to true particles
   *  @param matchedParticles the output matches between reconstructed and true particles
   *  @param matchedHits the output matches between reconstructed particles and hits
   */
  void GetRecoToTrueMatches(const lar_pandora::PFParticlesToHits &recoParticlesToHits, const lar_pandora::HitsToMCParticles &trueHitsToParticles,
       lar_pandora::MCParticlesToPFParticles &matchedParticles, lar_pandora::MCParticlesToHits &matchedHits) const;
   /**
   *  @brief Perform matching between true and reconstructed particles
   *
   *  @param recoParticlesToHits the mapping from reconstructed particles to hits
   *  @param trueHitsToParticles the mapping from hits to true particles
   *  @param matchedParticles the output matches between reconstructed and true particles
   *  @param matchedHits the output matches between reconstructed particles and hits
   *  @param recoVeto the veto list for reconstructed particles
   *  @param trueVeto the veto list for true particles
   */
  void GetRecoToTrueMatches(const lar_pandora::PFParticlesToHits &recoParticlesToHits, const lar_pandora::HitsToMCParticles &trueHitsToParticles,
               lar_pandora::MCParticlesToPFParticles &matchedParticles, lar_pandora::MCParticlesToHits &matchedHits, PFParticleSet &recoVeto, MCParticleSet &trueVeto) const;

  TTree* _tree1;
  int _run, _subrun, _event;
  int _ccnc, _pdg;
  int _nPFPtagged, _nuPFPwasTagged;
};


CosmicFlashTaggerAna::CosmicFlashTaggerAna(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p) 
{
  _pfp_producer            = p.get<std::string>("PFParticleProducer");
  _hitfinderLabel          = p.get<std::string>("HitProducer");
  _geantModuleLabel        = p.get<std::string>("GeantModule");
  _spacepointLabel         = p.get<std::string>("SpacePointProducer");
  _cosmic_tag_producer     = p.get<std::string>("CosmicTagProducer");

  art::ServiceHandle<art::TFileService> fs;
  _tree1 = fs->make<TTree>("tree","");
  _tree1->Branch("run",    &_run,    "run/I");
  _tree1->Branch("subrun", &_subrun, "subrun/I");
  _tree1->Branch("event",  &_event,  "event/I");
  _tree1->Branch("ccnc",   &_ccnc,   "ccnc/I");
  _tree1->Branch("pdg",    &_pdg,    "pdg/I");
  _tree1->Branch("nPFPtagged",     &_nPFPtagged,     "nPFPtagged/I");
  _tree1->Branch("nuPFPwasTagged", &_nuPFPwasTagged, "nuPFPwasTagged/I");
}

void CosmicFlashTaggerAna::analyze(art::Event const & e)
{
  _run    = e.id().run();
  _subrun = e.id().subRun();
  _event  = e.id().event();

  _ccnc   = -1;
  _pdg    = -1;

  art::ServiceHandle<cheat::BackTracker> bt;

  // *******************
  // Pandora MCParticle to PFParticle matching
  // *******************

  // --- Collect hits
  lar_pandora::HitVector hitVector;
  lar_pandora::LArPandoraHelper::CollectHits(e, _hitfinderLabel, hitVector);

  // --- Collect PFParticles and match Reco Particles to Hits
  lar_pandora::PFParticleVector  recoParticleVector;
  //lar_pandora::PFParticleVector  recoNeutrinoVector;
  lar_pandora::PFParticlesToHits recoParticlesToHits;
  lar_pandora::HitsToPFParticles recoHitsToParticles;

  lar_pandora::LArPandoraHelper::CollectPFParticles(e, _pfp_producer, recoParticleVector);
  //lar_pandora::LArPandoraHelper::SelectNeutrinoPFParticles(recoParticleVector, recoNeutrinoVector);
  lar_pandora::LArPandoraHelper::BuildPFParticleHitMaps(e, _pfp_producer, _spacepointLabel, recoParticlesToHits, recoHitsToParticles, lar_pandora::LArPandoraHelper::kAddDaughters);

  if (_debug)
    //std::cout << "  RecoNeutrinos: " << recoNeutrinoVector.size() << std::endl;

  if (_debug)
    std::cout << "  RecoParticles: " << recoParticleVector.size() << std::endl;

  // --- Collect MCParticles and match True Particles to Hits
  lar_pandora::MCParticleVector     trueParticleVector;
  lar_pandora::MCTruthToMCParticles truthToParticles;
  lar_pandora::MCParticlesToMCTruth particlesToTruth;
  lar_pandora::MCParticlesToHits    trueParticlesToHits;
  lar_pandora::HitsToMCParticles    trueHitsToParticles;

  if (!e.isRealData()) {
    lar_pandora::LArPandoraHelper::CollectMCParticles(e, _geantModuleLabel, trueParticleVector);
    lar_pandora::LArPandoraHelper::CollectMCParticles(e, _geantModuleLabel, truthToParticles, particlesToTruth);
    lar_pandora::LArPandoraHelper::BuildMCParticleHitMaps(e, _geantModuleLabel, hitVector, trueParticlesToHits, trueHitsToParticles, lar_pandora::LArPandoraHelper::kAddDaughters);
  }

  if (_debug)
    std::cout << "  TrueParticles: " << particlesToTruth.size() << std::endl;

  if (_debug)
    std::cout << "  TrueEvents: " << truthToParticles.size() << std::endl;

  lar_pandora::MCParticlesToPFParticles matchedParticles;    // This is a map: MCParticle to matched PFParticle
  lar_pandora::MCParticlesToHits        matchedParticleHits;

  // --- Do the matching
  this->GetRecoToTrueMatches(recoParticlesToHits, 
                             trueHitsToParticles, 
                             matchedParticles, 
                             matchedParticleHits);


  // *******************
  // CosmicTag analysis
  // *******************

  std::vector<art::Ptr<recob::PFParticle>> taggedPFP;
  std::vector<art::Ptr<recob::PFParticle>> neutrinoOriginPFP;

  // Loop over true particle and find the neutrino related ones
  for (lar_pandora::MCParticlesToPFParticles::const_iterator iter1 = matchedParticles.begin(), iterEnd1 = matchedParticles.end();
             iter1 != iterEnd1; ++iter1) {

     art::Ptr<simb::MCParticle>  mc_par = iter1->first;   // The MCParticle 
     art::Ptr<recob::PFParticle> pf_par = iter1->second;  // The matched PFParticle

     const art::Ptr<simb::MCTruth> mc_truth = bt->TrackIDToMCTruth(mc_par->TrackId());
     if (mc_truth->Origin() == NEUTRINO_ORIGIN) {
       if (_debug) {
         std::cout << "Neutrino related track found." << std::endl;
         std::cout << "Process (0==CC, 1==NC) " << mc_truth->GetNeutrino().CCNC()         << std::endl;
         std::cout << "Neutrino PDG           " << mc_truth->GetNeutrino().Nu().PdgCode() << std::endl;
         std::cout << "PDG  " << mc_par->PdgCode() << std::endl;
         std::cout << "Mass " << mc_par->Mass()    << std::endl;
         std::cout << "Proc " << mc_par->Process() << std::endl;
         std::cout << "Vx   " << mc_par->Vx()      << std::endl;
         std::cout << "Vy   " << mc_par->Vy()      << std::endl;
         std::cout << "Vz   " << mc_par->Vz()      << std::endl;
         std::cout << "T    " << mc_par->T()       << std::endl;
         double timeCorrection = 343.75;
         std::cout << "Remeber a time correction of " << timeCorrection << std::endl;
       }    
       if (_debug) {
         std::cout << "  The related PFP: " << std::endl;
         std::cout << "  has ID: " << pf_par->Self() << std::endl;
         //std::cout << "  Vx " << mc_par->Vx() << std::endl;
         //std::cout << "  Vy " << mc_par->Vy() << std::endl;
         //std::cout << "  Vz " << mc_par->Vz() << std::endl;
       }

       neutrinoOriginPFP.emplace_back(pf_par);

       _ccnc = mc_truth->GetNeutrino().CCNC();
       _pdg  = mc_truth->GetNeutrino().Nu().PdgCode();
     }
  }
  if (_debug) std::cout << "Neutrino related PFPs in this event: " << neutrinoOriginPFP.size() << std::endl;

  // Get the CosmicTag from the ART event
  art::Handle<std::vector<anab::CosmicTag>> cosmicTagHandle;
  e.getByLabel(_cosmic_tag_producer, cosmicTagHandle);

  if (!cosmicTagHandle.isValid() || cosmicTagHandle->empty()){
    std::cerr << "Cosmic tag is not valid or empty." << std::endl;
    return;
  }

  // Look up the associations to PFPs
  art::FindManyP<recob::PFParticle> cosmicPFPAssns(cosmicTagHandle, e, _cosmic_tag_producer);

  if (_debug) std::cout << "cosmicPFPAssns.size(): " << cosmicPFPAssns.size() << std::endl; 

  // Loop over the cosmic tags
  for (unsigned int ct = 0; ct < cosmicPFPAssns.size(); ct++) {

    // Get the cosmic tag
    art::Ptr<anab::CosmicTag> cosmicTag(cosmicTagHandle, ct);
    if(_debug) std::cout << "This cosmic tag (" << ct << ") has type: " << cosmicTag->CosmicType() << std::endl;

    // Get the PFP associated with this CT
    std::vector<art::Ptr<recob::PFParticle>> cosmicTagToPFP_v = cosmicPFPAssns.at(cosmicTag.key());
    if(_debug) std::cout << "Number of PFP associated with this Cosmic Tag: " << cosmicTagToPFP_v.size() << std::endl;

    if (cosmicTagToPFP_v.size() == 0){
      std::cerr << "No cosmic tags in this event." << std::endl;
      return;
    }

    if (cosmicTagToPFP_v.size() > 1){
      std::cerr << "More than 1 PFP associated with a single cosmic tag ?!" << std::endl;
      return;
    }

    if (cosmicTag->CosmicType() == anab::CosmicTagID_t::kFlash_BeamIncompatible) {
      taggedPFP.emplace_back(cosmicTagToPFP_v.at(0));
    }
  }

  if (_debug) std::cout << "Beam incompatible tagged PFPs in this event: " << taggedPFP.size() << std::endl;
  _nPFPtagged = -1;
  _nPFPtagged = taggedPFP.size();

  // Loop through the taggedPFP and see if there is a neutrino related one
  _nuPFPwasTagged = 0;
  for (unsigned int i = 0; i < taggedPFP.size(); i++) {
    for (unsigned int j = 0; j < neutrinoOriginPFP.size(); j++) {
      if(taggedPFP[i] == neutrinoOriginPFP[j]) {
        std::cout << ">>>>>>>>>>>>>>>>> A neutrino related PFP (with ID " << neutrinoOriginPFP[j]->Self() << ") was tagged as beam incompatible." << std::endl;
        _nuPFPwasTagged = 1;
      }
    }
  }

  _tree1->Fill();
}



//------------------------------------------------------------------------------------------------------------------------------------------

void CosmicFlashTaggerAna::GetRecoToTrueMatches(const lar_pandora::PFParticlesToHits &recoParticlesToHits, 
                                                const lar_pandora::HitsToMCParticles &trueHitsToParticles,
                                                lar_pandora::MCParticlesToPFParticles &matchedParticles, 
                                                lar_pandora::MCParticlesToHits &matchedHits) const
{   
  PFParticleSet recoVeto; MCParticleSet trueVeto;
    
  this->GetRecoToTrueMatches(recoParticlesToHits, trueHitsToParticles, matchedParticles, matchedHits, recoVeto, trueVeto);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void CosmicFlashTaggerAna::GetRecoToTrueMatches(const lar_pandora::PFParticlesToHits &recoParticlesToHits, 
                                                const lar_pandora::HitsToMCParticles &trueHitsToParticles,
                                                lar_pandora::MCParticlesToPFParticles &matchedParticles, 
                                                lar_pandora::MCParticlesToHits &matchedHits, 
                                                PFParticleSet &vetoReco, 
                                                MCParticleSet &vetoTrue) const
{
    bool foundMatches(false);

    for (lar_pandora::PFParticlesToHits::const_iterator iter1 = recoParticlesToHits.begin(), iterEnd1 = recoParticlesToHits.end();
        iter1 != iterEnd1; ++iter1)
    {
        const art::Ptr<recob::PFParticle> recoParticle = iter1->first;
        if (vetoReco.count(recoParticle) > 0)
            continue;

        const lar_pandora::HitVector &hitVector = iter1->second;

        lar_pandora::MCParticlesToHits truthContributionMap;

        for (lar_pandora::HitVector::const_iterator iter2 = hitVector.begin(), iterEnd2 = hitVector.end(); iter2 != iterEnd2; ++iter2)
        {
            const art::Ptr<recob::Hit> hit = *iter2;

            lar_pandora::HitsToMCParticles::const_iterator iter3 = trueHitsToParticles.find(hit);
            if (trueHitsToParticles.end() == iter3)
                continue;

            const art::Ptr<simb::MCParticle> trueParticle = iter3->second;
            if (vetoTrue.count(trueParticle) > 0)
                continue;

            truthContributionMap[trueParticle].push_back(hit);
        }

        lar_pandora::MCParticlesToHits::const_iterator mIter = truthContributionMap.end();

        for (lar_pandora::MCParticlesToHits::const_iterator iter4 = truthContributionMap.begin(), iterEnd4 = truthContributionMap.end();
            iter4 != iterEnd4; ++iter4)
        {
            if ((truthContributionMap.end() == mIter) || (iter4->second.size() > mIter->second.size()))
            {
                mIter = iter4;
            }
        }

        if (truthContributionMap.end() != mIter)
        {
            const art::Ptr<simb::MCParticle> trueParticle = mIter->first;

            lar_pandora::MCParticlesToHits::const_iterator iter5 = matchedHits.find(trueParticle);

            if ((matchedHits.end() == iter5) || (mIter->second.size() > iter5->second.size()))
            {
                matchedParticles[trueParticle] = recoParticle;
                matchedHits[trueParticle] = mIter->second;
                foundMatches = true;
            }
        }
    }

    if (!foundMatches)
        return;

    for (lar_pandora::MCParticlesToPFParticles::const_iterator pIter = matchedParticles.begin(), pIterEnd = matchedParticles.end();
        pIter != pIterEnd; ++pIter)
    {
        vetoTrue.insert(pIter->first);
        vetoReco.insert(pIter->second);
    }

    if (_recursiveMatching)
        this->GetRecoToTrueMatches(recoParticlesToHits, trueHitsToParticles, matchedParticles, matchedHits, vetoReco, vetoTrue);
}


DEFINE_ART_MODULE(CosmicFlashTaggerAna)

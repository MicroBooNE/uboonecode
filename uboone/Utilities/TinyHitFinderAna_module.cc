////////////////////////////////////////////////////////////////////////
//
// TinyHitFinderAna class
//
// echurch@fnal.gov
//
//  This algorithm is designed to analyze hits on wires after deconvolution and produce a handful of histograms.
////////////////////////////////////////////////////////////////////////
// Framework includes
#include "art/Framework/Core/ModuleMacros.h" 

#ifndef TINYHITFINDERANA_H
#define TINYHITFINDERANA_H

extern "C" {
#include <sys/types.h>
#include <sys/stat.h>
}

// ROOT includes
#include <TMath.h>
#include <TH1F.h>
#include <TGraph.h>
#include <TFile.h>

// C++ includes
#include <algorithm>
#include <sstream>
#include <fstream>
#include <bitset>

// Framework includes
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Persistency/Common/Ptr.h" 
#include "art/Persistency/Common/PtrVector.h" 
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// LArSoft includes
#include "Geometry/Geometry.h"
#include "Geometry/CryostatGeo.h"
#include "Geometry/TPCGeo.h"

#include "RecoBase/Hit.h"
#include "Utilities/LArFFT.h"
#include "RecoBase/Hit.h"

#include "TComplex.h"
#include "TString.h"
#include "TGraph.h"
#include "TH2.h"
#include "TTree.h"

#include "art/Framework/Core/EDAnalyzer.h"

#include <vector>
#include <string>

namespace geo { class Geometry; }

///Detector simulation of raw signals on wires
namespace hit {

  /// Base class for creation of raw signals on wires. 
  class TinyHitFinderAna : public art::EDAnalyzer {
    
  public:
        
    explicit TinyHitFinderAna(fhicl::ParameterSet const& pset); 
    virtual ~TinyHitFinderAna();
    
    /// read/write access to event
    void analyze (const art::Event& evt);
    void beginJob();
    void reconfigure(fhicl::ParameterSet const& p);

  private:

    std::string            fFFTTinyHitFinderModuleLabel;
    std::string            fLArG4ModuleLabel;
    

    TH1F* hfNp0;
    TH1F* hfNp1;
    TH1F* hfNp2;
    int fNp0;
    int fNp1;
    int fNp2;

  }; // class TinyHitFinderAna

}

namespace hit{

  //-------------------------------------------------
  TinyHitFinderAna::TinyHitFinderAna(fhicl::ParameterSet const& pset) 
    : EDAnalyzer(pset)
  {
    this->reconfigure(pset);
  }

  //-------------------------------------------------
  TinyHitFinderAna::~TinyHitFinderAna()
  {
  }

  void TinyHitFinderAna::reconfigure(fhicl::ParameterSet const& p)
  {
    fFFTTinyHitFinderModuleLabel = p.get< std::string >("HitsModuleLabel");
    fLArG4ModuleLabel        = p.get< std::string >("LArGeantModuleLabel");
    return;
  }
  //-------------------------------------------------
  void TinyHitFinderAna::beginJob() 
  {
    // get access to the TFile service
    art::ServiceHandle<art::TFileService> tfs;

    // 1000 max hits suffices for our isotropic ,1-0.5 GeV/c muons. 7000 is necessary for bnb+cosmics.
    hfNp0 = tfs->make<TH1F>("hitsPlane0", "hits on first plane" , 20, 0., 1000.);
    hfNp1 = tfs->make<TH1F>("hitsPlane1", "hits on second plane", 20, 0., 1000.);
    hfNp2 = tfs->make<TH1F>("hitsPlane2", "hits on third plane" , 20, 0., 1000.);

    return;

  }

  //-------------------------------------------------
  void TinyHitFinderAna::analyze(const art::Event& evt)
  {

    /*
    if (evt.isRealData()){
      throw cet::exception("TinyHitFinderAna: ") << "Not for use on Data yet...\n";
    }
    */


    art::Handle< std::vector<recob::Hit> > hitHandle;
    evt.getByLabel(fFFTTinyHitFinderModuleLabel,hitHandle);
    art::ServiceHandle<geo::Geometry> geom;  

    std::map<geo::PlaneID, std::vector< art::Ptr<recob::Hit> > > planeIDToHits;
    for(size_t h = 0; h < hitHandle->size(); ++h)
      planeIDToHits[hitHandle->at(h).WireID().planeID()].push_back(art::Ptr<recob::Hit>(hitHandle, h));
  
    fNp0 = 0;
    fNp1 = 0;
    fNp2 = 0;

    
    for(auto mapitr : planeIDToHits){
      
      geo::PlaneID pid = mapitr.first;
      auto itr = mapitr.second.begin();
      while(itr != mapitr.second.end()) {
	  
	//	fRun = evt.run();
	//	fEvt = evt.id().event();
	  
	
	if (pid.Plane == 0){
	  
	  ++fNp0;
	}
	  
	else if (pid.Plane == 1){
	  ++fNp1;
	}
	
	else if (pid.Plane == 2){
	  ++fNp2;
	}
	
	itr++;
      } // loop on Hits
    } // loop on map

    hfNp0->Fill(fNp0);
    hfNp1->Fill(fNp1);
    hfNp2->Fill(fNp2);
    //    std::cout << "TinyHitFinderAna_module: Plane 0,1,2 hits are " << fNp0 << ", " << fNp1 << ", " << fNp2 << "  ... if that can be believed." << std::endl;
    return;
  }//end analyze method
  
}//end namespace

namespace hit{

  DEFINE_ART_MODULE(TinyHitFinderAna)

} // end of hit namespace

#endif // TINYHITFINDERANA_H


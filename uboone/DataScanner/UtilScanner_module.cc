
#ifndef UtilScanner_H
#define UtilScanner_H

// Framework includes
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/FindManyP.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Registry/ServiceMacros.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "art/Persistency/Common/Ptr.h"
#include "art/Persistency/Common/PtrVector.h"

// LArSoft includes
#include "Geometry/Geometry.h"
#include "Geometry/CryostatGeo.h"
#include "Geometry/PlaneGeo.h"
#include "Geometry/OpDetGeo.h"
#include "Geometry/WireGeo.h"
#include "Geometry/Geometry.h"
#include "Utilities/LArProperties.h"
#include "Utilities/GeometryUtilities.h"
#include "Utilities/DetectorProperties.h"

#include "Base/Base-TypeDef.hh"


// ART includes.


#include <TTree.h>

namespace ana {
 
  class UtilScanner : public art::EDAnalyzer{
  public:
 
    UtilScanner(const fhicl::ParameterSet&);
    virtual ~UtilScanner();

    void beginJob();

    void analyze (const art::Event&); 

    /// Function to create utility TTrees
    void SaveUtilityData(fhicl::ParameterSet const& pset);
    void SaveGeometry(fhicl::ParameterSet const& pset);
    void SaveLArProperties(fhicl::ParameterSet const& pset);
    void SaveDetectorProperties(fhicl::ParameterSet const& pset);

  private:
    bool _util_saved;
    fhicl::ParameterSet _pset;
    TTree *_geom_tree;
    TTree *_detp_tree;
    TTree *_larp_tree;
  };

} 

#endif//  UtilScanner_H

// UtilScanner.cc

// Framework includes
#include "art/Framework/Core/ModuleMacros.h"

namespace ana {
  DEFINE_ART_MODULE(UtilScanner)
}

namespace ana {

  //-----------------------------------------------------------------------
  // Constructor
  UtilScanner::UtilScanner(fhicl::ParameterSet const& pset) : EDAnalyzer(pset)
  {
    _pset = pset;
    _util_saved=false;
    _larp_tree=0;
    _detp_tree=0;
    _geom_tree=0;
  }

  //-----------------------------------------------------------------------
  // Destructor
  UtilScanner::~UtilScanner(){}
   
  //-----------------------------------------------------------------------
  void UtilScanner::beginJob(){}
   

  //-----------------------------------------------------------------------
  void UtilScanner::analyze(const art::Event& evt) 
  {
    if(!_util_saved)
      SaveUtilityData(_pset);

    return;
  }

  void UtilScanner::SaveUtilityData(fhicl::ParameterSet const& pset)
  {
    SaveDetectorProperties(pset.get< fhicl::ParameterSet >("DetectorProperties"));
    SaveLArProperties(pset.get< fhicl::ParameterSet >("LArProperties"));
    SaveGeometry(pset.get< fhicl::ParameterSet >("Geometry"));
    _util_saved=true;
  }

  void UtilScanner::SaveGeometry(fhicl::ParameterSet const& pset)
  {
    if(_geom_tree) return;
    art::ServiceHandle<geo::Geometry> _geom;
    art::ServiceHandle<art::TFileService>  fileService;    
    TTree* _geom_tree = fileService->make<TTree>("Geometry","");

    //---Fill Variables ---//
    Double_t fDetLength = _geom->DetLength();
    Double_t fDetHalfWidth = _geom->DetHalfWidth();
    Double_t fDetHalfHeight = _geom->DetHalfHeight();

    Double_t fCryoLength = _geom->CryostatLength();
    Double_t fCryoHalfWidth = _geom->CryostatHalfWidth();
    Double_t fCryoHalfHeight = _geom->CryostatHalfHeight();

    Double_t start[3]={0.};
    Double_t end[3]={0.};
    std::vector<UChar_t>                fChannelToPlaneMap(_geom->Nchannels(),larlight::DATA::INVALID_UCHAR);
    std::vector<UShort_t>               fChannelToWireMap(_geom->Nchannels(),larlight::DATA::INVALID_USHORT);
    std::vector<std::vector<std::vector<Double_t> > > fWireStartVtx(_geom->Nplanes(),std::vector<std::vector<Double_t> >());
    std::vector<std::vector<std::vector<Double_t> > > fWireEndVtx(_geom->Nplanes(),std::vector<std::vector<Double_t> >());    
    for(size_t i=0; i<_geom->Nchannels(); ++i) {
      std::vector<geo::WireID> wids = _geom->ChannelToWire(i);
      fChannelToPlaneMap[i]=wids[0].Plane;
      fChannelToWireMap[i]=wids[0].Wire;
      if(!(fWireStartVtx.at(wids[0].Plane).size())) {
	fWireStartVtx.at(wids[0].Plane).resize(_geom->Nwires(wids[0].Plane),std::vector<double>(3,larlight::DATA::INVALID_DOUBLE));
	fWireEndVtx.at(wids[0].Plane).resize(_geom->Nwires(wids[0].Plane),std::vector<double>(3,larlight::DATA::INVALID_DOUBLE));
      }
      _geom->WireEndPoints(0,0,wids[0].Plane,wids[0].Wire,start,end);
      for(size_t coord =0; coord<3; ++coord) {
	fWireStartVtx.at(wids[0].Plane).at(wids[0].Wire).at(coord) = start[coord];
	fWireEndVtx.at(wids[0].Plane).at(wids[0].Wire).at(coord) = end[coord];
      }
    }

    // Vectors with length = # planes
    std::vector<std::vector<UShort_t> > fPlaneWireToChannelMap(_geom->Nplanes(),std::vector<UShort_t>());
    std::vector<larlight::GEO::SigType_t> fSignalType(_geom->Nplanes(),larlight::GEO::kMysteryType);
    std::vector<larlight::GEO::View_t> fViewType(_geom->Nplanes(),larlight::GEO::kUnknown);
    std::vector<Double_t> fPlanePitch(_geom->Nplanes(),-1.);

    for(size_t i=0; i<_geom->Nplanes(); ++i) {
      fSignalType[i] = (larlight::GEO::SigType_t)(_geom->Plane(i).SignalType());
      fViewType[i]   = (larlight::GEO::View_t)(_geom->Plane(i).View());
      fPlanePitch[i] = _geom->Plane(i).WirePitch();
      std::vector<UShort_t> wire_to_channel(_geom->Plane(i).Nwires(),larlight::DATA::INVALID_USHORT);
      for(size_t j=0; j<_geom->Plane(i).Nwires(); ++j)
	wire_to_channel[j]=_geom->PlaneWireToChannel(i,j);
      fPlaneWireToChannelMap[i]=wire_to_channel;
    }
  
    // Vectors with length = view
    // Find the maximum view type value
    std::set<geo::View_t> views = _geom->Views();
    size_t view_max = (*(views.rbegin()));
    std::vector<Double_t> fWirePitch(view_max+1,larlight::DATA::INVALID_DOUBLE);
    std::vector<Double_t> fWireAngle(view_max+1,larlight::DATA::INVALID_DOUBLE);
    for(auto iter=views.begin(); iter!=views.end(); ++iter) {
      fWirePitch[(size_t)((*iter))]=_geom->WirePitch((*iter));
      fWireAngle[(size_t)((*iter))]=_geom->WireAngleToVertical((*iter));
    }

    Double_t xyz[3]={0.};
    std::vector<std::vector<Float_t> > fOpChannelVtx(_geom->Cryostat().NOpDet(),std::vector<Float_t>(3,-1.));
    for(size_t i=0; i<_geom->Cryostat().NOpDet(); ++i) {

      _geom->Cryostat().OpDet(i).GetCenter(xyz);
      fOpChannelVtx[i][0]=xyz[0];
      fOpChannelVtx[i][1]=xyz[1];
      fOpChannelVtx[i][2]=xyz[2];
    }

    std::vector<std::vector<Double_t> > fPlaneOriginVtx(_geom->Nplanes(),std::vector<Double_t>(3,larlight::DATA::INVALID_DOUBLE));
    const Double_t orig[3]={0.};
    for(size_t i=0; i<_geom->Nplanes(); ++i) {
      _geom->Plane(i).LocalToWorld(orig,xyz);
      fPlaneOriginVtx[i][0] = xyz[0];
      fPlaneOriginVtx[i][1] = xyz[1];
      fPlaneOriginVtx[i][2] = xyz[2];
    }

    _geom_tree->Branch("fDetLength",&fDetLength,"fDetLength/D");
    _geom_tree->Branch("fDetHalfWidth",&fDetHalfWidth,"fDetHalfWidth/D");
    _geom_tree->Branch("fDetHalfHeight",&fDetHalfHeight,"fDetHalfHeight/D");

    _geom_tree->Branch("fCryoLength",&fCryoLength,"fCryoLength/D");
    _geom_tree->Branch("fCryoHalfWidth",&fCryoHalfWidth,"fCryoHalfWidth/D");
    _geom_tree->Branch("fCryoHalfHeight",&fCryoHalfHeight,"fCryoHalfHeight/D");

    _geom_tree->Branch("fChannelToPlaneMap","std::vector<UChar_t>",&fChannelToPlaneMap);
    _geom_tree->Branch("fChannelToWireMap","std::vector<UShort_t>",&fChannelToWireMap);
    _geom_tree->Branch("fPlaneWireToChannelMap","std::vector<std::vector<UShort_t> >",&fPlaneWireToChannelMap);
  
    _geom_tree->Branch("fSignalType","std::vector<larlight::GEO::SigType_t>",&fSignalType);
    _geom_tree->Branch("fViewType","std::vector<larlight::GEO::View_t>",&fViewType);
    _geom_tree->Branch("fPlanePitch","std::vector<Double_t>",&fPlanePitch);

    _geom_tree->Branch("fWireStartVtx","std::vector<std::vector<std::vector<Double_t> > >",&fWireStartVtx);
    _geom_tree->Branch("fWireEndVtx","std::vector<std::vector<std::vector<Double_t> > >",&fWireEndVtx);
    _geom_tree->Branch("fWirePitch","std::vector<Double_t>",&fWirePitch);
    _geom_tree->Branch("fWireAngle","std::vector<Double_t>",&fWireAngle);

    _geom_tree->Branch("fOpChannelVtx","std::vector<std::vector<Float_t> >",&fOpChannelVtx);

    _geom_tree->Branch("fPlaneOriginVtx","std::vector<std::vector<Double_t> >",&fPlaneOriginVtx);

    _geom_tree->Fill();
  }

  void UtilScanner::SaveDetectorProperties(fhicl::ParameterSet const& pset)
  {
    if(_detp_tree) return;
    art::ServiceHandle<art::TFileService>  fileService;    
    art::ServiceHandle<util::DetectorProperties> _detp;
    art::ServiceHandle<geo::Geometry> _geom;
    TTree* _detp_tree = fileService->make<TTree>("DetectorProperties","");

    //--- Fill Variables ---//
    Double_t fSamplingRate = _detp-> SamplingRate();          ///< in ns
    Int_t    fTriggerOffset = _detp->TriggerOffset();         ///< in # of clock ticks
    Double_t fElectronsToADC = _detp->ElectronsToADC();       ///< conversion factor for # of ionization electrons to 1 ADC count
    UInt_t   fNumberTimeSamples = _detp->NumberTimeSamples(); ///< number of clock ticks per event
    UInt_t   fReadOutWindowSize = _detp->ReadOutWindowSize(); ///< number of clock ticks per readout window
    Double_t fTimeOffsetU = _detp->TimeOffsetU();             ///< time offsets to convert spacepoint
    Double_t fTimeOffsetV = _detp->TimeOffsetV();             ///< coordinates to hit times on each
    Double_t fTimeOffsetZ = _detp->TimeOffsetZ();             ///< view
    Double_t fXTicksCoefficient = _detp->GetXTicksCoefficient(); ///< Parameters for x<-->ticks
    std::vector<Double_t> fXTicksOffsets(_geom->Nplanes(),3);
    for(unsigned int i=0; i<fXTicksOffsets.size(); ++i)
      fXTicksOffsets[i] = _detp->GetXTicksOffset(i,0,0);

    //--- Set TTree Branches ---//
    _detp_tree->Branch("fSamplingRate",&fSamplingRate,"fSamplingRate/D");
    _detp_tree->Branch("fTriggerOffset",&fTriggerOffset,"fTriggerOffset/I");
    _detp_tree->Branch("fElectronsToADC",&fElectronsToADC,"fElectronsToADC/D");
    _detp_tree->Branch("fNumberTimeSamples",&fNumberTimeSamples,"fNumberTimeSamples/i");
    _detp_tree->Branch("fReadOutWindowSize",&fReadOutWindowSize,"fReadOutWindowSize/i");
    _detp_tree->Branch("fTimeOffsetU",&fTimeOffsetU,"fTimeOffsetU/D");
    _detp_tree->Branch("fTimeOffsetV",&fTimeOffsetV,"fTimeOffsetV/D");
    _detp_tree->Branch("fTimeOffsetZ",&fTimeOffsetZ,"fTimeOffsetZ/D");
    _detp_tree->Branch("fXTicksCoefficient",&fXTicksCoefficient,"fXTicksCoefficient/D");
    _detp_tree->Branch("fXTicksOffsets","std::vector<Double_t>",&fXTicksOffsets);
  
    _detp_tree->Fill();
    return;
  }

  void UtilScanner::SaveLArProperties(fhicl::ParameterSet const& pset)
  {
    if(_larp_tree) return;
    art::ServiceHandle<util::LArProperties> _larp;
    art::ServiceHandle<geo::Geometry> _geom;

    //--- Fill Variables ---//
    std::vector< Double_t >          fEfield(_geom->Nplanes(),0);
    for(size_t i=0; i<fEfield.size(); ++i) { fEfield[i]=_larp->Efield(i);}
    Double_t                         fTemperature = _larp->Temperature();
    Double_t                         fElectronlifetime = _larp->ElectronLifetime(); ///< microseconds
    Double_t                         fRadiationLength = _larp->RadiationLength();  ///< g/cm^2

    Double_t                         fArgon39DecayRate = _larp->Argon39DecayRate(); ///<  decays per cm^3 per second
  
    // Following parameters are for use in Bethe-Bloch formula for dE/dx.
    Double_t fZ = pset.get<double>("AtomicNumber");       ///< Ar atomic number
    Double_t fA = pset.get<double>("AtomicMass");         ///< Ar atomic mass (g/mol)
    Double_t fI = pset.get<double>("ExcitationEnergy");   ///< Ar mean excitation energy (eV)
    Double_t fSa= pset.get<double>("SternheimerA");       ///< Sternheimer parameter a
    Double_t fSk= pset.get<double>("SternheimerK");       ///< Sternheimer parameter k
    Double_t fSx0 = pset.get<double>("SternheimerX0");    ///< Sternheimer parameter x0
    Double_t fSx1 = pset.get<double>("SternheimerX1");    ///< Sternheimer parameter x1
    Double_t fScbar = pset.get<double>("SternheimerCbar");///< Sternheimer parameter Cbar
  
    // Optical parameters for Dar 
    std::vector<Double_t> fFastScintEnergies = pset.get< std::vector<double> >("FastScintEnergies");
    std::vector<Double_t> fFastScintSpectrum = pset.get< std::vector<double> >("FastScintSpectrum");
    std::vector<Double_t> fSlowScintEnergies = pset.get< std::vector<double> >("SlowScintEnergies");
    std::vector<Double_t> fSlowScintSpectrum = pset.get< std::vector<double> >("SlowScintSpectrum");
    std::vector<Double_t> fAbsLengthEnergies = pset.get< std::vector<double> >("AbsLengthEnergies");
    std::vector<Double_t> fAbsLengthSpectrum = pset.get< std::vector<double> >("AbsLengthSpectrum");
    std::vector<Double_t> fRIndexEnergies    = pset.get< std::vector<double> >("RIndexEnergies"   );
    std::vector<Double_t> fRIndexSpectrum    = pset.get< std::vector<double> >("RIndexSpectrum"   );
    std::vector<Double_t> fRayleighEnergies  = pset.get< std::vector<double> >("RayleighEnergies" );
    std::vector<Double_t> fRayleighSpectrum  = pset.get< std::vector<double> >("RayleighSpectrum" );
  
    bool fScintByParticleType = pset.get<bool>("ScintByParticleType");

    Double_t fProtonScintYield        = pset.get<double>("ProtonScintYield"     );
    Double_t fProtonScintYieldRatio   = pset.get<double>("ProtonScintYieldRatio");
    Double_t fMuonScintYield          = pset.get<double>("MuonScintYield"       );
    Double_t fMuonScintYieldRatio     = pset.get<double>("MuonScintYieldRatio"  );
    Double_t fPionScintYield          = pset.get<double>("PionScintYield"       );
    Double_t fPionScintYieldRatio     = pset.get<double>("PionScintYieldRatio"  );
    Double_t fKaonScintYield          = pset.get<double>("KaonScintYield"       );
    Double_t fKaonScintYieldRatio     = pset.get<double>("KaonScintYieldRatio"  );
    Double_t fElectronScintYield      = pset.get<double>("ElectronScintYield"   );
    Double_t fElectronScintYieldRatio = pset.get<double>("ElectronScintYieldRatio");
    Double_t fAlphaScintYield         = pset.get<double>("AlphaScintYield"      );
    Double_t fAlphaScintYieldRatio    = pset.get<double>("AlphaScintYieldRatio" );  

    Double_t fScintResolutionScale = pset.get<double>("ScintResolutionScale");
    Double_t fScintFastTimeConst   = pset.get<double>("ScintFastTimeConst"  );
    Double_t fScintSlowTimeConst   = pset.get<double>("ScintSlowTimeConst"  );
    Double_t fScintBirksConstant   = pset.get<double>("ScintBirksConstant"  );
    Double_t fScintYield           = pset.get<double>("ScintYield"          );
    Double_t fScintYieldRatio      = pset.get<double>("ScintYieldRatio"     );
  
    bool fEnableCerenkovLight = pset.get<bool>("EnableCerenkovLight");

    std::vector<std::string> fReflectiveSurfaceNames = pset.get<std::vector<std::string> >("ReflectiveSurfaceNames");
    std::vector<Double_t> fReflectiveSurfaceEnergies = pset.get<std::vector<double> >("ReflectiveSurfaceEnergies");;
    std::vector<std::vector<Double_t> > fReflectiveSurfaceReflectances = pset.get<std::vector<std::vector<double> > >("ReflectiveSurfaceReflectances");
    std::vector<std::vector<Double_t> > fReflectiveSurfaceDiffuseFractions = pset.get<std::vector<std::vector<double> > >("ReflectiveSurfaceDiffuseFractions");

    //--- Set TTree Branches ---//
    art::ServiceHandle<art::TFileService>  fileService;    
    TTree* _larp_tree = fileService->make<TTree>("LArProperties","");

    _larp_tree->Branch("fEfield","std::vector<Double_t>", &fEfield);
    _larp_tree->Branch("fTemperature",&fTemperature,"fTemperature/D");
    _larp_tree->Branch("fElectronlifetime",&fElectronlifetime,"fElectronlifetime/D");
    _larp_tree->Branch("fRadiationLength",&fRadiationLength,"fRadiationLength/D");

    _larp_tree->Branch("fArgon39DecayRate",&fArgon39DecayRate,"fArgon39DecayRate/D");

    _larp_tree->Branch("fZ",&fZ,"fZ/D");
    _larp_tree->Branch("fA",&fA,"fA/D");
    _larp_tree->Branch("fI",&fI,"fI/D");
    _larp_tree->Branch("fSa",&fSa,"fSa/D");
    _larp_tree->Branch("fSk",&fSk,"fSk/D");
    _larp_tree->Branch("fSx0",&fSx0,"fSx0/D");
    _larp_tree->Branch("fSx1",&fSx1,"fSx1/D");
    _larp_tree->Branch("fScbar",&fScbar,"fScbar/D");

    _larp_tree->Branch("fFastScintSpectrum","std::vector<Double_t>",&fFastScintSpectrum);
    _larp_tree->Branch("fFastScintEnergies","std::vector<Double_t>",&fFastScintEnergies);
    _larp_tree->Branch("fSlowScintSpectrum","std::vector<Double_t>",&fSlowScintSpectrum);
    _larp_tree->Branch("fSlowScintEnergies","std::vector<Double_t>",&fSlowScintEnergies);
    _larp_tree->Branch("fRIndexSpectrum","std::vector<Double_t>",&fRIndexSpectrum);
    _larp_tree->Branch("fRIndexEnergies","std::vector<Double_t>",&fRIndexEnergies);
    _larp_tree->Branch("fAbsLengthSpectrum","std::vector<Double_t>",&fAbsLengthSpectrum);
    _larp_tree->Branch("fAbsLengthEnergies","std::vector<Double_t>",&fAbsLengthEnergies);
    _larp_tree->Branch("fRayleighSpectrum","std::vector<Double_t>",&fRayleighSpectrum);
    _larp_tree->Branch("fRayleighEnergies","std::vector<Double_t>",&fRayleighEnergies);

    _larp_tree->Branch("fScintByParticleType",&fScintByParticleType,"fScintByParticleType/O");

    _larp_tree->Branch("fProtonScintYield",&fProtonScintYield,"fProtonScintYield/D");
    _larp_tree->Branch("fProtonScintYieldRatio",&fProtonScintYieldRatio,"fProtonScintYieldRatio/D");
    _larp_tree->Branch("fMuonScintYield",&fMuonScintYield,"fMuonScintYield/D");
    _larp_tree->Branch("fMuonScintYieldRatio",&fMuonScintYieldRatio,"fMuonScintYieldRatio/D");
    _larp_tree->Branch("fPionScintYield",&fPionScintYield,"fPionScintYield/D");
    _larp_tree->Branch("fPionScintYieldRatio",&fPionScintYieldRatio,"fPionScintYieldRatio/D");
    _larp_tree->Branch("fKaonScintYield",&fKaonScintYield,"fKaonScintYield/D");
    _larp_tree->Branch("fKaonScintYieldRatio",&fKaonScintYieldRatio,"fKaonScintYieldRatio/D");
    _larp_tree->Branch("fElectronScintYield",&fElectronScintYield,"fElectronScintYield/D");
    _larp_tree->Branch("fElectronScintYieldRatio",&fElectronScintYieldRatio,"fElectronScintYieldRatio/D");
    _larp_tree->Branch("fAlphaScintYield",&fAlphaScintYield,"fAlphaScintYield/D");
    _larp_tree->Branch("fAlphaScintYieldRatio",&fAlphaScintYieldRatio,"fAlphaScintYieldRatio/D");

    _larp_tree->Branch("fScintYield",&fScintYield,"fScintYield/D");
    _larp_tree->Branch("fScintResolutionScale",&fScintResolutionScale,"fScintResolutionScale/D");
    _larp_tree->Branch("fScintFastTimeConst",&fScintFastTimeConst,"fScintFastTimeConst/D");
    _larp_tree->Branch("fScintSlowTimeConst",&fScintSlowTimeConst,"fScintSlowTimeConst/D");
    _larp_tree->Branch("fScintYieldRatio",&fScintYieldRatio,"fScintYieldRatio/D");  
    _larp_tree->Branch("fScintBirksConstant",&fScintBirksConstant,"fScintBirksConstant/D");

    _larp_tree->Branch("fEnableCerenkovLight",&fEnableCerenkovLight,"fEnableCerenkovLight/O");  

    _larp_tree->Branch("fReflectiveSurfaceNames","std::vector<std::string>",&fReflectiveSurfaceNames);
    _larp_tree->Branch("fReflectiveSurfaceEnergies","std::vector<Double_t>",&fReflectiveSurfaceEnergies);
    _larp_tree->Branch("fReflectiveSurfaceReflectances","std::vector<std::vector<Double_t> >",&fReflectiveSurfaceReflectances);
    _larp_tree->Branch("fReflectiveSurfaceDiffuseFractions","std::vector<std::vector<Double_t> >",&fReflectiveSurfaceDiffuseFractions);

    _larp_tree->Fill();

  }



} // namespace opdet



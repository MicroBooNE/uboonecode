////////////////////////////////////////////////////////////////////////
// Class:       FlashTrigger
// Module Type: filter
// File:        FlashTrigger_module.cc
//
// Generated at Tue Sep  8 12:34:50 2015 by Kazuhiro Terao using artmod
// from cetpkgsupport v1_08_06.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"

#include <memory>

#include "RecoBase/OpFlash.h"
#include "Geometry/Geometry.h"

#include <TH1D.h>
#include <string>
#include <iostream>
#include <vector>
class FlashTrigger;

class FlashTrigger : public art::EDFilter {
public:
  explicit FlashTrigger(fhicl::ParameterSet const & p);
  // The destructor generated by the compiler is fine for classes
  // without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  FlashTrigger(FlashTrigger const &) = delete;
  FlashTrigger(FlashTrigger &&) = delete;
  FlashTrigger & operator = (FlashTrigger const &) = delete;
  FlashTrigger & operator = (FlashTrigger &&) = delete;

  // Required functions.
  bool filter(art::Event & e) override;


private:

  // Declare member data here.
  TH1D *_hFlashEff;                     ///< Flash PE efficiency analysis histogram
  double _event_ctr;                    ///< A counter for the total number of events processed
  std::vector<double> _flash_ctr_v;     ///< A total number of flash above specific pe threshold set in _npe_threshold_v
  std::vector<double> _npe_threshold_v; ///< An array of npe threshold to fill a binned histogram (_hFlashEff)
  double _beam_time_diff_low;           ///< A lower boundary of time window with which we ask for a time coincident with beam gate
  double _beam_time_diff_high;          ///< A higher boundary of time window with which we ask for a time coincident with beam gate
  double _total_pe_threshold;           ///< The minimum number of photo-electrons threshold to pass an event
  unsigned int _multiplicity_threshold; ///< The minimum number of optical detectors considered to be "hit" to pass an event
  double _multiplicity_pe_threshold;    ///< A threshold value above which PMT is considered as "hit" for multiplicity condition
  std::string _flash_module;            ///< recob::OpFlash producer module label
  bool _verbose;                        ///< Verbosity flag
};


FlashTrigger::FlashTrigger(fhicl::ParameterSet const & p)
  : _hFlashEff (nullptr)
  , _event_ctr (0)
// Initialize member data here.
{
  // Call appropriate produces<>() functions here.
  _verbose                   = p.get< bool         > ( "Verbose"                 );
  _flash_module              = p.get< std::string  > ( "OpFlashModule"           );
  _beam_time_diff_low        = p.get< double       > ( "BeamDTLow"               );
  _beam_time_diff_high       = p.get< double       > ( "BeamDTHigh"              );
  _total_pe_threshold        = p.get< double       > ( "TotalPEThreshold"        );
  _multiplicity_pe_threshold = p.get< double       > ( "MultiplicityPEThreshold" );
  _multiplicity_threshold    = p.get< unsigned int > ( "MultiplicityThreshold"   );
  _npe_threshold_v           = p.get< std::vector<double> > ( "EffHistPEThresholdArray" );

  if(_flash_module.empty()) {
    std::cerr << "OpFlashModule is empty!" << std::endl;
    throw std::exception();
  }

  _flash_ctr_v.resize(_npe_threshold_v.size(),0);

  if(_verbose) {
    std::cout << "\033[93m" << __PRETTY_FUNCTION__ << "\033[00m" << std::endl
	      << "  recob::OpFlash module label: " << _flash_module.c_str() << std::endl
	      << "  DT window w.r.t. beam spill: " << _beam_time_diff_low << " => " << _beam_time_diff_high << std::endl
	      << "  \033[95mPass condition...\033[00m " << std::endl
	      << "  0) Flash within DT window with above " << _total_pe_threshold << " p.e." << std::endl
	      << "  1) Flash within DT window with multiplicity of " << _multiplicity_threshold << " PMT hits with > " 
	      << _multiplicity_pe_threshold << " p.e." << std::endl
	      << std::endl;
  }
  
}

bool FlashTrigger::filter(art::Event & e)
{
  _event_ctr += 1; // Increment total event ctr

  bool pass = false;

  // Retrieve OpFlash
  art::Handle< std::vector< recob::OpFlash > > flash_handle;
  e.getByLabel( _flash_module, flash_handle );

  // Geometry service
  art::ServiceHandle<geo::Geometry> geo;

  // If valid, perform
  if(flash_handle.isValid()) {
    for(auto const& flash : *flash_handle) {
      // Check timing
      if(flash.Time() < _beam_time_diff_low || _beam_time_diff_high < flash.Time()) {
	if(_verbose)
	  std::cout << "  Skipping a flash @ DT = " << flash.Time() << " [us] " << std::endl;
	continue;
      }

      // Fill analysis histogram
      auto const npe = flash.TotalPE();

      for(size_t thres_index=0; thres_index < _npe_threshold_v.size(); ++thres_index) {
	auto const& thres = _npe_threshold_v[thres_index];
	if(npe >= thres) _flash_ctr_v[thres_index] +=1;	
      }

      // Check Total PE level
      if( npe >= _total_pe_threshold )
	pass = true;
      
      // Check multiplicity
      unsigned int mult=0;
      for(size_t i=0; i < geo->NOpDets(); ++i )
	if(flash.PE(i) >= _multiplicity_pe_threshold) mult+=1;
      
      if(mult >= _multiplicity_threshold)
	pass = true;

      if(pass) {
	if(_verbose)
	  std::cout << "  Accepting a flash @ " << flash.Time() << " [us] with " << npe << " [p.e.]" << " ... multiplicity = " << mult << std::endl;
	break;
      }
    }
  }

  if(_flash_ctr_v.size()) {
    if(!_hFlashEff) {
      art::ServiceHandle<art::TFileService> tfs;
      _hFlashEff = tfs->make<TH1D>("hFlashEff","Flash Finding Efficiency; P.E. cut values; Efficiency",
				   _flash_ctr_v.size(), -0.5, _flash_ctr_v.size()-0.5);
    }
    for(size_t index=0; index < _flash_ctr_v.size(); ++index) {
      
      if(_flash_ctr_v[index] < 1) _hFlashEff->SetBinContent( index+1, 0);
      else _hFlashEff->SetBinContent( index+1, (double)(_flash_ctr_v[index]) / (double)(_event_ctr) );
      
    }
  }

  return pass;
}

DEFINE_ART_MODULE(FlashTrigger)
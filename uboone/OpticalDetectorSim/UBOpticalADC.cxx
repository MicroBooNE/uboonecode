#ifndef UBOPTICALADC_CXX
#define UBOPTICALADC_CXX

#include "UBOpticalADC.h"

namespace opdet {

  //----------------------------------------
  UBOpticalADC::UBOpticalADC() : UBADCBase()
  //----------------------------------------
  {
    Reset();
  }

  //------------------------
  void UBOpticalADC::Reset()
  //------------------------
  {
    UBADCBase::Reset();
    fSPE.Reset();
    fPED.Reset();
    fInputPhotonTime.clear();
    fDarkPhotonTime.clear();
    fPhotonTime.clear();
  }

  //--------------------------------------------------------------
  void UBOpticalADC::SetPhotons(const std::vector<double>& g4time)
  //--------------------------------------------------------------
  {
    fInputPhotonTime.clear();
    fInputPhotonTime.reserve(g4time.size());
    for(auto const &v : g4time) fInputPhotonTime.push_back(v);
  }

  //--------------------------------------------------------------------------
  //void UBOpticalADC::GenDarkNoise(double dark_rate, double period)
  void UBOpticalADC::GenDarkNoise(const unsigned int pmtid, const double g4start)
  //--------------------------------------------------------------------------
  {
    fDarkPhotonTime.clear();

    art::ServiceHandle<opdet::UBOpticalChConfig> ch_conf;
    art::ServiceHandle<geo::Geometry> geom;
    unsigned int ch = geom->OpChannel( pmtid, 0 ); // get channel reading out that PMT

    double dark_rate = ch_conf->GetParameter(kDarkRate,ch);

    unsigned int dark_count = RandomServer::GetME().Poisson(dark_rate * fDuration);

    fDarkPhotonTime.reserve(dark_count);

    for(size_t i=0; i<dark_count; ++i)

      fDarkPhotonTime.push_back(RandomServer::GetME().Uniform(fDuration*1.e3) + g4start);

  }

  //-------------------------------------------------------------------
  void UBOpticalADC::GenWaveform(const unsigned int ch, optdata::ChannelData& wf )
  //-------------------------------------------------------------------
  {
    //
    // Initialize, zero
    //
    size_t nticks = fTimeInfo.Ticks(fDuration);
    std::vector< float > wfm_tmp( nticks, 0.0 );
    if ( wf.size()!=nticks )
      wf.resize( nticks, 0 );
    else
      wf.assign( nticks, 0 );

    // services
    art::ServiceHandle<opdet::UBOpticalChConfig> ch_conf;

    //
    // Generate Signal & DarkNoise
    //
    // Configure to generate high gain SPE
    fSPE.Reset();

    fSPE.SetT0(ch_conf->GetParameter(kT0,ch),
	       ch_conf->GetParameter(kT0Spread,ch));
    
    fSPE.SetGain(ch_conf->GetParameter(kPMTGain,ch),
		 ch_conf->GetParameter(kGainSpread,ch));
    
    // Create combined photon time with QE applied on signal photons

    fPhotonTime.clear();
    fPhotonTime.reserve(fInputPhotonTime.size() + fDarkPhotonTime.size());
    const double qe = ch_conf->GetParameter(kQE,ch);
    for(auto const &v : fDarkPhotonTime) fPhotonTime.push_back(v);
    for(auto const &v : fInputPhotonTime)
      if(RandomServer::GetME().Uniform(1.) < qe) fPhotonTime.push_back(v);

    fSPE.SetPhotons(fPhotonTime);
    fSPE.Process(wfm_tmp,fTimeInfo);

    // convert from pe waveform to adc
    double gain_ratio = ch_conf->GetParameter(kSplitterGain,ch);
    for(auto &v : wfm_tmp) 
      v *= gain_ratio;

    //
    // Simulate pedestal
    //
    fPED.Reset();
    fPED.SetPedestal(ch_conf->GetParameter(kPedestalMean,ch),
		     ch_conf->GetParameter(kPedestalSpread,ch));
    fPED.Process(wfm_tmp,fTimeInfo);

    // Make sure algorithms did not alter the waveform size

    if(wf.size()!=nticks || wfm_tmp.size()!=nticks)
      throw UBOpticalException("Waveform length changed (prohibited)!");

    //
    // Digitize amplitude
    //
    wf.clear();
    Digitize(wfm_tmp,wf);
    
  }
  
}

#endif

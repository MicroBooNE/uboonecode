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
  void UBOpticalADC::GenDarkNoise(const unsigned int ch, const double g4start)
  //--------------------------------------------------------------------------
  {
    fDarkPhotonTime.clear();

    art::ServiceHandle<opdet::UBOpticalChConfig> ch_conf;

    double dark_rate = ch_conf->GetParameter(kDarkRate,ch);

    unsigned int dark_count = RandomServer::GetME().Poisson(dark_rate * fDuration);

    fDarkPhotonTime.reserve(dark_count);

    for(size_t i=0; i<dark_count; ++i)

      fDarkPhotonTime.push_back(RandomServer::GetME().Uniform(fDuration*1.e3) + g4start);

  }

  //-------------------------------------------------------------------
  void UBOpticalADC::GenWaveform(const unsigned int ch, 
				 std::vector<unsigned short>& high_wf,
				 std::vector<unsigned short>& low_wf  )
  //-------------------------------------------------------------------
  {
    //
    // Initialize
    //
    size_t nticks = fTimeInfo.Ticks(fDuration);
    std::vector<float> high_tmp_wf(nticks,0);
    std::vector<float> low_tmp_wf;
    low_tmp_wf.reserve(nticks);

    art::ServiceHandle<opdet::UBOpticalChConfig> ch_conf;
    //
    // Generate Signal & DarkNoise
    //
    // Configure to generate high gain SPE
    fSPE.Reset();

    fSPE.SetT0(ch_conf->GetParameter(kT0,ch),
	       ch_conf->GetParameter(kT0Spread,ch));

    fSPE.SetGain(ch_conf->GetParameter(kHighGain,ch),
		 ch_conf->GetParameter(kHighGain,ch) *
		 ch_conf->GetParameter(kGainSpread,ch));
    
    // Create combined photon time with QE applied on signal photons

    fPhotonTime.clear();
    fPhotonTime.reserve(fInputPhotonTime.size() + fDarkPhotonTime.size());
    const double qe = ch_conf->GetParameter(kQE,ch);
    for(auto const &v : fDarkPhotonTime) fPhotonTime.push_back(v);
    for(auto const &v : fInputPhotonTime)

      if(RandomServer::GetME().Uniform(1.) < qe) fPhotonTime.push_back(v);

    fSPE.SetPhotons(fPhotonTime);
    fSPE.Process(high_tmp_wf,fTimeInfo);

    // Copy waveform to low_tmp_wf with gain ratio
    double gain_ratio = (ch_conf->GetParameter(kLowGain,ch) / 
			 ch_conf->GetParameter(kHighGain,ch));
    for(auto const &v : high_tmp_wf) 

      low_tmp_wf.push_back(v * gain_ratio);

    //
    // Simulate pedestal
    //
    fPED.Reset();
    fPED.SetPedestal(ch_conf->GetParameter(kPedestalMean,ch),
		     ch_conf->GetParameter(kPedestalSpread,ch));
    fPED.Process(high_tmp_wf,fTimeInfo);
    fPED.Process(low_tmp_wf,fTimeInfo);

    // Make sure algorithms did not alter the waveform size

    if(high_tmp_wf.size()!=nticks || low_tmp_wf.size()!=nticks)

      throw UBOpticalException("Waveform length changed (prohibited)!");

    //
    // Digitize amplitude
    //
    high_wf.clear();
    low_wf.clear();
    Digitize(high_tmp_wf,high_wf);
    Digitize(low_tmp_wf,low_wf);
    
  }
  
}

#endif

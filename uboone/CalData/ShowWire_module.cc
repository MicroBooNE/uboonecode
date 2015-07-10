////////////////////////////////////////////////////////////////////////
// Class:       ShowWire
// Module Type: analyzer
// File:        ShowWire_module.cc
//
// Generated at Thu Jul 31 15:07:10 2014 by Wesley Ketchum using artmod
// from cetpkgsupport v1_05_04.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/FindOneP.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Registry/ServiceHandle.h" 
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "CalibrationDBI/Interface/IDetPedestalService.h"
#include "CalibrationDBI/Interface/IDetPedestalProvider.h"
#include "SimpleTypesAndConstants/geo_types.h" // geo::SigType_t
#include "SimpleTypesAndConstants/RawTypes.h" // raw::ChannelID_t
#include "Geometry/Geometry.h"
#include "RecoBase/Wire.h"
#include "RawData/RawDigit.h"

#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include <string>
#include <sstream>
#include <iostream>


namespace cal{ class ShowWire; }

class cal::ShowWire : public art::EDAnalyzer {
public:
  explicit ShowWire(fhicl::ParameterSet const & p);
  virtual ~ShowWire();

  void beginJob();
  void endJob();

  void analyze(art::Event const & e) override;
  void reconfigure(fhicl::ParameterSet const& p);

private:

  // Declare member data here.
  std::string fCalDataModuleLabel;
  unsigned int fPeakFitWindow;
  bool fDoSinglePulseChecks;
  unsigned int fSinglePulseLocation;
  unsigned int fSinglePulseBuffer;
  raw::ChannelID_t fChannel;
  bool fSaveWaveforms;

  std::string MakeHistName(const char*, const unsigned int, const unsigned int, const unsigned int, const unsigned int);
  std::string MakeHistTitle(const char*, const unsigned int, const unsigned int, const unsigned int, const unsigned int);
  void SetHistogram(TH1F*, const char*, const unsigned int, const unsigned int, const unsigned int, const unsigned int, const size_t, const Color_t);
  void FillWaveforms(recob::Wire const&, raw::RawDigit const&, TH1F*, TH1F*);

  float FindPeakTime(TH1F*);
  float FindPeak(TH1F*);
  float FindIntegral(TH1F*);

  TH1F* h_integral_wire;
  TH1F* h_integral_raw;
  TH2F* h_integral_2D;

  TH1F* h_peak_wire;
  TH1F* h_peak_raw;
  TH1F* h_peak_diff;
  TH2F* h_peak_2D;

  TH1F* h_peaktime_wire;
  TH1F* h_peaktime_raw;
  TH1F* h_peaktime_diff;
  TH2F* h_peaktime_2D;

  TH1F* h_wire;
  TH1F* h_raw;

};


cal::ShowWire::ShowWire(fhicl::ParameterSet const & p) 
  :
  EDAnalyzer(p)
 // More initializers here.
{
  this->reconfigure(p);
}

cal::ShowWire::~ShowWire()
{
  // Clean up dynamic memory and other resources here.
}
void cal::ShowWire::reconfigure(fhicl::ParameterSet const& p){
  fCalDataModuleLabel  = p.get<std::string>("CalDataModuleLabel");
  fPeakFitWindow       = p.get<unsigned int>("PeakFitWindow",3);
  fDoSinglePulseChecks = p.get<bool>("DoSinglePulseChecks",true);
  fSinglePulseLocation = p.get<unsigned int>("SinglePulseLocation",5000);
  fSinglePulseBuffer   = p.get<unsigned int>("SinglePulseBuffer",25);
  fChannel             = p.get<unsigned int>("Channel",7775);
  fSaveWaveforms       = p.get<bool>("SaveWaveforms",true);
  
}

void cal::ShowWire::beginJob(){

  if(!fDoSinglePulseChecks) return;

  art::ServiceHandle<art::TFileService> tfs;

  h_integral_wire = tfs->make<TH1F>("h_integral_wire","Integral of wire pulses",100,-1000,4000);
  h_integral_raw = tfs->make<TH1F>("h_integral_raw","Integral of rawdigit pulses",100,-1000,4000);
  h_integral_2D = tfs->make<TH2F>("h_integral_2D","Integral of wire and rawdigit pulses",100,-1000,4000,100,-1000,4000);
  
  h_peak_wire = tfs->make<TH1F>("h_peak_wire","Peak of wire pulses",100,-100,900);
  h_peak_raw = tfs->make<TH1F>("h_peak_raw","Peak of rawdigit pulses",100,-100,900);
  h_peak_diff = tfs->make<TH1F>("h_peak_diff","Difference in peak of wire and rawdigit pulses;Wire Peak Amplitude - Rawdigit Peak Amplitude",100,-10,10);
  h_peak_2D = tfs->make<TH2F>("h_peak_2D","Peak of wire and rawdigit pulses",100,-100,900,100,-100,900);
  
  h_peaktime_wire = tfs->make<TH1F>("h_peaktime_wire","Peak time of wire pulses",
				    100,fSinglePulseLocation-fSinglePulseBuffer,fSinglePulseLocation+fSinglePulseBuffer);
  h_peaktime_raw = tfs->make<TH1F>("h_peaktime_raw","Peak time of rawdigit pulses",
				   100,fSinglePulseLocation-fSinglePulseBuffer,fSinglePulseLocation+fSinglePulseBuffer);
  h_peaktime_diff = tfs->make<TH1F>("h_peaktime_diff","Difference in peak times of wire and rawdigit pulses;Wire Peak Time - Rawdigit Peak Time",200,-10,10);
  h_peaktime_2D = tfs->make<TH2F>("h_peaktime_2D","Peak time of wire and rawdigit pulses",
				  100,fSinglePulseLocation-fSinglePulseBuffer,fSinglePulseLocation+fSinglePulseBuffer,
				  100,fSinglePulseLocation-fSinglePulseBuffer,fSinglePulseLocation+fSinglePulseBuffer);
  

}

void cal::ShowWire::endJob(){}

float cal::ShowWire::FindPeakTime(TH1F* hist){
  //hist->GetXaxis()->SetRange((int)(fSinglePulseLocation-fSinglePulseBuffer),(int)(fSinglePulseLocation+fSinglePulseBuffer));

  int max_bin = hist->GetMaximumBin();
  float min_fit = hist->GetBinLowEdge(max_bin-fPeakFitWindow);
  float max_fit = hist->GetBinLowEdge(max_bin+fPeakFitWindow+1);
  hist->Fit("gaus","Q","",min_fit,max_fit);

  TF1 *thisfit = (TF1*) hist->GetFunction("gaus");
  return thisfit->GetParameter(1);
}

float cal::ShowWire::FindPeak(TH1F* hist){
  //hist->GetXaxis()->SetRange((int)(fSinglePulseLocation-fSinglePulseBuffer),(int)(fSinglePulseLocation+fSinglePulseBuffer));
  return hist->GetMaximum();
}

float cal::ShowWire::FindIntegral(TH1F* hist){
  //return hist->Integral((int)(fSinglePulseLocation-fSinglePulseBuffer),(int)(fSinglePulseLocation+fSinglePulseBuffer));
  return hist->Integral();
}

std::string cal::ShowWire::MakeHistName(const char* label, const unsigned int run, const unsigned int subrun, const unsigned int event, const unsigned int channel){
  std::stringstream ss_hist_name;
  ss_hist_name << "h" << label << "_Run" << run << "_SubRun" << subrun << "_Event" << event << "_Channel" << channel;
  return ss_hist_name.str();
}

std::string cal::ShowWire::MakeHistTitle(const char* label, const unsigned int run, const unsigned int subrun, const unsigned int event, const unsigned int channel){
  std::stringstream ss_hist_title;
  ss_hist_title << label << ", Run" << run << "_SubRun" << subrun << "_Event" << event << "_Channel" << channel << ";Ticks;ADCs";
  return ss_hist_title.str();
}

void cal::ShowWire::SetHistogram(TH1F* hist, 
				 const char* label, 
				 const unsigned int run, const unsigned int subrun, const unsigned int event, const unsigned int channel, 
				 const size_t nbins,
				 const Color_t color=kBlack){
  std::string hist_name = MakeHistName(label,run,subrun,event,channel);
  std::string hist_title = MakeHistTitle(label,run,subrun,event,channel);
  hist->SetNameTitle(hist_name.c_str(),hist_title.c_str());

  if(fDoSinglePulseChecks)
    hist->SetBins(2*fSinglePulseBuffer,fSinglePulseLocation-fSinglePulseBuffer,fSinglePulseLocation+fSinglePulseBuffer);
  else
    hist->SetBins(nbins,0,nbins);

  hist->SetLineColor(color);
}

void cal::ShowWire::FillWaveforms(recob::Wire const& wire, raw::RawDigit const& rawdigit, TH1F* h_wire, TH1F* h_raw){
  
  //get pedestal conditions
  const lariov::IDetPedestalProvider& pedestalRetrievalAlg = art::ServiceHandle<lariov::IDetPedestalService>()->GetPedestalProvider();
  float pedestal = pedestalRetrievalAlg.PedMean(rawdigit.Channel());

  size_t begin_iter=0;
  size_t end_iter = rawdigit.Samples();
  if(fDoSinglePulseChecks){
    begin_iter = fSinglePulseLocation-fSinglePulseBuffer;
    end_iter = fSinglePulseLocation+fSinglePulseBuffer;
  }

  for(size_t i=begin_iter; i<end_iter; i++){
    h_raw->SetBinContent(i-begin_iter+1,rawdigit.ADC(i) - pedestal);
    h_wire->SetBinContent(i-begin_iter+1,wire.Signal().at(i));

    //std::cout << "bin " << i-begin_iter+1 << " bigbin " << i << ", Raw=" << rawdigit.ADC(i)-pedestal << ", Wire=" << wire.Signal().at(i) << std::endl;

  }

}

void cal::ShowWire::analyze(art::Event const & e)
{
  // Implementation of required member function here.

  art::Handle< std::vector<recob::Wire> > wireVectorHandle;
  e.getByLabel(fCalDataModuleLabel,wireVectorHandle);
  std::vector<recob::Wire> const& wireVector(*wireVectorHandle);

  unsigned int run = e.run();
  unsigned int subrun = e.subRun();
  unsigned int event = e.event();

  art::ServiceHandle<art::TFileService> tfs;
  
  // get the raw::RawDigit associated by fCalDataModuleLabel to wires in
  // wireVectorHandle; RawDigitsFromWire.at(index) will be a
  // art::Ptr<raw::RawDigit>
  art::FindOneP<raw::RawDigit> RawDigitsFromWire
    (wireVectorHandle, e, fCalDataModuleLabel);

  art::ServiceHandle<geo::Geometry> geom;
  
  size_t iWire = 0;
  for(auto const& wire : wireVector){

    raw::ChannelID_t channel = wire.Channel();
    if(channel!=fChannel && fSaveWaveforms) continue;

    raw::RawDigit const& rawdigit(*(RawDigitsFromWire.at(iWire++)));
    size_t n_samples = rawdigit.Samples();

    //this is stupid I have to do this here and can't set it later...
    if(fSaveWaveforms){
      h_wire = tfs->make<TH1F>(MakeHistName("Wire",run,subrun,event,channel).c_str(),"",n_samples,0,n_samples);
      h_raw = tfs->make<TH1F>(MakeHistName("Raw",run,subrun,event,channel).c_str(),"",n_samples,0,n_samples);
    }
    else{
      h_wire = new TH1F(MakeHistName("Wire",run,subrun,event,channel).c_str(),"",n_samples,0,n_samples);
      h_raw = new TH1F(MakeHistName("Raw",run,subrun,event,channel).c_str(),"",n_samples,0,n_samples);
    }
    SetHistogram(h_wire,"Wire",run,subrun,event,channel,n_samples,kBlue);
    SetHistogram(h_raw,"Raw",run,subrun,event,channel,n_samples,kRed);
    FillWaveforms(wire,rawdigit,h_wire,h_raw);

    if(!fDoSinglePulseChecks) {
      if(!fSaveWaveforms) { delete h_wire; delete h_raw; }
      continue;
    }
    float integral_wire = FindIntegral(h_wire);
    float integral_raw = FindIntegral(h_raw);
    h_integral_wire->Fill( integral_wire );
    h_integral_raw->Fill( integral_raw );
    h_integral_2D->Fill( integral_wire, integral_raw );

    float peak_wire = FindPeak(h_wire);
    float peak_raw = FindPeak(h_raw);
    h_peak_wire->Fill( peak_wire );
    h_peak_raw->Fill( peak_raw );
    h_peak_diff->Fill( peak_wire - peak_raw );
    h_peak_2D->Fill( peak_wire, peak_raw );

    float peaktime_wire = FindPeakTime(h_wire);
    float peaktime_raw = FindPeakTime(h_raw);
    h_peaktime_wire->Fill( peaktime_wire );
    h_peaktime_raw->Fill( peaktime_raw );
    h_peaktime_diff->Fill( peaktime_wire - peaktime_raw );
    h_peaktime_2D->Fill( peaktime_wire, peaktime_raw );

    if(!fSaveWaveforms) { delete h_wire; delete h_raw; }

  }
    
}

DEFINE_ART_MODULE(cal::ShowWire)

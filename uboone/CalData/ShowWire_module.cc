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
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"

#include "RecoBase/Wire.h"
#include "RawData/RawDigit.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "TH1F.h"

#include <string>
#include <sstream>


namespace cal{ class ShowWire; }

class cal::ShowWire : public art::EDAnalyzer {
public:
  explicit ShowWire(fhicl::ParameterSet const & p);
  virtual ~ShowWire();

  void analyze(art::Event const & e) override;
  void reconfigure(fhicl::ParameterSet const& p);

private:

  // Declare member data here.
  std::string fCalDataModuleLabel;

  std::string MakeHistName(const char*, const unsigned int, const unsigned int, const unsigned int, const unsigned int);
  std::string MakeHistTitle(const char*, const unsigned int, const unsigned int, const unsigned int, const unsigned int);
  void SetHistogram(TH1F*, const char*, const unsigned int, const unsigned int, const unsigned int, const unsigned int, const size_t, const Color_t);
  void FillWaveforms(recob::Wire const&, raw::RawDigit const&, TH1F*, TH1F*);

};


cal::ShowWire::ShowWire(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p)  // ,
 // More initializers here.
{
  this->reconfigure(p);
}

cal::ShowWire::~ShowWire()
{
  // Clean up dynamic memory and other resources here.
}
void cal::ShowWire::reconfigure(fhicl::ParameterSet const& p){
  fCalDataModuleLabel = p.get<std::string>("CalDataModuleLabel");
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
  hist->SetBins(nbins,0,nbins);
  hist->SetLineColor(color);
}

void cal::ShowWire::FillWaveforms(recob::Wire const& wire, raw::RawDigit const& rawdigit, TH1F* h_wire, TH1F* h_raw){

  const float pedestal = rawdigit.GetPedestal();
  for(size_t i=0; i<rawdigit.Samples(); i++){
    h_raw->SetBinContent(i,rawdigit.ADC(i) - pedestal);
    h_wire->SetBinContent(i,wire.Signal().at(i));
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

  for(auto const& wire : wireVector){

    raw::RawDigit const& rawdigit( *(wire.RawDigit()) );
    unsigned int channel = wire.Channel();
    size_t n_samples = rawdigit.Samples();

    TH1F* h_wire = tfs->make<TH1F>();
    TH1F* h_raw = tfs->make<TH1F>();

    SetHistogram(h_wire,"Wire",run,subrun,event,channel,n_samples,kBlue);
    SetHistogram(h_raw,"Raw",run,subrun,event,channel,n_samples,kRed);
    FillWaveforms(wire,rawdigit,h_wire,h_raw);
  }
    
}

DEFINE_ART_MODULE(cal::ShowWire)

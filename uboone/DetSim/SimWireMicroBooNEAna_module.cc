////////////////////////////////////////////////////////////////////////
// Class:       SimWireMicroBooNEAna
// Module Type: analyzer
// File:        SimWireMicroBooNEAna_module.cc
//
// Generated at Wed May 21 14:57:20 2014 by Matthew Toups using artmod
// from cetpkgsupport v1_05_04.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Utilities/InputTag.h"
#include "RawData/RawDigit.h"
#include "RawData/raw.h"
#include "fhiclcpp/ParameterSet.h"
#include <iostream>

namespace detsim {
  class SimWireMicroBooNEAna;
}

class detsim::SimWireMicroBooNEAna : public art::EDAnalyzer {
public:
  explicit SimWireMicroBooNEAna(fhicl::ParameterSet const & p);
  virtual ~SimWireMicroBooNEAna();

  void analyze(art::Event const & evt) override;


private:

  std::string fDigitModuleLabel;
  // Declare member data here.

};


detsim::SimWireMicroBooNEAna::SimWireMicroBooNEAna(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p)  // ,
 // More initializers here.
{
  fDigitModuleLabel= p.get< std::string >("DigitModuleLabel");
}

detsim::SimWireMicroBooNEAna::~SimWireMicroBooNEAna() {
  // Clean up dynamic memory and other resources here.
}

void detsim::SimWireMicroBooNEAna::analyze(art::Event const & evt) {

   art::Handle< std::vector<raw::RawDigit> > digitVecHandle;

   evt.getByLabel(fDigitModuleLabel, digitVecHandle);
   if(!digitVecHandle.isValid()) throw cet::exception("") << "NO!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!";

   if (!digitVecHandle->size())  std::cout << "BLAH\n";
   
   mf::LogInfo("CalWireMicroBooNE") << "CalWireMicroBooNE:: digitVecHandle size is " << digitVecHandle->size();

   std::cout << "CalWireMicroBooNE:: digitVecHandle size is " << digitVecHandle->size() << std::endl;

   //std::cout << "DigSize: " << digitVecHandle->size() << std::endl;
   // How many raw digits are there? 8192?
   for(size_t rdIter = 0; rdIter < digitVecHandle->size(); ++rdIter){ // ++ move

     // get the reference to the current raw::RawDigit
     art::Ptr<raw::RawDigit> digitVec(digitVecHandle, rdIter);

     unsigned short mySize = digitVec->Samples();
     //std::cout << "Size: " << mySize << std::endl;
     std::vector<float> holder(mySize);
     std::vector<short> rawadc(mySize);

     uint32_t channel = digitVec->Channel();

     for(size_t i = 0; i<digitVec->NADC(); i++) {
       //std::cout << "i: " << i << "\tfADC[" << i << "]: " << digitVec->ADC(i) << std::endl;
     }
     //std::cout << "fADC size: " << digitVec->NADC() << "\tEntry 0: " << digitVec->ADC(1236) << std::endl;
     // uncompress the data
     raw::Uncompress(digitVec->ADCs(), rawadc, digitVec->Compression());

     // loop over all adc values and subtract the pedestal
     float pdstl = digitVec->GetPedestal();

     std::string adcinfo="";
     for(size_t bin = 0; bin < mySize; ++bin) {
       //if(rawadc[bin]>0) std::cout << "Bin: " << bin << "\trawadc[" << bin << "]: " << rawadc[bin] << std::endl;
       adcinfo += rawadc[bin];
       if(bin+1<mySize) 
	 adcinfo += ",";
     } 
     channel += pdstl;
     //std::cout << printf("Channel: %d\tPedestal: %f\tADCs:\n%s\n",channel,pdstl,adcinfo.c_str());
   }

   // Implementation of required member function here.
}

DEFINE_ART_MODULE(detsim::SimWireMicroBooNEAna)

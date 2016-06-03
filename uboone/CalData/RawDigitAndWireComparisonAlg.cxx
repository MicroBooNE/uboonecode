/*!
 * Title:   RawDigitAndWireComparisonAlg
 * Author:  wketchum@lanl.gov
 * Inputs:  raw::RawDigit and recob::Wire objects
 * Outputs: Tree comparing raw waveform and "calibrated" waveform
 *
 */

#include "RawDigitAndWireComparisonAlg.h"

#include "lardataobj/RawData/raw.h" // raw::Uncompress()

caldata::RawDigitAndWireComparisonAlg::RawDigitAndWireComparisonAlg(fhicl::ParameterSet const& p):
  fRawDigitPropertiesAlg(p.get<fhicl::ParameterSet>("RawDigitPropertiesAlg")),
  fRecoWirePropertiesAlg(p.get<fhicl::ParameterSet>("RecoWirePropertiesAlg"))
{
}

void caldata::RawDigitAndWireComparisonAlg::SetROIOutputTree(TTree* tree){
  fROICompareTree = tree;
  SetupROIOutputTree();
}

void caldata::RawDigitAndWireComparisonAlg::SetupROIOutputTree(){
  fROICompareTree->Branch("event",&fROICompare.event,"event/i");
  fROICompareTree->Branch("run",&fROICompare.run,"run/i");
  fROICompareTree->Branch("channel",&fROICompare.channel,"channel/i");
  fROICompareTree->Branch("plane",&fROICompare.plane,"plane/i");
  fROICompareTree->Branch("wireROI_index",&fROICompare.wireROI_index,"wireROI_index/i");
  fROICompareTree->Branch("wireROI_start",&fROICompare.wireROI_start,"wireROI_start/i");
  fROICompareTree->Branch("wireROI_size",&fROICompare.wireROI_size,"wireROI_size/i");
  fROICompareTree->Branch("wireROI_integral",&fROICompare.wireROI_integral,"wireROI_integral/D");
  fROICompareTree->Branch("wireROI_peak",&fROICompare.wireROI_peak,"wireROI_peak/F");
  fROICompareTree->Branch("wireROI_peaktime",&fROICompare.wireROI_peaktime,"wireROI_peaktime/i");
  fROICompareTree->Branch("digit_regionMaxTime",&fROICompare.digit_regionMaxTime,"digit_regionMaxTime/i");
  fROICompareTree->Branch("digit_regionMinTime",&fROICompare.digit_regionMinTime,"digit_regionMinTime/i");
  fROICompareTree->Branch("digit_localPed",&fROICompare.digit_localPed,"digit_localPed/F");
  fROICompareTree->Branch("digit_localNoise",&fROICompare.digit_localNoise,"digit_localNoise/F");
  fROICompareTree->Branch("digit_regionMax",&fROICompare.digit_regionMax,"digit_regionMax/S");
  fROICompareTree->Branch("digit_regionMin",&fROICompare.digit_regionMin,"digit_regionMin/S");
  fROICompareTree->Branch("digit_regionSum",&fROICompare.digit_regionSum,"digit_regionSum/D");
  fROICompareTree->Branch("digit_regionSize",&fROICompare.digit_regionSize,"digit_regionSize/i");
  fROICompareTree->Branch("digit_isSignal",&fROICompare.digit_isSignal,"digit_isSignal/O");
}

void caldata::RawDigitAndWireComparisonAlg::RunROICompare(std::vector<recob::Wire> const& wireVector,
							  std::vector<raw::RawDigit> const& digitVector,
							  std::vector<unsigned int> const& assocVector,
							  unsigned int run,
							  unsigned int event){
  if(assocVector.size() != wireVector.size())
    throw std::runtime_error("Error: Association vector not same size as wire vector");

  fROICompare.event = event; fROICompare.run = run;

  for(size_t i_wire=0; i_wire<wireVector.size(); i_wire++)
    RunROICompare(wireVector[i_wire],digitVector[ assocVector[i_wire] ]);
}

void caldata::RawDigitAndWireComparisonAlg::RunROICompare(recob::Wire   const& wire,
							  raw::RawDigit const& digit){

  //fill channel and plane info
  fROICompare.channel = digit.Channel();
  fROICompare.plane   = wire.View();

  //first, run the raw digit through the WaveformProperties/ROI alg
  std::vector<short> ADCs(digit.Samples());
  raw::Uncompress(digit.ADCs(), ADCs, digit.Compression());
  fRawDigitPropertiesAlg.ProcessWaveform(ADCs);

  //now loop over the ROIs
  unsigned int range_index=0;
  for(auto const& roi : wire.SignalROI().get_ranges()){

    //fill the wire bits
    fROICompare.wireROI_index = range_index;
    fROICompare.wireROI_start = roi.begin_index();
    fROICompare.wireROI_size = roi.size();

    fROICompare.wireROI_integral = fRecoWirePropertiesAlg.GetSum(roi.data());
    fROICompare.wireROI_peak     = fRecoWirePropertiesAlg.GetMax(roi.data());

    auto max_loc = fRecoWirePropertiesAlg.GetMaxLocation(roi.data());
    fROICompare.wireROI_peaktime = std::distance(roi.data().cbegin(),max_loc) + fROICompare.wireROI_start;
    
    //now. get the iterator to the location of the peak time in the digit
    std::vector<short>::const_iterator peaktick = 
      ADCs.cbegin()+fROICompare.wireROI_peaktime;
    
    //what are the local pedestal, noise, and region type?
    fROICompare.digit_localPed   = fRawDigitPropertiesAlg.GetLocalPedestal(ADCs,peaktick);
    fROICompare.digit_localNoise = fRawDigitPropertiesAlg.GetLocalNoise(ADCs,peaktick);
    fROICompare.digit_isSignal   = fRawDigitPropertiesAlg.IsSignalRegion(ADCs,peaktick);

    //get the region of that guy
    typename util::WaveformPropertiesAlg<short>::Region digit_region = 
      fRawDigitPropertiesAlg.GetRegion(ADCs,peaktick);

    //now we can fill the digit min/max/sum/etc.
    fROICompare.digit_regionSum = fRawDigitPropertiesAlg.GetSum(digit_region);

    fROICompare.digit_regionMax = fRawDigitPropertiesAlg.GetMax(digit_region);
    auto digit_max_loc = fRawDigitPropertiesAlg.GetMaxLocation(digit_region);
    fROICompare.digit_regionMaxTime = std::distance(ADCs.cbegin(),digit_max_loc);

    fROICompare.digit_regionMin = fRawDigitPropertiesAlg.GetMin(digit_region);
    auto digit_min_loc = fRawDigitPropertiesAlg.GetMinLocation(digit_region);
    fROICompare.digit_regionMinTime = std::distance(ADCs.cbegin(),digit_min_loc);

    fROICompare.digit_regionSize = std::distance(digit_region.Start(),digit_region.End());
    
    //now fill the tree
    fROICompareTree->Fill();

    range_index++;
  }
  
  
}

/*!
 * Title:   RawDigitAndWireComparisonAlg
 * Author:  wketchum@lanl.gov
 * Inputs:  raw::RawDigit and recob::Wire objects
 * Outputs: Tree comparing raw waveform and "calibrated" waveform
 *
 */

#include "RawDigitAndWireComparisonAlg.h"

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
  fROICompareTree->Branch("roicompare",&fROICompare,fROICompare.leaflist.c_str());
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
  fRawDigitPropertiesAlg.ProcessWaveform(digit.fADC);

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
    fROICompare.wireROI_peaktime = std::distance(roi.data().cbegin(),max_loc)+roi.begin_index();

    
    //now. get the iterator to the location of the peak time in the digit
    std::vector<short>::const_iterator peaktick = 
      digit.fADC.cbegin()+fROICompare.wireROI_peaktime;
    
    //what are the local pedestal, noise, and region type?
    fROICompare.digit_localPed   = fRawDigitPropertiesAlg.GetLocalPedestal(digit.fADC,peaktick);
    fROICompare.digit_localNoise = fRawDigitPropertiesAlg.GetLocalNoise(digit.fADC,peaktick);
    fROICompare.digit_isSignal   = fRawDigitPropertiesAlg.IsSignalRegion(digit.fADC,peaktick);

    //get the region of that guy
    typename util::WaveformPropertiesAlg<short>::Region digit_region = 
      fRawDigitPropertiesAlg.GetRegion(digit.fADC,peaktick);

    //now we can fill the digit min/max/sum/etc.
    fROICompare.digit_regionSum = fRawDigitPropertiesAlg.GetSum(digit_region);

    fROICompare.digit_regionMax = fRawDigitPropertiesAlg.GetMax(digit_region);
    auto digit_max_loc = fRawDigitPropertiesAlg.GetMaxLocation(digit_region);
    fROICompare.digit_regionMaxTime = std::distance(digit.fADC.cbegin(),digit_max_loc);

    fROICompare.digit_regionMin = fRawDigitPropertiesAlg.GetMin(digit_region);
    auto digit_min_loc = fRawDigitPropertiesAlg.GetMinLocation(digit_region);
    fROICompare.digit_regionMinTime = std::distance(digit.fADC.cbegin(),digit_min_loc);

    //now fill the tree
    fROICompareTree->Fill();

    range_index++;
  }
  
  
}

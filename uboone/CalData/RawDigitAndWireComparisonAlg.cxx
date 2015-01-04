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


    //some couts
    if(fROICompare.channel==310){
      std::cout << fROICompare.wireROI_size << " " << fROICompare.wireROI_start << std::endl;
      unsigned int row=0;
      std::cout << "wire roi is \n\t" << row++ << " ";
      for(size_t i_roi=0; i_roi<roi.size(); i_roi++){
	std::cout << roi.data()[i_roi] << " ";
	if(i_roi%10==9) std::cout << "\n\t" << row++ << " ";
      }
      std::cout << std::endl;
      std::cout << fROICompare.wireROI_integral << " " << fROICompare.wireROI_peak << " " << fROICompare.wireROI_peaktime << std::endl;

      row=0;
      std::cout << "digit is \n\t" << row++ << " ";
      for(size_t i_digit=0; i_digit<digit.fADC.size(); i_digit++){
	std::cout << digit.fADC[i_digit] << " ";
	if(i_digit%10==9) std::cout << "\n\t" << row++ << " ";
      }

    }
    
    
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

    /*
    if(fROICompare.channel==310){
      std::cout << "Digit corresponding peack value is " << std::distance(digit.fADC.cbegin(),peaktick) << " " << *peaktick << std::endl;
      std::cout << "Region start,end = (" << std::distance(digit.fADC.cbegin(),digit_region.Start()) << ","
		<< std::distance(digit.fADC.cbegin(),digit_region.End()) << ")" << std::endl;

      for(auto tick=digit_region.Start(); tick!=digit_region.End(); tick++){
	std::cout << *tick << " ";
      }
      std::cout << std::endl;

      auto const sbt = fRawDigitPropertiesAlg.GetSignalBaselineTrio(digit.fADC,peaktick);
      std::cout << "Region start,end = (" << std::distance(digit.fADC.cbegin(),sbt.SignalRegion.Start()) << ","
		<< std::distance(digit.fADC.cbegin(),sbt.SignalRegion.End()) << ")" << std::endl;
      for(auto tick=sbt.SignalRegion.Start(); tick!=sbt.SignalRegion.End(); tick++){
	std::cout << *tick << " ";
      }
      std::cout << std::endl;
      std::cout << "Region start,end = (" << std::distance(digit.fADC.cbegin(),sbt.BaselineRegion_Pre.Start()) << ","
		<< std::distance(digit.fADC.cbegin(),sbt.BaselineRegion_Pre.End()) << ")" << std::endl;
      for(auto tick=sbt.BaselineRegion_Pre.Start(); tick!=sbt.BaselineRegion_Pre.End(); tick++){
	std::cout << *tick << " ";
      }
      std::cout << std::endl;
      std::cout << "Region start,end = (" << std::distance(digit.fADC.cbegin(),sbt.BaselineRegion_Post.Start()) << ","
		<< std::distance(digit.fADC.cbegin(),sbt.BaselineRegion_Post.End()) << ")" << std::endl;
      for(auto tick=sbt.BaselineRegion_Post.Start(); tick!=sbt.BaselineRegion_Post.End(); tick++){
	std::cout << *tick << " ";
      }
      std::cout << std::endl;
      
      double sum1 = fRawDigitPropertiesAlg.GetSum(sbt.BaselineRegion_Pre);
      size_t n1 = std::distance(sbt.BaselineRegion_Pre.Start(),sbt.BaselineRegion_Pre.End());
      double sum2 = fRawDigitPropertiesAlg.GetSum(sbt.BaselineRegion_Post);
      size_t n2 = std::distance(sbt.BaselineRegion_Post.Start(),sbt.BaselineRegion_Post.End());
      std::cout << sum1 << "/" << n1 << " " << sum2 << "/" << n2 << " " << (double)(sum1+sum2)/(n1+n2) << std::endl;;

      std::cout << fROICompare.digit_localPed << " " << fROICompare.digit_localNoise << " " << fROICompare.digit_isSignal << std::endl;
      std::cout << fROICompare.digit_regionSum << std::endl;
      std::cout << fROICompare.digit_regionMaxTime << " " << fROICompare.digit_regionMax << std::endl;
      std::cout << fROICompare.digit_regionMinTime << " " << fROICompare.digit_regionMin << std::endl;

      std::cout << "Number of signal regions is " << fRawDigitPropertiesAlg.GetNSignalRegions(digit.fADC) << std::endl;
      
      for(auto const& range : fRawDigitPropertiesAlg.GetSignalRegions(digit.fADC)){
      std::cout << "\t Signal Region start,end = (" << std::distance(digit.fADC.cbegin(),range.Start()) << ","
		<< std::distance(digit.fADC.cbegin(),range.End()) << ")" << std::endl;
      }

      for(auto const& range : fRawDigitPropertiesAlg.GetBaselineRegions(digit.fADC)){
	std::cout << "\t BaselineRegion start,end = (" << std::distance(digit.fADC.cbegin(),range.Start()) << ","
		  << std::distance(digit.fADC.cbegin(),range.End()) << ")" << std::endl;
      }
      
    }
    */

    //now fill the tree
    fROICompareTree->Fill();

    range_index++;
  }
  
  
}

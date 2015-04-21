/*!
 * Title:   CalibrationTPC Algs
 * Author:  wketchum@lanl.gov, dcaratelli@nevis.columbia.edu
 * Inputs:  raw::RawDigit
 * Outputs: Histograms and other nice data
 *
 * Description:
 * These are the actual functions that do things for electronics calibration.
 */

#include "CalibrationTPC_Algs.h"
#include <cmath>
#include <iostream>

/*
#include "TMinuit.h"
//TMinuit instance
//TMinuit *gMinuit;

Double_t _arglist[2];

std::vector<double> _data;
std::vector<double> _voltage;
std::vector<double> _error;

double gainModel(double x,double par[]){
  
  // Model for gain curve (linear)
  Double_t val = par[0]+par[1]*x;
  return val;
}


void chiSquared(int& npar,double* gin,double& f,double par[],int iflag){
  
  Double_t chi = 0;
  
  for (size_t n=0; n < _data.size(); n++)
     chi += ((_data[n]-gainModel(_voltage[n],par))/_error[n])*((_data[n]-gainModel(_voltage[n],par))/_error[n]);
  
  f = chi;
}
*/

namespace calibration {
  

  CalibrationAlgs::CalibrationAlgs()
  {}

  void CalibrationAlgs::PrepareGainModel(){
    std::cout << "preparing gain model" << std::endl;
    // Set gain model
    gainModel = new TF1("gainModel","[0]+[1]*x",0,1000);
    gainModel->SetParName(0,"const");
    gainModel->SetParName(1,"slope");
    std::cout << "done preparing gain model" << std::endl;
    return;
  }  

  //-------------------------------------------------------------------------
  bool CalibrationAlgs::hasCompressedRawDigit( std::vector<raw::RawDigit> const& rawDigitVector){
    for(auto const& rawDigit : rawDigitVector)
      if(rawDigit.Compression() != raw::Compress_t::kNone) return true;

    return false;
  }


  //-------------------------------------------------------------------------
  void CalibrationAlgs::genChanMap( std::vector<raw::RawDigit> const& rawDigit,
				    std::map< unsigned int, uint32_t > & chanmap){
    
    const unsigned int n_channels = rawDigit.size();
    unsigned int ChanMin = 20000;

    for (unsigned int ich=0; ich < n_channels; ich++){
      chanmap[ich] = rawDigit.at(ich).Channel();
      if ( rawDigit.at(ich).Channel() < ChanMin )
	ChanMin = rawDigit.at(ich).Channel();
    }

    //shift channels so that they start at 0
    for (unsigned int ich=0; ich < n_channels; ich++){
      chanmap[ich] -= ChanMin;
    }

  }


  //-------------------------------------------------------------------------
  void CalibrationAlgs::calcPedestal_BaselineRegion( util::UniqueRangeSet<std::vector<short>::const_iterator > const& rawData,
						     float & pedestal){

    pedestal = 0;
    size_t nsamples = 0;
    for ( auto const& range : rawData ){
      nsamples += std::distance(range.Start(), range.End());
      for (std::vector<short>::const_iterator tick = range.Start(); tick != range.End(); tick++)
	pedestal += (float)*tick;
    }

    if (nsamples > 0)
      pedestal /= nsamples;

    return;
  }


  void CalibrationAlgs::calcSignal_SignalRegion( util::UniqueRangeSet<std::vector<short>::const_iterator> const& rawData,
						 float & pedestal, float & maxADC, float & area, float & minADC, float & time){
    
    maxADC   = 0;
    area     = 0;
    minADC   = 4095;
    time     = 0;
    
    for ( auto const& range : rawData ){
      for (std::vector<short>::const_iterator tick = range.Start(); tick != range.End(); tick++){
	float thisTick = (float)*tick;
	if (thisTick-pedestal > 0)
	  area += (thisTick-pedestal);
	if (thisTick > maxADC) { maxADC = thisTick; time = (float)(tick - range.Start()); }
	if (thisTick < minADC) { minADC = thisTick; }
      }
    }
    
    return;
  }
  
  //-------------------------------------------------------------------------
  void CalibrationAlgs::calcNoise_BaselineRegion( util::UniqueRangeSet<std::vector<short>::const_iterator> const& rawData,
						  float const& pedestal,
						  float & noise){

    noise = 0;
    size_t nsamples = 0;
    for ( auto const& range : rawData ){
      nsamples += std::distance(range.Start(), range.End());
      for (std::vector<short>::const_iterator tick = range.Start(); tick != range.End(); tick++)
        noise += (pedestal - (float)*tick) * (pedestal - (float)*tick);
    }

    if(nsamples < 2){
      std::cerr << "Error in CalibrationTPC_Algs::calcNoise_SingleChannel" 
		<< "\nNumber of samples (" << nsamples << ")"
		<< " must be greater than 1!";
      noise=0;
      return;
    }

    noise = sqrt(noise / (nsamples - 1));

    return;
  }


  void CalibrationAlgs::calcGain(std::vector<float> const& voltages, std::vector<float> const& data,
				 std::vector<float> const& errors, float& chi,
				 float& gain, float& gainErr, float& intercept, float& interceptErr,
				 float& residualSum, float& residualSumSquared)
  {
    
    // For now, since I can't get TMinuit to work here, use a TGraph to perform fitting with ROOT
    // Arrays should be the size of inputs.
    // Make sure all the same dimension!
    if ( (voltages.size() != data.size()) || (voltages.size() != errors.size()) )
      std::cerr << "Gain calculation cannot proceed. voltages, data and error vectors not the same size!" << std::endl;

    const int n = voltages.size();
    float Data[n];
    float Errors[n];
    float Voltages[n];
    // Error for voltages: none for now
    float VoltErr[n];
    for (int i=0; i < n; i++)
      VoltErr[i] = 0.;

    std::copy(data.begin(), data.end(), Data);
    std::copy(errors.begin(), errors.end(), Errors);
    std::copy(voltages.begin(), voltages.end(), Voltages);

    gr = new TGraphErrors(n,Voltages,Data,VoltErr,Errors);
    gainModel->SetRange(voltages.front(), voltages.back());
    // ease the fitting by setting the params
    float s = (data.back()-data.front())/(voltages.back()-voltages.front());
    float b = data.front()-s*voltages.front(); 
    gainModel->SetParameter(0,b);
    gainModel->SetParameter(1,s);
    // Perform fit
    gr->Fit("gainModel","QN");
    // And now get fit results and errors
    intercept    = gainModel->GetParameter(0);
    interceptErr = gainModel->GetParError(0);
    gain         = gainModel->GetParameter(1);
    gainErr      = gainModel->GetParError(1);
    chi          = gainModel->GetChisquare();

    // Now calculate residuals
    residualSum = 0;
    residualSumSquared = 0;
    float diff = 0;
    for(int j=0; j < n; j++){
      diff = gainModel->Eval(voltages[j])-data[j];
      residualSum += diff;
      residualSumSquared += (diff*diff);
    }
    
    return;
  }



  
}

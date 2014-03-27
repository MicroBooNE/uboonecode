/*!
 * Title:   CalibrationTPC Algs
 * Author:  wketchum@lanl.gov
 * Inputs:  raw::RawDigit
 * Outputs: Histograms and other nice data
 *
 * Description:
 * These are the actual functions that do things for electronics calibration.
 */

#include "CalibrationTPC_Algs.h"
#include <cmath>
#include <iostream>

namespace calibration{

  //-------------------------------------------------------------------------
  void analyzeEmptyEvent( std::vector<raw::RawDigit> const& rawDigit,
			  std::vector<float> & pedestal,
			  std::vector<float> & noise,
			  std::vector<std::vector<float> > & noise_spectrum){
    

    calcPedestal(rawDigit, pedestal);
    calcNoise(rawDigit, pedestal, noise, noise_spectrum);

  }


  //-------------------------------------------------------------------------
  void calcPedestal( std::vector<raw::RawDigit> const& rawDigit,
		     std::vector<float> & pedestal){

    const unsigned int n_channels = rawDigit.size();
    
    for(unsigned int ich=0; ich<n_channels; ich++)
      calcPedestal_SingleChannel(rawDigit.at(ich).fADC,
				 pedestal.at(ich));
 
  }


  //-------------------------------------------------------------------------
  void calcPedestal_SingleChannel( std::vector<short> const& rawData,
				   float & pedestal){

    const unsigned int n_samples = rawData.size();

    pedestal = 0;
    for(unsigned int it=0; it<n_samples; it++)
      pedestal += rawData.at(it);

    pedestal = pedestal / n_samples;

  }

  //-------------------------------------------------------------------------
  void calcNoise( std::vector<raw::RawDigit> const& rawDigit,
		  std::vector<float> const& pedestal,
		  std::vector<float> & noise,
		  std::vector< std::vector<float> > & noise_spectrum){


    const unsigned int n_channels = rawDigit.size();
    
    for(unsigned int ich=0; ich<n_channels; ich++)
      calcNoise_SingleChannel(rawDigit.at(ich).fADC,
			      pedestal.at(ich),
			      noise.at(ich),
			      noise_spectrum.at(ich));
    
  }

  //-------------------------------------------------------------------------
  void calcNoise_SingleChannel( std::vector<short> const& rawData,
				float const& pedestal,
				float & noise,
				std::vector<float> & noise_spectrum){

    const unsigned int n_samples = rawData.size();

    //should never happen, but let's be safe, not sorry
    if(n_samples < 2){
      std::cerr << "Error in CalibrationTPC_Algs::calcNoise_SingleChannel" 
		<< "\nNumber of samples (" << n_samples << ")"
		<< " must be greater than 1!";
      noise=0;
      return;
    }

    noise=0;
    for(unsigned int it=0; it<n_samples; it++){
      noise += (rawData.at(it)-pedestal)*(rawData.at(it)-pedestal);

      // STUFF HERE FOR NOISE SPECTRUM!!!!
    }

    noise = sqrt(noise / (n_samples - 1));
  }

}

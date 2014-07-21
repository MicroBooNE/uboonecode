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
#include "TMath.h"
#include "TComplex.h"

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
  void analyzeGainEvent( std::vector<raw::RawDigit> const& rawDigit,
			 std::vector<float> & pedestal,
			 std::vector<float> & noise,
			 std::vector<float> & maxADC,
			 std::vector<float> & Area,
			 std::vector<float> & minADC,
			 std::vector<float> & maxTime,
			 int const& prePulseTicks){
    
    calcGain(rawDigit, pedestal, noise, maxADC, Area, minADC, maxTime, prePulseTicks);

  }


  //-------------------------------------------------------------------------
  void genChanMap( std::vector<raw::RawDigit> const& rawDigit,
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
  void calcPedestal( std::vector<raw::RawDigit> const& rawDigit,
		     std::vector<float> & pedestal){

    const unsigned int n_channels = rawDigit.size();
    
    for(unsigned int ich=0; ich<n_channels; ich++)
      calcPedestal_SingleChannel(rawDigit.at(ich).fADC,
				 pedestal.at(ich)  );

 
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
  void calcGain( std::vector<raw::RawDigit> const& rawDigit,
		 std::vector<float> & pedestal,
		 std::vector<float> & noise,
		 std::vector<float> & maxADC,
		 std::vector<float> & Area,
		 std::vector<float> & minADC,
		 std::vector<float> & maxTime,
		 int const& prePulseTicks){

    const unsigned int n_channels = rawDigit.size();
    
    for(unsigned int ich=0; ich<n_channels; ich++)
      calcGain_SingleChannel(rawDigit.at(ich).fADC,
			     pedestal.at(ich),
			     noise.at(ich),
			     maxADC.at(ich),
			     Area.at(ich),
			     minADC.at(ich),
			     maxTime.at(ich),
			     prePulseTicks);
 
  }


  void calcGain_SingleChannel( std::vector<short> const& rawData,
			       float & pedestal,
			       float & noise,
			       float & maxADC,
			       float & Area,
			       float & minADC,
			       float & time,
			       int const& prePulseTicks){

    const unsigned int n_samples = rawData.size();

    pedestal = 0;
    noise    = 0;
    maxADC   = 0;
    Area     = 0;
    minADC   = 0;
    time     = 0;

    for(unsigned int it=0; it<n_samples; it++){
      short thisADC = rawData.at(it);
      if ( (int)it < prePulseTicks ) { pedestal += thisADC; }
      if ( thisADC > maxADC ){
	maxADC = thisADC;
	time = it;
      }
      if ( it < 70 ) { Area += thisADC; }
      if ( thisADC < minADC ) { minADC = thisADC; }
    }

    pedestal /= (float)prePulseTicks;

    //Subtract baseline
    //maxADC -= pedestal;
    //minADC -= pedestal;

    //calculate RMS on pre-pulse to make sure it is not too high (what threshold?)
    for(int it=0; it < prePulseTicks; it++)
      noise += (rawData.at(it)-pedestal)*(rawData.at(it)-pedestal);
    noise = sqrt( noise / (prePulseTicks -1) );

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
    //prepare vector to store baseline subtracted value for FFT
    std::vector<float> BaselineSubtracted(n_samples, 0.);
    //std::vector<TComplex> noiseFrequency( n_samples/2+1, 0.);
    //for ( unsigned int i=0; i < n_samples/2+1; i++)
    //  noiseFrequency.push_back(TC0.0);
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
      BaselineSubtracted.at(it) = (rawData.at(it)-pedestal);
    }

    //do FFT
    //    art::ServiceHandle<util::LArFFT> fFFT;
    //    fFFT->DoFFT(BaselineSubtracted,noiseFrequency);
    //BaselineSubtracted.clear();

    for (unsigned int it=0; it<(n_samples/2+1); it++)
      noise_spectrum.at(it) = 0.0;

    noise = sqrt(noise / (n_samples - 1));
  }


}

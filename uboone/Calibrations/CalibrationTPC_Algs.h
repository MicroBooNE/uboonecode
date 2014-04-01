#ifndef CALIBRATIONTPC_ALGS_H
#define CALIBRATIONTPC_ALGS_H
/*!
 * Title:   CalibrationTPC Algs
 * Author:  wketchum@lanl.gov
 * Inputs:  raw::RawDigit
 * Outputs: Histograms and other nice data
 *
 * Description:
 * These are the actual functions that do things for electronics calibration.
 */

#include "TVirtualFFT.h"
#include "TH1.h"
#include <vector>
#include "RawData/RawDigit.h"

namespace calibration{

  void analyzeEmptyEvent( std::vector<raw::RawDigit> const& rawDigit,
			  std::vector<float> & pedestal,
			  std::vector<float> & noise,
			  std::vector<std::vector<float> > & noise_spectrum);

  void calcPedestal( std::vector<raw::RawDigit> const& rawDigit,
		     std::vector<float> & pedestal);
  void calcPedestal_SingleChannel( std::vector<short> const& rawData,
				   float & pedestal);
  
  void calcNoise( std::vector<raw::RawDigit> const& rawDigit,
		  std::vector<float> const& pedestal,
		  std::vector<float> & noise,
		  std::vector< std::vector<float> > & noise_spectrum);
  void calcNoise_SingleChannel( std::vector<short> const& rawData,
				float const& pedestal,
				float & noise,
				std::vector<float> & noise_spectrum);

}
#endif //CALIBRATIONTPC_ALGS_H

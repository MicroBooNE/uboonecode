#ifndef CALIBRATIONTPC_ALGS_H
#define CALIBRATIONTPC_ALGS_H
/*!
 * Title:   CalibrationTPC_Algs
 * Author:  wketchum@lanl.gov
 * Inputs:  raw::RawDigit
 * Outputs: Histograms and other nice data
 *
 * Description:
 * This is the header file for the CalibrationTPC_Algs, which are a collection
 * algorithms for processing electronics calibrations runs.
 */

#include <vector>
#include "RawData/RawDigit.h"

namespace calibration {

  void analyzeEmptyEvent( std::vector<raw::RawDigit> const& rawDigit,
			  std::vector<float> & pedestal,
			  std::vector<float> & noise,
			  std::vector<std::vector<float> > & noise_spectrum);

} //end namespace calibration

#endif //CALIBRATIONTPC_ALGS_H

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

#include <vector>
#include "lardataobj/RawData/RawDigit.h"
#include "lardata/Utilities/UniqueRangeSet.h"
#include <map>


#include "TMath.h"
#include "TComplex.h"
#include "TF1.h"
#include "TGraphErrors.h"

#include "lardata/Utilities/LArFFT.h"

namespace calibration{

  class CalibrationAlgs {

  public:

    CalibrationAlgs();

    virtual ~CalibrationAlgs(){}

    void PrepareGainModel();

    bool hasCompressedRawDigit( std::vector<raw::RawDigit> const& rawDigitVector);
    
    void genChanMap( std::vector<raw::RawDigit> const& rawDigit,
		     std::map< unsigned int, uint32_t > & chanmap);
    
    void calcPedestal_BaselineRegion( util::UniqueRangeSet<std::vector<short>::const_iterator> const& rawData,
				      float & pedestal);
    
    void calcSignal_SignalRegion( util::UniqueRangeSet<std::vector<short>::const_iterator> const& rawData,
				  float & pedestal,
				  float & maxADC,
				  float & area,
				  float & minADC,
				  float & time);
    
    void calcNoise_BaselineRegion( util::UniqueRangeSet<std::vector<short>::const_iterator> const& rawData,
				   float const& pedestal,
				   float & noise);
    
    void calcGain(std::vector<float> const& voltages, std::vector<float> const& data,
		  std::vector<float> const& errors, float& chi,
		  float& gain, float& gainErr, float& intercept, float& interceptErr,
		  float& residualSum, float& residualSumSquared);
    
  private:

    TGraphErrors *gr;
    TF1 *gainModel;
    

  }; // end class CalibrationAlgs

}
#endif //CALIBRATIONTPC_ALGS_H

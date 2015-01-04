#ifndef RAWDIGITANDWIRECOMPARISONALG_H
#define RAWDIGITANDWIRECOMPARISONALG_H

/*!
 * Title:   RawDigitAndWireComparisonAlg
 * Author:  wketchum@lanl.gov
 * Inputs:  raw::RawDigit and recob::Wire objects
 * Outputs: Tree comparing raw waveform and "calibrated" waveform
 *
 */

#include <vector>
#include <string>
#include <exception>

#include "fhiclcpp/ParameterSet.h"
#include "RawData/RawDigit.h"
#include "RecoBase/Wire.h"
#include "WaveformPropertiesAlg.h"

#include "TTree.h"

namespace caldata{

  class RawDigitAndWireComparisonAlg{

  public:
    
    RawDigitAndWireComparisonAlg(fhicl::ParameterSet const& p);

    void SetROIOutputTree(TTree*);

    void RunROICompare(std::vector<recob::Wire> const&,
		       std::vector<raw::RawDigit> const&,
		       std::vector<unsigned int> const&,
		       unsigned int, unsigned int);

  private:
    
    util::WaveformPropertiesAlg<short> fRawDigitPropertiesAlg;
    util::WaveformPropertiesAlg<float> fRecoWirePropertiesAlg;
    
    void SetupROIOutputTree();

    void RunROICompare(recob::Wire   const&,
		       raw::RawDigit const&);

    struct ROITreeComp{
      unsigned int event;
      unsigned int run;
      unsigned int channel;
      unsigned int plane;

      unsigned int wireROI_index;
      unsigned int wireROI_start;
      size_t       wireROI_size;
      double       wireROI_integral;
      float        wireROI_peak;
      unsigned int wireROI_peaktime;

      unsigned int digit_regionMaxTime;
      unsigned int digit_regionMinTime;
      float        digit_localPed;
      float        digit_localNoise;
      short        digit_regionMax;
      short        digit_regionMin;
      double       digit_regionSum;
      bool         digit_isSignal;      
    };
    
    ROITreeComp fROICompare;

    TTree* fROICompareTree;

    
  };

}//end namespace caldata


#endif

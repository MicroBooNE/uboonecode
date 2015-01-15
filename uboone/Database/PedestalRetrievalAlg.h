#ifndef PEDESTALRETRIEVALALG_H
#define PEDESTALRETRIEVALALG_H

/*!
 * Title:   PedestalRetrievalAlg
 * Author:  Brandon Eberly (eberly@fnal.gov)
 * Inputs:  Channel, timestamp
 * Outputs: Pedestal Mean and RMS
 */
 
#include "fhiclcpp/ParameterSet.h"


namespace dtbse{

  class PedestalRetrievalAlg{
  
    public:

      PedestalRetrievalAlg(fhicl::ParameterSet const& p);
      virtual ~PedestalRetrievalAlg(){};

      void reconfigure(fhicl::ParameterSet const& p);

      void GetPedestalMean(unsigned int channel, float& pedmean, ulong timestamp=0 ) const;
      void GetPedestalRMS(unsigned int channel,  float& pedrms,  ulong timestamp=0 ) const;
      void GetPedestal(unsigned int channel, float& pedmean, float& pedrms, ulong timestamp=0 ) const;

    private:

      bool  RetrieveFromDB(unsigned int channel, float& pedmean, float& pedrms, ulong timestamp) const;

      float fDefaultCollectionPedMean;
      float fDefaultCollectionPedRMS;

      float fDefaultInductionPedMean; 
      float fDefaultInductionPedRMS;

      //insert here private data member handle to pedestal database   
  };
}//end namespace dtbse

#endif

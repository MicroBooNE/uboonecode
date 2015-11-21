#ifndef UBADCBASE_CXX
#define UBADCBASE_CXX

#include "UBADCBase.h"

#include "DetectorInfoServices/DetectorClocksService.h" // lardata

namespace opdet {

  //--------------------
  UBADCBase::UBADCBase()
  //--------------------
  {
    auto const* ts = lar::providerFrom<detinfo::DetectorClocksService>();
    fTimeInfo = ts->OpticalClock();
    fDuration = fTimeInfo.TickPeriod();
    Reset();
  }


  //------------------------------------------------------------
  void UBADCBase::SetTimeInfo(const util::ElecClock &start_freq,
				double duration)
  //------------------------------------------------------------
  {
    if(duration < start_freq.TickPeriod())

      throw UBOpticalException(Form("Waveform length (%g ns) smaller than the tick size (%g ns)",
				    duration,
				    start_freq.TickPeriod()
				    )
			       );
    
    fTimeInfo = start_freq;
    fDuration = duration;
  }

  //---------------------
  void UBADCBase::Reset()
  //---------------------
  {
    fDuration = 0;
    fTimeInfo.SetTime(0);
  }

  //---------------------------------------------------------------
  void UBADCBase::Digitize(const std::vector<float>& orig,
			   std::vector<unsigned short>& res) const 
  //---------------------------------------------------------------
  {
    
    res.clear();
    res.reserve(orig.size());

    for(size_t i=0; i<orig.size(); ++i) {

      float v = orig.at(i);

      if(v > ((float)kADC_MAX)) res.push_back(kADC_MAX);

      else if(v < 0) res.push_back(0);
      
      else res.push_back( (unsigned short)(v));
      
    }

  }

}

#endif

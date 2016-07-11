#include <vector>

namespace fememu {

  class BasicTriggerConfig {

  public:
    int fDiscr0delay;
    int fDiscr3delay;
    
    int fDiscr0threshold;
    int fDiscr3threshold;

    int fDiscr0precount;
    
    int fDiscr0deadtime;
    int fDiscr3deadtime;
    
    int fDiscr0width;
    int fDiscr3width;

    int fMinReadoutTicks;
    int fFrontBuffer;
  };

  void basicTrigger( int BeamWinSize, int NChannels, const BasicTriggerConfig& config, const std::vector< std::vector<int> >& chwfms, std::vector<int>& vmaxdiff, std::vector<int>& vmulti );

}

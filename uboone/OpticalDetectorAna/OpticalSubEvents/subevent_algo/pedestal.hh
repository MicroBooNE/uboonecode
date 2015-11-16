#include <vector>

namespace subevent {

  double calcPedestal( std::vector< double >& wfm, int samplelength, double variance_threshold, double default_ped );

  double removePedestal( std::vector< double >& wfm, int samplelength, double variance_threshold, double default_ped );
  
}

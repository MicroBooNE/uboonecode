#include <vector>

namespace cpysubevent {
  void runCFdiscriminatorCPP( std::vector< int >& t_fire, std::vector< int >& amp_fire, std::vector< int >& maxt_fire, std::vector< int >& diff_fire,
			      double* waveform, int delay, int threshold, int deadtime, int width, int arrlen );
}

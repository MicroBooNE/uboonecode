#include <vector>

namespace subevent {

  void calcScintResponseCPP( std::vector< double >& fexpectation, 
			     int tstart, int tend, int maxt, float sig, float maxamp, float fastconst, float slowconst, float nspertick );
}

#ifndef _BEAM_HTTPRESPONSE_H
#define _BEAM_HTTPRESPONSE_H

#include <string>

namespace gov {
namespace fnal {
namespace uboone {
namespace beam {
  
using namespace gov::fnal::uboone;

class httpResponse {

 public:
  httpResponse();
  ~httpResponse();

  std::string memory;      // The buffer from HTTP response

};
}  // end of namespace beam
}  // end of namespace uboone
}  // end of namespace fnal
}  // end of namespace gov

#endif /* #ifndef _BEAM_HTTPRESPONSE_H */

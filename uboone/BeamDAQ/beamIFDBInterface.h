#ifndef _UBOONETYPES_BEAMIFDBINTERFACE_H
#define _UBOONETYPES_BEAMIFDBINTERFACE_H

#include <stdlib.h>
#include <cstring>
#include <stdio.h>
#include <string>
#include "httpResponse.h"
#include <curl/curl.h>

namespace gov {
namespace fnal {
namespace uboone {
namespace beam {
  
using namespace gov::fnal::uboone;


class beamIFDBInterface {

public:
  beamIFDBInterface();
  ~beamIFDBInterface();

  void GetData(const char *url, httpResponse* response);

private:
  //curl stuff
  CURL *fCURLHandle;
  
  static size_t writeMemoryCallback(void *cont, size_t size, size_t nmemb, void *userp)
  {
    httpResponse *response = (httpResponse *)userp;
    for (unsigned int c = 0; c<size*nmemb; c++) {
      char* buf=(char*) cont;
      (response->memory).push_back(buf[c]);
    }
    return size*nmemb; //tell curl how many bytes we handled
  };

};
}  // end of namespace beam
}  // end of namespace uboone
}  // end of namespace fnal
}  // end of namespace gov

#endif /* #ifndef _UBOONETYPES_BEAMIFDBINTERFACE_H */

#ifndef _BEAMRUNHEADER_H
#define _BEAMRUNHEADER_H

#include <sys/types.h>
#include <inttypes.h>

#include <iostream>
#include <string>
#include <map>
#include "boost/date_time/posix_time/posix_time.hpp"

namespace gov {
namespace fnal {
namespace uboone {
namespace beam {

  using namespace gov::fnal::uboone::beam;

class beamRunHeader {

 public:
  beamRunHeader();

  int   fRun;
  int   fSubRun;
  boost::posix_time::ptime fRunStart;
  boost::posix_time::ptime fRunEnd;
  std::map<std::string, uint> fCounter; 

  bool fFinished;

inline bool operator<(const gov::fnal::uboone::beam::beamRunHeader& rhs) const {
  return (this->fRun<rhs.fRun) || (this->fRun==rhs.fRun && this->fSubRun<rhs.fSubRun);
}
 private:
 

};
}  // end of namespace datatypes
}  // end of namespace uboone
}  // end of namespace fnal
}  // end of namespace gov

std::ostream & operator<<(std::ostream &os, const gov::fnal::uboone::beam::beamRunHeader &brh);

inline bool operator>=(const gov::fnal::uboone::beam::beamRunHeader& lhs, const gov::fnal::uboone::beam::beamRunHeader& rhs) {
  return (lhs.fRun>rhs.fRun) || (lhs.fRun==rhs.fRun && lhs.fSubRun>=rhs.fSubRun);
}

inline bool operator==(const gov::fnal::uboone::beam::beamRunHeader& lhs, const gov::fnal::uboone::beam::beamRunHeader& rhs) {
  return (lhs.fRun==rhs.fRun && lhs.fSubRun==rhs.fSubRun);
}
#endif /* #ifndef BEAMRUNHEADER_H */

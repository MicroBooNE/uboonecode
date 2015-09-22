#include "beamRunHeader.h"
#include <climits>

using namespace gov::fnal::uboone::beam;
using namespace boost::posix_time;

std::ostream & operator<<(std::ostream &os, const beamRunHeader &brh)
{
  os <<"RUN:    "   << brh.fRun << std::endl
     <<"SUBRUN: "   << brh.fSubRun << std::endl
     <<"START:  " << to_iso_extended_string(brh.fRunStart) << std::endl
     <<"END:    " << to_iso_extended_string(brh.fRunEnd)   << std::endl;
  
  std::map<std::string, uint>::const_iterator cit=brh.fCounter.begin();
  while (cit!=brh.fCounter.end() ) {
    os << cit->first<<" events: "<<cit->second<<std::endl;
    cit++;
  }
  return os;
}

beamRunHeader::beamRunHeader()
{
  fRun      = -1;
  fSubRun      = -1;
  fRunStart = from_time_t(0);
  fRunEnd   = from_time_t(INT_MAX); //should work till ~2038 (at least on this machine :)

  fFinished = false;
                          
}

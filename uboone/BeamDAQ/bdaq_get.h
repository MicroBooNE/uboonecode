#ifndef _BEAMDAQ_H
#define _BEAMDAQ_H
#include <iostream>
#include <fstream>

#include "datatypes/ub_BeamHeader.h"
#include "datatypes/ub_BeamData.h"
#include "beamIFDBInterface.h"
#include "beamRun.h"
#include "beamRunHeader.h"
#include "beamDAQConfig.h"

#include <messagefacility/MessageLogger/MessageLogger.h>

#include <boost/thread/thread.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/algorithm/string.hpp>
#include "boost/date_time/gregorian/gregorian.hpp"
#include "boost/date_time/posix_time/posix_time.hpp"
#include <boost/date_time/c_local_time_adjustor.hpp>

static boost::posix_time::time_duration fZoneOffset = boost::posix_time::second_clock::local_time()-boost::posix_time::second_clock::universal_time();

#endif /* #ifndef _BEAMDAQ_H */

#ifndef SIMPLECHCONFIG_CXX
#define SIMPLECHCONFIG_CXX
#include "SimpleChConfig.h"

namespace opdet {

  const std::map<unsigned int, float>& SimpleChConfig::GetFloat (const ChConfigType_t type)
  {
    if ( !fInitialized ) {
      doInitialization();
      fInitialized = true;
    }
    
    auto par_chdata = fFloatParams.find( type );
    if ( par_chdata!=fFloatParams.end() )
      return (*par_chdata).second;
    
    // should not compile if mistake made
    throw UBOpticalException("Invalid parameter type!");
  }

  const std::map<unsigned int, int>& SimpleChConfig::GetInt (const ChConfigType_t type)
  {
    if ( !fInitialized ) {
      doInitialization();
      fInitialized = true;
    }
    
    auto par_chdata = fIntParams.find( type );
    if ( par_chdata!=fIntParams.end() )
      return (*par_chdata).second;
    
    // should not compile if mistake made
    throw UBOpticalException("Invalid parameter type!");
  }

  float SimpleChConfig::GetFloat( const ChConfigType_t type, const unsigned int ch)
  {
    
    if ( !fInitialized ) {
      doInitialization();
      fInitialized = true;
    }
    
    auto const& par_chdata = GetFloat( type );
    auto par_value = par_chdata.find( ch );
    if ( par_value!=par_chdata.end() )
      return (*par_value).second;
    
    char err[100];
    sprintf(err, "Invalid channel number=%d provided!",ch );
    throw UBOpticalException(err);
  }

  int SimpleChConfig::GetInt( const ChConfigType_t type, const unsigned int ch)
  {
    
    if ( !fInitialized ) {
      doInitialization();
      fInitialized = true;
    }
    
    auto const& par_chdata = GetInt( type );
    auto par_value = par_chdata.find( ch );
    if ( par_value!=par_chdata.end() )
      return (*par_value).second;
    
    char err[100];
    sprintf(err, "Invalid channel number=%d provided!",ch );
    throw UBOpticalException(err);
  }
}

#endif


#include "SimpleChConfig.h"

namespace opdet {
  const std::map<unsigned int, float>& SimpleChConfig::GetParameter(const ChConfigType_t type)
  {
    
    if ( !fInitialized ) {
      doInitialization();
      fInitialized = true;
    }
    
    auto par_chdata = fParams.find( (unsigned int)type );
    if ( par_chdata!=fParams.end() )
      return (*par_chdata).second;
    
    // should not compile if mistake made
    throw UBOpticalException("Invalid parameter type!");
  }
  
  
  float SimpleChConfig::GetParameter( const ChConfigType_t type, const unsigned int ch)
  {
    
    if ( !fInitialized ) {
      doInitialization();
      fInitialized = true;
    }
    
    auto par_chdata = GetParameter( type );
    auto par_value = par_chdata.find( ch );
    if ( par_value!=par_chdata.end() )
      return (*par_value).second;
   
    char err[100];
    sprintf(err, "Invalid channel number=%d provided!",ch );
    throw UBOpticalException(err);
  }
  
}

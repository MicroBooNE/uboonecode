
////////////////////////////////////////////////////////////////////////
//
//    trivia : Class to hold the Muon Counter System (MuCSDTOFFSET) data, 
//             a.k.a. the best class there's ever been, May 2015 


#include "MuCSDTOffset.h"

#include "cetlib/exception.h"

namespace MuCS {

MuCSDTOffset::MuCSDTOffset()
{}

MuCSDTOffset::~MuCSDTOffset()
{}
  
MuCSDTOffset::MuCSDTOffset( double inoffset ) 
{
  offset = inoffset;
 
}
  
double MuCSDTOffset::getoffset() const
{
  return offset;
  
}


}

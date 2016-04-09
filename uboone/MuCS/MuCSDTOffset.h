
////////////////////////////////////////////////////////////////////////
//
//    trivia : Class to hold the Muon Counter System (MuCSDTOFFSET) data, 
//             a.k.a. the best class there's ever been, May 2015 


#ifndef MUCSDTOFFSET_H
#define MUCSDTOFFSET_H

#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <iosfwd>

namespace MuCS {
 
class MuCSDTOffset 
{

 public:
  MuCSDTOffset();
  virtual ~MuCSDTOffset();  
  
  MuCSDTOffset( double offset ); 
  
  double getoffset() const;

    
 private:
  
  
  double offset;

    
};

}

#endif 


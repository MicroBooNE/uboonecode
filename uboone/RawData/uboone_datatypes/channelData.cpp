#include <vector>
#include <iostream>
#include "channelData.h"

using namespace gov::fnal::uboone::datatypes;

void channelData::decompress(){

  unsigned int i_bit,i_zero; // some counters we will use

  const size_t size16 = sizeof(uint16_t);
  std::unique_ptr<uint16_t> word(new uint16_t);
  std::vector<uint16_t> uncompressed_vector;

  for(unsigned int iword = 0; iword < channel_data_size; iword++){

    std::copy(getChannelDataPtr() + iword*size16,
	      getChannelDataPtr() + (iword+1)*size16,
	      (char*)word.get());

    //if it's not a compressed word, just put it on the uncompressed vector
    if( (*word & 0xf000)==0x0000 )
      uncompressed_vector.push_back(*word);
    
    else if( (*word & 0x8000)==0x8000 ){

      i_bit=0; //initialize bit counter

      // now it's time to go look for the codes
      while(i_bit<15){

	if( (((*word)>>i_bit)&0x1) == 0x1 ){ //here is our code

	  i_zero = 0; //initialize zero counter
	  while(i_bit<15){

	    if( (( (*word)>>(i_bit+1) ) & 0x1) == 0x0 ){
	      i_zero++;
	    }
	    else if( (( (*word)>>(i_bit+1) ) & 0x1) == 0x1 ){

	      if (i_zero==0){//difference: current-previous= 0
		uncompressed_vector.push_back( (*word)&0x7ff );
	      }
	      else if (i_zero==1){//difference: current-previous= -1
		uncompressed_vector.push_back( ((*word)&0x7ff) - 1 );
	      }
	      else if (i_zero==2){//difference +1
		uncompressed_vector.push_back( ((*word)&0x7ff) + 1 );
	      }
	      else if (i_zero==3){//difference -2
		uncompressed_vector.push_back( ((*word)&0x7ff) - 2 );
	      }
	      else if (i_zero==4){//difference +2
		uncompressed_vector.push_back( ((*word)&0x7ff) + 2 );
	      }
	      else if (i_zero==5){//difference -3
		uncompressed_vector.push_back( ((*word)&0x7ff) - 3 );
	      }
	      else if (i_zero==6){//difference +3
		uncompressed_vector.push_back( ((*word)&0x7ff) + 3 );
	      }
	      else {
		std::cout << "Something went wrong?" << std::endl;
		return;
	      }

	      i_zero=0;
	    }//end else if 0x1

	    i_bit++;
	  }//end inside while over bit counter
	}//end if we have a code

	else{
	  i_bit++;
	}

      }//end outside while

    }//end if huffman compressed word

    else{
      std::cout << "ERROR!!!!!" << std::endl;
      return;
    }

  }//end for loop over data words


  channel_data_size = uncompressed_vector.size()*size16;
  setChannelDataPtr((char*)(uncompressed_vector.data()));

}

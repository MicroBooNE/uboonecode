#include "cardDataPMT.h"

using namespace gov::fnal::uboone::datatypes;


char* cardDataPMT::getCardDataPtr(){
  
  if(cardData_IO_mode >= IO_GRANULARITY_CHANNEL){
    std::cout << "ERROR! Granularity is above card level." 
	      << "Cannot return pointer to card data!" << std::endl;
    return nullptr;
  }
  else{
    return card_data_ptr.get();
  }
}

void cardDataPMT::setCardDataPtr(char* ptr){

  if(cardData_IO_mode >= IO_GRANULARITY_CHANNEL){
    std::cout << "ERROR! Granularity is above card level." 
	      << "Cannot set pointer to card data!" << std::endl;
  }
  else{
    card_data_ptr.reset(ptr);
  }
}

void cardDataPMT::updateIOMode(uint8_t new_mode){

  //we are already at card granularity...so get out if that's the case
  if(new_mode <= IO_GRANULARITY_CARD)
    return;

  if(new_mode >= IO_GRANULARITY_CHANNEL && cardData_IO_mode < IO_GRANULARITY_CHANNEL){
    size_t total_data_read = 0;

    //get the channel header word
    std::unique_ptr<pmt_data_header_t> memblkDH(new pmt_data_header_t);
    std::copy(getCardDataPtr() + total_data_read,
	      getCardDataPtr() + total_data_read + sizeof(pmt_data_header_t),
	      (char*)memblkDH.get());
    pmt_data_header.setDataHeader(*memblkDH);
    total_data_read += sizeof(pmt_data_header_t);
    
    // std::cout << "Channel header " << std::hex << *channel_header << std::endl;
    
    FillPMTChannels(total_data_read);
        
    std::unique_ptr<pmt_data_trailer_t> memblkDT(new pmt_data_trailer_t);
    std::copy(getCardDataPtr() + total_data_read,
	      getCardDataPtr() + total_data_read + sizeof(pmt_data_trailer_t),
	      (char*)memblkDT.get());
    pmt_data_trailer.setDataTrailer(*memblkDT);
    total_data_read += sizeof(pmt_data_trailer_t);
    
    
    card_data_ptr.reset();
    cardData_IO_mode = IO_GRANULARITY_CHANNEL;
  }//endif channel granularity update

}

void cardDataPMT::insertChannel(int channel_number, channelDataPMT chD){
  channel_map.insert(std::pair<int,channelDataPMT>(channel_number,chD));
}

void cardDataPMT::insertWindow(windowHeaderPMT wH, windowDataPMT wD){
  int channel_number = wH.getChannelNumber();
  
  if(channel_map.find(channel_number) == channel_map.end()){
    channelDataPMT chD;
    chD.insertWindow(wH,wD);
    insertChannel(channel_number,chD);
  }
  else{
    (channel_map.find(channel_number)->second).insertWindow(wH,wD);
  }

}

void cardDataPMT::FillPMTChannels(size_t &total_data_read){

  bool full_header = true;
  std::unique_ptr<pmt_window_header_t> memblkWH(new pmt_window_header_t);

  const size_t size16 = sizeof(uint16_t);
  std::unique_ptr<uint16_t> word(new uint16_t);

  while(total_data_read < (card_data_size - sizeof(pmt_data_trailer_t))){

    if(full_header){
      std::copy(getCardDataPtr() + total_data_read,
		getCardDataPtr() + total_data_read + sizeof(pmt_window_header_t),
		(char*)memblkWH.get());
      total_data_read += sizeof(pmt_window_header_t);
    }
    else{
      std::copy(getCardDataPtr() + total_data_read,
		getCardDataPtr() + total_data_read + size16,
		(char*)word.get());
      total_data_read += size16;
      memblkWH->frame_and_sample1 = *word;
      
      std::copy(getCardDataPtr() + total_data_read,
		getCardDataPtr() + total_data_read + size16,
		(char*)word.get());
      total_data_read += size16;
      memblkWH->sample2 = *word;

      full_header=true;

    }
       
    size_t window_data_size=0;
    char* window_data_begin_ptr = getCardDataPtr() + total_data_read;
    while(1){
      std::copy(getCardDataPtr() + total_data_read,
		getCardDataPtr() + total_data_read + size16,
		(char*)word.get());
      total_data_read += size16;
      window_data_size += size16;

      //std::cout << std::hex << *word << std::endl;
      
      if( (*word & 0x3000)==0x3000) {


	std::copy(getCardDataPtr() + total_data_read,
		  getCardDataPtr() + total_data_read + size16,
		  (char*)word.get());

	//std::cout << std::hex << *word << " " << (*word & 0x3000) << std::endl;

	if( (*word & 0x3000)==0x2000){ //we have more adc words/secondary header words...
	  //std::cout << "We have a partial header coming up!" << std::endl;
	  full_header=false;
	}
	break;
      }
      
    }

    windowHeaderPMT windowH(*memblkWH);
    windowDataPMT windowD(window_data_size,window_data_begin_ptr);
    
    insertWindow(windowH,windowD);
    
  }
 
}


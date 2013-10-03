#include "crateDataPMT.h"

using namespace gov::fnal::uboone::datatypes;

char* crateDataPMT::getCrateDataPtr(){
  
  if(crateData_IO_mode >= IO_GRANULARITY_CARD){
    std::cout << "ERROR! Granularity is above crate level." 
	      << "Cannot return pointer to crate data!" << std::endl;
    return nullptr;
  }
  else {
    return crate_data_ptr.get();
  }
}

void crateDataPMT::setCrateDataPtr(char* ptr){

  if(crateData_IO_mode >= IO_GRANULARITY_CARD){
    std::cout << "ERROR! Granularity is above crate level." 
	      << "Cannot set pointer to crate data!" << std::endl;
  }
  else {
    crate_data_ptr.reset(ptr);
  }
}

void crateDataPMT::updateIOMode(uint8_t new_mode){

  //we are already at crate granularity...so get out if that's the case
  if(new_mode <= IO_GRANULARITY_CRATE)
    return;

  //  const size_t size16 = sizeof(uint16_t);

  if(new_mode >= IO_GRANULARITY_CARD && crateData_IO_mode < IO_GRANULARITY_CARD){

    size_t data_read = 0;
    std::unique_ptr<event_header_t> memblkEH(new event_header_t);
    std::unique_ptr<event_trailer_t> memblkET(new event_trailer_t);

    std::copy(getCrateDataPtr() + data_read,
	      getCrateDataPtr() + data_read + sizeof(event_header_t),
	      (char*)memblkEH.get());
    event_header.setEventHeader(*memblkEH);
    data_read += sizeof(event_header_t);
    // std::cout << "crateData.cpp read event_header: 0x" << std::hex << memblkEH->header << std::endl;
    
    int cards_read = 0;
    while(1){
            // std::cout << "crateData.cpp Reading card header at position " << data_read << std::endl;
      std::unique_ptr<pmt_card_header_t> memblkCardH(new pmt_card_header_t);
      std::copy(getCrateDataPtr() + data_read,
		getCrateDataPtr() + data_read + sizeof(pmt_card_header_t),
		(char*)memblkCardH.get());
      data_read += sizeof(pmt_card_header_t);
      
      cardHeaderPMT cardH(*memblkCardH);
      size_t cardDataSize = cardH.getCardDataSize();

        // std::cout << "Card header ...\n"
  // << std::hex << memblkCardH->id_and_module << " " << memblkCardH->word_count << " " 
  // << memblkCardH->event_number << " " << memblkCardH->frame_number<< " " << memblkCardH->checksum << std::dec << std::endl;


      std::shared_ptr<char> card_data(new char[cardDataSize]);
      std::copy(getCrateDataPtr() + data_read,
		getCrateDataPtr() + data_read + cardDataSize,
		(char*)card_data.get());
      //wait to increment data_read until after updating channel granularity

      cardDataPMT cardD(card_data,cardDataSize);
      if(new_mode == IO_GRANULARITY_CHANNEL)
     	cardD.updateIOMode(new_mode);

      //now increment the data_read variable
      data_read += cardDataSize;

      insertCard(cardH,cardD);
      cards_read++;
      std::copy(getCrateDataPtr() + data_read,
		getCrateDataPtr() + data_read + sizeof(event_trailer_t),
		(char*)memblkET.get());
      if(memblkET->trailer == 0xe0000000) break;
    }
    event_trailer.setEventTrailer(*memblkET);
    data_read += sizeof(event_trailer_t);
    // std::cout << "crateData.cpp read " << std::dec << cards_read << " cards with " << data_read << " bytes." << std::endl;    
    crate_data_ptr.reset();

    crateData_IO_mode = new_mode;
  } //endif on IO_GRANULARITY_CARD update

  if(new_mode == IO_GRANULARITY_CHANNEL && crateData_IO_mode < IO_GRANULARITY_CHANNEL){

    // this code activated when current granularity is card, wanted is channel
    std::map<cardHeaderPMT,cardDataPMT>::iterator card_it;
    for( card_it = card_map.begin(); card_it != card_map.end(); card_it++){
      (card_it->second).updateIOMode(new_mode);      
    }

    crateData_IO_mode = new_mode; //eventRecords io_mode

  }//endif on IO_GRANULARITY_CHANNEL update

}

void crateDataPMT::insertCard(cardHeaderPMT cH, cardDataPMT cD){
  card_map.insert(std::pair<cardHeaderPMT,cardDataPMT>(cH,cD));
}

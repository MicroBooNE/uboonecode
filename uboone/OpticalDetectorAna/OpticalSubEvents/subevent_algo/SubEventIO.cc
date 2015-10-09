#include "SubEventIO.hh"

namespace subevent {
  
  SubEventIO::SubEventIO( std::string filename, std::string mode_ )
    : mode(mode_) {
    
    eventid = -1;
    chmaxamp = 0.;
    nsubevents = 0;
    subevents = new SubEventList();

    if ( mode=='w' ) {
      fFile = new TFile( filename.c_str(), "RECREATE" );
      defineTree();
    }
  }

  SubEventIO::~SubEventIO() {
  }

  void SubEventIO::defineTree() {
    fTree = new TTree("subevents", "SubEvent Tree");
    //b_subeventlist = new TBranch( "subeventlist", subevents );
    fTree->Branch( "event", &eventid );
    fTree->Branch( "chmaxamp", &chmaxamp );
    fTree->Branch( "nsubevents", &nsubevents );
    fTree->Branch( "subeventlist", subevents );
  }

  void SubEventIO::transferSubEventList( SubEventList* subeventsrc ) {
    for ( SubEventListIter it=subeventsrc->begin(); it!=subeventsrc->end(); it++ ) {
      subevents->add( std::move( *it ) );
    }
  }

  void SubEventIO::fill() {
    nsubevents = subevents->size();
    fTree->Fill();
  }

  void SubEventIO::write() {
    fTree->Write();
  }


}

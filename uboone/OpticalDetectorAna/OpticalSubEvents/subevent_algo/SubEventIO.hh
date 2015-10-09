#ifndef __SUBEVENTIO__
#define __SUBEVENTIO__

#include <string>
#include "TTree.h"
#include "TFile.h"
#include "TBranch.h"
#include "SubEventList.hh"

namespace subevent {

  class SubEventIO {
    
  public:

    SubEventIO( std::string filename, std::string mode );
    ~SubEventIO();

    std::string mode; // r or w

    TFile* fFile;
    TTree* fTree;
    TBranch* b_subeventlist;
    SubEventList* subevents;
    int eventid;
    int nsubevents;
    double chmaxamp;

    void defineTree();
    void transferSubEventList( SubEventList* subevents );
    void write();
    void fill();
    void clearlist() { subevents->clear(); };
    TTree* getTree() { return fTree; };

    bool branchmade;
  };
}

#endif

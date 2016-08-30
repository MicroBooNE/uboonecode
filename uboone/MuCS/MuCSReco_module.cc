
////////////////////////////////////////////////////////////////////////
// Class:       MuCSReco
// Module Type: producer
////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////
//
//    trivia : Reco stage of the MuCSMerger process. 
//             Adds angle and position info to hit patterns based on MC.
//    author : Matt Bass
//    e-mail : Matthew.Bass@physics.ox.ac.uk
//
////////////////////////////////////////////////////////////////////////

#ifndef MUCSRECO_H
#define MUCSRECO_H 

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "TText.h"
#include "TTimeStamp.h"

#include <memory>
#include <iostream>
#include "vector"
#include "lardataobj/RawData/TriggerData.h"
#include "TH1F.h"
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include <sqlite3.h> 
#include "MuCSData.h"
#include "MuCSRecoData.h"

using namespace std;

class MuCSReco;

class MuCSReco : public art::EDProducer {
public:
  explicit MuCSReco( fhicl::ParameterSet const &pset );
  virtual ~MuCSReco();
  
  void reconfigure( fhicl::ParameterSet const &pset ); // override;
  void produce( art::Event &evt ) override;
      
private:
  void openDB();
  template<typename ... Args> string string_format( const std::string& format, Args ... args );
  sqlite3* db; //pointer to sqlite3 database object
  struct dbfields {
   Int_t matches;
   Int_t nentries;
   Float_t p;
   Float_t q;
   Float_t p_rms;
   Float_t q_rms;
  };
  void execFill(const std::vector<int> hitsa, const std::vector<int> hitsb, const std::string tablename, dbfields &ldbfields);

  //fcl parameters
  std::string fInputDB; //sqlite db containing hits info
  Float_t fTopBoxy; //y position to use for positions (top of top box)

};

void MuCSReco::reconfigure( fhicl::ParameterSet const &p ){
    fInputDB = p.get< std::string >( "InputDB" );
    fTopBoxy = p.get< Float_t >( "TopBoxy" );
    return;
}

MuCSReco::MuCSReco( fhicl::ParameterSet const &pset ){
  this->reconfigure( pset );
  
  produces< std::vector<MuCS::MuCSRecoData> >();  
  
  openDB();
}

void MuCSReco::openDB(){
  //open database file
  int res=sqlite3_open(fInputDB.c_str(),&db);
  if (res!= SQLITE_OK)
    throw cet::exception("MuCSReco") << "Error opening db: (" <<fInputDB<<") ("<<res<<"): " << sqlite3_errmsg(db) << "; memory used:<<"<<sqlite3_memory_used()<<"/"<<sqlite3_memory_highwater(0);
  else
    mf::LogInfo("MuCSReco")<<"Opened db "<< fInputDB<<"\n";
}
  
MuCSReco::~MuCSReco(){
  sqlite3_close(db);
}

//string version of sprintf for query formatting
template<typename ... Args> string MuCSReco::string_format( const std::string& format, Args ... args ){
    size_t size = snprintf( nullptr, 0, format.c_str(), args ... ) + 1; // Extra space for '\0'
    unique_ptr<char[]> buf( new char[ size ] ); 
    snprintf( buf.get(), size, format.c_str(), args ... );
    return string( buf.get(), buf.get() + size - 1 ); // We don't want the '\0' inside
}

void MuCSReco::execFill(const std::vector<int> hitsa, const std::vector<int> hitsb, const std::string tablename, dbfields &ldbfields){
  //build and execute query, populate query fields
  //this db query has three steps:
    //create temp table with all hit patterns that match the hits
    //select number of matches, nentries, position (p), theta (q), p_rms_squared, q_rms_squared from this table in aggregate
      //in aggregrate means a weighted avg/stddev is computed if multiple hits are returned
    //drop the temp table
  //Note this will always return one row if there are matches, no rows otherwise
  int res=0; //store results from sqlite_steps
  std::string kStatement1("CREATE TEMPORARY TABLE hitsel as select * from %s where (%s) AND (%s) AND (%s) AND (%s);");
  std::stringstream filters[4];
  for (unsigned int i =0;i<hitsa.size();i++){
    filters[0]<< (i>0 ?" OR ": "") << "hita1=" << hitsa[i];
    filters[1]<< (i>0 ?" OR ": "") << "hita2=" << hitsa[i];
  }
  for (unsigned int i =0;i<hitsb.size();i++){
    filters[2]<< (i>0 ?" OR ": "") << "hitb1=" << hitsb[i];
    filters[3]<< (i>0 ?" OR ": "") << "hitb2=" << hitsb[i];
  }
  sqlite3_stmt *statement1;
  std::string kthisStatement1=string_format(kStatement1,tablename.c_str(),filters[0].str().c_str(),filters[1].str().c_str(),filters[2].str().c_str(),filters[3].str().c_str());
  mf::LogInfo("MuCSReco")<<"Executing: "<<kthisStatement1<<"\n";
  if (sqlite3_prepare(db, kthisStatement1.c_str(), -1, &statement1, 0 )!= SQLITE_OK)
    throw cet::exception("MuCSReco") << "Error preparing statement!: " << sqlite3_errmsg(db) << "\n";
  if(sqlite3_step(statement1)!=SQLITE_DONE) 
    throw cet::exception("MuCSReco") << "Error stepping statement!: " << sqlite3_errmsg(db) << "\n";
  if(sqlite3_finalize(statement1)!=SQLITE_OK) 
    throw cet::exception("MuCSReco") << "Error finalizing statement!: " << sqlite3_errmsg(db) << "\n";
    
  std::string kStatement2("select count(nentries) as matches,"
    "sum(nentries) as nentries,"
    "sum(nentries*p/(select sum(nentries) from hitsel)) as p,"
    "sum(nentries*q/(select sum(nentries) from hitsel)) as q,"
    "case when (SELECT COUNT(*) from hitsel)>1 then sum(nentries*(p-(select sum(nentries*p/(select sum(nentries) from hitsel)) from hitsel))*(p-(select sum(nentries*p/(select sum(nentries) from hitsel)) from hitsel))/(select sum(nentries) from hitsel)) else p_rms*p_rms end as p_rms_squared,"
    "case when (SELECT COUNT(*) from hitsel)>1 then sum(nentries*(q-(select sum(nentries*q/(select sum(nentries) from hitsel)) from hitsel))*(q-(select sum(nentries*q/(select sum(nentries) from hitsel)) from hitsel))/(select sum(nentries) from hitsel)) else q_rms*q_rms end as q_rms_squared "
    "from hitsel;"); 
  sqlite3_stmt *statement2;
  mf::LogInfo("MuCSReco")<<"Executing: "<<kStatement2<<"\n";
  if ( sqlite3_prepare(db, kStatement2.c_str(), -1, &statement2, 0 )== SQLITE_OK   ){
    
    res = sqlite3_step(statement2);
    if ( res == SQLITE_ROW ){
      ldbfields.matches=sqlite3_column_int(statement2,0);
      ldbfields.nentries=sqlite3_column_double(statement2,1);
      ldbfields.p=sqlite3_column_double(statement2,2);
      ldbfields.q=sqlite3_column_double(statement2,3);
      ldbfields.p_rms=sqrt(sqlite3_column_double(statement2,4));
      ldbfields.q_rms=sqrt(sqlite3_column_double(statement2,5));
    }else{
      mf::LogInfo("MuCSReco")<<"No matches returned("<<res<<")!: " << sqlite3_errmsg(db) << "\n";
      ldbfields.matches=0;
      ldbfields.nentries=0;
      ldbfields.p=0.;
      ldbfields.q=0.;
      ldbfields.p_rms=0.;
      ldbfields.q_rms=0.;
    }
  }else{
    throw cet::exception("MuCSReco") << "Error preparing statement!: " << sqlite3_errmsg(db) << "\n";
  }
  if(sqlite3_finalize(statement2)!=SQLITE_OK) 
    throw cet::exception("MuCSReco") << "Error finalizing statement!: " << sqlite3_errmsg(db) << "\n";

  std::string kStatement3("DROP TABLE hitsel;");   
  sqlite3_stmt *statement3;
  mf::LogInfo("MuCSReco")<<"Executing: "<<kStatement3<<"\n";
  if (sqlite3_prepare(db, kStatement3.c_str(), -1, &statement3, 0 )!= SQLITE_OK){
    throw cet::exception("MuCSReco") << "Error preparing statement!: " << sqlite3_errmsg(db) << "\n";
  }
  if(sqlite3_step(statement3)!=SQLITE_DONE) 
    throw cet::exception("MuCSReco") << "Error stepping statement!: " << sqlite3_errmsg(db) << "\n";
  if(sqlite3_finalize(statement3)!=SQLITE_OK) 
    throw cet::exception("MuCSReco") << "Error finalizing statement!: " << sqlite3_errmsg(db) << "\n";
}

void MuCSReco::produce( art::Event &evt ){

  //get MuCSData object
  //art::Handle< std::vector<MuCS::MuCSData> > mucsdatahandle;
  //evt.getByLabel("merger",mucsdatahandle);
  //std::vector<MuCS::MuCSData> const& mucsdata = *mucsdatahandle;
  std::vector< art::Handle< std::vector<MuCS::MuCSData> > > mucslist;
  evt.getManyByType( mucslist );
  art::Handle< std::vector<MuCS::MuCSData> > mucs = mucslist[0]; 
  
  //cout<<"mucsdata size:"<<mucs[0].size()<<"\n";
  
  //for now, only concerned with the first entry in mucsdata
  /*cout<<"hits1 size:"<<mucs->at(0).Hits1().size()<<"\n"; 
  cout<<"hits2 size:"<<mucs->at(0).Hits2().size()<<"\n"; 
  cout<<"hits3 size:"<<mucs->at(0).Hits3().size()<<"\n"; 
  cout<<"hits7 size:"<<mucs->at(0).Hits7().size()<<"\n"; */
  
  dbfields xfields, zfields;
  Float_t xq=0.,xq_rms=0.,x=0.,x_rms=0.,zq=0.,zq_rms=0.,z=0.,z_rms=0.,y=0.;
  Int_t xmatches=0,zmatches=0;
  //only fill fields if there are hits
  if (mucs->at(0).Hits1().size()>0 && mucs->at(0).Hits2().size()>0 && mucs->at(0).Hits3().size()>0 && mucs->at(0).Hits7().size()>0){
    execFill(mucs->at(0).Hits3(),mucs->at(0).Hits1(),"hitmap13x",xfields); //assumes that PMT 3 is above PMT 1 spatially
    execFill(mucs->at(0).Hits7(),mucs->at(0).Hits2(),"hitmap27z",zfields);//assumes that PMT 7 is above PMT 2 spatially
    xq=xfields.q;
    xq_rms=xfields.q_rms;
    x=xfields.p;
    x_rms=xfields.p_rms;
    zq=zfields.q;
    zq_rms=zfields.q_rms;
    z=zfields.p;
    z_rms=zfields.p_rms;
    xmatches=xfields.matches;
    zmatches=zfields.matches;
    y=(xmatches>0 && zmatches>0) ? fTopBoxy : 0.; //only populate if there were database entries for this field

  }

  //now create and populate the MuCSReco object
  std::unique_ptr< std::vector<MuCS::MuCSRecoData> > mucsrecocol(new std::vector<MuCS::MuCSRecoData>);
  MuCS::MuCSRecoData mucsrecoevt( xq, xq_rms, x, x_rms, zq, zq_rms, z, z_rms, y, xmatches, zmatches );
  mucsrecocol->push_back( mucsrecoevt );
  evt.put( std::move( mucsrecocol ) );
  
}

DEFINE_ART_MODULE( MuCSReco )

#endif

////////////////////////////////////////////////////////////////////////
//
//    The end !
//
////////////////////////////////////////////////////////////////////////

/** ****************************************************************************
 * @file ArgoUrl.h
 * @brief Helper to create a URL for Argo event viewer
 * @author ntagg@otterbein.edu
 * @see  ArgoUrl.cpp
 *  
 * This little helper class will construct a well-formed URL that will jump to an event in Argo, 
 * including coordinate details to zoom in on a region of interest.
 * 
 * Syntax is 'method chaining'.  Order of operations is probably not important, but later operations take precidence.
 *
 * To use, simply include this file
 * #include "uboone/Utilities/ArgoUrl.h"
 *
 *
 *
 * Simple example: give URL to current event
 *     std::cout << "URL of current event is " << util::ArgoUrl() << std::endl();
 * 
 * Example: zoom in on a track vertex, store url in string, show me a region of +- 100 wires around that
 *      std::string url =  util::ArgoUrl().center(myTrk.Vertex()).wires(200);
 *
 * Example: zoom in on a specific wire and TDC count, with a total view width of at least 100 cm
 *      std::string url =  util::ArgoUrl().center(myWireID).center(tdc).width(100.);
 *
 * Example: Show an entire track.
 *      std::cout <<   util::ArgoUrl().frame(myTrk.Vertex(),myTrk.End());
 * 
 * By default, this helper class tries to figure out what input file you're looking at and which event in that file
 * by being devious and sneaky with ROOT internals.   
 * You can override this by giving explicit filenames and entries
 *
 * Contact: Nathaniel (ntagg@otterbein.edu)
 * ****************************************************************************/

#ifndef __ARGOURL_H_
#define __ARGOURL_H_

// For Argo thing
#include <sys/param.h>
#include <stdlib.h>
#include <string>
#include <map>
#include <sstream>
#include "TVector3.h"

#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>


// Framework includes
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "larcore/Geometry/Geometry.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"


namespace util {

class ArgoUrl {
  public:
    ArgoUrl();


    ArgoUrl& center(int tdc);                         ///< Center on this TDC.
    ArgoUrl& t1(int t);                               ///< Start TDC value (bottom side of view range)
    ArgoUrl& t2(int t);                               ///< End TDC value (bottom side of view range)

    ArgoUrl& wires(int w);                            ///< Width of view region in number of wires.         
    ArgoUrl& width(double w);                         ///< Select a region w cm wide

    ArgoUrl& center(const geo::WireID& w);             ///< Center on this wire (in the view of this wire)
    ArgoUrl& center(int plane, int wire);              ///< Same as call above, without casts.

    ArgoUrl& center(const TVector3& p);                ///< Center on this XYZ coordinate
    ArgoUrl& frame(const TVector3& p1, const TVector3& p2); ///< Try to frame both these coordinates

    ArgoUrl& base(const std::string& s) ; ///< In case you have your own Argo isntallation, overrride here.

    ArgoUrl& file(std::string filename);              ///< File location.
    ArgoUrl& entry(int entry);                        ///< event number within the fie


    operator std::string();                             ///< Deliver unto me the URL.

  protected:
    std::string _base;
    std::map<std::string,std::string> _hash;
    const int koffset = 2400; // Set for doing 2016 analysis, where the hit times have been shifted relative to what Argo expects.

  };

// Implmentation follows.
// This could all go in a .cxx library somehwere, but at present there exists no general Utility library, so screw it. Make it all inline.

inline ArgoUrl::ArgoUrl() :
_base("http://argo-microboone.fnal.gov/")
{  }




inline ArgoUrl& ArgoUrl::base(const std::string& s) // In case you have your own Argo isntallation, overrride here.
{
    _base = s; return *this;
}

 inline ArgoUrl& ArgoUrl::file(std::string filename)  
 {
    // Override the filename.
    _hash["filename"] = filename;
    return *this;
 }

 inline ArgoUrl& ArgoUrl::entry(int entry)
 {
    _hash["entry"] = std::to_string(entry);
    return *this;
 }


inline ArgoUrl& ArgoUrl::center(int tdc)           // Center on this TDC.
{ 
  tdc+=koffset; // temporary hack for timing region shift. 
  t1(tdc-100); 
  t2(tdc+100); 
  return *this; 
}

inline ArgoUrl& ArgoUrl::t1(int t)                     // Start TDC value (bottom side of view range)
{ t+=koffset; _hash["t1"] = std::to_string(t); return *this;}


inline ArgoUrl& ArgoUrl::t2(int t)                     // End TDC value (bottom side of view range)
{ t+=koffset; _hash["t2"] = std::to_string(t); return *this;}

inline ArgoUrl& ArgoUrl::wires(int w)              // Width of view region in number of wires.     
{ _hash["wires"] = std::to_string(w); return *this; }

inline ArgoUrl& ArgoUrl::width(double w)           // Select a region w cm wide
{ 
    _hash["wires"] = w/0.3; 
    return *this;
}

inline ArgoUrl& ArgoUrl::center(const geo::WireID& w)
{
   center(w.Plane,w.Wire); 
   return *this;
}

 inline ArgoUrl& ArgoUrl::center(int plane, int wire) // wire in this case is wire number on that plane, starting from 0 
{ 
    std::string key = "plane" + std::to_string(plane); 
    _hash[key] = std::to_string(wire); 
    return *this;
}

 inline ArgoUrl& ArgoUrl::center(const TVector3& p)
 {
    // Center on this XYZ coordinate
    auto const* detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();
    art::ServiceHandle<geo::Geometry> geo;



    // Loop planes.
    for (geo::PlaneID const& pID: geo->IteratePlaneIDs()) {
        try {       // Sometimes this throws out-of-bounds exceptions. We don't care.
            geo::WireID wire = geo->NearestWireID(p,pID);
            center(wire);
        } catch (...) {}; 

        int tdc = detprop->ConvertXToTicks(p.x(),pID);
        center(tdc-koffset); // do it all three planes, what the heck. Remove offset, since it doesn't apply to 3d coordinates.
     }
     return *this;
 }  
 
inline  ArgoUrl& ArgoUrl::frame(const TVector3& p1, const TVector3& p2) // Try to frame both these coordinates
{
    // Find the largest coordinate distance.
    auto const* detprop = lar::providerFrom<detinfo::DetectorPropertiesService>();

    TVector3 v= (p1+p2)*0.5;
    center(v);

    double dx = fabs(p1.x()-p2.x());
    double dy = fabs(p1.y()-p2.y());
    double dz = fabs(p1.z()-p2.z());
    if(dy>dx) dx = dy;
    if(dz>dx) dx = dz;
    width(dx);

    int ta = detprop->ConvertXToTicks(p1.x(),geo::PlaneID(0,0,2));
    int tb = detprop->ConvertXToTicks(p2.x(),geo::PlaneID(0,0,2));
    if(tb<ta) { double tt=ta; ta=tb; tb=tt; }
    _hash["t1"] = std::to_string(ta);
    _hash["t2"] = std::to_string(tb);


    return *this;
}



inline ArgoUrl::operator std::string()
{
    // Do EEET

    // We need a valid filename. Go find one if it hasn't been provided.  Ditto entry number.

    // FIXME: 
    // This currently looks for the input file, which is correct when user is only doing an analysis task.
    // It's NOT correct for an output file, where the user is doing reconstruction and wants to see her reco results.
    // Need a test bed to fix it, but it should work by looking for any file with an "Events" tree.

    if(_hash["filename"] == "" && _hash["entry"] == "") {

        // Try to get current file name.
        std::string filename="CHANGE_ME_TO_FILE_NAME.root";
        Long64_t entry =0;

        if(gROOT){
          TSeqCollection* filelist = gROOT->GetListOfFiles();
          for(int i=0;i<filelist->GetEntries();i++) {
            TFile* f = dynamic_cast<TFile*>(filelist->At(i));
            if(f) {
              // We have a file.  See if it's a READ-only.
              if(strcmp(f->GetOption(),"READ")==0) {
                // Look through for an 'Events' tree.
                TIter it2(f->GetList());
                TObject* o2;
                while((o2 = it2())) {
                  TTree* t = dynamic_cast<TTree*>(o2);
                  if(t) {
                    if(strcmp(t->GetName(),"Events")==0) {
                      filename = f->GetName();
                      entry = t->GetReadEntry(); // Get the current entry.
                    }
                  }
                }
              }
            }
          }
        }

        // Attempt to find fully-qualified pathname.
        char fullname[PATH_MAX];
        if(::realpath(filename.c_str(),fullname) != NULL) {
          filename = fullname;
        }
        if(_hash["filename"] == "") _hash["filename"] = filename;
        if(_hash["entry"] == "")    _hash["entry"] = std::to_string(entry);
    }

    if(_hash["t1"] == "") _hash["t1"] = "0"; // This element must be there to get coordinate info.

    std::string url = _base;
    for(std::map<std::string,std::string>::iterator it = _hash.begin(); it!= _hash.end(); it++){
      if(it==_hash.begin()) url.append("#");
      else                  url.append("&");
      url.append(it->first);
      url.append("=");
      url.append(it->second);
    }
    return url;
}

}
// Wrapper so it streams correctly as a string.
inline std::ostream & operator <<(std::ostream & os, util::ArgoUrl & t) { os << std::string(t); return os;}


#endif
////////////////////////////////////////////////////////////////////////
// Class:       AnodeCathodePMAlg
//
// Algorithm for quickly identifying anode-to-cathode tracks based on a
// simple pattern-matching approach.
//
// Input:  vector<recob::Hit>
// Output: true if anode-cathode track candidate in event (false if not)
// 
////////////////////////////////////////////////////////////////////////

#include <vector>
#include <unordered_set>
//#include <boost/functional/hash.hpp>
#include <tuple>
#include <chrono>

#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"

//#include "fhiclcpp/ParameterSet.h"
//#include "larcore/Geometry/GeometryCore.h"
//#include "lardata/DetectorInfo/DetectorProperties.h"

namespace pm      { class AnodeCathodePMAlg; }
namespace geo     { class GeometryCore; }
namespace detinfo { class DetectorProperties; }
//namespace fhicl   { class ParameterSet; }
namespace recob   { class Hit; }

class pm::AnodeCathodePMAlg{

  typedef std::pair<int,int> Indices_t;
  //typedef std::unordered_set< Indices_t, boost::hash<Indices_t> > Pattern_t;
  typedef std::vector< Indices_t > Pattern_t;
  
 public:
  AnodeCathodePMAlg();

  void Configure(//fhicl::ParameterSet         const&,
		 geo::GeometryCore           const&,
		 detinfo::DetectorProperties const&);

  void RunPatternMatching(std::vector<recob::Hit> const&,float &);
  
 private:

  void          ClearHitmap();
  unsigned int  GetBinTime(float);
  unsigned int  GetBinWire(float);
  unsigned int  GetBinWire(geo::WireID wid)
  { return GetBinWire((float)(wid.Wire)); }

  void    MakePatterns();
  void    FillHitmap( std::vector<recob::Hit> const&);
  float   GetFractionMatched();
  
  std::vector<Pattern_t> fPatterns;
  std::vector< std::vector<float> > fHitmap;

  geo::PlaneID    fPlaneID;
  
  unsigned int    fWiresPerBin;
  unsigned int    fNwires;
  unsigned int    fNbinsWires;
  geo::WireID     fStartWire;
  geo::WireID     fEndWire;
  
  float        fTimeTicksPerBin;
  float        fStartTimeTick;
  float        fEndTimeTick;
  unsigned int fNbinsTimeTicks;
  unsigned int fNtimeTicksPattern;
  
  float  fIntegralThreshold;
  float  fMatchFractionThreshold;
};

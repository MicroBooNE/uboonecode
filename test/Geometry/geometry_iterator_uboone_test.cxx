/**
 * @file   geometry_iterator_uboone_test.cxx
 * @brief  Unit test for geometry iterators on MicroBooNE detector
 * @date   May 11th, 2015
 * @author petrillo@fnal.gov
 * 
 * Usage: just run the executable.
 * Boost unit testing environment keeps the arguments secret anyway.
 */

// Boost test libraries; defining this symbol tells boost somehow to generate
// a main() function; Boost is pulled in by geometry_boost_unit_test_base.h
#define BOOST_TEST_MODULE GeometryIteratorTestMicroBooNE

// LArSoft libraries
#include "test/Geometry/geometry_unit_test_uboone.h"
#include "test/Geometry/boost_unit_test_base.h"
#include "test/Geometry/GeometryIteratorTestAlg.h"
#include "larcore/Geometry/ChannelMapStandardAlg.h"

//------------------------------------------------------------------------------
//---  The test environment
//---

// we define here all the configuration that is needed;
// in the specific, the type of the channel mapping and a proper test name,
// used for output only; we use MicroBooNEGeometryEnvironmentConfiguration
// as base class, that is already configured to use MicroBooNE geometry.
// We wrap it in testing::BoostCommandLineConfiguration<> so that it can learn
// the configuration file name from the command line.
struct MicroBooNEGeometryConfiguration:
  public testing::BoostCommandLineConfiguration<
    uboone::testing::MicroBooNEGeometryEnvironmentConfiguration
      <geo::ChannelMapStandardAlg>
    >
{
  /// Constructor: overrides the application name; ignores command line
  MicroBooNEGeometryConfiguration()
    { SetApplicationName("GeometryIteratorUnitTest"); }
}; // class MicroBooNEGeometryConfiguration

/*
 * Our fixture is based on GeometryTesterFixture, configured with the object
 * above.
 * It provides:
 * - `Tester`, a configured instance of the test algorithm.
 */
class MicroBooNEGeometryIteratorTestFixture:
  private testing::GeometryTesterEnvironment<MicroBooNEGeometryConfiguration>
{
    public:
  geo::GeometryIteratorTestAlg Tester;
  
  /// Constructor: initialize the tester with the Geometry from base class
  MicroBooNEGeometryIteratorTestFixture(): Tester(TesterParameters())
    { Tester.Setup(*Geometry()); }

}; // class MicroBooNEGeometryIteratorTestFixture



//------------------------------------------------------------------------------
//---  The tests
//---

BOOST_FIXTURE_TEST_SUITE
  (GeometryIteratorsMicroBooNE, MicroBooNEGeometryIteratorTestFixture)
// BOOST_GLOBAL_FIXTURE(MicroBooNEGeometryIteratorTestFixture)


BOOST_AUTO_TEST_CASE( AllTests )
{
  Tester.Run();
} // BOOST_AUTO_TEST_CASE( AllTests )

/*
BOOST_AUTO_TEST_CASE( CryostatIteratorsTest )
{
  Tester.CryostatIteratorsTest();
} // BOOST_AUTO_TEST_CASE( CryostatIteratorsTest )



BOOST_AUTO_TEST_CASE( TPCIteratorsTest )
{
  Tester.TPCIteratorsTest();
} // BOOST_AUTO_TEST_CASE( TPCIteratorsTest )



BOOST_AUTO_TEST_CASE( PlaneIteratorsTest )
{
  Tester.PlaneIteratorsTest();
} // BOOST_AUTO_TEST_CASE( PlaneIteratorsTest )



BOOST_AUTO_TEST_CASE( WireIteratorsTest )
{
  Tester.WireIteratorsTest();
} // BOOST_AUTO_TEST_CASE( WireIteratorsTest )
*/

BOOST_AUTO_TEST_SUITE_END()


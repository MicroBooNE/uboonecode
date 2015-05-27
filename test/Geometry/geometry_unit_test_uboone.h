/**
 * @file   geometry_unit_test_uboone.h
 * @brief  Class for objects initializing MicroBooNE geometry
 * @date   May 11th, 2015
 * @author petrillo@fnal.gov
 * 
 * Provides an environment for easy set up of MicroBooNE-aware tests.
 * Keep in mind that the channel mapping algorithm must be hard-coded and, if
 * using Boost unit test, the configuration file location must be hard coded too
 * (or you can use the provided configuration).
 * 
 * For an example of usage, see larcore/test/Geometry/geometry_test.cxx
 */

#ifndef TEST_GEOMETRY_UNIT_TEST_UBOONE_H
#define TEST_GEOMETRY_UNIT_TEST_UBOONE_H

// LArSoft libraries
#include "test/Geometry/geometry_unit_test_base.h"
// #include "Geometry/ChannelMapStandardAlg.h"

// C/C++ standard libraries
#include <string>


namespace geo {
  class ChannelMapStandardAlg;
} // namespace geo

namespace uboone {
  
  /// Namespace including MicroBooNE-specific testing
  namespace testing {
    
    /** ************************************************************************
     * @brief Class holding the configuration for a MicroBooNE fixture
     * @tparam CHANNELMAP the class used for channel mapping
     * @see BasicGeometryEnvironmentConfiguration
     *
     * This class needs to be fully constructed by the default constructor
     * in order to be useful as Boost unit test fixture.
     * It is supposed to be passed as a template parameter to another class
     * that can store an instance of it and extract configuration information
     * from it.
     * 
     * This class should be used with ChannelMapStandardAlg.
     * 
     * We reuse BasicGeometryEnvironmentConfiguration as base class and then we
     * fix its setup.
     */
    template <typename CHANNELMAP = geo::ChannelMapStandardAlg>
    struct MicroBooNEGeometryEnvironmentConfiguration:
      public ::testing::BasicGeometryEnvironmentConfiguration<CHANNELMAP>
    {
      // remember that BasicGeometryEnvironmentConfiguration is not polymorphic
      using base_t
        = ::testing::BasicGeometryEnvironmentConfiguration<CHANNELMAP>;
      
      /// Default constructor
      MicroBooNEGeometryEnvironmentConfiguration() { MicroBooNEdefaultInit(); }
      
      /// Constructor; accepts the name as parameter
      MicroBooNEGeometryEnvironmentConfiguration(std::string name):
        MicroBooNEGeometryEnvironmentConfiguration()
        { base_t::SetApplicationName(name); }
      
      
        private:
      void MicroBooNEdefaultInit()
        {
          // overwrite the configuration that happened in the base class:
          base_t::SetApplicationName("MicroBooNEGeometryTest");
          base_t::SetDefaultGeometryConfiguration(R"(
              services: {
                Geometry: {
                  SurfaceY: 6.9e2 # in cm, vertical distance to the surface
                  Name:     "microboonevX"
                  GDML:     "microboonevX.gdml"
                  ROOT:     "microboonevX.gdml"
                  SortingParameters: {}
                } # Geometry
              } # services
            )");
        }
    }; // class MicroBooNEGeometryEnvironmentConfiguration<>
    
    
  } // namespace testing
} // namespace uboone

#endif // TEST_GEOMETRY_UNIT_TEST_UBOONE_H

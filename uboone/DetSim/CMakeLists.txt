art_make( BASENAME_ONLY
          EXCLUDE SimWireMicroBooNE_module.cc
          MODULE_LIBRARIES larcore_Geometry
                        larcore_Geometry_Geometry_service
                        larsim_Simulation
                        lardata_Utilities
                        larevt_Filters
                        lardata_RawData
                        SignalShapingServiceMicroBooNE_service
			larevt_CalibrationDBI_Providers
                        larsim_RandomUtils_LArSeedService_service
                        ${ART_FRAMEWORK_CORE}
                        ${ART_FRAMEWORK_PRINCIPAL}
                        ${ART_FRAMEWORK_SERVICES_REGISTRY}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
                        ${ART_PERSISTENCY_COMMON}
                        ${ART_PERSISTENCY_PROVENANCE}
                        ${ART_UTILITIES}
                        ${MF_MESSAGELOGGER}
                        ${MF_UTILITIES}
                        ${FHICLCPP}
                        ${CETLIB}
                        ${CLHEP}
                        ${ROOT_GEOM}
                        ${ROOT_XMLIO}
                        ${ROOT_GDML}
                        ${ROOT_BASIC_LIB_LIST}
        )

simple_plugin(SimWireMicroBooNE "module"
          BASENAME_ONLY
              larcore_Geometry
              larcore_Geometry_Geometry_service
              larsim_Simulation
              lardata_Utilities
              larevt_Filters
              lardata_RawData
              lardata_DetectorInfoServices_DetectorClocksServiceStandard_service # BUG note that we link to a specific service implementation
              SignalShapingServiceMicroBooNE_service
              larevt_CalibrationDBI_Providers
              larsim_RandomUtils_LArSeedService_service
              ${ART_FRAMEWORK_CORE}
              ${ART_FRAMEWORK_PRINCIPAL}
              ${ART_FRAMEWORK_SERVICES_REGISTRY}
              ${ART_FRAMEWORK_SERVICES_OPTIONAL}
              ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
              ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
              ${ART_PERSISTENCY_COMMON}
              ${ART_PERSISTENCY_PROVENANCE}
              ${ART_UTILITIES}
              ${MF_MESSAGELOGGER}
              ${MF_UTILITIES}
              ${FHICLCPP}
              ${CETLIB}
              ${CLHEP}
              ${ROOT_GEOM}
              ${ROOT_XMLIO}
              ${ROOT_GDML}
              ${ROOT_BASIC_LIB_LIST}
        )

# install_headers()
install_fhicl()
install_source()

art_make( EXCLUDE UBOpticalChConfig_service.cc OpticalDRAMReadout_module.cc
	  #NO_PLUGINS
	  BASENAME_ONLY
	  LIBRARY_NAME  uboonecode_uboone_OpticalDetectorSim
          LIB_LIBRARIES UBOpticalChConfig_service
			RawData
			Utilities
			Geometry
                        ${SIMULATIONBASE}
                        ${ART_FRAMEWORK_CORE}
			${ART_FRAMEWORK_PRINCIPAL}
			${ART_FRAMEWORK_SERVICES_REGISTRY}
			${ART_FRAMEWORK_SERVICES_BASIC}
			${ART_FRAMEWORK_SERVICES_OPTIONAL}
			${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
			${ART_PERSISTENCY_COMMON}
			${ART_UTILITIES}
			${MF_MESSAGELOGGER}
			${MF_UTILITIES}
	                ${FHICLCPP}
			${CETLIB}
                        ${ROOT_BASIC_LIB_LIST}
          MODULE_LIBRARIES Simulation
			OpticalDetectorData
			UBOpticalChConfig_service
			UBOpReadoutMap_service
			uboonecode_uboone_OpticalDetectorSim
			Utilities
			UBooNEGeometryHelper_service
			uboone_Geometry
                        ${SIMULATIONBASE}
                        ${ART_FRAMEWORK_CORE}
			${ART_FRAMEWORK_PRINCIPAL}
			${ART_FRAMEWORK_SERVICES_REGISTRY}
			${ART_FRAMEWORK_SERVICES_BASIC}
			${ART_FRAMEWORK_SERVICES_OPTIONAL}
			${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
			${ART_PERSISTENCY_COMMON}
			${ART_UTILITIES}
			${MF_MESSAGELOGGER}
			${MF_UTILITIES}
	                ${FHICLCPP}
			${CETLIB}
                        ${ROOT_BASIC_LIB_LIST}

         )

simple_plugin(UBOpticalChConfig "service"
	      uboone_Geometry
	      UBOpReadoutMap_service	
	      larcore_Geometry_Geometry_service
	      Geometry
#	      uboonecode_uboone_OpticalDetectorSim
	      ${ART_UTILITIES}
	      ${MF_MESSAGELOGGER}
	      ${MF_UTILITIES}
	      ${FHICLCPP}
	      ${CETLIB}
	      ${ROOT_BASIC_LIB_LIST}
              BASENAME_ONLY
)

simple_plugin(OpticalDRAMReadout "module"
                  Simulation
                  lardata_DetectorInfoServices_DetectorClocksServiceStandard_service # BUG note that we link to a specific service implementation
                  OpticalDetectorData
                  UBOpticalChConfig_service
                  UBOpReadoutMap_service
                  uboonecode_uboone_OpticalDetectorSim
                  Utilities
                  UBooNEGeometryHelper_service
                  uboone_Geometry
                  ${SIMULATIONBASE}
                  ${ART_FRAMEWORK_CORE}
                  ${ART_FRAMEWORK_PRINCIPAL}
                  ${ART_FRAMEWORK_SERVICES_REGISTRY}
                  ${ART_FRAMEWORK_SERVICES_BASIC}
                  ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                  ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
                  ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
                  ${ART_PERSISTENCY_COMMON}
                  ${ART_UTILITIES}
                  ${MF_MESSAGELOGGER}
                  ${MF_UTILITIES}
                  ${FHICLCPP}
                  ${CETLIB}
                  ${ROOT_BASIC_LIB_LIST}
              BASENAME_ONLY
)

install_headers()
install_fhicl()
install_source()

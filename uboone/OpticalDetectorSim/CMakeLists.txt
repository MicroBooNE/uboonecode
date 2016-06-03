art_make( EXCLUDE UBOpticalChConfig_service.cc OpticalDRAMReadout_module.cc
	  #NO_PLUGINS
	  BASENAME_ONLY
	  LIBRARY_NAME  uboonecode_uboone_OpticalDetectorSim
          LIB_LIBRARIES UBOpticalChConfig_service
			lardataobj_RawData
			lardata_Utilities
			larcore_Geometry
                        nusimdata_SimulationBase
                        ${ART_FRAMEWORK_CORE}
			${ART_FRAMEWORK_PRINCIPAL}
			${ART_FRAMEWORK_SERVICES_REGISTRY}
			${ART_FRAMEWORK_SERVICES_BASIC}
			${ART_FRAMEWORK_SERVICES_OPTIONAL}
			${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
			art_Persistency_Common canvas_Persistency_Common
			art_Utilities canvas_Utilities
			${MF_MESSAGELOGGER}
			${MF_UTILITIES}
	                ${FHICLCPP}
			${CETLIB}
                        ${ROOT_BASIC_LIB_LIST}
          MODULE_LIBRARIES larsim_Simulation larsimobj_Simulation
			lardataobj_OpticalDetectorData
			UBOpticalChConfig_service
			UBOpReadoutMap_service
			uboonecode_uboone_OpticalDetectorSim
			lardata_Utilities
			UBooNEGeometryHelper_service
			uboone_Geometry
                        nusimdata_SimulationBase
                        ${ART_FRAMEWORK_CORE}
			${ART_FRAMEWORK_PRINCIPAL}
			${ART_FRAMEWORK_SERVICES_REGISTRY}
			${ART_FRAMEWORK_SERVICES_BASIC}
			${ART_FRAMEWORK_SERVICES_OPTIONAL}
			${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
			art_Persistency_Common canvas_Persistency_Common
			art_Utilities canvas_Utilities
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
	      larcore_Geometry
#	      uboonecode_uboone_OpticalDetectorSim
	      art_Utilities canvas_Utilities
	      ${MF_MESSAGELOGGER}
	      ${MF_UTILITIES}
	      ${FHICLCPP}
	      ${CETLIB}
	      ${ROOT_BASIC_LIB_LIST}
              BASENAME_ONLY
)

simple_plugin(OpticalDRAMReadout "module"
                  larsim_Simulation larsimobj_Simulation
                  lardata_DetectorInfoServices_DetectorClocksServiceStandard_service # BUG note that we link to a specific service implementation
                  lardataobj_OpticalDetectorData
                  UBOpticalChConfig_service
                  UBOpReadoutMap_service
                  uboonecode_uboone_OpticalDetectorSim
                  lardata_Utilities
                  UBooNEGeometryHelper_service
                  uboone_Geometry
                  nusimdata_SimulationBase
                  ${ART_FRAMEWORK_CORE}
                  ${ART_FRAMEWORK_PRINCIPAL}
                  ${ART_FRAMEWORK_SERVICES_REGISTRY}
                  ${ART_FRAMEWORK_SERVICES_BASIC}
                  ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                  ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
                  ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
                  art_Persistency_Common canvas_Persistency_Common
                  art_Utilities canvas_Utilities
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

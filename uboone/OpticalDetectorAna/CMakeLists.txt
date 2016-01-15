
art_make( BASENAME_ONLY
	  LIBRARY_NAME  uboone_OpticalDetectorAna
	  LIB_LIBRARIES
			larana_OpticalDetector
                        ${ROOT_BASIC_LIB_LIST}

	  MODULE_LIBRARIES larcore_Geometry
	  	        larcore_Geometry_Geometry_service
		        lardata_Utilities
	                uboone_Geometry
                        lardata_OpticalDetectorData
			larana_OpticalDetector
                        lardata_RawData
			lardata_RecoBase
			larsim_Simulation
			SimulationBase
			larsim_PhotonPropagation
			uboone_OpticalDetectorAna
                        ${ART_FRAMEWORK_CORE}
			${ART_FRAMEWORK_IO_SOURCES}
                        ${ART_FRAMEWORK_PRINCIPAL}
                        ${ART_FRAMEWORK_SERVICES_REGISTRY}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
                        ${ART_PERSISTENCY_COMMON}
                        ${ART_PERSISTENCY_PROVENANCE}
                        ${ART_UTILITIES}
                        ${MF_MESSAGELOGGER}
                        ${MF_UTILITIES}
                        ${FHICLCPP}
                        ${CETLIB}
                        ${ROOT_GEOM}
                        ${ROOT_XMLIO}
                        ${ROOT_GDML}
                        ${ROOT_BASIC_LIB_LIST}
			${SIMULATIONBASE}
        )							   

install_headers()
install_source()
install_fhicl()
add_subdirectory(job)
add_subdirectory(OpticalSubEvents)


add_subdirectory(Optical)

art_make( BASENAME_ONLY
	  LIBRARY_NAME  uboonecode_uboone_Calibrations
	  LIB_LIBRARIES
                        LArFFT_service
			${ART_FRAMEWORK_CORE}
			${ART_UTILITIES}
			${MF_MESSAGELOGGER}
			${MF_UTILITIES}
			${CETLIB}
      ${Boost_SYSTEM_LIBRARY}
			${ROOT_BASIC_LIB_LIST}
          MODULE_LIBRARIES 
			uboonecode_uboone_Calibrations
			RawData
			CalData
			Geometry
			Geometry_service
			BackTracker_service
			Simulation
			RecoBase
			RecoObjects  
			RecoAlg
			Utilities
			Filters
                        LArFFT_service
			SignalShapingServiceMicroBooNE_service
			${SIMULATIONBASE}
			${LARRECO_LIB}
			${LARDATA_LIB}
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

install_headers()
install_fhicl()
install_source()


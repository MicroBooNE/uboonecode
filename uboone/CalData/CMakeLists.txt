art_make( BASENAME_ONLY
          LIBRARY_NAME CalData
          LIB_LIBRARIES 
			lardata_RawData
                        lardata_RecoBase
	                lardata_Utilities
	                ${FHICLCPP}
			${CETLIB}
			${ROOT_BASIC_LIB_LIST}
	  MODULE_LIBRARIES larcore_Geometry
	                CalData
  	  		NoiseFilterAlgs
	  		larcore_Geometry_Geometry_service
	                lardata_Utilities
			lardata_Utilities_DetectorProperties_service
			larevt_Filters
			lardata_RawData
			larevt_CalibrationDBI_IOVData
			larevt_CalibrationDBI_Providers
			lardata_RecoBase
			lardata_RecoBaseArt
			SignalShapingServiceMicroBooNE_service
                        ${ART_FRAMEWORK_CORE}
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
			${ROOT_FFTW}
			${ROOT_BASIC_LIB_LIST}
        )

add_subdirectory(NoiseFilterAlgs)

install_headers()
install_fhicl()
install_source()

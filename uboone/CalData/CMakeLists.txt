art_make( BASENAME_ONLY
          LIBRARY_NAME CalData
          LIB_LIBRARIES 
			RawData
                        RecoBase
	                Utilities
	                ${FHICLCPP}
			${CETLIB}
			${ROOT_BASIC_LIB_LIST}
	  MODULE_LIBRARIES Geometry
  	  		CalData
	  		Geometry_service
	                Utilities
			Filters
			RawData
			CalibrationDBI_IOVData
			CalibrationDBI_WebDBI
			RecoBase
			RecoBaseArt
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
			${ROOT_BASIC_LIB_LIST}
        )

install_headers()
install_fhicl()
install_source()


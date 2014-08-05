

set( CalibrationTPC_Algs_sources CalibrationTPC_Algs.cc)
add_library(CalibrationTPC_Algs SHARED ${CalibrationTPC_Algs_sources})
install ( TARGETS CalibrationTPC_Algs DESTINATION ${flavorqual_dir}/lib )

art_make( BASENAME_ONLY
	  LIBRARY_NAME  uboonecode_uboone_Calibrations
          MODULE_LIBRARIES 
	  		CalibrationTPC_Algs
			RawData
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


cet_enable_asserts()

set( uboone_util_lib_list LArFFT_service
                        LArProperties_service
			Geometry
			Geometry_service
			Utilities
                        SignalShapingServiceMicroBooNE_service
  			${ART_FRAMEWORK_CORE}
			${ART_FRAMEWORK_PRINCIPAL}
			${ART_FRAMEWORK_SERVICES_REGISTRY}
			${ART_PERSISTENCY_COMMON}
			${ART_PERSISTENCY_PROVENANCE}
			${ART_UTILITIES}
			${MF_MESSAGELOGGER}
			${MF_UTILITIES}
			${CETLIB}
			${ROOT_BASIC_LIB_LIST}
    )


simple_plugin( SignalShapingServiceMicroBooNE  "service"
               ${uboone_util_lib_list}
	       BASENAME_ONLY
        )
             
 
install_headers()
install_fhicl()
install_source()

add_subdirectory(test)


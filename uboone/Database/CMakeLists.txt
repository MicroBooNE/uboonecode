art_make( BASENAME_ONLY
          LIBRARY_NAME Database
          LIB_LIBRARIES 
	                ${ART_FRAMEWORK_CORE}
			${ART_FRAMEWORK_PRINCIPAL}
			${ART_FRAMEWORK_SERVICES_REGISTRY}
			${ART_UTILITIES}
			${FHICLCPP}
			${CETLIB}
			${ROOT_BASIC_LIB_LIST}
                        Geometry
			Geometry_service
        )

install_headers()
install_fhicl()
install_source()


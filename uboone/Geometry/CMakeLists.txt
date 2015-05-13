
art_make( BASENAME_ONLY
	  LIBRARY_NAME uboone_Geometry
	  LIB_LIBRARIES Geometry
                        ${MF_MESSAGELOGGER}
                        ${CETLIB}
                        ${FHICLCPP}
                        ${ROOT_BASIC_LIB_LIST}
	  SERVICE_LIBRARIES Geometry
	  		uboone_Geometry
			${MF_UTILITIES}
	                ${FHICLCPP}
			${CETLIB}
                        ${ROOT_GEOM}
                        ${ROOT_XMLIO}
                        ${ROOT_GDML}
			${ROOT_BASIC_LIB_LIST}
        )							   

add_subdirectory(gdml)

install_headers()
install_fhicl()
install_source()



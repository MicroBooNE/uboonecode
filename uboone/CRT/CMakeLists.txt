include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )
include_directories( $ENV{NUTOOLS_INC} )

art_make( BASENAME_ONLY
	  EXCLUDE 
	  	  CRTChannelMapAlg.cc 
		  CRTGeometryHelper_service.cc
		  CRTDetSim_module.cc
          LIBRARY_NAME  
	  		CRTData
          LIB_LIBRARIES
                        ${ROOT_BASIC_LIB_LIST}
                        ${ROOT_EG}
        )


art_make_library( LIBRARY_NAME CRT
                  SOURCE
                            CRTChannelMapAlg.cc
                  LIBRARIES larcore_Geometry
                            ${MF_MESSAGELOGGER}
                            ${MF_UTILITIES}
                            ${FHICLCPP}
                            ${ROOT_BASIC_LIB_LIST}
                            ${ROOT_GEOM}
                            ${CETLIB}
)

simple_plugin( CRTGeometryHelper service
                larcore_Geometry
                CRT
                ${ART_FRAMEWORK_CORE}
                ${ART_FRAMEWORK_IO_SOURCES}
                ${ART_FRAMEWORK_PRINCIPAL}
                art_Persistency_Common canvas_Persistency_Common
                art_Persistency_Provenance canvas_Persistency_Provenance
                art_Utilities canvas_Utilities
                ${ART_FRAMEWORK_SERVICES_REGISTRY}
                ${MF_MESSAGELOGGER}
                ${MF_UTILITIES}
                ${ROOT_BASIC_LIB_LIST}
                ${CETLIB}
                BASENAME_ONLY
)
simple_plugin( CRTDetSim module
                larcore_Geometry
                CRT
		CRTData
                ${ART_FRAMEWORK_CORE}
                ${ART_FRAMEWORK_IO_SOURCES}
                ${ART_FRAMEWORK_PRINCIPAL}
                art_Persistency_Common canvas_Persistency_Common
                art_Persistency_Provenance canvas_Persistency_Provenance
                art_Utilities canvas_Utilities
                ${ART_FRAMEWORK_SERVICES_REGISTRY}
                ${MF_MESSAGELOGGER}
                ${MF_UTILITIES}
                ${ROOT_BASIC_LIB_LIST}
                ${CETLIB}
                BASENAME_ONLY
)

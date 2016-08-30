
#
simple_plugin( LArRawInputSourceUBooNE
               "source"
	       uboone_RawData
	       UBOpReadoutMap_service
               ${ART_FRAMEWORK_PERSISTENCY_COMMON}
               BASENAME_ONLY
             )

simple_plugin( RawDigitWriter "module"
	       UBOpReadoutMap_service
	       lardata_Utilities_DatabaseUtil_service
	       lardataobj_RawData
               ${ART_FRAMEWORK_CORE}
	       ${ART_FRAMEWORK_IO_SOURCES}
               ${ART_FRAMEWORK_PRINCIPAL}
               ${ART_FRAMEWORK_SERVICES_REGISTRY}
               ${ART_FRAMEWORK_SERVICES_OPTIONAL}
               ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
               art_Persistency_Common canvas_Persistency_Common
               art_Persistency_Provenance canvas_Persistency_Provenance
               art_Utilities canvas_Utilities
	       ${ROOT_BASIC_LIB_LIST}	       
               BASENAME_ONLY
             )

# need a special link list because of the explicit link to a specific service inplementation
simple_plugin( TriggerEnforcement "module"
                   lardata_DetectorInfoServices_DetectorClocksServiceStandard_service # BUG note that we link to a specific service implementation
                   lardataobj_RawData
                   uboone_RawData
                   ${ART_FRAMEWORK_CORE}
                   ${ART_FRAMEWORK_PRINCIPAL}
                   ${ART_FRAMEWORK_SERVICES_REGISTRY}
                   ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                   art_Persistency_Common canvas_Persistency_Common
                   art_Persistency_Provenance canvas_Persistency_Provenance
                   art_Utilities canvas_Utilities
                   ${MF_MESSAGELOGGER}
                   ${MF_UTILITIES}
                   ${FHICLCPP}
                   ${CETLIB}
                   ${ROOT_BASIC_LIB_LIST}
                   BASENAME_ONLY
             )


art_make( EXCLUDE LArRawInputSourceUBooNE_source.cc RawDigitWriter_module.cc TriggerEnforcement_module.cc
	  BASENAME_ONLY
	  LIBRARY_NAME uboone_RawData
	  LIB_LIBRARIES lardata_Utilities_DatabaseUtil_service
                        larcore_Geometry
	  	        larcore_Geometry_Geometry_service
		        lardata_Utilities
	                uboone_Geometry
			UBOpReadoutMap_service
                        lardataobj_OpticalDetectorData
                        lardataobj_RawData
                        larcoreobj_SummaryData
			${UBDut}
                        ${PQ}
			${Boost_SERIALIZATION_LIBRARY}
                        ${ART_FRAMEWORK_CORE}
			${ART_FRAMEWORK_IO_SOURCES}
                        ${ART_FRAMEWORK_PRINCIPAL}
                        ${ART_FRAMEWORK_SERVICES_REGISTRY}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
                        art_Persistency_Common canvas_Persistency_Common
                        art_Persistency_Provenance canvas_Persistency_Provenance
                        art_Utilities canvas_Utilities
                        ${MF_MESSAGELOGGER}
                        ${MF_UTILITIES}
                        ${FHICLCPP}
                        ${CETLIB}
                        ${ROOT_GEOM}
                        ${ROOT_XMLIO}
                        ${ROOT_GDML}
                        ${Boost_SYSTEM_LIBRARY}
                        ${ROOT_BASIC_LIB_LIST}

	  MODULE_LIBRARIES lardata_Utilities_DatabaseUtil_service
                        larcore_Geometry
	  	        larcore_Geometry_Geometry_service
		        lardata_Utilities
	                uboone_Geometry
			UBOpReadoutMap_service
                        lardataobj_OpticalDetectorData
                        lardataobj_RawData
                        larcoreobj_SummaryData
			uboone_RawData
			${UBDut}
                        ${PQ}
			${Boost_SERIALIZATION_LIBRARY}
                        ${ART_FRAMEWORK_CORE}
			${ART_FRAMEWORK_IO_SOURCES}
                        ${ART_FRAMEWORK_PRINCIPAL}
                        ${ART_FRAMEWORK_SERVICES_REGISTRY}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
                        art_Persistency_Common canvas_Persistency_Common
                        art_Persistency_Provenance canvas_Persistency_Provenance
                        art_Utilities canvas_Utilities
                        ${MF_MESSAGELOGGER}
                        ${MF_UTILITIES}
                        ${FHICLCPP}
                        ${CETLIB}
                        ${ROOT_GEOM}
                        ${ROOT_XMLIO}
                        ${ROOT_GDML}
                        ${Boost_SYSTEM_LIBRARY}
                        ${ROOT_BASIC_LIB_LIST}
        )							   

install_headers()
install_source()

add_subdirectory(job)

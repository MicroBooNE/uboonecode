
# find_ups_product is meant to be called from the top CMakeLists.txt file
find_ups_product( uboonedaq_datatypes v6_14_00 )

find_library( UBDut NAMES ubdata_types PATHS $ENV{UBOONEDAQ_DATATYPES_LIB} )

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
	       lardata_RawData
               ${ART_FRAMEWORK_CORE}
	       ${ART_FRAMEWORK_IO_SOURCES}
               ${ART_FRAMEWORK_PRINCIPAL}
               ${ART_FRAMEWORK_SERVICES_REGISTRY}
               ${ART_FRAMEWORK_SERVICES_OPTIONAL}
               ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
               ${ART_PERSISTENCY_COMMON}
               ${ART_PERSISTENCY_PROVENANCE}
               ${ART_UTILITIES}
	       ${ROOT_BASIC_LIB_LIST}	       
               BASENAME_ONLY
             )

# need a special link list because of the explicit link to a specific service inplementation
simple_plugin( TriggerEnforcement "module"
                   lardata_DetectorInfoServices_DetectorClocksServiceStandard_service # BUG note that we link to a specific service implementation
                   lardata_RawData
                   uboone_RawData
                   ${ART_FRAMEWORK_CORE}
                   ${ART_FRAMEWORK_PRINCIPAL}
                   ${ART_FRAMEWORK_SERVICES_REGISTRY}
                   ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                   ${ART_PERSISTENCY_COMMON}
                   ${ART_PERSISTENCY_PROVENANCE}
                   ${ART_UTILITIES}
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
                        lardata_OpticalDetectorData
                        lardata_RawData
                        larcore_SummaryData
			${UBDut}
                        ${PQ}
			${Boost_SERIALIZATION_LIBRARY}
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
                        ${Boost_SYSTEM_LIBRARY}
                        ${ROOT_BASIC_LIB_LIST}

	  MODULE_LIBRARIES lardata_Utilities_DatabaseUtil_service
                        larcore_Geometry
	  	        larcore_Geometry_Geometry_service
		        lardata_Utilities
	                uboone_Geometry
			UBOpReadoutMap_service
                        lardata_OpticalDetectorData
                        lardata_RawData
                        larcore_SummaryData
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
                        ${Boost_SYSTEM_LIBRARY}
                        ${ROOT_BASIC_LIB_LIST}
        )							   

install_headers()
install_source()

add_subdirectory(job)

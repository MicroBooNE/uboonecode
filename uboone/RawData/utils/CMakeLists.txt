
# find_ups_product is meant to be called from the top CMakeLists.txt file
find_ups_product( uboonedaq_datatypes v6_10_00 )

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
	       DatabaseUtil_service
	       RawData
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
 
               
art_make( EXCLUDE LArRawInputSourceUBooNE_source.cc RawDigitWriter_module.cc
	  BASENAME_ONLY
	  LIBRARY_NAME uboone_RawData
	  LIB_LIBRARIES DatabaseUtil_service
                  Geometry
	  	            Geometry_service
		              Utilities
	       	        TimeService_service
	                uboone_Geometry
			UBOpReadoutMap_service
                        OpticalDetectorData
                        RawData
                        SummaryData
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

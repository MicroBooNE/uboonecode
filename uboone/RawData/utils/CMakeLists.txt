
include_directories( $ENV{POSTGRESQL_INC} )
include_directories( $ENV{BOOST_INC} )
include_directories( $ENV{UBOONEDAQ_DATATYPES_INC} )

find_ups_product( postgresql v9_1_5 )
find_ups_product( uboonedaq_datatypes v6_08_00 )
find_ups_product( boost v1_57_0 )

cet_find_library( PQ NAMES pq PATHS ENV POSTGRESQL_LIBRARIES NO_DEFAULT_PATH )
find_library( UBDut NAMES ubdata_types PATHS $ENV{UBOONEDAQ_DATATYPES_LIB} )

cet_find_library( BOOSER NAMES boost_serialization PATHS ENV BOOST_LIB NO_DEFAULT_PATH )

simple_plugin( LArRawInputSourceUBooNE
               "source"
	       uboone_RawData
	       UBOpReadoutMap_service
               ${ART_FRAMEWORK_PERSISTENCY_COMMON}
               BASENAME_ONLY
             )
 
               
art_make( EXCLUDE LArRawInputSourceUBooNE_source.cc
	  BASENAME_ONLY
	  LIBRARY_NAME uboone_RawData
	  LIB_LIBRARIES Geometry
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
			${BOOSER}
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
                        ${ROOT_BASIC_LIB_LIST}
        )							   

install_headers()
install_source()

add_subdirectory(job)

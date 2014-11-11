
include_directories( $ENV{POSTGRESQL_INC} )
include_directories( $ENV{BOOST_INC} )
include_directories( $ENV{UBOONEDAQ_DATATYPES_INC} )

##<<<<<<< HEAD
find_ups_product( postgresql v9_1_5 )
find_ups_product( uboonedaq_datatypes v5_00_01 )
find_ups_product( boost v1_56_0 )

cet_find_library( PQ NAMES pq PATHS ENV POSTGRESQL_LIBRARIES NO_DEFAULT_PATH )
cet_find_library( UBDut NAMES udata_types PATHS ENV UBOONEDAQ_DATATYPES_LIB NO_DEFAULT_PATH )
cet_find_library( UBDuut NAMES udata_utilities PATHS ENV UBOONEDAQ_DATATYPES_LIB NO_DEFAULT_PATH )
cet_find_library( BOOSER NAMES boost_serialization PATHS ENV BOOST_LIB NO_DEFAULT_PATH )

##=======
##>>>>>>> v03_04_00
#simple_plugin( LArRawInputSourceUBooNE
#               "source"
#               Geometry
#               uboone_RawData
#               RawData
#               ${ART_FRAMEWORK_PERSISTENCY_COMMON}
#               BASENAME_ONLY
#             )
 
               
art_make( BASENAME_ONLY
	  LIBRARY_NAME uboone_RawData
	  LIB_LIBRARIES Geometry
                        OpticalDetectorData
                        RawData
                        SummaryData
			${UBDut}
			${UBDuut}
                        ${PQ}
			${BOOSER}
#                        daq_datatypes
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

	  SOURCE_LIBRARIES Geometry
	                OpticalDetectorData
			RawData
			uboone_RawData
			SummaryData
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

	  SERVICE_LIBRARIES Geometry
	  		lbne_Geometry
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
install_source()

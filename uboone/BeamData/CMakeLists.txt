find_ups_product( uboonedaq_datatypes v6_13_07 )
find_library( UB_DATA_TYPES NAMES ubdata_types PATHS $ENV{UBOONEDAQ_DATATYPES_LIB} )

art_make( BASENAME_ONLY
          MODULE_LIBRARIES
                        Utilities
                        TimeService_service
                        Filters
                        RawData
			${UB_DATA_TYPES}
			${Boost_SERIALIZATION_LIBRARY}
                        ${SEEDSERVICE_SERVICE} # artextensions
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
                        ${ROOT_BASIC_LIB_LIST}
        )

# install_headers()
install_fhicl()
install_source()


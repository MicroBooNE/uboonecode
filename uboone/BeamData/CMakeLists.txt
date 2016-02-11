find_ups_product( uboonedaq_datatypes v6_13_07 )
find_ups_product( ifdh_art v1_12_06 )

include_directories ( $ENV{IFDHC_FQ_DIR}/inc )

find_library( UB_DATA_TYPES NAMES ubdata_types PATHS $ENV{UBOONEDAQ_DATATYPES_LIB} )

art_make( BASENAME_ONLY
          MODULE_LIBRARIES
                        lardata_Utilities
                        lardata_Utilities_TimeService_service
                        larevt_Filters
                        lardata_RawData
			larcore_SummaryData
			BeamDAQ
			${UB_DATA_TYPES}
			${Boost_SERIALIZATION_LIBRARY}
                        ${SEEDSERVICE_SERVICE} # artextensions
                        ${ART_FRAMEWORK_CORE}
                        ${ART_FRAMEWORK_PRINCIPAL}
                        ${ART_FRAMEWORK_SERVICES_REGISTRY}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
                        ${ART_FRAMEWORK_SERVICES_SYSTEM_FILECATALOGMETADATA_SERVICE}
                        ${ART_PERSISTENCY_COMMON}
                        ${ART_PERSISTENCY_PROVENANCE}
                        ${ART_UTILITIES}
			ifdh
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


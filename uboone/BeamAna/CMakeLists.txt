art_make( BASENAME_ONLY
          MODULE_LIBRARIES
                        lardata_Utilities
                        larevt_Filters
                        lardata_RawData
			larcore_SummaryData
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

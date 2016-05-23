find_ups_product( uboonedaq_datatypes v6_13_07 )
find_ups_product( ifdh_art v1_12_06 )

include_directories ( $ENV{IFDHC_FQ_DIR}/inc )

find_library( UB_DATA_TYPES NAMES ubdata_types PATHS $ENV{UBOONEDAQ_DATATYPES_LIB} )
cet_make(LIBRARY_NAME BeamUtil LIBRARIES ${UB_DATA_TYPES} ${ROOT_SPECTRUM} ${ROOT_BASIC_LIB_LIST} EXCLUDE BeamData_module.cc BeamDataQualityFilter_module.cc)

simple_plugin ( BeamDataQualityFilter "module"
                   RawData
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
		   BASENAME_ONLY
              )

art_make( EXCLUDE BeamDataQualityFilter_module.cc
          BASENAME_ONLY
          MODULE_LIBRARIES
                        Utilities
                        TimeService_service
                        Filters
                        RawData
			SummaryData
			BeamUtil
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
		LIB_LIBRARIES
			${UB_DATA_TYPES}
			${ROOT_SPECTRUM}
                        ${ROOT_BASIC_LIB_LIST}
        )


# install_headers()
install_fhicl()
install_source()


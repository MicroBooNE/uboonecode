
include_directories ( $ENV{IFDHC_FQ_DIR}/inc )

find_library( UB_DATA_TYPES NAMES ubdata_types PATHS $ENV{UBOONEDAQ_DATATYPES_LIB} )
cet_make(LIBRARY_NAME BeamUtil LIBRARIES ${UB_DATA_TYPES} ${ROOT_SPECTRUM} ${ROOT_BASIC_LIB_LIST} EXCLUDE BeamData_module.cc BeamDataQualityFilter_module.cc)

simple_plugin ( BeamDataQualityFilter "module"
                   lardataobj_RawData
                   uboone_RawData
		   BeamUtil
                   ${ART_FRAMEWORK_CORE}
                   ${ART_FRAMEWORK_PRINCIPAL}
                   ${ART_FRAMEWORK_SERVICES_REGISTRY}
                   ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                   art_Persistency_Common canvas_Persistency_Common
                   art_Persistency_Provenance canvas_Persistency_Provenance
                   art_Utilities canvas_Utilities
                   ${MF_MESSAGELOGGER}
                   ${MF_UTILITIES}
		   BASENAME_ONLY
              )

art_make( EXCLUDE BeamDataQualityFilter_module.cc
          BASENAME_ONLY
          MODULE_LIBRARIES
                        lardata_Utilities
                        larevt_Filters
                        lardataobj_RawData
			larcoreobj_SummaryData
			BeamUtil
			BeamDAQ
			${UB_DATA_TYPES}
			${Boost_SERIALIZATION_LIBRARY}
                        nutools_RandomUtils_NuRandomService_service
                        ${ART_FRAMEWORK_CORE}
                        ${ART_FRAMEWORK_PRINCIPAL}
                        ${ART_FRAMEWORK_SERVICES_REGISTRY}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
                        ${ART_FRAMEWORK_SERVICES_SYSTEM_FILECATALOGMETADATA_SERVICE}
                        art_Persistency_Common canvas_Persistency_Common
                        art_Persistency_Provenance canvas_Persistency_Provenance
                        art_Utilities canvas_Utilities
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

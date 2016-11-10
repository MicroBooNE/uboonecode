# basic source code CMakeLists.txt

art_make(
	MODULE_LIBRARIES
			uboone_PatternFilter_PMAlgs
			larcore_Geometry
			lardata_DetectorInfo
			lardataobj_RecoBase
                        ${ART_FRAMEWORK_CORE}
			${ART_FRAMEWORK_IO_ROOT}
			${ART_FRAMEWORK_MODULES}
                        ${ART_FRAMEWORK_PRINCIPAL}
                        ${ART_FRAMEWORK_SERVICES_REGISTRY}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
                        art_Persistency_Common canvas_Persistency_Common
                        art_Utilities canvas_Utilities
			${MF_MESSAGELOGGER}
                        ${BOOST_LIB}
                        ${FHICLCPP}
                        ${CETLIB}
			${ROOT_BASIC_LIB_LIST}
)							   

install_headers()
install_source()
install_fhicl()

# Libraries I didn't need
#			ub_RawDigitAdderAna
#
#			${ART_FRAMEWORK_IO_SOURCES}
#                        art_Persistency_Provenance canvas_Persistency_Provenance
#                        ${MF_MESSAGELOGGER}
#                        ${MF_UTILITIES}
#                        ${ROOT_BASIC_LIB_LIST}

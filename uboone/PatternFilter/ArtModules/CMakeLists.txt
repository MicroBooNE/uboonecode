# basic source code CMakeLists.txt

art_make(
	LIBRARY_NAME AnodeCathodeTrackFilter_module
	MODULE_LIBRARIES
#			ub_RawDigitMixer
#			ub_OpDetWaveformMixer
#                        lardataobj_RawData
#			nusimdata_SimulationBase
#			lardataobj_MCBase
#			larevt_CalibrationDBI_IOVData
#			larevt_CalibrationDBI_Providers
                        ${ART_FRAMEWORK_CORE}
#			${ART_FRAMEWORK_IO_PRODUCTMIX}
			${ART_FRAMEWORK_IO_ROOT}
			${ART_FRAMEWORK_MODULES}
                        ${ART_FRAMEWORK_PRINCIPAL}
                        ${ART_FRAMEWORK_SERVICES_REGISTRY}
                        art_Persistency_Common canvas_Persistency_Common
                        art_Utilities canvas_Utilities
                        ${BOOST_LIB}
                        ${FHICLCPP}
                        ${CETLIB}
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

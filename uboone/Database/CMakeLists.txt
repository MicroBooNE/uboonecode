art_make( BASENAME_ONLY
	  SERVICE_LIBRARIES larevt_CalibrationDBI_Providers
	                    larevt_CalibrationDBI_IOVData
			    lardata_RawData
			    larcore_Geometry
			    larcore_Geometry_Geometry_service
			    ${ART_FRAMEWORK_SERVICES_REGISTRY}
			    ${ART_FRAMEWORK_CORE}
                            ${ART_FRAMEWORK_PRINCIPAL}
                            ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                            ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
                            ${ART_PERSISTENCY_COMMON}
                            ${ART_PERSISTENCY_PROVENANCE}
			    ${ART_UTILITIES}
			    ${FHICLCPP}
			    ${ROOT_BASIC_LIB_LIST}
)

install_headers()
install_fhicl()
install_source()


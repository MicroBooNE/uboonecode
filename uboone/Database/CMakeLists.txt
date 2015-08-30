art_make( BASENAME_ONLY
	  SERVICE_LIBRARIES CalibrationDBI_Providers
	                    CalibrationDBI_IOVData
			    RawData
			    Geometry
			    Geometry_service
			    ${ART_FRAMEWORK_SERVICES_REGISTRY}
			    ${ART_FRAMEWORK_CORE}
                            ${ART_FRAMEWORK_PRINCIPAL}
                            ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                            ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
                            ${ART_PERSISTENCY_COMMON}
                            ${ART_PERSISTENCY_PROVENANCE}
			    ${ART_UTILITIES}
			    ${FHICLCPP}
)

install_headers()
install_fhicl()
install_source()


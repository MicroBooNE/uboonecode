art_make( BASENAME_ONLY
	  SERVICE_LIBRARIES CalibrationDBI_Providers
                      ${ART_FRAMEWORK_PRINCIPAL}
                      ${ART_FRAMEWORK_SERVICES_REGISTRY}
                      ${ART_PERSISTENCY_COMMON}
                      ${ART_PERSISTENCY_PROVENANCE}
                      ${ART_FRAMEWORK_SERVICES_BASIC}
                      ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                      ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
                      ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
                      ${ART_FRAMEWORK_SERVICES_SYSTEM_FILECATALOGMETADATA_SERVICE}
                      ${ART_UTILITIES}
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
			    ${ROOT_BASIC_LIB_LIST}
)

install_headers()
install_fhicl()
install_source()


art_make( BASENAME_ONLY
	  SERVICE_LIBRARIES larevt_CalibrationDBI_Providers
                      ${ART_FRAMEWORK_PRINCIPAL}
                      ${ART_FRAMEWORK_SERVICES_REGISTRY}
                      art_Persistency_Common canvas_Persistency_Common
                      art_Persistency_Provenance canvas_Persistency_Provenance
                      ${ART_FRAMEWORK_SERVICES_BASIC}
                      ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                      ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
                      ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
                      ${ART_FRAMEWORK_SERVICES_SYSTEM_FILECATALOGMETADATA_SERVICE}
                      art_Utilities canvas_Utilities
	                    larevt_CalibrationDBI_IOVData
			    lardataobj_RawData
			    larcore_Geometry
			    larcore_Geometry_Geometry_service
			    ${ART_FRAMEWORK_SERVICES_REGISTRY}
			    ${ART_FRAMEWORK_CORE}
                            ${ART_FRAMEWORK_PRINCIPAL}
                            ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                            ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
                            art_Persistency_Common canvas_Persistency_Common
                            art_Persistency_Provenance canvas_Persistency_Provenance
			    art_Utilities canvas_Utilities
			    ${FHICLCPP}
			    ${ROOT_BASIC_LIB_LIST}
)

install_headers()
install_fhicl()
install_source()

art_make( BASENAME_ONLY
	  SERVICE_LIBRARIES CalibrationDBI_Providers
	                    CalibrationDBI_IOVData
			    ${ART_FRAMEWORK_SERVICES_REGISTRY}
			    ${FHICLCPP}
)

install_headers()
install_fhicl()
install_source()



art_make_library( LIBRARY_NAME uboone_cfdiscriminator
		  SOURCE CFDFire.cc CFDiscConfig.cc
		  LIBRARIES ${ART_FRAMEWORK_CORE}
                  	    ${ROOT_BASIC_LIB_LIST}
		)							   

install_headers()
install_source()

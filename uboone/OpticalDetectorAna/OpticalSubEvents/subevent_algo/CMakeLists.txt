
art_make_library( LIBRARY_NAME uboone_subeventalgo
		  SOURCE SubEvent.cc SubEventList.cc Flash.cc FlashList.cc SubEventModule.cc SubEventModConfig.cc
		         scintresponse.cc WaveformData.cc
		  LIBRARIES uboone_cfdiscriminator
		            ${ART_FRAMEWORK_CORE}
                  	    ${ROOT_BASIC_LIB_LIST}
		)							   

art_make( LIBRARY_NAME uboone_subeventalgodict
	  LIB_LIBRARIES uboone_subeventalgo
	  )

install_headers()
install_source()

# basic source code CMakeLists.txt

art_make_library( LIBRARY_NAME  ub_RawDigitAdderAna
	SOURCE RawDigitAdderAna.cxx
	LIBRARIES
                 lardataobj_RawData
                 ${ROOT_BASIC_LIB_LIST}
)

art_make_library( LIBRARY_NAME  ub_RawDigitAdder
	SOURCE RawDigitAdder.cxx
	             RawDigitAdder_HardSaturate.cxx
	LIBRARIES
                 ${ROOT_BASIC_LIB_LIST}
)

art_make_library( LIBRARY_NAME  ub_RawDigitMixer
	SOURCE RawDigitMixer.cxx
	LIBRARIES
		ub_RawDigitAdder
		lardataobj_RawData
		${ROOT_BASIC_LIB_LIST}
)

art_make_library( LIBRARY_NAME  ub_OpDetWaveformMixer
	SOURCE OpDetWaveformMixer.cxx
	LIBRARIES
		ub_RawDigitAdder
		lardataobj_RawData
		${ROOT_BASIC_LIB_LIST}
)

install_headers()
install_source()

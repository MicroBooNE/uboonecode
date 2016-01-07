# project name
project(datatypes)

# ----------------------------------------------------------------------
# Build this project's library:

set( LIB_NAME daq_datatypes)

set(LIB_SOURCES
evttypes.cpp
globalHeader.cpp
gps.cpp
triggerData.cpp
crateData.cpp
crateHeader.cpp
cardData.cpp
cardHeader.cpp
channelData.cpp
crateDataPMT.cpp
cardDataPMT.cpp
cardHeaderPMT.cpp
windowHeaderPMT.cpp
eventRecord.cpp
beamData.cpp
beamHeader.cpp
)

include_directories ( $ENV{WINDRIVER_INC} )
set(WINDRIVERLIB "-L$ENV{WINDRIVER_LIB} -lwdapi1100 -lnevispcie")
# rs:  removed 3 lines:
# find_ups_product( symmetricom        v2_00_00 )
# include_directories ( $ENV{SYMMETRICOM_INC} )
# set(SYMMETRICOMLIB "-L$ENV{SYMMETRICOM_LIB} -lbcsdk")

# The inclusion of gps/symm.h, necessary for stuffing the GPS time into the crate_Header, drags with it complaints
# about windriver.h functions that don't compile cleanly under gcc10+. We must add all these lines to
# get stuff to compile therefore.

set_source_files_properties(gps.cpp PROPERTIES COMPILE_FLAGS "-Wno-unused-variable -Wno-deprecated-declarations -Wno-unused-function -Wno-sign-compare -Wno-unused-parameter -Wno-format -Wno-format-y2k -Wno-unused-but-set-variable -Wno-return-type")
set_source_files_properties(crateHeader.cpp PROPERTIES COMPILE_FLAGS "-Wno-unused-variable -Wno-deprecated-declarations -Wno-unused-function -Wno-sign-compare -Wno-unused-parameter -Wno-format -Wno-format-y2k -Wno-unused-but-set-variable -Wno-return-type")
set_source_files_properties(crateData.cpp PROPERTIES COMPILE_FLAGS "-Wno-unused-variable -Wno-deprecated-declarations -Wno-unused-function -Wno-sign-compare -Wno-unused-parameter -Wno-format -Wno-format-y2k -Wno-unused-but-set-variable -Wno-return-type")
set_source_files_properties(cardHeader.cpp PROPERTIES COMPILE_FLAGS "-Wno-unused-variable -Wno-deprecated-declarations -Wno-unused-function -Wno-sign-compare -Wno-unused-parameter -Wno-format -Wno-format-y2k -Wno-unused-but-set-variable -Wno-return-type")
set_source_files_properties(cardData.cpp PROPERTIES COMPILE_FLAGS "-Wno-unused-variable -Wno-deprecated-declarations -Wno-unused-function -Wno-sign-compare -Wno-unused-parameter -Wno-format -Wno-format-y2k -Wno-unused-but-set-variable -Wno-return-type")
set_source_files_properties(channelData.cpp PROPERTIES COMPILE_FLAGS "-Wno-unused-variable -Wno-deprecated-declarations -Wno-unused-function -Wno-sign-compare -Wno-unused-parameter -Wno-format -Wno-format-y2k -Wno-unused-but-set-variable -Wno-return-type")
set_source_files_properties(crateDataPMT.cpp PROPERTIES COMPILE_FLAGS "-Wno-unused-variable -Wno-deprecated-declarations -Wno-unused-function -Wno-sign-compare -Wno-unused-parameter -Wno-format -Wno-format-y2k -Wno-unused-but-set-variable -Wno-return-type")
set_source_files_properties(cardHeaderPMT.cpp PROPERTIES COMPILE_FLAGS "-Wno-unused-variable -Wno-deprecated-declarations -Wno-unused-function -Wno-sign-compare -Wno-unused-parameter -Wno-format -Wno-format-y2k -Wno-unused-but-set-variable -Wno-return-type")
set_source_files_properties(cardDataPMT.cpp PROPERTIES COMPILE_FLAGS "-Wno-unused-variable -Wno-deprecated-declarations -Wno-unused-function -Wno-sign-compare -Wno-unused-parameter -Wno-format -Wno-format-y2k -Wno-unused-but-set-variable -Wno-return-type")
set_source_files_properties(windowHeaderPMT.cpp PROPERTIES COMPILE_FLAGS "-Wno-unused-variable -Wno-deprecated-declarations -Wno-unused-function -Wno-sign-compare -Wno-unused-parameter -Wno-format -Wno-format-y2k -Wno-unused-but-set-variable -Wno-return-type")
set_source_files_properties(triggerHeader.cpp PROPERTIES COMPILE_FLAGS "-Wno-unused-variable -Wno-deprecated-declarations -Wno-unused-function -Wno-sign-compare -Wno-unused-parameter -Wno-format -Wno-format-y2k -Wno-unused-but-set-variable -Wno-return-type")
set_source_files_properties(triggerData.cpp PROPERTIES COMPILE_FLAGS "-Wno-unused-variable -Wno-deprecated-declarations -Wno-unused-function -Wno-sign-compare -Wno-unused-parameter -Wno-format -Wno-format-y2k -Wno-unused-but-set-variable -Wno-return-type")

set_source_files_properties(eventRecord.cpp PROPERTIES COMPILE_FLAGS "-Wno-unused-variable -Wno-deprecated-declarations -Wno-unused-function -Wno-sign-compare -Wno-unused-parameter -Wno-format -Wno-format-y2k -Wno-unused-but-set-variable -Wno-return-type")

art_make_library(LIBRARY_NAME ${LIB_NAME} SOURCE ${LIB_SOURCES}
LIBRARIES # daq_share
${Boost_SERIALIZATION_LIBRARY}
${Boost_DATE_TIME_LIBRARY} 
${CMAKE_DL_LIBS})


#add_executable(write_read write_read.cpp)

#target_link_libraries(write_read
#    daq_share 
#    ${LIB_NAME}
#    ${COMMON_PROJECT_LIBS}
#    ${CMAKE_DL_LIBS} )

# ----------------------------------------------------------------------
# Install this project's headers, sources, library, and executables:

install_headers()
install_source()

install( TARGETS ${LIB_NAME}
         DESTINATION ${flavorqual_dir}/lib )


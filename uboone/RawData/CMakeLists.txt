# copied version of uboone_datatypes to build here. never used?
#add_subdirectory(uboone_datatypes)

# Include for Swizzler
set( UBOONEDAQ_DATATYPES_DIR $ENV{UBOONEDAQ_DATATYPES_DIR} )
if( UBOONEDAQ_DATATYPES_DIR )
   add_subdirectory(utils)
#   add_subdirectory(uboone_datatypes)
endif()

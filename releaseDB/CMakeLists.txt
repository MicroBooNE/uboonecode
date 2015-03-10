CONFIGURE_FILE ( ${CMAKE_CURRENT_SOURCE_DIR}/uboone.html.in
  ${CMAKE_CURRENT_BINARY_DIR}/uboone.html @ONLY )

INSTALL( FILES base_dependency_database
               README
         DESTINATION "${product}/${version}/releaseDB" )

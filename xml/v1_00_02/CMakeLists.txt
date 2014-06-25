FILE(GLOB xml_files *.xml)
install(FILES ${xml_files} DESTINATION ${product}/${version}/xml/v1_00_02)
install_source( EXTRAS ${xml_files} )

FILE(GLOB xml_files *.xml)
install(FILES ${xml_files} DESTINATION ${product}/${version}/xml/S2012.12.17)
install_source( EXTRAS ${xml_files} )

FILE(GLOB xml_files *.xml)
install(FILES ${xml_files} DESTINATION ${product}/${version}/xml/S2013.06.25)
install_source( EXTRAS ${xml_files} )

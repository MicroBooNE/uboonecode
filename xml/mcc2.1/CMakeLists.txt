FILE(GLOB xml_files *.xml)
install(FILES ${xml_files} DESTINATION ${product}/${version}/xml/mcc2.1)
